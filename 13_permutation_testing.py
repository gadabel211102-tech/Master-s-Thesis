#!/usr/bin/env python3
"""
Permutation testing for SNP carrier frequency differences between tumour and
healthy tissue samples. Complements the Fisher's exact test used in script 12
by providing a non-parametric alternative that makes no distributional
assumptions — important here because sample sizes are small and unequal across
cohorts.

SNP definition: variants with gnomAD NFE allele frequency > 1%, consistent
with scripts 11 and 12. Exome NFE AF is preferred over genome NFE AF where
both are available.

For each SNP and cohort, the script:
  1. Runs a permutation test (10,000 shuffles by default) and computes Cohen's h
  2. Runs Fisher's exact test on the same 2x2 carrier/non-carrier table
  3. Applies Benjamini-Hochberg FDR correction within each cohort independently
  4. Flags agreement and disagreement between the two test methods

Outputs: an Excel workbook with per-cohort result sheets and a summary, plus a
four-panel comparison figure.
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tqdm import tqdm
import warnings
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# --- Configuration -----------------------------------------------------------
BASE_PATH  = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
INPUT_FILE = os.path.join(BASE_PATH, "GSDMB_Annotated_Report_Fixed.xlsx")
OUTPUT_XLSX = os.path.join(BASE_PATH, "13_Permutation_Test_Results.xlsx")
OUTPUT_PLOT = os.path.join(BASE_PATH, "13_Permutation_Comparison.png")

N_PERMUTATIONS = 10000   # number of label shuffles per SNP; 10,000 gives resolution to p = 0.0001
RANDOM_SEED    = 42      # fixed seed for reproducibility across runs


def find_col(df: pd.DataFrame, target: str):
    """Case-insensitive column lookup — VEP field names vary in capitalisation between cache versions."""
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None


def permutation_test_snp(
    tumour_carriers:  int,
    tumour_total:     int,
    healthy_carriers: int,
    healthy_total:    int,
    n_permutations:   int = 10000,
    random_state:     int = 42
) -> dict:
    """
    Non-parametric permutation test for a difference in SNP carrier proportions
    between tumour and healthy groups.

    The null hypothesis is that carrier status is independent of tissue type.
    Under the null, randomly reassigning the 'tumour' or 'healthy' label to each
    individual should produce a distribution of proportion differences centred on
    zero. The observed difference is compared against this null distribution to
    obtain a two-tailed p-value.

    Procedure:
      1. Represent each individual as 1 (carrier) or 0 (non-carrier).
      2. Pool tumour and healthy arrays and shuffle N_PERMUTATIONS times.
      3. At each shuffle, split the pool at the original tumour group size and
         compute the difference in means (= difference in carrier proportions).
      4. p-value = fraction of permuted differences at least as extreme as observed.

    Effect size is reported as Cohen's h, the standard measure for comparing two
    proportions. h = 0.2, 0.5, 0.8 correspond to small, medium, and large effects.
    The arcsine transformation stabilises variance across the [0, 1] range.

    Returns a dict with keys: observed_diff, p_value, effect_size_h, perm_diffs,
    tumour_prop, healthy_prop.
    """
    np.random.seed(random_state)

    tumour_data  = np.array([1] * tumour_carriers  + [0] * (tumour_total  - tumour_carriers))
    healthy_data = np.array([1] * healthy_carriers + [0] * (healthy_total - healthy_carriers))

    obs_tumour_prop  = tumour_carriers  / tumour_total  if tumour_total  > 0 else 0
    obs_healthy_prop = healthy_carriers / healthy_total if healthy_total > 0 else 0
    observed_diff    = obs_tumour_prop - obs_healthy_prop

    # Pool all individuals and repeatedly draw a random partition of the same
    # sizes as the original groups; the mean of each partition equals its carrier proportion
    all_data = np.concatenate([tumour_data, healthy_data])
    n_tumour = len(tumour_data)

    perm_diffs = np.empty(n_permutations)
    for i in range(n_permutations):
        shuffled           = np.random.permutation(all_data)
        perm_diffs[i]      = np.mean(shuffled[:n_tumour]) - np.mean(shuffled[n_tumour:])

    # Two-tailed: count permuted differences at least as large in magnitude as observed
    p_value = np.mean(np.abs(perm_diffs) >= np.abs(observed_diff))

    # Cohen's h: arcsine-transformed difference in proportions
    h = 2 * (np.arcsin(np.sqrt(obs_tumour_prop)) - np.arcsin(np.sqrt(obs_healthy_prop)))

    return {
        "observed_diff": observed_diff,
        "p_value":       p_value,
        "effect_size_h": h,
        "perm_diffs":    perm_diffs,
        "tumour_prop":   obs_tumour_prop,
        "healthy_prop":  obs_healthy_prop
    }


def run_permutation_analysis():
    """
    Load the annotated variant report, filter for common SNPs, run permutation
    and Fisher tests for each SNP across three cohort groupings (Global, Breast,
    Endometrium), apply FDR correction, and export results.
    """
    print("=" * 70)
    print("SCRIPT 13: PERMUTATION TESTING FOR SNP ASSOCIATIONS")
    print("=" * 70)
    print(f"\n  Input:        {INPUT_FILE}")
    print(f"  Permutations: {N_PERMUTATIONS:,}")
    print(f"  Random seed:  {RANDOM_SEED}")
    print("  SNP filter:   gnomAD NFE AF > 1% (exome preferred, genome as fallback)")
    print("\n" + "=" * 70 + "\n")

    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: input file not found: {INPUT_FILE}")
        return

    df = pd.read_excel(INPUT_FILE, sheet_name="Biological_Annotations")

    # Locate required columns by case-insensitive name
    sym_c        = find_col(df, "SYMBOL")
    tis_c        = find_col(df, "TISSUE")
    sam_c        = find_col(df, "SAMPLE")
    var_c        = find_col(df, "Existing_variation")
    hgv_c        = find_col(df, "HGVSp")
    coh_c        = find_col(df, "COHORT")
    exome_nfe_c  = find_col(df, "gnomADe_NFE_AF")
    genome_nfe_c = find_col(df, "gnomADg_NFE_AF")
    imp_c        = find_col(df, "IMPACT")
    con_c        = find_col(df, "CONSEQUENCE")

    required = [sym_c, tis_c, sam_c, var_c, hgv_c, coh_c,
                exome_nfe_c, genome_nfe_c, imp_c, con_c]
    if any(c is None for c in required):
        print("ERROR: one or more required columns are missing from the input file.")
        print(f"  SYMBOL={sym_c}, TISSUE={tis_c}, SAMPLE={sam_c}, "
              f"Existing_variation={var_c}, HGVSp={hgv_c}, COHORT={coh_c}, "
              f"gnomADe_NFE_AF={exome_nfe_c}, gnomADg_NFE_AF={genome_nfe_c}, "
              f"IMPACT={imp_c}, CONSEQUENCE={con_c}")
        return

    # --- Build combined NFE AF ---------------------------------------------------
    # Exome and genome gnomAD NFE frequencies are stored in separate columns;
    # combine_first fills missing exome values with genome values so that every
    # variant has at most one AF estimate. Tracking the source column allows
    # downstream inspection of which database was used for each variant.
    df["gnomAD_NFE_AF_combined"] = df[exome_nfe_c].combine_first(df[genome_nfe_c])
    df["gnomAD_NFE_Source"] = np.where(
        df[exome_nfe_c].notna(), "Exome_NFE",
        np.where(df[genome_nfe_c].notna(), "Genome_NFE", "Missing")
    )

    # Retain only common SNPs: NFE AF > 1% excludes rare variants whose carrier
    # counts would be too small for meaningful frequency comparison
    df = df[df["gnomAD_NFE_AF_combined"] > 0.01].copy()
    if df.empty:
        print("No variants passed the gnomAD NFE AF > 1% filter.")
        return

    # Standardise tissue labels to Tumour / Healthy
    df[tis_c] = (df[tis_c].astype(str).str.strip()
                 .replace({"Tumor": "Tumour", "Tumour": "Tumour",
                           "Normal": "Healthy", "Healthy": "Healthy", "Control": "Healthy"}))
    df[coh_c] = df[coh_c].astype(str).str.strip()
    df = df[df[tis_c].isin(["Tumour", "Healthy"])].copy()

    # Build a stable variant identifier: rsID preferred; falls back to SYMBOL:HGVSp
    df["rsID"]       = df[var_c].astype(str).str.extract(r"(rs\d+)")
    df["Variant_ID"] = df["rsID"].fillna(
        df[sym_c].astype(str) + ":" + df[hgv_c].astype(str)
    )

    cohort_list = ["Global", "Breast", "Endometrium"]
    all_results = []

    for cohort_name in cohort_list:
        print(f"\n{'─' * 70}")
        print(f"Cohort: {cohort_name}")
        print(f"{'─' * 70}")

        cohort_df = df.copy() if cohort_name == "Global" else \
                    df[df[coh_c].str.contains(cohort_name, case=False, na=False)].copy()

        if cohort_df.empty:
            print(f"  No data for {cohort_name}"); continue

        # Sample counts: unique sample IDs within each tissue group
        n_tumour  = cohort_df[cohort_df[tis_c] == "Tumour"][sam_c].nunique()
        n_healthy = cohort_df[cohort_df[tis_c] == "Healthy"][sam_c].nunique()
        if n_tumour == 0 or n_healthy == 0:
            print(f"  Insufficient samples: Tumour={n_tumour}, Healthy={n_healthy}"); continue

        print(f"  Samples — Tumour: {n_tumour}, Healthy: {n_healthy}")
        unique_snps = cohort_df["Variant_ID"].dropna().unique()
        print(f"  Unique SNPs to test: {len(unique_snps)}")

        for var in tqdm(unique_snps, desc=f"  {cohort_name}"):
            var_data = cohort_df[cohort_df["Variant_ID"] == var]

            # Carrier count: number of distinct samples that carry this variant
            count_tumour  = var_data[var_data[tis_c] == "Tumour"][sam_c].nunique()
            count_healthy = var_data[var_data[tis_c] == "Healthy"][sam_c].nunique()

            perm_result = permutation_test_snp(
                tumour_carriers  = count_tumour,
                tumour_total     = n_tumour,
                healthy_carriers = count_healthy,
                healthy_total    = n_healthy,
                n_permutations   = N_PERMUTATIONS,
                random_state     = RANDOM_SEED
            )

            # Fisher's exact test on the same 2x2 contingency table:
            #   [ carriers_tumour,     non-carriers_tumour  ]
            #   [ carriers_healthy,    non-carriers_healthy ]
            a, b = count_tumour,  n_tumour  - count_tumour
            c, d = count_healthy, n_healthy - count_healthy
            fisher_or, fisher_p = stats.fisher_exact([[a, b], [c, d]])

            all_results.append({
                "Cohort":            cohort_name,
                "Variant_ID":        var,
                "Symbol":            var_data[sym_c].iloc[0],
                "Impact":            var_data[imp_c].iloc[0],
                "Consequence":       var_data[con_c].iloc[0],
                "gnomAD_NFE_AF":     var_data["gnomAD_NFE_AF_combined"].iloc[0],
                "gnomAD_NFE_Source": var_data["gnomAD_NFE_Source"].iloc[0],
                "Tumour_Carriers":   count_tumour,
                "Tumour_Total":      n_tumour,
                "Tumour_Freq_%":     perm_result["tumour_prop"] * 100,
                "Healthy_Carriers":  count_healthy,
                "Healthy_Total":     n_healthy,
                "Healthy_Freq_%":    perm_result["healthy_prop"] * 100,
                "Freq_Difference_%": perm_result["observed_diff"] * 100,
                "Permutation_P_Value": perm_result["p_value"],
                "Effect_Size_h":     perm_result["effect_size_h"],
                "Fisher_P_Value":    fisher_p,
                "Fisher_OR":         fisher_or,
                "N_Permutations":    N_PERMUTATIONS
            })

    if not all_results:
        print("No permutation results were generated.")
        return

    results_df = pd.DataFrame(all_results)

    # --- Multiple testing correction ----------------------------------------
    # FDR correction is applied within each cohort separately rather than across
    # all cohorts combined, because the tests within a cohort share a common
    # null population whereas cross-cohort comparisons are independent analyses
    corrected_results = []
    for cohort in results_df["Cohort"].unique():
        cohort_data = results_df[results_df["Cohort"] == cohort].copy()
        _, cohort_data["Permutation_FDR"], _, _ = multipletests(
            cohort_data["Permutation_P_Value"], method="fdr_bh"
        )
        _, cohort_data["Fisher_FDR"], _, _ = multipletests(
            cohort_data["Fisher_P_Value"], method="fdr_bh"
        )
        corrected_results.append(cohort_data)

    final_results = (pd.concat(corrected_results, ignore_index=True)
                     .sort_values(["Cohort", "Permutation_P_Value"]))

    # Significance flags at both uncorrected (raw) and FDR-adjusted thresholds
    final_results["Permutation_Significant_RAW"] = final_results["Permutation_P_Value"] < 0.05
    final_results["Fisher_Significant_RAW"]      = final_results["Fisher_P_Value"]      < 0.05
    final_results["Permutation_Significant_FDR"] = final_results["Permutation_FDR"]     < 0.05
    final_results["Fisher_Significant_FDR"]      = final_results["Fisher_FDR"]          < 0.05

    # Agreement columns indicate whether both methods reach the same significance
    # conclusion — high agreement validates the Fisher result for this sample size
    final_results["Agreement_RAW"] = (final_results["Permutation_Significant_RAW"] ==
                                      final_results["Fisher_Significant_RAW"])
    final_results["Agreement_FDR"] = (final_results["Permutation_Significant_FDR"] ==
                                      final_results["Fisher_Significant_FDR"])

    # --- Save results --------------------------------------------------------
    print(f"\n{'=' * 70}\nSaving results\n{'=' * 70}\n")

    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
        # Summary sheet: one row per cohort with counts of significant SNPs
        summary = final_results.groupby("Cohort").agg({
            "Variant_ID":                "count",
            "Permutation_Significant_RAW": "sum",
            "Fisher_Significant_RAW":      "sum",
            "Permutation_Significant_FDR": "sum",
            "Fisher_Significant_FDR":      "sum",
            "Agreement_RAW":               "sum",
            "Agreement_FDR":               "sum"
        }).reset_index()
        summary.columns = ["Cohort", "Total_SNPs",
                           "Perm_Significant_RAW", "Fisher_Significant_RAW",
                           "Perm_Significant_FDR", "Fisher_Significant_FDR",
                           "Test_Agreement_RAW",   "Test_Agreement_FDR"]
        summary.to_excel(writer, sheet_name="Summary", index=False)

        for cohort in cohort_list:
            cohort_data = final_results[final_results["Cohort"] == cohort]
            if not cohort_data.empty:
                cohort_data.to_excel(writer, sheet_name=cohort, index=False)

    print(f"Results saved: {OUTPUT_XLSX}\n")
    generate_comparison_plots(final_results)

    # --- Print summary statistics --------------------------------------------
    print(f"{'=' * 70}\nSummary statistics\n{'=' * 70}\n")
    for cohort in cohort_list:
        cohort_data = final_results[final_results["Cohort"] == cohort]
        if cohort_data.empty:
            continue
        print(f"{cohort} cohort:")
        print(f"  SNPs tested:                  {len(cohort_data)}")
        print(f"  Permutation sig. (raw p<0.05): {cohort_data['Permutation_Significant_RAW'].sum()}")
        print(f"  Fisher sig.      (raw p<0.05): {cohort_data['Fisher_Significant_RAW'].sum()}")
        print(f"  Permutation sig. (FDR<0.05):   {cohort_data['Permutation_Significant_FDR'].sum()}")
        print(f"  Fisher sig.      (FDR<0.05):   {cohort_data['Fisher_Significant_FDR'].sum()}")
        print(f"  Agreement (raw): {cohort_data['Agreement_RAW'].sum()}/{len(cohort_data)} "
              f"({cohort_data['Agreement_RAW'].mean()*100:.1f}%)")
        print(f"  Agreement (FDR): {cohort_data['Agreement_FDR'].sum()}/{len(cohort_data)} "
              f"({cohort_data['Agreement_FDR'].mean()*100:.1f}%)")

        top5 = cohort_data.nsmallest(5, "Permutation_P_Value")
        if not top5.empty:
            print("  Top 5 by permutation p-value:")
            for _, row in top5.iterrows():
                print(f"    {row['Variant_ID']}  "
                      f"p={row['Permutation_P_Value']:.4e}  "
                      f"FDR={row['Permutation_FDR']:.4e}  "
                      f"source={row['gnomAD_NFE_Source']}")
        print()

    print(f"{'=' * 70}\nPermutation analysis complete.\n{'=' * 70}\n")


def generate_comparison_plots(results_df: pd.DataFrame):
    """
    Four-panel figure comparing permutation and Fisher test results:

    Panel 1 (top-left):  scatter plot of -log10 raw p-values from both tests;
                         points on the diagonal indicate agreement, deviations
                         indicate cases where sample size affects one test more
    Panel 2 (top-right): Cohen's h vs permutation -log10 p; confirms that
                         significant results have meaningful effect sizes
    Panel 3 (bottom-left): histogram of carrier frequency differences
                         (tumour - healthy); centring on zero supports the null
    Panel 4 (bottom-right): stacked bar chart of FDR-based test agreement
                         categories per cohort
    """
    print("Generating comparison plots...")
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    cohorts = results_df["Cohort"].unique()
    colours = {"Global": "#3498db", "Breast": "#e74c3c", "Endometrium": "#2ecc71"}

    # Panel 1: raw p-value agreement scatter
    ax1 = axes[0, 0]
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        ax1.scatter(
            -np.log10(data["Fisher_P_Value"].replace(0, 1e-300)),
            -np.log10(data["Permutation_P_Value"].replace(0, 1e-300)),
            label=cohort, alpha=0.6, s=80, color=colours.get(cohort, "#95a5a6")
        )
    max_val = max(ax1.get_xlim()[1], ax1.get_ylim()[1])
    ax1.plot([0, max_val], [0, max_val], "k--", alpha=0.3, label="x=y")   # line of perfect agreement
    sig_line = -np.log10(0.05)
    ax1.axhline(sig_line, color="red", linestyle=":", alpha=0.5, label="p=0.05")
    ax1.axvline(sig_line, color="red", linestyle=":", alpha=0.5)
    ax1.set_xlabel("Fisher's Exact Test (-log10 raw p)", fontweight="bold")
    ax1.set_ylabel("Permutation Test (-log10 raw p)", fontweight="bold")
    ax1.set_title("Raw P-Value Comparison: Permutation vs Fisher", fontweight="bold", fontsize=12)
    ax1.legend(loc="best")
    ax1.grid(True, alpha=0.3)

    # Panel 2: effect size vs permutation significance
    ax2 = axes[0, 1]
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        ax2.scatter(
            data["Effect_Size_h"],
            -np.log10(data["Permutation_P_Value"].replace(0, 1e-300)),
            label=cohort, alpha=0.6, s=80, color=colours.get(cohort, "#95a5a6")
        )
    ax2.axhline(sig_line, color="red", linestyle=":", alpha=0.5, label="p=0.05")
    ax2.set_xlabel("Effect Size (Cohen's h)", fontweight="bold")
    ax2.set_ylabel("Permutation Test (-log10 raw p)", fontweight="bold")
    ax2.set_title("Effect Size vs Statistical Significance", fontweight="bold", fontsize=12)
    ax2.legend(loc="best")
    ax2.grid(True, alpha=0.3)

    # Panel 3: distribution of carrier frequency differences
    ax3 = axes[1, 0]
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        ax3.hist(data["Freq_Difference_%"], bins=30, alpha=0.5,
                 label=cohort, color=colours.get(cohort, "#95a5a6"), edgecolor="black")
    ax3.axvline(0, color="black", linestyle="--", alpha=0.5, label="No difference")
    ax3.set_xlabel("Frequency Difference (Tumour - Healthy %)", fontweight="bold")
    ax3.set_ylabel("Number of SNPs", fontweight="bold")
    ax3.set_title("Distribution of Carrier Frequency Differences", fontweight="bold", fontsize=12)
    ax3.legend(loc="best")
    ax3.grid(True, alpha=0.3, axis="y")

    # Panel 4: FDR-based agreement classification per cohort
    ax4 = axes[1, 1]
    agreement_data = []
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        agreement_data.append({
            "Cohort":           cohort,
            "Both Significant": ((data["Permutation_Significant_FDR"]) & (data["Fisher_Significant_FDR"])).sum(),
            "Perm Only":        ((data["Permutation_Significant_FDR"]) & (~data["Fisher_Significant_FDR"])).sum(),
            "Fisher Only":      ((~data["Permutation_Significant_FDR"]) & (data["Fisher_Significant_FDR"])).sum(),
            "Neither":          ((~data["Permutation_Significant_FDR"]) & (~data["Fisher_Significant_FDR"])).sum()
        })
    agreement_df = pd.DataFrame(agreement_data).set_index("Cohort")
    agreement_df.plot(kind="bar", stacked=True, ax=ax4,
                      color=["#2ecc71", "#f39c12", "#e74c3c", "#95a5a6"])
    ax4.set_xlabel("Cohort", fontweight="bold")
    ax4.set_ylabel("Number of SNPs", fontweight="bold")
    ax4.set_title("Test Agreement by FDR Significance", fontweight="bold", fontsize=12)
    ax4.legend(title="Significance category", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.setp(ax4.xaxis.get_majorticklabels(), rotation=0)
    ax4.grid(True, alpha=0.3, axis="y")

    plt.suptitle(
        f"Permutation Testing vs Fisher's Exact Test\n"
        f"({N_PERMUTATIONS:,} permutations per SNP, gnomAD NFE AF > 1%)",
        fontsize=16, fontweight="bold", y=0.995
    )
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Comparison plots saved: {OUTPUT_PLOT}\n")


if __name__ == "__main__":
    run_permutation_analysis()
