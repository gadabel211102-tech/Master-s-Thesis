#!/usr/bin/env python3
"""
Script 13: Permutation Testing for SNP Associations
===================================================

Purpose
-------
Re-evaluate tumour-versus-control SNP differences with a non-parametric
permutation framework that is less sensitive to small or imbalanced groups than
closed-form contingency testing alone.

Methodological rationale
------------------------
The script keeps the SNP definition aligned with scripts 11 and 12 by using
exome-or-genome gnomAD NFE AF > 1%. Comparison plots are expressed in
proportional terms wherever possible so they do not visually overstate patterns
from larger sample sets.

Date: 2026-02-12
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable
import warnings
from statsmodels.stats.multitest import multipletests

from pipeline_utils import (
    attach_amplicon_warning_columns,
    build_amplicon_warning_lookup,
    build_genomic_variant_id_series,
    combine_gnomad_nfe,
    common_nfe_variant_mask,
    ensure_directory,
    extract_rsid,
    find_col,
    get_paths,
    get_thresholds,
    pooled_control_frame,
    standardize_cohort_labels,
    standardize_tissue_labels,
    load_common_snp_ld_reference,
)
from pipeline_validation import print_validation_summary, validate_file_exists, validate_nonempty, validate_percentage_columns, validate_required_columns, validate_tissue_values

warnings.filterwarnings("ignore")

# --- CONFIGURATION ---
PATHS = get_paths()
THRESHOLDS = get_thresholds()
OUTPUT_DIR = ensure_directory(PATHS["permutation_testing_dir"])
INPUT_FILE = str(PATHS["annotated_report"])
OUTPUT_XLSX = str(OUTPUT_DIR / "GSDMB_SNP_Permutation_Test_Results.xlsx")
OUTPUT_PLOT = str(OUTPUT_DIR / "GSDMB_SNP_Permutation_Comparison.png")

# Permutation parameters
N_PERMUTATIONS = THRESHOLDS["n_permutations"]
RANDOM_SEED = THRESHOLDS["random_seed"]
LD_R2_THRESHOLD = 0.80


def permutation_test_snp(
    tumour_carriers,
    tumour_total,
    healthy_carriers,
    healthy_total,
    n_permutations=10000,
    random_state=42
):
    """
    Perform permutation test for SNP association.

    Null hypothesis:
        The SNP carrier status is independent of tissue type.
    """
    np.random.seed(random_state)

    # Binary arrays: 1 = carrier, 0 = non-carrier
    tumour_data = np.array([1] * tumour_carriers + [0] * (tumour_total - tumour_carriers))
    healthy_data = np.array([1] * healthy_carriers + [0] * (healthy_total - healthy_carriers))

    # Observed difference in proportions
    obs_tumour_prop = tumour_carriers / tumour_total if tumour_total > 0 else 0
    obs_healthy_prop = healthy_carriers / healthy_total if healthy_total > 0 else 0
    observed_diff = obs_tumour_prop - obs_healthy_prop

    # Combine and permute
    all_data = np.concatenate([tumour_data, healthy_data])
    n_tumour = len(tumour_data)

    perm_diffs = []
    for _ in range(n_permutations):
        shuffled = np.random.permutation(all_data)
        perm_tumour = shuffled[:n_tumour]
        perm_healthy = shuffled[n_tumour:]

        perm_tumour_prop = np.mean(perm_tumour)
        perm_healthy_prop = np.mean(perm_healthy)
        perm_diffs.append(perm_tumour_prop - perm_healthy_prop)

    perm_diffs = np.array(perm_diffs)

    # Two-tailed p-value
    p_value = np.mean(np.abs(perm_diffs) >= np.abs(observed_diff))

    # Cohen's h effect size for proportions
    h = 2 * (
        np.arcsin(np.sqrt(obs_tumour_prop)) -
        np.arcsin(np.sqrt(obs_healthy_prop))
    )

    return {
        "observed_diff": observed_diff,
        "p_value": p_value,
        "effect_size_h": h,
        "perm_diffs": perm_diffs,
        "tumour_prop": obs_tumour_prop,
        "healthy_prop": obs_healthy_prop
    }



def _singleton_block_id(variant_id):
    return f"Singleton::{variant_id}"


def _initialise_ld_block_columns(results_df, ld_reference):
    block_df = ld_reference["block_membership"].copy()
    if not block_df.empty:
        block_df["LD_Block_ID"] = block_df["Block"].astype(str)
        block_map = dict(zip(block_df["rsID"], block_df["LD_Block_ID"]))
        block_sizes = block_df.groupby("LD_Block_ID")["rsID"].size().to_dict()
        block_members = (
            block_df.groupby("LD_Block_ID")["rsID"]
            .apply(lambda series: ", ".join(series.astype(str)))
            .to_dict()
        )
    else:
        block_map = {}
        block_sizes = {}
        block_members = {}

    out = results_df.copy()
    out["LD_Threshold"] = ld_reference["threshold"]
    out["Threshold_Label"] = ld_reference["threshold_label"]
    out["LD_Block_ID"] = out["Variant_ID"].map(block_map)
    out["LD_Block_Assignment"] = np.where(out["LD_Block_ID"].notna(), "LD_Block", "Singleton")
    out.loc[out["LD_Block_ID"].isna(), "LD_Block_ID"] = out.loc[out["LD_Block_ID"].isna(), "Variant_ID"].map(_singleton_block_id)
    out["LD_Block_Size"] = out["LD_Block_ID"].map(block_sizes).fillna(1).astype(int)
    out["LD_Block_Members"] = out["LD_Block_ID"].map(block_members).fillna(out["Variant_ID"])
    return out, block_df


def _apply_block_fdr(results_df, cohort_col, block_col, specs):
    final_parts = []
    for cohort_value, cohort_df in results_df.groupby(cohort_col, sort=False):
        cohort_df = cohort_df.copy()
        for p_col, fdr_col, legacy_col in specs:
            _, cohort_df[legacy_col], _, _ = multipletests(cohort_df[p_col], method="fdr_bh")
        block_parts = []
        for _, block_df in cohort_df.groupby(block_col, sort=False):
            block_df = block_df.copy()
            block_df["LD_Block_Test_Count"] = int(len(block_df))
            for p_col, fdr_col, legacy_col in specs:
                if len(block_df) == 1:
                    block_df[fdr_col] = block_df[p_col]
                else:
                    _, block_df[fdr_col], _, _ = multipletests(block_df[p_col], method="fdr_bh")
            block_parts.append(block_df)
        cohort_df = pd.concat(block_parts, ignore_index=True)
        threshold_label = cohort_df["Threshold_Label"].iloc[0] if "Threshold_Label" in cohort_df.columns else "LD"
        cohort_df["FDR_Scope"] = f"Within_{cohort_value}_{threshold_label}_blocks"
        final_parts.append(cohort_df)
    return pd.concat(final_parts, ignore_index=True)

def run_permutation_analysis():
    """Main pipeline for permutation testing."""
    print("=" * 70)
    print("SCRIPT 13: PERMUTATION TESTING FOR SNP ASSOCIATIONS")
    print("=" * 70)
    print("\nConfiguration:")
    print(f"  - Input: {INPUT_FILE}")
    print(f"  - Permutations: {N_PERMUTATIONS:,}")
    print(f"  - Random seed: {RANDOM_SEED}")
    print("  - SNP threshold: >1% combined gnomAD NFE AF")
    print("    * gnomADe_NFE_AF if available")
    print("    * otherwise gnomADg_NFE_AF")
    print("\n" + "=" * 70 + "\n")

    # Load data
    validate_file_exists(INPUT_FILE, "Script 13 input")

    df = pd.read_excel(INPUT_FILE, sheet_name="Biological_Annotations")

    # Identify columns
    sym_c = find_col(df, "SYMBOL")
    tis_c = find_col(df, "TISSUE")
    sam_c = find_col(df, "SAMPLE")
    var_key_c = find_col(df, "Variant_Key")
    chrom_c = find_col(df, "CHROM")
    pos_c = find_col(df, "POS")
    ref_c = find_col(df, "REF")
    alt_c = find_col(df, "ALT")
    coh_c = find_col(df, "COHORT")
    exome_nfe_c = find_col(df, "gnomADe_NFE_AF")
    genome_nfe_c = find_col(df, "gnomADg_NFE_AF")
    imp_c = find_col(df, "IMPACT")
    con_c = find_col(df, "CONSEQUENCE")

    validate_required_columns(df, [sym_c, tis_c, sam_c, chrom_c, pos_c, ref_c, alt_c, coh_c, exome_nfe_c, genome_nfe_c, imp_c, con_c], "Script 13 Biological_Annotations")

    # --- Build combined NFE AF column ---
    df = combine_gnomad_nfe(df, exome_nfe_c, genome_nfe_c)

    # --- Filter for common SNPs (>1% exome or genome NFE AF) ---
    df = df[common_nfe_variant_mask(df, THRESHOLDS["min_nfe_af"], exome_nfe_c, genome_nfe_c)].copy()
    validate_nonempty(df, "Script 13 common SNP subset")

    # Standardise labels
    df[tis_c] = standardize_tissue_labels(df[tis_c])
    df[coh_c] = standardize_cohort_labels(df[coh_c])

    # Keep expected tissues only
    df = df[df[tis_c].isin(["Tumour", "Healthy"])].copy()
    validate_tissue_values(df, tis_c, "Script 13 filtered tissues")

    # Create genomic variant IDs
    df["Variant_ID"] = build_genomic_variant_id_series(
        variant_key_series=df[var_key_c] if var_key_c else None,
        chrom_series=df[chrom_c],
        pos_series=df[pos_c],
        ref_series=df[ref_c],
        alt_series=df[alt_c],
    )
    if "Existing_variation" in df.columns:
        df["rsID"] = df["Existing_variation"].apply(extract_rsid)
    else:
        df["rsID"] = pd.NA
    df["Variant_Display"] = df["rsID"].fillna(df["Variant_ID"])
    warning_lookup = build_amplicon_warning_lookup(df, variant_col="Variant_ID", chrom_col="CHROM", pos_col="POS")
    print_validation_summary(df, sam_c, "Script 13 filtered SNPs", [coh_c, tis_c])

    # Analysis groups
    cohort_list = ["Global", "Breast", "Endometrium"]
    cohort_filters = {"Breast": "Breast", "Endometrium": "Endometri"}
    all_results = []

    for cohort_name in cohort_list:
        print(f"\n{'─' * 70}")
        print(f"Processing: {cohort_name} Cohort")
        print(f"{'─' * 70}")

        if cohort_name == "Global":
            cohort_df = pooled_control_frame(df, coh_c, tis_c)
        else:
            cohort_df = pooled_control_frame(df, coh_c, tis_c, cohort_name, cohort_filters[cohort_name])

        if cohort_df.empty:
            print(f"  No data found for {cohort_name}")
            continue

        # Count total samples
        n_tumour = cohort_df[cohort_df[tis_c] == "Tumour"][sam_c].nunique()
        n_healthy = cohort_df[cohort_df[tis_c] == "Healthy"][sam_c].nunique()

        if n_tumour == 0 or n_healthy == 0:
            print(f"  Insufficient data: Tumour={n_tumour}, Control={n_healthy}")
            continue

        print(f"  Total samples: Tumour={n_tumour}, Control={n_healthy}")

        unique_snps = cohort_df["Variant_ID"].dropna().unique()
        print(f"  Unique SNPs to test: {len(unique_snps)}")

        for var in tqdm(unique_snps, desc=f"  Testing {cohort_name} SNPs"):
            var_data = cohort_df[cohort_df["Variant_ID"] == var]

            # Count carriers
            count_tumour = var_data[var_data[tis_c] == "Tumour"][sam_c].nunique()
            count_healthy = var_data[var_data[tis_c] == "Healthy"][sam_c].nunique()

            # Permutation test
            perm_result = permutation_test_snp(
                tumour_carriers=count_tumour,
                tumour_total=n_tumour,
                healthy_carriers=count_healthy,
                healthy_total=n_healthy,
                n_permutations=N_PERMUTATIONS,
                random_state=RANDOM_SEED
            )

            # Fisher exact test for comparison
            a, b = count_tumour, n_tumour - count_tumour
            c, d = count_healthy, n_healthy - count_healthy
            fisher_or, fisher_p = stats.fisher_exact([[a, b], [c, d]])

            all_results.append({
                "Cohort": cohort_name,
                "Variant_Display": var_data["Variant_Display"].iloc[0],
                "rsID": var_data["rsID"].iloc[0],
                "Variant_ID": var,
                "Symbol": var_data[sym_c].iloc[0],
                "Impact": var_data[imp_c].iloc[0],
                "Consequence": var_data[con_c].iloc[0],
                "gnomAD_NFE_AF": var_data["gnomAD_NFE_AF_combined"].iloc[0],
                "gnomAD_NFE_Source": var_data["gnomAD_NFE_Source"].iloc[0],
                "Tumour_Carriers": count_tumour,
                "Tumour_Total": n_tumour,
                "Tumour_Freq_%": perm_result["tumour_prop"] * 100,
                "Control_Carriers": count_healthy,
                "Control_Total": n_healthy,
                "Control_Freq_%": perm_result["healthy_prop"] * 100,
                "Freq_Difference_%": perm_result["observed_diff"] * 100,
                "Permutation_P_Value": perm_result["p_value"],
                "Effect_Size_h": perm_result["effect_size_h"],
                "Fisher_P_Value": fisher_p,
                "Fisher_OR": fisher_or,
                "N_Permutations": N_PERMUTATIONS
            })

    if not all_results:
        print("No permutation results were generated.")
        return

    # Convert to DataFrame and apply LD-aware block FDR
    results_df = pd.DataFrame(all_results)
    ld_reference = load_common_snp_ld_reference(
        annotated_report=PATHS["annotated_report"],
        phased_genotypes=PATHS["haplotype_phased"],
        min_nfe_af=THRESHOLDS["min_nfe_af"],
        common_snp_whitelist=PATHS["common_snps_dir"] / "GSDMB_Common_SNP_Frequency_Summary.xlsx",
        threshold=LD_R2_THRESHOLD,
    )
    results_df, ld_block_membership = _initialise_ld_block_columns(results_df, ld_reference)
    final_results = _apply_block_fdr(
        results_df,
        cohort_col="Cohort",
        block_col="LD_Block_ID",
        specs=[
            ("Permutation_P_Value", "Permutation_FDR", "Permutation_FDR_Global_Legacy"),
            ("Fisher_P_Value", "Fisher_FDR", "Fisher_FDR_Global_Legacy"),
        ],
    ).sort_values(["Cohort", "Permutation_P_Value"])

    # Significance flags
    validate_percentage_columns(final_results, ["Tumour_Freq_%", "Control_Freq_%"], "Script 13 results")
    final_results["Permutation_Significant_RAW"] = final_results["Permutation_P_Value"] < 0.05
    final_results["Fisher_Significant_RAW"] = final_results["Fisher_P_Value"] < 0.05
    final_results["Permutation_Significant_FDR"] = final_results["Permutation_FDR"] < 0.05
    final_results["Fisher_Significant_FDR"] = final_results["Fisher_FDR"] < 0.05
    final_results["Permutation_Significant_FDR_Global_Legacy"] = final_results["Permutation_FDR_Global_Legacy"] < 0.05
    final_results["Fisher_Significant_FDR_Global_Legacy"] = final_results["Fisher_FDR_Global_Legacy"] < 0.05
    final_results = attach_amplicon_warning_columns(final_results, warning_lookup)

    final_results["Agreement_RAW"] = (
        final_results["Permutation_Significant_RAW"] ==
        final_results["Fisher_Significant_RAW"]
    )
    final_results["Agreement_FDR"] = (
        final_results["Permutation_Significant_FDR"] ==
        final_results["Fisher_Significant_FDR"]
    )

    # Save results
    print(f"\n{'=' * 70}")
    print("SAVING RESULTS")
    print(f"{'=' * 70}\n")

    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
        summary = final_results.groupby("Cohort").agg({
            "Variant_ID": "count",
            "Permutation_Significant_RAW": "sum",
            "Fisher_Significant_RAW": "sum",
            "Permutation_Significant_FDR": "sum",
            "Fisher_Significant_FDR": "sum",
            "Agreement_RAW": "sum",
            "Agreement_FDR": "sum"
        }).reset_index()

        summary.columns = [
            "Cohort",
            "Total_SNPs",
            "Perm_Significant_RAW",
            "Fisher_Significant_RAW",
            "Perm_Significant_FDR",
            "Fisher_Significant_FDR",
            "Test_Agreement_RAW",
            "Test_Agreement_FDR"
        ]
        summary.to_excel(writer, sheet_name="Summary", index=False)

        for cohort in cohort_list:
            cohort_data = final_results[final_results["Cohort"] == cohort]
            if not cohort_data.empty:
                cohort_data.to_excel(writer, sheet_name=cohort, index=False)
        if not ld_block_membership.empty:
            ld_block_membership.to_excel(writer, sheet_name="LD_Block_Membership", index=False)

    print(f"Results saved to: {OUTPUT_XLSX}\n")

    # Generate comparison visualisation
    generate_comparison_plots(final_results)

    # Print summary statistics
    print(f"{'=' * 70}")
    print("SUMMARY STATISTICS")
    print(f"{'=' * 70}\n")

    for cohort in cohort_list:
        cohort_data = final_results[final_results["Cohort"] == cohort]
        if cohort_data.empty:
            continue

        print(f"{cohort} Cohort:")
        print(f"  Total SNPs tested: {len(cohort_data)}")
        print(f"  Permutation sig. (raw p<0.05): {cohort_data['Permutation_Significant_RAW'].sum()}")
        print(f"  Fisher sig. (raw p<0.05): {cohort_data['Fisher_Significant_RAW'].sum()}")
        print(f"  Permutation sig. (FDR<0.05): {cohort_data['Permutation_Significant_FDR'].sum()}")
        print(f"  Fisher sig. (FDR<0.05): {cohort_data['Fisher_Significant_FDR'].sum()}")
        print(
            f"  Test agreement (raw): {cohort_data['Agreement_RAW'].sum()}/{len(cohort_data)} "
            f"({cohort_data['Agreement_RAW'].mean() * 100:.1f}%)"
        )
        print(
            f"  Test agreement (FDR): {cohort_data['Agreement_FDR'].sum()}/{len(cohort_data)} "
            f"({cohort_data['Agreement_FDR'].mean() * 100:.1f}%)"
        )

        top5 = cohort_data.nsmallest(5, "Permutation_P_Value")
        if not top5.empty:
            print("\n  Top 5 SNPs by permutation test:")
            for _, row in top5.iterrows():
                print(
                    f"    {row['Variant_Display']} "
                    f"(p={row['Permutation_P_Value']:.4e}, "
                    f"perm_FDR={row['Permutation_FDR']:.4e}, "
                    f"source={row['gnomAD_NFE_Source']})"
                )
        print()

    print(f"{'=' * 70}")
    print("Permutation analysis complete!")
    print(f"{'=' * 70}\n")


def generate_comparison_plots(results_df):
    """Generate visualisation comparing permutation vs Fisher's exact test."""
    print("Generating comparison plots...")

    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))

    cohorts = results_df["Cohort"].unique()
    colors = {"Global": "#3498db", "Breast": "#e74c3c", "Endometrium": "#2ecc71"}

    # Plot 1: raw p-value comparison
    ax1 = axes[0, 0]
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        ax1.scatter(
            -np.log10(data["Fisher_P_Value"].replace(0, 1e-300)),
            -np.log10(data["Permutation_P_Value"].replace(0, 1e-300)),
            label=cohort,
            alpha=0.6,
            s=80,
            color=colors.get(cohort, "#95a5a6")
        )

    max_val = max(ax1.get_xlim()[1], ax1.get_ylim()[1])
    ax1.plot([0, max_val], [0, max_val], "k--", alpha=0.3, label="x = y")
    sig_line = -np.log10(0.05)
    ax1.axhline(sig_line, color="red", linestyle=":", alpha=0.5, label="P = 0.05")
    ax1.axvline(sig_line, color="red", linestyle=":", alpha=0.5)

    ax1.set_xlabel("Fisher exact test (-log10 raw P value)", fontweight="bold")
    ax1.set_ylabel("Permutation test (-log10 raw P value)", fontweight="bold")
    ax1.set_title("Raw P-value concordance", fontweight="bold", fontsize=12)
    ax1.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Plot 2: effect size vs permutation raw p-value
    ax2 = axes[0, 1]
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        ax2.scatter(
            data["Effect_Size_h"],
            -np.log10(data["Permutation_P_Value"].replace(0, 1e-300)),
            label=cohort,
            alpha=0.6,
            s=80,
            color=colors.get(cohort, "#95a5a6")
        )

    ax2.axhline(sig_line, color="red", linestyle=":", alpha=0.5, label="P = 0.05")
    ax2.set_xlabel("Effect size (Cohen's h)", fontweight="bold")
    ax2.set_ylabel("Permutation test (-log10 raw P value)", fontweight="bold")
    ax2.set_title("Effect size versus raw significance", fontweight="bold", fontsize=12)
    ax2.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Plot 3: frequency difference distribution (density-normalised)
    ax3 = axes[1, 0]
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        ax3.hist(
            data["Freq_Difference_%"],
            bins=30,
            density=True,
            alpha=0.45,
            label=cohort,
            color=colors.get(cohort, "#95a5a6"),
            edgecolor="black"
        )

    ax3.axvline(0, color="black", linestyle="--", alpha=0.5, label="No difference")
    ax3.set_xlabel("Carrier-frequency difference (tumour minus control, percentage points)", fontweight="bold")
    ax3.set_ylabel("Density", fontweight="bold")
    ax3.set_title("Tumour-control carrier-frequency shift", fontweight="bold", fontsize=12)
    ax3.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8)
    ax3.grid(True, alpha=0.3, axis="y")

    # Plot 4: agreement between tests using FDR significance (percentage of SNPs)
    ax4 = axes[1, 1]
    agreement_data = []
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        total_snps = len(data)
        if total_snps == 0:
            continue
        both_sig = ((data["Permutation_Significant_FDR"]) & (data["Fisher_Significant_FDR"])).sum()
        perm_only = ((data["Permutation_Significant_FDR"]) & (~data["Fisher_Significant_FDR"])).sum()
        fisher_only = ((~data["Permutation_Significant_FDR"]) & (data["Fisher_Significant_FDR"])).sum()
        neither = ((~data["Permutation_Significant_FDR"]) & (~data["Fisher_Significant_FDR"])).sum()

        agreement_data.append({
            "Cohort": cohort,
            "Both Significant": both_sig / total_snps * 100,
            "Perm Only": perm_only / total_snps * 100,
            "Fisher Only": fisher_only / total_snps * 100,
            "Neither": neither / total_snps * 100
        })

    agreement_df = pd.DataFrame(agreement_data).set_index("Cohort")
    agreement_df.plot(
        kind="bar",
        stacked=True,
        ax=ax4,
        color=["#2ecc71", "#f39c12", "#e74c3c", "#95a5a6"]
    )
    ax4.set_xlabel("Cohort", fontweight="bold")
    ax4.set_ylabel("Percentage of SNPs", fontweight="bold")
    ax4.set_title("FDR concordance by test", fontweight="bold", fontsize=12)
    ax4.legend(title="Significance", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.setp(ax4.xaxis.get_majorticklabels(), rotation=0)
    ax4.set_ylim(0, 100)
    ax4.grid(True, alpha=0.3, axis="y")

    plt.suptitle(
        f"Permutation vs Fisher Comparison\n"
        f"{N_PERMUTATIONS:,} permutations per SNP",
        fontsize=16,
        fontweight="bold",
        y=0.995
    )
    plt.tight_layout(rect=[0, 0, 0.92, 1])
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Comparison plots saved to: {OUTPUT_PLOT}\n")


if __name__ == "__main__":
    run_permutation_analysis()
