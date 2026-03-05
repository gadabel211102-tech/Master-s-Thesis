#!/usr/bin/env python3
"""
================================================================================
Script 13 — Permutation Testing for SNP Associations
================================================================================
Pipeline    : Step 13 of 17 — runs AFTER statistical enrichment analysis (12)
              and BEFORE the interactive dashboard (14)

--------------------------------------------------------------------------------
PURPOSE
--------------------------------------------------------------------------------
This script re-tests the SNP associations identified in script 12 using
permutation testing — a non-parametric statistical method that does not rely
on any assumptions about the distribution of the data. The results are compared
directly to the Fisher's exact test results from script 12 to assess robustness
and agreement between the two approaches.

The script addresses a key limitation of script 12: Fisher's exact test, while
valid for small samples, assumes that the observed cell counts follow a
hypergeometric distribution. When sample sizes are small or strongly imbalanced
between tissue groups, this assumption may not hold perfectly, and the p-values
may not be reliable. Permutation testing makes no distributional assumptions
and is therefore considered the more robust approach in these circumstances.

--------------------------------------------------------------------------------
WHAT IS PERMUTATION TESTING?
--------------------------------------------------------------------------------
A permutation test answers the following question:

  "If tissue type (tumour vs. healthy) had no effect on SNP carrier status,
   how often would we observe a difference as large as the one we measured,
   just by chance?"

The test works as follows (for one SNP):

  1. Observe the actual difference in carrier frequency between tumour and
     healthy samples (the "observed statistic").

  2. Pool all samples together (tumour and healthy), randomly shuffle the
     tissue labels, and re-split them into two groups of the same sizes as
     the original. Compute the difference in carrier frequency between these
     randomly shuffled groups (one "null" statistic).

  3. Repeat step 2 many times (here: 10,000 permutations), building up a
     distribution of null statistics that represents what differences would
     look like if tissue type had no effect.

  4. The permutation p-value is the proportion of null statistics that are as
     extreme as, or more extreme than, the observed statistic:

        p = (number of |null diffs| ≥ |observed diff|) / N_permutations

  5. This is a two-tailed test: we count both permutations where the null
     difference is ≥ the observed difference AND where it is ≤ its negative
     (i.e. equally extreme in the opposite direction).

A small p-value (< 0.05) means the observed difference is rarely seen under
the null hypothesis of no association, suggesting a real relationship.

Because permutation is computationally intensive (10,000 shuffles per SNP),
a progress bar (tqdm) is shown during execution.

--------------------------------------------------------------------------------
RANDOM SEED AND REPRODUCIBILITY
--------------------------------------------------------------------------------
The random number generator is seeded with RANDOM_SEED = 42 before each
permutation test. This ensures that re-running the script on the same data
produces identical p-values — a requirement for reproducible scientific
research. Any integer can be used as the seed; 42 is the conventional default.

--------------------------------------------------------------------------------
EFFECT SIZE: COHEN'S h
--------------------------------------------------------------------------------
Statistical significance (p-value) tells us whether an association is likely
to be real, but it does not tell us how large or meaningful the effect is.
A large study can produce a very small p-value for a trivially small difference.

Cohen's h is the standard effect size measure for comparisons between two
proportions. It is calculated as:

  h = 2 × (arcsin(√p₁) − arcsin(√p₂))

where p₁ is the tumour carrier frequency and p₂ is the healthy carrier
frequency. The arcsin transformation stabilises the variance of proportions.

Interpretation:
  |h| < 0.2   → small effect (likely not biologically meaningful)
  |h| ≈ 0.5   → medium effect
  |h| ≥ 0.8   → large effect (biologically meaningful)

A positive h indicates the SNP is more common in tumour; negative h indicates
it is less common (potentially protective).

--------------------------------------------------------------------------------
COMPARISON WITH FISHER'S EXACT TEST
--------------------------------------------------------------------------------
Fisher's exact test is also run for every SNP (replicating script 12) so that
both p-values appear side by side in the output. This enables direct comparison:

  Agreement    — both tests agree on significance (both p < 0.05 or both ≥ 0.05)
  Perm only    — permutation test significant but Fisher is not
  Fisher only  — Fisher significant but permutation test is not
  Neither      — both tests non-significant

Strong agreement between the two methods increases confidence in the findings.
Discordant results (where one test is significant and the other is not) may
indicate that one of the assumptions of Fisher's test is being violated, or
that the permutation distribution is irregular due to very small counts.

--------------------------------------------------------------------------------
MULTIPLE TESTING CORRECTION
--------------------------------------------------------------------------------
Benjamini-Hochberg FDR correction (as in script 12) is applied separately to
both the permutation p-values and the Fisher p-values, within each cohort.
Significance is assessed at FDR q < 0.05.

--------------------------------------------------------------------------------
FOUR-PANEL COMPARISON FIGURE
--------------------------------------------------------------------------------
Panel 1 (top-left):  P-value comparison scatter
  Each dot is a SNP. X-axis = −log10(Fisher p); Y-axis = −log10(permutation p).
  The diagonal line (x = y) shows perfect agreement between the two methods.
  Points above the diagonal = permutation more significant than Fisher.
  Points below the diagonal = Fisher more significant than permutation.

Panel 2 (top-right): Effect size vs permutation significance
  X-axis = Cohen's h (effect size); Y-axis = −log10(permutation p).
  Significant SNPs with large effect sizes are the most biologically important.

Panel 3 (bottom-left): Distribution of frequency differences
  Histogram of (Tumour% − Healthy%) across all tested SNPs.
  Centred near zero means most SNPs show no difference between tissue types.
  Skewed distributions or outliers may indicate biologically relevant variants.

Panel 4 (bottom-right): Test agreement stacked bar chart
  Per cohort, how many SNPs fall into each agreement category:
  Both significant / Permutation only / Fisher only / Neither.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        Required columns: SYMBOL, TISSUE, SAMPLE, Existing_variation, HGVSp,
                          COHORT, gnomADe_NFE_AF, IMPACT, CONSEQUENCE

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------
  13_Permutation_Test_Results.xlsx
    Sheets: Summary | Global | Breast | Endometrium
    Each SNP row includes: carrier counts, frequencies, permutation p-value,
    Fisher p-value, BH-FDR corrections for both, Cohen's h effect size,
    and test agreement flags.

  13_Permutation_Comparison.png
    Four-panel comparison figure (300 dpi).

--------------------------------------------------------------------------------
REFERENCES
--------------------------------------------------------------------------------
  Good PI (2005). Permutation, Parametric, and Bootstrap Tests of Hypotheses.
    3rd ed. Springer.

  Cohen J (1988). Statistical Power Analysis for the Behavioral Sciences.
    2nd ed. Lawrence Erlbaum Associates.
    (defines Cohen's h for comparing two proportions)

  Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
    practical and powerful approach to multiple testing.
    J R Stat Soc B, 57(1):289–300.

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python 13_permutation_testing.py

  NOTE: With N_PERMUTATIONS = 10,000 and many SNPs, this script may take
  several minutes to run. A progress bar (tqdm) tracks progress per cohort.

  (Paths and parameters are configured in the CONFIGURATION block below.)

--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------
  Python ≥ 3.8
  pandas, numpy, scipy, statsmodels, matplotlib, seaborn, tqdm, openpyxl

  Install: pip install pandas numpy scipy statsmodels matplotlib seaborn tqdm openpyxl
================================================================================
"""

# ──────────────────────────────────────────────────────────────────────────────
# STANDARD LIBRARY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import os       # File path construction and existence checks
import warnings # Suppress non-critical library warnings

# ──────────────────────────────────────────────────────────────────────────────
# THIRD-PARTY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt                        # Figure creation
import numpy as np                                     # Array operations, random shuffle
import pandas as pd                                    # Data loading and manipulation
import seaborn as sns                                  # Plot styling
from scipy import stats                                # Statistical functions
from scipy.stats import fisher_exact, permutation_test # Fisher's exact test
from statsmodels.stats.multitest import multipletests  # Benjamini-Hochberg FDR correction
from tqdm import tqdm                                  # Progress bar for permutation loops

warnings.filterwarnings("ignore")


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

BASE_PATH   = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
INPUT_FILE  = os.path.join(BASE_PATH, "GSDMB_Annotated_Report_Fixed.xlsx")
OUTPUT_XLSX = os.path.join(BASE_PATH, "13_Permutation_Test_Results.xlsx")
OUTPUT_PLOT = os.path.join(BASE_PATH, "13_Permutation_Comparison.png")

# Number of permutations per SNP. 10,000 is the standard recommendation for
# publication-quality results: it provides a minimum achievable p-value of
# 0.0001 (1/10,000), sufficient resolution for typical FDR thresholds.
# Increase to 100,000 for higher precision; decrease for faster testing runs.
N_PERMUTATIONS = 10_000

# Fixed random seed for reproducibility. Any integer is valid; 42 is the
# conventional default. Changing this value will produce slightly different
# p-values due to the stochastic nature of random shuffling.
RANDOM_SEED = 42

# Minimum gnomAD NFE allele frequency to classify a variant as a SNP.
# Consistent with scripts 11 and 12.
SNP_AF_THRESHOLD = 0.01

# FDR significance threshold applied to Benjamini-Hochberg corrected q-values.
FDR_THRESHOLD = 0.05


# ══════════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def find_col(df: pd.DataFrame, target: str):
    """
    Case-insensitive column name lookup, also stripping leading/trailing spaces.

    Parameters
    ----------
    df     : pd.DataFrame — the DataFrame to search
    target : str          — the desired column name (any case)

    Returns
    -------
    str or None
        The actual column name as it appears in the DataFrame, or None.
    """
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None


# ══════════════════════════════════════════════════════════════════════════════
# CORE STATISTICAL FUNCTION: PERMUTATION TEST
# ══════════════════════════════════════════════════════════════════════════════

def permutation_test_snp(
    tumour_carriers:  int,
    tumour_total:     int,
    healthy_carriers: int,
    healthy_total:    int,
    n_permutations:   int = N_PERMUTATIONS,
    random_state:     int = RANDOM_SEED,
) -> dict:
    """
    Perform a permutation test for a single SNP's association with tissue type.

    Algorithm
    ---------
    1. Represent the tumour group as a binary array of length tumour_total,
       where 1 = carries the SNP and 0 = does not carry the SNP.
       Do the same for the healthy group.

    2. Compute the observed test statistic: the difference in carrier
       proportions between tumour and healthy samples.
         observed_diff = (tumour_carriers / tumour_total)
                         − (healthy_carriers / healthy_total)

    3. Pool both groups into one combined array. Repeat n_permutations times:
         a. Randomly shuffle the pooled array (this destroys any real tissue-type
            signal, simulating the null hypothesis of no association).
         b. Re-split into two groups of the original sizes.
         c. Compute the proportion difference for this null permutation.

    4. The two-tailed permutation p-value is the proportion of null differences
       whose absolute value equals or exceeds the absolute observed difference:
         p = mean(|null_diffs| ≥ |observed_diff|)

    5. Compute Cohen's h as the effect size:
         h = 2 × (arcsin(√p_tumour) − arcsin(√p_healthy))

    Parameters
    ----------
    tumour_carriers  : int — number of tumour samples carrying the SNP
    tumour_total     : int — total number of tumour samples in this cohort
    healthy_carriers : int — number of healthy samples carrying the SNP
    healthy_total    : int — total number of healthy samples in this cohort
    n_permutations   : int — number of random shuffles to perform
    random_state     : int — random seed for reproducibility

    Returns
    -------
    dict with keys:
        observed_diff  (float) — tumour proportion minus healthy proportion
        p_value        (float) — two-tailed permutation p-value
        effect_size_h  (float) — Cohen's h effect size
        perm_diffs     (array) — the full null distribution (for diagnostics)
        tumour_prop    (float) — observed tumour carrier proportion
        healthy_prop   (float) — observed healthy carrier proportion
    """
    np.random.seed(random_state)

    # Build binary carrier arrays: 1 = carries SNP, 0 = does not carry SNP
    tumour_data  = np.array([1] * tumour_carriers  + [0] * (tumour_total  - tumour_carriers))
    healthy_data = np.array([1] * healthy_carriers + [0] * (healthy_total - healthy_carriers))

    # Observed carrier proportions and their difference
    obs_tumour_prop  = tumour_carriers  / tumour_total  if tumour_total  > 0 else 0.0
    obs_healthy_prop = healthy_carriers / healthy_total if healthy_total > 0 else 0.0
    observed_diff    = obs_tumour_prop - obs_healthy_prop

    # Pool all samples; record original group size for re-splitting
    all_data = np.concatenate([tumour_data, healthy_data])
    n_tumour = len(tumour_data)

    # Build the null distribution by repeated random shuffling
    null_diffs = np.empty(n_permutations)
    for i in range(n_permutations):
        shuffled         = np.random.permutation(all_data)
        perm_tumour      = shuffled[:n_tumour]
        perm_healthy     = shuffled[n_tumour:]
        null_diffs[i]    = np.mean(perm_tumour) - np.mean(perm_healthy)

    # Two-tailed p-value: fraction of null diffs as or more extreme than observed
    p_value = np.mean(np.abs(null_diffs) >= np.abs(observed_diff))

    # Cohen's h effect size for two proportions
    # arcsin transformation stabilises variance; factor of 2 scales to a
    # conventional effect size metric (Cohen, 1988)
    effect_h = 2.0 * (
        np.arcsin(np.sqrt(obs_tumour_prop)) -
        np.arcsin(np.sqrt(obs_healthy_prop))
    )

    return {
        "observed_diff": observed_diff,
        "p_value":       p_value,
        "effect_size_h": effect_h,
        "perm_diffs":    null_diffs,
        "tumour_prop":   obs_tumour_prop,
        "healthy_prop":  obs_healthy_prop,
    }


# ══════════════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def run_permutation_analysis():
    """
    Orchestrate the full permutation testing pipeline.

    Execution order
    ---------------
    1.  Load and filter the annotated variant table to SNPs (AF > 1%).
    2.  Standardise tissue/cohort labels; construct Variant IDs.
    3.  For each cohort (Global, Breast, Endometrium):
          For each unique SNP:
            a. Count tumour and healthy carriers.
            b. Run the permutation test (10,000 shuffles).
            c. Run Fisher's exact test for side-by-side comparison.
            d. Record all statistics.
    4.  Apply Benjamini-Hochberg FDR correction (separately for permutation
        and Fisher p-values, within each cohort).
    5.  Flag agreement/disagreement between the two test methods.
    6.  Save results to a multi-sheet Excel workbook.
    7.  Generate and save the four-panel comparison figure.
    8.  Print a summary of significant findings per cohort.
    """
    print("=" * 70)
    print("SCRIPT 13 — Permutation Testing for SNP Associations")
    print("=" * 70)
    print(f"\n  Input            : {INPUT_FILE}")
    print(f"  Permutations/SNP : {N_PERMUTATIONS:,}")
    print(f"  Random seed      : {RANDOM_SEED}")
    print(f"  SNP threshold    : gnomAD NFE AF > {SNP_AF_THRESHOLD}")
    print(f"  FDR threshold    : q < {FDR_THRESHOLD}")
    print("\n" + "=" * 70 + "\n")

    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: Input file not found: {INPUT_FILE}")
        print("       Please update 'BASE_PATH' in the CONFIGURATION block.")
        return

    # ── Step 1: Load data ─────────────────────────────────────────────────────
    df = pd.read_excel(INPUT_FILE, sheet_name="Biological_Annotations")
    print(f"Loaded {len(df)} rows from sheet 'Biological_Annotations'")

    # ── Step 2: Identify columns ──────────────────────────────────────────────
    sym_c = find_col(df, "SYMBOL")
    tis_c = find_col(df, "TISSUE")
    sam_c = find_col(df, "SAMPLE")
    var_c = find_col(df, "Existing_variation")
    hgv_c = find_col(df, "HGVSp")
    coh_c = find_col(df, "COHORT")
    nfe_c = find_col(df, "gnomADe_NFE_AF")
    imp_c = find_col(df, "IMPACT")
    con_c = find_col(df, "CONSEQUENCE")

    # ── Step 3: Filter and standardise ───────────────────────────────────────
    df = df[df[nfe_c] > SNP_AF_THRESHOLD].copy()
    print(f"Variants retained after SNP filter: {len(df)}\n")

    df[tis_c] = (
        df[tis_c].astype(str).str.strip()
        .replace({"Tumor": "Tumour", "Normal": "Healthy", "Control": "Healthy"})
    )
    df[coh_c] = df[coh_c].astype(str).str.strip()

    # Construct Variant IDs (rsID if available, else GENE:HGVSp)
    df["rsID"]       = df[var_c].astype(str).str.extract(r"(rs\d+)")
    df["Variant_ID"] = df["rsID"].fillna(
        df[sym_c].astype(str) + ":" + df[hgv_c].astype(str)
    )

    # ══════════════════════════════════════════════════════════════════════════
    # Step 4: Permutation test loop (one SNP × cohort combination at a time)
    # ══════════════════════════════════════════════════════════════════════════
    cohort_list = ["Global", "Breast", "Endometrium"]
    all_results = []

    for cohort_name in cohort_list:
        print(f"{'─' * 70}")
        print(f"Processing cohort: {cohort_name}")
        print(f"{'─' * 70}")

        cohort_df = (
            df.copy() if cohort_name == "Global"
            else df[df[coh_c].str.contains(cohort_name, case=False, na=False)].copy()
        )

        if cohort_df.empty:
            print(f"  No data found for cohort '{cohort_name}' — skipping.\n")
            continue

        n_tumour  = cohort_df[cohort_df[tis_c] == "Tumour"][sam_c].nunique()
        n_healthy = cohort_df[cohort_df[tis_c] == "Healthy"][sam_c].nunique()

        if n_tumour == 0 or n_healthy == 0:
            print(f"  Missing tissue group — Tumour={n_tumour}, Healthy={n_healthy}. Skipping.\n")
            continue

        unique_snps = cohort_df["Variant_ID"].unique()
        print(f"  Samples: {n_tumour} tumour, {n_healthy} healthy")
        print(f"  Unique SNPs to test: {len(unique_snps)}")

        # tqdm provides a real-time progress bar during the permutation loop,
        # which is important because 10,000 × N_snps iterations can be slow
        for variant_id in tqdm(unique_snps, desc=f"  Permuting {cohort_name} SNPs"):

            var_data      = cohort_df[cohort_df["Variant_ID"] == variant_id]
            count_tumour  = var_data[var_data[tis_c] == "Tumour"][sam_c].nunique()
            count_healthy = var_data[var_data[tis_c] == "Healthy"][sam_c].nunique()

            # ── Permutation test ──────────────────────────────────────────────
            perm_result = permutation_test_snp(
                tumour_carriers  = count_tumour,
                tumour_total     = n_tumour,
                healthy_carriers = count_healthy,
                healthy_total    = n_healthy,
                n_permutations   = N_PERMUTATIONS,
                random_state     = RANDOM_SEED,
            )

            # ── Fisher's exact test (for comparison with script 12) ───────────
            a, b = count_tumour,  max(0, n_tumour  - count_tumour)
            c, d = count_healthy, max(0, n_healthy - count_healthy)
            fisher_or, fisher_p = fisher_exact([[a, b], [c, d]])

            all_results.append({
                "Cohort":             cohort_name,
                "Variant_ID":         variant_id,
                "Symbol":             var_data[sym_c].iloc[0],
                "Impact":             var_data[imp_c].iloc[0],
                "Consequence":        var_data[con_c].iloc[0],
                "gnomAD_NFE_AF":      var_data[nfe_c].iloc[0],
                "Tumour_Carriers":    count_tumour,
                "Tumour_Total":       n_tumour,
                "Tumour_Freq_%":      perm_result["tumour_prop"]  * 100,
                "Healthy_Carriers":   count_healthy,
                "Healthy_Total":      n_healthy,
                "Healthy_Freq_%":     perm_result["healthy_prop"] * 100,
                "Freq_Difference_%":  perm_result["observed_diff"] * 100,
                "Permutation_P_Value": perm_result["p_value"],
                "Effect_Size_h":      perm_result["effect_size_h"],
                "Fisher_P_Value":     fisher_p,
                "Fisher_OR":          fisher_or,
                "N_Permutations":     N_PERMUTATIONS,
            })

    # ══════════════════════════════════════════════════════════════════════════
    # Step 5: Benjamini-Hochberg FDR correction (within each cohort)
    # ══════════════════════════════════════════════════════════════════════════
    # FDR correction is applied separately to permutation p-values and Fisher
    # p-values, and separately within each cohort (same rationale as script 12).
    results_df = pd.DataFrame(all_results)

    corrected_dfs = []
    for cohort in results_df["Cohort"].unique():
        sub = results_df[results_df["Cohort"] == cohort].copy()

        _, sub["Permutation_FDR"], _, _ = multipletests(
            sub["Permutation_P_Value"], method="fdr_bh"
        )
        _, sub["Fisher_FDR"], _, _ = multipletests(
            sub["Fisher_P_Value"], method="fdr_bh"
        )
        corrected_dfs.append(sub)

    final_results = (
        pd.concat(corrected_dfs)
        .sort_values(["Cohort", "Permutation_P_Value"])
        .reset_index(drop=True)
    )

    # ── Step 6: Flag significance and test agreement ──────────────────────────
    # Use raw p-values (not FDR) for the significance flags, consistent with
    # the exploratory nature of the agreement comparison.
    final_results["Permutation_Significant"] = (
        final_results["Permutation_P_Value"] < FDR_THRESHOLD
    )
    final_results["Fisher_Significant"] = (
        final_results["Fisher_P_Value"] < FDR_THRESHOLD
    )
    # Agreement = True if both tests reach the same significance conclusion
    final_results["Agreement"] = (
        final_results["Permutation_Significant"] == final_results["Fisher_Significant"]
    )

    # ── Step 7: Save results to Excel ────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print("SAVING RESULTS")
    print(f"{'=' * 70}\n")

    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:

        # Summary sheet: one row per cohort with aggregate counts
        summary = (
            final_results
            .groupby("Cohort")
            .agg(
                Total_SNPs            = ("Variant_ID", "count"),
                Perm_Significant      = ("Permutation_Significant", "sum"),
                Fisher_Significant    = ("Fisher_Significant", "sum"),
                Test_Agreement        = ("Agreement", "sum"),
            )
            .reset_index()
        )
        summary.to_excel(writer, sheet_name="Summary", index=False)

        # Individual cohort sheets
        for cohort in cohort_list:
            cohort_data = final_results[final_results["Cohort"] == cohort]
            if not cohort_data.empty:
                cohort_data.to_excel(writer, sheet_name=cohort, index=False)

    print(f"  Results saved: {OUTPUT_XLSX}\n")

    # ── Step 8: Generate comparison figure ───────────────────────────────────
    generate_comparison_plots(final_results)

    # ── Step 9: Print per-cohort summary ─────────────────────────────────────
    print(f"{'=' * 70}")
    print("SUMMARY STATISTICS")
    print(f"{'=' * 70}\n")

    for cohort in cohort_list:
        cohort_data = final_results[final_results["Cohort"] == cohort]
        if cohort_data.empty:
            continue

        n = len(cohort_data)
        n_perm_sig   = cohort_data["Permutation_Significant"].sum()
        n_fisher_sig = cohort_data["Fisher_Significant"].sum()
        n_agree      = cohort_data["Agreement"].sum()

        print(f"  {cohort} Cohort:")
        print(f"    SNPs tested               : {n}")
        print(f"    Permutation significant   : {n_perm_sig}  (p < {FDR_THRESHOLD})")
        print(f"    Fisher significant        : {n_fisher_sig}  (p < {FDR_THRESHOLD})")
        print(f"    Test agreement            : {n_agree}/{n} ({n_agree/n*100:.1f}%)")

        top5 = cohort_data.nsmallest(5, "Permutation_P_Value")
        if not top5.empty:
            print(f"    Top 5 SNPs (permutation p-value):")
            for _, row in top5.iterrows():
                print(f"      {row['Variant_ID']:30s}  "
                      f"p={row['Permutation_P_Value']:.4e}  "
                      f"FDR={row['Permutation_FDR']:.4e}  "
                      f"h={row['Effect_Size_h']:.3f}")
        print()

    print(f"{'=' * 70}")
    print("Script 13 complete.")
    print(f"{'=' * 70}\n")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE GENERATION
# ══════════════════════════════════════════════════════════════════════════════

def generate_comparison_plots(results_df: pd.DataFrame):
    """
    Generate a four-panel figure comparing permutation and Fisher's exact test
    results across all cohorts.

    Panel descriptions — see module docstring for full interpretation guide.

    Parameters
    ----------
    results_df : pd.DataFrame — full results table from run_permutation_analysis()
    """
    print("Generating comparison plots...")

    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    cohorts     = results_df["Cohort"].unique()
    sig_line    = -np.log10(FDR_THRESHOLD)  # Y/X position of the p=0.05 threshold line

    # Colour scheme: one colour per cohort, consistent across all four panels
    cohort_colours = {
        "Global":      "#3498db",
        "Breast":      "#e74c3c",
        "Endometrium": "#2ecc71",
    }

    # ── Panel 1: P-value comparison scatter ───────────────────────────────────
    ax1 = axes[0, 0]
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        ax1.scatter(
            -np.log10(data["Fisher_P_Value"]),
            -np.log10(data["Permutation_P_Value"]),
            label=cohort, alpha=0.6, s=80,
            color=cohort_colours.get(cohort, "#95a5a6"),
        )

    # Diagonal line of identity (x = y): perfect agreement between methods
    max_val = max(ax1.get_xlim()[1], ax1.get_ylim()[1])
    ax1.plot([0, max_val], [0, max_val], "k--", alpha=0.3, label="x = y (perfect agreement)")

    # Significance threshold lines on both axes
    ax1.axhline(sig_line, color="red", linestyle=":", alpha=0.5, label=f"p = {FDR_THRESHOLD}")
    ax1.axvline(sig_line, color="red", linestyle=":", alpha=0.5)

    ax1.set_xlabel("Fisher's Exact Test  (−log10 p)", fontweight="bold")
    ax1.set_ylabel("Permutation Test  (−log10 p)", fontweight="bold")
    ax1.set_title("P-Value Comparison:\nPermutation vs Fisher's Exact",
                  fontweight="bold", fontsize=12)
    ax1.legend(loc="best", fontsize=9)
    ax1.grid(True, alpha=0.3)

    # ── Panel 2: Effect size (Cohen's h) vs permutation significance ──────────
    ax2 = axes[0, 1]
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        ax2.scatter(
            data["Effect_Size_h"],
            -np.log10(data["Permutation_P_Value"]),
            label=cohort, alpha=0.6, s=80,
            color=cohort_colours.get(cohort, "#95a5a6"),
        )

    ax2.axhline(sig_line, color="red", linestyle=":", alpha=0.5, label=f"p = {FDR_THRESHOLD}")
    ax2.axvline(0, color="black", linestyle="--", alpha=0.3, label="h = 0 (no effect)")

    ax2.set_xlabel("Effect Size — Cohen's h", fontweight="bold")
    ax2.set_ylabel("Permutation Test  (−log10 p)", fontweight="bold")
    ax2.set_title("Effect Size vs Statistical Significance\n"
                  "(|h| < 0.2 small, |h| ≈ 0.5 medium, |h| ≥ 0.8 large)",
                  fontweight="bold", fontsize=11)
    ax2.legend(loc="best", fontsize=9)
    ax2.grid(True, alpha=0.3)

    # ── Panel 3: Distribution of frequency differences ────────────────────────
    ax3 = axes[1, 0]
    for cohort in cohorts:
        data = results_df[results_df["Cohort"] == cohort]
        ax3.hist(
            data["Freq_Difference_%"],
            bins=30, alpha=0.5, label=cohort,
            color=cohort_colours.get(cohort, "#95a5a6"),
            edgecolor="black",
        )

    # Vertical line at zero: no difference between tumour and healthy
    ax3.axvline(0, color="black", linestyle="--", alpha=0.5,
                label="No difference (0%)")
    ax3.set_xlabel("Frequency Difference: Tumour − Healthy (%)", fontweight="bold")
    ax3.set_ylabel("Number of SNPs", fontweight="bold")
    ax3.set_title("Distribution of Carrier Frequency Differences\n"
                  "(positive = enriched in tumour; negative = depleted)",
                  fontweight="bold", fontsize=11)
    ax3.legend(loc="best", fontsize=9)
    ax3.grid(True, alpha=0.3, axis="y")

    # ── Panel 4: Test agreement stacked bar chart ─────────────────────────────
    ax4 = axes[1, 1]
    agreement_rows = []
    for cohort in cohorts:
        data       = results_df[results_df["Cohort"] == cohort]
        both_sig   = ((data["Permutation_Significant"]) & ( data["Fisher_Significant"])).sum()
        perm_only  = ((data["Permutation_Significant"]) & (~data["Fisher_Significant"])).sum()
        fish_only  = ((~data["Permutation_Significant"]) & (data["Fisher_Significant"])).sum()
        neither    = ((~data["Permutation_Significant"]) & (~data["Fisher_Significant"])).sum()
        agreement_rows.append({
            "Cohort":            cohort,
            "Both Significant":  both_sig,
            "Permutation Only":  perm_only,
            "Fisher Only":       fish_only,
            "Neither":           neither,
        })

    agreement_df = pd.DataFrame(agreement_rows).set_index("Cohort")
    agreement_df.plot(
        kind="bar", stacked=True, ax=ax4,
        color=["#2ecc71", "#f39c12", "#e74c3c", "#95a5a6"],
    )
    ax4.set_xlabel("Cohort", fontweight="bold")
    ax4.set_ylabel("Number of SNPs", fontweight="bold")
    ax4.set_title("Test Agreement Between Methods\n"
                  "(Both Sig / Permutation Only / Fisher Only / Neither)",
                  fontweight="bold", fontsize=11)
    ax4.legend(title="Significance category",
               bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9)
    plt.setp(ax4.xaxis.get_majorticklabels(), rotation=0)
    ax4.grid(True, alpha=0.3, axis="y")

    # ── Save figure ───────────────────────────────────────────────────────────
    plt.suptitle(
        f"Permutation Testing vs Fisher's Exact Test — Comparison Summary\n"
        f"({N_PERMUTATIONS:,} permutations per SNP, seed={RANDOM_SEED})",
        fontsize=15, fontweight="bold", y=0.999,
    )
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"  Comparison figure saved: {OUTPUT_PLOT}\n")


# ── Script entry point ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    run_permutation_analysis()
