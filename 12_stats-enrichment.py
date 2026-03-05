#!/usr/bin/env python3
"""
================================================================================
Script 12 — SNP Statistical Enrichment Analysis
================================================================================
Pipeline    : Step 12 of 17 — runs AFTER SNP identification (11) and BEFORE
              permutation testing (13)

--------------------------------------------------------------------------------
PURPOSE
--------------------------------------------------------------------------------
This script tests whether any of the common SNPs identified in script 11 are
statistically enriched (more frequent) or depleted (less frequent) in tumour
tissue compared to healthy tissue within each cancer cohort.

The core question being asked is: for each SNP, is the probability of carrying
that variant significantly different between tumour samples and healthy samples?
A significant enrichment in tumour tissue could suggest that the variant
confers a cancer risk or plays a role in tumour development at the GSDMB locus.

Three parallel analyses are run:
  - Global     : all samples combined, regardless of cancer type
  - Breast     : breast cancer cohort only
  - Endometrium: endometrial cancer cohort only

--------------------------------------------------------------------------------
STATISTICAL METHOD: FISHER'S EXACT TEST
--------------------------------------------------------------------------------
For each SNP × cohort combination, a 2×2 contingency table is constructed:

                      Carries SNP    Does NOT carry SNP
  Tumour samples   |      a        |         b         |
  Healthy samples  |      c        |         d         |

Fisher's exact test calculates the probability of observing a table as extreme
as (or more extreme than) this one, assuming no true association between tissue
type and SNP carrier status (the null hypothesis). It is preferred over the
chi-squared test here because some cells in the contingency table may have
small counts, violating chi-squared's assumptions.

The result is a p-value: small p-values (typically < 0.05) indicate that the
observed difference in SNP frequency between tumour and healthy samples is
unlikely to have arisen by chance.

The Odds Ratio (OR) quantifies the direction and magnitude of any association:
  OR > 1  → SNP is more common in tumour than healthy (possible risk factor)
  OR = 1  → No difference (SNP equally common in both tissue types)
  OR < 1  → SNP is less common in tumour than healthy (possible protective effect)

--------------------------------------------------------------------------------
HALDANE-ANSCOMBE CORRECTION (+ 0.5)
--------------------------------------------------------------------------------
When any cell in the 2×2 table is zero (e.g. a SNP is completely absent from
all healthy samples, making c = 0), the standard Odds Ratio formula produces
a division by zero. The Haldane-Anscombe correction adds 0.5 to all four
cells before calculating the OR:

  OR_corrected = ((a + 0.5)(d + 0.5)) / ((b + 0.5)(c + 0.5))

This is a widely accepted small-sample correction that avoids infinite or
undefined odds ratios while introducing only minimal bias. All results
computed with this correction are flagged as "Corrected (+0.5)" in the output.

Fisher's exact test itself is always computed on the original (uncorrected)
counts — the correction applies only to the OR calculation.

--------------------------------------------------------------------------------
MULTIPLE TESTING CORRECTION: BENJAMINI-HOCHBERG FDR
--------------------------------------------------------------------------------
When testing many SNPs simultaneously (e.g. 50 SNPs per cohort), we expect
approximately 5% of tests to appear significant by chance even if no true
associations exist (false positives). Without correction, a study testing 100
SNPs at p < 0.05 would expect ~5 false positives.

The Benjamini-Hochberg (BH) procedure controls the False Discovery Rate (FDR):
the expected proportion of significant results that are false positives. It
ranks all p-values and adjusts them so that the FDR across all tests is
controlled at 5% (q < 0.05). This is less conservative than the Bonferroni
correction and better suited to exploratory genomic analyses where some false
positives are acceptable.

Results are reported with both the raw p-value and the BH-adjusted FDR
q-value. Statistical significance in this analysis requires FDR q < 0.05.

--------------------------------------------------------------------------------
VOLCANO PLOT
--------------------------------------------------------------------------------
The volcano plot is the standard visualisation for association analyses of
this kind. Each dot represents one SNP:
  - X-axis: Odds Ratio (log scale) — magnitude and direction of association
  - Y-axis: −log10(p-value) — statistical significance
            (higher = more significant; −log10(0.05) ≈ 1.3 is the threshold)

The plot is divided into four quadrants by reference lines:
  - Horizontal dashed red line at −log10(0.05): significance threshold
  - Vertical dashed blue line at OR = 1: line of no effect

Quadrant interpretation:
  Upper-right : significant enrichment in tumour (high OR, low p) → risk SNPs
  Upper-left  : significant depletion in tumour (low OR, low p)   → protective SNPs
  Lower half  : non-significant (above p = 0.05 threshold)

The top 5 most significant SNPs per cohort are labelled with their rsID or
gene:protein identifiers.

Two output figures are saved. Both currently use OR on the log scale; the
second (suffixed _log2OR.png) is reserved for potential future use of log2(OR)
on the X-axis if a linear-symmetric scale is preferred.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        Required columns: SYMBOL, IMPACT, CONSEQUENCE, TISSUE, SAMPLE,
                          Existing_variation, HGVSp, COHORT, gnomADe_NFE_AF

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------
  12_SNP_Enrichment_Results.xlsx
    Multi-sheet workbook, one sheet per analysis group (Global, Breast,
    Endometrium). Each row is one SNP with: frequencies, OR, raw p-value,
    BH-corrected FDR q-value, and correction method used.

  12_SNP_Volcano_Plots.png
    Three-panel volcano plot (one panel per cohort), OR on log scale.

  12_SNP_Volcano_Plots_log2OR.png
    Identical figure, reserved for log2(OR) X-axis variant.

--------------------------------------------------------------------------------
REFERENCES
--------------------------------------------------------------------------------
  Fisher RA (1922). On the interpretation of χ² from contingency tables, and
    the calculation of P. J R Stat Soc, 85(1):87–94.

  Haldane JBS (1956). The estimation and significance of the logarithm of a
    ratio of frequencies. Ann Hum Genet, 20(4):309–311.

  Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
    practical and powerful approach to multiple testing. J R Stat Soc B,
    57(1):289–300.

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python 12_stats-enrichment.py

  (Paths are configured in the CONFIGURATION block below.)

--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------
  Python ≥ 3.8
  pandas, numpy, scipy, statsmodels, matplotlib, seaborn, openpyxl

  Install: pip install pandas numpy scipy statsmodels matplotlib seaborn openpyxl
================================================================================
"""

# ──────────────────────────────────────────────────────────────────────────────
# STANDARD LIBRARY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import os  # File path construction and existence checks

# ──────────────────────────────────────────────────────────────────────────────
# THIRD-PARTY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt                        # Figure and axis creation
import numpy as np                                     # Numerical operations, log transforms
import pandas as pd                                    # Data loading and manipulation
import scipy.stats as stats                            # Fisher's exact test
import seaborn as sns                                  # FacetGrid volcano plots
from scipy.stats import chi2_contingency               # Available but not used (kept for reference)
from statsmodels.stats.multitest import multipletests  # Benjamini-Hochberg FDR correction


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

base_path    = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file   = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_xlsx  = os.path.join(base_path, "12_SNP_Enrichment_Results.xlsx")
output_plot  = os.path.join(base_path, "12_SNP_Volcano_Plots.png")

# The minimum gnomAD NFE allele frequency for a variant to be considered a SNP.
# Consistent with script 11.
SNP_AF_THRESHOLD = 0.01

# FDR significance threshold applied to BH-corrected q-values
FDR_THRESHOLD = 0.05

# Number of top significant SNPs to label per volcano plot panel
N_LABELS = 5


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
        The actual column name as it appears in the DataFrame, or None if
        not found.
    """
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None


# ══════════════════════════════════════════════════════════════════════════════
# VOLCANO PLOT DECORATION HELPER
# ══════════════════════════════════════════════════════════════════════════════

def add_volcano_decorations(
    ax,
    group_name: str,
    data: pd.DataFrame,
    x_col: str,
    use_log_scale: bool = True,
):
    """
    Add reference lines, shading, significance labels, and formatting to a
    single volcano plot axis.

    Elements added
    --------------
    Horizontal dashed red line  — p = 0.05 significance threshold
                                  (at y = −log10(0.05) ≈ 1.30)
    Vertical dashed blue line   — OR = 1, the line of no effect
    Blue shading (left half)    — OR < 1 region: SNPs less common in tumour
                                  (potentially protective)
    Red shading (right half)    — OR > 1 region: SNPs more common in tumour
                                  (potentially associated with cancer risk)
    Red horizontal shading      — significant region above the p = 0.05 line
    Annotation labels           — rsID or GENE:HGVSp labels for the top N
                                  most significant SNPs in this group

    Parameters
    ----------
    ax           : matplotlib Axes — the axis to decorate
    group_name   : str             — cohort label (e.g. "Breast") for filtering
    data         : pd.DataFrame    — full results table with 'Analysis_Group',
                                     '-log10_p', 'Label', and x_col columns
    x_col        : str             — column name for the X-axis (Odds_Ratio)
    use_log_scale: bool            — True if X-axis is on a log scale
    """
    # Horizontal threshold line at p = 0.05
    ax.axhline(
        -np.log10(0.05),
        color="red", linestyle="--", alpha=0.6, linewidth=2,
        label="p = 0.05 threshold",
    )

    # Vertical reference line at OR = 1 (no effect)
    ax.axvline(
        1, color="blue", linestyle="--", alpha=0.6, linewidth=2,
        label="OR = 1 (no effect)",
    )

    # Background shading to indicate direction of effect
    xlim = ax.get_xlim()
    ax.axvspan(xlim[0], 1,      alpha=0.05, color="blue",
               label="Protective (OR < 1)")
    ax.axvspan(1,       xlim[1], alpha=0.05, color="red",
               label="Risk factor (OR > 1)")

    # Horizontal shading for the significant region (above the p threshold)
    ax.axhspan(
        -np.log10(FDR_THRESHOLD), ax.get_ylim()[1],
        alpha=0.10, color="red", label=f"Significant (FDR < {FDR_THRESHOLD})",
    )

    # Label top N significant SNPs with their Variant_ID
    group_data = data[data["Analysis_Group"] == group_name]
    labeled_rows = group_data[group_data["Label"] != ""]
    for _, row in labeled_rows.iterrows():
        ax.annotate(
            row["Label"],
            xy=(row[x_col], row["-log10_p"]),
            xytext=(10, 10),
            textcoords="offset points",
            fontsize=8, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7),
            arrowprops=dict(
                arrowstyle="->", connectionstyle="arc3,rad=0",
                color="black", lw=1,
            ),
        )

    ax.set_ylabel("−log10(p-value)", fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="best", fontsize=8)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def run_snp_association_analysis():
    """
    Orchestrate the full SNP enrichment analysis pipeline.

    Execution order
    ---------------
    1.  Load and filter the annotated variant table to SNPs (AF > 1%).
    2.  Standardise tissue and cohort labels; construct Variant IDs.
    3.  For each cohort (Global, Breast, Endometrium):
          a. Count tumour and healthy carrier samples per SNP.
          b. Build the 2×2 contingency table.
          c. Run Fisher's exact test to obtain a p-value.
          d. Calculate the Odds Ratio, applying the Haldane-Anscombe
             correction (+0.5) when any cell count is zero.
    4.  Apply Benjamini-Hochberg FDR correction within each cohort.
    5.  Save all results to a multi-sheet Excel workbook.
    6.  Compute plot-ready columns (−log10 p, log2 OR, significance flags).
    7.  Label the top N most significant SNPs per cohort for annotation.
    8.  Generate and save two volcano plot figures.
    """
    print("=" * 65)
    print("SCRIPT 12 — SNP Statistical Enrichment Analysis")
    print("=" * 65)
    print(f"Input : {input_file}")
    print(f"Output: {output_xlsx}")
    print(f"        {output_plot}\n")

    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        print("       Please update 'base_path' in the CONFIGURATION block.")
        return

    # ── Step 1: Load data ─────────────────────────────────────────────────────
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")
    print(f"Loaded {len(df)} rows from sheet 'Biological_Annotations'")

    # ── Step 2: Identify columns (case-insensitive) ───────────────────────────
    sym_c = find_col(df, "SYMBOL")
    imp_c = find_col(df, "IMPACT")
    con_c = find_col(df, "CONSEQUENCE")
    tis_c = find_col(df, "TISSUE")
    sam_c = find_col(df, "SAMPLE")
    var_c = find_col(df, "Existing_variation")
    hgv_c = find_col(df, "HGVSp")
    coh_c = find_col(df, "COHORT")
    nfe_c = find_col(df, "gnomADe_NFE_AF")

    # ── Step 3: Filter to SNPs and standardise labels ─────────────────────────
    # Retain only common variants (gnomAD NFE AF > 1%), consistent with script 11
    df = df[df[nfe_c] > SNP_AF_THRESHOLD].copy()
    print(f"Variants retained after SNP filter (AF > {SNP_AF_THRESHOLD}): {len(df)}")

    df[tis_c] = (
        df[tis_c].astype(str).str.strip()
        .replace({"Tumor": "Tumour", "Normal": "Healthy", "Control": "Healthy"})
    )
    df[coh_c] = df[coh_c].astype(str).str.strip()

    # ── Step 4: Construct Variant IDs ─────────────────────────────────────────
    # Use rsID where available; fall back to GENE:HGVSp for unregistered variants
    df["rsID"]       = df[var_c].astype(str).str.extract(r"(rs\d+)")
    df["Variant_ID"] = df["rsID"].fillna(
        df[sym_c].astype(str) + ":" + df[hgv_c].astype(str)
    )

    # ══════════════════════════════════════════════════════════════════════════
    # Step 5: Fisher's exact test with Haldane-Anscombe correction
    # ══════════════════════════════════════════════════════════════════════════
    cohort_list = ["Global", "Breast", "Endometrium"]
    all_results = []

    for cohort_name in cohort_list:

        # Subset data for the current cohort
        # "Global" includes all samples; named cohorts filter by cohort label
        if cohort_name == "Global":
            cohort_df = df.copy()
        else:
            cohort_df = df[
                df[coh_c].str.contains(cohort_name, case=False, na=False)
            ].copy()

        if cohort_df.empty:
            print(f"  Cohort '{cohort_name}': no data found, skipping.")
            continue

        # Count total unique samples in each tissue group for this cohort.
        # These are the denominators for frequency calculations.
        n_tumour  = cohort_df[cohort_df[tis_c] == "Tumour"][sam_c].nunique()
        n_healthy = cohort_df[cohort_df[tis_c] == "Healthy"][sam_c].nunique()

        if n_tumour == 0 or n_healthy == 0:
            print(f"  Cohort '{cohort_name}': missing tumour or healthy samples, skipping.")
            continue

        print(f"\n  Cohort '{cohort_name}': "
              f"{n_tumour} tumour samples, {n_healthy} healthy samples, "
              f"{cohort_df['Variant_ID'].nunique()} unique SNPs")

        for variant_id in cohort_df["Variant_ID"].unique():

            # Get all rows for this SNP in this cohort
            var_data = cohort_df[cohort_df["Variant_ID"] == variant_id]

            # Count how many unique samples carry this SNP in each tissue
            count_tumour  = var_data[var_data[tis_c] == "Tumour"][sam_c].nunique()
            count_healthy = var_data[var_data[tis_c] == "Healthy"][sam_c].nunique()

            # Build the 2×2 contingency table:
            #   a = tumour carriers      b = tumour non-carriers
            #   c = healthy carriers     d = healthy non-carriers
            a = count_tumour
            b = max(0, n_tumour  - count_tumour)
            c = count_healthy
            d = max(0, n_healthy - count_healthy)

            # Fisher's exact test — always uses raw counts (no correction)
            _, p_value = stats.fisher_exact([[int(a), int(b)], [int(c), int(d)]])

            # Odds Ratio calculation:
            # If any cell is 0, apply the Haldane-Anscombe correction (+0.5)
            # to avoid division by zero or infinite OR values.
            if a == 0 or b == 0 or c == 0 or d == 0:
                or_val = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                note   = "Corrected (+0.5)"
            else:
                or_val = (a * d) / (b * c) if (b * c) != 0 else np.inf
                note   = "Standard"

            all_results.append({
                "Analysis_Group":     cohort_name,
                "SNP_ID":             variant_id,
                "Symbol":             var_data[sym_c].iloc[0],
                "Impact":             var_data[imp_c].iloc[0],
                "Tumour_Freq_%":      (count_tumour  / n_tumour)  * 100,
                "Healthy_Freq_%":     (count_healthy / n_healthy) * 100,
                "Odds_Ratio":         or_val,
                "P_Value":            p_value,
                "Calculation_Method": note,
            })

    # ══════════════════════════════════════════════════════════════════════════
    # Step 6: Benjamini-Hochberg FDR correction (per cohort)
    # ══════════════════════════════════════════════════════════════════════════
    # FDR correction is applied within each cohort separately, not globally.
    # This is the standard practice: the multiple testing burden is defined
    # by the number of SNPs tested within a single analysis, not across all
    # cohorts combined (which would be an overly conservative correction since
    # the same SNP tested in Breast and Endometrium is an independent question).
    res_df = pd.DataFrame(all_results)

    corrected_dfs = []
    for group in res_df["Analysis_Group"].unique():
        sub = res_df[res_df["Analysis_Group"] == group].copy()

        # multipletests returns: reject array, corrected p-values, alphacSidak, alphaBonf
        _, sub["FDR_P_Value"], _, _ = multipletests(sub["P_Value"], method="fdr_bh")
        corrected_dfs.append(sub)

    final_res = (
        pd.concat(corrected_dfs)
        .sort_values(["Analysis_Group", "P_Value"])
        .reset_index(drop=True)
    )

    # ── Step 7: Save results to Excel ────────────────────────────────────────
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        for group in cohort_list:
            if group in final_res["Analysis_Group"].unique():
                final_res[final_res["Analysis_Group"] == group].to_excel(
                    writer, sheet_name=group, index=False
                )
    print(f"\nResults saved: {output_xlsx}")

    # ══════════════════════════════════════════════════════════════════════════
    # Step 8: Prepare columns for plotting
    # ══════════════════════════════════════════════════════════════════════════

    # −log10(p-value): transforms p-values so that more significant results
    # plot higher on the Y-axis. p = 0 is replaced with a small positive value
    # (1e-10) to avoid log(0) = −∞.
    final_res["-log10_p"] = -np.log10(
        final_res["P_Value"].replace(0, 1e-10)
    )

    # log2(OR): not used on the current X-axis but pre-computed for reference.
    # Zero and infinite OR values are replaced with small/large finite values
    # before the log transform to avoid NaN or ±∞.
    final_res["log2_OR"] = np.log2(
        final_res["Odds_Ratio"].replace([0, np.inf], [0.001, 1000])
    )

    # Flag SNPs that pass the FDR significance threshold
    final_res["Is_Significant"] = final_res["FDR_P_Value"] < FDR_THRESHOLD

    # Prepare label column: by default empty; filled for the top N significant
    # SNPs per cohort so they are annotated on the volcano plot.
    final_res["Label"] = ""
    for group in final_res["Analysis_Group"].unique():
        group_mask = final_res["Analysis_Group"] == group
        sig_mask   = (final_res["FDR_P_Value"] < FDR_THRESHOLD) & group_mask

        if sig_mask.any():
            # Label only the N most significant (smallest FDR q-value) SNPs
            top_indices = final_res[sig_mask].nsmallest(N_LABELS, "FDR_P_Value").index
            final_res.loc[top_indices, "Label"] = final_res.loc[top_indices, "SNP_ID"]

    # ══════════════════════════════════════════════════════════════════════════
    # Step 9: Volcano Plot 1 — OR on log scale
    # ══════════════════════════════════════════════════════════════════════════
    # One panel per cohort (Global / Breast / Endometrium).
    # Dots coloured by VEP Impact category; X-axis on logarithmic scale so
    # OR = 0.1 and OR = 10 are equidistant from OR = 1.
    plt.figure(figsize=(18, 10))
    g = sns.FacetGrid(
        final_res,
        col="Analysis_Group", hue="Impact", palette="Set1",
        height=6, aspect=1.3, despine=False,
    )
    g.map(sns.scatterplot, "Odds_Ratio", "-log10_p",
          s=150, edgecolor="black", alpha=0.8)

    for ax, group_name in zip(g.axes.flat, final_res["Analysis_Group"].unique()):
        ax.set_xscale("log")  # Log scale: symmetric around OR = 1
        add_volcano_decorations(ax, group_name, final_res, "Odds_Ratio",
                                use_log_scale=True)
        ax.set_xlabel("Odds Ratio (log scale)", fontweight="bold")

    g.add_legend(title="VEP Impact", bbox_to_anchor=(1.02, 0.5), loc="center left")
    plt.subplots_adjust(top=0.90, right=0.95)
    g.fig.suptitle(
        "SNP Association Analysis: Tumour vs. Healthy\n"
        "Fisher's Exact Test with Benjamini-Hochberg FDR Correction",
        fontsize=16, fontweight="bold",
    )
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Volcano plot saved: {output_plot}")

    # ══════════════════════════════════════════════════════════════════════════
    # Step 10: Volcano Plot 2 — second version (currently identical to plot 1;
    # X-axis uses OR on log scale, reserved for log2(OR) if preferred)
    # ══════════════════════════════════════════════════════════════════════════
    output_plot_log2 = output_plot.replace(".png", "_log2OR.png")
    plt.figure(figsize=(18, 10))
    g2 = sns.FacetGrid(
        final_res,
        col="Analysis_Group", hue="Impact", palette="Set1",
        height=6, aspect=1.3, despine=False,
    )
    g2.map(sns.scatterplot, "Odds_Ratio", "-log10_p",
           s=150, edgecolor="black", alpha=0.8)

    for ax, group_name in zip(g2.axes.flat, final_res["Analysis_Group"].unique()):
        ax.set_xscale("log")
        add_volcano_decorations(ax, group_name, final_res, "Odds_Ratio",
                                use_log_scale=True)
        ax.set_xlabel("Odds Ratio (log scale)", fontweight="bold")

    g2.add_legend(title="VEP Impact", bbox_to_anchor=(1.02, 0.5), loc="center left")
    plt.subplots_adjust(top=0.90, right=0.95)
    g2.fig.suptitle(
        "SNP Association Analysis: Tumour vs. Healthy\n"
        "Fisher's Exact Test with Benjamini-Hochberg FDR Correction",
        fontsize=16, fontweight="bold",
    )
    plt.savefig(output_plot_log2, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Volcano plot (log2OR variant) saved: {output_plot_log2}")

    # ── Completion summary ────────────────────────────────────────────────────
    n_sig = final_res["Is_Significant"].sum()
    print(f"\n{'=' * 65}")
    print("Script 12 complete.")
    print(f"  Total SNP × cohort tests run : {len(final_res)}")
    print(f"  Significant after FDR correction: {n_sig} "
          f"(FDR < {FDR_THRESHOLD})")
    print(f"{'=' * 65}")


# ── Script entry point ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    run_snp_association_analysis()
