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
  - Global      : all samples combined, regardless of cancer type
  - Breast      : breast cancer cohort only
  - Endometrium : endometrial cancer cohort only

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

FDR correction is applied within each cohort separately — the multiple testing
burden is defined by the number of SNPs tested within a single analysis, not
across all cohorts combined.

--------------------------------------------------------------------------------
TWO VOLCANO PLOTS: RAW vs FDR-ADJUSTED
--------------------------------------------------------------------------------
This script produces two volcano plots using the same layout but different
Y-axis values:

  Plot 1 — Raw p-value volcano  (12_SNP_Volcano_Plots_RAW.png)
    Y-axis = −log10(raw Fisher p-value)
    Threshold line at −log10(0.05) ≈ 1.30
    Shows all nominally significant SNPs before multiple testing correction.
    Useful for exploratory screening — more sensitive but less specific.

  Plot 2 — FDR-adjusted volcano  (12_SNP_Volcano_Plots_FDR.png)
    Y-axis = −log10(BH-corrected FDR q-value)
    Threshold line at −log10(0.05) ≈ 1.30
    Shows only SNPs that remain significant after correcting for the number
    of tests performed. More stringent — the primary result for reporting.

Presenting both plots allows the reader to assess how many nominally
significant results survive multiple testing correction — an important
transparency measure in genomic association studies.

Within each plot, panels are arranged in the order: Breast | Endometrium |
Global. Both X and Y axes are shared across panels to facilitate direct
visual comparison between cohorts.

--------------------------------------------------------------------------------
VOLCANO PLOT INTERPRETATION
--------------------------------------------------------------------------------
Each dot in the volcano plot represents one SNP:
  - X-axis : Odds Ratio on a log scale — direction and magnitude of association
  - Y-axis : −log10(p-value) — statistical significance (higher = more significant)

Reference elements:
  - Horizontal dashed red line  : significance threshold (p = 0.05)
  - Vertical dashed blue line   : OR = 1 (no association)
  - Blue left shading           : OR < 1 region (SNP depleted in tumour → protective)
  - Red right shading           : OR > 1 region (SNP enriched in tumour → risk)
  - Red top shading             : significant region (above threshold line)

The top 5 most significant SNPs per panel are labelled with their rsID or
GENE:HGVSp identifier.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        Required columns: SYMBOL, IMPACT, TISSUE, SAMPLE,
                          Existing_variation, HGVSp, COHORT, gnomADe_NFE_AF

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------
  12_SNP_Enrichment_Results.xlsx
    Multi-sheet workbook (Global | Breast | Endometrium). Each row is one SNP
    with: carrier counts, frequencies, OR, raw p-value, BH FDR q-value,
    and correction method used.

  12_SNP_Volcano_Plots_RAW.png
    Three-panel volcano plot using raw Fisher p-values (300 dpi).

  12_SNP_Volcano_Plots_FDR.png
    Three-panel volcano plot using BH-adjusted FDR q-values (300 dpi).

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

import matplotlib.pyplot as plt                        # Figure creation
import numpy as np                                     # Numerical operations
import pandas as pd                                    # Data loading and manipulation
import scipy.stats as stats                            # Fisher's exact test
import seaborn as sns                                  # FacetGrid volcano plots
from statsmodels.stats.multitest import multipletests  # Benjamini-Hochberg FDR correction


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

base_path   = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file  = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_xlsx = os.path.join(base_path, "12_SNP_Enrichment_Results.xlsx")

# Two separate volcano plot outputs: one for raw p-values, one for FDR q-values
output_plot_raw = os.path.join(base_path, "12_SNP_Volcano_Plots_RAW.png")
output_plot_fdr = os.path.join(base_path, "12_SNP_Volcano_Plots_FDR.png")

# Minimum gnomAD NFE allele frequency to classify a variant as a SNP.
# Consistent with scripts 11 and 13.
SNP_AF_THRESHOLD = 0.01

# Significance threshold applied to both raw and FDR-adjusted p-values.
SIG_THRESHOLD = 0.05

# Number of top significant SNPs to label per volcano panel.
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
# VOLCANO PLOT HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def add_volcano_decorations(
    ax,
    group_name:  str,
    data:        pd.DataFrame,
    y_col:       str,
    label_col:   str,
    threshold_p: float,
    y_label:     str,
):
    """
    Add reference lines, shading, SNP labels, and axis formatting to one
    volcano plot panel.

    Elements added
    --------------
    Horizontal dashed red line  — significance threshold at −log10(threshold_p)
    Vertical dashed blue line   — OR = 1 (no association)
    Blue left shading           — OR < 1 region (protective direction)
    Red right shading           — OR > 1 region (risk direction)
    Red top shading             — significant region above the threshold line
    SNP labels                  — annotated text boxes for pre-selected top SNPs

    Parameters
    ----------
    ax          : matplotlib Axes — the panel to decorate
    group_name  : str             — cohort name (e.g. "Breast") for data filtering
    data        : pd.DataFrame    — full results table
    y_col       : str             — column name for the Y-axis values
    label_col   : str             — column containing label text (empty = unlabelled)
    threshold_p : float           — p-value threshold (e.g. 0.05)
    y_label     : str             — Y-axis label string
    """
    y_thresh = -np.log10(threshold_p)

    # Significance threshold line (horizontal)
    ax.axhline(
        y_thresh,
        color="red", linestyle="--", alpha=0.7, linewidth=2,
        label=f"threshold = {threshold_p}",
    )

    # Line of no effect: OR = 1 means tumour and healthy frequencies are equal
    ax.axvline(
        1,
        color="blue", linestyle="--", alpha=0.7, linewidth=2,
        label="OR = 1 (no effect)",
    )

    # Background shading indicating direction of association
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.axvspan(xlim[0], 1,       alpha=0.05, color="blue", label="Protective (OR < 1)")
    ax.axvspan(1,       xlim[1], alpha=0.05, color="red",  label="Risk factor (OR > 1)")

    # Horizontal shading for the region above the significance threshold
    ax.axhspan(y_thresh, ylim[1], alpha=0.08, color="red", label="Above threshold")

    # Annotate pre-selected significant SNPs with labelled text boxes.
    # Labels were assigned before plotting (in the main function) to the top
    # N_LABELS most significant SNPs per cohort.
    group_data = data[data["Analysis_Group"] == group_name]
    labelled   = group_data[group_data[label_col] != ""]
    for _, row in labelled.iterrows():
        ax.annotate(
            row[label_col],
            xy=(row["Odds_Ratio_Plot"], row[y_col]),
            xytext=(8, 8),
            textcoords="offset points",
            fontsize=8, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="yellow", alpha=0.7),
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
        )

    ax.set_ylabel(y_label, fontweight="bold")
    ax.set_xlabel("Odds Ratio (log scale)", fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="lower left", fontsize=8, frameon=True)


def make_volcano_plot(
    data:        pd.DataFrame,
    y_col:       str,
    label_col:   str,
    threshold_p: float,
    y_label:     str,
    title:       str,
    output_path: str,
):
    """
    Build and save a three-panel volcano plot (one panel per cohort).

    Layout
    ------
    Panels are ordered: Breast | Endometrium | Global — placing the two
    disease-specific cohorts first and the combined global analysis last.
    Both X and Y axes are shared across panels (sharex=True, sharey=True)
    so that relative effect sizes and significance levels can be visually
    compared directly between cohorts.

    The X-axis uses a log scale so that OR = 0.1 and OR = 10 are equidistant
    from OR = 1, providing a symmetric visual representation of protective and
    risk associations.

    This function is called twice: once for the raw p-value plot and once for
    the FDR-adjusted plot. The y_col and label_col parameters control which
    Y-axis values and SNP labels are used in each case.

    Parameters
    ----------
    data        : pd.DataFrame — full results table
    y_col       : str          — Y-axis column ("neglog10_raw_p" or "neglog10_fdr_p")
    label_col   : str          — label column ("Label_RAW" or "Label_FDR")
    threshold_p : float        — significance threshold for the horizontal line
    y_label     : str          — Y-axis label string
    title       : str          — overall figure title
    output_path : str          — full path for the saved PNG file
    """
    group_order = ["Breast", "Endometrium", "Global"]

    # Restrict to the three expected cohorts and enforce panel ordering
    plot_data = data[data["Analysis_Group"].isin(group_order)].copy()
    plot_data["Analysis_Group"] = pd.Categorical(
        plot_data["Analysis_Group"],
        categories=group_order,
        ordered=True,
    )
    plot_data = plot_data.sort_values("Analysis_Group")

    # FacetGrid creates one subplot per Analysis_Group value.
    # hue=Impact colours dots by VEP impact category (HIGH/MODERATE/MODIFIER/LOW).
    g = sns.FacetGrid(
        plot_data,
        col      = "Analysis_Group",
        hue      = "Impact",
        col_order= group_order,
        palette  = "Set1",
        height   = 5.8,
        aspect   = 1.15,
        despine  = False,
        sharex   = True,   # Shared X-axis: OR scale consistent across cohorts
        sharey   = True,   # Shared Y-axis: significance scale consistent across cohorts
    )

    g.map_dataframe(
        sns.scatterplot,
        x         = "Odds_Ratio_Plot",
        y         = y_col,
        s         = 120,
        edgecolor = "black",
        alpha     = 0.8,
    )

    # Add decorations and log scale to each panel
    for ax, group_name in zip(g.axes.flat, group_order):
        ax.set_xscale("log")  # Log scale: OR = 0.1 and OR = 10 equidistant from OR = 1
        add_volcano_decorations(
            ax          = ax,
            group_name  = group_name,
            data        = plot_data,
            y_col       = y_col,
            label_col   = label_col,
            threshold_p = threshold_p,
            y_label     = y_label,
        )

    # Place the VEP Impact colour legend to the right of the figure
    g.add_legend(title="VEP Impact")
    if g._legend is not None:
        g._legend.set_bbox_to_anchor((1.02, 0.5))
        g._legend._loc = 6  # Centre-left anchor point

    g.fig.subplots_adjust(top=0.82, right=0.86, wspace=0.08)
    g.fig.suptitle(title, fontsize=16, fontweight="bold", y=0.98)

    g.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    print(f"  Volcano plot saved: {output_path}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def run_snp_association_analysis():
    """
    Orchestrate the full SNP enrichment analysis pipeline.

    Execution order
    ---------------
    1.  Load and validate the annotated variant table.
    2.  Filter to SNPs (gnomAD NFE AF > 1%) and standardise tissue labels.
    3.  Remove rows with unexpected tissue labels (anything other than
        "Tumour" or "Healthy" after standardisation).
    4.  Construct Variant IDs (rsID where available, else GENE:HGVSp).
    5.  For each cohort (Global, Breast, Endometrium):
          For each unique SNP:
            a. Build the 2×2 contingency table.
            b. Run Fisher's exact test (raw p-value).
            c. Compute the Odds Ratio with Haldane-Anscombe correction if
               any cell count is zero.
    6.  Apply Benjamini-Hochberg FDR correction within each cohort.
    7.  Compute plot-ready columns (−log10 p for raw and FDR, safe OR for
        log scale, significance flags, SNP labels).
    8.  Save results to a multi-sheet Excel workbook.
    9.  Generate and save the raw p-value volcano plot.
    10. Generate and save the FDR-adjusted volcano plot.
    """
    print("=" * 65)
    print("SCRIPT 12 — SNP Statistical Enrichment Analysis")
    print("=" * 65)
    print(f"Input : {input_file}")
    print(f"Output: {output_xlsx}")
    print(f"        {output_plot_raw}")
    print(f"        {output_plot_fdr}\n")

    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        print("       Please update 'base_path' in the CONFIGURATION block.")
        return

    # ── Step 1: Load data ─────────────────────────────────────────────────────
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")
    print(f"Loaded {len(df)} rows from sheet 'Biological_Annotations'")

    # ── Step 2: Identify required columns ────────────────────────────────────
    sym_c = find_col(df, "SYMBOL")
    imp_c = find_col(df, "IMPACT")
    tis_c = find_col(df, "TISSUE")
    sam_c = find_col(df, "SAMPLE")
    var_c = find_col(df, "Existing_variation")
    hgv_c = find_col(df, "HGVSp")
    coh_c = find_col(df, "COHORT")
    nfe_c = find_col(df, "gnomADe_NFE_AF")

    # Validate that all required columns were found before proceeding.
    # Printing the full mapping helps diagnose which column name is missing
    # rather than crashing with an unhelpful AttributeError later.
    required = [sym_c, imp_c, tis_c, sam_c, var_c, hgv_c, coh_c, nfe_c]
    if any(c is None for c in required):
        print("ERROR: One or more required columns could not be found.")
        print("       Column mapping found:")
        print(f"         SYMBOL={sym_c}, IMPACT={imp_c}, TISSUE={tis_c}, SAMPLE={sam_c}")
        print(f"         Existing_variation={var_c}, HGVSp={hgv_c}, "
              f"COHORT={coh_c}, gnomADe_NFE_AF={nfe_c}")
        return

    # ── Step 3: Filter to SNPs and standardise labels ─────────────────────────
    df = df[df[nfe_c] > SNP_AF_THRESHOLD].copy()
    print(f"Variants after SNP filter (AF > {SNP_AF_THRESHOLD}): {len(df)}")

    # Harmonise tissue label variants to the two canonical labels used
    # throughout this pipeline: "Tumour" and "Healthy"
    df[tis_c] = (
        df[tis_c].astype(str).str.strip()
        .replace({
            "Tumor":   "Tumour",
            "Tumour":  "Tumour",
            "Normal":  "Healthy",
            "Healthy": "Healthy",
            "Control": "Healthy",
        })
    )
    df[coh_c] = df[coh_c].astype(str).str.strip()

    # Remove any rows that did not resolve to one of the two expected labels.
    # This prevents unexpected tissue values from silently distorting counts.
    df = df[df[tis_c].isin(["Tumour", "Healthy"])].copy()

    # ── Step 4: Construct Variant IDs ─────────────────────────────────────────
    df["rsID"]       = df[var_c].astype(str).str.extract(r"(rs\d+)")
    df["Variant_ID"] = df["rsID"].fillna(
        df[sym_c].astype(str) + ":" + df[hgv_c].astype(str)
    )

    # ══════════════════════════════════════════════════════════════════════════
    # Step 5: Fisher's exact test with Haldane-Anscombe OR correction
    # ══════════════════════════════════════════════════════════════════════════
    cohort_list = ["Global", "Breast", "Endometrium"]
    all_results = []

    for cohort_name in cohort_list:

        cohort_df = (
            df.copy() if cohort_name == "Global"
            else df[df[coh_c].str.contains(cohort_name, case=False, na=False)].copy()
        )

        if cohort_df.empty:
            print(f"  Cohort '{cohort_name}': no data, skipping.")
            continue

        n_tumour  = cohort_df[cohort_df[tis_c] == "Tumour"][sam_c].nunique()
        n_healthy = cohort_df[cohort_df[tis_c] == "Healthy"][sam_c].nunique()

        if n_tumour == 0 or n_healthy == 0:
            print(f"  Cohort '{cohort_name}': missing tissue group, skipping.")
            continue

        unique_snps = cohort_df["Variant_ID"].dropna().unique()
        print(f"\n  Cohort '{cohort_name}': "
              f"{n_tumour} tumour samples, {n_healthy} healthy samples, "
              f"{len(unique_snps)} unique SNPs")

        for variant_id in unique_snps:

            var_data      = cohort_df[cohort_df["Variant_ID"] == variant_id]
            count_tumour  = var_data[var_data[tis_c] == "Tumour"][sam_c].nunique()
            count_healthy = var_data[var_data[tis_c] == "Healthy"][sam_c].nunique()

            # 2×2 contingency table cells (cast to int for scipy)
            a = int(count_tumour)
            b = int(max(0, n_tumour  - count_tumour))   # tumour non-carriers
            c = int(count_healthy)
            d = int(max(0, n_healthy - count_healthy))  # healthy non-carriers

            # Fisher's exact test on raw (uncorrected) counts
            _, p_value = stats.fisher_exact([[a, b], [c, d]])

            # Odds Ratio — apply Haldane-Anscombe +0.5 if any cell is zero
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
                "Tumour_Carriers":    count_tumour,
                "Healthy_Carriers":   count_healthy,
                "Tumour_Total":       n_tumour,
                "Healthy_Total":      n_healthy,
                "Tumour_Freq_%":      (count_tumour  / n_tumour)  * 100,
                "Healthy_Freq_%":     (count_healthy / n_healthy) * 100,
                "Odds_Ratio":         or_val,
                "P_Value":            p_value,
                "Calculation_Method": note,
            })

    if not all_results:
        print("No SNP association results generated — check input data and filters.")
        return

    # ══════════════════════════════════════════════════════════════════════════
    # Step 6: Benjamini-Hochberg FDR correction (within each cohort)
    # ══════════════════════════════════════════════════════════════════════════
    res_df = pd.DataFrame(all_results)
    corrected_dfs = []

    for group in res_df["Analysis_Group"].unique():
        sub = res_df[res_df["Analysis_Group"] == group].copy()
        _, sub["FDR_P_Value"], _, _ = multipletests(sub["P_Value"], method="fdr_bh")
        corrected_dfs.append(sub)

    final_res = (
        pd.concat(corrected_dfs, ignore_index=True)
        .sort_values(["Analysis_Group", "P_Value"])
    )

    # ══════════════════════════════════════════════════════════════════════════
    # Step 7: Compute plot-ready derived columns
    # ══════════════════════════════════════════════════════════════════════════

    # Replace exact zero p-values with a very small positive value before the
    # log transform to avoid log(0) = −∞. 1e-300 is effectively zero for any
    # practical interpretation but avoids a numerical error.
    final_res["P_Value_Safe"]     = final_res["P_Value"].replace(0, 1e-300)
    final_res["FDR_P_Value_Safe"] = final_res["FDR_P_Value"].replace(0, 1e-300)

    final_res["neglog10_raw_p"] = -np.log10(final_res["P_Value_Safe"])
    final_res["neglog10_fdr_p"] = -np.log10(final_res["FDR_P_Value_Safe"])

    # Replace zero and infinite OR values with finite substitutes for log-scale
    # plotting. These edge cases arise where the Haldane correction was applied.
    final_res["Odds_Ratio_Plot"] = (
        final_res["Odds_Ratio"]
        .replace(0,      1e-6)   # Near-zero OR → plots at far left of log scale
        .replace(np.inf, 1e6)    # Infinite OR  → plots at far right of log scale
    )

    # Significance flags — separate for raw p and FDR q
    final_res["Raw_Significant"] = final_res["P_Value"]     < SIG_THRESHOLD
    final_res["FDR_Significant"] = final_res["FDR_P_Value"] < SIG_THRESHOLD

    # Pre-assign SNP labels for the top N_LABELS significant variants per cohort.
    # Separate label columns are used for the two volcano plots so that each
    # figure labels only its own relevant significant SNPs.
    final_res["Label_RAW"] = ""
    final_res["Label_FDR"] = ""

    for group in final_res["Analysis_Group"].unique():
        group_mask = final_res["Analysis_Group"] == group

        raw_sig = final_res[group_mask & final_res["Raw_Significant"]]
        if not raw_sig.empty:
            top_raw = raw_sig.nsmallest(N_LABELS, "P_Value").index
            final_res.loc[top_raw, "Label_RAW"] = final_res.loc[top_raw, "SNP_ID"]

        fdr_sig = final_res[group_mask & final_res["FDR_Significant"]]
        if not fdr_sig.empty:
            top_fdr = fdr_sig.nsmallest(N_LABELS, "FDR_P_Value").index
            final_res.loc[top_fdr, "Label_FDR"] = final_res.loc[top_fdr, "SNP_ID"]

    # ── Step 8: Save results to Excel ────────────────────────────────────────
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        for group in cohort_list:
            if group in final_res["Analysis_Group"].unique():
                final_res[final_res["Analysis_Group"] == group].to_excel(
                    writer, sheet_name=group, index=False
                )
    print(f"\nResults saved: {output_xlsx}")

    # ── Step 9: Raw p-value volcano plot ─────────────────────────────────────
    make_volcano_plot(
        data        = final_res,
        y_col       = "neglog10_raw_p",
        label_col   = "Label_RAW",
        threshold_p = SIG_THRESHOLD,
        y_label     = "−log10(raw p-value)",
        title       = ("SNP Association Analysis: Tumour vs. Healthy\n"
                       "Raw Fisher p-values  |  Threshold: p < 0.05"),
        output_path = output_plot_raw,
    )

    # ── Step 10: FDR-adjusted volcano plot ───────────────────────────────────
    make_volcano_plot(
        data        = final_res,
        y_col       = "neglog10_fdr_p",
        label_col   = "Label_FDR",
        threshold_p = SIG_THRESHOLD,
        y_label     = "−log10(FDR-adjusted q-value)",
        title       = ("SNP Association Analysis: Tumour vs. Healthy\n"
                       "Benjamini-Hochberg FDR-adjusted q-values  |  Threshold: q < 0.05"),
        output_path = output_plot_fdr,
    )

    # ── Completion summary ────────────────────────────────────────────────────
    n_raw_sig = final_res["Raw_Significant"].sum()
    n_fdr_sig = final_res["FDR_Significant"].sum()

    print(f"\n{'=' * 65}")
    print("Script 12 complete.")
    print(f"  Total SNP × cohort tests run : {len(final_res)}")
    print(f"  Nominally significant        : {n_raw_sig}  (raw p < {SIG_THRESHOLD})")
    print(f"  After FDR correction         : {n_fdr_sig}  (FDR q < {SIG_THRESHOLD})")
    print(f"{'=' * 65}")


# ── Script entry point ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    run_snp_association_analysis()
