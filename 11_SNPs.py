#!/usr/bin/env python3
"""
================================================================================
Script 11 — SNP Identification and Population Frequency Analysis
================================================================================
Pipeline    : Step 11 of 17 — runs AFTER variant statistics (10) and BEFORE
              statistical enrichment analysis (12)

--------------------------------------------------------------------------------
PURPOSE
--------------------------------------------------------------------------------
This script identifies Single Nucleotide Polymorphisms (SNPs) within the
GSDMB locus and compares their observed frequencies in this study cohort
against published population reference frequencies from gnomAD.

A SNP, in the context used here, is a common genetic variant present at a
frequency greater than 1% in the general population. This distinguishes SNPs
from rare pathogenic mutations, which typically occur at much lower frequencies.
SNPs are important in this study because:

  1. They help characterise the baseline genetic diversity of the GSDMB locus
     across the cancer cohorts studied.
  2. SNPs shared with known population databases (gnomAD) confirm that the
     detected variants are genuine germline polymorphisms rather than
     sequencing artefacts.
  3. Differences in SNP frequencies between tumour and healthy tissue, or
     between cancer types, may point to population stratification effects or
     locus-specific selection pressures.

--------------------------------------------------------------------------------
WHAT IS gnomAD NFE FREQUENCY?
--------------------------------------------------------------------------------
gnomAD (Genome Aggregation Database) is the largest publicly available
reference database of human genetic variation, aggregating sequencing data
from over 125,000 exomes and 15,000 genomes across diverse populations
(Karczewski et al., 2020).

The "NFE" subpopulation refers to Non-Finnish Europeans — a large, relatively
homogeneous reference group commonly used as a comparator in European-ancestry
cancer cohort studies.

The gnomAD NFE allele frequency (gnomADe_NFE_AF) for a variant is the
proportion of NFE chromosomes in the gnomAD database that carry that variant.
A value of 0.01 (1%) is the conventional threshold separating:
  - Common variants / SNPs  (AF ≥ 0.01) — the focus of this script
  - Rare variants           (AF < 0.01) — analysed elsewhere in the pipeline

Comparing study frequency (how often we observe a SNP in our cohort) to
gnomAD NFE frequency (the population baseline) reveals whether our cohort
is enriched or depleted for that variant relative to the general population.

--------------------------------------------------------------------------------
POPULATION BENCHMARKING PLOT (Identity Plot)
--------------------------------------------------------------------------------
The "identity plot" (scatter plot) for each cohort places each SNP as a dot:
  - X-axis: gnomAD NFE population frequency (the expected value)
  - Y-axis: observed frequency in our study cohort (the measured value)
  - Diagonal reference line: the line of identity (y = x × 100), representing
    perfect agreement between our cohort and the population baseline

Points falling on or near the diagonal indicate SNPs that behave as expected
for a random sample of the European population. Points above the diagonal
indicate variants more common in our cohort than in the general population —
potentially of biological or clinical interest.

--------------------------------------------------------------------------------
VARIANT ID CONSTRUCTION
--------------------------------------------------------------------------------
Where a variant has a known rsID (the standard identifier format for SNPs in
public databases, e.g. rs2305480), that rsID is used as the Variant_ID for
clarity and cross-referencing with external databases.

Where no rsID is available (novel or unregistered variants), a fallback
identifier is constructed as "GENE:HGVSp" (e.g. "GSDMB:p.Arg12Gln"),
providing a human-readable unique identifier from information already in
the table.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        Required columns: gnomADe_NFE_AF, Existing_variation, SYMBOL, HGVSp,
                          Cohort, Tissue, Sample, Consequence, IMPACT

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------
  11_Master_Unique_SNP_Summary.xlsx
    One row per unique SNP, with observed frequency columns for each
    cohort × tissue combination, gnomAD NFE reference frequency,
    consequence type, and impact category.

  11_BarPlot_SNPs_Per_Gene.png
    Bar chart showing the number of unique SNPs identified per gene.

  11_IdentityPlot_<Cohort>.png  (one per cohort)
    Scatter plot benchmarking observed SNP frequencies against gnomAD NFE
    population baseline, coloured by tissue type.

  11_SNP_Consequences_Distribution.png
    Horizontal bar chart of the top 10 most common SNP consequence types.

--------------------------------------------------------------------------------
REFERENCES
--------------------------------------------------------------------------------
  Karczewski KJ, et al. (2020). The mutational constraint spectrum quantified
    from variation in 141,456 humans. Nature, 581:434–443.
    https://doi.org/10.1038/s41586-020-2308-7

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python 11_SNPs.py

  (Paths are configured in the CONFIGURATION block below.)

--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------
  Python ≥ 3.8
  pandas, matplotlib, seaborn, numpy, openpyxl

  Install: pip install pandas matplotlib seaborn numpy openpyxl
================================================================================
"""

# ──────────────────────────────────────────────────────────────────────────────
# STANDARD LIBRARY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import os  # File path construction and existence checks

# ──────────────────────────────────────────────────────────────────────────────
# THIRD-PARTY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt  # Figure creation and formatting
import numpy as np               # Numerical operations (available for downstream use)
import pandas as pd              # Data loading, grouping, pivot tables
import seaborn as sns            # Statistical bar and scatter plots


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
# Update base_path to match your working directory before running.

base_path  = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_dir = base_path

# The minimum population allele frequency threshold for a variant to be
# classified as a SNP in this analysis. 0.01 = 1%, the conventional boundary
# between common variants (SNPs) and rare variants.
SNP_AF_THRESHOLD = 0.01


# ══════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def identify_snps_pipeline():
    """
    Identify SNPs, calculate cohort frequencies, benchmark against gnomAD,
    and generate summary tables and figures.

    Execution steps
    ---------------
    1.  Load the annotated variant table.
    2.  Filter to common variants: gnomADe_NFE_AF > 1% (SNP threshold).
    3.  Construct human-readable Variant IDs (rsID where available, else
        GENE:HGVSp as fallback).
    4.  Standardise cohort and tissue labels.
    5.  Calculate observed carrier frequency for each SNP within each
        cohort × tissue group.
    6.  Build a pivot table with one row per SNP and frequency columns for
        each cohort × tissue combination.
    7.  Save the pivot table as an Excel workbook.
    8.  Generate three figure types:
          a. Bar chart: unique SNP count per gene
          b. Identity plot per cohort: observed vs gnomAD frequency
          c. Consequence distribution bar chart
    """
    print("=" * 65)
    print("SCRIPT 11 — SNP Identification and Population Frequency Analysis")
    print("=" * 65)
    print(f"Input : {input_file}")
    print(f"Output: {output_dir}\n")

    # ── Guard: check input file exists ───────────────────────────────────────
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        print("       Please update 'base_path' in the CONFIGURATION block.")
        return

    # ── Step 1: Load annotated variant table ─────────────────────────────────
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")
    print(f"Loaded {len(df)} rows from sheet 'Biological_Annotations'")

    # ── Step 2: Filter to SNPs (gnomAD NFE AF > 1%) ───────────────────────────
    # Only variants with a population allele frequency above the SNP threshold
    # are retained. Variants below this threshold are rare and handled
    # separately elsewhere in the pipeline.
    df_snps = df[df["gnomADe_NFE_AF"] > SNP_AF_THRESHOLD].copy()
    print(f"Variants passing SNP threshold (AF > {SNP_AF_THRESHOLD}): {len(df_snps)} rows")

    # ── Step 3: Construct Variant IDs ─────────────────────────────────────────
    # Primary ID: extract rsID from the Existing_variation column using a
    # regular expression. rsIDs follow the format "rs" followed by digits
    # (e.g. rs2305480). A single cell may contain multiple IDs separated by
    # commas or semicolons; the regex captures only the first rsID found.
    df_snps["rsID"] = (
        df_snps["Existing_variation"]
        .astype(str)
        .str.extract(r"(rs\d+)")
    )

    # Fallback ID: for variants without a registered rsID, construct a
    # descriptive identifier from the gene symbol and protein-level change.
    # Example: "GSDMB:p.Arg12Gln"
    df_snps["Variant_ID"] = df_snps["rsID"].fillna(
        df_snps["SYMBOL"].astype(str) + ":" + df_snps["HGVSp"].astype(str)
    )

    # ── Step 4: Standardise cohort and tissue labels ──────────────────────────
    # Normalise any alternative spellings to the two canonical tissue labels
    # used throughout this pipeline: "Tumour" and "Healthy".
    df_snps["Cohort"] = df_snps["Cohort"].astype(str).str.strip()
    df_snps["Tissue"] = (
        df_snps["Tissue"]
        .astype(str)
        .str.strip()
        .replace({"Tumor": "Tumour", "Normal": "Healthy", "Control": "Healthy"})
    )

    # ── Step 5: Calculate observed carrier frequency per group ────────────────
    # For each cohort × tissue combination, count the total number of unique
    # samples. This denominator is needed to express SNP counts as percentages.
    total_samples_dict = (
        df_snps
        .groupby(["Cohort", "Tissue"])["Sample"]
        .nunique()
        .to_dict()
    )

    # Aggregate: for each (SNP × cohort × tissue) combination, count how many
    # unique samples carry that SNP. Including Consequence and IMPACT in the
    # groupby key ensures these annotations are preserved in the output table.
    summary = (
        df_snps
        .groupby([
            "Variant_ID", "Cohort", "Tissue",
            "SYMBOL", "Consequence", "IMPACT", "gnomADe_NFE_AF",
        ])
        .agg(Carrier_Count=("Sample", "nunique"))
        .reset_index()
    )

    # Add the total sample count for each group and compute frequency
    summary["Total_Group_Samples"] = summary.apply(
        lambda row: total_samples_dict.get((row["Cohort"], row["Tissue"]), 0),
        axis=1,
    )
    summary["Frequency_%"] = (
        summary["Carrier_Count"] / summary["Total_Group_Samples"]
    ) * 100

    print(f"Unique SNPs identified: "
          f"{summary['Variant_ID'].nunique()} across "
          f"{summary['SYMBOL'].nunique()} genes")

    # ── Step 6: Build pivot table ─────────────────────────────────────────────
    # Reshape the summary table so that each SNP occupies one row, and the
    # observed frequency in each cohort × tissue group becomes its own column.
    # This format is more readable and suitable for the Excel output.
    master_pivot = summary.pivot_table(
        index=["Variant_ID", "SYMBOL", "Consequence", "IMPACT", "gnomADe_NFE_AF"],
        columns=["Cohort", "Tissue"],
        values="Frequency_%",
        aggfunc="first",
    )

    # Flatten the multi-level column index into descriptive string labels
    # e.g. ("Breast", "Tumour") → "Breast_Tumour_Frequency_%"
    master_pivot.columns = [
        f"{col[0]}_{col[1]}_Frequency_%" for col in master_pivot.columns
    ]

    # Fill missing frequency values with 0.
    # A missing value means no sample in that group carried the SNP, so 0% is
    # the correct value (not "unknown").
    master_pivot = master_pivot.fillna(0).reset_index()

    # Rename the gnomAD column to a cleaner label for the Excel output
    master_pivot = master_pivot.rename(columns={"gnomADe_NFE_AF": "gnomAD_NFE_AF"})

    # Reorder columns to place Consequence immediately after SYMBOL
    # for logical reading order in the Excel output
    cols = master_pivot.columns.tolist()
    if "SYMBOL" in cols and "Consequence" in cols:
        symbol_idx = cols.index("SYMBOL")
        cols.remove("Consequence")
        cols.insert(symbol_idx + 1, "Consequence")
        master_pivot = master_pivot[cols]

    # ── Step 7: Save Excel summary ────────────────────────────────────────────
    excel_path = os.path.join(output_dir, "11_Master_Unique_SNP_Summary.xlsx")
    master_pivot.to_excel(excel_path, index=False)
    print(f"SNP summary table saved: {excel_path}")

    # ══════════════════════════════════════════════════════════════════════════
    # FIGURE A: Unique SNPs per gene (bar chart)
    # ══════════════════════════════════════════════════════════════════════════
    # Counts distinct Variant_IDs per gene symbol to show which genes in the
    # locus harbour the most common polymorphisms.
    gene_counts = (
        master_pivot
        .groupby("SYMBOL")["Variant_ID"]
        .nunique()
        .sort_values(ascending=False)
        .reset_index()
    )

    plt.figure(figsize=(10, 6))
    sns.barplot(data=gene_counts, x="SYMBOL", y="Variant_ID", palette="viridis")
    plt.title("Total Unique SNPs Identified per Gene", fontsize=14, fontweight="bold")
    plt.xlabel("Gene Symbol")
    plt.ylabel("Unique SNP Count")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    fig_a_path = os.path.join(output_dir, "11_BarPlot_SNPs_Per_Gene.png")
    plt.savefig(fig_a_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {fig_a_path}")

    # ══════════════════════════════════════════════════════════════════════════
    # FIGURE B: Identity plots — observed frequency vs gnomAD NFE baseline
    # (one plot per cohort)
    # ══════════════════════════════════════════════════════════════════════════
    # For each cancer cohort, plot observed study frequency (Y) against
    # gnomAD NFE reference frequency (X). The diagonal line represents the
    # line of identity — where our cohort perfectly matches the population
    # baseline. Deviations from the diagonal indicate enrichment or depletion
    # relative to the European reference population.
    for cohort in summary["Cohort"].unique():
        df_cohort = summary[summary["Cohort"] == cohort]

        plt.figure(figsize=(10, 8))
        sns.scatterplot(
            data    = df_cohort,
            x       = "gnomADe_NFE_AF",
            y       = "Frequency_%",
            hue     = "Tissue",
            s       = 130,
            palette = {"Tumour": "#e74c3c", "Healthy": "#3498db"},
            alpha   = 0.7,
        )

        # Line of identity: if our study frequency matched gnomAD exactly,
        # all points would fall on this line (gnomAD AF × 100 = study %)
        plt.plot(
            [0, 1], [0, 100],
            linestyle="--", color="grey", alpha=0.4,
            label="Line of identity (gnomAD baseline)",
        )

        plt.title(
            f"SNP Frequency Benchmarking: {cohort} Cohort\n"
            f"(Observed study frequency vs. gnomAD NFE population)",
            fontsize=14, fontweight="bold",
        )
        plt.xlabel("gnomAD NFE Population Frequency")
        plt.ylabel("Observed Study Frequency (%)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()

        fig_b_path = os.path.join(output_dir, f"11_IdentityPlot_{cohort}.png")
        plt.savefig(fig_b_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Figure saved: {fig_b_path}")

    # ══════════════════════════════════════════════════════════════════════════
    # FIGURE C: SNP consequence type distribution (horizontal bar chart)
    # ══════════════════════════════════════════════════════════════════════════
    # Shows the top 10 most common VEP consequence types among the identified
    # SNPs, providing a quick overview of the molecular nature of common
    # variation at this locus (e.g. predominantly intronic, synonymous, etc.)
    consequence_counts = master_pivot["Consequence"].value_counts().head(10)

    plt.figure(figsize=(12, 6))
    sns.barplot(
        x       = consequence_counts.values,
        y       = consequence_counts.index,
        palette = "mako",
    )
    plt.title("Top 10 Most Common SNP Consequence Types",
              fontsize=14, fontweight="bold")
    plt.xlabel("Number of Unique SNPs")
    plt.ylabel("VEP Consequence Type")
    plt.tight_layout()
    fig_c_path = os.path.join(output_dir, "11_SNP_Consequences_Distribution.png")
    plt.savefig(fig_c_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {fig_c_path}")

    # ── Completion summary ────────────────────────────────────────────────────
    print(f"\n{'=' * 65}")
    print("Script 11 complete.")
    print(f"  Unique SNPs identified : {len(master_pivot)}")
    print(f"  Consequence types found: {master_pivot['Consequence'].nunique()}")
    print(f"  All outputs saved to   : {output_dir}")
    print(f"{'=' * 65}")


# ── Script entry point ─────────────────────────────────────────────────────────
# Ensures identify_snps_pipeline() is only called when this script is run
# directly, not when imported as a module by another script.
if __name__ == "__main__":
    identify_snps_pipeline()
