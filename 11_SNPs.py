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
This script identifies which of the detected variants are common single
nucleotide polymorphisms (SNPs) — genetic variants that are present at
meaningful frequency in the general population — and calculates how often
each SNP is observed (as a carrier) in each cohort × tissue group in this
study.

The results serve two purposes:
  1. Descriptive: characterise the SNP landscape at the GSDMB locus across
     breast and endometrial cancer cohorts compared to healthy controls.
  2. Preparatory: produce the SNP summary table used by scripts 12, 13, and
     14 for statistical enrichment analysis, permutation testing, and the
     interactive dashboard.

--------------------------------------------------------------------------------
WHAT IS A SNP?
--------------------------------------------------------------------------------
A Single Nucleotide Polymorphism (SNP) is a position in the genome where a
single DNA base differs between individuals in a population. By convention, a
variant is classified as a SNP when it is present in ≥ 1% of the population
(allele frequency ≥ 0.01). Variants below this threshold are considered rare
mutations rather than common polymorphisms.

This 1% threshold is applied here using gnomAD NFE allele frequencies (see
below). A higher allele frequency means the variant is more common in the
general population and is therefore more likely to represent a neutral
polymorphism rather than a disease-causing rare mutation — though population
frequency alone does not determine pathogenicity.

--------------------------------------------------------------------------------
GNOMAD AND THE NFE SUBPOPULATION
--------------------------------------------------------------------------------
gnomAD (Genome Aggregation Database) is a large-scale public database
aggregating exome and genome sequencing data from tens of thousands of
individuals. It provides allele frequency estimates for virtually every known
variant, stratified by continental ancestry group.

NFE (Non-Finnish European) is the ancestry subpopulation used in this analysis.
The breast and endometrial cancer cohorts in this study are predominantly of
European ancestry, so comparing to the NFE subpopulation provides the most
relevant baseline. Using a matched ancestry reference avoids confounding by
known allele frequency differences between continental populations.

gnomAD provides two separate databases:
  gnomADe — Exome database (protein-coding regions only; largest sample size)
  gnomADg — Genome database (whole genome; includes non-coding regions)

--------------------------------------------------------------------------------
COMBINED NFE AF: EXOME + GENOME FALLBACK LOGIC
--------------------------------------------------------------------------------
This script builds a combined allele frequency column using both gnomAD sources:

  Priority 1: gnomADe_NFE_AF (exome database)
    Use if available. The exome database has the largest sample size for
    coding variants and is the primary reference for protein-coding SNPs.

  Priority 2: gnomADg_NFE_AF (genome database)
    Use if exome frequency is missing (NaN). The genome database covers
    introns, UTRs, and intergenic regions that are absent from exome capture.
    Without this fallback, common non-coding variants would be incorrectly
    excluded because their exome frequency is blank — not because they are rare.

  Missing: if both sources are NaN, the variant is excluded.

This two-source strategy ensures that non-coding SNPs (e.g. intronic variants,
UTR variants) are retained when the genome database has a frequency entry but
the exome database does not. The source used for each variant is recorded in
the gnomAD_NFE_Source column of the output for full transparency.

This is an update from earlier pipeline versions that used only
gnomADe_NFE_AF, which caused some genuine non-coding SNPs to be missed.

--------------------------------------------------------------------------------
CARRIER FREQUENCY CALCULATION
--------------------------------------------------------------------------------
For each SNP × cohort × tissue group, the script counts:
  - Carrier_Count       : unique samples carrying at least one copy of the SNP
  - Total_Group_Samples : total unique samples in that cohort × tissue group
  - Frequency_%         : Carrier_Count / Total_Group_Samples × 100

This is CARRIER frequency (presence/absence per sample), not allele frequency.
A sample is counted as a carrier regardless of whether it carries one copy
(heterozygous) or two copies (homozygous) of the SNP.

The frequency is calculated per cohort × tissue group so that, for example,
the frequency in "Breast Tumour" can be compared with "Breast Healthy" and
with the gnomAD population reference in scripts 12–13.

--------------------------------------------------------------------------------
PIVOT TABLE OUTPUT
--------------------------------------------------------------------------------
The long-format summary (one row per SNP × group) is reshaped into a wide-
format pivot table where each row is one unique SNP and each column is the
frequency in one group (e.g. "Breast_Tumour_Frequency_%"). Gaps are filled
with 0 rather than NaN because an undetected SNP has a frequency of 0%,
not "unknown". This wide-format table is the input to scripts 12, 13, and 14.

--------------------------------------------------------------------------------
FOUR OUTPUT FIGURES
--------------------------------------------------------------------------------
Figure 1 — Bar plot: unique SNPs per gene  (11_BarPlot_SNPs_Per_Gene.png)
  Count of unique SNPs per gene. Shows which gene in the locus contributes
  the most common polymorphisms.

Figure 2 — Identity/benchmarking scatter  (11_IdentityPlot_{Cohort}.png)
  One figure per cancer cohort. X = gnomAD NFE population frequency;
  Y = observed carrier frequency in this study (%).
  Diagonal reference line = y = x × 100 (perfect agreement with population).
  Points above the line are enriched; points below are depleted.
  Colour = tissue (tumour/healthy); shape = gnomAD annotation source.

Figure 3 — SNP consequence distribution  (11_SNP_Consequences_Distribution.png)
  Top 10 VEP consequence types among identified SNPs.

Figure 4 — gnomAD source QC bar chart  (11_SNP_gnomAD_Source_Distribution.png)
  How many SNPs used the exome vs. genome source. A substantial "Genome_NFE"
  bar confirms the fallback logic recovered variants that would otherwise
  have been discarded.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        Required columns: Existing_variation, SYMBOL, HGVSp, Cohort, Tissue,
                          Sample, Consequence, IMPACT,
                          gnomADe_NFE_AF, gnomADg_NFE_AF

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------
  11_Master_Unique_SNP_Summary.xlsx
    Wide-format pivot: one row per unique SNP; columns include Variant_ID,
    SYMBOL, Consequence, IMPACT, gnomAD_NFE_AF, gnomAD_NFE_Source, and one
    Frequency_% column per cohort × tissue group.

  11_BarPlot_SNPs_Per_Gene.png
  11_IdentityPlot_{Cohort}.png
  11_SNP_Consequences_Distribution.png
  11_SNP_gnomAD_Source_Distribution.png

--------------------------------------------------------------------------------
REFERENCES
--------------------------------------------------------------------------------
  Karczewski KJ et al. (2020). The mutational constraint spectrum quantified
    from variation in 141,456 humans. Nature, 581:434–443.
    (gnomAD v2 exome database)

  Chen S et al. (2024). A genomic mutational constraint map using variation in
    76,156 human genomes. Nature, 625:92–100.
    (gnomAD v4 genome database)

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python 11_SNPs.py

  (Paths are configured in the CONFIGURATION block below.)

--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------
  Python ≥ 3.8
  pandas, numpy, matplotlib, seaborn, openpyxl

  Install: pip install pandas numpy matplotlib seaborn openpyxl
================================================================================
"""

# ──────────────────────────────────────────────────────────────────────────────
# STANDARD LIBRARY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import os  # File path construction and existence checks

# ──────────────────────────────────────────────────────────────────────────────
# THIRD-PARTY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt  # Figure creation and saving
import numpy as np               # np.where() for conditional column assignment
import pandas as pd              # Data loading, groupby, pivot
import seaborn as sns            # Statistical plot styling


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

base_path  = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_dir = base_path

# Minimum combined gnomAD NFE allele frequency to classify a variant as a SNP.
# The conventional population genetics threshold: variants present in ≥ 1% of
# the NFE population are considered common polymorphisms.
SNP_AF_THRESHOLD = 0.01

# Top N consequence types to display in the consequence distribution bar chart
TOP_N_CONSEQUENCES = 10


# ══════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def identify_snps_pipeline():
    """
    Identify common SNPs and calculate per-cohort carrier frequencies.

    Execution order
    ---------------
    1.  Load the annotated variant table and validate required columns.
    2.  Build the combined gnomAD NFE AF column (exome first, genome fallback).
    3.  Record the annotation source (Exome_NFE / Genome_NFE / Missing) for
        each variant for downstream transparency.
    4.  Filter to SNPs: variants with combined NFE AF > SNP_AF_THRESHOLD.
    5.  Construct Variant IDs (rsID where available, else GENE:HGVSp).
    6.  Standardise tissue and cohort labels.
    7.  Calculate per-SNP carrier counts and carrier frequencies for each
        cohort × tissue group.
    8.  Pivot to wide format and save the master SNP summary to Excel.
    9.  Generate four figures:
          Figure 1 — unique SNP count per gene (bar chart)
          Figure 2 — benchmarking scatter per cohort (identity plots)
          Figure 3 — top consequence types (horizontal bar chart)
          Figure 4 — gnomAD annotation source QC (bar chart)
    """
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        print("       Please update 'base_path' in the CONFIGURATION block.")
        return

    print("=" * 65)
    print("SCRIPT 11 — SNP Identification and Population Frequency Analysis")
    print("=" * 65)

    # ── Step 1: Load and validate ─────────────────────────────────────────────
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")
    print(f"Loaded {len(df)} rows from sheet 'Biological_Annotations'")

    required_cols = [
        "Existing_variation", "SYMBOL", "HGVSp", "Cohort", "Tissue",
        "Sample", "Consequence", "IMPACT",
        "gnomADe_NFE_AF",   # gnomAD exome NFE allele frequency (primary source)
        "gnomADg_NFE_AF",   # gnomAD genome NFE allele frequency (fallback source)
    ]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}")
        print("       Verify column names in the input Excel file.")
        return

    # ══════════════════════════════════════════════════════════════════════════
    # Step 2: Build the combined gnomAD NFE AF column
    # ══════════════════════════════════════════════════════════════════════════
    # pandas combine_first() fills NaN values in the primary series with the
    # corresponding value from the secondary series. This implements the rule:
    #   "use exome AF if available; otherwise use genome AF".
    # The result is never NaN when at least one source has a value.
    df["gnomAD_NFE_AF_combined"] = df["gnomADe_NFE_AF"].combine_first(
        df["gnomADg_NFE_AF"]
    )

    # ── Step 3: Record annotation source per variant ──────────────────────────
    # np.where() evaluates conditions in order, assigning the first matching
    # label. This creates a human-readable audit trail carried through to the
    # output Excel file so that reviewers can verify the source logic.
    df["gnomAD_NFE_Source"] = np.where(
        df["gnomADe_NFE_AF"].notna(),   "Exome_NFE",
        np.where(
            df["gnomADg_NFE_AF"].notna(), "Genome_NFE",
            "Missing"   # Both sources absent: variant will fail the AF filter below
        )
    )

    source_summary = df["gnomAD_NFE_Source"].value_counts()
    print(f"\n  gnomAD NFE annotation source breakdown (all variants):")
    for src, n in source_summary.items():
        print(f"    {src}: {n:,} rows")

    # ══════════════════════════════════════════════════════════════════════════
    # Step 4: Filter to SNPs (combined NFE AF > 1%)
    # ══════════════════════════════════════════════════════════════════════════
    df_snps = df[df["gnomAD_NFE_AF_combined"] > SNP_AF_THRESHOLD].copy()

    if df_snps.empty:
        print(f"\nWARNING: No variants passed the combined gnomAD NFE AF > "
              f"{SNP_AF_THRESHOLD} filter.")
        print("         Check that gnomADe_NFE_AF / gnomADg_NFE_AF columns "
              "are populated in the input file.")
        return

    print(f"\n  Variants retained after SNP filter (combined AF > {SNP_AF_THRESHOLD}): "
          f"{len(df_snps):,}")

    # ── Step 5: Construct Variant IDs ─────────────────────────────────────────
    # rsID is the standardised identifier for known SNPs in public databases.
    # The regex extracts the first "rs" accession from the Existing_variation
    # field, which may contain multiple pipe- or semicolon-delimited IDs
    # (e.g. "rs12345;COSV98765"). For variants without an rsID (e.g. novel
    # variants not yet in dbSNP), the fallback GENE:HGVSp notation is used.
    df_snps["rsID"] = (
        df_snps["Existing_variation"].astype(str)
        .str.extract(r"(rs\d+)")
    )
    df_snps["Variant_ID"] = df_snps["rsID"].fillna(
        df_snps["SYMBOL"].astype(str) + ":" + df_snps["HGVSp"].astype(str)
    )

    # ── Step 6: Standardise labels ────────────────────────────────────────────
    df_snps["Cohort"] = df_snps["Cohort"].astype(str).str.strip()
    df_snps["Tissue"] = (
        df_snps["Tissue"].astype(str).str.strip()
        .replace({"Tumor": "Tumour", "Normal": "Healthy", "Control": "Healthy"})
    )

    # ══════════════════════════════════════════════════════════════════════════
    # Step 7: Calculate carrier frequencies
    # ══════════════════════════════════════════════════════════════════════════

    # Denominator: total unique samples per cohort × tissue group.
    # Stored as a dict for efficient row-wise lookup.
    total_samples_dict = (
        df_snps.groupby(["Cohort", "Tissue"])["Sample"]
        .nunique()
        .to_dict()
    )

    # Numerator: unique samples carrying each SNP per cohort × tissue group.
    # nunique() on Sample avoids double-counting a sample that appears in
    # multiple rows due to VEP transcript expansion (one row per transcript
    # per variant — see script 10).
    summary = (
        df_snps.groupby([
            "Variant_ID", "Cohort", "Tissue", "SYMBOL",
            "Consequence", "IMPACT",
            "gnomAD_NFE_AF_combined",
            "gnomAD_NFE_Source",
        ])
        .agg({"Sample": "nunique"})
        .reset_index()
        .rename(columns={"Sample": "Carrier_Count"})
    )

    # Look up the group denominator for each row
    summary["Total_Group_Samples"] = summary.apply(
        lambda x: total_samples_dict.get((x["Cohort"], x["Tissue"]), 0), axis=1
    )

    summary["Frequency_%"] = (
        summary["Carrier_Count"] / summary["Total_Group_Samples"]
    ) * 100

    # ══════════════════════════════════════════════════════════════════════════
    # Step 8: Pivot to wide format and save to Excel
    # ══════════════════════════════════════════════════════════════════════════
    # Reshape from long (one row per SNP × group) to wide (one row per SNP).
    # Each cohort × tissue group becomes a separate frequency column.
    master_pivot = summary.pivot_table(
        index=[
            "Variant_ID", "SYMBOL", "Consequence", "IMPACT",
            "gnomAD_NFE_AF_combined", "gnomAD_NFE_Source",
        ],
        columns  = ["Cohort", "Tissue"],
        values   = "Frequency_%",
        aggfunc  = "first",   # One unique value per SNP × group — no aggregation needed
    )

    # Flatten the multi-level column headers produced by pivot_table:
    # ("Breast", "Tumour") → "Breast_Tumour_Frequency_%"
    master_pivot.columns = [
        f"{col[0]}_{col[1]}_Frequency_%" for col in master_pivot.columns.values
    ]

    # Fill NaN with 0: a SNP absent from a group has 0% frequency, not unknown
    master_pivot = master_pivot.fillna(0)
    master_pivot.reset_index(inplace=True)

    # Rename combined AF column to a cleaner output name
    master_pivot.rename(columns={"gnomAD_NFE_AF_combined": "gnomAD_NFE_AF"}, inplace=True)

    # Place metadata columns first, then frequency columns
    preferred_order = [
        "Variant_ID", "SYMBOL", "Consequence", "IMPACT",
        "gnomAD_NFE_AF", "gnomAD_NFE_Source",
    ]
    freq_cols    = [c for c in master_pivot.columns if c not in preferred_order]
    master_pivot = master_pivot[preferred_order + freq_cols]

    output_excel = os.path.join(output_dir, "11_Master_Unique_SNP_Summary.xlsx")
    master_pivot.to_excel(output_excel, index=False)
    print(f"\n  Master SNP summary saved: {output_excel}")

    # ══════════════════════════════════════════════════════════════════════════
    # Step 9: Four output figures
    # ══════════════════════════════════════════════════════════════════════════
    sns.set_style("whitegrid")

    # ── Figure 1: Unique SNPs per gene ────────────────────────────────────────
    gene_counts = (
        master_pivot.groupby("SYMBOL")["Variant_ID"]
        .nunique()
        .sort_values(ascending=False)
        .reset_index()
    )

    plt.figure(figsize=(10, 6))
    sns.barplot(
        data    = gene_counts,
        x       = "SYMBOL",
        y       = "Variant_ID",
        hue     = "SYMBOL",
        palette = "viridis",
        legend  = False,
    )
    plt.title("Total Unique SNPs Identified per Gene",
              fontsize=14, fontweight="bold")
    plt.ylabel("Unique SNP Count")
    plt.xlabel("Gene Symbol")
    plt.xticks(rotation=45)
    plt.tight_layout()
    fig1_path = os.path.join(output_dir, "11_BarPlot_SNPs_Per_Gene.png")
    plt.savefig(fig1_path, dpi=300)
    plt.close()
    print(f"  Figure 1 saved: {fig1_path}")

    # ── Figure 2: Benchmarking / identity scatter (one per cohort) ────────────
    # X = gnomAD NFE population allele frequency
    # Y = observed carrier frequency in this study (%)
    # Diagonal line = population baseline (y = x × 100)
    #   Points above: enriched in our cohort relative to the general population
    #   Points below: depleted in our cohort relative to the general population
    # Marker shape encodes annotation source (exome vs. genome) so any
    # systematic differences introduced by the fallback logic are visible.
    for cohort in summary["Cohort"].dropna().unique():
        plt.figure(figsize=(10, 8))
        df_c = summary[summary["Cohort"] == cohort].copy()

        sns.scatterplot(
            data    = df_c,
            x       = "gnomAD_NFE_AF_combined",
            y       = "Frequency_%",
            hue     = "Tissue",
            style   = "gnomAD_NFE_Source",    # Shape = exome vs. genome source
            s       = 130,
            palette = {"Tumour": "#e74c3c", "Healthy": "#3498db"},
            alpha   = 0.75,
        )

        # Reference line: perfect agreement with population baseline
        # X is a proportion [0, 1]; Y is percentage [0, 100] → slope = 100
        plt.plot([0, 1], [0, 100], "--", color="grey", alpha=0.4,
                 label="Population baseline (y = x × 100)")
        plt.title(f"SNP Frequency Benchmarking: {cohort} Cohort",
                  fontsize=15, fontweight="bold")
        plt.xlabel("Combined gnomAD NFE Allele Frequency (population reference)")
        plt.ylabel("Study Carrier Frequency (%)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()
        fig2_path = os.path.join(output_dir, f"11_IdentityPlot_{cohort}.png")
        plt.savefig(fig2_path, dpi=300)
        plt.close()
        print(f"  Figure 2 saved: {fig2_path}")

    # ── Figure 3: SNP consequence type distribution ───────────────────────────
    # Horizontal bar chart of the top 10 consequence types.
    # Reveals whether the SNP landscape at this locus is dominated by intronic
    # (non-coding), synonymous, or missense (protein-altering) variants.
    plt.figure(figsize=(12, 6))
    consequence_counts = master_pivot["Consequence"].value_counts().head(TOP_N_CONSEQUENCES)
    sns.barplot(
        x       = consequence_counts.values,
        y       = consequence_counts.index,
        hue     = consequence_counts.index,
        palette = "mako",
        legend  = False,
    )
    plt.title(f"Top {TOP_N_CONSEQUENCES} Most Common SNP Consequence Types",
              fontsize=14, fontweight="bold")
    plt.xlabel("Number of Unique SNPs")
    plt.ylabel("VEP Consequence Type")
    plt.tight_layout()
    fig3_path = os.path.join(output_dir, "11_SNP_Consequences_Distribution.png")
    plt.savefig(fig3_path, dpi=300)
    plt.close()
    print(f"  Figure 3 saved: {fig3_path}")

    # ── Figure 4: gnomAD annotation source QC ────────────────────────────────
    # Transparency bar chart showing how many unique SNPs in the final output
    # used the exome database vs. the genome database as their AF source.
    # A substantial "Genome_NFE" bar is the methodological justification for
    # the combined AF approach: it proves that the fallback logic recovered
    # real SNPs that would have been discarded by an exome-only approach.
    plt.figure(figsize=(6, 5))
    source_counts = master_pivot["gnomAD_NFE_Source"].value_counts()
    sns.barplot(
        x       = source_counts.index,
        y       = source_counts.values,
        hue     = source_counts.index,
        palette = "Set2",
        legend  = False,
    )
    plt.title("gnomAD Annotation Source Used for Combined NFE AF\n"
              "(QC — confirms exome + genome fallback logic)",
              fontsize=13, fontweight="bold")
    plt.xlabel("gnomAD Database Source")
    plt.ylabel("Number of Unique SNPs")
    plt.tight_layout()
    fig4_path = os.path.join(output_dir, "11_SNP_gnomAD_Source_Distribution.png")
    plt.savefig(fig4_path, dpi=300)
    plt.close()
    print(f"  Figure 4 saved: {fig4_path}")

    # ── Completion summary ────────────────────────────────────────────────────
    print(f"\n{'=' * 65}")
    print("Script 11 complete.")
    print(f"  Unique SNPs identified       : {len(master_pivot)}")
    print(f"  Unique consequence types     : {master_pivot['Consequence'].nunique()}")
    print(f"  gnomAD NFE AF logic:")
    print(f"    Primary source  : gnomADe_NFE_AF (exome database)")
    print(f"    Fallback source : gnomADg_NFE_AF (genome database)")
    print(f"  Master summary saved to      : {output_excel}")
    print(f"{'=' * 65}")


# ── Script entry point ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    identify_snps_pipeline()
