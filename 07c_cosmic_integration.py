#!/usr/bin/env python3
"""
================================================================================
Script 07c — COSMIC Cancer Gene Census Integration
================================================================================
Author      : [Your Name]
Institution : [Your University / Department]
Thesis      : Characterisation of GSDMB genetic variation across cancer cohorts
Pipeline    : Step 07c of 17 — runs AFTER variant QC (07b) and BEFORE
              genomic mapping (08)
Last updated: 2025

--------------------------------------------------------------------------------
PURPOSE
--------------------------------------------------------------------------------
This script compares the detected variants against the COSMIC (Catalogue Of
Somatic Mutations In Cancer) Cancer Gene Census — the most comprehensive
curated reference of genes with established roles in human cancer
(Sondka et al., 2018; https://cancer.sanger.ac.uk/census).

The goal is to contextualise the variants found in this study within the
broader landscape of cancer genetics: which variants fall in known cancer
driver genes? Which positions are recurrently mutated across patients
(mutational hotspots)? Do variants cluster within specific cancer types?

The analysis has four analytical modules:

  1. COSMIC Census Gene Annotation
     ───────────────────────────────
     Each variant is cross-referenced against a curated list of COSMIC Cancer
     Gene Census entries relevant to this study (GSDMB locus and key cancer
     drivers). Each gene in the Census is classified by:
       - Tier: a confidence rating for the gene's involvement in cancer
           Tier 1 = strong evidence from multiple independent studies
           Tier 2 = moderate evidence
           Tier 3 = limited or indirect evidence
       - Role: functional role in cancer biology (oncogene, tumour suppressor,
               pyroptosis mediator, transcription factor, etc.)
     Variants are then assigned a Cancer Relevance score (High / Medium /
     Low / Unknown) based on the combination of gene Tier and predicted
     functional impact (HIGH / MODERATE / LOW / MODIFIER from VEP).

  2. Mutational Hotspot Identification
     ─────────────────────────────────
     A "mutational hotspot" is a genomic position that is independently
     mutated in multiple patients. Recurrence across patients is important
     because random sequencing errors or passenger mutations are unlikely
     to recur at the same position. A variant found in ≥3 independent samples
     at the same position is flagged as a hotspot. The recurrence rate is
     expressed as a percentage of all samples in that tissue group.

  3. Variant Impact by Tissue Type
     ────────────────────────────────
     Variant functional impact categories (HIGH, MODERATE, MODIFIER, LOW)
     are compared between tumour and healthy tissue samples. This replaces
     somatic/germline classification, which is not possible here because the
     study design uses unpaired samples (tumour and healthy tissue come from
     different patients — see note below).

  4. Cancer Type Specificity
     ────────────────────────
     Variants in tumour samples are examined to determine whether they are
     specific to one cancer type (breast or endometrium) or shared across
     cancer types (pan-cancer). Pan-cancer variants in known driver genes are
     of particular biological interest.

--------------------------------------------------------------------------------
CRITICAL NOTE ON STUDY DESIGN: UNPAIRED SAMPLES
--------------------------------------------------------------------------------
Standard somatic mutation calling compares a tumour sample directly to a
matched normal (healthy) sample from the SAME patient (a "matched tumour-
normal pair"). This paired design allows confident identification of mutations
that arose in the tumour and were absent in that patient's germline.

This cohort does NOT have matched pairs. The tumour and healthy tissue samples
come from different individuals. Therefore:
  - We CANNOT classify individual variants as somatic or germline.
  - We CANNOT filter germline polymorphisms using a matched normal.
  - Instead, we compare variant patterns ACROSS the two tissue groups at the
    POPULATION level: are certain variants or positions enriched in tumour
    samples relative to healthy controls across the cohort?

This is a valid and widely used approach in population-level cancer genomics,
but the limitation must be clearly stated when interpreting results.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        One row per variant × transcript (VEP format).
        Required columns: SYMBOL, Pos, TISSUE, SAMPLE, IMPACT,
                          CONSEQUENCE, COHORT, HGVSp, Existing_variation

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------
  07c_COSMIC_Integration_Results.xlsx
    ├── Sheet: COSMIC_Gene_Summary        — per-gene variant counts and roles
    ├── Sheet: All_COSMIC_Variants        — all rows matching Census genes
    ├── Sheet: High_Cancer_Relevance      — variants with High relevance score
    ├── Sheet: Mutational_Hotspots        — recurrent positions (≥3 samples)
    ├── Sheet: Cancer_Specificity         — breast vs endometrium specificity
    └── Sheet: Tissue_Impact_Summary      — impact counts by tissue type

  07c_COSMIC_Analysis_Plots.png
    Six-panel figure (300 dpi) summarising all four analytical modules.

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python 07c_cosmic_integration.py

  (Paths are configured in the CONFIGURATION block below.)

--------------------------------------------------------------------------------
REFERENCES
--------------------------------------------------------------------------------
  Sondka Z, et al. (2018). The COSMIC Cancer Gene Census: describing genetic
    dysfunction across all human cancers. Nat Rev Cancer, 18(11):696–705.
    https://doi.org/10.1038/s41568-018-0060-1

  Tate JG, et al. (2019). COSMIC: the Catalogue Of Somatic Mutations In Cancer.
    Nucleic Acids Res, 47(D1):D941–D947.
    https://doi.org/10.1093/nar/gky1015

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

import os        # File path construction and existence checks
import warnings  # Suppress non-critical library warnings during execution

# ──────────────────────────────────────────────────────────────────────────────
# THIRD-PARTY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt          # Base plotting library
import numpy as np                       # Numerical operations and NaN handling
import pandas as pd                      # DataFrame manipulation and Excel I/O
import seaborn as sns                    # Statistical plot styling
from matplotlib.patches import Patch    # Custom legend elements for plots

warnings.filterwarnings("ignore")  # Suppress openpyxl and seaborn deprecation warnings


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
# All file paths are defined here. Update BASE_PATH to match your working
# directory before running.

BASE_PATH    = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
INPUT_FILE   = os.path.join(BASE_PATH, "GSDMB_Annotated_Report_Fixed.xlsx")
OUTPUT_XLSX  = os.path.join(BASE_PATH, "07c_COSMIC_Integration_Results.xlsx")
OUTPUT_PLOT  = os.path.join(BASE_PATH, "07c_COSMIC_Analysis_Plots.png")


# ══════════════════════════════════════════════════════════════════════════════
# COSMIC CANCER GENE CENSUS — REFERENCE DATA
# ══════════════════════════════════════════════════════════════════════════════
#
# This dictionary contains manually curated entries for genes relevant to the
# GSDMB locus and the cancer types studied here (breast and endometrial cancer).
#
# Each entry reflects the gene's classification in the COSMIC Cancer Gene
# Census (https://cancer.sanger.ac.uk/census) as of the date of this analysis.
#
# Fields:
#   Role     — functional role of the gene in cancer biology
#   Tier     — confidence of cancer gene status:
#                 1 = strong, well-replicated evidence across multiple studies
#                 2 = moderate evidence (fewer studies or indirect mechanisms)
#                 3 = limited or context-dependent evidence
#   Somatic  — whether somatic mutations in this gene are documented in COSMIC
#   Germline — whether germline mutations in this gene confer cancer risk
#
# Gene notes:
#   GSDMB   — the primary gene of interest; encodes Gasdermin B, a pore-forming
#              protein involved in pyroptosis (inflammatory cell death); located
#              at chromosome 17q21.32 within the ORMDL3/GSDMB asthma risk locus
#   ORMDL3  — ER transmembrane protein; located adjacent to GSDMB; implicated
#              in ER stress responses and potentially in cancer signalling
#   IKZF3   — AIOLOS transcription factor; Tier 1 due to strong evidence in
#              haematological malignancies; also at the 17q21.32 locus
#   GSDMA   — paralogue of GSDMB; also a pyroptosis effector; included for
#              comparison within the gasdermin family
#   ERBB2   — HER2 oncogene; classically amplified/overexpressed in ~20% of
#              breast cancers; located at 17q12, adjacent to the GSDMB region
#   TP53    — the most commonly mutated gene in human cancer; included as a
#              reference Tier 1 tumour suppressor
#   BRCA1/2 — high-penetrance breast and ovarian cancer susceptibility genes;
#              included as reference germline cancer risk genes
# ──────────────────────────────────────────────────────────────────────────────

COSMIC_CENSUS_GENES = {
    "GSDMB":  {"Role": "Pyroptosis",           "Tier": 2, "Somatic": True,  "Germline": False},
    "ORMDL3": {"Role": "ER stress",             "Tier": 2, "Somatic": True,  "Germline": False},
    "IKZF3":  {"Role": "Transcription factor",  "Tier": 1, "Somatic": True,  "Germline": True},
    "ZPBP2":  {"Role": "Unknown",               "Tier": 3, "Somatic": False, "Germline": False},
    "GSDMA":  {"Role": "Pyroptosis",            "Tier": 2, "Somatic": True,  "Germline": False},
    "PGAP3":  {"Role": "GPI anchor",            "Tier": 3, "Somatic": False, "Germline": False},
    "ERBB2":  {"Role": "Oncogene",              "Tier": 1, "Somatic": True,  "Germline": True},
    "TP53":   {"Role": "Tumour suppressor",     "Tier": 1, "Somatic": True,  "Germline": True},
    "BRCA1":  {"Role": "Tumour suppressor",     "Tier": 1, "Somatic": True,  "Germline": True},
    "BRCA2":  {"Role": "Tumour suppressor",     "Tier": 1, "Somatic": True,  "Germline": True},
}


# ══════════════════════════════════════════════════════════════════════════════
# IMPACT SCORE MAPPING
# ══════════════════════════════════════════════════════════════════════════════
#
# VEP (Variant Effect Predictor) assigns each variant an IMPACT category
# reflecting the predicted severity of its functional consequence on the
# protein product:
#
#   HIGH     — likely loss of protein function (e.g. frameshift, stop gain,
#              splice site disruption). Score = 4.
#   MODERATE — likely change in protein function but not loss (e.g. missense
#              substitution, in-frame indel). Score = 3.
#   MODIFIER — non-coding or unlikely to affect protein function (e.g.
#              intronic, upstream gene variant). Score = 2.
#   LOW      — likely benign change (e.g. synonymous substitution). Score = 1.
#
# These integer scores are used in the Cancer Relevance classification to
# combine gene-level (Tier) and variant-level (Impact) evidence.
# ──────────────────────────────────────────────────────────────────────────────

IMPACT_SCORES = {
    "HIGH":     4,
    "MODERATE": 3,
    "MODIFIER": 2,
    "LOW":      1,
}


# ══════════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def find_col(df: pd.DataFrame, target: str):
    """
    Case-insensitive column name lookup.

    Excel files exported from different tools may capitalise column names
    inconsistently (e.g. "SYMBOL" vs "Symbol" vs "symbol"). This function
    normalises the comparison so the script works regardless of capitalisation.

    Parameters
    ----------
    df     : pd.DataFrame — the DataFrame to search
    target : str          — the desired column name (any case)

    Returns
    -------
    str or None
        The actual column name as it appears in the DataFrame, or None if the
        column does not exist.
    """
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None


# ══════════════════════════════════════════════════════════════════════════════
# MODULE 1: CANCER RELEVANCE CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════════

def classify_cosmic_relevance(row, sym_c: str = "SYMBOL", imp_c: str = "IMPACT") -> str:
    """
    Assign a Cancer Relevance category to a single variant.

    Classification logic
    --------------------
    The relevance score combines two independent pieces of evidence:

      Gene-level evidence (COSMIC Tier):
        How well-established is this gene's role in cancer?
        Tier 1 = strong, Tier 2 = moderate, Tier 3 = limited.

      Variant-level evidence (VEP IMPACT score):
        How functionally damaging is this specific variant predicted to be?
        HIGH (4) > MODERATE (3) > MODIFIER (2) > LOW (1).

    Classification rules (applied in order):
      - Gene not in COSMIC Census                       → "Unknown"
      - Tier 1 gene AND HIGH impact                     → "High"
      - Tier 1 or 2 gene AND MODERATE or higher impact  → "Medium"
      - Tier 1, 2, or 3 gene (any impact)               → "Low"

    This approach mirrors the tiered evidence framework used in clinical
    variant interpretation guidelines (e.g. ACMG/AMP).

    Parameters
    ----------
    row   : pd.Series — one row of the variant DataFrame
    sym_c : str       — name of the gene symbol column
    imp_c : str       — name of the VEP impact column

    Returns
    -------
    str : "High" | "Medium" | "Low" | "Unknown"
    """
    gene   = row[sym_c]
    impact = row[imp_c]

    # If the gene is not in our curated COSMIC reference, relevance is unknown
    if gene not in COSMIC_CENSUS_GENES:
        return "Unknown"

    gene_info    = COSMIC_CENSUS_GENES[gene]
    tier         = gene_info["Tier"]
    impact_score = IMPACT_SCORES.get(impact, 0)

    if tier == 1 and impact_score >= 4:
        return "High"
    elif tier <= 2 and impact_score >= 3:
        return "Medium"
    elif tier <= 3:
        return "Low"
    else:
        return "Unknown"


# ══════════════════════════════════════════════════════════════════════════════
# MODULE 2: MUTATIONAL HOTSPOT IDENTIFICATION
# ══════════════════════════════════════════════════════════════════════════════

def identify_hotspot_positions(
    df: pd.DataFrame,
    sym_c: str,
    pos_c: str,
    tissue_c: str,
    sam_c: str = "Sample",
) -> pd.DataFrame:
    """
    Identify genomic positions that are recurrently mutated across patients.

    Background
    ----------
    A mutational hotspot is a position in the genome that is independently
    mutated in multiple unrelated patients. The recurrence is biologically
    significant because:
      - Random sequencing errors are not reproducible between patients.
      - Passenger mutations (neutral mutations that arise by chance) are
        unlikely to recur at precisely the same position across many patients.
      - Recurrent mutations in cancer are therefore likely to confer a
        selective growth advantage — i.e. they are likely driver mutations.

    Famous examples of cancer hotspots include TP53 codon 248, KRAS codon 12,
    and PIK3CA H1047. This function searches for equivalent hotspot patterns
    within the GSDMB locus.

    Hotspot definition used here
    ----------------------------
    A position is classified as a hotspot if it is mutated in ≥3 independent
    samples within the same tissue group (tumour or healthy). This threshold
    is commonly used in targeted panel hotspot analyses and represents a
    conservative balance between sensitivity and specificity given typical
    cohort sizes.

    Recurrence rate calculation
    ---------------------------
    Recurrence_Rate_% = (samples_with_mutation / total_samples_in_tissue) × 100

    Parameters
    ----------
    df       : pd.DataFrame — variant table (one row per variant × transcript)
    sym_c    : str          — gene symbol column name
    pos_c    : str          — genomic position column name
    tissue_c : str          — tissue type column name
    sam_c    : str          — sample identifier column name

    Returns
    -------
    pd.DataFrame
        Hotspot table (sorted by sample count descending), one row per
        gene × position × tissue combination that meets the hotspot threshold.
        Empty DataFrame if no hotspots are found.
    """
    # Count unique samples per (gene, position, tissue) combination
    hotspots = (
        df.groupby([sym_c, pos_c, tissue_c])
        .agg({sam_c: "nunique"})
        .reset_index()
    )
    hotspots.columns = ["Gene", "Position", "Tissue", "Sample_Count"]

    # Flag positions meeting the recurrence threshold (≥3 independent samples)
    hotspots["Is_Hotspot"] = hotspots["Sample_Count"] >= 3

    # Calculate recurrence rate relative to the total samples in each tissue group
    total_samples_by_tissue = df.groupby(tissue_c)[sam_c].nunique().to_dict()
    hotspots["Recurrence_Rate_%"] = hotspots.apply(
        lambda row: (row["Sample_Count"] / total_samples_by_tissue.get(row["Tissue"], 1)) * 100,
        axis=1,
    )

    # Return only confirmed hotspots, sorted by recurrence
    return (
        hotspots[hotspots["Is_Hotspot"]]
        .sort_values("Sample_Count", ascending=False)
        .reset_index(drop=True)
    )


# ══════════════════════════════════════════════════════════════════════════════
# MODULE 3: VARIANT IMPACT BY TISSUE TYPE
# ══════════════════════════════════════════════════════════════════════════════

def summarise_tissue_impact(df: pd.DataFrame, tis_c: str, imp_c: str) -> pd.DataFrame:
    """
    Cross-tabulate variant functional impact categories by tissue type.

    Background
    ----------
    In this unpaired study design, direct somatic mutation calling (comparing
    a patient's tumour to their own matched normal) is not possible because
    tumour and healthy samples come from different patients. Instead, we
    compare the distribution of predicted functional impacts across the two
    tissue groups (Tumour vs. Healthy) at the population level.

    If high-impact variants (e.g. frameshift mutations, stop gains) are
    significantly more prevalent in tumour tissue, this provides indirect
    evidence for tumour-associated mutation enrichment at the population level.

    Impact categories (from VEP):
      HIGH     — predicted loss of protein function
      MODERATE — predicted change in protein function (e.g. missense)
      LOW      — predicted benign change (e.g. synonymous)
      MODIFIER — non-coding or negligible protein effect

    Parameters
    ----------
    df    : pd.DataFrame — variant table with tissue and impact columns
    tis_c : str          — tissue type column name
    imp_c : str          — VEP IMPACT column name

    Returns
    -------
    pd.DataFrame
        Cross-tabulation: rows = tissue types, columns = impact categories,
        values = variant counts. Columns are ordered HIGH → MODERATE →
        MODIFIER → LOW for intuitive reading.
    """
    tissue_impact = (
        df.groupby([tis_c, imp_c])
        .size()
        .unstack(fill_value=0)
    )

    # Reorder columns from most to least severe impact (where columns exist)
    ordered_cols = [c for c in ["HIGH", "MODERATE", "MODIFIER", "LOW"]
                    if c in tissue_impact.columns]
    return tissue_impact[ordered_cols]


# ══════════════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def analyze_cosmic_context():
    """
    Orchestrate the full COSMIC integration analysis.

    Execution order
    ---------------
    1. Load the annotated variant table.
    2. Standardise tissue type labels (e.g. "Tumor" → "Tumour").
    3. Module 1: Annotate variants with COSMIC Census gene information and
       assign Cancer Relevance scores.
    4. Module 2: Identify recurrent mutational hotspots (≥3 samples).
    5. Module 3: Summarise variant impact distribution by tissue type.
    6. Module 4: Determine cancer type specificity (breast vs endometrium
       vs pan-cancer) for tumour-tissue variants.
    7. Save all results to a multi-sheet Excel workbook.
    8. Generate and save the six-panel summary figure.
    9. Print a consolidated summary to the console.
    """
    print("=" * 70)
    print("SCRIPT 07c: COSMIC CANCER MUTATION DATABASE INTEGRATION")
    print("=" * 70)
    print(f"\nConfiguration:")
    print(f"  Input file        : {INPUT_FILE}")
    print(f"  COSMIC genes loaded: {len(COSMIC_CENSUS_GENES)}")
    print(f"  Analysis focus    : Cancer-relevant mutations")
    print(f"  Study design      : Unpaired (tumour vs. healthy, different patients)")
    print("\n" + "=" * 70 + "\n")

    # ── Load data ────────────────────────────────────────────────────────────
    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: Input file not found: {INPUT_FILE}")
        print("       Please update BASE_PATH in the CONFIGURATION block.")
        return

    df = pd.read_excel(INPUT_FILE, sheet_name="Biological_Annotations")
    print(f"Loaded {len(df)} rows from sheet 'Biological_Annotations'\n")

    # ── Identify columns (case-insensitive) ──────────────────────────────────
    sym_c = find_col(df, "SYMBOL")           # Gene symbol (e.g. "GSDMB")
    pos_c = find_col(df, "Pos")              # Genomic position
    tis_c = find_col(df, "TISSUE")           # Tissue type (Tumour / Healthy)
    sam_c = find_col(df, "SAMPLE")           # Sample identifier
    imp_c = find_col(df, "IMPACT")           # VEP impact category
    con_c = find_col(df, "CONSEQUENCE")      # VEP consequence (e.g. missense_variant)
    coh_c = find_col(df, "COHORT")           # Cancer cohort (e.g. Breast, Endometrium)
    hgv_c = find_col(df, "HGVSp")           # Protein-level HGVS notation
    var_c = find_col(df, "Existing_variation")  # Known variant IDs (e.g. rsIDs)

    # ── Standardise tissue labels ─────────────────────────────────────────────
    # Harmonise common alternative spellings to the two canonical labels used
    # throughout this pipeline: "Tumour" and "Healthy"
    df[tis_c] = (
        df[tis_c]
        .astype(str)
        .str.strip()
        .replace({"Tumor": "Tumour", "Normal": "Healthy", "Control": "Healthy"})
    )

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 1: COSMIC Census Gene Annotation & Cancer Relevance Scoring
    # ══════════════════════════════════════════════════════════════════════════
    print("STEP 1: Annotating Variants with COSMIC Cancer Gene Census")
    print("─" * 70)

    # Flag whether each variant's gene appears in our COSMIC reference
    df["In_COSMIC_Census"] = df[sym_c].isin(COSMIC_CENSUS_GENES.keys())

    # Look up the biological role for each gene from the census dictionary
    df["COSMIC_Role"] = df[sym_c].map(
        lambda x: COSMIC_CENSUS_GENES.get(x, {}).get("Role", "Not in census")
    )

    # Look up the evidence tier for each gene
    df["COSMIC_Tier"] = df[sym_c].map(
        lambda x: COSMIC_CENSUS_GENES.get(x, {}).get("Tier", np.nan)
    )

    # Assign Cancer Relevance score based on gene tier + variant impact
    df["Cancer_Relevance"] = df.apply(
        classify_cosmic_relevance, axis=1, sym_c=sym_c, imp_c=imp_c
    )

    # Report findings
    census_variants = df[df["In_COSMIC_Census"]]
    print(f"  Variants in COSMIC Census genes : {len(census_variants)} "
          f"across {census_variants[sym_c].nunique()} genes")
    print(f"  Cancer relevance breakdown:")
    print(f"    High    : {(df['Cancer_Relevance'] == 'High').sum()}")
    print(f"    Medium  : {(df['Cancer_Relevance'] == 'Medium').sum()}")
    print(f"    Low     : {(df['Cancer_Relevance'] == 'Low').sum()}")
    print(f"    Unknown : {(df['Cancer_Relevance'] == 'Unknown').sum()}\n")

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 2: Mutational Hotspot Identification
    # ══════════════════════════════════════════════════════════════════════════
    print("STEP 2: Identifying Mutational Hotspots (≥3 samples per position)")
    print("─" * 70)

    hotspots = identify_hotspot_positions(df, sym_c, pos_c, tis_c, sam_c)
    print(f"  Hotspot positions identified: {len(hotspots)}")

    if not hotspots.empty:
        print(f"\n  Top 5 hotspots:")
        for _, row in hotspots.head(5).iterrows():
            print(f"    {row['Gene']}:{row['Position']}  —  "
                  f"{row['Sample_Count']} samples "
                  f"({row['Recurrence_Rate_%']:.1f}% of {row['Tissue']} samples)")
    print()

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 3: Variant Impact by Tissue Type
    # ══════════════════════════════════════════════════════════════════════════
    print("STEP 3: Variant Impact Distribution — Tumour vs. Healthy")
    print("─" * 70)
    print("  NOTE: Somatic/germline classification is not performed.")
    print("  Tumour and healthy samples come from different patients (unpaired")
    print("  design). Comparison is made at the population level by tissue type.\n")

    tissue_impact = summarise_tissue_impact(df, tis_c, imp_c)

    for tissue, row in tissue_impact.iterrows():
        print(f"  {tissue}:")
        for impact_cat, count in row.items():
            print(f"    {impact_cat:10s}: {count}")
    print()

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 4: Cancer Type Specificity (Breast vs Endometrium vs Pan-Cancer)
    # ══════════════════════════════════════════════════════════════════════════
    print("STEP 4: Cancer Type Specificity Analysis")
    print("─" * 70)

    # Restrict to tumour samples only for this analysis
    # (healthy tissue does not carry cancer-type labels in the same way)
    tumour_variants = df[df[tis_c] == "Tumour"]

    # Build a pivot table: rows = gene × position, columns = cohort,
    # values = variant count. This shows which variants appear in which
    # cancer type cohorts.
    cancer_specific = (
        tumour_variants
        .groupby([sym_c, pos_c, coh_c])
        .size()
        .reset_index(name="Count")
    )
    cancer_specific = cancer_specific.pivot_table(
        index=[sym_c, pos_c],
        columns=coh_c,
        values="Count",
        fill_value=0,
    ).reset_index()

    # Classify each variant position by its cancer type distribution
    # (only if both breast and endometrium cohorts are present in the data)
    if "Breast" in cancer_specific.columns and "Endometrium" in cancer_specific.columns:
        cancer_specific["Breast_Specific"]      = (
            (cancer_specific["Breast"] > 0) & (cancer_specific["Endometrium"] == 0)
        )
        cancer_specific["Endometrium_Specific"] = (
            (cancer_specific["Endometrium"] > 0) & (cancer_specific["Breast"] == 0)
        )
        cancer_specific["Pan_Cancer"]           = (
            (cancer_specific["Breast"] > 0) & (cancer_specific["Endometrium"] > 0)
        )

        print(f"  Cancer specificity (tumour samples):")
        print(f"    Breast-only variants      : {cancer_specific['Breast_Specific'].sum()}")
        print(f"    Endometrium-only variants : {cancer_specific['Endometrium_Specific'].sum()}")
        print(f"    Pan-cancer variants       : {cancer_specific['Pan_Cancer'].sum()}\n")
    else:
        print("  WARNING: 'Breast' and/or 'Endometrium' cohort not found in data.")
        print("           Cancer specificity analysis skipped.\n")

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 5: Per-Gene Summary and Save Results
    # ══════════════════════════════════════════════════════════════════════════
    print("STEP 5: Generating Summary Reports and Saving Outputs")
    print("─" * 70)

    # Aggregate per gene: how many samples affected, total variants, high-relevance variants
    summary_by_gene = (
        df[df["In_COSMIC_Census"]]
        .groupby(sym_c)
        .agg(
            Unique_Samples        = (sam_c, "nunique"),
            Total_Variants        = (var_c, "count"),
            High_Relevance_Variants = ("Cancer_Relevance", lambda x: (x == "High").sum()),
        )
        .reset_index()
    )
    summary_by_gene.columns = [
        "Gene", "Unique_Samples", "Total_Variants", "High_Relevance_Variants"
    ]

    # Re-attach gene metadata from the census dictionary
    summary_by_gene["COSMIC_Role"] = summary_by_gene["Gene"].map(
        lambda x: COSMIC_CENSUS_GENES.get(x, {}).get("Role", "")
    )
    summary_by_gene["COSMIC_Tier"] = summary_by_gene["Gene"].map(
        lambda x: COSMIC_CENSUS_GENES.get(x, {}).get("Tier", np.nan)
    )
    summary_by_gene = summary_by_gene.sort_values("Total_Variants", ascending=False)

    # ── Save all results to a multi-sheet Excel workbook ─────────────────────
    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:

        # Sheet 1: Per-gene summary
        summary_by_gene.to_excel(
            writer, sheet_name="COSMIC_Gene_Summary", index=False
        )

        # Sheet 2: All variant rows that fall within COSMIC Census genes
        df[df["In_COSMIC_Census"]].copy().to_excel(
            writer, sheet_name="All_COSMIC_Variants", index=False
        )

        # Sheet 3: Only High Cancer Relevance variants (if any)
        high_relevance = df[df["Cancer_Relevance"] == "High"].copy()
        if not high_relevance.empty:
            high_relevance.to_excel(
                writer, sheet_name="High_Cancer_Relevance", index=False
            )

        # Sheet 4: Mutational hotspots
        if not hotspots.empty:
            hotspots.to_excel(
                writer, sheet_name="Mutational_Hotspots", index=False
            )

        # Sheet 5: Cancer type specificity (if cohort columns were present)
        if "Breast_Specific" in cancer_specific.columns:
            cancer_specific.to_excel(
                writer, sheet_name="Cancer_Specificity", index=False
            )

        # Sheet 6: Impact distribution by tissue type
        tissue_impact.to_excel(writer, sheet_name="Tissue_Impact_Summary")

    print(f"  Results saved: {OUTPUT_XLSX}\n")

    # ── Generate plots ────────────────────────────────────────────────────────
    generate_cosmic_plots(
        df, summary_by_gene, hotspots, cancer_specific,
        tissue_impact, sym_c, tis_c, coh_c, imp_c
    )

    # ── Print consolidated summary ────────────────────────────────────────────
    print("=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"\n  Total variant rows analysed       : {len(df)}")
    print(f"  Variants in COSMIC census genes   : {len(census_variants)} "
          f"({len(census_variants)/len(df)*100:.1f}%)")

    print(f"\n  Cancer relevance breakdown:")
    for level in ["High", "Medium", "Low", "Unknown"]:
        count = (df["Cancer_Relevance"] == level).sum()
        print(f"    {level:8s}: {count:5d}  ({count/len(df)*100:.1f}%)")

    print(f"\n  Top 5 most mutated COSMIC genes:")
    for _, row in summary_by_gene.head(5).iterrows():
        print(f"    {row['Gene']} ({row['COSMIC_Role']}): "
              f"{row['Total_Variants']} variants in {row['Unique_Samples']} samples")

    print("\n" + "=" * 70)
    print("COSMIC integration analysis complete.")
    print("=" * 70 + "\n")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE GENERATION
# ══════════════════════════════════════════════════════════════════════════════

def generate_cosmic_plots(
    df: pd.DataFrame,
    summary_by_gene: pd.DataFrame,
    hotspots: pd.DataFrame,
    cancer_specific: pd.DataFrame,
    tissue_impact: pd.DataFrame,
    sym_c: str,
    tis_c: str,
    coh_c: str,
    imp_c: str = "IMPACT",
) -> None:
    """
    Generate a six-panel COSMIC integration summary figure (300 dpi PNG).

    Panel layout (3 rows × 3 columns grid):
    ─────────────────────────────────────────
    Row 1, cols 0–1  (wide): Mutation burden per COSMIC Census gene (top 10)
                             Bars colour-coded by COSMIC Tier (red/amber/grey).

    Row 1, col 2     (small): Cancer Relevance pie chart
                              Proportion of variants classified as High /
                              Medium / Low / Unknown relevance.

    Row 2, col 0     : Variant impact stacked bar by tissue type
                       (Tumour vs. Healthy) — replaces somatic/germline plot
                       for unpaired study design.

    Row 2, col 1     : Top 10 mutational hotspots by recurrence rate (%).
                       Horizontal bar chart ordered by recurrence.

    Row 2, col 2     : Cancer type specificity bar chart
                       (Breast-only, Endometrium-only, Pan-cancer).

    Row 3, cols 0–2  (full width): COSMIC Census genes across cohorts
                       Stacked horizontal bar chart showing variant counts
                       per gene split by cancer type cohort.

    Parameters
    ----------
    df              : pd.DataFrame  — full variant table (with Cancer_Relevance column)
    summary_by_gene : pd.DataFrame  — per-gene summary from Step 5
    hotspots        : pd.DataFrame  — hotspot table from Module 2
    cancer_specific : pd.DataFrame  — cancer specificity pivot table from Step 4
    tissue_impact   : pd.DataFrame  — tissue × impact cross-tabulation from Module 3
    sym_c           : str           — gene symbol column name
    tis_c           : str           — tissue type column name
    coh_c           : str           — cohort column name
    imp_c           : str           — IMPACT column name
    """
    print("Generating COSMIC analysis plots...")

    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(18, 12))
    gs  = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)

    # Colour scheme for COSMIC Tiers
    # Tier 1 (red) = strongest evidence; Tier 2 (amber); Tier 3 (grey)
    tier_colours = {1: "#e74c3c", 2: "#f39c12", 3: "#95a5a6"}

    # ── Panel 1: Mutation burden per COSMIC Census gene ──────────────────────
    ax1 = fig.add_subplot(gs[0, :2])   # Spans first two columns of row 1

    top_genes  = summary_by_gene.head(10)
    bar_colors = [tier_colours.get(tier, "#95a5a6") for tier in top_genes["COSMIC_Tier"]]

    ax1.barh(range(len(top_genes)), top_genes["Total_Variants"], color=bar_colors)
    ax1.set_yticks(range(len(top_genes)))
    ax1.set_yticklabels(
        [f"{gene} ({role})" for gene, role in zip(top_genes["Gene"], top_genes["COSMIC_Role"])]
    )
    ax1.set_xlabel("Total Variants", fontweight="bold")
    ax1.set_title("Top 10 COSMIC Census Genes by Mutation Burden",
                  fontweight="bold", fontsize=12)
    ax1.grid(True, alpha=0.3, axis="x")

    # Custom legend explaining the tier colour coding
    legend_elements = [
        Patch(facecolor="#e74c3c", label="Tier 1 — Strong evidence"),
        Patch(facecolor="#f39c12", label="Tier 2 — Moderate evidence"),
        Patch(facecolor="#95a5a6", label="Tier 3 — Limited evidence"),
    ]
    ax1.legend(handles=legend_elements, loc="lower right", fontsize=9)

    # ── Panel 2: Cancer Relevance pie chart ───────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 2])

    relevance_counts = df["Cancer_Relevance"].value_counts()
    relevance_colours = {
        "High": "#e74c3c", "Medium": "#f39c12",
        "Low": "#3498db",  "Unknown": "#95a5a6",
    }
    pie_colours = [relevance_colours.get(cat, "#95a5a6") for cat in relevance_counts.index]

    ax2.pie(
        relevance_counts.values,
        labels=relevance_counts.index,
        autopct="%1.1f%%",
        colors=pie_colours,
        startangle=90,
    )
    ax2.set_title("Cancer Relevance\nClassification", fontweight="bold", fontsize=11)

    # ── Panel 3: Variant impact stacked bar by tissue type ────────────────────
    ax3 = fig.add_subplot(gs[1, 0])

    if not tissue_impact.empty:
        tissue_impact.plot(
            kind="bar", stacked=True, ax=ax3,
            color=["#d9534f", "#5bc0de", "#f0ad4e", "#5cb85c"],
        )
        ax3.set_xlabel("Tissue Type", fontweight="bold")
        ax3.set_ylabel("Variant Count", fontweight="bold")
        ax3.set_title(
            "Variant Impact by Tissue Type\n(Tumour vs. Healthy — unpaired design)",
            fontweight="bold", fontsize=10,
        )
        plt.setp(ax3.xaxis.get_majorticklabels(), rotation=0)
        ax3.legend(title="VEP Impact", loc="upper right", fontsize=8)
    else:
        ax3.text(0.5, 0.5, "No tissue impact data available",
                 ha="center", va="center", fontsize=11)
        ax3.axis("off")

    # ── Panel 4: Mutational hotspot recurrence ────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 1])

    if not hotspots.empty:
        top_hotspots = (
            hotspots.nlargest(10, "Recurrence_Rate_%")
            .sort_values("Recurrence_Rate_%", ascending=True)
        )
        ax4.barh(range(len(top_hotspots)), top_hotspots["Recurrence_Rate_%"],
                 color="#e74c3c", alpha=0.7)
        ax4.set_yticks(range(len(top_hotspots)))
        ax4.set_yticklabels(
            [f"{row['Gene']}:{row['Position']}" for _, row in top_hotspots.iterrows()],
            fontsize=8,
        )
        ax4.set_xlabel("Recurrence Rate (%)", fontweight="bold")
        ax4.set_title("Top Mutational Hotspots\n(% of samples in tissue group)",
                      fontweight="bold", fontsize=11)
        ax4.grid(True, alpha=0.3, axis="x")
    else:
        ax4.text(0.5, 0.5, "No hotspots identified\n(threshold: ≥3 samples)",
                 ha="center", va="center", fontsize=11)
        ax4.axis("off")

    # ── Panel 5: Cancer type specificity ──────────────────────────────────────
    ax5 = fig.add_subplot(gs[1, 2])

    if "Breast_Specific" in cancer_specific.columns:
        specificity_counts = {
            "Breast-specific":       cancer_specific["Breast_Specific"].sum(),
            "Endometrium-specific":  cancer_specific["Endometrium_Specific"].sum(),
            "Pan-cancer":            cancer_specific["Pan_Cancer"].sum(),
        }
        ax5.bar(
            specificity_counts.keys(),
            specificity_counts.values(),
            color=["#e74c3c", "#2ecc71", "#3498db"],
            alpha=0.75,
            edgecolor="black",
        )
        ax5.set_ylabel("Number of Variant Positions", fontweight="bold")
        ax5.set_title("Cancer Type Specificity\n(tumour samples only)",
                      fontweight="bold", fontsize=11)
        plt.setp(ax5.xaxis.get_majorticklabels(), rotation=15)
        ax5.grid(True, alpha=0.3, axis="y")
    else:
        ax5.text(0.5, 0.5, "Cancer specificity data\nnot available",
                 ha="center", va="center", fontsize=11)
        ax5.axis("off")

    # ── Panel 6: COSMIC Census genes across cancer cohorts (full-width row) ───
    ax6 = fig.add_subplot(gs[2, :])   # Spans all three columns of row 3

    # Pivot: rows = genes, columns = cohorts, values = variant counts
    cosmic_by_cohort = (
        df[df["In_COSMIC_Census"]]
        .groupby([sym_c, coh_c])
        .size()
        .unstack(fill_value=0)
    )

    if len(cosmic_by_cohort) > 0:
        # Show top 15 genes by total variant count
        top_cosmic_genes = cosmic_by_cohort.sum(axis=1).nlargest(15)
        plot_data = cosmic_by_cohort.loc[top_cosmic_genes.index]

        plot_data.plot(
            kind="barh", stacked=True, ax=ax6,
            color=["#e74c3c", "#2ecc71", "#3498db"],
        )
        ax6.set_xlabel("Variant Count", fontweight="bold")
        ax6.set_ylabel("Gene", fontweight="bold")
        ax6.set_title("COSMIC Census Genes Across Cancer Types (top 15 genes)",
                      fontweight="bold", fontsize=12)
        ax6.legend(title="Cohort", bbox_to_anchor=(1.02, 1), loc="upper left")
        ax6.grid(True, alpha=0.3, axis="x")
    else:
        ax6.text(0.5, 0.5, "No COSMIC Census variants found in data",
                 ha="center", va="center", fontsize=12)
        ax6.axis("off")

    # ── Save figure ───────────────────────────────────────────────────────────
    plt.suptitle(
        "COSMIC Cancer Gene Census Integration Analysis",
        fontsize=16, fontweight="bold", y=1.0,
    )
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"  COSMIC plots saved: {OUTPUT_PLOT}\n")


# ── Script entry point ─────────────────────────────────────────────────────────
# This block ensures that analyze_cosmic_context() is only called when the
# script is executed directly, not when it is imported as a module.
if __name__ == "__main__":
    analyze_cosmic_context()
