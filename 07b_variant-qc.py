#!/usr/bin/env python3
"""
================================================================================
Script 07b — Variant-Level Quality Control
================================================================================
Pipeline    : Step 07b of 17 — runs AFTER annotation merge (07) and BEFORE
              genomic mapping (08)

--------------------------------------------------------------------------------
PURPOSE
--------------------------------------------------------------------------------
This script performs quality control (QC) at the level of individual genetic
variants. It takes as input the fully annotated variant table produced by
script 07 and applies three independent, complementary checks to determine
whether the variant calls are biologically plausible and technically reliable.

The three checks are:

  1. Ti/Tv ratio (Transition-to-Transversion ratio)
     ─────────────────────────────────────────────
     DNA mutations can be classified into two categories:
       - Transitions (Ti): mutations between chemically similar bases
         (purine ↔ purine, or pyrimidine ↔ pyrimidine).
         Example: Adenine→Guanine, Cytosine→Thymine.
       - Transversions (Tv): mutations between chemically dissimilar bases
         (purine ↔ pyrimidine).
         Example: Adenine→Cytosine, Guanine→Thymine.

     In the human genome, transitions are biologically more likely than
     transversions, so real variant data produces a Ti/Tv ratio of roughly
     2.0–2.1 for whole-exome or targeted panel sequencing. A ratio that
     deviates substantially from this range can indicate:
       - Too low  (< 1.5): poor-quality base calls, sequencing artefacts,
                           or insufficient variants to calculate reliably.
       - Too high (> 3.0): over-filtering of transversions, enrichment of
                           CpG sites (which mutate as C→T preferentially).

  2. Allelic Balance (AB)
     ─────────────────────
     A heterozygous variant (one mutated copy, one normal copy) should, in
     theory, have reads supporting each allele at roughly equal frequencies —
     i.e. an allele fraction (AF) of approximately 0.5 (50%). Strong
     deviations from 0.5 in heterozygous calls can indicate:
       - Sample contamination (another sample's DNA artificially elevating
         one allele)
       - Loss of heterozygosity (LOH): one allele is deleted in tumour tissue
       - Copy number variations (CNV): gene amplification skewing read counts
       - Strand bias or PCR amplification artefacts

     Each heterozygous call is flagged as PASS, ACCEPTABLE, or IMBALANCED
     depending on how far its AF deviates from 0.5.

  3. Batch Effects
     ──────────────
     "Batch effects" refer to technical differences in variant counts that
     arise not from true biological variation, but from artefacts introduced
     by processing samples in different experimental batches (e.g. different
     sequencing runs, reagent lots, or library preparation dates). This check
     computes the coefficient of variation (CV%) of variant counts across
     samples within each cohort/tissue group. A high CV% (> 30%) suggests
     that some samples in the group have many more or fewer variants than
     their peers — a potential sign of a technical problem requiring
     investigation.

--------------------------------------------------------------------------------
IMPORTANT NOTE ON AMPLICON COVERAGE QC
--------------------------------------------------------------------------------
This script does NOT perform per-amplicon coverage quality control.
Coverage QC requires the complete depth track for every amplicon (mosdepth
output), not just the positions where a variant was called. Amplicons with
no variants would be entirely invisible here and would be silently omitted,
producing a misleading picture of panel performance. For amplicon-level
coverage assessment, refer to scripts 02, 03, and 04.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        One row per variant × transcript combination (VEP annotation format).
        Required columns: CHROM, POS, REF, ALT, Sample, GT, AF
        Optional columns: Cohort, Tissue (used for batch effect grouping)

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------
  <output_dir>/
    ├── Variant_QC_Report.xlsx       Full annotated table with AB flags appended
    ├── TiTv_Summary.csv             Transition/transversion counts and ratio
    ├── Batch_Effects_Summary.csv    Per-group variant count statistics
    └── Variant_QC_Plots.png         Three-panel summary figure (300 dpi)

--------------------------------------------------------------------------------
USAGE EXAMPLES
--------------------------------------------------------------------------------
  # Minimal usage (default thresholds):
  python 07b_variant-qc.py --input /path/to/GSDMB_Annotated_Report_Fixed.xlsx

  # Custom allelic balance thresholds and explicit output directory:
  python 07b_variant-qc.py \\
      --input annotations.xlsx \\
      --output ./qc_results \\
      --ab-lower 0.20 \\
      --ab-upper 0.80

  # Verbose mode (prints debug-level messages):
  python 07b_variant-qc.py --input annotations.xlsx --verbose

--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------
  Python ≥ 3.8
  pandas        (tabular data manipulation)
  numpy         (numerical operations, NaN handling)
  matplotlib    (figure generation)
  openpyxl      (Excel file reading/writing, used by pandas)

  Install: pip install pandas numpy matplotlib openpyxl
================================================================================
"""

# ──────────────────────────────────────────────────────────────────────────────
# STANDARD LIBRARY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import argparse    # Command-line argument parsing
import logging     # Structured log messages (timestamps, levels: INFO/WARNING/ERROR)
import os          # File and directory path operations
import sys         # System-level operations (e.g. exit on error)

# ──────────────────────────────────────────────────────────────────────────────
# THIRD-PARTY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt  # Generating the QC summary figure
import numpy as np               # Numerical operations (NaN, array math)
import pandas as pd              # DataFrame operations (loading, filtering, grouping)


# ══════════════════════════════════════════════════════════════════════════════
# QUALITY CONTROL THRESHOLDS
# ══════════════════════════════════════════════════════════════════════════════
#
# All numeric thresholds are centralised here so that they are easy to find,
# understand, and adjust without hunting through the code.
#
# Justification for threshold values:
#
#   Ti/Tv:
#     Optimal range 2.0–2.1 reflects the expected ratio for human coding
#     sequence under neutral evolution (Ts/Tv ~2.0; Ajay et al., 2011;
#     DePristo et al., 2011 GATK paper). Acceptable range 1.5–3.0 is the
#     wider window used in clinical and research WES pipelines.
#
#   Allelic balance:
#     A true germline heterozygous variant should have AF ≈ 0.5 ± noise.
#     Lower bound 0.25 and upper bound 0.75 are standard in GATK best
#     practices for flagging potential artefacts. The optimal window 0.40–0.60
#     reflects minimal sampling noise.
#
#   Batch CV%:
#     Coefficient of variation thresholds are empirical. CV < 15% is
#     considered low variation; CV 15–30% is acceptable for heterogeneous
#     cancer cohorts; CV > 50% suggests outlier samples or batch problems.
# ──────────────────────────────────────────────────────────────────────────────

THRESHOLDS = {
    "titv": {
        "optimal_min":   2.0,   # Ideal minimum Ti/Tv for targeted sequencing
        "optimal_max":   2.1,   # Ideal maximum Ti/Tv for targeted sequencing
        "acceptable_min": 1.5,  # Below this → likely artefacts or too few variants
        "acceptable_max": 3.0,  # Above this → possible CpG enrichment or over-filtering
    },
    "allelic_balance": {
        "lower":         0.25,  # AF below this in a het call → flagged as imbalanced
        "upper":         0.75,  # AF above this in a het call → flagged as imbalanced
        "optimal_lower": 0.40,  # AF range considered truly well-balanced
        "optimal_upper": 0.60,
    },
    "batch_cv": {
        "excellent":  15,  # CV% < 15 → very consistent variant counts across samples
        "acceptable": 30,  # CV% 15–30 → acceptable for mixed cancer cohorts
        "poor":       50,  # CV% > 50 → investigate for batch artefacts or outliers
    },
}


# ══════════════════════════════════════════════════════════════════════════════
# LOGGING SETUP
# ══════════════════════════════════════════════════════════════════════════════

def setup_logging(verbose: bool = False) -> logging.Logger:
    """
    Configure the Python logging module for this script.

    Parameters
    ----------
    verbose : bool
        If True, DEBUG-level messages are also printed (more detail).
        If False (default), only INFO, WARNING, and ERROR messages appear.

    Returns
    -------
    logging.Logger
        Configured logger instance used throughout the script.

    Log format example:
        14:32:05 [INFO] Loaded 4821 rows from Biological_Annotations
    """
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════════════════════
# UTILITY / HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def find_col(df: pd.DataFrame, name: str):
    """
    Case-insensitive column name lookup.

    Excel files exported from different tools sometimes capitalise column names
    inconsistently (e.g. "CHROM" vs "Chrom" vs "chrom"). This function
    normalises the comparison so the script handles any capitalisation.

    Parameters
    ----------
    df   : pd.DataFrame  — the DataFrame to search
    name : str           — the desired column name (any case)

    Returns
    -------
    str or None
        The actual column name as it appears in the DataFrame, or None if not
        found.
    """
    for col in df.columns:
        if col.upper() == name.upper():
            return col
    return None


def deduplicate_to_variants(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Collapse transcript-expanded rows to one row per unique variant call.

    Background
    ----------
    Script 07 (annotation merge) uses VEP (Variant Effect Predictor) to
    annotate variants. VEP expands the CSQ (consequence) field so that each
    variant × transcript combination becomes its own row. For example, a
    single SNP affecting 3 transcripts will appear as 3 rows in the table.

    Why this matters for QC
    -----------------------
    If we count transitions and transversions, or flag allelic imbalance,
    without first collapsing these duplicates, every variant will be counted
    once per transcript it affects. This would inflate the Ti/Tv counts and
    generate duplicate allelic balance flags for the same underlying call.

    Deduplication strategy
    ----------------------
    The combination (CHROM, POS, REF, ALT, Sample) uniquely identifies a
    variant call. We keep one row per such combination, discarding the
    additional transcript-level rows. All downstream QC operates on this
    deduplicated table.

    Parameters
    ----------
    df     : pd.DataFrame   — full annotated table from script 07
    logger : logging.Logger — for progress messages

    Returns
    -------
    pd.DataFrame
        Deduplicated table with one row per unique variant call.
    """
    # Identify which of the five key columns are present in this table
    key_cols = [c for c in ["CHROM", "POS", "REF", "ALT", "Sample"]
                if find_col(df, c)]

    if not key_cols:
        logger.warning(
            "Cannot identify variant key columns (CHROM, POS, REF, ALT, Sample) — "
            "using full table. Counts may be inflated by transcript duplicates."
        )
        return df

    # Standardise column names to uppercase so the rest of the script can
    # refer to them consistently regardless of the source file capitalisation
    rename_map = {find_col(df, c): c for c in key_cols}
    df = df.rename(columns=rename_map)

    n_before = len(df)
    df = df.drop_duplicates(subset=key_cols)
    n_after = len(df)

    logger.info(
        f"Deduplicated {n_before} transcript-level rows → "
        f"{n_after} unique variant calls "
        f"(removed {n_before - n_after} duplicate transcript rows)"
    )
    return df


# ══════════════════════════════════════════════════════════════════════════════
# QC MODULE 1: TRANSITION / TRANSVERSION (Ti/Tv) RATIO
# ══════════════════════════════════════════════════════════════════════════════

# These sets define which single-nucleotide substitutions are transitions and
# which are transversions.
#
# TRANSITIONS (chemically similar base exchange):
#   Purines:     A ↔ G   (both have double-ring structure)
#   Pyrimidines: C ↔ T   (both have single-ring structure)
#
# TRANSVERSIONS (chemically dissimilar base exchange):
#   Purine ↔ Pyrimidine: A↔C, A↔T, G↔C, G↔T (and their reverses)
#
# The notation "A>G" means reference allele A was changed to alternate allele G.

TRANSITIONS   = {"A>G", "G>A", "C>T", "T>C"}
TRANSVERSIONS = {"A>C", "C>A", "A>T", "T>A", "G>C", "C>G", "G>T", "T>G"}


def calculate_titv(df: pd.DataFrame, logger: logging.Logger) -> dict:
    """
    Calculate the Transition-to-Transversion (Ti/Tv) ratio for all SNVs.

    Background
    ----------
    Single-nucleotide variants (SNVs) fall into two mutational classes:
    transitions and transversions (see TRANSITIONS / TRANSVERSIONS sets above).
    Because transitions are chemically more probable, real human sequencing data
    consistently shows Ti/Tv ≈ 2.0–2.1 for coding regions (whole-exome or
    targeted panel). This ratio is widely used as a sequencing quality metric:

      Ti/Tv < 1.5 → Too low — likely dominated by sequencing errors, which
                    tend to produce random base changes (enriching Tv) rather
                    than the biologically expected transitions. Could also
                    indicate too few SNVs to calculate meaningfully.

      Ti/Tv 1.5–2.0 or 2.1–3.0 → Acceptable but worth noting. Some panels
                    and targeted regions show values in this range due to
                    CpG site enrichment (C→T transitions at CpG dinucleotides
                    are exceptionally frequent) or specific filtering choices.

      Ti/Tv 2.0–2.1 → Optimal for typical exonic/targeted sequencing.

      Ti/Tv > 3.0 → High — possible over-filtering of transversions, or a
                    dataset strongly enriched for CpG mutations.

    Indels (insertions or deletions) are skipped — they are not SNVs and do
    not contribute to the Ti/Tv calculation.

    Parameters
    ----------
    df     : pd.DataFrame   — deduplicated variant table
    logger : logging.Logger — for progress and warning messages

    Returns
    -------
    dict with keys:
        ti        (int)   — count of transitions
        tv        (int)   — count of transversions
        snv_total (int)   — total SNVs counted (ti + tv)
        ratio     (float) — Ti/Tv ratio (NaN if tv == 0)
        status    (str)   — human-readable QC verdict
        colour    (str)   — hex colour code for plot formatting
    """
    ref_col = find_col(df, "REF")
    alt_col = find_col(df, "ALT")

    # Guard: both REF and ALT columns must be present
    if not ref_col or not alt_col:
        logger.warning("REF or ALT column missing — cannot calculate Ti/Tv ratio.")
        return {
            "ti": 0, "tv": 0, "snv_total": 0,
            "ratio": float("nan"),
            "status": "FAIL — REF or ALT column missing",
            "colour": "#e74c3c",  # Red
        }

    ti = tv = skipped = 0

    for ref, alt in zip(
        df[ref_col].astype(str).str.upper(),
        df[alt_col].astype(str).str.upper()
    ):
        # Skip indels: an SNV must have exactly one base for both REF and ALT
        if len(ref) != 1 or len(alt) != 1:
            skipped += 1
            continue

        # Skip ambiguous bases (anything outside standard A/C/G/T)
        if ref not in "ACGT" or alt not in "ACGT":
            skipped += 1
            continue

        # Classify the substitution
        mutation = f"{ref}>{alt}"
        if mutation in TRANSITIONS:
            ti += 1
        elif mutation in TRANSVERSIONS:
            tv += 1

    snv_total = ti + tv
    logger.info(
        f"SNV classification: transitions={ti}, transversions={tv}, "
        f"total SNVs={snv_total}, indels/other skipped={skipped}"
    )

    # Guard: cannot divide by zero if no transversions found
    if tv == 0:
        logger.warning(
            "No transversions found — Ti/Tv ratio cannot be calculated. "
            "This may indicate too few variants or an extreme filter applied upstream."
        )
        return {
            "ti": ti, "tv": tv, "snv_total": snv_total,
            "ratio": float("nan"),
            "status": "FAIL — no transversions detected",
            "colour": "#e74c3c",
        }

    ratio = ti / tv
    t = THRESHOLDS["titv"]

    # Assign QC status based on the ratio value
    if ratio < t["acceptable_min"]:
        status = "FAIL — ratio too low (likely artefacts or insufficient variants)"
        colour = "#e74c3c"  # Red
    elif ratio > t["acceptable_max"]:
        status = "WARN — ratio high (CpG site enrichment or transversion over-filtering?)"
        colour = "#f39c12"  # Amber
    elif t["optimal_min"] <= ratio <= t["optimal_max"]:
        status = "OPTIMAL"
        colour = "#2ecc71"  # Green
    else:
        status = "ACCEPTABLE"
        colour = "#3498db"  # Blue

    logger.info(f"Ti/Tv ratio = {ratio:.3f}  →  {status}")
    return {
        "ti": ti, "tv": tv, "snv_total": snv_total,
        "ratio": ratio, "status": status, "colour": colour,
    }


# ══════════════════════════════════════════════════════════════════════════════
# QC MODULE 2: ALLELIC BALANCE
# ══════════════════════════════════════════════════════════════════════════════

def check_allelic_balance(
    df: pd.DataFrame,
    logger: logging.Logger,
    ab_lower: float,
    ab_upper: float,
) -> pd.DataFrame:
    """
    Flag heterozygous variant calls where the allele fraction deviates
    unexpectedly from the expected value of 0.5.

    Background
    ----------
    A heterozygous genotype (written as "0/1" in VCF format) means the
    individual has one reference allele and one alternate allele at that
    position. In a pure diploid sample with balanced sequencing, we expect
    approximately half the reads to support each allele, giving an allele
    fraction (AF) of ~0.5.

    If a heterozygous call shows a strongly skewed AF (e.g. 0.10 or 0.90),
    this is called allelic imbalance and can be caused by:

      Contamination:
        DNA from another sample artificially inflates one allele's read count.

      Loss of heterozygosity (LOH):
        One allele is physically deleted in the tumour sample, so only one
        allele remains — the surviving allele appears at high frequency.

      Copy number variation (CNV):
        Gene amplification can create many copies of one allele, skewing
        the read fraction away from 0.5.

      Technical artefacts:
        Strand bias, PCR amplification errors, or library preparation
        problems can all produce apparent imbalance in sequencing data.

    How this function works
    -----------------------
    1. Identify heterozygous calls: GT column contains "0/1" or "0|1"
       (the pipe "|" notation is used for phased genotypes and is equivalent
       here).
    2. Retrieve the AF value (allele fraction) for each heterozygous call.
    3. Compare AF to the configured thresholds and assign a quality flag.

    Three output columns are added to the DataFrame:
      AB_AF       (float)  — the allele fraction used (NaN for non-het calls)
      AB_Flag     (str)    — PASS | IMBALANCED (AF=X.XXX) | NO_DATA
      AB_Quality  (str)    — OPTIMAL | ACCEPTABLE | POOR | HOM/MISSING

    Parameters
    ----------
    df        : pd.DataFrame   — deduplicated variant table
    logger    : logging.Logger — for progress and warning messages
    ab_lower  : float          — lower AF threshold (from CLI or THRESHOLDS)
    ab_upper  : float          — upper AF threshold (from CLI or THRESHOLDS)

    Returns
    -------
    pd.DataFrame
        Original table with three new columns appended (AB_AF, AB_Flag,
        AB_Quality).
    """
    t = THRESHOLDS["allelic_balance"]
    df = df.copy()

    # Initialise the three new output columns with default values
    # These defaults apply to homozygous calls and rows where data are missing
    df["AB_AF"]      = np.nan    # NaN = not applicable (not a heterozygous call)
    df["AB_Flag"]    = "NO_DATA"
    df["AB_Quality"] = "HOM/MISSING"

    gt_col = find_col(df, "GT")   # Genotype column (e.g. "0/1", "1/1")
    af_col = find_col(df, "AF")   # Allele fraction column (e.g. 0.487)

    if not gt_col:
        logger.warning(
            "GT (genotype) column not found — skipping allelic balance check. "
            "This column is required to identify heterozygous calls."
        )
        return df

    # ── Identify heterozygous calls ──────────────────────────────────────────
    # A heterozygous call has two different alleles separated by "/" or "|"
    # Examples of heterozygous GT values:  "0/1", "0|1", "1/2"
    # Examples of homozygous GT values:    "0/0", "1/1"
    # We skip any GT containing a "." (missing data, e.g. "./.")
    het_indices = []
    for row_index, gt_raw in df[gt_col].items():
        gt = str(gt_raw).replace("|", "/")   # Treat phased (|) same as unphased (/)
        if "/" not in gt:
            continue
        alleles = gt.split("/")
        if len(alleles) == 2 and alleles[0] != alleles[1] and "." not in alleles:
            het_indices.append(row_index)

    logger.info(
        f"Heterozygous calls identified: {len(het_indices)} out of {len(df)} total variant rows"
    )

    if not het_indices:
        logger.warning(
            "No heterozygous calls found — all calls are homozygous or have missing "
            "genotype data. Allelic balance check will not produce any flags."
        )
        return df

    if not af_col:
        logger.warning(
            "AF (allele fraction) column not found — cannot compute allelic balance. "
            "The AF column is required to check the balance of heterozygous calls."
        )
        return df

    # ── Evaluate each heterozygous call ─────────────────────────────────────
    n_flagged = n_optimal = n_acceptable = 0

    for idx in het_indices:
        try:
            af = float(df.at[idx, af_col])
        except (ValueError, TypeError):
            # AF value is missing or non-numeric for this row; skip silently
            continue

        # Store the computed AF for this call
        df.at[idx, "AB_AF"] = af

        # Classify based on how far AF deviates from the ideal 0.5
        if af < ab_lower or af > ab_upper:
            # Strongly imbalanced — outside the acceptable window
            df.at[idx, "AB_Flag"]    = f"IMBALANCED (AF={af:.3f})"
            df.at[idx, "AB_Quality"] = "POOR"
            n_flagged += 1

        elif t["optimal_lower"] <= af <= t["optimal_upper"]:
            # Well-balanced — within the optimal window around 0.5
            df.at[idx, "AB_Flag"]    = "PASS"
            df.at[idx, "AB_Quality"] = "OPTIMAL"
            n_optimal += 1

        else:
            # Slightly off-centre but still within the acceptable window
            df.at[idx, "AB_Flag"]    = "PASS"
            df.at[idx, "AB_Quality"] = "ACCEPTABLE"
            n_acceptable += 1

    # ── Report summary statistics ────────────────────────────────────────────
    n_total   = len(het_indices)
    pct_flagged = n_flagged / n_total * 100

    logger.info(
        f"Allelic balance results: "
        f"optimal={n_optimal} ({n_optimal/n_total*100:.1f}%), "
        f"acceptable={n_acceptable} ({n_acceptable/n_total*100:.1f}%), "
        f"imbalanced={n_flagged} ({pct_flagged:.1f}%)"
    )

    # Warn if a large proportion of heterozygous calls are imbalanced
    if pct_flagged > 20:
        logger.warning(
            f"{pct_flagged:.0f}% of heterozygous calls are imbalanced — "
            "this exceeds the 20% warning threshold. "
            "Consider investigating for sample contamination, loss of heterozygosity (LOH), "
            "or copy number variations (CNV) in affected samples."
        )

    return df


# ══════════════════════════════════════════════════════════════════════════════
# QC MODULE 3: BATCH EFFECT DETECTION
# ══════════════════════════════════════════════════════════════════════════════

def detect_batch_effects(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Detect potential batch effects by examining the variability of variant
    counts across samples within each cohort and tissue group.

    Background
    ----------
    In a well-controlled experiment, samples processed under the same conditions
    should yield broadly comparable numbers of detected variants. When variant
    counts vary greatly between samples within the same biological group, this
    is suspicious and may indicate:

      Batch effects:
        Samples sequenced in different runs, with different reagent lots, or
        by different operators may have systematically different call rates
        unrelated to biology.

      Genuine biological heterogeneity:
        Cancer samples can differ substantially in their mutation burden (the
        total number of somatic mutations). A high variant count may reflect a
        hypermutated tumour rather than a technical problem.

      Sample quality outliers:
        Low-quality samples (degraded DNA, insufficient input) may yield
        fewer variants. Such outliers should be reviewed and potentially
        excluded before downstream analysis.

    Metric: Coefficient of Variation (CV%)
    ─────────────────────────────────────
    CV% = (standard deviation / mean) × 100

    This metric expresses variability as a percentage of the mean, allowing
    fair comparison between groups with different absolute variant counts.

      CV% < 15  →  EXCELLENT: very consistent counts across samples
      CV% 15–30 →  ACCEPTABLE: some variation, typical for cancer cohorts
      CV% 30–50 →  HIGH VARIATION: worth investigating
      CV% > 50  →  VERY HIGH VARIATION: likely batch problem or outlier

    Parameters
    ----------
    df     : pd.DataFrame   — deduplicated variant table
    logger : logging.Logger — for progress and warning messages

    Returns
    -------
    pd.DataFrame
        Summary table, one row per group (cohort × tissue combination), with
        columns: Group, N_Samples, Mean_Variants, Median_Variants,
                 Std_Variants, CV_%, Min_Variants, Max_Variants, QC_Flag.
        Returns an empty DataFrame if the Sample column is missing.
    """
    cohort_col = find_col(df, "Cohort")
    tissue_col = find_col(df, "Tissue")
    sample_col = find_col(df, "Sample")

    if not sample_col:
        logger.warning(
            "Sample column not found — skipping batch effect detection. "
            "A 'Sample' column is required to count variants per sample."
        )
        return pd.DataFrame()

    # Determine grouping columns (cohort and/or tissue if available)
    group_cols = [c for c in [cohort_col, tissue_col] if c]

    if not group_cols:
        # No grouping metadata available — treat all samples as one global group
        logger.warning(
            "Neither 'Cohort' nor 'Tissue' column found — "
            "batch analysis will treat all samples as a single group."
        )
        df = df.copy()
        df["_group"] = "All samples"
        group_cols   = ["_group"]

    # Count the number of variants per sample within each group
    # This gives one row per (group, sample) with the variant count
    per_sample_counts = (
        df.groupby(group_cols + [sample_col])
        .size()
        .reset_index(name="N_Variants")
    )

    t = THRESHOLDS["batch_cv"]
    summary_rows = []

    for group_keys, group_data in per_sample_counts.groupby(group_cols):

        # Normalise single-key groups to a tuple for uniform handling
        if not isinstance(group_keys, tuple):
            group_keys = (group_keys,)

        # Human-readable label for this group (e.g. "CancerCohort / Lung")
        group_label = " / ".join(str(k) for k in group_keys)

        mean_count   = group_data["N_Variants"].mean()
        std_count    = group_data["N_Variants"].std(ddof=1) if len(group_data) > 1 else 0.0
        cv_pct       = std_count / mean_count * 100 if mean_count > 0 else 0.0

        # Assign QC flag and colour based on CV%
        if cv_pct < t["excellent"]:
            flag   = "EXCELLENT"
            colour = "#2ecc71"  # Green
        elif cv_pct < t["acceptable"]:
            flag   = "ACCEPTABLE"
            colour = "#3498db"  # Blue
        elif cv_pct < t["poor"]:
            flag   = "HIGH VARIATION"
            colour = "#f39c12"  # Amber
        else:
            flag   = "VERY HIGH VARIATION"
            colour = "#e74c3c"  # Red

        summary_rows.append({
            "Group":            group_label,
            "N_Samples":        len(group_data),
            "Mean_Variants":    round(mean_count, 1),
            "Median_Variants":  round(group_data["N_Variants"].median(), 1),
            "Std_Variants":     round(std_count, 1),
            "CV_%":             round(cv_pct, 1),
            "Min_Variants":     group_data["N_Variants"].min(),
            "Max_Variants":     group_data["N_Variants"].max(),
            "QC_Flag":          flag,
            "_colour":          colour,   # Internal use — removed before CSV export
        })

    batch_summary = pd.DataFrame(summary_rows)

    # Print a clean summary table to the log
    logger.info(
        "\nBatch effect summary (variant count CV% per group):\n" +
        batch_summary[["Group", "N_Samples", "Mean_Variants", "CV_%", "QC_Flag"]]
        .to_string(index=False)
    )

    high_variation_groups = batch_summary[batch_summary["CV_%"] > t["acceptable"]]
    if not high_variation_groups.empty:
        logger.warning(
            f"{len(high_variation_groups)} group(s) have CV% > {t['acceptable']}% — "
            "possible batch effects or sample outliers. Review these groups before "
            "downstream association analyses."
        )

    return batch_summary


# ══════════════════════════════════════════════════════════════════════════════
# QC MODULE 4: SUMMARY FIGURE
# ══════════════════════════════════════════════════════════════════════════════

def plot_qc_summary(
    titv: dict,
    df: pd.DataFrame,
    batch_df: pd.DataFrame,
    out_dir: str,
    logger: logging.Logger,
) -> None:
    """
    Generate a three-panel quality control summary figure (PNG, 300 dpi).

    The figure provides a visual overview of all three QC modules at once,
    intended for inclusion in thesis reports or supplementary materials.

    Panel layout
    ────────────
    Panel 1 (left): Ti/Tv bar chart
        Two bars showing the absolute counts of transitions and transversions.
        The panel title displays the computed ratio and QC status.
        Reference thresholds are shown in an annotation box.

    Panel 2 (centre): Allelic balance histogram
        Distribution of allele fractions (AF) for all heterozygous calls.
        Colour-coded background zones show the optimal, acceptable, and poor
        AF ranges. A dashed line marks the expected AF of 0.5. The title
        reports the number and percentage of flagged (imbalanced) calls.

    Panel 3 (right): Batch effects horizontal bar chart
        One bar per cohort/tissue group, showing CV% of variant counts.
        Vertical reference lines mark the excellent, acceptable, and poor
        thresholds. Bars are colour-coded green/blue/amber/red.

    Parameters
    ----------
    titv     : dict            — output from calculate_titv()
    df       : pd.DataFrame    — variant table (must contain AB_AF, AB_Flag columns)
    batch_df : pd.DataFrame    — output from detect_batch_effects()
    out_dir  : str             — directory where the PNG will be saved
    logger   : logging.Logger  — for progress messages
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(
        "Variant-Level QC Summary",
        fontsize=14, fontweight="bold", y=1.01
    )

    # ── Panel 1: Ti/Tv bar chart ─────────────────────────────────────────────
    ax = axes[0]

    if not np.isnan(titv["ratio"]):
        # Draw one bar for transitions (Ti) and one for transversions (Tv)
        ax.bar(
            ["Transitions (Ti)", "Transversions (Tv)"],
            [titv["ti"], titv["tv"]],
            color=["#3498db", "#e74c3c"],
            edgecolor="black", linewidth=1.2
        )
        t = THRESHOLDS["titv"]

        # Panel title coloured to match QC status (green/amber/red/blue)
        ax.set_title(
            f"Ti/Tv Ratio: {titv['ratio']:.3f}\n{titv['status']}",
            fontweight="bold", color=titv["colour"], fontsize=11
        )

        # Small annotation box showing the reference thresholds
        ax.text(
            0.5, 0.02,
            f"Optimal {t['optimal_min']}–{t['optimal_max']}  |  "
            f"Acceptable {t['acceptable_min']}–{t['acceptable_max']}",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=8,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5)
        )
    else:
        # If ratio could not be computed, display the failure reason
        ax.text(
            0.5, 0.5, f"Ti/Tv\n{titv['status']}",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=12, color="red"
        )

    ax.set_ylabel("Count")
    ax.grid(axis="y", alpha=0.3)

    # ── Panel 2: Allelic balance histogram ───────────────────────────────────
    ax = axes[1]
    t_ab = THRESHOLDS["allelic_balance"]

    if "AB_AF" in df.columns and df["AB_AF"].notna().any():
        het_af = df["AB_AF"].dropna()

        # Draw colour-coded background zones to show QC regions
        ax.axvspan(0,                    t_ab["lower"],         alpha=0.15, color="#e74c3c", label="Poor")
        ax.axvspan(t_ab["lower"],        t_ab["optimal_lower"], alpha=0.15, color="#3498db", label="Acceptable")
        ax.axvspan(t_ab["optimal_lower"],t_ab["optimal_upper"], alpha=0.20, color="#2ecc71", label="Optimal")
        ax.axvspan(t_ab["optimal_upper"],t_ab["upper"],         alpha=0.15, color="#3498db")
        ax.axvspan(t_ab["upper"],        1.0,                   alpha=0.15, color="#e74c3c")

        ax.hist(het_af, bins=30, color="#9b59b6", edgecolor="black",
                alpha=0.75, zorder=3)
        ax.axvline(0.5, color="darkgreen", linestyle="--", linewidth=1.5,
                   label="Expected (0.5)", zorder=4)

        # Count flagged calls for the panel title
        n_flagged = (df["AB_Flag"].str.startswith("IMBALANCED") == True).sum()
        ax.set_title(
            f"Allelic Balance (heterozygous calls)\n"
            f"n={len(het_af)}  |  flagged={n_flagged} ({n_flagged/len(het_af)*100:.1f}%)",
            fontweight="bold", fontsize=11
        )
        ax.legend(loc="upper right", fontsize=7)
    else:
        ax.text(
            0.5, 0.5,
            "No allelic fraction data\n(AF column missing or all calls are homozygous)",
            transform=ax.transAxes, ha="center", va="center", fontsize=11
        )
        ax.set_title("Allelic Balance", fontweight="bold")

    ax.set_xlabel("Allele Fraction (AF)")
    ax.set_ylabel("Count")
    ax.set_xlim(0, 1)
    ax.grid(axis="y", alpha=0.3)

    # ── Panel 3: Batch effects horizontal bar chart ──────────────────────────
    ax = axes[2]

    if not batch_df.empty:
        t_cv    = THRESHOLDS["batch_cv"]
        colours = batch_df["_colour"].tolist()
        y_positions = range(len(batch_df))

        ax.barh(
            list(y_positions), batch_df["CV_%"],
            color=colours, edgecolor="black", linewidth=1.2, alpha=0.75
        )
        ax.set_yticks(list(y_positions))
        ax.set_yticklabels(batch_df["Group"], fontsize=9)

        # Vertical reference lines for each threshold level
        ax.axvline(t_cv["excellent"],  color="#2ecc71", linestyle="--",
                   linewidth=1.5, label=f"Excellent (< {t_cv['excellent']}%)")
        ax.axvline(t_cv["acceptable"], color="#f39c12", linestyle="--",
                   linewidth=1.5, label=f"Acceptable (< {t_cv['acceptable']}%)")
        ax.axvline(t_cv["poor"],       color="#e74c3c", linestyle="--",
                   linewidth=1.5, label=f"Poor (< {t_cv['poor']}%)")

        ax.set_xlabel("CV% of variant counts across samples")
        ax.set_title(
            "Batch Effects\n(variant count variability per group)",
            fontweight="bold", fontsize=11
        )
        ax.legend(loc="lower right", fontsize=7)
    else:
        ax.text(0.5, 0.5, "No batch data available",
                transform=ax.transAxes, ha="center", va="center", fontsize=12)
        ax.set_title("Batch Effects", fontweight="bold")

    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "Variant_QC_Plots.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"QC summary figure saved: {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """
    Orchestrate the full variant-level QC pipeline.

    Execution order
    ---------------
    1. Parse command-line arguments.
    2. Load the annotated variant table (Biological_Annotations sheet).
    3. Deduplicate transcript-expanded rows to one row per variant call.
    4. Calculate Ti/Tv ratio (Module 1).
    5. Check allelic balance for heterozygous calls (Module 2).
    6. Detect batch effects via CV% analysis (Module 3).
    7. Generate and save the three-panel QC figure (Module 4).
    8. Save all output files (Excel with AB flags, CSV summaries).
    9. Print a consolidated QC summary to the log.
    """

    # ── Command-line argument definition ────────────────────────────────────
    parser = argparse.ArgumentParser(
        description=(
            "Variant-level QC for GSDMB targeted sequencing data.\n"
            "Checks Ti/Tv ratio, allelic balance, and batch effects.\n"
            "Designed as step 07b in the GSDMB variant characterisation pipeline."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help=(
            "Path to the annotated Excel file produced by script 07 "
            "(typically named GSDMB_Annotated_Report_Fixed.xlsx). "
            "Must contain a sheet named 'Biological_Annotations'."
        ),
    )
    parser.add_argument(
        "--output", "-o",
        default=".",
        help="Directory where all output files will be saved (default: current directory).",
    )
    parser.add_argument(
        "--ab-lower",
        type=float,
        default=THRESHOLDS["allelic_balance"]["lower"],
        help=(
            f"Lower allele fraction threshold for allelic balance flagging. "
            f"Heterozygous calls with AF below this value are flagged as IMBALANCED. "
            f"Default: {THRESHOLDS['allelic_balance']['lower']}"
        ),
    )
    parser.add_argument(
        "--ab-upper",
        type=float,
        default=THRESHOLDS["allelic_balance"]["upper"],
        help=(
            f"Upper allele fraction threshold for allelic balance flagging. "
            f"Heterozygous calls with AF above this value are flagged as IMBALANCED. "
            f"Default: {THRESHOLDS['allelic_balance']['upper']}"
        ),
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose (DEBUG-level) logging for additional diagnostic detail.",
    )
    args = parser.parse_args()

    # ── Initialise logger and output directory ───────────────────────────────
    logger = setup_logging(args.verbose)
    os.makedirs(args.output, exist_ok=True)

    logger.info("=" * 65)
    logger.info("VARIANT-LEVEL QC  |  Script 07b  |  GSDMB pipeline")
    logger.info("=" * 65)
    logger.info(f"Input file   : {args.input}")
    logger.info(f"Output dir   : {args.output}")
    logger.info(f"AB thresholds: lower={args.ab_lower}, upper={args.ab_upper}")
    logger.info("=" * 65)

    # ── Step 1: Load annotated variant table ─────────────────────────────────
    try:
        df_raw = pd.read_excel(args.input, sheet_name="Biological_Annotations")
        logger.info(f"Loaded {len(df_raw)} rows from sheet 'Biological_Annotations'")
    except FileNotFoundError:
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Failed to load input file: {e}")
        sys.exit(1)

    # ── Step 2: Deduplicate to one row per variant call ───────────────────────
    # (Collapses VEP's per-transcript expansion back to per-variant)
    df = deduplicate_to_variants(df_raw, logger)

    # ── Step 3: Ti/Tv ratio ──────────────────────────────────────────────────
    logger.info("\n── MODULE 1: Ti/Tv RATIO ───────────────────────────────────")
    titv = calculate_titv(df, logger)

    # ── Step 4: Allelic balance ───────────────────────────────────────────────
    logger.info("\n── MODULE 2: ALLELIC BALANCE ───────────────────────────────")
    df = check_allelic_balance(df, logger, args.ab_lower, args.ab_upper)

    # ── Step 5: Batch effects ─────────────────────────────────────────────────
    logger.info("\n── MODULE 3: BATCH EFFECTS ─────────────────────────────────")
    batch_df = detect_batch_effects(df, logger)

    # ── Step 6: Generate QC figure ────────────────────────────────────────────
    logger.info("\n── MODULE 4: GENERATING PLOTS ──────────────────────────────")
    plot_qc_summary(titv, df, batch_df, args.output, logger)

    # ── Step 7: Save output files ─────────────────────────────────────────────
    logger.info("\n── SAVING OUTPUT FILES ─────────────────────────────────────")

    # Re-attach AB flag columns to the full (pre-deduplication) table.
    # This ensures that every transcript row for a flagged variant also carries
    # the allelic balance flag — important for downstream filtering in script 08.
    key_cols = [c for c in ["CHROM", "POS", "REF", "ALT", "Sample"]
                if c in df.columns]

    if key_cols and all(c in df_raw.columns for c in key_cols):
        # Extract only the AB flag columns plus variant keys, then merge
        ab_cols_to_merge = df[key_cols + ["AB_AF", "AB_Flag", "AB_Quality"]].drop_duplicates(
            subset=key_cols
        )
        df_output = df_raw.merge(ab_cols_to_merge, on=key_cols, how="left")
        logger.info(
            f"AB flags merged back to full table: {len(df_output)} rows "
            f"(one row per variant × transcript)"
        )
    else:
        # Fallback: use the deduplicated table directly
        df_output = df
        logger.warning(
            "Could not merge AB flags back to full table — "
            "key columns missing. Output will contain deduplicated rows only."
        )

    # Save annotated table as Excel
    excel_path = os.path.join(args.output, "Variant_QC_Report.xlsx")
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df_output.to_excel(writer, sheet_name="Biological_Annotations", index=False)
    logger.info(f"Annotated variant table saved: {excel_path}")

    # Save Ti/Tv summary as CSV
    titv_path = os.path.join(args.output, "TiTv_Summary.csv")
    pd.DataFrame([{k: v for k, v in titv.items() if k != "colour"}]).to_csv(
        titv_path, index=False
    )
    logger.info(f"Ti/Tv summary CSV saved: {titv_path}")

    # Save batch effects summary as CSV (remove internal colour column first)
    if not batch_df.empty:
        batch_path = os.path.join(args.output, "Batch_Effects_Summary.csv")
        batch_df.drop(columns=["_colour"], errors="ignore").to_csv(
            batch_path, index=False
        )
        logger.info(f"Batch effects summary CSV saved: {batch_path}")

    # ── Step 8: Print consolidated QC summary ────────────────────────────────
    logger.info("\n" + "=" * 65)
    logger.info("QC COMPLETE — Consolidated Summary")
    logger.info("=" * 65)

    ratio_str = f"{titv['ratio']:.3f}" if not np.isnan(titv["ratio"]) else "N/A"
    logger.info(f"  Ti/Tv ratio      : {ratio_str}  [{titv['status']}]")

    if "AB_Flag" in df.columns:
        n_het     = df["AB_AF"].notna().sum()
        n_flagged = df["AB_Flag"].str.startswith("IMBALANCED").sum()
        logger.info(
            f"  Allelic balance  : {n_flagged} / {n_het} heterozygous calls "
            f"flagged as imbalanced ({n_flagged/n_het*100:.1f}% of hets)"
            if n_het > 0 else
            "  Allelic balance  : no heterozygous calls found"
        )

    if not batch_df.empty:
        n_high_cv = (batch_df["CV_%"] > THRESHOLDS["batch_cv"]["acceptable"]).sum()
        logger.info(
            f"  Batch effects    : {n_high_cv} group(s) with CV% > "
            f"{THRESHOLDS['batch_cv']['acceptable']}% (acceptable threshold)"
        )

    logger.info("=" * 65)
    logger.info("Output files:")
    logger.info(f"  {os.path.join(args.output, 'Variant_QC_Report.xlsx')}")
    logger.info(f"  {os.path.join(args.output, 'TiTv_Summary.csv')}")
    logger.info(f"  {os.path.join(args.output, 'Batch_Effects_Summary.csv')}")
    logger.info(f"  {os.path.join(args.output, 'Variant_QC_Plots.png')}")
    logger.info("=" * 65)


# ── Script entry point ────────────────────────────────────────────────────────
# This block ensures that main() is only called when the script is run
# directly (e.g. `python 07b_variant-qc.py`), not when it is imported as
# a module by another script.
if __name__ == "__main__":
    main()
