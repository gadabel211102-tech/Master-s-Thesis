#!/usr/bin/env python3
"""
Variant-Level Quality Control (07b)
=====================================
Runs after script 07 (annotation merge) and before script 08 (mapping).

Checks three things that are genuinely answerable from the annotated variant table:

  1. Ti/Tv ratio        — is the SNV spectrum biologically plausible?
  2. Allelic balance    — are heterozygous calls real or artefactual?
  3. Batch effects      — do variant counts vary suspiciously across cohort/tissue groups?

NOTE on amplicon QC: this script intentionally does NOT attempt per-amplicon
coverage QC. Coverage QC requires the full depth track (mosdepth output), not
just the rows where a variant was called. Amplicons with no variants are
invisible here and would be silently omitted, giving a false picture of panel
performance. Use scripts 02/03/04 for amplicon-level coverage assessment.

Input:  GSDMB_Annotated_Report_Fixed.xlsx  (sheet: Biological_Annotations)
Output: <output_dir>/
          Variant_QC_Report.xlsx       — full table with AB flags appended
          TiTv_Summary.csv             — transition / transversion counts and ratio
          Batch_Effects_Summary.csv    — per-group variant count statistics
          Variant_QC_Plots.png         — three-panel summary figure

Usage:
  python 07b_variant_qc.py --input /path/to/GSDMB_Annotated_Report_Fixed.xlsx
  python 07b_variant_qc.py --input annotations.xlsx --output ./qc --ab-lower 0.2 --ab-upper 0.8
"""

import argparse
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# THRESHOLDS
# ─────────────────────────────────────────────────────────────────────────────

THRESHOLDS = {
    "titv": {
        "optimal_min": 2.0,
        "optimal_max": 2.1,
        "acceptable_min": 1.5,
        "acceptable_max": 3.0,
    },
    "allelic_balance": {
        "lower": 0.25,        # below this → flag as imbalanced
        "upper": 0.75,        # above this → flag as imbalanced
        "optimal_lower": 0.40,
        "optimal_upper": 0.60,
    },
    "batch_cv": {
        "excellent": 15,      # CV% below this is fine
        "acceptable": 30,
        "poor": 50,
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────────────────────────────────────

def setup_logging(verbose: bool = False) -> logging.Logger:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def find_col(df: pd.DataFrame, name: str):
    """Case-insensitive column lookup. Returns actual column name or None."""
    for c in df.columns:
        if c.upper() == name.upper():
            return c
    return None


def deduplicate_to_variants(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Script 07 produces one row per variant × transcript (VEP expands CSQ).
    For QC purposes we want one row per unique variant call (CHROM/POS/REF/ALT/Sample).
    Drop transcript duplicates so we don't inflate Ti/Tv counts or flag the same
    variant multiple times for allelic imbalance.
    """
    key_cols = [c for c in ["CHROM", "POS", "REF", "ALT", "Sample"] if find_col(df, c)]
    if not key_cols:
        logger.warning("Cannot identify variant key columns — using full table (may inflate counts)")
        return df
    # Normalise column names for the key
    rename = {find_col(df, c): c for c in key_cols}
    df = df.rename(columns=rename)
    before = len(df)
    df = df.drop_duplicates(subset=key_cols)
    logger.info(f"Deduplicated {before} transcript rows → {len(df)} unique variant calls")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 1. Ti/Tv RATIO
# ─────────────────────────────────────────────────────────────────────────────

TRANSITIONS  = {"A>G", "G>A", "C>T", "T>C"}
TRANSVERSIONS = {"A>C", "C>A", "A>T", "T>A", "G>C", "C>G", "G>T", "T>G"}


def calculate_titv(df: pd.DataFrame, logger: logging.Logger) -> dict:
    """
    Count transitions and transversions across all SNVs.
    Only single-nucleotide substitutions are included; indels are skipped.

    A healthy targeted panel should produce Ti/Tv ≈ 2.0–2.1.
    Lower values suggest poor-quality calls or technical artefacts.
    Higher values can indicate over-filtering or CpG site enrichment.
    """
    ref_col = find_col(df, "REF")
    alt_col = find_col(df, "ALT")

    if not ref_col or not alt_col:
        logger.warning("REF or ALT column missing — cannot calculate Ti/Tv")
        return {"ti": 0, "tv": 0, "snv_total": 0, "ratio": float("nan"),
                "status": "FAIL_MISSING_COLUMNS", "colour": "#e74c3c"}

    ti = tv = skipped = 0
    for ref, alt in zip(df[ref_col].astype(str).str.upper(),
                        df[alt_col].astype(str).str.upper()):
        if len(ref) != 1 or len(alt) != 1:
            skipped += 1   # indel
            continue
        if ref not in "ACGT" or alt not in "ACGT":
            skipped += 1
            continue
        mut = f"{ref}>{alt}"
        if mut in TRANSITIONS:
            ti += 1
        elif mut in TRANSVERSIONS:
            tv += 1

    snv_total = ti + tv
    logger.info(f"SNVs: {snv_total} (transitions={ti}, transversions={tv}, indels/other skipped={skipped})")

    if tv == 0:
        logger.warning("No transversions found — Ti/Tv cannot be calculated (too few variants?)")
        return {"ti": ti, "tv": tv, "snv_total": snv_total,
                "ratio": float("nan"), "status": "FAIL_NO_TV", "colour": "#e74c3c"}

    ratio = ti / tv
    t = THRESHOLDS["titv"]

    if ratio < t["acceptable_min"]:
        status, colour = "FAIL — ratio too low (likely artefacts or very few variants)", "#e74c3c"
    elif ratio > t["acceptable_max"]:
        status, colour = "WARN — ratio high (CpG enrichment or over-filtering?)", "#f39c12"
    elif t["optimal_min"] <= ratio <= t["optimal_max"]:
        status, colour = "OPTIMAL", "#2ecc71"
    else:
        status, colour = "ACCEPTABLE", "#3498db"

    logger.info(f"Ti/Tv ratio: {ratio:.3f}  →  {status}")
    return {"ti": ti, "tv": tv, "snv_total": snv_total,
            "ratio": ratio, "status": status, "colour": colour}


# ─────────────────────────────────────────────────────────────────────────────
# 2. ALLELIC BALANCE
# ─────────────────────────────────────────────────────────────────────────────

def check_allelic_balance(df: pd.DataFrame, logger: logging.Logger,
                          ab_lower: float, ab_upper: float) -> pd.DataFrame:
    """
    For every heterozygous call, check that the alternate-allele fraction
    sits near 0.5 (expected for a true germline het or a ~50% VAF somatic).

    Strong deviation flags:
      - Contamination (another sample's reads boosting one allele)
      - Loss of heterozygosity (one allele lost in tumour)
      - Copy number changes inflating one allele count
      - Strand/amplification artefacts

    Adds columns:
      AB_AF       — computed allele fraction (NaN for hom/missing)
      AB_Flag     — PASS | IMBALANCED | NO_DATA
      AB_Quality  — OPTIMAL | ACCEPTABLE | POOR | HOM/MISSING
    """
    t = THRESHOLDS["allelic_balance"]
    df = df.copy()
    df["AB_AF"]      = np.nan
    df["AB_Flag"]    = "NO_DATA"
    df["AB_Quality"] = "HOM/MISSING"

    gt_col = find_col(df, "GT")
    af_col = find_col(df, "AF")

    if not gt_col:
        logger.warning("GT column not found — skipping allelic balance check")
        return df

    het_idx = []
    for idx, gt_raw in df[gt_col].items():
        gt = str(gt_raw).replace("|", "/")
        if "/" not in gt:
            continue
        parts = gt.split("/")
        if len(parts) == 2 and parts[0] != parts[1] and "." not in parts:
            het_idx.append(idx)

    logger.info(f"Heterozygous calls to check: {len(het_idx)} / {len(df)}")

    if not het_idx:
        logger.warning("No heterozygous calls found — all homozygous or missing?")
        return df

    if not af_col:
        logger.warning("AF column not found — cannot compute allelic fractions")
        return df

    flagged = optimal = acceptable = 0
    for idx in het_idx:
        try:
            af = float(df.at[idx, af_col])
        except (ValueError, TypeError):
            continue

        df.at[idx, "AB_AF"] = af

        if af < ab_lower or af > ab_upper:
            df.at[idx, "AB_Flag"]    = f"IMBALANCED (AF={af:.3f})"
            df.at[idx, "AB_Quality"] = "POOR"
            flagged += 1
        elif t["optimal_lower"] <= af <= t["optimal_upper"]:
            df.at[idx, "AB_Flag"]    = "PASS"
            df.at[idx, "AB_Quality"] = "OPTIMAL"
            optimal += 1
        else:
            df.at[idx, "AB_Flag"]    = "PASS"
            df.at[idx, "AB_Quality"] = "ACCEPTABLE"
            acceptable += 1

    n = len(het_idx)
    pct_flagged = flagged / n * 100
    logger.info(f"Allelic balance — optimal: {optimal} ({optimal/n*100:.1f}%), "
                f"acceptable: {acceptable} ({acceptable/n*100:.1f}%), "
                f"imbalanced: {flagged} ({pct_flagged:.1f}%)")
    if pct_flagged > 20:
        logger.warning(f">{pct_flagged:.0f}% of hets are imbalanced — "
                       "check for contamination, LOH, or CNV in these samples")

    return df


# ─────────────────────────────────────────────────────────────────────────────
# 3. BATCH EFFECTS
# ─────────────────────────────────────────────────────────────────────────────

def detect_batch_effects(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Count variants per sample, then compute the coefficient of variation (CV%)
    within each cohort × tissue group.

    High CV suggests one of:
      - Batch effects (samples sequenced in different runs / reagent lots)
      - Genuine biology (some tumours highly mutated)
      - Sample quality outliers that should be excluded

    Returns a summary DataFrame, one row per cohort × tissue group.
    """
    cohort_col  = find_col(df, "Cohort")
    tissue_col  = find_col(df, "Tissue")
    sample_col  = find_col(df, "Sample")

    if not sample_col:
        logger.warning("Sample column not found — skipping batch effect detection")
        return pd.DataFrame()

    group_cols = [c for c in [cohort_col, tissue_col] if c]
    if not group_cols:
        logger.warning("No Cohort or Tissue column found — using single global group")
        df = df.copy()
        df["_group"] = "All"
        group_cols = ["_group"]

    per_sample = df.groupby(group_cols + [sample_col]).size().reset_index(name="N_Variants")

    t = THRESHOLDS["batch_cv"]
    rows = []
    for keys, grp in per_sample.groupby(group_cols):
        if not isinstance(keys, tuple):
            keys = (keys,)
        label = " / ".join(str(k) for k in keys)
        mean = grp["N_Variants"].mean()
        std  = grp["N_Variants"].std(ddof=1) if len(grp) > 1 else 0.0
        cv   = std / mean * 100 if mean > 0 else 0.0

        if cv < t["excellent"]:
            flag, colour = "EXCELLENT",      "#2ecc71"
        elif cv < t["acceptable"]:
            flag, colour = "ACCEPTABLE",     "#3498db"
        elif cv < t["poor"]:
            flag, colour = "HIGH VARIATION", "#f39c12"
        else:
            flag, colour = "VERY HIGH VAR",  "#e74c3c"

        rows.append({
            "Group":           label,
            "N_Samples":       len(grp),
            "Mean_Variants":   round(mean, 1),
            "Median_Variants": round(grp["N_Variants"].median(), 1),
            "Std_Variants":    round(std, 1),
            "CV_%":            round(cv, 1),
            "Min_Variants":    grp["N_Variants"].min(),
            "Max_Variants":    grp["N_Variants"].max(),
            "QC_Flag":         flag,
            "_colour":         colour,
        })

    batch_df = pd.DataFrame(rows)
    logger.info("\nBatch effect summary:\n" +
                batch_df[["Group", "N_Samples", "Mean_Variants", "CV_%", "QC_Flag"]]
                .to_string(index=False))

    high_var = batch_df[batch_df["CV_%"] > t["acceptable"]]
    if not high_var.empty:
        logger.warning(f"{len(high_var)} group(s) with CV > {t['acceptable']}% — "
                       "possible batch effects or sample outliers")
    return batch_df


# ─────────────────────────────────────────────────────────────────────────────
# 4. PLOTS
# ─────────────────────────────────────────────────────────────────────────────

def plot_qc_summary(titv: dict, df: pd.DataFrame,
                    batch_df: pd.DataFrame, out_dir: str,
                    logger: logging.Logger):
    """Three-panel QC figure: Ti/Tv bar, allelic balance histogram, batch CV."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Variant-Level QC Summary", fontsize=14, fontweight="bold", y=1.01)

    # ── Panel 1: Ti/Tv ──────────────────────────────────────────────────────
    ax = axes[0]
    if not np.isnan(titv["ratio"]):
        ax.bar(["Transitions", "Transversions"],
               [titv["ti"], titv["tv"]],
               color=["#3498db", "#e74c3c"], edgecolor="black", linewidth=1.2)
        t = THRESHOLDS["titv"]
        ax.set_title(
            f"Ti/Tv Ratio: {titv['ratio']:.3f}\n{titv['status']}",
            fontweight="bold", color=titv["colour"], fontsize=11)
        ax.text(0.5, 0.02,
                f"Optimal {t['optimal_min']}–{t['optimal_max']}  |  "
                f"Acceptable {t['acceptable_min']}–{t['acceptable_max']}",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=8,
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    else:
        ax.text(0.5, 0.5, f"Ti/Tv\n{titv['status']}",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=12, color="red")
    ax.set_ylabel("Count")
    ax.grid(axis="y", alpha=0.3)

    # ── Panel 2: Allelic balance histogram ──────────────────────────────────
    ax = axes[1]
    t_ab = THRESHOLDS["allelic_balance"]
    if "AB_AF" in df.columns and df["AB_AF"].notna().any():
        het_af = df["AB_AF"].dropna()

        # Colour-coded background zones
        ax.axvspan(0,             t_ab["lower"],         alpha=0.15, color="#e74c3c", label="Poor")
        ax.axvspan(t_ab["lower"], t_ab["optimal_lower"], alpha=0.15, color="#3498db", label="Acceptable")
        ax.axvspan(t_ab["optimal_lower"], t_ab["optimal_upper"], alpha=0.2, color="#2ecc71", label="Optimal")
        ax.axvspan(t_ab["optimal_upper"], t_ab["upper"], alpha=0.15, color="#3498db")
        ax.axvspan(t_ab["upper"], 1.0,                   alpha=0.15, color="#e74c3c")

        ax.hist(het_af, bins=30, color="#9b59b6", edgecolor="black", alpha=0.75, zorder=3)
        ax.axvline(0.5, color="darkgreen", linestyle="--", linewidth=1.5,
                   label="Expected (0.5)", zorder=4)

        n_flagged = (df["AB_Flag"].str.startswith("IMBALANCED") == True).sum()
        ax.set_title(
            f"Allelic Balance (heterozygous calls)\n"
            f"n={len(het_af)}  flagged={n_flagged} ({n_flagged/len(het_af)*100:.1f}%)",
            fontweight="bold", fontsize=11)
        ax.legend(loc="upper right", fontsize=7)
    else:
        ax.text(0.5, 0.5, "No allelic fraction data\n(AF column missing or all homozygous)",
                transform=ax.transAxes, ha="center", va="center", fontsize=11)
        ax.set_title("Allelic Balance", fontweight="bold")
    ax.set_xlabel("Allele Fraction (AF)")
    ax.set_ylabel("Count")
    ax.set_xlim(0, 1)
    ax.grid(axis="y", alpha=0.3)

    # ── Panel 3: Batch CV% ───────────────────────────────────────────────────
    ax = axes[2]
    if not batch_df.empty:
        t_cv = THRESHOLDS["batch_cv"]
        colours = batch_df["_colour"].tolist()
        y_pos = range(len(batch_df))
        ax.barh(list(y_pos), batch_df["CV_%"], color=colours,
                edgecolor="black", linewidth=1.2, alpha=0.75)
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels(batch_df["Group"], fontsize=9)
        ax.axvline(t_cv["excellent"],  color="#2ecc71", linestyle="--",
                   linewidth=1.5, label=f"Excellent (<{t_cv['excellent']}%)")
        ax.axvline(t_cv["acceptable"], color="#f39c12", linestyle="--",
                   linewidth=1.5, label=f"Acceptable (<{t_cv['acceptable']}%)")
        ax.axvline(t_cv["poor"],       color="#e74c3c", linestyle="--",
                   linewidth=1.5, label=f"Poor (<{t_cv['poor']}%)")
        ax.set_xlabel("CV% of variant counts across samples")
        ax.set_title("Batch Effects\n(variant count variability per group)",
                     fontweight="bold", fontsize=11)
        ax.legend(loc="lower right", fontsize=7)
    else:
        ax.text(0.5, 0.5, "No batch data", transform=ax.transAxes,
                ha="center", va="center", fontsize=12)
        ax.set_title("Batch Effects", fontweight="bold")
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "Variant_QC_Plots.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"QC plot saved: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Variant-level QC: Ti/Tv, allelic balance, batch effects",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input",  "-i", required=True,
                        help="Annotated Excel file from script 07 "
                             "(GSDMB_Annotated_Report_Fixed.xlsx)")
    parser.add_argument("--output", "-o", default=".",
                        help="Output directory (default: current directory)")
    parser.add_argument("--ab-lower", type=float,
                        default=THRESHOLDS["allelic_balance"]["lower"],
                        help="Lower AF threshold for allelic balance "
                             f"(default: {THRESHOLDS['allelic_balance']['lower']})")
    parser.add_argument("--ab-upper", type=float,
                        default=THRESHOLDS["allelic_balance"]["upper"],
                        help="Upper AF threshold for allelic balance "
                             f"(default: {THRESHOLDS['allelic_balance']['upper']})")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logger = setup_logging(args.verbose)
    os.makedirs(args.output, exist_ok=True)

    logger.info("=" * 65)
    logger.info("VARIANT-LEVEL QC  (script 07b)")
    logger.info("=" * 65)
    logger.info(f"Input : {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"AB thresholds: {args.ab_lower} – {args.ab_upper}")

    # ── Load ────────────────────────────────────────────────────────────────
    try:
        df_raw = pd.read_excel(args.input, sheet_name="Biological_Annotations")
        logger.info(f"Loaded {len(df_raw)} rows from Biological_Annotations")
    except Exception as e:
        logger.error(f"Failed to load input: {e}")
        sys.exit(1)

    # ── Deduplicate to one row per variant call ──────────────────────────────
    df = deduplicate_to_variants(df_raw, logger)

    # ── 1. Ti/Tv ────────────────────────────────────────────────────────────
    logger.info("\n── 1. Ti/Tv RATIO ──────────────────────────────────────────")
    titv = calculate_titv(df, logger)

    # ── 2. Allelic balance ───────────────────────────────────────────────────
    logger.info("\n── 2. ALLELIC BALANCE ──────────────────────────────────────")
    df = check_allelic_balance(df, logger, args.ab_lower, args.ab_upper)

    # ── 3. Batch effects ─────────────────────────────────────────────────────
    logger.info("\n── 3. BATCH EFFECTS ────────────────────────────────────────")
    batch_df = detect_batch_effects(df, logger)

    # ── 4. Plots ─────────────────────────────────────────────────────────────
    logger.info("\n── 4. PLOTS ────────────────────────────────────────────────")
    plot_qc_summary(titv, df, batch_df, args.output, logger)

    # ── 5. Save outputs ──────────────────────────────────────────────────────
    logger.info("\n── 5. SAVING OUTPUTS ───────────────────────────────────────")

    # Annotated table with AB flags
    # Re-attach AB columns to the full (non-deduplicated) table so every transcript
    # row for a flagged variant is also flagged.
    key_cols = [c for c in ["CHROM", "POS", "REF", "ALT", "Sample"] if c in df.columns]
    if key_cols and all(c in df_raw.columns for c in key_cols):
        ab_cols = df[key_cols + ["AB_AF", "AB_Flag", "AB_Quality"]].drop_duplicates(subset=key_cols)
        df_out = df_raw.merge(ab_cols, on=key_cols, how="left")
    else:
        df_out = df  # fallback

    excel_path = os.path.join(args.output, "Variant_QC_Report.xlsx")
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df_out.to_excel(writer, sheet_name="Biological_Annotations", index=False)
    logger.info(f"Annotated table saved: {excel_path}")

    # Ti/Tv summary
    titv_path = os.path.join(args.output, "TiTv_Summary.csv")
    pd.DataFrame([{k: v for k, v in titv.items() if k != "colour"}]).to_csv(
        titv_path, index=False)
    logger.info(f"Ti/Tv summary saved: {titv_path}")

    # Batch effects summary
    if not batch_df.empty:
        batch_path = os.path.join(args.output, "Batch_Effects_Summary.csv")
        batch_df.drop(columns=["_colour"], errors="ignore").to_csv(batch_path, index=False)
        logger.info(f"Batch effects summary saved: {batch_path}")

    # ── Final summary ────────────────────────────────────────────────────────
    logger.info("\n" + "=" * 65)
    logger.info("QC COMPLETE — Summary")
    logger.info("=" * 65)
    ratio_str = f"{titv['ratio']:.3f}" if not np.isnan(titv["ratio"]) else "N/A"
    logger.info(f"  Ti/Tv ratio  : {ratio_str}  [{titv['status']}]")
    if "AB_Flag" in df.columns:
        n_het     = df["AB_AF"].notna().sum()
        n_flagged = df["AB_Flag"].str.startswith("IMBALANCED").sum()
        logger.info(f"  Allelic bal. : {n_flagged}/{n_het} hets flagged as imbalanced")
    if not batch_df.empty:
        high_cv = (batch_df["CV_%"] > THRESHOLDS["batch_cv"]["acceptable"]).sum()
        logger.info(f"  Batch effects: {high_cv} group(s) with CV > "
                    f"{THRESHOLDS['batch_cv']['acceptable']}%")
    logger.info("=" * 65)


if __name__ == "__main__":
    main()