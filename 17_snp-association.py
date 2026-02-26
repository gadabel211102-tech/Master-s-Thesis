#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
17_snp-association.py  (v2 — column-verified against HARMONISED_B_v3)
=======================================================================
SNP–Clinical Variable Association Analysis

Joins the GSDMB variant report (GSDMB_Annotated_Report_Fixed.xlsx) with the
harmonised clinical master (MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx)
and tests each SNP for association with clinical variables.

SAMPLE FILTERING — SEQUENCED SAMPLES ONLY
-------------------------------------------
Analysis is restricted to samples confirmed as sequenced, derived from the
manifest sheets in the harmonised master.  Exclusion rules (verified from data):

  AT=AUs (AU endometrial):
    EXCLUDE extraction_flag == 'NO HAY'   — no biological material available
    EXCLUDE pd_status == 'NO HACER'       — explicitly flagged "do not process"
    → 114 / 118 DNA samples retained

  MT-T_N (Spanish breast tumour + paired normal):
    No exclusion flags present → all 71 DNA samples retained

  EN (Spanish endometrial healthy):
    No exclusion flags present → all 80 DNA samples retained

  MN (Spanish breast healthy):
    No exclusion flags present → all 100 DNA samples retained

The GSDMB variant report (inner join on snp_code) provides a second-level
filter: only snp_codes that produced variant calls will be present in that
file.  Together, these two filters ensure analysis is limited to samples
that have both been (a) approved for sequencing and (b) successfully run.

ANALYSES
--------
1. Tumour vs Healthy       — Fisher's exact test per cohort + globally
2. SNP × Clinical variable — Fisher's exact (binary/categorical) or
                              Mann–Whitney U (continuous), per tumour cohort
3. Genotype-dose           — Kruskal-Wallis / chi-square trend for top hits
4. Survival (AU endo only) — KM curves + log-rank + Cox PH (OS, PFS)

VARIABLE COVERAGE (verified from data)
---------------------------------------
BREAST (SP, MT-T_N sheet, n=71 tumour samples)
  ✓ Grade             clin_dcs__GRADO              (1/2/3, n=65)
  ✓ ER status         clin_dcs__RE                 (POSITIVO/NEGATIVO, n=65)
  ✓ PR status         clin_dcs__RP                 (POSITIVO/NEGATIVO, n=65)
  ✓ Recurrence        clin_dcs__Recaida/Progresión (NO vs other, n=65)
  ✓ Metastasis        clin_dcs__MTxDISTANCIA       (NO/NO-LOCAL/SI)
  ✓ HER2 copies       canon__her2_copies           (continuous, n=68)
  ✓ Age               canon__age                   (continuous)
  ✓ KI67              clin_dcs__KI67               (messy — coerced numeric)
  ✓ Exitus            clin_dcs__Exitus             (NO/SI)
  ✓ HER2+ subtype     clin_her2__DX                (Ca mama HER2+ / TN)
  ~ Survival (derived) from Fecha_dx + Ultima_fecha + Exitus

ENDOMETRIAL AU (AT=AUs sheet, n=114 sequenced DNA samples)
  ✓ FIGO stage        clin_au_endo__FIGO_STAGE     (IA/IB/II/IIIA/IIIC1/IIIC2/IVB)
  ✓ Grade             clin_au_endo__GRADE          (G1/G2/G3)
  ✓ Histology group   clin_au_endo__HISTOLOGY_GROUP (NEEC/EEC)
  ✓ LVSI              clin_au_endo__LVSI            (YES/NO)
  ✓ Myometrial inv.   clin_au_endo__MYOMETRIAL_INFILTRATION (<50%/>50%)
  ✓ MSI status        canon__msi_status            (Stable/UNSTABLE)
  ✓ Molecular class   canon__molecular_class       (P53/NSMP/MMRd/POLE)
  ✓ ER status         canon__er_status             (Positive/Negative)
  ✓ PR status         canon__pr_status             (Positive/Negative)
  ✓ Progression       clin_au_endo__PD_STATUS      (PD/NO PD)
  ✓ Exitus            clin_au_endo__EXITUS         (YES/NO)
  ✓ Risk recurrence   clin_au_endo__RISK_OF_RECURRENCE
  ✓ OS (months)       canon__os_months
  ✓ PFS (months)      canon__pfs_months
  ✓ Age               canon__age

OUTPUTS
-------
  17_SNP_Clinical_Association_Results.xlsx
    • tumour_vs_healthy
    • breast_clinical_assoc
    • endo_clinical_assoc
    • genotype_dose
    • survival_au_endo
    • summary_significant
    • sample_manifest

  17_SNP_Volcano_Plots_raw_p.png
  17_SNP_Volcano_Plots_FDR.png
  17_SNP_Heatmap_Breast_raw_p.png
  17_SNP_Heatmap_Breast_FDR.png
  17_SNP_Heatmap_Endometrial_raw_p.png
  17_SNP_Heatmap_Endometrial_FDR.png
  17_KM_Curves_AU_Endo.pdf

RUN
---
python3 17_snp-association.py \\
  --gsdmb   /path/to/GSDMB_Annotated_Report_Fixed.xlsx \\
  --master  /path/to/MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx \\
  --out_dir /path/to/output/
"""

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.backends.backend_pdf as pdf_backend
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import false_discovery_control, fisher_exact, mannwhitneyu

# lifelines optional — only needed for survival plots
try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    _HAS_LIFELINES = True
except ImportError:
    _HAS_LIFELINES = False
    print("WARNING: lifelines not installed — survival plots will be skipped.")

warnings.filterwarnings("ignore")

# ── DEFAULT PATHS ─────────────────────────────────────────────────────────────
DEFAULT_GSDMB  = Path("/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_Annotated_Report_Fixed.xlsx")
DEFAULT_MASTER = Path("/home/gadeaalonsoj/tfm/MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx")
DEFAULT_OUT    = Path("/home/gadeaalonsoj/tfm/gsdmb_final_results/")

# Pass-BAM manifests — one file per cohort, listing QC-passed BAMs
# These are the ground truth for which samples have actually been sequenced.
# Set to None to skip manifest filtering (falls back to extraction_flag logic).
DEFAULT_MANIFESTS = {
    "endometrium-tumour": Path("/home/gadeaalonsoj/tfm/manifests/endometrium-tumour-pass_manifest.txt"),
    "endometrium-normal": Path("/home/gadeaalonsoj/tfm/manifests/endometrium-normal-pass_manifest.txt"),
    "breast-tumour":      Path("/home/gadeaalonsoj/tfm/manifests/breast-tumour-pass_manifest.txt"),
    "breast-normal":      Path("/home/gadeaalonsoj/tfm/manifests/breast-normal-pass_manifest.txt"),
}

FDR_THRESHOLD  = 0.10
MIN_CARRIERS   = 3
MIN_NFE_AF     = 0.01


# ── MANIFEST PARSING ──────────────────────────────────────────────────────────

def parse_manifests(manifest_paths: Dict[str, Path]) -> Tuple[set, set]:
    """
    Parse pass-BAM manifest files and return:
      sequenced_codes : set of snp_codes confirmed sequenced (excl. replicates)
      replicate_codes : set of snp_codes that are technical Repeticion replicates

    Manifest format: one BAM path per line, e.g.:
      /path/to/pass_bams/DNA_SNP_MT-T_10_IonCode_0108_MamaHER2_LIB1y2_874.bam
      /path/to/pass_bams/DNA_MT_T_17_Repeticion_IonCode_0119.bam

    SNP code extraction mirrors the logic in _extract_snp_code().
    All four manifests are parsed and codes pooled across cohorts.
    """
    import os

    sequenced: set = set()
    replicates: set = set()

    for cohort_label, manifest_path in manifest_paths.items():
        if manifest_path is None or not Path(manifest_path).exists():
            print(f"      WARNING: manifest not found for {cohort_label}: {manifest_path}")
            continue

        with open(manifest_path) as fh:
            lines = [l.strip() for l in fh if l.strip()]

        for bam_path in lines:
            basename = os.path.basename(bam_path).replace(".bam", "")
            code = _extract_snp_code(basename)
            if code is None:
                print(f"      WARNING: could not parse snp_code from: {bam_path}")
                continue
            is_rep = bool(re.search(r"repeticion", basename, re.IGNORECASE))
            if is_rep:
                replicates.add(code)
            else:
                sequenced.add(code)

    print(f"      Manifest: {len(sequenced)} unique sequenced codes, "
          f"{len(replicates)} technical replicate codes")
    return sequenced, replicates


# ── SNP CODE EXTRACTION ───────────────────────────────────────────────────────

def _extract_snp_code(sample_name: str) -> Optional[str]:
    """
    Derive the SNP_XX_NN lab code from the sequencer sample name.
    Examples:
      SNP_AT_2__BL_039__IonCode_0118     → SNP_AT_2
      DNA_SNP_EN_10_IonCode_0101_...     → SNP_EN_10
      DNA_MT-T_37_IonCode_0121           → SNP_MT-T_37
    """
    s = str(sample_name)
    m = re.match(r"^(SNP_AT_\d+)", s)
    if m: return m.group(1)
    m = re.match(r"^DNA_SNP_(EN|MN|MT-T)_(\d+)_", s)
    if m: return f"SNP_{m.group(1)}_{m.group(2)}"
    m = re.match(r"^SNP_DNA_(MN|EN)_(\d+)_", s)
    if m: return f"SNP_{m.group(1)}_{m.group(2)}"
    m = re.match(r"^DNA_MT[-_]T_(\d+)_", s)
    if m: return f"SNP_MT-T_{m.group(1)}"
    return None


# ── CLINICAL VARIABLE DERIVATIONS ─────────────────────────────────────────────

def _derive_breast_recurrence(val) -> Optional[int]:
    """1 if any recurrence/progression text, 0 if 'NO', else NaN."""
    if pd.isna(val): return None
    return 0 if str(val).strip().upper() == "NO" else 1


def _derive_breast_exitus(val) -> Optional[int]:
    """1 if SI/YES, 0 if NO."""
    if pd.isna(val): return None
    v = str(val).strip().upper()
    if v.startswith("SI") or v.startswith("YES"): return 1
    if v.startswith("NO"): return 0
    return None


def _derive_breast_metastasis(val) -> Optional[int]:
    """1 = distant metastasis (SI), 0 = NO or local only."""
    if pd.isna(val): return None
    v = str(val).strip().upper()
    if v == "SI": return 1
    if v in ("NO", "NO-LOCAL"): return 0
    return None


def _derive_her2_subtype(val) -> Optional[str]:
    """Collapse DX field to HER2+, TN, or Other."""
    if pd.isna(val): return None
    v = str(val).strip().upper()
    if "HER2+" in v or "HER2 +" in v: return "HER2+"
    if "TN" in v or "TRIPLE" in v:     return "TN"
    return "Other"


def _derive_figo_numeric(val) -> Optional[float]:
    """Convert FIGO stage string to numeric for ordinal trend tests."""
    if pd.isna(val): return None
    mapping = {
        "IA": 1.0, "IB": 1.5, "I": 1.5,
        "II": 2.0,
        "IIIA": 3.1, "IIIB": 3.2, "IIIC1": 3.3, "IIIC2": 3.4, "III": 3.0,
        "IVA": 4.1, "IVB": 4.2, "IV": 4.0,
    }
    return mapping.get(str(val).strip().upper(), None)


def _derive_endo_grade(val) -> Optional[float]:
    """G1→1, G2→2, G3→3."""
    if pd.isna(val): return None
    m = re.search(r"(\d)", str(val))
    return float(m.group(1)) if m else None


def _derive_risk_ordinal(val) -> Optional[float]:
    """LOW=1, INTERMEDIATE=2, INTERMEDIATE-HIGH=3, HIGH=4."""
    if pd.isna(val): return None
    v = str(val).strip().upper()
    mapping = {
        "LOW": 1.0, "INTERMEDIATE": 2.0,
        "INTERMEDIATE-HIGH": 3.0, "HIGH": 4.0,
    }
    return mapping.get(v, None)


def _derive_ki67_numeric(val) -> Optional[float]:
    """Coerce messy KI67 to float (handles 0.3, '20-25%', 40 etc.)."""
    if pd.isna(val): return None
    s = str(val).replace("%", "").strip()
    # Range: take midpoint
    m = re.match(r"([\d.]+)\s*[-–]\s*([\d.]+)", s)
    if m:
        mid = (float(m.group(1)) + float(m.group(2))) / 2
        # If >1, assume percentage; convert to fraction
        return mid / 100 if mid > 1 else mid
    try:
        v = float(s)
        return v / 100 if v > 1 else v
    except (ValueError, TypeError):
        return None


def _derive_breast_os(row) -> Optional[float]:
    """Derive OS in months from fecha_dx + ultima_fecha + exitus."""
    try:
        dx  = pd.to_datetime(row["clin_dcs__Fecha dx"],    errors="coerce")
        ult = pd.to_datetime(row["clin_dcs__Última fecha disponible"], errors="coerce")
        if pd.isna(dx) or pd.isna(ult): return None
        return max((ult - dx).days / 30.44, 0)
    except Exception:
        return None


# ── HELPERS ───────────────────────────────────────────────────────────────────

def _haldane_or(a, b, c, d) -> Tuple[float, str]:
    if 0 in (a, b, c, d):
        a, b, c, d = a+0.5, b+0.5, c+0.5, d+0.5
        note = "Corrected (+0.5)"
    else:
        note = "Standard"
    denom = b * c
    return ((a * d) / denom if denom else np.inf), note


def _apply_fdr(df: pd.DataFrame, p_col: str = "P_Value",
               out_col: str = "FDR_P_Value") -> pd.DataFrame:
    df = df.copy()
    valid = df[p_col].notna()
    if valid.sum() > 1:
        df.loc[valid, out_col] = false_discovery_control(
            df.loc[valid, p_col].values, method="bh"
        )
    else:
        df[out_col] = df[p_col]
    return df


def _normalise_gt(gt) -> str:
    return "0/1" if str(gt) == "1/0" else str(gt)


# ── DATA LOADING & MERGING ────────────────────────────────────────────────────

def load_and_merge(gsdmb_path: Path, master_path: Path,
                   manifest_paths: Optional[Dict[str, Path]] = None) -> pd.DataFrame:
    """
    Load GSDMB variant report + clinical master; join on snp_code.
    Adds derived clinical columns with _DERIVED suffix.

    manifest_paths : dict mapping cohort label → path to pass-BAM manifest txt.
                     When provided, the clinical master is pre-filtered to only
                     snp_codes confirmed sequenced in those manifests.
                     Replicate BAMs (Repeticion) are retained but flagged.
    """
    print("[1/4] Loading GSDMB variant report …")
    gsdmb = pd.read_excel(gsdmb_path, sheet_name="Biological_Annotations")
    gsdmb["snp_code"] = gsdmb["Sample"].apply(_extract_snp_code)
    n_fail = gsdmb["snp_code"].isna().sum()
    if n_fail: print(f"      WARNING: could not extract snp_code for {n_fail} rows")

    if "CANONICAL" in gsdmb.columns:
        canon = gsdmb[gsdmb["CANONICAL"].astype(str).str.upper() == "YES"]
        if len(canon): gsdmb = canon.copy()

    if "GT" in gsdmb.columns:
        gsdmb["GT"] = gsdmb["GT"].apply(_normalise_gt)

    gsdmb["Tissue"] = gsdmb["Tissue"].astype(str).str.strip().replace(
        {"Normal": "Healthy", "Control": "Healthy", "Tumor": "Tumour"}
    )
    gsdmb["rsID"] = gsdmb["Existing_variation"].astype(str).str.extract(r"(rs\d+)")
    gsdmb["Variant_ID"] = gsdmb["rsID"].fillna(
        gsdmb["SYMBOL"].astype(str) + ":" + gsdmb["HGVSp"].astype(str)
    )
    gsdmb["gnomADe_NFE_AF"] = pd.to_numeric(gsdmb["gnomADe_NFE_AF"], errors="coerce")
    gsdmb = gsdmb[gsdmb["gnomADe_NFE_AF"] > MIN_NFE_AF].copy()
    print(f"      {len(gsdmb)} rows after gnomAD NFE AF > {MIN_NFE_AF} filter")

    print("[2/4] Loading harmonised clinical master …")
    master = pd.read_excel(master_path, sheet_name="harmonised_plus_canon")
    master_dna = master[master["nucleic_acid"] == "DNA"].copy()

    # ── MANIFEST FILTER: restrict to QC-passed sequenced samples ────────────
    if manifest_paths:
        print("      Applying pass-BAM manifest filter …")
        sequenced_codes, replicate_codes = parse_manifests(manifest_paths)

        # Flag replicates before filtering (kept but marked)
        master_dna["is_replicate_manifest"] = master_dna["snp_code"].isin(replicate_codes)

        n_before = len(master_dna)
        # Keep only codes confirmed sequenced (include replicates — flagged separately)
        all_manifest_codes = sequenced_codes | replicate_codes
        master_dna = master_dna[master_dna["snp_code"].isin(all_manifest_codes)].copy()
        n_after = len(master_dna)

        print(f"      Master rows: {n_before} → {n_after} (kept only manifest-confirmed codes)")
        print(f"      Sequenced samples per cohort:")
        for sheet, label in [("AT=AUs", "endo-tumour AU"), ("MT-T_N", "breast-tumour SP"),
                               ("EN", "endo-normal SP"), ("MN", "breast-normal SP")]:
            n = (master_dna["sheet"] == sheet).sum()
            n_rep = master_dna[(master_dna["sheet"] == sheet) & master_dna["is_replicate_manifest"]].shape[0]
            print(f"        {sheet:8s} ({label}): {n - n_rep} sequenced"
                  + (f" + {n_rep} replicates" if n_rep else ""))
    else:
        # Fallback: use extraction_flag from harmonised master (less precise)
        print("      No manifests provided — falling back to extraction_flag filter …")
        master_dna["is_replicate_manifest"] = False
        at_aus_excluded = (
            (master_dna["sheet"] == "AT=AUs")
            & (master_dna["extraction_flag"].isin(["NO HAY"]) | (master_dna["pd_status"] == "NO HACER"))
        )
        master_dna = master_dna[~at_aus_excluded].copy()

    # Deduplicate to one row per patient
    master_dna = master_dna.drop_duplicates("case_id").copy()

    # ── Derived columns ──────────────────────────────────────────────────────
    # Breast
    master_dna["BREAST_RECURRENCE_DERIVED"] = master_dna["clin_dcs__Recaida/Progresión"].apply(_derive_breast_recurrence)
    master_dna["BREAST_EXITUS_DERIVED"]     = master_dna["clin_dcs__Exitus"].apply(_derive_breast_exitus)
    master_dna["BREAST_METASTASIS_DERIVED"] = master_dna["clin_dcs__MTxDISTANCIA"].apply(_derive_breast_metastasis)
    master_dna["BREAST_HER2_SUBTYPE"]       = master_dna["clin_her2__DX"].apply(_derive_her2_subtype)
    master_dna["BREAST_ER_BIN"]             = master_dna["clin_dcs__RE"].map({"POSITIVO": 1, "NEGATIVO": 0})
    master_dna["BREAST_PR_BIN"]             = master_dna["clin_dcs__RP"].map({"POSITIVO": 1, "NEGATIVO": 0})
    master_dna["BREAST_KI67_NUMERIC"]       = master_dna["clin_dcs__KI67"].apply(_derive_ki67_numeric)
    master_dna["BREAST_OS_MONTHS_DERIVED"]  = master_dna.apply(_derive_breast_os, axis=1)

    # Endometrial AU
    master_dna["ENDO_FIGO_NUMERIC"]         = master_dna["clin_au_endo__FIGO_STAGE"].apply(_derive_figo_numeric)
    master_dna["ENDO_GRADE_NUMERIC"]        = master_dna["clin_au_endo__GRADE"].apply(_derive_endo_grade)
    master_dna["ENDO_LVSI_BIN"]             = master_dna["clin_au_endo__LVSI"].map({"YES": 1, "NO": 0})
    master_dna["ENDO_MYOINV_BIN"]           = master_dna["clin_au_endo__MYOMETRIAL_INFILTRATION"].map({"<50%": 0, ">50%": 1})
    master_dna["ENDO_MSI_BIN"]              = master_dna["canon__msi_status"].map({"UNSTABLE": 1, "Stable": 0})
    master_dna["ENDO_NEEC_BIN"]             = master_dna["clin_au_endo__HISTOLOGY_GROUP"].map({"NEEC": 1, "EEC": 0})
    master_dna["ENDO_PD_BIN"]               = master_dna["clin_au_endo__PD_STATUS"].map({"PD": 1, "NO PD": 0})
    master_dna["ENDO_EXITUS_BIN"]           = master_dna["clin_au_endo__EXITUS"].map({"YES": 1, "NO": 0})
    master_dna["ENDO_RISK_ORDINAL"]         = master_dna["clin_au_endo__RISK_OF_RECURRENCE"].apply(_derive_risk_ordinal)
    master_dna["ENDO_ER_BIN"]               = master_dna["canon__er_status"].map({"Positive": 1, "Negative": 0})
    master_dna["ENDO_PR_BIN"]               = master_dna["canon__pr_status"].map({"Positive": 1, "Negative": 0})

    # Keep only relevant columns
    canon_cols = [c for c in master_dna.columns
                  if c.startswith("canon__") and not c.endswith("__source")]
    derived_cols = [c for c in master_dna.columns if c.endswith("_DERIVED") or c.endswith("_BIN")
                    or c.endswith("_NUMERIC") or c.endswith("_ORDINAL") or c.endswith("_SUBTYPE")]
    raw_breast_cols = ["clin_dcs__GRADO", "clin_dcs__RE", "clin_dcs__RP",
                       "clin_her2__DX", "clin_dcs__Exitus",
                       "clin_dcs__Recaida/Progresión", "clin_dcs__MTxDISTANCIA",
                       "clin_dcs__KI67", "clin_dcs__Fecha dx",
                       "clin_dcs__Última fecha disponible"]
    raw_endo_cols  = ["clin_au_endo__FIGO_STAGE", "clin_au_endo__GRADE",
                       "clin_au_endo__HISTOLOGY_GROUP", "clin_au_endo__LVSI",
                       "clin_au_endo__MYOMETRIAL_INFILTRATION",
                       "clin_au_endo__MSI_STATUS_IHC", "clin_au_endo__RISK_OF_RECURRENCE",
                       "clin_au_endo__PD_STATUS", "clin_au_endo__EXITUS",
                       "clin_au_endo__MOLECULAR CLASSIFICATION_according to IHC and/or NGS profile"]
    keep = (["snp_code", "case_id", "sample_id", "sheet", "tissue", "tumour_normal",
              "histology", "pd_status"]
            + canon_cols + derived_cols
            + [c for c in raw_breast_cols + raw_endo_cols if c in master_dna.columns])
    master_slim = master_dna[[c for c in dict.fromkeys(keep) if c in master_dna.columns]].copy()

    print("[3/4] Merging on snp_code …")
    merged = gsdmb.merge(master_slim, on="snp_code", how="inner")
    merged["is_replicate"] = (
        merged["Sample"].str.contains("Repeticion", case=False, na=False)
        | merged.get("is_replicate_manifest", pd.Series(False, index=merged.index))
    )
    print(f"      {len(merged)} rows | {merged['Sample'].nunique()} samples | {merged['Variant_ID'].nunique()} variants")
    print("[4/4] Done.\n")
    return merged


# ── ANALYSIS 1: TUMOUR VS HEALTHY ────────────────────────────────────────────

def tumour_vs_healthy(df: pd.DataFrame) -> pd.DataFrame:
    print("=== Analysis 1: Tumour vs Healthy ===")
    df = df[~df["is_replicate"]].copy()
    cohort_groups = {
        "Global":        df,
        "Breast":        df[df["Cohort"].str.contains("Breast",    case=False, na=False)],
        "Endometrium":   df[df["Cohort"].str.contains("Endometri", case=False, na=False)],
    }
    rows = []
    for grp_name, gdf in cohort_groups.items():
        if gdf.empty: continue
        n_t = gdf[gdf["Tissue"] == "Tumour"]["Sample"].nunique()
        n_h = gdf[gdf["Tissue"] == "Healthy"]["Sample"].nunique()
        if not n_t or not n_h: continue
        for var_id, vdf in gdf.groupby("Variant_ID"):
            a = vdf[vdf["Tissue"] == "Tumour"]["Sample"].nunique()
            c = vdf[vdf["Tissue"] == "Healthy"]["Sample"].nunique()
            b, d = n_t - a, n_h - c
            if a + c < MIN_CARRIERS: continue
            _, p = fisher_exact([[a, b], [c, d]])
            or_v, meth = _haldane_or(a, b, c, d)
            sym = vdf["SYMBOL"].dropna().iloc[0] if vdf["SYMBOL"].notna().any() else ""
            rows.append({
                "Analysis_Group": grp_name,
                "Variant_ID": var_id, "Gene": sym,
                "Consequence": vdf["Consequence"].dropna().iloc[0] if "Consequence" in vdf and vdf["Consequence"].notna().any() else "",
                "IMPACT": vdf["IMPACT"].dropna().iloc[0] if vdf["IMPACT"].notna().any() else "",
                "gnomAD_NFE_AF": vdf["gnomADe_NFE_AF"].dropna().iloc[0] if vdf["gnomADe_NFE_AF"].notna().any() else np.nan,
                "N_Tumour": n_t, "N_Healthy": n_h,
                "Carriers_Tumour": a, "Carriers_Healthy": c,
                "Freq_Tumour_%": round(a / n_t * 100, 2),
                "Freq_Healthy_%": round(c / n_h * 100, 2),
                "Odds_Ratio": round(or_v, 4), "P_Value": p, "OR_Method": meth,
            })
    res = pd.DataFrame(rows)
    if res.empty: print("  No results.\n"); return res
    parts = [_apply_fdr(res[res["Analysis_Group"] == g].copy()) for g in res["Analysis_Group"].unique()]
    res = pd.concat(parts, ignore_index=True).sort_values(["Analysis_Group", "P_Value"])
    res["Nominal_Sig"] = res["P_Value"] < 0.05
    res["FDR_Sig"]     = res["FDR_P_Value"] < FDR_THRESHOLD
    print(f"  {len(res)} tests | {res['Nominal_Sig'].sum()} nominal | {res['FDR_Sig'].sum()} FDR\n")
    return res


# ── ANALYSIS 2: SNP × CLINICAL VARIABLE ─────────────────────────────────────
#
# Each cohort has completely different outcome variables from separate clinical systems:
#
#   BREAST (MT-T_N) — Spanish HER2 cohort (DCs MAMA HER2 clinical sheet)
#     Outcomes: grade, ER/PR, recurrence, metastasis, exitus, HER2 copies, KI67, OS, subtype
#     Confounder for adjustment: age at diagnosis
#
#   ENDOMETRIAL AU (AT=AUs) — Australian cohort (Clinical DATA AU_Endometrial)
#     Outcomes: FIGO stage, grade, myometrial invasion, LVSI, MSI, molecular class,
#               histology, ER/PR, progression, exitus, OS/PFS, risk of recurrence
#     Confounder for adjustment: age at surgery
#
# For BINARY outcomes: both unadjusted (Fisher exact) and age-adjusted
#   (logistic regression: outcome ~ snp_carrier + age_z) are run.
#   The adjusted model is skipped automatically if fewer than MIN_EVENTS_FOR_LOGISTIC
#   events exist in the smaller group.
# For CONTINUOUS outcomes: Mann-Whitney U (unadjusted only).
# For NOMINAL outcomes: Chi-square (unadjusted only).

MIN_EVENTS_FOR_LOGISTIC = 5   # minimum events in smaller group for adjusted model

CLINICAL_VARS_BREAST: Dict[str, Dict] = {
    # ── Continuous (Mann-Whitney U) ───────────────────────────────────────
    "BREAST_KI67_NUMERIC":      {"type": "continuous", "label": "KI67",
                                  "note": "Proliferation index (coerced to fraction)"},
    "canon__her2_copies":       {"type": "continuous", "label": "HER2 copies (FISH)",
                                  "note": "Continuous copy number"},
    "clin_dcs__GRADO":          {"type": "continuous", "label": "Tumour grade",
                                  "note": "Ordinal 1/2/3"},
    "BREAST_OS_MONTHS_DERIVED": {"type": "continuous", "label": "Overall survival (months)",
                                  "note": "Derived from fecha_dx + ultima_fecha"},
    "canon__age":               {"type": "continuous", "label": "Age at diagnosis"},
    # ── Binary (Fisher exact + age-adjusted logistic) ─────────────────────
    "BREAST_ER_BIN":            {"type": "binary", "label": "ER positive",
                                  "note": "POSITIVO=1 vs NEGATIVO=0"},
    "BREAST_PR_BIN":            {"type": "binary", "label": "PR positive",
                                  "note": "POSITIVO=1 vs NEGATIVO=0"},
    "BREAST_RECURRENCE_DERIVED":{"type": "binary", "label": "Recurrence / progression",
                                  "note": "Any recurrence=1 vs NO=0"},
    "BREAST_METASTASIS_DERIVED":{"type": "binary", "label": "Distant metastasis",
                                  "note": "SI=1 vs NO/NO-LOCAL=0"},
    "BREAST_EXITUS_DERIVED":    {"type": "binary", "label": "Exitus",
                                  "note": "SI=1 vs NO=0"},
    # ── Nominal (chi-square) ─────────────────────────────────────────────
    "BREAST_HER2_SUBTYPE":      {"type": "nominal", "label": "HER2 subtype",
                                  "note": "HER2+ / TN / Other"},
}

CLINICAL_VARS_ENDO: Dict[str, Dict] = {
    # ── Continuous (Mann-Whitney U) ───────────────────────────────────────
    "ENDO_FIGO_NUMERIC":        {"type": "continuous", "label": "FIGO stage",
                                  "note": "Ordinal: IA=1, IB=1.5, II=2, IIIA=3.1 … IVB=4.2"},
    "ENDO_GRADE_NUMERIC":       {"type": "continuous", "label": "Tumour grade",
                                  "note": "G1=1, G2=2, G3=3"},
    "ENDO_RISK_ORDINAL":        {"type": "continuous", "label": "Risk of recurrence",
                                  "note": "LOW=1, INT=2, INT-HIGH=3, HIGH=4"},
    "canon__os_months":         {"type": "continuous", "label": "Overall survival (months)"},
    "canon__pfs_months":        {"type": "continuous", "label": "Progression-free survival (months)"},
    "canon__age":               {"type": "continuous", "label": "Age at surgery"},
    # ── Binary (Fisher exact + age-adjusted logistic) ─────────────────────
    "ENDO_LVSI_BIN":            {"type": "binary", "label": "LVSI",
                                  "note": "YES=1 vs NO=0"},
    "ENDO_MYOINV_BIN":          {"type": "binary", "label": "Myometrial invasion ≥50%",
                                  "note": ">50%=1 vs <50%=0"},
    "ENDO_MSI_BIN":             {"type": "binary", "label": "MSI-H",
                                  "note": "UNSTABLE=1 vs Stable=0"},
    "ENDO_NEEC_BIN":            {"type": "binary", "label": "Non-endometrioid histology",
                                  "note": "NEEC=1 vs EEC=0"},
    "ENDO_PD_BIN":              {"type": "binary", "label": "Disease progression",
                                  "note": "PD=1 vs NO PD=0"},
    "ENDO_EXITUS_BIN":          {"type": "binary", "label": "Exitus",
                                  "note": "YES=1 vs NO=0"},
    "ENDO_ER_BIN":              {"type": "binary", "label": "ER positive",
                                  "note": "Positive=1 vs Negative=0"},
    "ENDO_PR_BIN":              {"type": "binary", "label": "PR positive",
                                  "note": "Positive=1 vs Negative=0"},
    # ── Nominal (chi-square) ─────────────────────────────────────────────
    "canon__molecular_class":   {"type": "nominal", "label": "Molecular classification",
                                  "note": "POLE / MMRd / NSMP / P53"},
    "clin_au_endo__FIGO_STAGE": {"type": "nominal", "label": "FIGO stage (categorical)",
                                  "note": "IA/IB/II/IIIA/IIIC1/IIIC2/IVB"},
}


def _test_continuous(c_vals, nc_vals) -> Optional[Dict]:
    """Mann-Whitney U between carriers and non-carriers. No age adjustment."""
    c  = pd.to_numeric(c_vals,  errors="coerce").dropna()
    nc = pd.to_numeric(nc_vals, errors="coerce").dropna()
    if len(c) < MIN_CARRIERS or len(nc) < MIN_CARRIERS:
        return None
    stat, p = mannwhitneyu(c.values, nc.values, alternative="two-sided")
    return {
        "Test_Unadj":          "Mann-Whitney U",
        "N_Carriers":          len(c),
        "N_NonCarriers":       len(nc),
        "Median_Carriers":     round(c.median(), 3),
        "Median_NonCarriers":  round(nc.median(), 3),
        "Stat_Unadj":          round(stat, 3),
        "P_Unadj":             p,
        "Test_Adj":            "N/A",
        "OR_Adj":              np.nan,
        "OR_Adj_CI95":         np.nan,
        "P_Adj":               np.nan,
        "N_Adj":               np.nan,
        "Adj_Note":            "Not applicable for continuous outcome",
    }


def _test_binary(c_df, nc_df, outcome_col, age_col="canon__age") -> Optional[Dict]:
    """
    Unadjusted Fisher exact test + age-adjusted logistic regression.
    outcome_col must already be encoded as 0/1.
    """
    try:
        from statsmodels.formula.api import logit as sm_logit
    except ImportError:
        sm_logit = None

    combined = pd.concat([
        c_df[[outcome_col, age_col]].assign(snp_carrier=1),
        nc_df[[outcome_col, age_col]].assign(snp_carrier=0),
    ])
    combined[outcome_col] = pd.to_numeric(combined[outcome_col], errors="coerce")
    combined = combined.dropna(subset=[outcome_col])

    n_c  = int((combined["snp_carrier"] == 1).sum())
    n_nc = int((combined["snp_carrier"] == 0).sum())
    if n_c < MIN_CARRIERS or n_nc < MIN_CARRIERS:
        return None

    # Counts for contingency table
    a = int(((combined["snp_carrier"]==1) & (combined[outcome_col]==1)).sum())
    b = int(((combined["snp_carrier"]==1) & (combined[outcome_col]==0)).sum())
    c = int(((combined["snp_carrier"]==0) & (combined[outcome_col]==1)).sum())
    d = int(((combined["snp_carrier"]==0) & (combined[outcome_col]==0)).sum())
    if a + b == 0 or c + d == 0:
        return None

    # ── 1. Unadjusted Fisher exact ────────────────────────────────────────
    _, p_unadj = fisher_exact([[a, b], [c, d]])
    or_unadj, or_meth = _haldane_or(a, b, c, d)

    result = {
        "Test_Unadj":         "Fisher exact",
        "N_Carriers":         n_c,
        "N_NonCarriers":      n_nc,
        "Events_Carriers":    a,
        "Events_NonCarriers": c,
        "OR_Unadj":           round(or_unadj, 4),
        "OR_Unadj_Method":    or_meth,
        "P_Unadj":            p_unadj,
    }

    # ── 2. Age-adjusted logistic regression ──────────────────────────────
    events_min   = min(a + c, b + d)
    model_df     = combined[[outcome_col, age_col, "snp_carrier"]].dropna()
    n_with_age   = len(model_df)

    if events_min < MIN_EVENTS_FOR_LOGISTIC or n_with_age < 10 or sm_logit is None:
        reason = (f"Too few events ({events_min} < {MIN_EVENTS_FOR_LOGISTIC})"
                  if events_min < MIN_EVENTS_FOR_LOGISTIC
                  else f"Too few with age data ({n_with_age})")
        result.update({
            "Test_Adj": "Logistic (age-adjusted)", "OR_Adj": np.nan,
            "OR_Adj_CI95": np.nan, "P_Adj": np.nan,
            "N_Adj": n_with_age, "Adj_Note": f"Skipped: {reason}",
        })
        return result

    try:
        model_df = model_df.copy()
        model_df.columns = ["outcome", "age", "snp"]
        model_df["age_z"] = (model_df["age"] - model_df["age"].mean()) / model_df["age"].std()
        fit = sm_logit("outcome ~ snp + age_z", data=model_df).fit(disp=0, maxiter=200)
        or_adj = round(float(np.exp(fit.params["snp"])), 4)
        ci     = np.exp(fit.conf_int().loc["snp"])
        p_adj  = round(float(fit.pvalues["snp"]), 4)
        result.update({
            "Test_Adj":    "Logistic (age-adjusted)",
            "OR_Adj":      or_adj,
            "OR_Adj_CI95": f"[{ci.iloc[0]:.3f}, {ci.iloc[1]:.3f}]",
            "P_Adj":       p_adj,
            "N_Adj":       len(model_df),
            "Adj_Note":    "Converged OK",
        })
    except Exception as e:
        result.update({
            "Test_Adj": "Logistic (age-adjusted)", "OR_Adj": np.nan,
            "OR_Adj_CI95": np.nan, "P_Adj": np.nan,
            "N_Adj": n_with_age, "Adj_Note": f"Failed: {str(e)[:80]}",
        })
    return result


def _test_nominal(c_df, nc_df, col) -> Optional[Dict]:
    """Chi-square for multi-category variables. Unadjusted only."""
    n_c  = c_df[col].notna().sum()
    n_nc = nc_df[col].notna().sum()
    if n_c < MIN_CARRIERS or n_nc < MIN_CARRIERS:
        return None
    ct = pd.crosstab(
        pd.concat([pd.Series(["C"]  * len(c_df),  index=c_df.index),
                   pd.Series(["NC"] * len(nc_df), index=nc_df.index)]),
        pd.concat([c_df[col], nc_df[col]])
    ).dropna(axis=1)
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return None
    chi2, p, dof, _ = stats.chi2_contingency(ct)
    return {
        "Test_Unadj":    "Chi-square",
        "N_Carriers":    int(n_c),
        "N_NonCarriers": int(n_nc),
        "Stat_Unadj":    round(chi2, 3),
        "P_Unadj":       p,
        "Test_Adj":      "N/A",
        "OR_Adj":        np.nan,
        "OR_Adj_CI95":   np.nan,
        "P_Adj":         np.nan,
        "N_Adj":         np.nan,
        "Adj_Note":      "Not applicable for nominal outcome",
    }


def _run_clin_for_cohort(df_t: pd.DataFrame, var_dict: Dict, cohort_label: str) -> pd.DataFrame:
    """
    Run all clinical variable tests for one cancer cohort (tumour samples only).
    Healthy controls (MN/EN) are never passed here — they are excluded upstream
    by the Tissue == 'Tumour' filter in clinical_associations().
    """
    rows = []
    for var_id, vdf in df_t.groupby("Variant_ID"):
        sym = vdf["SYMBOL"].dropna().iloc[0] if vdf["SYMBOL"].notna().any() else ""
        carrier_samples = set(vdf["Sample"].unique())
        sample_data     = df_t.drop_duplicates("Sample").set_index("Sample")
        c_df  = sample_data[sample_data.index.isin(carrier_samples)]
        nc_df = sample_data[~sample_data.index.isin(carrier_samples)]
        if len(c_df) < MIN_CARRIERS or len(nc_df) < MIN_CARRIERS:
            continue

        for col, meta in var_dict.items():
            if col not in sample_data.columns:
                continue
            if meta["type"] == "continuous":
                res = _test_continuous(c_df[col], nc_df[col])
            elif meta["type"] == "binary":
                res = _test_binary(c_df, nc_df, col, age_col="canon__age")
            elif meta["type"] == "nominal":
                res = _test_nominal(c_df, nc_df, col)
            else:
                res = None
            if res is None:
                continue
            row = {
                "Cohort":       cohort_label,
                "Variant_ID":   var_id,
                "Gene":         sym,
                "Clinical_Var": col,
                "Clin_Label":   meta["label"],
                "Clin_Type":    meta["type"],
                "Note":         meta.get("note", ""),
            }
            row.update(res)
            rows.append(row)

    res_df = pd.DataFrame(rows)
    if res_df.empty:
        return res_df

    # FDR correction separately for unadjusted and adjusted p-values,
    # applied within each clinical variable across all variants
    parts = []
    for clin_col in res_df["Clinical_Var"].unique():
        sub = res_df[res_df["Clinical_Var"] == clin_col].copy()
        sub = _apply_fdr(sub, p_col="P_Unadj", out_col="FDR_Unadj")
        if sub["P_Adj"].notna().any():
            sub = _apply_fdr(sub, p_col="P_Adj", out_col="FDR_Adj")
        else:
            sub["FDR_Adj"] = np.nan
        parts.append(sub)

    res_df = pd.concat(parts, ignore_index=True).sort_values(["Clinical_Var", "P_Unadj"])
    res_df["Nominal_Sig_Unadj"] = res_df["P_Unadj"] < 0.05
    res_df["FDR_Sig_Unadj"]     = res_df["FDR_Unadj"] < FDR_THRESHOLD
    res_df["Nominal_Sig_Adj"]   = res_df["P_Adj"].notna() & (res_df["P_Adj"] < 0.05)
    res_df["FDR_Sig_Adj"]       = res_df["FDR_Adj"].notna() & (res_df["FDR_Adj"] < FDR_THRESHOLD)
    return res_df


def clinical_associations(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    print("=== Analysis 2: SNP × Clinical Variable ===")
    df = df[~df["is_replicate"] & (df["Tissue"] == "Tumour")].copy()
    breast_df = df[df["Cohort"].str.contains("Breast",    case=False, na=False)]
    endo_df   = df[df["Cohort"].str.contains("Endometri", case=False, na=False)]
    print(f"  Breast tumour:      {breast_df['Sample'].nunique()} samples")
    print(f"  Endometrial tumour: {endo_df['Sample'].nunique()} samples")
    breast_res = _run_clin_for_cohort(breast_df, CLINICAL_VARS_BREAST, "Breast_Tumour")
    endo_res   = _run_clin_for_cohort(endo_df,   CLINICAL_VARS_ENDO,   "Endometrial_Tumour")
    for label, res in [("Breast", breast_res), ("Endometrial", endo_res)]:
        if res.empty: print(f"  {label}: no results")
        else:
            print(f"  {label}: {len(res)} tests | {res['Nominal_Sig_Unadj'].sum()} nominal (unadj) | {res['FDR_Sig_Unadj'].sum()} FDR (unadj) | {res['Nominal_Sig_Adj'].sum()} nominal (adj)")
    print()
    return breast_res, endo_res


# ── ANALYSIS 3: GENOTYPE DOSE ─────────────────────────────────────────────────

def genotype_dose_analysis(df: pd.DataFrame, top_variants: List[str],
                            var_dict: Dict, cohort_label: str) -> pd.DataFrame:
    print(f"=== Analysis 3: Genotype-dose — {cohort_label} ===")
    df = df[~df["is_replicate"] & (df["Tissue"] == "Tumour")].copy()
    cohort_filter = "Breast" if "Breast" in cohort_label else "Endometri"
    df = df[df["Cohort"].str.contains(cohort_filter, case=False, na=False)]
    if "GT" not in df.columns:
        print("  No GT column — skipping.\n"); return pd.DataFrame()

    rows = []
    for var_id in top_variants:
        vdf = df[df["Variant_ID"] == var_id]
        if vdf.empty: continue
        sym = vdf["SYMBOL"].dropna().iloc[0] if vdf["SYMBOL"].notna().any() else ""
        sd  = vdf.drop_duplicates("Sample").set_index("Sample").copy()
        sd["dose"] = sd["GT"].map({"0/0": 0, "0/1": 1, "1/1": 2})

        for col, meta in var_dict.items():
            if col not in sd.columns: continue
            groups = [sd[sd["dose"] == g][col].dropna() for g in [0, 1, 2]]
            groups = [g for g in groups if len(g) >= 2]
            if len(groups) < 2: continue
            if meta["type"] == "continuous":
                stat, p = stats.kruskal(*[pd.to_numeric(g, errors="coerce").dropna().values for g in groups])
                test = "Kruskal-Wallis"
            else:
                combined = pd.concat([sd[sd["dose"] == g][[col]].assign(dose_group=g)
                                       for g in [0,1,2] if len(sd[sd["dose"]==g]) >= 2])
                ct = pd.crosstab(combined["dose_group"], combined[col])
                if ct.shape[0] < 2 or ct.shape[1] < 2: continue
                stat, p, _, _ = stats.chi2_contingency(ct)
                test = "Chi-square (trend)"
            rows.append({
                "Cohort": cohort_label, "Variant_ID": var_id, "Gene": sym,
                "Clinical_Var": col, "Clin_Label": meta["label"], "Test": test,
                "N_0/0": len(groups[0]) if len(groups) > 0 else 0,
                "N_0/1": len(groups[1]) if len(groups) > 1 else 0,
                "N_1/1": len(groups[2]) if len(groups) > 2 else 0,
                "Stat": round(stat, 3), "P_Value": p,
            })
    res = pd.DataFrame(rows)
    if not res.empty:
        res = _apply_fdr(res)
        res["Nominal_Sig"] = res["P_Value"] < 0.05
        res["FDR_Sig"]     = res["FDR_P_Value"] < FDR_THRESHOLD
        print(f"  {len(res)} tests | {res['Nominal_Sig'].sum()} nominal | {res['FDR_Sig'].sum()} FDR\n")
    else:
        print("  No results.\n")
    return res


# ── ANALYSIS 4: SURVIVAL (AU ENDO) ───────────────────────────────────────────

def survival_analysis(df: pd.DataFrame) -> pd.DataFrame:
    if not _HAS_LIFELINES:
        print("=== Analysis 4: Survival — SKIPPED (lifelines not installed) ===\n")
        return pd.DataFrame()

    print("=== Analysis 4: Survival — AU Endometrial ===")
    df = df[~df["is_replicate"] & (df["Tissue"] == "Tumour")
            & df["Cohort"].str.contains("Endometri", case=False, na=False)].copy()

    rows = []
    km_pages = []

    for endpoint, ev_col in [("OS", "ENDO_EXITUS_BIN"), ("PFS", "ENDO_PD_BIN")]:
        t_col = "canon__os_months" if endpoint == "OS" else "canon__pfs_months"
        for var_id, vdf in df.groupby("Variant_ID"):
            sym = vdf["SYMBOL"].dropna().iloc[0] if vdf["SYMBOL"].notna().any() else ""
            carrier_set = set(vdf["Sample"].unique())
            sd = df.drop_duplicates("Sample").set_index("Sample")
            if t_col not in sd.columns or ev_col not in sd.columns: continue
            valid = sd[[t_col, ev_col]].dropna()
            if len(valid) < 6: continue
            T = valid[t_col]
            E = valid[ev_col].astype(float)
            x = pd.Series(sd.index.isin(carrier_set).astype(int), index=sd.index)
            x = x.loc[valid.index]
            if x.nunique() < 2: continue

            # Log-rank
            try:
                lr = logrank_test(T[x==1], T[x==0], E[x==1], E[x==0])
                lr_p = round(lr.p_value, 4)
            except Exception: lr_p = np.nan

            # Cox PH
            try:
                cox_df = pd.DataFrame({"T": T, "E": E, "snp": x})
                cph = CoxPHFitter()
                cph.fit(cox_df, duration_col="T", event_col="E", show_progress=False)
                hr  = round(np.exp(cph.params_["snp"]), 3)
                ci  = np.exp(cph.confidence_intervals_.loc["snp"])
                p_cox = round(cph.summary.loc["snp", "p"], 4)
                hr_ci = f"[{ci.iloc[0]:.3f}, {ci.iloc[1]:.3f}]"
            except Exception:
                hr = np.nan; hr_ci = np.nan; p_cox = np.nan

            rows.append({
                "Endpoint": endpoint, "Variant_ID": var_id, "Gene": sym,
                "N_total": len(valid), "N_carriers": int(x.sum()),
                "N_events": int(E.sum()),
                "Logrank_P": lr_p, "HR": hr, "HR_CI95": hr_ci, "P_Cox": p_cox,
            })

            # KM plot
            fig, ax = plt.subplots(figsize=(6, 4))
            kmf = KaplanMeierFitter()
            for g, label, col in [(1, "Carrier", "#D32F2F"), (0, "Non-carrier", "#1976D2")]:
                m = x == g
                if m.sum() < 2: continue
                kmf.fit(T[m], E[m], label=f"{label} (n={m.sum()})")
                kmf.plot_survival_function(ax=ax, ci_show=True, color=col)
            ax.set_title(f"{var_id} — {sym} | {endpoint} | p={lr_p}", fontsize=9)
            ax.set_xlabel(f"{endpoint} (months)")
            ax.set_ylabel("Survival probability")
            ax.legend(fontsize=8)
            plt.tight_layout()
            km_pages.append(fig)

    res = pd.DataFrame(rows)
    if not res.empty:
        for ep in res["Endpoint"].unique():
            sub = res[res["Endpoint"] == ep].copy()
            res.loc[sub.index, "FDR_P_Cox"] = _apply_fdr(sub, "P_Cox", "FDR_P_Cox")["FDR_P_Cox"]
        print(f"  {len(res)} tests across {res['Endpoint'].nunique()} endpoints\n")
    else:
        print("  No results (likely insufficient OS/PFS data).\n")

    return res, km_pages


# ── VISUALISATION ─────────────────────────────────────────────────────────────

def _draw_volcano(tvh, p_col, sig_col, thresh, title_suffix, out_path):
    groups = tvh["Analysis_Group"].unique()
    fig, axes = plt.subplots(1, len(groups), figsize=(7*len(groups), 6), squeeze=False)
    axes = axes[0]
    impact_colours = {"HIGH": "#d32f2f", "MODERATE": "#f57c00",
                      "LOW": "#388e3c", "MODIFIER": "#90a4ae"}
    for ax, grp in zip(axes, groups):
        sub = tvh[tvh["Analysis_Group"] == grp].copy()
        sub["-log10p"] = -np.log10(pd.to_numeric(sub[p_col], errors="coerce").clip(lower=1e-10))
        sub["log2OR"]  = np.log2(sub["Odds_Ratio"].replace({0: 0.001, np.inf: 1000}).clip(0.001, 1000))
        for imp, col in impact_colours.items():
            m = sub["IMPACT"] == imp
            ax.scatter(sub.loc[m, "log2OR"], sub.loc[m, "-log10p"],
                       c=col, label=imp, alpha=0.8, s=70, edgecolors="k", linewidths=0.4)
        ax.axhline(-np.log10(thresh), color="red", linestyle="--", linewidth=1.2, alpha=0.8)
        ax.axvline(0, color="grey", linestyle="--", linewidth=1, alpha=0.5)
        for _, row in sub[sub[sig_col]].nsmallest(8, p_col).iterrows():
            ax.annotate(row["Variant_ID"], xy=(row["log2OR"], row["-log10p"]),
                        xytext=(8, 4), textcoords="offset points", fontsize=7, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.2", fc="yellow", alpha=0.75),
                        arrowprops=dict(arrowstyle="->", lw=0.8))
        ax.set_title(grp, fontsize=13, fontweight="bold")
        ax.set_xlabel("log₂(OR)"); ax.set_ylabel(f"-log₁₀(p)")
        ax.grid(True, linestyle=":", alpha=0.35)
    handles = [mpatches.Patch(color=c, label=i) for i, c in impact_colours.items()]
    fig.legend(handles=handles, title="VEP Impact", bbox_to_anchor=(1.01, 0.5), loc="center left")
    fig.suptitle(f"Tumour vs Healthy — {title_suffix}", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def make_volcano_plots(tvh, out_dir):
    if tvh.empty: return
    _draw_volcano(tvh, "P_Value",     "Nominal_Sig", 0.05,          "Raw p-value",     out_dir/"17_SNP_Volcano_raw_p.png")
    _draw_volcano(tvh, "FDR_P_Value", "FDR_Sig",     FDR_THRESHOLD, "FDR-corrected",   out_dir/"17_SNP_Volcano_FDR.png")


def _draw_heatmap(clin_res, p_col, title_suffix, out_path, filter_thresh=0.20):
    pivot = clin_res.pivot_table(index="Variant_ID", columns="Clin_Label",
                                  values=p_col, aggfunc="min").dropna(how="all")
    pivot = pivot[(pivot < filter_thresh).any(axis=1)]
    if pivot.empty:
        print(f"  Heatmap ({title_suffix}): no hits below {filter_thresh} — skipped")
        return
    log_piv = -np.log10(pivot.fillna(1).clip(lower=1e-10))
    fig, ax = plt.subplots(figsize=(max(7, len(pivot.columns)*1.3), max(5, len(pivot)*0.45)))
    sns.heatmap(log_piv, ax=ax, cmap="YlOrRd", linewidths=0.5, linecolor="lightgrey",
                annot=True, fmt=".2f", cbar_kws={"label": f"-log₁₀(p)"})
    sig_t = -np.log10(0.05 if "Raw" in title_suffix else FDR_THRESHOLD)
    for i, var in enumerate(log_piv.index):
        for j, clin in enumerate(log_piv.columns):
            if log_piv.loc[var, clin] >= sig_t:
                ax.text(j+0.5, i+0.15, "★", ha="center", va="top",
                        fontsize=9, color="white", fontweight="bold")
    ax.set_title(f"SNP × Clinical — {title_suffix}", fontsize=12, fontweight="bold")
    ax.set_xlabel("Clinical Variable"); ax.set_ylabel("Variant")
    plt.xticks(rotation=40, ha="right", fontsize=9); plt.yticks(fontsize=8)
    plt.tight_layout(); plt.savefig(out_path, dpi=300, bbox_inches="tight"); plt.close()
    print(f"  Saved: {out_path}")


def make_heatmaps(breast_clin, endo_clin, out_dir):
    for clin_res, label in [(breast_clin, "Breast"), (endo_clin, "Endometrial")]:
        if clin_res.empty: continue
        _draw_heatmap(clin_res, "P_Unadj",   f"Raw p-value (unadjusted) — {label}",
                      out_dir/f"17_SNP_Heatmap_{label}_raw_p.png",  filter_thresh=0.10)
        _draw_heatmap(clin_res, "FDR_Unadj", f"FDR-corrected (unadjusted) — {label}",
                      out_dir/f"17_SNP_Heatmap_{label}_FDR.png",    filter_thresh=0.30)
        # Age-adjusted heatmap (only rows where adjusted test was run)
        adj_res = clin_res[clin_res["P_Adj"].notna()].copy()
        if not adj_res.empty:
            _draw_heatmap(adj_res, "P_Adj",   f"Raw p-value (age-adjusted) — {label}",
                          out_dir/f"17_SNP_Heatmap_{label}_adj_raw_p.png", filter_thresh=0.10)
            _draw_heatmap(adj_res, "FDR_Adj", f"FDR-corrected (age-adjusted) — {label}",
                          out_dir/f"17_SNP_Heatmap_{label}_adj_FDR.png",   filter_thresh=0.30)


def save_km_pdf(km_pages, out_path):
    if not km_pages: return
    with pdf_backend.PdfPages(out_path) as pdf:
        for fig in km_pages:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"  Saved: {out_path} ({len(km_pages)} KM curves)")


# ── MAIN ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="SNP–clinical association analysis (script 17 v2)")
    p.add_argument("--gsdmb",   default=str(DEFAULT_GSDMB))
    p.add_argument("--master",  default=str(DEFAULT_MASTER))
    p.add_argument("--out_dir", default=str(DEFAULT_OUT))
    p.add_argument("--manifest_endo_tumour", default=str(DEFAULT_MANIFESTS["endometrium-tumour"]),
                   help="Pass-BAM manifest for AU endometrial tumour")
    p.add_argument("--manifest_endo_normal", default=str(DEFAULT_MANIFESTS["endometrium-normal"]),
                   help="Pass-BAM manifest for SP endometrial normal")
    p.add_argument("--manifest_breast_tumour", default=str(DEFAULT_MANIFESTS["breast-tumour"]),
                   help="Pass-BAM manifest for SP breast tumour")
    p.add_argument("--manifest_breast_normal", default=str(DEFAULT_MANIFESTS["breast-normal"]),
                   help="Pass-BAM manifest for SP breast normal")
    p.add_argument("--no_manifests", action="store_true",
                   help="Skip manifest filter (use extraction_flag fallback instead)")
    return p.parse_args()


def main():
    args    = parse_args()
    gsdmb   = Path(args.gsdmb)
    master  = Path(args.master)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_xlsx = out_dir / "17_SNP_Clinical_Association_Results.xlsx"

    # Build manifest dict (None if --no_manifests flag set)
    if args.no_manifests:
        manifest_paths = None
    else:
        manifest_paths = {
            "endometrium-tumour": Path(args.manifest_endo_tumour),
            "endometrium-normal": Path(args.manifest_endo_normal),
            "breast-tumour":      Path(args.manifest_breast_tumour),
            "breast-normal":      Path(args.manifest_breast_normal),
        }

    merged = load_and_merge(gsdmb, master, manifest_paths=manifest_paths)

    # Analysis 1
    tvh = tumour_vs_healthy(merged)

    # Analysis 2
    breast_clin, endo_clin = clinical_associations(merged)

    # Analysis 3
    def _top_vars(res):
        return [] if res.empty else res[res["Nominal_Sig_Unadj"]]["Variant_ID"].unique().tolist()
    breast_dose = genotype_dose_analysis(merged, _top_vars(breast_clin), CLINICAL_VARS_BREAST, "Breast_Tumour")
    endo_dose   = genotype_dose_analysis(merged, _top_vars(endo_clin),   CLINICAL_VARS_ENDO,   "Endometrial_Tumour")

    # Analysis 4 — survival
    surv_res = pd.DataFrame()
    km_pages = []
    surv_out = survival_analysis(merged)
    if isinstance(surv_out, tuple):
        surv_res, km_pages = surv_out

    # Summary
    sig_parts = []
    if not tvh.empty:
        s = tvh[tvh["Nominal_Sig"]].copy(); s["Analysis_Type"] = "Tumour_vs_Healthy"
        sig_parts.append(s[["Analysis_Type", "Analysis_Group", "Variant_ID", "Gene",
                             "Consequence", "IMPACT", "Freq_Tumour_%", "Freq_Healthy_%",
                             "Odds_Ratio", "P_Value", "FDR_P_Value", "Nominal_Sig", "FDR_Sig"]])
    for res, label in [(breast_clin, "Breast_Tumour"), (endo_clin, "Endometrial_Tumour")]:
        if not res.empty:
            # Include rows nominal in either unadjusted OR adjusted test
            sig_mask = res["Nominal_Sig_Unadj"] | res["Nominal_Sig_Adj"]
            s = res[sig_mask].copy(); s["Analysis_Type"] = f"Clinical_{label}"
            keep_cols = ["Analysis_Type", "Cohort", "Variant_ID", "Gene", "Clin_Label",
                         "Clin_Type", "N_Carriers", "N_NonCarriers",
                         "Test_Unadj", "OR_Unadj", "P_Unadj", "FDR_Unadj",
                         "Nominal_Sig_Unadj", "FDR_Sig_Unadj",
                         "Test_Adj", "OR_Adj", "OR_Adj_CI95", "P_Adj", "FDR_Adj",
                         "Nominal_Sig_Adj", "FDR_Sig_Adj", "Adj_Note"]
            s = s[[c for c in keep_cols if c in s.columns]]
            sig_parts.append(s)
    summary = (pd.concat(sig_parts, ignore_index=True).sort_values("P_Unadj")
               if sig_parts else pd.DataFrame({"Note": ["No nominally significant results."]}))

    # Sample manifest
    keep_clin = ["canon__age", "canon__grade", "clin_au_endo__FIGO_STAGE",
                 "canon__molecular_class", "canon__msi_status", "canon__er_status",
                 "canon__pr_status", "canon__her2_copies", "canon__os_months",
                 "canon__pfs_months", "ENDO_EXITUS_BIN", "BREAST_EXITUS_DERIVED",
                 "BREAST_RECURRENCE_DERIVED", "ENDO_PD_BIN"]
    manifest = merged.drop_duplicates("Sample")[[
        c for c in ["Sample", "snp_code", "sample_id", "Cohort", "Tissue", "sheet",
                    "tumour_normal", "is_replicate"] + keep_clin
        if c in merged.columns
    ]].sort_values(["Cohort", "Tissue", "Sample"]).reset_index(drop=True)

    # Write Excel
    print("=== Writing output Excel ===")
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        if not tvh.empty:       tvh.to_excel(xw,         sheet_name="tumour_vs_healthy",     index=False)
        if not breast_clin.empty: breast_clin.to_excel(xw,sheet_name="breast_clinical_assoc", index=False)
        if not endo_clin.empty:   endo_clin.to_excel(xw,  sheet_name="endo_clinical_assoc",   index=False)
        dose_all = pd.concat([breast_dose, endo_dose], ignore_index=True)
        if not dose_all.empty:  dose_all.to_excel(xw,    sheet_name="genotype_dose",          index=False)
        if not surv_res.empty:  surv_res.to_excel(xw,    sheet_name="survival_au_endo",       index=False)
        summary.to_excel(xw,                             sheet_name="summary_significant",    index=False)
        manifest.to_excel(xw,                            sheet_name="sample_manifest",        index=False)
    print(f"  Saved: {out_xlsx}\n")

    # Plots
    print("=== Generating plots ===")
    make_volcano_plots(tvh, out_dir)
    make_heatmaps(breast_clin, endo_clin, out_dir)
    save_km_pdf(km_pages, out_dir / "17_KM_Curves_AU_Endo.pdf")

    # Final summary
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"  Samples in manifest (sequenced): {merged[~merged['is_replicate']]['Sample'].nunique()}")
    print(f"  Variants tested  : {merged['Variant_ID'].nunique()}")
    if not tvh.empty:
        print(f"  Tumour vs healthy   — nominal: {tvh['Nominal_Sig'].sum()} | FDR<{FDR_THRESHOLD}: {tvh['FDR_Sig'].sum()}")
    if not breast_clin.empty:
        print(f"  Breast clinical     — nominal unadj: {breast_clin['Nominal_Sig_Unadj'].sum()} | FDR: {breast_clin['FDR_Sig_Unadj'].sum()} | nominal adj: {breast_clin['Nominal_Sig_Adj'].sum()}")
    if not endo_clin.empty:
        print(f"  Endometrial clinical— nominal unadj: {endo_clin['Nominal_Sig_Unadj'].sum()} | FDR: {endo_clin['FDR_Sig_Unadj'].sum()} | nominal adj: {endo_clin['Nominal_Sig_Adj'].sum()}")
    if not surv_res.empty:
        print(f"  Survival (AU endo)  — {len(surv_res)} tests")
    print(f"\n  Results : {out_xlsx}")
    print(f"  Plots   : {out_dir}")


if __name__ == "__main__":
    main()