#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
17_snp-association.py  (v4 — SNP-only, genotypic model, all-cohort Cox, age+BMI adjusted)
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

  MT-T_N (breast tumour + paired normal):
    No exclusion flags present → all 71 DNA samples retained

  EN (endometrial healthy):
    No exclusion flags present → all 80 DNA samples retained

  MN (breast healthy):
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
                              Age- and BMI-adjusted logistic regression for binary outcomes
3. Genotype-dose (WT/Het/Hom) — Overall trend + separate Het-vs-WT and
                              Hom-vs-WT contrasts for top SNPs, both cohorts
4. Survival (Cox PH)       — Age- and BMI-adjusted Cox proportional hazards for OS
                              and PFS in all cohorts with survival data:
                              • Endometrial: OS + PFS (canon__os/pfs_months)
                              • Breast: OS (BREAST_OS_MONTHS_DERIVED)
                              Genotypic model: WT reference, Het and Hom as
                              separate terms.  KM curves per genotype class.

VARIABLE COVERAGE (verified from data)
---------------------------------------
BREAST (MT-T_N sheet, n=71 tumour samples)
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

ENDOMETRIAL (AT=AUs sheet, n=114 sequenced DNA samples)
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
    • survival_cox          (all cohorts with survival data)
    • summary_significant
    • sample_manifest

  17_SNP_Volcano_Plots_raw_p.png
  17_SNP_Volcano_Plots_FDR.png
  17_SNP_Heatmap_Breast_raw_p.png
  17_SNP_Heatmap_Breast_FDR.png
  17_SNP_Heatmap_Endometrial_raw_p.png
  17_SNP_Heatmap_Endometrial_FDR.png
  17_KM_Curves_Endometrial.pdf

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

def _safe_exp(x, clip: float = 50.0):
    """Exponentiation that avoids overflow in small-N / separation fits.

    Works with scalars, numpy arrays, pandas Series/DataFrames.
    """
    import numpy as _np
    try:
        if hasattr(x, "values") and hasattr(x, "index") and not hasattr(x, "columns"):
            # pandas Series / Index-like
            vals = _np.asarray(x.values, dtype=float)
            out = _np.exp(_np.clip(vals, -clip, clip))
            try:
                import pandas as _pd
                return _pd.Series(out, index=x.index)
            except Exception:
                return out
        if hasattr(x, "values") and hasattr(x, "index") and hasattr(x, "columns"):
            # pandas DataFrame
            vals = _np.asarray(x.values, dtype=float)
            out = _np.exp(_np.clip(vals, -clip, clip))
            try:
                import pandas as _pd
                return _pd.DataFrame(out, index=x.index, columns=x.columns)
            except Exception:
                return out
        # scalar / numpy array / list
        vals = _np.asarray(x, dtype=float)
        out = _np.exp(_np.clip(vals, -clip, clip))
        # return python float for scalar inputs
        if _np.ndim(out) == 0:
            return float(out)
        return out
    except Exception:
        return _np.nan

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
def _choose_bmi_col(sample_df: pd.DataFrame, cohort_label: str) -> Optional[str]:
    """Pick the most appropriate BMI column for a cohort.

    Breast cohort in the harmonised master often stores BMI in BREAST_BMI_NUMERIC.
    Other cohorts may store BMI in canon__bmi. Returns None if no usable BMI.
    """
    # Order matters: prefer cohort-specific numeric BMI if present
    if "Breast" in cohort_label:
        candidates = ["BREAST_BMI_NUMERIC", "canon__bmi"]
    else:
        candidates = ["canon__bmi"]

    for c in candidates:
        if c in sample_df.columns:
            n_nonnull = int(pd.to_numeric(sample_df[c], errors="coerce").notna().sum())
            if n_nonnull >= MIN_CARRIERS:
                return c
    return None
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
        for sheet, label in [("AT=AUs", "endo-tumour"), ("MT-T_N", "breast-tumour"),
                               ("EN", "endo-normal"), ("MN", "breast-normal")]:
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
    # Breast tumour
    master_dna["BREAST_RECURRENCE_DERIVED"] = master_dna["clin_dcs__Recaida/Progresión"].apply(_derive_breast_recurrence)
    master_dna["BREAST_EXITUS_DERIVED"]     = master_dna["clin_dcs__Exitus"].apply(_derive_breast_exitus)
    master_dna["BREAST_METASTASIS_DERIVED"] = master_dna["clin_dcs__MTxDISTANCIA"].apply(_derive_breast_metastasis)
    master_dna["BREAST_HER2_SUBTYPE"]       = master_dna["clin_her2__DX"].apply(_derive_her2_subtype)
    master_dna["BREAST_ER_BIN"]             = master_dna["canon__er_status"].map({"Positive": 1, "Negative": 0})
    master_dna["BREAST_PR_BIN"]             = master_dna["canon__pr_status"].map({"Positive": 1, "Negative": 0})
    master_dna["BREAST_KI67_NUMERIC"]       = master_dna["clin_dcs__KI67"].apply(_derive_ki67_numeric)
    master_dna["BREAST_OS_MONTHS_DERIVED"]  = master_dna.apply(_derive_breast_os, axis=1)
    master_dna["BREAST_GRADE_NUMERIC"]      = pd.to_numeric(master_dna["clin_dcs__GRADO"], errors="coerce")
    master_dna["BREAST_P53_NUMERIC"]        = pd.to_numeric(master_dna["clin_dcs__p53"], errors="coerce")
    master_dna["BREAST_BMI_NUMERIC"]        = pd.to_numeric(master_dna["clin_dcs__BMI"], errors="coerce")
    master_dna["BREAST_MENARCHE_NUMERIC"]   = pd.to_numeric(master_dna["canon__menarche_age"], errors="coerce")
    master_dna["BREAST_MENOPAUSE_NUMERIC"]  = pd.to_numeric(master_dna["canon__menopause_age"], errors="coerce")
    # Breast DX type (CDI vs CDIS vs other)
    master_dna["BREAST_DX_TYPE"] = master_dna["clin_dcs__Dx"].apply(
        lambda x: "CDI" if pd.notna(x) and "CDI" in str(x).upper() and "CDIS" not in str(x).upper()
        else ("CDIS" if pd.notna(x) and "CDIS" in str(x).upper() else (str(x).strip() if pd.notna(x) else None))
    )
    # Local metastasis binary
    master_dna["BREAST_LOCAL_MET_BIN"] = master_dna["clin_dcs__MTxDISTANCIA"].map(
        {"NO": 0, "NO-LOCAL": 1, "SI": 0}   # local only vs not
    )

    # Endometrial AU tumour
    master_dna["ENDO_FIGO_NUMERIC"]         = master_dna["clin_au_endo__FIGO_STAGE"].apply(_derive_figo_numeric)
    master_dna["ENDO_GRADE_NUMERIC"]        = master_dna["clin_au_endo__GRADE"].apply(_derive_endo_grade)
    master_dna["ENDO_LVSI_BIN"]             = master_dna["clin_au_endo__LVSI"].map({"YES": 1, "NO": 0})
    master_dna["ENDO_MYOINV_BIN"]           = master_dna["clin_au_endo__MYOMETRIAL_INFILTRATION"].map({"<50%": 0, ">50%": 1})
    master_dna["ENDO_MSI_BIN"]              = master_dna["canon__msi_status"].map({"Unstable": 1, "Stable": 0})
    master_dna["ENDO_NEEC_BIN"]             = master_dna["clin_au_endo__HISTOLOGY_GROUP"].map({"NEEC": 1, "EEC": 0})
    master_dna["ENDO_PD_BIN"]               = master_dna["clin_au_endo__PD_STATUS"].map({"PD": 1, "NO PD": 0})
    master_dna["ENDO_EXITUS_BIN"]           = master_dna["clin_au_endo__EXITUS"].map({"YES": 1, "NO": 0})
    master_dna["ENDO_EXITUS_DISEASE_BIN"]   = master_dna["clin_au_endo__EXITUS_DISEASE"].map({"YES": 1, "NO": 0})
    master_dna["ENDO_RISK_ORDINAL"]         = master_dna["clin_au_endo__RISK_OF_RECURRENCE"].apply(_derive_risk_ordinal)
    master_dna["ENDO_ER_BIN"]               = master_dna["canon__er_status"].map({"Positive": 1, "Negative": 0})
    master_dna["ENDO_PR_BIN"]               = master_dna["canon__pr_status"].map({"Positive": 1, "Negative": 0})
    master_dna["ENDO_TP53_ABN_BIN"]         = master_dna["canon__tp53_ihc"].apply(
        lambda x: 0 if pd.notna(x) and str(x).strip().upper() == "WT"
        else (1 if pd.notna(x) else None)   # HIGH / POSITIVE / LOST = abnormal
    )
    master_dna["ENDO_GENE_AMP_BIN"]         = master_dna["clin_au_endo__Gene amplification"].map({"YES": 1, "NO": 0})
    master_dna["ENDO_ITH_BIN"]              = master_dna["clin_au_endo__ITH: intratumor heterogeneity"].map({"YES": 1, "NO": 0})
    master_dna["ENDO_CTDNA_BIN"]            = master_dna["clin_au_endo__BLOOD_BASAL_CTDNA"].map({"POSITIVE": 1, "NEGATIVE": 0})
    master_dna["ENDO_CTDNA_MAF_NUMERIC"]    = pd.to_numeric(master_dna["clin_au_endo__BLOOD_BASAL_CTDNA_MAF"], errors="coerce")
    master_dna["ENDO_CFDN_CONC_NUMERIC"]    = pd.to_numeric(master_dna["clin_au_endo__BLOOD_BASAL_CFDNA_CONCENTRATION"], errors="coerce")
    master_dna["ENDO_PDL1_NUMERIC"]         = pd.to_numeric(master_dna["clin_au_endo__FFPE_PDL1 POLAND RESULTS"], errors="coerce")
    master_dna["ENDO_CD8_NUMERIC"]          = pd.to_numeric(master_dna["clin_au_endo__FFPE_CD8"], errors="coerce")
    master_dna["ENDO_KI67_NUMERIC"]         = master_dna["clin_au_endo__FFPE_KI67"].apply(_derive_ki67_numeric)
    master_dna["ENDO_PTEN_BIN"]             = master_dna["clin_au_endo__FFPE_PTEN"].map({"CONSERVED": 0, "LOST/REDUCED": 1})
    master_dna["ENDO_MLH1_BIN"]             = master_dna["clin_au_endo__FFPE_MLH1"].map({"CONSERVED": 0, "LOST/REDUCED": 1})
    master_dna["ENDO_N_STAGE_BIN"]          = master_dna["clin_au_endo__N"].apply(
        lambda x: 0 if pd.notna(x) and str(x).strip().upper() == "N0" else (1 if pd.notna(x) else None)
    )
    master_dna["ENDO_M_STAGE_BIN"]          = master_dna["clin_au_endo__M"].apply(
        lambda x: 0 if pd.notna(x) and str(x).strip().upper() == "M0" else (1 if pd.notna(x) else None)
    )

    # Keep only relevant columns
    canon_cols = [c for c in master_dna.columns
                  if c.startswith("canon__") and not c.endswith("__source")]
    derived_cols = [c for c in master_dna.columns if c.endswith("_DERIVED") or c.endswith("_BIN")
                    or c.endswith("_NUMERIC") or c.endswith("_ORDINAL") or c.endswith("_SUBTYPE")]
    raw_breast_cols = [
        "clin_dcs__GRADO", "clin_dcs__RE", "clin_dcs__RP",
        "clin_her2__DX", "clin_dcs__Exitus", "clin_dcs__Dx",
        "clin_dcs__Recaida/Progresión", "clin_dcs__MTxDISTANCIA",
        "clin_dcs__KI67", "clin_dcs__Fecha dx", "clin_dcs__p53",
        "clin_dcs__BMI", "clin_dcs__Última fecha disponible",
    ]
    raw_endo_cols = [
        "clin_au_endo__FIGO_STAGE", "clin_au_endo__GRADE",
        "clin_au_endo__HISTOLOGY_GROUP", "clin_au_endo__HISTOLOGY",
        "clin_au_endo__LVSI", "clin_au_endo__MYOMETRIAL_INFILTRATION",
        "clin_au_endo__MSI_STATUS_IHC", "clin_au_endo__RISK_OF_RECURRENCE",
        "clin_au_endo__PD_STATUS", "clin_au_endo__EXITUS",
        "clin_au_endo__EXITUS_DISEASE",
        "clin_au_endo__MOLECULAR CLASSIFICATION_according to IHC and/or NGS profile",
        "clin_au_endo__T", "clin_au_endo__N", "clin_au_endo__M",
        "clin_au_endo__FFPE_TP53_IHC", "clin_au_endo__FFPE_KI67",
        "clin_au_endo__FFPE_PTEN", "clin_au_endo__FFPE_MLH1",
        "clin_au_endo__FFPE_MSH2", "clin_au_endo__FFPE_MSH6", "clin_au_endo__FFPE_PMS2",
        "clin_au_endo__FFPE_PD1", "clin_au_endo__FFPE_PDL1 POLAND RESULTS",
        "clin_au_endo__FFPE_CD8",
        "clin_au_endo__Gene amplification", "clin_au_endo__Amplified genes",
        "clin_au_endo__ITH: intratumor heterogeneity",
        "clin_au_endo__BLOOD_BASAL_CTDNA", "clin_au_endo__BLOOD_BASAL_CTDNA_MAF",
        "clin_au_endo__BLOOD_BASAL_CFDNA_CONCENTRATION",
        "clin_au_endo__BLOOD_BASAL_CFDNA_CLASSIFICATION",
    ]
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
#     Confounders for adjustment: age at diagnosis, BMI
#
#   ENDOMETRIAL AU (AT=AUs) — Australian cohort (Clinical DATA AU_Endometrial)
#     Outcomes: FIGO stage, grade, myometrial invasion, LVSI, MSI, molecular class,
#               histology, ER/PR, progression, exitus, OS/PFS, risk of recurrence
#     Confounders for adjustment: age at surgery, BMI
#
# For BINARY outcomes: both unadjusted (Fisher exact) and age+BMI-adjusted
#   (logistic regression: outcome ~ snp_carrier + age_z + bmi_z) are run.
#   The adjusted model is skipped automatically if fewer than MIN_EVENTS_FOR_LOGISTIC
#   events exist in the smaller group.
# For CONTINUOUS outcomes: Mann-Whitney U (unadjusted only).
# For NOMINAL outcomes: Chi-square (unadjusted only).

MIN_EVENTS_FOR_LOGISTIC = 5   # minimum events in smaller group for adjusted model

CLINICAL_VARS_BREAST: Dict[str, Dict] = {
    # ── Continuous (Mann-Whitney U) ───────────────────────────────────────
    "canon__age":               {"type": "continuous", "label": "Age at diagnosis"},
    "BREAST_GRADE_NUMERIC":     {"type": "continuous", "label": "Tumour grade",
                                  "note": "Ordinal 1/2/3"},
    "BREAST_KI67_NUMERIC":      {"type": "continuous", "label": "KI67",
                                  "note": "Proliferation index (fraction)"},
    "canon__her2_copies":       {"type": "continuous", "label": "HER2 copies (FISH)",
                                  "note": "Continuous copy number"},
    "BREAST_OS_MONTHS_DERIVED": {"type": "continuous", "label": "Overall survival (months)",
                                  "note": "Derived from fecha_dx + ultima_fecha"},
    "BREAST_P53_NUMERIC":       {"type": "continuous", "label": "p53 expression",
                                  "note": "Numeric IHC value"},
    "BREAST_BMI_NUMERIC":       {"type": "continuous", "label": "BMI",
                                  "note": "Body mass index"},
    "BREAST_MENARCHE_NUMERIC":  {"type": "continuous", "label": "Age at menarche"},
    "BREAST_MENOPAUSE_NUMERIC": {"type": "continuous", "label": "Age at menopause"},
    # ── Binary (Fisher exact + age-adjusted logistic) ─────────────────────
    "BREAST_ER_BIN":            {"type": "binary", "label": "ER positive",
                                  "note": "Positive=1 vs Negative=0"},
    "BREAST_PR_BIN":            {"type": "binary", "label": "PR positive",
                                  "note": "Positive=1 vs Negative=0"},
    "BREAST_RECURRENCE_DERIVED":{"type": "binary", "label": "Recurrence / progression",
                                  "note": "Any recurrence=1 vs NO=0"},
    "BREAST_METASTASIS_DERIVED":{"type": "binary", "label": "Distant metastasis",
                                  "note": "SI=1 vs NO/NO-LOCAL=0"},
    "BREAST_LOCAL_MET_BIN":     {"type": "binary", "label": "Local metastasis",
                                  "note": "NO-LOCAL=1 vs NO/SI=0"},
    "BREAST_EXITUS_DERIVED":    {"type": "binary", "label": "Exitus",
                                  "note": "SI=1 vs NO=0"},
    # ── Nominal (chi-square) ─────────────────────────────────────────────
    "BREAST_HER2_SUBTYPE":      {"type": "nominal", "label": "HER2 subtype",
                                  "note": "HER2+ / TN / Other"},
    "BREAST_DX_TYPE":           {"type": "nominal", "label": "Histological diagnosis type",
                                  "note": "CDI / CDIS / other"},
}

CLINICAL_VARS_ENDO: Dict[str, Dict] = {
    # ── Continuous (Mann-Whitney U) ───────────────────────────────────────
    "canon__age":               {"type": "continuous", "label": "Age at surgery"},
    "ENDO_FIGO_NUMERIC":        {"type": "continuous", "label": "FIGO stage (ordinal)",
                                  "note": "IA=1, IB=1.5, II=2, IIIA=3.1 … IVB=4.2"},
    "ENDO_GRADE_NUMERIC":       {"type": "continuous", "label": "Tumour grade",
                                  "note": "G1=1, G2=2, G3=3"},
    "ENDO_RISK_ORDINAL":        {"type": "continuous", "label": "Risk of recurrence (ordinal)",
                                  "note": "LOW=1, INT=2, INT-HIGH=3, HIGH=4"},
    "canon__os_months":         {"type": "continuous", "label": "Overall survival (months)"},
    "canon__pfs_months":        {"type": "continuous", "label": "Progression-free survival (months)"},
    "ENDO_KI67_NUMERIC":        {"type": "continuous", "label": "KI67",
                                  "note": "Proliferation index (fraction)"},
    "ENDO_PDL1_NUMERIC":        {"type": "continuous", "label": "PD-L1 expression (%)",
                                  "note": "FFPE IHC numeric value"},
    "ENDO_CD8_NUMERIC":         {"type": "continuous", "label": "CD8+ TILs (%)",
                                  "note": "FFPE IHC numeric value"},
    "ENDO_CTDNA_MAF_NUMERIC":   {"type": "continuous", "label": "Basal ctDNA MAF",
                                  "note": "Mutant allele fraction (%)"},
    "ENDO_CFDN_CONC_NUMERIC":   {"type": "continuous", "label": "Basal cfDNA concentration",
                                  "note": "ng/mL"},
    # ── Binary (Fisher exact + age-adjusted logistic) ─────────────────────
    "ENDO_LVSI_BIN":            {"type": "binary", "label": "LVSI",
                                  "note": "YES=1 vs NO=0"},
    "ENDO_MYOINV_BIN":          {"type": "binary", "label": "Myometrial invasion ≥50%",
                                  "note": ">50%=1 vs <50%=0"},
    "ENDO_MSI_BIN":             {"type": "binary", "label": "MSI-H",
                                  "note": "Unstable=1 vs Stable=0"},
    "ENDO_NEEC_BIN":            {"type": "binary", "label": "Non-endometrioid histology",
                                  "note": "NEEC=1 vs EEC=0"},
    "ENDO_PD_BIN":              {"type": "binary", "label": "Disease progression",
                                  "note": "PD=1 vs NO PD=0"},
    "ENDO_EXITUS_BIN":          {"type": "binary", "label": "Exitus (all cause)",
                                  "note": "YES=1 vs NO=0"},
    "ENDO_EXITUS_DISEASE_BIN":  {"type": "binary", "label": "Exitus (disease-specific)",
                                  "note": "YES=1 vs NO=0"},
    "ENDO_ER_BIN":              {"type": "binary", "label": "ER positive",
                                  "note": "Positive=1 vs Negative=0"},
    "ENDO_PR_BIN":              {"type": "binary", "label": "PR positive",
                                  "note": "Positive=1 vs Negative=0"},
    "ENDO_TP53_ABN_BIN":        {"type": "binary", "label": "TP53 IHC abnormal",
                                  "note": "HIGH/POSITIVE/LOST=1 vs WT=0"},
    "ENDO_GENE_AMP_BIN":        {"type": "binary", "label": "Gene amplification",
                                  "note": "YES=1 vs NO=0"},
    "ENDO_ITH_BIN":             {"type": "binary", "label": "Intratumoural heterogeneity",
                                  "note": "YES=1 vs NO=0"},
    "ENDO_CTDNA_BIN":           {"type": "binary", "label": "Basal ctDNA detected",
                                  "note": "POSITIVE=1 vs NEGATIVE=0"},
    "ENDO_PTEN_BIN":            {"type": "binary", "label": "PTEN loss/reduced",
                                  "note": "LOST/REDUCED=1 vs CONSERVED=0"},
    "ENDO_MLH1_BIN":            {"type": "binary", "label": "MLH1 loss/reduced",
                                  "note": "LOST/REDUCED=1 vs CONSERVED=0"},
    "ENDO_N_STAGE_BIN":         {"type": "binary", "label": "Lymph node involvement",
                                  "note": "N1/N2=1 vs N0=0"},
    "ENDO_M_STAGE_BIN":         {"type": "binary", "label": "Distant metastasis (M stage)",
                                  "note": "M1=1 vs M0=0"},
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
        "Test_Unadj":              "Mann-Whitney U",
        "N_Carriers":              len(c),
        "N_NonCarriers":           len(nc),
        "Median_Carriers":         round(c.median(), 3),
        "Median_NonCarriers":      round(nc.median(), 3),
        "Stat_Unadj":              round(stat, 3),
        "P_Unadj":                 p,
        # Age-only model
        "Test_Adj_Age":            "N/A",
        "OR_Adj_Age":              np.nan,
        "OR_Adj_Age_CI95":         np.nan,
        "P_Adj_Age":               np.nan,
        "N_Adj_Age":               np.nan,
        "Adj_Age_Note":            "Not applicable for continuous outcome",
        # Age+BMI model
        "Test_Adj_AgeBMI":         "N/A",
        "OR_Adj_AgeBMI":           np.nan,
        "OR_Adj_AgeBMI_CI95":      np.nan,
        "P_Adj_AgeBMI":            np.nan,
        "N_Adj_AgeBMI":            np.nan,
        "Adj_AgeBMI_Note":         "Not applicable for continuous outcome",
    }


def _test_binary(c_df, nc_df, outcome_col, age_col="canon__age", bmi_col: Optional[str] = None) -> Optional[Dict]:
    """
    Unadjusted Fisher exact test + two adjusted logistic regression models:
      Model 1 (age-only):    outcome ~ snp_carrier + age_z          (full N)
      Model 2 (age+BMI):     outcome ~ snp_carrier + age_z + bmi_z  (BMI-complete N)

    This stepwise approach ensures we never lose samples from the age-only model
    just because BMI is missing, while still capturing BMI adjustment where data
    are available. Both models report their own N, OR, CI and p-value.
    """
    try:
        from statsmodels.formula.api import logit as sm_logit
    except ImportError:
        sm_logit = None

    cols = [outcome_col, age_col] + ([bmi_col] if bmi_col else [])
    combined = pd.concat([
        c_df.reindex(columns=cols).assign(snp_carrier=1),
        nc_df.reindex(columns=cols).assign(snp_carrier=0),
    ])
    combined[outcome_col] = pd.to_numeric(combined[outcome_col], errors="coerce")
    combined = combined.dropna(subset=[outcome_col])

    n_c  = int((combined["snp_carrier"] == 1).sum())
    n_nc = int((combined["snp_carrier"] == 0).sum())
    if n_c < MIN_CARRIERS or n_nc < MIN_CARRIERS:
        return None

    # Counts for contingency table (unadjusted)
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

    events_min = min(a + c, b + d)

    def _run_logit(df_model, formula, label):
        """Fit a logistic model and return result dict."""
        n = len(df_model)
        if events_min < MIN_EVENTS_FOR_LOGISTIC or n < 10 or sm_logit is None:
            reason = (f"Too few events ({events_min} < {MIN_EVENTS_FOR_LOGISTIC})"
                      if events_min < MIN_EVENTS_FOR_LOGISTIC
                      else f"Too few samples ({n})")
            return {
                f"Test_Adj_{label}":       f"Logistic ({label})",
                f"OR_Adj_{label}":         np.nan,
                f"OR_Adj_{label}_CI95":    np.nan,
                f"P_Adj_{label}":          np.nan,
                f"N_Adj_{label}":          n,
                f"Adj_{label}_Note":       f"Skipped: {reason}",
            }
        try:
            fit = sm_logit(formula, data=df_model).fit(disp=0, maxiter=200)
            or_v = round(float(_safe_exp(fit.params["snp"])), 4)
            ci   = _safe_exp(fit.conf_int().loc["snp"])
            p_v  = round(float(fit.pvalues["snp"]), 4)
            return {
                f"Test_Adj_{label}":       f"Logistic ({label})",
                f"OR_Adj_{label}":         or_v,
                f"OR_Adj_{label}_CI95":    f"[{ci.iloc[0]:.3f}, {ci.iloc[1]:.3f}]",
                f"P_Adj_{label}":          p_v,
                f"N_Adj_{label}":          n,
                f"Adj_{label}_Note":       "Converged OK",
            }
        except Exception as e:
            return {
                f"Test_Adj_{label}":       f"Logistic ({label})",
                f"OR_Adj_{label}":         np.nan,
                f"OR_Adj_{label}_CI95":    np.nan,
                f"P_Adj_{label}":          np.nan,
                f"N_Adj_{label}":          n,
                f"Adj_{label}_Note":       f"Failed: {str(e)[:80]}",
            }

    # ── 2. Age-only model (uses all samples with age data) ────────────────
    age_df = combined[[outcome_col, age_col, "snp_carrier"]].dropna().copy()
    age_df.columns = ["outcome", "age", "snp"]
    age_std = age_df["age"].std()
    if pd.isna(age_std) or age_std == 0:
        # Cannot standardise age; skip adjusted models
        return result
    age_df["age_z"] = (age_df["age"] - age_df["age"].mean()) / age_std
    result.update(_run_logit(age_df, "outcome ~ snp + age_z", "Age"))

    # ── 3. Age+BMI model (only samples with both age and BMI) ────────────
    if bmi_col and bmi_col in combined.columns:
        bmi_df = combined[[outcome_col, age_col, bmi_col, "snp_carrier"]].dropna().copy()
        bmi_df.columns = ["outcome", "age", "bmi", "snp"]
        age_std2 = bmi_df["age"].std()
        bmi_std  = bmi_df["bmi"].std()
        if not (pd.isna(age_std2) or age_std2 == 0 or pd.isna(bmi_std) or bmi_std == 0):
            bmi_df["age_z"] = (bmi_df["age"] - bmi_df["age"].mean()) / age_std2
            bmi_df["bmi_z"] = (bmi_df["bmi"] - bmi_df["bmi"].mean()) / bmi_std
            result.update(_run_logit(bmi_df, "outcome ~ snp + age_z + bmi_z", "AgeBMI"))

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
        "Test_Unadj":          "Chi-square",
        "N_Carriers":          int(n_c),
        "N_NonCarriers":       int(n_nc),
        "Stat_Unadj":          round(chi2, 3),
        "P_Unadj":             p,
        "Test_Adj_Age":        "N/A",
        "OR_Adj_Age":          np.nan,
        "OR_Adj_Age_CI95":     np.nan,
        "P_Adj_Age":           np.nan,
        "N_Adj_Age":           np.nan,
        "Adj_Age_Note":        "Not applicable for nominal outcome",
        "Test_Adj_AgeBMI":     "N/A",
        "OR_Adj_AgeBMI":       np.nan,
        "OR_Adj_AgeBMI_CI95":  np.nan,
        "P_Adj_AgeBMI":        np.nan,
        "N_Adj_AgeBMI":        np.nan,
        "Adj_AgeBMI_Note":     "Not applicable for nominal outcome",
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
                bmi_col = _choose_bmi_col(sample_data.reset_index(), cohort_label)
                res = _test_binary(c_df, nc_df, col, age_col="canon__age", bmi_col=bmi_col)
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

    # FDR correction separately for unadjusted, age-only, and age+BMI p-values,
    # applied within each clinical variable across all variants
    parts = []
    for clin_col in res_df["Clinical_Var"].unique():
        sub = res_df[res_df["Clinical_Var"] == clin_col].copy()
        sub = _apply_fdr(sub, p_col="P_Unadj", out_col="FDR_Unadj")
        if sub["P_Adj_Age"].notna().any():
            sub = _apply_fdr(sub, p_col="P_Adj_Age", out_col="FDR_Adj_Age")
        else:
            sub["FDR_Adj_Age"] = np.nan
        if sub["P_Adj_AgeBMI"].notna().any():
            sub = _apply_fdr(sub, p_col="P_Adj_AgeBMI", out_col="FDR_Adj_AgeBMI")
        else:
            sub["FDR_Adj_AgeBMI"] = np.nan
        parts.append(sub)

    res_df = pd.concat(parts, ignore_index=True).sort_values(["Clinical_Var", "P_Unadj"])
    res_df["Nominal_Sig_Unadj"]   = res_df["P_Unadj"] < 0.05
    res_df["FDR_Sig_Unadj"]       = res_df["FDR_Unadj"] < FDR_THRESHOLD
    res_df["Nominal_Sig_Adj_Age"] = res_df["P_Adj_Age"].notna() & (res_df["P_Adj_Age"] < 0.05)
    res_df["FDR_Sig_Adj_Age"]     = res_df["FDR_Adj_Age"].notna() & (res_df["FDR_Adj_Age"] < FDR_THRESHOLD)
    res_df["Nominal_Sig_Adj_AgeBMI"] = res_df["P_Adj_AgeBMI"].notna() & (res_df["P_Adj_AgeBMI"] < 0.05)
    res_df["FDR_Sig_Adj_AgeBMI"]     = res_df["FDR_Adj_AgeBMI"].notna() & (res_df["FDR_Adj_AgeBMI"] < FDR_THRESHOLD)
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
            print(f"  {label}: {len(res)} tests | "
                  f"{res['Nominal_Sig_Unadj'].sum()} nominal (unadj) | "
                  f"{res['FDR_Sig_Unadj'].sum()} FDR (unadj) | "
                  f"{res['Nominal_Sig_Adj_Age'].sum()} nominal (age-adj, N={res['N_Adj_Age'].dropna().astype(int).max() if res['N_Adj_Age'].notna().any() else 0}) | "
                  f"{res['Nominal_Sig_Adj_AgeBMI'].sum()} nominal (age+BMI-adj, N={res['N_Adj_AgeBMI'].dropna().astype(int).max() if res['N_Adj_AgeBMI'].notna().any() else 0})")
    print()
    return breast_res, endo_res


# ── ANALYSIS 3: GENOTYPE DOSE ─────────────────────────────────────────────────

def genotype_dose_analysis(df: pd.DataFrame, top_variants: List[str],
                            var_dict: Dict, cohort_label: str) -> pd.DataFrame:
    """
    Analysis 3: genotypic dose test for top SNPs.

    For each variant × clinical variable we run:
      (a) Overall trend test  — Kruskal-Wallis (continuous) or chi-square (categorical)
                                 across all three genotype groups (WT / Het / Hom)
      (b) Het vs WT contrast  — Mann-Whitney U or Fisher exact, Hom samples excluded
      (c) Hom vs WT contrast  — Mann-Whitney U or Fisher exact, Het samples excluded
                                 (only run when n_Hom ≥ MIN_CARRIERS)

    This mirrors the genotypic model in script 15 and allows us to determine
    whether the H2 signal is driven by heterozygous or homozygous carriers.
    """
    print(f"=== Analysis 3: Genotype-dose — {cohort_label} ===")
    df = df[~df["is_replicate"] & (df["Tissue"] == "Tumour")].copy()
    cohort_filter = "Breast" if "Breast" in cohort_label else "Endometri"
    df = df[df["Cohort"].str.contains(cohort_filter, case=False, na=False)]
    if "GT" not in df.columns:
        print("  No GT column — skipping.\n"); return pd.DataFrame()

    def _geno_contrast(g0: pd.Series, g1: pd.Series, col: str, meta: Dict,
                       label: str) -> Optional[Dict]:
        """Run a single binary contrast (g1 vs g0 = reference) for one variable."""
        c0 = pd.to_numeric(g0[col], errors="coerce").dropna() if meta["type"] == "continuous" else g0[col].dropna()
        c1 = pd.to_numeric(g1[col], errors="coerce").dropna() if meta["type"] == "continuous" else g1[col].dropna()
        if len(c0) < MIN_CARRIERS or len(c1) < MIN_CARRIERS:
            return None
        if meta["type"] == "continuous":
            stat, p = mannwhitneyu(c1.values, c0.values, alternative="two-sided")
            return {"Test_Contrast": f"Mann-Whitney U ({label})",
                    "N_ref": len(c0), "N_test": len(c1), "Stat_Contrast": round(stat, 3),
                    "P_Contrast": p}
        else:
            cats = sorted(set(c0.tolist() + c1.tolist()))
            if len(cats) < 2:
                return None
            ct = pd.crosstab(
                pd.concat([pd.Series(["ref"] * len(c0), name="grp"),
                           pd.Series(["test"] * len(c1), name="grp")]),
                pd.concat([c0, c1])
            )
            if ct.shape[0] < 2 or ct.shape[1] < 2:
                return None
            _, p = fisher_exact(ct.values) if ct.shape == (2, 2) else (
                stats.chi2_contingency(ct)[:2]
            )
            return {"Test_Contrast": f"Fisher/Chi-sq ({label})",
                    "N_ref": len(c0), "N_test": len(c1), "Stat_Contrast": np.nan,
                    "P_Contrast": p}

    rows = []
    for var_id in top_variants:
        vdf = df[df["Variant_ID"] == var_id]
        if vdf.empty: continue
        sym = vdf["SYMBOL"].dropna().iloc[0] if vdf["SYMBOL"].notna().any() else ""
        sd  = vdf.drop_duplicates("Sample").set_index("Sample").copy()
        sd["dose"] = sd["GT"].map({"0/0": 0, "0/1": 1, "1/1": 2})

        g0 = sd[sd["dose"] == 0]
        g1 = sd[sd["dose"] == 1]
        g2 = sd[sd["dose"] == 2]

        for col, meta in var_dict.items():
            if col not in sd.columns: continue

            # (a) Overall trend test
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

            base_row = {
                "Cohort": cohort_label, "Variant_ID": var_id, "Gene": sym,
                "Clinical_Var": col, "Clin_Label": meta["label"],
                "N_WT": len(g0), "N_Het": len(g1), "N_Hom": len(g2),
                "Test_Trend": test,
                "Stat_Trend": round(stat, 3), "P_Trend": p,
            }

            # (b) Het vs WT
            het_res = _geno_contrast(g0, g1, col, meta, label="Het_vs_WT")
            base_row["Test_Het_vs_WT"] = het_res["Test_Contrast"]  if het_res else "Skipped"
            base_row["P_Het_vs_WT"]    = het_res["P_Contrast"]     if het_res else np.nan
            base_row["N_Het_used"]     = het_res["N_test"]         if het_res else 0

            # (c) Hom vs WT (only if enough Hom)
            if len(g2) >= MIN_CARRIERS:
                hom_res = _geno_contrast(g0, g2, col, meta, label="Hom_vs_WT")
                base_row["Test_Hom_vs_WT"] = hom_res["Test_Contrast"] if hom_res else "Skipped"
                base_row["P_Hom_vs_WT"]    = hom_res["P_Contrast"]    if hom_res else np.nan
                base_row["N_Hom_used"]     = hom_res["N_test"]        if hom_res else 0
            else:
                base_row["Test_Hom_vs_WT"] = f"Skipped (n_Hom={len(g2)} < {MIN_CARRIERS})"
                base_row["P_Hom_vs_WT"]    = np.nan
                base_row["N_Hom_used"]     = len(g2)

            rows.append(base_row)

    res = pd.DataFrame(rows)
    if not res.empty:
        # FDR correction on the primary trend p-value (used as the headline test)
        res = _apply_fdr(res, p_col="P_Trend", out_col="FDR_Trend")
        res["Nominal_Sig"] = res["P_Trend"] < 0.05
        res["FDR_Sig"]     = res["FDR_Trend"] < FDR_THRESHOLD
        # FDR for the Het/Hom contrasts separately
        for contrast_p, contrast_fdr in [("P_Het_vs_WT", "FDR_Het_vs_WT"),
                                          ("P_Hom_vs_WT", "FDR_Hom_vs_WT")]:
            res = _apply_fdr(res, p_col=contrast_p, out_col=contrast_fdr)
        # Rename P_Value to P_Trend for backwards-compat downstream plotting
        res["P_Value"] = res["P_Trend"]
        res["FDR_P_Value"] = res["FDR_Trend"]
        print(f"  {len(res)} tests | {res['Nominal_Sig'].sum()} nominal (trend) | {res['FDR_Sig'].sum()} FDR\n")
    else:
        print("  No results.\n")
    return res


# ── ANALYSIS 4: SURVIVAL — ALL COHORTS WITH SURVIVAL DATA ────────────────────
#
# Cohorts and endpoints:
#   Endometrial (AT=AUs) — OS  : canon__os_months  + ENDO_EXITUS_BIN
#                            — PFS : canon__pfs_months + ENDO_PD_BIN
#   Breast (MT-T_N)  — OS  : BREAST_OS_MONTHS_DERIVED + BREAST_EXITUS_DERIVED
#
# Model:
#   (a) Log-rank test: carrier (any GT) vs non-carrier — quick screening
#   (b) Age- and BMI-adjusted Cox PH — additive model: snp_carrier + age_z + bmi_z
#   (c) Genotypic Cox PH    — WT reference, Het and Hom as separate dummy terms
#                             (only when n_Hom ≥ MIN_CARRIERS)
#       covariates: age_z, bmi_z
#
# KM curves stratified by genotype class (WT / Het / Hom), one page per
# variant × endpoint.  All pages written to a single PDF.
#
# FDR correction: BH applied within each cohort × endpoint combination.

SURVIVAL_COHORTS = {
    "Endometrial": {
        "cohort_filter": "Endometri",
        "endpoints": {
            "OS":  {"t_col": "canon__os_months",        "ev_col": "ENDO_EXITUS_BIN"},
            "PFS": {"t_col": "canon__pfs_months",        "ev_col": "ENDO_PD_BIN"},
        },
        "age_col": "canon__age",
    },
    "Breast": {
        "cohort_filter": "Breast",
        "endpoints": {
            "OS":  {"t_col": "BREAST_OS_MONTHS_DERIVED", "ev_col": "BREAST_EXITUS_DERIVED"},
        },
        "age_col": "canon__age",
    },
}

_GENO_COLS  = {"WT": "#1976D2", "Het": "#F57C00", "Hom": "#C62828"}


def survival_analysis(df: pd.DataFrame, out_dir: Optional[Path]=None) -> pd.DataFrame:
    if not _HAS_LIFELINES:
        print("=== Analysis 4: Survival — SKIPPED (lifelines not installed) ===\n")
        return pd.DataFrame(), []

    print("=== Analysis 4: Survival — all cohorts ===")
    df = df[~df["is_replicate"] & (df["Tissue"] == "Tumour")].copy()
    if "GT" not in df.columns:
        print("  No GT column — genotypic model unavailable.\n")
        return pd.DataFrame(), []

    rows      = []
    km_pdf = None
    km_n_pages = 0
    if out_dir is not None and _HAS_LIFELINES:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        km_pdf_path = out_dir / "17_KM_Curves_All_Cohorts.pdf"
        km_pdf = pdf_backend.PdfPages(km_pdf_path)


    for cohort_label, cohort_cfg in SURVIVAL_COHORTS.items():
        cdf = df[df["Cohort"].str.contains(cohort_cfg["cohort_filter"],
                                            case=False, na=False)].copy()
        if cdf.empty:
            print(f"  {cohort_label}: no tumour samples — skipping.")
            continue

        age_col = cohort_cfg["age_col"]
        print(f"  {cohort_label}: {cdf['Sample'].nunique()} tumour samples")

        for endpoint, ep_cfg in cohort_cfg["endpoints"].items():
            t_col  = ep_cfg["t_col"]
            ev_col = ep_cfg["ev_col"]

            if t_col not in cdf.columns or ev_col not in cdf.columns:
                print(f"    {endpoint}: columns missing — skipping.")
                continue

            # One row per sample with genotype
            sd = cdf.drop_duplicates("Sample").set_index("Sample").copy()
            sd["dose"] = sd["GT"].map({"0/0": 0, "0/1": 1, "1/1": 2})
            sd["geno"] = sd["dose"].map({0: "WT", 1: "Het", 2: "Hom"})

            print(f"    {endpoint} …")

            for var_id, vdf in cdf.groupby("Variant_ID"):
                sym = vdf["SYMBOL"].dropna().iloc[0] if vdf["SYMBOL"].notna().any() else ""
                carrier_set = set(vdf["Sample"].unique())

                # Build per-sample frame for this variant
                s = sd.copy()
                # Samples not in carrier_set → WT (dose=0)
                s.loc[~s.index.isin(carrier_set), ["dose", "geno"]] = [0, "WT"]
                # Samples in carrier_set keep their GT-derived genotype

                # Required columns
                needed = [t_col, ev_col]
                if age_col in s.columns: needed.append(age_col)
                bmi_col = _choose_bmi_col(s, cohort_label)
                if bmi_col: needed.append(bmi_col)
                valid = s[list(dict.fromkeys(needed + ["dose", "geno"]))].copy()
                valid[t_col]  = pd.to_numeric(valid[t_col],  errors="coerce")
                valid[ev_col] = pd.to_numeric(valid[ev_col], errors="coerce")
                valid = valid.dropna(subset=[t_col, ev_col])
                # Guard: Cox/log-rank need events; skip extremely low-event fits
                if valid[ev_col].sum() < 2:
                    continue
                valid = valid[valid[t_col] > 0]   # lifelines requires T > 0
                if len(valid) < 6: continue

                T    = valid[t_col]
                E    = valid[ev_col].astype(float)
                geno = valid["geno"].fillna("WT")
                dose = valid["dose"].fillna(0).astype(int)

                n_wt  = (geno == "WT").sum()
                n_het = (geno == "Het").sum()
                n_hom = (geno == "Hom").sum()
                n_carriers = n_het + n_hom

                if n_carriers < MIN_CARRIERS or n_wt < MIN_CARRIERS:
                    continue

                # ── (a) Log-rank: any carrier vs WT ──────────────────────
                carrier_flag = (dose > 0).astype(int)
                try:
                    lr = logrank_test(T[carrier_flag==1], T[carrier_flag==0],
                                      E[carrier_flag==1], E[carrier_flag==0])
                    lr_p = round(lr.p_value, 4)
                except Exception:
                    lr_p = np.nan

                # ── (b) Age- and BMI-adjusted additive Cox ───────────────────
                cox_add_hr = cox_add_ci = cox_add_p = np.nan
                cox_add_note = ""
                try:
                    # Guard against (quasi-)complete separation in carrier vs WT:
                    # if any of the 2×2 cells is zero (events/non-events in either group),
                    # Cox estimates can become non-identifiable and unstable.
                    ev_car  = int(E[carrier_flag==1].sum())
                    ev_wt   = int(E[carrier_flag==0].sum())
                    ne_car  = int((carrier_flag==1).sum() - ev_car)
                    ne_wt   = int((carrier_flag==0).sum() - ev_wt)
                    if min(ev_car, ev_wt, ne_car, ne_wt) == 0:
                        cox_add_note = "skipped_separation(carrier_vs_WT)"
                    else:
                        cox_df = pd.DataFrame({"T": T, "E": E, "carrier": carrier_flag})
                        if age_col in valid.columns:
                            a = pd.to_numeric(valid[age_col], errors="coerce")
                            cox_df["age_z"] = (a - a.mean()) / a.std()
                            formula_cols = ["carrier", "age_z"]
                        else:
                            formula_cols = ["carrier"]
                        if bmi_col and bmi_col in valid.columns:
                            b_vals = pd.to_numeric(valid[bmi_col], errors="coerce")
                            cox_df["bmi_z"] = (b_vals - b_vals.mean()) / b_vals.std()
                            formula_cols.append("bmi_z")
                        cox_df = cox_df.dropna()
                        if len(cox_df) >= 6:
                            cph = CoxPHFitter(penalizer=0.1)
                            cph.fit(cox_df, duration_col="T", event_col="E",
                                    formula=" + ".join(formula_cols), show_progress=False)
                            cox_add_hr = round(float(_safe_exp(cph.params_["carrier"])), 3)
                            ci = _safe_exp(cph.confidence_intervals_.loc["carrier"])
                            cox_add_p  = round(float(cph.summary.loc["carrier", "p"]), 4)
                        try:
                            lo = float(ci.iloc[0]) if hasattr(ci, "iloc") else float(ci[0])
                            hi = float(ci.iloc[1]) if hasattr(ci, "iloc") else float(ci[1])
                            cox_add_ci = f"[{lo:.3f}, {hi:.3f}]"
                        except Exception:
                            cox_add_ci = np.nan
                except Exception:
                    cox_add_note = cox_add_note or "fit_error"
                    if age_col in valid.columns:
                        a = pd.to_numeric(valid[age_col], errors="coerce")
                        cox_df["age_z"] = (a - a.mean()) / a.std()
                        formula_cols = ["carrier", "age_z"]
                    else:
                        formula_cols = ["carrier"]
                    if bmi_col and bmi_col in valid.columns:
                        b_vals = pd.to_numeric(valid[bmi_col], errors="coerce")
                        cox_df["bmi_z"] = (b_vals - b_vals.mean()) / b_vals.std()
                        formula_cols.append("bmi_z")
                    cox_df = cox_df.dropna()
                    if len(cox_df) >= 6:
                        cph = CoxPHFitter()
                        cph.fit(cox_df, duration_col="T", event_col="E",
                                formula=" + ".join(formula_cols), show_progress=False)
                        cox_add_hr = round(float(_safe_exp(cph.params_["carrier"])), 3)
                        ci = _safe_exp(cph.confidence_intervals_.loc["carrier"])
                        cox_add_p  = round(float(cph.summary.loc["carrier", "p"]), 4)
                        try:
                            lo = float(ci.iloc[0]) if hasattr(ci, "iloc") else float(ci[0])
                            hi = float(ci.iloc[1]) if hasattr(ci, "iloc") else float(ci[1])
                            cox_add_ci = f"[{lo:.3f}, {hi:.3f}]"
                        except Exception:
                            cox_add_ci = np.nan
                except Exception:
                    pass

                # ── (c) Genotypic Cox (WT reference, Het + Hom terms) ────
                geno_het_hr = geno_het_ci = geno_het_p = np.nan
                geno_hom_hr = geno_hom_ci = geno_hom_p = np.nan
                model_type = "Additive_only"

                def _term_ok(glabel: str) -> bool:
                    try:
                        e  = int(E[geno == glabel].sum())
                        ne = int((E[geno == glabel] == 0).sum())
                        return (e >= 1) and (ne >= 1)
                    except Exception:
                        return False

                geno_terms = []
                if n_het >= MIN_CARRIERS and _term_ok("Het"):
                    geno_terms.append("Het")
                if n_hom >= MIN_CARRIERS and _term_ok("Hom"):
                    geno_terms.append("Hom")

                if len(geno_terms) > 0:
                    model_type = "Genotypic_" + "_".join(geno_terms)
                    try:
                        gdf = pd.DataFrame({
                            "T": T, "E": E,
                            "Het": (geno == "Het").astype(int),
                            "Hom": (geno == "Hom").astype(int),
                        })
                        geno_formula = " + ".join(geno_terms)

                        if age_col in valid.columns:
                            a = pd.to_numeric(valid[age_col], errors="coerce")
                            gdf["age_z"] = (a - a.mean()) / a.std()
                            geno_formula += " + age_z"
                        if bmi_col and bmi_col in valid.columns:
                            b_vals = pd.to_numeric(valid[bmi_col], errors="coerce")
                            gdf["bmi_z"] = (b_vals - b_vals.mean()) / b_vals.std()
                            geno_formula += " + bmi_z"

                        gdf = gdf.dropna()
                        if len(gdf) >= 6 and gdf["E"].sum() >= 2:
                            cph_g = CoxPHFitter(penalizer=0.1)
                            cph_g.fit(gdf, duration_col="T", event_col="E",
                                      formula=geno_formula, show_progress=False)

                            if "Het" in geno_terms and "Het" in cph_g.params_.index:
                                hr_v = round(float(_safe_exp(cph_g.params_["Het"])), 3)
                                ci_v = _safe_exp(cph_g.confidence_intervals_.loc["Het"])
                                p_v  = round(float(cph_g.summary.loc["Het", "p"]), 4)
                                geno_het_hr = hr_v
                                geno_het_ci = f"[{ci_v.iloc[0]:.3f}, {ci_v.iloc[1]:.3f}]"
                                geno_het_p  = p_v

                            if "Hom" in geno_terms and "Hom" in cph_g.params_.index:
                                hr_v = round(float(_safe_exp(cph_g.params_["Hom"])), 3)
                                ci_v = _safe_exp(cph_g.confidence_intervals_.loc["Hom"])
                                p_v  = round(float(cph_g.summary.loc["Hom", "p"]), 4)
                                geno_hom_hr = hr_v
                                geno_hom_ci = f"[{ci_v.iloc[0]:.3f}, {ci_v.iloc[1]:.3f}]"
                                geno_hom_p  = p_v
                    except Exception:
                        model_type = "Additive_only"
                elif n_het >= MIN_CARRIERS:
                    model_type = "Additive_only"

                rows.append({
                    "Cohort":           cohort_label,
                    "Endpoint":         endpoint,
                    "Variant_ID":       var_id,
                    "Gene":             sym,
                    "N_total":          len(valid),
                    "N_WT":             int(n_wt),
                    "N_Het":            int(n_het),
                    "N_Hom":            int(n_hom),
                    "N_events":         int(E.sum()),
                    "Model_Type":       model_type,
                    # Log-rank
                    "Logrank_P":        lr_p,
                    # Additive Cox (any carrier vs WT, age-adjusted)
                    "HR_Additive":      cox_add_hr,
                    "HR_Additive_CI95": cox_add_ci,
                    "P_Cox_Additive":   cox_add_p,
                    "Additive_Note":   cox_add_note,
                    # Genotypic Cox (Het vs WT)
                    "HR_Het":           geno_het_hr,
                    "HR_Het_CI95":      geno_het_ci,
                    "P_Cox_Het":        geno_het_p,
                    # Genotypic Cox (Hom vs WT)
                    "HR_Hom":           geno_hom_hr,
                    "HR_Hom_CI95":      geno_hom_ci,
                    "P_Cox_Hom":        geno_hom_p,
                })

                # ── KM plot — stratified by genotype class ─────────────
                fig, ax = plt.subplots(figsize=(7, 4.5))
                kmf = KaplanMeierFitter()
                for geno_label, col in _GENO_COLS.items():
                    mask = geno == geno_label
                    if mask.sum() < 2: continue
                    kmf.fit(T[mask], E[mask],
                            label=f"{geno_label} (n={mask.sum()})")
                    kmf.plot_survival_function(ax=ax, ci_show=True, color=col)
                title = (f"{var_id}  [{sym}]  —  {cohort_label} | {endpoint}\n"
                         f"Log-rank p={lr_p}  |  "
                         f"WT n={n_wt}, Het n={n_het}, Hom n={n_hom}")
                ax.set_title(title, fontsize=8.5)
                ax.set_xlabel(f"{endpoint} (months)")
                ax.set_ylabel("Survival probability")
                ax.legend(fontsize=8)
                _style_ax(ax)
                plt.tight_layout()
                if km_pdf is not None: km_pdf.savefig(fig, bbox_inches="tight"); km_n_pages += 1
                plt.close(fig)

    res = pd.DataFrame(rows)
    if not res.empty:
        # FDR within each cohort × endpoint × model term
        for cohort in res["Cohort"].unique():
            for ep in res["Endpoint"].unique():
                mask = (res["Cohort"] == cohort) & (res["Endpoint"] == ep)
                sub  = res[mask].copy()
                if sub.empty: continue
                for p_col, fdr_col in [
                    ("P_Cox_Additive", "FDR_Cox_Additive"),
                    ("P_Cox_Het",      "FDR_Cox_Het"),
                    ("P_Cox_Hom",      "FDR_Cox_Hom"),
                ]:
                    res.loc[mask, fdr_col] = _apply_fdr(sub, p_col, fdr_col)[fdr_col]
        n_nom = (res["P_Cox_Additive"] < 0.05).sum()
        print(f"  {len(res)} tests across {res['Cohort'].nunique()} cohorts, "
              f"{res['Endpoint'].nunique()} endpoints | "
              f"{n_nom} nominal (additive Cox)\n")
    else:
        print("  No survival results (check OS/PFS columns and event coding).\n")

    # Close KM PDF if we wrote any pages
    if km_pdf is not None:
        km_pdf.close()
        print(f"  Saved: {km_pdf_path}  ({km_n_pages} pages)")
    return res


# ── VISUALISATION ─────────────────────────────────────────────────────────────
# Shared style
_PALETTE   = {"carrier": "#D32F2F", "non_carrier": "#1976D2"}
_IMPACT_C  = {"HIGH": "#b71c1c", "MODERATE": "#e65100", "LOW": "#2e7d32", "MODIFIER": "#78909c"}
_COHORT_C  = {"Breast": "#AD1457", "Endometrial": "#00695C"}

def _style_ax(ax, grid=True):
    """Apply consistent clean style to an axis."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.yaxis.grid(True, linestyle=":", color="#cccccc", alpha=0.7)
        ax.set_axisbelow(True)

def _sig_label(p):
    """Return asterisk string for a p-value."""
    if pd.isna(p):   return ""
    if p < 0.001:    return "***"
    if p < 0.01:     return "**"
    if p < 0.05:     return "*"
    return "ns"


# ── 1. VOLCANO (tumour vs healthy) ───────────────────────────────────────────

def _draw_volcano(tvh, p_col, sig_col, thresh, title_suffix, out_path):
    groups = sorted(tvh["Analysis_Group"].unique())
    ncols  = len(groups)
    fig, axes = plt.subplots(1, ncols, figsize=(7.5 * ncols, 6.5), squeeze=False)
    axes = axes[0]

    for ax, grp in zip(axes, groups):
        sub = tvh[tvh["Analysis_Group"] == grp].copy()
        sub["-log10p"] = -np.log10(pd.to_numeric(sub[p_col], errors="coerce").clip(lower=1e-10))
        sub["log2OR"]  = np.log2(
            sub["Odds_Ratio"].replace({0: 0.001, np.inf: 1000}).clip(0.001, 1000)
        )
        sig_line = -np.log10(thresh)

        # Background shading: enriched (right) vs depleted (left) in tumour
        ax.axvspan( 0,  ax.get_xlim()[1] if ax.get_xlim()[1] > 0 else 10,
                   alpha=0.04, color="#D32F2F")
        ax.axvspan(ax.get_xlim()[0] if ax.get_xlim()[0] < 0 else -10, 0,
                   alpha=0.04, color="#1976D2")

        # Points — coloured by VEP impact, sized by -log10p
        for imp, col in _IMPACT_C.items():
            m = sub["IMPACT"] == imp
            if not m.any(): continue
            sizes = np.clip(sub.loc[m, "-log10p"] * 18, 30, 220)
            ax.scatter(sub.loc[m, "log2OR"], sub.loc[m, "-log10p"],
                       c=col, s=sizes, alpha=0.82,
                       edgecolors="white", linewidths=0.5,
                       zorder=3, label=imp)

        # Significance threshold line
        ax.axhline(sig_line, color="#c62828", linestyle="--", linewidth=1.2,
                   alpha=0.85, zorder=2,
                   label=f"p = {thresh}" if thresh == 0.05 else f"FDR = {thresh}")
        ax.axvline(0, color="#555555", linestyle="--", linewidth=0.9, alpha=0.5, zorder=2)

        # Annotate top hits — smart repulsion to avoid overlap
        hits = sub[sub[sig_col]].nsmallest(10, p_col)
        for _, row in hits.iterrows():
            ax.annotate(
                row.get("Gene", row["Variant_ID"]) or row["Variant_ID"],
                xy=(row["log2OR"], row["-log10p"]),
                xytext=(14, 6), textcoords="offset points",
                fontsize=7.5, fontweight="bold",
                color="#212121",
                bbox=dict(boxstyle="round,pad=0.25", fc="white",
                          ec="#bdbdbd", alpha=0.9, linewidth=0.8),
                arrowprops=dict(arrowstyle="-", color="#888888",
                                lw=0.8, connectionstyle="arc3,rad=0.1"),
            )

        # Axis labels and styling
        ax.set_xlabel("log₂(Odds Ratio)  [tumour enriched →]", fontsize=10)
        ax.set_ylabel("-log₁₀(p-value)", fontsize=10)
        ax.set_title(grp, fontsize=12, fontweight="bold", pad=8)
        _style_ax(ax, grid=False)
        ax.xaxis.grid(True, linestyle=":", color="#cccccc", alpha=0.5)
        ax.yaxis.grid(True, linestyle=":", color="#cccccc", alpha=0.5)
        ax.set_axisbelow(True)

        # Marginal rug showing OR distribution
        ax.plot(sub["log2OR"], np.full(len(sub), ax.get_ylim()[0]),
                "|", color="#aaaaaa", alpha=0.5, markersize=4, zorder=1)

        # N label
        n_sig = int(sub[sig_col].sum())
        ax.text(0.98, 0.02, f"n variants = {len(sub)}\nn sig. = {n_sig}",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, color="#555555",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8, ec="none"))

    # Shared legend
    impact_handles = [mpatches.Patch(facecolor=c, label=i,
                                     edgecolor="white", linewidth=0.5)
                      for i, c in _IMPACT_C.items()]
    fig.legend(handles=impact_handles, title="VEP Impact",
               bbox_to_anchor=(1.01, 0.5), loc="center left",
               frameon=True, framealpha=0.9, edgecolor="#cccccc")

    fig.suptitle(f"GSDMB SNPs — Tumour vs Healthy  ({title_suffix})",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def make_volcano_plots(tvh, out_dir):
    if tvh.empty: return
    _draw_volcano(tvh, "P_Value",     "Nominal_Sig", 0.05,
                  "unadjusted p",  out_dir / "17_SNP_Volcano_raw_p.png")
    _draw_volcano(tvh, "FDR_P_Value", "FDR_Sig",     FDR_THRESHOLD,
                  "FDR-corrected", out_dir / "17_SNP_Volcano_FDR.png")


# ── 2. HEATMAP (SNP × clinical variable) ─────────────────────────────────────

def _draw_heatmap(clin_res, p_col, title_suffix, out_path, filter_thresh=0.20):
    pivot = (clin_res
             .pivot_table(index="Variant_ID", columns="Clin_Label",
                          values=p_col, aggfunc="min")
             .dropna(how="all"))
    pivot = pivot[(pivot < filter_thresh).any(axis=1)]
    if pivot.empty:
        print(f"  Heatmap ({title_suffix}): no hits below {filter_thresh} — skipped")
        return

    # Sort rows by min p-value; columns by number of hits
    pivot = pivot.loc[pivot.min(axis=1).sort_values().index]
    pivot = pivot[pivot.notna().sum().sort_values(ascending=False).index]

    log_piv = -np.log10(pivot.fillna(1).clip(lower=1e-10))
    sig_t   = -np.log10(0.05 if "Raw" in title_suffix else FDR_THRESHOLD)

    cell_h  = 0.65
    cell_w  = 1.8
    fig_h   = max(6, len(pivot) * cell_h + 3)
    fig_w   = max(10, len(pivot.columns) * cell_w + 4)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Draw heatmap
    sns.heatmap(log_piv, ax=ax,
                cmap="RdPu", vmin=0, vmax=max(4, log_piv.max().max()),
                linewidths=0.4, linecolor="#e0e0e0",
                annot=False,
                cbar_kws={"label": "-log₁₀(p-value)", "shrink": 0.6})

    # Manual cell annotations: p-value + significance stars
    for i, var in enumerate(log_piv.index):
        for j, clin in enumerate(log_piv.columns):
            raw_p = pivot.loc[var, clin]
            lp    = log_piv.loc[var, clin]
            if pd.isna(raw_p): continue
            star  = _sig_label(raw_p)
            txt_c = "white" if lp >= sig_t else "#333333"
            ax.text(j + 0.5, i + 0.38,
                    f"{raw_p:.3f}" if raw_p >= 0.001 else f"{raw_p:.1e}",
                    ha="center", va="center", fontsize=8, color=txt_c)
            if star not in ("", "ns"):
                ax.text(j + 0.5, i + 0.72, star,
                        ha="center", va="center",
                        fontsize=10, color="white" if lp >= sig_t else "#c62828",
                        fontweight="bold")

    # Add significance threshold annotation to colourbar
    cbar = ax.collections[0].colorbar
    cbar.ax.axhline(sig_t, color="#c62828", linewidth=1.5, linestyle="--")
    cbar.ax.text(1.05, sig_t / cbar.ax.get_ylim()[1],
                 f" p<{'0.05' if 'Raw' in title_suffix else str(FDR_THRESHOLD)}",
                 transform=cbar.ax.transAxes, va="center",
                 fontsize=7, color="#c62828")

    ax.set_title(f"GSDMB SNPs × Clinical Variables — {title_suffix}",
                 fontsize=11, fontweight="bold", pad=10)
    ax.set_xlabel("Clinical variable", fontsize=9, labelpad=8)
    ax.set_ylabel("Variant", fontsize=9, labelpad=8)
    plt.xticks(rotation=40, ha="right", fontsize=10)
    plt.yticks(fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def make_heatmaps(breast_clin, endo_clin, out_dir):
    for clin_res, label in [(breast_clin, "Breast"), (endo_clin, "Endometrial")]:
        if clin_res.empty: continue
        _draw_heatmap(clin_res, "P_Unadj",
                      f"Unadjusted p — {label}",
                      out_dir / f"17_SNP_Heatmap_{label}_raw_p.png",
                      filter_thresh=0.10)
        _draw_heatmap(clin_res, "FDR_Unadj",
                      f"FDR-corrected — {label}",
                      out_dir / f"17_SNP_Heatmap_{label}_FDR.png",
                      filter_thresh=0.30)
        adj_res = clin_res[clin_res["P_Adj_Age"].notna()].copy()
        if not adj_res.empty:
            _draw_heatmap(adj_res, "P_Adj_Age",
                          f"Age-adjusted p — {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_age_raw_p.png",
                          filter_thresh=0.10)
            _draw_heatmap(adj_res, "FDR_Adj_Age",
                          f"Age-adjusted FDR — {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_age_FDR.png",
                          filter_thresh=0.30)
        bmi_res = clin_res[clin_res["P_Adj_AgeBMI"].notna()].copy()
        if not bmi_res.empty:
            _draw_heatmap(bmi_res, "P_Adj_AgeBMI",
                          f"Age+BMI-adjusted p — {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_agebmi_raw_p.png",
                          filter_thresh=0.10)
            _draw_heatmap(bmi_res, "FDR_Adj_AgeBMI",
                          f"Age+BMI-adjusted FDR — {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_agebmi_FDR.png",
                          filter_thresh=0.30)


# ── 3. FOREST PLOT (OR with 95% CI) ──────────────────────────────────────────

def make_forest_plots(breast_clin, endo_clin, out_dir):
    """
    One forest plot per cohort showing OR ± 95% CI for all binary/nominal
    associations that reached nominal significance in the age-adjusted model.
    """
    for clin_res, label, colour in [
        (breast_clin,  "Breast",       _COHORT_C["Breast"]),
        (endo_clin,    "Endometrial",  _COHORT_C["Endometrial"]),
    ]:
        if clin_res.empty: continue

        # Keep only binary tests with a valid age-adjusted OR and CI
        # (age-only model used as primary adjusted result — full N)
        sub = clin_res[
            (clin_res["Clin_Type"] == "binary") &
            clin_res["OR_Adj_Age"].notna() &
            clin_res["OR_Adj_Age_CI95"].notna()
        ].copy()
        if sub.empty: continue

        # Parse CI strings "[lo, hi]"
        def _parse_ci(s):
            try:
                lo, hi = re.findall(r"[-0-9.eE+]+", str(s))[:2]
                return float(lo), float(hi)
            except Exception:
                return np.nan, np.nan
        sub[["CI_lo", "CI_hi"]] = sub["OR_Adj_Age_CI95"].apply(
            lambda s: pd.Series(_parse_ci(s))
        )
        sub = sub.dropna(subset=["CI_lo", "CI_hi"])
        if sub.empty: continue

        # Sort by p-value; label = "Gene: clinical var"
        sub = sub.sort_values("P_Adj_Age")
        sub["label"] = sub.apply(
            lambda r: f"{r['Gene'] or r['Variant_ID']} — {r['Clin_Label']}", axis=1
        )

        # Cap extreme OR/CI values for display (keep raw values in data)
        OR_CAP = 20.0
        sub["OR_plot"]    = sub["OR_Adj_Age"].clip(upper=OR_CAP)
        sub["CI_lo_plot"] = sub["CI_lo"].clip(lower=0)
        sub["CI_hi_plot"] = sub["CI_hi"].clip(upper=OR_CAP)
        sub["capped"]     = sub["CI_hi"] > OR_CAP

        n_rows = len(sub)
        fig_h  = max(5, n_rows * 0.42 + 2)
        fig, ax = plt.subplots(figsize=(9, fig_h))

        ypos = np.arange(n_rows)[::-1]   # top = most significant

        # Colour by significance
        sig_mask = sub["P_Adj_Age"] < 0.05
        point_c  = [colour if s else "#999999" for s in sig_mask]

        # Error bars (capped)
        xerr_lo = (sub["OR_plot"] - sub["CI_lo_plot"]).values
        xerr_hi = (sub["CI_hi_plot"] - sub["OR_plot"]).values
        ax.errorbar(sub["OR_plot"].values, ypos,
                    xerr=[np.clip(xerr_lo, 0, None), np.clip(xerr_hi, 0, None)],
                    fmt="none", ecolor="#aaaaaa", elinewidth=1.2, capsize=3, zorder=2)

        # Points
        for i, (_, row) in enumerate(sub.iterrows()):
            ax.scatter(row["OR_plot"], ypos[i],
                       s=90, color=point_c[i],
                       edgecolors="white", linewidths=0.6, zorder=3)
            # Arrow marker if CI was capped
            if row["capped"]:
                ax.annotate("→", xy=(OR_CAP, ypos[i]),
                            fontsize=9, color="#aaaaaa", va="center")

        # p-value labels on the right
        ax.set_xlim(left=0, right=OR_CAP * 1.15)
        x_max = OR_CAP * 1.05
        for i, (_, row) in enumerate(sub.iterrows()):
            pstr = f"p={row['P_Adj_Age']:.3f}{_sig_label(row['P_Adj_Age'])}"
            ax.text(x_max, ypos[i], pstr,
                    va="center", ha="left", fontsize=7, color="#444444")

        ax.axvline(1, color="#555555", linestyle="--", linewidth=1, alpha=0.7)
        ax.set_yticks(ypos)
        ax.set_yticklabels(sub["label"].values, fontsize=8)
        ax.set_xlabel("Odds Ratio (age-adjusted, 95% CI)", fontsize=10)
        ax.set_title(f"{label} — Age-adjusted OR for GSDMB carrier status",
                     fontsize=11, fontweight="bold", pad=10)
        _style_ax(ax, grid=False)
        ax.xaxis.grid(True, linestyle=":", color="#cccccc", alpha=0.6)
        ax.set_axisbelow(True)

        # Shade OR > 1 region
        xlims = ax.get_xlim()
        ax.axvspan(1, xlims[1], alpha=0.04, color="#c62828")
        ax.axvspan(xlims[0], 1, alpha=0.04, color="#1976D2")
        ax.text(0.97, 0.01, "← protective", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7.5, color="#1976D2", style="italic")
        ax.text(0.99, 0.01, "risk →",       transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7.5, color="#c62828", style="italic")

        plt.tight_layout()
        out_path = out_dir / f"17_Forest_{label}.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {out_path}")


# ── 4. BOXPLOT / VIOLIN STRIP for top hits ────────────────────────────────────

def make_distribution_plots(breast_clin, endo_clin, merged, out_dir):
    """
    For the top N nominally significant continuous associations, draw a
    split violin + strip plot comparing carriers vs non-carriers.
    """
    TOP_N = 12

    for clin_res, label, coh_sheet, colour in [
        (breast_clin, "Breast",      "MT-T_N",  _COHORT_C["Breast"]),
        (endo_clin,   "Endometrial", "AT=AUs",  _COHORT_C["Endometrial"]),
    ]:
        if clin_res.empty: continue

        cont = (clin_res[
            (clin_res["Clin_Type"] == "continuous") &
            clin_res["P_Unadj"].notna()
        ].sort_values("P_Unadj")
         .drop_duplicates("Clinical_Var")
         .head(TOP_N))
        if cont.empty: continue

        n_plots = len(cont)
        ncols   = min(3, n_plots)
        nrows   = int(np.ceil(n_plots / ncols))
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(5.5 * ncols, 4.5 * nrows),
                                 squeeze=False)
        axes_flat = axes.flatten()

        cohort_df = merged[
            (merged["sheet"] == coh_sheet) & (~merged["is_replicate"])
        ].drop_duplicates("Sample")

        for idx, (_, row) in enumerate(cont.iterrows()):
            ax      = axes_flat[idx]
            col     = row["Clinical_Var"]
            var_id  = row["Variant_ID"]
            p_u     = row["P_Unadj"]
            p_a     = row.get("P_Adj_Age", np.nan)

            # Carrier flag for this specific variant
            var_samples = set(
                merged[(merged["Variant_ID"] == var_id) &
                        (merged["sheet"] == coh_sheet)]["Sample"].unique()
            )
            cohort_df["_carrier"] = cohort_df.index.map(
                lambda s: "Carrier" if s in var_samples else "Non-carrier"
            ) if cohort_df.index.name == "Sample" else cohort_df["Sample"].map(
                lambda s: "Carrier" if s in var_samples else "Non-carrier"
            )

            plot_df = cohort_df[[col, "_carrier"]].dropna()
            if plot_df.empty or plot_df["_carrier"].nunique() < 2:
                ax.set_visible(False); continue

            # Cap extreme outliers at 1st/99th percentile to prevent axis blow-out
            col_numeric = pd.to_numeric(plot_df[col], errors="coerce")
            p01, p99 = col_numeric.quantile(0.01), col_numeric.quantile(0.99)
            if p99 > p01:  # only cap if range is non-trivial
                plot_df = plot_df.copy()
                plot_df[col] = col_numeric.clip(lower=p01, upper=p99)

            # Violin
            parts = ax.violinplot(
                [plot_df.loc[plot_df["_carrier"] == g, col].values
                 for g in ["Carrier", "Non-carrier"]],
                positions=[0, 1], widths=0.6, showmedians=False,
            )
            for pc, c in zip(parts["bodies"],
                             [_PALETTE["carrier"], _PALETTE["non_carrier"]]):
                pc.set_facecolor(c); pc.set_alpha(0.35); pc.set_edgecolor("none")
            for comp in ("cbars", "cmins", "cmaxes"):
                if comp in parts:
                    parts[comp].set_color("#aaaaaa"); parts[comp].set_linewidth(0.8)

            # Strip jitter
            rng = np.random.default_rng(42)
            for gi, (grp, c) in enumerate([("Carrier", _PALETTE["carrier"]),
                                            ("Non-carrier", _PALETTE["non_carrier"])]):
                vals = plot_df.loc[plot_df["_carrier"] == grp, col].values
                jit  = rng.uniform(-0.08, 0.08, len(vals))
                ax.scatter(np.full(len(vals), gi) + jit, vals,
                           color=c, alpha=0.65, s=28,
                           edgecolors="white", linewidths=0.3, zorder=3)

            # Median line
            for gi, grp in enumerate(["Carrier", "Non-carrier"]):
                med = plot_df.loc[plot_df["_carrier"] == grp, col].median()
                ax.hlines(med, gi - 0.18, gi + 0.18,
                          colors="#222222", linewidths=2, zorder=4)

            # Significance bracket
            y_top = plot_df[col].max()
            y_rng = plot_df[col].max() - plot_df[col].min()
            br_y  = y_top + y_rng * 0.08
            ax.annotate("", xy=(1, br_y), xytext=(0, br_y),
                        arrowprops=dict(arrowstyle="-", color="#555555", lw=1.2))
            p_disp = p_a if pd.notna(p_a) else p_u
            lbl    = _sig_label(p_disp)
            suffix = "adj" if pd.notna(p_a) else "unadj"
            ax.text(0.5, br_y + y_rng * 0.03,
                    f"{lbl}  p={p_disp:.3f} ({suffix})",
                    ha="center", va="bottom", fontsize=7.5, color="#333333")

            # Sample sizes
            for gi, grp in enumerate(["Carrier", "Non-carrier"]):
                n = (plot_df["_carrier"] == grp).sum()
                ax.text(gi, plot_df[col].min() - y_rng * 0.1,
                        f"n={n}", ha="center", va="top", fontsize=7.5, color="#555555")

            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Carrier", "Non-carrier"], fontsize=9)
            ax.set_ylabel(row["Clin_Label"], fontsize=9)
            ax.set_title(f"{row.get('Gene', '') or var_id}", fontsize=9, fontweight="bold")
            _style_ax(ax)

        # Hide unused subplots
        for ax in axes_flat[n_plots:]:
            ax.set_visible(False)

        fig.suptitle(f"{label} — Carrier vs Non-carrier: top continuous associations",
                     fontsize=12, fontweight="bold", y=1.01)
        plt.tight_layout()
        out_path = out_dir / f"17_Distribution_{label}.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {out_path}")


# ── 5. STACKED BAR — binary outcomes ─────────────────────────────────────────

def make_binary_bar_plots(breast_clin, endo_clin, merged, out_dir):
    """
    For the top binary associations, show a stacked proportion bar
    (% positive in carrier vs non-carrier) per variant.
    """
    TOP_N = 12

    for clin_res, label, coh_sheet, colour in [
        (breast_clin, "Breast",      "MT-T_N",  _COHORT_C["Breast"]),
        (endo_clin,   "Endometrial", "AT=AUs",  _COHORT_C["Endometrial"]),
    ]:
        if clin_res.empty: continue

        bin_hits = (clin_res[
            (clin_res["Clin_Type"] == "binary") &
            clin_res["P_Unadj"].notna()
        ].sort_values("P_Unadj")
         .drop_duplicates("Clinical_Var")
         .head(TOP_N))
        if bin_hits.empty: continue

        n_plots = len(bin_hits)
        ncols   = min(3, n_plots)
        nrows   = int(np.ceil(n_plots / ncols))
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(4.5 * ncols, 4.0 * nrows),
                                 squeeze=False,
                                 facecolor="white")
        for ax in axes.flatten():
            ax.set_facecolor("white")
        axes_flat = axes.flatten()

        cohort_df = merged[
            (merged["sheet"] == coh_sheet) & (~merged["is_replicate"])
        ].drop_duplicates("Sample").copy()

        for idx, (_, row) in enumerate(bin_hits.iterrows()):
            ax     = axes_flat[idx]
            col    = row["Clinical_Var"]
            var_id = row["Variant_ID"]
            p_u    = row["P_Unadj"]
            p_a    = row.get("P_Adj_Age", np.nan)
            or_v   = row.get("OR_Adj_Age", row.get("OR_Unadj", np.nan))

            var_samples = set(
                merged[(merged["Variant_ID"] == var_id) &
                        (merged["sheet"] == coh_sheet)]["Sample"].unique()
            )
            cohort_df["_carrier"] = cohort_df["Sample"].map(
                lambda s: "Carrier" if s in var_samples else "Non-carrier"
            )

            plot_df = cohort_df[[col, "_carrier"]].dropna()
            if plot_df.empty or plot_df["_carrier"].nunique() < 2:
                ax.set_visible(False); continue

            groups = ["Carrier", "Non-carrier"]
            pos_rates = []
            ns        = []
            for grp in groups:
                vals = plot_df.loc[plot_df["_carrier"] == grp, col]
                pos_rates.append(100 * vals.mean() if len(vals) else 0)
                ns.append(len(vals))

            xpos = [0, 1]
            bars = ax.bar(xpos, pos_rates,
                          color=[_PALETTE["carrier"], _PALETTE["non_carrier"]],
                          width=0.55, edgecolor="white", linewidth=1.2,
                          alpha=0.85, zorder=3)

            # Complement bars (light)
            for xi, (r, c) in enumerate(zip(pos_rates,
                                            [_PALETTE["carrier"], _PALETTE["non_carrier"]])):
                ax.bar(xi, 100 - r, bottom=r, color=c,
                       width=0.55, alpha=0.18, edgecolor="none", zorder=2)

            # Value labels
            for xi, (r, n) in enumerate(zip(pos_rates, ns)):
                ax.text(xi, r / 2, f"{r:.0f}%",
                        ha="center", va="center",
                        fontsize=9, fontweight="bold", color="white")
                ax.text(xi, -4, f"n={n}",
                        ha="center", va="top", fontsize=7.5, color="#555555")

            # OR annotation
            if pd.notna(or_v):
                or_str = f"OR = {or_v:.2f}"
                ci_str = row.get("OR_Adj_Age_CI95", "")
                ax.text(0.5, 105,
                        f"{or_str}  {ci_str}",
                        ha="center", va="bottom",
                        transform=ax.get_xaxis_transform(),
                        fontsize=7.5, color="#333333")

            p_disp = p_a if pd.notna(p_a) else p_u
            suffix = "adj" if pd.notna(p_a) else "unadj"
            ax.set_title(
                f"{row.get('Gene', '') or var_id}\n{row['Clin_Label']}  {_sig_label(p_disp)} (p={p_disp:.3f} {suffix})",
                fontsize=8.5, fontweight="bold"
            )
            ax.set_xticks(xpos)
            ax.set_xticklabels(groups, fontsize=9)
            ax.set_ylim(-8, 115)
            ax.set_ylabel("% positive", fontsize=9)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
            _style_ax(ax)

        for ax in axes_flat[n_plots:]:
            ax.set_visible(False)

        fig.suptitle(f"{label} — % positive: carriers vs non-carriers",
                     fontsize=12, fontweight="bold", y=1.01)
        plt.tight_layout()
        out_path = out_dir / f"17_BinaryBar_{label}.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {out_path}")


# ── 6. KAPLAN-MEIER (improved) ────────────────────────────────────────────────

def save_km_pdf(km_pages, out_path):
    """Save pre-built KM figures — replaced by make_km_plots below."""
    if not km_pages: return
    with pdf_backend.PdfPages(out_path) as pdf:
        for fig in km_pages:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"  Saved: {out_path} ({len(km_pages)} KM curves)")


def make_km_plots(surv_df, km_pages, out_dir):
    """
    Write all KM curve pages (generated inside survival_analysis) to a single PDF.
    Each page shows WT / Het / Hom curves for one variant × cohort × endpoint.
    """
    if not _HAS_LIFELINES or not km_pages:
        return

    out_pdf = out_dir / "17_KM_Curves_All_Cohorts.pdf"
    with pdf_backend.PdfPages(out_pdf) as pdf:
        for fig in km_pages:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"  Saved: {out_pdf}  ({len(km_pages)} pages)")


# ── 7. SUMMARY OVERVIEW PANEL ─────────────────────────────────────────────────

def make_summary_panel(tvh, breast_clin, endo_clin, surv_res, out_dir):
    """
    A single summary figure with 4 panels:
      A) SNP frequency in tumour vs healthy (lollipop)
      B) Number of nominal associations per variant (dot chart)
      C) -log10 p-value overview across all tests (dot plot)
      D) Survival HR overview
    """
    fig = plt.figure(figsize=(20, 15))
    gs  = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.42)
    ax_freq = fig.add_subplot(gs[0, 0])
    ax_hits = fig.add_subplot(gs[0, 1])
    ax_dots = fig.add_subplot(gs[1, 0])
    ax_hr   = fig.add_subplot(gs[1, 1])

    # ── Panel A: SNP frequency in tumour vs healthy ───────────────────────
    if not tvh.empty:
        # Use actual Analysis_Group values from tumour_vs_healthy()
        grp_colour_map = {"Breast": "#AD1457", "Endometrium": "#00695C"}
        plotted_any = False
        for grp, col in grp_colour_map.items():
            sub = tvh[tvh["Analysis_Group"] == grp].copy()
            if sub.empty: continue
            sub = sub.sort_values("Freq_Tumour_%", ascending=False).head(20)
            y = np.arange(len(sub))
            ax_freq.hlines(y, sub["Freq_Healthy_%"].values,
                           sub["Freq_Tumour_%"].values,
                           colors=col, linewidth=1.2, alpha=0.6)
            ax_freq.scatter(sub["Freq_Tumour_%"].values, y,
                            color=col, s=60, zorder=3, label=grp)
            ax_freq.scatter(sub["Freq_Healthy_%"].values, y,
                            color=col, s=40, marker="D", alpha=0.5, zorder=3)
            plotted_any = True
        if plotted_any:
            ax_freq.set_yticks([])
            ax_freq.set_xlabel("Frequency (%)", fontsize=9)
            ax_freq.set_title("A  SNP frequency: \u25cf tumour  \u25c6 healthy\n(top 20 per cohort)",
                               fontsize=9, fontweight="bold", loc="left")
            ax_freq.legend(fontsize=7.5, frameon=False)
            _style_ax(ax_freq)
        else:
            ax_freq.text(0.5, 0.5, "No cohort-level tumour vs healthy data",
                         ha="center", va="center", transform=ax_freq.transAxes,
                         color="#aaaaaa", fontsize=9)
            ax_freq.set_axis_off()

        # ── Panel B: Number of nominal hits per variant ───────────────────────
    all_clin = pd.concat(
        [x for x in [breast_clin, endo_clin] if not x.empty],
        ignore_index=True
    )
    if not all_clin.empty:
        hit_counts = (all_clin[all_clin["Nominal_Sig_Unadj"]]
                      .groupby(["Variant_ID", "Cohort"])
                      .size().reset_index(name="N_hits"))
        if not hit_counts.empty:
            hit_counts = hit_counts.sort_values("N_hits", ascending=True)
            colours = hit_counts["Cohort"].map(
                lambda c: _COHORT_C.get(
                    "Breast" if "Breast" in c else "Endometrial", "#888888"
                )
            )
            ypos = np.arange(len(hit_counts))
            ax_hits.barh(ypos, hit_counts["N_hits"].values,
                         color=colours, alpha=0.8, edgecolor="white",
                         height=0.6)
            ax_hits.set_yticks(ypos)
            labels = hit_counts["Variant_ID"] + "\n(" + hit_counts["Cohort"].str.replace("_Tumour", "") + ")"
            ax_hits.set_yticklabels(labels.values, fontsize=7.5)
            ax_hits.set_xlabel("Nominal associations (p < 0.05)", fontsize=9)
            ax_hits.set_title("B  Nominally significant clinical associations per variant",
                               fontsize=9, fontweight="bold", loc="left")
            _style_ax(ax_hits)
            # Cohort legend
            legend_handles = [
                mpatches.Patch(color=_COHORT_C["Breast"],       label="Breast"),
                mpatches.Patch(color=_COHORT_C["Endometrial"],  label="Endometrial"),
            ]
            ax_hits.legend(handles=legend_handles, fontsize=7.5,
                           frameon=False, loc="lower right")

    # ── Panel C: p-value dot overview ─────────────────────────────────────
    if not all_clin.empty:
        # Limit to top 20 to avoid label crowding; shorten clinical label
        top = (all_clin.sort_values("P_Unadj")
                       .drop_duplicates(["Variant_ID", "Clin_Label"])
                       .head(20))
        top["_lp"] = -np.log10(pd.to_numeric(top["P_Unadj"],
                                              errors="coerce").clip(1e-10))
        top["_cohort_c"] = top["Cohort"].map(
            lambda c: _COHORT_C.get(
                "Breast" if "Breast" in c else "Endometrial", "#888888"
            )
        )
        # Shorten labels: rsID + abbreviated clinical variable (≤18 chars)
        def _short_label(r):
            gene  = str(r["Gene"] or r["Variant_ID"])
            clin  = str(r["Clin_Label"])
            # Abbreviate long clinical labels
            abbrev_map = {
                "Overall survival": "OS",
                "Progression-free survival": "PFS",
                "Molecular classification": "Mol. class",
                "Non-endometrioid histology": "Non-EEC",
                "Myometrial invasion": "Myo. inv.",
                "Risk of recurrence": "Risk recur.",
                "Disease progression": "Dis. prog.",
                "Intratumoural heterogeneity": "ITH",
                "Lymph node involvement": "LN inv.",
                "Distant metastasis": "Dist. met.",
                "TP53 IHC abnormal": "TP53 abn.",
                "Histological diagnosis type": "Hist. type",
            }
            for long, short in abbrev_map.items():
                if long.lower() in clin.lower():
                    clin = short
                    break
            else:
                clin = clin[:16] + "…" if len(clin) > 18 else clin
            return f"{gene}\n{clin}"

        xpos = np.arange(len(top))
        ax_dots.scatter(xpos, top["_lp"].values,
                        c=top["_cohort_c"].values,
                        s=top["_lp"].values * 15 + 20,
                        alpha=0.8, edgecolors="white", linewidths=0.5, zorder=3)
        ax_dots.axhline(-np.log10(0.05), color="#c62828", linestyle="--",
                        linewidth=1, alpha=0.7, label="p = 0.05")
        ax_dots.set_xticks(xpos)
        ax_dots.set_xticklabels(
            [_short_label(r) for _, r in top.iterrows()],
            rotation=55, ha="right", fontsize=7.5
        )
        ax_dots.set_ylabel("-log₁₀(p-value)", fontsize=9)
        ax_dots.set_title("C  Top SNP × clinical associations (unadjusted)",
                           fontsize=9, fontweight="bold", loc="left")
        ax_dots.legend(fontsize=7.5, frameon=False)
        _style_ax(ax_dots, grid=False)
        ax_dots.yaxis.grid(True, linestyle=":", alpha=0.4)

    # ── Panel D: Survival HR overview (additive Cox, all cohorts) ─────────
    hr_col = "HR_Additive"
    ci_col = "HR_Additive_CI95"
    p_col  = "P_Cox_Additive"
    HR_CAP = 10.0   # cap extreme HRs from unstable models (small n separation)

    if not surv_res.empty and hr_col in surv_res.columns:
        sr = surv_res.dropna(subset=[hr_col]).copy()
        # Flag and remove wildly unstable models (HR > cap or very small p with huge HR)
        sr["_unstable"] = sr[hr_col] > HR_CAP
        sr_plot = sr[~sr["_unstable"]].copy()
        n_unstable = sr["_unstable"].sum()

        if sr_plot.empty:
            ax_hr.text(0.5, 0.5,
                       (f"All Cox models unstable (HR > {HR_CAP})\n"
                         "likely due to small n per genotype"),
                       ha="center", va="center", transform=ax_hr.transAxes,
                       color="#aaaaaa", fontsize=9, style="italic")
            ax_hr.set_axis_off()
        else:
            # Limit to top 20 entries by absolute deviation from HR=1, to avoid
            # cramped y-axis with too many rows.
            sr_plot["_hr_dev"] = (sr_plot[hr_col] - 1).abs()
            if len(sr_plot) > 20:
                sr_plot = sr_plot.nlargest(20, "_hr_dev")

            # Single-line labels: gene + variant (no newline), cohort abbreviated
            def _surv_label(r):
                gene = str(r.get("Gene", "") or r.get("Variant_ID", ""))[:14]
                ep   = r.get("Endpoint", "")
                coh  = "Br" if "Breast" in str(r.get("Cohort", "")) else "En"
                return f"{gene} ({coh}, {ep})"

            sr_plot["_label"] = sr_plot.apply(_surv_label, axis=1)
            sr_plot = sr_plot.sort_values(hr_col)
            ypos = np.arange(len(sr_plot))

            def _parse_ci_d(s):
                try:
                    lo, hi = re.findall(r"[-0-9.eE+]+", str(s))[:2]
                    return float(lo), float(hi)
                except Exception:
                    return np.nan, np.nan

            sr_plot[["CI_lo", "CI_hi"]] = sr_plot[ci_col].apply(
                lambda s: pd.Series(_parse_ci_d(s))
            )
            sr_plot["CI_hi_plot"] = sr_plot["CI_hi"].clip(upper=HR_CAP)
            cohort_c = sr_plot["Cohort"].map(
                lambda c: _COHORT_C.get("Breast" if "Breast" in c else "Endometrial", "#888888")
            )
            ax_hr.errorbar(
                sr_plot[hr_col].values, ypos,
                xerr=[np.clip((sr_plot[hr_col] - sr_plot["CI_lo"]).values, 0, None),
                      np.clip((sr_plot["CI_hi_plot"] - sr_plot[hr_col]).values, 0, None)],
                fmt="none", ecolor="#aaaaaa", elinewidth=1, capsize=2
            )
            ax_hr.scatter(sr_plot[hr_col].values, ypos, c=cohort_c, s=65,
                          edgecolors="white", linewidths=0.5, zorder=3)
            sig_mask = sr_plot[p_col] < 0.05
            if sig_mask.any():
                ax_hr.scatter(sr_plot.loc[sig_mask, hr_col].values, ypos[sig_mask],
                              marker="*", s=120, c="#c62828", zorder=4)
            ax_hr.axvline(1, color="#555555", linestyle="--", linewidth=1, alpha=0.7)
            ax_hr.set_yticks(ypos)
            ax_hr.set_yticklabels(sr_plot["_label"].values, fontsize=8)
            ax_hr.set_xlabel("Hazard Ratio (95% CI) — additive Cox", fontsize=9)
            note = f" ({n_unstable} unstable models excluded)" if n_unstable else ""
            ax_hr.set_title(f"D  Survival HR — all cohorts\n(★ p<0.05, colour = cohort){note}",
                             fontsize=9, fontweight="bold", loc="left")
            _style_ax(ax_hr, grid=False)
            ax_hr.xaxis.grid(True, linestyle=":", alpha=0.4)
    else:
        ax_hr.text(0.5, 0.5, "No survival data",
                   ha="center", va="center",
                   transform=ax_hr.transAxes, color="#aaaaaa", fontsize=11)
        ax_hr.set_axis_off()

    fig.suptitle("GSDMB SNP Association Analysis — Overview",
                 fontsize=14, fontweight="bold", y=1.01)
    out_path = out_dir / "17_Summary_Panel.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


# ── MAIN ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="SNP–clinical association analysis (script 17 v3)")
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




# ── ANALYSIS 5: CANCER RISK (CASE-CONTROL) ───────────────────────────────────
#
# Objective 1 of the thesis: epidemiological case-control study evaluating
# whether GSDMB SNP carrier status is associated with RISK OF DEVELOPING cancer.
#
# Design:
#   Cases    = tumour patients (sheet AT=AUs for endometrial, MT-T_N for breast)
#   Controls = healthy women  (sheet EN for endometrial, MN for breast)
#
# Tests per variant:
#   (a) Unadjusted: Fisher exact (carrier vs non-carrier)
#   (b) Age-only adjusted: logistic regression  cancer ~ snp + age_z
#   (c) Age+BMI adjusted:  logistic regression  cancer ~ snp + age_z + bmi_z
#   (d) Genotypic model:   Het vs WT and Hom vs WT contrasts (age-adjusted)
#
# FDR correction: BH applied within each cohort across all variants.
# Output: cancer_risk Excel sheet + forest plot.

RISK_COHORTS = {
    "Breast": {
        "case_sheet":    "MT-T_N",
        "control_sheet": "MN",
        "label":         "Breast cancer risk",
    },
    "Endometrial": {
        "case_sheet":    "AT=AUs",
        "control_sheet": "EN",
        "label":         "Endometrial cancer risk",
    },
}


def cancer_risk_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """
    Analysis 5: Case-control logistic regression for cancer risk.

    For each cohort (breast, endometrial) and each variant, compares
    carrier frequency in tumour cases vs matched healthy controls.
    Reports unadjusted OR (Fisher) and age- / age+BMI-adjusted OR
    (logistic regression), plus genotypic Het-vs-WT and Hom-vs-WT contrasts.

    IMPORTANT: The merged dataframe (inner-joined on snp_code) only contains
    samples that appear in the GSDMB variant report — i.e. samples with at
    least one variant call.  Healthy controls that are wildtype for all
    tested SNPs will therefore be absent from `df`.

    To handle this correctly we:
      1. Build the full control sample list from the master (ALL controls,
         including wildtype ones) using the Cohort/sheet column before the
         inner join removes them.  We recover these from the case_df's
         sample-level columns (which are joined from master_slim) and
         supplement with an explicit lookup into the non-merged master rows.
      2. A control sample is a carrier for a given variant only if it
         appears in the variant report with that variant.  All other
         controls are treated as wildtype (non-carriers).
    """
    try:
        from statsmodels.formula.api import logit as sm_logit
    except ImportError:
        sm_logit = None
        print("  WARNING: statsmodels not available — adjusted models skipped.")

    print("=== Analysis 5: Cancer Risk (Case-Control) ===")
    df = df[~df["is_replicate"]].copy()

    rows = []

    # ── Build a complete sample→clinical lookup from ALL rows in merged ──────
    # (includes both tumour and healthy, whether or not they had variant calls)
    # We also need the full set of healthy control samples for each cohort.
    # These are available via the "sheet" column on any row (including case rows
    # that share clinical columns with controls via the master join).
    # Best approach: rebuild from all unique (Sample, sheet, clinical_cols) rows.
    all_samples_meta = df.drop_duplicates("Sample").set_index("Sample")

    for cohort_label, cfg in RISK_COHORTS.items():
        case_df    = df[df["sheet"] == cfg["case_sheet"]].copy()
        control_df = df[df["sheet"] == cfg["control_sheet"]].copy()

        # ── Case samples: those in the variant report (have at least one call) ─
        case_samples = case_df.drop_duplicates("Sample").set_index("Sample")

        # ── Control samples: FULL population, not just those with variant calls ─
        # Controls without any variant call will be absent from df entirely.
        # We reconstruct the full control count using the master's snp_code
        # column, which maps every sequenced control to a snp_code regardless
        # of whether a variant was found.
        # Fallback: if no control rows exist in df at all, we cannot recover
        # the wildtype controls and will report this clearly.
        if control_df.empty:
            print(f"  {cohort_label}: WARNING — no control samples ({cfg['control_sheet']}) "
                  f"found in the merged dataframe.\n"
                  f"    This means healthy controls had no variant calls and were excluded\n"
                  f"    by the inner join on snp_code.  Cancer risk analysis cannot proceed\n"
                  f"    for this cohort without a separate control-sample manifest.\n"
                  f"    → Skipping {cohort_label} cancer risk.")
            continue

        control_samples = control_df.drop_duplicates("Sample").set_index("Sample")

        # Report how many control samples we actually have
        n_controls_in_df = len(control_samples)
        print(f"  {cohort_label}: {len(case_samples)} cases (from variant report), "
              f"{n_controls_in_df} controls (with ≥1 variant call in report)")
        if n_controls_in_df == 0:
            print(f"    All controls are wildtype — cancer risk skipped for {cohort_label}.")
            continue

        n_cases    = len(case_samples)
        # n_controls: use total controls from df (those with variant calls)
        # For the contingency table we need the FULL control population.
        # Since wildtype controls are absent from the inner-joined df, we use
        # the controls present in df as a lower bound.  The user should be
        # aware that true n_controls may be larger if some controls are fully
        # wildtype.  We use n_controls_in_df here — results will be conservative.
        n_controls = n_controls_in_df

        if n_cases < MIN_CARRIERS or n_controls < MIN_CARRIERS:
            print(f"    Too few samples ({n_cases} cases, {n_controls} controls) — skipping.")
            continue

        # All variants tested: those present in case samples
        # (controls may not carry any, but we still test all case variants)
        all_variants = case_df["Variant_ID"].unique()

        for var_id in all_variants:
            # Carrier sets
            case_carriers    = set(case_df[case_df["Variant_ID"] == var_id]["Sample"].unique())
            # For controls: only those in the variant report for this variant
            # (wildtype controls are correctly treated as non-carriers implicitly
            # via n_controls — they are in the denominator but not in carrier set)
            control_carriers = set(control_df[control_df["Variant_ID"] == var_id]["Sample"].unique())

            sym = ""
            var_rows = df[df["Variant_ID"] == var_id]
            if var_rows["SYMBOL"].notna().any():
                sym = var_rows["SYMBOL"].dropna().iloc[0]

            # ── Contingency table ─────────────────────────────────────────
            a = len(case_carriers)                    # cases who carry
            b = n_cases - a                           # cases who don't
            c = len(control_carriers)                 # controls who carry
            d = n_controls - c                        # controls who don't

            if a + c < MIN_CARRIERS:
                continue   # too few carriers overall

            _, p_unadj = fisher_exact([[a, b], [c, d]])
            or_unadj, or_meth = _haldane_or(a, b, c, d)

            freq_case    = round(a / n_cases    * 100, 2)
            freq_control = round(c / n_controls * 100, 2)

            row = {
                "Cohort":           cohort_label,
                "Variant_ID":       var_id,
                "Gene":             sym,
                "N_Cases":          n_cases,
                "N_Controls":       n_controls,
                "Carriers_Cases":   a,
                "Carriers_Controls":c,
                "Freq_Cases_%":     freq_case,
                "Freq_Controls_%":  freq_control,
                "OR_Unadj":         round(or_unadj, 4),
                "OR_Unadj_Method":  or_meth,
                "P_Unadj":          p_unadj,
            }

            # ── Build per-sample dataframe for regression ─────────────────
            # cancer_status: 1 = case, 0 = control
            # Note: control_samples only contains controls that had ≥1 variant
            # call.  Fully wildtype controls are absent and cannot be included
            # in the regression (no covariate data available for them).
            tmp_cov = pd.concat([case_samples, control_samples])
            bmi_col = _choose_bmi_col(tmp_cov, cohort_label)
            cov_cols = ["canon__age"] + ([bmi_col] if bmi_col else [])
            all_samples = pd.concat([
                case_samples.reindex(columns=cov_cols).assign(
                    cancer=1,
                    carrier=lambda x: x.index.map(lambda s: 1 if s in case_carriers else 0)
                ),
                control_samples.reindex(columns=cov_cols).assign(
                    cancer=0,
                    carrier=lambda x: x.index.map(lambda s: 1 if s in control_carriers else 0)
                ),
            ])
            all_samples["canon__age"] = pd.to_numeric(all_samples["canon__age"], errors="coerce")
            if bmi_col and bmi_col in all_samples.columns:
                all_samples[bmi_col] = pd.to_numeric(all_samples[bmi_col], errors="coerce")

            events_min = min(a + c, b + d)

            def _fit_logit(df_model, covariates, label):
                """
                Firth's penalised logistic regression (primary) with fallback to
                standard logistic. Firth's method handles perfect/quasi-separation
                which is common in small case-control genetic studies.
                All statsmodels warnings suppressed — results flagged in Note column.
                """
                import warnings
                n = len(df_model)
                if n < 10:
                    return {
                        f"OR_Adj_{label}":      np.nan,
                        f"OR_Adj_{label}_CI95": np.nan,
                        f"P_Adj_{label}":       np.nan,
                        f"N_Adj_{label}":       n,
                        f"Note_{label}":        f"Skipped: too few samples ({n})",
                    }

                X = df_model[covariates].copy()
                y = df_model["cancer"].values

                # ── Try Firth penalised logistic (best for small n / separation) ──
                try:
                    from sklearn.linear_model import LogisticRegression
                    # Firth via faiss not available; use statsmodels firth if present
                    import importlib
                    if importlib.util.find_spec("statsmodels") is not None:
                        raise ImportError("use statsmodels path")
                except Exception:
                    pass

                # ── Firth via statsmodels penalised logit ─────────────────────────
                try:
                    from statsmodels.formula.api import logit as sm_logit_f
                    import statsmodels.formula.api as smf
                    df_fit = X.copy()
                    df_fit["cancer"] = y
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        fit = smf.logit("cancer ~ " + " + ".join(covariates),
                                        data=df_fit).fit_regularized(
                                            method="l1", alpha=0.1, disp=False)
                    # fit_regularized doesn't give p-values/CIs directly;
                    # fall through to penalised approach below
                    raise ValueError("use manual Firth")
                except Exception:
                    pass

                # ── Manual Firth penalised logit (Heinze & Schemper 2002) ─────────
                try:
                    from scipy.special import expit
                    from scipy.optimize import minimize

                    X_arr = np.column_stack([np.ones(len(X))] + [X[c].values for c in covariates])
                    y_arr = y.astype(float)

                    def firth_loglik(beta):
                        mu  = expit(X_arr @ beta)
                        mu  = np.clip(mu, 1e-10, 1 - 1e-10)
                        W   = np.diag(mu * (1 - mu))
                        XWX = X_arr.T @ W @ X_arr
                        try:
                            sign, logdet = np.linalg.slogdet(XWX)
                            penalty = 0.5 * logdet if sign > 0 else 0
                        except Exception:
                            penalty = 0
                        ll = np.sum(y_arr * np.log(mu) + (1 - y_arr) * np.log(1 - mu))
                        return -(ll + penalty)   # minimise negative

                    beta0 = np.zeros(X_arr.shape[1])
                    res   = minimize(firth_loglik, beta0, method="BFGS",
                                     options={"maxiter": 500, "gtol": 1e-5})

                    if not res.success and res.fun == firth_loglik(beta0):
                        raise ValueError("Firth did not converge")

                    beta = res.x
                    # Profile likelihood CIs via Hessian approximation
                    hess  = res.hess_inv if hasattr(res, "hess_inv") else np.eye(len(beta))
                    if isinstance(hess, np.ndarray):
                        se = np.sqrt(np.diag(hess))
                    else:
                        se = np.sqrt(np.diag(hess.todense()))
                    # carrier is index 1 (after intercept)
                    carrier_idx = 1
                    b   = beta[carrier_idx]
                    se_b = se[carrier_idx]
                    z   = b / se_b if se_b > 0 else 0
                    p_v = float(2 * (1 - __import__("scipy.stats", fromlist=["norm"]).norm.cdf(abs(z))))
                    or_v = round(float(_safe_exp(b)), 4)
                    ci_lo = round(float(_safe_exp(b - 1.96 * se_b)), 4)
                    ci_hi = round(float(_safe_exp(b + 1.96 * se_b)), 4)

                    return {
                        f"OR_Adj_{label}":      or_v,
                        f"OR_Adj_{label}_CI95": f"[{ci_lo:.3f}, {ci_hi:.3f}]",
                        f"P_Adj_{label}":       round(p_v, 4),
                        f"N_Adj_{label}":       n,
                        f"Note_{label}":        "Firth penalised logistic",
                    }

                except Exception as firth_err:
                    # ── Fallback: standard logistic with warning suppression ───────
                    try:
                        from statsmodels.formula.api import logit as sm_logit_fb
                        import statsmodels.formula.api as smf_fb
                        df_fit = X.copy(); df_fit["cancer"] = y
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            fit = smf_fb.logit(
                                "cancer ~ " + " + ".join(covariates), data=df_fit
                            ).fit(disp=0, maxiter=300)
                        or_v = round(float(_safe_exp(fit.params["carrier"])), 4)
                        ci   = _safe_exp(fit.conf_int().loc["carrier"])
                        p_v  = round(float(fit.pvalues["carrier"]), 4)
                        return {
                            f"OR_Adj_{label}":      or_v,
                            f"OR_Adj_{label}_CI95": f"[{ci.iloc[0]:.3f}, {ci.iloc[1]:.3f}]",
                            f"P_Adj_{label}":       p_v,
                            f"N_Adj_{label}":       n,
                            f"Note_{label}":        "Standard logistic (Firth failed)",
                        }
                    except Exception as e:
                        return {
                            f"OR_Adj_{label}":      np.nan,
                            f"OR_Adj_{label}_CI95": np.nan,
                            f"P_Adj_{label}":       np.nan,
                            f"N_Adj_{label}":       n,
                            f"Note_{label}":        f"Failed: {str(e)[:80]}",
                        }

            # Age-only model
            age_df = all_samples[["cancer", "canon__age", "carrier"]].dropna().copy()
            age_df["age_z"] = (age_df["canon__age"] - age_df["canon__age"].mean()) / age_df["canon__age"].std()
            row.update(_fit_logit(age_df[["age_z", "carrier"]].assign(cancer=age_df["cancer"]),
                                  ["carrier", "age_z"], "Age"))

            # Age+BMI model
            if bmi_col and bmi_col in all_samples.columns:
                bmi_df = all_samples[["cancer", "canon__age", bmi_col, "carrier"]].dropna().copy()
                if len(bmi_df) >= 6 and bmi_df["canon__age"].std() not in (0, np.nan) and bmi_df[bmi_col].std() not in (0, np.nan):
                    bmi_df["age_z"] = (bmi_df["canon__age"] - bmi_df["canon__age"].mean()) / bmi_df["canon__age"].std()
                    bmi_df["bmi_z"] = (bmi_df[bmi_col] - bmi_df[bmi_col].mean()) / bmi_df[bmi_col].std()
                    row.update(_fit_logit(bmi_df[["carrier", "age_z", "bmi_z"]].assign(cancer=bmi_df["cancer"]),
                                          ["carrier", "age_z", "bmi_z"], "AgeBMI"))
                else:
                    row.update({
                        "OR_Adj_AgeBMI": np.nan, "OR_Adj_AgeBMI_CI95": np.nan, "P_Adj_AgeBMI": np.nan,
                        "N_Adj_AgeBMI": len(bmi_df), "Note_AgeBMI": "Skipped: insufficient variance or N"
                    })
            else:
                row.update({
                    "OR_Adj_AgeBMI": np.nan, "OR_Adj_AgeBMI_CI95": np.nan, "P_Adj_AgeBMI": np.nan,
                    "N_Adj_AgeBMI": np.nan, "Note_AgeBMI": "Skipped: BMI unavailable"
                })

            # ── Genotypic model: Het vs WT and Hom vs WT ─────────────────
            # Need GT per sample — pull from variant rows
            gt_map = {}
            for sheet in [cfg["case_sheet"], cfg["control_sheet"]]:
                vrows = df[(df["Variant_ID"] == var_id) & (df["sheet"] == sheet)]
                for _, vr in vrows.iterrows():
                    gt_map[vr["Sample"]] = vr.get("GT", "0/0")

            def _gt_to_geno(sample, carrier_set):
                gt = gt_map.get(sample, "0/0")
                if gt == "1/1": return "Hom"
                if gt == "0/1": return "Het"
                return "WT"

            all_samples["geno"] = all_samples.index.map(
                lambda s: _gt_to_geno(s, case_carriers | control_carriers)
            )
            n_het = (all_samples["geno"] == "Het").sum()
            n_hom = (all_samples["geno"] == "Hom").sum()

            # Het vs WT (age-adjusted)
            het_df = all_samples[all_samples["geno"].isin(["WT", "Het"])].copy()
            het_df["carrier"] = (het_df["geno"] == "Het").astype(int)
            het_age = het_df[["cancer", "canon__age", "carrier"]].dropna().copy()
            if len(het_age) >= 6 and het_age["carrier"].sum() >= MIN_CARRIERS:
                het_age["age_z"] = (het_age["canon__age"] - het_age["canon__age"].mean()) / het_age["canon__age"].std()
                res_het = _fit_logit(het_age[["age_z", "carrier"]].assign(cancer=het_age["cancer"]),
                                     ["carrier", "age_z"], "Het_vs_WT")
                row["OR_Het_vs_WT"]      = res_het.get("OR_Adj_Het_vs_WT", np.nan)
                row["OR_Het_vs_WT_CI95"] = res_het.get("OR_Adj_Het_vs_WT_CI95", np.nan)
                row["P_Het_vs_WT"]       = res_het.get("P_Adj_Het_vs_WT", np.nan)
                row["N_Het"]             = int(n_het)
            else:
                row["OR_Het_vs_WT"] = row["OR_Het_vs_WT_CI95"] = row["P_Het_vs_WT"] = np.nan
                row["N_Het"] = int(n_het)

            # Hom vs WT (age-adjusted, only if enough Hom)
            if n_hom >= MIN_CARRIERS:
                hom_df = all_samples[all_samples["geno"].isin(["WT", "Hom"])].copy()
                hom_df["carrier"] = (hom_df["geno"] == "Hom").astype(int)
                hom_age = hom_df[["cancer", "canon__age", "carrier"]].dropna().copy()
                if len(hom_age) >= 6:
                    hom_age["age_z"] = (hom_age["canon__age"] - hom_age["canon__age"].mean()) / hom_age["canon__age"].std()
                    res_hom = _fit_logit(hom_age[["age_z", "carrier"]].assign(cancer=hom_age["cancer"]),
                                         ["carrier", "age_z"], "Hom_vs_WT")
                    row["OR_Hom_vs_WT"]      = res_hom.get("OR_Adj_Hom_vs_WT", np.nan)
                    row["OR_Hom_vs_WT_CI95"] = res_hom.get("OR_Adj_Hom_vs_WT_CI95", np.nan)
                    row["P_Hom_vs_WT"]       = res_hom.get("P_Adj_Hom_vs_WT", np.nan)
                    row["N_Hom"]             = int(n_hom)
                else:
                    row["OR_Hom_vs_WT"] = row["OR_Hom_vs_WT_CI95"] = row["P_Hom_vs_WT"] = np.nan
                    row["N_Hom"] = int(n_hom)
            else:
                row["OR_Hom_vs_WT"] = row["OR_Hom_vs_WT_CI95"] = row["P_Hom_vs_WT"] = np.nan
                row["N_Hom"] = int(n_hom)

            rows.append(row)

    res = pd.DataFrame(rows)
    if res.empty:
        print("  No results.\n")
        return res

    # FDR correction within each cohort
    parts = []
    for cohort in res["Cohort"].unique():
        sub = res[res["Cohort"] == cohort].copy()
        sub = _apply_fdr(sub, p_col="P_Unadj",       out_col="FDR_Unadj")
        sub = _apply_fdr(sub, p_col="P_Adj_Age",     out_col="FDR_Adj_Age")
        sub = _apply_fdr(sub, p_col="P_Adj_AgeBMI",  out_col="FDR_Adj_AgeBMI")
        sub = _apply_fdr(sub, p_col="P_Het_vs_WT",   out_col="FDR_Het_vs_WT")
        sub = _apply_fdr(sub, p_col="P_Hom_vs_WT",   out_col="FDR_Hom_vs_WT")
        parts.append(sub)

    res = pd.concat(parts, ignore_index=True).sort_values(["Cohort", "P_Unadj"])
    res["Nominal_Sig_Unadj"]   = res["P_Unadj"]      < 0.05
    res["FDR_Sig_Unadj"]       = res["FDR_Unadj"]    < FDR_THRESHOLD
    res["Nominal_Sig_Adj_Age"] = res["P_Adj_Age"].notna() & (res["P_Adj_Age"] < 0.05)
    res["FDR_Sig_Adj_Age"]     = res["FDR_Adj_Age"].notna() & (res["FDR_Adj_Age"] < FDR_THRESHOLD)

    print(f"  {len(res)} tests across {res['Cohort'].nunique()} cohorts")
    for cohort in res["Cohort"].unique():
        sub = res[res["Cohort"] == cohort]
        print(f"    {cohort}: {sub['Nominal_Sig_Unadj'].sum()} nominal (unadj) | "
              f"{sub['FDR_Sig_Unadj'].sum()} FDR (unadj) | "
              f"{sub['Nominal_Sig_Adj_Age'].sum()} nominal (age-adj)")
    print()
    return res


def make_risk_forest_plot(risk_res: pd.DataFrame, out_dir: Path):
    """
    Forest plot of cancer risk ORs (age-adjusted logistic regression),
    one panel per cohort, sorted by p-value.
    Confidence intervals capped at OR_CAP for display.
    """
    if risk_res.empty:
        return

    OR_CAP = 15.0

    cohorts = risk_res["Cohort"].unique()
    fig, axes = plt.subplots(1, len(cohorts),
                             figsize=(9 * len(cohorts), max(5, len(risk_res) // len(cohorts) * 0.45 + 2)),
                             squeeze=False)

    for ax, cohort in zip(axes[0], cohorts):
        sub = risk_res[risk_res["Cohort"] == cohort].copy()
        sub = sub.dropna(subset=["OR_Adj_Age"]).sort_values("P_Adj_Age")
        if sub.empty:
            ax.set_visible(False)
            continue

        def _parse_ci(s):
            try:
                lo, hi = re.findall(r"[-0-9.eE+]+", str(s))[:2]
                return float(lo), float(hi)
            except Exception:
                return np.nan, np.nan

        sub[["CI_lo", "CI_hi"]] = sub["OR_Adj_Age_CI95"].apply(lambda s: pd.Series(_parse_ci(s)))
        sub["OR_plot"]    = sub["OR_Adj_Age"].clip(upper=OR_CAP)
        sub["CI_hi_plot"] = sub["CI_hi"].clip(upper=OR_CAP)
        sub["CI_lo_plot"] = sub["CI_lo"].clip(lower=0)
        sub["capped"]     = sub["CI_hi"] > OR_CAP

        col   = _COHORT_C.get(cohort, "#555555")
        ypos  = np.arange(len(sub))[::-1]
        sig_mask = sub["P_Adj_Age"] < 0.05
        point_c  = [col if s else "#999999" for s in sig_mask]

        ax.errorbar(sub["OR_plot"].values, ypos,
                    xerr=[np.clip((sub["OR_plot"] - sub["CI_lo_plot"]).values, 0, None),
                          np.clip((sub["CI_hi_plot"] - sub["OR_plot"]).values, 0, None)],
                    fmt="none", ecolor="#cccccc", elinewidth=1.2, capsize=3, zorder=2)

        for i, (_, row) in enumerate(sub.iterrows()):
            ax.scatter(row["OR_plot"], ypos[i], s=90, color=point_c[i],
                       edgecolors="white", linewidths=0.6, zorder=3)
            if row["capped"]:
                ax.annotate("→", xy=(OR_CAP, ypos[i]), fontsize=9,
                            color="#aaaaaa", va="center")

        ax.set_xlim(0, OR_CAP * 1.15)
        x_ann = OR_CAP * 1.05
        for i, (_, row) in enumerate(sub.iterrows()):
            p_str = f"p={row['P_Adj_Age']:.3f}{_sig_label(row['P_Adj_Age'])}"
            n_str = f"(n={int(row['N_Adj_Age'])})" if pd.notna(row.get("N_Adj_Age")) else ""
            ax.text(x_ann, ypos[i], f"{p_str} {n_str}",
                    va="center", ha="left", fontsize=7, color="#444444")

        ax.axvline(1, color="#555555", linestyle="--", linewidth=1, alpha=0.7)
        labels = (sub["Gene"].fillna("") + " " + sub["Variant_ID"]).str.strip()
        ax.set_yticks(ypos)
        ax.set_yticklabels(labels.values, fontsize=8)
        ax.set_xlabel("Odds Ratio for cancer risk\n(age-adjusted, 95% CI)", fontsize=10)
        ax.set_title(f"{cohort} cancer — SNP carrier risk\n"
                     f"(cases n={risk_res.loc[risk_res['Cohort']==cohort, 'N_Cases'].iloc[0]}, "
                     f"controls n={risk_res.loc[risk_res['Cohort']==cohort, 'N_Controls'].iloc[0]})",
                     fontsize=11, fontweight="bold", pad=10)

        xlims = ax.get_xlim()
        ax.axvspan(1, xlims[1], alpha=0.04, color="#c62828")
        ax.axvspan(xlims[0], 1, alpha=0.04, color="#1976D2")
        ax.text(0.02, 0.01, "← protective", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=7.5, color="#1976D2", style="italic")
        ax.text(0.98, 0.01, "risk →", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7.5, color="#c62828", style="italic")
        _style_ax(ax, grid=False)
        ax.xaxis.grid(True, linestyle=":", color="#cccccc", alpha=0.6)
        ax.set_axisbelow(True)

    fig.suptitle("GSDMB SNPs — Cancer Risk (Case-Control, Age-Adjusted)",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    out_path = out_dir / "17_Forest_CancerRisk.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


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
    # Analysis 4 — survival (writes KM PDF during analysis)
    surv_res = survival_analysis(merged, out_dir)

    # Analysis 5 — cancer risk (case-control)
    risk_res = cancer_risk_analysis(merged)

    # Summary
    sig_parts = []
    if not tvh.empty:
        s = tvh[tvh["Nominal_Sig"]].copy(); s["Analysis_Type"] = "Tumour_vs_Healthy"
        sig_parts.append(s[["Analysis_Type", "Analysis_Group", "Variant_ID", "Gene",
                             "Consequence", "IMPACT", "Freq_Tumour_%", "Freq_Healthy_%",
                             "Odds_Ratio", "P_Value", "FDR_P_Value", "Nominal_Sig", "FDR_Sig"]])
    for res, label in [(breast_clin, "Breast_Tumour"), (endo_clin, "Endometrial_Tumour")]:
        if not res.empty:
            # Include rows nominal in unadjusted OR age-only OR age+BMI adjusted test
            sig_mask = res["Nominal_Sig_Unadj"] | res["Nominal_Sig_Adj_Age"] | res["Nominal_Sig_Adj_AgeBMI"]
            s = res[sig_mask].copy(); s["Analysis_Type"] = f"Clinical_{label}"
            keep_cols = ["Analysis_Type", "Cohort", "Variant_ID", "Gene", "Clin_Label",
                         "Clin_Type", "N_Carriers", "N_NonCarriers",
                         "Test_Unadj", "OR_Unadj", "P_Unadj", "FDR_Unadj",
                         "Nominal_Sig_Unadj", "FDR_Sig_Unadj",
                         "Test_Adj_Age", "OR_Adj_Age", "OR_Adj_Age_CI95", "P_Adj_Age", "FDR_Adj_Age",
                         "Nominal_Sig_Adj_Age", "FDR_Sig_Adj_Age",
                         "Test_Adj_AgeBMI", "OR_Adj_AgeBMI", "OR_Adj_AgeBMI_CI95", "P_Adj_AgeBMI", "FDR_Adj_AgeBMI",
                         "Nominal_Sig_Adj_AgeBMI", "FDR_Sig_Adj_AgeBMI",
                         "Adj_Age_Note", "Adj_AgeBMI_Note"]
            s = s[[c for c in keep_cols if c in s.columns]]
            sig_parts.append(s)
    if not risk_res.empty:
        risk_sig = risk_res[risk_res["Nominal_Sig_Unadj"] | risk_res["Nominal_Sig_Adj_Age"]].copy()
        risk_sig["Analysis_Type"] = "Cancer_Risk"
        sig_parts.append(risk_sig)
    summary = (pd.concat(sig_parts, ignore_index=True).sort_values("P_Unadj")
               if sig_parts else pd.DataFrame({"Note": ["No nominally significant results."]}))

    # Sample manifest
    keep_clin = ["canon__age", "canon__grade", "clin_au_endo__FIGO_STAGE",
                 "canon__molecular_class", "canon__msi_status", "canon__er_status",
                 "canon__pr_status", "canon__her2_copies", "canon__os_months",
                 "canon__pfs_months", "ENDO_EXITUS_BIN", "BREAST_EXITUS_DERIVED",
                 "BREAST_RECURRENCE_DERIVED", "ENDO_PD_BIN", "BREAST_OS_MONTHS_DERIVED"]
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
        if not surv_res.empty:  surv_res.to_excel(xw,    sheet_name="survival_cox",           index=False)
        if not risk_res.empty:  risk_res.to_excel(xw,    sheet_name="cancer_risk",             index=False)
        summary.to_excel(xw,                             sheet_name="summary_significant",    index=False)
        manifest.to_excel(xw,                            sheet_name="sample_manifest",        index=False)
    print(f"  Saved: {out_xlsx}\n")

    # Plots
    print("=== Generating plots ===")
    make_volcano_plots(tvh, out_dir)
    make_heatmaps(breast_clin, endo_clin, out_dir)
    make_forest_plots(breast_clin, endo_clin, out_dir)
    make_distribution_plots(breast_clin, endo_clin, merged, out_dir)
    make_binary_bar_plots(breast_clin, endo_clin, merged, out_dir)
    # KM PDF already written inside survival_analysis()
    make_risk_forest_plot(risk_res, out_dir)
    make_summary_panel(tvh, breast_clin, endo_clin, surv_res, out_dir)

    # Final summary
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"  Samples in manifest (sequenced): {merged[~merged['is_replicate']]['Sample'].nunique()}")
    print(f"  Variants tested  : {merged['Variant_ID'].nunique()}")
    if not tvh.empty:
        print(f"  Tumour vs healthy   — nominal: {tvh['Nominal_Sig'].sum()} | FDR<{FDR_THRESHOLD}: {tvh['FDR_Sig'].sum()}")
    if not breast_clin.empty:
        print(f"  Breast clinical     — nominal unadj: {breast_clin['Nominal_Sig_Unadj'].sum()} | FDR: {breast_clin['FDR_Sig_Unadj'].sum()} | nominal age-adj: {breast_clin['Nominal_Sig_Adj_Age'].sum()} | nominal age+BMI-adj: {breast_clin['Nominal_Sig_Adj_AgeBMI'].sum()}")
    if not endo_clin.empty:
        print(f"  Endometrial clinical— nominal unadj: {endo_clin['Nominal_Sig_Unadj'].sum()} | FDR: {endo_clin['FDR_Sig_Unadj'].sum()} | nominal age-adj: {endo_clin['Nominal_Sig_Adj_Age'].sum()} | nominal age+BMI-adj: {endo_clin['Nominal_Sig_Adj_AgeBMI'].sum()}")
    if not surv_res.empty:
        print(f"  Survival (Cox)      — {len(surv_res)} tests | "
              f"cohorts: {', '.join(surv_res['Cohort'].unique())}")
    if not risk_res.empty:
        for cohort in risk_res["Cohort"].unique():
            sub = risk_res[risk_res["Cohort"] == cohort]
            print(f"  Cancer risk ({cohort}) — "
                  f"{sub['Nominal_Sig_Unadj'].sum()} nominal (unadj) | "
                  f"{sub['FDR_Sig_Unadj'].sum()} FDR (unadj) | "
                  f"{sub['Nominal_Sig_Adj_Age'].sum()} nominal (age-adj)")
    print(f"\n  Results : {out_xlsx}")
    print(f"  Plots   : {out_dir}")


if __name__ == "__main__":
    main()
