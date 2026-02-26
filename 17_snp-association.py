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

    cell_h  = 0.52
    cell_w  = 1.4
    fig_h   = max(5, len(pivot) * cell_h + 2.5)
    fig_w   = max(8, len(pivot.columns) * cell_w + 3)

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
                    ha="center", va="center", fontsize=6.5, color=txt_c)
            if star not in ("", "ns"):
                ax.text(j + 0.5, i + 0.72, star,
                        ha="center", va="center",
                        fontsize=8, color="white" if lp >= sig_t else "#c62828",
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
    plt.xticks(rotation=40, ha="right", fontsize=8)
    plt.yticks(fontsize=7.5)
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
        adj_res = clin_res[clin_res["P_Adj"].notna()].copy()
        if not adj_res.empty:
            _draw_heatmap(adj_res, "P_Adj",
                          f"Age-adjusted p — {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_raw_p.png",
                          filter_thresh=0.10)
            _draw_heatmap(adj_res, "FDR_Adj",
                          f"Age-adjusted FDR — {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_FDR.png",
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

        # Keep only binary tests with a valid adjusted OR and CI
        sub = clin_res[
            (clin_res["Clin_Type"] == "binary") &
            clin_res["OR_Adj"].notna() &
            clin_res["OR_Adj_CI95"].notna()
        ].copy()
        if sub.empty: continue

        # Parse CI strings "[lo, hi]"
        def _parse_ci(s):
            try:
                lo, hi = re.findall(r"[-0-9.eE+]+", str(s))[:2]
                return float(lo), float(hi)
            except Exception:
                return np.nan, np.nan
        sub[["CI_lo", "CI_hi"]] = sub["OR_Adj_CI95"].apply(
            lambda s: pd.Series(_parse_ci(s))
        )
        sub = sub.dropna(subset=["CI_lo", "CI_hi"])
        if sub.empty: continue

        # Sort by p-value; label = "Gene: clinical var"
        sub = sub.sort_values("P_Adj")
        sub["label"] = sub.apply(
            lambda r: f"{r['Gene'] or r['Variant_ID']} — {r['Clin_Label']}", axis=1
        )

        n_rows = len(sub)
        fig_h  = max(5, n_rows * 0.42 + 2)
        fig, ax = plt.subplots(figsize=(9, fig_h))

        ypos = np.arange(n_rows)[::-1]   # top = most significant

        # Colour by significance
        sig_mask = sub["P_Adj"] < 0.05
        point_c  = [colour if s else "#999999" for s in sig_mask]

        # Error bars
        xerr_lo = (sub["OR_Adj"] - sub["CI_lo"]).values
        xerr_hi = (sub["CI_hi"] - sub["OR_Adj"]).values
        ax.errorbar(sub["OR_Adj"].values, ypos,
                    xerr=[np.clip(xerr_lo, 0, None), np.clip(xerr_hi, 0, None)],
                    fmt="none", ecolor="#aaaaaa", elinewidth=1.2, capsize=3, zorder=2)

        # Points
        for i, (_, row) in enumerate(sub.iterrows()):
            ax.scatter(row["OR_Adj"], ypos[i],
                       s=90, color=point_c[i],
                       edgecolors="white", linewidths=0.6, zorder=3)

        # p-value and OR labels on the right
        x_max = ax.get_xlim()[1]
        for i, (_, row) in enumerate(sub.iterrows()):
            pstr = f"p={row['P_Adj']:.3f}{_sig_label(row['P_Adj'])}"
            ax.text(x_max * 1.01, ypos[i], pstr,
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
            p_a     = row.get("P_Adj", np.nan)

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
                                 squeeze=False)
        axes_flat = axes.flatten()

        cohort_df = merged[
            (merged["sheet"] == coh_sheet) & (~merged["is_replicate"])
        ].drop_duplicates("Sample").copy()

        for idx, (_, row) in enumerate(bin_hits.iterrows()):
            ax     = axes_flat[idx]
            col    = row["Clinical_Var"]
            var_id = row["Variant_ID"]
            p_u    = row["P_Unadj"]
            p_a    = row.get("P_Adj", np.nan)
            or_v   = row.get("OR_Adj", row.get("OR_Unadj", np.nan))

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
                ci_str = row.get("OR_Adj_CI95", "")
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


def make_km_plots(surv_df, merged, out_dir):
    """
    Improved KM plots: at-risk table, shaded CI, clean styling.
    One page per variant × endpoint combination.
    """
    if not _HAS_LIFELINES or surv_df.empty: return

    endpoint_meta = {
        "OS":  ("canon__os_months",  "ENDO_EXITUS_BIN",  "Overall Survival"),
        "PFS": ("canon__pfs_months", "ENDO_PD_BIN",      "Progression-free Survival"),
    }

    endo_df = (merged[
        (~merged["is_replicate"]) &
        (merged["sheet"] == "AT=AUs")
    ].drop_duplicates("Sample")
     .set_index("Sample"))

    out_pdf = out_dir / "17_KM_Curves_AU_Endo.pdf"
    with pdf_backend.PdfPages(out_pdf) as pdf:
        for _, srow in surv_df.iterrows():
            ep      = srow["Endpoint"]
            var_id  = srow["Variant_ID"]
            gene    = srow.get("Gene", "") or var_id
            t_col, ev_col, ep_label = endpoint_meta[ep]

            if t_col not in endo_df.columns or ev_col not in endo_df.columns:
                continue

            # Carrier flag
            var_samples = set(
                merged[(merged["Variant_ID"] == var_id) &
                        (merged["sheet"] == "AT=AUs")]["Sample"].unique()
            )
            T = endo_df[t_col]
            E = endo_df[ev_col].astype(float)
            x = endo_df.index.map(lambda s: "Carrier" if s in var_samples
                                  else "Non-carrier")
            valid = pd.concat([T, E, x.rename("group")], axis=1).dropna()
            valid.columns = ["T", "E", "group"]
            if valid["group"].nunique() < 2: continue

            # Build figure with at-risk table
            fig = plt.figure(figsize=(8, 6))
            gs  = fig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0.05)
            ax_km   = fig.add_subplot(gs[0])
            ax_risk = fig.add_subplot(gs[1], sharex=ax_km)

            kmf = KaplanMeierFitter()
            time_points = np.linspace(0, valid["T"].max(), 8)
            risk_table  = {}

            for grp, c in [("Carrier", _PALETTE["carrier"]),
                            ("Non-carrier", _PALETTE["non_carrier"])]:
                m = valid["group"] == grp
                if m.sum() < 2: continue
                kmf.fit(valid.loc[m, "T"], valid.loc[m, "E"],
                        label=f"{grp} (n={m.sum()})")
                kmf.plot_survival_function(
                    ax=ax_km, ci_show=True,
                    color=c, ci_alpha=0.12, linewidth=2.2,
                )
                # At-risk counts at each time point
                risk_table[grp] = [
                    int((valid.loc[m, "T"] >= t).sum()) for t in time_points
                ]

            # Significance + HR annotation
            lr_p = srow.get("Logrank_P", np.nan)
            hr   = srow.get("HR", np.nan)
            hr_ci= srow.get("HR_CI95", "")
            info = (f"Log-rank p = {lr_p:.4f}" if pd.notna(lr_p) else "")
            if pd.notna(hr):
                info += f"\nHR = {hr:.2f} {hr_ci}"
            ax_km.text(0.97, 0.97, info,
                       transform=ax_km.transAxes,
                       ha="right", va="top", fontsize=8.5,
                       bbox=dict(boxstyle="round,pad=0.4", fc="white",
                                 ec="#cccccc", alpha=0.9))

            ax_km.set_ylabel("Survival probability", fontsize=10)
            ax_km.set_xlabel("")
            ax_km.set_ylim(-0.03, 1.08)
            ax_km.set_title(
                f"{gene} ({var_id})  —  {ep_label}",
                fontsize=11, fontweight="bold"
            )
            _style_ax(ax_km, grid=False)
            ax_km.yaxis.grid(True, linestyle=":", alpha=0.4)
            ax_km.axhline(0.5, color="#aaaaaa", linestyle=":", linewidth=1)

            # At-risk table
            ax_risk.set_xlim(ax_km.get_xlim())
            ax_risk.set_ylim(-0.5, len(risk_table) - 0.5)
            ax_risk.set_yticks(range(len(risk_table)))
            ax_risk.set_yticklabels(list(risk_table.keys()), fontsize=8)
            for yi, (grp, counts) in enumerate(risk_table.items()):
                for xi, (tp, cnt) in enumerate(zip(time_points, counts)):
                    ax_risk.text(tp, yi, str(cnt),
                                 ha="center", va="center", fontsize=7.5,
                                 color=(_PALETTE["carrier"] if grp == "Carrier"
                                        else _PALETTE["non_carrier"]))
            ax_risk.set_xlabel("Time (months)", fontsize=10)
            ax_risk.set_ylabel("At risk", fontsize=8, labelpad=4)
            ax_risk.spines["top"].set_visible(False)
            ax_risk.spines["right"].set_visible(False)
            ax_risk.spines["left"].set_visible(False)
            ax_risk.yaxis.set_tick_params(length=0)
            ax_risk.xaxis.grid(False)
            ax_risk.yaxis.grid(False)
            plt.setp(ax_km.get_xticklabels(), visible=False)

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"  Saved: {out_pdf}")


# ── 7. SUMMARY OVERVIEW PANEL ─────────────────────────────────────────────────

def make_summary_panel(tvh, breast_clin, endo_clin, surv_res, out_dir):
    """
    A single summary figure with 4 panels:
      A) SNP frequency in tumour vs healthy (lollipop)
      B) Number of nominal associations per variant (dot chart)
      C) -log10 p-value overview across all tests (dot plot)
      D) Survival HR overview
    """
    fig = plt.figure(figsize=(18, 13))
    gs  = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.4)
    ax_freq = fig.add_subplot(gs[0, 0])
    ax_hits = fig.add_subplot(gs[0, 1])
    ax_dots = fig.add_subplot(gs[1, 0])
    ax_hr   = fig.add_subplot(gs[1, 1])

    # ── Panel A: SNP frequency in tumour vs healthy ───────────────────────
    if not tvh.empty:
        grp_cols = {"Breast_Tumour_vs_Healthy": "#AD1457",
                    "Endo_Tumour_vs_Healthy":   "#00695C"}
        for grp, col in grp_cols.items():
            sub = tvh[tvh["Analysis_Group"] == grp].copy()
            if sub.empty: continue
            sub = sub.sort_values("Freq_Tumour_%", ascending=False).head(20)
            y = np.arange(len(sub))
            ax_freq.hlines(y, sub["Freq_Healthy_%"].values,
                           sub["Freq_Tumour_%"].values,
                           colors=col, linewidth=1.2, alpha=0.6)
            ax_freq.scatter(sub["Freq_Tumour_%"].values, y,
                            color=col, s=60, zorder=3,
                            label=grp.replace("_", " "))
            ax_freq.scatter(sub["Freq_Healthy_%"].values, y,
                            color=col, s=40, marker="D",
                            alpha=0.5, zorder=3)
        ax_freq.set_yticks([])
        ax_freq.set_xlabel("Frequency (%)", fontsize=9)
        ax_freq.set_title("A  SNP frequency: \u25cf tumour  \u25c6 healthy\n(top 20 per cohort)",
                           fontsize=9, fontweight="bold", loc="left")
        ax_freq.legend(fontsize=7.5, frameon=False)
        _style_ax(ax_freq)

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
        top = (all_clin.sort_values("P_Unadj")
                       .drop_duplicates(["Variant_ID", "Clin_Label"])
                       .head(30))
        top["_lp"] = -np.log10(pd.to_numeric(top["P_Unadj"],
                                              errors="coerce").clip(1e-10))
        top["_cohort_c"] = top["Cohort"].map(
            lambda c: _COHORT_C.get(
                "Breast" if "Breast" in c else "Endometrial", "#888888"
            )
        )
        xpos = np.arange(len(top))
        ax_dots.scatter(xpos, top["_lp"].values,
                        c=top["_cohort_c"].values,
                        s=top["_lp"].values * 15 + 20,
                        alpha=0.8, edgecolors="white", linewidths=0.5, zorder=3)
        ax_dots.axhline(-np.log10(0.05), color="#c62828", linestyle="--",
                        linewidth=1, alpha=0.7, label="p = 0.05")
        ax_dots.set_xticks(xpos)
        ax_dots.set_xticklabels(
            [f"{r['Gene'] or r['Variant_ID']}\n{r['Clin_Label']}" for _, r in top.iterrows()],
            rotation=45, ha="right", fontsize=6.5
        )
        ax_dots.set_ylabel("-log₁₀(p-value)", fontsize=9)
        ax_dots.set_title("C  Top SNP × clinical associations (unadjusted)",
                           fontsize=9, fontweight="bold", loc="left")
        ax_dots.legend(fontsize=7.5, frameon=False)
        _style_ax(ax_dots, grid=False)
        ax_dots.yaxis.grid(True, linestyle=":", alpha=0.4)

    # ── Panel D: Survival HR overview ─────────────────────────────────────
    if not surv_res.empty and "HR" in surv_res.columns:
        sr = surv_res.dropna(subset=["HR"]).copy()
        sr["_label"] = sr["Gene"].fillna("") + " " + sr["Variant_ID"] + "\n(" + sr["Endpoint"] + ")"
        sr = sr.sort_values("HR")
        ypos = np.arange(len(sr))

        def _parse_ci(s):
            try:
                lo, hi = re.findall(r"[-0-9.eE+]+", str(s))[:2]
                return float(lo), float(hi)
            except Exception:
                return np.nan, np.nan

        sr[["CI_lo", "CI_hi"]] = sr["HR_CI95"].apply(
            lambda s: pd.Series(_parse_ci(s))
        )
        ax_hr.errorbar(
            sr["HR"].values, ypos,
            xerr=[np.clip((sr["HR"] - sr["CI_lo"]).values, 0, None),
                  np.clip((sr["CI_hi"] - sr["HR"]).values, 0, None)],
            fmt="none", ecolor="#aaaaaa", elinewidth=1, capsize=2
        )
        c_hr = sr["Logrank_P"].apply(
            lambda p: "#c62828" if pd.notna(p) and p < 0.05 else "#999999"
        )
        ax_hr.scatter(sr["HR"].values, ypos, c=c_hr, s=65,
                      edgecolors="white", linewidths=0.5, zorder=3)
        ax_hr.axvline(1, color="#555555", linestyle="--", linewidth=1, alpha=0.7)
        ax_hr.set_yticks(ypos)
        ax_hr.set_yticklabels(sr["_label"].values, fontsize=7.5)
        ax_hr.set_xlabel("Hazard Ratio (95% CI)", fontsize=9)
        ax_hr.set_title("D  Survival HR \u2014 AU endometrial cohort\n(\u25cf p<0.05  \u25cb n.s.)",
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
    make_forest_plots(breast_clin, endo_clin, out_dir)
    make_distribution_plots(breast_clin, endo_clin, merged, out_dir)
    make_binary_bar_plots(breast_clin, endo_clin, merged, out_dir)
    make_km_plots(surv_res, merged, out_dir)
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
        print(f"  Breast clinical     — nominal unadj: {breast_clin['Nominal_Sig_Unadj'].sum()} | FDR: {breast_clin['FDR_Sig_Unadj'].sum()} | nominal adj: {breast_clin['Nominal_Sig_Adj'].sum()}")
    if not endo_clin.empty:
        print(f"  Endometrial clinical— nominal unadj: {endo_clin['Nominal_Sig_Unadj'].sum()} | FDR: {endo_clin['FDR_Sig_Unadj'].sum()} | nominal adj: {endo_clin['Nominal_Sig_Adj'].sum()}")
    if not surv_res.empty:
        print(f"  Survival (AU endo)  — {len(surv_res)} tests")
    print(f"\n  Results : {out_xlsx}")
    print(f"  Plots   : {out_dir}")


if __name__ == "__main__":
    main()
