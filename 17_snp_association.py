#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
17_snp_association.py  (v4 — SNP-only, genotypic model, all-cohort Cox, age+BMI adjusted)
=======================================================================
SNP–Clinical Variable Association Analysis

Joins the GSDMB variant report (GSDMB_Annotated_Report.xlsx) with the
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
1. Tumour vs Control       — Fisher's exact test per cohort + globally
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
  ? Death             clin_dcs__Exitus             (No/Yes)
  ✓ HER2+ subtype     clin_her2__DX                (Ca mama HER2+ / TN)
  ~ Survival (derived) from Fecha_dx + Ultima_fecha + death-status source fields

ENDOMETRIAL (AT=AUs sheet, n=114 sequenced DNA samples)
  ✓ FIGO stage        clin_au_endo__FIGO_STAGE     (IA/IB/II/IIIA/IIIC1/IIIC2/IVB)
  ✓ Grade             clin_au_endo__GRADE          (G1/G2/G3)
  ✓ Histology group   clin_au_endo__HISTOLOGY_GROUP (NEEC/EEC)
  ✓ LVSI              clin_au_endo__LVSI            (YES/NO)
  ✓ Myometrial inv.   clin_au_endo__MYOMETRIAL_INFILTRATION (<50%/>50%)
  ✓ MSI status        canon__msi_status            (Stable/Unstable)
  ✓ Molecular class   canon__molecular_class       (P53/NSMP/MMRd/POLE)
  ✓ TP53 IHC (abnml) canon__p53_status            (Aberrant/Normal — produced by 16b)
  ✓ ER status         canon__er_status             (Positive/Negative)
  ✓ PR status         canon__pr_status             (Positive/Negative)
  ✓ Progression       clin_au_endo__PD_STATUS      (PD / No PD)
  ? Death             clin_au_endo__EXITUS         (Yes/No)
  ✓ Risk recurrence   clin_au_endo__RISK_OF_RECURRENCE
  ✓ OS (months)       canon__os_months
  ✓ PFS (months)      canon__pfs_months
  ✓ Age               canon__age

OUTPUTS
-------
  17_SNP_Clinical_Association_Results.xlsx
    • tumour_vs_control
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
python3 17_snp_association.py \\
  --gsdmb   /path/to/GSDMB_Annotated_Report.xlsx \\
  --master  /path/to/MASTER_SNP_plus_clinical_HARMONISED.xlsx \\
  --out_dir /path/to/output/
"""

from __future__ import annotations

import argparse
import os
import re
import warnings
import textwrap
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

from association_runtime import script17_defaults
from figure_style import COMPARATIVE_TAG, COHORT_COLORS, GENOTYPE_COLORS, IMPACT_COLORS, arm_color, cohort_color, tagged_title
from pipeline_utils import attach_amplicon_warning_columns, build_amplicon_warning_lookup
from pipeline_validation import print_validation_summary, validate_file_exists, validate_percentage_columns

# lifelines optional — only needed for survival plots
try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    _HAS_LIFELINES = True
except ImportError:
    _HAS_LIFELINES = False
    print("WARNING: lifelines not installed — survival plots will be skipped.")

# Suppress only known noisy third-party warnings; convergence warnings from
# statsmodels and lifelines are intentionally left visible so that failed or
# unreliable model fits are not silently swallowed.
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=".*tight_layout.*")
warnings.filterwarnings("ignore", message=".*More than 20 figures.*")

# ── DEFAULT PATHS ─────────────────────────────────────────────────────────────
DEFAULTS = script17_defaults()
DEFAULT_GSDMB  = DEFAULTS["gsdmb"]
DEFAULT_MASTER = DEFAULTS["master"]
DEFAULT_OUT    = DEFAULTS["out_dir"]
DEFAULT_VARIANT_WHITELIST = DEFAULTS.get("variant_whitelist")

# Pass-BAM manifests - one file per cohort, listing QC-passed BAMs.
# These are the ground truth for which samples have actually been sequenced.
# Set to None to skip manifest filtering (falls back to extraction_flag logic).
DEFAULT_MANIFESTS = DEFAULTS.get("manifests", {
    "endometrium-tumour": Path("/home/gadeaalonsoj/tfm/manifests/endometrium-tumour-pass_manifest.txt"),
    "endometrium-normal": Path("/home/gadeaalonsoj/tfm/manifests/endometrium-normal-pass_manifest.txt"),
    "breast-tumour":      Path("/home/gadeaalonsoj/tfm/manifests/breast-tumour-pass_manifest.txt"),
    "breast-normal":      Path("/home/gadeaalonsoj/tfm/manifests/breast-normal-pass_manifest.txt"),
})

FDR_THRESHOLD  = float(DEFAULTS.get("fdr_threshold", 0.10))
MIN_CARRIERS   = int(DEFAULTS.get("min_carriers", 5))   # EPV-informed minimum; n=3 gives unreliable estimates
MIN_NFE_AF     = float(DEFAULTS.get("min_nfe_af", 0.01))

MIN_COMPARISON_CARRIERS = 3  # minimum cases/controls in each cell for stable comparison
MIN_EVENTS_FOR_LOGISTIC = 10  # EPV rule: ~10 events per predictor for stable logistic


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
      SNP_AT_2__BL_039__IonCode_0118     ? SNP_AT_2
      DNA_SNP_EN_10_IonCode_0101_...     ? SNP_EN_10
      DNA_MT-T_37_IonCode_0121           ? SNP_MT-T_37
    """
    s = str(sample_name)
    m = re.match(r"^(SNP_(?:AT|EN|MN|MT-T)_\d+)", s, re.IGNORECASE)
    if m: return m.group(1).upper()
    m = re.match(r"^DNA_SNP_(EN|MN)_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_{m.group(1).upper()}_{m.group(2)}"
    m = re.match(r"^DNA_SNP_AT_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_AT_{m.group(1)}"
    m = re.match(r"^DNA_AT_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_AT_{m.group(1)}"
    m = re.match(r"^SNP_DNA_AT_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_AT_{m.group(1)}"
    m = re.match(r"^SNP_DNA_(MN|EN)_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_{m.group(1).upper()}_{m.group(2)}"
    m = re.match(r"^DNA_SNP_MT-T_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_MT-T_{m.group(1)}"
    m = re.match(r"^DNA_MT[-_]T_(\d+)_", s, re.IGNORECASE)
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
    """
    Coerce messy KI67 to a fraction in [0, 1].

    Handles inputs such as:
      ? Integer percentage:  40      ? 0.40
      ? String percentage:   '40%'   ? 0.40
      ? Decimal fraction:    0.40    ? 0.40  (already a fraction, left as-is)
      ? Range (percent):     '20-25%'? midpoint 22.5 ? 0.225
      ? Range (fraction):    '0.2-0.25' ? midpoint 0.225 (already a fraction)

    The heuristic: if the raw numeric value (or midpoint) is > 1 it is assumed
    to be a percentage and is divided by 100.  Values ? 1 are treated as
    already being a fraction and are returned unchanged.
    """
    if pd.isna(val):
        return None
    s = str(val).replace("%", "").strip()
    m = re.match(r"([\d.]+)\s*[-?]\s*([\d.]+)", s)
    if m:
        mid = (float(m.group(1)) + float(m.group(2))) / 2
        return mid / 100 if mid > 1 else mid
    try:
        v = float(s)
        return v / 100 if v > 1 else v
    except (ValueError, TypeError):
        return None


def _derive_breast_os(row) -> Optional[float]:
    """Fallback OS in months from diagnosis to last available date."""
    try:
        dx = pd.to_datetime(row.get("clin_dcs__Fecha dx"), errors="coerce")
        ult = pd.to_datetime(
            row.get("clin_dcs__?ltima fecha disponible",
                    row.get("clin_dcs__?ltima fecha disponible")),
            errors="coerce",
        )
        if pd.isna(dx) or pd.isna(ult):
            return None
        return max((ult - dx).days / 30.44, 0)
    except Exception:
        return None



def _ki67_fraction_from_pct(val) -> Optional[float]:
    """Convert a 0-100 KI67 percentage into a 0-1 fraction."""
    if pd.isna(val):
        return None
    try:
        v = float(val)
    except (TypeError, ValueError):
        return None
    if 0 <= v <= 100:
        return round(v / 100.0, 4)
    return None


def _extract_variant_id(row: pd.Series) -> str:
    """
    Prefer rsIDs where available; otherwise fall back to a genomic label.
    This keeps downstream plots readable while remaining deterministic for
    novel or non-rs variants.
    """
    existing = str(row.get("Existing_variation", "") or "")
    for token in re.split(r"[,&;\s]+", existing):
        token = token.strip()
        if re.fullmatch(r"rs\d+", token, flags=re.IGNORECASE):
            return token

    chrom = str(row.get("CHROM", "") or "").strip()
    pos = pd.to_numeric(row.get("POS"), errors="coerce")
    ref = str(row.get("REF", "") or "").strip()
    alt = str(row.get("ALT", "") or "").strip()
    if chrom and pd.notna(pos) and ref and alt:
        return f"{chrom}:{int(pos)}_{ref}>{alt}"
    if chrom and pd.notna(pos):
        return f"{chrom}:{int(pos)}"
    return existing.split("&")[0].strip() if existing else "UNKNOWN_VARIANT"


def _load_variant_whitelist_filters(whitelist_path: Path) -> Dict[str, set]:
    """
    Read a flexible Excel whitelist and return identifiers that can be used to
    filter the annotated variant table. The whitelist can contain rsIDs, a
    Variant_ID column, or genomic coordinates (CHROM/POS/REF/ALT).
    """
    wl = pd.read_excel(whitelist_path)
    wl.columns = [str(c).strip() for c in wl.columns]

    ids: set = set()
    rsids: set = set()
    coord_keys: set = set()
    positions: set = set()

    for col in wl.columns:
        series = wl[col].dropna().astype(str).str.strip()
        ids.update(v for v in series if v)
        rsids.update(v for v in series if re.fullmatch(r"rs\d+", v, flags=re.IGNORECASE))

    pos_col = next((c for c in wl.columns if re.fullmatch(r"POS|Position", c, flags=re.IGNORECASE)), None)
    ref_col = next((c for c in wl.columns if c.upper() == "REF"), None)
    alt_col = next((c for c in wl.columns if c.upper() == "ALT"), None)
    chrom_col = next((c for c in wl.columns if c.upper() == "CHROM"), None)

    if pos_col is not None:
        pos_vals = pd.to_numeric(wl[pos_col], errors="coerce").dropna().astype(int)
        positions.update(pos_vals.tolist())

    if all(c is not None for c in [chrom_col, pos_col, ref_col, alt_col]):
        tmp = wl[[chrom_col, pos_col, ref_col, alt_col]].copy()
        tmp[pos_col] = pd.to_numeric(tmp[pos_col], errors="coerce")
        tmp = tmp.dropna(subset=[pos_col, ref_col, alt_col])
        tmp[chrom_col] = tmp[chrom_col].astype(str).str.strip()
        tmp[ref_col] = tmp[ref_col].astype(str).str.strip().str.upper()
        tmp[alt_col] = tmp[alt_col].astype(str).str.strip().str.upper()
        coord_keys.update(
            f"{r[chrom_col]}|{int(r[pos_col])}|{r[ref_col]}|{r[alt_col]}"
            for _, r in tmp.iterrows()
        )

    return {
        "ids": {v for v in ids if v},
        "rsids": {v.lower() for v in rsids if v},
        "coord_keys": coord_keys,
        "positions": positions,
    }


CLINICAL_VARS_BREAST: Dict[str, Dict] = {
    "BREAST_GRADE_NUMERIC":      {"type": "continuous", "label": "Grade", "note": "Ordinal 1-3 breast grade"},
    "BREAST_ER_BIN":             {"type": "binary",     "label": "ER status", "note": "Positive vs negative"},
    "BREAST_PR_BIN":             {"type": "binary",     "label": "PR status", "note": "Positive vs negative"},
    "BREAST_RECURRENCE_DERIVED": {"type": "binary",     "label": "Recurrence / progression", "note": "Derived from clinical recurrence/progression status"},
    "BREAST_METASTASIS_DERIVED": {"type": "binary",     "label": "Distant metastasis", "note": "Derived from metastasis status"},
    "canon__her2_copies":        {"type": "continuous", "label": "HER2 copies", "note": "Continuous HER2 copy number"},
    "canon__age":                {"type": "continuous", "label": "Age", "note": "Age at diagnosis or surgery"},
    "BREAST_KI67_NUMERIC":       {"type": "continuous", "label": "KI67", "note": "Fractional KI67 value"},
    "BREAST_EXITUS_DERIVED":     {"type": "binary",     "label": "Death", "note": "Derived overall death status"},
    "BREAST_HER2_SUBTYPE":       {"type": "nominal",    "label": "HER2 subtype", "note": "HER2+, TN, or Other"},
    "BREAST_OS_MONTHS_DERIVED":  {"type": "continuous", "label": "Overall survival (months)", "note": "Derived from canon or follow-up dates"},
}

CLINICAL_VARS_ENDO: Dict[str, Dict] = {
    "ENDO_FIGO_NUMERIC":         {"type": "continuous", "label": "FIGO stage", "note": "Ordered FIGO mapping"},
    "ENDO_GRADE_NUMERIC":        {"type": "continuous", "label": "Grade", "note": "Ordinal 1-3 endometrial grade"},
    "ENDO_NEEC_BIN":             {"type": "binary",     "label": "Histology group (NEEC)", "note": "NEEC = 1 vs EEC = 0"},
    "ENDO_LVSI_BIN":             {"type": "binary",     "label": "LVSI", "note": "Present vs absent"},
    "ENDO_MYOINV_BIN":           {"type": "binary",     "label": "Myometrial invasion >50%", "note": ">50% = 1 vs <50% = 0"},
    "ENDO_MSI_BIN":              {"type": "binary",     "label": "MSI status", "note": "Unstable = 1 vs stable = 0"},
    "canon__molecular_class":    {"type": "nominal",    "label": "Molecular class", "note": "Harmonised endometrial molecular class"},
    "ENDO_TP53_ABN_BIN":         {"type": "binary",     "label": "TP53 IHC aberrant", "note": "Aberrant = 1 vs normal = 0"},
    "ENDO_ER_BIN":               {"type": "binary",     "label": "ER status", "note": "Positive vs negative"},
    "ENDO_PR_BIN":               {"type": "binary",     "label": "PR status", "note": "Positive vs negative"},
    "ENDO_PD_BIN":               {"type": "binary",     "label": "Progression", "note": "Progressive disease yes = 1 vs no = 0"},
    "ENDO_EXITUS_BIN":           {"type": "binary",     "label": "Death", "note": "Overall death status"},
    "ENDO_EXITUS_DISEASE_BIN":   {"type": "binary",     "label": "Disease-specific death", "note": "Death due to disease"},
    "ENDO_RISK_ORDINAL":         {"type": "continuous", "label": "Risk of recurrence", "note": "Ordered recurrence-risk grouping"},
    "canon__os_months":          {"type": "continuous", "label": "Overall survival (months)", "note": "Harmonised OS months"},
    "canon__pfs_months":         {"type": "continuous", "label": "Progression-free survival (months)", "note": "Harmonised PFS months"},
    "canon__age":                {"type": "continuous", "label": "Age", "note": "Age at diagnosis or surgery"},
    "ENDO_M_STAGE_BIN":          {"type": "binary",     "label": "Metastatic stage (M1)", "note": "M1 = 1 vs M0 = 0"},
}


def load_clinical_master(master_path: Path,
                         manifest_paths: Optional[Dict[str, Path]] = None) -> pd.DataFrame:
    """
    Load the harmonised master, restrict to DNA samples, optionally apply the
    pass-BAM manifest filter, and derive the clinical columns required by the
    downstream association analyses.
    """
    print("  Loading harmonised clinical master ?")
    try:
        master = pd.read_excel(master_path, sheet_name="harmonised_plus_canon")
    except Exception:
        print("  WARNING: 'harmonised_plus_canon' not found ? trying sheet 0")
        master = pd.read_excel(master_path, sheet_name=0)

    master.columns = [str(c).strip() for c in master.columns]

    if "nucleic_acid" in master.columns:
        n_before = len(master)
        master = master[master["nucleic_acid"].astype(str).str.upper() == "DNA"].copy()
        print(f"  DNA rows retained: {len(master)} / {n_before}")

    if "snp_code" not in master.columns:
        for alt in ["SNP_code", "Sample", "SAMPLE", "sample_id"]:
            if alt in master.columns:
                master = master.rename(columns={alt: "snp_code"})
                break
    if "snp_code" not in master.columns:
        raise ValueError("Cannot find snp_code column in harmonised master workbook")

    master["snp_code"] = (
        master["snp_code"]
        .astype(str)
        .str.strip()
        .str.replace(r"\s+", "", regex=True)
    )

    replicate_codes: set = set()
    if manifest_paths:
        sequenced_codes, replicate_codes = parse_manifests(manifest_paths)
        before = master["snp_code"].nunique()
        master = master[master["snp_code"].isin(sequenced_codes)].copy()
        after = master["snp_code"].nunique()
        print(f"  Manifest filter retained {after} / {before} DNA samples")
    else:
        if "extraction_flag" in master.columns and "sheet" in master.columns:
            excl = (
                (master["sheet"] == "AT=AUs")
                & (
                    master["extraction_flag"].astype(str).str.upper().eq("NO HAY")
                    | master.get("pd_status", pd.Series("", index=master.index)).astype(str).str.upper().eq("NO HACER")
                )
            )
            dropped = int(excl.sum())
            if dropped:
                print(f"  Fallback filter removed {dropped} AU rows flagged NO HAY / NO HACER")
            master = master[~excl].copy()

    if "is_replicate" in master.columns:
        master["is_replicate"] = master["is_replicate"].fillna(False).astype(bool)
    else:
        master["is_replicate"] = False
    if replicate_codes:
        master["is_replicate"] = master["is_replicate"] | master["snp_code"].isin(replicate_codes)

    sheet_to_tissue = {
        "MT-T_N": "Tumour",
        "AT=AUs": "Tumour",
        "MN": "Healthy",
        "EN": "Healthy",
    }
    sheet_to_cohort = {
        "MT-T_N": "Breast_Tumour",
        "MN": "Breast_Healthy",
        "AT=AUs": "Endometrial_Tumour",
        "EN": "Endometrial_Healthy",
    }
    if "sheet" in master.columns:
        master["Tissue"] = master["sheet"].map(sheet_to_tissue).fillna(master.get("tissue"))
        master["Cohort"] = master["sheet"].map(sheet_to_cohort).fillna(master["sheet"])
    else:
        tissue_raw = master.get("tissue", pd.Series("", index=master.index)).astype(str).str.strip()
        master["Tissue"] = tissue_raw.replace({"Normal": "Healthy", "Tumour": "Tumour", "Tumor": "Tumour"})
        master["Cohort"] = master.get("cohort", "Unknown")

    if "canon__bmi" not in master.columns:
        master["canon__bmi"] = np.nan

    if "canon__her2_copies" not in master.columns:
        master["canon__her2_copies"] = pd.to_numeric(
            master.get("clin_dcs__COPIAS HER2", master.get("clin_her2__COPIAS HER2", pd.Series(np.nan, index=master.index))),
            errors="coerce",
        )

    for src, dst, fn in [
        ("clin_dcs__Recaida/Progresi?n", "BREAST_RECURRENCE_DERIVED", _derive_breast_recurrence),
        ("clin_dcs__Exitus",             "BREAST_EXITUS_DERIVED",     _derive_breast_exitus),
        ("clin_dcs__MTxDISTANCIA",       "BREAST_METASTASIS_DERIVED", _derive_breast_metastasis),
        ("clin_her2__DX",                "BREAST_HER2_SUBTYPE",       _derive_her2_subtype),
    ]:
        if src in master.columns:
            master[dst] = master[src].apply(fn)

    if "canon__er_status" in master.columns:
        master["BREAST_ER_BIN"] = master["canon__er_status"].map({"Positive": 1, "Negative": 0})
        master["ENDO_ER_BIN"] = master["canon__er_status"].map({"Positive": 1, "Negative": 0})
    elif "clin_dcs__RE" in master.columns:
        master["BREAST_ER_BIN"] = (
            master["clin_dcs__RE"].astype(str).str.strip().str.upper().map({"POSITIVO": 1, "NEGATIVO": 0})
        )

    if "canon__pr_status" in master.columns:
        master["BREAST_PR_BIN"] = master["canon__pr_status"].map({"Positive": 1, "Negative": 0})
        master["ENDO_PR_BIN"] = master["canon__pr_status"].map({"Positive": 1, "Negative": 0})
    elif "clin_dcs__RP" in master.columns:
        master["BREAST_PR_BIN"] = (
            master["clin_dcs__RP"].astype(str).str.strip().str.upper().map({"POSITIVO": 1, "NEGATIVO": 0})
        )

    ki67_pct = master.get("canon__ki67_pct", pd.Series(np.nan, index=master.index))
    grade = master.get("canon__grade", pd.Series(np.nan, index=master.index))
    figo = master.get("canon__figo_stage", pd.Series(np.nan, index=master.index))
    risk = master.get("canon__risk_group", pd.Series(np.nan, index=master.index))
    lvsi = master.get("canon__lvsi", pd.Series(pd.NA, index=master.index))
    myoinv = master.get("canon__myometrial_invasion", pd.Series(pd.NA, index=master.index))
    pd_flag = master.get("canon__pd_flag", pd.Series(pd.NA, index=master.index))
    exitus_flag = master.get("canon__exitus_flag", pd.Series(pd.NA, index=master.index))
    breast_os_raw = master.get("canon__os_months", pd.Series(np.nan, index=master.index))

    master["BREAST_KI67_NUMERIC"] = ki67_pct.apply(_ki67_fraction_from_pct)
    master["BREAST_GRADE_NUMERIC"] = pd.to_numeric(grade, errors="coerce")
    master["BREAST_OS_MONTHS_DERIVED"] = pd.to_numeric(breast_os_raw, errors="coerce").combine_first(
        master.apply(_derive_breast_os, axis=1)
    )

    master["ENDO_FIGO_NUMERIC"] = figo.apply(_derive_figo_numeric)
    master["ENDO_GRADE_NUMERIC"] = pd.to_numeric(grade, errors="coerce")
    master["ENDO_RISK_ORDINAL"] = risk.apply(_derive_risk_ordinal)
    master["ENDO_LVSI_BIN"] = lvsi.map({"Yes": 1, "No": 0, "YES": 1, "NO": 0})
    master["ENDO_MYOINV_BIN"] = myoinv.map({"<50%": 0, ">50%": 1})
    master["ENDO_MSI_BIN"] = master.get("canon__msi_status", pd.Series(pd.NA, index=master.index)).map({"Unstable": 1, "Stable": 0})
    master["ENDO_NEEC_BIN"] = master.get("clin_au_endo__HISTOLOGY_GROUP", pd.Series(pd.NA, index=master.index)).map({"NEEC": 1, "EEC": 0})
    master["ENDO_PD_BIN"] = pd.to_numeric(pd_flag, errors="coerce")
    master["ENDO_EXITUS_BIN"] = pd.to_numeric(exitus_flag, errors="coerce")
    master["ENDO_EXITUS_DISEASE_BIN"] = (
        master.get("clin_au_endo__EXITUS_DISEASE", pd.Series(pd.NA, index=master.index))
        .astype(str).str.strip().str.upper().map({"YES": 1, "NO": 0})
    )
    master["ENDO_TP53_ABN_BIN"] = master.get("canon__p53_status", pd.Series(pd.NA, index=master.index)).map({"Aberrant": 1, "Normal": 0})
    if "clin_au_endo__M" in master.columns:
        master["ENDO_M_STAGE_BIN"] = master["clin_au_endo__M"].apply(
            lambda x: 0 if pd.notna(x) and str(x).strip().upper() == "M0" else (1 if pd.notna(x) else None)
        )

    master = master.drop_duplicates("snp_code").copy()
    print(f"  Clinical master ready: {len(master)} unique DNA samples")
    if "sheet" in master.columns:
        print(f"  Sheet breakdown: {master['sheet'].value_counts().to_dict()}")
    return master


def load_and_merge(gsdmb_path: Path,
                   master_path: Path,
                   manifest_paths: Optional[Dict[str, Path]] = None,
                   whitelist_path: Optional[Path] = None) -> pd.DataFrame:
    """
    Join the annotated GSDMB callset to the harmonised clinical master.

    The merge is intentionally a sample-level inner join on `snp_code`: only
    samples present in the clinical master and represented in the annotated
    callset are retained, matching the documented script-17 behaviour.
    """
    print("=== Loading GSDMB SNP calls and harmonised clinical data ===")
    master = load_clinical_master(master_path, manifest_paths=manifest_paths)

    print("  Loading annotated GSDMB variants ?")
    variants = pd.read_excel(gsdmb_path, sheet_name="Biological_Annotations")
    variants.columns = [str(c).strip() for c in variants.columns]

    if "Sample" not in variants.columns:
        raise ValueError("Annotated GSDMB workbook is missing the Sample column")

    variants["Sample"] = variants["Sample"].astype(str).str.strip()
    variants["snp_code"] = variants["Sample"].apply(_extract_snp_code)
    failed_rows = int(variants["snp_code"].isna().sum())
    if failed_rows:
        examples = variants.loc[variants["snp_code"].isna(), "Sample"].drop_duplicates().head(8).tolist()
        print(f"  WARNING: could not derive snp_code for {failed_rows} rows; examples: {examples}")
    variants = variants.dropna(subset=["snp_code"]).copy()

    if "GT" in variants.columns:
        variants["GT"] = (
            variants["GT"].astype(str).str.strip()
            .replace({"0|0": "0/0", "0|1": "0/1", "1|0": "1/0", "1|1": "1/1"})
        )

    if "Cohort" in variants.columns:
        variants = variants.rename(columns={"Cohort": "Cohort_callset"})
    if "Tissue" in variants.columns:
        variants = variants.rename(columns={"Tissue": "Tissue_callset"})

    variants["Variant_ID"] = variants.apply(_extract_variant_id, axis=1)
    if all(c in variants.columns for c in ["CHROM", "POS", "REF", "ALT"]):
        variants["_coord_key"] = (
            variants["CHROM"].astype(str).str.strip() + "|"
            + pd.to_numeric(variants["POS"], errors="coerce").fillna(-1).astype(int).astype(str) + "|"
            + variants["REF"].astype(str).str.strip().str.upper() + "|"
            + variants["ALT"].astype(str).str.strip().str.upper()
        )
    else:
        variants["_coord_key"] = ""

    if whitelist_path is not None:
        filters = _load_variant_whitelist_filters(whitelist_path)
        rsid_mask = variants["Variant_ID"].astype(str).str.lower().isin(filters["rsids"])
        id_mask = variants["Variant_ID"].astype(str).isin(filters["ids"])
        coord_mask = variants["_coord_key"].isin(filters["coord_keys"])
        pos_mask = pd.to_numeric(variants.get("POS"), errors="coerce").fillna(-1).astype(int).isin(filters["positions"])
        keep_mask = id_mask | rsid_mask | coord_mask | pos_mask
        before = variants["Variant_ID"].nunique()
        variants = variants[keep_mask].copy()
        after = variants["Variant_ID"].nunique()
        print(f"  Variant whitelist retained {after} / {before} unique variants")

    dedup_cols = [c for c in ["Sample", "Variant_ID", "CHROM", "POS", "REF", "ALT"] if c in variants.columns]
    if dedup_cols:
        variants = variants.drop_duplicates(subset=dedup_cols).copy()

    merged = variants.merge(master, on="snp_code", how="inner")

    print(
        f"  Annotated samples: {variants['snp_code'].nunique()} | "
        f"clinical samples: {master['snp_code'].nunique()} | "
        f"matched samples: {merged['snp_code'].nunique()}"
    )
    print(f"  Matched rows: {len(merged)} | unique variants: {merged['Variant_ID'].nunique()}")
    if not merged.empty:
        print(f"  Tissue breakdown: {merged.drop_duplicates('Sample')['Tissue'].value_counts().to_dict()}")
        print(f"  Cohort breakdown: {merged.drop_duplicates('Sample')['Cohort'].value_counts().to_dict()}")

    return merged.drop(columns=["_coord_key"], errors="ignore")


def tumour_vs_control(merged: pd.DataFrame) -> pd.DataFrame:
    """
    Analysis 1: compare tumour vs pooled healthy controls for each variant.

    Cohort-specific comparisons use the relevant tumour arm (Breast or
    Endometrium) against the pooled healthy control arm; a Global comparison
    uses all tumours vs all healthy controls.
    """
    print("=== Analysis 1: Tumour vs Control ===")
    df = merged[~merged["is_replicate"]].copy()
    if df.empty:
        print("  No merged rows - skipping.\n")
        return pd.DataFrame()

    sample_manifest = df.drop_duplicates("Sample").copy()
    pooled_controls = sample_manifest[sample_manifest["Tissue"] == "Healthy"].copy()

    comparison_defs = {
        "Breast": {
            "mask": (sample_manifest["Tissue"] == "Tumour")
                    & sample_manifest["Cohort"].astype(str).str.contains("Breast", case=False, na=False),
        },
        "Endometrium": {
            "mask": (sample_manifest["Tissue"] == "Tumour")
                    & sample_manifest["Cohort"].astype(str).str.contains("Endometri", case=False, na=False),
        },
        "Global": {
            "mask": sample_manifest["Tissue"] == "Tumour",
        },
    }

    rows = []
    for analysis_group, cfg in comparison_defs.items():
        tumour_samples = sample_manifest.loc[cfg["mask"]].copy()
        if tumour_samples.empty or pooled_controls.empty:
            continue

        n_tumour = len(tumour_samples)
        n_control = len(pooled_controls)
        print(f"  {analysis_group}: tumour n={n_tumour}, pooled control n={n_control}")
        if n_tumour < MIN_CARRIERS or n_control < MIN_CARRIERS:
            continue

        tumour_df = df[df["Sample"].isin(tumour_samples["Sample"])]
        control_df = df[df["Sample"].isin(pooled_controls["Sample"])]
        all_variants = pd.concat(
            [tumour_df[["Variant_ID"]], control_df[["Variant_ID"]]],
            ignore_index=True,
        )["Variant_ID"].dropna().unique()

        for var_id in all_variants:
            tum_carriers = set(tumour_df.loc[tumour_df["Variant_ID"] == var_id, "Sample"].unique())
            ctl_carriers = set(control_df.loc[control_df["Variant_ID"] == var_id, "Sample"].unique())

            a = len(tum_carriers)
            b = n_tumour - a
            c = len(ctl_carriers)
            d = n_control - c

            if a + c < MIN_COMPARISON_CARRIERS:
                continue

            var_rows = df[df["Variant_ID"] == var_id]
            gene = var_rows["SYMBOL"].dropna().iloc[0] if "SYMBOL" in var_rows.columns and var_rows["SYMBOL"].notna().any() else ""
            consequence = var_rows["Consequence"].dropna().iloc[0] if "Consequence" in var_rows.columns and var_rows["Consequence"].notna().any() else np.nan
            impact = var_rows["IMPACT"].dropna().iloc[0] if "IMPACT" in var_rows.columns and var_rows["IMPACT"].notna().any() else np.nan

            _, p_val = fisher_exact([[a, b], [c, d]])
            or_val, or_method = _haldane_or(a, b, c, d)

            rows.append({
                "Analysis_Group": analysis_group,
                "Variant_ID": var_id,
                "Gene": gene,
                "Consequence": consequence,
                "IMPACT": impact,
                "N_Tumour": n_tumour,
                "N_Control": n_control,
                "Carriers_Tumour": a,
                "Carriers_Control": c,
                "Freq_Tumour_%": round(a / n_tumour * 100, 2),
                "Freq_Control_%": round(c / n_control * 100, 2),
                "Odds_Ratio": round(or_val, 4),
                "OR_Method": or_method,
                "P_Value": p_val,
            })

    res = pd.DataFrame(rows)
    if res.empty:
        print("  No tumour-vs-control results.\n")
        return res

    parts = []
    for grp in res["Analysis_Group"].unique():
        sub = res[res["Analysis_Group"] == grp].copy()
        sub = _apply_fdr(sub, p_col="P_Value", out_col="FDR_P_Value")
        parts.append(sub)
    res = pd.concat(parts, ignore_index=True).sort_values(["Analysis_Group", "P_Value", "Variant_ID"])
    res["Nominal_Sig"] = res["P_Value"] < 0.05
    res["FDR_Sig"] = res["FDR_P_Value"] < FDR_THRESHOLD

    print(f"  {len(res)} comparisons across {res['Analysis_Group'].nunique()} groups | "
          f"{res['Nominal_Sig'].sum()} nominal | {res['FDR_Sig'].sum()} FDR<{FDR_THRESHOLD}\n")
    return res

def _apply_fdr(df: pd.DataFrame, p_col: str, out_col: str) -> pd.DataFrame:
    df = df.copy()
    valid = df[p_col].notna()
    if valid.sum() > 1:
        df.loc[valid, out_col] = false_discovery_control(df.loc[valid, p_col].values, method="bh")
    else:
        df[out_col] = df[p_col]
    return df


def _zscore(series: pd.Series) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce")
    sd = vals.std()
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.zeros(len(vals)), index=vals.index)
    return (vals - vals.mean()) / sd


def _haldane_or(a: int, b: int, c: int, d: int) -> Tuple[float, str]:
    """Odds ratio with Haldane-Anscombe correction when needed."""
    if min(a, b, c, d) == 0:
        a, b, c, d = [x + 0.5 for x in (a, b, c, d)]
        return (a * d) / (b * c), "Haldane-Anscombe"
    return (a * d) / (b * c), "Standard"


# ?? STATISTICAL TESTS ?????????????????????????????????????????????????????????

def _test_continuous(c_vals, nc_vals) -> Optional[Dict]:
    """Mann-Whitney U between carriers and non-carriers. No age adjustment."""
    c  = pd.to_numeric(c_vals,  errors="coerce").dropna()
    nc = pd.to_numeric(nc_vals, errors="coerce").dropna()
    if len(c) < MIN_CARRIERS or len(nc) < MIN_CARRIERS:
        return None
    stat, p = mannwhitneyu(c.values, nc.values, alternative="two-sided")
    r_rb = round(1 - (2 * stat) / (len(c) * len(nc)), 4)
    return {
        "Test_Unadj":              "Mann-Whitney U",
        "N_Carriers":              len(c),
        "N_NonCarriers":           len(nc),
        "Median_Carriers":         round(c.median(), 3),
        "IQR_Carriers":            f"{round(c.quantile(0.25),3)}?{round(c.quantile(0.75),3)}",
        "Median_NonCarriers":      round(nc.median(), 3),
        "IQR_NonCarriers":         f"{round(nc.quantile(0.25),3)}?{round(nc.quantile(0.75),3)}",
        "Effect_RankBiserial":     r_rb,
        "Stat_Unadj":              round(stat, 3),
        "P_Unadj":                 p,
        "Test_Adj_Age":            "N/A",
        "OR_Adj_Age":              np.nan,
        "OR_Adj_Age_CI95":         np.nan,
        "P_Adj_Age":               np.nan,
        "N_Adj_Age":               np.nan,
        "Adj_Age_Note":            "Not applicable for continuous outcome",
        "Test_Adj_BMI":            "N/A",
        "OR_Adj_BMI":              np.nan,
        "OR_Adj_BMI_CI95":         np.nan,
        "P_Adj_BMI":               np.nan,
        "N_Adj_BMI":               np.nan,
        "Adj_BMI_Note":            "Not applicable for continuous outcome",
        "Test_Adj_AgeBMI":         "N/A",
        "OR_Adj_AgeBMI":           np.nan,
        "OR_Adj_AgeBMI_CI95":      np.nan,
        "P_Adj_AgeBMI":            np.nan,
        "N_Adj_AgeBMI":            np.nan,
        "Adj_AgeBMI_Note":         "Not applicable for continuous outcome",
    }


def _test_binary(c_df, nc_df, outcome_col, age_col="canon__age", bmi_col="canon__bmi", include_bmi_model: bool = True) -> Optional[Dict]:
    """
    Unadjusted Fisher exact test + three adjusted logistic regression models:
      Model 1 (age-only):    outcome ~ snp_carrier + age_z
      Model 2 (BMI-only):    outcome ~ snp_carrier + bmi_z
      Model 3 (age+BMI):     outcome ~ snp_carrier + age_z + bmi_z

    Each model reports its own N, OR, CI and p-value so missing BMI does not
    silently shrink the age-adjusted analysis.
    """
    try:
        from statsmodels.formula.api import logit as sm_logit
    except ImportError:
        sm_logit = None

    combined = pd.concat([
        c_df[[outcome_col, age_col, bmi_col]].assign(snp_carrier=1),
        nc_df[[outcome_col, age_col, bmi_col]].assign(snp_carrier=0),
    ])
    combined[outcome_col] = pd.to_numeric(combined[outcome_col], errors="coerce")
    combined = combined.dropna(subset=[outcome_col])

    n_c  = int((combined["snp_carrier"] == 1).sum())
    n_nc = int((combined["snp_carrier"] == 0).sum())
    if n_c < MIN_CARRIERS or n_nc < MIN_CARRIERS:
        return None

    a = int(((combined["snp_carrier"]==1) & (combined[outcome_col]==1)).sum())
    b = int(((combined["snp_carrier"]==1) & (combined[outcome_col]==0)).sum())
    c = int(((combined["snp_carrier"]==0) & (combined[outcome_col]==1)).sum())
    d = int(((combined["snp_carrier"]==0) & (combined[outcome_col]==0)).sum())
    if a + b == 0 or c + d == 0:
        return None

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
            or_v = round(float(np.exp(fit.params["snp"])), 4)
            ci   = np.exp(fit.conf_int().loc["snp"])
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

    age_df = combined[[outcome_col, age_col, "snp_carrier"]].dropna().copy()
    age_df.columns = ["outcome", "age", "snp"]
    age_df["age_z"] = _zscore(age_df["age"])
    result.update(_run_logit(age_df, "outcome ~ snp + age_z", "Age"))

    if include_bmi_model:
        bmi_only_df = combined[[outcome_col, bmi_col, "snp_carrier"]].dropna().copy()
        bmi_only_df.columns = ["outcome", "bmi", "snp"]
        bmi_only_df["bmi_z"] = _zscore(bmi_only_df["bmi"])
        result.update(_run_logit(bmi_only_df, "outcome ~ snp + bmi_z", "BMI"))

        bmi_df = combined[[outcome_col, age_col, bmi_col, "snp_carrier"]].dropna().copy()
        bmi_df.columns = ["outcome", "age", "bmi", "snp"]
        bmi_df["age_z"] = _zscore(bmi_df["age"])
        bmi_df["bmi_z"] = _zscore(bmi_df["bmi"])
        result.update(_run_logit(bmi_df, "outcome ~ snp + age_z + bmi_z", "AgeBMI"))
    else:
        result.update({
            "Test_Adj_BMI":        "N/A",
            "OR_Adj_BMI":          np.nan,
            "OR_Adj_BMI_CI95":     np.nan,
            "P_Adj_BMI":           np.nan,
            "N_Adj_BMI":           np.nan,
            "Adj_BMI_Note":        "BMI-only model not requested for this analysis",
            "Test_Adj_AgeBMI":     "N/A",
            "OR_Adj_AgeBMI":       np.nan,
            "OR_Adj_AgeBMI_CI95":  np.nan,
            "P_Adj_AgeBMI":        np.nan,
            "N_Adj_AgeBMI":        np.nan,
            "Adj_AgeBMI_Note":     "Age+BMI model not requested for this analysis",
        })

    return result


def _test_nominal(c_df, nc_df, col) -> Optional[Dict]:
    """Chi-square for multi-category variables. Unadjusted only."""
    n_c  = c_df[col].notna().sum()
    n_nc = nc_df[col].notna().sum()
    if n_c < MIN_CARRIERS or n_nc < MIN_CARRIERS:
        return None
    grp_s = pd.Series(["C"] * len(c_df) + ["NC"] * len(nc_df), name="grp")
    val_s = pd.concat(
        [c_df[col].reset_index(drop=True), nc_df[col].reset_index(drop=True)]
    ).reset_index(drop=True)
    ct = pd.crosstab(grp_s, val_s).dropna(axis=1)
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
        "Test_Adj_BMI":        "N/A",
        "OR_Adj_BMI":          np.nan,
        "OR_Adj_BMI_CI95":     np.nan,
        "P_Adj_BMI":           np.nan,
        "N_Adj_BMI":           np.nan,
        "Adj_BMI_Note":        "Not applicable for nominal outcome",
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

    PRIMARY MODEL — Genotypic (Het vs WT, Hom vs WT):
      Each variant × clinical variable produces up to two rows:
        Contrast = "Het_vs_WT": heterozygous carriers vs wildtype
        Contrast = "Hom_vs_WT": homozygous carriers vs wildtype
                                 (only when n_Hom >= MIN_CARRIERS)
      This is the biologically correct model for GSDMB: isoform expression
      is dose-dependent (Hom minor allele completely loses pyroptotic isoforms;
      Het has partial expression). Lumping Het+Hom as "carrier" masks this.

    SECONDARY MODEL — Additive (any carrier vs WT):
      Het + Hom pooled vs WT. Stored in *_Additive columns of the same row
      for completeness and comparability with published literature.

    GT column must be present and normalised (0/0, 0/1, 1/1).
    If GT is absent, falls back to additive-only (original behaviour).

    FDR correction: BH applied within each clinical variable × contrast
    across all variants (correct unit — not pooled across variables).
    """
    GT_MAP = {"0/0": "WT", "0/1": "Het", "1/0": "Het", "1/1": "Hom"}
    has_gt = "GT" in df_t.columns

    rows = []
    for var_id, vdf in df_t.groupby("Variant_ID"):
        sym         = vdf["SYMBOL"].dropna().iloc[0] if vdf["SYMBOL"].notna().any() else ""
        sample_data = df_t.drop_duplicates("Sample").set_index("Sample").copy()

        # ── Assign genotype groups ────────────────────────────────────────
        if has_gt:
            gt_lookup = (
                vdf.drop_duplicates("Sample")
                   .set_index("Sample")["GT"]
                   .apply(lambda g: GT_MAP.get(str(g), "WT"))
            )
            sample_data["_geno"] = sample_data.index.map(gt_lookup).fillna("WT")
        else:
            # Fallback: anyone in vdf rows is a carrier (additive only)
            carrier_set = set(vdf["Sample"].unique())
            sample_data["_geno"] = sample_data.index.map(
                lambda s: "Het" if s in carrier_set else "WT"
            )

        wt_df  = sample_data[sample_data["_geno"] == "WT"]
        het_df = sample_data[sample_data["_geno"] == "Het"]
        hom_df = sample_data[sample_data["_geno"] == "Hom"]
        # Additive: Het + Hom pooled
        carrier_df = sample_data[sample_data["_geno"].isin(["Het", "Hom"])]

        n_wt  = len(wt_df)
        n_het = len(het_df)
        n_hom = len(hom_df)

        if n_wt < MIN_CARRIERS:
            continue  # need a reference group

        # Contrasts to run: (label, test_df, ref_df)
        contrasts = []
        if n_het >= MIN_CARRIERS:
            contrasts.append(("Het_vs_WT", het_df, wt_df))
        if n_hom >= MIN_CARRIERS:
            contrasts.append(("Hom_vs_WT", hom_df, wt_df))
        if not contrasts:
            continue

        for col, meta in var_dict.items():
            if col not in sample_data.columns:
                continue

            # ── Additive (secondary) ──────────────────────────────────────
            if len(carrier_df) >= MIN_CARRIERS:
                if meta["type"] == "continuous":
                    add_res = _test_continuous(carrier_df[col], wt_df[col])
                elif meta["type"] == "binary":
                    add_res = _test_binary(carrier_df, wt_df, col,
                                           age_col="canon__age", bmi_col="canon__bmi",
                                           include_bmi_model=True)
                elif meta["type"] == "nominal":
                    add_res = _test_nominal(carrier_df, wt_df, col)
                else:
                    add_res = None
            else:
                add_res = None

            # Rename additive result keys with _Additive suffix
            additive_cols = {}
            if add_res:
                for k, v in add_res.items():
                    additive_cols[f"{k}_Additive"] = v

            # ── Genotypic contrasts (primary) ─────────────────────────────
            for contrast_label, test_df, ref_df in contrasts:
                if meta["type"] == "continuous":
                    res = _test_continuous(test_df[col], ref_df[col])
                elif meta["type"] == "binary":
                    res = _test_binary(test_df, ref_df, col,
                                       age_col="canon__age", bmi_col="canon__bmi",
                                       include_bmi_model=True)
                elif meta["type"] == "nominal":
                    res = _test_nominal(test_df, ref_df, col)
                else:
                    res = None

                if res is None and add_res is None:
                    continue

                # Use genotypic result if available, else NaN placeholders
                if res is None:
                    res = {k: np.nan for k in (add_res or {}).keys()}

                row = {
                    "Cohort":        cohort_label,
                    "Variant_ID":    var_id,
                    "Gene":          sym,
                    "Contrast":      contrast_label,
                    "N_WT":          n_wt,
                    "N_Het":         n_het,
                    "N_Hom":         n_hom,
                    "Clinical_Var":  col,
                    "Clin_Label":    meta["label"],
                    "Clin_Type":     meta["type"],
                    "Note":          meta.get("note", ""),
                }
                row.update(res)
                row.update(additive_cols)
                rows.append(row)

    res_df = pd.DataFrame(rows)
    if res_df.empty:
        return res_df

    # ── FDR: within each clinical variable × contrast ─────────────────────
    # This is the correct unit: across all variants for a given outcome + contrast.
    parts = []
    for (clin_col, contrast) in res_df[["Clinical_Var", "Contrast"]].drop_duplicates().itertuples(index=False):
        mask = (res_df["Clinical_Var"] == clin_col) & (res_df["Contrast"] == contrast)
        sub  = res_df[mask].copy()
        sub  = _apply_fdr(sub, p_col="P_Unadj", out_col="FDR_Unadj")
        if "P_Adj_Age" in sub.columns and sub["P_Adj_Age"].notna().any():
            sub = _apply_fdr(sub, p_col="P_Adj_Age", out_col="FDR_Adj_Age")
        else:
            sub["FDR_Adj_Age"] = np.nan
        if "P_Adj_BMI" in sub.columns and sub["P_Adj_BMI"].notna().any():
            sub = _apply_fdr(sub, p_col="P_Adj_BMI", out_col="FDR_Adj_BMI")
        else:
            sub["FDR_Adj_BMI"] = np.nan
        if "P_Adj_AgeBMI" in sub.columns and sub["P_Adj_AgeBMI"].notna().any():
            sub = _apply_fdr(sub, p_col="P_Adj_AgeBMI", out_col="FDR_Adj_AgeBMI")
        else:
            sub["FDR_Adj_AgeBMI"] = np.nan
        parts.append(sub)

    res_df = pd.concat(parts, ignore_index=True).sort_values(
        ["Clinical_Var", "Contrast", "P_Unadj"]
    )
    res_df["Nominal_Sig_Unadj"]      = res_df["P_Unadj"] < 0.05
    res_df["FDR_Sig_Unadj"]          = res_df["FDR_Unadj"] < FDR_THRESHOLD
    res_df["Nominal_Sig_Adj_Age"]    = res_df["P_Adj_Age"].notna() & (res_df["P_Adj_Age"] < 0.05)
    res_df["FDR_Sig_Adj_Age"]        = res_df["FDR_Adj_Age"].notna() & (res_df["FDR_Adj_Age"] < FDR_THRESHOLD)
    res_df["Nominal_Sig_Adj_BMI"]    = res_df["P_Adj_BMI"].notna() & (res_df["P_Adj_BMI"] < 0.05)
    res_df["FDR_Sig_Adj_BMI"]        = res_df["FDR_Adj_BMI"].notna() & (res_df["FDR_Adj_BMI"] < FDR_THRESHOLD)
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
        if res.empty:
            print(f"  {label}: no results")
        else:
            for contrast in sorted(res["Contrast"].unique()):
                sub  = res[res["Contrast"] == contrast]
                n_age = int(sub["N_Adj_Age"].dropna().astype(int).max()) if "N_Adj_Age" in sub.columns and sub["N_Adj_Age"].notna().any() else 0
                n_bmi = int(sub["N_Adj_BMI"].dropna().astype(int).max()) if "N_Adj_BMI" in sub.columns and sub["N_Adj_BMI"].notna().any() else 0
                n_age_bmi = int(sub["N_Adj_AgeBMI"].dropna().astype(int).max()) if "N_Adj_AgeBMI" in sub.columns and sub["N_Adj_AgeBMI"].notna().any() else 0
                print(f"  {label} [{contrast}]: {len(sub)} tests | "
                      f"{sub['Nominal_Sig_Unadj'].sum()} nominal | "
                      f"{sub['FDR_Sig_Unadj'].sum()} FDR | "
                      f"{sub['Nominal_Sig_Adj_Age'].sum()} nominal age-adj (max N={n_age}) | "
                      f"{sub['Nominal_Sig_Adj_BMI'].sum()} nominal BMI-adj (max N={n_bmi}) | "
                      f"{sub['Nominal_Sig_Adj_AgeBMI'].sum()} nominal age+BMI-adj (max N={n_age_bmi})")
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
            grp_series = pd.Series(
                ["ref"] * len(c0) + ["test"] * len(c1), name="grp"
            )
            val_series = pd.concat(
                [c0.reset_index(drop=True), c1.reset_index(drop=True)]
            ).reset_index(drop=True)
            ct = pd.crosstab(grp_series, val_series)
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
                # Spearman rho of dose (0/1/2) vs outcome — tests monotonic
                # ordered alternative (dose-response), correct for genotypic trend.
                dose_arr = np.concatenate([
                    np.full(len(sd[sd["dose"]==g][col].dropna()), g)
                    for g in [0, 1, 2]
                ])
                out_arr = np.concatenate([
                    pd.to_numeric(sd[sd["dose"]==g][col], errors="coerce").dropna().values
                    for g in [0, 1, 2]
                ])
                if len(dose_arr) < 4:
                    continue
                stat, p = stats.spearmanr(dose_arr, out_arr)
                test = "Spearman trend (dose-response)"
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
#   Endometrial (AT=AUs) — OS  : canon__os_months  + ENDO_EXITUS_DISEASE_BIN (disease-specific)
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
            "OS":  {"t_col": "canon__os_months",        "ev_col": "ENDO_EXITUS_DISEASE_BIN"},  # disease-specific is primary
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

_GENO_COLS  = GENOTYPE_COLORS.copy()


def survival_analysis(df: pd.DataFrame) -> Tuple[pd.DataFrame, List]:
    if not _HAS_LIFELINES:
        print("=== Analysis 4: Survival — SKIPPED (lifelines not installed) ===\n")
        return pd.DataFrame(), []

    print("=== Analysis 4: Survival — all cohorts ===")
    df = df[~df["is_replicate"] & (df["Tissue"] == "Tumour")].copy()
    if "GT" not in df.columns:
        print("  No GT column — genotypic model unavailable.\n")
        return pd.DataFrame(), []

    rows      = []
    km_pages  = []

    def _prepare_cox_frame(base_df: pd.DataFrame, model_terms: List[str]) -> Tuple[pd.DataFrame, List[str], str, str]:
        cox_df = base_df.copy()
        total_n = len(cox_df)
        age_complete = int(cox_df["age"].notna().sum()) if "age" in cox_df.columns else 0
        bmi_complete = int(cox_df["bmi"].notna().sum()) if "bmi" in cox_df.columns else 0

        if "age" in cox_df.columns and age_complete >= 6:
            if "bmi" in cox_df.columns and bmi_complete >= 6:
                cox_df = cox_df.dropna(subset=["age", "bmi"]).copy()
                if len(cox_df) >= 6:
                    cox_df["age_z"] = _zscore(cox_df["age"])
                    cox_df["bmi_z"] = _zscore(cox_df["bmi"])
                    note = (
                        f"Age+BMI-adjusted Cox on BMI-complete subset ({len(cox_df)}/{total_n} samples)."
                        if len(cox_df) < total_n
                        else "Age+BMI-adjusted Cox on full analyzable set."
                    )
                    return cox_df, model_terms + ["age_z", "bmi_z"], "Age+BMI", note
            cox_df = cox_df.dropna(subset=["age"]).copy()
            if len(cox_df) >= 6:
                cox_df["age_z"] = _zscore(cox_df["age"])
                note = (
                    f"BMI adjustment not feasible: {bmi_complete}/{total_n} BMI-complete samples; age-adjusted Cox used."
                    if "bmi" in base_df.columns
                    else "BMI column unavailable; age-adjusted Cox used."
                )
                return cox_df, model_terms + ["age_z"], "Age", note

        note = (
            f"Age adjustment not feasible: {age_complete}/{total_n} age-complete samples; unadjusted Cox used."
            if "age" in base_df.columns
            else "Age column unavailable; unadjusted Cox used."
        )
        return cox_df.copy(), model_terms.copy(), "Unadjusted", note

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
                if "canon__bmi" in s.columns: needed.append("canon__bmi")
                valid = s[needed].copy()
                valid[t_col]  = pd.to_numeric(valid[t_col],  errors="coerce")
                valid[ev_col] = pd.to_numeric(valid[ev_col], errors="coerce")
                valid = valid.dropna(subset=[t_col, ev_col])
                valid = valid[valid[t_col] > 0]   # lifelines requires T > 0
                if len(valid) < 6: continue

                T    = valid[t_col]
                E    = valid[ev_col].astype(float)
                geno = s.loc[valid.index, "geno"].fillna("WT")
                dose = s.loc[valid.index, "dose"].fillna(0).astype(int)

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
                cox_ph_flag = "PH_check_unavailable"
                cox_adjustment = "Unadjusted"
                cox_adjustment_note = "Cox model not fitted"
                cox_n = np.nan
                try:
                    cox_base = pd.DataFrame({"T": T, "E": E, "carrier": carrier_flag})
                    if age_col in valid.columns:
                        cox_base["age"] = pd.to_numeric(valid[age_col], errors="coerce")
                    if "canon__bmi" in valid.columns:
                        cox_base["bmi"] = pd.to_numeric(valid["canon__bmi"], errors="coerce")
                    cox_df, formula_cols, cox_adjustment, cox_adjustment_note = _prepare_cox_frame(cox_base, ["carrier"])
                    cox_n = len(cox_df)
                    if len(cox_df) >= 6 and cox_df["carrier"].nunique() > 1:
                        cph = CoxPHFitter()
                        cph.fit(cox_df[["T", "E"] + formula_cols], duration_col="T", event_col="E",
                                formula=" + ".join(formula_cols), show_progress=False)
                        cox_add_hr = round(float(np.exp(cph.params_["carrier"])), 3)
                        ci = np.exp(cph.confidence_intervals_.loc["carrier"])
                        cox_add_p  = round(float(cph.summary.loc["carrier", "p"]), 4)
                        cox_add_ci = f"[{ci.iloc[0]:.3f}, {ci.iloc[1]:.3f}]"
                        try:
                            ph_results = cph.check_assumptions(cox_df[["T", "E"] + formula_cols], p_value_threshold=0.05,
                                                               show_plots=False)
                            ph_pvals = [r.p_value for r in ph_results if hasattr(r, "p_value")]
                            cox_ph_flag = (f"PH_WARNING(min_p={min(ph_pvals):.3f})"
                                           if ph_pvals and min(ph_pvals) < 0.05
                                           else "PH_OK")
                        except Exception:
                            cox_ph_flag = "PH_check_unavailable"
                except Exception:
                    pass

                # ?? (c) Genotypic Cox (WT reference, Het + Hom terms) ????
                geno_het_hr = geno_het_ci = geno_het_p = np.nan
                geno_hom_hr = geno_hom_ci = geno_hom_p = np.nan
                model_type = "Additive_only"

                if n_hom >= MIN_CARRIERS:
                    model_type = "Genotypic"
                    try:
                        gdf_base = pd.DataFrame({
                            "T": T, "E": E,
                            "Het": (geno == "Het").astype(int),
                            "Hom": (geno == "Hom").astype(int),
                        })
                        if age_col in valid.columns:
                            gdf_base["age"] = pd.to_numeric(valid[age_col], errors="coerce")
                        if "canon__bmi" in valid.columns:
                            gdf_base["bmi"] = pd.to_numeric(valid["canon__bmi"], errors="coerce")
                        gdf, geno_terms, _, _ = _prepare_cox_frame(gdf_base, ["Het", "Hom"])
                        if len(gdf) >= 6:
                            cph_g = CoxPHFitter()
                            cph_g.fit(gdf[["T", "E"] + geno_terms], duration_col="T", event_col="E",
                                      formula=" + ".join(geno_terms), show_progress=False)
                            for term, hr_var, ci_var, p_var in [
                                ("Het", "geno_het_hr", "geno_het_ci", "geno_het_p"),
                                ("Hom", "geno_hom_hr", "geno_hom_ci", "geno_hom_p"),
                            ]:
                                if term in cph_g.params_.index:
                                    hr_v = round(float(np.exp(cph_g.params_[term])), 3)
                                    ci_v = np.exp(cph_g.confidence_intervals_.loc[term])
                                    p_v  = round(float(cph_g.summary.loc[term, "p"]), 4)
                                    ci_s = f"[{ci_v.iloc[0]:.3f}, {ci_v.iloc[1]:.3f}]"
                                    if hr_var == "geno_het_hr":
                                        geno_het_hr = hr_v; geno_het_ci = ci_s; geno_het_p = p_v
                                    else:
                                        geno_hom_hr = hr_v; geno_hom_ci = ci_s; geno_hom_p = p_v
                    except Exception:
                        pass
                elif n_het >= MIN_CARRIERS:
                    model_type = "Het_vs_WT"

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
                    "N_Cox":            cox_n,
                    "Cox_Adjustment":   cox_adjustment,
                    "Cox_Adjustment_Note": cox_adjustment_note,
                    # Log-rank
                    "Logrank_P":        lr_p,
                    # Additive Cox (any carrier vs WT, explicit adjustment)
                    "HR_Additive":      cox_add_hr,
                    "HR_Additive_CI95": cox_add_ci,
                    "P_Cox_Additive":   cox_add_p,
                    "PH_Assumption_Flag": cox_ph_flag,
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
                plt.tight_layout(rect=[0, 0, 0.9, 1])
                km_pages.append(fig)

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

    return res, km_pages


# ── VISUALISATION ─────────────────────────────────────────────────────────────
# Shared style
_PALETTE   = {"carrier": arm_color('Tumour'), "non_carrier": arm_color('Control')}
_IMPACT_C  = IMPACT_COLORS.copy()
_COHORT_C  = {"Breast": cohort_color('Breast'), "Endometrial": cohort_color('Endometrial')}

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


def _short_label(value, max_len=34):
    if pd.isna(value):
        return ""
    text = re.sub(r"\s+", " ", str(value)).strip()
    return text if len(text) <= max_len else text[: max_len - 1].rstrip() + "..."


def _wrap_label(value, width=20, max_lines=3):
    text = re.sub(r"\s+", " ", str(value)).strip()
    wrapped = textwrap.wrap(text, width=width, break_long_words=False, break_on_hyphens=False)
    if not wrapped:
        return text
    if len(wrapped) > max_lines:
        tail = " ".join(wrapped[max_lines - 1:])
        wrapped = wrapped[: max_lines - 1] + [_short_label(tail, max_len=width)]
    return "\n".join(wrapped)


def _set_wrapped_ticklabels(ax, axis="x", width=18, max_lines=3, rotation=0, ha="center", fontsize=9):
    labels = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    wrapped = [_wrap_label(label.get_text(), width=width, max_lines=max_lines) for label in labels]
    if axis == "x":
        ax.set_xticklabels(wrapped, rotation=rotation, ha=ha, fontsize=fontsize)
    else:
        ax.set_yticklabels(wrapped, rotation=rotation, ha=ha, fontsize=fontsize)


def _add_margin_labels(ax, points, side="right", max_labels=5, fontsize=7.5):
    if not points:
        return
    points = sorted(points, key=lambda item: item[1], reverse=True)[:max_labels]
    y_min, y_max = ax.get_ylim()
    pad = (y_max - y_min) * 0.06 if y_max != y_min else 0.5
    y_positions = np.linspace(y_max - pad, y_min + pad, len(points))
    x_text = 1.03 if side == "right" else -0.03
    ha = "left" if side == "right" else "right"
    for (x_val, y_val, label), y_text in zip(points, y_positions):
        ax.annotate(
            _short_label(label, 28),
            xy=(x_val, y_val),
            xycoords="data",
            xytext=(x_text, y_text),
            textcoords=ax.get_yaxis_transform(),
            ha=ha,
            va="center",
            fontsize=fontsize,
            clip_on=False,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#bdbdbd", alpha=0.92, linewidth=0.8),
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.8),
        )


# ── 1. VOLCANO (tumour vs control) ───────────────────────────────────────────

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

        # Annotate the strongest hits in the plot margins to avoid label pile-ups
        hits = sub[sub[sig_col]].nsmallest(10, p_col)
        left_hits = []
        right_hits = []
        for _, row in hits.iterrows():
            label = row.get("Gene", row["Variant_ID"]) or row["Variant_ID"]
            point = (float(row["log2OR"]), float(row["-log10p"]), label)
            if point[0] < 0:
                left_hits.append(point)
            else:
                right_hits.append(point)
        _add_margin_labels(ax, left_hits, side="left", max_labels=5)
        _add_margin_labels(ax, right_hits, side="right", max_labels=5)

        # Axis labels and styling
        ax.set_xlabel("log₂(Odds Ratio)  [tumour enriched →]", fontsize=10)
        ax.set_ylabel("-log₁₀(p-value)", fontsize=10)
        if not sub.empty:
            n_tumour = int(sub['N_Tumour'].iloc[0])
            n_control = int(sub['N_Control'].iloc[0])
            ax.set_title(f"{grp} | Tumour n={n_tumour}, Pooled control n={n_control}", fontsize=12, fontweight="bold", pad=8)
        else:
            ax.set_title(f"{grp}", fontsize=12, fontweight="bold", pad=8)
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

    fig.suptitle(
        f"GSDMB SNP Volcano ({title_suffix})",
        fontsize=13,
        fontweight="bold",
        y=1.01,
    )
    plt.tight_layout(rect=[0.06, 0, 0.94, 1])
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def make_volcano_plots(tvh, out_dir):
    if tvh.empty: return
    _draw_volcano(tvh, "P_Value",     "Nominal_Sig", 0.05,
                  "unadjusted p",  out_dir / "17_SNP_Volcano_raw_p.png")
    _draw_volcano(tvh, "FDR_P_Value", "FDR_Sig",     FDR_THRESHOLD,
                  "FDR-corrected", out_dir / "17_SNP_Volcano_FDR.png")


def _filter_significant_clinical_rows(clin_res: pd.DataFrame, p_col: str | None = None) -> pd.DataFrame:
    """Keep only statistically significant rows for clinical-association figures."""
    if clin_res is None or clin_res.empty:
        return pd.DataFrame()
    sub = clin_res.copy()
    if p_col and p_col in sub.columns:
        cutoff = FDR_THRESHOLD if "FDR" in p_col else 0.05
        return sub[sub[p_col].notna() & (sub[p_col] < cutoff)].copy()

    mask = pd.Series(False, index=sub.index)
    for col in ["Nominal_Sig_Adj_AgeBMI", "Nominal_Sig_Adj_BMI", "Nominal_Sig_Adj_Age", "Nominal_Sig_Unadj"]:
        if col in sub.columns:
            mask |= sub[col].fillna(False).astype(bool)
    return sub[mask].copy()


# ── 2. HEATMAP (SNP × clinical variable) ─────────────────────────────────────

def _draw_heatmap(clin_res, p_col, title_suffix, out_path, filter_thresh=0.20):
    # Label rows as "VariantID [Contrast]" so Het and Hom are separate rows
    clin_res = clin_res.copy()
    if "Contrast" in clin_res.columns:
        clin_res["_row_label"] = clin_res["Variant_ID"] + " [" + clin_res["Contrast"] + "]"
    else:
        clin_res["_row_label"] = clin_res["Variant_ID"]

    pivot = (clin_res
             .pivot_table(index="_row_label", columns="Clin_Label",
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

    ax.set_title(
        f"SNP-Clinical Heatmap ({title_suffix})",
        fontsize=11,
        fontweight="bold",
        pad=10,
    )
    ax.set_xlabel("Clinical variable", fontsize=9, labelpad=8)
    ax.set_ylabel("Variant", fontsize=9, labelpad=8)
    _set_wrapped_ticklabels(ax, axis="x", width=18, max_lines=3, rotation=0, ha="center", fontsize=9)
    _set_wrapped_ticklabels(ax, axis="y", width=26, max_lines=2, rotation=0, ha="right", fontsize=8.5)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def make_heatmaps(breast_clin, endo_clin, out_dir, significant_only=False):
    suffix = " (significant only)" if significant_only else ""
    raw_thresh = 0.05 if significant_only else 1.01
    fdr_thresh = FDR_THRESHOLD if significant_only else 1.01
    for clin_res, label in [(breast_clin, "Breast"), (endo_clin, "Endometrial")]:
        if clin_res.empty:
            continue
        raw_res = _filter_significant_clinical_rows(clin_res, "P_Unadj") if significant_only else clin_res[clin_res["P_Unadj"].notna()].copy()
        if not raw_res.empty:
            _draw_heatmap(raw_res, "P_Unadj",
                          f"Unadjusted p{suffix} ? {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_raw_p.png",
                          filter_thresh=raw_thresh)
        fdr_res = _filter_significant_clinical_rows(clin_res, "FDR_Unadj") if significant_only else clin_res[clin_res["FDR_Unadj"].notna()].copy()
        if not fdr_res.empty:
            _draw_heatmap(fdr_res, "FDR_Unadj",
                          f"FDR-corrected{suffix} ? {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_FDR.png",
                          filter_thresh=fdr_thresh)
        adj_res = _filter_significant_clinical_rows(clin_res, "P_Adj_Age") if significant_only else clin_res[clin_res["P_Adj_Age"].notna()].copy()
        if not adj_res.empty:
            _draw_heatmap(adj_res, "P_Adj_Age",
                          f"Age-adjusted p{suffix} ? {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_age_raw_p.png",
                          filter_thresh=raw_thresh)
        adj_fdr = _filter_significant_clinical_rows(clin_res, "FDR_Adj_Age") if significant_only else clin_res[clin_res["FDR_Adj_Age"].notna()].copy()
        if not adj_fdr.empty:
            _draw_heatmap(adj_fdr, "FDR_Adj_Age",
                          f"Age-adjusted FDR{suffix} ? {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_age_FDR.png",
                          filter_thresh=fdr_thresh)
        bmi_only_res = _filter_significant_clinical_rows(clin_res, "P_Adj_BMI") if significant_only else clin_res[clin_res["P_Adj_BMI"].notna()].copy()
        if not bmi_only_res.empty:
            _draw_heatmap(bmi_only_res, "P_Adj_BMI",
                          f"BMI-adjusted p{suffix} ? {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_bmi_raw_p.png",
                          filter_thresh=raw_thresh)
        bmi_only_fdr = _filter_significant_clinical_rows(clin_res, "FDR_Adj_BMI") if significant_only else clin_res[clin_res["FDR_Adj_BMI"].notna()].copy()
        if not bmi_only_fdr.empty:
            _draw_heatmap(bmi_only_fdr, "FDR_Adj_BMI",
                          f"BMI-adjusted FDR{suffix} ? {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_bmi_FDR.png",
                          filter_thresh=fdr_thresh)
        bmi_res = _filter_significant_clinical_rows(clin_res, "P_Adj_AgeBMI") if significant_only else clin_res[clin_res["P_Adj_AgeBMI"].notna()].copy()
        if not bmi_res.empty:
            _draw_heatmap(bmi_res, "P_Adj_AgeBMI",
                          f"Age+BMI-adjusted p{suffix} ? {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_agebmi_raw_p.png",
                          filter_thresh=raw_thresh)
        bmi_fdr = _filter_significant_clinical_rows(clin_res, "FDR_Adj_AgeBMI") if significant_only else clin_res[clin_res["FDR_Adj_AgeBMI"].notna()].copy()
        if not bmi_fdr.empty:
            _draw_heatmap(bmi_fdr, "FDR_Adj_AgeBMI",
                          f"Age+BMI-adjusted FDR{suffix} ? {label}",
                          out_dir / f"17_SNP_Heatmap_{label}_adj_agebmi_FDR.png",
                          filter_thresh=fdr_thresh)


# ?? 3. FOREST PLOT (OR with 95% CI) ??????????????????????????????????????????

def _select_heatmap_focus_rows(clin_res: pd.DataFrame, p_col: str, max_rows: int = 12) -> pd.DataFrame:
    """Return a cleaner subset of heatmap rows for main-text presentation."""
    if clin_res.empty or p_col not in clin_res.columns:
        return clin_res

    focus = clin_res[clin_res[p_col].notna()].copy()
    if focus.empty:
        return clin_res

    if "Contrast" in focus.columns:
        focus["_focus_row"] = focus["Variant_ID"].astype(str) + " [" + focus["Contrast"].astype(str) + "]"
    else:
        focus["_focus_row"] = focus["Variant_ID"].astype(str)

    rank_df = (focus.groupby("_focus_row", as_index=False)[p_col]
               .min()
               .sort_values(p_col))

    if "FDR_Unadj" in focus.columns and (focus["FDR_Unadj"] < FDR_THRESHOLD).any():
        keep = rank_df[rank_df["_focus_row"].isin(focus.loc[focus["FDR_Unadj"] < FDR_THRESHOLD, "_focus_row"])].head(max_rows)
    elif "Nominal_Sig_Unadj" in focus.columns and focus["Nominal_Sig_Unadj"].any():
        keep = rank_df[rank_df["_focus_row"].isin(focus.loc[focus["Nominal_Sig_Unadj"], "_focus_row"])].head(max_rows)
    else:
        keep = rank_df.head(max_rows)

    keep_rows = set(keep["_focus_row"])
    out = focus[focus["_focus_row"].isin(keep_rows)].copy()
    return out.drop(columns=["_focus_row"], errors="ignore")


def make_main_text_heatmaps(breast_clin, endo_clin, out_dir):
    """Generate cleaner focused heatmaps for main-text presentation."""
    for clin_res, label in [(breast_clin, "Breast"), (endo_clin, "Endometrial")]:
        if clin_res.empty:
            continue
        significant = _filter_significant_clinical_rows(clin_res)
        if significant.empty:
            continue
        focused = _select_heatmap_focus_rows(significant, "P_Unadj", max_rows=12)
        _draw_heatmap(
            focused,
            "P_Unadj",
            f"Main-text overview ? {label}",
            out_dir / f"17_SNP_Heatmap_{label}_MainText.png",
            filter_thresh=0.05,
        )


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
            clin_res["OR_Adj_Age_CI95"].notna() &
            clin_res["P_Adj_Age"].notna() &
            (clin_res["P_Adj_Age"] < 0.05)
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
        fig, ax = plt.subplots(figsize=(10.5, fig_h))

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

        # p-value labels on the right margin stay readable even when CIs are wide
        ax.set_xlim(left=0, right=OR_CAP * 1.10)
        for i, (_, row) in enumerate(sub.iterrows()):
            pstr = f"p={row['P_Adj_Age']:.3f}{_sig_label(row['P_Adj_Age'])}"
            ax.text(1.02, ypos[i], pstr, transform=ax.get_yaxis_transform(),
                    va="center", ha="left", fontsize=6.8, color="#444444", clip_on=False)

        ax.axvline(1, color="#555555", linestyle="--", linewidth=1, alpha=0.7)
        ax.set_yticks(ypos)
        ax.set_yticklabels([_wrap_label(v, width=34, max_lines=2) for v in sub["label"].values], fontsize=8)
        ax.set_xlabel("Odds Ratio (age-adjusted, 95% CI)", fontsize=10)
        ax.set_title(
            f"{label}: Age-Adjusted Associations of GSDMB SNP Carrier Status With Binary Clinical Outcomes",
            fontsize=11,
            fontweight="bold",
            pad=10,
        )
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

        plt.tight_layout(rect=[0, 0, 0.9, 1])
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
            clin_res["P_Unadj"].notna() &
            (clin_res["P_Unadj"] < 0.05)
        ].sort_values("P_Unadj")
         .drop_duplicates(["Clinical_Var", "Contrast"] if "Contrast" in clin_res.columns else ["Clinical_Var"])
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
            suffix = "age-adj" if pd.notna(p_a) else "unadj"
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
            gene  = row.get("Gene", "") or ""
            snp   = var_id
            contrast_lbl = f" [{row['Contrast']}]" if "Contrast" in row and pd.notna(row.get("Contrast")) else ""
            title = f"{gene}  [{snp}]{contrast_lbl}" if gene and gene != snp else f"{snp}{contrast_lbl}"
            ax.set_title(f"Variant-Level Distribution Plot\n{title}", fontsize=9, fontweight="bold")
            _style_ax(ax)

        # Hide unused subplots
        for ax in axes_flat[n_plots:]:
            ax.set_visible(False)

        fig.suptitle(
            f"Continuous Clinical Variables According to GSDMB SNP Carrier Status in {label}",
            fontsize=12,
            fontweight="bold",
            y=1.01,
        )
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

    Fix applied (v4.1): three bugs caused the canvas to render ~10-20x
    taller than expected, squashing the chart at the very bottom:
      1. ax.text(xi, -4, ...) and ax.text(0.5, 105, ...) placed text
         outside the data limits; bbox_inches='tight' then expanded
         the canvas to capture them.
      2. suptitle(y=1.01) sat above the figure boundary, adding further
         vertical expansion via bbox_inches='tight'.
    Fix: all annotations now use transform=ax.transAxes or
         transform=ax.get_xaxis_transform() so coordinates are always
         0-1 and never escape the axes/figure bounds.
         suptitle moved to y=0.98 with tight_layout(rect=[0,0,1,0.96]).
    """
    TOP_N = 12

    for clin_res, label, coh_sheet, colour in [
        (breast_clin, "Breast",      "MT-T_N",  _COHORT_C["Breast"]),
        (endo_clin,   "Endometrial", "AT=AUs",  _COHORT_C["Endometrial"]),
    ]:
        if clin_res.empty:
            continue

        sig_clin = _filter_significant_clinical_rows(clin_res)
        bin_hits = (
            sig_clin[
                (sig_clin["Clin_Type"] == "binary") &
                sig_clin["P_Unadj"].notna()
            ]
            .sort_values("P_Unadj")
            .drop_duplicates(["Clinical_Var", "Contrast"] if "Contrast" in clin_res.columns else ["Clinical_Var"])
            .head(TOP_N)
        )
        if bin_hits.empty:
            continue

        n_plots = len(bin_hits)
        ncols   = min(3, n_plots)
        nrows   = int(np.ceil(n_plots / ncols))

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(4.5 * ncols, 3.8 * nrows),
            squeeze=False,
            facecolor="white",
        )
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
                merged[
                    (merged["Variant_ID"] == var_id) &
                    (merged["sheet"] == coh_sheet)
                ]["Sample"].unique()
            )
            cohort_df["_carrier"] = cohort_df["Sample"].map(
                lambda s: "Carrier" if s in var_samples else "Non-carrier"
            )

            plot_df = cohort_df[[col, "_carrier"]].dropna()
            if plot_df.empty or plot_df["_carrier"].nunique() < 2:
                ax.set_visible(False)
                continue

            groups    = ["Carrier", "Non-carrier"]
            pos_rates = []
            ns        = []
            for grp in groups:
                vals = plot_df.loc[plot_df["_carrier"] == grp, col]
                pos_rates.append(100 * vals.mean() if len(vals) else 0)
                ns.append(len(vals))

            xpos = [0, 1]
            ax.bar(
                xpos, pos_rates,
                color=[_PALETTE["carrier"], _PALETTE["non_carrier"]],
                width=0.55, edgecolor="white", linewidth=1.2,
                alpha=0.85, zorder=3,
            )

            # Complement bars (light, showing 100 - positive%)
            for xi, (r, c) in enumerate(
                zip(pos_rates, [_PALETTE["carrier"], _PALETTE["non_carrier"]])
            ):
                ax.bar(
                    xi, 100 - r, bottom=r, color=c,
                    width=0.55, alpha=0.18, edgecolor="none", zorder=2,
                )

            # Value labels inside bars — only when bar is tall enough
            for xi, r in enumerate(pos_rates):
                if r > 8:
                    ax.text(
                        xi, r / 2, f"{r:.0f}%",
                        ha="center", va="center",
                        fontsize=9, fontweight="bold", color="white",
                    )

            # Sample size labels — y in AXES fraction via get_xaxis_transform()
            # so they never push the data-limit downward
            for xi, n in enumerate(ns):
                ax.text(
                    xi, -0.08, f"n={n}",
                    ha="center", va="top",
                    transform=ax.get_xaxis_transform(),
                    fontsize=7.5, color="#555555",
                )

            # OR + p annotation as an in-axes text box (axes coordinates 0-1)
            # — never extends outside the axes, so bbox_inches='tight' is safe
            p_disp = p_a if pd.notna(p_a) else p_u
            suffix = "age-adj" if pd.notna(p_a) else "unadj"
            sig    = _sig_label(p_disp)

            ann_lines = [f"{sig}  p={p_disp:.3f} ({suffix})"]
            if pd.notna(or_v):
                ci_str = row.get("OR_Adj_Age_CI95", "")
                ann_lines.insert(0, f"OR={or_v:.2f}  {ci_str}".strip())

            ax.text(
                0.5, 0.97,
                "\n".join(ann_lines),
                ha="center", va="top",
                transform=ax.transAxes,
                fontsize=7.5, color="#333333",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor="#cccccc", alpha=0.85),
            )

            # Fixed y-axis 0-100 — no text escapes this range
            ax.set_ylim(0, 100)
            ax.set_xticks(xpos)
            ax.set_xticklabels(groups, fontsize=9)
            ax.set_ylabel("% positive", fontsize=9)
            ax.yaxis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f"{v:.0f}%")
            )
            gene  = row.get("Gene", "") or ""
            snp   = var_id if var_id != gene else var_id
            title = (f"{gene}  [{snp}]\n{row['Clin_Label']}"
                     if gene and gene != snp else f"{snp}\n{row['Clin_Label']}")
            ax.set_title(f"Carrier-Status Binary Outcome Plot\n{title}", fontsize=8.5, fontweight="bold")
            _style_ax(ax)

        # Hide unused subplots
        for ax in axes_flat[n_plots:]:
            ax.set_visible(False)

        # suptitle at y=0.98 (inside the figure) — prevents bbox expansion
        fig.suptitle(
            f"{label} — % positive: carriers vs non-carriers",
            fontsize=12, fontweight="bold",
            y=0.98,
        )
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        out_path = out_dir / f"17_BinaryBar_{label}.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.1)
        plt.close()
        print(f"  Saved: {out_path}")




# ── 5b. VIOLIN STRIP — continuous outcomes by genotype (WT / Het / Hom) ───────

_GENO_PLOT_C = GENOTYPE_COLORS.copy()
_GENO_PLOT_HATCH = {"WT": "", "Het": "//", "Hom": "xx"}
_GENO_PLOT_TEXT = {"WT": "#222222", "Het": "#222222", "Hom": "white"}
_GENO_PLOT_MARKER = {"WT": "o", "Het": "s", "Hom": "^"}

def make_genotype_distribution_plots(breast_clin, endo_clin, merged, out_dir):
    """
    For the top N nominally significant continuous associations, draw a
    violin + strip plot with three groups: WT (0/0), Het (0/1), Hom (1/1).
    """
    TOP_N = 12
    GT_MAP = {"0/0": "WT", "0/1": "Het", "1/0": "Het", "1/1": "Hom"}
    GENO_ORDER = ["WT", "Het", "Hom"]

    for clin_res, label, coh_sheet in [
        (breast_clin, "Breast",      "MT-T_N"),
        (endo_clin,   "Endometrial", "AT=AUs"),
    ]:
        if clin_res.empty:
            continue

        cont = (
            clin_res[
                (clin_res["Clin_Type"] == "continuous") &
                clin_res["P_Unadj"].notna() &
                (clin_res["P_Unadj"] < 0.05)
            ]
            .sort_values("P_Unadj")
            .drop_duplicates("Clinical_Var")
            .head(TOP_N)
        )
        if cont.empty:
            continue

        if "GT" not in merged.columns:
            print(f"  Genotype distribution plots ({label}): no GT column — skipped")
            continue

        n_plots = len(cont)
        ncols   = min(3, n_plots)
        nrows   = int(np.ceil(n_plots / ncols))
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(5.5 * ncols, 4.5 * nrows),
            squeeze=False,
        )
        axes_flat = axes.flatten()

        cohort_df = merged[
            (merged["sheet"] == coh_sheet) & (~merged["is_replicate"])
        ].drop_duplicates("Sample").copy()

        for idx, (_, row) in enumerate(cont.iterrows()):
            ax     = axes_flat[idx]
            col    = row["Clinical_Var"]
            var_id = row["Variant_ID"]
            p_u    = row["P_Unadj"]
            p_a    = row.get("P_Adj_Age", np.nan)

            gt_lookup = (
                merged[
                    (merged["Variant_ID"] == var_id) &
                    (merged["sheet"] == coh_sheet)
                ]
                .drop_duplicates("Sample")
                .set_index("Sample")["GT"]
                .apply(lambda g: GT_MAP.get(str(g), None))
            )
            cohort_df["_geno"] = cohort_df["Sample"].map(gt_lookup).fillna("WT")

            plot_df = cohort_df[[col, "_geno"]].dropna()
            plot_df = plot_df[plot_df["_geno"].isin(GENO_ORDER)]

            present = [g for g in GENO_ORDER
                       if (plot_df["_geno"] == g).sum() >= MIN_CARRIERS]
            if len(present) < 2:
                ax.set_visible(False)
                continue

            col_numeric = pd.to_numeric(plot_df[col], errors="coerce")
            p01, p99 = col_numeric.quantile(0.01), col_numeric.quantile(0.99)
            if p99 > p01:
                plot_df = plot_df.copy()
                plot_df[col] = col_numeric.clip(lower=p01, upper=p99)

            parts = ax.violinplot(
                [plot_df.loc[plot_df["_geno"] == g, col].values for g in present],
                positions=list(range(len(present))),
                widths=0.6, showmedians=False,
            )
            for pc, g in zip(parts["bodies"], present):
                pc.set_facecolor(_GENO_PLOT_C[g])
                pc.set_alpha(0.35)
                pc.set_edgecolor("#4A4A4A")
                pc.set_linewidth(0.9)
                pc.set_hatch(_GENO_PLOT_HATCH[g])
            for comp in ("cbars", "cmins", "cmaxes"):
                if comp in parts:
                    parts[comp].set_color("#aaaaaa")
                    parts[comp].set_linewidth(0.8)

            rng = np.random.default_rng(42)
            for gi, g in enumerate(present):
                vals = plot_df.loc[plot_df["_geno"] == g, col].values
                jit  = rng.uniform(-0.08, 0.08, len(vals))
                ax.scatter(
                    np.full(len(vals), gi) + jit, vals,
                    color=_GENO_PLOT_C[g], alpha=0.75, s=34,
                    marker=_GENO_PLOT_MARKER[g],
                    edgecolors="#4A4A4A", linewidths=0.4, zorder=3,
                )

            for gi, g in enumerate(present):
                med = plot_df.loc[plot_df["_geno"] == g, col].median()
                ax.hlines(med, gi - 0.18, gi + 0.18,
                          colors="#222222", linewidths=2, zorder=4)

            kw_groups = [
                pd.to_numeric(plot_df.loc[plot_df["_geno"] == g, col],
                              errors="coerce").dropna().values
                for g in present
            ]
            try:
                # Use Spearman trend (consistent with Analysis 3 dose test)
                dose_plot = np.concatenate([
                    np.full(len(g), gi) for gi, g in enumerate(kw_groups)
                ])
                vals_plot = np.concatenate(kw_groups)
                sp_r, sp_p = stats.spearmanr(dose_plot, vals_plot)
                kw_lbl = f"Spearman ρ={sp_r:.2f}, p={sp_p:.3f} {_sig_label(sp_p)}"
            except Exception:
                sp_p, kw_lbl = np.nan, ""

            ax.text(
                0.03, 0.97,
                kw_lbl,
                ha="left", va="top",
                transform=ax.transAxes,
                fontsize=7.4, color="#333333",
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#d0d0d0", alpha=0.85),
            )

            for gi, g in enumerate(present):
                n = (plot_df["_geno"] == g).sum()
                ax.text(
                    gi, -0.08, f"n={n}",
                    ha="center", va="top",
                    transform=ax.get_xaxis_transform(),
                    fontsize=7.5, color="#555555",
                )

            ax.set_xticks(list(range(len(present))))
            ax.set_xticklabels(present, fontsize=9)
            ax.set_ylabel(_wrap_label(row["Clin_Label"], width=22, max_lines=2), fontsize=9)
            gene  = row.get("Gene", "") or ""
            title = f"{gene} [{var_id}]" if gene and gene != var_id else var_id
            ax.set_title(_wrap_label(title, width=28, max_lines=2), fontsize=9, fontweight="bold", pad=8)
            _style_ax(ax)

        for ax in axes_flat[n_plots:]:
            ax.set_visible(False)

        fig.suptitle(
            f"Continuous Clinical Variables According to GSDMB SNP Genotype in {label}",
            fontsize=12,
            fontweight="bold",
            y=0.985,
        )
        plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=2.0, w_pad=1.4)
        out_path = out_dir / f"17_GenoDistribution_{label}.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {out_path}")


# ── 5c. STACKED BAR — binary outcomes by genotype (WT / Het / Hom) ────────────

def make_genotype_binary_bar_plots(breast_clin, endo_clin, merged, out_dir):
    """
    For the top N binary associations, show stacked proportion bars
    (% positive) for WT, Het, and Hom separately.
    """
    TOP_N = 12
    GT_MAP = {"0/0": "WT", "0/1": "Het", "1/0": "Het", "1/1": "Hom"}
    GENO_ORDER = ["WT", "Het", "Hom"]

    for clin_res, label, coh_sheet in [
        (breast_clin, "Breast",      "MT-T_N"),
        (endo_clin,   "Endometrial", "AT=AUs"),
    ]:
        if clin_res.empty:
            continue

        sig_clin = _filter_significant_clinical_rows(clin_res)
        bin_hits = (
            sig_clin[
                (sig_clin["Clin_Type"] == "binary") &
                sig_clin["P_Unadj"].notna()
            ]
            .sort_values("P_Unadj")
            .drop_duplicates(["Clinical_Var", "Contrast"] if "Contrast" in clin_res.columns else ["Clinical_Var"])
            .head(TOP_N)
        )
        if bin_hits.empty:
            continue

        if "GT" not in merged.columns:
            print(f"  Genotype binary bar plots ({label}): no GT column — skipped")
            continue

        n_plots = len(bin_hits)
        ncols   = min(3, n_plots)
        nrows   = int(np.ceil(n_plots / ncols))

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(4.5 * ncols, 3.8 * nrows),
            squeeze=False,
            facecolor="white",
        )
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

            gt_lookup = (
                merged[
                    (merged["Variant_ID"] == var_id) &
                    (merged["sheet"] == coh_sheet)
                ]
                .drop_duplicates("Sample")
                .set_index("Sample")["GT"]
                .apply(lambda g: GT_MAP.get(str(g), None))
            )
            cohort_df["_geno"] = cohort_df["Sample"].map(gt_lookup).fillna("WT")

            plot_df = cohort_df[[col, "_geno"]].dropna()
            plot_df = plot_df[plot_df["_geno"].isin(GENO_ORDER)]

            present = [g for g in GENO_ORDER
                       if (plot_df["_geno"] == g).sum() >= MIN_CARRIERS]
            if len(present) < 2:
                ax.set_visible(False)
                continue

            pos_rates = []
            ns        = []
            for g in present:
                vals = plot_df.loc[plot_df["_geno"] == g, col]
                vals_num = pd.to_numeric(vals, errors="coerce").dropna()
                pos_rates.append(100 * vals_num.mean() if len(vals_num) else 0)
                ns.append(len(vals_num))

            colours = [_GENO_PLOT_C[g] for g in present]
            xpos    = list(range(len(present)))

            for xi, (r, c, g) in enumerate(zip(pos_rates, colours, present)):
                ax.bar(xi, r, color=c,
                       width=0.55, edgecolor="#4A4A4A", linewidth=1.0,
                       hatch=_GENO_PLOT_HATCH[g], alpha=0.9, zorder=3)
                ax.bar(xi, 100 - r, bottom=r, color=c,
                       width=0.55, alpha=0.16, edgecolor="none", zorder=2)

            for xi, (r, g) in enumerate(zip(pos_rates, present)):
                if r > 8:
                    ax.text(xi, r / 2, f"{g}\n{r:.0f}%",
                            ha="center", va="center",
                            fontsize=8.5, fontweight="bold", color=_GENO_PLOT_TEXT[g])

            for xi, n in enumerate(ns):
                ax.text(xi, -0.08, f"n={n}",
                        ha="center", va="top",
                        transform=ax.get_xaxis_transform(),
                        fontsize=7.5, color="#555555")

            try:
                ct_vals = []
                for g in present:
                    vals_num = pd.to_numeric(
                        plot_df.loc[plot_df["_geno"] == g, col], errors="coerce"
                    ).dropna()
                    ct_vals.append([int((vals_num == 1).sum()),
                                    int((vals_num == 0).sum())])
                ct = np.array(ct_vals)
                if ct.shape[0] >= 2 and ct.shape[1] >= 2 and ct.min() >= 0:
                    if ct.shape == (2, 2):
                        _, chi_p = fisher_exact(ct)
                    else:
                        _, chi_p, _, _ = stats.chi2_contingency(ct)
                else:
                    chi_p = np.nan
            except Exception:
                chi_p = np.nan

            p_disp = p_a if pd.notna(p_a) else p_u
            suffix = "age-adj" if pd.notna(p_a) else "unadj"
            sig    = _sig_label(p_disp)
            ann_lines = [f"{sig}  p={p_disp:.3f} ({suffix})"]
            if pd.notna(chi_p):
                ann_lines.append(f"Geno χ² p={chi_p:.3f} {_sig_label(chi_p)}")

            ax.text(
                0.5, 0.97,
                "\n".join(ann_lines),
                ha="center", va="top",
                transform=ax.transAxes,
                fontsize=7.5, color="#333333",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor="#cccccc", alpha=0.85),
            )

            ax.set_ylim(0, 100)
            ax.set_xticks(xpos)
            ax.set_xticklabels(present, fontsize=9)
            ax.set_ylabel("% positive", fontsize=9)
            ax.yaxis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f"{v:.0f}%")
            )
            gene  = row.get("Gene", "") or ""
            snp_label = f"{gene}  [{var_id}]" if gene and gene != var_id else var_id
            ax.set_title(
                f"{snp_label}\nBinary Outcome: {row['Clin_Label']}",
                fontsize=8.5,
                fontweight="bold",
            )
            _style_ax(ax)

        for ax in axes_flat[n_plots:]:
            ax.set_visible(False)

        fig.suptitle(
            f"Binary Clinical Outcomes According to GSDMB SNP Genotype in {label}",
            fontsize=12,
            fontweight="bold",
            y=0.98,
        )
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        out_path = out_dir / f"17_GenoBinaryBar_{label}.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.1)
        plt.close()
        print(f"  Saved: {out_path}")


# ── 6. KAPLAN-MEIER (improved) ────────────────────────────────────────────────

# ?? 5d. GENOTYPE COMPOSITION ? case/control comparison using WT / Het / Hom ??

def make_genotype_composition_plots(tvh, merged, out_dir):
    """
    For the top SNPs in each tumour-vs-control comparison, show the sample
    composition of WT, heterozygous, and homozygous genotypes within each arm.
    Bars are normalised to percentages so the visual comparison is not biased by
    unequal sample sizes.
    """
    if tvh.empty:
        return
    if "GT" not in merged.columns:
        print("  Genotype composition plots: no GT column ? skipped")
        return

    TOP_N = 8
    GT_MAP = {"0/0": "WT", "0/1": "Het", "1/0": "Het", "1/1": "Hom"}
    GENO_ORDER = ["WT", "Het", "Hom"]

    sample_manifest = (
        merged[~merged["is_replicate"]]
        .drop_duplicates("Sample")
        .copy()
    )
    pooled_controls = sample_manifest[sample_manifest["Tissue"] == "Healthy"].copy()

    comparison_defs = {
        "Breast": {
            "tumour_mask": (sample_manifest["Tissue"] == "Tumour")
                           & sample_manifest["Cohort"].astype(str).str.contains("Breast", case=False, na=False),
            "tumour_label": "Breast tumour",
            "control_label": "Pooled control",
        },
        "Endometrium": {
            "tumour_mask": (sample_manifest["Tissue"] == "Tumour")
                           & sample_manifest["Cohort"].astype(str).str.contains("Endometri", case=False, na=False),
            "tumour_label": "Endometrium tumour",
            "control_label": "Pooled control",
        },
        "Global": {
            "tumour_mask": sample_manifest["Tissue"] == "Tumour",
            "tumour_label": "All tumours",
            "control_label": "All controls",
        },
    }

    for analysis_group, cfg in comparison_defs.items():
        sub = (
            tvh[tvh["Analysis_Group"] == analysis_group]
            .sort_values(["FDR_Sig", "Nominal_Sig", "P_Value"], ascending=[False, False, True])
            .drop_duplicates("Variant_ID")
            .copy()
        )
        if sub.empty:
            continue

        if sub["Nominal_Sig"].any() or sub["FDR_Sig"].any():
            chosen = sub[sub["Nominal_Sig"] | sub["FDR_Sig"]].head(TOP_N).copy()
        else:
            chosen = sub.head(TOP_N).copy()
        if chosen.empty:
            continue

        tumour_samples = sample_manifest.loc[cfg["tumour_mask"]].copy()
        control_samples = pooled_controls.copy()
        if tumour_samples.empty or control_samples.empty:
            continue

        arm_df = pd.concat(
            [
                control_samples.assign(_arm=cfg["control_label"]),
                tumour_samples.assign(_arm=cfg["tumour_label"]),
            ],
            ignore_index=True,
        )[["Sample", "_arm"]]

        n_plots = len(chosen)
        ncols = min(3, n_plots)
        nrows = int(np.ceil(n_plots / ncols))
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(5.0 * ncols, 4.2 * nrows),
            squeeze=False,
            facecolor="white",
        )
        axes_flat = axes.flatten()

        for ax in axes_flat:
            ax.set_facecolor("white")

        for idx, (_, row) in enumerate(chosen.iterrows()):
            ax = axes_flat[idx]
            var_id = row["Variant_ID"]
            gene = row.get("Gene", "") or ""

            gt_lookup = (
                merged[(merged["Variant_ID"] == var_id) & (~merged["is_replicate"])]
                .drop_duplicates("Sample")
                .set_index("Sample")["GT"]
                .apply(lambda g: GT_MAP.get(str(g), "WT"))
            )

            plot_df = arm_df.copy()
            plot_df["Genotype"] = plot_df["Sample"].map(gt_lookup).fillna("WT")

            counts = (
                plot_df.groupby(["_arm", "Genotype"])
                .size()
                .unstack(fill_value=0)
                .reindex(index=[cfg["control_label"], cfg["tumour_label"]], columns=GENO_ORDER, fill_value=0)
            )
            totals = counts.sum(axis=1)
            if (totals == 0).any():
                ax.set_visible(False)
                continue
            perc = counts.div(totals, axis=0) * 100

            xpos = np.arange(len(counts.index))
            bottoms = np.zeros(len(counts.index))
            for geno in GENO_ORDER:
                vals = perc[geno].to_numpy(dtype=float)
                ax.bar(
                    xpos,
                    vals,
                    bottom=bottoms,
                    color=_GENO_PLOT_C[geno],
                    edgecolor="#4A4A4A",
                    linewidth=0.9,
                    hatch=_GENO_PLOT_HATCH[geno],
                    width=0.58,
                    label=geno,
                )
                for xi, val, bottom in zip(xpos, vals, bottoms):
                    if val >= 9:
                        ax.text(
                            xi,
                            bottom + val / 2,
                            f"{geno}\n{val:.0f}%",
                            ha="center",
                            va="center",
                            fontsize=8.2,
                            color=_GENO_PLOT_TEXT[geno],
                            fontweight="bold",
                        )
                bottoms += vals

            ax.set_ylim(0, 100)
            ax.set_xticks(xpos)
            ax.set_xticklabels([
                f"{label}\n(n={int(totals.loc[label])})" for label in counts.index
            ], fontsize=8.5)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
            ax.set_ylabel("Genotype composition (%)", fontsize=9)
            snp_label = f"{gene} [{var_id}]" if gene and gene != var_id else var_id
            sig = "FDR" if bool(row.get("FDR_Sig", False)) else ("p<0.05" if bool(row.get("Nominal_Sig", False)) else "top hit")
            p_val = row.get("FDR_P_Value") if bool(row.get("FDR_Sig", False)) else row.get("P_Value")
            or_text = row.get("Odds_Ratio", np.nan)
            ax.set_title(
                f"Genotype Composition for {snp_label}\nOR={or_text:.2f} | p={p_val:.3g} ({sig})",
                fontsize=8.5,
                fontweight="bold",
            )
            _style_ax(ax)

        for ax in axes_flat[n_plots:]:
            ax.set_visible(False)

        handles = [
            plt.Rectangle((0, 0), 1, 1, facecolor=_GENO_PLOT_C[g], edgecolor="#4A4A4A", hatch=_GENO_PLOT_HATCH[g])
            for g in GENO_ORDER
        ]
        fig.legend(handles, ["WT", "Het", "Hom"], title="Genotype state (WT / Het / Hom)",
                   loc="upper center", ncol=3, frameon=False,
                   bbox_to_anchor=(0.5, 1.02))
        fig.suptitle(
            tagged_title(f"Genotype Composition in Tumour and Pooled Control Samples: {analysis_group} Comparison", COMPARATIVE_TAG),
            fontsize=12,
            fontweight="bold",
            y=1.04,
        )
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        out_path = out_dir / f"17_GenotypeComposition_{analysis_group}.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {out_path}")


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

    out_pdf = out_dir / "GSDMB_SNP_KM_Curves_All_Cohorts.pdf"
    with pdf_backend.PdfPages(out_pdf) as pdf:
        for fig in km_pages:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"  Saved: {out_pdf}  ({len(km_pages)} pages)")


# ── 7. SUMMARY OVERVIEW PANEL ─────────────────────────────────────────────────

def make_summary_panel(tvh, breast_clin, endo_clin, surv_res, out_dir):
    """
    A single summary figure with 4 panels:
      A) SNP frequency in tumour vs control (lollipop)
      B) Nominal association rate per variant (bar chart)
      C) -log10 p-value overview across all tests (dot plot)
      D) Survival HR overview
    """
    fig = plt.figure(figsize=(19.5, 12.8))
    gs  = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.42)
    ax_freq = fig.add_subplot(gs[0, 0])
    ax_hits = fig.add_subplot(gs[0, 1])
    ax_dots = fig.add_subplot(gs[1, 0])
    ax_hr   = fig.add_subplot(gs[1, 1])

    # ── Panel A: SNP frequency in tumour vs control ───────────────────────
    if not tvh.empty:
        # Use actual Analysis_Group values from tumour_vs_control()
        grp_colour_map = {"Breast": cohort_color("Breast"), "Endometrium": cohort_color("Endometrium")}
        plotted_any = False
        for grp, col in grp_colour_map.items():
            sub = tvh[tvh["Analysis_Group"] == grp].copy()
            if sub.empty: continue
            sub = sub.sort_values("Freq_Tumour_%", ascending=False).head(20)
            y = np.arange(len(sub))
            ax_freq.hlines(y, sub["Freq_Control_%"].values,
                           sub["Freq_Tumour_%"].values,
                           colors=col, linewidth=1.2, alpha=0.6)
            ax_freq.scatter(sub["Freq_Tumour_%"].values, y,
                            color=col, s=60, zorder=3, label=grp)
            ax_freq.scatter(sub["Freq_Control_%"].values, y,
                            color=col, s=40, marker="D", alpha=0.5, zorder=3)
            plotted_any = True
        if plotted_any:
            ax_freq.set_yticks([])
            ax_freq.set_xlabel("Frequency (%)", fontsize=9)
            ax_freq.set_xlim(0, 100)
            ax_freq.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
            ax_freq.set_title(
                "A  Distribution of GSDMB SNP Frequencies in Tumour and Pooled Control Samples\n(top 20 variants per cohort)",
                fontsize=9,
                fontweight="bold",
                loc="left",
            )
            ax_freq.legend(fontsize=7.5, frameon=False)
            _style_ax(ax_freq)
        else:
            ax_freq.text(0.5, 0.5, "No cohort-level tumour vs control data",
                         ha="center", va="center", transform=ax_freq.transAxes,
                         color="#aaaaaa", fontsize=9)
            ax_freq.set_axis_off()

        # ── Panel B: Number of nominal hits per variant ───────────────────────
    all_clin = pd.concat(
        [x for x in [breast_clin, endo_clin] if not x.empty],
        ignore_index=True
    )
    if not all_clin.empty:
        test_counts = (
            all_clin.groupby(["Variant_ID", "Cohort"])
            .size().reset_index(name="N_tests")
        )
        hit_counts = (
            all_clin[all_clin["Nominal_Sig_Unadj"]]
            .groupby(["Variant_ID", "Cohort"])
            .size().reset_index(name="N_hits")
        )
        hit_rates = test_counts.merge(hit_counts, on=["Variant_ID", "Cohort"], how="left")
        hit_rates["N_hits"] = hit_rates["N_hits"].fillna(0)
        hit_rates = hit_rates[hit_rates["N_hits"] > 0].copy()
        if not hit_rates.empty:
            hit_rates["Hit_Rate_%"] = hit_rates["N_hits"] / hit_rates["N_tests"] * 100
            hit_rates = hit_rates.sort_values(["Hit_Rate_%", "N_hits"], ascending=False).head(14).sort_values(["Hit_Rate_%", "N_hits"], ascending=True)
            colours = hit_rates["Cohort"].map(
                lambda c: _COHORT_C.get(
                    "Breast" if "Breast" in c else "Endometrial", "#888888"
                )
            )
            ypos = np.arange(len(hit_rates))
            ax_hits.barh(ypos, hit_rates["Hit_Rate_%"].values,
                         color=colours, alpha=0.8, edgecolor="white",
                         height=0.6)
            ax_hits.set_yticks(ypos)
            labels = hit_rates["Variant_ID"] + "\n(" + hit_rates["Cohort"].str.replace("_Tumour", "") + ")"
            ax_hits.set_yticklabels([_wrap_label(v, width=20, max_lines=2) for v in labels.values], fontsize=7.1)
            ax_hits.set_xlabel("Nominal association rate (% of tests, p < 0.05)", fontsize=9)
            ax_hits.set_title("B  Rate of nominally significant SNP-clinical associations\n(top 14 cohort-variant pairs)",
                               fontsize=9, fontweight="bold", loc="left")
            for y, rate, hits, tests in zip(ypos, hit_rates["Hit_Rate_%"], hit_rates["N_hits"], hit_rates["N_tests"]):
                ax_hits.text(rate + 0.6, y, f"{int(hits)}/{int(tests)}", va="center", fontsize=7.2, color="#555555")
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
                       .head(18))
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
            [_wrap_label(f"{r['Gene'] or r['Variant_ID']} | {r['Clin_Label']}", width=18, max_lines=3) for _, r in top.iterrows()],
            rotation=45, ha="right", fontsize=6.5
        )
        ax_dots.set_ylabel("-log₁₀(p-value)", fontsize=9)
        ax_dots.set_title(
            "C  Most Significant Unadjusted Associations Between GSDMB SNPs and Clinical Variables\n(top 18 variant-variable pairs)",
            fontsize=9,
            fontweight="bold",
            loc="left",
        )
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
            sr_plot["_label"] = (sr_plot["Gene"].fillna("") + " " + sr_plot["Variant_ID"]
                            + "\n(" + sr_plot["Cohort"].str.replace("_", " ")
                            + ", " + sr_plot["Endpoint"] + ")")
            sr_plot = sr_plot.sort_values([p_col, hr_col]).head(18).sort_values(hr_col)
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
            ax_hr.set_yticklabels([_wrap_label(v, width=24, max_lines=2) for v in sr_plot["_label"].values], fontsize=6.8)
            ax_hr.set_xlabel("Hazard Ratio (95% CI) — additive Cox", fontsize=9)
            note = f" ({n_unstable} unstable models excluded)" if n_unstable else ""
            ax_hr.set_title(
                f"D  Survival Associations of GSDMB SNPs Across Cohorts\n(top 18 stable models; * denotes p<0.05; colour = cohort){note}",
                fontsize=9,
                fontweight="bold",
                loc="left",
            )
            _style_ax(ax_hr, grid=False)
            ax_hr.xaxis.grid(True, linestyle=":", alpha=0.4)
    else:
        ax_hr.text(0.5, 0.5, "No survival data",
                   ha="center", va="center",
                   transform=ax_hr.transAxes, color="#aaaaaa", fontsize=11)
        ax_hr.set_axis_off()

    fig.suptitle(
        "Overview of GSDMB SNP Association Analyses",
        fontsize=14,
        fontweight="bold",
        y=1.01,
    )
    out_path = out_dir / "GSDMB_SNP_Association_Summary_Panel.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


# ── MAIN ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="SNP–clinical association analysis (script 17 v3)")
    p.add_argument("--gsdmb",   default=str(DEFAULT_GSDMB))
    p.add_argument("--master",  default=str(DEFAULT_MASTER))
    p.add_argument("--out_dir", default=str(DEFAULT_OUT))
    p.add_argument("--variant_whitelist", default=(str(DEFAULT_VARIANT_WHITELIST) if DEFAULT_VARIANT_WHITELIST else None),
                   help="Excel whitelist from script 11 defining the thesis SNP backbone")
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
    carrier frequency in tumour cases vs pooled healthy controls (breast + endometrial).
    Reports unadjusted OR (Fisher) and age- / age+BMI-adjusted OR
    (logistic regression), plus genotypic Het-vs-WT and Hom-vs-WT contrasts.
    """
    try:
        from statsmodels.formula.api import logit as sm_logit
    except ImportError:
        sm_logit = None
        print("  WARNING: statsmodels not available — adjusted models skipped.")

    print("=== Analysis 5: Cancer Risk (Case-Control) ===")
    df = df[~df["is_replicate"]].copy()

    rows = []

    pooled_controls = df[df["Tissue"] == "Healthy"].copy()

    for cohort_label, cfg in RISK_COHORTS.items():
        case_df    = df[df["sheet"] == cfg["case_sheet"]].copy()
        control_df = pooled_controls.copy()

        # One row per sample — use sample-level data (not variant-level rows)
        case_samples    = case_df.drop_duplicates("Sample").set_index("Sample")
        control_samples = control_df.drop_duplicates("Sample").set_index("Sample")

        n_cases    = len(case_samples)
        n_controls = len(control_samples)

        print(f"  {cohort_label}: {n_cases} cases, {n_controls} pooled controls")
        if n_cases < MIN_CARRIERS or n_controls < MIN_CARRIERS:
            print(f"    Too few samples — skipping.")
            continue

        # All variants present in either the tumour cohort or the pooled healthy arm
        all_variants = pd.concat([
            case_df[["Variant_ID"]],
            control_df[["Variant_ID"]],
        ], ignore_index=True)["Variant_ID"].dropna().unique()

        for var_id in all_variants:
            # Carrier sets
            case_carriers    = set(case_df[case_df["Variant_ID"] == var_id]["Sample"].unique())
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
            all_samples = pd.concat([
                case_samples[["canon__age", "canon__bmi"]].assign(
                    cancer=1,
                    carrier=lambda x: x.index.map(lambda s: 1 if s in case_carriers else 0)
                ),
                control_samples[["canon__age", "canon__bmi"]].assign(
                    cancer=0,
                    carrier=lambda x: x.index.map(lambda s: 1 if s in control_carriers else 0)
                ),
            ])
            all_samples["canon__age"] = pd.to_numeric(all_samples["canon__age"], errors="coerce")
            all_samples["canon__bmi"] = pd.to_numeric(all_samples["canon__bmi"], errors="coerce")

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
                    or_v = round(float(np.exp(b)), 4)
                    ci_lo = round(float(np.exp(b - 1.96 * se_b)), 4)
                    ci_hi = round(float(np.exp(b + 1.96 * se_b)), 4)

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
                        or_v = round(float(np.exp(fit.params["carrier"])), 4)
                        ci   = np.exp(fit.conf_int().loc["carrier"])
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
            age_df["age_z"] = _zscore(age_df["canon__age"])
            row.update(_fit_logit(age_df[["age_z", "carrier"]].assign(cancer=age_df["cancer"]),
                                  ["carrier", "age_z"], "Age"))

            # Age+BMI model
            bmi_df = all_samples[["cancer", "canon__age", "canon__bmi", "carrier"]].dropna().copy()
            bmi_df["age_z"] = _zscore(bmi_df["canon__age"])
            bmi_df["bmi_z"] = _zscore(bmi_df["canon__bmi"])
            row.update(_fit_logit(bmi_df[["age_z", "bmi_z", "carrier"]].assign(cancer=bmi_df["cancer"]),
                                  ["carrier", "age_z", "bmi_z"], "AgeBMI"))

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

        ax.set_xlim(0, OR_CAP * 1.10)
        for i, (_, row) in enumerate(sub.iterrows()):
            p_str = f"p={row['P_Adj_Age']:.3f}{_sig_label(row['P_Adj_Age'])}"
            n_str = f"n={int(row['N_Adj_Age'])}" if pd.notna(row.get("N_Adj_Age")) else ""
            label_txt = f"{p_str}  {n_str}".strip()
            ax.text(1.02, ypos[i], label_txt, transform=ax.get_yaxis_transform(),
                    va="center", ha="left", fontsize=6.8, color="#444444", clip_on=False)

        ax.axvline(1, color="#555555", linestyle="--", linewidth=1, alpha=0.7)
        labels = (sub["Gene"].fillna("") + " " + sub["Variant_ID"]).str.strip()
        ax.set_yticks(ypos)
        ax.set_yticklabels([_wrap_label(v, width=30, max_lines=2) for v in labels.values], fontsize=8)
        ax.set_xlabel("Odds Ratio for cancer risk\n(age-adjusted, 95% CI)", fontsize=10)
        ax.set_title(
            f"Age-Adjusted Association Between GSDMB SNP Carrier Status and {cohort} Cancer Risk\n"
            f"(cases n={risk_res.loc[risk_res['Cohort']==cohort, 'N_Cases'].iloc[0]}, "
            f"controls n={risk_res.loc[risk_res['Cohort']==cohort, 'N_Controls'].iloc[0]})",
            fontsize=11,
            fontweight="bold",
            pad=10,
        )

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

    fig.suptitle(
        "Age-Adjusted Associations Between GSDMB SNPs and Cancer Risk",
        fontsize=13,
        fontweight="bold",
        y=1.01,
    )
    plt.tight_layout(rect=[0, 0, 0.9, 1])
    out_path = out_dir / "17_Forest_CancerRisk.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def main():
    args    = parse_args()
    gsdmb   = Path(args.gsdmb)
    master  = Path(args.master)
    out_dir = Path(args.out_dir)
    validate_file_exists(gsdmb, "Script 17 GSDMB input")
    validate_file_exists(master, "Script 17 master input")
    variant_whitelist = None if args.variant_whitelist in {None, "", "None"} else Path(args.variant_whitelist)
    if variant_whitelist is not None:
        validate_file_exists(variant_whitelist, "Script 17 variant whitelist")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_xlsx = out_dir / "GSDMB_SNP_Clinical_Association_Results.xlsx"

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

    merged = load_and_merge(gsdmb, master, manifest_paths=manifest_paths, whitelist_path=variant_whitelist)
    print_validation_summary(merged, "Sample", "Script 17 merged analysis input", ["Cohort", "Tissue"])
    warning_lookup = build_amplicon_warning_lookup(merged, variant_col="Variant_ID", chrom_col="CHROM", pos_col="POS")

    # Analysis 1
    tvh = tumour_vs_control(merged)
    if not tvh.empty:
        validate_percentage_columns(tvh, ["Freq_Tumour_%", "Freq_Control_%"], "Script 17 tumour vs control")
    tvh = attach_amplicon_warning_columns(tvh, warning_lookup)

    # Analysis 2
    breast_clin, endo_clin = clinical_associations(merged)
    breast_clin = attach_amplicon_warning_columns(breast_clin, warning_lookup)
    endo_clin = attach_amplicon_warning_columns(endo_clin, warning_lookup)

    # Analysis 3 — genotype-dose
    # Run on ALL variants, not just nominal hits, to avoid selection bias.
    # FDR correction within genotype_dose_analysis handles the multiple-testing
    # burden; the FDR_Sig column in the output flags the variants that survive.
    all_breast_vars = breast_clin["Variant_ID"].unique().tolist() if not breast_clin.empty else []
    all_endo_vars   = endo_clin["Variant_ID"].unique().tolist()   if not endo_clin.empty   else []
    breast_dose = genotype_dose_analysis(merged, all_breast_vars, CLINICAL_VARS_BREAST, "Breast_Tumour")
    endo_dose   = genotype_dose_analysis(merged, all_endo_vars,   CLINICAL_VARS_ENDO,   "Endometrial_Tumour")
    breast_dose = attach_amplicon_warning_columns(breast_dose, warning_lookup)
    endo_dose = attach_amplicon_warning_columns(endo_dose, warning_lookup)

    # Analysis 4 — survival
    surv_res = pd.DataFrame()
    km_pages = []
    surv_out = survival_analysis(merged)
    if isinstance(surv_out, tuple):
        surv_res, km_pages = surv_out
    surv_res = attach_amplicon_warning_columns(surv_res, warning_lookup)

    # Analysis 5 — cancer risk (case-control)
    risk_res = cancer_risk_analysis(merged)
    risk_res = attach_amplicon_warning_columns(risk_res, warning_lookup)

    # Summary
    sig_parts = []
    if not tvh.empty:
        s = tvh[tvh["Nominal_Sig"]].copy()
        s["Analysis_Type"] = "Tumour_vs_Control"
        # Alias P_Value → P_Unadj so the shared sort key works across all result types
        s["P_Unadj"] = s["P_Value"]
        sig_parts.append(s[["Analysis_Type", "Analysis_Group", "Variant_ID", "Gene",
                             "Coverage_Risk_Flag", "Coverage_Risk_Amplicon",
                             "Coverage_Risk_Region", "Coverage_Risk_Note",
                             "Consequence", "IMPACT", "Freq_Tumour_%", "Freq_Control_%",
                             "Odds_Ratio", "P_Value", "P_Unadj", "FDR_P_Value",
                             "Nominal_Sig", "FDR_Sig"]])
    for res, label in [(breast_clin, "Breast_Tumour"), (endo_clin, "Endometrial_Tumour")]:
        if not res.empty:
            sig_mask = res["Nominal_Sig_Unadj"] | res["Nominal_Sig_Adj_Age"] | res["Nominal_Sig_Adj_BMI"] | res["Nominal_Sig_Adj_AgeBMI"]
            s = res[sig_mask].copy(); s["Analysis_Type"] = f"Clinical_{label}"
            keep_cols = ["Analysis_Type", "Cohort", "Variant_ID", "Gene",
                         "Coverage_Risk_Flag", "Coverage_Risk_Amplicon",
                         "Coverage_Risk_Region", "Coverage_Risk_Note",
                         "Contrast", "N_WT", "N_Het", "N_Hom",
                         "Clin_Label", "Clin_Type",
                         "Test_Unadj", "OR_Unadj", "P_Unadj", "FDR_Unadj",
                         "Nominal_Sig_Unadj", "FDR_Sig_Unadj",
                         "Test_Adj_Age", "OR_Adj_Age", "OR_Adj_Age_CI95", "P_Adj_Age", "FDR_Adj_Age",
                         "Nominal_Sig_Adj_Age", "FDR_Sig_Adj_Age", "Adj_Age_Note",
                         "Test_Adj_BMI", "OR_Adj_BMI", "OR_Adj_BMI_CI95", "P_Adj_BMI", "FDR_Adj_BMI",
                         "Nominal_Sig_Adj_BMI", "FDR_Sig_Adj_BMI", "Adj_BMI_Note",
                         "Test_Adj_AgeBMI", "OR_Adj_AgeBMI", "OR_Adj_AgeBMI_CI95", "P_Adj_AgeBMI", "FDR_Adj_AgeBMI",
                         "Nominal_Sig_Adj_AgeBMI", "FDR_Sig_Adj_AgeBMI", "Adj_AgeBMI_Note",
                         # Additive secondary model
                         "P_Unadj_Additive", "OR_Unadj_Additive",
                         "OR_Adj_Age_Additive", "P_Adj_Age_Additive",
                         "OR_Adj_BMI_Additive", "P_Adj_BMI_Additive",
                         "OR_Adj_AgeBMI_Additive", "P_Adj_AgeBMI_Additive"]
            s = s[[c for c in keep_cols if c in s.columns]]
            sig_parts.append(s)
    if not risk_res.empty:
        risk_sig = risk_res[risk_res["Nominal_Sig_Unadj"] | risk_res["Nominal_Sig_Adj_Age"]].copy()
        risk_sig["Analysis_Type"] = "Cancer_Risk"
        sig_parts.append(risk_sig)
    summary = (pd.concat(sig_parts, ignore_index=True).sort_values("P_Unadj")
               if sig_parts else pd.DataFrame({"Note": ["No nominally significant results."]}))
    summary = attach_amplicon_warning_columns(summary, warning_lookup)

    # Sample manifest
    keep_clin = ["BREAST_RECURRENCE_DERIVED", "BREAST_METASTASIS_DERIVED",
                 "BREAST_EXITUS_DERIVED", "BREAST_OS_MONTHS_DERIVED",
                 "ENDO_PD_BIN", "ENDO_M_STAGE_BIN", "ENDO_EXITUS_BIN",
                 "ENDO_EXITUS_DISEASE_BIN", "canon__os_months", "canon__pfs_months",
                 "ENDO_RISK_ORDINAL"]
    manifest = merged.drop_duplicates("Sample")[[
        c for c in ["Sample", "snp_code", "sample_id", "Cohort", "Tissue", "sheet",
                    "tumour_normal", "is_replicate"] + keep_clin
        if c in merged.columns
    ]].sort_values(["Cohort", "Tissue", "Sample"]).reset_index(drop=True)

    # Write Excel
    print("=== Writing output Excel ===")
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        if not tvh.empty:       tvh.to_excel(xw,         sheet_name="tumour_vs_control",     index=False)
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
    significant_dir = out_dir / "significant_only_figures"
    significant_dir.mkdir(parents=True, exist_ok=True)
    make_volcano_plots(tvh, out_dir)
    make_heatmaps(breast_clin, endo_clin, out_dir, significant_only=False)
    make_heatmaps(breast_clin, endo_clin, significant_dir, significant_only=True)
    make_main_text_heatmaps(breast_clin, endo_clin, out_dir)
    make_main_text_heatmaps(breast_clin, endo_clin, significant_dir)
    make_forest_plots(breast_clin, endo_clin, out_dir)
    make_distribution_plots(breast_clin, endo_clin, merged, out_dir)
    make_binary_bar_plots(breast_clin, endo_clin, merged, out_dir)
    make_km_plots(surv_res, km_pages, out_dir)
    make_risk_forest_plot(risk_res, out_dir)
    make_summary_panel(tvh, breast_clin, endo_clin, surv_res, out_dir)
    make_genotype_distribution_plots(breast_clin, endo_clin, merged, out_dir)
    make_genotype_binary_bar_plots(breast_clin, endo_clin, merged, out_dir)
    make_genotype_composition_plots(tvh, merged, out_dir)

    # Final summary
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"  Samples in manifest (sequenced): {merged[~merged['is_replicate']]['Sample'].nunique()}")
    print(f"  Variants tested  : {merged['Variant_ID'].nunique()}")
    if not tvh.empty:
        print(f"  Tumour vs control   — nominal: {tvh['Nominal_Sig'].sum()} | FDR<{FDR_THRESHOLD}: {tvh['FDR_Sig'].sum()}")
    if not breast_clin.empty:
        print(f"  Breast clinical     — nominal unadj: {breast_clin['Nominal_Sig_Unadj'].sum()} | FDR: {breast_clin['FDR_Sig_Unadj'].sum()} | nominal age-adj: {breast_clin['Nominal_Sig_Adj_Age'].sum()} | nominal BMI-adj: {breast_clin['Nominal_Sig_Adj_BMI'].sum()} | nominal age+BMI-adj: {breast_clin['Nominal_Sig_Adj_AgeBMI'].sum()}")
    if not endo_clin.empty:
        print(f"  Endometrial clinical— nominal unadj: {endo_clin['Nominal_Sig_Unadj'].sum()} | FDR: {endo_clin['FDR_Sig_Unadj'].sum()} | nominal age-adj: {endo_clin['Nominal_Sig_Adj_Age'].sum()} | nominal BMI-adj: {endo_clin['Nominal_Sig_Adj_BMI'].sum()} | nominal age+BMI-adj: {endo_clin['Nominal_Sig_Adj_AgeBMI'].sum()}")
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


