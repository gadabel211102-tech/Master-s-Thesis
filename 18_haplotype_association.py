#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
18_haplo_clinical_association.py
=======================================================================
Haplotype–Clinical Variable Association Analysis

Mirrors script 17 (SNP–clinical associations) but uses BEAGLE-phased
haplotype carrier status as the exposure variable instead of individual
SNP genotypes.

DATA SOURCES
------------
1. phased_genotypes.tsv  — BEAGLE output from script 17_haplotypes.sh
   Rows = SNPs, Cols = CHROM, POS, REF, ALT, <sample1>, <sample2>, ...
   GT format: "0|1", "1|0", "0|0", "1|1"

2. MASTER_SNP_plus_clinical__HARMONISED_B_v4.xlsx  — harmonised clinical master
   (same file used by script 17; provides all clinical variables and the
   snp_code → sample mapping)

3. GSDMB_Annotated_Report_Fixed.xlsx  — used only to identify the 14 established
   SNP positions (gnomAD NFE AF > 1%) that define the haplotype backbone

HAPLOTYPE CONSTRUCTION
----------------------
For each sample, a haplotype string is built from the phased alleles at the
established SNP positions (e.g. "01011001" = alt allele present at positions
2, 4, 5, 8).  Each sample has two haplotype strings (one per chromosome).
A sample is a "carrier" of a given haplotype if it appears on at least one
chromosome (dosage ≥ 1).

ANALYSES (parallel to script 17)
---------------------------------
1. Tumour vs Control       — Fisher's exact test per haplotype
2. Haplotype × Clinical    — Mann–Whitney U (continuous) or Fisher exact
                             (binary) or Chi-square (nominal), per tumour
                             cohort.  Age- and BMI-adjusted logistic
                             regression for binary outcomes.
3. Haplotype-dose          — Trend test across dosage 0/1/2 (Kruskal-Wallis
                             or chi-square) for nominally significant hits
4. Survival (Cox PH)       — Age-adjusted Cox for OS/PFS, KM curves per
                             haplotype carrier status
5. Cancer Risk             — Case-control logistic regression (cases = tumour,
                             controls = healthy)

OUTPUTS
-------
  18_Haplo_Clinical_Association_Results.xlsx
    • haplotype_frequencies
    • tumour_vs_control
    • breast_clinical_assoc
    • endo_clinical_assoc
    • haplotype_dose
    • survival_cox
    • cancer_risk
    • summary_significant
    • sample_manifest

  18_Haplo_Volcano_raw_p.png
  18_Haplo_Heatmap_Breast.png
  18_Haplo_Heatmap_Endometrial.png
  18_Haplo_Forest_Breast.png
  18_Haplo_Forest_Endometrial.png
  18_Haplo_Forest_CancerRisk.png
  18_KM_Curves_All_Cohorts.pdf

RUN
---
python3 18_haplo_clinical_association.py \\
  --phased   /path/to/19_haplotype_phased/phased_genotypes.tsv \\
  --master   /path/to/MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx \\
  --annot    /path/to/GSDMB_Annotated_Report_Fixed.xlsx \\
  --out_dir  /path/to/output/

NOTES
-----
• Only haplotypes with global frequency ≥ MIN_HAP_FREQ are tested.
• Sample names in the phased TSV must match snp_code values in the master.
  The script normalises both (strip spaces, uppercase) before joining.
• Replicates (is_replicate = True in the master) are excluded from all
  analyses, consistent with script 17.
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

from association_runtime import script18_defaults
from figure_style import COMPARATIVE_TAG, COHORT_COLORS, HAPLOTYPE_STATUS_COLORS, arm_color, cohort_color
from pipeline_validation import print_validation_summary, validate_file_exists, validate_percentage_columns

try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    _HAS_LIFELINES = True
except ImportError:
    _HAS_LIFELINES = False
    print("WARNING: lifelines not installed — survival plots will be skipped.")

warnings.filterwarnings("ignore")

# ── DEFAULT PATHS ─────────────────────────────────────────────────────────────
DEFAULT_PHASED  = Path("/home/gadeaalonsoj/tfm/gsdmb_final_results/19_haplotype_phased/phased_genotypes.tsv")
DEFAULT_MASTER  = Path("/home/gadeaalonsoj/tfm/MASTER_SNP_plus_clinical__HARMONISED_B_v4.xlsx")
DEFAULT_ANNOT   = Path("/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_Annotated_Report_Fixed.xlsx")
DEFAULT_OUT          = Path("/home/gadeaalonsoj/tfm/gsdmb_final_results/")
DEFAULT_HAPLO_RESULTS = Path("/home/gadeaalonsoj/tfm/gsdmb_final_results/19_haplo_stats_results/19_Haplotype_Results_v7_blocks_and_genes.xlsx")

# ── CONSTANTS ─────────────────────────────────────────────────────────────────
MIN_HAP_FREQ        = 0.02   # global frequency threshold ? rare haplotypes skipped
MIN_CARRIERS        = 5      # minimum carriers (or non-carriers) to run a test
MIN_COMPARISON_CARRIERS = 3  # minimum callable carriers within a specific case-control comparison
MIN_EVENTS_LOGISTIC = 10     # minimum events per predictor (EPV rule) for stable logistic regression
FDR_THRESHOLD       = 0.10

# ── CLINICAL VARIABLE DEFINITIONS ────────────────────────────────────────────
# Identical to script 17 — haplotype carrier status replaces SNP carrier status
# as the exposure, but the outcome variables and test types are unchanged.

CLINICAL_VARS_BREAST: Dict[str, Dict] = {
    # Prognostic/predictive variables — matched to script 17.
    # REMOVED: age, BMI, menarche, menopause (not tumour outcomes),
    #   OS_MONTHS as continuous (Cox in Analysis 4),
    #   BREAST_LOCAL_MET_BIN (broken) → replaced by ANY_METASTASIS_BIN.
    "BREAST_GRADE_NUMERIC":      {"type": "continuous", "label": "Tumour grade (1/2/3)", "note": "Ordinal 1/2/3", "bmi_adjust": True},
    "BREAST_KI67_NUMERIC":       {"type": "continuous", "label": "KI67 (proliferation index)", "bmi_adjust": False},
    "canon__her2_copies":        {"type": "continuous", "label": "HER2 copies (FISH)", "bmi_adjust": False},
    "BREAST_P53_NUMERIC":        {"type": "continuous", "label": "p53 expression (IHC)", "bmi_adjust": False},
    "BREAST_ER_BIN":             {"type": "binary",     "label": "ER positive", "note": "Positive=1 vs Negative=0", "bmi_adjust": False},
    "BREAST_PR_BIN":             {"type": "binary",     "label": "PR positive", "note": "Positive=1 vs Negative=0", "bmi_adjust": False},
    "BREAST_RECURRENCE_DERIVED": {"type": "binary",     "label": "Recurrence / progression", "note": "Any recurrence=1 vs NO=0", "bmi_adjust": True},
    "BREAST_ANY_METASTASIS_BIN": {"type": "binary",     "label": "Any metastasis",
                                   "note": "Local or distant=1 vs none=0", "bmi_adjust": True},
    "BREAST_EXITUS_DERIVED":     {"type": "binary",     "label": "Exitus", "note": "SI=1 vs NO=0", "bmi_adjust": True},
    "BREAST_HER2_SUBTYPE":       {"type": "nominal",    "label": "HER2 subtype", "note": "HER2+ / TN / Other — core GSDMB-relevant subtype", "bmi_adjust": False},
    "BREAST_DX_TYPE":            {"type": "nominal",    "label": "Histological diagnosis type", "note": "CDI / CDIS / other", "bmi_adjust": False},
}
CLINICAL_VARS_ENDO: Dict[str, Dict] = {
    # Prognostic/predictive variables — matched to script 17.
    # REMOVED: age, OS/PFS as continuous (Analysis 4 Cox),
    #   ctDNA/cfDNA + ITH (too sparse), M_STAGE (near-zero events),
    #   duplicate nominal FIGO and risk (ordinal numeric retained for power).
    # EXITUS: disease-specific is primary; all-cause supplementary (correlated).
    "ENDO_FIGO_NUMERIC":         {"type": "continuous", "label": "FIGO stage (ordinal)", "note": "IA=1, IB=1.5, II=2, IIIA=3.1 … IVB=4.2", "bmi_adjust": True},
    "ENDO_GRADE_NUMERIC":        {"type": "continuous", "label": "Tumour grade (G1/G2/G3)", "note": "G1=1, G2=2, G3=3", "bmi_adjust": True},
    "ENDO_RISK_ORDINAL":         {"type": "continuous", "label": "Risk of recurrence (ordinal)", "note": "LOW=1, INT=2, INT-HIGH=3, HIGH=4", "bmi_adjust": True},
    "ENDO_KI67_NUMERIC":         {"type": "continuous", "label": "KI67 (proliferation index)", "note": "Fraction 0–1", "bmi_adjust": False},
    "ENDO_PDL1_NUMERIC":         {"type": "continuous", "label": "PD-L1 expression (%)",
                                   "note": "GSDMB pyroptosis → immune activation → PD-L1", "bmi_adjust": False},
    "ENDO_CD8_NUMERIC":          {"type": "continuous", "label": "CD8+ TILs (%)",
                                   "note": "Immune infiltration — pyroptosis hypothesis", "bmi_adjust": False},
    "ENDO_EXITUS_DISEASE_BIN":   {"type": "binary",     "label": "Disease-specific exitus [PRIMARY]", "note": "YES=1 vs NO=0 — primary survival endpoint", "bmi_adjust": True},
    "ENDO_EXITUS_BIN":           {"type": "binary",     "label": "All-cause exitus [SUPPLEMENTARY]",
                                   "note": "Correlated with disease-specific; interpret together", "bmi_adjust": True},
    "ENDO_PD_BIN":               {"type": "binary",     "label": "Disease progression", "note": "PD=1 vs NO PD=0", "bmi_adjust": True},
    "ENDO_PTEN_BIN":             {"type": "binary",     "label": "PTEN loss/reduced", "note": "LOST/REDUCED=1 vs CONSERVED=0", "bmi_adjust": False},
    "ENDO_MLH1_BIN":             {"type": "binary",     "label": "MLH1 loss/reduced", "note": "LOST/REDUCED=1 vs CONSERVED=0 — mismatch repair marker", "bmi_adjust": False},
    "ENDO_N_STAGE_BIN":          {"type": "binary",     "label": "Lymph node involvement", "note": "N1/N2=1 vs N0=0", "bmi_adjust": True},
    "ENDO_LVSI_BIN":             {"type": "binary",     "label": "LVSI", "note": "Lymphovascular space invasion — YES=1 vs NO=0", "bmi_adjust": False},
    "ENDO_MYOINV_BIN":           {"type": "binary", "label": "Myometrial invasion ≥50%",
                                   "note": ">50%=1 vs <50%=0", "bmi_adjust": True},
    "ENDO_MSI_BIN":              {"type": "binary",     "label": "MSI-H", "note": "Unstable=1 vs Stable=0", "bmi_adjust": False},
    "ENDO_NEEC_BIN":             {"type": "binary",     "label": "Non-endometrioid histology", "note": "NEEC=1 vs EEC=0", "bmi_adjust": False},
    "ENDO_ER_BIN":               {"type": "binary",     "label": "ER positive", "note": "Positive=1 vs Negative=0", "bmi_adjust": False},
    "ENDO_PR_BIN":               {"type": "binary",     "label": "PR positive", "note": "Positive=1 vs Negative=0", "bmi_adjust": False},
    "ENDO_TP53_ABN_BIN":         {"type": "binary",     "label": "TP53 IHC abnormal", "note": "Aberrant=1 vs WT=0 — defines p53-abn molecular subtype", "bmi_adjust": False},
    "ENDO_GENE_AMP_BIN":         {"type": "binary",     "label": "Gene amplification", "note": "YES=1 vs NO=0 — relevant to GSDMB locus amplification", "bmi_adjust": False},
    "canon__molecular_class":    {"type": "nominal",    "label": "Molecular classification",
                                   "note": "POLE / MMRd / NSMP / P53", "bmi_adjust": False},
}
SURVIVAL_COHORTS = {
    "Endometrial": {
        "cohort_filter": "Endometri",
        "endpoints": {
            "OS":  {"t_col": "canon__os_months",        "ev_col": "ENDO_EXITUS_DISEASE_BIN"},  # disease-specific is primary
            "PFS": {"t_col": "canon__pfs_months",       "ev_col": "ENDO_PD_BIN"},
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

RISK_COHORTS = {
    "Breast":      {"case_sheet": "MT-T_N",  "control_sheet": "MN"},
    "Endometrial": {"case_sheet": "AT=AUs",  "control_sheet": "EN"},
}

# ── STYLE ─────────────────────────────────────────────────────────────────────
_COHORT_C = {"Breast": cohort_color("Breast"), "Endometrial": cohort_color("Endometrial")}

def _style_ax(ax, grid=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.yaxis.grid(True, linestyle=":", color="#cccccc", alpha=0.7)
        ax.set_axisbelow(True)

def _sig_label(p):
    if pd.isna(p):  return ""
    if p < 0.001:   return "***"
    if p < 0.01:    return "**"
    if p < 0.05:    return "*"
    return "ns"

# ── HELPERS ───────────────────────────────────────────────────────────────────

def _haldane_or(a, b, c, d) -> Tuple[float, str]:
    if 0 in (a, b, c, d):
        a, b, c, d = a+0.5, b+0.5, c+0.5, d+0.5
        note = "Corrected (+0.5)"
    else:
        note = "Standard"
    denom = b * c
    return ((a * d) / denom if denom else np.inf), note

def _apply_fdr(df: pd.DataFrame, p_col: str, out_col: str) -> pd.DataFrame:
    df = df.copy()
    valid = df[p_col].notna()
    if valid.sum() > 1:
        df.loc[valid, out_col] = false_discovery_control(
            df.loc[valid, p_col].values, method="bh"
        )
    else:
        df[out_col] = df[p_col]
    return df

def _run_logit(df_model, formula, label_prefix):
    """Fit logistic regression; return OR, CI, p keyed by label_prefix."""
    try:
        from statsmodels.formula.api import logit as sm_logit
        res = sm_logit(formula, data=df_model).fit(disp=0, method="bfgs", maxiter=200)
        coef = res.params.get("carrier", np.nan)
        ci   = res.conf_int()
        lo   = ci.loc["carrier", 0] if "carrier" in ci.index else np.nan
        hi   = ci.loc["carrier", 1] if "carrier" in ci.index else np.nan
        p    = res.pvalues.get("carrier", np.nan)
        return {
            f"OR_Adj_{label_prefix}":      round(np.exp(coef), 4),
            f"OR_Adj_{label_prefix}_CI95": f"[{round(np.exp(lo),3)}, {round(np.exp(hi),3)}]",
            f"P_Adj_{label_prefix}":       p,
            f"N_Adj_{label_prefix}":       len(df_model),
            f"Adj_{label_prefix}_Note":    "converged" if res.mle_retvals.get("converged") else "not converged",
        }
    except Exception as e:
        return {
            f"OR_Adj_{label_prefix}":      np.nan,
            f"OR_Adj_{label_prefix}_CI95": np.nan,
            f"P_Adj_{label_prefix}":       np.nan,
            f"N_Adj_{label_prefix}":       np.nan,
            f"Adj_{label_prefix}_Note":    str(e)[:80],
        }

# ── STATISTICAL TESTS ─────────────────────────────────────────────────────────

def _test_continuous(c_vals, nc_vals) -> Optional[Dict]:
    c  = pd.to_numeric(c_vals,  errors="coerce").dropna()
    nc = pd.to_numeric(nc_vals, errors="coerce").dropna()
    if len(c) < MIN_CARRIERS or len(nc) < MIN_CARRIERS:
        return None
    stat, p = mannwhitneyu(c, nc, alternative="two-sided")
    # Rank-biserial correlation: standard effect size for Mann-Whitney U
    # r = 1 - (2U) / (n1 * n2), ranges -1 to +1
    r_rb = round(1 - (2 * stat) / (len(c) * len(nc)), 4)
    return {
        "Test_Unadj":          "Mann-Whitney U",
        "N_Carriers":          len(c),
        "N_NonCarriers":       len(nc),
        "Stat_Unadj":          round(stat, 3),
        "P_Unadj":             p,
        "Effect_RankBiserial": r_rb,
        "Median_Carriers":     round(c.median(), 3),
        "Median_NonCarriers":  round(nc.median(), 3),
        "IQR_Carriers":        f"{round(c.quantile(0.25),3)}–{round(c.quantile(0.75),3)}",
        "IQR_NonCarriers":     f"{round(nc.quantile(0.25),3)}–{round(nc.quantile(0.75),3)}",
        "OR_Adj_Age": np.nan, "OR_Adj_Age_CI95": np.nan,
        "P_Adj_Age":  np.nan, "N_Adj_Age": np.nan,
        "OR_Adj_AgeBMI": np.nan, "OR_Adj_AgeBMI_CI95": np.nan,
        "P_Adj_AgeBMI":  np.nan, "N_Adj_AgeBMI": np.nan,
    }

def _test_binary(c_df, nc_df, col, age_col, bmi_col, include_bmi_model: bool = True) -> Optional[Dict]:
    c_vals  = pd.to_numeric(c_df[col],  errors="coerce")
    nc_vals = pd.to_numeric(nc_df[col], errors="coerce")
    a = int(c_vals.sum());  b = int(c_vals.notna().sum()) - a
    c = int(nc_vals.sum()); d = int(nc_vals.notna().sum()) - c
    if min(a+c, b+d) < MIN_CARRIERS:
        return None
    _, p_unadj = fisher_exact([[a, b], [c, d]])
    or_unadj, _ = _haldane_or(a, b, c, d)

    result = {
        "Test_Unadj":    "Fisher exact",
        "N_Carriers":    a + b,
        "N_NonCarriers": c + d,
        "Stat_Unadj":    round(or_unadj, 4),
        "P_Unadj":       p_unadj,
        "Median_Carriers":    np.nan,
        "Median_NonCarriers": np.nan,
    }

    # Build combined frame for adjusted models
    combined = pd.concat([
        c_df[[col, age_col, bmi_col]].assign(carrier=1),
        nc_df[[col, age_col, bmi_col]].assign(carrier=0),
    ]).rename(columns={col: "outcome", age_col: "age", bmi_col: "bmi"})
    combined["outcome"] = pd.to_numeric(combined["outcome"], errors="coerce")
    combined["age"]     = pd.to_numeric(combined["age"],     errors="coerce")
    combined["bmi"]     = pd.to_numeric(combined["bmi"],     errors="coerce")

    events_min = min(a + c, b + d)
    if events_min >= MIN_EVENTS_LOGISTIC:
        age_df = combined[["outcome", "age", "carrier"]].dropna().copy()
        if len(age_df) >= 6:
            age_df["age_z"] = (age_df["age"] - age_df["age"].mean()) / age_df["age"].std()
            result.update(_run_logit(age_df.rename(columns={"outcome": "outcome"}),
                                     "outcome ~ carrier + age_z", "Age"))
        else:
            result.update({k: np.nan for k in
                           ["OR_Adj_Age","OR_Adj_Age_CI95","P_Adj_Age","N_Adj_Age","Adj_Age_Note"]})

        if include_bmi_model:
            bmi_df = combined[["outcome", "age", "bmi", "carrier"]].dropna().copy()
            if len(bmi_df) >= 6:
                bmi_df["age_z"] = (bmi_df["age"] - bmi_df["age"].mean()) / bmi_df["age"].std()
                bmi_df["bmi_z"] = (bmi_df["bmi"] - bmi_df["bmi"].mean()) / bmi_df["bmi"].std()
                result.update(_run_logit(bmi_df.rename(columns={"outcome": "outcome"}),
                                         "outcome ~ carrier + age_z + bmi_z", "AgeBMI"))
            else:
                result.update({k: np.nan for k in
                               ["OR_Adj_AgeBMI","OR_Adj_AgeBMI_CI95","P_Adj_AgeBMI","N_Adj_AgeBMI","Adj_AgeBMI_Note"]})
        else:
            result.update({k: np.nan for k in
                           ["OR_Adj_AgeBMI","OR_Adj_AgeBMI_CI95","P_Adj_AgeBMI","N_Adj_AgeBMI","Adj_AgeBMI_Note"]})
    else:
        for sfx in ["Age", "AgeBMI"]:
            result.update({k: np.nan for k in
                           [f"OR_Adj_{sfx}", f"OR_Adj_{sfx}_CI95", f"P_Adj_{sfx}",
                            f"N_Adj_{sfx}", f"Adj_{sfx}_Note"]})
    return result

def _test_nominal(c_df, nc_df, col) -> Optional[Dict]:
    n_c  = c_df[col].notna().sum()
    n_nc = nc_df[col].notna().sum()
    if n_c < MIN_CARRIERS or n_nc < MIN_CARRIERS:
        return None
    ct = pd.crosstab(
        pd.concat([pd.Series(["C"] * len(c_df), index=c_df.index),
                   pd.Series(["NC"]* len(nc_df),index=nc_df.index)]),
        pd.concat([c_df[col], nc_df[col]])
    ).dropna(axis=1)
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return None
    # Drop columns where all groups have zero count (would give zero expected)
    ct = ct.loc[:, ct.sum(axis=0) > 0]
    if ct.shape[1] < 2:
        return None
    try:
        chi2, p, _, _ = stats.chi2_contingency(ct)
    except Exception:
        return None
    return {
        "Test_Unadj":    "Chi-square",
        "N_Carriers":    int(n_c),
        "N_NonCarriers": int(n_nc),
        "Stat_Unadj":    round(chi2, 3),
        "P_Unadj":       p,
        "Median_Carriers":    np.nan,
        "Median_NonCarriers": np.nan,
        "OR_Adj_Age":    np.nan, "OR_Adj_Age_CI95":    np.nan,
        "P_Adj_Age":     np.nan, "N_Adj_Age":          np.nan,
        "OR_Adj_AgeBMI": np.nan, "OR_Adj_AgeBMI_CI95": np.nan,
        "P_Adj_AgeBMI":  np.nan, "N_Adj_AgeBMI":       np.nan,
    }

# ── DATA LOADING ──────────────────────────────────────────────────────────────

def load_established_snp_positions(annot_path: Path) -> List[int]:
    """Read annotated report; return POS values for gnomAD NFE AF > 1% variants."""
    print("  Reading established SNP positions from annotated report …")
    df = pd.read_excel(annot_path, sheet_name="Biological_Annotations")

    # Case-insensitive column search
    def _find(df, pattern):
        for c in df.columns:
            if re.search(pattern, c, re.IGNORECASE):
                return c
        return None

    pos_col = _find(df, r"^pos$|^position$")
    nfe_col = _find(df, r"gnomad.*nfe.*af|nfe.*af")
    if not pos_col:
        raise ValueError("Cannot find POS column in annotated report")
    if not nfe_col:
        raise ValueError("Cannot find gnomAD NFE AF column in annotated report")

    df[nfe_col] = pd.to_numeric(df[nfe_col], errors="coerce")
    pos_list = (
        df[df[nfe_col] > 0.01][pos_col]
        .dropna()
        .astype(int)
        .unique()
        .tolist()
    )
    print(f"  Found {len(pos_list)} established SNP positions")
    return sorted(pos_list)


def load_phased_genotypes(phased_path: Path, keep_positions: List[int], ref_snp_labels: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Load BEAGLE phased genotype TSV; filter to established SNP positions.

    Returns a DataFrame indexed by sample name with columns:
      hap1  — phased allele string, chromosome 1 (e.g. "01011")
      hap2  — phased allele string, chromosome 2
    Each character corresponds to one SNP position (0=ref, 1=alt), in
    ascending positional order.
    """
    print("  Loading phased genotype TSV …")
    geno = pd.read_csv(phased_path, sep="\t", dtype=str)

    geno["POS"] = pd.to_numeric(geno["POS"], errors="coerce")
    geno["REF"] = geno["REF"].astype(str).str.strip().str.upper()
    geno["ALT"] = geno["ALT"].astype(str).str.strip().str.upper()
    geno["SNP_Label"] = (
        geno["POS"].astype("Int64").astype(str) + "_" +
        geno["REF"] + ">" + geno["ALT"]
    )

    # Prefer exact backbone matching by POS+REF+ALT when reference labels are available.
    if ref_snp_labels:
        ref_label_set = {str(lbl).strip() for lbl in ref_snp_labels}
        geno = geno[geno["SNP_Label"].isin(ref_label_set)].copy()
    else:
        geno = geno[geno["POS"].isin(keep_positions)].copy()

    geno = geno.sort_values("POS").drop_duplicates(subset=["SNP_Label"]).reset_index(drop=True)

    n_snps = len(geno)
    print(f"  {n_snps} established SNPs retained in phased matrix")
    if n_snps == 0:
        raise ValueError(
            "No SNPs remain after filtering phased genotypes to the requested haplotype backbone.\n"
            "Check that POS/REF/ALT labels in the phased TSV match the reference haplotype definition."
        )

    snp_labels = geno["SNP_Label"].tolist()

    sample_cols = [c for c in geno.columns if c not in ("CHROM", "POS", "ID", "REF", "ALT", "SNP_Label")]

    # Parse phased GTs: "0|1" -> allele1=0, allele2=1
    hap1_rows, hap2_rows = [], []
    for _, row in geno.iterrows():
        a1_list, a2_list = [], []
        for samp in sample_cols:
            gt = str(row[samp]).strip()
            if gt in (".", ".|.", "./.", ""):
                a1_list.append("N"); a2_list.append("N")
            elif "|" in gt:
                parts = gt.split("|")
                a1_list.append(parts[0]); a2_list.append(parts[1])
            elif "/" in gt:
                # Unphased fallback — treat as unordered
                parts = gt.split("/")
                a1_list.append(parts[0]); a2_list.append(parts[1])
            else:
                a1_list.append("N"); a2_list.append("N")
        hap1_rows.append(a1_list)
        hap2_rows.append(a2_list)

    # Build per-sample haplotype strings (length = n_snps)
    records = []
    for i, samp in enumerate(sample_cols):
        h1 = "".join(row[i] for row in hap1_rows)
        h2 = "".join(row[i] for row in hap2_rows)
        records.append({"Sample_phased": samp, "hap1": h1, "hap2": h2})

    hap_df = pd.DataFrame(records).set_index("Sample_phased")
    print(f"  Haplotype strings built for {len(hap_df)} samples "
          f"({n_snps} SNPs per string)")
    return hap_df, snp_labels


def enumerate_haplotypes(hap_df: pd.DataFrame,
                          min_freq: float = MIN_HAP_FREQ
                         ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Count all observed haplotypes; return:
      freq_df  — haplotype × frequency table (global and per haplotype string)
      carrier_df — sample × haplotype carrier matrix (1 = carrier on ≥1 chrom)
    """
    all_haps: Dict[str, int] = {}
    for h in ("hap1", "hap2"):
        for s in hap_df[h]:
            if "N" not in s:
                all_haps[s] = all_haps.get(s, 0) + 1

    total_chroms = sum(all_haps.values())
    freq_df = pd.DataFrame([
        {"Haplotype": hap, "Count": cnt,
         "Global_Freq": round(cnt / total_chroms, 4)}
        for hap, cnt in sorted(all_haps.items(), key=lambda x: -x[1])
    ])
    freq_df = freq_df[freq_df["Global_Freq"] >= min_freq].reset_index(drop=True)
    freq_df.index = [f"H{i+1}" for i in range(len(freq_df))]
    freq_df.index.name = "Haplotype_ID"
    freq_df = freq_df.reset_index()

    print(f"  {len(all_haps)} distinct haplotypes observed; "
          f"{len(freq_df)} with freq ≥ {min_freq}")

    # Build carrier matrix
    rows = []
    for _, freq_row in freq_df.iterrows():
        hap_str = freq_row["Haplotype"]
        hap_id  = freq_row["Haplotype_ID"]
        for samp in hap_df.index:
            h1 = hap_df.loc[samp, "hap1"]
            h2 = hap_df.loc[samp, "hap2"]
            callable_hap = int("N" not in h1 and "N" not in h2)
            if callable_hap:
                dosage = int(h1 == hap_str) + int(h2 == hap_str)
                carrier = int(dosage >= 1)
            else:
                dosage = np.nan
                carrier = np.nan
            rows.append({
                "Sample_phased": samp,
                "Haplotype_ID":  hap_id,
                "Haplotype":     hap_str,
                "Global_Freq":   freq_row["Global_Freq"],
                "Dosage":        dosage,
                "Carrier":       carrier,
                "Callable":      callable_hap,
            })
    carrier_df = pd.DataFrame(rows)
    return freq_df, carrier_df


def _extract_snp_code(sample_name: str) -> Optional[str]:
    s = str(sample_name)
    m = re.match(r"^(SNP_(?:AT|EN|MN|MT-T)_\d+)", s, re.IGNORECASE)
    if m: return m.group(1)
    m = re.match(r"^DNA_SNP_(EN|MN)_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_{m.group(1).upper()}_{m.group(2)}"
    m = re.match(r"^DNA_SNP_AT_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_AT_{m.group(1)}"
    m = re.match(r"^DNA_MT[-_]T_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_MT-T_{m.group(1)}"
    m = re.match(r"^SNP_DNA_(MN|EN)_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_{m.group(1).upper()}_{m.group(2)}"
    # New: DNA_SNP_MT-T_N_... with long library suffix
    m = re.match(r"^DNA_SNP_MT-T_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_MT-T_{m.group(1)}"
    return None


def _derive_breast_recurrence(val):
    if pd.isna(val): return None
    return 0 if str(val).strip().upper() == "NO" else 1

def _derive_breast_exitus(val):
    if pd.isna(val): return None
    v = str(val).strip().upper()
    if v.startswith("SI") or v.startswith("YES"): return 1
    if v.startswith("NO"): return 0
    return None

def _derive_breast_metastasis(val):
    if pd.isna(val): return None
    v = str(val).strip().upper()
    return 1 if v == "SI" else (0 if v in ("NO", "NO-LOCAL") else None)

def _derive_her2_subtype(val):
    if pd.isna(val): return None
    v = str(val).strip().upper()
    if "HER2+" in v or "HER2 +" in v: return "HER2+"
    if "TN" in v or "TRIPLE" in v: return "TN"
    return "Other"

def _derive_figo_numeric(val):
    if pd.isna(val): return None
    mapping = {"IA": 1.0, "IB": 1.5, "I": 1.5, "II": 2.0,
               "IIIA": 3.1, "IIIB": 3.2, "IIIC1": 3.3, "IIIC2": 3.4, "III": 3.0,
               "IVA": 4.1, "IVB": 4.2, "IV": 4.0}
    return mapping.get(str(val).strip().upper(), None)

def _derive_endo_grade(val):
    if pd.isna(val): return None
    m = re.search(r"(\d)", str(val))
    return float(m.group(1)) if m else None

def _derive_risk_ordinal(val):
    if pd.isna(val): return None
    return {"LOW": 1.0, "INTERMEDIATE": 2.0, "INTERMEDIATE-HIGH": 3.0, "HIGH": 4.0}.get(str(val).strip().upper(), None)

def _derive_ki67_numeric(val):
    if pd.isna(val): return None
    s = str(val).replace("%", "").strip()
    m = re.match(r"([\d.]+)\s*[-–]\s*([\d.]+)", s)
    if m:
        mid = (float(m.group(1)) + float(m.group(2))) / 2
        return mid / 100 if mid > 1 else mid
    try:
        v = float(s); return v / 100 if v > 1 else v
    except (ValueError, TypeError): return None

def _derive_breast_os(row):
    try:
        dx  = pd.to_datetime(row.get("clin_dcs__Fecha dx"),               errors="coerce")
        ult = pd.to_datetime(row.get("clin_dcs__Última fecha disponible"), errors="coerce")
        if pd.isna(dx) or pd.isna(ult): return None
        return max((ult - dx).days / 30.44, 0)
    except Exception: return None


def _ki67_fraction_from_pct(val):
    if pd.isna(val):
        return None
    try:
        v = float(val)
    except (TypeError, ValueError):
        return None
    return round(v / 100.0, 4) if 0.0 <= v <= 100.0 else None


def load_clinical_master(master_path: Path) -> pd.DataFrame:
    """
    Load 'harmonised_plus_canon' sheet, restrict to DNA rows, build Tissue,
    Cohort, and all derived clinical columns (mirrors script 17 logic).
    """
    print("  Loading harmonised clinical master …")
    try:
        master = pd.read_excel(master_path, sheet_name="harmonised_plus_canon")
    except Exception:
        print("  WARNING: 'harmonised_plus_canon' not found — trying sheet 0")
        master = pd.read_excel(master_path, sheet_name=0)

    master.columns = [str(c).strip() for c in master.columns]
    if "nucleic_acid" in master.columns:
        master = master[master["nucleic_acid"] == "DNA"].copy()

    if "snp_code" not in master.columns:
        for c in ["SNP_code", "Sample", "SAMPLE", "sample_id"]:
            if c in master.columns:
                master = master.rename(columns={c: "snp_code"}); break
    if "snp_code" not in master.columns:
        raise ValueError(f"Cannot find snp_code. Columns: {list(master.columns[:20])}")
    master["snp_code"] = master["snp_code"].astype(str).str.strip().str.replace(r"\s+", "", regex=True)

    if "extraction_flag" in master.columns and "sheet" in master.columns:
        excl = (master["sheet"] == "AT=AUs") & (
            master["extraction_flag"].isin(["NO HAY"]) |
            (master.get("pd_status", pd.Series("", index=master.index)) == "NO HACER")
        )
        master = master[~excl].copy()

    if "is_replicate" not in master.columns:
        master["is_replicate"] = master["snp_code"].str.contains("Repeticion", case=False, na=False)

    # Derive Tissue and Cohort from sheet column (most reliable)
    sheet_to_tissue = {"MT-T_N": "Tumour", "AT=AUs": "Tumour", "MN": "Healthy", "EN": "Healthy"}
    sheet_to_cohort = {"MT-T_N": "Breast_Tumour", "MN": "Breast_Healthy",
                       "AT=AUs": "Endometrial_Tumour", "EN": "Endometrial_Healthy"}
    if "sheet" in master.columns:
        master["Tissue"] = master["sheet"].map(sheet_to_tissue)
        master["Cohort"] = master["sheet"].map(sheet_to_cohort).fillna(master["sheet"])
    else:
        master["Tissue"] = None
        master["Cohort"] = "Unknown"

    # Breast derived columns
    for src, dst, fn in [
        ("clin_dcs__Recaida/Progresi\u00f3n", "BREAST_RECURRENCE_DERIVED", _derive_breast_recurrence),
        ("clin_dcs__Exitus",             "BREAST_EXITUS_DERIVED",     _derive_breast_exitus),
        ("clin_dcs__MTxDISTANCIA",       "BREAST_METASTASIS_DERIVED", _derive_breast_metastasis),
        ("clin_her2__DX",                "BREAST_HER2_SUBTYPE",       _derive_her2_subtype),
    ]:
        if src in master.columns:
            master[dst] = master[src].apply(fn)
    if "clin_dcs__MTxDISTANCIA" in master.columns:
        master["BREAST_ANY_METASTASIS_BIN"] = master["clin_dcs__MTxDISTANCIA"].map(
            {"NO": 0, "NO-LOCAL": 1, "SI": 1}
        )
    if "canon__er_status" in master.columns:
        master["BREAST_ER_BIN"] = master["canon__er_status"].map({"Positive": 1, "Negative": 0})
        master["ENDO_ER_BIN"]   = master["canon__er_status"].map({"Positive": 1, "Negative": 0})
    if "canon__pr_status" in master.columns:
        master["BREAST_PR_BIN"] = master["canon__pr_status"].map({"Positive": 1, "Negative": 0})
        master["ENDO_PR_BIN"]   = master["canon__pr_status"].map({"Positive": 1, "Negative": 0})

    ki67_pct = master["canon__ki67_pct"] if "canon__ki67_pct" in master.columns else pd.Series(np.nan, index=master.index)
    grade = master["canon__grade"] if "canon__grade" in master.columns else pd.Series(np.nan, index=master.index)
    bmi = master["canon__bmi"] if "canon__bmi" in master.columns else pd.Series(np.nan, index=master.index)
    breast_os_raw = master["canon__os_months"] if "canon__os_months" in master.columns else pd.Series(np.nan, index=master.index)
    figo = master["canon__figo_stage"] if "canon__figo_stage" in master.columns else pd.Series(np.nan, index=master.index)
    lvsi = master["canon__lvsi"] if "canon__lvsi" in master.columns else pd.Series(pd.NA, index=master.index)
    myoinv = master["canon__myometrial_invasion"] if "canon__myometrial_invasion" in master.columns else pd.Series(pd.NA, index=master.index)
    risk = master["canon__risk_group"] if "canon__risk_group" in master.columns else pd.Series(pd.NA, index=master.index)
    pd_flag = master["canon__pd_flag"] if "canon__pd_flag" in master.columns else pd.Series(pd.NA, index=master.index)
    exitus_flag = master["canon__exitus_flag"] if "canon__exitus_flag" in master.columns else pd.Series(pd.NA, index=master.index)

    master["BREAST_KI67_NUMERIC"]      = ki67_pct.apply(_ki67_fraction_from_pct)
    master["BREAST_GRADE_NUMERIC"]     = pd.to_numeric(grade, errors="coerce")
    master["BREAST_P53_NUMERIC"]       = pd.to_numeric(master.get("clin_dcs__p53", pd.Series(np.nan, index=master.index)), errors="coerce")
    master["BREAST_BMI_NUMERIC"]       = pd.to_numeric(bmi, errors="coerce")
    master["BREAST_MENARCHE_NUMERIC"]  = pd.to_numeric(master.get("canon__menarche_age", pd.Series(np.nan, index=master.index)), errors="coerce")
    master["BREAST_MENOPAUSE_NUMERIC"] = pd.to_numeric(master.get("canon__menopause_age", pd.Series(np.nan, index=master.index)), errors="coerce")
    master["BREAST_OS_MONTHS_DERIVED"] = pd.to_numeric(breast_os_raw, errors="coerce").combine_first(master.apply(_derive_breast_os, axis=1))
    if "clin_dcs__Dx" in master.columns:
        master["BREAST_DX_TYPE"] = master["clin_dcs__Dx"].apply(
            lambda x: "CDI" if pd.notna(x) and "CDI" in str(x).upper() and "CDIS" not in str(x).upper()
            else ("CDIS" if pd.notna(x) and "CDIS" in str(x).upper()
                  else (str(x).strip() if pd.notna(x) else None)))

    # Endometrial derived columns
    master["ENDO_FIGO_NUMERIC"]       = figo.apply(_derive_figo_numeric)
    master["ENDO_GRADE_NUMERIC"]      = pd.to_numeric(grade, errors="coerce")
    master["ENDO_RISK_ORDINAL"]       = risk.apply(_derive_risk_ordinal)
    master["ENDO_KI67_NUMERIC"]       = ki67_pct.apply(_ki67_fraction_from_pct)
    master["ENDO_LVSI_BIN"]           = lvsi.map({"Yes": 1, "No": 0})
    master["ENDO_MYOINV_BIN"]         = myoinv.map({"<50%": 0, ">50%": 1})
    master["ENDO_MSI_BIN"]            = master.get("canon__msi_status", pd.Series(pd.NA, index=master.index)).map({"Unstable": 1, "Stable": 0})
    master["ENDO_NEEC_BIN"]           = master.get("clin_au_endo__HISTOLOGY_GROUP", pd.Series(pd.NA, index=master.index)).map({"NEEC": 1, "EEC": 0})
    master["ENDO_PD_BIN"]             = pd.to_numeric(pd_flag, errors="coerce")
    master["ENDO_EXITUS_BIN"]         = pd.to_numeric(exitus_flag, errors="coerce")
    master["ENDO_EXITUS_DISEASE_BIN"] = master.get("clin_au_endo__EXITUS_DISEASE", pd.Series(pd.NA, index=master.index)).map({"YES": 1, "NO": 0})
    master["ENDO_TP53_ABN_BIN"]       = master.get("canon__p53_status", pd.Series(pd.NA, index=master.index)).map({"Aberrant": 1, "Normal": 0})
    master["ENDO_GENE_AMP_BIN"]       = master.get("clin_au_endo__Gene amplification", pd.Series(pd.NA, index=master.index)).map({"YES": 1, "NO": 0})
    master["ENDO_ITH_BIN"]            = master.get("clin_au_endo__ITH: intratumor heterogeneity", pd.Series(pd.NA, index=master.index)).map({"YES": 1, "NO": 0})
    master["ENDO_CTDNA_BIN"]          = master.get("clin_au_endo__BLOOD_BASAL_CTDNA", pd.Series(pd.NA, index=master.index)).map({"POSITIVE": 1, "NEGATIVE": 0})
    master["ENDO_CTDNA_MAF_NUMERIC"]  = pd.to_numeric(master.get("clin_au_endo__BLOOD_BASAL_CTDNA_MAF", pd.Series(np.nan, index=master.index)), errors="coerce")
    master["ENDO_CFDN_CONC_NUMERIC"]  = pd.to_numeric(master.get("clin_au_endo__BLOOD_BASAL_CFDNA_CONCENTRATION", pd.Series(np.nan, index=master.index)), errors="coerce")
    master["ENDO_PDL1_NUMERIC"]       = pd.to_numeric(master.get("clin_au_endo__FFPE_PDL1 POLAND RESULTS", pd.Series(np.nan, index=master.index)), errors="coerce")
    master["ENDO_CD8_NUMERIC"]        = pd.to_numeric(master.get("clin_au_endo__FFPE_CD8", pd.Series(np.nan, index=master.index)), errors="coerce")
    master["ENDO_PTEN_BIN"]           = master.get("clin_au_endo__FFPE_PTEN", pd.Series(pd.NA, index=master.index)).map({"CONSERVED": 0, "LOST/REDUCED": 1})
    master["ENDO_MLH1_BIN"]           = master.get("clin_au_endo__FFPE_MLH1", pd.Series(pd.NA, index=master.index)).map({"CONSERVED": 0, "LOST/REDUCED": 1})
    for src, dst, ref in [("clin_au_endo__N", "ENDO_N_STAGE_BIN", "N0"), ("clin_au_endo__M", "ENDO_M_STAGE_BIN", "M0")]:
        if src in master.columns:
            master[dst] = master[src].apply(
                lambda x, r=ref: 0 if pd.notna(x) and str(x).strip().upper() == r else (1 if pd.notna(x) else None))
    if "canon__bmi" not in master.columns:
        master["canon__bmi"] = np.nan

    master = master.drop_duplicates("snp_code").copy()
    print(f"  Master: {len(master)} unique samples")
    return master


def merge_haplotypes_with_clinical(carrier_df: pd.DataFrame,
                                   master: pd.DataFrame) -> pd.DataFrame:
    """
    Join carrier_df to clinical master. Converts BAM-style phased TSV column
    names (DNA_MT-T_37_IonCode_0121) to snp_codes (SNP_MT-T_37) before joining.
    """
    print("  Joining haplotype carrier data to clinical master …")
    carrier_df = carrier_df.copy()
    carrier_df["snp_code"] = carrier_df["Sample_phased"].apply(_extract_snp_code)

    n_failed = carrier_df["snp_code"].isna().sum()
    if n_failed:
        examples = carrier_df.loc[carrier_df["snp_code"].isna(), "Sample_phased"].unique()[:5].tolist()
        print(f"  WARNING: could not extract snp_code for {n_failed} rows (e.g. {examples})")
    carrier_df = carrier_df.dropna(subset=["snp_code"])

    merged = carrier_df.merge(master, on="snp_code", how="inner")
    n_phased  = carrier_df["snp_code"].nunique()
    n_matched = merged["snp_code"].nunique()
    print(f"  Phased samples: {n_phased} | matched to master: {n_matched} "
          f"| unmatched: {n_phased - n_matched}")

    if n_matched == 0:
        print(f"  Extracted codes (first 8): {sorted(carrier_df['snp_code'].dropna().unique())[:8]}")
        print(f"  Master codes    (first 8): {sorted(master['snp_code'].unique())[:8]}")
    else:
        unique_samp = merged.drop_duplicates("snp_code")
        print(f"  Tissue breakdown : {unique_samp['Tissue'].value_counts().to_dict()}")
        print(f"  Cohort breakdown : {unique_samp['Cohort'].value_counts().to_dict()}")
    return merged


# ── ANALYSIS 1: TUMOUR VS HEALTHY ─────────────────────────────────────────────

def tumour_vs_control(merged: pd.DataFrame) -> pd.DataFrame:
    print("=== Analysis 1: Tumour vs Control ===")
    df = merged[~merged["is_replicate"]].copy()

    rows = []
    for hap_id, hdf in df.groupby("Haplotype_ID"):
        hap_str  = hdf["Haplotype"].iloc[0]
        hap_freq = hdf["Global_Freq"].iloc[0]
        sample_hdf = hdf.drop_duplicates("snp_code").copy()
        if "Callable" in sample_hdf.columns:
            sample_hdf = sample_hdf[sample_hdf["Callable"] == 1]
        sample_hdf = sample_hdf.dropna(subset=["Carrier"])
        if sample_hdf.empty:
            continue
        for cohort_label in ["Breast", "Endometrial"]:
            coh_filter = "Breast" if cohort_label == "Breast" else "Endometri"
            tum_df = sample_hdf[
                sample_hdf["Cohort"].str.contains(coh_filter, case=False, na=False)
                & (sample_hdf["Tissue"] == "Tumour")
            ]
            hlt_df = sample_hdf[sample_hdf["Tissue"] == "Healthy"]
            if tum_df.empty or hlt_df.empty:
                continue

            tum_samp = tum_df.drop_duplicates("snp_code")
            hlt_samp = hlt_df.drop_duplicates("snp_code")
            if len(tum_samp) < MIN_CARRIERS or len(hlt_samp) < MIN_CARRIERS:
                continue
            a = int(tum_samp["Carrier"].sum())
            b = len(tum_samp) - a
            c = int(hlt_samp["Carrier"].sum())
            d = len(hlt_samp) - c

            if (a + c) < MIN_COMPARISON_CARRIERS or (b + d) < MIN_COMPARISON_CARRIERS:
                continue

            _, p = fisher_exact([[a, b], [c, d]])
            or_val, or_meth = _haldane_or(a, b, c, d)
            rows.append({
                "Cohort":          cohort_label,
                "Haplotype_ID":    hap_id,
                "Haplotype":       hap_str,
                "Global_Freq":     hap_freq,
                "N_Tumour":        len(tum_samp),
                "N_Control":       len(hlt_samp),
                "Carriers_Tumour": int(a),
                "Carriers_Control":int(c),
                "Freq_Tumour_%":   round(a / len(tum_samp) * 100, 2),
                "Freq_Control_%":  round(c / len(hlt_samp) * 100, 2),
                "OR":              round(or_val, 4),
                "OR_Method":       or_meth,
                "P_Value":         p,
                "Nominal_Sig":     p < 0.05,
            })

    res = pd.DataFrame(rows)
    if not res.empty:
        res = _apply_fdr(res, "P_Value", "FDR_P_Value")
        res["FDR_Sig"] = res["FDR_P_Value"] < FDR_THRESHOLD
        print(f"  {len(res)} tests | {res['Nominal_Sig'].sum()} nominal | "
              f"{res['FDR_Sig'].sum()} FDR<{FDR_THRESHOLD}\n")
    else:
        print("  No results.\n")
    return res


# ── ANALYSIS 2: HAPLOTYPE × CLINICAL VARIABLE ─────────────────────────────────

def _run_clin_for_cohort(df_t: pd.DataFrame,
                          var_dict: Dict,
                          cohort_label: str) -> pd.DataFrame:
    rows = []
    for hap_id, hdf in df_t.groupby("Haplotype_ID"):
        hap_str  = hdf["Haplotype"].iloc[0]
        hap_freq = hdf["Global_Freq"].iloc[0]
        sample_data = hdf.drop_duplicates("snp_code").set_index("snp_code")
        if "Callable" in sample_data.columns:
            sample_data = sample_data[sample_data["Callable"] == 1]
        sample_data = sample_data.dropna(subset=["Carrier"])
        c_df  = sample_data[sample_data["Carrier"] == 1]
        nc_df = sample_data[sample_data["Carrier"] == 0]
        if len(c_df) < MIN_CARRIERS or len(nc_df) < MIN_CARRIERS:
            continue
        bmi_col = next((c for c in ["canon__bmi", "BREAST_BMI_NUMERIC", "bmi", "BMI"]
                         if c in sample_data.columns), "canon__bmi")
        for col, meta in var_dict.items():
            if col not in sample_data.columns:
                continue
            if meta["type"] == "continuous":
                res = _test_continuous(c_df[col], nc_df[col])
            elif meta["type"] == "binary":
                # BMI adjustment only for cancer risk (Analysis 5); here we
                # adjust for age only to avoid over-adjustment within tumour cohorts
                res = _test_binary(c_df, nc_df, col,
                                   age_col="canon__age", bmi_col=bmi_col,
                                   include_bmi_model=meta.get("bmi_adjust", False))
            elif meta["type"] == "nominal":
                res = _test_nominal(c_df, nc_df, col)
            else:
                res = None
            if res is None:
                continue
            row = {
                "Cohort":        cohort_label,
                "Haplotype_ID":  hap_id,
                "Haplotype":     hap_str,
                "Global_Freq":   hap_freq,
                "Clinical_Var":  col,
                "Clin_Label":    meta["label"],
                "Clin_Type":     meta["type"],
            }
            row.update(res)
            rows.append(row)

    res_df = pd.DataFrame(rows)
    if res_df.empty:
        return res_df

    # FDR within each clinical variable across all haplotypes — correct unit.
    # Pooling all haplotype × variable combinations is too conservative.
    parts = []
    for clin_col in res_df["Clinical_Var"].unique():
        sub = res_df[res_df["Clinical_Var"] == clin_col].copy()
        sub = _apply_fdr(sub, "P_Unadj", "FDR_Unadj")
        sub["Nominal_Sig_Unadj"] = sub["P_Unadj"] < 0.05
        sub["FDR_Sig_Unadj"]     = sub["FDR_Unadj"] < FDR_THRESHOLD
        if "P_Adj_Age" in sub.columns and sub["P_Adj_Age"].notna().any():
            sub = _apply_fdr(sub, "P_Adj_Age", "FDR_Adj_Age")
            sub["Nominal_Sig_Adj_Age"] = sub["P_Adj_Age"] < 0.05
            sub["FDR_Sig_Adj_Age"]     = sub["FDR_Adj_Age"] < FDR_THRESHOLD
        else:
            sub["FDR_Adj_Age"] = np.nan
            sub["Nominal_Sig_Adj_Age"] = False
            sub["FDR_Sig_Adj_Age"]     = False
        parts.append(sub)

    res_df = pd.concat(parts, ignore_index=True).sort_values(["Clinical_Var", "P_Unadj"])
    return res_df


def clinical_associations(merged: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    print("=== Analysis 2: Haplotype × Clinical Variable ===")
    df = merged[~merged["is_replicate"] & (merged["Tissue"] == "Tumour")].copy()
    breast_df = df[df["Cohort"].str.contains("Breast",    case=False, na=False)]
    endo_df   = df[df["Cohort"].str.contains("Endometri", case=False, na=False)]
    print(f"  Breast tumour:      {breast_df['snp_code'].nunique()} samples")
    print(f"  Endometrial tumour: {endo_df['snp_code'].nunique()} samples")

    breast_res = _run_clin_for_cohort(breast_df, CLINICAL_VARS_BREAST, "Breast_Tumour")
    endo_res   = _run_clin_for_cohort(endo_df,   CLINICAL_VARS_ENDO,   "Endometrial_Tumour")

    for label, res in [("Breast", breast_res), ("Endometrial", endo_res)]:
        if res.empty:
            print(f"  {label}: no results")
        else:
            nom    = res["Nominal_Sig_Unadj"].sum()
            fdr    = res["FDR_Sig_Unadj"].sum()
            age_n  = res["Nominal_Sig_Adj_Age"].sum()    if "Nominal_Sig_Adj_Age"    in res.columns else 0
            bmi_n  = res["Nominal_Sig_Adj_AgeBMI"].sum() if "Nominal_Sig_Adj_AgeBMI" in res.columns else 0
            print(f"  {label}: {len(res)} tests | {nom} nominal | {fdr} FDR | "
                  f"{age_n} nominal (age-adj) | {bmi_n} nominal (age+BMI-adj)")
    print()
    return breast_res, endo_res


# ── ANALYSIS 3: HAPLOTYPE DOSE ────────────────────────────────────────────────

def haplotype_dose_analysis(merged: pd.DataFrame,
                             top_haplotypes: List[str],
                             var_dict: Dict,
                             cohort_label: str) -> pd.DataFrame:
    """
    For nominally significant haplotypes from analysis 2, test dosage
    (0 / 1 / 2 copies) against clinical variables using a trend test.
    """
    print(f"=== Analysis 3: Haplotype-dose — {cohort_label} ===")
    if not top_haplotypes:
        print("  No nominally significant haplotypes — skipping.\n")
        return pd.DataFrame()

    df = merged[~merged["is_replicate"] & (merged["Tissue"] == "Tumour")].copy()
    coh_filter = "Breast" if "Breast" in cohort_label else "Endometri"
    df = df[df["Cohort"].str.contains(coh_filter, case=False, na=False)]

    rows = []
    for hap_id in top_haplotypes:
        hdf = df[df["Haplotype_ID"] == hap_id].drop_duplicates("snp_code").set_index("snp_code")
        if "Callable" in hdf.columns:
            hdf = hdf[hdf["Callable"] == 1]
        hdf = hdf.dropna(subset=["Dosage"])
        if hdf.empty:
            continue
        hap_str = hdf["Haplotype"].iloc[0]

        for col, meta in var_dict.items():
            if col not in hdf.columns:
                continue
            g0 = hdf[hdf["Dosage"] == 0]
            g1 = hdf[hdf["Dosage"] == 1]
            g2 = hdf[hdf["Dosage"] == 2]

            if meta["type"] == "continuous":
                vals = {k: pd.to_numeric(grp[col], errors="coerce").dropna()
                        for k, grp in [("WT", g0), ("Het", g1), ("Hom", g2)]}
                non_empty = [v for v in vals.values() if len(v) >= MIN_CARRIERS]
                if len(non_empty) < 2:
                    continue
                # Jonckheere-Terpstra trend test via Spearman correlation
                # of dosage (0/1/2) vs. outcome — tests monotonic ordered
                # alternative, which is scientifically correct for dose-response
                dose_vals = pd.concat([
                    pd.Series([d]*len(v), name="dose")
                    for d, v in zip([0,1,2],[
                        pd.to_numeric(g0[col], errors="coerce").dropna(),
                        pd.to_numeric(g1[col], errors="coerce").dropna(),
                        pd.to_numeric(g2[col], errors="coerce").dropna(),
                    ])
                ])
                out_vals = pd.concat([
                    pd.to_numeric(g0[col], errors="coerce").dropna(),
                    pd.to_numeric(g1[col], errors="coerce").dropna(),
                    pd.to_numeric(g2[col], errors="coerce").dropna(),
                ])
                if len(dose_vals) < 4:
                    continue
                stat, p_trend = stats.spearmanr(dose_vals.values, out_vals.values)
            elif meta["type"] == "binary":
                ct_rows = []
                for grp in [g0, g1, g2]:
                    v = pd.to_numeric(grp[col], errors="coerce").dropna()
                    if len(v) < MIN_CARRIERS:
                        continue
                    ct_rows.append([int(v.sum()), len(v) - int(v.sum())])
                if len(ct_rows) < 2:
                    continue
                ct_arr = np.array(ct_rows)
                # Skip if any column sums to zero (chi2 would crash on zero expected)
                if (ct_arr.sum(axis=0) == 0).any():
                    continue
                try:
                    if ct_arr.shape == (2, 2):
                        _, p_trend = fisher_exact(ct_arr.tolist())
                        stat = p_trend  # no chi2 stat for Fisher
                    else:
                        chi2, p_trend, _, _ = stats.chi2_contingency(ct_arr)
                        stat = chi2
                except Exception:
                    continue
            else:
                continue

            rows.append({
                "Cohort":       cohort_label,
                "Haplotype_ID": hap_id,
                "Haplotype":    hap_str,
                "Clinical_Var": col,
                "Clin_Label":   meta["label"],
                "N_WT":   len(g0), "N_Het": len(g1), "N_Hom": len(g2),
                "Stat_Trend":   round(stat, 3),
                "P_Trend":      p_trend,
            })

    res = pd.DataFrame(rows)
    if not res.empty:
        res = _apply_fdr(res, "P_Trend", "FDR_Trend")
        res["Nominal_Sig"] = res["P_Trend"] < 0.05
        res["FDR_Sig"]     = res["FDR_Trend"] < FDR_THRESHOLD
        print(f"  {len(res)} tests | {res['Nominal_Sig'].sum()} nominal | "
              f"{res['FDR_Sig'].sum()} FDR\n")
    else:
        print("  No results.\n")
    return res


# ── ANALYSIS 4: SURVIVAL ──────────────────────────────────────────────────────

def survival_analysis(merged: pd.DataFrame) -> Tuple[pd.DataFrame, List]:
    if not _HAS_LIFELINES:
        print("=== Analysis 4: Survival — SKIPPED (lifelines not installed) ===\n")
        return pd.DataFrame(), []

    print("=== Analysis 4: Survival — all cohorts ===")
    df = merged[~merged["is_replicate"] & (merged["Tissue"] == "Tumour")].copy()

    rows, km_pages = [], []

    for cohort_label, cfg in SURVIVAL_COHORTS.items():
        cdf = df[df["Cohort"].str.contains(cfg["cohort_filter"], case=False, na=False)]
        if cdf.empty:
            continue
        age_col = cfg["age_col"]
        print(f"  {cohort_label}: {cdf['snp_code'].nunique()} tumour samples")

        for endpoint, ep_cfg in cfg["endpoints"].items():
            t_col  = ep_cfg["t_col"]
            ev_col = ep_cfg["ev_col"]
            if t_col not in cdf.columns or ev_col not in cdf.columns:
                print(f"    {endpoint}: columns missing — skipping.")
                continue
            print(f"    {endpoint} …")

            for hap_id, hdf in cdf.groupby("Haplotype_ID"):
                hap_str = hdf["Haplotype"].iloc[0]
                s = hdf.drop_duplicates("snp_code").set_index("snp_code").copy()
                needed = [t_col, ev_col]
                if age_col in s.columns: needed.append(age_col)
                if "canon__bmi" in s.columns: needed.append("canon__bmi")
                valid = s[needed].copy()
                valid[t_col]  = pd.to_numeric(valid[t_col],  errors="coerce")
                valid[ev_col] = pd.to_numeric(valid[ev_col], errors="coerce")
                valid = valid.dropna(subset=[t_col, ev_col])
                valid = valid[valid[t_col] > 0]
                if len(valid) < 6:
                    continue

                if "Callable" in s.columns:
                    s = s[s["Callable"] == 1].copy()
                carrier = pd.to_numeric(s.loc[valid.index, "Carrier"], errors="coerce")
                callable_idx = carrier.dropna().index
                if len(callable_idx) < 6:
                    continue
                valid = valid.loc[callable_idx].copy()
                carrier = carrier.loc[callable_idx].astype(int)
                n_c   = int(carrier.sum())
                n_nc  = int((carrier == 0).sum())
                if n_c < MIN_CARRIERS or n_nc < MIN_CARRIERS:
                    continue

                T = valid[t_col]
                E = valid[ev_col].astype(float)

                # Log-rank
                lr = logrank_test(T[carrier == 0], T[carrier == 1],
                                  E[carrier == 0], E[carrier == 1])
                lr_p = lr.p_value

                # Cox additive
                cox_df = pd.DataFrame({"T": T, "E": E, "carrier": carrier})
                if age_col in valid.columns:
                    cox_df["age"] = pd.to_numeric(valid[age_col], errors="coerce")
                    cox_df["age_z"] = (cox_df["age"] - cox_df["age"].mean()) / cox_df["age"].std()
                cox_df = cox_df.dropna()
                p_cox  = np.nan
                hr     = np.nan
                hr_lo  = np.nan
                hr_hi  = np.nan
                if len(cox_df) >= 6 and cox_df["carrier"].nunique() > 1:
                    try:
                        cph = CoxPHFitter()
                        cov = ["carrier"] + (["age_z"] if "age_z" in cox_df.columns else [])
                        cph.fit(cox_df[["T", "E"] + cov], duration_col="T", event_col="E")
                        p_cox  = cph.summary.loc["carrier", "p"]
                        hr     = cph.summary.loc["carrier", "exp(coef)"]
                        hr_lo  = cph.summary.loc["carrier", "exp(coef) lower 95%"]
                        hr_hi  = cph.summary.loc["carrier", "exp(coef) upper 95%"]
                    except Exception:
                        pass

                rows.append({
                    "Cohort":        cohort_label,
                    "Haplotype_ID":  hap_id,
                    "Haplotype":     hap_str,
                    "Endpoint":      endpoint,
                    "N_Total":       len(valid),
                    "N_Carriers":    int(n_c),
                    "N_NonCarriers": int(n_nc),
                    "P_LogRank":     lr_p,
                    "HR_Cox":        round(hr, 4) if not np.isnan(hr) else np.nan,
                    "HR_CI95_Lo":    round(hr_lo, 4) if not np.isnan(hr_lo) else np.nan,
                    "HR_CI95_Hi":    round(hr_hi, 4) if not np.isnan(hr_hi) else np.nan,
                    "P_Cox":         p_cox,
                })

                # KM plot
                fig, ax = plt.subplots(figsize=(6, 4))
                for grp_val, grp_label, col in [(0, "Non-carrier", HAPLOTYPE_STATUS_COLORS["Non-carrier"]),
                                                 (1, "Carrier",     HAPLOTYPE_STATUS_COLORS["Heterozygous carrier"])]:
                    mask = carrier == grp_val
                    if mask.sum() < 2:
                        continue
                    kmf = KaplanMeierFitter()
                    kmf.fit(T[mask], E[mask],
                            label=f"{grp_label} (n={mask.sum()})")
                    kmf.plot_survival_function(ax=ax, ci_show=True, color=col)
                ax.set_title(
                    f"Kaplan-Meier Analysis of {endpoint} According to {hap_id} Carrier Status\n"
                    f"{cohort_label} | Log-rank p={lr_p:.4f}",
                    fontsize=9,
                )
                ax.set_xlabel(f"{endpoint} (months)")
                ax.set_ylabel("Survival probability")
                ax.legend(fontsize=8)
                _style_ax(ax)
                plt.tight_layout()
                km_pages.append(fig)

    res = pd.DataFrame(rows)
    if not res.empty:
        fdr_parts = []
        for (cohort, ep), grp in res.groupby(["Cohort", "Endpoint"]):
            grp = _apply_fdr(grp, "P_Cox", "FDR_Cox") if grp["P_Cox"].notna().sum() > 1 else grp.assign(FDR_Cox=grp["P_Cox"])
            fdr_parts.append(grp)
        res = pd.concat(fdr_parts, ignore_index=True)
        res["Nominal_Sig"] = res["P_Cox"] < 0.05
        n_nom = res["Nominal_Sig"].sum()
        print(f"  {len(res)} tests | {n_nom} nominal (Cox)\n")
    else:
        print("  No survival results.\n")
    return res, km_pages


# ── ANALYSIS 5: CANCER RISK ───────────────────────────────────────────────────

def cancer_risk_analysis(merged: pd.DataFrame) -> pd.DataFrame:
    print("=== Analysis 5: Cancer Risk (Case-Control) ===")
    df = merged[~merged["is_replicate"]].copy()

    rows = []
    for cohort_label, cfg in RISK_COHORTS.items():
        case_df    = df[df["sheet"] == cfg["case_sheet"]].copy()
        control_df = df[df["sheet"] == cfg["control_sheet"]].copy()
        if "Callable" in case_df.columns:
            case_df = case_df[case_df["Callable"] == 1].copy()
        if "Callable" in control_df.columns:
            control_df = control_df[control_df["Callable"] == 1].copy()
        case_df = case_df.dropna(subset=["Carrier"])
        control_df = control_df.dropna(subset=["Carrier"])

        for hap_i, hap_id in enumerate(df["Haplotype_ID"].unique()):
            hap_str  = df[df["Haplotype_ID"] == hap_id]["Haplotype"].iloc[0]
            hap_freq = df[df["Haplotype_ID"] == hap_id]["Global_Freq"].iloc[0]

            hap_case = case_df[case_df["Haplotype_ID"] == hap_id].drop_duplicates("snp_code").set_index("snp_code")
            hap_control = control_df[control_df["Haplotype_ID"] == hap_id].drop_duplicates("snp_code").set_index("snp_code")
            n_cases = len(hap_case)
            n_controls = len(hap_control)
            if hap_i == 0:
                print(f"  {cohort_label}: {n_cases} callable cases, {n_controls} callable controls")
            if n_cases < MIN_CARRIERS or n_controls < MIN_CARRIERS:
                continue

            a = int(hap_case["Carrier"].sum())
            b = n_cases - a
            c = int(hap_control["Carrier"].sum())
            d = n_controls - c

            if (a + c) < MIN_COMPARISON_CARRIERS or (b + d) < MIN_COMPARISON_CARRIERS:
                continue

            _, p_unadj = fisher_exact([[a, b], [c, d]])
            or_unadj, or_meth = _haldane_or(a, b, c, d)

            # 95 % CI via Woolf / log method (with Haldane correction)
            a_c, b_c, c_c, d_c = (a+0.5, b+0.5, c+0.5, d+0.5) if 0 in (a,b,c,d) else (a,b,c,d)
            try:
                log_or   = np.log(a_c * d_c / (b_c * c_c))
                se_log   = np.sqrt(1/a_c + 1/b_c + 1/c_c + 1/d_c)
                or_ci_lo = round(np.exp(log_or - 1.96 * se_log), 4)
                or_ci_hi = round(np.exp(log_or + 1.96 * se_log), 4)
            except Exception:
                or_ci_lo, or_ci_hi = np.nan, np.nan

            row = {
                "Cohort":            cohort_label,
                "Haplotype_ID":      hap_id,
                "Haplotype":         hap_str,
                "Global_Freq":       hap_freq,
                "N_Cases":           n_cases,
                "N_Controls":        n_controls,
                "Carriers_Cases":    a,
                "Carriers_Controls": c,
                "Freq_Cases_%":      round(a / n_cases    * 100, 2),
                "Freq_Controls_%":   round(c / n_controls * 100, 2),
                "OR_Unadj":          round(or_unadj, 4),
                "OR_CI95_Lo":        or_ci_lo,
                "OR_CI95_Hi":        or_ci_hi,
                "OR_Method":         or_meth,
                "P_Unadj":           p_unadj,
            }

            # Age-adjusted logistic
            all_samp = pd.concat([
                hap_case[["canon__age", "canon__bmi", "Carrier"]].rename(columns={"Carrier": "carrier"}).assign(cancer=1),
                hap_control[["canon__age", "canon__bmi", "Carrier"]].rename(columns={"Carrier": "carrier"}).assign(cancer=0),
            ])
            all_samp["canon__age"] = pd.to_numeric(all_samp["canon__age"], errors="coerce")
            all_samp["canon__bmi"] = pd.to_numeric(all_samp["canon__bmi"], errors="coerce")

            if min(a + c, b + d) >= MIN_EVENTS_LOGISTIC:
                age_df = all_samp[["cancer", "canon__age", "carrier"]].dropna().copy()
                if len(age_df) >= 6:
                    age_df["age_z"] = (age_df["canon__age"] - age_df["canon__age"].mean()) / \
                                       age_df["canon__age"].std()
                    adj = _run_logit(age_df.rename(columns={"cancer": "outcome"}),
                                     "outcome ~ carrier + age_z", "Age")
                    row["OR_Adj_Age"]      = adj.get("OR_Adj_Age", np.nan)
                    row["OR_Adj_Age_CI95"] = adj.get("OR_Adj_Age_CI95", np.nan)
                    row["P_Adj_Age"]       = adj.get("P_Adj_Age", np.nan)
                    row["N_Adj_Age"]       = adj.get("N_Adj_Age", np.nan)
                else:
                    row.update({k: np.nan for k in
                                ["OR_Adj_Age","OR_Adj_Age_CI95","P_Adj_Age","N_Adj_Age"]})
            else:
                row.update({k: np.nan for k in
                            ["OR_Adj_Age","OR_Adj_Age_CI95","P_Adj_Age","N_Adj_Age"]})
            rows.append(row)

    res = pd.DataFrame(rows)
    if not res.empty:
        for cohort in res["Cohort"].unique():
            mask = res["Cohort"] == cohort
            res = _apply_fdr(res[mask], "P_Unadj", "FDR_Unadj").combine_first(res)
        res["Nominal_Sig_Unadj"] = res["P_Unadj"] < 0.05
        res["FDR_Sig_Unadj"]     = res["FDR_Unadj"] < FDR_THRESHOLD
        if "P_Adj_Age" in res.columns:
            res["Nominal_Sig_Adj_Age"] = res["P_Adj_Age"] < 0.05
        for cohort in res["Cohort"].unique():
            sub = res[res["Cohort"] == cohort]
            print(f"  {cohort}: {sub['Nominal_Sig_Unadj'].sum()} nominal (unadj) | "
                  f"{sub['FDR_Sig_Unadj'].sum()} FDR (unadj) | "
                  f"{sub.get('Nominal_Sig_Adj_Age', pd.Series([False]*len(sub))).sum()} nominal (age-adj)")
    print()
    return res


# ── VISUALISATION ─────────────────────────────────────────────────────────────

def make_volcano(tvh: pd.DataFrame, out_dir: Path):
    if tvh.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    for cohort, col in _COHORT_C.items():
        sub = tvh[tvh["Cohort"] == cohort]
        if sub.empty:
            continue
        ax.scatter(sub["OR"].apply(lambda x: np.log2(x)),
                   -np.log10(sub["P_Value"]),
                   c=col, label=f"{cohort} (tumour n={int(sub['N_Tumour'].iloc[0])}, control n={int(sub['N_Control'].iloc[0])})", alpha=0.7, s=60, zorder=3)
        for _, row in sub[sub["Nominal_Sig"]].iterrows():
            ax.annotate(row["Haplotype_ID"],
                        (np.log2(row["OR"]), -np.log10(row["P_Value"])),
                        fontsize=7, ha="left", va="bottom")
    ax.axhline(-np.log10(0.05), ls="--", c="#aaaaaa", lw=1, label="p=0.05")
    ax.axvline(0, ls=":", c="#cccccc", lw=0.8)
    ax.set_xlabel("log₂(OR)  —  tumour vs control")
    ax.set_ylabel("-log₁₀(p)")
    ax.set_title(f"Association Between Haplotype Carrier Status and Tumour-Control Status [{COMPARATIVE_TAG}]", fontweight="bold")
    ax.legend(fontsize=9)
    _style_ax(ax)
    plt.tight_layout()
    out = out_dir / "18_Haplo_Volcano_raw_p.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def make_heatmap(clin_res: pd.DataFrame, cohort_label: str,
                 p_col: str, title_suffix: str, out_path: Path):
    """
    Dual-panel heatmap showing ALL haplotypes with global frequency >= MIN_HAP_FREQ.

    Panel A  -- -log10(p) for each haplotype x clinical variable combination.
                Cells with p < 0.05 annotated with "*"; p < 0.01 with "**".
    Panel B  -- Effect size: log2(OR) for binary/nominal variables, or
                median(carriers) - median(non-carriers) for continuous variables.

    A frequency bar to the right contextualises haplotype prevalence.
    X-axis labels are rotated and sized to avoid overlap.
    """
    if clin_res.empty:
        print(f"  Heatmap ({cohort_label} | {title_suffix}): no data -- skipped")
        return

    sub = clin_res.copy()
    if sub.empty:
        print(f"  Heatmap ({cohort_label} | {title_suffix}): empty -- skipped")
        return

    # Build pivot tables
    pivot_p = sub.pivot_table(index="Haplotype_ID", columns="Clin_Label",
                               values=p_col, aggfunc="min")

    # Effect size: median diff for continuous, OR for binary/nominal
    if "Clin_Type" in sub.columns and "Median_Carriers" in sub.columns:
        sub = sub.copy()
        sub["_ES"] = np.where(
            sub["Clin_Type"] == "continuous",
            pd.to_numeric(sub["Median_Carriers"],    errors="coerce") -
            pd.to_numeric(sub["Median_NonCarriers"], errors="coerce"),
            pd.to_numeric(sub["Stat_Unadj"],         errors="coerce"),
        )
    else:
        sub["_ES"] = pd.to_numeric(sub.get("Stat_Unadj", np.nan), errors="coerce")

    pivot_es = sub.pivot_table(index="Haplotype_ID", columns="Clin_Label",
                                values="_ES", aggfunc="mean")

    all_haps = pivot_p.index.union(pivot_es.index)
    all_vars = pivot_p.columns.union(pivot_es.columns)
    pivot_p  = pivot_p.reindex(index=all_haps, columns=all_vars)
    pivot_es = pivot_es.reindex(index=all_haps, columns=all_vars)

    # Sort haplotypes H1, H2, ...
    import re as _re
    def _hk(h):
        m = _re.search(r"(\d+)", str(h))
        return int(m.group(1)) if m else 999
    sorted_haps = sorted(all_haps, key=_hk)
    pivot_p  = pivot_p.loc[sorted_haps]
    pivot_es = pivot_es.loc[sorted_haps]

    # Per-haplotype global frequency
    freq_map = {}
    if "Haplotype_ID" in sub.columns and "Global_Freq" in sub.columns:
        freq_map = (sub.drop_duplicates("Haplotype_ID")
                       .set_index("Haplotype_ID")["Global_Freq"].to_dict())
    freqs = np.array([freq_map.get(h, np.nan) for h in sorted_haps])

    n_haps = len(sorted_haps)
    n_vars  = len(all_vars)

    # Significance annotation matrix
    def _sa(p):
        if pd.isna(p): return ""
        if p < 0.01:   return "**"
        if p < 0.05:   return "*"
        return ""
    annot_mat = pivot_p.applymap(_sa)

    # Figure sizing: extra width to prevent x-label overlap
    # Allow at least 1.1 inches per variable; row height 0.65 per haplotype
    col_w  = max(1.0, 12 / max(n_vars, 1))
    row_h  = max(0.65, 6  / max(n_haps, 1))
    # Bottom margin for rotated labels -- scale with longest label
    max_lbl = max((len(str(c)) for c in all_vars), default=10)
    bottom_margin = min(0.40, max(0.22, max_lbl * 0.012))

    fig_w = max(12, n_vars * col_w * 2 + 4)   # x2 for two panels + space
    fig_h = max(5,  n_haps * row_h + 2.5)

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = fig.add_gridspec(1, 3,
                           width_ratios=[n_vars, n_vars, max(1, n_vars // 8)],
                           wspace=0.30,
                           left=0.10, right=0.95,
                           top=0.88, bottom=bottom_margin)
    ax_p   = fig.add_subplot(gs[0])
    ax_es  = fig.add_subplot(gs[1])
    ax_bar = fig.add_subplot(gs[2])

    # --- Panel A: -log10(p) ---
    log_p_mat = -np.log10(pivot_p.fillna(1).clip(lower=1e-300))
    vmax_p = max(3.0, float(log_p_mat.max().max()))

    # Font size for annotations scales with grid size
    ann_fs = max(5, min(8, int(60 / max(n_vars * n_haps, 1) ** 0.5)))

    sns.heatmap(log_p_mat, ax=ax_p, cmap="YlOrRd",
                annot=annot_mat, fmt="s",
                annot_kws={"size": ann_fs + 1, "weight": "bold", "color": "#333333"},
                linewidths=0.25, linecolor="#e0e0e0",
                vmin=0, vmax=vmax_p,
                cbar_kws={"label": "-log10(p)", "shrink": 0.55},
                xticklabels=True, yticklabels=True)
    ax_p.set_title(f"Statistical significance [{title_suffix}]", fontsize=9, fontweight="bold", pad=6)
    ax_p.set_xlabel("")
    ax_p.set_ylabel("Haplotype", fontsize=9)
    # Rotate x labels; auto-size font to avoid overlap
    x_fs = max(6, min(9, int(120 / max(n_vars, 1))))
    ax_p.tick_params(axis="x", rotation=45, labelsize=x_fs)
    ax_p.tick_params(axis="y", labelsize=8.5)
    ax_p.set_xticklabels(ax_p.get_xticklabels(), ha="right", rotation_mode="anchor")

    # --- Panel B: effect size ---
    es_vals = pivot_es.values.astype(float)
    has_non_cont = ("Clin_Type" in sub.columns and
                    (sub["Clin_Type"] != "continuous").any())
    if has_non_cont:
        # log2(OR) for OR-type columns; leave continuous median diffs as-is
        # Identify which columns are continuous
        cont_cols = set(sub[sub["Clin_Type"] == "continuous"]["Clin_Label"].unique())
        es_display = np.copy(es_vals)
        for ci, col in enumerate(pivot_es.columns):
            if col not in cont_cols:
                col_vals = es_vals[:, ci]
                mask = np.isfinite(col_vals) & (col_vals > 0)
                es_display[mask, ci] = np.log2(col_vals[mask])
                es_display[~mask & np.isfinite(col_vals), ci] = np.nan
        cbar_label = "log2(OR)  /  Δmedian"
    else:
        es_display = es_vals
        cbar_label = "Δmedian (carriers − non-carriers)"

    es_df = pd.DataFrame(es_display, index=pivot_es.index, columns=pivot_es.columns)
    finite = es_display[np.isfinite(es_display)]
    abs_max = max(float(np.nanmax(np.abs(finite))) if len(finite) > 0 else 1.0, 0.5)

    sns.heatmap(es_df, ax=ax_es, cmap="coolwarm",
                annot=True, fmt=".2f",
                annot_kws={"size": ann_fs},
                linewidths=0.25, linecolor="#e0e0e0",
                center=0, vmin=-abs_max, vmax=abs_max,
                cbar_kws={"label": cbar_label, "shrink": 0.55},
                xticklabels=True, yticklabels=False)
    ax_es.set_title("Effect size estimate", fontsize=9, fontweight="bold", pad=6)
    ax_es.set_xlabel("")
    ax_es.set_ylabel("")
    ax_es.tick_params(axis="x", rotation=45, labelsize=x_fs)
    ax_es.set_xticklabels(ax_es.get_xticklabels(), ha="right", rotation_mode="anchor")

    # --- Frequency bar ---
    y_pos   = np.arange(n_haps) + 0.5
    blues   = plt.cm.Blues(np.linspace(0.35, 0.85, n_haps))  # type: ignore
    ax_bar.barh(y_pos, np.where(np.isnan(freqs), 0, freqs) * 100,
                height=0.7, color=blues, edgecolor="#555555", linewidth=0.4)
    finite_f = freqs[np.isfinite(freqs)]
    xmax_bar = float(np.nanmax(finite_f)) * 120 if len(finite_f) > 0 else 10
    ax_bar.set_xlim(0, max(xmax_bar, 5))
    ax_bar.set_ylim(0, n_haps)
    ax_bar.set_yticks([])
    ax_bar.set_xlabel("Global\nfreq (%)", fontsize=7.5)
    ax_bar.set_title("Global frequency", fontsize=9, fontweight="bold")
    ax_bar.tick_params(axis="x", labelsize=7)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    for i, f in enumerate(freqs):
        if np.isfinite(f):
            ax_bar.text(f * 100 + 0.3, y_pos[i], f"{f*100:.1f}%",
                        va="center", fontsize=6.5, color="#333333")

    fig.suptitle(
        f"Associations Between GSDMB Haplotypes and Clinical Variables in {cohort_label} ({title_suffix})\n"
        f"All haplotypes >= {MIN_HAP_FREQ*100:.0f}% global frequency  |  "
        f"* p<0.05   ** p<0.01",
        fontsize=10,
        fontweight="bold",
        y=0.97,
    )

    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def make_forest(clin_res: pd.DataFrame, cohort_label: str, out_path: Path):
    if clin_res.empty or "Stat_Unadj" not in clin_res.columns:
        return
    sub = clin_res[clin_res["Nominal_Sig_Unadj"] &
                   (clin_res["Clin_Type"] == "binary")].copy()
    if sub.empty:
        print(f"  Forest ({cohort_label}): no nominal binary results — skipped")
        return
    sub = sub.sort_values("P_Unadj")
    fig, ax = plt.subplots(figsize=(7, max(4, len(sub) * 0.45)))
    ys = range(len(sub))
    col = _COHORT_C.get(cohort_label.split("_")[0], "#333333")
    ax.scatter(sub["Stat_Unadj"], ys, color=col, zorder=3, s=50)
    ax.axvline(1, ls="--", c="#aaaaaa", lw=1)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([f"{r['Haplotype_ID']} | {r['Clin_Label']}"
                        for _, r in sub.iterrows()], fontsize=8)
    ax.set_xlabel("Odds Ratio (unadjusted)")
    ax.set_title(
        f"Associations Between GSDMB Haplotypes and Binary Clinical Outcomes in {cohort_label}",
        fontweight="bold",
    )
    _style_ax(ax, grid=False)
    ax.xaxis.grid(True, linestyle=":", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def make_risk_forest(risk_res: pd.DataFrame, out_path: Path):
    """
    Forest plot for cancer risk analysis.
    Shows OR point estimate with 95% CI whiskers, carrier counts, and p-value
    annotations. Haplotypes are sorted by cohort then by OR.
    """
    if risk_res.empty or "OR_Unadj" not in risk_res.columns:
        return

    # Sort: Breast first, then Endometrial; within each cohort sort by OR
    df = risk_res.copy().reset_index(drop=True)
    cohort_order = {"Breast": 0, "Endometrial": 1}
    df["_cohort_rank"] = df["Cohort"].map(cohort_order).fillna(2)
    df = df.sort_values(["_cohort_rank", "OR_Unadj"], ascending=[True, True]).reset_index(drop=True)

    n   = len(df)
    ys  = np.arange(n)
    colours = [_COHORT_C.get(r["Cohort"], "#555555") for _, r in df.iterrows()]

    fig, ax = plt.subplots(figsize=(9, max(5, n * 0.5 + 1.5)))

    for i, (_, row) in enumerate(df.iterrows()):
        or_v  = row["OR_Unadj"]
        lo    = row.get("OR_CI95_Lo", np.nan)
        hi    = row.get("OR_CI95_Hi", np.nan)
        col   = colours[i]
        p_val = row.get("P_Unadj", np.nan)

        # CI whisker
        if pd.notna(lo) and pd.notna(hi):
            ax.plot([lo, hi], [i, i], color=col, lw=1.5, alpha=0.7, zorder=2)
            ax.plot([lo, lo], [i - 0.12, i + 0.12], color=col, lw=1.5, alpha=0.7)
            ax.plot([hi, hi], [i - 0.12, i + 0.12], color=col, lw=1.5, alpha=0.7)

        # Point estimate
        ax.scatter([or_v], [i], color=col, s=70, zorder=4, edgecolors="white", linewidth=0.6)

        # p-value and carrier count annotation on the right
        n_c    = int(row.get("Carriers_Cases",    0))
        n_ctrl = int(row.get("Carriers_Controls", 0))
        sig    = _sig_label(p_val)
        p_txt  = f"p={p_val:.3f}" if pd.notna(p_val) else ""
        ann    = f"{p_txt} {sig}  (cases: {n_c}, ctrl: {n_ctrl})".strip()

        # Position annotation just beyond the right edge of the CI (or OR if no CI)
        x_ann = hi if pd.notna(hi) else or_v
        ax.text(x_ann, i, f"  {ann}", va="center", fontsize=7, color="#333333")

    # Reference line at OR = 1
    ax.axvline(1, ls="--", c="#aaaaaa", lw=1, zorder=1)

    # Y-axis labels: haplotype ID + cohort
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{r['Haplotype_ID']}  ({r['Cohort']})"
                        for _, r in df.iterrows()], fontsize=8.5)
    ax.set_xlabel("Odds Ratio — unadjusted (95% CI)", fontsize=10)
    ax.set_title(
        "Associations Between GSDMB Haplotypes and Cancer Risk",
        fontweight="bold",
        fontsize=11,
    )

    # Log scale x-axis helps when CIs are very wide
    all_vals = pd.concat([df["OR_Unadj"],
                          df.get("OR_CI95_Lo", pd.Series(dtype=float)),
                          df.get("OR_CI95_Hi", pd.Series(dtype=float))]).dropna()
    if all_vals.min() > 0:
        ax.set_xscale("log")
        ax.set_xlabel("Odds Ratio — unadjusted (95% CI, log scale)", fontsize=10)
        from matplotlib.ticker import ScalarFormatter
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.get_major_formatter().set_scientific(False)

    _style_ax(ax, grid=False)
    ax.xaxis.grid(True, linestyle=":", alpha=0.4)
    handles = [mpatches.Patch(color=c, label=l) for l, c in _COHORT_C.items()]
    ax.legend(handles=handles, fontsize=9, loc="lower right")

    # Add a note about small-n haplotypes
    small_n = df[(df.get("Carriers_Cases", pd.Series(0, index=df.index)) +
                  df.get("Carriers_Controls", pd.Series(0, index=df.index))) <= 5]
    if not small_n.empty:
        note = "Note: wide CIs reflect small carrier numbers (n carriers + controls ≤ 5)"
        fig.text(0.5, 0.01, note, ha="center", fontsize=7.5,
                 color="#888888", style="italic")

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def _select_haplotype_heatmap_focus_rows(clin_res: pd.DataFrame, p_col: str, max_rows: int = 10) -> pd.DataFrame:
    """Return a cleaner subset of haplotypes for main-text heatmap presentation."""
    if clin_res.empty or p_col not in clin_res.columns:
        return clin_res
    focus = clin_res[clin_res[p_col].notna()].copy()
    if focus.empty:
        return clin_res
    rank_df = (focus.groupby('Haplotype_ID', as_index=False)[p_col]
               .min()
               .sort_values(p_col))
    if 'FDR_Sig_Unadj' in focus.columns and focus['FDR_Sig_Unadj'].any():
        keep = rank_df[rank_df['Haplotype_ID'].isin(focus.loc[focus['FDR_Sig_Unadj'], 'Haplotype_ID'])].head(max_rows)
    elif 'Nominal_Sig_Unadj' in focus.columns and focus['Nominal_Sig_Unadj'].any():
        keep = rank_df[rank_df['Haplotype_ID'].isin(focus.loc[focus['Nominal_Sig_Unadj'], 'Haplotype_ID'])].head(max_rows)
    else:
        keep = rank_df.head(max_rows)
    return focus[focus['Haplotype_ID'].isin(set(keep['Haplotype_ID']))].copy()


def make_main_text_heatmaps(b_clin: pd.DataFrame, e_clin: pd.DataFrame, out_dir: Path):
    """Generate focused haplotype heatmaps suitable for main-text presentation."""
    for clin_res, cohort_label in [(b_clin, 'Breast'), (e_clin, 'Endometrial')]:
        if clin_res.empty:
            continue
        focused = _select_haplotype_heatmap_focus_rows(clin_res, 'P_Unadj', max_rows=10)
        make_heatmap(
            focused,
            cohort_label,
            'P_Unadj',
            'Focused main-text view',
            out_dir / f'18_Haplo_Heatmap_{cohort_label}_MainText.png',
        )


def make_haplotype_composition_plots(tvh: pd.DataFrame, merged: pd.DataFrame, out_dir: Path):
    """Plot non-carrier / heterozygous / homozygous haplotype composition by comparison arm."""
    if tvh.empty:
        return

    top_n = 8
    status_order = ['Non-carrier', 'Heterozygous carrier', 'Homozygous carrier']
    dosage_labels = {0: 'Non-carrier', 1: 'Heterozygous carrier', 2: 'Homozygous carrier'}

    sample_manifest = merged[~merged['is_replicate']].drop_duplicates('snp_code').copy()
    pooled_controls = sample_manifest[sample_manifest['Tissue'] == 'Healthy'].copy()

    comparison_defs = {
        'Breast': {
            'tumour_mask': (sample_manifest['Tissue'] == 'Tumour') & sample_manifest['Cohort'].astype(str).str.contains('Breast', case=False, na=False),
            'tumour_label': 'Breast tumour',
            'control_label': 'Pooled control',
        },
        'Endometrial': {
            'tumour_mask': (sample_manifest['Tissue'] == 'Tumour') & sample_manifest['Cohort'].astype(str).str.contains('Endometri', case=False, na=False),
            'tumour_label': 'Endometrium tumour',
            'control_label': 'Pooled control',
        },
    }

    for cohort, cfg in comparison_defs.items():
        sub = (tvh[tvh['Cohort'] == cohort]
               .sort_values(['FDR_Sig', 'Nominal_Sig', 'P_Value'], ascending=[False, False, True])
               .drop_duplicates('Haplotype_ID'))
        if sub.empty:
            continue
        chosen = sub[sub['Nominal_Sig'] | sub['FDR_Sig']].head(top_n).copy() if (sub['Nominal_Sig'].any() or sub['FDR_Sig'].any()) else sub.head(top_n).copy()
        if chosen.empty:
            continue

        tumour_samples = sample_manifest.loc[cfg['tumour_mask'], ['snp_code']].copy()
        control_samples = pooled_controls[['snp_code']].copy()
        if tumour_samples.empty or control_samples.empty:
            continue

        arm_df = pd.concat([
            control_samples.assign(_arm=cfg['control_label']),
            tumour_samples.assign(_arm=cfg['tumour_label']),
        ], ignore_index=True)

        n_plots = len(chosen)
        ncols = min(3, n_plots)
        nrows = int(np.ceil(n_plots / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5.1 * ncols, 4.3 * nrows), squeeze=False, facecolor='white')
        axes_flat = axes.flatten()

        for ax in axes_flat:
            ax.set_facecolor('white')

        for idx, (_, row) in enumerate(chosen.iterrows()):
            ax = axes_flat[idx]
            hap_id = row['Haplotype_ID']
            hap_lookup = (merged[(~merged['is_replicate']) & (merged['Haplotype_ID'] == hap_id)]
                          .drop_duplicates('snp_code')
                          .set_index('snp_code')['Dosage'])
            plot_df = arm_df.copy()
            plot_df['_dosage'] = plot_df['snp_code'].map(hap_lookup)
            plot_df = plot_df[plot_df['_dosage'].notna()].copy()
            plot_df['Status'] = plot_df['_dosage'].astype(int).map(dosage_labels)

            counts = (plot_df.groupby(['_arm', 'Status']).size().unstack(fill_value=0)
                      .reindex(index=[cfg['control_label'], cfg['tumour_label']], columns=status_order, fill_value=0))
            totals = counts.sum(axis=1)
            perc = counts.div(totals, axis=0) * 100

            xpos = np.arange(len(counts.index))
            bottoms = np.zeros(len(counts.index))
            for status in status_order:
                vals = perc[status].to_numpy(dtype=float)
                ax.bar(xpos, vals, bottom=bottoms, color=HAPLOTYPE_STATUS_COLORS[status], edgecolor='white', linewidth=1.0, width=0.58)
                for xi, val, bottom in zip(xpos, vals, bottoms):
                    if val >= 9:
                        ax.text(xi, bottom + val / 2, f'{val:.0f}%', ha='center', va='center', fontsize=8.5, color='white', fontweight='bold')
                bottoms += vals

            ax.set_ylim(0, 100)
            ax.set_xticks(xpos)
            ax.set_xticklabels([f'{label}\n(n={int(totals.loc[label])})' for label in counts.index], fontsize=8.5)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v:.0f}%'))
            ax.set_ylabel('Haplotype composition (%)', fontsize=9)
            p_val = row['FDR_P_Value'] if bool(row.get('FDR_Sig', False)) else row['P_Value']
            sig_label = 'FDR' if bool(row.get('FDR_Sig', False)) else ('p<0.05' if bool(row.get('Nominal_Sig', False)) else 'top hit')
            ax.set_title(f'Haplotype {hap_id}\nOR={row["OR"]:.2f} | p={p_val:.3g} ({sig_label})', fontsize=8.5, fontweight='bold')
            _style_ax(ax)

        for ax in axes_flat[n_plots:]:
            ax.set_visible(False)

        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=HAPLOTYPE_STATUS_COLORS[s], edgecolor='white') for s in status_order]
        fig.legend(handles, status_order, title='Haplotype status', loc='upper center', ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
        fig.suptitle(f'Haplotype Composition in Tumour and Control Samples: {cohort} Comparison [{COMPARATIVE_TAG}]', fontsize=12, fontweight='bold', y=1.04)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        out_path = out_dir / f'18_HaplotypeComposition_{cohort}.png'
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f'  Saved: {out_path}')


def make_km_pdf(km_pages: List, out_path: Path):
    if not km_pages:
        return
    with pdf_backend.PdfPages(out_path) as pdf:
        for fig in km_pages:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"  Saved: {out_path}  ({len(km_pages)} pages)")


# ── MAIN ─────────────────────────────────────────────────────────────────────




def _order_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Move key identifier columns to the front; leave the rest in original order."""
    if df is None or df.empty:
        return df
    priority = [
        "Region_Name", "Region_Type", "Region_Sheet",
        "Haplotype_ID", "Haplotype", "Haplotype_Source",
        "Global_Freq", "Count", "Num_SNPs", "SNP_Positions", "SNP_Labels",
        "Cohort", "Tissue", "snp_code",
        "Clin_Var", "Clin_Label", "Variable", "Analysis_Type",
    ]
    front = [c for c in priority if c in df.columns]
    rest  = [c for c in df.columns if c not in front]
    return df[front + rest]


def _normalise_haplotype_meta(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure ref_freq_df has the required columns that downstream functions expect:
      Haplotype_ID, Haplotype, Global_Freq, Count, Num_SNPs, SNP_Positions, SNP_Labels.
    The Haplotype string (0/1 per SNP) and Haplotype_ID are already built by
    discover_reference_regions; this function just standardises types and fills gaps.
    """
    df = df.copy()

    # Haplotype_ID: prefer existing column, else use index
    if "Haplotype_ID" not in df.columns:
        if df.index.name == "Haplotype_ID":
            df = df.reset_index()
        else:
            df["Haplotype_ID"] = [f"H{i+1}" for i in range(len(df))]
    df["Haplotype_ID"] = df["Haplotype_ID"].astype(str).str.strip()

    # Haplotype string (0/1 per allele, length = Num_SNPs)
    if "Haplotype" not in df.columns:
        df["Haplotype"] = ""
    df["Haplotype"] = df["Haplotype"].astype(str)

    # Frequency / Count
    for col, default in [("Global_Freq", np.nan), ("Count", np.nan),
                         ("Num_SNPs", 0), ("SNP_Positions", ""), ("SNP_Labels", "")]:
        if col not in df.columns:
            df[col] = default
    df["Global_Freq"] = pd.to_numeric(df["Global_Freq"], errors="coerce")
    df["Count"]       = pd.to_numeric(df["Count"],       errors="coerce")
    df["Num_SNPs"]    = pd.to_numeric(df["Num_SNPs"],    errors="coerce").fillna(0).astype(int)

    return df


def reconcile_reference_to_phased_backbone(
    ref_positions: List[int],
    ref_labels:    List[str],
    ref_freq_df:   pd.DataFrame,
    phased_labels: List[str],
) -> Tuple[List[int], List[str], pd.DataFrame]:
    """
    When the SNP labels in the phased TSV don't exactly match those in the
    reference haplotype table (e.g. a SNP present in the reference is missing
    from the phased data, or vice-versa), restrict both to the common set and
    rebuild the Haplotype strings in ref_freq_df accordingly.

    Returns updated (positions, labels, ref_freq_df).
    """
    ref_set    = set(ref_labels)
    phased_set = set(phased_labels)
    common     = [lbl for lbl in ref_labels if lbl in phased_set]  # preserve ref order

    if not common:
        print("  WARNING: No overlapping SNP labels between reference and phased data. "
              "Returning reference as-is; carrier matching may fail.")
        return ref_positions, ref_labels, ref_freq_df

    dropped_ref    = ref_set    - phased_set
    dropped_phased = phased_set - ref_set
    if dropped_ref:
        print(f"  Reconcile: {len(dropped_ref)} reference SNPs absent from phased data "
              f"(dropped from haplotype strings): {sorted(dropped_ref)[:5]}{'…' if len(dropped_ref)>5 else ''}")
    if dropped_phased:
        print(f"  Reconcile: {len(dropped_phased)} phased SNPs absent from reference "
              f"(ignored): {sorted(dropped_phased)[:5]}{'…' if len(dropped_phased)>5 else ''}")

    # Indices of kept positions within the original reference label list
    keep_idx   = [i for i, lbl in enumerate(ref_labels) if lbl in phased_set]
    new_pos    = [ref_positions[i] for i in keep_idx]
    new_labels = [ref_labels[i]    for i in keep_idx]

    # Rebuild Haplotype strings by extracting the kept character positions
    df = ref_freq_df.copy()
    if "Haplotype" in df.columns and df["Haplotype"].str.len().gt(0).any():
        def _trim_hap(hap_str):
            if len(hap_str) == len(ref_labels):
                return "".join(hap_str[i] for i in keep_idx)
            return hap_str  # already trimmed or unexpected length — leave alone
        df["Haplotype"] = df["Haplotype"].apply(_trim_hap)

    df["Num_SNPs"]     = len(new_pos)
    df["SNP_Positions"] = ";".join(map(str, new_pos))
    df["SNP_Labels"]    = ";".join(new_labels)

    return new_pos, new_labels, df


def build_carrier_matrix_from_reference(
    hap_df:         pd.DataFrame,
    freq_df:        pd.DataFrame,
    ref_snp_labels: Optional[List[str]] = None,
    phased_snp_labels: Optional[List[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Given phased haplotype strings (hap_df: index=Sample_phased, cols hap1/hap2)
    and a reference frequency table (freq_df with Haplotype and Haplotype_ID cols),
    compute per-sample carrier status for each reference haplotype.

    ref_snp_labels / phased_snp_labels: when provided (and lengths differ),
    the phased haplotype strings are subsetted to only the positions that appear
    in the reference, so that exact string matching works correctly after
    reconciliation.

    Returns (freq_df_unchanged, carrier_df).
    """
    # Build index mapping: which characters in phased hap strings correspond to
    # each reference SNP label (handles the case where phased has extra SNPs).
    keep_idx: Optional[List[int]] = None
    if ref_snp_labels and phased_snp_labels and (ref_snp_labels != phased_snp_labels):
        phased_pos = {lbl: i for i, lbl in enumerate(phased_snp_labels)}
        keep_idx = [phased_pos[lbl] for lbl in ref_snp_labels if lbl in phased_pos]
        if len(keep_idx) != len(ref_snp_labels):
            # Some reference labels not found in phased — warn but continue
            missing = [l for l in ref_snp_labels if l not in phased_pos]
            print(f"  WARNING: {len(missing)} ref SNP labels not found in phased data: {missing[:3]}")
            keep_idx = None  # fall back to full-string comparison

    def _subset_hap(h: str) -> str:
        """Extract only the reference-position characters from a phased hap string."""
        if keep_idx is None or len(h) != len(phased_snp_labels or []):
            return h
        return "".join(h[i] for i in keep_idx)

    meta_cols = [c for c in freq_df.columns
                 if c not in ("Haplotype", "Haplotype_ID", "Global_Freq", "Count")]

    # Pre-compute subsetted haplotype strings for all samples (faster than re-slicing per ref hap)
    hap1_sub = {samp: _subset_hap(str(hap_df.loc[samp, "hap1"])) for samp in hap_df.index}
    hap2_sub = {samp: _subset_hap(str(hap_df.loc[samp, "hap2"])) for samp in hap_df.index}

    rows = []
    for _, freq_row in freq_df.iterrows():
        hap_str = str(freq_row["Haplotype"])
        hap_id  = str(freq_row["Haplotype_ID"])
        gfreq   = freq_row.get("Global_Freq", np.nan)

        for samp in hap_df.index:
            h1 = hap1_sub[samp]
            h2 = hap2_sub[samp]
            callable_hap = int("N" not in h1 and "N" not in h2)
            if callable_hap:
                dosage = int(h1 == hap_str) + int(h2 == hap_str)
                carrier = int(dosage >= 1)
            else:
                dosage = np.nan
                carrier = np.nan
            rec = {
                "Sample_phased": samp,
                "Haplotype_ID":  hap_id,
                "Haplotype":     hap_str,
                "Global_Freq":   gfreq,
                "Dosage":        dosage,
                "Carrier":       carrier,
                "Callable":      callable_hap,
            }
            for mc in meta_cols:
                rec[mc] = freq_row.get(mc, np.nan)
            rows.append(rec)

    carrier_df = pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["Sample_phased","Haplotype_ID","Haplotype","Global_Freq","Dosage","Carrier","Callable"]
    )

    n_carriers = carrier_df[carrier_df["Carrier"] == 1]["Sample_phased"].nunique() if not carrier_df.empty else 0
    print(f"  Carrier matrix: {len(hap_df)} samples × {len(freq_df)} haplotypes "
          f"({n_carriers} samples carry ≥1 reference haplotype)")
    return freq_df, carrier_df


META_HAP_COLS = {"Frequency", "Count", "Haplotype_ID", "Allele_String", "Alt_Alleles", "Variant_Content"}

def discover_reference_regions(haplo_results_path: Path) -> List[Dict]:
    print("  Scanning script 15 workbook for region-level haplotype sheets …")
    xl = pd.ExcelFile(haplo_results_path)
    sheets = xl.sheet_names
    if "SNPs_Used" not in sheets:
        raise ValueError("Script 15 workbook is missing SNPs_Used sheet.")
    snps_used = pd.read_excel(haplo_results_path, sheet_name="SNPs_Used")
    required = {"POS_int", "rsID_clean", "REF", "ALT"}
    missing = required - set(snps_used.columns)
    if missing:
        raise ValueError(f"SNPs_Used is missing required columns: {sorted(missing)}")

    snps_used["POS_int"] = pd.to_numeric(snps_used["POS_int"], errors="coerce")
    snps_used = snps_used.dropna(subset=["POS_int"]).copy()
    snps_used["POS_int"] = snps_used["POS_int"].astype(int)
    snps_used["rsID_clean"] = snps_used["rsID_clean"].astype(str).str.strip()
    snps_used["snp_label"] = (
        snps_used["POS_int"].astype(str) + "_" +
        snps_used["REF"].astype(str).str.strip() + ">" +
        snps_used["ALT"].astype(str).str.strip()
    )
    rsid_map = snps_used.set_index("rsID_clean")[["POS_int", "snp_label"]].to_dict("index")

    regions = []
    for sheet in sheets:
        if not sheet.endswith("_Freqs"):
            continue
        df = pd.read_excel(haplo_results_path, sheet_name=sheet)
        rsid_cols = [c for c in df.columns if c in rsid_map]
        if not rsid_cols:
            continue

        if sheet == "Full_Global_Freqs":
            region_name = "Full_Region"
            region_type = "Full"
        elif sheet.startswith("LD_Block"):
            region_name = sheet.replace("_Freqs", "")
            region_type = "LD_Block"
        elif sheet.startswith("Gene_"):
            region_name = sheet.replace("_Freqs", "")
            region_type = "Gene"
        else:
            region_name = sheet.replace("_Freqs", "")
            region_type = "Other"

        ref_freq_df = df.copy()
        hap_matrix = ref_freq_df[rsid_cols].apply(pd.to_numeric, errors="coerce").astype("Int64")
        if hap_matrix.isna().any().any():
            raise ValueError(f"Non-numeric haplotype allele codes found in {sheet}")
        ref_freq_df["Haplotype"] = (hap_matrix.astype(int) - 1).astype(str).agg("".join, axis=1)
        ref_freq_df["Global_Freq"] = pd.to_numeric(ref_freq_df.get("Frequency", np.nan), errors="coerce")
        ref_freq_df["Count"] = pd.to_numeric(ref_freq_df.get("Count", np.nan), errors="coerce")
        positions = [int(rsid_map[r]["POS_int"]) for r in rsid_cols]
        labels = [str(rsid_map[r]["snp_label"]) for r in rsid_cols]
        ref_freq_df["Num_SNPs"] = len(positions)
        ref_freq_df["SNP_Positions"] = ";".join(map(str, positions))
        ref_freq_df["SNP_Labels"] = ";".join(labels)
        ref_freq_df["Haplotype_Source"] = region_name
        ref_freq_df["Region_Name"] = region_name
        ref_freq_df["Region_Type"] = region_type
        ref_freq_df["Region_Sheet"] = sheet
        ref_freq_df = _normalise_haplotype_meta(ref_freq_df)

        regions.append({
            "region_name": region_name,
            "region_type": region_type,
            "sheet": sheet,
            "snp_positions": positions,
            "ref_snp_labels": labels,
            "ref_freq_df": ref_freq_df,
        })

    regions = sorted(regions, key=lambda r: (0 if r["region_type"]=="Full" else 1 if r["region_type"]=="LD_Block" else 2, r["region_name"]))
    print(f"  Detected {len(regions)} analysable regions:")
    for r in regions:
        print(f'    {r["region_name"]}  [{r["region_type"]}]  - {len(r["snp_positions"])} SNPs from {r["sheet"]}')
    return regions


def _prefix_region_columns(df: pd.DataFrame, region: Dict) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    out = df.copy()
    out["Region_Name"] = region["region_name"]
    out["Region_Type"] = region["region_type"]
    out["Region_Sheet"] = region["sheet"]
    return out


def run_region_analysis(region: Dict, args, master: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    print("\n" + "-" * 80)
    print(f'REGION: {region["region_name"]}  [{region["region_type"]}]')
    print("-" * 80 + "\n")

    print("[2/5] Loading and processing phased genotypes …")
    hap_df, snp_labels = load_phased_genotypes(Path(args.phased), region["snp_positions"], region["ref_snp_labels"])

    snp_positions = list(region["snp_positions"])
    ref_snp_labels = list(region["ref_snp_labels"])
    ref_freq_df = region["ref_freq_df"].copy()
    if snp_labels != ref_snp_labels:
        snp_positions, ref_snp_labels, ref_freq_df = reconcile_reference_to_phased_backbone(
            snp_positions, ref_snp_labels, ref_freq_df, snp_labels
        )

    print("\n[3/5] Defining haplotypes to test …")
    freq_df = ref_freq_df.copy()
    if args.min_freq is not None:
        before = len(freq_df)
        freq_df = freq_df[freq_df["Global_Freq"] >= args.min_freq].reset_index(drop=True)
        after = len(freq_df)
        if after != before:
            print(f"  Applied min_freq filter to reference haplotypes: {before} -> {after}")

    # prefix IDs to avoid collisions across regions
    freq_df["Haplotype_ID"] = freq_df["Haplotype_ID"].astype(str).apply(lambda x: f'{region["region_name"]}__{x}')
    freq_df["Haplotype_Source"] = region["region_name"]
    freq_df["Region_Name"] = region["region_name"]
    freq_df["Region_Type"] = region["region_type"]
    freq_df["Region_Sheet"] = region["sheet"]

    freq_df, carrier_df = build_carrier_matrix_from_reference(hap_df, freq_df, ref_snp_labels, snp_labels)

    print("\n[5/5] Merging haplotype data with clinical master …")
    merged = merge_haplotypes_with_clinical(carrier_df, master)
    if merged.empty:
        print("  No matched samples for this region, skipping analyses.")
        return {k: pd.DataFrame() for k in ["freq_df","merged","tvh","b_clin","e_clin","b_dose","e_dose","surv_res","risk_res","manifest"]}

    n_samples = merged["snp_code"].nunique()
    n_haplotypes = merged["Haplotype_ID"].nunique()
    print(f"\n  Merged: {n_samples} samples × {n_haplotypes} haplotypes\n")

    tvh = tumour_vs_control(merged)
    if not tvh.empty:
        validate_percentage_columns(tvh, ["Freq_Tumour_%", "Freq_Control_%"], "Script 18 tumour vs control")
    b_clin, e_clin = clinical_associations(merged)
    # Dose analysis on ALL haplotypes (not pre-filtered by Analysis 2 p-values,
    # which would introduce circular selection bias)
    all_hap_ids = merged["Haplotype_ID"].unique().tolist()
    b_dose = haplotype_dose_analysis(merged, all_hap_ids, CLINICAL_VARS_BREAST, "Breast_Tumour")
    e_dose = haplotype_dose_analysis(merged, all_hap_ids, CLINICAL_VARS_ENDO, "Endometrial_Tumour")
    surv_res, km_pages = survival_analysis(merged)
    risk_res = cancer_risk_analysis(merged)

    manifest_cols = [c for c in [
        "snp_code", "Cohort", "Tissue", "sheet", "is_replicate",
        "Haplotype_ID", "Haplotype", "Global_Freq", "Count", "Num_SNPs",
        "SNP_Positions", "SNP_Labels", "Haplotype_Source", "Dosage", "Carrier"
    ] if c in merged.columns]
    manifest = (
        merged[manifest_cols]
        .drop_duplicates(subset=["snp_code", "Haplotype_ID"])
        .sort_values(["Cohort", "Tissue", "snp_code"])
        .reset_index(drop=True)
    )

    return {
        "freq_df": _prefix_region_columns(freq_df, region),
        "merged": _prefix_region_columns(merged, region),
        "tvh": _prefix_region_columns(tvh, region),
        "b_clin": _prefix_region_columns(b_clin, region),
        "e_clin": _prefix_region_columns(e_clin, region),
        "b_dose": _prefix_region_columns(b_dose, region),
        "e_dose": _prefix_region_columns(e_dose, region),
        "surv_res": _prefix_region_columns(surv_res, region),
        "risk_res": _prefix_region_columns(risk_res, region),
        "manifest": _prefix_region_columns(manifest, region),
        "km_pages": km_pages,
    }


def parse_args():
    p = argparse.ArgumentParser(
        description="Haplotype–clinical association analysis (script 18)"
    )
    p.add_argument("--phased",   default=str(DEFAULT_PHASED),
                   help="BEAGLE phased genotype TSV")
    p.add_argument("--master",   default=str(DEFAULT_MASTER),
                   help="Harmonised clinical master Excel")
    p.add_argument("--annot",    default=str(DEFAULT_ANNOT),
                   help="GSDMB annotated report (for SNP position filter)")
    p.add_argument("--out_dir",  default=str(DEFAULT_OUT))
    p.add_argument("--min_freq", default=MIN_HAP_FREQ, type=float,
                   help="Minimum global haplotype frequency to test")
    p.add_argument("--haplo_results", default=str(DEFAULT_HAPLO_RESULTS),
                   help="Path to script 15 haplotype results Excel (optional; "
                        "enables reference-region-aware analysis)")
    return p.parse_args()



def main():
    args = parse_args()
    validate_file_exists(args.phased, "Script 18 phased haplotype input")
    validate_file_exists(args.master, "Script 18 master input")
    validate_file_exists(args.annot, "Script 18 annotation input")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_xlsx = out_dir / "18_Haplo_Clinical_Association_Results.xlsx"

    print("=" * 60)
    print("SCRIPT 18: HAPLOTYPE-CLINICAL ASSOCIATION ANALYSIS")
    print("=" * 60 + "\n")


    haplo_results_path = Path(args.haplo_results) if args.haplo_results else None

    # Hard-fail early with a clear message if the file doesn't exist
    if haplo_results_path is None or not haplo_results_path.exists():
        print(f"\nERROR: --haplo_results file not found: {haplo_results_path}")
        print("  This file is required to run the haplotype association analysis.")
        print("  It is the Excel output produced by script 15 (15_haplo_stats.r).")
        print(f"  Expected location: {DEFAULT_HAPLO_RESULTS}")
        print("\n  If the file is elsewhere, re-run with:")
        print("    python3 18_haplotype_association.py \\")
        print("      --haplo_results /actual/path/to/19_Haplotype_Results_v7_blocks_and_genes.xlsx")
        raise SystemExit(1)

    print("[1/5] Loading reference haplotypes from script 15 …")
    print(f"  Using: {haplo_results_path}")
    regions = discover_reference_regions(haplo_results_path)

    print("\n[4/5] Loading clinical master …")

    master = load_clinical_master(Path(args.master))
    print_validation_summary(master, "snp_code" if "snp_code" in master.columns else None, "Script 18 clinical master")

    all_res = []
    for region in regions:
        all_res.append(run_region_analysis(region, args, master))

    def _cat(key):
        dfs = [r[key] for r in all_res if key in r and r[key] is not None and not r[key].empty]
        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    freq_df  = _cat("freq_df")
    merged   = _cat("merged")
    tvh      = _cat("tvh")
    b_clin   = _cat("b_clin")
    e_clin   = _cat("e_clin")
    b_dose   = _cat("b_dose")
    e_dose   = _cat("e_dose")
    dose_all = pd.concat([b_dose, e_dose], ignore_index=True) if (not b_dose.empty or not e_dose.empty) else pd.DataFrame()
    surv_res = _cat("surv_res")
    risk_res = _cat("risk_res")
    manifest = _cat("manifest")
    all_km_pages = [pg for r in all_res for pg in r.get("km_pages", [])]

    print("\n=== Writing output Excel ===")

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        _any_sheet = False

        def _write_sheet(df, name):
            nonlocal _any_sheet
            if not df.empty:
                _order_columns(df).to_excel(xw, sheet_name=name, index=False)
                _any_sheet = True

        _write_sheet(freq_df,  "haplotype_frequencies")
        _write_sheet(tvh,      "tumour_vs_control")
        _write_sheet(b_clin,   "breast_clinical_assoc")
        _write_sheet(e_clin,   "endo_clinical_assoc")
        _write_sheet(dose_all, "haplotype_dose")
        _write_sheet(surv_res, "survival_cox")
        _write_sheet(risk_res, "cancer_risk")
        _write_sheet(manifest, "sample_manifest")

        sig_rows = []
        for res_df, analysis in [
            (tvh,      "Tumour_vs_Control"),
            (b_clin,   "Breast_Clinical"),
            (e_clin,   "Endo_Clinical"),
            (surv_res, "Survival"),
            (risk_res, "Cancer_Risk"),
        ]:
            if res_df.empty:
                continue
            p_col = (
                "P_Value" if "P_Value" in res_df.columns else
                "P_Unadj" if "P_Unadj" in res_df.columns else
                "P_Cox"   if "P_Cox"   in res_df.columns else None
            )
            if p_col:
                n_tests = res_df[p_col].notna().sum()
                bonf = 0.05 / n_tests if n_tests > 0 else 0.05
                sig = res_df[res_df[p_col] < 0.05].copy()
                sig["Analysis_Type"] = analysis
                sig["N_Tests_In_Set"] = n_tests
                sig["Bonferroni_Threshold"] = round(bonf, 6)
                sig["Passes_Bonferroni"] = sig[p_col] < bonf
                sig_rows.append(sig)
        summary = pd.concat(sig_rows, ignore_index=True) if sig_rows else pd.DataFrame()
        _write_sheet(summary, "summary_significant")

        # openpyxl requires ≥1 visible sheet — write a diagnostic sheet when
        # no analysis produced results (e.g. --haplo_results not supplied).
        if not _any_sheet:
            hint = (
                "No analysis results were produced.\n\n"
                "Most likely cause: --haplo_results was not provided, so no\n"
                "reference haplotype regions were available to analyse.\n\n"
                "Re-run with:\n"
                "  python3 18_haplotype_association.py \\\n"
                "    --haplo_results /path/to/19_Haplotype_Results.xlsx\n\n"
                f"Regions found    : {len(regions)}\n"
                f"use_reference    : {use_reference}\n"
                f"Regions analysed : {len(all_res)}\n"
            )
            print(f"\nWARNING: No result sheets written.\n{hint}")
            pd.DataFrame({"Diagnostic_Message": [hint]}).to_excel(
                xw, sheet_name="run_diagnostics", index=False
            )

    print(f"  Saved: {out_xlsx}\n")

    # ── FIGURE GENERATION ────────────────────────────────────────────────────
    print("=== Generating figures ===")

    # Volcano: tumour vs control
    make_volcano(tvh, out_dir)

    # Heatmaps: unadjusted p-values
    make_heatmap(
        b_clin, "Breast", "P_Unadj", "Unadjusted (full view)",
        out_dir / "18_Haplo_Heatmap_Breast.png",
    )
    make_heatmap(
        e_clin, "Endometrial", "P_Unadj", "Unadjusted (full view)",
        out_dir / "18_Haplo_Heatmap_Endometrial.png",
    )
    make_main_text_heatmaps(b_clin, e_clin, out_dir)

    # Heatmaps: age-adjusted p-values (if available)
    if not b_clin.empty and "P_Adj_Age" in b_clin.columns and b_clin["P_Adj_Age"].notna().any():
        make_heatmap(
            b_clin, "Breast", "P_Adj_Age", "Age-adjusted (full view)",
            out_dir / "18_Haplo_Heatmap_Breast_AgeAdj.png",
        )
    if not e_clin.empty and "P_Adj_Age" in e_clin.columns and e_clin["P_Adj_Age"].notna().any():
        make_heatmap(
            e_clin, "Endometrial", "P_Adj_Age", "Age-adjusted (full view)",
            out_dir / "18_Haplo_Heatmap_Endometrial_AgeAdj.png",
        )

    make_haplotype_composition_plots(tvh, merged, out_dir)

    # Forest plots: binary clinical associations
    make_forest(b_clin, "Breast",      out_dir / "18_Haplo_Forest_Breast.png")
    make_forest(e_clin, "Endometrial", out_dir / "18_Haplo_Forest_Endometrial.png")

    # Forest plot: cancer risk (case-control)
    make_risk_forest(risk_res, out_dir / "18_Haplo_Forest_CancerRisk.png")

    # Kaplan-Meier PDF
    make_km_pdf(all_km_pages, out_dir / "18_KM_Curves_All_Cohorts.pdf")

    print()
    print("=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"  Regions analysed : {len(regions)}")
    print(f"  Haplotypes tested: {freq_df['Haplotype_ID'].nunique() if not freq_df.empty else 0}")
    print(f"  Results          : {out_xlsx}")
    print(f"  Output folder    : {out_dir}")

if __name__ == "__main__":
    main()
