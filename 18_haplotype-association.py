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

2. MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx  — harmonised clinical master
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
1. Tumour vs Healthy       — Fisher's exact test per haplotype
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
    • tumour_vs_healthy
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
DEFAULT_MASTER  = Path("/home/gadeaalonsoj/tfm/MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx")
DEFAULT_ANNOT   = Path("/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_Annotated_Report_Fixed.xlsx")
DEFAULT_OUT     = Path("/home/gadeaalonsoj/tfm/gsdmb_final_results/")

# ── CONSTANTS ─────────────────────────────────────────────────────────────────
MIN_HAP_FREQ        = 0.02   # global frequency threshold — rare haplotypes skipped
MIN_CARRIERS        = 3      # minimum carriers (or non-carriers) to run a test
MIN_EVENTS_LOGISTIC = 5      # minimum events for adjusted logistic regression
FDR_THRESHOLD       = 0.10

# ── CLINICAL VARIABLE DEFINITIONS ────────────────────────────────────────────
# Identical to script 17 — haplotype carrier status replaces SNP carrier status
# as the exposure, but the outcome variables and test types are unchanged.

CLINICAL_VARS_BREAST: Dict[str, Dict] = {
    "canon__age":               {"type": "continuous", "label": "Age at diagnosis"},
    "BREAST_GRADE_NUMERIC":     {"type": "continuous", "label": "Tumour grade"},
    "BREAST_KI67_NUMERIC":      {"type": "continuous", "label": "KI67"},
    "canon__her2_copies":       {"type": "continuous", "label": "HER2 copies (FISH)"},
    "BREAST_OS_MONTHS_DERIVED": {"type": "continuous", "label": "Overall survival (months)"},
    "BREAST_P53_NUMERIC":       {"type": "continuous", "label": "p53 expression"},
    "BREAST_BMI_NUMERIC":       {"type": "continuous", "label": "BMI"},
    "BREAST_MENARCHE_NUMERIC":  {"type": "continuous", "label": "Age at menarche"},
    "BREAST_MENOPAUSE_NUMERIC": {"type": "continuous", "label": "Age at menopause"},
    "BREAST_ER_BIN":            {"type": "binary",     "label": "ER positive"},
    "BREAST_PR_BIN":            {"type": "binary",     "label": "PR positive"},
    "BREAST_RECURRENCE_DERIVED":{"type": "binary",     "label": "Recurrence / progression"},
    "BREAST_METASTASIS_DERIVED":{"type": "binary",     "label": "Distant metastasis"},
    "BREAST_LOCAL_MET_BIN":     {"type": "binary",     "label": "Local metastasis"},
    "BREAST_EXITUS_DERIVED":    {"type": "binary",     "label": "Exitus"},
    "BREAST_HER2_SUBTYPE":      {"type": "nominal",    "label": "HER2 subtype"},
    "BREAST_DX_TYPE":           {"type": "nominal",    "label": "Histological diagnosis type"},
}

CLINICAL_VARS_ENDO: Dict[str, Dict] = {
    "canon__age":                       {"type": "continuous", "label": "Age at surgery"},
    "ENDO_FIGO_NUMERIC":                {"type": "continuous", "label": "FIGO stage (ordinal)"},
    "ENDO_GRADE_NUMERIC":               {"type": "continuous", "label": "Tumour grade"},
    "ENDO_RISK_ORDINAL":                {"type": "continuous", "label": "Risk of recurrence (ordinal)"},
    "canon__os_months":                 {"type": "continuous", "label": "Overall survival (months)"},
    "canon__pfs_months":                {"type": "continuous", "label": "PFS (months)"},
    "ENDO_KI67_NUMERIC":                {"type": "continuous", "label": "KI67"},
    "ENDO_CTDNA_MAF_NUMERIC":           {"type": "continuous", "label": "ctDNA MAF"},
    "ENDO_CFDN_CONC_NUMERIC":           {"type": "continuous", "label": "cfDNA concentration"},
    "ENDO_PDL1_NUMERIC":                {"type": "continuous", "label": "PDL1"},
    "ENDO_CD8_NUMERIC":                 {"type": "continuous", "label": "CD8"},
    "ENDO_EXITUS_BIN":                  {"type": "binary",     "label": "Exitus"},
    "ENDO_PD_BIN":                      {"type": "binary",     "label": "Progression / PD"},
    "ENDO_PTEN_BIN":                    {"type": "binary",     "label": "PTEN loss"},
    "ENDO_MLH1_BIN":                    {"type": "binary",     "label": "MLH1 loss"},
    "ENDO_N_STAGE_BIN":                 {"type": "binary",     "label": "Node positive"},
    "ENDO_M_STAGE_BIN":                 {"type": "binary",     "label": "Distant metastasis"},
    "ENDO_LVSI_BIN":                    {"type": "binary",     "label": "LVSI"},
    "ENDO_MYOINVASION_BIN":             {"type": "binary",     "label": "Myometrial invasion >50%"},
    "canon__msi_status":                {"type": "nominal",    "label": "MSI status"},
    "canon__molecular_class":           {"type": "nominal",    "label": "Molecular classification"},
    "clin_au_endo__HISTOLOGY_GROUP":    {"type": "nominal",    "label": "Histology group"},
    "clin_au_endo__RISK_OF_RECURRENCE": {"type": "nominal",    "label": "Risk of recurrence (category)"},
}

SURVIVAL_COHORTS = {
    "Endometrial": {
        "cohort_filter": "Endometri",
        "endpoints": {
            "OS":  {"t_col": "canon__os_months",        "ev_col": "ENDO_EXITUS_BIN"},
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
_COHORT_C = {"Breast": "#AD1457", "Endometrial": "#00695C"}

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
    return {
        "Test_Unadj":   "Mann-Whitney U",
        "N_Carriers":   len(c),
        "N_NonCarriers":len(nc),
        "Stat_Unadj":   round(stat, 3),
        "P_Unadj":      p,
        "Median_Carriers":    round(c.median(), 3),
        "Median_NonCarriers": round(nc.median(), 3),
        "OR_Adj_Age": np.nan, "OR_Adj_Age_CI95": np.nan,
        "P_Adj_Age":  np.nan, "N_Adj_Age": np.nan,
        "OR_Adj_AgeBMI": np.nan, "OR_Adj_AgeBMI_CI95": np.nan,
        "P_Adj_AgeBMI":  np.nan, "N_Adj_AgeBMI": np.nan,
    }

def _test_binary(c_df, nc_df, col, age_col, bmi_col) -> Optional[Dict]:
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
    chi2, p, _, _ = stats.chi2_contingency(ct)
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


def load_phased_genotypes(phased_path: Path, keep_positions: List[int]) -> pd.DataFrame:
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

    # Keep only established SNP positions
    geno["POS"] = pd.to_numeric(geno["POS"], errors="coerce")
    geno = geno[geno["POS"].isin(keep_positions)].copy()
    geno = geno.sort_values("POS").reset_index(drop=True)

    n_snps = len(geno)
    print(f"  {n_snps} established SNPs retained in phased matrix")
    if n_snps == 0:
        raise ValueError(
            "No SNPs remain after filtering phased genotypes to established positions.\n"
            "Check that POS values in the TSV match those in the annotated report."
        )

    snp_labels = (geno["POS"].astype(str) + "_" +
                  geno["REF"].astype(str) + ">" +
                  geno["ALT"].astype(str)).tolist()

    sample_cols = [c for c in geno.columns if c not in ("CHROM", "POS", "REF", "ALT")]

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
    carrier_records = {}
    for hap_id, row in freq_df.iterrows():
        hap_str = row["Haplotype"]
        hap_id_label = row["Haplotype_ID"]
        carriers = []
        for samp in hap_df.index:
            h1, h2 = hap_df.loc[samp, "hap1"], hap_df.loc[samp, "hap2"]
            dosage = int(h1 == hap_str and "N" not in h1) + \
                     int(h2 == hap_str and "N" not in h2)
            carriers.append({"Sample_phased": samp,
                              "Haplotype_ID":  hap_id_label,
                              "Haplotype":     hap_str,
                              "Dosage":        dosage,
                              "Carrier":       int(dosage >= 1)})
    carrier_df = pd.DataFrame(carriers if carriers else
                               [{"Sample_phased": s, "Haplotype_ID": "H1",
                                 "Haplotype": "", "Dosage": 0, "Carrier": 0}
                                for s in hap_df.index])

    # Rebuild properly
    rows = []
    for _, freq_row in freq_df.iterrows():
        hap_str = freq_row["Haplotype"]
        hap_id  = freq_row["Haplotype_ID"]
        for samp in hap_df.index:
            h1 = hap_df.loc[samp, "hap1"]
            h2 = hap_df.loc[samp, "hap2"]
            dosage = (int(h1 == hap_str) if "N" not in h1 else 0) + \
                     (int(h2 == hap_str) if "N" not in h2 else 0)
            rows.append({
                "Sample_phased": samp,
                "Haplotype_ID":  hap_id,
                "Haplotype":     hap_str,
                "Global_Freq":   freq_row["Global_Freq"],
                "Dosage":        dosage,
                "Carrier":       int(dosage >= 1),
            })
    carrier_df = pd.DataFrame(rows)
    return freq_df, carrier_df


def load_clinical_master(master_path: Path) -> pd.DataFrame:
    """
    Load the harmonised clinical master.  This is the same Excel loaded by
    script 17 — it contains snp_code, sheet, Cohort, Tissue, is_replicate,
    and all derived clinical columns.

    Returns one row per unique sample (deduplicated on snp_code).
    """
    print("  Loading harmonised clinical master …")
    # The master may be multi-sheet; load the first (or 'Sheet1') sheet
    try:
        master = pd.read_excel(master_path, sheet_name=0)
    except Exception:
        master = pd.read_excel(master_path)

    # Normalise key column names (strip whitespace)
    master.columns = [str(c).strip() for c in master.columns]

    # Identify the snp_code column (may also be called 'Sample' or 'snp_code')
    snp_col = None
    for candidate in ["snp_code", "SNP_code", "Sample", "SAMPLE"]:
        if candidate in master.columns:
            snp_col = candidate
            break
    if snp_col is None:
        raise ValueError(
            "Cannot find snp_code / Sample column in clinical master.\n"
            f"Available columns: {list(master.columns[:20])}"
        )
    if snp_col != "snp_code":
        master = master.rename(columns={snp_col: "snp_code"})

    master["snp_code"] = master["snp_code"].astype(str).str.strip().str.replace(" ", "")

    # Ensure is_replicate exists
    if "is_replicate" not in master.columns:
        master["is_replicate"] = False

    # Deduplicate to one row per snp_code
    master = master.drop_duplicates("snp_code").copy()
    print(f"  Master: {len(master)} unique samples")
    return master


def merge_haplotypes_with_clinical(carrier_df: pd.DataFrame,
                                   master: pd.DataFrame) -> pd.DataFrame:
    """
    Join carrier_df (long format: Sample_phased × Haplotype_ID) to the
    clinical master on snp_code.

    The sample names in the phased TSV come from the VCF pipeline (folder
    names). These should match the snp_code values in the master after
    normalisation (strip spaces, uppercase).

    Returns a long-format DataFrame with one row per sample × haplotype,
    containing all clinical variables from the master.
    """
    print("  Joining haplotype carrier data to clinical master …")

    # Normalise join key on both sides
    carrier_df = carrier_df.copy()
    carrier_df["snp_code"] = (carrier_df["Sample_phased"]
                               .astype(str).str.strip().str.replace(" ", "").str.upper())
    master_join = master.copy()
    master_join["snp_code"] = master_join["snp_code"].str.upper()

    merged = carrier_df.merge(master_join, on="snp_code", how="inner")

    n_phased  = carrier_df["snp_code"].nunique()
    n_matched = merged["snp_code"].nunique()
    print(f"  Phased samples: {n_phased} | matched to master: {n_matched} "
          f"| unmatched: {n_phased - n_matched}")

    if n_matched == 0:
        print("\n  WARNING: No samples matched. Check that sample names in the phased")
        print("  TSV correspond to snp_code values in the master (after normalisation).")

    return merged


# ── ANALYSIS 1: TUMOUR VS HEALTHY ─────────────────────────────────────────────

def tumour_vs_healthy(merged: pd.DataFrame) -> pd.DataFrame:
    print("=== Analysis 1: Tumour vs Healthy ===")
    df = merged[~merged["is_replicate"]].copy()

    rows = []
    for hap_id, hdf in df.groupby("Haplotype_ID"):
        hap_str  = hdf["Haplotype"].iloc[0]
        hap_freq = hdf["Global_Freq"].iloc[0]
        for cohort_label in ["Breast", "Endometrial"]:
            coh_filter = "Breast" if cohort_label == "Breast" else "Endometri"
            coh_df = hdf[hdf["Cohort"].str.contains(coh_filter, case=False, na=False)]
            tum_df = coh_df[coh_df["Tissue"] == "Tumour"]
            hlt_df = coh_df[coh_df["Tissue"] == "Healthy"]
            if tum_df.empty or hlt_df.empty:
                continue

            tum_samp = tum_df.drop_duplicates("snp_code")
            hlt_samp = hlt_df.drop_duplicates("snp_code")
            a = tum_samp["Carrier"].sum()
            b = len(tum_samp) - a
            c = hlt_samp["Carrier"].sum()
            d = len(hlt_samp) - c

            if a + c < MIN_CARRIERS:
                continue

            _, p = fisher_exact([[a, b], [c, d]])
            or_val, or_meth = _haldane_or(a, b, c, d)
            rows.append({
                "Cohort":          cohort_label,
                "Haplotype_ID":    hap_id,
                "Haplotype":       hap_str,
                "Global_Freq":     hap_freq,
                "N_Tumour":        len(tum_samp),
                "N_Healthy":       len(hlt_samp),
                "Carriers_Tumour": int(a),
                "Carriers_Healthy":int(c),
                "Freq_Tumour_%":   round(a / len(tum_samp) * 100, 2),
                "Freq_Healthy_%":  round(c / len(hlt_samp) * 100, 2),
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
        c_df  = sample_data[sample_data["Carrier"] == 1]
        nc_df = sample_data[sample_data["Carrier"] == 0]
        if len(c_df) < MIN_CARRIERS or len(nc_df) < MIN_CARRIERS:
            continue
        for col, meta in var_dict.items():
            if col not in sample_data.columns:
                continue
            if meta["type"] == "continuous":
                res = _test_continuous(c_df[col], nc_df[col])
            elif meta["type"] == "binary":
                res = _test_binary(c_df, nc_df, col,
                                   age_col="canon__age", bmi_col="canon__bmi")
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

    res_df = _apply_fdr(res_df, "P_Unadj", "FDR_Unadj")
    res_df["Nominal_Sig_Unadj"]     = res_df["P_Unadj"]      < 0.05
    res_df["FDR_Sig_Unadj"]         = res_df["FDR_Unadj"]    < FDR_THRESHOLD
    if "P_Adj_Age" in res_df.columns:
        res_df = _apply_fdr(res_df, "P_Adj_Age", "FDR_Adj_Age")
        res_df["Nominal_Sig_Adj_Age"] = res_df["P_Adj_Age"] < 0.05
    if "P_Adj_AgeBMI" in res_df.columns:
        res_df = _apply_fdr(res_df, "P_Adj_AgeBMI", "FDR_Adj_AgeBMI")
        res_df["Nominal_Sig_Adj_AgeBMI"] = res_df["P_Adj_AgeBMI"] < 0.05
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
                stat, p_trend = stats.kruskal(*non_empty)
            elif meta["type"] == "binary":
                ct_rows = []
                for grp in [g0, g1, g2]:
                    v = pd.to_numeric(grp[col], errors="coerce").dropna()
                    if len(v) < MIN_CARRIERS:
                        continue
                    ct_rows.append([int(v.sum()), len(v) - int(v.sum())])
                if len(ct_rows) < 2:
                    continue
                chi2, p_trend, _, _ = stats.chi2_contingency(ct_rows)
                stat = chi2
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

                carrier = s.loc[valid.index, "Carrier"].fillna(0).astype(int)
                n_c   = carrier.sum()
                n_nc  = (carrier == 0).sum()
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
                for grp_val, grp_label, col in [(0, "Non-carrier", "#1976D2"),
                                                 (1, "Carrier",     "#D32F2F")]:
                    mask = carrier == grp_val
                    if mask.sum() < 2:
                        continue
                    kmf = KaplanMeierFitter()
                    kmf.fit(T[mask], E[mask],
                            label=f"{grp_label} (n={mask.sum()})")
                    kmf.plot_survival_function(ax=ax, ci_show=True, color=col)
                ax.set_title(f"{hap_id} [{hap_str}]  —  {cohort_label} | {endpoint}\n"
                             f"Log-rank p={lr_p:.4f}",
                             fontsize=9)
                ax.set_xlabel(f"{endpoint} (months)")
                ax.set_ylabel("Survival probability")
                ax.legend(fontsize=8)
                _style_ax(ax)
                plt.tight_layout()
                km_pages.append(fig)

    res = pd.DataFrame(rows)
    if not res.empty:
        for cohort in res["Cohort"].unique():
            for ep in res["Endpoint"].unique():
                mask = (res["Cohort"] == cohort) & (res["Endpoint"] == ep)
                if mask.sum() > 1:
                    res = _apply_fdr(res[mask], "P_Cox", "FDR_Cox").combine_first(res)
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

        case_samp    = case_df.drop_duplicates("snp_code").set_index("snp_code")
        control_samp = control_df.drop_duplicates("snp_code").set_index("snp_code")
        n_cases    = len(case_samp)
        n_controls = len(control_samp)
        print(f"  {cohort_label}: {n_cases} cases, {n_controls} controls")

        if n_cases < MIN_CARRIERS or n_controls < MIN_CARRIERS:
            print("    Too few samples — skipping.")
            continue

        for hap_id in df["Haplotype_ID"].unique():
            hap_str  = df[df["Haplotype_ID"] == hap_id]["Haplotype"].iloc[0]
            hap_freq = df[df["Haplotype_ID"] == hap_id]["Global_Freq"].iloc[0]

            c_case    = case_samp[case_samp.index.isin(
                case_df[case_df["Haplotype_ID"] == hap_id]["snp_code"])]
            c_control = control_samp[control_samp.index.isin(
                control_df[control_df["Haplotype_ID"] == hap_id]["snp_code"])]

            a = int(case_df[(case_df["Haplotype_ID"] == hap_id) &
                            (case_df["Carrier"] == 1)]["snp_code"].nunique())
            b = n_cases - a
            c = int(control_df[(control_df["Haplotype_ID"] == hap_id) &
                               (control_df["Carrier"] == 1)]["snp_code"].nunique())
            d = n_controls - c

            if a + c < MIN_CARRIERS:
                continue

            _, p_unadj = fisher_exact([[a, b], [c, d]])
            or_unadj, or_meth = _haldane_or(a, b, c, d)

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
                "OR_Method":         or_meth,
                "P_Unadj":           p_unadj,
            }

            # Age-adjusted logistic
            all_samp = pd.concat([
                case_samp[["canon__age", "canon__bmi"]].assign(
                    cancer=1,
                    carrier=case_samp.index.map(
                        lambda s: 1 if s in
                        set(case_df[case_df["Haplotype_ID"] == hap_id]["snp_code"])
                        else 0)),
                control_samp[["canon__age", "canon__bmi"]].assign(
                    cancer=0,
                    carrier=control_samp.index.map(
                        lambda s: 1 if s in
                        set(control_df[control_df["Haplotype_ID"] == hap_id]["snp_code"])
                        else 0)),
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
                   c=col, label=cohort, alpha=0.7, s=60, zorder=3)
        for _, row in sub[sub["Nominal_Sig"]].iterrows():
            ax.annotate(row["Haplotype_ID"],
                        (np.log2(row["OR"]), -np.log10(row["P_Value"])),
                        fontsize=7, ha="left", va="bottom")
    ax.axhline(-np.log10(0.05), ls="--", c="#aaaaaa", lw=1, label="p=0.05")
    ax.axvline(0, ls=":", c="#cccccc", lw=0.8)
    ax.set_xlabel("log₂(OR)  —  tumour vs healthy")
    ax.set_ylabel("-log₁₀(p)")
    ax.set_title("Haplotype: Tumour vs Healthy", fontweight="bold")
    ax.legend(fontsize=9)
    _style_ax(ax)
    plt.tight_layout()
    out = out_dir / "18_Haplo_Volcano_raw_p.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def make_heatmap(clin_res: pd.DataFrame, cohort_label: str,
                 p_col: str, title_suffix: str, out_path: Path):
    if clin_res.empty:
        return
    sub = clin_res[clin_res[p_col] < 0.1].copy()
    if sub.empty:
        print(f"  Heatmap ({title_suffix}): no hits below 0.1 — skipped")
        return
    pivot = sub.pivot_table(index="Haplotype_ID", columns="Clin_Label",
                            values=p_col, aggfunc="min")
    fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns) * 0.8),
                                    max(4, len(pivot) * 0.5)))
    sns.heatmap(-np.log10(pivot.fillna(1)), ax=ax, cmap="YlOrRd",
                annot=True, fmt=".1f", linewidths=0.4, cbar_kws={"label": "-log₁₀(p)"})
    ax.set_title(f"Haplotype × Clinical  —  {cohort_label}  ({title_suffix})",
                 fontweight="bold")
    ax.set_xlabel(""); ax.set_ylabel("Haplotype")
    plt.tight_layout()
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
    ax.set_title(f"Haplotype Forest  —  {cohort_label}", fontweight="bold")
    _style_ax(ax, grid=False)
    ax.xaxis.grid(True, linestyle=":", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def make_risk_forest(risk_res: pd.DataFrame, out_path: Path):
    if risk_res.empty or "OR_Unadj" not in risk_res.columns:
        return
    fig, ax = plt.subplots(figsize=(7, max(4, len(risk_res) * 0.45)))
    ys = range(len(risk_res))
    colours = [_COHORT_C.get(r["Cohort"], "#555") for _, r in risk_res.iterrows()]
    ax.scatter(risk_res["OR_Unadj"], ys, color=colours, zorder=3, s=50)
    ax.axvline(1, ls="--", c="#aaaaaa", lw=1)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([f"{r['Haplotype_ID']} ({r['Cohort']})"
                        for _, r in risk_res.iterrows()], fontsize=8)
    ax.set_xlabel("Odds Ratio (unadjusted)")
    ax.set_title("Haplotype Cancer Risk  —  Case-Control", fontweight="bold")
    _style_ax(ax, grid=False)
    ax.xaxis.grid(True, linestyle=":", alpha=0.4)
    handles = [mpatches.Patch(color=c, label=l) for l, c in _COHORT_C.items()]
    ax.legend(handles=handles, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def make_km_pdf(km_pages: List, out_path: Path):
    if not km_pages:
        return
    with pdf_backend.PdfPages(out_path) as pdf:
        for fig in km_pages:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"  Saved: {out_path}  ({len(km_pages)} pages)")


# ── MAIN ─────────────────────────────────────────────────────────────────────

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
    return p.parse_args()


def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_xlsx = out_dir / "18_Haplo_Clinical_Association_Results.xlsx"

    print("=" * 60)
    print("SCRIPT 18: HAPLOTYPE–CLINICAL ASSOCIATION ANALYSIS")
    print("=" * 60 + "\n")

    # ── Load data ────────────────────────────────────────────────────────
    print("[1/5] Loading established SNP positions …")
    snp_positions = load_established_snp_positions(Path(args.annot))

    print("\n[2/5] Loading and processing phased genotypes …")
    hap_df, snp_labels = load_phased_genotypes(Path(args.phased), snp_positions)

    print(f"\n  SNP positions used in haplotype backbone:")
    for i, lbl in enumerate(snp_labels, 1):
        print(f"    Position {i}: {lbl}")

    print("\n[3/5] Enumerating haplotypes …")
    freq_df, carrier_df = enumerate_haplotypes(hap_df, min_freq=args.min_freq)

    print(f"\n  Haplotypes to be tested (freq ≥ {args.min_freq}):")
    for _, row in freq_df.iterrows():
        print(f"    {row['Haplotype_ID']}: {row['Haplotype']}  "
              f"(freq={row['Global_Freq']:.3f})")

    print("\n[4/5] Loading clinical master …")
    master = load_clinical_master(Path(args.master))

    print("\n[5/5] Merging haplotype data with clinical master …")
    merged = merge_haplotypes_with_clinical(carrier_df, master)

    if merged.empty:
        print("\nERROR: Merge produced no rows. Cannot proceed.")
        print("Common causes:")
        print("  • Sample names in phased TSV don't match snp_code in master")
        print("  • The master sheet loaded is not the variant-level sheet")
        return

    n_samples   = merged["snp_code"].nunique()
    n_haplotypes = merged["Haplotype_ID"].nunique()
    print(f"\n  Merged: {n_samples} samples × {n_haplotypes} haplotypes\n")

    # ── Analyses ────────────────────────────────────────────────────────
    tvh       = tumour_vs_healthy(merged)
    b_clin, e_clin = clinical_associations(merged)

    def _top_haps(res):
        if res.empty or "Nominal_Sig_Unadj" not in res.columns:
            return []
        return res[res["Nominal_Sig_Unadj"]]["Haplotype_ID"].unique().tolist()

    b_dose = haplotype_dose_analysis(merged, _top_haps(b_clin),
                                     CLINICAL_VARS_BREAST, "Breast_Tumour")
    e_dose = haplotype_dose_analysis(merged, _top_haps(e_clin),
                                     CLINICAL_VARS_ENDO,   "Endometrial_Tumour")

    surv_res, km_pages = survival_analysis(merged)
    risk_res = cancer_risk_analysis(merged)

    # ── Write Excel ──────────────────────────────────────────────────────
    print("=== Writing output Excel ===")
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        freq_df.to_excel(xw,      sheet_name="haplotype_frequencies",   index=False)
        if not tvh.empty:
            tvh.to_excel(xw,      sheet_name="tumour_vs_healthy",        index=False)
        if not b_clin.empty:
            b_clin.to_excel(xw,   sheet_name="breast_clinical_assoc",   index=False)
        if not e_clin.empty:
            e_clin.to_excel(xw,   sheet_name="endo_clinical_assoc",     index=False)
        dose_all = pd.concat([b_dose, e_dose], ignore_index=True)
        if not dose_all.empty:
            dose_all.to_excel(xw, sheet_name="haplotype_dose",          index=False)
        if not surv_res.empty:
            surv_res.to_excel(xw, sheet_name="survival_cox",            index=False)
        if not risk_res.empty:
            risk_res.to_excel(xw, sheet_name="cancer_risk",             index=False)

        # Summary of significant results
        sig_rows = []
        for res_df, analysis in [(tvh, "Tumour_vs_Healthy"),
                                  (b_clin, "Breast_Clinical"),
                                  (e_clin, "Endo_Clinical"),
                                  (surv_res, "Survival"),
                                  (risk_res, "Cancer_Risk")]:
            if res_df.empty:
                continue
            p_col = ("P_Value" if "P_Value" in res_df.columns else
                     "P_Unadj" if "P_Unadj" in res_df.columns else
                     "P_Cox"   if "P_Cox"   in res_df.columns else None)
            if p_col:
                sig = res_df[res_df[p_col] < 0.05].copy()
                sig["Analysis_Type"] = analysis
                sig_rows.append(sig)
        summary = pd.concat(sig_rows, ignore_index=True) if sig_rows else pd.DataFrame()
        if not summary.empty:
            summary.to_excel(xw, sheet_name="summary_significant", index=False)

        # Sample manifest
        manifest_cols = [c for c in ["snp_code", "Cohort", "Tissue", "sheet",
                                      "is_replicate", "Haplotype_ID", "Haplotype",
                                      "Dosage", "Carrier", "canon__age",
                                      "canon__os_months", "canon__pfs_months",
                                      "ENDO_EXITUS_BIN", "BREAST_OS_MONTHS_DERIVED",
                                      "BREAST_EXITUS_DERIVED"]
                         if c in merged.columns]
        manifest = (merged[manifest_cols]
                    .drop_duplicates(subset=["snp_code", "Haplotype_ID"])
                    .sort_values(["Cohort", "Tissue", "snp_code"])
                    .reset_index(drop=True))
        manifest.to_excel(xw, sheet_name="sample_manifest", index=False)

    print(f"  Saved: {out_xlsx}\n")

    # ── Plots ────────────────────────────────────────────────────────────
    print("=== Generating plots ===")
    make_volcano(tvh, out_dir)
    make_heatmap(b_clin, "Breast",      "P_Unadj",    "raw p",
                 out_dir / "18_Haplo_Heatmap_Breast.png")
    make_heatmap(e_clin, "Endometrial", "P_Unadj",    "raw p",
                 out_dir / "18_Haplo_Heatmap_Endometrial.png")
    make_heatmap(e_clin, "Endometrial", "FDR_Unadj",  "FDR",
                 out_dir / "18_Haplo_Heatmap_Endometrial_FDR.png")
    make_forest(b_clin, "Breast_Tumour",
                out_dir / "18_Haplo_Forest_Breast.png")
    make_forest(e_clin, "Endometrial_Tumour",
                out_dir / "18_Haplo_Forest_Endometrial.png")
    make_risk_forest(risk_res,
                     out_dir / "18_Haplo_Forest_CancerRisk.png")
    make_km_pdf(km_pages, out_dir / "18_KM_Curves_All_Cohorts.pdf")

    # ── Final summary ────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"  Samples analysed  : {n_samples}")
    print(f"  Haplotypes tested : {n_haplotypes}")
    if not tvh.empty:
        print(f"  Tumour vs healthy — nominal: {tvh['Nominal_Sig'].sum()} | "
              f"FDR<{FDR_THRESHOLD}: {tvh['FDR_Sig'].sum()}")
    if not b_clin.empty:
        print(f"  Breast clinical   — nominal unadj: {b_clin['Nominal_Sig_Unadj'].sum()} | "
              f"FDR: {b_clin['FDR_Sig_Unadj'].sum()}")
    if not e_clin.empty:
        print(f"  Endo clinical     — nominal unadj: {e_clin['Nominal_Sig_Unadj'].sum()} | "
              f"FDR: {e_clin['FDR_Sig_Unadj'].sum()}")
    if not surv_res.empty:
        print(f"  Survival (Cox)    — {len(surv_res)} tests | "
              f"{surv_res['Nominal_Sig'].sum()} nominal")
    if not risk_res.empty:
        print(f"  Cancer risk       — "
              f"{risk_res.get('Nominal_Sig_Unadj', pd.Series()).sum()} nominal (unadj)")
    print(f"\n  Results : {out_xlsx}")
    print(f"  Plots   : {out_dir}")


if __name__ == "__main__":
    main()