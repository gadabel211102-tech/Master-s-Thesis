#!/usr/bin/env python3
"""Stage 23: integrate RNA measurements with the significant DNA findings.

This script is the final synthesis layer for thesis objective 2. It collects:

- the cleaned Excel-based isoform measurements from stage 20
- BAM-derived panel-gene ratios from stage 20b
- the significant SNP backbone from stage 12
- the significant rare-variant / mutation exposures from stage 17b
- the significant haplotype signals from stage 22, or stage 15 as a fallback

The output is an integration workbook that shows which RNA variables can be
tested against which significant SNP or haplotype exposures, together with
summary plots and manifest tables that document sample matching and QC status.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sample_identity_utils import build_analysis_sample_map
from pipeline_utils import extract_rsid, find_col

try:
    import seaborn as sns
except ImportError:
    sns = None

ROOT = Path(__file__).resolve().parent
RES = ROOT / "analysis_results"
OUT_DIR = RES / "23_rna_integration"

# Centralising path definitions at the top makes the final integration stage
# easier to audit, because every upstream dependency is visible in one place.
PATHS = {
    "master": ROOT / "MASTER_SNP_plus_clinical_HARMONISED.xlsx",
    "rna_qc_manifest": RES / "19b_objective2_rna_qc" / "RNA_QC_Manifest.tsv",
    "stage20": RES / "20_isoform_expression_associations" / "GSDMB_Objective2_Isoform_Results.xlsx",
    "stage20b": RES / "20b_panel_gene_expression" / "20b_Panel_Gene_Expression.xlsx",
    "stage12": RES / "12_snp_enrichment" / "GSDMB_SNP_Enrichment_Results.xlsx",
    "stage17b_breast": RES / "17b_rare_variant_associations" / "17b_breast_mutation_assoc.tsv",
    "stage17b_endo": RES / "17b_rare_variant_associations" / "17b_endo_mutation_assoc.tsv",
    "stage17b_exposure_catalogue": RES / "17b_rare_variant_associations" / "17b_exposure_catalogue.tsv",
    "stage15": RES / "15_haplotype_statistics" / "GSDMB_Haplotype_Results.xlsx",
    "phased": RES / "15_haplotype_phasing" / "phased_genotypes.tsv",
    "stage22": RES / "22_haplotype_first_interpretation" / "22_Haplotype_First_Interpretation.xlsx",
}

ISOFORM_COLS = [
    "G1", "G2", "G3", "G3b", "G4", "GSDMB",
    "GSDMB1_expression", "GSDMB2_expression", "GSDMB3_expression", "GSDMB4_expression",
    "exon6_positive_fraction", "exon6_positive_to_negative_log2",
    "exon7_containing_fraction", "exon7_deficient_fraction", "exon7_containing_to_deficient_log2",
    "GSDMB2_fraction", "GSDMB_to_ERBB2_ratio", "exon6_positive_GSDMB_to_ERBB2_ratio",
    "pyroptotic_isoform_fraction", "non_pyroptotic_isoform_fraction",
    "G4_vs_total", "G2_vs_total",
    "exon6_containing_fraction", "exon6_lacking_fraction", "exon6_balance_log2",
    "exon7_containing_fraction", "exon7_lacking_fraction", "exon7_balance_log2",
]
ISOFORM_BASE_COLS = ("G1", "G2", "G3", "G3b", "G4")
ISOFORM_PROFILE_VARS = tuple(f"excel_iso_ratio__{iso}_to_isoform_total" for iso in ISOFORM_BASE_COLS)
EXCEL_LEGACY_ENDPOINTS = {
    "pyroptotic_isoform_fraction", "non_pyroptotic_isoform_fraction",
    "G4_vs_total", "G2_vs_total",
    "exon6_containing_fraction", "exon6_lacking_fraction", "exon6_balance_log2",
    "exon7_lacking_fraction", "exon7_balance_log2",
}
QC_COLS = [
    "RNA_QC_Analysis_Ready", "RNA_QC_Exploratory_Ready", "RNA_QC_Final_Status",
    "RNA_QC_Exclusion_Reason", "RNA_QC_Eligibility_Note",
]
IMMUNE_SIGNATURE_GENES = ("PTPRC", "CD14", "CD68", "PDCD1", "CD274", "CD8A", "CD8B", "GZMB", "CXCL9", "CXCL10", "LAG3", "CTLA4", "TIGIT", "IFNG")
META_PREFIXES = ("canon__", "BREAST_", "ENDO_", "clin_")
META_COLS = {
    "snp_code", "sample_id", "case_id", "sheet", "analysis_group", "analysis_role",
    "cohort", "tumour_normal", "Tissue", "match_status", "match_method",
    "match_candidates", "master_join_success", "analysis_include_primary",
    "NOMBRE DE LA MUESTRA", "CODIGO JC", "raw1", "raw2", "sx1", "sx2",
    "jc_comments", "SNP28", "rs869402_x", "observaciones EVA", "TIENEN RNA",
    "FALTA MUESTRA", "isoform_total_scale", "isoform_total_sum", "isoform_total_flag",
    "endpoint_missing_count", "DNA_BAM_Group", "DNA_BAM_Path", "DNA_BAM_Name",
    "DNA_BAM_Size_MB", "DNA_Baseline_Present", "RNA_BAM_Group", "RNA_BAM_Path",
    "RNA_BAM_Name", "RNA_BAM_Size_MB", "RNA_BAM_Present", "RNA_BAM_Has_Index",
    "RNA_BAM_Group_Matched", "RNA_Machine_Total_Reads", "RNA_Machine_Mapped_Reads",
    "RNA_Machine_Mean_Read_Length", "RNA_Machine_Q20_Percent", "RNA_Machine_Path",
    "RNA_Machine_Name", "RNA_Machine_QC_Status", "RNA_Machine_QC_Note",
    "RNA_Machine_Source", "RNA_Target_Coverage_Available", "RNA_Target_Coverage_Targets",
    "RNA_Target_Coverage_Median_Depth", "RNA_Target_Coverage_Targets_20x_Pct",
    "RNA_Target_Coverage_Targets_100x_Pct", "RNA_Target_Coverage_GSDMB_Mean_Depth",
    "RNA_Target_Coverage_GSDMB_Fraction_Covered", "RNA_Target_Coverage_Waning",
    "RNA_Target_Coverage_QC_Status", "RNA_Gene_Coverage_Targets_75Reads_Pct",
    "RNA_Gene_Coverage_Targets_100Reads_Pct", "RNA_Gene_Coverage_Median_Reads_Per_Gene",
    "RNA_Gene_Coverage_GSDMB_Reads", "RNA_Gene_Coverage_Housekeeping_Reads_Median",
    "RNA_Gene_Coverage_Waning",
} | set(QC_COLS)
FDR_THRESHOLD = 0.10
MIN_N = 8
MIN_CARRIERS = 3


# Small coercion helpers are used throughout because the upstream RNA and
# clinical workbooks contain mixed formatting and legacy sample labels.
def s(v):
    return "" if pd.isna(v) else str(v).strip()


def sx(v):
    t = s(v)
    patterns = [
        (r"^(SNP_(?:AT|EN|MN|MT-T)_\d+)", lambda m: m.group(1).upper()),
        (r"^SNP_DNA_(EN|MN|AT)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_DNA_MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^DNA_SNP_(EN|MN|AT)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^DNA_(AT|EN|MN)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^DNA_(?:SNP_)?MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^MAMAH2_MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^MAMAH2_MT[-_]N_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^RNA_SNP_(AT|EN|MN|MT-T|MT-N)[-_]?(\d+)", lambda m: f"SNP_{'MT-T' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_RNA_(AT|EN|MN|MT-T|MT-N)_(\d+)", lambda m: f"SNP_{'MT-T' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_(AT|EN|MN|MT-T|MT-N)_RNA_(\d+)", lambda m: f"SNP_{'MT-T' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^RNA_(AT|EN|MN)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^RNA_MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^RNA_MT[-_]N_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
    ]
    for pat, fmt in patterns:
        m = re.match(pat, t, re.I)
        if m:
            return fmt(m)
    return None


def first_available_column(df, candidates):
    for candidate in candidates:
        col = find_col(df, candidate)
        if col is not None:
            return col
    return None


def gt_dose(gt):
    t = s(gt)
    if t in {"", ".", "./.", ".|."}:
        return np.nan
    parts = t.split("|") if "|" in t else t.split("/") if "/" in t else None
    try:
        return float(sum(int(x) for x in parts)) if parts else np.nan
    except Exception:
        return np.nan


def gt_haps(gt):
    t = s(gt)
    if t in {"", ".", "./.", ".|."}:
        return "N", "N"
    if "|" in t:
        return tuple(t.split("|", 1))
    if "/" in t:
        return tuple(t.split("/", 1))
    return "N", "N"


def genomic_snp_id(chrom, pos, ref, alt):
    chrom_s = s(chrom)
    if chrom_s and not chrom_s.lower().startswith("chr"):
        chrom_s = f"chr{chrom_s}"
    pos_s = s(pos)
    ref_s = s(ref).upper()
    alt_s = s(alt).upper()
    if not chrom_s or not pos_s or not ref_s or not alt_s:
        return np.nan
    return f"{chrom_s}:{pos_s}:{ref_s}:{alt_s}"


def analysis_id_to_snp_code(value):
    sample = s(value)
    if not sample:
        return None
    sample = re.sub(r"__RUN\d+$", "", sample)
    code = sx(sample)
    return code or sample


def first_valid(series):
    for v in series:
        if pd.notna(v) and str(v).strip() not in {"", "nan", "None"}:
            return v
    return np.nan


def bh(vals):
    arr = np.asarray(vals, dtype=float)
    out = np.full(arr.shape, np.nan)
    mask = np.isfinite(arr)
    if not mask.any():
        return out
    p = arr[mask]
    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)
    adj = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, ranked[i] * n / (i + 1))
        adj[i] = min(prev, 1.0)
    back = np.empty(n)
    back[order] = adj
    out[mask] = back
    return out


def add_fdr(df, p_col, group_cols):
    if df.empty:
        return df
    df = df.copy()
    df["FDR"] = np.nan
    for _, idx in df.groupby(group_cols, dropna=False).groups.items():
        df.loc[idx, "FDR"] = bh(pd.to_numeric(df.loc[idx, p_col], errors="coerce"))
    df["Nominal_Sig"] = pd.to_numeric(df[p_col], errors="coerce") < 0.05
    df["FDR_Sig"] = pd.to_numeric(df["FDR"], errors="coerce") < FDR_THRESHOLD
    return df


def clean_name(name):
    name = re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")
    return re.sub(r"_+", "_", name)


def categorise_excel_column(col):
    if col in ISOFORM_COLS:
        if col == "GSDMB":
            return "excel_total", f"excel_total__{clean_name(col)}", "overall_gene_expression"
        if col in {"GSDMB1_expression", "GSDMB2_expression", "GSDMB3_expression", "GSDMB4_expression"}:
            return "excel_isoform", f"excel_iso__{clean_name(col)}", "isoform_expression"
        if col in {"exon6_positive_fraction", "exon7_containing_fraction", "exon7_deficient_fraction", "GSDMB2_fraction", "exon6_containing_fraction", "exon6_lacking_fraction", "G2_vs_total", "G4_vs_total", "pyroptotic_isoform_fraction", "non_pyroptotic_isoform_fraction"}:
            return "excel_iso_ratio", f"excel_iso_ratio__{clean_name(col)}", "isoform_composition"
        if col in {"exon6_positive_to_negative_log2", "exon7_containing_to_deficient_log2", "exon6_balance_log2", "exon7_balance_log2"}:
            return "excel_iso_ratio", f"excel_iso_ratio__{clean_name(col)}", "isoform_log_ratio"
        if col in {"GSDMB_to_ERBB2_ratio", "exon6_positive_GSDMB_to_ERBB2_ratio"}:
            return "excel_erbb2_ratio", f"excel_erbb2_ratio__{clean_name(col)}", "erbb2_normalised_expression"
        return "excel_isoform", f"excel_iso__{clean_name(col)}", "isoform_expression"
    if col == "GSDMB EXPRESSION IN PANEL":
        return "excel_panel", f"excel_panel__{clean_name(col)}", "other_rna_quant"
    return "excel_panel", f"excel_panel__{clean_name(col)}", "other_rna_quant"


def expression_nonmissing(df, cols):
    if not cols:
        return pd.Series(0, index=df.index)
    return df[cols].apply(pd.to_numeric, errors="coerce").notna().sum(axis=1)


def _safe_log2_ratio(numerator, denominator, pseudocount=1e-6):
    numerator = pd.to_numeric(numerator, errors="coerce")
    denominator = pd.to_numeric(denominator, errors="coerce")
    return np.log2((numerator + pseudocount) / (denominator + pseudocount))


def add_isoform_ratio_features(collapsed, var_rows):
    iso_cols = {iso: f"excel_iso__{clean_name(iso)}" for iso in ISOFORM_BASE_COLS}
    available = {iso: col for iso, col in iso_cols.items() if col in collapsed.columns}
    if len(available) < 2:
        return collapsed, var_rows

    iso_values = pd.DataFrame({
        iso: pd.to_numeric(collapsed[col], errors="coerce")
        for iso, col in available.items()
    }, index=collapsed.index)
    total = iso_values.sum(axis=1, min_count=1)
    valid_total = total > 0
    existing_vars = {row.get("Integrated_Var") for row in var_rows}

    def record_var(new_col, label, source, level):
        if new_col in existing_vars:
            return
        var_rows.append({
            "Integrated_Var": new_col,
            "Original_Var": label,
            "RNA_Source": source,
            "RNA_Level": level,
            "Nonmissing_Collapsed": int(pd.to_numeric(collapsed[new_col], errors="coerce").notna().sum()),
        })
        existing_vars.add(new_col)

    def add_numeric_feature(new_col, label, values, source, level):
        if new_col not in collapsed.columns:
            collapsed[new_col] = values
        record_var(new_col, label, source, level)

    iso_aliases = {
        "GSDMB1_expression": iso_values["G1"] if "G1" in iso_values.columns else pd.Series(np.nan, index=collapsed.index),
        "GSDMB2_expression": iso_values["G2"] if "G2" in iso_values.columns else pd.Series(np.nan, index=collapsed.index),
        "GSDMB3_expression": iso_values[[c for c in ("G3", "G3b") if c in iso_values.columns]].sum(axis=1, min_count=1),
        "GSDMB4_expression": iso_values["G4"] if "G4" in iso_values.columns else pd.Series(np.nan, index=collapsed.index),
    }
    for label, values in iso_aliases.items():
        add_numeric_feature(
            f"excel_iso__{label}",
            label.replace("_", " "),
            values,
            "excel_isoform",
            "isoform_expression",
        )

    for iso in ISOFORM_BASE_COLS:
        if iso not in iso_values.columns:
            continue
        new_col = f"excel_iso_ratio__{iso}_to_isoform_total"
        collapsed[new_col] = np.where(valid_total, iso_values[iso] / total, np.nan)

    def add_log_ratio(new_col, label, numerator, denominator):
        add_numeric_feature(new_col, label, _safe_log2_ratio(numerator, denominator), "excel_iso_ratio", "isoform_log_ratio")

    pyro = iso_values[[c for c in ("G3", "G3b", "G4") if c in iso_values.columns]].sum(axis=1, min_count=1)
    non_pyro = iso_values[[c for c in ("G1", "G2") if c in iso_values.columns]].sum(axis=1, min_count=1)
    if pyro.notna().any():
        add_numeric_feature(
            "excel_iso_ratio__exon6_positive_fraction",
            "(G3+G3b+G4) / sum(G1,G2,G3,G3b,G4)",
            np.where(valid_total, pyro / total, np.nan),
            "excel_iso_ratio",
            "isoform_composition",
        )
        if "G2" in iso_values.columns:
            add_numeric_feature(
                "excel_iso_ratio__GSDMB2_fraction",
                "GSDMB2 / sum(G1,G2,G3,G3b,G4)",
                np.where(valid_total, iso_values["G2"] / total, np.nan),
                "excel_iso_ratio",
                "isoform_composition",
            )
    if pyro.notna().any() and non_pyro.notna().any():
        add_log_ratio("excel_iso_ratio__exon6_positive_to_negative_log2", "log2((G3+G3b+G4)/(G1+G2))", pyro, non_pyro)
    exon7_containing = iso_values[[c for c in ("G1", "G3", "G3b") if c in iso_values.columns]].sum(axis=1, min_count=1)
    exon7_deficient = iso_values[[c for c in ("G2", "G4") if c in iso_values.columns]].sum(axis=1, min_count=1)
    if exon7_containing.notna().any():
        add_numeric_feature(
            "excel_iso_ratio__exon7_containing_fraction",
            "(G1+G3+G3b) / sum(G1,G2,G3,G3b,G4)",
            np.where(valid_total, exon7_containing / total, np.nan),
            "excel_iso_ratio",
            "isoform_composition",
        )
    if exon7_deficient.notna().any():
        add_numeric_feature(
            "excel_iso_ratio__exon7_deficient_fraction",
            "(G2+G4) / sum(G1,G2,G3,G3b,G4)",
            np.where(valid_total, exon7_deficient / total, np.nan),
            "excel_iso_ratio",
            "isoform_composition",
        )
    if exon7_containing.notna().any() and exon7_deficient.notna().any():
        add_log_ratio("excel_iso_ratio__exon7_containing_to_deficient_log2", "log2((G1+G3+G3b)/(G2+G4))", exon7_containing, exon7_deficient)
    erbb2_col = "excel_panel__ERBB2" if "excel_panel__ERBB2" in collapsed.columns else None
    gsdmb_col = "excel_total__GSDMB" if "excel_total__GSDMB" in collapsed.columns else None
    if erbb2_col:
        erbb2 = pd.to_numeric(collapsed[erbb2_col], errors="coerce")
        if gsdmb_col:
            add_numeric_feature(
                "excel_erbb2_ratio__GSDMB_to_ERBB2_ratio",
                "Total GSDMB / ERBB2",
                np.where(erbb2 > 0, pd.to_numeric(collapsed[gsdmb_col], errors="coerce") / erbb2, np.nan),
                "excel_erbb2_ratio",
                "erbb2_normalised_expression",
            )
        if pyro.notna().any():
            add_numeric_feature(
                "excel_erbb2_ratio__exon6_positive_GSDMB_to_ERBB2_ratio",
                "(G3+G3b+G4) / ERBB2",
                np.where(erbb2 > 0, pyro / erbb2, np.nan),
                "excel_erbb2_ratio",
                "erbb2_normalised_expression",
            )
    return collapsed, var_rows


def load_excel_na(stage20_path):
    # Stage 20 can contain repeated rows per sample. We collapse to one
    # best-supported row per `snp_code` while keeping provenance columns that
    # reveal whether multiple source rows were encountered.
    raw = pd.read_excel(stage20_path, sheet_name="expression_cleaned")
    raw["snp_code"] = raw["snp_code"].astype(str).str.strip()
    raw = raw[raw["snp_code"].notna() & (raw["snp_code"] != "") & (raw["snp_code"] != "nan")].copy()
    value_cols = []
    for col in raw.columns:
        if col in EXCEL_LEGACY_ENDPOINTS:
            continue
        if col in META_COLS or col.startswith(META_PREFIXES):
            continue
        if re.match(r"^rs\d+(_[A-Za-z]+)?$", str(col), re.I):
            continue
        if re.match(r"^\d+_[ACGT]+>[ACGT]+$", str(col), re.I):
            continue
        if re.match(r"^chr[\w.]+:\d+:[^:]+:[^:]+$", str(col), re.I):
            continue
        vals = pd.to_numeric(raw[col], errors="coerce")
        if vals.notna().sum() > 0:
            value_cols.append(col)
    raw["_expr_nonmissing"] = expression_nonmissing(raw, value_cols)
    raw["_matched_rank"] = (raw.get("match_status", "") == "matched").astype(int)
    raw["_join_rank"] = raw.get("master_join_success", False).fillna(False).astype(int)
    raw["_qc_rank"] = raw.get("RNA_QC_Exploratory_Ready", raw.get("RNA_QC_Analysis_Ready", False)).fillna(False).astype(int)
    raw = raw.sort_values(["snp_code", "_matched_rank", "_join_rank", "_qc_rank", "_expr_nonmissing"], ascending=[True, False, False, False, False])

    group_sizes = raw.groupby("snp_code").size().rename("excel_source_row_count")
    mixed_groups = raw.groupby("snp_code")["analysis_group"].nunique(dropna=True).rename("excel_source_group_count")
    mixed_names = raw.groupby("snp_code")["NOMBRE DE LA MUESTRA"].nunique(dropna=True).rename("excel_source_name_count")
    agg = {col: first_valid for col in raw.columns if not col.startswith("_")}
    collapsed = raw.groupby("snp_code", dropna=False).agg(agg).reset_index(drop=True)
    collapsed = collapsed.merge(group_sizes.reset_index(), on="snp_code", how="left")
    collapsed = collapsed.merge(mixed_groups.reset_index(), on="snp_code", how="left")
    collapsed = collapsed.merge(mixed_names.reset_index(), on="snp_code", how="left")
    collapsed["excel_mixed_group_sources"] = collapsed["excel_source_group_count"].fillna(0) > 1
    collapsed["excel_mixed_name_sources"] = collapsed["excel_source_name_count"].fillna(0) > 1

    rename_map = {}
    var_rows = []
    for col in value_cols:
        source, new_name, level = categorise_excel_column(col)
        rename_map[col] = new_name
        var_rows.append({
            "Integrated_Var": new_name,
            "Original_Var": col,
            "RNA_Source": source,
            "RNA_Level": level,
            "Nonmissing_Collapsed": int(pd.to_numeric(collapsed[col], errors="coerce").notna().sum()),
        })
    collapsed = collapsed.rename(columns=rename_map)
    collapsed, var_rows = add_isoform_ratio_features(collapsed, var_rows)
    vars_df = pd.DataFrame(var_rows).drop_duplicates("Integrated_Var", keep="first").sort_values(["RNA_Source", "Integrated_Var"]).reset_index(drop=True)
    return collapsed, vars_df, raw


def load_bam_na(stage20b_path):
    # BAM-derived measurements arrive already normalised; this loader mainly
    # standardises naming so the Excel and BAM layers can be merged cleanly.
    ratio = pd.read_excel(stage20b_path, sheet_name="normalised_ratio")
    ratio["snp_code"] = ratio["snp_code"].astype(str).str.strip()
    meta = [c for c in ["snp_code", "bam_path", "bam_group", "cohort", "tissue", "analysis_role"] if c in ratio.columns]
    gene_cols = [c for c in ratio.columns if c not in meta and pd.to_numeric(ratio[c], errors="coerce").notna().sum() > 0]
    rename_map = {c: f"bam_ratio__{clean_name(c)}" for c in gene_cols}
    out = ratio[meta + gene_cols].rename(columns=rename_map).copy()
    vars_df = pd.DataFrame([
        {
            "Integrated_Var": rename_map[c],
            "Original_Var": c,
            "RNA_Source": "bam_ratio",
            "RNA_Level": "overall_gene_expression" if c == "GSDMB" else "other_rna_quant",
            "Nonmissing_Collapsed": int(pd.to_numeric(ratio[c], errors="coerce").notna().sum()),
        }
        for c in gene_cols
    ]).sort_values(["RNA_Source", "Integrated_Var"]).reset_index(drop=True)
    return out, vars_df

def build_na_manifest(excel_df, bam_df, qc_df):
    # The manifest is a human-readable audit table: it shows which samples have
    # Excel RNA data, BAM-derived ratios, and QC approval.
    meta_cols = [c for c in [
        "snp_code", "sample_id", "case_id", "analysis_group", "analysis_role", "cohort",
        "tumour_normal", "Tissue", "match_status", "match_method", "master_join_success",
        "RNA_QC_Analysis_Ready", "RNA_QC_Exploratory_Ready", "RNA_QC_Final_Status",
        "RNA_QC_Exclusion_Reason", "RNA_QC_Eligibility_Note", "NOMBRE DE LA MUESTRA",
        "CODIGO JC", "excel_source_row_count", "excel_mixed_group_sources",
        "excel_mixed_name_sources",
    ] if c in excel_df.columns]
    manifest = excel_df[meta_cols].copy()
    if not qc_df.empty:
        keep = [c for c in ["snp_code"] + QC_COLS if c in qc_df.columns]
        gate = qc_df[keep].copy()
        gate["snp_code"] = gate["snp_code"].astype(str).str.strip()
        gate = gate[gate["snp_code"].notna() & (gate["snp_code"] != "") & (gate["snp_code"] != "nan")]
        gate = gate.sort_values([c for c in ["RNA_QC_Exploratory_Ready", "RNA_QC_Analysis_Ready"] if c in gate.columns], ascending=False).drop_duplicates("snp_code")
        manifest = manifest.drop(columns=[c for c in QC_COLS if c in manifest.columns], errors="ignore").merge(gate, on="snp_code", how="outer")
    if not bam_df.empty:
        bam_meta = bam_df[[c for c in ["snp_code", "bam_path", "bam_group", "cohort", "tissue", "analysis_role"] if c in bam_df.columns]].drop_duplicates("snp_code")
        manifest = manifest.merge(bam_meta, on="snp_code", how="outer", suffixes=("", "_bam"))
    manifest["has_excel_row"] = manifest["analysis_group"].notna()
    manifest["has_bam_ratio"] = manifest.get("bam_path", pd.Series(index=manifest.index, dtype=object)).notna()
    return manifest


def merge_na_layers(excel_df, bam_df):
    # Integration is keyed on `snp_code` so every downstream association table
    # can point back to the same sample identifier used in the thesis workbooks.
    merged = excel_df.copy()
    if not bam_df.empty:
        merged = merged.merge(bam_df, on="snp_code", how="outer", suffixes=("", "_bam"))
        for col in ["cohort", "analysis_role"]:
            bam_col = f"{col}_bam"
            if bam_col in merged.columns:
                merged[col] = merged[col].fillna(merged[bam_col])
        if "tissue" in merged.columns:
            merged["tumour_normal"] = merged.get("tumour_normal").fillna(merged["tissue"])
    return merged


def load_significant_snps(stage12_path):
    # Only the statistically supported SNPs are carried into stage 23. This
    # keeps the RNA layer interpretive rather than exploratory.
    parts = []
    for sheet in ["Global", "Breast", "Endometrium"]:
        df = pd.read_excel(stage12_path, sheet_name=sheet)
        column_map = {
            "FDR_P_Value": ["FDR_P_Value", "Carrier_FDR_P_Value", "Primary_Genotype_FDR_P_Value", "Best_Genotype_FDR"],
            "P_Value": ["P_Value", "Carrier_P_Value", "Primary_Genotype_P_Value", "Best_Genotype_P"],
            "LD_Block_ID": ["LD_Block_ID", "LD_Block_Assignment"],
            "Odds_Ratio": ["Odds_Ratio", "Carrier_Odds_Ratio", "Het_vs_WT_Odds_Ratio", "Hom_vs_WT_Odds_Ratio"],
        }
        for target, candidates in column_map.items():
            source = first_available_column(df, candidates)
            if source is not None and source != target:
                df[target] = df[source]
        missing = [name for name in ["FDR_P_Value", "P_Value"] if name not in df.columns]
        if missing or first_available_column(df, ["rsID", "SNP_rsID", "SNP_ID"]) is None:
            raise ValueError(
                f"Stage-12 {sheet} sheet is missing required columns after schema normalisation: {', '.join(missing or ['SNP_ID/rsID'])}"
            )
        if "LD_Block_ID" not in df.columns:
            df["LD_Block_ID"] = pd.NA
        if "Odds_Ratio" not in df.columns:
            df["Odds_Ratio"] = np.nan
        rsid_col = first_available_column(df, ["rsID", "SNP_rsID"])
        snp_id_col = first_available_column(df, ["SNP_ID"])
        if rsid_col is not None:
            df["SNP_rsID"] = df[rsid_col].map(extract_rsid).fillna(df[rsid_col])
        else:
            df["SNP_rsID"] = pd.NA
        if snp_id_col is not None:
            df["SNP_Technical_ID"] = df[snp_id_col]
        else:
            df["SNP_Technical_ID"] = pd.NA
        df["SNP_ID"] = df["SNP_rsID"].fillna(df["SNP_Technical_ID"])
        df["SNP_Match_Key"] = df["SNP_ID"]
        df["FDR_P_Value"] = pd.to_numeric(df["FDR_P_Value"], errors="coerce")
        df["P_Value"] = pd.to_numeric(df["P_Value"], errors="coerce")
        df["Odds_Ratio"] = pd.to_numeric(df["Odds_Ratio"], errors="coerce")
        if df.empty:
            continue
        df = df[df["FDR_P_Value"].fillna(1) < FDR_THRESHOLD].copy()
        if df.empty:
            continue
        df["Significant_In"] = sheet
        parts.append(df)
    if not parts:
        return pd.DataFrame(columns=["SNP_ID", "SNP_rsID", "SNP_Technical_ID", "SNP_Match_Key", "Source_Groups", "Min_FDR", "Min_P", "LD_Block_ID", "Min_OR"])
    all_sig = pd.concat(parts, ignore_index=True)
    return (
        all_sig.groupby("SNP_ID", dropna=False)
        .agg(
            SNP_rsID=("SNP_rsID", first_valid),
            SNP_Technical_ID=("SNP_Technical_ID", lambda x: "; ".join(sorted({s(v) for v in x if s(v)}))),
            SNP_Match_Key=("SNP_Match_Key", first_valid),
            Source_Groups=("Significant_In", lambda x: "; ".join(sorted(set(map(str, x))))),
            Min_FDR=("FDR_P_Value", "min"),
            Min_P=("P_Value", "min"),
            LD_Block_ID=("LD_Block_ID", lambda x: "; ".join(sorted(set(map(str, x))))),
            Min_OR=("Odds_Ratio", "min"),
        )
        .reset_index()
        .sort_values(["Min_FDR", "SNP_ID"])
        .reset_index(drop=True)
    )


def empty_significant_mutation_table():
    return pd.DataFrame(
        columns=[
            "Exposure_ID", "Base_Exposure_ID", "Comparison", "Exposure_Label",
            "Exposure_Type", "Gene", "Consequence", "IMPACT",
            "Component_Genes", "Component_Variants", "Component_Variant_Count",
            "Significant_Clinical_Vars", "Significant_Models", "Min_FDR",
            "N_Carriers_Cohort", "Carrier_Samples", "RNA_Testable",
            "RNA_Test_Note",
        ]
    )


def empty_significant_haplotype_table():
    return pd.DataFrame(
        columns=[
            "Region", "Comparison", "Haplotype_String", "Haplotype_Def", "Haplotype_Copies",
            "Global_Freq", "Freq_Tumour", "Freq_Healthy", "OR", "OR_95CI_Lo", "OR_95CI_Hi",
            "P_Fisher", "P_Model", "N_Tumour", "N_Healthy", "FDR", "Significant_Fisher",
            "Significant_FDR", "Haplotype_ID", "Exposure_ID", "Interpretation_Source",
        ]
    )


def load_significant_haplotypes_from_stage15(stage15_path):
    xl = pd.ExcelFile(stage15_path)
    frames = []
    pattern = re.compile(r"^LD_Block_r2_080_Block\d+_Assoc$")
    for sheet in xl.sheet_names:
        if sheet != "Full_Associations" and not pattern.match(sheet):
            continue
        df = pd.read_excel(stage15_path, sheet_name=sheet)
        if "Significant_FDR" not in df.columns:
            continue
        df = df[df["Significant_FDR"].fillna(False)].copy()
        if df.empty:
            continue
        df["Interpretation_Source"] = "Stage15_Association_Fallback"
        frames.append(df)
    if not frames:
        return empty_significant_haplotype_table()
    df = pd.concat(frames, ignore_index=True)
    df["Exposure_ID"] = df["Region"].astype(str) + "::" + df["Haplotype_ID"].astype(str)
    return df.sort_values(["FDR", "Region", "Haplotype_ID"]).reset_index(drop=True)


def load_significant_haplotypes(stage22_path, stage15_path):
    # Stage 22 is preferred because it contains the curated interpretation
    # layer, but stage 23 can still run from stage 15 if stage 22 is absent.
    if stage22_path.exists():
        df = pd.read_excel(stage22_path, sheet_name="Primary_Haplotypes")
        df = df[df["FDR"].fillna(1) < FDR_THRESHOLD].copy()
        df["Exposure_ID"] = df["Region"].astype(str) + "::" + df["Haplotype_ID"].astype(str)
        df["Interpretation_Source"] = "Stage22_Primary_Haplotypes"
        return df.sort_values(["FDR", "Region", "Haplotype_ID"]).reset_index(drop=True)
    return load_significant_haplotypes_from_stage15(stage15_path)


def simplify_comparison(value):
    text = s(value).lower()
    if "breast" in text:
        return "Breast"
    if "endo" in text:
        return "Endometrial"
    return "All_Primary"


def load_significant_mutations(breast_path, endo_path, exposure_catalogue_path):
    if not (breast_path.exists() and endo_path.exists() and exposure_catalogue_path.exists()):
        return empty_significant_mutation_table()

    frames = []
    for path in [breast_path, endo_path]:
        df = pd.read_csv(path, sep="	")
        sig_cols = [c for c in df.columns if c.startswith("FDR_Sig_")]
        if not sig_cols:
            continue
        df = df[df[sig_cols].fillna(False).any(axis=1)].copy()
        if df.empty:
            continue
        df["Comparison"] = df["Cohort"].map(simplify_comparison)
        model_labels = [c.replace("FDR_Sig_", "") for c in sig_cols]
        df["Significant_Models"] = df.apply(
            lambda row: "; ".join(
                label for label, sig_col in zip(model_labels, sig_cols) if bool(row.get(sig_col, False))
            ),
            axis=1,
        )

        def _min_sig_fdr(row):
            vals = []
            for label, sig_col in zip(model_labels, sig_cols):
                if not bool(row.get(sig_col, False)):
                    continue
                val = pd.to_numeric(row.get(f"FDR_{label}"), errors="coerce")
                if pd.notna(val):
                    vals.append(float(val))
            return min(vals) if vals else np.nan

        df["Min_FDR"] = df.apply(_min_sig_fdr, axis=1)
        frames.append(df)

    if not frames:
        return empty_significant_mutation_table()

    assoc = pd.concat(frames, ignore_index=True)
    expo = pd.read_csv(exposure_catalogue_path, sep="	")
    expo_cols = [c for c in [
        "Exposure_ID", "Exposure_Label", "Exposure_Type", "Gene", "Consequence",
        "IMPACT", "Component_Genes", "Component_Variants", "Component_Variant_Count",
        "Carrier_Samples",
    ] if c in expo.columns]
    assoc = assoc.merge(expo[expo_cols].drop_duplicates("Exposure_ID"), on="Exposure_ID", how="left", suffixes=("", "_catalogue"))
    for col in ["Exposure_Label", "Exposure_Type", "Gene", "Consequence", "IMPACT", "Component_Genes", "Component_Variants", "Component_Variant_Count", "Carrier_Samples"]:
        cat_col = f"{col}_catalogue"
        if cat_col in assoc.columns:
            assoc[col] = assoc[col].fillna(assoc[cat_col])
        if col not in assoc.columns:
            assoc[col] = np.nan

    assoc["Base_Exposure_ID"] = assoc["Exposure_ID"].astype(str)
    assoc["Exposure_ID"] = assoc["Comparison"].astype(str) + "::" + assoc["Base_Exposure_ID"].astype(str)
    grouped = (
        assoc.groupby(["Exposure_ID", "Comparison", "Base_Exposure_ID"], dropna=False)
        .agg(
            Exposure_Label=("Exposure_Label", first_valid),
            Exposure_Type=("Exposure_Type", first_valid),
            Gene=("Gene", first_valid),
            Consequence=("Consequence", first_valid),
            IMPACT=("IMPACT", first_valid),
            Component_Genes=("Component_Genes", first_valid),
            Component_Variants=("Component_Variants", first_valid),
            Component_Variant_Count=("Component_Variant_Count", "max"),
            Significant_Clinical_Vars=("Clinical_Var", lambda x: "; ".join(sorted({s(v) for v in x if s(v)}))),
            Significant_Models=("Significant_Models", lambda x: "; ".join(sorted({part.strip() for v in x for part in s(v).split(";") if part.strip()}))),
            Min_FDR=("Min_FDR", "min"),
            N_Carriers_Cohort=("N_Carriers_Cohort", "max"),
            Carrier_Samples=("Carrier_Samples", first_valid),
        )
        .reset_index()
        .sort_values(["Comparison", "Min_FDR", "Base_Exposure_ID"], na_position="last")
        .reset_index(drop=True)
    )
    return grouped


def load_phased_matrix(phased_path, stage15_path):
    # The phased matrix lets us calculate SNP dosage and haplotype exposure for
    # the same samples represented in the RNA workbook.
    phased = pd.read_csv(phased_path, sep="	", dtype=str)
    phased["POS"] = pd.to_numeric(phased["POS"], errors="coerce")
    phased["REF"] = phased["REF"].astype(str).str.upper()
    phased["ALT"] = phased["ALT"].astype(str).str.upper()
    phased["SNP_ID"] = phased.apply(lambda row: genomic_snp_id(row.get("CHROM"), row.get("POS"), row.get("REF"), row.get("ALT")), axis=1)
    snps_used = pd.read_excel(stage15_path, sheet_name="SNPs_Used").rename(columns={"POS_int": "POS", "rsID_clean": "Backbone_rsID"})
    snps_used["POS"] = pd.to_numeric(snps_used["POS"], errors="coerce")
    snps_used["REF"] = snps_used["REF"].astype(str).str.upper()
    snps_used["ALT"] = snps_used["ALT"].astype(str).str.upper()
    snps_used["SNP_ID"] = snps_used.apply(lambda row: genomic_snp_id("chr17", row.get("POS"), row.get("REF"), row.get("ALT")), axis=1)
    phased = phased.merge(snps_used[["POS", "REF", "ALT", "Backbone_rsID", "Gene"]], on=["POS", "REF", "ALT"], how="left")
    phased["SNP_rsID"] = phased["Backbone_rsID"].astype(str).str.strip()
    phased.loc[phased["SNP_rsID"].isin({"", "nan", "None"}), "SNP_rsID"] = pd.NA
    phased["SNP_Match_Key"] = phased["SNP_rsID"].fillna(phased["SNP_ID"])
    sample_cols = [c for c in phased.columns if c not in {"CHROM", "POS", "ID", "REF", "ALT", "SNP_ID", "SNP_rsID", "SNP_Match_Key", "Backbone_rsID", "Gene"}]
    return phased, sample_cols, snps_used


def build_significant_snp_dosage(phased, sample_cols, sig_snp_df):
    targets = set(sig_snp_df["SNP_Match_Key"].dropna().astype(str))
    available = phased[phased["SNP_Match_Key"].isin(targets)].copy()
    identity = build_analysis_sample_map(sample_cols, PATHS["phased"], sx)
    rep = identity.sort_values(["analysis_sample_id", "raw_sample_name"]).drop_duplicates("analysis_sample_id")
    rows = []
    for _, row in available.iterrows():
        display_id = s(row.get("SNP_rsID")) or s(row.get("SNP_Match_Key")) or s(row.get("SNP_ID"))
        for _, id_row in rep.iterrows():
            col = id_row["raw_sample_name"]
            rows.append({
                "analysis_sample_id": id_row["analysis_sample_id"],
                "snp_code": id_row["snp_code"],
                "SNP_ID": display_id,
                "SNP_Technical_ID": row.get("SNP_ID", np.nan),
                "SNP_rsID": row.get("SNP_rsID", np.nan),
                "SNP_Match_Key": row.get("SNP_Match_Key", np.nan),
                "Gene": row.get("Gene", np.nan),
                "Dosage": gt_dose(row[col]),
                "GT": row[col],
            })
    dosage = pd.DataFrame(rows)
    tested = set(dosage["SNP_Match_Key"].dropna().astype(str)) if not dosage.empty else set()
    summary = sig_snp_df.copy()
    summary["RNA_Testable"] = summary["SNP_Match_Key"].astype(str).isin(tested)
    summary["RNA_Test_Note"] = np.where(summary["RNA_Testable"], "Phased/common backbone dosage available", "Significant in stage 12 but absent from phased/common backbone")
    return dosage, summary


def region_membership(stage15_path, snps_used):
    members = {"Full_Region": list(snps_used["Backbone_rsID"].dropna().astype(str))}
    blocks = pd.read_excel(stage15_path, sheet_name="LD_Blocks")
    for (thr, block), sub in blocks.groupby(["Threshold_Label", "Block"], dropna=False):
        members[f"LD_Block_{thr}_{block}"] = sub.sort_values("SNP_Order")["rsID"].dropna().astype(str).tolist()
    return members


def build_significant_haplotype_dosage(phased, sample_cols, stage15_path, hap_df, snps_used):
    members = region_membership(stage15_path, snps_used)
    row_lookup = {str(r["Backbone_rsID"]): r for _, r in phased.dropna(subset=["Backbone_rsID"]).iterrows()}
    identity = build_analysis_sample_map(sample_cols, PATHS["phased"], sx)
    rep = identity.sort_values(["analysis_sample_id", "raw_sample_name"]).drop_duplicates("analysis_sample_id")
    rows = []
    checks = []
    for _, hap in hap_df.iterrows():
        region = str(hap["Region"])
        target = str(hap["Haplotype_String"])
        exposure = str(hap["Exposure_ID"])
        snp_ids = members.get(region, [])
        callable_n = 0
        alt_count = 0
        for _, id_row in rep.iterrows():
            col = id_row["raw_sample_name"]
            code = id_row["snp_code"]
            analysis_sample_id = id_row["analysis_sample_id"]
            hap1, hap2 = [], []
            for snp_id in snp_ids:
                src = row_lookup.get(snp_id)
                if src is None:
                    hap1.append("N")
                    hap2.append("N")
                    continue
                a, b = gt_haps(src[col])
                hap1.append(a)
                hap2.append(b)
            h1, h2 = "".join(hap1), "".join(hap2)
            if "N" in h1 or "N" in h2:
                dose = np.nan
                carrier = np.nan
                callable_flag = 0
            else:
                dose = int(h1 == target) + int(h2 == target)
                carrier = int(dose >= 1)
                callable_flag = 1
                callable_n += 1
                alt_count += dose
            rows.append({
                "analysis_sample_id": analysis_sample_id,
                "snp_code": code,
                "Exposure_ID": exposure,
                "Region": region,
                "Haplotype_ID": hap["Haplotype_ID"],
                "Haplotype_String": target,
                "Dosage": dose,
                "Carrier": carrier,
                "Callable": callable_flag,
            })
        checks.append({
            "Exposure_ID": exposure,
            "Region": region,
            "Callable_Samples": callable_n,
            "RNA_Matched_Global_Freq": (alt_count / (2 * callable_n)) if callable_n else np.nan,
            "Stage22_Global_Freq": hap.get("Global_Freq", np.nan),
            "Abs_Delta_vs_Stage22_Global": abs(((alt_count / (2 * callable_n)) if callable_n else np.nan) - float(hap.get("Global_Freq", np.nan))) if callable_n and pd.notna(hap.get("Global_Freq", np.nan)) else np.nan,
            "Frequency_Check_Note": "Recomputed from phased dosage across all callable study samples",
        })
    return pd.DataFrame(rows), pd.DataFrame(checks)


def context_for_sources(source_groups):
    groups = {x.strip() for x in s(source_groups).split(";") if x.strip()}
    out = []
    if "Global" in groups:
        out.append("All_Primary")
    if "Breast" in groups:
        out.append("Breast")
    if "Endometrium" in groups:
        out.append("Endometrial")
    return out or ["All_Primary"]


def build_na_variable_catalog(excel_vars, bam_vars):
    return pd.concat([excel_vars, bam_vars], ignore_index=True).sort_values(["RNA_Source", "Integrated_Var"]).reset_index(drop=True)


def primary_subset(na_df):
    gate = na_df.get("RNA_QC_Exploratory_Ready", na_df.get("RNA_QC_Analysis_Ready", False)).fillna(False).astype(bool)
    role = na_df.get("analysis_role", "").astype(str)
    return na_df[gate & (role == "primary_tumour")].copy()


def subset_for_context(df, context):
    return df.copy() if context == "All_Primary" else df[df["cohort"] == context].copy()


def build_significant_mutation_dosage(sig_mut_df, sample_codes):
    if sig_mut_df.empty:
        return pd.DataFrame(), empty_significant_mutation_table()

    sample_codes = sorted({s(code) for code in sample_codes if s(code)})
    rows = []
    summary_rows = []
    for _, mut in sig_mut_df.iterrows():
        carrier_ids = [part.strip() for part in s(mut.get("Carrier_Samples")).split(";") if part.strip()]
        carrier_codes = {analysis_id_to_snp_code(part) for part in carrier_ids}
        carrier_codes = {code for code in carrier_codes if code}
        matched_carriers = sorted(code for code in sample_codes if code in carrier_codes)
        matched_noncarriers = [code for code in sample_codes if code not in carrier_codes]
        for code in sample_codes:
            rows.append({
                "analysis_sample_id": code,
                "snp_code": code,
                "Exposure_ID": mut["Exposure_ID"],
                "Base_Exposure_ID": mut.get("Base_Exposure_ID", np.nan),
                "Comparison": mut.get("Comparison", np.nan),
                "Exposure_Type": mut.get("Exposure_Type", np.nan),
                "Gene": mut.get("Gene", np.nan),
                "Component_Genes": mut.get("Component_Genes", np.nan),
                "Component_Variants": mut.get("Component_Variants", np.nan),
                "Component_Variant_Count": mut.get("Component_Variant_Count", np.nan),
                "Dosage": 1.0 if code in carrier_codes else 0.0,
            })
        note = (
            f"{len(matched_carriers)} RNA-matched carrier snp_code(s); "
            f"{len(matched_noncarriers)} RNA-matched non-carrier snp_code(s)"
        )
        summary_row = mut.to_dict()
        summary_row["RNA_Testable"] = bool(matched_carriers and matched_noncarriers)
        summary_row["RNA_Test_Note"] = note
        summary_rows.append(summary_row)
    return pd.DataFrame(rows), pd.DataFrame(summary_rows)


def regression_rows(na_df, dosage_df, exposures_df, exposure_type, var_catalog):
    value_cols = var_catalog["Integrated_Var"].tolist()
    var_level_map = var_catalog.set_index("Integrated_Var")["RNA_Level"].to_dict()
    rows = []
    exposure_meta = exposures_df.copy()
    if exposure_type == "SNP":
        id_col = "SNP_ID"
        if exposure_meta.empty or id_col not in exposure_meta.columns or id_col not in dosage_df.columns:
            return pd.DataFrame()
        cohort_map = exposure_meta.set_index(id_col)["Source_Groups"].to_dict()
        meta_map = exposure_meta.set_index(id_col).to_dict(orient="index")
    else:
        id_col = "Exposure_ID"
        if exposure_meta.empty or id_col not in exposure_meta.columns or id_col not in dosage_df.columns:
            return pd.DataFrame()
        cohort_map = {r["Exposure_ID"]: r["Comparison"] for _, r in exposure_meta.iterrows()}
        meta_map = exposure_meta.set_index(id_col).to_dict(orient="index")

    for exposure_id, g in dosage_df.groupby(id_col, dropna=False):
        if exposure_type == "SNP":
            contexts = context_for_sources(cohort_map.get(exposure_id, ""))
        else:
            comp = simplify_comparison(cohort_map.get(exposure_id, "All_Primary"))
            contexts = [comp]
        for context in contexts:
            base = subset_for_context(na_df, context)
            merge_cols = [c for c in ["analysis_sample_id", "snp_code", "Dosage"] if c in g.columns]
            dedup_keys = [c for c in ["analysis_sample_id"] if c in merge_cols]
            merged = base.merge(g[merge_cols].drop_duplicates(dedup_keys if dedup_keys else None), on="snp_code", how="inner")
            if merged.empty:
                continue
            for var in value_cols:
                if var not in merged.columns:
                    continue
                t = merged[[var, "Dosage"]].copy()
                t[var] = pd.to_numeric(t[var], errors="coerce")
                t["Dosage"] = pd.to_numeric(t["Dosage"], errors="coerce")
                t = t.dropna()
                if len(t) < MIN_N or t["Dosage"].nunique() < 2 or int((t["Dosage"] > 0).sum()) < MIN_CARRIERS:
                    continue
                if t[var].nunique() < 2:
                    continue
                slope, intercept, r_value, p_value, se = stats.linregress(t["Dosage"], t[var])
                dose0 = t.loc[t["Dosage"] == 0, var].dropna()
                dose1 = t.loc[t["Dosage"] == 1, var].dropna()
                dose2 = t.loc[t["Dosage"] == 2, var].dropna()
                carrier = t.loc[t["Dosage"] > 0, var].dropna()
                noncarrier = t.loc[t["Dosage"] == 0, var].dropna()
                p_carrier = np.nan
                if len(carrier) >= MIN_CARRIERS and len(noncarrier) >= MIN_CARRIERS:
                    p_carrier = stats.mannwhitneyu(carrier, noncarrier, alternative="two-sided")[1]
                meta = meta_map.get(exposure_id, {})
                row = {
                    "Exposure_Type": exposure_type,
                    id_col: exposure_id,
                    "Context": context,
                    "RNA_Var": var,
                    "RNA_Source": var.split("__", 1)[0],
                    "RNA_Level": var_level_map.get(var, "other_rna_quant"),
                    "N": len(t),
                    "Carrier_Count": int((t["Dosage"] > 0).sum()),
                    "N_Dose0": int((t["Dosage"] == 0).sum()),
                    "N_Dose1": int((t["Dosage"] == 1).sum()),
                    "N_Dose2": int((t["Dosage"] == 2).sum()),
                    "Slope": slope,
                    "Intercept": intercept,
                    "R": r_value,
                    "SE": se,
                    "P_Value": p_value,
                    "P_Carrier": p_carrier,
                    "Median_Dose0": dose0.median() if len(dose0) else np.nan,
                    "Median_Dose1": dose1.median() if len(dose1) else np.nan,
                    "Median_Dose2": dose2.median() if len(dose2) else np.nan,
                    "Median_Carrier_Minus_Noncarrier": carrier.median() - noncarrier.median() if len(carrier) and len(noncarrier) else np.nan,
                    "Effect_Direction": "Higher_with_dosage" if slope > 0 else "Lower_with_dosage",
                }
                for key in ["Source_Groups", "LD_Block_ID", "Min_FDR", "Min_P", "Region", "Haplotype_ID", "Comparison", "Global_Freq", "Overlapping_Significant_SNPs", "Consistency_Summary", "Base_Exposure_ID", "Exposure_Label", "Gene", "Consequence", "IMPACT", "Component_Genes", "Component_Variants", "Component_Variant_Count", "Significant_Clinical_Vars", "Significant_Models"]:
                    if key in meta:
                        row[key] = meta[key]
                rows.append(row)
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return add_fdr(result, "P_Value", ["Exposure_Type", id_col, "Context", "RNA_Source"])


def _clr_profile_frame(df, profile_cols):
    work = df[["Dosage"] + list(profile_cols)].copy()
    for col in work.columns:
        work[col] = pd.to_numeric(work[col], errors="coerce")
    work = work.dropna()
    if work.empty:
        return np.empty((0, len(profile_cols))), work
    valid = (work[list(profile_cols)] >= 0).all(axis=1) & (work[list(profile_cols)].sum(axis=1) > 0)
    work = work[valid].copy()
    if work.empty:
        return np.empty((0, len(profile_cols))), work
    props = work[list(profile_cols)].div(work[list(profile_cols)].sum(axis=1), axis=0)
    logs = np.log(props.to_numpy(dtype=float) + 1e-6)
    clr = logs - logs.mean(axis=1, keepdims=True)
    return clr, work


def _profile_pseudo_f(clr, groups):
    groups = pd.Series(groups).reset_index(drop=True)
    if clr.shape[0] < MIN_N:
        return np.nan, np.nan, np.nan, np.nan
    levels = [level for level in sorted(groups.dropna().unique()) if (groups == level).sum() > 0]
    if len(levels) < 2 or clr.shape[0] <= len(levels):
        return np.nan, np.nan, np.nan, np.nan
    centroid = clr.mean(axis=0)
    total_ss = float(((clr - centroid) ** 2).sum())
    within_ss = 0.0
    for level in levels:
        x = clr[(groups == level).to_numpy()]
        within_ss += float(((x - x.mean(axis=0)) ** 2).sum())
    between_ss = max(total_ss - within_ss, 0.0)
    df_between = len(levels) - 1
    df_within = clr.shape[0] - len(levels)
    if df_between <= 0 or df_within <= 0 or within_ss <= 0:
        return np.nan, between_ss, within_ss, total_ss
    pseudo_f = (between_ss / df_between) / (within_ss / df_within)
    return float(pseudo_f), between_ss, within_ss, total_ss


def _permutation_profile_p(clr, groups, n_perm=499, rng=None):
    rng = rng or np.random.default_rng(42)
    observed, between_ss, within_ss, total_ss = _profile_pseudo_f(clr, groups)
    if not np.isfinite(observed):
        return observed, np.nan, between_ss, within_ss, total_ss
    group_values = pd.Series(groups).to_numpy()
    hits = 0
    for _ in range(n_perm):
        permuted = rng.permutation(group_values)
        perm_f, _, _, _ = _profile_pseudo_f(clr, permuted)
        if np.isfinite(perm_f) and perm_f >= observed:
            hits += 1
    p_value = (hits + 1) / (n_perm + 1)
    return observed, float(p_value), between_ss, within_ss, total_ss


def isoform_profile_rows(na_df, dosage_df, exposures_df, exposure_type, n_perm=499):
    profile_cols = [col for col in ISOFORM_PROFILE_VARS if col in na_df.columns]
    if len(profile_cols) < 3:
        return pd.DataFrame()
    rows = []
    exposure_meta = exposures_df.copy()
    if exposure_type == "SNP":
        id_col = "SNP_ID"
        if exposure_meta.empty or id_col not in exposure_meta.columns or id_col not in dosage_df.columns:
            return pd.DataFrame()
        cohort_map = exposure_meta.set_index(id_col)["Source_Groups"].to_dict()
        meta_map = exposure_meta.set_index(id_col).to_dict(orient="index")
    else:
        id_col = "Exposure_ID"
        if exposure_meta.empty or id_col not in exposure_meta.columns or id_col not in dosage_df.columns:
            return pd.DataFrame()
        cohort_map = {r["Exposure_ID"]: r["Comparison"] for _, r in exposure_meta.iterrows()}
        meta_map = exposure_meta.set_index(id_col).to_dict(orient="index")

    rng = np.random.default_rng(42)
    for exposure_id, g in dosage_df.groupby(id_col, dropna=False):
        contexts = context_for_sources(cohort_map.get(exposure_id, "")) if exposure_type == "SNP" else [simplify_comparison(cohort_map.get(exposure_id, "All_Primary"))]
        for context in contexts:
            base = subset_for_context(na_df, context)
            merge_cols = [c for c in ["analysis_sample_id", "snp_code", "Dosage"] if c in g.columns]
            dedup_keys = [c for c in ["analysis_sample_id"] if c in merge_cols]
            merged = base.merge(g[merge_cols].drop_duplicates(dedup_keys if dedup_keys else None), on="snp_code", how="inner")
            if merged.empty:
                continue
            clr, work = _clr_profile_frame(merged, profile_cols)
            if work.empty:
                continue
            carrier_count = int((work["Dosage"] > 0).sum())
            dose_group_sizes = work.groupby("Dosage", dropna=True).size()
            if len(work) < MIN_N or work["Dosage"].nunique() < 2 or carrier_count < MIN_CARRIERS or int((dose_group_sizes >= MIN_CARRIERS).sum()) < 2:
                continue
            pseudo_f, p_value, between_ss, within_ss, total_ss = _permutation_profile_p(clr, work["Dosage"], n_perm=n_perm, rng=rng)
            if not np.isfinite(p_value):
                continue
            slopes = {}
            for col in profile_cols:
                iso = col.replace("excel_iso_ratio__", "").replace("_to_isoform_total", "")
                y = pd.to_numeric(work[col], errors="coerce")
                x = pd.to_numeric(work["Dosage"], errors="coerce")
                t = pd.concat([x, y], axis=1).dropna()
                if len(t) >= MIN_N and t.iloc[:, 0].nunique() >= 2 and t.iloc[:, 1].nunique() >= 2:
                    slopes[iso] = float(stats.linregress(t.iloc[:, 0], t.iloc[:, 1]).slope)
            ordered = sorted(slopes.items(), key=lambda item: abs(item[1]), reverse=True)
            meta = meta_map.get(exposure_id, {})
            row = {
                "Exposure_Type": exposure_type,
                id_col: exposure_id,
                "Context": context,
                "Profile_Test": "PERMANOVA_Aitchison_CLR_by_dosage",
                "N": len(work),
                "Carrier_Count": carrier_count,
                "N_Dose0": int((work["Dosage"] == 0).sum()),
                "N_Dose1": int((work["Dosage"] == 1).sum()),
                "N_Dose2": int((work["Dosage"] == 2).sum()),
                "Pseudo_F": pseudo_f,
                "P_Value": p_value,
                "Between_SS": between_ss,
                "Within_SS": within_ss,
                "Total_SS": total_ss,
                "Top_Isoform_Shift": ordered[0][0] if ordered else np.nan,
                "Top_Isoform_Slope": ordered[0][1] if ordered else np.nan,
                "Positive_Isoform_Shifts": "; ".join(f"{iso}:{slope:.4g}" for iso, slope in ordered if slope > 0),
                "Negative_Isoform_Shifts": "; ".join(f"{iso}:{slope:.4g}" for iso, slope in ordered if slope < 0),
                "Profile_Variables": "; ".join(profile_cols),
            }
            for key in ["Source_Groups", "LD_Block_ID", "Min_FDR", "Min_P", "Region", "Haplotype_ID", "Comparison", "Global_Freq", "Overlapping_Significant_SNPs", "Consistency_Summary", "Base_Exposure_ID", "Exposure_Label", "Gene", "Consequence", "IMPACT", "Component_Genes", "Component_Variants", "Component_Variant_Count", "Significant_Clinical_Vars", "Significant_Models"]:
                if key in meta:
                    row[key] = meta[key]
            rows.append(row)
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return add_fdr(result, "P_Value", ["Exposure_Type", "Context"])


def plot_heatmap(df, row_col, out_path, title):
    if df.empty:
        return
    work = df.copy()
    work["_score"] = -np.log10(pd.to_numeric(work["P_Value"], errors="coerce").clip(lower=1e-300))
    work["_signed"] = work["_score"] * np.sign(pd.to_numeric(work["Slope"], errors="coerce").fillna(0))
    top = work.groupby(row_col)["_score"].max().sort_values(ascending=False).head(20).index
    work = work[work[row_col].isin(top)].copy()
    work["RNA_Label"] = work["RNA_Var"].astype(str).str.replace(r"^[^_]+__", "", regex=True)
    piv = work.pivot_table(index=row_col, columns="RNA_Label", values="_signed", aggfunc="max")
    if piv.empty:
        return
    fig_w = max(12, 0.6 * len(piv.columns) + 4)
    fig_h = max(4, 0.4 * len(piv.index) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    if sns is not None:
        sns.heatmap(
            piv,
            cmap="coolwarm",
            center=0,
            linewidths=0.4,
            linecolor="white",
            ax=ax,
            cbar_kws={"label": "signed -log10(P value)"},
        )
    else:
        im = ax.imshow(piv.values, aspect="auto", cmap="coolwarm")
        ax.set_xticks(range(len(piv.columns)))
        ax.set_xticklabels(piv.columns, rotation=40, ha="right", fontsize=8)
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels(piv.index, fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="signed -log10(P value)")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def summarise_hits(snp_assoc, mutation_assoc, hap_assoc):
    rows = []
    for _, df, id_col in [("SNP", snp_assoc, "SNP_ID"), ("Mutation", mutation_assoc, "Exposure_ID"), ("Haplotype", hap_assoc, "Exposure_ID")]:
        if df.empty:
            continue
        sig = df[df["FDR_Sig"].fillna(False)].copy()
        if sig.empty:
            top = df.sort_values(["P_Value", "N"], ascending=[True, False]).groupby(id_col, dropna=False).head(3)
            top = top.assign(Summary_Status="Top_nominal_no_FDR_hit")
            rows.append(top)
        else:
            sig = sig.assign(Summary_Status="FDR_significant")
            rows.append(sig)
    if not rows:
        return pd.DataFrame({"Note": ["No RNA associations met the reporting thresholds."]})
    out = pd.concat(rows, ignore_index=True)
    keep = [c for c in [
        "Summary_Status", "Exposure_Type", "SNP_ID", "Exposure_ID", "Context", "RNA_Var", "RNA_Source", "RNA_Level",
        "P_Value", "FDR", "Slope", "R", "N", "Carrier_Count", "Effect_Direction",
        "Source_Groups", "Base_Exposure_ID", "Exposure_Label", "Gene", "Significant_Clinical_Vars", "Significant_Models",
        "Component_Genes", "Component_Variants", "Component_Variant_Count",
        "Region", "Haplotype_ID", "Comparison", "Overlapping_Significant_SNPs",
    ] if c in out.columns]
    return out[keep].sort_values(["Summary_Status", "P_Value", "N"], ascending=[True, True, False]).reset_index(drop=True)


def build_gsdmb_isoform_answer(snp_assoc, mutation_assoc, hap_assoc):
    rows = []
    requested_levels = {
        "overall_gene_expression",
        "isoform_expression",
        "isoform_composition",
        "isoform_log_ratio",
        "erbb2_normalised_expression",
    }
    for exposure_type, df, id_col in [("SNP", snp_assoc, "SNP_ID"), ("Mutation", mutation_assoc, "Exposure_ID"), ("Haplotype", hap_assoc, "Exposure_ID")]:
        if df.empty:
            continue
        work = df[df["RNA_Level"].astype(str).isin(requested_levels)].copy()
        if work.empty:
            continue
        work["Exposure_Label_Final"] = work[id_col].astype(str)
        work["RNA_Label"] = work["RNA_Var"].astype(str).str.replace(r"^[^_]+__", "", regex=True).str.replace("_", " ", regex=False)
        work["Summary_Status"] = np.where(
            work["FDR_Sig"].fillna(False),
            "FDR_significant",
            np.where(pd.to_numeric(work["P_Value"], errors="coerce") < 0.05, "Nominal_only", "Top_ranked_context"),
        )
        work["Interpretation"] = work.apply(
            lambda row: f"{row['RNA_Label']} is {'higher' if float(row['Slope']) > 0 else 'lower'} with higher {exposure_type.lower()} dosage in {row['Context']}",
            axis=1,
        )
        keep = [c for c in [
            "Summary_Status", "Exposure_Type", id_col, "Exposure_Label_Final", "Context", "RNA_Label",
            "P_Value", "FDR", "Slope", "R", "N", "Carrier_Count", "N_Dose0", "N_Dose1", "N_Dose2",
            "Effect_Direction", "Interpretation", "Region", "Haplotype_ID", "Comparison",
            "Base_Exposure_ID", "Exposure_Label", "Gene", "Significant_Clinical_Vars", "Significant_Models",
            "Component_Genes", "Component_Variants", "Component_Variant_Count",
        ] if c in work.columns]
        rows.append(work[keep])
    if not rows:
        return pd.DataFrame({"Note": ["No SNP, mutation, or haplotype associations reached the GSDMB isoform-expression reporting threshold."]})
    out = pd.concat(rows, ignore_index=True)
    out = out.sort_values(["Summary_Status", "P_Value", "N"], ascending=[True, True, False]).reset_index(drop=True)
    return out


def build_gsdmb_isoform_summary(gsdmb_isoform_answer):
    """Collapse repeated LD-equivalent exposure rows into a supervisor-facing summary."""
    if gsdmb_isoform_answer.empty or "Note" in gsdmb_isoform_answer.columns:
        return gsdmb_isoform_answer.copy()

    work = gsdmb_isoform_answer.copy()
    work["Exposure_Label_Final"] = work["Exposure_Label_Final"].astype(str)
    group_cols = [
        "Summary_Status",
        "Exposure_Type",
        "Context",
        "RNA_Label",
        "Effect_Direction",
        "P_Value",
        "FDR",
        "Slope",
        "R",
        "N",
        "Carrier_Count",
        "N_Dose0",
        "N_Dose1",
        "N_Dose2",
        "Interpretation",
    ]
    rows = []
    for keys, sub in work.groupby(group_cols, dropna=False, sort=False):
        exposure_list = list(dict.fromkeys(sub["Exposure_Label_Final"].astype(str).tolist()))
        representative = exposure_list[0] if exposure_list else ""
        note = ""
        if str(keys[1]) == "SNP" and len(exposure_list) > 1:
            note = f"{len(exposure_list)} linked SNPs share this isoform-association pattern."
        elif len(exposure_list) > 1:
            note = f"{len(exposure_list)} exposures share this isoform-association pattern."
        row = {
            "Summary Status": keys[0],
            "Exposure Type": keys[1],
            "Context": keys[2],
            "RNA Label": keys[3],
            "Effect Direction": keys[4],
            "P-Value": keys[5],
            "FDR": keys[6],
            "Slope": keys[7],
            "R": keys[8],
            "N": keys[9],
            "Carrier Count": keys[10],
            "N Dose 0": keys[11],
            "N Dose 1": keys[12],
            "N Dose 2": keys[13],
            "Representative Exposure": representative,
            "Exposure Count": len(exposure_list),
            "Exposure List": "; ".join(exposure_list),
            "Interpretation": keys[14],
            "Summary Note": note,
        }
        if "Component_Genes" in sub.columns:
            row["Component Genes"] = "; ".join(sorted({part.strip() for value in sub["Component_Genes"].dropna().astype(str) for part in value.split(";") if part.strip()}))
        if "Component_Variants" in sub.columns:
            row["Component Variants"] = "; ".join(sorted({part.strip() for value in sub["Component_Variants"].dropna().astype(str) for part in value.split(";") if part.strip()}))
        if "Component_Variant_Count" in sub.columns:
            row["Component Variant Count"] = pd.to_numeric(sub["Component_Variant_Count"], errors="coerce").max()
        rows.append(row)
    order = {"FDR_significant": 0, "Nominal_only": 1, "Top_ranked_context": 2}
    out = pd.DataFrame(rows)
    out["__order"] = out["Summary Status"].map(order).fillna(9)
    out = out.sort_values(["__order", "P-Value", "Exposure Type", "Context", "RNA Label"], kind="stable").drop(columns="__order")
    return out.reset_index(drop=True)


def source_inventory(manifest, na_catalog):
    return pd.DataFrame([
        {
            "Source": "Excel_matched_rows",
            "Samples": int(manifest["has_excel_row"].fillna(False).sum()),
            "Variables": int((na_catalog["RNA_Source"].astype(str).str.startswith("excel")).sum()),
            "Note": "Matched and collapsed from stage-20 expression workbook; these values are treated as source-normalised assay outputs and are preserved in their original workbook scale.",
        },
        {
            "Source": "BAM_ratio_rows",
            "Samples": int(manifest["has_bam_ratio"].fillna(False).sum()),
            "Variables": int((na_catalog["RNA_Source"] == "bam_ratio").sum()),
            "Note": "Ratio-normalised gene-level panel expression from stage 20b using both strict and exploratory BAMs when available.",
        },
    ])


def gate_usage_summary(manifest):
    if manifest.empty:
        return pd.DataFrame({"Note": ["No RNA manifest rows were available."]})
    rows = []
    for label, sub in [("All_RNA", manifest), ("Primary_Tumour", manifest[manifest.get("analysis_role", "").astype(str) == "primary_tumour"].copy())]:
        if sub.empty:
            continue
        rows.append({
            "Scope": label,
            "Unique_Samples": int(sub["snp_code"].dropna().astype(str).nunique()),
            "Strict_Ready": int(sub.get("RNA_QC_Analysis_Ready", False).fillna(False).sum()),
            "Exploratory_Ready": int(sub.get("RNA_QC_Exploratory_Ready", sub.get("RNA_QC_Analysis_Ready", False)).fillna(False).sum()),
            "Exploratory_Only": int((sub.get("RNA_QC_Exploratory_Ready", False).fillna(False) & ~sub.get("RNA_QC_Analysis_Ready", False).fillna(False)).sum()),
            "With_BAM_Ratio": int(sub.get("has_bam_ratio", False).fillna(False).sum()),
            "Excel_Only": int((sub.get("has_excel_row", False).fillna(False) & ~sub.get("has_bam_ratio", False).fillna(False)).sum()),
            "Notes": "Stage-23 associations use the exploratory-ready primary subset; BAM ratio availability can still lag behind manifest eligibility if stage 20b is stale.",
        })
    return pd.DataFrame(rows)


def build_excel_bam_concordance(na_df):
    if na_df.empty:
        return pd.DataFrame({"Note": ["No RNA rows were available for Excel-vs-BAM concordance checks."]})
    work = primary_subset(na_df)
    if work.empty:
        return pd.DataFrame({"Note": ["No exploratory-ready primary RNA rows were available for concordance checks."]})
    pairs = []
    cols = set(work.columns)
    if "excel_total__GSDMB" in cols and "bam_ratio__GSDMB" in cols:
        pairs.append(("GSDMB", "excel_total__GSDMB", "bam_ratio__GSDMB", "Excel total GSDMB vs BAM ratio-normalised GSDMB"))
    if "excel_panel__GSDMB_EXPRESSION_IN_PANEL" in cols and "bam_ratio__GSDMB" in cols:
        pairs.append(("GSDMB_PANEL", "excel_panel__GSDMB_EXPRESSION_IN_PANEL", "bam_ratio__GSDMB", "Excel panel GSDMB vs BAM ratio-normalised GSDMB"))
    for col in sorted(cols):
        if not str(col).startswith("excel_panel__"):
            continue
        gene = str(col).split("__", 1)[1]
        if gene == "GSDMB_EXPRESSION_IN_PANEL":
            continue
        bam_col = f"bam_ratio__{gene}"
        if bam_col in cols:
            pairs.append((gene, col, bam_col, "Excel panel gene vs BAM ratio-normalised gene"))
    rows = []
    for gene, excel_col, bam_col, comparison_note in pairs:
        for context, sub in [("All_Primary", work), ("Breast", work[work["cohort"] == "Breast"].copy()), ("Endometrial", work[work["cohort"] == "Endometrial"].copy())]:
            if sub.empty:
                continue
            t = sub[[excel_col, bam_col]].copy()
            t[excel_col] = pd.to_numeric(t[excel_col], errors="coerce")
            t[bam_col] = pd.to_numeric(t[bam_col], errors="coerce")
            t = t.dropna()
            if len(t) < 5:
                continue
            if t[excel_col].nunique() < 2 or t[bam_col].nunique() < 2:
                continue
            spearman_rho, spearman_p = stats.spearmanr(t[excel_col], t[bam_col], nan_policy="omit")
            try:
                pearson_r = np.corrcoef(t[excel_col], t[bam_col])[0, 1]
            except Exception:
                pearson_r = np.nan
            rows.append({
                "Context": context,
                "Gene": gene,
                "Excel_Var": excel_col,
                "BAM_Var": bam_col,
                "N": len(t),
                "Spearman_Rho": spearman_rho,
                "Spearman_P": spearman_p,
                "Pearson_R": pearson_r,
                "Excel_Median": t[excel_col].median(),
                "BAM_Median": t[bam_col].median(),
                "Comparison_Note": comparison_note,
                "Scale_Note": "Excel values are treated as source-normalised workbook measurements, whereas BAM values come from stage-20b ratio normalisation. Unless the paired Excel field is known to share a directly compatible scale, these checks are cross-platform concordance assessments rather than direct unit-matched comparisons.",
            })
    if not rows:
        return pd.DataFrame({"Note": ["No overlapping Excel/BAM gene pairs had enough matched primary RNA samples for concordance testing."]})
    return pd.DataFrame(rows).sort_values(["Gene", "Context"]).reset_index(drop=True)


def load_optional_workbook_sheet(path: Path, sheet_name: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        xl = pd.ExcelFile(path)
    except Exception:
        return pd.DataFrame()
    if sheet_name not in xl.sheet_names:
        return pd.DataFrame()
    return pd.read_excel(path, sheet_name=sheet_name)


def compute_immune_response_score(df: pd.DataFrame) -> tuple[pd.Series, list[str]]:
    genes = []
    score_cols = {}
    for gene in IMMUNE_SIGNATURE_GENES:
        candidates = [
            gene,
            f"bam_ratio__{gene}",
            f"excel_panel__{gene}",
            f"excel_total__{gene}",
            f"excel_iso__{gene}",
        ]
        source = next((col for col in candidates if col in df.columns), None)
        if source is not None:
            genes.append(gene)
            score_cols[gene] = source
    score = pd.Series(np.nan, index=df.index, dtype=float)
    if len(genes) < 2:
        return score, genes
    zscores = pd.DataFrame({gene: pd.to_numeric(df[source], errors="coerce") for gene, source in score_cols.items()})
    zscores = pd.DataFrame({
        gene: (
            (col - col.mean()) / col.std(ddof=0)
            if pd.notna(col.std(ddof=0)) and col.std(ddof=0) != 0
            else pd.Series(np.nan, index=df.index)
        )
        for gene, col in zscores.items()
    })
    score = zscores.mean(axis=1, skipna=True)
    score[zscores.notna().sum(axis=1) < 2] = np.nan
    return score, genes


def _analysis_frame_for_exposure(base_df: pd.DataFrame, dosage_df: pd.DataFrame, id_col: str, exposure_id: object) -> pd.DataFrame:
    if dosage_df.empty or id_col not in dosage_df.columns:
        return pd.DataFrame()
    g = dosage_df[dosage_df[id_col].astype(str) == str(exposure_id)].copy()
    if g.empty:
        return pd.DataFrame()
    merge_cols = [c for c in ["analysis_sample_id", "snp_code", "Dosage"] if c in g.columns]
    if not merge_cols:
        return pd.DataFrame()
    dedup_keys = [c for c in ["analysis_sample_id"] if c in merge_cols]
    g = g[merge_cols].drop_duplicates(dedup_keys if dedup_keys else None)
    merged = base_df.merge(g, on="snp_code", how="inner")
    if "analysis_sample_id_x" in merged.columns and "analysis_sample_id_y" in merged.columns:
        merged["analysis_sample_id"] = merged["analysis_sample_id_x"].fillna(merged["analysis_sample_id_y"])
    elif "analysis_sample_id_x" in merged.columns:
        merged["analysis_sample_id"] = merged["analysis_sample_id_x"]
    elif "analysis_sample_id_y" in merged.columns:
        merged["analysis_sample_id"] = merged["analysis_sample_id_y"]
    return merged


def _fit_ols(data: pd.DataFrame, y_col: str, x_cols: list[str]) -> tuple[dict[str, float] | None, pd.DataFrame]:
    cols = [y_col] + x_cols
    work = data[cols].copy()
    for col in cols:
        work[col] = pd.to_numeric(work[col], errors="coerce")
    work = work.dropna()
    if work.empty or len(work) < 8 or work[y_col].nunique() < 2:
        return None, work
    design = np.column_stack([np.ones(len(work))] + [work[col].to_numpy(dtype=float) for col in x_cols])
    response = work[y_col].to_numpy(dtype=float)
    try:
        beta, _, _, _ = np.linalg.lstsq(design, response, rcond=None)
    except Exception:
        return None, work
    keys = ["Intercept"] + x_cols
    return dict(zip(keys, map(float, beta))), work


def build_immune_response_summary(na_df: pd.DataFrame) -> pd.DataFrame:
    if na_df.empty or "Immune_Response_Score" not in na_df.columns:
        return pd.DataFrame({"Note": ["No immune-response score could be derived from the RNA expression matrix."]})
    rows = []
    for scope, sub in [
        ("All_Primary", na_df),
        ("Breast", na_df[na_df["cohort"] == "Breast"].copy()),
        ("Endometrial", na_df[na_df["cohort"] == "Endometrial"].copy()),
    ]:
        work = sub[pd.to_numeric(sub["Immune_Response_Score"], errors="coerce").notna()].copy()
        if work.empty:
            continue
        rows.append({
            "Scope": scope,
            "N": len(work),
            "Median_Score": work["Immune_Response_Score"].median(),
            "IQR_Lo": work["Immune_Response_Score"].quantile(0.25),
            "IQR_Hi": work["Immune_Response_Score"].quantile(0.75),
            "Available_Genes": "; ".join([
                gene for gene in IMMUNE_SIGNATURE_GENES
                if any(col in work.columns for col in [gene, f"bam_ratio__{gene}", f"excel_panel__{gene}", f"excel_total__{gene}", f"excel_iso__{gene}"])
            ]),
        })
        for label, var in [("GSDMB", "GSDMB"), ("Pyroptotic_Fraction", "pyroptotic_isoform_fraction"), ("Non_Pyroptotic_Fraction", "non_pyroptotic_isoform_fraction")]:
            if var not in work.columns:
                continue
            x = pd.to_numeric(work["Immune_Response_Score"], errors="coerce")
            y = pd.to_numeric(work[var], errors="coerce")
            t = pd.concat([x, y], axis=1).dropna()
            if len(t) < 5 or t.iloc[:, 0].nunique() < 2 or t.iloc[:, 1].nunique() < 2:
                continue
            rho, p = stats.spearmanr(t.iloc[:, 0], t.iloc[:, 1], nan_policy="omit")
            rows.append({
                "Scope": scope,
                "N": len(t),
                "Median_Score": np.nan,
                "IQR_Lo": np.nan,
                "IQR_Hi": np.nan,
                "Available_Genes": f"Correlation_vs_{label}",
                "Spearman_Rho": rho,
                "Spearman_P": p,
            })
    return pd.DataFrame(rows)


def build_bootstrap_stability_summary(na_df: pd.DataFrame, dosage_df: pd.DataFrame, assoc_df: pd.DataFrame, id_col: str, exposure_type: str, top_n: int = 20, n_boot: int = 200, random_state: int = 42) -> pd.DataFrame:
    if assoc_df.empty or "P_Value" not in assoc_df.columns or "RNA_Var" not in assoc_df.columns:
        return pd.DataFrame({"Note": [f"No {exposure_type} association rows were available for bootstrap stability checks."]})
    candidates = assoc_df.copy()
    if "FDR_Sig" in candidates.columns and candidates["FDR_Sig"].fillna(False).any():
        candidates = candidates[candidates["FDR_Sig"].fillna(False)].copy()
    else:
        candidates = candidates[candidates["Nominal_Sig"].fillna(False)].copy() if "Nominal_Sig" in candidates.columns else candidates.copy()
    if candidates.empty:
        return pd.DataFrame({"Note": [f"No {exposure_type} associations met the nominal/FDR filter for bootstrap stability."]})
    candidates = candidates.sort_values(["P_Value", "N"], ascending=[True, False]).head(top_n).copy()
    rng = np.random.default_rng(random_state)
    rows = []
    for _, row in candidates.iterrows():
        base = subset_for_context(na_df, row.get("Context", "All_Primary"))
        frame = _analysis_frame_for_exposure(base, dosage_df, id_col, row[id_col] if id_col in row.index else row.get("SNP_ID", row.get("Exposure_ID")))
        if frame.empty or row["RNA_Var"] not in frame.columns:
            continue
        t = frame[[row["RNA_Var"], "Dosage"]].copy()
        t[row["RNA_Var"]] = pd.to_numeric(t[row["RNA_Var"]], errors="coerce")
        t["Dosage"] = pd.to_numeric(t["Dosage"], errors="coerce")
        t = t.dropna()
        if len(t) < 8 or t["Dosage"].nunique() < 2 or t[row["RNA_Var"]].nunique() < 2:
            continue
        orig_slope = stats.linregress(t["Dosage"], t[row["RNA_Var"]]).slope
        boot_slopes = []
        for _ in range(n_boot):
            sample = t.iloc[rng.integers(0, len(t), len(t))].copy()
            if sample["Dosage"].nunique() < 2 or sample[row["RNA_Var"]].nunique() < 2:
                continue
            try:
                boot_slopes.append(stats.linregress(sample["Dosage"], sample[row["RNA_Var"]]).slope)
            except Exception:
                continue
        if not boot_slopes:
            continue
        boot_arr = np.asarray(boot_slopes, dtype=float)
        rows.append({
            "Exposure_Type": exposure_type,
            id_col: row[id_col] if id_col in row.index else row.get("SNP_ID", row.get("Exposure_ID")),
            "Context": row.get("Context", np.nan),
            "RNA_Var": row["RNA_Var"],
            "N": len(t),
            "Original_Slope": orig_slope,
            "Bootstrap_Median_Slope": float(np.nanmedian(boot_arr)),
            "Bootstrap_Lo": float(np.nanpercentile(boot_arr, 2.5)),
            "Bootstrap_Hi": float(np.nanpercentile(boot_arr, 97.5)),
            "Bootstrap_Sign_Concordance": float(np.mean(np.sign(boot_arr) == np.sign(orig_slope))),
            "Bootstrap_N": len(boot_arr),
            "CI_Crosses_Zero": not (np.nanpercentile(boot_arr, 2.5) > 0 or np.nanpercentile(boot_arr, 97.5) < 0),
            "Source_P_Value": row.get("P_Value", np.nan),
            "Source_FDR": row.get("FDR", np.nan),
        })
    if not rows:
        return pd.DataFrame({"Note": [f"No {exposure_type} association rows had enough data for bootstrap stability estimation."]})
    return pd.DataFrame(rows).sort_values(["Source_FDR", "Source_P_Value", "Bootstrap_Sign_Concordance"], ascending=[True, True, False]).reset_index(drop=True)


def build_mediation_summary(na_df: pd.DataFrame, dosage_df: pd.DataFrame, assoc_df: pd.DataFrame, id_col: str, exposure_type: str, top_n: int = 10, n_boot: int = 200, random_state: int = 42) -> pd.DataFrame:
    if assoc_df.empty or "RNA_Var" not in assoc_df.columns:
        return pd.DataFrame({"Note": [f"No {exposure_type} association rows were available for mediation analysis."]})
    candidates = assoc_df.copy()
    if "FDR_Sig" in candidates.columns and candidates["FDR_Sig"].fillna(False).any():
        candidates = candidates[candidates["FDR_Sig"].fillna(False)].copy()
    else:
        candidates = candidates[candidates["Nominal_Sig"].fillna(False)].copy() if "Nominal_Sig" in candidates.columns else candidates.copy()
    if candidates.empty:
        return pd.DataFrame({"Note": [f"No {exposure_type} associations met the nominal/FDR filter for mediation analysis."]})
    candidates = candidates.sort_values(["P_Value", "N"], ascending=[True, False]).head(top_n).copy()
    mediator_map = {
        "Breast": ["GSDMB", "pyroptotic_isoform_fraction", "non_pyroptotic_isoform_fraction", "Immune_Response_Score"],
        "Endometrial": ["GSDMB", "pyroptotic_isoform_fraction", "non_pyroptotic_isoform_fraction", "Immune_Response_Score"],
        "All_Primary": ["GSDMB", "pyroptotic_isoform_fraction", "non_pyroptotic_isoform_fraction", "Immune_Response_Score"],
    }
    outcome_map = {
        "Breast": ["BREAST_RECURRENCE_DERIVED", "BREAST_METASTASIS_DERIVED", "BREAST_EXITUS_DERIVED", "BREAST_OS_MONTHS_DERIVED"],
        "Endometrial": ["ENDO_PD_BIN", "ENDO_EXITUS_BIN", "ENDO_OS_MONTHS_DERIVED", "ENDO_PFS_MONTHS_DERIVED"],
        "All_Primary": ["BREAST_RECURRENCE_DERIVED", "BREAST_METASTASIS_DERIVED", "BREAST_EXITUS_DERIVED", "BREAST_OS_MONTHS_DERIVED", "ENDO_PD_BIN", "ENDO_EXITUS_BIN", "ENDO_OS_MONTHS_DERIVED", "ENDO_PFS_MONTHS_DERIVED"],
    }
    rng = np.random.default_rng(random_state)
    rows = []
    for _, row in candidates.iterrows():
        context = row.get("Context", "All_Primary")
        base = subset_for_context(na_df, context)
        frame = _analysis_frame_for_exposure(base, dosage_df, id_col, row[id_col] if id_col in row.index else row.get("SNP_ID", row.get("Exposure_ID")))
        if frame.empty:
            continue
        covars = [c for c in ["canon__age", "canon__bmi"] if c in frame.columns]
        if row["RNA_Var"] not in frame.columns:
            continue
        mediators = [m for m in mediator_map.get(context, mediator_map["All_Primary"]) if m in frame.columns]
        outcomes = [o for o in outcome_map.get(context, outcome_map["All_Primary"]) if o in frame.columns]
        for mediator in mediators:
            for outcome in outcomes:
                work = frame[[row["RNA_Var"], mediator, outcome, "Dosage"] + covars].copy()
                for col in work.columns:
                    work[col] = pd.to_numeric(work[col], errors="coerce")
                work = work.dropna()
                if len(work) < 8 or work["Dosage"].nunique() < 2 or work[mediator].nunique() < 2 or work[outcome].nunique() < 2:
                    continue
                a_fit, a_work = _fit_ols(work, mediator, ["Dosage"] + covars)
                b_fit, b_work = _fit_ols(work, outcome, ["Dosage", mediator] + covars)
                c_fit, _ = _fit_ols(work, outcome, ["Dosage"] + covars)
                if a_fit is None or b_fit is None or c_fit is None:
                    continue
                indirect = a_fit["Dosage"] * b_fit[mediator]
                boot_indirect = []
                for _ in range(n_boot):
                    sample = work.iloc[rng.integers(0, len(work), len(work))].copy()
                    if sample["Dosage"].nunique() < 2 or sample[mediator].nunique() < 2 or sample[outcome].nunique() < 2:
                        continue
                    a_b, _ = _fit_ols(sample, mediator, ["Dosage"] + covars)
                    b_b, _ = _fit_ols(sample, outcome, ["Dosage", mediator] + covars)
                    if a_b is None or b_b is None:
                        continue
                    boot_indirect.append(a_b["Dosage"] * b_b[mediator])
                if not boot_indirect:
                    continue
                boot_arr = np.asarray(boot_indirect, dtype=float)
                rows.append({
                    "Exposure_Type": exposure_type,
                    id_col: row[id_col] if id_col in row.index else row.get("SNP_ID", row.get("Exposure_ID")),
                    "Context": context,
                    "RNA_Var": row["RNA_Var"],
                    "Mediator": mediator,
                    "Outcome": outcome,
                    "N": len(work),
                    "Direct_Effect": b_fit["Dosage"],
                    "Mediator_a": a_fit["Dosage"],
                    "Mediator_b": b_fit[mediator],
                    "Total_Effect": c_fit["Dosage"],
                    "Indirect_Effect": indirect,
                    "Indirect_Lo": float(np.nanpercentile(boot_arr, 2.5)),
                    "Indirect_Hi": float(np.nanpercentile(boot_arr, 97.5)),
                    "Bootstrap_N": len(boot_arr),
                    "Indirect_Sign_Concordance": float(np.mean(np.sign(boot_arr) == np.sign(indirect))),
                    "Source_P_Value": row.get("P_Value", np.nan),
                    "Source_FDR": row.get("FDR", np.nan),
                    "Covariates": "age + bmi" if covars else "none",
                })
    if not rows:
        return pd.DataFrame({"Note": [f"No {exposure_type} rows yielded a stable mediation fit."]})
    return pd.DataFrame(rows).sort_values(["Source_FDR", "Source_P_Value", "Indirect_Sign_Concordance"], ascending=[True, True, False]).reset_index(drop=True)


def main():
    # Stage 23 is a workbook-first integration step, so it fails early when an
    # expected upstream result is missing rather than silently producing a
    # partial workbook.
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    optional_inputs = {"stage22", "stage17b_breast", "stage17b_endo", "stage17b_exposure_catalogue"}
    for name, path in PATHS.items():
        if name in optional_inputs:
            continue
        if not path.exists():
            raise FileNotFoundError(f"Missing required input for stage 23: {name} -> {path}")

    qc_manifest = pd.read_csv(PATHS["rna_qc_manifest"], sep="\t")
    excel_na, excel_vars, excel_raw = load_excel_na(PATHS["stage20"])
    bam_na, bam_vars = load_bam_na(PATHS["stage20b"])
    manifest = build_na_manifest(excel_na, bam_na, qc_manifest)
    na = merge_na_layers(excel_na, bam_na)
    if "Immune_Response_Score" not in na.columns:
        na["Immune_Response_Score"], immune_genes = compute_immune_response_score(na)
    else:
        immune_genes = [gene for gene in IMMUNE_SIGNATURE_GENES if gene in na.columns]
    na_catalog = build_na_variable_catalog(excel_vars, bam_vars)
    sig_snps = load_significant_snps(PATHS["stage12"])
    sig_mutations = load_significant_mutations(PATHS["stage17b_breast"], PATHS["stage17b_endo"], PATHS["stage17b_exposure_catalogue"])
    sig_haps = load_significant_haplotypes(PATHS["stage22"], PATHS["stage15"])
    therapy_response_assoc = load_optional_workbook_sheet(PATHS["stage20"], "therapy_response_assoc")
    phased, sample_cols, snps_used = load_phased_matrix(PATHS["phased"], PATHS["stage15"])
    snp_dosage, sig_snps_summary = build_significant_snp_dosage(phased, sample_cols, sig_snps)
    mutation_dosage, sig_mutations_summary = build_significant_mutation_dosage(sig_mutations, na.get("snp_code", pd.Series(dtype=object)).dropna().astype(str).tolist())
    hap_dosage, hap_freq_checks = build_significant_haplotype_dosage(phased, sample_cols, PATHS["stage15"], sig_haps, snps_used)

    tested_snp_codes = set(snp_dosage["snp_code"].dropna()) if not snp_dosage.empty else set()
    tested_mutation_codes = set(mutation_dosage.loc[pd.to_numeric(mutation_dosage.get("Dosage"), errors="coerce") > 0, "snp_code"].dropna()) if not mutation_dosage.empty else set()
    tested_hap_codes = set(hap_dosage["snp_code"].dropna()) if not hap_dosage.empty else set()
    manifest["has_sig_snp_dosage"] = manifest["snp_code"].isin(tested_snp_codes)
    manifest["has_sig_mutation_exposure"] = manifest["snp_code"].isin(tested_mutation_codes)
    manifest["has_sig_haplotype_dosage"] = manifest["snp_code"].isin(tested_hap_codes)
    manifest["sample_match_logic"] = "Exact canonical snp_code join after stage-20 RNA matching and BAM filename canonicalisation"

    primary_na = primary_subset(na)
    snp_assoc = regression_rows(primary_na, snp_dosage, sig_snps_summary, "SNP", na_catalog)
    mutation_assoc = regression_rows(primary_na, mutation_dosage, sig_mutations_summary, "Mutation", na_catalog)
    hap_assoc = regression_rows(primary_na, hap_dosage, sig_haps, "Haplotype", na_catalog)
    snp_isoform_profile_assoc = isoform_profile_rows(primary_na, snp_dosage, sig_snps_summary, "SNP")
    mutation_isoform_profile_assoc = isoform_profile_rows(primary_na, mutation_dosage, sig_mutations_summary, "Mutation")
    hap_isoform_profile_assoc = isoform_profile_rows(primary_na, hap_dosage, sig_haps, "Haplotype")
    immune_summary = build_immune_response_summary(na)
    bootstrap_stability_summary = pd.concat([
        build_bootstrap_stability_summary(primary_na, snp_dosage, snp_assoc, "SNP_ID", "SNP", top_n=15, n_boot=200),
        build_bootstrap_stability_summary(primary_na, hap_dosage, hap_assoc, "Exposure_ID", "Haplotype", top_n=15, n_boot=200),
    ], ignore_index=True)
    mediation_summary = pd.concat([
        build_mediation_summary(primary_na, snp_dosage, snp_assoc, "SNP_ID", "SNP", top_n=10, n_boot=200),
        build_mediation_summary(primary_na, hap_dosage, hap_assoc, "Exposure_ID", "Haplotype", top_n=10, n_boot=200),
    ], ignore_index=True)
    summary = summarise_hits(snp_assoc, mutation_assoc, hap_assoc)
    gsdmb_isoform_answer = build_gsdmb_isoform_answer(snp_assoc, mutation_assoc, hap_assoc)
    gsdmb_isoform_summary = build_gsdmb_isoform_summary(gsdmb_isoform_answer)
    unmatched_excel = excel_raw[excel_raw.get("match_status", "") != "matched"].copy()
    inventory = source_inventory(manifest, na_catalog)
    gate_summary = gate_usage_summary(manifest)
    excel_bam_concordance = build_excel_bam_concordance(na)

    with pd.ExcelWriter(OUT_DIR / "23_RNA_Genetic_Integration.xlsx", engine="openpyxl") as writer:
        pd.DataFrame({
            "Notes": [
                "Stage 23 links cleaned RNA data to current significant SNP, mutation, and haplotype signals.",
                "Excel-derived variables preserve isoform and panel-expression provenance with explicit prefixes.",
                "BAM-derived variables use stage-20b ratio-normalised panel expression with bam_ratio__ prefixes.",
                "Mutation and rare-variant burden rows retain Component_Genes, Component_Variants, and Component_Variant_Count from stage 17b.",
                "If stage 22 has not been run, significant haplotypes are recovered directly from the stage-15 association workbook so stage 23 remains part of the core pipeline.",
                "Associations are restricted to primary tumour samples passing the exploratory RNA QC gate, while the manifest still preserves which rows were strict-ready versus exploratory-only.",
                "Use the GSDMB Isoform Summary sheet for the collapsed supervisor-facing answer, and gsdmb_isoform_answer for the fully expanded row-level output.",
                "Therapy-response associations are carried through from stage 20 so the final integration workbook keeps treatment-related RNA signals in one place.",
                "Immune-response scoring is derived from the panel-gene layer using the available immune markers; mediation and bootstrap summaries use age and BMI as the only covariates.",
                "Excel RNA values are treated as source-normalised assay outputs and preserved in their original workbook units; no extra repo-level rescaling is applied to those columns.",
                "Excel-vs-BAM comparison sheets therefore report cross-platform concordance unless directly compatible units are explicitly known for the paired variables.",
                "The primary RNA endpoint set is total GSDMB, GSDMB1-4 expression, exon-6-positive fraction/balance, exon-7-containing and exon-7-deficient fraction/balance, GSDMB2 fraction, total GSDMB/ERBB2, exon-6-positive GSDMB/ERBB2, isoform-profile correlations with measured genes, and their SNP/haplotype/clinical associations. GSDMB5 is documented but not tested unless a measured G5/GSDMB5 column is available.",
                "Isoform-profile association sheets use a permutation test on Aitchison CLR-transformed isoform proportions to ask whether an exposure shifts several isoforms jointly.",
                "FDR is applied within each exposure x context x RNA source family.",
            ]
        }).to_excel(writer, sheet_name="README", index=False)
        na.to_excel(writer, sheet_name="rna_cleaned", index=False)
        manifest.sort_values(["cohort", "analysis_group", "snp_code"], na_position="last").to_excel(writer, sheet_name="rna_match_manifest", index=False)
        na_catalog.to_excel(writer, sheet_name="rna_variable_catalog", index=False)
        inventory.to_excel(writer, sheet_name="rna_source_inventory", index=False)
        gate_summary.to_excel(writer, sheet_name="rna_gate_summary", index=False)
        excel_bam_concordance.to_excel(writer, sheet_name="excel_bam_concordance", index=False)
        sig_snps_summary.to_excel(writer, sheet_name="significant_snps", index=False)
        sig_mutations_summary.to_excel(writer, sheet_name="significant_mutations", index=False)
        sig_haps.to_excel(writer, sheet_name="significant_haplotypes", index=False)
        hap_freq_checks.to_excel(writer, sheet_name="haplotype_freq_checks", index=False)
        (therapy_response_assoc if not therapy_response_assoc.empty else pd.DataFrame({"Note": ["No therapy-response RNA associations were available in stage 20."]})).to_excel(writer, sheet_name="therapy_response_assoc", index=False)
        (snp_assoc if not snp_assoc.empty else pd.DataFrame({"Note": ["No significant-SNP RNA tests met the minimum thresholds."]})).to_excel(writer, sheet_name="snp_rna_assoc", index=False)
        (mutation_assoc if not mutation_assoc.empty else pd.DataFrame({"Note": ["No significant-mutation RNA tests met the minimum thresholds."]})).to_excel(writer, sheet_name="mutation_rna_assoc", index=False)
        (hap_assoc if not hap_assoc.empty else pd.DataFrame({"Note": ["No significant-haplotype RNA tests met the minimum thresholds."]})).to_excel(writer, sheet_name="haplotype_rna_assoc", index=False)
        (snp_isoform_profile_assoc if not snp_isoform_profile_assoc.empty else pd.DataFrame({"Note": ["No significant-SNP isoform-profile tests met the minimum thresholds."]})).to_excel(writer, sheet_name="snp_isoform_profile", index=False)
        (mutation_isoform_profile_assoc if not mutation_isoform_profile_assoc.empty else pd.DataFrame({"Note": ["No significant-mutation isoform-profile tests met the minimum thresholds."]})).to_excel(writer, sheet_name="mutation_isoform_profile", index=False)
        (hap_isoform_profile_assoc if not hap_isoform_profile_assoc.empty else pd.DataFrame({"Note": ["No significant-haplotype isoform-profile tests met the minimum thresholds."]})).to_excel(writer, sheet_name="haplotype_isoform_profile", index=False)
        immune_summary.to_excel(writer, sheet_name="immune_response_summary", index=False)
        (bootstrap_stability_summary if not bootstrap_stability_summary.empty else pd.DataFrame({"Note": ["No bootstrap stability rows could be estimated."]})).to_excel(writer, sheet_name="bootstrap_stability_summary", index=False)
        (mediation_summary if not mediation_summary.empty else pd.DataFrame({"Note": ["No mediation rows could be estimated."]})).to_excel(writer, sheet_name="mediation_summary", index=False)
        gsdmb_isoform_summary.to_excel(writer, sheet_name="GSDMB Isoform Summary", index=False)
        gsdmb_isoform_answer.to_excel(writer, sheet_name="gsdmb_isoform_answer", index=False)
        summary.to_excel(writer, sheet_name="expression_effect_summary", index=False)
        unmatched_excel.to_excel(writer, sheet_name="unmatched_excel_rows", index=False)

    if not snp_assoc.empty:
        plot_heatmap(snp_assoc.sort_values("P_Value"), "SNP_ID", OUT_DIR / "23_SNP_RNA_Heatmap.png", "Significant SNP vs RNA Associations")
    if not mutation_assoc.empty:
        plot_heatmap(mutation_assoc.sort_values("P_Value"), "Exposure_ID", OUT_DIR / "23_Mutation_RNA_Heatmap.png", "Significant Mutation vs RNA Associations")
    if not hap_assoc.empty:
        plot_heatmap(hap_assoc.sort_values("P_Value"), "Exposure_ID", OUT_DIR / "23_Haplotype_RNA_Heatmap.png", "Significant Haplotype vs RNA Associations")

    print("=== Stage 23 RNA integration summary ===")
    print(f"Collapsed Excel RNA samples      : {excel_na['snp_code'].nunique()}")
    print(f"BAM RNA samples                 : {bam_na['snp_code'].nunique()}")
    print(f"Primary RNA samples analysed    : {primary_na['snp_code'].nunique()}")
    print(f"Primary strict-ready rows       : {int(primary_na.get('RNA_QC_Analysis_Ready', False).fillna(False).sum())}")
    print(f"Primary exploratory-only rows   : {int((primary_na.get('RNA_QC_Exploratory_Ready', False).fillna(False) & ~primary_na.get('RNA_QC_Analysis_Ready', False).fillna(False)).sum())}")
    print(f"Significant SNPs retained       : {len(sig_snps_summary)}")
    print(f"Significant SNPs testable       : {int(sig_snps_summary['RNA_Testable'].sum())}")
    print(f"Significant mutations retained  : {len(sig_mutations_summary)}")
    print(f"Significant mutations testable  : {int(sig_mutations_summary.get('RNA_Testable', pd.Series(dtype=bool)).fillna(False).sum()) if not sig_mutations_summary.empty else 0}")
    print(f"Significant haplotypes retained : {len(sig_haps)}")
    print(f"SNP vs RNA tests                : {len(snp_assoc)}")
    print(f"Mutation vs RNA tests           : {len(mutation_assoc)}")
    print(f"Haplotype vs RNA tests          : {len(hap_assoc)}")
    print(f"SNP isoform-profile tests       : {len(snp_isoform_profile_assoc)}")
    print(f"Mutation isoform-profile tests  : {len(mutation_isoform_profile_assoc)}")
    print(f"Haplotype isoform-profile tests : {len(hap_isoform_profile_assoc)}")
    print(f"Therapy-response rows           : {len(therapy_response_assoc)}")
    print(f"Immune markers used             : {len(immune_genes)}")
    print(f"Immune-score summary rows       : {len(immune_summary)}")
    print(f"Bootstrap stability rows        : {len(bootstrap_stability_summary)}")
    print(f"Mediation summary rows          : {len(mediation_summary)}")
    print(f"Results workbook                : {OUT_DIR / '23_RNA_Genetic_Integration.xlsx'}")


if __name__ == "__main__":
    main()
