#!/usr/bin/env python3
"""Stage 23: integrate RNA measurements with the significant DNA findings.

This script is the final synthesis layer for thesis objective 2. It collects:

- the cleaned Excel-based isoform measurements from stage 20
- BAM-derived panel-gene ratios from stage 20b
- the significant SNP backbone from stage 12
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

try:
    import seaborn as sns
except ImportError:
    sns = None

ROOT = Path("/home/gadeaalonsoj/tfm")
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
    "stage15": RES / "15_haplotype_statistics" / "GSDMB_Haplotype_Results.xlsx",
    "phased": RES / "15_haplotype_phasing" / "phased_genotypes.tsv",
    "stage22": RES / "22_haplotype_first_interpretation" / "22_Haplotype_First_Interpretation.xlsx",
}

ISOFORM_COLS = [
    "G1", "G2", "G3", "G3b", "G4", "GSDMB",
    "pyroptotic_isoform_fraction", "non_pyroptotic_isoform_fraction",
    "G4_vs_total", "G2_vs_total",
]
QC_COLS = [
    "RNA_QC_Analysis_Ready", "RNA_QC_Exploratory_Ready", "RNA_QC_Final_Status",
    "RNA_QC_Exclusion_Reason", "RNA_QC_Eligibility_Note",
]
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
        (r"^RNA_SNP_(AT|EN|MN|MT-T|MT-N)[-_]?(\d+)", lambda m: f"SNP_{'MN' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_RNA_(AT|EN|MN|MT-T|MT-N)_(\d+)", lambda m: f"SNP_{'MN' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^SNP_(AT|EN|MN|MT-T|MT-N)_RNA_(\d+)", lambda m: f"SNP_{'MN' if m.group(1).upper() == 'MT-N' else m.group(1).upper()}_{m.group(2)}"),
        (r"^RNA_(AT|EN|MN)_(\d+)", lambda m: f"SNP_{m.group(1).upper()}_{m.group(2)}"),
        (r"^RNA_MT[-_]T_(\d+)", lambda m: f"SNP_MT-T_{m.group(1)}"),
        (r"^RNA_MT[-_]N_(\d+)", lambda m: f"SNP_MN_{m.group(1)}"),
    ]
    for pat, fmt in patterns:
        m = re.match(pat, t, re.I)
        if m:
            return fmt(m)
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
        return "excel_isoform", f"excel_iso__{clean_name(col)}", "isoform_expression"
    if col == "GSDMB EXPRESSION IN PANEL":
        return "excel_panel", f"excel_panel__{clean_name(col)}", "other_rna_quant"
    return "excel_panel", f"excel_panel__{clean_name(col)}", "other_rna_quant"


def expression_nonmissing(df, cols):
    if not cols:
        return pd.Series(0, index=df.index)
    return df[cols].apply(pd.to_numeric, errors="coerce").notna().sum(axis=1)


def load_excel_na(stage20_path):
    # Stage 20 can contain repeated rows per sample. We collapse to one
    # best-supported row per `snp_code` while keeping provenance columns that
    # reveal whether multiple source rows were encountered.
    raw = pd.read_excel(stage20_path, sheet_name="expression_cleaned")
    raw["snp_code"] = raw["snp_code"].astype(str).str.strip()
    raw = raw[raw["snp_code"].notna() & (raw["snp_code"] != "") & (raw["snp_code"] != "nan")].copy()
    value_cols = []
    for col in raw.columns:
        if col in META_COLS or col.startswith(META_PREFIXES):
            continue
        if re.match(r"^rs\d+(_[A-Za-z]+)?$", str(col), re.I):
            continue
        if re.match(r"^\d+_[ACGT]+>[ACGT]+$", str(col), re.I):
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
    vars_df = pd.DataFrame(var_rows).sort_values(["RNA_Source", "Integrated_Var"]).reset_index(drop=True)
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
        df = df[df["FDR_P_Value"].fillna(1) < FDR_THRESHOLD].copy()
        if df.empty:
            continue
        df["Significant_In"] = sheet
        parts.append(df)
    all_sig = pd.concat(parts, ignore_index=True)
    return (
        all_sig.groupby("SNP_ID", dropna=False)
        .agg(
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


def load_phased_matrix(phased_path, stage15_path):
    # The phased matrix lets us calculate SNP dosage and haplotype exposure for
    # the same samples represented in the RNA workbook.
    phased = pd.read_csv(phased_path, sep="\t", dtype=str)
    phased["POS"] = pd.to_numeric(phased["POS"], errors="coerce")
    phased["REF"] = phased["REF"].astype(str).str.upper()
    phased["ALT"] = phased["ALT"].astype(str).str.upper()
    snps_used = pd.read_excel(stage15_path, sheet_name="SNPs_Used").rename(columns={"POS_int": "POS", "rsID_clean": "SNP_ID"})
    snps_used["POS"] = pd.to_numeric(snps_used["POS"], errors="coerce")
    snps_used["REF"] = snps_used["REF"].astype(str).str.upper()
    snps_used["ALT"] = snps_used["ALT"].astype(str).str.upper()
    phased = phased.merge(snps_used[["POS", "REF", "ALT", "SNP_ID", "Gene"]], on=["POS", "REF", "ALT"], how="left")
    sample_cols = [c for c in phased.columns if c not in {"CHROM", "POS", "ID", "REF", "ALT", "SNP_ID", "Gene"}]
    return phased, sample_cols, snps_used


def build_significant_snp_dosage(phased, sample_cols, sig_snp_df):
    targets = set(sig_snp_df["SNP_ID"].dropna().astype(str))
    available = phased[phased["SNP_ID"].isin(targets)].copy()
    rows = []
    for _, row in available.iterrows():
        for col in sample_cols:
            code = sx(col)
            if not code:
                continue
            rows.append({
                "snp_code": code,
                "SNP_ID": row["SNP_ID"],
                "Gene": row.get("Gene", np.nan),
                "Dosage": gt_dose(row[col]),
                "GT": row[col],
            })
    dosage = pd.DataFrame(rows)
    tested = set(dosage["SNP_ID"].dropna().astype(str)) if not dosage.empty else set()
    summary = sig_snp_df.copy()
    summary["RNA_Testable"] = summary["SNP_ID"].astype(str).isin(tested)
    summary["RNA_Test_Note"] = np.where(summary["RNA_Testable"], "Phased/common backbone dosage available", "Significant in stage 12 but absent from phased/common backbone")
    return dosage, summary


def region_membership(stage15_path, snps_used):
    members = {"Full_Region": list(snps_used["SNP_ID"].dropna().astype(str))}
    blocks = pd.read_excel(stage15_path, sheet_name="LD_Blocks")
    for (thr, block), sub in blocks.groupby(["Threshold_Label", "Block"], dropna=False):
        members[f"LD_Block_{thr}_{block}"] = sub.sort_values("SNP_Order")["rsID"].dropna().astype(str).tolist()
    return members


def build_significant_haplotype_dosage(phased, sample_cols, stage15_path, hap_df, snps_used):
    members = region_membership(stage15_path, snps_used)
    row_lookup = {str(r["SNP_ID"]): r for _, r in phased.dropna(subset=["SNP_ID"]).iterrows()}
    rows = []
    checks = []
    for _, hap in hap_df.iterrows():
        region = str(hap["Region"])
        target = str(hap["Haplotype_String"])
        exposure = str(hap["Exposure_ID"])
        snp_ids = members.get(region, [])
        callable_n = 0
        alt_count = 0
        for col in sample_cols:
            code = sx(col)
            if not code:
                continue
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

def regression_rows(na_df, dosage_df, exposures_df, exposure_type, var_catalog):
    value_cols = var_catalog["Integrated_Var"].tolist()
    var_level_map = var_catalog.set_index("Integrated_Var")["RNA_Level"].to_dict()
    rows = []
    exposure_meta = exposures_df.copy()
    if exposure_type == "SNP":
        id_col = "SNP_ID"
        cohort_map = exposure_meta.set_index(id_col)["Source_Groups"].to_dict()
        meta_map = exposure_meta.set_index(id_col).to_dict(orient="index")
    else:
        id_col = "Exposure_ID"
        cohort_map = {r["Exposure_ID"]: r["Comparison"] for _, r in exposure_meta.iterrows()}
        meta_map = exposure_meta.set_index(id_col).to_dict(orient="index")

    for exposure_id, g in dosage_df.groupby(id_col, dropna=False):
        if exposure_type == "SNP":
            contexts = context_for_sources(cohort_map.get(exposure_id, ""))
        else:
            comp = s(cohort_map.get(exposure_id, "Breast"))
            contexts = ["Breast"] if comp == "Breast" else ["All_Primary"]
        for context in contexts:
            base = subset_for_context(na_df, context)
            merged = base.merge(g[["snp_code", "Dosage"]].drop_duplicates(), on="snp_code", how="inner")
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
                for key in ["Source_Groups", "LD_Block_ID", "Min_FDR", "Min_P", "Region", "Haplotype_ID", "Comparison", "Global_Freq", "Overlapping_Significant_SNPs", "Consistency_Summary"]:
                    if key in meta:
                        row[key] = meta[key]
                rows.append(row)
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return add_fdr(result, "P_Value", ["Exposure_Type", id_col, "Context", "RNA_Source"])


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
        sns.heatmap(piv, cmap="coolwarm", center=0, linewidths=0.4, linecolor="white", ax=ax)
    else:
        im = ax.imshow(piv.values, aspect="auto", cmap="coolwarm")
        ax.set_xticks(range(len(piv.columns)))
        ax.set_xticklabels(piv.columns, rotation=40, ha="right", fontsize=8)
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels(piv.index, fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def summarise_hits(snp_assoc, hap_assoc):
    rows = []
    for _, df, id_col in [("SNP", snp_assoc, "SNP_ID"), ("Haplotype", hap_assoc, "Exposure_ID")]:
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
        "Source_Groups", "Region", "Haplotype_ID", "Comparison", "Overlapping_Significant_SNPs",
    ] if c in out.columns]
    return out[keep].sort_values(["Summary_Status", "P_Value", "N"], ascending=[True, True, False]).reset_index(drop=True)


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


def main():
    # Stage 23 is a workbook-first integration step, so it fails early when an
    # expected upstream result is missing rather than silently producing a
    # partial workbook.
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, path in PATHS.items():
        if name == "stage22":
            continue
        if not path.exists():
            raise FileNotFoundError(f"Missing required input for stage 23: {name} -> {path}")

    qc_manifest = pd.read_csv(PATHS["rna_qc_manifest"], sep="\t")
    excel_na, excel_vars, excel_raw = load_excel_na(PATHS["stage20"])
    bam_na, bam_vars = load_bam_na(PATHS["stage20b"])
    manifest = build_na_manifest(excel_na, bam_na, qc_manifest)
    na = merge_na_layers(excel_na, bam_na)
    na_catalog = build_na_variable_catalog(excel_vars, bam_vars)
    sig_snps = load_significant_snps(PATHS["stage12"])
    sig_haps = load_significant_haplotypes(PATHS["stage22"], PATHS["stage15"])
    phased, sample_cols, snps_used = load_phased_matrix(PATHS["phased"], PATHS["stage15"])
    snp_dosage, sig_snps_summary = build_significant_snp_dosage(phased, sample_cols, sig_snps)
    hap_dosage, hap_freq_checks = build_significant_haplotype_dosage(phased, sample_cols, PATHS["stage15"], sig_haps, snps_used)

    tested_snp_codes = set(snp_dosage["snp_code"].dropna()) if not snp_dosage.empty else set()
    tested_hap_codes = set(hap_dosage["snp_code"].dropna()) if not hap_dosage.empty else set()
    manifest["has_sig_snp_dosage"] = manifest["snp_code"].isin(tested_snp_codes)
    manifest["has_sig_haplotype_dosage"] = manifest["snp_code"].isin(tested_hap_codes)
    manifest["sample_match_logic"] = "Exact canonical snp_code join after stage-20 RNA matching and BAM filename canonicalisation"

    primary_na = primary_subset(na)
    snp_assoc = regression_rows(primary_na, snp_dosage, sig_snps_summary, "SNP", na_catalog)
    hap_assoc = regression_rows(primary_na, hap_dosage, sig_haps, "Haplotype", na_catalog)
    summary = summarise_hits(snp_assoc, hap_assoc)
    unmatched_excel = excel_raw[excel_raw.get("match_status", "") != "matched"].copy()
    inventory = source_inventory(manifest, na_catalog)
    gate_summary = gate_usage_summary(manifest)
    excel_bam_concordance = build_excel_bam_concordance(na)

    with pd.ExcelWriter(OUT_DIR / "23_RNA_Genetic_Integration.xlsx", engine="openpyxl") as writer:
        pd.DataFrame({
            "Notes": [
                "Stage 23 links cleaned RNA data to current significant SNP and haplotype signals.",
                "Excel-derived variables preserve isoform and panel-expression provenance with explicit prefixes.",
                "BAM-derived variables use stage-20b ratio-normalised panel expression with bam_ratio__ prefixes.",
                "If stage 22 has not been run, significant haplotypes are recovered directly from the stage-15 association workbook so stage 23 remains part of the core pipeline.",
                "Associations are restricted to primary tumour samples passing the exploratory RNA QC gate, while the manifest still preserves which rows were strict-ready versus exploratory-only.",
                "Excel RNA values are treated as source-normalised assay outputs and preserved in their original workbook units; no extra repo-level rescaling is applied to those columns.",
                "Excel-vs-BAM comparison sheets therefore report cross-platform concordance unless directly compatible units are explicitly known for the paired variables.",
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
        sig_haps.to_excel(writer, sheet_name="significant_haplotypes", index=False)
        hap_freq_checks.to_excel(writer, sheet_name="haplotype_freq_checks", index=False)
        (snp_assoc if not snp_assoc.empty else pd.DataFrame({"Note": ["No significant-SNP RNA tests met the minimum thresholds."]})).to_excel(writer, sheet_name="snp_rna_assoc", index=False)
        (hap_assoc if not hap_assoc.empty else pd.DataFrame({"Note": ["No significant-haplotype RNA tests met the minimum thresholds."]})).to_excel(writer, sheet_name="haplotype_rna_assoc", index=False)
        summary.to_excel(writer, sheet_name="expression_effect_summary", index=False)
        unmatched_excel.to_excel(writer, sheet_name="unmatched_excel_rows", index=False)

    if not snp_assoc.empty:
        plot_heatmap(snp_assoc.sort_values("P_Value"), "SNP_ID", OUT_DIR / "23_SNP_RNA_Heatmap.png", "Significant SNP vs RNA associations")
    if not hap_assoc.empty:
        plot_heatmap(hap_assoc.sort_values("P_Value"), "Exposure_ID", OUT_DIR / "23_Haplotype_RNA_Heatmap.png", "Significant haplotype vs RNA associations")

    print("=== Stage 23 RNA integration summary ===")
    print(f"Collapsed Excel RNA samples      : {excel_na['snp_code'].nunique()}")
    print(f"BAM RNA samples                 : {bam_na['snp_code'].nunique()}")
    print(f"Primary RNA samples analysed    : {primary_na['snp_code'].nunique()}")
    print(f"Primary strict-ready rows       : {int(primary_na.get('RNA_QC_Analysis_Ready', False).fillna(False).sum())}")
    print(f"Primary exploratory-only rows   : {int((primary_na.get('RNA_QC_Exploratory_Ready', False).fillna(False) & ~primary_na.get('RNA_QC_Analysis_Ready', False).fillna(False)).sum())}")
    print(f"Significant SNPs retained       : {len(sig_snps_summary)}")
    print(f"Significant SNPs testable       : {int(sig_snps_summary['RNA_Testable'].sum())}")
    print(f"Significant haplotypes retained : {len(sig_haps)}")
    print(f"SNP vs RNA tests                : {len(snp_assoc)}")
    print(f"Haplotype vs RNA tests          : {len(hap_assoc)}")
    print(f"Results workbook                : {OUT_DIR / '23_RNA_Genetic_Integration.xlsx'}")


if __name__ == "__main__":
    main()





