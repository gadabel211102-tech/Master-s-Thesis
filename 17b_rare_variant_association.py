#!/usr/bin/env python3
"""Rare-variant / mutation association analysis for the GSDMB thesis workflow.

This stage complements script 17 by focusing on variants that are not common
SNPs in gnomAD NFE. The default exposure definition is:

  rare / mutation = not common in gnomAD NFE exome or genome
                    OR missing both NFE frequencies

The script keeps the conservative, common-SNP association workflow intact and
adds a parallel mutation-focused workbook with four main views:

1. Single rare-variant frequencies and tumour-vs-control enrichment.
2. Gene-level rare-variant burden frequencies.
3. Recurrent rare-variant pair frequencies.
4. Tumour-only clinical associations for single variants, gene burdens and
   simple burden flags (any rare variant, multiple rare variants).

The analytical unit is the sequenced sample manifest from the harmonised master,
so non-carriers are defined explicitly from all sequenced samples rather than
only from rows present in the annotated variant workbook.
"""

from __future__ import annotations

import argparse
import importlib.util
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, false_discovery_control, fisher_exact, mannwhitneyu

from association_runtime import script17b_defaults
from pipeline_utils import combine_gnomad_nfe, common_nfe_variant_mask
from sample_identity_utils import attach_analysis_sample_ids, build_analysis_sample_map
from pipeline_validation import print_validation_summary, validate_file_exists
from variant_label_utils import build_variant_display_table as shared_build_variant_display_table


def _load_script17_module():
    """Load script 17 dynamically so we can reuse validated helpers."""
    path = Path(__file__).with_name("17_snp_association.py")
    spec = importlib.util.spec_from_file_location("script17_module", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load helper module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S17 = _load_script17_module()

DEFAULTS = script17b_defaults()
DEFAULT_ANNOT = DEFAULTS["annot"]
DEFAULT_MASTER = DEFAULTS["master"]
DEFAULT_PHASED = DEFAULTS["phased"]
DEFAULT_OUT = DEFAULTS["out_dir"]
DEFAULT_MANIFESTS = DEFAULTS.get("manifests", {})
DEFAULT_FDR = float(DEFAULTS.get("fdr_threshold", 0.10))
DEFAULT_MIN_CARRIERS = int(DEFAULTS.get("min_carriers", 5))
DEFAULT_RARE_AF = float(DEFAULTS.get("rare_af_threshold", 0.01))
DEFAULT_MIN_PAIR_CARRIERS = int(DEFAULTS.get("min_pair_carriers", 3))

MIN_COMPARISON_CARRIERS = 3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rare-variant association analysis alongside script 17")
    parser.add_argument("--annot", default=str(DEFAULT_ANNOT), help="Annotated workbook from stage 07")
    parser.add_argument("--master", default=str(DEFAULT_MASTER), help="Harmonised clinical master workbook")
    parser.add_argument("--out_dir", default=str(DEFAULT_OUT), help="Output directory for stage 17b")
    parser.add_argument("--phased", default=str(DEFAULT_PHASED), help="Phased genotype backbone for exact sample-identity collapsing")
    parser.add_argument("--sheet", default="Biological_Annotations", help="Annotated workbook sheet name")
    parser.add_argument("--rare-af-threshold", type=float, default=DEFAULT_RARE_AF, help="Rare threshold for combined gnomAD NFE AF")
    parser.add_argument("--fdr-threshold", type=float, default=DEFAULT_FDR, help="Benjamini-Hochberg FDR threshold")
    parser.add_argument("--min-carriers", type=int, default=DEFAULT_MIN_CARRIERS, help="Minimum carriers/non-carriers for clinical tests")
    parser.add_argument("--min-pair-carriers", type=int, default=DEFAULT_MIN_PAIR_CARRIERS, help="Minimum carriers to retain a recurrent rare-variant pair")
    parser.add_argument("--exclude-missing-af", action="store_true", help="Exclude variants with missing gnomAD NFE AF from the rare set")
    parser.add_argument("--manifest-endo-tumour", default=str(DEFAULT_MANIFESTS.get("endometrium-tumour", "")))
    parser.add_argument("--manifest-endo-normal", default=str(DEFAULT_MANIFESTS.get("endometrium-normal", "")))
    parser.add_argument("--manifest-breast-tumour", default=str(DEFAULT_MANIFESTS.get("breast-tumour", "")))
    parser.add_argument("--manifest-breast-normal", default=str(DEFAULT_MANIFESTS.get("breast-normal", "")))
    parser.add_argument("--no-manifests", action="store_true", help="Skip manifest filtering and rely on the master workbook only")
    return parser.parse_args()


def _normalise_gt(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .replace({"0|0": "0/0", "0|1": "0/1", "1|0": "1/0", "1|1": "1/1"})
    )


def _first_nonempty(series: pd.Series, fallback: str = "") -> str:
    for value in series:
        text = str(value).strip()
        if text and text.lower() != "nan":
            return text
    return fallback


def _safe_join(values: Iterable[object], sep: str = "; ") -> str:
    cleaned = sorted({str(v).strip() for v in values if str(v).strip() and str(v).strip().lower() != "nan"})
    return sep.join(cleaned)


def _join_component_lists(values: Iterable[object]) -> str:
    parts = []
    for value in values:
        parts.extend(part.strip() for part in str(value).split(";") if part.strip())
    return _safe_join(parts)


def _component_count(values: Iterable[object]) -> int:
    joined = _join_component_lists(values)
    return len([part for part in joined.split("; ") if part])


def _build_variant_key(df: pd.DataFrame) -> pd.Series:
    chrom = df["CHROM"].astype(str).str.strip()
    pos = pd.to_numeric(df["POS"], errors="coerce").fillna(-1).astype(int).astype(str)
    ref = df["REF"].astype(str).str.strip().str.upper()
    alt = df["ALT"].astype(str).str.strip().str.upper()
    return chrom + ":" + pos + "_" + ref + ">" + alt


def _pair_label(a: str, b: str, label_lookup: dict[str, str]) -> str:
    left = label_lookup.get(a, a)
    right = label_lookup.get(b, b)
    return f"{left} + {right}"


def _exposure_meta_defaults(frame: pd.DataFrame) -> dict[str, object]:
    return {
        "Variant_ID": _first_nonempty(frame.get("Variant_ID", pd.Series([], dtype=object))),
        "rsID": _first_nonempty(frame.get("rsID", pd.Series([], dtype=object))),
        "Variant_Label": _first_nonempty(frame.get("Variant_Label", pd.Series([], dtype=object))),
        "Genomic_Label": _first_nonempty(frame.get("Genomic_Label", pd.Series([], dtype=object))),
        "Label_Source": _first_nonempty(frame.get("Label_Source", pd.Series([], dtype=object))),
    }


EXPOSURE_COLUMNS = [
    "analysis_sample_id", "snp_code", "Exposure_ID", "Exposure_Label", "Exposure_Type",
    "Variant_ID", "rsID", "Variant_Label", "Genomic_Label", "Label_Source", "Gene",
    "Consequence", "IMPACT", "HGVSp", "VARIANT_CLASS",
    "Component_Genes", "Component_Variants", "Component_Variant_Count",
]


def _apply_fdr(df: pd.DataFrame, p_col: str, out_col: str, group_cols: list[str]) -> pd.DataFrame:
    if df.empty or p_col not in df.columns:
        return df
    out = []
    for _, sub in df.groupby(group_cols, dropna=False):
        sub = sub.copy()
        valid = sub[p_col].notna()
        if valid.sum() > 1:
            sub.loc[valid, out_col] = false_discovery_control(sub.loc[valid, p_col].to_numpy(), method="bh")
        else:
            sub[out_col] = sub[p_col]
        out.append(sub)
    return pd.concat(out, ignore_index=True) if out else df.copy()


def _test_continuous_local(carriers: pd.Series, noncarriers: pd.Series) -> dict[str, object] | None:
    carrier_vals = pd.to_numeric(carriers, errors="coerce").dropna()
    noncarrier_vals = pd.to_numeric(noncarriers, errors="coerce").dropna()
    if carrier_vals.empty or noncarrier_vals.empty:
        return None
    try:
        stat, p_val = mannwhitneyu(carrier_vals, noncarrier_vals, alternative="two-sided")
    except ValueError:
        return None
    return {
        "Test_Unadj": "Mann-Whitney U",
        "Stat_Unadj": float(stat),
        "P_Unadj": float(p_val),
        "Median_Carriers": float(carrier_vals.median()),
        "Median_NonCarriers": float(noncarrier_vals.median()),
        "Mean_Carriers": float(carrier_vals.mean()),
        "Mean_NonCarriers": float(noncarrier_vals.mean()),
        "N_Carrier_Valid": int(carrier_vals.shape[0]),
        "N_NonCarrier_Valid": int(noncarrier_vals.shape[0]),
        "Effect_Direction": "Higher in carriers" if carrier_vals.median() > noncarrier_vals.median() else ("Lower in carriers" if carrier_vals.median() < noncarrier_vals.median() else "Similar medians"),
        "Test_Adj_Age": "",
        "P_Adj_Age": np.nan,
        "Test_Adj_BMI": "",
        "P_Adj_BMI": np.nan,
        "Test_Adj_AgeBMI": "",
        "P_Adj_AgeBMI": np.nan,
    }


def _test_nominal_local(c_df: pd.DataFrame, nc_df: pd.DataFrame, clin_col: str) -> dict[str, object] | None:
    carrier = c_df[[clin_col]].copy()
    carrier["Carrier_Group"] = "Carrier"
    noncarrier = nc_df[[clin_col]].copy()
    noncarrier["Carrier_Group"] = "NonCarrier"
    combined = pd.concat([carrier, noncarrier], ignore_index=True)
    combined[clin_col] = combined[clin_col].astype(str).str.strip()
    combined = combined[combined[clin_col] != ""].copy()
    combined = combined[~combined[clin_col].str.lower().isin({"nan", "none", "<na>"})].copy()
    if combined.empty:
        return None
    table = pd.crosstab(combined["Carrier_Group"], combined[clin_col])
    if table.shape[0] < 2 or table.shape[1] < 2:
        return None
    try:
        stat, p_val, dof, expected = chi2_contingency(table, correction=False)
    except ValueError:
        return None
    expected_min = float(np.min(expected)) if expected.size else np.nan
    return {
        "Test_Unadj": "Chi-square",
        "Stat_Unadj": float(stat),
        "P_Unadj": float(p_val),
        "Chi2_DOF": int(dof),
        "Expected_Min_Cell": expected_min,
        "Nominal_Table": " | ".join(f"{col}: carrier={int(table.loc['Carrier', col])}, noncarrier={int(table.loc['NonCarrier', col])}" for col in table.columns),
        "Test_Adj_Age": "",
        "P_Adj_Age": np.nan,
        "Test_Adj_BMI": "",
        "P_Adj_BMI": np.nan,
        "Test_Adj_AgeBMI": "",
        "P_Adj_AgeBMI": np.nan,
        "Stat_Note": "Interpret cautiously when expected cell counts are small." if np.isfinite(expected_min) and expected_min < 5 else "",
    }


def _test_binary_local(
    c_df: pd.DataFrame,
    nc_df: pd.DataFrame,
    clin_col: str,
    age_col: str = "canon__age",
    bmi_col: str = "canon__bmi",
    include_bmi_model: bool = True,
) -> dict[str, object] | None:
    carriers = c_df.copy()
    carriers["carrier"] = 1
    noncarriers = nc_df.copy()
    noncarriers["carrier"] = 0
    combined = pd.concat([carriers, noncarriers], ignore_index=True).copy()
    combined["y"] = S17._binary_yes_no(combined[clin_col])
    valid = combined.dropna(subset=["y"]).copy()
    if valid.empty or valid["y"].nunique() < 2:
        return None

    a = int(((valid["carrier"] == 1) & (valid["y"] == 1)).sum())
    b = int(((valid["carrier"] == 1) & (valid["y"] == 0)).sum())
    c = int(((valid["carrier"] == 0) & (valid["y"] == 1)).sum())
    d = int(((valid["carrier"] == 0) & (valid["y"] == 0)).sum())
    or_fisher, p_fisher = fisher_exact([[a, b], [c, d]])
    or_haldane, or_method = S17._haldane_or(a, b, c, d)

    or_unadj, lo_unadj, hi_unadj, p_unadj = S17._fit_glm_binomial(valid, "y", ["carrier"])
    if pd.isna(or_unadj):
        or_unadj = float(or_haldane)
        lo_unadj, hi_unadj = np.nan, np.nan
        p_unadj = float(p_fisher)
        unadj_method = f"Fisher exact ({or_method})"
    else:
        unadj_method = "GLM binomial"

    age_valid = valid.assign(age=pd.to_numeric(valid.get(age_col), errors="coerce"))
    or_age, lo_age, hi_age, p_age = S17._fit_glm_binomial(age_valid, "y", ["carrier", "age"]) if age_col in valid.columns else (np.nan, np.nan, np.nan, np.nan)

    if include_bmi_model and bmi_col in valid.columns:
        bmi_valid = valid.assign(bmi=pd.to_numeric(valid.get(bmi_col), errors="coerce"))
        or_bmi, lo_bmi, hi_bmi, p_bmi = S17._fit_glm_binomial(bmi_valid, "y", ["carrier", "bmi"])
    else:
        or_bmi, lo_bmi, hi_bmi, p_bmi = np.nan, np.nan, np.nan, np.nan

    if include_bmi_model and age_col in valid.columns and bmi_col in valid.columns:
        ab_valid = valid.assign(
            age=pd.to_numeric(valid.get(age_col), errors="coerce"),
            bmi=pd.to_numeric(valid.get(bmi_col), errors="coerce"),
        )
        or_ab, lo_ab, hi_ab, p_ab = S17._fit_glm_binomial(ab_valid, "y", ["carrier", "age", "bmi"])
    else:
        or_ab, lo_ab, hi_ab, p_ab = np.nan, np.nan, np.nan, np.nan

    return {
        "Test_Unadj": unadj_method,
        "P_Unadj": float(p_unadj) if pd.notna(p_unadj) else np.nan,
        "OR_Unadj": float(or_unadj) if pd.notna(or_unadj) else np.nan,
        "OR_Unadj_CI95": f"[{lo_unadj:.3g}, {hi_unadj:.3g}]" if pd.notna(lo_unadj) and pd.notna(hi_unadj) else np.nan,
        "OR_Unadj_Method": unadj_method,
        "Fisher_P_Value": float(p_fisher),
        "Fisher_OR": float(or_fisher) if np.isfinite(or_fisher) else np.nan,
        "Carrier_Events": a,
        "Carrier_NonEvents": b,
        "NonCarrier_Events": c,
        "NonCarrier_NonEvents": d,
        "Test_Adj_Age": "GLM binomial" if pd.notna(p_age) else "",
        "OR_Adj_Age": float(or_age) if pd.notna(or_age) else np.nan,
        "OR_Adj_Age_CI95": f"[{lo_age:.3g}, {hi_age:.3g}]" if pd.notna(lo_age) and pd.notna(hi_age) else np.nan,
        "P_Adj_Age": float(p_age) if pd.notna(p_age) else np.nan,
        "Test_Adj_BMI": "GLM binomial" if pd.notna(p_bmi) else "",
        "OR_Adj_BMI": float(or_bmi) if pd.notna(or_bmi) else np.nan,
        "OR_Adj_BMI_CI95": f"[{lo_bmi:.3g}, {hi_bmi:.3g}]" if pd.notna(lo_bmi) and pd.notna(hi_bmi) else np.nan,
        "P_Adj_BMI": float(p_bmi) if pd.notna(p_bmi) else np.nan,
        "Test_Adj_AgeBMI": "GLM binomial" if pd.notna(p_ab) else "",
        "OR_Adj_AgeBMI": float(or_ab) if pd.notna(or_ab) else np.nan,
        "OR_Adj_AgeBMI_CI95": f"[{lo_ab:.3g}, {hi_ab:.3g}]" if pd.notna(lo_ab) and pd.notna(hi_ab) else np.nan,
        "P_Adj_AgeBMI": float(p_ab) if pd.notna(p_ab) else np.nan,
        "N_Adj_Age": int(age_valid[["carrier", "age", "y"]].dropna().shape[0]) if age_col in valid.columns else np.nan,
        "N_Adj_BMI": int(bmi_valid[["carrier", "bmi", "y"]].dropna().shape[0]) if include_bmi_model and bmi_col in valid.columns else np.nan,
        "N_Adj_AgeBMI": int(ab_valid[["carrier", "age", "bmi", "y"]].dropna().shape[0]) if include_bmi_model and age_col in valid.columns and bmi_col in valid.columns else np.nan,
        "Stat_Note": "Perfect separation or sparse outcome cells can make adjusted ORs unstable." if min(a, b, c, d) == 0 else "",
    }

def load_master(master_path: Path, manifest_paths: Optional[Dict[str, Path]], phased_path: Path) -> pd.DataFrame:
    master = S17.load_clinical_master(master_path, manifest_paths=manifest_paths)
    master = master[~master["is_replicate"]].copy()
    master["snp_code"] = master["snp_code"].astype(str).str.strip()
    master = master.drop_duplicates("snp_code").copy()

    phased_header = pd.read_csv(phased_path, sep="	", nrows=0)
    raw_sample_names = [c for c in phased_header.columns if c not in {"CHROM", "POS", "REF", "ALT", "ID"}]
    identity_map = build_analysis_sample_map(raw_sample_names, phased_path, S17._extract_snp_code)
    expanded = identity_map.merge(master, on="snp_code", how="inner")
    expanded["Sample"] = expanded["analysis_sample_id"]
    return expanded.drop_duplicates("analysis_sample_id").copy()


def load_rare_variants(
    annot_path: Path,
    sheet_name: str,
    master: pd.DataFrame,
    rare_af_threshold: float,
    include_missing_af: bool,
) -> pd.DataFrame:
    print("=== Loading annotated rare variants ===")
    annot = pd.read_excel(annot_path, sheet_name=sheet_name, engine="openpyxl")
    annot.columns = [str(c).strip() for c in annot.columns]

    if "Sample" not in annot.columns:
        raise ValueError("Annotated workbook is missing the Sample column")
    needed = ["CHROM", "POS", "REF", "ALT"]
    missing = [col for col in needed if col not in annot.columns]
    if missing:
        raise ValueError(f"Annotated workbook is missing required columns: {missing}")

    annot["Sample"] = annot["Sample"].astype(str).str.strip()
    annot = attach_analysis_sample_ids(
        annot,
        raw_col="Sample",
        phased_path=DEFAULT_PHASED,
        extract_snp_code=S17._extract_snp_code,
    )
    annot = annot.dropna(subset=["snp_code", "analysis_sample_id"]).copy()
    annot["snp_code"] = annot["snp_code"].astype(str).str.strip()

    if "GT" in annot.columns:
        annot["GT"] = _normalise_gt(annot["GT"])
    else:
        annot["GT"] = pd.NA

    for col in ["POS", "DP", "AF", "gnomADe_NFE_AF", "gnomADg_NFE_AF"]:
        if col in annot.columns:
            annot[col] = pd.to_numeric(annot[col], errors="coerce")
        else:
            annot[col] = np.nan

    annot["Variant_Key"] = _build_variant_key(annot)
    annot["Gene_Label"] = (
        annot.get("SYMBOL", pd.Series("", index=annot.index))
        .fillna(annot.get("Gene", pd.Series("", index=annot.index)))
        .astype(str)
        .str.strip()
        .replace({"nan": "", "None": ""})
    )
    display_lookup = shared_build_variant_display_table(
        annot,
        variant_col="Variant_Key",
        existing_col="Existing_variation",
        symbol_cols=["SYMBOL", "Gene"],
        hgvs_col="HGVSp" if "HGVSp" in annot.columns else ("HGVSc" if "HGVSc" in annot.columns else None),
        chrom_col="CHROM",
        pos_col="POS",
        ref_col="REF",
        alt_col="ALT",
        output_id_col="Variant_Key",
    )
    annot = annot.merge(display_lookup, on="Variant_Key", how="left", suffixes=("", "__display"))
    for col in ["rsID", "Display_Label", "Variant_Label", "Genomic_Label", "Label_Source", "HGVS_Short"]:
        display_col = f"{col}__display"
        if display_col not in annot.columns:
            continue
        if col in annot.columns:
            annot[col] = annot[col].replace({"": pd.NA}).fillna(annot[display_col])
            annot = annot.drop(columns=[display_col], errors="ignore")
        else:
            annot = annot.rename(columns={display_col: col})
    annot["Variant_ID"] = annot["Variant_Label"].fillna(annot["Variant_Key"])
    annot["rsID"] = annot["rsID"].replace({"": pd.NA})

    annot = combine_gnomad_nfe(annot)
    common_mask = common_nfe_variant_mask(annot, rare_af_threshold)
    rare_mask = ~common_mask
    if not include_missing_af:
        af_missing_mask = annot[["gnomADe_NFE_AF", "gnomADg_NFE_AF", "gnomAD_NFE_AF_combined"]].isna().all(axis=1)
        rare_mask = rare_mask & ~af_missing_mask
    annot = annot[rare_mask].copy()

    dedup_cols = ["analysis_sample_id", "Variant_Key"]
    annot["_missing_score"] = annot[["HGVSp", "Consequence", "IMPACT", "Existing_variation"]].isna().sum(axis=1)
    annot = annot.sort_values(["snp_code", "Variant_Key", "_missing_score"]).drop_duplicates(dedup_cols, keep="first")
    annot = annot.drop(columns=["_missing_score"], errors="ignore")

    merged = annot.merge(master, on=["analysis_sample_id", "snp_code"], how="inner", suffixes=("", "_master"))
    merged["Sample"] = merged["analysis_sample_id"]
    print(
        f"  Rare rows retained: {len(merged)} | unique variants: {merged['Variant_Key'].nunique()} | "
        f"matched analysis samples: {merged['analysis_sample_id'].nunique()}"
    )
    return merged


def build_variant_catalogue(rare_df: pd.DataFrame, sample_manifest: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "Variant_Key", "Variant_Label", "Variant_ID", "rsID", "Genomic_Label", "Label_Source", "Gene", "CHROM", "POS", "REF", "ALT",
        "VARIANT_CLASS", "Consequence", "IMPACT", "HGVSp", "Existing_variation",
        "gnomAD_NFE_AF_combined", "gnomAD_NFE_Source", "N_Carriers_Total",
        "N_Carriers_Tumour", "N_Carriers_Control", "Carrier_Cohorts", "Carrier_Tissues", "Carrier_Samples",
    ]
    if rare_df.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    control_codes = set(sample_manifest.loc[sample_manifest["Tissue"] == "Healthy", "analysis_sample_id"])
    tumour_codes = set(sample_manifest.loc[sample_manifest["Tissue"] == "Tumour", "analysis_sample_id"])
    for variant_key, sub in rare_df.groupby("Variant_Key"):
        carriers = set(sub["analysis_sample_id"].astype(str))
        rows.append({
            "Variant_Key": variant_key,
            "Variant_Label": _first_nonempty(sub["Variant_Label"], variant_key),
            "Variant_ID": _first_nonempty(sub["Variant_ID"], variant_key),
            "rsID": _first_nonempty(sub.get("rsID", pd.Series([], dtype=object))),
            "Genomic_Label": _first_nonempty(sub.get("Genomic_Label", pd.Series([], dtype=object)), variant_key),
            "Label_Source": _first_nonempty(sub.get("Label_Source", pd.Series([], dtype=object))),
            "Gene": _first_nonempty(sub["Gene_Label"]),
            "CHROM": _first_nonempty(sub["CHROM"]),
            "POS": pd.to_numeric(sub["POS"], errors="coerce").dropna().astype(int).min() if sub["POS"].notna().any() else np.nan,
            "REF": _first_nonempty(sub["REF"]),
            "ALT": _first_nonempty(sub["ALT"]),
            "VARIANT_CLASS": _first_nonempty(sub.get("VARIANT_CLASS", pd.Series([], dtype=object))),
            "Consequence": _first_nonempty(sub.get("Consequence", pd.Series([], dtype=object))),
            "IMPACT": _first_nonempty(sub.get("IMPACT", pd.Series([], dtype=object))),
            "HGVSp": _first_nonempty(sub.get("HGVSp", pd.Series([], dtype=object))),
            "Existing_variation": _first_nonempty(sub.get("Existing_variation", pd.Series([], dtype=object))),
            "gnomAD_NFE_AF_combined": pd.to_numeric(sub["gnomAD_NFE_AF_combined"], errors="coerce").dropna().min() if sub["gnomAD_NFE_AF_combined"].notna().any() else np.nan,
            "gnomAD_NFE_Source": _first_nonempty(sub.get("gnomAD_NFE_Source", pd.Series([], dtype=object))),
            "N_Carriers_Total": len(carriers),
            "N_Carriers_Tumour": len(carriers & tumour_codes),
            "N_Carriers_Control": len(carriers & control_codes),
            "Carrier_Cohorts": _safe_join(sub["Cohort"]),
            "Carrier_Tissues": _safe_join(sub["Tissue"]),
            "Carrier_Samples": _safe_join(sorted(carriers)),
        })
    return pd.DataFrame(rows, columns=columns).sort_values(["Gene", "Variant_Label", "Variant_Key"]).reset_index(drop=True)


def build_sample_burden(rare_df: pd.DataFrame, sample_manifest: pd.DataFrame) -> pd.DataFrame:
    burden = (
        rare_df.groupby("analysis_sample_id", as_index=False)
        .agg(
            N_Rare_Variants=("Variant_Key", pd.Series.nunique),
            N_Rare_Genes=("Gene_Label", lambda s: len({v for v in s if str(v).strip()})),
            Rare_Variant_List=("Variant_Label", _safe_join),
            Rare_Genes=("Gene_Label", _safe_join),
        )
    )
    out = sample_manifest.merge(burden, on="analysis_sample_id", how="left")
    out["N_Rare_Variants"] = out["N_Rare_Variants"].fillna(0).astype(int)
    out["Rare_Genes"] = out["Rare_Genes"].fillna("")
    out["Rare_Variant_List"] = out["Rare_Variant_List"].fillna("")
    out["Has_Any_Rare_Variant"] = np.where(out["N_Rare_Variants"] >= 1, "YES", "NO")
    out["Has_Multiple_Rare_Variants"] = np.where(out["N_Rare_Variants"] >= 2, "YES", "NO")
    keep_cols = [
        col for col in [
            "analysis_sample_id", "snp_code", "sample_id", "sheet", "Cohort", "Tissue", "canon__age", "canon__bmi",
            "N_Rare_Variants", "N_Rare_Genes", "Has_Any_Rare_Variant", "Has_Multiple_Rare_Variants",
            "Rare_Genes", "Rare_Variant_List",
        ] if col in out.columns
    ]
    return out[keep_cols].sort_values(["Cohort", "Tissue", "analysis_sample_id"]).reset_index(drop=True)


def build_single_variant_exposures(rare_df: pd.DataFrame) -> pd.DataFrame:
    cols = ["analysis_sample_id", "snp_code", "Variant_Key", "Variant_ID", "rsID", "Variant_Label", "Genomic_Label", "Label_Source", "Gene_Label", "Consequence", "IMPACT", "HGVSp", "VARIANT_CLASS"]
    exp = rare_df[cols].drop_duplicates().copy()
    exp["Exposure_ID"] = exp["Variant_Key"]
    exp["Exposure_Label"] = exp["Variant_Label"]
    exp["Gene"] = exp["Gene_Label"]
    exp["Exposure_Type"] = "Rare_Variant"
    exp["Component_Genes"] = exp["Gene_Label"]
    exp["Component_Variants"] = exp["Variant_Label"]
    exp["Component_Variant_Count"] = 1
    return exp[EXPOSURE_COLUMNS]


def build_gene_burden_exposures(rare_df: pd.DataFrame) -> pd.DataFrame:
    sub = rare_df[rare_df["Gene_Label"].astype(str).str.strip() != ""].copy()
    exp = (
        sub.groupby(["analysis_sample_id", "snp_code", "Gene_Label"], as_index=False)
        .agg(
            Component_Genes=("Gene_Label", _safe_join),
            Component_Variants=("Variant_Label", _safe_join),
            Component_Variant_Count=("Variant_Key", pd.Series.nunique),
        )
        .rename(columns={"Gene_Label": "Gene"})
    )
    exp["Exposure_ID"] = "GeneBurden:" + exp["Gene"].astype(str)
    exp["Exposure_Label"] = exp["Gene"].astype(str) + " rare burden"
    exp["Exposure_Type"] = "Gene_Burden"
    exp["Variant_ID"] = ""
    exp["rsID"] = pd.NA
    exp["Variant_Label"] = exp["Exposure_Label"]
    exp["Genomic_Label"] = ""
    exp["Label_Source"] = "Gene_Burden"
    exp["Consequence"] = "Any rare variant in gene"
    exp["IMPACT"] = "Mixed"
    exp["HGVSp"] = ""
    exp["VARIANT_CLASS"] = "Burden"
    return exp[EXPOSURE_COLUMNS]


def build_burden_flag_exposures(sample_burden: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in sample_burden.iterrows():
        component_genes = row.get("Rare_Genes", "")
        component_variants = row.get("Rare_Variant_List", "")
        component_count = row.get("N_Rare_Variants", 0)
        if row["N_Rare_Variants"] >= 1:
            rows.append({"analysis_sample_id": row["analysis_sample_id"], "snp_code": row["snp_code"], "Exposure_ID": "Any_Rare_Variant", "Exposure_Label": "Any rare variant", "Exposure_Type": "Burden_Flag", "Variant_ID": "", "rsID": pd.NA, "Variant_Label": "Any rare variant", "Genomic_Label": "", "Label_Source": "Burden_Flag", "Gene": component_genes, "Consequence": "Any rare variant in GSDMB-region panel", "IMPACT": "Mixed", "HGVSp": "", "VARIANT_CLASS": "Burden", "Component_Genes": component_genes, "Component_Variants": component_variants, "Component_Variant_Count": component_count})
        if row["N_Rare_Variants"] >= 2:
            rows.append({"analysis_sample_id": row["analysis_sample_id"], "snp_code": row["snp_code"], "Exposure_ID": "Multiple_Rare_Variants", "Exposure_Label": "Multiple rare variants", "Exposure_Type": "Burden_Flag", "Variant_ID": "", "rsID": pd.NA, "Variant_Label": "Multiple rare variants", "Genomic_Label": "", "Label_Source": "Burden_Flag", "Gene": component_genes, "Consequence": "At least two rare variants in GSDMB-region panel", "IMPACT": "Mixed", "HGVSp": "", "VARIANT_CLASS": "Burden", "Component_Genes": component_genes, "Component_Variants": component_variants, "Component_Variant_Count": component_count})
    return pd.DataFrame(rows, columns=EXPOSURE_COLUMNS)

def build_pair_exposures(single_variant_exp: pd.DataFrame, min_pair_carriers: int) -> pd.DataFrame:
    label_lookup = (
        single_variant_exp.drop_duplicates("Exposure_ID")
        .set_index("Exposure_ID")["Exposure_Label"]
        .to_dict()
    )
    gene_lookup = (
        single_variant_exp.drop_duplicates("Exposure_ID")
        .set_index("Exposure_ID")["Gene"]
        .to_dict()
    )
    per_sample = single_variant_exp.groupby("analysis_sample_id")["Exposure_ID"].apply(lambda s: sorted(set(s))).to_dict()
    sample_code_lookup = single_variant_exp.groupby("analysis_sample_id")["snp_code"].agg(_first_nonempty).to_dict()
    pair_rows = []
    for analysis_sample_id, variant_keys in per_sample.items():
        if len(variant_keys) < 2:
            continue
        for left, right in combinations(variant_keys, 2):
            pair_rows.append({"analysis_sample_id": analysis_sample_id, "snp_code": sample_code_lookup.get(analysis_sample_id, ""), "Exposure_ID": f"Pair:{left}||{right}", "Left_Variant": left, "Right_Variant": right})
    if not pair_rows:
        return pd.DataFrame(columns=EXPOSURE_COLUMNS)
    pair_df = pd.DataFrame(pair_rows)
    carrier_counts = pair_df.groupby("Exposure_ID")["analysis_sample_id"].nunique()
    keep_ids = carrier_counts[carrier_counts >= min_pair_carriers].index
    pair_df = pair_df[pair_df["Exposure_ID"].isin(keep_ids)].copy()
    if pair_df.empty:
        return pd.DataFrame(columns=EXPOSURE_COLUMNS)
    pair_df["Exposure_Label"] = pair_df.apply(lambda r: _pair_label(r["Left_Variant"], r["Right_Variant"], label_lookup), axis=1)
    pair_df["Exposure_Type"] = "Variant_Pair"
    pair_df["Variant_ID"] = ""
    pair_df["rsID"] = pd.NA
    pair_df["Variant_Label"] = pair_df["Exposure_Label"]
    pair_df["Genomic_Label"] = ""
    pair_df["Label_Source"] = "Variant_Pair"
    pair_df["Gene"] = "GSDMB"
    pair_df["Consequence"] = "Recurrent rare-variant pair"
    pair_df["IMPACT"] = "Mixed"
    pair_df["HGVSp"] = ""
    pair_df["VARIANT_CLASS"] = "Pair"
    pair_df["Component_Genes"] = pair_df.apply(lambda r: _safe_join([gene_lookup.get(r["Left_Variant"], ""), gene_lookup.get(r["Right_Variant"], "")]), axis=1)
    pair_df["Component_Variants"] = pair_df["Exposure_Label"]
    pair_df["Component_Variant_Count"] = 2
    return pair_df[EXPOSURE_COLUMNS].drop_duplicates()


def build_exposure_catalogue(exposure_df: pd.DataFrame, sample_manifest: pd.DataFrame) -> pd.DataFrame:
    if exposure_df.empty:
        return pd.DataFrame()
    controls = set(sample_manifest.loc[sample_manifest["Tissue"] == "Healthy", "analysis_sample_id"])
    tumours = set(sample_manifest.loc[sample_manifest["Tissue"] == "Tumour", "analysis_sample_id"])
    rows = []
    for exp_id, sub in exposure_df.groupby("Exposure_ID"):
        carriers = set(sub["analysis_sample_id"])
        meta_defaults = _exposure_meta_defaults(sub)
        rows.append({
            "Exposure_ID": exp_id,
            "Exposure_Label": _first_nonempty(sub["Exposure_Label"], exp_id),
            "Exposure_Type": _first_nonempty(sub["Exposure_Type"]),
            "Variant_ID": meta_defaults["Variant_ID"],
            "rsID": meta_defaults["rsID"],
            "Variant_Label": meta_defaults["Variant_Label"],
            "Genomic_Label": meta_defaults["Genomic_Label"],
            "Label_Source": meta_defaults["Label_Source"],
            "Exposure_rsIDs": _safe_join(sub.get("rsID", pd.Series([], dtype=object))),
            "Gene": _first_nonempty(sub["Gene"]),
            "Consequence": _first_nonempty(sub["Consequence"]),
            "IMPACT": _first_nonempty(sub["IMPACT"]),
            "Component_Genes": _join_component_lists(sub.get("Component_Genes", pd.Series([], dtype=object))),
            "Component_Variants": _join_component_lists(sub.get("Component_Variants", pd.Series([], dtype=object))),
            "Component_Variant_Count": _component_count(sub.get("Component_Variants", pd.Series([], dtype=object))),
            "N_Carriers_Total": len(carriers),
            "N_Carriers_Tumour": len(carriers & tumours),
            "N_Carriers_Control": len(carriers & controls),
            "Carrier_Samples": _safe_join(sorted(carriers)),
        })
    return pd.DataFrame(rows).sort_values(["Exposure_Type", "N_Carriers_Total", "Exposure_Label"], ascending=[True, False, True]).reset_index(drop=True)


def tumour_vs_control(exposure_df: pd.DataFrame, sample_manifest: pd.DataFrame, fdr_threshold: float) -> pd.DataFrame:
    if exposure_df.empty:
        return pd.DataFrame()
    sample_manifest = sample_manifest.drop_duplicates("analysis_sample_id").copy()
    pooled_controls = set(sample_manifest.loc[sample_manifest["Tissue"] == "Healthy", "analysis_sample_id"])
    comparisons = {
        "Breast": set(sample_manifest.loc[(sample_manifest["Tissue"] == "Tumour") & (sample_manifest["Cohort"].astype(str).str.contains("Breast", case=False, na=False)), "analysis_sample_id"]),
        "Endometrium": set(sample_manifest.loc[(sample_manifest["Tissue"] == "Tumour") & (sample_manifest["Cohort"].astype(str).str.contains("Endometri", case=False, na=False)), "analysis_sample_id"]),
        "Global": set(sample_manifest.loc[sample_manifest["Tissue"] == "Tumour", "analysis_sample_id"]),
    }
    carrier_lookup = exposure_df.groupby("Exposure_ID")["analysis_sample_id"].agg(lambda s: set(s)).to_dict()
    meta = exposure_df.drop_duplicates("Exposure_ID").set_index("Exposure_ID")
    rows = []
    for analysis_group, tumour_codes in comparisons.items():
        n_tumour = len(tumour_codes)
        n_control = len(pooled_controls)
        if n_tumour < MIN_COMPARISON_CARRIERS or n_control < MIN_COMPARISON_CARRIERS:
            continue
        for exposure_id, carriers in carrier_lookup.items():
            a = len(carriers & tumour_codes)
            c = len(carriers & pooled_controls)
            if a + c < MIN_COMPARISON_CARRIERS:
                continue
            b = n_tumour - a
            d = n_control - c
            _, p_val = fisher_exact([[a, b], [c, d]])
            odds_ratio, or_method = S17._haldane_or(a, b, c, d)
            meta_row = meta.loc[exposure_id]
            meta_defaults = _exposure_meta_defaults(pd.DataFrame([meta_row]))
            comparison_carriers = (carriers & tumour_codes) | (carriers & pooled_controls)
            component_sub = exposure_df[
                (exposure_df["Exposure_ID"] == exposure_id) &
                (exposure_df["analysis_sample_id"].isin(comparison_carriers))
            ]
            rows.append({
                "Analysis_Group": analysis_group,
                "Exposure_ID": exposure_id,
                "Exposure_Label": meta_row["Exposure_Label"],
                "Exposure_Type": meta_row["Exposure_Type"],
                "Variant_ID": meta_defaults["Variant_ID"],
                "rsID": meta_defaults["rsID"],
                "Variant_Label": meta_defaults["Variant_Label"],
                "Genomic_Label": meta_defaults["Genomic_Label"],
                "Label_Source": meta_defaults["Label_Source"],
                "Gene": meta_row.get("Gene", ""),
                "Consequence": meta_row.get("Consequence", ""),
                "IMPACT": meta_row.get("IMPACT", ""),
                "Component_Genes": _join_component_lists(component_sub.get("Component_Genes", pd.Series([], dtype=object))),
                "Component_Variants": _join_component_lists(component_sub.get("Component_Variants", pd.Series([], dtype=object))),
                "Component_Variant_Count": _component_count(component_sub.get("Component_Variants", pd.Series([], dtype=object))),
                "N_Tumour": n_tumour,
                "N_Control": n_control,
                "Carriers_Tumour": a,
                "Carriers_Control": c,
                "Freq_Tumour_%": round(a / n_tumour * 100, 2),
                "Freq_Control_%": round(c / n_control * 100, 2),
                "Odds_Ratio": round(float(odds_ratio), 4),
                "OR_Method": or_method,
                "P_Value": p_val,
            })
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    res = _apply_fdr(res, "P_Value", "FDR_P_Value", ["Analysis_Group", "Exposure_Type"])
    res["Nominal_Sig"] = res["P_Value"] < 0.05
    res["FDR_Sig"] = res["FDR_P_Value"] < fdr_threshold
    return res.sort_values(["Exposure_Type", "Analysis_Group", "P_Value", "Exposure_Label"]).reset_index(drop=True)


def clinical_associations(
    exposure_df: pd.DataFrame,
    sample_manifest: pd.DataFrame,
    clinical_vars: Dict[str, Dict],
    cohort_filter: str,
    cohort_label: str,
    min_carriers: int,
    fdr_threshold: float,
) -> pd.DataFrame:
    if exposure_df.empty:
        return pd.DataFrame()
    cohort_df = sample_manifest[(sample_manifest["Tissue"] == "Tumour") & (sample_manifest["Cohort"].astype(str).str.contains(cohort_filter, case=False, na=False))].copy()
    if cohort_df.empty:
        return pd.DataFrame()
    carrier_lookup = exposure_df.groupby("Exposure_ID")["analysis_sample_id"].agg(lambda s: set(s)).to_dict()
    meta_lookup = exposure_df.drop_duplicates("Exposure_ID").set_index("Exposure_ID")
    rows = []
    for exposure_id, carriers in carrier_lookup.items():
        carrier_mask = cohort_df["analysis_sample_id"].isin(carriers)
        c_df = cohort_df[carrier_mask].copy()
        nc_df = cohort_df[~carrier_mask].copy()
        if len(c_df) < min_carriers or len(nc_df) < min_carriers:
            continue
        meta_row = meta_lookup.loc[exposure_id]
        meta_defaults = _exposure_meta_defaults(pd.DataFrame([meta_row]))
        component_sub = exposure_df[
            (exposure_df["Exposure_ID"] == exposure_id) &
            (exposure_df["analysis_sample_id"].isin(c_df["analysis_sample_id"]))
        ]
        for clin_col, clin_meta in clinical_vars.items():
            if clin_col not in cohort_df.columns:
                continue
            if clin_meta["type"] == "continuous":
                res = _test_continuous_local(c_df[clin_col], nc_df[clin_col])
            elif clin_meta["type"] == "binary":
                res = _test_binary_local(
                    c_df,
                    nc_df,
                    clin_col,
                    age_col="canon__age",
                    bmi_col="canon__bmi",
                    include_bmi_model=bool(clin_meta.get("bmi_adjust", True)),
                )
            elif clin_meta["type"] == "nominal":
                res = _test_nominal_local(c_df, nc_df, clin_col)
            else:
                res = None
            if res is None:
                continue
            row = {
                "Cohort": cohort_label,
                "Exposure_ID": exposure_id,
                "Exposure_Label": meta_row["Exposure_Label"],
                "Exposure_Type": meta_row["Exposure_Type"],
                "Variant_ID": meta_defaults["Variant_ID"],
                "rsID": meta_defaults["rsID"],
                "Variant_Label": meta_defaults["Variant_Label"],
                "Genomic_Label": meta_defaults["Genomic_Label"],
                "Label_Source": meta_defaults["Label_Source"],
                "Gene": meta_row.get("Gene", ""),
                "Consequence": meta_row.get("Consequence", ""),
                "IMPACT": meta_row.get("IMPACT", ""),
                "Component_Genes": _join_component_lists(component_sub.get("Component_Genes", pd.Series([], dtype=object))),
                "Component_Variants": _join_component_lists(component_sub.get("Component_Variants", pd.Series([], dtype=object))),
                "Component_Variant_Count": _component_count(component_sub.get("Component_Variants", pd.Series([], dtype=object))),
                "Clinical_Var": clin_col,
                "Clin_Label": clin_meta["label"],
                "Clin_Type": clin_meta["type"],
                "Note": clin_meta.get("note", ""),
                "N_Cohort_Samples": len(cohort_df),
                "N_Carriers_Cohort": len(c_df),
                "N_NonCarriers_Cohort": len(nc_df),
                "Carrier_Prevalence_%": round(len(c_df) / len(cohort_df) * 100, 2),
            }
            row.update(res)
            rows.append(row)
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    for p_col, out_col in [("P_Unadj", "FDR_Unadj"), ("P_Adj_Age", "FDR_Adj_Age"), ("P_Adj_BMI", "FDR_Adj_BMI"), ("P_Adj_AgeBMI", "FDR_Adj_AgeBMI")]:
        if p_col in res.columns:
            res = _apply_fdr(res, p_col, out_col, ["Clinical_Var", "Exposure_Type"])
    res["Nominal_Sig_Unadj"] = res.get("P_Unadj", pd.Series(np.nan, index=res.index)).lt(0.05)
    res["FDR_Sig_Unadj"] = res.get("FDR_Unadj", pd.Series(np.nan, index=res.index)).lt(fdr_threshold)
    res["Nominal_Sig_Adj_Age"] = res.get("P_Adj_Age", pd.Series(np.nan, index=res.index)).lt(0.05)
    res["FDR_Sig_Adj_Age"] = res.get("FDR_Adj_Age", pd.Series(np.nan, index=res.index)).lt(fdr_threshold)
    res["Nominal_Sig_Adj_BMI"] = res.get("P_Adj_BMI", pd.Series(np.nan, index=res.index)).lt(0.05)
    res["FDR_Sig_Adj_BMI"] = res.get("FDR_Adj_BMI", pd.Series(np.nan, index=res.index)).lt(fdr_threshold)
    res["Nominal_Sig_Adj_AgeBMI"] = res.get("P_Adj_AgeBMI", pd.Series(np.nan, index=res.index)).lt(0.05)
    res["FDR_Sig_Adj_AgeBMI"] = res.get("FDR_Adj_AgeBMI", pd.Series(np.nan, index=res.index)).lt(fdr_threshold)
    return res.sort_values(["Exposure_Type", "Clinical_Var", "P_Unadj", "Exposure_Label"]).reset_index(drop=True)


def build_summary_significant(*tables: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for df in tables:
        if df is None or df.empty:
            continue
        if "P_Value" in df.columns:
            sub = df[df["Nominal_Sig"] | df["FDR_Sig"]].copy()
            if not sub.empty:
                sub["Analysis"] = "Tumour_vs_Control"
                sub["P_Primary"] = sub["P_Value"]
                sub["FDR_Primary"] = sub["FDR_P_Value"]
                parts.append(sub)
        elif "P_Unadj" in df.columns:
            sig_mask = (
                df["Nominal_Sig_Unadj"] | df["FDR_Sig_Unadj"] |
                df["Nominal_Sig_Adj_Age"] | df["FDR_Sig_Adj_Age"] |
                df["Nominal_Sig_Adj_BMI"] | df["FDR_Sig_Adj_BMI"] |
                df["Nominal_Sig_Adj_AgeBMI"] | df["FDR_Sig_Adj_AgeBMI"]
            )
            sub = df[sig_mask].copy()
            if not sub.empty:
                sub["Analysis"] = "Clinical_Association"
                sub["P_Primary"] = sub["P_Unadj"]
                sub["FDR_Primary"] = sub.get("FDR_Unadj", np.nan)
                parts.append(sub)
    if not parts:
        return pd.DataFrame({"Note": ["No nominally significant rare-variant or mutation associations were detected."]})
    return pd.concat(parts, ignore_index=True).sort_values(["P_Primary", "Analysis", "Exposure_Type", "Exposure_Label"]).reset_index(drop=True)

def build_readme(
    rare_af_threshold: float,
    include_missing_af: bool,
    min_carriers: int,
    min_pair_carriers: int,
    sample_manifest: pd.DataFrame,
    rare_df: pd.DataFrame,
    exposure_catalogue: pd.DataFrame,
) -> pd.DataFrame:
    text = [
        "17b rare-variant / mutation association workbook",
        "",
        "Purpose:",
        "This workbook complements the common-SNP analysis by focusing on mutations and other rare variants in GSDMB.",
        "",
        "Rare-variant definition:",
        f"Combined gnomAD NFE AF < {rare_af_threshold:.2%}.",
        "Variants with missing gnomAD NFE AF are " + ("included" if include_missing_af else "excluded") + ".",
        "",
        "What is tested:",
        "1. Single rare variants / mutations.",
        "2. Gene-level rare burden (any rare variant in the gene).",
        "3. Recurrent rare-variant pairs shared by multiple samples.",
        "4. Simple burden flags (any rare variant, multiple rare variants).",
        "5. Exploratory clinical association testing for recurrent rare-variant pairs.",
        "Burden and pair rows include Component_Genes, Component_Variants, and Component_Variant_Count so the exact burden composition is visible downstream.",
        "",
        "Interpretation notes:",
        "Tumour-vs-control tables compare tumour samples against pooled healthy controls from the sequenced manifest.",
        "Clinical association tables use tumour samples only and compare carriers versus non-carriers.",
        "The non-carrier group is defined from all sequenced samples in the harmonised master, not only from samples with called variants.",
        f"Clinical tests require at least {min_carriers} carriers and {min_carriers} non-carriers.",
        f"Rare-variant pairs require at least {min_pair_carriers} carrier samples to enter the catalogue.",
        "Pair-clinical association sheets are exploratory and should be interpreted cautiously because pair carrier counts are typically small.",
        "",
        "Quick counts:",
        f"Sequenced non-replicate analysis samples: {sample_manifest['analysis_sample_id'].nunique()}",
        f"Rare-variant rows after transcript collapsing: {len(rare_df)}",
        f"Unique rare variants: {rare_df['Variant_Key'].nunique() if not rare_df.empty else 0}",
        f"Samples carrying at least one rare variant: {rare_df['analysis_sample_id'].nunique() if not rare_df.empty else 0}",
        f"Distinct exposures tested: {len(exposure_catalogue)}",
    ]
    return pd.DataFrame({"README": text})


def main() -> None:
    args = parse_args()
    annot_path = Path(args.annot)
    master_path = Path(args.master)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_xlsx = out_dir / "GSDMB_Rare_Variant_Association_Results.xlsx"

    validate_file_exists(annot_path, "Script 17b annotated workbook")
    validate_file_exists(master_path, "Script 17b master workbook")

    if args.no_manifests:
        manifest_paths = None
    else:
        manifest_paths = {
            "endometrium-tumour": Path(args.manifest_endo_tumour),
            "endometrium-normal": Path(args.manifest_endo_normal),
            "breast-tumour": Path(args.manifest_breast_tumour),
            "breast-normal": Path(args.manifest_breast_normal),
        }

    sample_manifest = load_master(master_path, manifest_paths, Path(args.phased))
    rare_df = load_rare_variants(
        annot_path,
        args.sheet,
        sample_manifest,
        rare_af_threshold=args.rare_af_threshold,
        include_missing_af=not args.exclude_missing_af,
    )
    print_validation_summary(rare_df, "Sample", "Script 17b rare merged input", ["Cohort", "Tissue"])

    variant_catalogue = build_variant_catalogue(rare_df, sample_manifest)
    sample_burden = build_sample_burden(rare_df, sample_manifest)
    single_variant_exp = build_single_variant_exposures(rare_df)
    gene_burden_exp = build_gene_burden_exposures(rare_df)
    burden_flag_exp = build_burden_flag_exposures(sample_burden)
    pair_exp = build_pair_exposures(single_variant_exp, args.min_pair_carriers)

    exposure_df = pd.concat([single_variant_exp, gene_burden_exp, burden_flag_exp, pair_exp], ignore_index=True)
    exposure_catalogue = build_exposure_catalogue(exposure_df, sample_manifest)

    tvc_all = tumour_vs_control(exposure_df, sample_manifest, args.fdr_threshold)
    tvc_single = tvc_all[tvc_all["Exposure_Type"] == "Rare_Variant"].copy() if not tvc_all.empty else pd.DataFrame()
    tvc_gene = tvc_all[tvc_all["Exposure_Type"] == "Gene_Burden"].copy() if not tvc_all.empty else pd.DataFrame()
    tvc_pair = tvc_all[tvc_all["Exposure_Type"] == "Variant_Pair"].copy() if not tvc_all.empty else pd.DataFrame()
    tvc_flags = tvc_all[tvc_all["Exposure_Type"] == "Burden_Flag"].copy() if not tvc_all.empty else pd.DataFrame()

    clinical_exp = exposure_df[exposure_df["Exposure_Type"].isin(["Rare_Variant", "Gene_Burden", "Burden_Flag"])].copy()
    breast_assoc = clinical_associations(
        clinical_exp,
        sample_manifest,
        S17.CLINICAL_VARS_BREAST,
        cohort_filter="Breast",
        cohort_label="Breast_Tumour",
        min_carriers=args.min_carriers,
        fdr_threshold=args.fdr_threshold,
    )
    endo_assoc = clinical_associations(
        clinical_exp,
        sample_manifest,
        S17.CLINICAL_VARS_ENDO,
        cohort_filter="Endometri",
        cohort_label="Endometrial_Tumour",
        min_carriers=args.min_carriers,
        fdr_threshold=args.fdr_threshold,
    )
    breast_pair_assoc = clinical_associations(
        pair_exp,
        sample_manifest,
        S17.CLINICAL_VARS_BREAST,
        cohort_filter="Breast",
        cohort_label="Breast_Tumour",
        min_carriers=args.min_carriers,
        fdr_threshold=args.fdr_threshold,
    )
    endo_pair_assoc = clinical_associations(
        pair_exp,
        sample_manifest,
        S17.CLINICAL_VARS_ENDO,
        cohort_filter="Endometri",
        cohort_label="Endometrial_Tumour",
        min_carriers=args.min_carriers,
        fdr_threshold=args.fdr_threshold,
    )
    if not breast_pair_assoc.empty:
        breast_pair_assoc["Exploratory_Note"] = "Exploratory pair-level clinical association; interpret cautiously due to small carrier counts."
    if not endo_pair_assoc.empty:
        endo_pair_assoc["Exploratory_Note"] = "Exploratory pair-level clinical association; interpret cautiously due to small carrier counts."
    summary = build_summary_significant(tvc_all, breast_assoc, endo_assoc, breast_pair_assoc, endo_pair_assoc)
    readme = build_readme(
        rare_af_threshold=args.rare_af_threshold,
        include_missing_af=not args.exclude_missing_af,
        min_carriers=args.min_carriers,
        min_pair_carriers=args.min_pair_carriers,
        sample_manifest=sample_manifest,
        rare_df=rare_df,
        exposure_catalogue=exposure_catalogue,
    )

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        variant_catalogue.to_excel(writer, sheet_name="rare_variant_catalogue", index=False)
        sample_burden.to_excel(writer, sheet_name="sample_burden", index=False)
        exposure_catalogue.to_excel(writer, sheet_name="exposure_catalogue", index=False)
        if not tvc_single.empty:
            tvc_single.to_excel(writer, sheet_name="tvc_single_variants", index=False)
        if not tvc_gene.empty:
            tvc_gene.to_excel(writer, sheet_name="tvc_gene_burdens", index=False)
        if not tvc_pair.empty:
            tvc_pair.to_excel(writer, sheet_name="tvc_variant_pairs", index=False)
        if not tvc_flags.empty:
            tvc_flags.to_excel(writer, sheet_name="tvc_burden_flags", index=False)
        if not breast_assoc.empty:
            breast_assoc.to_excel(writer, sheet_name="breast_mutation_assoc", index=False)
        if not endo_assoc.empty:
            endo_assoc.to_excel(writer, sheet_name="endo_mutation_assoc", index=False)
        if not breast_pair_assoc.empty:
            breast_pair_assoc.to_excel(writer, sheet_name="breast_pair_assoc", index=False)
        if not endo_pair_assoc.empty:
            endo_pair_assoc.to_excel(writer, sheet_name="endo_pair_assoc", index=False)
        summary.to_excel(writer, sheet_name="summary_significant", index=False)
        manifest_export = sample_manifest[[
            col for col in [
                "analysis_sample_id", "snp_code", "sample_id", "sheet", "Cohort", "Tissue", "canon__age", "canon__bmi", "is_replicate"
            ] if col in sample_manifest.columns
        ]].sort_values(["Cohort", "Tissue", "snp_code"])
        manifest_export.to_excel(writer, sheet_name="sample_manifest", index=False)

    variant_catalogue.to_csv(out_dir / "17b_rare_variant_catalogue.tsv", sep="\t", index=False)
    exposure_catalogue.to_csv(out_dir / "17b_exposure_catalogue.tsv", sep="\t", index=False)
    if not tvc_all.empty:
        tvc_all.to_csv(out_dir / "17b_tumour_vs_control.tsv", sep="\t", index=False)
    if not breast_assoc.empty:
        breast_assoc.to_csv(out_dir / "17b_breast_mutation_assoc.tsv", sep="\t", index=False)
    if not endo_assoc.empty:
        endo_assoc.to_csv(out_dir / "17b_endo_mutation_assoc.tsv", sep="\t", index=False)

    print("\n" + "=" * 60)
    print("RARE-VARIANT ASSOCIATION ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"  Sequenced analysis samples: {sample_manifest['analysis_sample_id'].nunique()}")
    print(f"  Rare variants tested: {variant_catalogue['Variant_Key'].nunique() if not variant_catalogue.empty else 0}")
    print(f"  Exposure catalogue size: {len(exposure_catalogue)}")
    if not tvc_all.empty:
        print(f"  Tumour-vs-control results: {len(tvc_all)} rows | nominal={int(tvc_all['Nominal_Sig'].sum())} | FDR={int(tvc_all['FDR_Sig'].sum())}")
    if not breast_assoc.empty:
        print(f"  Breast clinical rows: {len(breast_assoc)} | nominal unadj={int(breast_assoc['Nominal_Sig_Unadj'].sum())} | FDR unadj={int(breast_assoc['FDR_Sig_Unadj'].sum())}")
    if not endo_assoc.empty:
        print(f"  Endometrial clinical rows: {len(endo_assoc)} | nominal unadj={int(endo_assoc['Nominal_Sig_Unadj'].sum())} | FDR unadj={int(endo_assoc['FDR_Sig_Unadj'].sum())}")
    if not breast_pair_assoc.empty:
        print(f"  Breast pair-clinical rows: {len(breast_pair_assoc)} | nominal unadj={int(breast_pair_assoc['Nominal_Sig_Unadj'].sum())} | FDR unadj={int(breast_pair_assoc['FDR_Sig_Unadj'].sum())}")
    if not endo_pair_assoc.empty:
        print(f"  Endometrial pair-clinical rows: {len(endo_pair_assoc)} | nominal unadj={int(endo_pair_assoc['Nominal_Sig_Unadj'].sum())} | FDR unadj={int(endo_pair_assoc['FDR_Sig_Unadj'].sum())}")
    print(f"  Workbook: {out_xlsx}")


if __name__ == "__main__":
    main()
