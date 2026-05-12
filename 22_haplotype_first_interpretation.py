#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
22_haplotype_first_interpretation.py
====================================

Build a haplotype-first interpretation layer on top of the existing
frequency/association outputs without changing the underlying tests.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from pipeline_utils import find_col, get_paths

FDR_THRESHOLD = 0.10
PRIMARY_COMPARISONS = ("Breast", "All", "Endometrium")
HAP_EFFECT_COLORS = {
    "Tumour-enriched": "#b2182b",
    "Healthy-enriched": "#2166ac",
    "Balanced": "#6c757d",
}
COMPARISON_MARKERS = {
    "Breast": "s",
    "All": "o",
    "Endometrium": "^",
}


def harmonised_defaults() -> dict[str, Path]:
    paths = get_paths()
    hap_stats = paths["haplotype_results"]
    return {
        "hap_stats": hap_stats,
        "hap_stats_female": hap_stats.with_name("GSDMB_Haplotype_Results_FemaleOnly_StudySamples.xlsx"),
        "snp_enrichment": paths["snp_enrichment_dir"] / "GSDMB_SNP_Enrichment_Results.xlsx",
        "hap_compare": paths["thousand_genomes_results"] / "19_1000G_Haplotype_Comparison.xlsx",
        "phased": paths["haplotype_phased"],
        "sample_meta": paths["haplotype_phased_dir"] / "sample_metadata.tsv",
        "out_dir": paths["results_dir"] / "22_haplotype_first_interpretation",
    }


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_excel(path: Path, sheet_name: str) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet_name)


def first_available_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    for candidate in candidates:
        col = find_col(df, candidate)
        if col is not None:
            return col
    return None


def normalise_stage12_snp_sheet(df: pd.DataFrame, sheet_name: str) -> pd.DataFrame:
    out = df.copy()
    required = first_available_column(out, ["FDR_Significant", "Primary_Genotype_FDR_Significant"])
    if required is None:
        raise ValueError(f"Stage-12 {sheet_name} sheet is missing an FDR significance column.")
    out = out[out[required].eq(True)].copy()
    if out.empty:
        return out

    schema_map = {
        "SNP_ID": ["rsID", "SNP_ID"],
        "Symbol": ["Symbol", "Gene"],
        "LD_Block_ID": ["LD_Block_ID", "LD_Block_Assignment"],
        "Odds_Ratio": ["Odds_Ratio", "Carrier_Odds_Ratio", "Het_vs_WT_Odds_Ratio", "Hom_vs_WT_Odds_Ratio"],
        "P_Value": ["P_Value", "Carrier_P_Value", "Primary_Genotype_P_Value", "Best_Genotype_P"],
        "FDR_P_Value": ["FDR_P_Value", "Carrier_FDR_P_Value", "Primary_Genotype_FDR_P_Value", "Best_Genotype_FDR"],
        "POS": ["POS"],
        "Dose_Difference": ["Dose_Difference"],
    }
    for target, candidates in schema_map.items():
        source = first_available_column(out, candidates)
        if source is not None and source != target:
            out[target] = out[source]

    missing = [name for name in ["SNP_ID", "Odds_Ratio", "P_Value", "FDR_P_Value"] if name not in out.columns]
    if missing:
        raise ValueError(
            f"Stage-12 {sheet_name} sheet is missing required columns after schema normalisation: {', '.join(missing)}"
        )

    for numeric_col in ["Odds_Ratio", "P_Value", "FDR_P_Value", "POS", "Dose_Difference"]:
        if numeric_col in out.columns:
            out[numeric_col] = pd.to_numeric(out[numeric_col], errors="coerce")
    return out


def inspect_female_subset_workbook(path: Path) -> dict[str, object]:
    info: dict[str, object] = {
        "available": False,
        "verified": False,
        "sample_filter_mode": "",
        "sample_filter_file": "",
        "analysis_label": "",
        "note": "Female-only study subset workbook not found; study-only female frequency checks will be left missing.",
    }
    if not path.exists():
        return info

    info["available"] = True
    try:
        xl = pd.ExcelFile(path)
    except Exception as exc:
        info["note"] = f"Female-only study subset workbook could not be opened: {exc}"
        return info

    if "Analysis_Metadata" not in xl.sheet_names:
        info["note"] = "Female-only study subset workbook is present but has no Analysis_Metadata sheet, so the sample filter cannot be verified."
        return info

    meta = pd.read_excel(path, sheet_name="Analysis_Metadata")
    if {"Setting", "Value"}.issubset(meta.columns):
        meta_map = dict(zip(meta["Setting"].astype(str), meta["Value"]))
    else:
        meta_map = {}
    info["sample_filter_mode"] = str(meta_map.get("Sample_Filter_Mode", "") or "")
    info["sample_filter_file"] = str(meta_map.get("Sample_Filter_File", "") or "")
    info["analysis_label"] = str(meta_map.get("Analysis_Label", "") or "")
    info["verified"] = info["sample_filter_mode"] == "ExplicitSampleList"
    if info["verified"]:
        info["note"] = "Female-only study frequency checks use the explicit study-sample subset workbook generated by stage 15R."
    else:
        info["note"] = "Female-only study subset workbook is present but not marked as an explicit sample-list rerun; female study frequencies will be treated as unverified."
    return info

def load_significant_haplotypes(hap_stats_path: Path) -> pd.DataFrame:
    xl = pd.ExcelFile(hap_stats_path)
    frames = []
    pattern = re.compile(r"^LD_Block_r2_080_Block\d+_Assoc$")
    for sheet in xl.sheet_names:
        if sheet != "Full_Associations" and not pattern.match(sheet):
            continue
        df = read_excel(hap_stats_path, sheet)
        if "Significant_FDR" not in df.columns:
            continue
        df = df[df["Significant_FDR"].eq(True) & df["Comparison"].isin(PRIMARY_COMPARISONS)].copy()
        if df.empty:
            continue
        df["Assoc_Sheet"] = sheet
        frames.append(df)
    if not frames:
        return pd.DataFrame()

    sig = pd.concat(frames, ignore_index=True)
    sig = sig[sig["Region"].eq("Full_Region") | sig["Region"].str.fullmatch(r"LD_Block_r2_080_Block\d+", na=False)].copy()
    sig["Haplotype_Effect"] = np.where(
        sig["OR"] > 1,
        "Tumour-enriched",
        np.where(sig["OR"] < 1, "Healthy-enriched", "Balanced"),
    )
    sig["Primary_Comparison_Rank"] = sig["Comparison"].map({"Breast": 0, "All": 1, "Endometrium": 2})
    sig = sig.sort_values(["FDR", "P_Fisher", "Primary_Comparison_Rank", "Region", "Haplotype_ID"], kind="stable")
    return sig.reset_index(drop=True)


def load_significant_snps(snp_results_path: Path) -> pd.DataFrame:
    comparison_map = {"Global": "All", "Breast": "Breast", "Endometrium": "Endometrium"}
    frames = []
    for sheet, comparison in comparison_map.items():
        df = read_excel(snp_results_path, sheet)
        df = normalise_stage12_snp_sheet(df, sheet)
        if df.empty:
            continue
        df["Comparison"] = comparison
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    sig = pd.concat(frames, ignore_index=True)
    odds = pd.to_numeric(sig["Odds_Ratio"], errors="coerce")
    dose = pd.to_numeric(sig.get("Dose_Difference", pd.Series(np.nan, index=sig.index)), errors="coerce")
    sig["SNP_Effect"] = np.select(
        [
            odds > 1,
            odds < 1,
            dose > 0,
            dose < 0,
        ],
        [
            "Tumour-enriched",
            "Healthy-enriched",
            "Tumour-enriched",
            "Healthy-enriched",
        ],
        default="Balanced",
    )
    return sig.sort_values(["Comparison", "FDR_P_Value", "P_Value", "SNP_ID"], kind="stable").reset_index(drop=True)


def load_region_lookup(
    hap_stats_path: Path,
    compare_path: Path,
) -> tuple[pd.DataFrame, dict[str, list[str]], dict[str, list[int]], dict[str, tuple[int, int]]]:
    snps_used = read_excel(hap_stats_path, "SNPs_Used").rename(columns={"POS_int": "POS", "rsID_clean": "rsID"})
    full_region_rsids = snps_used["rsID"].dropna().astype(str).tolist()
    full_region_pos = snps_used["POS"].dropna().astype(int).tolist()

    ld_blocks = read_excel(compare_path, "Study_LD_Blocks")
    ld_blocks["POS"] = ld_blocks["POS"].astype(int)
    ld_blocks["Region"] = ld_blocks["Region"].astype(str)

    region_to_rsids = {"Full_Region": full_region_rsids}
    region_to_pos = {"Full_Region": full_region_pos}
    region_to_span = {"Full_Region": (min(full_region_pos), max(full_region_pos))}
    for region, sub in ld_blocks.groupby("Region", sort=False):
        ordered = sub.sort_values("POS")
        rsids = ordered["rsID"].astype(str).tolist()
        pos = ordered["POS"].astype(int).tolist()
        region_to_rsids[region] = rsids
        region_to_pos[region] = pos
        region_to_span[region] = (min(pos), max(pos))
    return snps_used, region_to_rsids, region_to_pos, region_to_span


def load_freq_sheet_definitions(
    hap_stats_path: Path,
    region_to_rsids: dict[str, list[str]],
    regions: Iterable[str],
    allow_missing: bool = False,
) -> dict[str, pd.DataFrame]:
    out: dict[str, pd.DataFrame] = {}
    xl = pd.ExcelFile(hap_stats_path)
    available_sheets = set(xl.sheet_names)
    for region in regions:
        sheet = "Full_Global_Freqs" if region == "Full_Region" else f"{region}_Freqs"
        if sheet not in available_sheets:
            if allow_missing:
                continue
            raise ValueError(f"Worksheet named '{sheet}' not found in {hap_stats_path}")
        df = pd.read_excel(xl, sheet_name=sheet)
        rsid_cols = [col for col in region_to_rsids[region] if col in df.columns]
        keep = ["Haplotype_ID", "Frequency", "Count", "Variant_Content"] + rsid_cols
        out[region] = df[keep].copy()
    return out


def extract_snp_position(snp_id: str) -> int | None:
    text = str(snp_id or "").strip()
    if not text:
        return None
    match = re.search(r"\|(\d+)_", text)
    if match:
        return int(match.group(1))
    return None


def display_snp_label(snp_id: str, pos: int | None, pos_to_rsid: dict[int, str]) -> str:
    if pos is not None:
        rsid = pos_to_rsid.get(int(pos))
        if rsid and str(rsid).strip().lower() != "nan":
            return str(rsid)
    text = str(snp_id or "").strip()
    return text if text else "SNP"


def block_label_for_position(pos: int | None, region_to_span: dict[str, tuple[int, int]]) -> str:
    if pos is None:
        return "Unmapped"
    candidates = []
    for region, span in region_to_span.items():
        if region == "Full_Region":
            continue
        start, end = span
        if start <= pos <= end:
            candidates.append((end - start, region))
    if not candidates:
        full_span = region_to_span.get("Full_Region")
        if full_span is not None and full_span[0] <= pos <= full_span[1]:
            return "Outside LD blocks"
        return "Unmapped"
    _, region = sorted(candidates)[0]
    return region_display_label(region)


def region_key_for_position(pos: int | None, region_to_span: dict[str, tuple[int, int]]) -> str | None:
    if pos is None:
        return None
    candidates = []
    for region, span in region_to_span.items():
        if region == "Full_Region":
            continue
        start, end = span
        if start <= pos <= end:
            candidates.append((end - start, start, region))
    if candidates:
        return sorted(candidates)[0][2]
    full_span = region_to_span.get("Full_Region")
    if full_span is not None and full_span[0] <= pos <= full_span[1]:
        return "Full_Region"
    return None


def enrich_significant_snps(
    sig_snps: pd.DataFrame,
    rsid_to_pos: dict[str, int],
    pos_to_rsid: dict[int, str],
    region_to_span: dict[str, tuple[int, int]],
) -> pd.DataFrame:
    if sig_snps.empty:
        return sig_snps.copy()
    out = sig_snps.copy()
    rsid_series = out.get("rsID", pd.Series(pd.NA, index=out.index, dtype=object)).astype(str).str.strip()
    rsid_series = rsid_series.mask(rsid_series.eq("")).mask(rsid_series.str.lower().eq("nan"))
    pos_series = pd.to_numeric(out.get("POS", pd.Series(np.nan, index=out.index)), errors="coerce")
    pos_series = pos_series.combine_first(out["SNP_ID"].map(extract_snp_position))
    rsid_pos_series = pd.Series(
        [rsid_to_pos.get(value) if pd.notna(value) else np.nan for value in rsid_series.tolist()],
        index=out.index,
        dtype=float,
    )
    pos_series = pos_series.where(pos_series.notna(), rsid_pos_series)
    out["SNP_rsID"] = rsid_series
    out["SNP_POS"] = pos_series
    out["SNP_Display"] = [
        (
            str(rsid)
            if pd.notna(rsid) and str(rsid).strip().lower() != "nan"
            else display_snp_label(snp_id, int(pos) if pd.notna(pos) else None, pos_to_rsid)
        )
        for snp_id, rsid, pos in zip(out["SNP_ID"].astype(str).tolist(), out["SNP_rsID"].tolist(), out["SNP_POS"].tolist())
    ]
    out["Mapped_Study_Region"] = [
        region_key_for_position(int(pos) if pd.notna(pos) else None, region_to_span)
        for pos in out["SNP_POS"].tolist()
    ]
    out["Mapped_Study_Block_Label"] = [
        block_label_for_position(int(pos) if pd.notna(pos) else None, region_to_span)
        for pos in out["SNP_POS"].tolist()
    ]
    out["Mapped_Backbone_SNP"] = out["SNP_POS"].isin(list(pos_to_rsid.keys()))
    return out


def parse_gt(gt: str) -> tuple[int, int] | None:
    if pd.isna(gt):
        return None
    text = str(gt).strip()
    if text in {"", ".", "./.", ".|.", "nan"}:
        return None
    delim = "|" if "|" in text else "/" if "/" in text else None
    if delim is None:
        return None
    left, right = text.split(delim)
    if left not in {"0", "1"} or right not in {"0", "1"}:
        return None
    return int(left), int(right)


def get_sample_sets(meta: pd.DataFrame, comparison: str) -> tuple[list[str], list[str]]:
    if comparison == "Breast":
        sub = meta[meta["Cohort"].eq("Breast")].copy()
    elif comparison == "Endometrium":
        sub = meta[meta["Cohort"].str.contains("Endomet", case=False, na=False)].copy()
    elif comparison == "All":
        sub = meta.copy()
    else:
        raise ValueError(f"Unsupported comparison: {comparison}")
    tumour = sub[sub["Tissue"].eq("Tumour")]["Sample"].astype(str).tolist()
    healthy = sub[sub["Tissue"].eq("Healthy")]["Sample"].astype(str).tolist()
    return tumour, healthy


def region_target_bits(
    freq_defs: dict[str, pd.DataFrame],
    region: str,
    haplotype_id: str,
    region_rsids: list[str],
) -> tuple[str, dict[str, str]]:
    row = freq_defs[region].loc[freq_defs[region]["Haplotype_ID"].eq(haplotype_id)]
    if row.empty:
        raise KeyError(f"Haplotype {haplotype_id} not found in {region}")
    row = row.iloc[0]
    bits = []
    allele_state: dict[str, str] = {}
    for rsid in region_rsids:
        value = row.get(rsid)
        is_alt = pd.notna(value) and int(float(value)) == 2
        bits.append("1" if is_alt else "0")
        allele_state[rsid] = "ALT" if is_alt else "REF"
    return "".join(bits), allele_state

def fetch_gt_value(phased_by_pos: pd.DataFrame, pos: int, sample: str):
    value = phased_by_pos.loc[pos, sample]
    if isinstance(value, pd.Series):
        return value.iloc[0]
    return value


def compute_haplotype_dosages(
    phased_by_pos: pd.DataFrame,
    region_positions: list[int],
    sample_names: Iterable[str],
    target_bits: str,
) -> pd.Series:
    dosage = {}
    for sample in sample_names:
        if sample not in phased_by_pos.columns:
            continue
        bits_a = []
        bits_b = []
        callable_sample = True
        for pos in region_positions:
            if pos not in phased_by_pos.index:
                callable_sample = False
                break
            parsed = parse_gt(fetch_gt_value(phased_by_pos, pos, sample))
            if parsed is None:
                callable_sample = False
                break
            a, b = parsed
            bits_a.append(str(a))
            bits_b.append(str(b))
        if not callable_sample:
            continue
        dosage[sample] = int("".join(bits_a) == target_bits) + int("".join(bits_b) == target_bits)
    return pd.Series(dosage, name="Dosage")


def compute_snp_alt_dosage(phased_by_pos: pd.DataFrame, pos: int, sample_names: Iterable[str]) -> pd.Series:
    dosage = {}
    if pos not in phased_by_pos.index:
        return pd.Series(dtype=float, name="Alt_Dosage")
    for sample in sample_names:
        if sample not in phased_by_pos.columns:
            continue
        parsed = parse_gt(fetch_gt_value(phased_by_pos, pos, sample))
        if parsed is None:
            continue
        dosage[sample] = sum(parsed)
    return pd.Series(dosage, name="Alt_Dosage")


def dosage_shift_label(case_pct_1: float, ctrl_pct_1: float, case_pct_2: float, ctrl_pct_2: float) -> str:
    het_delta = abs(case_pct_1 - ctrl_pct_1)
    hom_delta = abs(case_pct_2 - ctrl_pct_2)
    if hom_delta > het_delta + 5:
        return "Homozygous-dominant shift"
    if het_delta > hom_delta + 5:
        return "Heterozygous-dominant shift"
    return "Mixed/balanced shift"


def summarise_haplotype_dosage(
    hap_rows: pd.DataFrame,
    freq_defs: dict[str, pd.DataFrame],
    region_to_rsids: dict[str, list[str]],
    region_to_pos: dict[str, list[int]],
    phased_by_pos: pd.DataFrame,
    meta: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for row in hap_rows.itertuples(index=False):
        tumour_samples, healthy_samples = get_sample_sets(meta, row.Comparison)
        target_bits, _ = region_target_bits(freq_defs, row.Region, row.Haplotype_ID, region_to_rsids[row.Region])
        tumour_dosage = compute_haplotype_dosages(phased_by_pos, region_to_pos[row.Region], tumour_samples, target_bits)
        healthy_dosage = compute_haplotype_dosages(phased_by_pos, region_to_pos[row.Region], healthy_samples, target_bits)
        if tumour_dosage.empty or healthy_dosage.empty:
            continue
        tumour_counts = tumour_dosage.value_counts().reindex([0, 1, 2], fill_value=0)
        healthy_counts = healthy_dosage.value_counts().reindex([0, 1, 2], fill_value=0)
        tumour_pct = (tumour_counts / tumour_counts.sum()) * 100.0
        healthy_pct = (healthy_counts / healthy_counts.sum()) * 100.0
        rows.append({
            "Region": row.Region,
            "Comparison": row.Comparison,
            "Haplotype_ID": row.Haplotype_ID,
            "Haplotype_Effect": row.Haplotype_Effect,
            "Target_Haplotype_Bits": target_bits,
            "N_Tumour_Callable": int(tumour_counts.sum()),
            "N_Healthy_Callable": int(healthy_counts.sum()),
            "Tumour_Dose0": int(tumour_counts[0]),
            "Tumour_Dose1": int(tumour_counts[1]),
            "Tumour_Dose2": int(tumour_counts[2]),
            "Healthy_Dose0": int(healthy_counts[0]),
            "Healthy_Dose1": int(healthy_counts[1]),
            "Healthy_Dose2": int(healthy_counts[2]),
            "Tumour_Dose0_%": float(tumour_pct[0]),
            "Tumour_Dose1_%": float(tumour_pct[1]),
            "Tumour_Dose2_%": float(tumour_pct[2]),
            "Healthy_Dose0_%": float(healthy_pct[0]),
            "Healthy_Dose1_%": float(healthy_pct[1]),
            "Healthy_Dose2_%": float(healthy_pct[2]),
            "Dominant_Dosage_Shift": dosage_shift_label(float(tumour_pct[1]), float(healthy_pct[1]), float(tumour_pct[2]), float(healthy_pct[2])),
        })
    return pd.DataFrame(rows)


def summarise_snp_genotypes(
    sig_snps: pd.DataFrame,
    rsid_to_pos: dict[str, int],
    phased_by_pos: pd.DataFrame,
    meta: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for snp in sig_snps.itertuples(index=False):
        pos = rsid_to_pos.get(snp.SNP_ID)
        if pos is None:
            continue
        tumour_samples, healthy_samples = get_sample_sets(meta, snp.Comparison)
        tumour_dosage = compute_snp_alt_dosage(phased_by_pos, pos, tumour_samples)
        healthy_dosage = compute_snp_alt_dosage(phased_by_pos, pos, healthy_samples)
        if tumour_dosage.empty or healthy_dosage.empty:
            continue
        tumour_counts = tumour_dosage.value_counts().reindex([0, 1, 2], fill_value=0)
        healthy_counts = healthy_dosage.value_counts().reindex([0, 1, 2], fill_value=0)
        tumour_pct = (tumour_counts / tumour_counts.sum()) * 100.0
        healthy_pct = (healthy_counts / healthy_counts.sum()) * 100.0
        tumour_mean = float(tumour_dosage.mean())
        healthy_mean = float(healthy_dosage.mean())
        phased_effect = "Balanced"
        if tumour_mean > healthy_mean + 0.02:
            phased_effect = "Tumour-enriched"
        elif tumour_mean < healthy_mean - 0.02:
            phased_effect = "Healthy-enriched"
        rows.append({
            "Comparison": snp.Comparison,
            "SNP_ID": snp.SNP_ID,
            "SNP_rsID": getattr(snp, "SNP_rsID", None),
            "SNP_Display": getattr(snp, "SNP_Display", None),
            "Gene": snp.Symbol,
            "LD_Block_ID": snp.LD_Block_ID,
            "Stage12_OR": float(snp.Odds_Ratio),
            "Stage12_FDR": float(snp.FDR_P_Value),
            "Stage12_Effect": snp.SNP_Effect,
            "Phased_GT_Effect": phased_effect,
            "Tumour_WT": int(tumour_counts[0]),
            "Tumour_Het": int(tumour_counts[1]),
            "Tumour_Hom": int(tumour_counts[2]),
            "Healthy_WT": int(healthy_counts[0]),
            "Healthy_Het": int(healthy_counts[1]),
            "Healthy_Hom": int(healthy_counts[2]),
            "Tumour_WT_%": float(tumour_pct[0]),
            "Tumour_Het_%": float(tumour_pct[1]),
            "Tumour_Hom_%": float(tumour_pct[2]),
            "Healthy_WT_%": float(healthy_pct[0]),
            "Healthy_Het_%": float(healthy_pct[1]),
            "Healthy_Hom_%": float(healthy_pct[2]),
            "Mean_Alt_Dosage_Tumour": tumour_mean,
            "Mean_Alt_Dosage_Healthy": healthy_mean,
            "Dominant_Genotype_Shift": dosage_shift_label(float(tumour_pct[1]), float(healthy_pct[1]), float(tumour_pct[2]), float(healthy_pct[2])),
        })
    return pd.DataFrame(rows)


def expected_snp_effect_from_haplotype(haplotype_effect: str, allele_state: str) -> str:
    if haplotype_effect not in {"Tumour-enriched", "Healthy-enriched"}:
        return "Balanced"
    if allele_state == "ALT":
        return haplotype_effect
    return "Healthy-enriched" if haplotype_effect == "Tumour-enriched" else "Tumour-enriched"


def build_haplotype_snp_links(
    primary_haps: pd.DataFrame,
    sig_snps: pd.DataFrame,
    snp_genotypes: pd.DataFrame,
    freq_defs: dict[str, pd.DataFrame],
    region_to_rsids: dict[str, list[str]],
    region_to_pos: dict[str, list[int]],
) -> pd.DataFrame:
    geno_lookup = snp_genotypes.set_index(["Comparison", "SNP_ID"], drop=False)
    base_columns = [
        "Region", "Comparison", "Haplotype_ID", "Haplotype_Effect", "Haplotype_FDR",
        "Linked_SNP_ID", "Linked_SNP_Display", "Linked_SNP_rsID", "Linked_SNP_POS",
        "Linked_SNP_Gene", "Linked_SNP_LD_Block_Source", "Linked_SNP_Mapped_Study_Region",
        "Linked_SNP_Mapped_Block_Label", "Link_Scope", "Exact_Region_Match",
        "Linked_SNP_Stage12_FDR", "Allele_State_In_Haplotype",
        "Expected_SNP_Effect_From_Haplotype", "Observed_SNP_Effect_From_PhasedGT",
        "Observed_SNP_Effect_From_Stage12", "Consistency_Call",
        "Tumour_WT", "Tumour_Het", "Tumour_Hom",
        "Healthy_WT", "Healthy_Het", "Healthy_Hom",
        "Dominant_Genotype_Shift",
    ]
    rows = []
    for hap in primary_haps.itertuples(index=False):
        region_rsids = region_to_rsids[hap.Region]
        region_pos = region_to_pos[hap.Region]
        region_rsids_set = {str(value).strip() for value in region_rsids if str(value).strip()}
        region_pos_set = {int(value) for value in region_pos}
        region_pos_to_rsid = {int(pos): str(rsid) for rsid, pos in zip(region_rsids, region_pos)}
        _, allele_states = region_target_bits(freq_defs, hap.Region, hap.Haplotype_ID, region_rsids)
        overlap = sig_snps[
            sig_snps["Comparison"].eq(hap.Comparison)
            & (
                sig_snps["SNP_rsID"].astype(str).isin(region_rsids_set)
                | pd.to_numeric(sig_snps["SNP_POS"], errors="coerce").isin(region_pos_set)
            )
        ].copy()
        if overlap.empty:
            rows.append({
                "Region": hap.Region,
                "Comparison": hap.Comparison,
                "Haplotype_ID": hap.Haplotype_ID,
                "Haplotype_Effect": hap.Haplotype_Effect,
                "Haplotype_FDR": hap.FDR,
                "Linked_SNP_ID": None,
                "Linked_SNP_Display": None,
                "Linked_SNP_rsID": None,
                "Linked_SNP_POS": np.nan,
                "Linked_SNP_Gene": None,
                "Linked_SNP_LD_Block_Source": None,
                "Linked_SNP_Mapped_Study_Region": None,
                "Linked_SNP_Mapped_Block_Label": None,
                "Link_Scope": "No overlapping significant SNP",
                "Exact_Region_Match": False,
                "Linked_SNP_Stage12_FDR": np.nan,
                "Allele_State_In_Haplotype": None,
                "Expected_SNP_Effect_From_Haplotype": None,
                "Observed_SNP_Effect_From_PhasedGT": None,
                "Observed_SNP_Effect_From_Stage12": None,
                "Consistency_Call": "No overlapping significant SNP",
                "Tumour_WT": np.nan,
                "Tumour_Het": np.nan,
                "Tumour_Hom": np.nan,
                "Healthy_WT": np.nan,
                "Healthy_Het": np.nan,
                "Healthy_Hom": np.nan,
                "Dominant_Genotype_Shift": None,
            })
            continue
        for snp in overlap.itertuples(index=False):
            geno_row = None
            if (hap.Comparison, snp.SNP_ID) in geno_lookup.index:
                maybe = geno_lookup.loc[(hap.Comparison, snp.SNP_ID)]
                geno_row = maybe.iloc[0] if isinstance(maybe, pd.DataFrame) else maybe
            snp_pos = int(snp.SNP_POS) if pd.notna(snp.SNP_POS) else None
            member_rsid = None
            if pd.notna(snp.SNP_rsID) and str(snp.SNP_rsID).strip() in allele_states:
                member_rsid = str(snp.SNP_rsID).strip()
            elif snp_pos is not None:
                member_rsid = region_pos_to_rsid.get(snp_pos)
            allele_state = allele_states.get(member_rsid) if member_rsid is not None else None
            if allele_state in {"ALT", "REF"}:
                expected = expected_snp_effect_from_haplotype(hap.Haplotype_Effect, allele_state)
            else:
                expected = "Ambiguous"
            observed = geno_row["Phased_GT_Effect"] if geno_row is not None else None
            if observed in {None, "Balanced"}:
                consistency = "Ambiguous"
            elif observed == expected:
                consistency = "Consistent"
            else:
                consistency = "Discordant"
            exact_region_match = str(snp.Mapped_Study_Region) == str(hap.Region)
            if hap.Region == "Full_Region":
                link_scope = "Full-region umbrella"
            elif exact_region_match:
                link_scope = "Exact LD-block member"
            else:
                link_scope = "Region member"
            rows.append({
                "Region": hap.Region,
                "Comparison": hap.Comparison,
                "Haplotype_ID": hap.Haplotype_ID,
                "Haplotype_Effect": hap.Haplotype_Effect,
                "Haplotype_FDR": hap.FDR,
                "Linked_SNP_ID": snp.SNP_ID,
                "Linked_SNP_Display": snp.SNP_Display,
                "Linked_SNP_rsID": snp.SNP_rsID,
                "Linked_SNP_POS": snp_pos,
                "Linked_SNP_Gene": snp.Symbol,
                "Linked_SNP_LD_Block_Source": snp.LD_Block_ID,
                "Linked_SNP_Mapped_Study_Region": snp.Mapped_Study_Region,
                "Linked_SNP_Mapped_Block_Label": snp.Mapped_Study_Block_Label,
                "Link_Scope": link_scope,
                "Exact_Region_Match": bool(exact_region_match),
                "Linked_SNP_Stage12_FDR": float(snp.FDR_P_Value),
                "Allele_State_In_Haplotype": allele_state,
                "Expected_SNP_Effect_From_Haplotype": expected,
                "Observed_SNP_Effect_From_PhasedGT": observed,
                "Observed_SNP_Effect_From_Stage12": snp.SNP_Effect,
                "Consistency_Call": consistency,
                "Tumour_WT": geno_row["Tumour_WT"] if geno_row is not None else np.nan,
                "Tumour_Het": geno_row["Tumour_Het"] if geno_row is not None else np.nan,
                "Tumour_Hom": geno_row["Tumour_Hom"] if geno_row is not None else np.nan,
                "Healthy_WT": geno_row["Healthy_WT"] if geno_row is not None else np.nan,
                "Healthy_Het": geno_row["Healthy_Het"] if geno_row is not None else np.nan,
                "Healthy_Hom": geno_row["Healthy_Hom"] if geno_row is not None else np.nan,
                "Dominant_Genotype_Shift": geno_row["Dominant_Genotype_Shift"] if geno_row is not None else None,
            })
    return pd.DataFrame(rows, columns=base_columns)


def build_significant_snp_haplotype_map(sig_snps: pd.DataFrame, hap_links: pd.DataFrame) -> pd.DataFrame:
    if sig_snps.empty:
        return pd.DataFrame()
    base = sig_snps.drop_duplicates(["Comparison", "SNP_ID"]).copy()
    if hap_links.empty or "Linked_SNP_ID" not in hap_links.columns:
        base["Primary_Haplotype_Link_Count"] = 0
        base["Exact_Block_Primary_Haplotype_Count"] = 0
        base["Linked_Primary_Haplotypes"] = ""
        base["Exact_Block_Primary_Haplotypes"] = ""
        base["Full_Region_Umbrella_Haplotypes"] = ""
        base["Link_Consistency_Calls"] = ""
        base["Mapping_Status"] = "No primary significant haplotype link"
        return base

    link_only = hap_links.dropna(subset=["Linked_SNP_ID"]).copy()
    if link_only.empty:
        base["Primary_Haplotype_Link_Count"] = 0
        base["Exact_Block_Primary_Haplotype_Count"] = 0
        base["Linked_Primary_Haplotypes"] = ""
        base["Exact_Block_Primary_Haplotypes"] = ""
        base["Full_Region_Umbrella_Haplotypes"] = ""
        base["Link_Consistency_Calls"] = ""
        base["Mapping_Status"] = "No primary significant haplotype link"
        return base

    def _join_unique(values: pd.Series) -> str:
        seen = []
        for value in values.astype(str):
            value = value.strip()
            if value and value.lower() != "nan" and value not in seen:
                seen.append(value)
        return "; ".join(seen)

    link_only["Haplotype_Label"] = link_only.apply(
        lambda row: f"{region_display_label(str(row['Region']))} {row['Haplotype_ID']} [{row['Comparison']}]",
        axis=1,
    )
    summary_rows = []
    for (comparison, snp_id), sub in link_only.groupby(["Comparison", "Linked_SNP_ID"], sort=False):
        summary_rows.append({
            "Comparison": comparison,
            "SNP_ID": snp_id,
            "Primary_Haplotype_Link_Count": int(sub["Haplotype_Label"].nunique()),
            "Exact_Block_Primary_Haplotype_Count": int(sub["Exact_Region_Match"].fillna(False).astype(bool).sum()),
            "Linked_Primary_Haplotypes": _join_unique(sub["Haplotype_Label"]),
            "Exact_Block_Primary_Haplotypes": _join_unique(sub.loc[sub["Exact_Region_Match"].fillna(False).astype(bool), "Haplotype_Label"]),
            "Full_Region_Umbrella_Haplotypes": _join_unique(sub.loc[sub["Region"].astype(str).eq("Full_Region"), "Haplotype_Label"]),
            "Link_Consistency_Calls": _join_unique(sub["Consistency_Call"]),
        })
    summary = pd.DataFrame(summary_rows)
    out = base.merge(summary, on=["Comparison", "SNP_ID"], how="left")
    for col in [
        "Primary_Haplotype_Link_Count",
        "Exact_Block_Primary_Haplotype_Count",
    ]:
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0).astype(int)
    for col in [
        "Linked_Primary_Haplotypes",
        "Exact_Block_Primary_Haplotypes",
        "Full_Region_Umbrella_Haplotypes",
        "Link_Consistency_Calls",
    ]:
        out[col] = out[col].fillna("")
    out["Mapping_Status"] = np.select(
        [
            out["Exact_Block_Primary_Haplotype_Count"] > 0,
            out["Primary_Haplotype_Link_Count"] > 0,
            out["Mapped_Study_Region"].astype(str).str.startswith("LD_Block_"),
        ],
        [
            "Mapped to exact primary LD-block haplotype",
            "Mapped only through full-region umbrella haplotype",
            "Significant SNP in backbone LD block but no primary haplotype link",
        ],
        default="Significant SNP outside labelled LD blocks or unmapped to backbone",
    )
    return out.sort_values(["Comparison", "FDR_P_Value", "P_Value", "SNP_ID"], kind="stable").reset_index(drop=True)


def best_primary_haplotype_rows(sig_haps: pd.DataFrame) -> pd.DataFrame:
    if sig_haps.empty:
        return sig_haps
    best = sig_haps.sort_values(["FDR", "P_Fisher", "Primary_Comparison_Rank"], kind="stable").groupby(["Region", "Haplotype_ID"], as_index=False).first()
    return best.sort_values(["FDR", "P_Fisher", "Region"], kind="stable").reset_index(drop=True)


def phased_haplotype_frequency(
    phased_by_pos: pd.DataFrame,
    region_positions: list[int],
    sample_names: Iterable[str],
    target_bits: str,
) -> tuple[float, int]:
    dosage = compute_haplotype_dosages(phased_by_pos, region_positions, sample_names, target_bits)
    if dosage.empty:
        return (math.nan, 0)
    n = int(dosage.shape[0])
    return (float(dosage.sum() / (2 * n)), n)

def build_frequency_checks(
    primary_haps: pd.DataFrame,
    freq_defs: dict[str, pd.DataFrame],
    region_to_rsids: dict[str, list[str]],
    region_to_pos: dict[str, list[int]],
    phased_by_pos: pd.DataFrame,
    meta: pd.DataFrame,
    hap_stats_female_path: Path,
    compare_path: Path,
    female_subset_info: dict[str, object],
) -> pd.DataFrame:
    fixed_compare = read_excel(compare_path, "Fixed_Region_Compare")
    conc = read_excel(compare_path, "Concordance_Summary")
    female_verified = bool(female_subset_info.get("verified"))
    female_freq_defs = (
        load_freq_sheet_definitions(
            hap_stats_female_path,
            region_to_rsids,
            primary_haps["Region"].dropna().unique(),
            allow_missing=True,
        )
        if female_verified else {}
    )
    all_samples = meta["Sample"].astype(str).tolist()
    rows = []
    for hap in primary_haps.itertuples(index=False):
        target_bits, _ = region_target_bits(freq_defs, hap.Region, hap.Haplotype_ID, region_to_rsids[hap.Region])
        freq_row = freq_defs[hap.Region].loc[freq_defs[hap.Region]["Haplotype_ID"].eq(hap.Haplotype_ID)].iloc[0]
        if female_verified and hap.Region in female_freq_defs:
            female_row = female_freq_defs[hap.Region].loc[female_freq_defs[hap.Region]["Haplotype_ID"].eq(hap.Haplotype_ID)]
            female_freq = float(female_row.iloc[0]["Frequency"]) if not female_row.empty else math.nan
        else:
            female_freq = math.nan

        tumour_samples, healthy_samples = get_sample_sets(meta, hap.Comparison)
        overall_freq, n_all = phased_haplotype_frequency(phased_by_pos, region_to_pos[hap.Region], all_samples, target_bits)
        tumour_freq, n_tumour = phased_haplotype_frequency(phased_by_pos, region_to_pos[hap.Region], tumour_samples, target_bits)
        healthy_freq, n_healthy = phased_haplotype_frequency(phased_by_pos, region_to_pos[hap.Region], healthy_samples, target_bits)

        eur_row = fixed_compare[
            fixed_compare["Region"].eq(hap.Region)
            & fixed_compare["Population"].eq("EUR")
            & fixed_compare["Study_Haplotype_ID"].eq(hap.Haplotype_ID)
        ]
        eur_ref_freq = float(eur_row.iloc[0]["Reference_Frequency"]) if not eur_row.empty else math.nan
        eur_abs_delta = float(eur_row.iloc[0]["Abs_Frequency_Difference"]) if not eur_row.empty else math.nan

        conc_row = conc[conc["Region"].eq(hap.Region) & conc["Population"].eq("EUR")]
        shared_mass = float(conc_row.iloc[0]["Shared_Frequency_Mass"]) if not conc_row.empty and pd.notna(conc_row.iloc[0]["Shared_Frequency_Mass"]) else math.nan
        jsd = float(conc_row.iloc[0]["Jensen_Shannon_Distance"]) if not conc_row.empty and pd.notna(conc_row.iloc[0]["Jensen_Shannon_Distance"]) else math.nan
        exact_match = bool(conc_row.iloc[0]["Exact_Block_Match"]) if not conc_row.empty and pd.notna(conc_row.iloc[0]["Exact_Block_Match"]) else math.nan

        rows.append({
            "Region": hap.Region,
            "Comparison": hap.Comparison,
            "Haplotype_ID": hap.Haplotype_ID,
            "Assoc_Global_Freq": float(hap.Global_Freq),
            "FreqSheet_Global_Freq": float(freq_row["Frequency"]),
            "Phased_Global_Freq_Recomputed": overall_freq,
            "Abs_Delta_Assoc_vs_Recomputed": abs(float(hap.Global_Freq) - overall_freq),
            "Assoc_Tumour_Freq": float(hap.Freq_Tumour),
            "Phased_Tumour_Freq_Recomputed": tumour_freq,
            "Abs_Delta_Tumour": abs(float(hap.Freq_Tumour) - tumour_freq),
            "Assoc_Healthy_Freq": float(hap.Freq_Healthy),
            "Phased_Healthy_Freq_Recomputed": healthy_freq,
            "Abs_Delta_Healthy": abs(float(hap.Freq_Healthy) - healthy_freq),
            "Callable_All_Samples_Recomputed": n_all,
            "Callable_Tumour_Samples_Recomputed": n_tumour,
            "Callable_Healthy_Samples_Recomputed": n_healthy,
            "FemaleOnly_Study_Global_Freq": female_freq,
            "FemaleOnly_Study_Available": bool(female_subset_info.get("available")),
            "FemaleOnly_Study_Verified": female_verified,
            "FemaleOnly_Study_Filter_Mode": str(female_subset_info.get("sample_filter_mode", "")),
            "FemaleOnly_Study_Filter_File": str(female_subset_info.get("sample_filter_file", "")),
            "FemaleOnly_Study_Note": str(female_subset_info.get("note", "")),
            "Abs_Delta_Main_vs_FemaleOnly": abs(float(freq_row["Frequency"]) - female_freq) if pd.notna(female_freq) else math.nan,
            "EUR_Female_Reference_Freq": eur_ref_freq,
            "Abs_Delta_Study_vs_EUR_Female": eur_abs_delta,
            "EUR_Shared_Frequency_Mass": shared_mass,
            "EUR_Jensen_Shannon_Distance": jsd,
            "EUR_Exact_Block_Match": exact_match,
        })
    return pd.DataFrame(rows)


def build_dataset_consistency(sig_haps: pd.DataFrame, frequency_checks: pd.DataFrame) -> pd.DataFrame:
    if sig_haps.empty:
        return pd.DataFrame()
    pivot = sig_haps.pivot_table(index=["Region", "Haplotype_ID"], columns="Comparison", values=["OR", "FDR", "Freq_Tumour", "Freq_Healthy"], aggfunc="first")
    pivot.columns = ["_".join([str(part) for part in col if str(part) != ""]).strip("_") for col in pivot.columns]
    out = pivot.reset_index().merge(
        frequency_checks[[
            "Region", "Haplotype_ID", "FemaleOnly_Study_Global_Freq", "Abs_Delta_Main_vs_FemaleOnly",
            "EUR_Female_Reference_Freq", "Abs_Delta_Study_vs_EUR_Female", "EUR_Shared_Frequency_Mass",
            "EUR_Jensen_Shannon_Distance", "EUR_Exact_Block_Match",
        ]],
        on=["Region", "Haplotype_ID"],
        how="left",
    )

    def classify(row: pd.Series) -> str:
        breast_sig = pd.notna(row.get("FDR_Breast")) and row.get("FDR_Breast", 1.0) <= FDR_THRESHOLD
        endo_sig = pd.notna(row.get("FDR_Endometrium")) and row.get("FDR_Endometrium", 1.0) <= FDR_THRESHOLD
        all_sig = pd.notna(row.get("FDR_All")) and row.get("FDR_All", 1.0) <= FDR_THRESHOLD
        if breast_sig and endo_sig:
            return "Shared across tumour cohorts"
        if breast_sig:
            return "Breast-driven"
        if endo_sig:
            return "Endometrium-driven"
        if all_sig:
            return "Global only"
        return "Not significant after deduplication"

    out["Study_Cohort_Consistency"] = out.apply(classify, axis=1)
    return out.sort_values(["Region", "Haplotype_ID"], kind="stable")


def build_primary_haplotype_table(
    primary_haps: pd.DataFrame,
    hap_links: pd.DataFrame,
    frequency_checks: pd.DataFrame,
    region_to_rsids: dict[str, list[str]],
) -> pd.DataFrame:
    linked = hap_links.dropna(subset=["Linked_SNP_ID"]).copy()
    summary_rows = []
    for (region, comparison, haplotype_id), sub in linked.groupby(["Region", "Comparison", "Haplotype_ID"], sort=False):
        summary_rows.append({
            "Region": region,
            "Comparison": comparison,
            "Haplotype_ID": haplotype_id,
            "Overlapping_Significant_SNPs": "; ".join(dict.fromkeys(sub["Linked_SNP_Display"].astype(str))),
            "Exact_Block_Linked_SNPs": "; ".join(dict.fromkeys(sub.loc[sub["Exact_Region_Match"].fillna(False).astype(bool), "Linked_SNP_Display"].astype(str))),
            "Full_Region_Umbrella_SNPs": "; ".join(dict.fromkeys(sub.loc[sub["Region"].astype(str).eq("Full_Region"), "Linked_SNP_Display"].astype(str))),
            "Consistency_Summary": "; ".join(sorted(set(sub["Consistency_Call"].astype(str)))),
        })
    link_summary = pd.DataFrame(summary_rows)
    out = primary_haps.merge(link_summary, on=["Region", "Comparison", "Haplotype_ID"], how="left")
    out = out.merge(
        frequency_checks[[
            "Region", "Comparison", "Haplotype_ID", "Phased_Global_Freq_Recomputed",
            "Abs_Delta_Assoc_vs_Recomputed", "EUR_Female_Reference_Freq",
            "EUR_Shared_Frequency_Mass", "EUR_Jensen_Shannon_Distance",
        ]],
        on=["Region", "Comparison", "Haplotype_ID"],
        how="left",
    )
    out["Region_SNPs"] = out["Region"].map(lambda region: "; ".join(region_to_rsids.get(region, [])))
    return out.sort_values(["FDR", "P_Fisher"], kind="stable")


def region_display_label(region: str) -> str:
    if region == "Full_Region":
        return "Full region"
    match = re.fullmatch(r"LD_Block_r2_080_Block(\d+)", str(region))
    if match:
        return f"LD block {match.group(1)}"
    return str(region).replace("_", " ")


def place_staggered_labels(
    ax: plt.Axes,
    items: list[tuple[float, str]],
    base_y: float,
    direction: int,
    row_step: float,
    *,
    fontsize: float = 7.5,
    max_rows: int = 4,
    text_color: str = "#303030",
    bbox_fc: str = "white",
    bbox_ec: str = "#d0d0d0",
    rotation: float = 0.0,
) -> None:
    if not items:
        return
    rows_end: list[float] = []
    for x, label in sorted(items, key=lambda pair: (pair[0], pair[1])):
        label_width_mb = max(0.010, 0.00115 * len(str(label)))
        half_width = label_width_mb / 2.0
        row_idx = None
        for idx, row_end in enumerate(rows_end):
            if x - half_width > row_end + 0.0015:
                row_idx = idx
                rows_end[idx] = x + half_width
                break
        if row_idx is None:
            row_idx = len(rows_end)
            rows_end.append(x + half_width)
        row_idx = min(row_idx, max_rows - 1)
        y = base_y + (direction * row_step * row_idx)
        ax.text(
            x,
            y,
            str(label),
            ha="center",
            va="bottom" if direction > 0 else "top",
            fontsize=fontsize,
            rotation=rotation,
            color=text_color,
            bbox=dict(boxstyle="round,pad=0.16", fc=bbox_fc, ec=bbox_ec, alpha=0.92),
            zorder=8,
            clip_on=False,
        )


def build_ld_block_map_audit(
    primary_haps: pd.DataFrame,
    region_to_span: dict[str, tuple[int, int]],
    region_to_pos: dict[str, list[int]],
    sig_snps: pd.DataFrame,
    hap_links: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    snp_region_summary = (
        sig_snps.dropna(subset=["Mapped_Study_Region"])
        .groupby("Mapped_Study_Region", as_index=False)
        .agg(
            Significant_SNP_Count=("SNP_Display", "nunique"),
            Significant_SNPs=("SNP_Display", lambda s: "; ".join(dict.fromkeys(s.astype(str)))),
        )
        .rename(columns={"Mapped_Study_Region": "Region"})
    ) if not sig_snps.empty else pd.DataFrame(columns=["Region", "Significant_SNP_Count", "Significant_SNPs"])
    linked_region_summary = (
        hap_links.dropna(subset=["Linked_SNP_ID"])
        .groupby("Region", as_index=False)
        .agg(
            Primary_Linked_SNP_Count=("Linked_SNP_Display", "nunique"),
            Primary_Linked_SNPs=("Linked_SNP_Display", lambda s: "; ".join(dict.fromkeys(s.astype(str)))),
        )
    ) if not hap_links.empty else pd.DataFrame(columns=["Region", "Primary_Linked_SNP_Count", "Primary_Linked_SNPs"])
    region_order = sorted(region_to_span, key=lambda region: region_to_span[region][0])
    for region in region_order:
        sub = primary_haps[primary_haps["Region"].eq(region)].copy()
        start, end = region_to_span[region]
        row_count = len(sub)
        row_labels = "; ".join(
            f"{hap_id} ({comparison})"
            for hap_id, comparison in sub[["Haplotype_ID", "Comparison"]].itertuples(index=False)
        )
        if row_count:
            status = "Shown as labelled haplotype row(s) in the genomic map"
            note = "This region has at least one primary significant haplotype row in the map."
        elif str(region).startswith("LD_Block_"):
            status = "Shown on the shared LD-block track only"
            note = "The block is valid in the study LD backbone but has no primary significant haplotype row after stage-22 filtering."
        else:
            status = "Reference span only"
            note = "Full-region backbone reference span."
        rows.append(
            {
                "Region": region,
                "Display_Label": region_display_label(region),
                "Start_Mb": round(start / 1e6, 6),
                "End_Mb": round(end / 1e6, 6),
                "SNP_Count": len(region_to_pos.get(region, [])),
                "Primary_Significant_Row_Count": row_count,
                "Primary_Rows": row_labels,
                "Significant_SNP_Count": 0,
                "Significant_SNPs": "",
                "Primary_Linked_SNP_Count": 0,
                "Primary_Linked_SNPs": "",
                "Map_Status": status,
                "Map_Note": note,
            }
        )
    out = pd.DataFrame(rows)
    out = out.merge(snp_region_summary, on="Region", how="left", suffixes=("", "_sig"))
    out = out.merge(linked_region_summary, on="Region", how="left", suffixes=("", "_link"))
    for target, fallback in [
        ("Significant_SNP_Count", "Significant_SNP_Count_sig"),
        ("Significant_SNPs", "Significant_SNPs_sig"),
        ("Primary_Linked_SNP_Count", "Primary_Linked_SNP_Count_link"),
        ("Primary_Linked_SNPs", "Primary_Linked_SNPs_link"),
    ]:
        if fallback in out.columns:
            out[target] = out[target].where(out[target].notna(), out[fallback])
            out = out.drop(columns=[fallback])
    out["Significant_SNP_Count"] = pd.to_numeric(out["Significant_SNP_Count"], errors="coerce").fillna(0).astype(int)
    out["Primary_Linked_SNP_Count"] = pd.to_numeric(out["Primary_Linked_SNP_Count"], errors="coerce").fillna(0).astype(int)
    out["Significant_SNPs"] = out["Significant_SNPs"].fillna("")
    out["Primary_Linked_SNPs"] = out["Primary_Linked_SNPs"].fillna("")
    return out


def build_snp_genomic_mapping_audit(
    sig_snps: pd.DataFrame,
    hap_links: pd.DataFrame,
    region_to_pos: dict[str, list[int]],
) -> pd.DataFrame:
    if sig_snps.empty:
        return pd.DataFrame()
    snp_map = build_significant_snp_haplotype_map(sig_snps, hap_links)
    out = snp_map.copy()
    out["Mapped_Region_SNP_Count"] = out["Mapped_Study_Region"].map(lambda region: len(region_to_pos.get(region, [])) if pd.notna(region) else 0)
    out["Mapping_Note"] = np.where(
        out["Primary_Haplotype_Link_Count"] > 0,
        "Mapped to at least one primary significant haplotype row in the stage-22 interpretation layer.",
        np.where(
            out["Mapped_Study_Region"].astype(str).str.startswith("LD_Block_"),
            "Mapped to a backbone LD block, but that block has no linked primary significant haplotype row in the same comparison.",
            "Mapped outside labelled LD blocks or absent from the phased/backbone study region.",
        ),
    )
    keep = [
        "Comparison", "SNP_ID", "SNP_Display", "SNP_rsID", "SNP_POS", "Gene",
        "LD_Block_ID", "Mapped_Study_Region", "Mapped_Study_Block_Label", "Mapped_Backbone_SNP",
        "Mapped_Region_SNP_Count", "Primary_Haplotype_Link_Count", "Exact_Block_Primary_Haplotype_Count",
        "Linked_Primary_Haplotypes", "Exact_Block_Primary_Haplotypes", "Full_Region_Umbrella_Haplotypes",
        "Odds_Ratio", "P_Value", "FDR_P_Value", "Mapping_Status", "Mapping_Note",
    ]
    available = [col for col in keep if col in out.columns]
    return out[available].sort_values(["Comparison", "FDR_P_Value", "P_Value", "SNP_ID"], kind="stable").reset_index(drop=True)


def plot_shared_genomic_map(
    primary_haps: pd.DataFrame,
    sig_snps: pd.DataFrame,
    hap_links: pd.DataFrame,
    region_to_span: dict[str, tuple[int, int]],
    region_to_pos: dict[str, list[int]],
    pos_to_rsid: dict[int, str],
    out_path: Path,
) -> None:
    if primary_haps.empty:
        return

    plot_haps = primary_haps.copy()
    plot_haps["Region_Start"] = plot_haps["Region"].map(lambda region: region_to_span[region][0])
    plot_haps = plot_haps.sort_values(["Region_Start", "FDR"], kind="stable").reset_index(drop=True)
    plot_haps["y"] = list(range(len(plot_haps), 0, -1))
    plot_haps["Row_Label"] = plot_haps.apply(
        lambda row: f"{region_display_label(row['Region'])} {row['Haplotype_ID']} ({row['Comparison']})",
        axis=1,
    )

    link_plot = hap_links.dropna(subset=["Linked_SNP_ID"]).copy()
    if not link_plot.empty:
        link_plot["SNP_POS"] = pd.to_numeric(link_plot["Linked_SNP_POS"], errors="coerce")
        link_plot = link_plot.dropna(subset=["SNP_POS"]).copy()
        if not link_plot.empty:
            link_plot["SNP_POS"] = link_plot["SNP_POS"].astype(int)
    if not link_plot.empty and "Linked_SNP_Display" in link_plot.columns:
        plot_haps = plot_haps.merge(
            link_plot.groupby(["Region", "Comparison", "Haplotype_ID"], as_index=False).agg(
                Linked_SNP_Label=("Linked_SNP_Display", lambda s: ", ".join(s.astype(str))),
            ),
            on=["Region", "Comparison", "Haplotype_ID"],
            how="left",
        )
    else:
        plot_haps["Linked_SNP_Label"] = pd.NA

    fig = plt.figure(figsize=(24, 19.0))
    gs = fig.add_gridspec(3, 1, height_ratios=[3.0, 1.15, 1.75], hspace=0.08)
    ax_top = fig.add_subplot(gs[0])
    ax_mid = fig.add_subplot(gs[1], sharex=ax_top)
    ax_bottom = fig.add_subplot(gs[2], sharex=ax_top)

    ld_track_y_top = 0.8
    ld_track_y_bottom = 0.85
    ld_regions = [region for region in region_to_span if str(region).startswith("LD_Block_")]
    ld_regions = sorted(ld_regions, key=lambda region: region_to_span[region][0])
    block_colors = sns.color_palette("tab20", n_colors=max(len(ld_regions), 1))
    block_label_items: list[tuple[float, str]] = []
    for idx, region in enumerate(ld_regions):
        start, end = region_to_span[region]
        track_color = block_colors[idx % len(block_colors)]
        ax_top.hlines(ld_track_y_top, start / 1e6, end / 1e6, color=track_color, lw=10, alpha=0.95, zorder=0)
        for pos in region_to_pos[region]:
            ax_top.vlines(pos / 1e6, ld_track_y_top - 0.10, ld_track_y_top + 0.10, color="white", lw=1.0, alpha=0.95, zorder=1)
        block_label_items.append(((start + end) / 2 / 1e6, region_display_label(region)))

    if not link_plot.empty:
        hap_y_lookup = plot_haps.set_index(["Region", "Comparison", "Haplotype_ID"])["y"].to_dict()
        for link in link_plot.itertuples(index=False):
            key = (link.Region, link.Comparison, link.Haplotype_ID)
            if key not in hap_y_lookup:
                continue
            x = link.SNP_POS / 1e6
            y = hap_y_lookup[key]
            color = HAP_EFFECT_COLORS.get(link.Haplotype_Effect, "#6c757d")
            exact = bool(link.Exact_Region_Match)
            ax_top.plot(
                [x, x],
                [0.24, y - 0.10],
                color=color,
                lw=1.15 if exact else 0.8,
                linestyle="-" if exact else "--",
                alpha=0.28 if exact else 0.14,
                zorder=1.6,
            )

    for row in plot_haps.itertuples(index=False):
        start, end = region_to_span[row.Region]
        y = row.y
        color = HAP_EFFECT_COLORS.get(row.Haplotype_Effect, "#6c757d")
        linewidth = 14 if row.Region == "Full_Region" else 8
        ax_top.hlines(y, start / 1e6, end / 1e6, color=color, lw=linewidth, alpha=0.9, zorder=1)
        for pos in region_to_pos[row.Region]:
            ax_top.vlines(pos / 1e6, y - 0.16, y + 0.16, color="white", lw=1.2, alpha=0.95, zorder=2)

        linked = link_plot[
            link_plot["Region"].eq(row.Region)
            & link_plot["Comparison"].eq(row.Comparison)
            & link_plot["Haplotype_ID"].eq(row.Haplotype_ID)
        ].drop_duplicates(["Linked_SNP_ID", "SNP_POS"])
        if not linked.empty:
            ax_top.scatter(
                linked["SNP_POS"] / 1e6,
                np.repeat(y, len(linked)),
                s=52,
                color="white",
                edgecolor="black",
                linewidth=0.8,
                zorder=4,
            )

    snp_plot = sig_snps.drop_duplicates(["Comparison", "SNP_ID"]).copy()
    snp_plot["SNP_POS"] = pd.to_numeric(snp_plot.get("SNP_POS", pd.Series(np.nan, index=snp_plot.index)), errors="coerce")
    snp_plot = snp_plot.dropna(subset=["SNP_POS"]).copy()
    snp_plot = snp_plot.loc[:, ~snp_plot.columns.duplicated()].copy()
    snp_plot["SNP_POS"] = snp_plot["SNP_POS"].astype(int)
    snp_plot["SNP_Display"] = snp_plot.get("SNP_Display", pd.Series(index=snp_plot.index, dtype=object)).fillna("")
    snp_plot["SNP_Display"] = np.where(
        snp_plot["SNP_Display"].astype(str).str.strip().ne(""),
        snp_plot["SNP_Display"],
        [
            display_snp_label(snp_id, int(pos), pos_to_rsid)
            for snp_id, pos in zip(snp_plot["SNP_ID"].astype(str).tolist(), snp_plot["SNP_POS"].tolist())
        ],
    )
    snp_plot["Block_Label"] = snp_plot.get("Mapped_Study_Block_Label", pd.Series(index=snp_plot.index, dtype=object)).fillna("")
    snp_plot["Block_Label"] = np.where(
        snp_plot["Block_Label"].astype(str).str.strip().ne(""),
        snp_plot["Block_Label"],
        [block_label_for_position(int(pos), region_to_span) for pos in snp_plot["SNP_POS"].tolist()],
    )
    top_snps = snp_plot.sort_values(["FDR_P_Value", "P_Value", "SNP_ID"], kind="stable").head(6).copy()
    top_snps = top_snps.sort_values("SNP_POS", kind="stable").reset_index(drop=True)
    top_snp_y = float(plot_haps["y"].max()) + 0.95

    if not top_snps.empty:
        for i, row in enumerate(top_snps.itertuples(index=False)):
            marker = COMPARISON_MARKERS.get(row.Comparison, "o")
            color = HAP_EFFECT_COLORS.get(row.SNP_Effect, "#6c757d")
            ax_top.scatter(row.SNP_POS / 1e6, top_snp_y, s=84, marker=marker, color=color, edgecolor="black", linewidth=0.8, zorder=6)
            top_snps.loc[top_snps.index[i], "Top_Label"] = f"{row.SNP_Display}\n{row.Block_Label}"

    top_label_items = [(row.SNP_POS / 1e6, f"{row.SNP_Display}\n{row.Block_Label}") for row in top_snps.itertuples(index=False)]

    ax_top.set_yticks([top_snp_y, ld_track_y_top] + plot_haps["y"].tolist())
    ax_top.set_yticklabels(["Top SNPs", "All study LD blocks"] + plot_haps["Row_Label"].tolist())
    ax_top.set_title(
        "Significant SNP-to-Haplotype Associations in the GSDMB Locus with LD Context",
        fontsize=16,
        pad=12,
    )
    ax_top.set_ylim(0.0, top_snp_y + 1.55)
    ax_top.grid(axis="x", linestyle=":", alpha=0.35)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.set_axisbelow(True)
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_top.tick_params(axis="y", labelsize=11)

    place_staggered_labels(
        ax_top,
        top_label_items,
        base_y=top_snp_y + 0.08,
        direction=1,
        row_step=0.38,
        fontsize=7.6,
        max_rows=3,
        bbox_ec="#bdbdbd",
        rotation=0.0,
    )

    snp_plot["y"] = 0.2
    bottom_label_items: list[tuple[float, str]] = []
    for row in snp_plot.itertuples(index=False):
        marker = COMPARISON_MARKERS.get(row.Comparison, "o")
        color = HAP_EFFECT_COLORS.get(row.SNP_Effect, "#6c757d")
        size = max(60, -math.log10(max(float(row.FDR_P_Value), 1e-12)) * 42)
        ax_bottom.scatter(row.SNP_POS / 1e6, row.y, marker=marker, s=size, color=color, edgecolor="black", linewidth=0.6, zorder=3)
        bottom_label_items.append((row.SNP_POS / 1e6, str(row.SNP_Display)))

    place_staggered_labels(
        ax_bottom,
        bottom_label_items,
        base_y=-0.14,
        direction=-1,
        row_step=0.18,
        fontsize=6.2,
        max_rows=6,
        rotation=0.0,
    )

    ax_bottom.set_xlabel("Chr17 position (Mb)")
    ax_bottom.set_yticks([0.2])
    ax_bottom.set_yticklabels(["Significant SNPs"])
    ax_bottom.set_ylim(-1.0, 0.95)
    ax_bottom.grid(axis="x", linestyle=":", alpha=0.35)
    ax_bottom.spines["top"].set_visible(False)
    ax_bottom.spines["right"].set_visible(False)
    ax_bottom.set_axisbelow(True)
    ax_bottom.tick_params(axis="y", labelsize=10)

    for idx, region in enumerate(ld_regions):
        start, end = region_to_span[region]
        track_color = block_colors[idx % len(block_colors)]
        ax_mid.hlines(ld_track_y_bottom, start / 1e6, end / 1e6, color=track_color, lw=10, alpha=0.95, zorder=0)
        for pos in region_to_pos[region]:
            ax_mid.vlines(pos / 1e6, ld_track_y_bottom - 0.10, ld_track_y_bottom + 0.10, color="white", lw=1.0, alpha=0.95, zorder=1)

    place_staggered_labels(
        ax_mid,
        block_label_items,
        base_y=ld_track_y_bottom + 0.12,
        direction=1,
        row_step=0.25,
        fontsize=7.0,
        max_rows=7,
        bbox_ec="#d0d0d0",
        rotation=35.0,
    )

    ax_mid.set_yticks([ld_track_y_bottom])
    ax_mid.set_yticklabels(["All study LD blocks"])
    ax_mid.set_ylim(0.2, 1.75)
    ax_mid.grid(axis="x", linestyle=":", alpha=0.35)
    ax_mid.spines["top"].set_visible(False)
    ax_mid.spines["right"].set_visible(False)
    ax_mid.set_axisbelow(True)
    ax_mid.tick_params(axis="x", labelbottom=False)
    ax_mid.tick_params(axis="y", labelsize=10)

    hap_handles = [plt.Line2D([0], [0], color=color, lw=8, label=label) for label, color in HAP_EFFECT_COLORS.items()]
    snp_handles = [plt.Line2D([0], [0], marker=marker, linestyle="", markerfacecolor="#666666", markeredgecolor="black", label=label) for label, marker in COMPARISON_MARKERS.items()]
    linked_handle = plt.Line2D([0], [0], marker="o", linestyle="", markerfacecolor="white", markeredgecolor="black", label="Linked SNP on haplotype row")
    exact_link_handle = plt.Line2D([0], [0], color="#666666", lw=1.2, linestyle="-", label="Exact SNP-to-haplotype block link")
    umbrella_link_handle = plt.Line2D([0], [0], color="#666666", lw=1.0, linestyle="--", label="Full-region umbrella link")
    block_handle = plt.Line2D([0], [0], color="#adb5bd", lw=8, label="All study LD blocks")
    ax_top.legend(
        handles=hap_handles + snp_handles + [linked_handle, exact_link_handle, umbrella_link_handle, block_handle],
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
        ncol=1,
        handletextpad=0.5,
        fontsize=11,
    )
    fig.subplots_adjust(left=0.06, right=0.82, top=0.95, bottom=0.08, hspace=0.10)
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)

def plot_dosage_genotype_heatmap(
    hap_dosage: pd.DataFrame,
    snp_genotypes: pd.DataFrame,
    primary_haps: pd.DataFrame,
    out_path: Path,
) -> None:
    top_haps = primary_haps.head(8)[["Region", "Comparison", "Haplotype_ID"]].copy()
    hap_heat = hap_dosage.merge(top_haps, on=["Region", "Comparison", "Haplotype_ID"], how="inner")
    rows = []
    for _, row in hap_heat.iterrows():
        rows.append({
            "Label": f"HAP {row['Region']} {row['Haplotype_ID']} [{row['Comparison']}]",
            "WT": row["Tumour_Dose0_%"] - row["Healthy_Dose0_%"],
            "Het": row["Tumour_Dose1_%"] - row["Healthy_Dose1_%"],
            "Hom": row["Tumour_Dose2_%"] - row["Healthy_Dose2_%"],
        })
    for _, row in snp_genotypes.sort_values(["Stage12_FDR", "Comparison", "SNP_ID"], kind="stable").head(12).iterrows():
        snp_label = row.get("SNP_rsID")
        if pd.isna(snp_label) or str(snp_label).strip().lower() == "nan" or str(snp_label).strip() == "":
            snp_label = row.get("SNP_Display")
        if pd.isna(snp_label) or str(snp_label).strip().lower() == "nan" or str(snp_label).strip() == "":
            snp_label = row["SNP_ID"]
        rows.append({
            "Label": f"SNP {snp_label} [{row['Comparison']}]",
            "WT": row["Tumour_WT_%"] - row["Healthy_WT_%"],
            "Het": row["Tumour_Het_%"] - row["Healthy_Het_%"],
            "Hom": row["Tumour_Hom_%"] - row["Healthy_Hom_%"],
        })
    heat_df = pd.DataFrame(rows)
    if heat_df.empty:
        return
    matrix = heat_df.set_index("Label")
    fig_height = max(6.0, 0.32 * len(matrix) + 2.2)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    sns.heatmap(
        matrix,
        cmap="RdBu_r",
        center=0,
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"label": "Tumour minus healthy (percentage points)"},
        ax=ax,
    )
    ax.set_title("Dosage and genotype separation for primary haplotypes and significant SNPs")
    ax.set_xlabel("Genotype class")
    ax.set_ylabel("")
    ax.set_xticklabels(["WT", "Het", "Hom ALT"], rotation=0)
    fig.tight_layout()
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def dataframe_from_lines(lines: list[str]) -> pd.DataFrame:
    return pd.DataFrame({"Notes": lines})


def main() -> None:
    paths = harmonised_defaults()
    out_dir = ensure_dir(paths["out_dir"])
    out_xlsx = out_dir / "22_Haplotype_First_Interpretation.xlsx"
    out_map = out_dir / "22_Haplotype_First_Genomic_Map.png"
    out_heatmap = out_dir / "22_Haplotype_SNP_Dose_Genotype_Shifts.png"

    sig_haps = load_significant_haplotypes(paths["hap_stats"])
    sig_snps = load_significant_snps(paths["snp_enrichment"])
    snps_used, region_to_rsids, region_to_pos, region_to_span = load_region_lookup(paths["hap_stats"], paths["hap_compare"])
    freq_defs = load_freq_sheet_definitions(paths["hap_stats"], region_to_rsids, sig_haps["Region"].dropna().unique())

    phased = pd.read_csv(paths["phased"], sep="\t", dtype=str)
    phased["POS"] = phased["POS"].astype(int)
    phased_by_pos = phased.set_index("POS")
    meta = pd.read_csv(paths["sample_meta"], sep="\t", dtype=str)
    meta = meta[meta["Sample"].isin(phased.columns[5:])].copy()
    rsid_to_pos = dict(zip(snps_used["rsID"].astype(str), snps_used["POS"].astype(int)))

    # Script 12 uses build_genomic_variant_id_series() which produces positional
    # IDs like "17|39909987_C_T" rather than rsIDs.  Bridge the gap by loading
    # the common-SNP frequency summary written by script 11, which carries both
    # Variant_ID and rsID side-by-side, so that either form resolves to a POS.
    from pipeline_utils import get_paths as _gp
    _snp_freq_summary = _gp().get("common_snps_dir", Path()) / "GSDMB_Common_SNP_Frequency_Summary.xlsx"
    if _snp_freq_summary.exists():
        _freq_bridge = pd.read_excel(_snp_freq_summary, usecols=["Variant_ID", "rsID"]).dropna(subset=["rsID"])
        _rsid_to_vid: dict[str, str] = dict(zip(_freq_bridge["rsID"].astype(str), _freq_bridge["Variant_ID"].astype(str)))
        _variantid_to_pos: dict[str, int] = {
            vid: pos
            for rsid, pos in rsid_to_pos.items()
            for vid in [_rsid_to_vid.get(rsid)]
            if vid is not None
        }
        if _variantid_to_pos:
            rsid_to_pos.update(_variantid_to_pos)
            print(f"  INFO: Bridged {len(_variantid_to_pos)} positional Variant_IDs into rsid_to_pos "
                  "via script-11 frequency summary.")
        else:
            print("  WARNING: Script-11 frequency summary found but no Variant_IDs could be bridged. "
                  "Check rsID overlap between SNPs_Used and GSDMB_Common_SNP_Frequency_Summary.xlsx.")
    else:
        print(f"  WARNING: {_snp_freq_summary} not found — rsid_to_pos contains rsIDs only. "
              "SNPs stored as positional Variant_IDs will not resolve and snp_genotypes will be empty.")

    primary_haps = best_primary_haplotype_rows(sig_haps)
    pos_to_rsid = {int(pos): str(rsid) for pos, rsid in zip(snps_used["POS"].astype(int), snps_used["rsID"].astype(str)) if pd.notna(pos)}
    sig_snps = enrich_significant_snps(sig_snps, rsid_to_pos, pos_to_rsid, region_to_span)
    ld_block_map_audit = build_ld_block_map_audit(primary_haps, region_to_span, region_to_pos, sig_snps, pd.DataFrame())
    hap_dosage = summarise_haplotype_dosage(primary_haps, freq_defs, region_to_rsids, region_to_pos, phased_by_pos, meta)
    snp_genotypes = summarise_snp_genotypes(sig_snps, rsid_to_pos, phased_by_pos, meta)
    hap_links = build_haplotype_snp_links(primary_haps, sig_snps, snp_genotypes, freq_defs, region_to_rsids, region_to_pos)
    snp_haplotype_map = build_significant_snp_haplotype_map(sig_snps, hap_links)
    snp_mapping_audit = build_snp_genomic_mapping_audit(sig_snps, hap_links, region_to_pos)
    ld_block_map_audit = build_ld_block_map_audit(primary_haps, region_to_span, region_to_pos, sig_snps, hap_links)
    female_subset_info = inspect_female_subset_workbook(paths["hap_stats_female"])
    frequency_checks = build_frequency_checks(primary_haps, freq_defs, region_to_rsids, region_to_pos, phased_by_pos, meta, paths["hap_stats_female"], paths["hap_compare"], female_subset_info)
    dataset_consistency = build_dataset_consistency(sig_haps, frequency_checks)
    primary_table = build_primary_haplotype_table(primary_haps, hap_links, frequency_checks, region_to_rsids)
    plot_shared_genomic_map(primary_haps, sig_snps, hap_links, region_to_span, region_to_pos, pos_to_rsid, out_map)
    plot_dosage_genotype_heatmap(hap_dosage, snp_genotypes, primary_haps, out_heatmap)

    notes = [
        "Primary significance source: stage-15 haplotype tumour-vs-healthy results restricted to Full_Region and r2_080 LD blocks.",
        "Supporting SNP significance source: stage-12 LD-block-aware SNP enrichment results (Global/Breast/Endometrium sheets).",
        "Genotype labels use the same convention as stage 17: 0/0 = WT, 0/1 or 1/0 = Het, 1/1 = Hom.",
        "Haplotype dosage is the haplotype analogue of WT/Het/Hom: dosage 0 / 1 / 2.",
        "Genotype and dosage summaries were recomputed from phased study genotypes so haplotype and SNP summaries use the same backbone.",
        "Significant SNP-to-haplotype links are mapped by bridged rsID/position backbone coordinates, not only by raw stage-12 SNP_ID labels.",
        "1000 Genomes reference context comes from the female-filtered comparison output generated in stage 19.",
        str(female_subset_info.get("note", "")),
        "The genomic map includes explicit SNP-to-LD-block guides plus SNP-to-haplotype link lines so exact block mappings are visually separated from full-region umbrella links.",
        "Full_Region links are intentionally broad umbrella mappings; exact LD-block links are reported separately in the SNP and haplotype mapping tables.",
    ]

    with pd.ExcelWriter(out_xlsx) as writer:
        dataframe_from_lines(notes).to_excel(writer, sheet_name="README", index=False)
        primary_table.to_excel(writer, sheet_name="Primary_Haplotypes", index=False)
        snp_haplotype_map.to_excel(writer, sheet_name="Significant_SNP_Hap_Map", index=False)
        hap_links.to_excel(writer, sheet_name="Haplotype_SNP_Comparison", index=False)
        snp_mapping_audit.to_excel(writer, sheet_name="SNP_Genomic_Mapping_Audit", index=False)
        hap_dosage.to_excel(writer, sheet_name="Haplotype_Dosage_Summary", index=False)
        snp_genotypes.to_excel(writer, sheet_name="SNP_Genotype_Summary", index=False)
        frequency_checks.to_excel(writer, sheet_name="Frequency_Checks", index=False)
        dataset_consistency.to_excel(writer, sheet_name="Dataset_Consistency", index=False)
        ld_block_map_audit.to_excel(writer, sheet_name="LD_Block_Map_Audit", index=False)
        snps_used.to_excel(writer, sheet_name="Backbone_SNPs", index=False)

    print(f"Wrote workbook: {out_xlsx}")
    print(f"Wrote map:      {out_map}")
    print(f"Wrote heatmap:  {out_heatmap}")


if __name__ == "__main__":
    main()
