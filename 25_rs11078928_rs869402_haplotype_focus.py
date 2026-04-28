#!/usr/bin/env python3
"""
Focused phased-haplotype and case-control analysis for rs11078928 and rs869402.

This script reuses the existing stage-15 phased outputs to answer a specific
haplotype question while also generating explicit SNP-level case-control
comparisons for presentation-ready reporting.

Outputs are written to:
  analysis_results/25_rs11078928_rs869402_haplotype_focus/
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from openpyxl import load_workbook
from openpyxl.drawing.image import Image as XLImage
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy.stats import fisher_exact


BASE = Path(__file__).resolve().parent
RESULTS = BASE / "analysis_results"
HAPLO_STATS_XLSX = RESULTS / "15_haplotype_statistics" / "GSDMB_Haplotype_Results.xlsx"
PHASED_TSV = RESULTS / "15_haplotype_phasing" / "phased_genotypes.tsv"
SAMPLE_METADATA_TSV = RESULTS / "15_haplotype_phasing" / "sample_metadata.tsv"
FORCED_GENOTYPE_XLSX = RESULTS / "05b_forced_genotypes" / "GSDMB_Forced_Genotypes_Union_Sites.xlsx"
OUT_DIR = RESULTS / "25_rs11078928_rs869402_haplotype_focus"

TARGET_RSIDS = ("rs11078928", "rs869402")
GENOTYPE_ORDER = ["WT", "Het", "Hom"]
GROUP_ORDER = [
    "Breast_Healthy",
    "Breast_Tumour",
    "Endometrium_Healthy",
    "Endometrium_Tumour",
    "Healthy_All",
    "Tumour_All",
    "All_Samples",
]
REPORT_GROUP_ORDER = [
    "Breast_Healthy",
    "Breast_Tumour",
    "Endometrium_Healthy",
    "Endometrium_Tumour",
    "Healthy_All",
    "Tumour_All",
]
GROUP_LABELS = {
    "Breast_Healthy": "Breast normal",
    "Breast_Tumour": "Breast tumour",
    "Endometrium_Healthy": "Endometrium normal",
    "Endometrium_Tumour": "Endometrium tumour",
    "Healthy_All": "Pooled normal controls",
    "Tumour_All": "Pooled tumours",
    "All_Samples": "All samples",
}
GROUP_SHORT_LABELS = {
    "Breast_Healthy": "Breast normal",
    "Breast_Tumour": "Breast tumour",
    "Endometrium_Healthy": "Endo normal",
    "Endometrium_Tumour": "Endo tumour",
    "Healthy_All": "Pooled normal",
    "Tumour_All": "Pooled tumour",
}
MINI_LABELS = {
    "00": "Neither_ALT",
    "01": "rs869402_only",
    "10": "rs11078928_only",
    "11": "Both_ALT",
}
# Earlier versions of the focused mini-haplotype summaries re-declared their
# own comparison pairs locally. The arm definitions happened to match the
# intended logic, but that duplication made the pooled comparison easy to
# misread as an arbitrary breast-versus-endometrium contrast. All tumour-vs-
# normal summaries now reuse this single explicit comparison specification:
# within-cohort comparisons for Breast and Endometrium, plus one pooled tumour
# versus pooled normal comparison across both cohorts.
COMPARISON_SPECS = [
    ("Breast", "Breast tumour vs pooled normal", "Healthy_All", "Breast_Tumour"),
    ("Endometrium", "Endometrium tumour vs pooled normal", "Healthy_All", "Endometrium_Tumour"),
    ("All", "Pooled tumour vs pooled normal", "Healthy_All", "Tumour_All"),
]
HEATMAP_SPECS = [
    ("Breast_Healthy", "Breast_Normal", "Breast normal controls"),
    ("Breast_Tumour", "Breast_Tumour", "Breast tumours"),
    ("Endometrium_Healthy", "Endometrium_Normal", "Endometrium normal controls"),
    ("Endometrium_Tumour", "Endometrium_Tumour", "Endometrial tumours"),
    ("Healthy_All", "Pooled_Normal", "Pooled normal controls"),
    ("Tumour_All", "Pooled_Tumour", "Pooled tumours"),
]
MATRIX_SHEET_MAP = {
    "Breast_Healthy": "Breast_Normal_GT_Matrix",
    "Breast_Tumour": "Breast_GT_Matrix",
    "Endometrium_Healthy": "Endometrium_Normal_GT_Matrix",
    "Endometrium_Tumour": "Endometrium_GT_Matrix",
    "Healthy_All": "Pooled_Normal_GT_Matrix",
    "Tumour_All": "Pooled_Tumour_GT_Matrix",
}


def ordered_group(df: pd.DataFrame, column: str) -> pd.DataFrame:
    out = df.copy()
    out[column] = pd.Categorical(out[column], categories=GROUP_ORDER, ordered=True)
    return out.sort_values(column, kind="stable").reset_index(drop=True)


def safe_fraction(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if denominator else float("nan")


def describe_comparison_logic(normal_group: str, tumour_group: str) -> str:
    return (
        f"Explicit tumour-versus-normal contrast: {GROUP_LABELS[tumour_group]} "
        f"versus {GROUP_LABELS[normal_group]}."
    )


def build_genotype_shift_rows(normal: pd.Series, tumour: pd.Series) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    for gt in GENOTYPE_ORDER:
        normal_pct = float(normal[f"{gt}_Pct_Callable"])
        tumour_pct = float(tumour[f"{gt}_Pct_Callable"])
        rows.append(
            {
                "Genotype": gt,
                "Normal_Pct": normal_pct,
                "Tumour_Pct": tumour_pct,
                "Delta": tumour_pct - normal_pct,
            }
        )
    return rows


def describe_dominant_genotype_shift(shift_rows: List[Dict[str, float]]) -> str:
    non_zero = [row for row in shift_rows if not math.isclose(abs(row["Delta"]), 0.0, abs_tol=1e-12)]
    if not non_zero:
        return "WT, Het, and Hom proportions were identical in tumours and normals."
    max_abs = max(abs(row["Delta"]) for row in non_zero)
    leaders = [
        row
        for row in non_zero
        if math.isclose(abs(row["Delta"]), max_abs, rel_tol=1e-9, abs_tol=1e-12)
    ]
    parts = []
    for row in leaders:
        direction = "higher" if row["Delta"] > 0 else "lower"
        gt_label = "Homozygous ALT" if row["Genotype"] == "Hom" else row["Genotype"]
        parts.append(
            f"{gt_label} proportion was {direction} in tumours "
            f"({row['Tumour_Pct']:.1%} vs {row['Normal_Pct']:.1%})"
        )
    return " and ".join(parts)


def describe_full_genotype_balance(shift_rows: List[Dict[str, float]]) -> str:
    parts = []
    for row in shift_rows:
        if math.isclose(abs(row["Delta"]), 0.0, abs_tol=1e-12):
            continue
        direction = "higher" if row["Delta"] > 0 else "lower"
        gt_label = "Homozygous ALT" if row["Genotype"] == "Hom" else row["Genotype"]
        parts.append(
            f"{gt_label}: {direction} in tumours by {abs(row['Delta']):.1%} "
            f"({row['Tumour_Pct']:.1%} vs {row['Normal_Pct']:.1%})"
        )
    return "; ".join(parts) if parts else "No genotype-balance difference detected."


def parse_phased_gt(value: object) -> Dict[str, object]:
    raw = "." if pd.isna(value) else str(value).strip()
    if raw in {"", ".", "./.", ".|."}:
        return {"raw": raw, "callable": False, "dose": pd.NA, "gt_label": "NoCall", "bits": None}

    separator = "|" if "|" in raw else "/" if "/" in raw else None
    if separator is None:
        return {"raw": raw, "callable": False, "dose": pd.NA, "gt_label": "NoCall", "bits": None}

    bits = raw.split(separator)
    if len(bits) != 2 or any(bit not in {"0", "1"} for bit in bits):
        return {"raw": raw, "callable": False, "dose": pd.NA, "gt_label": "NoCall", "bits": None}

    dose = sum(int(bit) for bit in bits)
    gt_label = {0: "WT", 1: "Het", 2: "Hom"}[dose]
    phased_bits = bits if separator == "|" else None
    return {"raw": raw, "callable": True, "dose": dose, "gt_label": gt_label, "bits": phased_bits}


def load_required_inputs() -> Tuple[pd.ExcelFile, pd.DataFrame, pd.DataFrame]:
    missing = [str(p) for p in [HAPLO_STATS_XLSX, PHASED_TSV, SAMPLE_METADATA_TSV] if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required input(s): " + ", ".join(missing))
    workbook = pd.ExcelFile(HAPLO_STATS_XLSX)
    phased = pd.read_csv(PHASED_TSV, sep="\t")
    sample_metadata = pd.read_csv(SAMPLE_METADATA_TSV, sep="\t")
    return workbook, phased, sample_metadata


def get_target_positions(workbook: pd.ExcelFile) -> pd.DataFrame:
    snps = pd.read_excel(workbook, sheet_name="SNPs_Used")
    target = snps[snps["rsID_clean"].isin(TARGET_RSIDS)].copy()
    if set(target["rsID_clean"]) != set(TARGET_RSIDS):
        found = sorted(target["rsID_clean"].dropna().astype(str).unique().tolist())
        raise ValueError(f"Could not find all target rsIDs in SNPs_Used. Found: {found}")
    return target[["rsID_clean", "POS_int", "Gene", "REF", "ALT"]].sort_values("POS_int", kind="stable")


def get_ld_context(workbook: pd.ExcelFile) -> Tuple[pd.DataFrame, pd.DataFrame]:
    blocks = pd.read_excel(workbook, sheet_name="LD_Blocks")
    rows: List[Dict[str, object]] = []
    for rsid in TARGET_RSIDS:
        hit = blocks[blocks["rsID"].astype(str).eq(rsid)].copy()
        if hit.empty:
            rows.append(
                {
                    "rsID": rsid,
                    "LD_Block_Status": "No explicit LD block at r2=0.80",
                    "LD_Block_Name": pd.NA,
                    "LD_Block_SNPs": pd.NA,
                    "Joint_Assessment_Within_LD_Block": "Not feasible",
                    "Notes": "This SNP is absent from the explicit r2=0.80 LD-block table in stage 15.",
                }
            )
            continue
        block_row = hit.iloc[0]
        block_name = f"LD_Block_{block_row['Threshold_Label']}_{block_row['Block']}"
        block_members = blocks[
            blocks["Threshold_Label"].astype(str).eq(str(block_row["Threshold_Label"]))
            & blocks["Block"].astype(str).eq(str(block_row["Block"]))
        ].sort_values("SNP_Order", kind="stable")
        rows.append(
            {
                "rsID": rsid,
                "LD_Block_Status": "Present in explicit LD block",
                "LD_Block_Name": block_name,
                "LD_Block_SNPs": "; ".join(block_members["rsID"].astype(str)),
                "Joint_Assessment_Within_LD_Block": "Feasible"
                if all(r in set(block_members["rsID"].astype(str)) for r in TARGET_RSIDS)
                else "Not feasible",
                "Notes": "Joint LD-block analysis requires both target SNPs to be in the same explicit block.",
            }
        )
    return pd.DataFrame(rows), blocks


def with_aggregate_groups(df: pd.DataFrame) -> pd.DataFrame:
    extra_frames = [
        df[df["Tissue"].astype(str).eq("Healthy")].assign(Group="Healthy_All"),
        df[df["Tissue"].astype(str).eq("Tumour")].assign(Group="Tumour_All"),
        df.assign(Group="All_Samples"),
    ]
    return pd.concat([df] + extra_frames, ignore_index=True)


def load_forced_genotype_inputs() -> Tuple[pd.ExcelFile, pd.DataFrame, pd.DataFrame]:
    if not FORCED_GENOTYPE_XLSX.exists():
        raise FileNotFoundError(f"Missing forced-genotype workbook: {FORCED_GENOTYPE_XLSX}")
    workbook = pd.ExcelFile(FORCED_GENOTYPE_XLSX)
    per_sample = pd.read_excel(workbook, sheet_name="Per_Sample_Genotypes")
    locus_summary = pd.read_excel(workbook, sheet_name="Locus_Summary")
    return workbook, per_sample, locus_summary


def build_forced_target_sample_table(per_sample: pd.DataFrame) -> pd.DataFrame:
    lookup = {
        "rs11078928": "chr17:39908216:T",
        "rs869402": "chr17:39911790:T",
    }
    tissue_map = {"normal": "Healthy", "tumour": "Tumour", "tumor": "Tumour"}
    keep = per_sample[per_sample["Locus_Key"].astype(str).isin(lookup.values())].copy()
    keep["Tissue"] = keep["Tissue"].astype(str).str.strip().str.lower().map(tissue_map)
    keep = keep[keep["Tissue"].notna()].copy()
    keep["Cohort"] = keep["Cohort"].astype(str).str.strip()
    keep["Group"] = keep["Cohort"] + "_" + keep["Tissue"]
    keep["Callable"] = keep["Meets_Callable_DP"].astype(bool)
    return keep


def make_forced_target_snp_status_table(target_positions: pd.DataFrame, locus_summary: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    lookup = target_positions.set_index("rsID_clean")
    for rsid in TARGET_RSIDS:
        locus_key = {
            "rs11078928": "chr17:39908216:T",
            "rs869402": "chr17:39911790:T",
        }[rsid]
        locus_matches = locus_summary[locus_summary["Locus_Key"].astype(str).eq(locus_key)].copy()
        if locus_matches.empty:
            raise KeyError(f"{locus_key} was not found in the forced-genotype locus summary.")
        # Some workbooks carry duplicate summary rows for the same locus; use
        # the widest all-sample row rather than assuming the key is unique.
        locus_row = (
            locus_matches.assign(N_Samples_Num=pd.to_numeric(locus_matches["N_Samples"], errors="coerce"))
            .sort_values(["N_Samples_Num"], ascending=False, kind="stable")
            .iloc[0]
        )
        callable_n = int(locus_row["N_WT"] + locus_row["N_Het_ALT"] + locus_row["N_Hom_ALT"])
        noncallable_n = int(locus_row["N_Below_DP_Threshold"] + locus_row["N_No_Call"] + locus_row["N_Other_ALT_Model"])
        rows.append(
            {
                "SNP": rsid,
                "POS": int(lookup.at[rsid, "POS_int"]),
                "Gene": str(lookup.at[rsid, "Gene"]),
                "REF": str(lookup.at[rsid, "REF"]),
                "ALT": str(lookup.at[rsid, "ALT"]),
                "Present_in_SNPs_Used": True,
                "Present_in_Forced_Genotype_Workbook": True,
                "Callable_Samples": callable_n,
                "NonCallable_Samples": noncallable_n,
                "Included_In_Group_Summary": True,
                "Included_In_Comparison_Summary": callable_n > 0,
                "Inclusion_Note": "Forced-genotype counts are used directly for the genotype distribution and carrier summary in this report.",
            }
        )
    return pd.DataFrame(rows)


def make_forced_target_snp_group_summary(per_sample: pd.DataFrame) -> pd.DataFrame:
    sample_table = build_forced_target_sample_table(per_sample)
    rows: List[Dict[str, object]] = []
    for rsid in TARGET_RSIDS:
        locus_key = {
            "rs11078928": "chr17:39908216:T",
            "rs869402": "chr17:39911790:T",
        }[rsid]
        locus_all = sample_table[sample_table["Locus_Key"].astype(str).eq(locus_key)].copy()
        locus_callable = locus_all[locus_all["Callable"]].copy()
        all_expanded = with_aggregate_groups(locus_all)
        callable_expanded = with_aggregate_groups(locus_callable)
        for group in GROUP_ORDER:
            sub_all = all_expanded[all_expanded["Group"].astype(str).eq(group)].copy()
            sub = callable_expanded[callable_expanded["Group"].astype(str).eq(group)].copy()
            wt = int((sub["Genotype_Status"].astype(str).eq("WT")).sum())
            het = int((sub["Genotype_Status"].astype(str).eq("Het_ALT")).sum())
            hom = int((sub["Genotype_Status"].astype(str).eq("Hom_ALT")).sum())
            callable_n = int(len(sub))
            total_n = int(len(sub_all))
            carriers = het + hom
            alt_alleles = het + 2 * hom
            rows.append(
                {
                    "SNP": rsid,
                    "Group": group,
                    "Group_Label": GROUP_LABELS.get(group, group),
                    "Total_Samples": total_n,
                    "Callable_Samples": callable_n,
                    "NonCallable_Samples": total_n - callable_n,
                    "WT_Count": wt,
                    "Het_Count": het,
                    "Hom_Count": hom,
                    "WT_Pct_Callable": safe_fraction(wt, callable_n),
                    "Het_Pct_Callable": safe_fraction(het, callable_n),
                    "Hom_Pct_Callable": safe_fraction(hom, callable_n),
                    "Carrier_Count": carriers,
                    "NonCarrier_Count": callable_n - carriers,
                    "Carrier_Pct_Callable": safe_fraction(carriers, callable_n),
                    "ALT_Allele_Count": alt_alleles,
                    "ALT_Allele_Denominator": int(2 * callable_n),
                    "ALT_Allele_Freq_Callable": safe_fraction(alt_alleles, 2 * callable_n),
                    "Percentages_Denominator": "Callable samples for this SNP",
                }
            )
    out = pd.DataFrame(rows)
    out["Group"] = pd.Categorical(out["Group"], categories=GROUP_ORDER, ordered=True)
    out = out.sort_values(["SNP", "Group"], kind="stable").reset_index(drop=True)
    out["Group"] = out["Group"].astype(str)
    return out


def extract_sample_genotypes(phased: pd.DataFrame, sample_metadata: pd.DataFrame, target_positions: pd.DataFrame) -> pd.DataFrame:
    target_pos = dict(zip(target_positions["rsID_clean"], target_positions["POS_int"]))
    phased_subset = phased[phased["POS"].isin(target_pos.values())].copy()
    if phased_subset["POS"].nunique() != len(TARGET_RSIDS):
        raise ValueError("Could not recover all target SNP rows from the phased genotype matrix.")

    phased_subset = phased_subset.set_index("POS")
    sample_cols = [c for c in phased_subset.columns if c not in {"CHROM", "ID", "REF", "ALT"}]
    rows: List[Dict[str, object]] = []
    for sample in sample_cols:
        row: Dict[str, object] = {"Sample": sample}
        both_callable = True
        both_phased = True
        for rsid, pos in target_pos.items():
            parsed = parse_phased_gt(phased_subset.at[int(pos), sample])
            row[f"{rsid}_GT_Raw"] = parsed["raw"]
            row[f"{rsid}_Callable"] = int(parsed["callable"])
            row[f"{rsid}_Dose"] = parsed["dose"] if parsed["callable"] else pd.NA
            row[f"{rsid}_GT"] = parsed["gt_label"]
            row[f"{rsid}_Carrier"] = int(parsed["dose"] > 0) if parsed["callable"] else pd.NA
            row[f"{rsid}_ALT_Allele_Count"] = int(parsed["dose"]) if parsed["callable"] else pd.NA
            both_callable = both_callable and bool(parsed["callable"])
            both_phased = both_phased and parsed["bits"] is not None
        row["Both_Target_SNPs_Callable"] = int(both_callable)
        row["Both_Target_SNPs_Phased"] = int(both_phased)
        rows.append(row)

    out = pd.DataFrame(rows).merge(sample_metadata, on="Sample", how="left")
    out["Cohort"] = out["Cohort"].fillna("Unknown")
    out["Tissue"] = out["Tissue"].fillna("Unknown")
    out["Group"] = out["Group"].fillna("Unknown")
    return out


def extract_haplotype_copies(sample_genotypes: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for record in sample_genotypes.itertuples(index=False):
        parsed_110 = parse_phased_gt(record.rs11078928_GT_Raw)
        parsed_869 = parse_phased_gt(record.rs869402_GT_Raw)
        if parsed_110["bits"] is None or parsed_869["bits"] is None:
            continue
        for copy_idx, (bit_110, bit_869) in enumerate(zip(parsed_110["bits"], parsed_869["bits"]), start=1):
            mini_code = f"{bit_110}{bit_869}"
            rows.append(
                {
                    "Sample": record.Sample,
                    "Cohort": record.Cohort,
                    "Tissue": record.Tissue,
                    "Group": record.Group,
                    "Haplotype_Copy": copy_idx,
                    "rs11078928_ALT": int(bit_110),
                    "rs869402_ALT": int(bit_869),
                    "Mini_Haplotype_Code": mini_code,
                    "Mini_Haplotype_Label": MINI_LABELS[mini_code],
                }
            )
    return pd.DataFrame(rows)

def make_mini_haplotype_frequency_table(hap_copies: pd.DataFrame) -> pd.DataFrame:
    expanded = with_aggregate_groups(hap_copies)
    freq = (
        expanded.groupby(["Group", "Mini_Haplotype_Code", "Mini_Haplotype_Label"], as_index=False)
        .size()
        .rename(columns={"size": "Haplotype_Copy_Count"})
    )
    totals = freq.groupby("Group")["Haplotype_Copy_Count"].transform("sum")
    freq["Haplotype_Copy_Frequency"] = freq["Haplotype_Copy_Count"] / totals
    freq["Callable_Sample_Count"] = (totals / 2).astype(int)
    freq["Frequency_Definition"] = (
        "Frequency = phased haplotype copies with this two-SNP pattern divided by all callable phased "
        "haplotype copies in the same group."
    )
    freq["Callable_Filter"] = (
        "Only samples with both target SNPs phased were included; each callable sample contributes two "
        "haplotype copies."
    )
    return ordered_group(freq, "Group")


def make_sample_level_tables(sample_genotypes: pd.DataFrame, hap_copies: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pattern_rows: List[Dict[str, object]] = []
    for sample, sub in hap_copies.groupby("Sample", sort=False):
        ordered_sub = sub.sort_values("Haplotype_Copy", kind="stable")
        labels = ordered_sub["Mini_Haplotype_Label"].astype(str).tolist()
        pattern_rows.append(
            {
                "Sample": sample,
                "Sample_Haplotype_Pattern": " + ".join(labels),
                "Carries_Both_ALT_Haplotype": int((ordered_sub["Mini_Haplotype_Code"] == "11").any()),
                "Carries_rs11078928_only_Haplotype": int((ordered_sub["Mini_Haplotype_Code"] == "10").any()),
                "Carries_rs869402_only_Haplotype": int((ordered_sub["Mini_Haplotype_Code"] == "01").any()),
            }
        )

    pattern_df = pd.DataFrame(pattern_rows)
    sample_level = sample_genotypes[
        [
            "Sample",
            "Cohort",
            "Tissue",
            "Group",
            "Both_Target_SNPs_Callable",
            "rs11078928_Dose",
            "rs869402_Dose",
            "rs11078928_GT",
            "rs869402_GT",
            "rs11078928_Carrier",
            "rs869402_Carrier",
        ]
    ].copy()
    sample_level = sample_level.rename(
        columns={
            "rs11078928_Carrier": "Carries_rs11078928_ALT",
            "rs869402_Carrier": "Carries_rs869402_ALT",
        }
    )
    sample_level = sample_level.merge(pattern_df, on="Sample", how="left")

    fill_zero_cols = [
        "Carries_rs11078928_ALT",
        "Carries_rs869402_ALT",
        "Carries_Both_ALT_Haplotype",
        "Carries_rs11078928_only_Haplotype",
        "Carries_rs869402_only_Haplotype",
    ]
    for col in fill_zero_cols:
        sample_level[col] = sample_level[col].fillna(0).astype(int)
    sample_level["Sample_Haplotype_Pattern"] = sample_level["Sample_Haplotype_Pattern"].fillna("No phased haplotype available")
    sample_level["Carries_Both_Variants_Across_Any_Haplotype"] = (
        (sample_level["Carries_rs11078928_ALT"] == 1) & (sample_level["Carries_rs869402_ALT"] == 1)
    ).astype(int)

    expanded = with_aggregate_groups(sample_level)
    carrier_cols = [
        "Carries_rs11078928_ALT",
        "Carries_rs869402_ALT",
        "Carries_Both_ALT_Haplotype",
        "Carries_rs11078928_only_Haplotype",
        "Carries_rs869402_only_Haplotype",
        "Carries_Both_Variants_Across_Any_Haplotype",
    ]
    summary = expanded.groupby("Group", as_index=False)[carrier_cols].sum()
    summary["Total_Sample_Count"] = expanded.groupby("Group").size().values
    summary["Callable_Sample_Count"] = expanded.groupby("Group")["Both_Target_SNPs_Callable"].sum().values
    summary["NonCallable_Sample_Count"] = summary["Total_Sample_Count"] - summary["Callable_Sample_Count"]
    for col in carrier_cols:
        summary[f"{col}_Proportion_of_Callable"] = [
            safe_fraction(v, d) for v, d in zip(summary[col], summary["Callable_Sample_Count"])
        ]
    summary = ordered_group(summary, "Group")

    patterns = (
        expanded[expanded["Both_Target_SNPs_Callable"] == 1]
        .groupby(["Group", "Sample_Haplotype_Pattern"], as_index=False)
        .size()
        .rename(columns={"size": "Sample_Count"})
    )
    if not patterns.empty:
        patterns["Callable_Sample_Denominator"] = patterns.groupby("Group")["Sample_Count"].transform("sum")
        patterns["Sample_Proportion_of_Callable"] = patterns["Sample_Count"] / patterns["Callable_Sample_Denominator"]
        patterns = ordered_group(patterns, "Group")
    return summary, patterns, expanded


def make_genotype_concordance_tables(expanded_samples: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    count_rows: List[pd.DataFrame] = []
    summary_rows: List[Dict[str, object]] = []
    for group in GROUP_ORDER:
        sub_all = expanded_samples[expanded_samples["Group"].astype(str).eq(group)].copy()
        if sub_all.empty:
            continue
        sub = sub_all[sub_all["Both_Target_SNPs_Callable"] == 1].copy()
        if not sub.empty:
            counts = (
                sub.groupby(["rs11078928_GT", "rs869402_GT"], as_index=False)
                .size()
                .rename(columns={"size": "Sample_Count"})
            )
            counts["Group"] = group
            counts["Row_Proportion_within_rs11078928_GT"] = counts.groupby("rs11078928_GT")["Sample_Count"].transform(
                lambda s: s / s.sum()
            )
            count_rows.append(
                counts[["Group", "rs11078928_GT", "rs869402_GT", "Sample_Count", "Row_Proportion_within_rs11078928_GT"]]
            )

        cond_110 = sub[sub["rs11078928_Dose"] > 0]
        cond_869 = sub[sub["rs869402_Dose"] > 0]

        def summarise_partner(df: pd.DataFrame, focal_gt: str) -> str:
            vc = df.loc[df["rs11078928_GT"] == focal_gt, "rs869402_GT"].value_counts()
            if vc.empty:
                return "n/a"
            order = [g for g in GENOTYPE_ORDER if g in vc.index]
            return ", ".join(f"{g}:{int(vc[g])}" for g in order)

        summary_rows.append(
            {
                "Group": group,
                "Group_Label": GROUP_LABELS.get(group, group),
                "Total_Sample_Count": int(len(sub_all)),
                "Callable_Sample_Count": int(len(sub)),
                "NonCallable_Sample_Count": int(len(sub_all) - len(sub)),
                "rs11078928_ALT_Carriers": int(len(cond_110)),
                "rs869402_ALT_Carriers": int(len(cond_869)),
                "P_rs869402_ALT_given_rs11078928_ALT": safe_fraction((cond_110["rs869402_Dose"] > 0).sum(), len(cond_110)),
                "P_rs11078928_ALT_given_rs869402_ALT": safe_fraction((cond_869["rs11078928_Dose"] > 0).sum(), len(cond_869)),
                "rs869402_GT_if_rs11078928_Het": summarise_partner(sub, "Het"),
                "rs869402_GT_if_rs11078928_Hom": summarise_partner(sub, "Hom"),
            }
        )

    counts_out = ordered_group(pd.concat(count_rows, ignore_index=True), "Group") if count_rows else pd.DataFrame()
    summary_out = ordered_group(pd.DataFrame(summary_rows), "Group") if summary_rows else pd.DataFrame()
    return counts_out, summary_out


def make_group_genotype_heatmap(genotype_counts: pd.DataFrame, group: str, filename_stem: str, title_label: str) -> Path:
    out_path = OUT_DIR / f"25_{filename_stem}_Genotype_Heatmap.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sub = genotype_counts[genotype_counts["Group"].astype(str).eq(group)].copy()
    if sub.empty:
        raise ValueError(f"{group} genotype counts were not available for heatmap generation.")

    count_matrix = (
        sub.pivot(index="rs11078928_GT", columns="rs869402_GT", values="Sample_Count")
        .reindex(index=GENOTYPE_ORDER, columns=GENOTYPE_ORDER)
        .fillna(0)
    )
    prop_matrix = (
        sub.pivot(index="rs11078928_GT", columns="rs869402_GT", values="Row_Proportion_within_rs11078928_GT")
        .reindex(index=GENOTYPE_ORDER, columns=GENOTYPE_ORDER)
        .fillna(0.0)
    )

    fig, ax = plt.subplots(figsize=(5.4, 4.2), dpi=220)
    im = ax.imshow(count_matrix.values, cmap="Blues", vmin=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(f"{title_label} callable sample count", fontsize=9)

    ax.set_xticks(range(len(GENOTYPE_ORDER)), GENOTYPE_ORDER)
    ax.set_yticks(range(len(GENOTYPE_ORDER)), GENOTYPE_ORDER)
    ax.set_xlabel("rs869402 genotype", fontsize=10)
    ax.set_ylabel("rs11078928 genotype", fontsize=10)
    ax.set_title(
        f"{title_label}: genotype pairing between rs11078928 and rs869402\nCounts with row percentages within rs11078928 genotype",
        fontsize=11,
        pad=10,
    )

    max_count = max(float(count_matrix.to_numpy().max()), 1.0)
    for i, row_label in enumerate(GENOTYPE_ORDER):
        for j, col_label in enumerate(GENOTYPE_ORDER):
            count = int(count_matrix.loc[row_label, col_label])
            prop = float(prop_matrix.loc[row_label, col_label])
            label = f"{count}\n({prop * 100:.0f}%)" if count > 0 else "0"
            color = "white" if count >= 0.45 * max_count else "#1F2A37"
            ax.text(j, i, label, ha="center", va="center", fontsize=10, fontweight="bold", color=color)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([x - 0.5 for x in range(1, len(GENOTYPE_ORDER))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(GENOTYPE_ORDER))], minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return out_path


def make_requested_heatmaps(genotype_counts: pd.DataFrame) -> Dict[str, Path]:
    return {
        group: make_group_genotype_heatmap(genotype_counts, group, stem, title_label)
        for group, stem, title_label in HEATMAP_SPECS
    }

def make_target_snp_status_table(target_positions: pd.DataFrame, sample_genotypes: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    total_samples = int(len(sample_genotypes))
    lookup = target_positions.set_index("rsID_clean")
    for rsid in TARGET_RSIDS:
        callable_n = int(sample_genotypes[f"{rsid}_Callable"].sum())
        noncallable_n = total_samples - callable_n
        rows.append(
            {
                "SNP": rsid,
                "POS": int(lookup.at[rsid, "POS_int"]),
                "Gene": str(lookup.at[rsid, "Gene"]),
                "REF": str(lookup.at[rsid, "REF"]),
                "ALT": str(lookup.at[rsid, "ALT"]),
                "Present_in_SNPs_Used": True,
                "Present_in_Phased_Matrix": True,
                "Callable_Samples": callable_n,
                "NonCallable_Samples": noncallable_n,
                "Included_In_Group_Summary": True,
                "Included_In_Comparison_Summary": callable_n > 0,
                "Inclusion_Note": "Explicitly targeted in TARGET_RSIDS and retained irrespective of generic ranking or threshold rules."
                if callable_n > 0
                else "Explicitly targeted in TARGET_RSIDS, but all samples were non-callable in the phased genotype matrix.",
            }
        )
    return pd.DataFrame(rows)


def make_target_snp_group_summary(sample_genotypes: pd.DataFrame) -> pd.DataFrame:
    expanded = with_aggregate_groups(sample_genotypes)
    rows: List[Dict[str, object]] = []
    for rsid in TARGET_RSIDS:
        for group in GROUP_ORDER:
            sub_all = expanded[expanded["Group"].astype(str).eq(group)].copy()
            if sub_all.empty:
                continue
            sub = sub_all[sub_all[f"{rsid}_Callable"] == 1].copy()
            wt = int((sub[f"{rsid}_Dose"] == 0).sum())
            het = int((sub[f"{rsid}_Dose"] == 1).sum())
            hom = int((sub[f"{rsid}_Dose"] == 2).sum())
            callable_n = int(len(sub))
            total_n = int(len(sub_all))
            carriers = int((sub[f"{rsid}_Dose"] > 0).sum())
            alt_alleles = int(sub[f"{rsid}_Dose"].fillna(0).sum()) if callable_n else 0
            rows.append(
                {
                    "SNP": rsid,
                    "Group": group,
                    "Group_Label": GROUP_LABELS.get(group, group),
                    "Total_Samples": total_n,
                    "Callable_Samples": callable_n,
                    "NonCallable_Samples": total_n - callable_n,
                    "WT_Count": wt,
                    "Het_Count": het,
                    "Hom_Count": hom,
                    "WT_Pct_Callable": safe_fraction(wt, callable_n),
                    "Het_Pct_Callable": safe_fraction(het, callable_n),
                    "Hom_Pct_Callable": safe_fraction(hom, callable_n),
                    "Carrier_Count": carriers,
                    "NonCarrier_Count": callable_n - carriers,
                    "Carrier_Pct_Callable": safe_fraction(carriers, callable_n),
                    "ALT_Allele_Count": alt_alleles,
                    "ALT_Allele_Denominator": int(2 * callable_n),
                    "ALT_Allele_Freq_Callable": safe_fraction(alt_alleles, 2 * callable_n),
                    "Percentages_Denominator": "Callable samples for this SNP",
                }
            )
    out = pd.DataFrame(rows)
    out["Group"] = pd.Categorical(out["Group"], categories=GROUP_ORDER, ordered=True)
    out = out.sort_values(["SNP", "Group"], kind="stable").reset_index(drop=True)
    out["Group"] = out["Group"].astype(str)
    return out


def make_forced_heatmap_concordance_tables(per_sample: pd.DataFrame) -> pd.DataFrame:
    lookup = {
        "rs11078928": "chr17:39908216:T",
        "rs869402": "chr17:39911790:T",
    }
    tissue_map = {"normal": "Healthy", "tumour": "Tumour", "tumor": "Tumour"}
    gt_map = {"WT": "WT", "Het_ALT": "Het", "Hom_ALT": "Hom"}
    keep = per_sample[per_sample["Locus_Key"].astype(str).isin(lookup.values())].copy()
    keep["Cohort"] = keep["Cohort"].astype(str).str.strip()
    keep["Tissue"] = keep["Tissue"].astype(str).str.strip().str.lower().map(tissue_map)
    keep = keep[keep["Tissue"].notna()].copy()
    keep["Group"] = keep["Cohort"] + "_" + keep["Tissue"]
    keep["Callable"] = keep["Meets_Callable_DP"].astype(bool)

    sample_rows: List[Dict[str, object]] = []
    for (sample, cohort, tissue, group), sub in keep.groupby(["Sample", "Cohort", "Tissue", "Group"], sort=False):
        row: Dict[str, object] = {
            "Sample": sample,
            "Cohort": cohort,
            "Tissue": tissue,
            "Group": group,
        }
        callable_flags = []
        for rsid, locus_key in lookup.items():
            locus_row = sub[sub["Locus_Key"].astype(str).eq(locus_key)]
            if locus_row.empty:
                row[f"{rsid}_GT"] = pd.NA
                row[f"{rsid}_Callable"] = False
                callable_flags.append(False)
                continue
            locus_row = locus_row.iloc[0]
            gt_status = str(locus_row["Genotype_Status"])
            row[f"{rsid}_GT"] = gt_map.get(gt_status, pd.NA)
            row[f"{rsid}_Callable"] = bool(locus_row["Callable"])
            callable_flags.append(bool(locus_row["Callable"]))
        row["Both_Target_SNPs_Callable"] = int(all(callable_flags))
        sample_rows.append(row)

    sample_df = pd.DataFrame(sample_rows)
    if sample_df.empty:
        return pd.DataFrame(columns=["Group", "rs11078928_GT", "rs869402_GT", "Sample_Count", "Row_Proportion_within_rs11078928_GT"])

    expanded = with_aggregate_groups(sample_df)
    count_rows: List[pd.DataFrame] = []
    for group in GROUP_ORDER:
        sub = expanded[(expanded["Group"].astype(str).eq(group)) & (expanded["Both_Target_SNPs_Callable"] == 1)].copy()
        if sub.empty:
            continue
        counts = (
            sub.groupby(["rs11078928_GT", "rs869402_GT"], as_index=False)
            .size()
            .rename(columns={"size": "Sample_Count"})
        )
        counts["Group"] = group
        counts["Row_Proportion_within_rs11078928_GT"] = counts.groupby("rs11078928_GT")["Sample_Count"].transform(
            lambda s: s / s.sum()
        )
        count_rows.append(
            counts[["Group", "rs11078928_GT", "rs869402_GT", "Sample_Count", "Row_Proportion_within_rs11078928_GT"]]
        )

    return ordered_group(pd.concat(count_rows, ignore_index=True), "Group") if count_rows else pd.DataFrame()


def make_target_snp_comparison_table(target_snp_group_summary: pd.DataFrame) -> pd.DataFrame:
    lookup = target_snp_group_summary.set_index(["SNP", "Group"])
    rows: List[Dict[str, object]] = []
    for rsid in TARGET_RSIDS:
        for cohort_key, comparison_label, normal_group, tumour_group in COMPARISON_SPECS:
            normal = lookup.loc[(rsid, normal_group)]
            tumour = lookup.loc[(rsid, tumour_group)]
            odds_ratio, p_value = fisher_exact(
                [
                    [int(tumour["Carrier_Count"]), int(tumour["NonCarrier_Count"])],
                    [int(normal["Carrier_Count"]), int(normal["NonCarrier_Count"])],
                ]
            )

            shift_rows = build_genotype_shift_rows(normal, tumour)
            genotype_deltas = {row["Genotype"]: row["Delta"] for row in shift_rows}
            dominant_desc = describe_dominant_genotype_shift(shift_rows)
            genotype_balance = describe_full_genotype_balance(shift_rows)
            callability_note = (
                "All samples were callable for this SNP in both groups."
                if int(normal["NonCallable_Samples"]) + int(tumour["NonCallable_Samples"]) == 0
                else f"Non-callable samples: {GROUP_LABELS[normal_group]} {int(normal['NonCallable_Samples'])}; {GROUP_LABELS[tumour_group]} {int(tumour['NonCallable_Samples'])}."
            )
            rows.append(
                {
                    "SNP": rsid,
                    "Comparison_Key": cohort_key,
                    "Comparison": comparison_label,
                    "Normal_Group": normal_group,
                    "Tumour_Group": tumour_group,
                    "Normal_Label": GROUP_LABELS[normal_group],
                    "Tumour_Label": GROUP_LABELS[tumour_group],
                    "Comparison_Logic": describe_comparison_logic(normal_group, tumour_group),
                    "Carrier_Test_Model": "Any ALT carrier versus WT/non-carrier among callable samples",
                    "Genotype_Frequency_Model": "WT, Het, and Hom percentages among callable samples",
                    "Normal_Total_Samples": int(normal["Total_Samples"]),
                    "Normal_Callable_Samples": int(normal["Callable_Samples"]),
                    "Normal_NonCallable_Samples": int(normal["NonCallable_Samples"]),
                    "Tumour_Total_Samples": int(tumour["Total_Samples"]),
                    "Tumour_Callable_Samples": int(tumour["Callable_Samples"]),
                    "Tumour_NonCallable_Samples": int(tumour["NonCallable_Samples"]),
                    "Normal_WT_Count": int(normal["WT_Count"]),
                    "Normal_Het_Count": int(normal["Het_Count"]),
                    "Normal_Hom_Count": int(normal["Hom_Count"]),
                    "Tumour_WT_Count": int(tumour["WT_Count"]),
                    "Tumour_Het_Count": int(tumour["Het_Count"]),
                    "Tumour_Hom_Count": int(tumour["Hom_Count"]),
                    "Normal_WT_Pct_Callable": float(normal["WT_Pct_Callable"]),
                    "Normal_Het_Pct_Callable": float(normal["Het_Pct_Callable"]),
                    "Normal_Hom_Pct_Callable": float(normal["Hom_Pct_Callable"]),
                    "Tumour_WT_Pct_Callable": float(tumour["WT_Pct_Callable"]),
                    "Tumour_Het_Pct_Callable": float(tumour["Het_Pct_Callable"]),
                    "Tumour_Hom_Pct_Callable": float(tumour["Hom_Pct_Callable"]),
                    "WT_Pct_Delta_Tumour_minus_Normal": genotype_deltas["WT"],
                    "Het_Pct_Delta_Tumour_minus_Normal": genotype_deltas["Het"],
                    "Hom_Pct_Delta_Tumour_minus_Normal": genotype_deltas["Hom"],
                    "Normal_Carrier_Count": int(normal["Carrier_Count"]),
                    "Normal_NonCarrier_Count": int(normal["NonCarrier_Count"]),
                    "Tumour_Carrier_Count": int(tumour["Carrier_Count"]),
                    "Tumour_NonCarrier_Count": int(tumour["NonCarrier_Count"]),
                    "Normal_Carrier_Pct_Callable": float(normal["Carrier_Pct_Callable"]),
                    "Tumour_Carrier_Pct_Callable": float(tumour["Carrier_Pct_Callable"]),
                    "Carrier_Pct_Delta_Tumour_minus_Normal": float(tumour["Carrier_Pct_Callable"] - normal["Carrier_Pct_Callable"]),
                    "Normal_ALT_Allele_Freq_Callable": float(normal["ALT_Allele_Freq_Callable"]),
                    "Tumour_ALT_Allele_Freq_Callable": float(tumour["ALT_Allele_Freq_Callable"]),
                    "ALT_Allele_Freq_Delta_Tumour_minus_Normal": float(tumour["ALT_Allele_Freq_Callable"] - normal["ALT_Allele_Freq_Callable"]),
                    "Carrier_Odds_Ratio_Fisher": float(odds_ratio),
                    "Carrier_P_Fisher": float(p_value),
                    "Dominant_Genotype_Shift": dominant_desc,
                    "Genotype_Balance_Summary": genotype_balance,
                    "Callability_Note": callability_note,
                    "Percentages_Denominator": "Callable samples for the relevant SNP within each group",
                }
            )
    out = pd.DataFrame(rows)
    comparison_order = [label for _, label, _, _ in COMPARISON_SPECS]
    out["Comparison"] = pd.Categorical(out["Comparison"], categories=comparison_order, ordered=True)
    out = out.sort_values(["SNP", "Comparison"], kind="stable").reset_index(drop=True)
    out["Comparison"] = out["Comparison"].astype(str)
    return out


def make_target_snp_interpretation_table(target_snp_comparisons: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for rsid in TARGET_RSIDS:
        sub = target_snp_comparisons[target_snp_comparisons["SNP"].eq(rsid)].copy()
        for record in sub.itertuples(index=False):
            carrier_saturated = min(record.Normal_Carrier_Pct_Callable, record.Tumour_Carrier_Pct_Callable) >= 0.95
            if carrier_saturated:
                lead = (
                    f"Carrier frequency was near-saturated in both groups ({record.Normal_Carrier_Count}/{record.Normal_Callable_Samples}, {record.Normal_Carrier_Pct_Callable:.1%} normal; "
                    f"{record.Tumour_Carrier_Count}/{record.Tumour_Callable_Samples}, {record.Tumour_Carrier_Pct_Callable:.1%} tumour), so genotype balance is more informative than carrier status."
                )
            else:
                if record.Carrier_Pct_Delta_Tumour_minus_Normal > 0:
                    direction = "higher"
                elif record.Carrier_Pct_Delta_Tumour_minus_Normal < 0:
                    direction = "lower"
                else:
                    direction = "similar"
                lead = (
                    f"Carrier frequency was {direction} in tumour than normal ({record.Tumour_Carrier_Count}/{record.Tumour_Callable_Samples}, {record.Tumour_Carrier_Pct_Callable:.1%} vs "
                    f"{record.Normal_Carrier_Count}/{record.Normal_Callable_Samples}, {record.Normal_Carrier_Pct_Callable:.1%})."
                )

            geno_sentence = record.Genotype_Balance_Summary + "."

            if record.Carrier_P_Fisher < 0.05:
                stat_sentence = f"Carrier Fisher exact P = {record.Carrier_P_Fisher:.3g}, so the carrier contrast reached nominal significance."
            else:
                stat_sentence = f"Carrier Fisher exact P = {record.Carrier_P_Fisher:.3g}; this carrier contrast did not reach nominal significance."

            caveat_parts = [record.Callability_Note]
            if carrier_saturated:
                caveat_parts.append("Because carrier frequency is close to ceiling, ALT allele frequency and the WT/Het/Hom split are more informative than carrier status alone.")
            caveat = " ".join(caveat_parts)
            rows.append(
                {
                    "SNP": rsid,
                    "Comparison_Key": record.Comparison_Key,
                    "Comparison": record.Comparison,
                    "Summary": " ".join([lead, geno_sentence, stat_sentence]),
                    "Caveat": caveat,
                    "Comparison_Logic": record.Comparison_Logic,
                }
            )
    out = pd.DataFrame(rows)
    comparison_order = [label for _, label, _, _ in COMPARISON_SPECS]
    out["Comparison"] = pd.Categorical(out["Comparison"], categories=comparison_order, ordered=True)
    out = out.sort_values(["SNP", "Comparison"], kind="stable").reset_index(drop=True)
    out["Comparison"] = out["Comparison"].astype(str)
    return out

def make_target_snp_case_control_figure(target_snp_group_summary: pd.DataFrame, rsid: str) -> Path:
    out_path = OUT_DIR / f"25_{rsid}_Case_Control_Comparison.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sub = target_snp_group_summary[
        target_snp_group_summary["SNP"].eq(rsid) & target_snp_group_summary["Group"].isin(REPORT_GROUP_ORDER)
    ].copy()
    sub["Group"] = pd.Categorical(sub["Group"], categories=REPORT_GROUP_ORDER, ordered=True)
    sub = sub.sort_values("Group", kind="stable").reset_index(drop=True)

    x = list(range(len(sub)))
    labels = []
    for row in sub.itertuples(index=False):
        label = f"{GROUP_SHORT_LABELS[row.Group]}\ncallable {row.Callable_Samples}/{row.Total_Samples}"
        if row.NonCallable_Samples:
            label += f"\nno-call {row.NonCallable_Samples}"
        labels.append(label)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.6), dpi=220, gridspec_kw={"width_ratios": [1.2, 1.0]})

    bottom = [0.0] * len(sub)
    genotype_colors = {"WT": "#D6E4F0", "Het": "#6FA8DC", "Hom": "#0B5394"}
    for gt_col, label in [("WT_Pct_Callable", "WT"), ("Het_Pct_Callable", "Het"), ("Hom_Pct_Callable", "Hom ALT")]:
        heights = sub[gt_col].fillna(0.0).tolist()
        ax1.bar(x, heights, bottom=bottom, color=genotype_colors[label.split()[0]], edgecolor="white", linewidth=0.7, label=label)
        bottom = [b + h for b, h in zip(bottom, heights)]
    ax1.set_xticks(x, labels)
    ax1.set_ylim(0, 1.05)
    ax1.set_ylabel("Genotype percentage among callable samples")
    ax1.set_title(f"{rsid}: Genotype distribution (WT / Het / Hom)", fontsize=11)
    ax1.legend(frameon=False, fontsize=9, loc="upper left")
    ax1.tick_params(axis="x", labelrotation=15)

    width = 0.35
    carrier = sub["Carrier_Pct_Callable"].fillna(0.0).tolist()
    alt_af = sub["ALT_Allele_Freq_Callable"].fillna(0.0).tolist()
    ax2.bar([i - width / 2 for i in x], carrier, width=width, color="#F4A261", label="Carrier frequency")
    ax2.bar([i + width / 2 for i in x], alt_af, width=width, color="#2A9D8F", label="ALT allele frequency")
    ax2.set_xticks(x, labels)
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("Frequency among callable samples")
    ax2.set_title(f"{rsid}: Carrier frequency (ALT presence) and ALT allele frequency", fontsize=11)
    ax2.legend(frameon=False, fontsize=9, loc="upper left")
    ax2.tick_params(axis="x", labelrotation=15)
    for xpos, value in zip([i - width / 2 for i in x], carrier):
        ax2.text(xpos, value + 0.02, f"{value:.0%}", ha="center", va="bottom", fontsize=8)
    for xpos, value in zip([i + width / 2 for i in x], alt_af):
        ax2.text(xpos, value + 0.02, f"{value:.0%}", ha="center", va="bottom", fontsize=8)

    fig.suptitle(
        f"{rsid}: explicit tumour-vs-normal comparison with pooled normal controls included",
        fontsize=12,
        y=1.02,
    )
    fig.text(
        0.5,
        -0.02,
        "Genotype distribution uses callable samples for the relevant SNP. Carrier frequency is ALT presence among callable samples; ALT allele frequency denominator = 2 x callable samples.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return out_path


def make_target_snp_case_control_figures(target_snp_group_summary: pd.DataFrame) -> Dict[str, Path]:
    return {rsid: make_target_snp_case_control_figure(target_snp_group_summary, rsid) for rsid in TARGET_RSIDS}


def build_results_summary_table(
    target_snp_comparisons: pd.DataFrame,
    target_snp_interpretation: pd.DataFrame,
    direct: pd.DataFrame,
    significance: pd.DataFrame,
) -> pd.DataFrame:
    comparison_lookup = target_snp_comparisons.set_index(["SNP", "Comparison"])
    interpretation_lookup = target_snp_interpretation.set_index(["SNP", "Comparison"])
    direct_lookup = direct.set_index("Comparison_Key")

    def summarise_significance(cohort: str) -> str:
        sub = significance[significance["Comparison_Key"].astype(str).eq(cohort)].copy()
        if sub.empty:
            return "No mini-haplotype significance table available."
        fdr_hits = sub[sub["Significant_FDR_BH"] == True].copy()
        nominal_hits = sub[(sub["Significant_P_Fisher"] == True) & (sub["Significant_FDR_BH"] != True)].copy()
        parts = []
        if not fdr_hits.empty:
            desc = "; ".join(
                f"{row['Mini_Haplotype_Label']} ({row['Healthy_Freq']:.1%} healthy vs {row['Tumour_Freq']:.1%} tumour; FDR {row['FDR_BH']:.3g})"
                for _, row in fdr_hits.sort_values("FDR_BH", kind="stable").iterrows()
            )
            parts.append(f"FDR-significant: {desc}")
        if not nominal_hits.empty:
            desc = "; ".join(
                f"{row['Mini_Haplotype_Label']} (P {row['P_Fisher']:.3g}, FDR {row['FDR_BH']:.3g})"
                for _, row in nominal_hits.sort_values("P_Fisher", kind="stable").iterrows()
            )
            parts.append(f"Nominal only: {desc}")
        if not parts:
            return "No mini-haplotype shifts reached nominal or FDR significance."
        return " | ".join(parts)

    rows = []
    for cohort_key, comparison_label, _, _ in COMPARISON_SPECS:
        comp_110 = comparison_lookup.loc[("rs11078928", comparison_label)]
        comp_869 = comparison_lookup.loc[("rs869402", comparison_label)]
        direct_row = direct_lookup.loc[cohort_key]
        rows.append(
            {
                "Comparison": comparison_label,
                "Comparison logic": comp_110["Comparison_Logic"],
                "rs11078928 callable normal/tumour": f"{int(comp_110['Normal_Callable_Samples'])}/{int(comp_110['Normal_Total_Samples'])} vs {int(comp_110['Tumour_Callable_Samples'])}/{int(comp_110['Tumour_Total_Samples'])}",
                "rs869402 callable normal/tumour": f"{int(comp_869['Normal_Callable_Samples'])}/{int(comp_869['Normal_Total_Samples'])} vs {int(comp_869['Tumour_Callable_Samples'])}/{int(comp_869['Tumour_Total_Samples'])}",
                "rs11078928 summary": interpretation_lookup.loc[("rs11078928", comparison_label), "Summary"],
                "rs869402 summary": interpretation_lookup.loc[("rs869402", comparison_label), "Summary"],
                "Mini-haplotype readout": (
                    f"Both_ALT {float(direct_row['Healthy_Both_ALT_Freq']):.1%}->{float(direct_row['Tumour_Both_ALT_Freq']):.1%}; "
                    f"rs11078928-only delta {float(direct_row['Delta_rs11078928_only']):+.1%}; "
                    f"rs869402-only delta {float(direct_row['Delta_rs869402_only']):+.1%}"
                ),
                "Mini-haplotype support": summarise_significance(cohort_key),
            }
        )
    return pd.DataFrame(rows)


def build_genotype_matrix_exports(genotype_counts: pd.DataFrame) -> Dict[str, Dict[str, pd.DataFrame]]:
    exports: Dict[str, Dict[str, pd.DataFrame]] = {}
    for group, _, title_label in HEATMAP_SPECS:
        sub = genotype_counts[genotype_counts["Group"].astype(str).eq(group)].copy()
        count_matrix = (
            sub.pivot(index="rs11078928_GT", columns="rs869402_GT", values="Sample_Count")
            .reindex(index=GENOTYPE_ORDER, columns=GENOTYPE_ORDER)
            .fillna(0)
            .astype(int)
        )
        count_matrix.index.name = "rs11078928 genotype"
        count_matrix.columns = [f"rs869402 {col}" for col in count_matrix.columns]

        row_pct_matrix = (
            sub.pivot(index="rs11078928_GT", columns="rs869402_GT", values="Row_Proportion_within_rs11078928_GT")
            .reindex(index=GENOTYPE_ORDER, columns=GENOTYPE_ORDER)
            .fillna(0.0)
        )
        row_pct_matrix.index.name = "rs11078928 genotype"
        row_pct_matrix.columns = [f"rs869402 {col}" for col in row_pct_matrix.columns]
        row_pct_display = row_pct_matrix.map(lambda x: f"{float(x):.1%}")

        exports[group] = {
            "title": pd.DataFrame(
                {
                    "Item": ["Focus", "Read this as"],
                    "Value": [
                        title_label,
                        "Rows are rs11078928 genotype; columns are rs869402 genotype. Row percentages stay within the rs11078928 genotype row.",
                    ],
                }
            ),
            "counts": count_matrix.reset_index(),
            "row_pct": row_pct_display.reset_index(),
        }
    return exports

def make_full_region_tables(workbook: pd.ExcelFile) -> Tuple[pd.DataFrame, pd.DataFrame]:
    full_freq = pd.read_excel(workbook, sheet_name="Full_Global_Freqs")
    keep = full_freq[(full_freq["rs11078928"] == 2) | (full_freq["rs869402"] == 2)].copy()
    keep["Pattern_Label"] = keep.apply(
        lambda r: (
            "Both_ALT"
            if r["rs11078928"] == 2 and r["rs869402"] == 2
            else "rs11078928_only"
            if r["rs11078928"] == 2
            else "rs869402_only"
        ),
        axis=1,
    )
    keep = keep[
        [
            "Haplotype_ID",
            "Pattern_Label",
            "Frequency",
            "Count",
            "rs11078928",
            "rs869402",
            "Variant_Content",
        ]
    ].sort_values(["Pattern_Label", "Frequency"], ascending=[True, False], kind="stable")

    full_assoc = pd.read_excel(workbook, sheet_name="Full_Associations").copy()
    full_assoc["has_rs11078928"] = full_assoc["Haplotype_Def"].astype(str).str.contains("rs11078928", regex=False)
    full_assoc["has_rs869402"] = full_assoc["Haplotype_Def"].astype(str).str.contains("rs869402", regex=False)
    full_assoc = full_assoc[full_assoc["has_rs11078928"] | full_assoc["has_rs869402"]].copy()
    full_assoc["Pattern_Label"] = full_assoc.apply(
        lambda r: (
            "Both_ALT"
            if r["has_rs11078928"] and r["has_rs869402"]
            else "rs11078928_only"
            if r["has_rs11078928"]
            else "rs869402_only"
        ),
        axis=1,
    )
    summary = (
        full_assoc.groupby(["Comparison", "Pattern_Label"], as_index=False)[["Freq_Tumour", "Freq_Healthy"]]
        .sum()
        .sort_values(["Comparison", "Pattern_Label"], kind="stable")
    )
    summary["Tumour_minus_Healthy"] = summary["Freq_Tumour"] - summary["Freq_Healthy"]
    return keep.reset_index(drop=True), summary.reset_index(drop=True)


def load_ld_block7_tables(workbook: pd.ExcelFile) -> Tuple[pd.DataFrame, pd.DataFrame]:
    block_freqs = pd.read_excel(workbook, sheet_name="LD_Block_r2_080_Block7_Freqs").copy()
    block_assoc = pd.read_excel(workbook, sheet_name="LD_Block_r2_080_Block7_Assoc").copy()
    return block_freqs, block_assoc


def make_direct_comparison_table(mini_freqs: pd.DataFrame) -> pd.DataFrame:
    subset = mini_freqs[
        mini_freqs["Group"].isin(
            ["Breast_Healthy", "Breast_Tumour", "Endometrium_Healthy", "Endometrium_Tumour", "Healthy_All", "Tumour_All"]
        )
    ].copy()
    freq_pivot = subset.pivot_table(
        index="Group", columns="Mini_Haplotype_Label", values="Haplotype_Copy_Frequency", aggfunc="first", observed=False
    ).fillna(0.0)
    count_pivot = subset.pivot_table(
        index="Group", columns="Mini_Haplotype_Label", values="Haplotype_Copy_Count", aggfunc="sum", observed=False
    ).fillna(0)
    rows = []
    for cohort, comparison_label, healthy_group, tumour_group in COMPARISON_SPECS:
        h = freq_pivot.loc[healthy_group]
        t = freq_pivot.loc[tumour_group]
        healthy_total = int(count_pivot.loc[healthy_group].sum())
        tumour_total = int(count_pivot.loc[tumour_group].sum())
        delta_110 = float(t.get("rs11078928_only", 0.0) - h.get("rs11078928_only", 0.0))
        delta_869 = float(t.get("rs869402_only", 0.0) - h.get("rs869402_only", 0.0))
        delta_both = float(t.get("Both_ALT", 0.0) - h.get("Both_ALT", 0.0))
        inverse_supported = (delta_110 * delta_869) < 0 and abs(delta_both) < max(abs(delta_110), abs(delta_869))
        rows.append(
            {
                "Comparison_Key": cohort,
                "Comparison": comparison_label,
                "Healthy_Label": GROUP_LABELS[healthy_group],
                "Tumour_Label": GROUP_LABELS[tumour_group],
                "Comparison_Logic": describe_comparison_logic(healthy_group, tumour_group),
                "Healthy_Callable_Samples": int(healthy_total / 2),
                "Tumour_Callable_Samples": int(tumour_total / 2),
                "Healthy_Haplotype_Copy_Total": healthy_total,
                "Tumour_Haplotype_Copy_Total": tumour_total,
                "Healthy_rs11078928_only_Freq": float(h.get("rs11078928_only", 0.0)),
                "Tumour_rs11078928_only_Freq": float(t.get("rs11078928_only", 0.0)),
                "Healthy_rs869402_only_Freq": float(h.get("rs869402_only", 0.0)),
                "Tumour_rs869402_only_Freq": float(t.get("rs869402_only", 0.0)),
                "Healthy_Both_ALT_Freq": float(h.get("Both_ALT", 0.0)),
                "Tumour_Both_ALT_Freq": float(t.get("Both_ALT", 0.0)),
                "Delta_rs11078928_only": delta_110,
                "Delta_rs869402_only": delta_869,
                "Delta_Both_ALT": delta_both,
                "Inverse_Pattern_Supported": inverse_supported,
                "Frequency_Definition": (
                    "Each frequency is the proportion of callable phased haplotype copies in that arm; "
                    "Breast and endometrium rows compare tumours against pooled normal controls, while the pooled row compares pooled tumours with pooled normal controls."
                ),
            }
        )
    return pd.DataFrame(rows)


def make_mini_haplotype_significance_table(mini_freqs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort, comparison_label, healthy_group, tumour_group in COMPARISON_SPECS:
        sub = mini_freqs[mini_freqs["Group"].isin([healthy_group, tumour_group])].copy()
        healthy_total = int(sub.loc[sub["Group"].eq(healthy_group), "Haplotype_Copy_Count"].sum())
        tumour_total = int(sub.loc[sub["Group"].eq(tumour_group), "Haplotype_Copy_Count"].sum())
        lookup = {
            (str(group), str(label)): int(count)
            for group, label, count in sub[["Group", "Mini_Haplotype_Label", "Haplotype_Copy_Count"]].itertuples(index=False)
        }
        for label in ["Both_ALT", "rs869402_only", "rs11078928_only", "Neither_ALT"]:
            healthy_count = lookup.get((healthy_group, label), 0)
            tumour_count = lookup.get((tumour_group, label), 0)
            odds_ratio, p_value = fisher_exact(
                [[tumour_count, tumour_total - tumour_count], [healthy_count, healthy_total - healthy_count]]
            )
            rows.append(
                {
                    "Comparison_Key": cohort,
                    "Comparison": comparison_label,
                    "Healthy_Label": GROUP_LABELS[healthy_group],
                    "Tumour_Label": GROUP_LABELS[tumour_group],
                    "Mini_Haplotype_Label": label,
                    "Tumour_Count": tumour_count,
                    "Tumour_Total": tumour_total,
                    "Healthy_Count": healthy_count,
                    "Healthy_Total": healthy_total,
                    "Tumour_Freq": safe_fraction(tumour_count, tumour_total),
                    "Healthy_Freq": safe_fraction(healthy_count, healthy_total),
                    "Odds_Ratio": odds_ratio,
                    "P_Fisher": p_value,
                    "Comparison_Logic": describe_comparison_logic(healthy_group, tumour_group),
                    "Frequency_Definition": (
                        "Frequencies are based on callable phased haplotype copies; each callable sample contributes two copies."
                    ),
                }
            )
    out = pd.DataFrame(rows)
    adjusted = []
    for _, sub in out.groupby("Comparison_Key", sort=False):
        sub = sub.sort_values("P_Fisher", kind="stable").reset_index(drop=True)
        m = len(sub)
        bh = [min(float(sub.loc[i, "P_Fisher"]) * m / (i + 1), 1.0) for i in range(m)]
        for i in range(m - 2, -1, -1):
            bh[i] = min(bh[i], bh[i + 1])
        sub["FDR_BH"] = bh
        sub["Significant_P_Fisher"] = sub["P_Fisher"] < 0.05
        sub["Significant_FDR_BH"] = sub["FDR_BH"] < 0.05
        adjusted.append(sub)
    return pd.concat(adjusted, ignore_index=True)


def build_interpretation(
    target_status: pd.DataFrame,
    target_snp_interpretation: pd.DataFrame,
    direct: pd.DataFrame,
    ld_context: pd.DataFrame,
    full_relevant: pd.DataFrame,
) -> str:
    status_lookup = target_status.set_index("SNP")
    interpretation_lookup = target_snp_interpretation.set_index(["SNP", "Comparison"])
    direct_lookup = direct.set_index("Comparison_Key")
    rs869402_ld = ld_context.loc[ld_context["rsID"].eq("rs869402"), "LD_Block_Status"].iloc[0]
    global_both = full_relevant.loc[full_relevant["Pattern_Label"].eq("Both_ALT"), "Frequency"].sum()
    global_110_only = full_relevant.loc[full_relevant["Pattern_Label"].eq("rs11078928_only"), "Frequency"].sum()
    global_869_only = full_relevant.loc[full_relevant["Pattern_Label"].eq("rs869402_only"), "Frequency"].sum()

    lines = [
        "Focused phased-haplotype and case-control interpretation for rs11078928 and rs869402",
        "",
        "Target inclusion and genotype source",
    ]
    for rsid in TARGET_RSIDS:
        status = status_lookup.loc[rsid]
        lines.append(
            f"- {rsid}: included explicitly; {int(status['Callable_Samples'])}/{int(status['Callable_Samples']) + int(status['NonCallable_Samples'])} samples callable; {status['Inclusion_Note']}"
        )
    lines.extend(
        [
            "- WT, Het, and Hom in this report come directly from the forced-genotype workbook (analysis_results/05b_forced_genotypes/GSDMB_Forced_Genotypes_Union_Sites.xlsx); they are not inferred from the phased haplotype matrix or from carrier subtraction.",
            "- Pooled controls refer to Healthy_All (breast normal + endometrium normal), and pooled tumours refer to Tumour_All.",
            "",
            "Per-comparison SNP summaries",
        ]
    )

    for _, comparison_label, _, _ in COMPARISON_SPECS:
        lines.append(f"- {comparison_label} | rs11078928: {interpretation_lookup.loc[('rs11078928', comparison_label), 'Summary']}")
        lines.append(f"- {comparison_label} | rs869402: {interpretation_lookup.loc[('rs869402', comparison_label), 'Summary']}")

    breast = direct_lookup.loc["Breast"]
    endo = direct_lookup.loc["Endometrium"]
    pooled = direct_lookup.loc["All"]
    lines.extend(
        [
            "",
            "Joint haplotype context",
            (
                f"- Across the full-region phased backbone, the dominant relevant background carried both ALT alleles (global frequency {global_both:.3f}), whereas rs11078928-only ({global_110_only:.3f}) and rs869402-only ({global_869_only:.3f}) backgrounds were less common."
            ),
            (
                f"- In breast, the Both_ALT mini-haplotype increased from {breast['Healthy_Both_ALT_Freq']:.1%} in pooled normal controls to {breast['Tumour_Both_ALT_Freq']:.1%} in tumours; in endometrium the increase was smaller ({endo['Healthy_Both_ALT_Freq']:.1%} to {endo['Tumour_Both_ALT_Freq']:.1%})."
            ),
            (
                f"- In pooled controls versus pooled tumours, Both_ALT rose from {pooled['Healthy_Both_ALT_Freq']:.1%} to {pooled['Tumour_Both_ALT_Freq']:.1%}, while both single-variant backgrounds were lower in tumours."
            ),
            (
                f"- rs869402 is {str(rs869402_ld).lower()}, so the joint relationship between rs11078928 and rs869402 is clearer in the full-region phased analysis than in the explicit LD-block framework."
            ),
        ]
    )
    return "\n".join(lines) + "\n"

def build_supervisor_summary_table(results_summary: pd.DataFrame) -> pd.DataFrame:
    """Return a clean, supervisor-facing tabular summary without note rows."""
    return results_summary.rename(
        columns={
            "Comparison logic": "Comparison Logic",
            "rs11078928 callable normal/tumour": "rs11078928 Callable Normal/Tumour",
            "rs869402 callable normal/tumour": "rs869402 Callable Normal/Tumour",
            "rs11078928 summary": "rs11078928 Summary",
            "rs869402 summary": "rs869402 Summary",
            "Mini-haplotype readout": "Mini-Haplotype Readout",
            "Mini-haplotype support": "Mini-Haplotype Support",
        }
    ).copy()


def polish_workbook(workbook_path: Path, heatmap_paths: Dict[str, Path], snp_figure_paths: Dict[str, Path]) -> None:
    wb = load_workbook(workbook_path)
    header_fill = PatternFill(fill_type="solid", fgColor="1F7C98")
    header_font = Font(color="FFFFFF", bold=True)
    title_font = Font(bold=True, size=12)
    wrap_alignment = Alignment(wrap_text=True, vertical="top")
    strong_fill = PatternFill(fill_type="solid", fgColor="E2F0D9")
    nominal_fill = PatternFill(fill_type="solid", fgColor="FFF2CC")
    caution_fill = PatternFill(fill_type="solid", fgColor="FCE4D6")

    def style_header_row(ws, row_idx: int) -> None:
        if row_idx > ws.max_row:
            return
        for cell in ws[row_idx]:
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

    heatmap_sheet = wb.create_sheet("Heatmaps", 2)
    heatmap_sheet["A1"] = "Genotype heatmaps for interpretation"
    heatmap_sheet["A1"].font = Font(bold=True, size=14)
    heatmap_sheet["A2"] = "Rows = rs11078928 genotype; columns = rs869402 genotype; labels show callable sample counts and row percentages within each rs11078928 genotype row."
    heatmap_sheet["A2"].alignment = wrap_alignment
    row_cursor = 4
    for group, _, title_label in HEATMAP_SPECS:
        heatmap_sheet[f"A{row_cursor}"] = title_label
        heatmap_sheet[f"A{row_cursor}"].font = title_font
        img = XLImage(str(heatmap_paths[group]))
        img.width = 620
        img.height = 455
        heatmap_sheet.add_image(img, f"A{row_cursor + 1}")
        row_cursor += 24

    if "README" in wb.sheetnames:
        ws = wb["README"]
        ws["A1"] = "Workbook guide for interpretation"
        ws["A1"].font = Font(bold=True, size=12)
        ws["A2"] = "Begin with Results_Summary, then review rs11078928 and rs869402 sheets, SNP_Interpretation, and the heatmaps. Genotype distribution (WT / Het / Hom) is shown directly from phased genotypes, while carrier frequency (ALT presence) is reported separately."
        ws["A2"].alignment = wrap_alignment

    if "Results_Summary" in wb.sheetnames:
        ws = wb["Results_Summary"]
        ws["A1"] = "Results summary for rs11078928 and rs869402"
        ws["A1"].font = Font(bold=True, size=12)
        ws["A2"] = "Breast and endometrium rows compare tumours against pooled normal controls. The pooled row compares pooled tumours with pooled normal controls across both cohorts."
        ws["A3"] = "Use the SNP-specific sheets for genotype distribution (WT / Het / Hom), carrier frequency (ALT presence), ALT allele frequency, and callable denominators."

    if "MiniHap_Signif" in wb.sheetnames:
        ws = wb["MiniHap_Signif"]
        ws["A1"] = "Mini-haplotype significance results"
        ws["A1"].font = Font(bold=True, size=12)
        ws["A2"] = "Green rows remained significant after Fisher exact testing with BH-FDR correction on callable phased haplotype copies."
        ws["A3"] = "Yellow rows reached nominal significance only and should be treated as supportive rather than definitive. The pooled row is pooled tumours versus pooled normal controls; it is not a breast-versus-endometrium comparison."

    if "SNP_Interpretation" in wb.sheetnames:
        ws = wb["SNP_Interpretation"]
        ws["A1"] = "SNP-specific interpretation notes"
        ws["A1"].font = Font(bold=True, size=12)
        ws["A2"] = "These summaries are comparison-focused and explicitly separate genotype distribution (WT / Het / Hom) from carrier frequency (ALT presence). Read genotype balance first when the carrier rate is close to saturation."

    if "Mini_Haplotype_Freqs" in wb.sheetnames:
        ws = wb["Mini_Haplotype_Freqs"]
        ws["A1"] = "Mini-haplotype frequencies"
        ws["A1"].font = Font(bold=True, size=12)
        ws["A2"] = "Each frequency is the proportion of callable phased haplotype copies with that two-SNP pattern in the stated group. Only samples with both target SNPs phased are included, and each callable sample contributes two haplotype copies."
        ws["A2"].alignment = wrap_alignment

    if "Direct_Comparison" in wb.sheetnames:
        ws = wb["Direct_Comparison"]
        ws["A1"] = "Direct mini-haplotype comparison summary"
        ws["A1"].font = Font(bold=True, size=12)
        ws["A2"] = "Breast and endometrium rows compare tumours against pooled normal controls. The pooled row is pooled tumours versus pooled normal controls across both cohorts."
        ws["A2"].alignment = wrap_alignment

    for rsid in TARGET_RSIDS:
        if rsid in wb.sheetnames:
            ws = wb[rsid]
            ws["A1"] = f"{rsid}: Genotype distribution (WT / Het / Hom) and carrier frequency (ALT presence)"
            ws["A1"].font = Font(bold=True, size=12)
            ws["A2"] = "Upper table: genotype distribution (WT / Het / Hom) by group, plus carrier frequency (ALT presence) and ALT allele frequency. Lower table: breast, endometrium, and pooled tumour-vs-normal comparison summary with callable denominators and explicit arm labels."
            ws["A11"] = "Case-control comparison summary"
            ws["A11"].font = title_font
            img = XLImage(str(snp_figure_paths[rsid]))
            img.width = 760
            img.height = 300
            ws.add_image(img, "Y2")

    for group, _, title_label in HEATMAP_SPECS:
        sheet_name = MATRIX_SHEET_MAP[group]
        if sheet_name in wb.sheetnames:
            ws = wb[sheet_name]
            ws["A1"] = f"{title_label}: sample-count matrix"
            ws["A1"].font = Font(bold=True, size=12)
            ws["A2"] = "Rows = rs11078928 genotype; columns = rs869402 genotype."
            ws["A9"] = f"{title_label}: row percentages within rs11078928 genotype"
            ws["A9"].font = Font(bold=True, size=12)

    header_rows_map = {
        "README": [3],
        "Results_Summary": [4],
        "MiniHap_Signif": [4],
        "SNP_Interpretation": [3],
    }
    for rsid in TARGET_RSIDS:
        header_rows_map[rsid] = [3, 13]
    for sheet_name in MATRIX_SHEET_MAP.values():
        header_rows_map[sheet_name] = [3, 11]

    for ws in wb.worksheets:
        for col_idx in range(1, ws.max_column + 1):
            width = 0
            col_letter = get_column_letter(col_idx)
            for cell in ws[col_letter]:
                value = "" if cell.value is None else str(cell.value)
                width = max(width, min(len(value) + 2, 55))
                if cell.row > 1:
                    cell.alignment = wrap_alignment
            ws.column_dimensions[col_letter].width = max(12, width)

        for row_idx in header_rows_map.get(ws.title, [1] if ws.title != "Heatmaps" else []):
            style_header_row(ws, row_idx)

        if ws.title == "README":
            ws.freeze_panes = "A4"
        elif ws.title in {"Results_Summary", "MiniHap_Signif"}:
            ws.freeze_panes = "A5"
        elif ws.title in TARGET_RSIDS or ws.title in set(MATRIX_SHEET_MAP.values()) or ws.title == "SNP_Interpretation":
            ws.freeze_panes = "A4"
        elif ws.title != "Heatmaps":
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions

    if "MiniHap_Signif" in wb.sheetnames:
        ws = wb["MiniHap_Signif"]
        headers = {str(cell.value): cell.column for cell in ws[4]}
        fdr_col = headers.get("Significant_FDR_BH")
        nominal_col = headers.get("Significant_P_Fisher")
        p_col = headers.get("P_Fisher")
        fdrbh_col = headers.get("FDR_BH")
        for row_idx in range(5, ws.max_row + 1):
            fdr_value = ws.cell(row=row_idx, column=fdr_col).value if fdr_col else None
            nominal_value = ws.cell(row=row_idx, column=nominal_col).value if nominal_col else None
            if fdr_value is True:
                for col_idx in range(1, ws.max_column + 1):
                    ws.cell(row=row_idx, column=col_idx).fill = strong_fill
            elif nominal_value is True:
                for col_idx in range(1, ws.max_column + 1):
                    ws.cell(row=row_idx, column=col_idx).fill = nominal_fill
        for col_idx in [p_col, fdrbh_col]:
            if col_idx:
                for row_idx in range(5, ws.max_row + 1):
                    ws.cell(row=row_idx, column=col_idx).number_format = "0.000E+00"

    if "Results_Summary" in wb.sheetnames:
        ws = wb["Results_Summary"]
        headers = {str(cell.value): cell.column for cell in ws[4]}
        stat_col = headers.get("Mini-haplotype support")
        if stat_col:
            for row_idx in range(5, ws.max_row + 1):
                text_value = str(ws.cell(row=row_idx, column=stat_col).value or "")
                if "FDR-significant" in text_value:
                    ws.cell(row=row_idx, column=stat_col).fill = strong_fill
                elif "Nominal only" in text_value:
                    ws.cell(row=row_idx, column=stat_col).fill = nominal_fill
                else:
                    ws.cell(row=row_idx, column=stat_col).fill = caution_fill

    wb.save(workbook_path)


def write_outputs(
    target_positions: pd.DataFrame,
    target_status: pd.DataFrame,
    ld_context: pd.DataFrame,
    mini_freqs: pd.DataFrame,
    sample_carriers: pd.DataFrame,
    sample_patterns: pd.DataFrame,
    genotype_counts: pd.DataFrame,
    genotype_summary: pd.DataFrame,
    heatmap_counts: pd.DataFrame,
    target_snp_group_summary: pd.DataFrame,
    target_snp_comparisons: pd.DataFrame,
    target_snp_interpretation: pd.DataFrame,
    full_relevant: pd.DataFrame,
    full_summary: pd.DataFrame,
    block7_freqs: pd.DataFrame,
    block7_assoc: pd.DataFrame,
    direct: pd.DataFrame,
    significance: pd.DataFrame,
    interpretation: str,
) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    heatmap_paths = make_requested_heatmaps(heatmap_counts)
    snp_figure_paths = make_target_snp_case_control_figures(target_snp_group_summary)
    results_summary = build_results_summary_table(target_snp_comparisons, target_snp_interpretation, direct, significance)
    supervisor_summary = build_supervisor_summary_table(results_summary)
    matrix_exports = build_genotype_matrix_exports(heatmap_counts)

    tsvs = {
        "25_Target_SNP_Metadata.tsv": target_positions,
        "25_Target_SNP_Status.tsv": target_status,
        "25_Target_SNP_Group_Summary.tsv": target_snp_group_summary,
        "25_Target_SNP_Comparison_Summary.tsv": target_snp_comparisons,
        "25_Target_SNP_Interpretation.tsv": target_snp_interpretation,
        "25_LD_Context.tsv": ld_context,
        "25_Mini_Haplotype_Frequencies.tsv": mini_freqs,
        "25_Sample_Carrier_Summary.tsv": sample_carriers,
        "25_Sample_Haplotype_Patterns.tsv": sample_patterns,
        "25_Genotype_Concordance_Counts.tsv": genotype_counts,
        "25_Genotype_Conditional_Summary.tsv": genotype_summary,
        "25_Results_Summary.tsv": results_summary,
        "25_Supervisor_Summary.tsv": supervisor_summary,
        "25_Full_Region_Relevant_Haplotypes.tsv": full_relevant,
        "25_Full_Region_Comparison_Summary.tsv": full_summary,
        "25_LD_Block7_Haplotypes.tsv": block7_freqs,
        "25_LD_Block7_Associations.tsv": block7_assoc,
        "25_Direct_Comparison_Summary.tsv": direct,
        "25_Mini_Haplotype_Significance.tsv": significance,
    }
    for group, exports in matrix_exports.items():
        tsvs[f"25_{group}_Genotype_Count_Matrix.tsv"] = exports["counts"]
        tsvs[f"25_{group}_Genotype_RowPct_Matrix.tsv"] = exports["row_pct"]
    for name, df in tsvs.items():
        df.to_csv(OUT_DIR / name, sep="\t", index=False)

    (OUT_DIR / "25_Interpretation.txt").write_text(interpretation, encoding="utf-8")

    workbook_path = OUT_DIR / "25_rs11078928_rs869402_haplotype_focus.xlsx"
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        pd.DataFrame(
            {
                "Item": [
                    "Analytical question",
                    "Recommended reading order",
                    "Primary genotype source",
                    "Primary haplotype workbook",
                    "Genotype denominator rule",
                    "Mini-haplotype denominator rule",
                    "What pooled means here",
                    "Target inclusion rule",
                    "Key limitation",
                ],
                "Value": [
                    "Focused phased-haplotype relationship plus explicit tumour-vs-normal SNP comparison for rs11078928 and rs869402.",
                    "Supervisor Summary, Results_Summary, rs11078928, rs869402, SNP_Interpretation, MiniHap_Signif, then Heatmaps.",
                    str(FORCED_GENOTYPE_XLSX),
                    str(HAPLO_STATS_XLSX),
                    "WT, Het, Hom counts use callable samples for the relevant SNP; carrier counts are ALT presence among callable samples.",
                    "Mini-haplotype frequencies use callable phased haplotype copies as the denominator, equivalent to 2 x callable samples when both target SNPs are phased.",
                    "Pooled normal controls = Healthy_All; pooled tumours = Tumour_All.",
                    "rs11078928 and rs869402 are hard-coded targets in TARGET_RSIDS and are not filtered out by ranking or generic threshold rules.",
                    "rs869402 is not part of any explicit r2=0.80 LD block in the current stage-15 workbook, so joint LD-block interpretation remains limited.",
                ],
            }
        ).to_excel(writer, sheet_name="README", index=False, startrow=2)
        supervisor_summary.to_excel(writer, sheet_name="Supervisor Summary", index=False)
        results_summary.to_excel(writer, sheet_name="Results_Summary", index=False, startrow=3)
        target_snp_interpretation.to_excel(writer, sheet_name="SNP_Interpretation", index=False, startrow=2)
        for rsid in TARGET_RSIDS:
            target_snp_group_summary[
                target_snp_group_summary["SNP"].eq(rsid) & target_snp_group_summary["Group"].isin(REPORT_GROUP_ORDER)
            ].to_excel(writer, sheet_name=rsid, index=False, startrow=2)
            target_snp_comparisons[target_snp_comparisons["SNP"].eq(rsid)].to_excel(
                writer, sheet_name=rsid, index=False, startrow=12
            )
        for group, exports in matrix_exports.items():
            sheet_name = MATRIX_SHEET_MAP[group]
            exports["counts"].to_excel(writer, sheet_name=sheet_name, index=False, startrow=2)
            exports["row_pct"].to_excel(writer, sheet_name=sheet_name, index=False, startrow=10)
        target_positions.to_excel(writer, sheet_name="Target_SNPs", index=False)
        target_status.to_excel(writer, sheet_name="Target_SNP_Status", index=False)
        ld_context.to_excel(writer, sheet_name="LD_Context", index=False)
        mini_freqs.to_excel(writer, sheet_name="Mini_Haplotype_Freqs", index=False)
        sample_carriers.to_excel(writer, sheet_name="Carrier_Summary", index=False)
        sample_patterns.to_excel(writer, sheet_name="Sample_Patterns", index=False)
        genotype_counts.to_excel(writer, sheet_name="Genotype_Concordance", index=False)
        genotype_summary.to_excel(writer, sheet_name="Genotype_Conditional", index=False)
        full_relevant.to_excel(writer, sheet_name="Full_Region_Haps", index=False)
        full_summary.to_excel(writer, sheet_name="Full_Region_Summary", index=False)
        block7_freqs.to_excel(writer, sheet_name="LD_Block7_Haps", index=False)
        block7_assoc.to_excel(writer, sheet_name="LD_Block7_Assoc", index=False)
        direct.to_excel(writer, sheet_name="Direct_Comparison", index=False)
        significance.to_excel(writer, sheet_name="MiniHap_Signif", index=False, startrow=3)
        pd.DataFrame({"Interpretation": interpretation.splitlines()}).to_excel(writer, sheet_name="Interpretation", index=False)

    polish_workbook(workbook_path, heatmap_paths, snp_figure_paths)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Focused phased-haplotype analysis for rs11078928 and rs869402")
    return parser.parse_args()


def main() -> None:
    parse_args()
    workbook, phased, sample_metadata = load_required_inputs()
    forced_workbook, forced_per_sample, forced_locus_summary = load_forced_genotype_inputs()
    target_positions = get_target_positions(workbook)
    ld_context, _ = get_ld_context(workbook)
    sample_genotypes = extract_sample_genotypes(phased, sample_metadata, target_positions)
    hap_copies = extract_haplotype_copies(sample_genotypes)
    mini_freqs = make_mini_haplotype_frequency_table(hap_copies)
    sample_carriers, sample_patterns, expanded_samples = make_sample_level_tables(sample_genotypes, hap_copies)
    genotype_counts, genotype_summary = make_genotype_concordance_tables(expanded_samples)
    target_status = make_forced_target_snp_status_table(target_positions, forced_locus_summary)
    target_snp_group_summary = make_forced_target_snp_group_summary(forced_per_sample)
    heatmap_counts = make_forced_heatmap_concordance_tables(forced_per_sample)
    target_snp_comparisons = make_target_snp_comparison_table(target_snp_group_summary)
    target_snp_interpretation = make_target_snp_interpretation_table(target_snp_comparisons)
    full_relevant, full_summary = make_full_region_tables(workbook)
    block7_freqs, block7_assoc = load_ld_block7_tables(workbook)
    direct = make_direct_comparison_table(mini_freqs)
    significance = make_mini_haplotype_significance_table(mini_freqs)
    interpretation = build_interpretation(target_status, target_snp_interpretation, direct, ld_context, full_relevant)
    write_outputs(
        target_positions=target_positions,
        target_status=target_status,
        ld_context=ld_context,
        mini_freqs=mini_freqs,
        sample_carriers=sample_carriers,
        sample_patterns=sample_patterns,
        genotype_counts=genotype_counts,
        genotype_summary=genotype_summary,
        heatmap_counts=heatmap_counts,
        target_snp_group_summary=target_snp_group_summary,
        target_snp_comparisons=target_snp_comparisons,
        target_snp_interpretation=target_snp_interpretation,
        full_relevant=full_relevant,
        full_summary=full_summary,
        block7_freqs=block7_freqs,
        block7_assoc=block7_assoc,
        direct=direct,
        significance=significance,
        interpretation=interpretation,
    )
    print("Focused haplotype analysis complete.")
    print(f"Output directory: {OUT_DIR}")


if __name__ == "__main__":
    main()
