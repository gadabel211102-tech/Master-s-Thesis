#!/usr/bin/env python3
"""
Focused phased-haplotype analysis for rs11078928 and rs869402.

This script reuses the existing stage-15 phased outputs to answer a specific
haplotype question:

- which phased haplotypes carry rs11078928 and/or rs869402
- whether the two variants occupy alternative haplotypic backgrounds
- whether any inverse pattern is visible across breast and endometrial cohorts
- whether the signal is clearer in full-region or LD-block definitions

Outputs are written to:
  analysis_results/25_rs11078928_rs869402_haplotype_focus/
"""

from __future__ import annotations

import argparse
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
OUT_DIR = RESULTS / "25_rs11078928_rs869402_haplotype_focus"

TARGET_RSIDS = ("rs11078928", "rs869402")
GROUP_ORDER = [
    "Breast_Healthy",
    "Breast_Tumour",
    "Endometrium_Healthy",
    "Endometrium_Tumour",
    "Healthy_All",
    "Tumour_All",
    "All_Samples",
]
MINI_LABELS = {
    "00": "Neither_ALT",
    "01": "rs869402_only",
    "10": "rs11078928_only",
    "11": "Both_ALT",
}
HEATMAP_SPECS = [
    ("Breast_Tumour", "Breast_Tumour", "Breast tumours"),
    ("Endometrium_Tumour", "Endometrium_Tumour", "Endometrial tumours"),
    ("Tumour_All", "Tumour_All", "All tumours pooled"),
]
DIRECT_COHORT_MAP = {
    "Breast_Tumour": "Breast",
    "Endometrium_Tumour": "Endometrium",
    "Tumour_All": "All",
}
MATRIX_SHEET_MAP = {
    "Breast_Tumour": "Breast_GT_Matrix",
    "Endometrium_Tumour": "Endometrium_GT_Matrix",
    "Tumour_All": "Global_GT_Matrix",
}


def ordered_group(df: pd.DataFrame, column: str) -> pd.DataFrame:
    out = df.copy()
    out[column] = pd.Categorical(out[column], categories=GROUP_ORDER, ordered=True)
    return out.sort_values(column, kind="stable").reset_index(drop=True)


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
                "Joint_Assessment_Within_LD_Block": "Feasible" if all(r in set(block_members["rsID"].astype(str)) for r in TARGET_RSIDS) else "Not feasible",
                "Notes": "Joint LD-block analysis requires both target SNPs to be in the same explicit block.",
            }
        )
    return pd.DataFrame(rows), blocks


def extract_haplotype_copies(phased: pd.DataFrame, sample_metadata: pd.DataFrame, target_positions: pd.DataFrame) -> pd.DataFrame:
    target_pos = dict(zip(target_positions["rsID_clean"], target_positions["POS_int"]))
    phased_subset = phased[phased["POS"].isin(target_pos.values())].copy()
    if phased_subset["POS"].nunique() != 2:
        raise ValueError("Could not recover exactly two phased rows for the target SNPs.")
    phased_subset = phased_subset.set_index("POS")
    sample_cols = [c for c in phased_subset.columns if c not in {"CHROM", "ID", "REF", "ALT"}]
    rows: List[Dict[str, object]] = []
    pos_110 = int(target_pos["rs11078928"])
    pos_869 = int(target_pos["rs869402"])
    for sample in sample_cols:
        gt_110 = str(phased_subset.at[pos_110, sample])
        gt_869 = str(phased_subset.at[pos_869, sample])
        if "|" not in gt_110 or "|" not in gt_869:
            continue
        bits_110 = gt_110.split("|")
        bits_869 = gt_869.split("|")
        for copy_idx, (bit_110, bit_869) in enumerate(zip(bits_110, bits_869), start=1):
            mini_code = f"{bit_110}{bit_869}"
            rows.append(
                {
                    "Sample": sample,
                    "Haplotype_Copy": copy_idx,
                    "rs11078928_ALT": int(bit_110),
                    "rs869402_ALT": int(bit_869),
                    "Mini_Haplotype_Code": mini_code,
                    "Mini_Haplotype_Label": MINI_LABELS[mini_code],
                }
            )
    out = pd.DataFrame(rows).merge(sample_metadata, on="Sample", how="left")
    out["Group"] = out["Group"].fillna("Unknown")
    return out


def with_aggregate_groups(hap_copies: pd.DataFrame) -> pd.DataFrame:
    extra_frames = [
        hap_copies[hap_copies["Tissue"].astype(str).eq("Healthy")].assign(Group="Healthy_All"),
        hap_copies[hap_copies["Tissue"].astype(str).eq("Tumour")].assign(Group="Tumour_All"),
        hap_copies.assign(Group="All_Samples"),
    ]
    return pd.concat([hap_copies] + extra_frames, ignore_index=True)


def make_mini_haplotype_frequency_table(hap_copies: pd.DataFrame) -> pd.DataFrame:
    expanded = with_aggregate_groups(hap_copies)
    freq = (
        expanded.groupby(["Group", "Mini_Haplotype_Code", "Mini_Haplotype_Label"], as_index=False)
        .size()
        .rename(columns={"size": "Haplotype_Copy_Count"})
    )
    totals = freq.groupby("Group")["Haplotype_Copy_Count"].transform("sum")
    freq["Haplotype_Copy_Frequency"] = freq["Haplotype_Copy_Count"] / totals
    return ordered_group(freq, "Group")


def make_sample_level_tables(hap_copies: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_rows = []
    for sample, sub in hap_copies.groupby("Sample", sort=False):
        ordered_sub = sub.sort_values("Haplotype_Copy", kind="stable")
        labels = ordered_sub["Mini_Haplotype_Label"].astype(str).tolist()
        dose_110 = int(ordered_sub["rs11078928_ALT"].sum())
        dose_869 = int(ordered_sub["rs869402_ALT"].sum())
        gt_map = {0: "WT", 1: "Het", 2: "Hom"}
        sample_rows.append(
            {
                "Sample": sample,
                "Cohort": ordered_sub["Cohort"].iloc[0],
                "Tissue": ordered_sub["Tissue"].iloc[0],
                "Group": ordered_sub["Group"].iloc[0],
                "Sample_Haplotype_Pattern": " + ".join(labels),
                "Carries_rs11078928_ALT": int((ordered_sub["rs11078928_ALT"] == 1).any()),
                "Carries_rs869402_ALT": int((ordered_sub["rs869402_ALT"] == 1).any()),
                "Carries_Both_ALT_Haplotype": int((ordered_sub["Mini_Haplotype_Code"] == "11").any()),
                "Carries_rs11078928_only_Haplotype": int((ordered_sub["Mini_Haplotype_Code"] == "10").any()),
                "Carries_rs869402_only_Haplotype": int((ordered_sub["Mini_Haplotype_Code"] == "01").any()),
                "rs11078928_Dose": dose_110,
                "rs869402_Dose": dose_869,
                "rs11078928_GT": gt_map[dose_110],
                "rs869402_GT": gt_map[dose_869],
            }
        )
    sample_level = pd.DataFrame(sample_rows)
    sample_level["Carries_Both_Variants_Across_Any_Haplotype"] = (
        (sample_level["Carries_rs11078928_ALT"] == 1) & (sample_level["Carries_rs869402_ALT"] == 1)
    ).astype(int)

    expanded = pd.concat(
        [
            sample_level,
            sample_level[sample_level["Tissue"].astype(str).eq("Healthy")].assign(Group="Healthy_All"),
            sample_level[sample_level["Tissue"].astype(str).eq("Tumour")].assign(Group="Tumour_All"),
            sample_level.assign(Group="All_Samples"),
        ],
        ignore_index=True,
    )

    carrier_cols = [
        "Carries_rs11078928_ALT",
        "Carries_rs869402_ALT",
        "Carries_Both_ALT_Haplotype",
        "Carries_rs11078928_only_Haplotype",
        "Carries_rs869402_only_Haplotype",
        "Carries_Both_Variants_Across_Any_Haplotype",
    ]
    summary = expanded.groupby("Group", as_index=False)[carrier_cols].sum()
    summary["Sample_Count"] = expanded.groupby("Group").size().values
    for col in carrier_cols:
        summary[f"{col}_Proportion"] = summary[col] / summary["Sample_Count"]
    summary = ordered_group(summary, "Group")

    patterns = (
        expanded.groupby(["Group", "Sample_Haplotype_Pattern"], as_index=False)
        .size()
        .rename(columns={"size": "Sample_Count"})
    )
    patterns["Sample_Proportion"] = patterns.groupby("Group")["Sample_Count"].transform(lambda s: s / s.sum())
    patterns = ordered_group(patterns, "Group")
    return summary, patterns, expanded


def make_genotype_concordance_tables(expanded_samples: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    count_rows = []
    summary_rows = []
    for group in GROUP_ORDER:
        sub = expanded_samples[expanded_samples["Group"].astype(str).eq(group)].copy()
        if sub.empty:
            continue
        counts = (
            sub.groupby(["rs11078928_GT", "rs869402_GT"], as_index=False)
            .size()
            .rename(columns={"size": "Sample_Count"})
        )
        counts["Group"] = group
        counts["Row_Proportion_within_rs11078928_GT"] = counts.groupby("rs11078928_GT")["Sample_Count"].transform(lambda s: s / s.sum())
        count_rows.append(counts[["Group", "rs11078928_GT", "rs869402_GT", "Sample_Count", "Row_Proportion_within_rs11078928_GT"]])

        cond_110 = sub[sub["rs11078928_Dose"] > 0]
        cond_869 = sub[sub["rs869402_Dose"] > 0]

        def summarise_partner(df: pd.DataFrame, focal_gt: str) -> str:
            vc = df.loc[df["rs11078928_GT"] == focal_gt, "rs869402_GT"].value_counts()
            if vc.empty:
                return "n/a"
            order = [g for g in ["WT", "Het", "Hom"] if g in vc.index]
            return ", ".join(f"{g}:{int(vc[g])}" for g in order)

        summary_rows.append(
            {
                "Group": group,
                "Sample_Count": int(len(sub)),
                "rs11078928_ALT_Carriers": int(len(cond_110)),
                "rs869402_ALT_Carriers": int(len(cond_869)),
                "P_rs869402_ALT_given_rs11078928_ALT": float((cond_110["rs869402_Dose"] > 0).mean()) if len(cond_110) else pd.NA,
                "P_rs11078928_ALT_given_rs869402_ALT": float((cond_869["rs11078928_Dose"] > 0).mean()) if len(cond_869) else pd.NA,
                "rs869402_GT_if_rs11078928_Het": summarise_partner(sub, "Het"),
                "rs869402_GT_if_rs11078928_Hom": summarise_partner(sub, "Hom"),
            }
        )

    counts_out = ordered_group(pd.concat(count_rows, ignore_index=True), "Group")
    summary_out = ordered_group(pd.DataFrame(summary_rows), "Group")
    return counts_out, summary_out


def make_group_genotype_heatmap(genotype_counts: pd.DataFrame, group: str, filename_stem: str, title_label: str) -> Path:
    out_path = OUT_DIR / f"25_{filename_stem}_Genotype_Heatmap.png"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    order = ["WT", "Het", "Hom"]
    sub = genotype_counts[genotype_counts["Group"].astype(str).eq(group)].copy()
    if sub.empty:
        raise ValueError(f"{group} genotype counts were not available for heatmap generation.")

    count_matrix = (
        sub.pivot(index="rs11078928_GT", columns="rs869402_GT", values="Sample_Count")
        .reindex(index=order, columns=order)
        .fillna(0)
    )
    prop_matrix = (
        sub.pivot(index="rs11078928_GT", columns="rs869402_GT", values="Row_Proportion_within_rs11078928_GT")
        .reindex(index=order, columns=order)
        .fillna(0.0)
    )

    fig, ax = plt.subplots(figsize=(5.4, 4.2), dpi=220)
    im = ax.imshow(count_matrix.values, cmap="Blues", vmin=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(f"{title_label} sample count", fontsize=9)

    ax.set_xticks(range(len(order)), order)
    ax.set_yticks(range(len(order)), order)
    ax.set_xlabel("rs869402 genotype", fontsize=10)
    ax.set_ylabel("rs11078928 genotype", fontsize=10)
    ax.set_title(
        f"{title_label}: genotype pairing between rs11078928 and rs869402\nCounts with row percentages within rs11078928 genotype",
        fontsize=11,
        pad=10,
    )

    max_count = max(float(count_matrix.to_numpy().max()), 1.0)
    for i, row_label in enumerate(order):
        for j, col_label in enumerate(order):
            count = int(count_matrix.loc[row_label, col_label])
            prop = float(prop_matrix.loc[row_label, col_label])
            label = f"{count}\n({prop * 100:.0f}%)" if count > 0 else "0"
            color = "white" if count >= 0.45 * max_count else "#1F2A37"
            ax.text(j, i, label, ha="center", va="center", fontsize=10, fontweight="bold", color=color)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([x - 0.5 for x in range(1, len(order))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(order))], minor=True)
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


def build_supervisor_summary_table(genotype_summary: pd.DataFrame, direct: pd.DataFrame) -> pd.DataFrame:
    rows = []
    direct_lookup = direct.set_index("Cohort")
    interpretation_map = {
        "Breast_Tumour": "Very strong co-inheritance; once rs11078928 is homozygous, rs869402 is usually homozygous too.",
        "Endometrium_Tumour": "Co-inheritance remains strong, but dosage pairing is looser and includes more WT/Het partner states.",
        "Tumour_All": "Across pooled tumours, concordance remains high and the dominant paired state is rs11078928 Hom with rs869402 Hom.",
    }
    for group, _, title_label in HEATMAP_SPECS:
        row = genotype_summary.loc[genotype_summary["Group"].eq(group)].iloc[0]
        direct_row = direct_lookup.loc[DIRECT_COHORT_MAP[group]]
        rows.append(
            {
                "Focus": title_label,
                "Samples": int(row["Sample_Count"]),
                "rs11078928 ALT carriers": int(row["rs11078928_ALT_Carriers"]),
                "rs869402 ALT carriers": int(row["rs869402_ALT_Carriers"]),
                "P(rs869402 ALT | rs11078928 ALT)": f"{float(row['P_rs869402_ALT_given_rs11078928_ALT']):.1%}",
                "P(rs11078928 ALT | rs869402 ALT)": f"{float(row['P_rs11078928_ALT_given_rs869402_ALT']):.1%}",
                "rs869402 if rs11078928 Het": str(row["rs869402_GT_if_rs11078928_Het"]),
                "rs869402 if rs11078928 Hom": str(row["rs869402_GT_if_rs11078928_Hom"]),
                "Mini-haplotype frequency readout": (
                    f"Both_ALT {float(direct_row['Healthy_Both_ALT_Freq']):.1%}->{float(direct_row['Tumour_Both_ALT_Freq']):.1%}; "
                    f"rs11078928-only delta {float(direct_row['Delta_rs11078928_only']):+.1%}; "
                    f"rs869402-only delta {float(direct_row['Delta_rs869402_only']):+.1%}"
                ),
                "Supervisor interpretation": interpretation_map[group],
            }
        )
    return pd.DataFrame(rows)


def build_genotype_matrix_exports(genotype_counts: pd.DataFrame) -> Dict[str, Dict[str, pd.DataFrame]]:
    order = ["WT", "Het", "Hom"]
    exports: Dict[str, Dict[str, pd.DataFrame]] = {}
    for group, _, title_label in HEATMAP_SPECS:
        sub = genotype_counts[genotype_counts["Group"].astype(str).eq(group)].copy()
        count_matrix = (
            sub.pivot(index="rs11078928_GT", columns="rs869402_GT", values="Sample_Count")
            .reindex(index=order, columns=order)
            .fillna(0)
            .astype(int)
        )
        count_matrix.index.name = "rs11078928 genotype"
        count_matrix.columns = [f"rs869402 {col}" for col in count_matrix.columns]

        row_pct_matrix = (
            sub.pivot(index="rs11078928_GT", columns="rs869402_GT", values="Row_Proportion_within_rs11078928_GT")
            .reindex(index=order, columns=order)
            .fillna(0.0)
        )
        row_pct_matrix.index.name = "rs11078928 genotype"
        row_pct_matrix.columns = [f"rs869402 {col}" for col in row_pct_matrix.columns]
        row_pct_display = row_pct_matrix.map(lambda x: f"{float(x):.1%}")

        exports[group] = {
            "title": pd.DataFrame({"Item": ["Focus", "Read this as"], "Value": [title_label, "Rows are rs11078928 genotype; columns are rs869402 genotype. Row percentages stay within the rs11078928 genotype row."]}),
            "counts": count_matrix.reset_index(),
            "row_pct": row_pct_display.reset_index(),
        }
    return exports


def polish_workbook(workbook_path: Path, heatmap_paths: Dict[str, Path]) -> None:
    wb = load_workbook(workbook_path)
    header_fill = PatternFill(fill_type="solid", fgColor="1F7C98")
    header_font = Font(color="FFFFFF", bold=True)
    title_font = Font(bold=True, size=12)
    wrap_alignment = Alignment(wrap_text=True, vertical="top")

    heatmap_sheet = wb.create_sheet("Heatmaps", 2)
    heatmap_sheet["A1"] = "Supervisor-facing genotype heatmaps"
    heatmap_sheet["A1"].font = Font(bold=True, size=14)
    heatmap_sheet["A2"] = "Rows = rs11078928 genotype; columns = rs869402 genotype; labels show sample counts and row percentages within the rs11078928 genotype."
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

    for ws in wb.worksheets:
        max_col = ws.max_column
        max_row = ws.max_row
        if max_row > 1 and ws.title != "Heatmaps":
            for cell in ws[1]:
                cell.fill = header_fill
                cell.font = header_font
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
        for col_idx in range(1, max_col + 1):
            width = 0
            col_letter = get_column_letter(col_idx)
            for cell in ws[col_letter]:
                value = "" if cell.value is None else str(cell.value)
                width = max(width, min(len(value) + 2, 45))
                if cell.row > 1:
                    cell.alignment = wrap_alignment
            ws.column_dimensions[col_letter].width = max(12, width)

    wb.save(workbook_path)


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
    subset = mini_freqs[mini_freqs["Group"].isin(["Breast_Healthy", "Breast_Tumour", "Endometrium_Healthy", "Endometrium_Tumour", "Healthy_All", "Tumour_All"])].copy()
    pivot = subset.pivot_table(index="Group", columns="Mini_Haplotype_Label", values="Haplotype_Copy_Frequency", aggfunc="first", observed=False).fillna(0.0)
    rows = []
    for cohort, healthy_group, tumour_group in [
        ("Breast", "Breast_Healthy", "Breast_Tumour"),
        ("Endometrium", "Endometrium_Healthy", "Endometrium_Tumour"),
        ("All", "Healthy_All", "Tumour_All"),
    ]:
        h = pivot.loc[healthy_group]
        t = pivot.loc[tumour_group]
        delta_110 = float(t.get("rs11078928_only", 0.0) - h.get("rs11078928_only", 0.0))
        delta_869 = float(t.get("rs869402_only", 0.0) - h.get("rs869402_only", 0.0))
        delta_both = float(t.get("Both_ALT", 0.0) - h.get("Both_ALT", 0.0))
        inverse_supported = (delta_110 * delta_869) < 0 and abs(delta_both) < max(abs(delta_110), abs(delta_869))
        rows.append(
            {
                "Cohort": cohort,
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
            }
        )
    return pd.DataFrame(rows)


def make_mini_haplotype_significance_table(mini_freqs: pd.DataFrame) -> pd.DataFrame:
    comparisons = [
        ("Breast", "Breast_Healthy", "Breast_Tumour"),
        ("Endometrium", "Endometrium_Healthy", "Endometrium_Tumour"),
        ("All", "Healthy_All", "Tumour_All"),
    ]
    rows = []
    for cohort, healthy_group, tumour_group in comparisons:
        sub = mini_freqs[mini_freqs["Group"].isin([healthy_group, tumour_group])].copy()
        healthy_total = int(sub.loc[sub["Group"].eq(healthy_group), "Haplotype_Copy_Count"].sum())
        tumour_total = int(sub.loc[sub["Group"].eq(tumour_group), "Haplotype_Copy_Count"].sum())
        for label in ["Both_ALT", "rs869402_only", "rs11078928_only", "Neither_ALT"]:
            healthy_count = int(sub.loc[(sub["Group"].eq(healthy_group)) & (sub["Mini_Haplotype_Label"].eq(label)), "Haplotype_Copy_Count"].iloc[0])
            tumour_count = int(sub.loc[(sub["Group"].eq(tumour_group)) & (sub["Mini_Haplotype_Label"].eq(label)), "Haplotype_Copy_Count"].iloc[0])
            odds_ratio, p_value = fisher_exact([[tumour_count, tumour_total - tumour_count], [healthy_count, healthy_total - healthy_count]])
            rows.append(
                {
                    "Cohort": cohort,
                    "Mini_Haplotype_Label": label,
                    "Tumour_Count": tumour_count,
                    "Tumour_Total": tumour_total,
                    "Healthy_Count": healthy_count,
                    "Healthy_Total": healthy_total,
                    "Tumour_Freq": tumour_count / tumour_total,
                    "Healthy_Freq": healthy_count / healthy_total,
                    "Odds_Ratio": odds_ratio,
                    "P_Fisher": p_value,
                }
            )
    out = pd.DataFrame(rows)
    adjusted = []
    for _, sub in out.groupby("Cohort", sort=False):
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

def build_interpretation(direct: pd.DataFrame, ld_context: pd.DataFrame, full_relevant: pd.DataFrame, genotype_summary: pd.DataFrame) -> str:
    global_both = full_relevant.loc[full_relevant["Pattern_Label"].eq("Both_ALT"), "Frequency"].sum()
    global_110_only = full_relevant.loc[full_relevant["Pattern_Label"].eq("rs11078928_only"), "Frequency"].sum()
    global_869_only = full_relevant.loc[full_relevant["Pattern_Label"].eq("rs869402_only"), "Frequency"].sum()
    breast = direct.loc[direct["Cohort"].eq("Breast")].iloc[0]
    endo = direct.loc[direct["Cohort"].eq("Endometrium")].iloc[0]
    pooled = genotype_summary.loc[genotype_summary["Group"].eq("Tumour_All")].iloc[0]
    breast_gt = genotype_summary.loc[genotype_summary["Group"].eq("Breast_Tumour")].iloc[0]
    endo_gt = genotype_summary.loc[genotype_summary["Group"].eq("Endometrium_Tumour")].iloc[0]
    rs869402_ld = ld_context.loc[ld_context["rsID"].eq("rs869402"), "LD_Block_Status"].iloc[0]

    lines = [
        "Focused phased-haplotype and genotype-concordance interpretation for rs11078928 and rs869402",
        "",
        (
            f"Across the full-region phased backbone, the dominant relevant background carried both alt alleles "
            f"(global frequency {global_both:.3f}), whereas rs869402-only ({global_869_only:.3f}) and "
            f"rs11078928-only ({global_110_only:.3f}) backgrounds were less frequent."
        ),
        (
            f"In breast tumours, co-inheritance is strongest: {breast_gt['P_rs869402_ALT_given_rs11078928_ALT']:.3f} of rs11078928-alt carriers also carried rs869402 alt, and {breast_gt['P_rs11078928_ALT_given_rs869402_ALT']:.3f} of rs869402-alt carriers also carried rs11078928 alt."
        ),
        (
            f"Breast dosage pairing is the clearest supervisor-facing result: rs11078928 Het maps to {breast_gt['rs869402_GT_if_rs11078928_Het']}, whereas rs11078928 Hom maps mostly to {breast_gt['rs869402_GT_if_rs11078928_Hom']}."
        ),
        (
            f"Endometrial tumours still show strong concordance but a looser dosage pattern: {endo_gt['P_rs869402_ALT_given_rs11078928_ALT']:.3f} of rs11078928-alt carriers also carried rs869402 alt, with heterozygous rs11078928 mapping to {endo_gt['rs869402_GT_if_rs11078928_Het']} and homozygous rs11078928 mapping to {endo_gt['rs869402_GT_if_rs11078928_Hom']}."
        ),
        (
            f"Across pooled tumours, concordance remains high overall: {pooled['P_rs869402_ALT_given_rs11078928_ALT']:.3f} of rs11078928-alt carriers also carried rs869402 alt, and {pooled['P_rs11078928_ALT_given_rs869402_ALT']:.3f} of rs869402-alt carriers also carried rs11078928 alt."
        ),
        (
            f"In parallel, the haplotype-frequency view still shows the strongest breast shift as an increase in the Both_ALT background ({breast['Healthy_Both_ALT_Freq']:.3f} in healthy vs {breast['Tumour_Both_ALT_Freq']:.3f} in tumour), while both rs869402-only and rs11078928-only backgrounds decreased."
        ),
        (
            f"In endometrial tumours, the same qualitative pattern was weaker: Both_ALT increased from {endo['Healthy_Both_ALT_Freq']:.3f} to {endo['Tumour_Both_ALT_Freq']:.3f}, while both single-variant backgrounds were only slightly lower in tumours than in healthy samples."
        ),
        (
            "This means the key interpretation is concordance and dosage pairing, not a clean reciprocal model in which one variant background simply replaces the other."
        ),
        (
            f"The relationship is clearer in full-region haplotypes than in LD-block haplotypes, because rs11078928 is captured within LD_Block_r2_080_Block7 whereas rs869402 is {rs869402_ld.lower()} and therefore cannot be jointly assessed with rs11078928 inside the current explicit LD-block framework."
        ),
        (
            "The workbook now includes breast, endometrial, and pooled-tumour genotype heatmaps plus cohort-specific count and row-percentage matrices so the supervisor can inspect the pairing directly."
        ),
    ]
    return "\n".join(lines) + "\n"


def write_outputs(
    target_positions: pd.DataFrame,
    ld_context: pd.DataFrame,
    mini_freqs: pd.DataFrame,
    sample_carriers: pd.DataFrame,
    sample_patterns: pd.DataFrame,
    genotype_counts: pd.DataFrame,
    genotype_summary: pd.DataFrame,
    full_relevant: pd.DataFrame,
    full_summary: pd.DataFrame,
    block7_freqs: pd.DataFrame,
    block7_assoc: pd.DataFrame,
    direct: pd.DataFrame,
    significance: pd.DataFrame,
    interpretation: str,
) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    heatmap_paths = make_requested_heatmaps(genotype_counts)
    supervisor_summary = build_supervisor_summary_table(genotype_summary, direct)
    matrix_exports = build_genotype_matrix_exports(genotype_counts)

    tsvs = {
        "25_Target_SNP_Metadata.tsv": target_positions,
        "25_LD_Context.tsv": ld_context,
        "25_Mini_Haplotype_Frequencies.tsv": mini_freqs,
        "25_Sample_Carrier_Summary.tsv": sample_carriers,
        "25_Sample_Haplotype_Patterns.tsv": sample_patterns,
        "25_Genotype_Concordance_Counts.tsv": genotype_counts,
        "25_Genotype_Conditional_Summary.tsv": genotype_summary,
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
        df.to_csv(OUT_DIR / name, sep="	", index=False)

    (OUT_DIR / "25_Interpretation.txt").write_text(interpretation, encoding="utf-8")

    workbook_path = OUT_DIR / "25_rs11078928_rs869402_haplotype_focus.xlsx"
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        pd.DataFrame(
            {
                "Item": [
                    "Question",
                    "Best sheets for supervisor review",
                    "Primary input workbook",
                    "Primary phased matrix",
                    "Interpretation rule",
                    "What global means here",
                    "Key limitation",
                    "Statistical interpretation",
                    "Where to find the formal tests",
                ],
                "Value": [
                    "Focused phased-haplotype relationship between rs11078928 and rs869402",
                    "Supervisor_Summary, MiniHap_Signif, Heatmaps, Breast_GT_Matrix, Endometrium_GT_Matrix, Global_GT_Matrix",
                    str(HAPLO_STATS_XLSX),
                    str(PHASED_TSV),
                    "Treat as a haplotype-focused question, not isolated SNP correlation",
                    "Global refers to pooled tumour samples across breast and endometrium.",
                    "rs869402 is not part of any explicit r2=0.80 LD block in the current stage-15 workbook",
                    "Exact tests in this workbook are based on phased haplotype-copy frequencies, not independent per-sample carrier counts.",
                    "MiniHap_Signif reports Fisher exact p-values plus BH FDR for Both_ALT, rs869402_only, rs11078928_only, and Neither_ALT within each cohort.",
                ],
            }
        ).to_excel(writer, sheet_name="README", index=False)
        supervisor_summary.to_excel(writer, sheet_name="Supervisor_Summary", index=False)
        for group, exports in matrix_exports.items():
            sheet_name = MATRIX_SHEET_MAP[group]
            exports["counts"].to_excel(writer, sheet_name=sheet_name, index=False, startrow=2)
            exports["row_pct"].to_excel(writer, sheet_name=sheet_name, index=False, startrow=10)
        target_positions.to_excel(writer, sheet_name="Target_SNPs", index=False)
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
        significance.to_excel(writer, sheet_name="MiniHap_Signif", index=False)
        pd.DataFrame({"Interpretation": interpretation.splitlines()}).to_excel(writer, sheet_name="Interpretation", index=False)

    wb = load_workbook(workbook_path)
    for group, exports in matrix_exports.items():
        ws = wb[MATRIX_SHEET_MAP[group]]
        title_label = exports["title"].iloc[0, 1]
        ws["A1"] = f"{title_label}: sample-count matrix"
        ws["A1"].font = Font(bold=True, size=12)
        ws["A2"] = "Rows = rs11078928 genotype; columns = rs869402 genotype."
        ws["A9"] = f"{title_label}: row percentages within rs11078928 genotype"
        ws["A9"].font = Font(bold=True, size=12)
    wb.save(workbook_path)
    polish_workbook(workbook_path, heatmap_paths)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Focused phased-haplotype analysis for rs11078928 and rs869402")
    return parser.parse_args()


def main() -> None:
    parse_args()
    workbook, phased, sample_metadata = load_required_inputs()
    target_positions = get_target_positions(workbook)
    ld_context, _ = get_ld_context(workbook)
    hap_copies = extract_haplotype_copies(phased, sample_metadata, target_positions)
    mini_freqs = make_mini_haplotype_frequency_table(hap_copies)
    sample_carriers, sample_patterns, expanded_samples = make_sample_level_tables(hap_copies)
    genotype_counts, genotype_summary = make_genotype_concordance_tables(expanded_samples)
    full_relevant, full_summary = make_full_region_tables(workbook)
    block7_freqs, block7_assoc = load_ld_block7_tables(workbook)
    direct = make_direct_comparison_table(mini_freqs)
    significance = make_mini_haplotype_significance_table(mini_freqs)
    interpretation = build_interpretation(direct, ld_context, full_relevant, genotype_summary)
    write_outputs(
        target_positions=target_positions,
        ld_context=ld_context,
        mini_freqs=mini_freqs,
        sample_carriers=sample_carriers,
        sample_patterns=sample_patterns,
        genotype_counts=genotype_counts,
        genotype_summary=genotype_summary,
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





