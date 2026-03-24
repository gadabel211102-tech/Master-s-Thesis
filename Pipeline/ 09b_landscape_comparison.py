"""Script 09b: Normal-versus-Tumour Landscape Comparison
=========================================================

Purpose
-------
Formalise the follow-up review of the landscape plots from scripts 08 and 09.
The script compares normal and tumour carrier frequencies, identifies
normal-only and tumour-only variants, and checks whether the visual patterns
remain after restricting the comparison to higher-quality samples.

Inputs
------
- ``Biological_Annotations`` sheet from the consolidated annotated report
- coverage-audit CSV files produced by ``03_qc_visualisation.py``

Outputs
-------
- ``Landscape_Normal_vs_Tumour_with_QC_and_Significance.xlsx``
- companion CSV summaries for downstream review and manuscript support
"""

from __future__ import annotations

from math import comb
import re
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from pipeline_utils import ensure_directory, find_col, get_paths
from pipeline_validation import (
    print_validation_summary,
    validate_file_exists,
    validate_required_columns,
)

PATHS = get_paths()
INPUT_FILE = PATHS["annotated_report"]
QC_DIR = PATHS["qc_visualisation_dir"]
OUTPUT_DIR = ensure_directory(
    PATHS.get("landscape_comparison_dir", PATHS["results_dir"] / "09b_landscape_comparison")
)
WORKBOOK_PATH = OUTPUT_DIR / "Landscape_Normal_vs_Tumour_with_QC_Significance_Genes_UPDATED.xlsx"
SUMMARY_CSV = OUTPUT_DIR / "Landscape_Normal_vs_Tumour_Summary.csv"
GSDMB_TESTS_CSV = OUTPUT_DIR / "GSDMB_Frequency_Tests.csv"
GLOBAL_TESTS_CSV = OUTPUT_DIR / "Global_Frequency_Tests.csv"
SIGNIFICANCE_SUMMARY_CSV = OUTPUT_DIR / "Landscape_Frequency_Significance_Summary.csv"

HEADER_FILL = PatternFill("solid", fgColor="1F4E78")
HEADER_FONT = Font(color="FFFFFF", bold=True)
SECTION_FILL = PatternFill("solid", fgColor="D9EAF7")


def normalize_sample_name(value: object) -> str:
    """Reduce sample names to a stable comparison key across QC and annotations."""
    text = str(value).strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")


def extract_rsid(value: object) -> str | None:
    """Return the first rsID embedded in an annotation field, if present."""
    text = str(value).strip()
    match = re.search(r"(rs\d+)", text, flags=re.IGNORECASE)
    return match.group(1).lower() if match else None


def classify_location(consequence: object) -> str:
    """Collapse VEP consequence strings into an easier-to-read location class."""
    text = str(consequence).strip().lower()
    if not text or text == 'nan':
        return 'Unknown'
    if 'intergenic_variant' in text:
        return 'Intergenic'
    if 'upstream_gene_variant' in text:
        return 'Upstream'
    if 'downstream_gene_variant' in text:
        return 'Downstream'
    if 'intron_variant' in text:
        return 'Intronic'
    if 'splice' in text:
        return 'Splicing'
    if '5_prime_utr_variant' in text or '3_prime_utr_variant' in text:
        return 'UTR'
    if any(term in text for term in [
        'missense_variant', 'synonymous_variant', 'stop_gained', 'stop_lost',
        'frameshift_variant', 'start_lost', 'protein_altering_variant',
        'coding_sequence_variant', 'inframe_insertion', 'inframe_deletion'
    ]):
        return 'Exonic/Coding'
    if 'non_coding_transcript_exon_variant' in text:
        return 'Exonic/Non-coding'
    return 'Other'


def fisher_two_sided(normal_carriers: int, normal_noncarriers: int, tumour_carriers: int, tumour_noncarriers: int) -> float:
    """Compute a two-sided Fisher exact p-value for a 2x2 carrier table."""
    row1 = normal_carriers + normal_noncarriers
    row2 = tumour_carriers + tumour_noncarriers
    col1 = normal_carriers + tumour_carriers
    total = row1 + row2

    def hypergeom_prob(x_val: int) -> float:
        return comb(row1, x_val) * comb(row2, col1 - x_val) / comb(total, col1)

    min_x = max(0, col1 - row2)
    max_x = min(row1, col1)
    observed = hypergeom_prob(normal_carriers)
    p_value = 0.0
    for x_val in range(min_x, max_x + 1):
        prob = hypergeom_prob(x_val)
        if prob <= observed + 1e-12:
            p_value += prob
    return min(p_value, 1.0)


def benjamini_hochberg(p_values: list[float | None]) -> list[float | None]:
    """Apply Benjamini-Hochberg correction to a list of p-values."""
    n_tests = len(p_values)
    if n_tests == 0:
        return []

    order = sorted(
        range(n_tests),
        key=lambda idx: float("inf") if pd.isna(p_values[idx]) else p_values[idx],
    )
    ranked = [p_values[idx] for idx in order]
    adjusted: list[float | None] = [None] * n_tests
    previous = 1.0

    for rank in range(n_tests, 0, -1):
        p_value = ranked[rank - 1]
        current = 1.0 if pd.isna(p_value) else min(previous, p_value * n_tests / rank)
        previous = current
        adjusted[order[rank - 1]] = current
    return adjusted


def load_annotations() -> pd.DataFrame:
    """Load and standardise the annotated variant workbook."""
    validate_file_exists(INPUT_FILE, "Script 09b input")
    df = pd.read_excel(INPUT_FILE, sheet_name="Biological_Annotations")

    sample_col = find_col(df, "Sample")
    cohort_col = find_col(df, "Cohort")
    tissue_col = find_col(df, "Tissue")
    symbol_col = find_col(df, "Symbol")
    impact_col = find_col(df, "Impact")
    pos_col = find_col(df, "Pos")
    existing_col = find_col(df, "Existing_variation")
    gene_col = find_col(df, "Gene")
    consequence_col = find_col(df, "Consequence")
    feature_type_col = find_col(df, "Feature_type")
    exon_col = find_col(df, "EXON")
    intron_col = find_col(df, "INTRON")

    required = [
        sample_col, cohort_col, tissue_col, symbol_col, impact_col, pos_col,
        existing_col, gene_col, consequence_col, feature_type_col, exon_col, intron_col,
    ]
    validate_required_columns(df, required, "Script 09b Biological_Annotations")
    print_validation_summary(df, sample_col, "Script 09b raw annotations", [cohort_col, tissue_col])

    out = df[[
        sample_col, cohort_col, tissue_col, symbol_col, impact_col, pos_col,
        existing_col, gene_col, consequence_col, feature_type_col, exon_col, intron_col,
    ]].copy()
    out.columns = [
        "Sample", "Cohort", "Tissue", "Symbol", "Impact", "Pos",
        "Existing_variation", "Gene_ID", "Consequence", "Feature_type", "EXON", "INTRON",
    ]
    out["Cohort"] = out["Cohort"].astype(str).str.strip().str.title()
    out["Tissue"] = out["Tissue"].astype(str).str.strip().str.title()
    out["Impact"] = out["Impact"].astype(str).str.strip().str.upper()
    out["sample_key"] = out["Sample"].map(normalize_sample_name)
    out["rsID"] = out["Existing_variation"].map(extract_rsid)
    out["Location_Class"] = out["Consequence"].map(classify_location)
    return out[out["Tissue"].isin(["Normal", "Tumour"])].copy()


def load_qc_audits() -> pd.DataFrame:
    """Load the coverage-audit summaries produced by script 03."""
    audit_map = [
        ("Breast", "Normal", QC_DIR / "CoverageAudit_Breast_Control.csv"),
        ("Breast", "Tumour", QC_DIR / "CoverageAudit_Breast_Tumour.csv"),
        ("Endometrium", "Normal", QC_DIR / "CoverageAudit_Endometrium_Control.csv"),
        ("Endometrium", "Tumour", QC_DIR / "CoverageAudit_Endometrium_Tumour.csv"),
    ]

    frames: list[pd.DataFrame] = []
    for cohort, tissue, path in audit_map:
        validate_file_exists(path, f"Script 09b QC audit {cohort}/{tissue}")
        frame = pd.read_csv(path)
        sample_col = find_col(frame, "sample")
        zero_col = find_col(frame, "zeros")
        validate_required_columns(frame, [sample_col, zero_col], f"Script 09b QC audit {cohort}/{tissue}")

        sub = frame[[sample_col, zero_col]].copy()
        sub.columns = ["sample_raw", "zero_amplicons"]
        sub["sample_key"] = sub["sample_raw"].map(normalize_sample_name)
        sub["Cohort"] = cohort
        sub["Tissue"] = tissue
        sub["high_qc"] = sub["zero_amplicons"].fillna(999).eq(0)
        frames.append(sub)
    return pd.concat(frames, ignore_index=True)


def prepare_scope_table(annotations: pd.DataFrame, scope: str) -> pd.DataFrame:
    """Return the per-sample carrier table used for one scope."""
    if scope == "GSDMB":
        subset = annotations[annotations["Symbol"].astype(str).str.contains("GSDMB", case=False, na=False)].copy()
        subset = subset.drop_duplicates(["Sample", "sample_key", "Cohort", "Tissue", "Pos", "Impact"])
    else:
        subset = annotations.drop_duplicates(["Sample", "sample_key", "Cohort", "Tissue", "Symbol", "Pos", "Impact"])
    return subset


def build_frequency_tests(scope_df: pd.DataFrame, key_cols: list[str]) -> pd.DataFrame:
    """Summarise per-variant carrier frequencies and significance tests."""
    denominators = scope_df.groupby(["Cohort", "Tissue"])["Sample"].nunique().to_dict()
    carriers = (
        scope_df.groupby(key_cols + ["Tissue"])
        .agg(
            Carrier_Count=("Sample", "nunique"),
            rsID=("rsID", lambda values: next((x for x in values if pd.notna(x)), None)),
            Gene_Symbol=("Symbol", lambda values: next((x for x in values if pd.notna(x)), None)),
            Gene_ID=("Gene_ID", lambda values: next((x for x in values if pd.notna(x)), None)),
            Consequence=("Consequence", lambda values: next((x for x in values if pd.notna(x)), None)),
            Feature_type=("Feature_type", lambda values: next((x for x in values if pd.notna(x)), None)),
            EXON=("EXON", lambda values: next((x for x in values if pd.notna(x)), None)),
            INTRON=("INTRON", lambda values: next((x for x in values if pd.notna(x)), None)),
            Location_Class=("Location_Class", lambda values: next((x for x in values if pd.notna(x)), None)),
        )
        .reset_index()
    )

    rows: list[dict[str, object]] = []
    for keys, group in carriers.groupby(key_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(key_cols, keys))
        row["Normal_Carriers"] = int(group.loc[group["Tissue"] == "Normal", "Carrier_Count"].sum())
        row["Tumour_Carriers"] = int(group.loc[group["Tissue"] == "Tumour", "Carrier_Count"].sum())
        row["Normal_N"] = int(denominators.get((row["Cohort"], "Normal"), 0))
        row["Tumour_N"] = int(denominators.get((row["Cohort"], "Tumour"), 0))
        row["rsID"] = next((x for x in group["rsID"] if pd.notna(x)), None)
        row["Gene_Symbol"] = row.get("Symbol") or next((x for x in group["Gene_Symbol"] if pd.notna(x)), None)
        row["Gene_ID"] = next((x for x in group["Gene_ID"] if pd.notna(x)), None)
        row["Consequence"] = next((x for x in group["Consequence"] if pd.notna(x)), None)
        row["Feature_type"] = next((x for x in group["Feature_type"] if pd.notna(x)), None)
        row["EXON"] = next((x for x in group["EXON"] if pd.notna(x)), None)
        row["INTRON"] = next((x for x in group["INTRON"] if pd.notna(x)), None)
        row["Location_Class"] = next((x for x in group["Location_Class"] if pd.notna(x)), None)
        row["Normal_Frequency_pct"] = (row["Normal_Carriers"] / row["Normal_N"] * 100) if row["Normal_N"] else None
        row["Tumour_Frequency_pct"] = (row["Tumour_Carriers"] / row["Tumour_N"] * 100) if row["Tumour_N"] else None
        row["Tumour_minus_Normal_pct"] = row["Tumour_Frequency_pct"] - row["Normal_Frequency_pct"]
        row["Variant_Status"] = (
            "shared"
            if row["Normal_Carriers"] > 0 and row["Tumour_Carriers"] > 0
            else ("tumour_only" if row["Tumour_Carriers"] > 0 else "normal_only")
        )
        row["Direction"] = (
            "Lower in tumour"
            if row["Tumour_minus_Normal_pct"] < -1e-9
            else ("Higher in tumour" if row["Tumour_minus_Normal_pct"] > 1e-9 else "Same frequency")
        )
        row["p_value"] = fisher_two_sided(
            int(row["Normal_Carriers"]),
            int(row["Normal_N"] - row["Normal_Carriers"]),
            int(row["Tumour_Carriers"]),
            int(row["Tumour_N"] - row["Tumour_Carriers"]),
        )
        rows.append(row)

    tests = pd.DataFrame(rows)
    tests["FDR_BH"] = None
    for _, indices in tests.groupby("Cohort").groups.items():
        tests.loc[list(indices), "FDR_BH"] = benjamini_hochberg(tests.loc[list(indices), "p_value"].tolist())
    tests["p_value"] = tests["p_value"].astype(float)
    tests["FDR_BH"] = tests["FDR_BH"].astype(float)
    tests["Raw_p_lt_0.05"] = tests["p_value"] < 0.05
    tests["FDR_lt_0.05"] = tests["FDR_BH"] < 0.05
    preferred_order = [
        "Cohort", "Gene_Symbol", "Gene_ID", "Pos", "Impact", "Consequence",
        "Location_Class", "Feature_type", "EXON", "INTRON", "Normal_Carriers",
        "Tumour_Carriers", "Normal_N", "Tumour_N", "rsID", "Normal_Frequency_pct",
        "Tumour_Frequency_pct", "Tumour_minus_Normal_pct", "Variant_Status",
        "Direction", "p_value", "FDR_BH", "Raw_p_lt_0.05", "FDR_lt_0.05",
    ]
    tail_cols = [col for col in tests.columns if col not in preferred_order]
    return tests[preferred_order + tail_cols]


def build_variant_status_summary(tests_df: pd.DataFrame) -> pd.DataFrame:
    """Count shared, normal-only, and tumour-only variants by impact class."""
    return (
        tests_df.groupby(["Cohort", "Impact", "Variant_Status"])
        .size()
        .reset_index(name="n_variants")
        .sort_values(["Cohort", "Impact", "Variant_Status"])
    )


def build_direction_summary(tests_df: pd.DataFrame, direction_col: str, out_col: str) -> pd.DataFrame:
    """Count how many shared variants move up or down in tumour."""
    shared = tests_df[tests_df["Variant_Status"] == "shared"].copy()
    return (
        shared.groupby(["Cohort", "Impact", direction_col])
        .size()
        .reset_index(name=out_col)
        .sort_values(["Cohort", "Impact", direction_col])
    )


def build_qc_sensitivity(scope_df: pd.DataFrame, qc_df: pd.DataFrame, key_cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Repeat the comparison on aligned high-QC samples only."""
    aligned = scope_df.merge(
        qc_df[["sample_key", "Cohort", "Tissue", "zero_amplicons", "high_qc"]],
        on=["sample_key", "Cohort", "Tissue"],
        how="inner",
    )

    aligned_counts = (
        aligned[aligned["high_qc"]]
        .groupby(key_cols + ["Tissue"])
        .agg(Carrier_Count=("Sample", "nunique"))
        .reset_index()
    )
    high_qc_denoms = aligned[aligned["high_qc"]].groupby(["Cohort", "Tissue"])["Sample"].nunique().to_dict()

    rows: list[dict[str, object]] = []
    for keys, group in aligned_counts.groupby(key_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(key_cols, keys))
        normal_n = int(high_qc_denoms.get((row["Cohort"], "Normal"), 0))
        tumour_n = int(high_qc_denoms.get((row["Cohort"], "Tumour"), 0))
        normal_carriers = int(group.loc[group["Tissue"] == "Normal", "Carrier_Count"].sum())
        tumour_carriers = int(group.loc[group["Tissue"] == "Tumour", "Carrier_Count"].sum())
        row["Normal_Frequency_pct_HighQC"] = (normal_carriers / normal_n * 100) if normal_n else None
        row["Tumour_Frequency_pct_HighQC"] = (tumour_carriers / tumour_n * 100) if tumour_n else None
        if row["Normal_Frequency_pct_HighQC"] is None or row["Tumour_Frequency_pct_HighQC"] is None:
            row["Tumour_minus_Normal_pct_HighQC"] = None
            row["Direction_HighQC"] = "Not estimable"
        else:
            delta = row["Tumour_Frequency_pct_HighQC"] - row["Normal_Frequency_pct_HighQC"]
            row["Tumour_minus_Normal_pct_HighQC"] = delta
            row["Direction_HighQC"] = (
                "Lower in tumour" if delta < -1e-9 else ("Higher in tumour" if delta > 1e-9 else "Same frequency")
            )
        rows.append(row)

    high_qc_table = pd.DataFrame(rows)
    qc_sample_counts = (
        aligned.groupby(["Cohort", "Tissue"])
        .apply(
            lambda group: pd.Series(
                {
                    "aligned_samples": group["sample_key"].nunique(),
                    "high_qc_samples": group.loc[group["high_qc"], "sample_key"].nunique(),
                    "mean_zero_amplicons": group[["sample_key", "zero_amplicons"]].drop_duplicates()["zero_amplicons"].mean(),
                }
            ),
            include_groups=False,
        )
        .reset_index()
    )
    return high_qc_table, qc_sample_counts


def build_summary(scope_label: str, status_df: pd.DataFrame, tests_df: pd.DataFrame, qc_counts_df: pd.DataFrame) -> pd.DataFrame:
    """Create a readable summary table for manuscript or email support."""
    rows: list[dict[str, object]] = []
    for cohort in sorted(tests_df["Cohort"].dropna().unique()):
        cohort_status = status_df[status_df["Cohort"] == cohort].copy()
        cohort_tests = tests_df[tests_df["Cohort"] == cohort].copy()
        pivot = cohort_status.pivot_table(index="Impact", columns="Variant_Status", values="n_variants", fill_value=0)
        shared = cohort_tests[cohort_tests["Variant_Status"] == "shared"]
        lower_count = int((shared["Direction"] == "Lower in tumour").sum())
        higher_count = int((shared["Direction"] == "Higher in tumour").sum())
        raw_hits = cohort_tests[cohort_tests["Raw_p_lt_0.05"]].sort_values("p_value")
        hit_text = ", ".join(
            f"{row.rsID or 'no_rsID'} (p={row.p_value:.4g})"
            for row in raw_hits.head(5).itertuples(index=False)
        ) or "None"
        qc_rows = qc_counts_df[qc_counts_df["Cohort"] == cohort]
        qc_text = "; ".join(
            f"{row.Tissue}: aligned n={int(row.aligned_samples)}, high-QC n={int(row.high_qc_samples)}"
            for row in qc_rows.itertuples(index=False)
        )
        best_row = cohort_tests.sort_values("p_value").iloc[0]
        rows.append(
            {
                "Map": scope_label,
                "Cohort": cohort,
                "Tumour_only_variants": int((cohort_tests["Variant_Status"] == "tumour_only").sum()),
                "Normal_only_variants": int((cohort_tests["Variant_Status"] == "normal_only").sum()),
                "Shared_variants": int((cohort_tests["Variant_Status"] == "shared").sum()),
                "Shared_lower_in_tumour": lower_count,
                "Shared_higher_in_tumour": higher_count,
                "Raw_p_lt_0.05": int(cohort_tests["Raw_p_lt_0.05"].sum()),
                "FDR_lt_0.05": int(cohort_tests["FDR_lt_0.05"].sum()),
                "Best_rsID": best_row["rsID"] if pd.notna(best_row["rsID"]) else "no_rsID",
                "Best_variant_status": best_row["Variant_Status"],
                "Best_direction": best_row["Direction"],
                "Best_p_value": float(best_row["p_value"]),
                "Best_FDR_BH": float(best_row["FDR_BH"]),
                "Example_nominal_hits": hit_text,
                "Tumour_only_breakdown": "; ".join(
                    f"{impact}: {int(pivot.loc[impact, 'tumour_only'])}"
                    for impact in pivot.index
                    if "tumour_only" in pivot.columns and int(pivot.loc[impact, "tumour_only"]) > 0
                ) or "None",
                "Normal_only_breakdown": "; ".join(
                    f"{impact}: {int(pivot.loc[impact, 'normal_only'])}"
                    for impact in pivot.index
                    if "normal_only" in pivot.columns and int(pivot.loc[impact, "normal_only"]) > 0
                ) or "None",
                "QC_alignment": qc_text or "Not available",
            }
        )
    return pd.DataFrame(rows)


def style_workbook(path: Path) -> None:
    """Apply lightweight formatting so the workbook is easy to share."""
    workbook = load_workbook(path)
    for sheet in workbook.worksheets:
        if sheet.title == "Guide":
            sheet.column_dimensions["A"].width = 22
            sheet.column_dimensions["B"].width = 100
            for row in sheet.iter_rows(min_row=1, max_row=sheet.max_row, min_col=1, max_col=2):
                row[0].font = Font(bold=True)
                row[0].fill = SECTION_FILL
                row[0].alignment = Alignment(vertical="top", wrap_text=True)
                row[1].alignment = Alignment(vertical="top", wrap_text=True)
            continue

        sheet.freeze_panes = "A2"
        for cell in sheet[1]:
            cell.fill = HEADER_FILL
            cell.font = HEADER_FONT
            cell.alignment = Alignment(wrap_text=True, horizontal="center", vertical="center")

        for col_cells in sheet.columns:
            width = max(len(str(cell.value)) if cell.value is not None else 0 for cell in col_cells) + 2
            sheet.column_dimensions[get_column_letter(col_cells[0].column)].width = min(max(width, 12), 36)

        for row in sheet.iter_rows(min_row=2):
            for cell in row:
                cell.alignment = Alignment(wrap_text=True, vertical="top")
    workbook.save(path)


def main() -> None:
    """Run the landscape follow-up analysis and write the deliverables."""
    print("======================================================================")
    print("SCRIPT 09B: NORMAL-VERSUS-TUMOUR LANDSCAPE COMPARISON")
    print("======================================================================")
    print(f"Input annotations: {INPUT_FILE}")
    print(f"QC audit directory: {QC_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print("======================================================================")

    annotations = load_annotations()
    qc_audits = load_qc_audits()

    gsdmb_scope = prepare_scope_table(annotations, "GSDMB")
    global_scope = prepare_scope_table(annotations, "GLOBAL")

    gsdmb_tests = build_frequency_tests(gsdmb_scope, ["Cohort", "Pos", "Impact"])
    global_tests = build_frequency_tests(global_scope, ["Cohort", "Symbol", "Pos", "Impact"])

    gsdmb_status = build_variant_status_summary(gsdmb_tests)
    global_status = build_variant_status_summary(global_tests)
    gsdmb_direction = build_direction_summary(gsdmb_tests, "Direction", "n_shared_variants")
    global_direction = build_direction_summary(global_tests, "Direction", "n_shared_variants")

    gsdmb_high_qc, gsdmb_qc_counts = build_qc_sensitivity(gsdmb_scope, qc_audits, ["Cohort", "Pos", "Impact"])
    global_high_qc, global_qc_counts = build_qc_sensitivity(global_scope, qc_audits, ["Cohort", "Symbol", "Pos", "Impact"])
    gsdmb_tests = gsdmb_tests.merge(gsdmb_high_qc, on=["Cohort", "Pos", "Impact"], how="left")
    global_tests = global_tests.merge(global_high_qc, on=["Cohort", "Symbol", "Pos", "Impact"], how="left")

    if "Symbol" in global_tests.columns:
        global_tests = global_tests.drop(columns=["Symbol"])

    gsdmb_qc_source = gsdmb_tests[gsdmb_tests["Direction_HighQC"].notna()].copy()
    gsdmb_qc_source["Direction_QC"] = gsdmb_qc_source["Direction_HighQC"]
    global_qc_source = global_tests[global_tests["Direction_HighQC"].notna()].copy()
    global_qc_source["Direction_QC"] = global_qc_source["Direction_HighQC"]

    gsdmb_qc_direction = build_direction_summary(
        gsdmb_qc_source,
        "Direction_QC",
        "n_shared_variants_high_qc",
    )
    global_qc_direction = build_direction_summary(
        global_qc_source,
        "Direction_QC",
        "n_shared_variants_high_qc",
    )

    summary = pd.concat(
        [
            build_summary("GSDMB (script 9)", gsdmb_status, gsdmb_tests, gsdmb_qc_counts),
            build_summary("Global (script 8)", global_status, global_tests, global_qc_counts),
        ],
        ignore_index=True,
    )

    def label_significant_row(row: pd.Series) -> str:
        gene_part = f"{row['Gene_Symbol']} " if 'Gene_Symbol' in row.index and pd.notna(row.get('Gene_Symbol')) else ''
        rs_part = row["rsID"] if pd.notna(row.get("rsID")) else f"chr17:{int(row['Pos'])}"
        return f"{gene_part}{rs_part} (p={row['p_value']:.4g}, FDR={row['FDR_BH']:.4g})"

    significance_counts_rows: list[dict[str, object]] = []
    significance_detail_rows: list[dict[str, object]] = []
    detail_columns = [
        "Map", "Cohort", "Variants_tested", "Raw_p_lt_0.05", "FDR_lt_0.05",
        "Significance_Level", "Also_FDR_significant", "Gene_Symbol", "Gene_ID", "rsID",
        "Pos", "Impact", "Consequence", "Location_Class", "Feature_type",
        "EXON", "INTRON", "Variant_Status", "Direction", "Normal_Carriers",
        "Tumour_Carriers", "Normal_N", "Tumour_N", "Normal_Frequency_pct",
        "Tumour_Frequency_pct", "Tumour_minus_Normal_pct", "p_value", "FDR_BH",
    ]

    for map_label, frame in [("GSDMB (script 9)", gsdmb_tests), ("Global (script 8)", global_tests)]:
        for cohort, sub in frame.groupby("Cohort"):
            nominal_hits = sub[sub["Raw_p_lt_0.05"]].sort_values(["p_value", "FDR_BH", "Pos"])
            fdr_hits = sub[sub["FDR_lt_0.05"]].sort_values(["FDR_BH", "p_value", "Pos"])
            counts = {
                "Map": map_label,
                "Cohort": cohort,
                "Variants_tested": len(sub),
                "Raw_p_lt_0.05": int(sub["Raw_p_lt_0.05"].sum()),
                "FDR_lt_0.05": int(sub["FDR_lt_0.05"].sum()),
                "Shared_variants": int((sub["Variant_Status"] == "shared").sum()),
                "Tumour_only_variants": int((sub["Variant_Status"] == "tumour_only").sum()),
                "Normal_only_variants": int((sub["Variant_Status"] == "normal_only").sum()),
                "Nominal_significant_variants": "; ".join(label_significant_row(row) for _, row in nominal_hits.iterrows()) or "None",
                "FDR_significant_variants": "; ".join(label_significant_row(row) for _, row in fdr_hits.iterrows()) or "None",
            }
            significance_counts_rows.append(counts)

            if nominal_hits.empty:
                significance_detail_rows.append(
                    {
                        **{key: counts[key] for key in ["Map", "Cohort", "Variants_tested", "Raw_p_lt_0.05", "FDR_lt_0.05"]},
                        "Significance_Level": "None",
                        "Also_FDR_significant": False,
                        "Gene_Symbol": None,
                        "Gene_ID": None,
                        "rsID": None,
                        "Pos": None,
                        "Impact": None,
                        "Consequence": None,
                        "Location_Class": None,
                        "Feature_type": None,
                        "EXON": None,
                        "INTRON": None,
                        "Variant_Status": None,
                        "Direction": None,
                        "Normal_Carriers": None,
                        "Tumour_Carriers": None,
                        "Normal_N": None,
                        "Tumour_N": None,
                        "Normal_Frequency_pct": None,
                        "Tumour_Frequency_pct": None,
                        "Tumour_minus_Normal_pct": None,
                        "p_value": None,
                        "FDR_BH": None,
                    }
                )
                continue

            for _, row in nominal_hits.iterrows():
                detail = {
                    **{key: counts[key] for key in ["Map", "Cohort", "Variants_tested", "Raw_p_lt_0.05", "FDR_lt_0.05"]},
                    "Significance_Level": "Nominal (p < 0.05)",
                    "Also_FDR_significant": bool(row["FDR_lt_0.05"]),
                }
                for col in detail_columns:
                    if col in detail:
                        continue
                    detail[col] = row[col] if col in row.index else None
                significance_detail_rows.append(detail)

    significance_counts = pd.DataFrame(significance_counts_rows)
    significance_summary = pd.DataFrame(significance_detail_rows)
    significance_summary = significance_summary[detail_columns]

    qc_counts = pd.concat(
        [gsdmb_qc_counts.assign(Map="GSDMB (script 9)"), global_qc_counts.assign(Map="Global (script 8)")],
        ignore_index=True,
    )

    guide = pd.DataFrame(
        [
            ["Purpose", "Summary of the visual patterns seen in scripts 8 and 9, with supporting variant counts, frequency tests, QC sensitivity checks, and a detailed significant-variant overview."],
            ["How to read frequencies", "Frequencies are carrier frequencies (% of samples carrying the variant) within each cohort/tissue panel."],
            ["Variant status", "shared = present in both tissues; tumour_only = seen only in tumour; normal_only = seen only in normal."],
            ["Direction", "Tumour_minus_Normal_pct < 0 means the variant is less frequent in tumour than in normal."],
            ["QC sensitivity", "High-QC comparisons use only aligned samples with zero_amplicons = 0 in the coverage audit generated by script 03."],
            ["Where to find p/FDR", "Open Significance_Overview for the per-variant significant results with gene/consequence/location details, Significance_Counts for the map-level counts, and the *_with_p_FDR sheets for the full tested tables."],
            ["Frequency significance", "p_value is a two-sided Fisher exact test on carrier vs non-carrier counts within each cohort. FDR_BH is Benjamini-Hochberg correction within each cohort and map."],
        ]
    )

    gsdmb_key = gsdmb_tests.sort_values(["Cohort", "p_value", "Tumour_minus_Normal_pct"]).groupby("Cohort", group_keys=False).head(12)
    global_key = global_tests.sort_values(["Cohort", "p_value", "Tumour_minus_Normal_pct"]).groupby("Cohort", group_keys=False).head(12)

    with pd.ExcelWriter(WORKBOOK_PATH, engine="openpyxl") as writer:
        guide.to_excel(writer, sheet_name="Guide", index=False, header=False)
        summary.to_excel(writer, sheet_name="Summary", index=False)
        significance_summary.to_excel(writer, sheet_name="Significance_Overview", index=False)
        significance_counts.to_excel(writer, sheet_name="Significance_Counts", index=False)
        qc_counts.to_excel(writer, sheet_name="QC_Sample_Counts", index=False)
        gsdmb_status.to_excel(writer, sheet_name="GSDMB_Variant_Status", index=False)
        gsdmb_direction.to_excel(writer, sheet_name="GSDMB_Shared_Direction", index=False)
        gsdmb_key.to_excel(writer, sheet_name="GSDMB_Key_with_p_FDR", index=False)
        gsdmb_tests.to_excel(writer, sheet_name="GSDMB_All_with_p_FDR", index=False)
        gsdmb_qc_direction.to_excel(writer, sheet_name="GSDMB_QC_Sensitivity", index=False)
        global_status.to_excel(writer, sheet_name="Global_Variant_Status", index=False)
        global_direction.to_excel(writer, sheet_name="Global_Shared_Direction", index=False)
        global_key.to_excel(writer, sheet_name="Global_Key_with_p_FDR", index=False)
        global_tests.to_excel(writer, sheet_name="Global_All_with_p_FDR", index=False)
        global_qc_direction.to_excel(writer, sheet_name="Global_QC_Sensitivity", index=False)

    style_workbook(WORKBOOK_PATH)

    summary.to_csv(SUMMARY_CSV, index=False)
    gsdmb_tests.to_csv(GSDMB_TESTS_CSV, index=False)
    global_tests.to_csv(GLOBAL_TESTS_CSV, index=False)
    significance_summary.to_csv(SIGNIFICANCE_SUMMARY_CSV, index=False)

    print("Saved workbook:", WORKBOOK_PATH)
    print("Saved summary CSV:", SUMMARY_CSV)
    print("Saved GSDMB tests CSV:", GSDMB_TESTS_CSV)
    print("Saved global tests CSV:", GLOBAL_TESTS_CSV)
    print("Saved significance summary CSV:", SIGNIFICANCE_SUMMARY_CSV)


if __name__ == "__main__":
    main()
