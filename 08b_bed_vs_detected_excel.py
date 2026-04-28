"""
Script 08b: Designed BED Versus Detected Variant Comparison
===========================================================

Purpose
-------
Create a supervisor-friendly workbook that compares the original panel design
BED against the variants actually detected downstream, while separating true
design-linked findings from extra calls not explicitly represented in the BED
annotations.

Methodological rationale
------------------------
The primary truth for "in design" is interval overlap with the BED rather than
rsID alone. rsIDs are still used when present because they improve traceability
and can rescue design-linked variants when BED naming is more informative than
strict interval annotation.

Inputs
------
- Designed BED file, with optional annotation columns such as amplicon ID,
  rsID and gene.
- Annotated Excel workbook, typically the ``Biological_Annotations`` sheet
  generated downstream of variant calling and VEP annotation.

Outputs
-------
- Multi-sheet Excel workbook with:
  - ``Detected_vs_BED``
  - ``BED_targets_summary``
  - ``Extra_detected_only``
  - ``Summary``
- Optional CSV export of the main variant-centric sheet.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from pipeline_utils import common_nfe_variant_mask, ensure_directory, get_paths, get_thresholds, get_warning_amplicons, normalise_chromosome_label
from pipeline_validation import validate_file_exists

PATHS = get_paths()
THRESHOLDS = get_thresholds()
MIN_NFE_AF = float(THRESHOLDS.get("min_nfe_af", 0.01))
DEFAULT_ANNOTATED_XLSX = PATHS.get(
    "annotated_report",
    Path("/home/gadeaalonsoj/tfm/analysis_results/07_annotated_variants/GSDMB_Annotated_Variants.xlsx"),
)
DEFAULT_OUTPUT_DIR = ensure_directory(PATHS.get("results_dir", Path.cwd() / "analysis_results") / "08b_bed_vs_detected")
DEFAULT_OUTPUT_XLSX = DEFAULT_OUTPUT_DIR / "GSDMB_BED_vs_Detected_Comparison.xlsx"
DEFAULT_TECHNICAL_AUDIT_CSV = PATHS.get("technical_audit_dir", Path("/home/gadeaalonsoj/tfm/analysis_results/04_technical_audit")) / "GSDMB_Technical_Audit.csv"
DEFAULT_SHEET = "Biological_Annotations"

HEADER_FILL = PatternFill(fill_type="solid", fgColor="1F4E78")
HEADER_FONT = Font(color="FFFFFF", bold=True)
GREEN_FILL = PatternFill(fill_type="solid", fgColor="E2F0D9")
ORANGE_FILL = PatternFill(fill_type="solid", fgColor="FCE4D6")
GREY_FILL = PatternFill(fill_type="solid", fgColor="E7E6E6")
SUMMARY_LABEL_FILL = PatternFill(fill_type="solid", fgColor="D9EAF7")
MISSING_TOKENS = {"", ".", "na", "n/a", "nan", "none", "null", "-"}
ANNOTATION_SOURCE_NOTE = "Stage-07 Biological_Annotations sheet (merged VEP annotations from script 07)."
IMPACT_PRIORITY = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
CONSEQUENCE_PRIORITY = {
    "transcript_ablation": 0,
    "splice_acceptor_variant": 1,
    "splice_donor_variant": 2,
    "stop_gained": 3,
    "frameshift_variant": 4,
    "stop_lost": 5,
    "start_lost": 6,
    "transcript_amplification": 7,
    "inframe_insertion": 8,
    "inframe_deletion": 9,
    "missense_variant": 10,
    "protein_altering_variant": 11,
    "splice_donor_5th_base_variant": 12,
    "splice_region_variant": 13,
    "splice_polypyrimidine_tract_variant": 14,
    "synonymous_variant": 15,
    "start_retained_variant": 16,
    "stop_retained_variant": 17,
    "5_prime_UTR_variant": 18,
    "3_prime_UTR_variant": 19,
    "non_coding_transcript_exon_variant": 20,
    "intron_variant": 21,
    "upstream_gene_variant": 22,
    "downstream_gene_variant": 23,
    "regulatory_region_variant": 24,
    "intergenic_variant": 25,
}
FUNCTIONAL_CLASS_MAP = {
    "synonymous_variant": "synonymous",
    "missense_variant": "missense",
    "stop_gained": "nonsense",
    "frameshift_variant": "frameshift",
    "splice_acceptor_variant": "splice_site",
    "splice_donor_variant": "splice_site",
    "splice_region_variant": "splice_site",
    "splice_donor_5th_base_variant": "splice_site",
    "splice_polypyrimidine_tract_variant": "splice_site",
    "intron_variant": "intronic",
    "upstream_gene_variant": "regulatory",
    "downstream_gene_variant": "regulatory",
    "regulatory_region_variant": "regulatory",
    "intergenic_variant": "intergenic",
    "5_prime_UTR_variant": "utr",
    "3_prime_UTR_variant": "utr",
}
INFERRED_IMPACT_MAP = {
    "transcript_ablation": "HIGH",
    "splice_acceptor_variant": "HIGH",
    "splice_donor_variant": "HIGH",
    "stop_gained": "HIGH",
    "frameshift_variant": "HIGH",
    "stop_lost": "HIGH",
    "start_lost": "HIGH",
    "transcript_amplification": "HIGH",
    "inframe_insertion": "MODERATE",
    "inframe_deletion": "MODERATE",
    "missense_variant": "MODERATE",
    "protein_altering_variant": "MODERATE",
    "splice_donor_5th_base_variant": "LOW",
    "splice_region_variant": "LOW",
    "splice_polypyrimidine_tract_variant": "LOW",
    "synonymous_variant": "LOW",
    "start_retained_variant": "LOW",
    "stop_retained_variant": "LOW",
    "5_prime_UTR_variant": "MODIFIER",
    "3_prime_UTR_variant": "MODIFIER",
    "non_coding_transcript_exon_variant": "MODIFIER",
    "intron_variant": "MODIFIER",
    "upstream_gene_variant": "MODIFIER",
    "downstream_gene_variant": "MODIFIER",
    "regulatory_region_variant": "MODIFIER",
    "intergenic_variant": "MODIFIER",
}


def normalise_yes_flag(value: Any) -> bool:
    """Interpret common YES/true-like transcript flags consistently."""
    return clean_text(value).upper() in {"YES", "Y", "TRUE", "1"}


def split_consequence_terms(value: Any) -> list[str]:
    """Split compound VEP consequence strings while preserving term order."""
    text = clean_text(value)
    if not text:
        return []
    terms: list[str] = []
    seen: set[str] = set()
    for raw_term in re.split(r"[,&]", text):
        term = raw_term.strip()
        if term and term not in seen:
            seen.add(term)
            terms.append(term)
    return terms


def select_worst_consequence_term(value: Any) -> str:
    """Return the highest-priority VEP consequence term from a compound annotation."""
    terms = split_consequence_terms(value)
    if not terms:
        return ""
    return min(terms, key=lambda term: CONSEQUENCE_PRIORITY.get(term, 999))


def map_functional_class(consequence_term: str) -> str:
    """Map detailed VEP consequences to a compact, report-friendly functional class."""
    term = clean_text(consequence_term)
    if not term:
        return ""
    return FUNCTIONAL_CLASS_MAP.get(term, "other")


def infer_impact_from_consequence(consequence_term: str) -> str:
    """Infer a VEP-style impact tier when IMPACT is missing from the input row."""
    term = clean_text(consequence_term)
    if not term:
        return ""
    return INFERRED_IMPACT_MAP.get(term, "MODIFIER")


def parse_sift_prediction(raw_value: Any, fallback_value: Any) -> str:
    """Normalise SIFT strings such as deleterious(0.02) or single-letter fallback codes."""
    raw_text = clean_text(raw_value).lower()
    if raw_text:
        match = re.search(r"(deleterious|tolerated)", raw_text)
        if match:
            return match.group(1)
    fallback_tokens = re.findall(r"[A-Za-z]+", clean_text(fallback_value).upper())
    for token in fallback_tokens:
        if token == "D":
            return "deleterious"
        if token == "T":
            return "tolerated"
    return ""


def parse_polyphen_prediction(raw_value: Any, fallback_value: Any) -> str:
    """Normalise PolyPhen strings such as probably_damaging(0.92) or B/P/D fallback codes."""
    raw_text = clean_text(raw_value).lower()
    if raw_text:
        for label in ["probably_damaging", "possibly_damaging", "benign"]:
            if label in raw_text:
                return label
    fallback_tokens = re.findall(r"[A-Za-z]+", clean_text(fallback_value).upper())
    for token in fallback_tokens:
        if token == "D":
            return "probably_damaging"
        if token == "P":
            return "possibly_damaging"
        if token == "B":
            return "benign"
    return ""


def describe_annotation_basis(row: pd.Series) -> str:
    """Explain why a transcript row became the representative row for one variant."""
    if normalise_yes_flag(row.get("PICK", "")):
        return "picked_transcript"
    if normalise_yes_flag(row.get("CANONICAL", "")):
        return "canonical_transcript"
    if clean_text(row.get("MANE_SELECT", "")):
        return "mane_transcript"
    return "worst_available_transcript"


def choose_representative_annotation_row(group: pd.DataFrame) -> pd.Series:
    """Choose one transcript row per variant using explicit, reproducible rules.

    The stage-07 workbook is transcript-expanded. For the new detailed consequence
    fields we select one representative row by prioritising VEP/PICK when present,
    then CANONICAL, then MANE, and finally worst available IMPACT/consequence.
    This keeps the rule deterministic while still favouring biologically severe
    annotations when several transcripts remain tied.
    """
    ranked = group.copy()
    ranked["_pick_rank"] = ranked.get("PICK", pd.Series(index=ranked.index, dtype=object)).map(normalise_yes_flag).map({True: 0, False: 1})
    ranked["_canonical_rank"] = ranked.get("CANONICAL", pd.Series(index=ranked.index, dtype=object)).map(normalise_yes_flag).map({True: 0, False: 1})
    ranked["_mane_rank"] = ranked.get("MANE_SELECT", pd.Series(index=ranked.index, dtype=object)).map(lambda value: 0 if clean_text(value) else 1)
    ranked["_impact_rank"] = ranked.get("IMPACT", pd.Series(index=ranked.index, dtype=object)).map(lambda value: IMPACT_PRIORITY.get(clean_text(value).upper(), 99))
    ranked["_worst_consequence_term"] = ranked.get("Consequence", pd.Series(index=ranked.index, dtype=object)).map(select_worst_consequence_term)
    ranked["_consequence_rank"] = ranked["_worst_consequence_term"].map(lambda value: CONSEQUENCE_PRIORITY.get(clean_text(value), 999))
    ranked = ranked.sort_values(
        ["_pick_rank", "_canonical_rank", "_mane_rank", "_impact_rank", "_consequence_rank"],
        ascending=[True, True, True, True, True],
        kind="stable",
    )
    return ranked.iloc[0]


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the command-line interface for the reporting step."""
    parser = argparse.ArgumentParser(
        description="Compare designed BED targets against detected annotated variants and export a formatted Excel workbook."
    )
    parser.add_argument("--bed", required=True, help="Path to the designed BED file.")
    parser.add_argument(
        "--annotated_xlsx",
        default=str(DEFAULT_ANNOTATED_XLSX),
        help=f"Path to the annotated Excel workbook (default: {DEFAULT_ANNOTATED_XLSX}).",
    )
    parser.add_argument(
        "--sheet",
        default=DEFAULT_SHEET,
        help=f"Sheet name containing detected variants (default: {DEFAULT_SHEET}).",
    )
    parser.add_argument(
        "--technical_audit_csv",
        default=str(DEFAULT_TECHNICAL_AUDIT_CSV),
        help=f"Optional technical audit CSV used to annotate amplicon quality (default: {DEFAULT_TECHNICAL_AUDIT_CSV}).",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_XLSX),
        help=f"Output Excel workbook path (default: {DEFAULT_OUTPUT_XLSX}).",
    )
    parser.add_argument(
        "--csv_output",
        default="",
        help="Optional CSV export path for the main Detected_vs_BED sheet.",
    )
    return parser


def is_missing_like(value: Any) -> bool:
    """Return True when a scalar is effectively blank for reporting purposes."""
    if pd.isna(value):
        return True
    text = str(value).strip()
    return text.lower() in MISSING_TOKENS


def clean_text(value: Any) -> str:
    """Normalise free-text cells to stripped strings with blanks for missing values."""
    return "" if is_missing_like(value) else str(value).strip()


def normalise_header_key(value: Any) -> str:
    """Build a forgiving header-comparison key."""
    return re.sub(r"[^A-Z0-9]+", "", str(value).strip().upper())


def find_column_case_insensitive(
    df: pd.DataFrame,
    candidates: Iterable[str],
    *,
    required: bool = False,
    context: str = "",
) -> str | None:
    """Find the first matching column using punctuation-insensitive comparison."""
    lookup = {normalise_header_key(col): col for col in df.columns}
    for candidate in candidates:
        resolved = lookup.get(normalise_header_key(candidate))
        if resolved is not None:
            return resolved
    if required:
        candidate_text = ", ".join(candidates)
        raise ValueError(f"{context}: could not find any of the expected columns -> {candidate_text}")
    return None


def first_nonempty(series: pd.Series) -> Any:
    """Return the first non-empty value from a series, else blank."""
    for value in series:
        if not is_missing_like(value):
            return value
    return ""


def extract_all_rsids(value: Any) -> list[str]:
    """Extract all rsIDs from a semi-structured annotation field."""
    text = clean_text(value)
    if not text:
        return []
    return sorted({match.lower() for match in re.findall(r"(rs\d+)", text, flags=re.IGNORECASE)})


def extract_first_rsid(value: Any) -> str:
    """Extract the first rsID from a scalar value, if present."""
    rsids = extract_all_rsids(value)
    return rsids[0] if rsids else ""


def join_unique(values: Iterable[Any], *, max_items: int | None = None) -> str:
    """Collapse values into a unique, semicolon-delimited string."""
    unique_values = sorted({clean_text(value) for value in values if not is_missing_like(value)})
    if max_items is not None:
        unique_values = unique_values[:max_items]
    return "; ".join(unique_values)


def flatten_unique(list_values: Iterable[Iterable[Any]]) -> list[str]:
    """Flatten a list of iterables and return unique non-empty text values."""
    flattened: set[str] = set()
    for values in list_values:
        for value in values:
            text = clean_text(value)
            if text:
                flattened.add(text)
    return sorted(flattened)


def summarise_examples(values: Iterable[Any], max_items: int = 5) -> str:
    """Return a compact list of example values for Excel reporting."""
    unique_values = sorted({clean_text(value) for value in values if not is_missing_like(value)})
    if not unique_values:
        return ""
    if len(unique_values) <= max_items:
        return "; ".join(unique_values)
    shown = "; ".join(unique_values[:max_items])
    return f"{shown}; +{len(unique_values) - max_items} more"


def chrom_sort_key(value: Any) -> tuple[int, str]:
    """Provide a stable human-friendly chromosome sort order."""
    chrom = clean_text(value)
    if not chrom:
        return (99, "")
    chrom = normalise_chromosome_label(chrom)
    suffix = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if suffix.isdigit():
        return (int(suffix), chrom)
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (special.get(suffix.upper(), 98), chrom)


def build_variant_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    """Build a stable variant identifier for deduplication and matching."""
    return f"{clean_text(chrom)}:{int(pos)}:{clean_text(ref)}>{clean_text(alt)}"


def looks_like_bed_header(fields: list[str]) -> bool:
    """Heuristically detect a BED header row after skipping track/browser lines."""
    normalised = [normalise_header_key(value) for value in fields]
    header_terms = {
        "CHR", "CHROM", "CHROMOSOME", "START", "END", "AMPLICONID",
        "AMPLICON", "RSID", "GENE", "NAME", "TARGET",
    }
    hits = sum(field in header_terms for field in normalised)
    start_is_header = len(normalised) > 1 and normalised[1] == "START"
    end_is_header = len(normalised) > 2 and normalised[2] == "END"
    return hits >= 2 or start_is_header or end_is_header


def build_default_bed_columns(n_columns: int) -> list[str]:
    """Create a default BED column layout when the file has no explicit header."""
    base = ["chr", "start", "end", "amplicon_id", "annotation_name", "gene"]
    if n_columns <= len(base):
        return base[:n_columns]
    extras = [f"extra_col_{idx}" for idx in range(1, n_columns - len(base) + 1)]
    return base + extras


def create_bed_label(row: pd.Series) -> str:
    """Create a clean human-readable BED annotation label."""
    for key in ["rsid", "gene", "annotation_name", "amplicon_id"]:
        value = clean_text(row.get(key, ""))
        if value:
            return value
    return f"{row['chr']}:{row['start']}-{row['end']}"


def load_bed_file(bed_path: Path) -> pd.DataFrame:
    """Read a BED file while tolerating track/browser/header lines and mixed annotations."""
    rows: list[list[str]] = []
    with bed_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped or stripped.startswith(("track", "browser", "#")):
                continue
            rows.append(stripped.split("\t"))

    if not rows:
        raise ValueError(f"BED file contains no usable target rows -> {bed_path}")

    max_columns = max(len(row) for row in rows)
    if looks_like_bed_header(rows[0]):
        header = [clean_text(value) or f"col_{idx + 1}" for idx, value in enumerate(rows[0])]
        if len(header) < max_columns:
            header.extend([f"extra_col_{idx}" for idx in range(1, max_columns - len(header) + 1)])
        data_rows = [row + [""] * (max_columns - len(row)) for row in rows[1:]]
    else:
        header = build_default_bed_columns(max_columns)
        data_rows = [row + [""] * (max_columns - len(row)) for row in rows]

    bed_raw = pd.DataFrame(data_rows, columns=header)
    if bed_raw.empty:
        raise ValueError(f"BED file contains a header but no target rows -> {bed_path}")

    chrom_col = find_column_case_insensitive(bed_raw, ["chr", "chrom", "chromosome"], required=True, context="BED")
    start_col = find_column_case_insensitive(bed_raw, ["start"], required=True, context="BED")
    end_col = find_column_case_insensitive(bed_raw, ["end", "stop"], required=True, context="BED")
    amplicon_col = find_column_case_insensitive(
        bed_raw,
        ["amplicon_id", "amplicon", "ampliconid", "target_id", "id", "name"],
        context="BED",
    )
    annotation_col = find_column_case_insensitive(
        bed_raw,
        ["rsid", "annotation_name", "target_name", "variant_name", "name"],
        context="BED",
    )
    gene_col = find_column_case_insensitive(bed_raw, ["gene", "symbol", "gene_name"], context="BED")

    bed_df = pd.DataFrame({
        "chr": bed_raw[chrom_col].map(normalise_chromosome_label),
        "start": pd.to_numeric(bed_raw[start_col], errors="coerce"),
        "end": pd.to_numeric(bed_raw[end_col], errors="coerce"),
        "amplicon_id": bed_raw[amplicon_col].map(clean_text) if amplicon_col else "",
        "annotation_name": bed_raw[annotation_col].map(clean_text) if annotation_col else "",
        "gene": bed_raw[gene_col].map(clean_text) if gene_col else "",
    })
    bed_df = bed_df.dropna(subset=["start", "end"]).copy()
    bed_df["start"] = bed_df["start"].astype(int)
    bed_df["end"] = bed_df["end"].astype(int)
    bed_df = bed_df[(bed_df["chr"] != "") & (bed_df["end"] >= bed_df["start"])].copy()
    if bed_df.empty:
        raise ValueError(f"No valid BED intervals remained after parsing -> {bed_path}")

    bed_df["rsid"] = bed_df["annotation_name"].apply(extract_first_rsid)
    bed_df["BED_annotation_label"] = bed_df.apply(create_bed_label, axis=1)
    bed_df["bed_start_1based"] = bed_df["start"] + 1
    bed_df["bed_end_1based"] = bed_df["end"]
    bed_df["BED_Target_ID"] = range(1, len(bed_df) + 1)
    bed_df = bed_df.sort_values(["chr", "start", "end", "BED_Target_ID"], key=lambda col: col.map(chrom_sort_key) if col.name == "chr" else col)
    return bed_df.reset_index(drop=True)



def load_technical_audit(technical_audit_csv: Path | None) -> pd.DataFrame:
    """Load amplicon-level technical audit metrics when available."""
    out_cols = [
        "BED_Target_ID", "Selected_Worst_Amplicon", "Global_Failure_%",
        "Global_Mean_Depth", "Global_Median_Depth",
        "Breast_Control_fail_%", "Breast_Tumour_fail_%",
        "Endometrium_Control_fail_%", "Endometrium_Tumour_fail_%",
    ]
    if technical_audit_csv is None or not technical_audit_csv.exists():
        return pd.DataFrame(columns=out_cols)

    audit_raw = pd.read_csv(technical_audit_csv)
    chrom_col = find_column_case_insensitive(audit_raw, ["chr", "chrom", "chromosome"], required=True, context="Technical audit")
    start_col = find_column_case_insensitive(audit_raw, ["start"], required=True, context="Technical audit")
    end_col = find_column_case_insensitive(audit_raw, ["end"], required=True, context="Technical audit")
    numeric_id_col = find_column_case_insensitive(audit_raw, ["numeric_id", "amplicon_id", "id"], context="Technical audit")

    audit = pd.DataFrame({
        "chr": audit_raw[chrom_col].map(normalise_chromosome_label),
        "start": pd.to_numeric(audit_raw[start_col], errors="coerce").astype("Int64"),
        "end": pd.to_numeric(audit_raw[end_col], errors="coerce").astype("Int64"),
        "amplicon_id": audit_raw[numeric_id_col].map(clean_text) if numeric_id_col else "",
        "Global_Failure_%": pd.to_numeric(audit_raw.get("Global_Failure_%"), errors="coerce"),
        "Global_Mean_Depth": pd.to_numeric(audit_raw.get("Global_Mean_Depth"), errors="coerce"),
        "Global_Median_Depth": pd.to_numeric(audit_raw.get("Global_Median_Depth"), errors="coerce"),
        "Breast_Control_fail_%": pd.to_numeric(audit_raw.get("Breast_Control_fail_%"), errors="coerce"),
        "Breast_Tumour_fail_%": pd.to_numeric(audit_raw.get("Breast_Tumour_fail_%"), errors="coerce"),
        "Endometrium_Control_fail_%": pd.to_numeric(audit_raw.get("Endometrium_Control_fail_%"), errors="coerce"),
        "Endometrium_Tumour_fail_%": pd.to_numeric(audit_raw.get("Endometrium_Tumour_fail_%"), errors="coerce"),
    }).dropna(subset=["start", "end"]).copy()
    audit["start"] = audit["start"].astype(int)
    audit["end"] = audit["end"].astype(int)

    warning_amplicons = get_warning_amplicons()
    warning_ids = set()
    for row in warning_amplicons.itertuples(index=False):
        overlaps = audit[
            (audit["chr"] == row.chrom)
            & (audit["start"] <= int(row.end))
            & (audit["end"] >= int(row.start))
        ]
        warning_ids.update(overlaps["amplicon_id"].astype(str).tolist())

    audit["Selected_Worst_Amplicon"] = audit["amplicon_id"].astype(str).isin(warning_ids).map({True: "YES", False: "NO"})
    return audit[[
        "chr", "start", "end", "amplicon_id", "Selected_Worst_Amplicon",
        "Global_Failure_%", "Global_Mean_Depth", "Global_Median_Depth",
        "Breast_Control_fail_%", "Breast_Tumour_fail_%",
        "Endometrium_Control_fail_%", "Endometrium_Tumour_fail_%",
    ]].drop_duplicates()

def load_detected_annotations(annotated_xlsx: Path, sheet_name: str) -> pd.DataFrame:
    """Load and standardise the stage-07 VEP-derived annotation workbook.

    Script 08b does not re-run VEP itself. Instead it consumes the consolidated
    ``Biological_Annotations`` sheet produced by script 07, which already carries
    transcript-level VEP fields such as ``Consequence``, ``IMPACT``, ``CANONICAL``,
    ``MANE_SELECT``, ``SIFT`` and ``PolyPhen``.
    """
    detected_raw = pd.read_excel(annotated_xlsx, sheet_name=sheet_name)
    chrom_col = find_column_case_insensitive(detected_raw, ["CHROM", "chr", "chrom"], required=True, context="Detected variants")
    pos_col = find_column_case_insensitive(detected_raw, ["POS", "position", "start"], required=True, context="Detected variants")
    ref_col = find_column_case_insensitive(detected_raw, ["REF", "reference"], required=True, context="Detected variants")
    alt_col = find_column_case_insensitive(detected_raw, ["ALT", "alternate"], required=True, context="Detected variants")

    optional_cols = {
        "Existing_variation": find_column_case_insensitive(detected_raw, ["Existing_variation", "existing_variation", "rsid"], context="Detected variants"),
        "SYMBOL": find_column_case_insensitive(detected_raw, ["SYMBOL", "gene", "symbol"], context="Detected variants"),
        "Consequence": find_column_case_insensitive(detected_raw, ["Consequence"], context="Detected variants"),
        "IMPACT": find_column_case_insensitive(detected_raw, ["IMPACT", "Impact"], context="Detected variants"),
        "HGVSp": find_column_case_insensitive(detected_raw, ["HGVSp", "Protein_change"], context="Detected variants"),
        "Feature": find_column_case_insensitive(detected_raw, ["Feature", "Transcript", "feature"], context="Detected variants"),
        "CANONICAL": find_column_case_insensitive(detected_raw, ["CANONICAL"], context="Detected variants"),
        "MANE_SELECT": find_column_case_insensitive(detected_raw, ["MANE_SELECT"], context="Detected variants"),
        "PICK": find_column_case_insensitive(detected_raw, ["PICK"], context="Detected variants"),
        "SIFT": find_column_case_insensitive(detected_raw, ["SIFT"], context="Detected variants"),
        "PolyPhen": find_column_case_insensitive(detected_raw, ["PolyPhen"], context="Detected variants"),
        "SIFT_pred": find_column_case_insensitive(detected_raw, ["SIFT_pred"], context="Detected variants"),
        "Polyphen2_HDIV_pred": find_column_case_insensitive(detected_raw, ["Polyphen2_HDIV_pred"], context="Detected variants"),
        "Sample": find_column_case_insensitive(detected_raw, ["Sample", "sample_id"], context="Detected variants"),
        "Cohort": find_column_case_insensitive(detected_raw, ["Cohort", "cohort"], context="Detected variants"),
        "Tissue": find_column_case_insensitive(detected_raw, ["Tissue", "tissue"], context="Detected variants"),
        "DP": find_column_case_insensitive(detected_raw, ["DP", "Depth"], context="Detected variants"),
        "AF": find_column_case_insensitive(detected_raw, ["AF", "VAF", "Allele_fraction"], context="Detected variants"),
        "gnomADe_NFE_AF": find_column_case_insensitive(detected_raw, ["gnomADe_NFE_AF"], context="Detected variants"),
        "gnomADg_NFE_AF": find_column_case_insensitive(detected_raw, ["gnomADg_NFE_AF"], context="Detected variants"),
    }

    detected = pd.DataFrame({
        "CHROM": detected_raw[chrom_col].map(normalise_chromosome_label),
        "POS": pd.to_numeric(detected_raw[pos_col], errors="coerce"),
        "REF": detected_raw[ref_col].map(lambda value: clean_text(value).upper()),
        "ALT": detected_raw[alt_col].map(lambda value: clean_text(value).upper()),
    })
    for out_col, source_col in optional_cols.items():
        if source_col is None:
            detected[out_col] = pd.NA
        else:
            detected[out_col] = detected_raw[source_col]

    detected = detected.dropna(subset=["POS"]).copy()
    detected["POS"] = detected["POS"].astype(int)
    detected = detected[(detected["CHROM"] != "") & (detected["REF"] != "") & (detected["ALT"] != "")].copy()
    if detected.empty:
        raise ValueError("No valid detected variant rows remained after standardising the annotated workbook.")

    detected["Variant_Key"] = detected.apply(
        lambda row: build_variant_key(row["CHROM"], row["POS"], row["REF"], row["ALT"]),
        axis=1,
    )
    detected["Detected_rsIDs"] = detected["Existing_variation"].apply(extract_all_rsids)
    detected["Sample"] = detected["Sample"].map(clean_text) if "Sample" in detected.columns else ""
    detected["Cohort"] = detected["Cohort"].map(clean_text) if "Cohort" in detected.columns else ""
    detected["Tissue"] = detected["Tissue"].map(clean_text) if "Tissue" in detected.columns else ""
    detected["DP"] = pd.to_numeric(detected["DP"], errors="coerce")
    detected["AF"] = pd.to_numeric(detected["AF"], errors="coerce")
    detected["gnomADe_NFE_AF"] = pd.to_numeric(detected["gnomADe_NFE_AF"], errors="coerce")
    detected["gnomADg_NFE_AF"] = pd.to_numeric(detected["gnomADg_NFE_AF"], errors="coerce")
    detected["gnomAD_NFE_AF_combined"] = detected["gnomADe_NFE_AF"].combine_first(detected["gnomADg_NFE_AF"])
    detected["Common_SNP_By_NFE"] = common_nfe_variant_mask(detected, MIN_NFE_AF)
    return detected

def collapse_variant_annotations(group: pd.DataFrame) -> dict[str, Any]:
    """Collapse transcript-expanded rows into one biologically useful variant record.

    Existing script-08b columns still summarise all transcript values across the
    grouped rows. The new appended fields instead describe one representative
    transcript chosen by ``choose_representative_annotation_row`` so each variant
    gets a deterministic detailed consequence, simplified functional class,
    impact tier, and parsed predictor labels.
    """
    has_sample_col = "Sample" in group.columns
    sample_level = group.copy()
    if has_sample_col and sample_level["Sample"].replace("", pd.NA).notna().any():
        sample_level = (
            sample_level.groupby(["Variant_Key", "Sample"], dropna=False, sort=False)
            .agg({
                "Cohort": lambda values: join_unique(values),
                "Tissue": lambda values: join_unique(values),
                "DP": "mean",
                "AF": "mean",
            })
            .reset_index()
        )
        sample_values = [value for value in sample_level["Sample"] if clean_text(value)]
        n_samples = len(set(sample_values))
        sample_examples = summarise_examples(sample_values)
        cohorts = join_unique(sample_level["Cohort"])
        tissues = join_unique(sample_level["Tissue"])
        mean_dp = sample_level["DP"].mean() if sample_level["DP"].notna().any() else pd.NA
        mean_af = sample_level["AF"].mean() if sample_level["AF"].notna().any() else pd.NA
    else:
        n_samples = 0
        sample_examples = ""
        cohorts = join_unique(group["Cohort"]) if "Cohort" in group.columns else ""
        tissues = join_unique(group["Tissue"]) if "Tissue" in group.columns else ""
        mean_dp = group["DP"].mean() if group["DP"].notna().any() else pd.NA
        mean_af = group["AF"].mean() if group["AF"].notna().any() else pd.NA

    representative = choose_representative_annotation_row(group)
    detailed_consequence = clean_text(representative.get("Consequence", ""))
    worst_consequence_term = select_worst_consequence_term(detailed_consequence)
    impact_severity = clean_text(representative.get("IMPACT", "")).upper()
    if not impact_severity:
        impact_severity = infer_impact_from_consequence(worst_consequence_term)

    sift_raw = clean_text(representative.get("SIFT", ""))
    polyphen_raw = clean_text(representative.get("PolyPhen", ""))
    sift_prediction = parse_sift_prediction(sift_raw, representative.get("SIFT_pred", ""))
    polyphen_prediction = parse_polyphen_prediction(polyphen_raw, representative.get("Polyphen2_HDIV_pred", ""))

    detected_rsids = flatten_unique(group["Detected_rsIDs"])
    gnomad_nfe = group["gnomAD_NFE_AF_combined"].dropna().max() if group["gnomAD_NFE_AF_combined"].notna().any() else pd.NA
    common_snp_by_nfe = bool(pd.notna(gnomad_nfe) and float(gnomad_nfe) >= MIN_NFE_AF)
    return {
        "CHROM": first_nonempty(group["CHROM"]),
        "POS": int(first_nonempty(group["POS"])),
        "REF": first_nonempty(group["REF"]),
        "ALT": first_nonempty(group["ALT"]),
        "Variant_Key": first_nonempty(group["Variant_Key"]),
        "rsID_detected": "; ".join(detected_rsids),
        "SYMBOL": join_unique(group["SYMBOL"]),
        "Consequence": join_unique(group["Consequence"]),
        "IMPACT": join_unique(group["IMPACT"]),
        "HGVSp": join_unique(group["HGVSp"]),
        "N_samples": n_samples,
        "Sample_examples": sample_examples,
        "Cohorts": cohorts,
        "Tissues": tissues,
        "Mean_DP": mean_dp,
        "Mean_AF": mean_af,
        "gnomAD_NFE_AF": gnomad_nfe,
        "Common_SNP_By_NFE": common_snp_by_nfe,
        "Detected_rsIDs": detected_rsids,
        "detailed_consequence": detailed_consequence,
        "functional_class": map_functional_class(worst_consequence_term),
        "impact_severity": impact_severity,
        "representative_transcript": clean_text(representative.get("Feature", "")),
        "annotation_basis": describe_annotation_basis(representative),
        "sift_raw": sift_raw,
        "sift_prediction": sift_prediction,
        "polyphen_raw": polyphen_raw,
        "polyphen_prediction": polyphen_prediction,
    }

def deduplicate_detected_variants(detected: pd.DataFrame) -> pd.DataFrame:
    """Collapse the annotated workbook to one row per unique biological variant."""
    variant_rows = [
        collapse_variant_annotations(group)
        for _, group in detected.groupby("Variant_Key", sort=False, dropna=False)
    ]
    variant_df = pd.DataFrame(variant_rows)
    if variant_df.empty:
        raise ValueError("No unique variants were available after collapsing transcript-expanded rows.")
    variant_df = variant_df.sort_values(["CHROM", "POS", "REF", "ALT"], key=lambda col: col.map(chrom_sort_key) if col.name == "CHROM" else col)
    return variant_df.reset_index(drop=True)


def technical_interpretation(audit_row: dict[str, Any]) -> tuple[float | None, str]:
    """Summarise amplicon technical context in plain language."""
    fail_values = [
        audit_row.get("Global_Failure_%"),
        audit_row.get("Breast_Control_fail_%"),
        audit_row.get("Breast_Tumour_fail_%"),
        audit_row.get("Endometrium_Control_fail_%"),
        audit_row.get("Endometrium_Tumour_fail_%"),
    ]
    numeric_fail_values = [float(value) for value in fail_values if pd.notna(value)]
    max_failure_pct = max(numeric_fail_values) if numeric_fail_values else None
    mean_depth = audit_row.get("Global_Mean_Depth")
    worst_flag = clean_text(audit_row.get("Selected_Worst_Amplicon", "NO")) == "YES"

    if worst_flag:
        return max_failure_pct, "Caution: selected worst-performing amplicon"
    if max_failure_pct is not None and max_failure_pct >= 5:
        return max_failure_pct, "Caution: elevated failure rate"
    if (max_failure_pct is not None and max_failure_pct >= 2) or (pd.notna(mean_depth) and float(mean_depth) < 1500):
        return max_failure_pct, "Watch: moderate technical concern"
    if max_failure_pct is None and pd.isna(mean_depth):
        return max_failure_pct, "No technical audit data"
    return max_failure_pct, "Technically acceptable"


def build_interval_lookup(bed_df: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], dict[str, set[int]]]:
    """Prepare simple chromosome and rsID lookups for interval-first matching."""
    chrom_lookup = {
        chrom: sub_df.copy()
        for chrom, sub_df in bed_df.groupby("chr", sort=False)
    }
    rsid_lookup: dict[str, set[int]] = defaultdict(set)
    for row in bed_df.itertuples(index=False):
        if clean_text(row.rsid):
            rsid_lookup[row.rsid].add(int(row.BED_Target_ID))
    return chrom_lookup, rsid_lookup


def variant_overlap_bounds(row: pd.Series) -> tuple[int, int]:
    """Calculate a simple 1-based variant span for overlap checks."""
    ref = clean_text(row["REF"])
    ref_span = max(1, len(ref.replace("-", "")))
    start = int(row["POS"])
    end = start + ref_span - 1
    return start, end


def format_variant_example(row: pd.Series) -> str:
    """Build a compact label used in the BED-target summary sheet."""
    rsid = clean_text(row.get("rsID_detected", ""))
    if rsid:
        return rsid
    return row["Variant_Key"]



def match_variants_to_bed(variant_df: pd.DataFrame, bed_df: pd.DataFrame) -> tuple[pd.DataFrame, dict[int, list[str]]]:
    """Link variants to BED targets while separating explicit design hits from extras."""
    chrom_lookup, rsid_lookup = build_interval_lookup(bed_df)
    bed_index = bed_df.set_index("BED_Target_ID", drop=False)
    matched_rows: list[dict[str, Any]] = []
    bed_to_variants: dict[int, list[str]] = defaultdict(list)

    for row in variant_df.to_dict(orient="records"):
        variant = pd.Series(row)
        interval_ids: set[int] = set()
        rsid_ids: set[int] = set()
        explicit_non_rsid_ids: set[int] = set()
        chrom = variant["CHROM"]
        var_start, var_end = variant_overlap_bounds(variant)

        chrom_targets = chrom_lookup.get(chrom)
        if chrom_targets is not None:
            overlaps = chrom_targets[
                (chrom_targets["bed_start_1based"] <= var_end)
                & (chrom_targets["bed_end_1based"] >= var_start)
            ]
            interval_ids = set(overlaps["BED_Target_ID"].tolist())
            explicit_non_rsid_ids = set(overlaps[overlaps["rsid"].fillna("").astype(str).str.strip() == ""]["BED_Target_ID"].tolist())

        for rsid in variant["Detected_rsIDs"]:
            rsid_ids.update(rsid_lookup.get(rsid, set()))

        matched_ids = sorted(interval_ids | rsid_ids)
        targeted_region_overlap = "YES" if interval_ids else "NO"
        explicit_bed_annotation = "YES" if (rsid_ids or explicit_non_rsid_ids) else "NO"
        in_bed_design = explicit_bed_annotation
        extra_finding = "NO" if explicit_bed_annotation == "YES" else "YES"

        if interval_ids and rsid_ids:
            match_type = "both"
        elif interval_ids:
            match_type = "interval"
        elif rsid_ids:
            match_type = "rsID"
        else:
            match_type = "none"

        if explicit_bed_annotation == "YES":
            match_status = "Explicitly represented in BED"
        elif targeted_region_overlap == "YES":
            match_status = "Within targeted BED interval only"
        else:
            match_status = "Outside BED design"

        if matched_ids:
            matched_bed = bed_index.loc[matched_ids]
            bed_annotation_label = join_unique(matched_bed["BED_annotation_label"])
            bed_amplicon_id = join_unique(matched_bed["amplicon_id"])
            bed_rsid = join_unique(matched_bed["rsid"])
            bed_gene = join_unique(matched_bed["gene"])
            for bed_id in matched_ids:
                bed_to_variants[int(bed_id)].append(variant["Variant_Key"])
        else:
            bed_annotation_label = ""
            bed_amplicon_id = ""
            bed_rsid = ""
            bed_gene = ""

        matched_rows.append({
            "CHROM": variant["CHROM"],
            "POS": int(variant["POS"]),
            "REF": variant["REF"],
            "ALT": variant["ALT"],
            "Variant_Key": variant["Variant_Key"],
            "rsID_detected": variant["rsID_detected"],
            "SYMBOL": variant["SYMBOL"],
            "Consequence": variant["Consequence"],
            "IMPACT": variant["IMPACT"],
            "HGVSp": variant["HGVSp"],
            "BED_match_status": match_status,
            "BED_match_type": match_type,
            "BED_display_label": bed_annotation_label,
            "BED_amplicon_id": bed_amplicon_id,
            "BED_gene": bed_gene,
            "gnomAD_NFE_AF": variant["gnomAD_NFE_AF"],
            "Common_SNP_By_NFE": "YES" if variant["Common_SNP_By_NFE"] else "NO",
            "Targeted_Region_Overlap": targeted_region_overlap,
            "Explicit_BED_Annotation": explicit_bed_annotation,
            "In_BED_Design": in_bed_design,
            "Extra_Finding": extra_finding,
            "N_samples": variant["N_samples"],
            "Sample_examples": variant["Sample_examples"],
            "Cohorts": variant["Cohorts"],
            "Tissues": variant["Tissues"],
            "Mean_DP": variant["Mean_DP"],
            "Mean_AF": variant["Mean_AF"],
            "detailed_consequence": variant["detailed_consequence"],
            "functional_class": variant["functional_class"],
            "impact_severity": variant["impact_severity"],
            "representative_transcript": variant["representative_transcript"],
            "annotation_basis": variant["annotation_basis"],
            "sift_raw": variant["sift_raw"],
            "sift_prediction": variant["sift_prediction"],
            "polyphen_raw": variant["polyphen_raw"],
            "polyphen_prediction": variant["polyphen_prediction"],
            "Variant_Example_Label": format_variant_example(variant),
        })

    matched_df = pd.DataFrame(matched_rows)
    matched_df = matched_df.sort_values(["CHROM", "POS", "REF", "ALT"], key=lambda col: col.map(chrom_sort_key) if col.name == "CHROM" else col)
    return matched_df.reset_index(drop=True), bed_to_variants


def build_bed_targets_summary(
    bed_df: pd.DataFrame,
    matched_variants: pd.DataFrame,
    bed_to_variants: dict[int, list[str]],
    technical_audit_df: pd.DataFrame,
) -> pd.DataFrame:
    """Build one summary row per BED target."""
    example_lookup = matched_variants.set_index("Variant_Key")["Variant_Example_Label"].to_dict()
    summary_rows: list[dict[str, Any]] = []
    audit_lookup = {}
    if not technical_audit_df.empty:
        for row in technical_audit_df.to_dict(orient="records"):
            key = (row["chr"], int(row["start"]), int(row["end"]), clean_text(row["amplicon_id"]))
            audit_lookup[key] = row

    for row in bed_df.to_dict(orient="records"):
        bed_id = int(row["BED_Target_ID"])
        variant_keys = sorted(set(bed_to_variants.get(bed_id, [])))
        examples = [example_lookup.get(key, key) for key in variant_keys]
        audit_row = audit_lookup.get((row["chr"], int(row["start"]), int(row["end"]), clean_text(row["amplicon_id"])), {})
        max_failure_pct, technical_note = technical_interpretation(audit_row)
        summary_rows.append({
            "chr": row["chr"],
            "start": row["start"],
            "end": row["end"],
            "amplicon_id": row["amplicon_id"],
            "rsid": row["rsid"],
            "gene": row["gene"],
            "BED_display_label": row["BED_annotation_label"],
            "any_detected_variant": "YES" if variant_keys else "NO",
            "n_detected_variants": len(variant_keys),
            "detected_variant_examples": summarise_examples(examples, max_items=5),
            "Selected_Worst_Amplicon": audit_row.get("Selected_Worst_Amplicon", "NO"),
            "Max_Group_Failure_%": max_failure_pct if max_failure_pct is not None else pd.NA,
            "Technical_Interpretation": technical_note,
            "Global_Failure_%": audit_row.get("Global_Failure_%", pd.NA),
            "Global_Mean_Depth": audit_row.get("Global_Mean_Depth", pd.NA),
            "Global_Median_Depth": audit_row.get("Global_Median_Depth", pd.NA),
            "Breast_Control_fail_%": audit_row.get("Breast_Control_fail_%", pd.NA),
            "Breast_Tumour_fail_%": audit_row.get("Breast_Tumour_fail_%", pd.NA),
            "Endometrium_Control_fail_%": audit_row.get("Endometrium_Control_fail_%", pd.NA),
            "Endometrium_Tumour_fail_%": audit_row.get("Endometrium_Tumour_fail_%", pd.NA),
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.sort_values(["chr", "start", "end"], key=lambda col: col.map(chrom_sort_key) if col.name == "chr" else col)
    return summary_df.reset_index(drop=True)


def build_summary_sheet(bed_summary: pd.DataFrame, matched_variants: pd.DataFrame) -> pd.DataFrame:
    """Create a compact counts-and-percentages summary sheet."""
    total_bed_targets = len(bed_summary)
    bed_targets_with_variant = int((bed_summary["any_detected_variant"] == "YES").sum())
    bed_targets_without_variant = total_bed_targets - bed_targets_with_variant
    total_detected_variants = len(matched_variants)
    targeted_total = int((matched_variants["Targeted_Region_Overlap"] == "YES").sum())
    targeted_region_only = int(((matched_variants["Targeted_Region_Overlap"] == "YES") & (matched_variants["Explicit_BED_Annotation"] == "NO")).sum())
    outside_bed = int((matched_variants["Targeted_Region_Overlap"] == "NO").sum())
    explicit_bed = int((matched_variants["Explicit_BED_Annotation"] == "YES").sum())
    detected_extra = int((matched_variants["Extra_Finding"] == "YES").sum())
    targeted_common_snps = int(((matched_variants["Targeted_Region_Overlap"] == "YES") & (matched_variants["Common_SNP_By_NFE"] == "YES")).sum())
    targeted_other_variants = int(((matched_variants["Targeted_Region_Overlap"] == "YES") & (matched_variants["Common_SNP_By_NFE"] == "NO")).sum())
    extra_targeted_common_snps = int(((matched_variants["Targeted_Region_Overlap"] == "YES") & (matched_variants["Extra_Finding"] == "YES") & (matched_variants["Common_SNP_By_NFE"] == "YES")).sum())
    extra_targeted_other_variants = int(((matched_variants["Targeted_Region_Overlap"] == "YES") & (matched_variants["Extra_Finding"] == "YES") & (matched_variants["Common_SNP_By_NFE"] == "NO")).sum())
    sift_available = int(matched_variants["sift_prediction"].astype(str).str.strip().ne("").sum())
    polyphen_available = int(matched_variants["polyphen_prediction"].astype(str).str.strip().ne("").sum())

    def pct(numerator: int, denominator: int) -> float:
        return (numerator / denominator * 100.0) if denominator else 0.0

    summary_rows = [
        ("total BED targets", total_bed_targets),
        ("BED targets with at least one detected variant", bed_targets_with_variant),
        ("BED targets with no detected variant", bed_targets_without_variant),
        ("BED target detection percentage", pct(bed_targets_with_variant, total_bed_targets)),
        ("total unique detected variants", total_detected_variants),
        ("unique designed-panel variants found", explicit_bed),
        ("unique variants found within designed amplicons", targeted_total),
        (f"unique common SNPs within designed amplicons (gnomAD NFE >= {MIN_NFE_AF:.2%})", targeted_common_snps),
        (f"unique other variants within designed amplicons (gnomAD NFE < {MIN_NFE_AF:.2%} or missing)", targeted_other_variants),
        (f"unique extra common SNPs within designed amplicons (gnomAD NFE >= {MIN_NFE_AF:.2%})", extra_targeted_common_snps),
        (f"unique extra other variants within designed amplicons (gnomAD NFE < {MIN_NFE_AF:.2%} or missing)", extra_targeted_other_variants),
        ("detected variants within targeted BED intervals only", targeted_region_only),
        ("detected variants outside BED intervals", outside_bed),
        ("detected variants classified as extra", detected_extra),
        ("percentage explicitly represented in BED", pct(explicit_bed, total_detected_variants)),
        ("percentage extra", pct(detected_extra, total_detected_variants)),
        ("variants with SIFT prediction", sift_available),
        ("variants with PolyPhen prediction", polyphen_available),
        ("annotation source", ANNOTATION_SOURCE_NOTE),
        ("representative transcript rule", "PICK if available, then CANONICAL, then MANE_SELECT, then worst IMPACT/consequence."),
        ("functional_class mapping note", "Uses representative transcript worst consequence term: synonymous, missense, nonsense, frameshift, splice_site, intronic, regulatory, intergenic, utr, other."),
        ("impact_severity note", "Uses VEP IMPACT when present; otherwise inferred from consequence with explicit rule table in script 08b."),
        ("deleteriousness note", "SIFT/PolyPhen are transcript-level VEP fields. Values are blank where predictors are not available or not applicable."),
        ("counting note", "All counts are unique variant-level counts using CHROM+POS+REF+ALT."),
        ("SNP definition note", f"Common SNP = gnomAD_NFE_AF_combined >= {MIN_NFE_AF:.2%}."),
        ("technical interpretation note 1", "Technical interpretation uses the amplicon audit to flag weak regions."),
        ("technical interpretation note 2", "Caution = selected worst-performing amplicon or clearly elevated failure."),
        ("technical interpretation note 3", "Watch = moderate concern based on failure percentage or lower mean depth."),
        ("column guide 1", "BED_display_label = best human-readable BED label for that target or match."),
        ("column guide 2", "In_BED_Design = variant explicitly represented by the BED annotation itself."),
        ("column guide 3", "Extra_Finding = variant detected but not explicitly represented by the BED annotation."),
        ("column guide 4", "Targeted_Region_Overlap = variant falls inside a targeted BED interval, even if extra."),
        ("column guide 5", "Common_SNP_By_NFE = gnomAD NFE frequency at or above the common-SNP threshold."),
        ("column guide 6", "Technical_Interpretation = quick amplicon-quality context from the technical audit."),
        ("column guide 7", "detailed_consequence = full VEP consequence string from the representative transcript."),
        ("column guide 8", "functional_class = reproducible simplified class derived from the worst representative consequence term."),
        ("column guide 9", "impact_severity = representative VEP IMPACT or inferred severity when IMPACT is blank."),
        ("column guide 10", "representative_transcript = transcript feature ID chosen by the deterministic prioritisation rule."),
        ("column guide 11", "annotation_basis = whether the representative transcript was selected via PICK, canonical, MANE, or worst available transcript."),
        ("column guide 12", "sift_prediction / polyphen_prediction = parsed qualitative labels for quick review."),
        ("interpretation note 1", "Extra can still mean the variant falls inside a targeted amplicon."),
        ("interpretation note 2", "Here, extra means not explicitly represented by the BED annotation itself."),
    ]
    return pd.DataFrame(summary_rows, columns=["Metric", "Value"])

def autosize_worksheet(ws) -> None:
    """Set readable column widths without producing excessively wide sheets."""
    for column_cells in ws.iter_cols(min_row=1, max_row=ws.max_row, min_col=1, max_col=ws.max_column):
        max_length = 0
        column_letter = get_column_letter(column_cells[0].column)
        for cell in column_cells:
            if cell.value is None:
                continue
            max_length = max(max_length, len(str(cell.value)))
        ws.column_dimensions[column_letter].width = min(max(max_length + 2, 12), 45)


def apply_main_sheet_highlighting(ws) -> None:
    """Highlight explicit BED hits and extra findings for fast supervisor review."""
    header_lookup = {clean_text(cell.value): cell.column for cell in ws[1]}
    in_bed_col = header_lookup.get("In_BED_Design")
    extra_col = header_lookup.get("Extra_Finding")
    targeted_col = header_lookup.get("Targeted_Region_Overlap")
    explicit_col = header_lookup.get("Explicit_BED_Annotation")
    status_col = header_lookup.get("BED_match_status")
    type_col = header_lookup.get("BED_match_type")

    for row_idx in range(2, ws.max_row + 1):
        in_bed_value = clean_text(ws.cell(row=row_idx, column=in_bed_col).value) if in_bed_col else ""
        extra_value = clean_text(ws.cell(row=row_idx, column=extra_col).value) if extra_col else ""
        targeted_value = clean_text(ws.cell(row=row_idx, column=targeted_col).value) if targeted_col else ""
        explicit_value = clean_text(ws.cell(row=row_idx, column=explicit_col).value) if explicit_col else ""

        if extra_value == "YES":
            for col_idx in [extra_col, targeted_col, status_col, type_col]:
                if col_idx:
                    ws.cell(row=row_idx, column=col_idx).fill = ORANGE_FILL
        elif in_bed_value == "YES" or explicit_value == "YES":
            for col_idx in [in_bed_col, explicit_col, status_col, type_col]:
                if col_idx:
                    ws.cell(row=row_idx, column=col_idx).fill = GREEN_FILL
        elif targeted_value == "YES":
            for col_idx in [targeted_col, status_col, type_col]:
                if col_idx:
                    ws.cell(row=row_idx, column=col_idx).fill = GREY_FILL
        else:
            for col_idx in [status_col, type_col]:
                if col_idx:
                    ws.cell(row=row_idx, column=col_idx).fill = GREY_FILL

def apply_bed_summary_highlighting(ws) -> None:
    """Highlight BED targets with and without detected variants and flag weak amplicons."""
    header_lookup = {clean_text(cell.value): cell.column for cell in ws[1]}
    detected_col = header_lookup.get("any_detected_variant")
    worst_col = header_lookup.get("Selected_Worst_Amplicon")
    if detected_col:
        for row_idx in range(2, ws.max_row + 1):
            cell = ws.cell(row=row_idx, column=detected_col)
            value = clean_text(cell.value)
            if value == "YES":
                cell.fill = GREEN_FILL
            elif value == "NO":
                cell.fill = GREY_FILL
    if worst_col:
        for row_idx in range(2, ws.max_row + 1):
            cell = ws.cell(row=row_idx, column=worst_col)
            value = clean_text(cell.value)
            if value == "YES":
                cell.fill = ORANGE_FILL
            elif value == "NO":
                cell.fill = GREEN_FILL
    interpretation_col = header_lookup.get("Technical_Interpretation")
    if interpretation_col:
        for row_idx in range(2, ws.max_row + 1):
            cell = ws.cell(row=row_idx, column=interpretation_col)
            value = clean_text(cell.value)
            if value.startswith("Caution"):
                cell.fill = ORANGE_FILL
            elif value.startswith("Watch"):
                cell.fill = GREY_FILL
            elif value == "Technically acceptable":
                cell.fill = GREEN_FILL


def apply_summary_formatting(ws) -> None:
    """Apply a simple formatted layout to the Summary sheet."""
    for row_idx in range(2, ws.max_row + 1):
        ws.cell(row=row_idx, column=1).fill = SUMMARY_LABEL_FILL
        ws.cell(row=row_idx, column=1).font = Font(bold=True)
    header_lookup = {clean_text(cell.value): cell.column for cell in ws[1]}
    value_col = header_lookup.get("Value")
    metric_lookup = {clean_text(ws.cell(row=row_idx, column=1).value): row_idx for row_idx in range(2, ws.max_row + 1)}
    if value_col:
        for metric_name in [
            "BED target detection percentage",
            "percentage explicitly represented in BED",
            "percentage extra",
        ]:
            row_idx = metric_lookup.get(metric_name)
            if row_idx:
                ws.cell(row=row_idx, column=value_col).number_format = "0.0"


def format_workbook(output_path: Path) -> None:
    """Apply clean, thesis-friendly formatting to all workbook sheets."""
    workbook = load_workbook(output_path)
    for ws in workbook.worksheets:
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
        for cell in ws[1]:
            cell.fill = HEADER_FILL
            cell.font = HEADER_FONT
            cell.alignment = Alignment(horizontal="center", vertical="center")
        autosize_worksheet(ws)

        if ws.title in {"Detected_vs_BED", "Extra_detected_only"}:
            apply_main_sheet_highlighting(ws)
            header_lookup = {clean_text(cell.value): cell.column for cell in ws[1]}
            for col_name in ["Mean_DP", "Mean_AF"]:
                col_idx = header_lookup.get(col_name)
                if col_idx:
                    for row_idx in range(2, ws.max_row + 1):
                        ws.cell(row=row_idx, column=col_idx).number_format = "0.0000"
        elif ws.title == "BED_targets_summary":
            apply_bed_summary_highlighting(ws)
            header_lookup = {clean_text(cell.value): cell.column for cell in ws[1]}
            for col_name in [
                "Max_Group_Failure_%", "Global_Failure_%", "Breast_Control_fail_%", "Breast_Tumour_fail_%",
                "Endometrium_Control_fail_%", "Endometrium_Tumour_fail_%"
            ]:
                col_idx = header_lookup.get(col_name)
                if col_idx:
                    for row_idx in range(2, ws.max_row + 1):
                        ws.cell(row=row_idx, column=col_idx).number_format = "0.00"
            for col_name in ["Global_Mean_Depth", "Global_Median_Depth"]:
                col_idx = header_lookup.get(col_name)
                if col_idx:
                    for row_idx in range(2, ws.max_row + 1):
                        ws.cell(row=row_idx, column=col_idx).number_format = "0.0"
        elif ws.title == "Summary":
            apply_summary_formatting(ws)

    workbook.save(output_path)


def export_outputs(
    matched_variants: pd.DataFrame,
    bed_summary: pd.DataFrame,
    summary_df: pd.DataFrame,
    output_path: Path,
    csv_output: Path | None,
) -> None:
    """Write the Excel workbook and optional CSV export."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    extra_only = matched_variants[matched_variants["Extra_Finding"] == "YES"].copy()

    for frame in [matched_variants, bed_summary, extra_only, summary_df]:
        if "Variant_Example_Label" in frame.columns:
            frame.drop(columns=["Variant_Example_Label"], inplace=True)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        matched_variants.to_excel(writer, sheet_name="Detected_vs_BED", index=False)
        bed_summary.to_excel(writer, sheet_name="BED_targets_summary", index=False)
        extra_only.to_excel(writer, sheet_name="Extra_detected_only", index=False)
        summary_df.to_excel(writer, sheet_name="Summary", index=False)

    format_workbook(output_path)
    if csv_output is not None:
        csv_output.parent.mkdir(parents=True, exist_ok=True)
        matched_variants.to_csv(csv_output, index=False)


def run_report(
    bed_path: Path,
    annotated_xlsx: Path,
    sheet_name: str,
    technical_audit_csv: Path | None,
    output_path: Path,
    csv_output: Path | None,
) -> None:
    """Orchestrate the BED-versus-detected comparison workflow."""
    print(">>> Script 08b: Comparing designed BED targets against detected variants...")
    validate_file_exists(bed_path, "Script 08b BED input")
    validate_file_exists(annotated_xlsx, "Script 08b annotated workbook input")

    bed_df = load_bed_file(bed_path)
    technical_audit_df = load_technical_audit(technical_audit_csv)
    detected_df = load_detected_annotations(annotated_xlsx, sheet_name)
    variant_df = deduplicate_detected_variants(detected_df)
    matched_variants, bed_to_variants = match_variants_to_bed(variant_df, bed_df)
    bed_summary = build_bed_targets_summary(bed_df, matched_variants, bed_to_variants, technical_audit_df)
    summary_df = build_summary_sheet(bed_summary, matched_variants)
    export_outputs(matched_variants, bed_summary, summary_df, output_path, csv_output)

    designed_count = int((matched_variants["In_BED_Design"] == "YES").sum())
    extra_count = int((matched_variants["Extra_Finding"] == "YES").sum())
    print(f"    BED targets parsed: {len(bed_df)}")
    print(f"    Unique detected variants: {len(matched_variants)}")
    print(f"    Variants linked to BED: {designed_count}")
    print(f"    Extra detected variants: {extra_count}")
    print(f"    Excel workbook saved to: {output_path}")
    if csv_output is not None:
        print(f"    CSV export saved to: {csv_output}")


def main() -> None:
    """Parse CLI arguments and run the comparison."""
    args = build_argument_parser().parse_args()
    bed_path = Path(args.bed).expanduser().resolve()
    annotated_xlsx = Path(args.annotated_xlsx).expanduser().resolve()
    technical_audit_csv = Path(args.technical_audit_csv).expanduser().resolve() if clean_text(args.technical_audit_csv) else None
    output_path = Path(args.output).expanduser().resolve()
    csv_output = Path(args.csv_output).expanduser().resolve() if clean_text(args.csv_output) else None
    run_report(bed_path, annotated_xlsx, args.sheet, technical_audit_csv, output_path, csv_output)


if __name__ == "__main__":
    main()
