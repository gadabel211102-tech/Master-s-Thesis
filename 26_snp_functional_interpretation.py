#!/usr/bin/env python3
"""Stage 26: SNP functional interpretation and biological-context integration."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from association_runtime import script26_defaults
from pipeline_utils import (
    attach_amplicon_warning_columns,
    build_amplicon_warning_lookup,
    combine_gnomad_nfe,
    ensure_directory,
    extract_rsid,
)


LOCUS_GENES = {"ERBB2", "GSDMB", "IKZF3", "ORMDL3", "GSDMA", "PGAP3", "PSMD3", "ZPBP2", "LRRC3C"}
IMPACT_ORDER = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
HER2_LABELS = {"HER2 subtype", "HER2 copies"}
OUTPUT_WORKBOOK = "26_SNP_Functional_Interpretation.xlsx"
OUTPUT_TABLE_TSV = "26_SNP_Interpretation_Table.tsv"
OUTPUT_SUMMARY_TXT = "26_Thesis_Interpretation_Summary.txt"
OUTPUT_INTRONIC_ALAMUT_XLSX = "26_Intronic_Variant_Alamut_Shortlist.xlsx"
OUTPUT_INTRONIC_ALAMUT_TSV = "26_Intronic_Variant_Alamut_Shortlist.tsv"
OUTPUT_SUPERVISOR_SUMMARY_XLSX = "26_Supervisor_Shareable_Summary.xlsx"


DEFAULTS = script26_defaults()
OUT_DIR = ensure_directory(DEFAULTS["out_dir"])
FDR_THRESHOLD = float(DEFAULTS["fdr_threshold"])
RESULTS_ROOT = OUT_DIR.parent
PATHS = {
    key: Path(value) if not isinstance(value, Path) else value
    for key, value in DEFAULTS.items()
    if key not in {"fdr_threshold", "grouping"}
}
PATHS.update({
    "stage17b": RESULTS_ROOT / "17b_rare_variant_associations" / "GSDMB_Rare_Variant_Association_Results.xlsx",
    "stage23": RESULTS_ROOT / "23_rna_integration" / "23_RNA_Genetic_Integration.xlsx",
    "stage25": RESULTS_ROOT / "25_rs11078928_rs869402_haplotype_focus" / "25_Supervisor_Summary.tsv",
})


def s(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip()


def safe_float(value: object) -> float:
    try:
        return float(value)
    except Exception:
        return math.nan


def unique_keep_order(values) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = s(value)
        if not text or text.lower() in {"nan", "none", ".", "na", "n/a", "<na>"}:
            continue
        if text not in seen:
            seen.add(text)
            out.append(text)
    return out


def join_unique(values, sep: str = "; ", limit: int | None = None) -> str:
    items = unique_keep_order(values)
    if limit is not None and len(items) > limit:
        return sep.join(items[:limit]) + f"{sep}..."
    return sep.join(items)


def _first_nonempty_text(series) -> str:
    for value in series:
        text = s(value)
        if text:
            return text
    return ""


def _min_numeric_across(row: pd.Series, columns: list[str]) -> float:
    values = [pd.to_numeric(row.get(col), errors="coerce") for col in columns if col in row.index]
    values = [float(value) for value in values if pd.notna(value)]
    return min(values) if values else math.nan


def split_pubmed_ids(value: object) -> list[str]:
    return re.findall(r"\d{5,}", s(value))


def variant_id_from_row(row: pd.Series) -> str:
    existing = s(row.get("Existing_variation"))
    match = re.search(r"(rs\d+)", existing, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    chrom = s(row.get("CHROM"))
    pos = pd.to_numeric(row.get("POS"), errors="coerce")
    ref = s(row.get("REF")).upper()
    alt = s(row.get("ALT")).upper()
    if chrom and pd.notna(pos) and ref and alt:
        return f"{chrom}:{int(pos)}_{ref}>{alt}"
    if chrom and pd.notna(pos):
        return f"{chrom}:{int(pos)}"
    return existing or "UNKNOWN_VARIANT"


def truthy(value: object) -> bool:
    text = s(value).upper()
    return text in {"TRUE", "YES", "Y", "1"}


def fmt_num(value: object, digits: int = 3) -> str:
    num = safe_float(value)
    if math.isnan(num):
        return "NA"
    return f"{num:.{digits}g}"


def genes_from_text(text: str) -> set[str]:
    return {item.strip().upper() for item in re.split(r"[;,]", text) if item.strip()}


def _split_gene_tokens(value: object) -> list[str]:
    return [item.strip() for item in re.split(r"[;,]", s(value)) if item.strip()]


def _looks_like_ensembl_gene_text(value: object) -> bool:
    tokens = _split_gene_tokens(value)
    return bool(tokens) and all(re.fullmatch(r"ENSG\d+(?:\.\d+)?", token, flags=re.IGNORECASE) for token in tokens)


def _display_gene_name(primary_gene: object, annotation_gene: object) -> str:
    primary = s(primary_gene)
    annotation = s(annotation_gene)
    if primary and not _looks_like_ensembl_gene_text(primary):
        return primary
    if annotation and not _looks_like_ensembl_gene_text(annotation):
        return annotation
    return primary or annotation


def representative_annotation_row(sub: pd.DataFrame) -> pd.Series:
    ranked = sub.copy()
    ranked["_canonical"] = ranked.get("CANONICAL", pd.Series(index=ranked.index, dtype=object)).astype(str).str.upper().eq("YES").astype(int)
    ranked["_mane"] = (
        ranked.get("MANE_SELECT", pd.Series(index=ranked.index, dtype=object)).astype(str).str.upper().ne("")
        & ranked.get("MANE_SELECT", pd.Series(index=ranked.index, dtype=object)).astype(str).str.upper().ne("NAN")
        & ranked.get("MANE_SELECT", pd.Series(index=ranked.index, dtype=object)).astype(str).str.upper().ne(".")
    ).astype(int)
    ranked["_impact_rank"] = ranked.get("IMPACT", pd.Series(index=ranked.index, dtype=object)).astype(str).str.upper().map(IMPACT_ORDER).fillna(99)
    ranked["_cadd"] = pd.to_numeric(
        ranked.get("CADD_PHRED", ranked.get("CADD_phred", pd.Series(index=ranked.index, dtype=float))),
        errors="coerce",
    ).fillna(-1)
    ranked = ranked.sort_values(["_canonical", "_mane", "_impact_rank", "_cadd"], ascending=[False, False, True, False], kind="stable")
    return ranked.iloc[0]


def load_excel_sheet(path: Path, sheet_name: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_excel(path, sheet_name=sheet_name)


def load_curated_evidence(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["SNP_ID", "Evidence_Kind", "Evidence_Scope", "Target_Gene", "Effect_Direction", "Evidence_Strength", "Statement", "Source_Label", "Citation", "PMID", "URL"])
    df = pd.read_csv(path, sep="\t")
    df["SNP_ID"] = df["SNP_ID"].astype(str).str.strip()
    return df


def load_annotation_backbone(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Biological_Annotations")
    df["SNP_ID"] = df.apply(variant_id_from_row, axis=1)
    df = combine_gnomad_nfe(df)
    warning_lookup = build_amplicon_warning_lookup(df, variant_col="SNP_ID", chrom_col="CHROM", pos_col="POS")
    df = attach_amplicon_warning_columns(df, warning_lookup, variant_col="SNP_ID")

    rows: list[dict[str, object]] = []
    for snp_id, sub in df.groupby("SNP_ID", sort=False):
        rep = representative_annotation_row(sub)
        gnomad = pd.to_numeric(sub["gnomAD_NFE_AF_combined"], errors="coerce")
        max_af = pd.to_numeric(sub.get("MAX_AF", pd.Series(index=sub.index, dtype=float)), errors="coerce")
        pubmed_ids = unique_keep_order([pmid for value in sub.get("PUBMED", pd.Series(index=sub.index, dtype=object)) for pmid in split_pubmed_ids(value)])
        gene_text = join_unique(sub.get("SYMBOL", pd.Series(index=sub.index, dtype=object))) or join_unique(sub.get("Gene", pd.Series(index=sub.index, dtype=object)))
        cons_text = join_unique(sub.get("Consequence", pd.Series(index=sub.index, dtype=object)))
        motif_names = join_unique(sub.get("MOTIF_NAME", pd.Series(index=sub.index, dtype=object)))
        tfs = join_unique(sub.get("TRANSCRIPTION_FACTORS", pd.Series(index=sub.index, dtype=object)))
        cadd = pd.to_numeric(sub.get("CADD_PHRED", sub.get("CADD_phred", pd.Series(index=sub.index, dtype=float))), errors="coerce")
        cadd_best = cadd.max(skipna=True) if cadd.notna().any() else math.nan
        sift = join_unique(sub.get("SIFT", pd.Series(index=sub.index, dtype=object))) or join_unique(sub.get("SIFT_pred", pd.Series(index=sub.index, dtype=object)))
        polyphen = join_unique(sub.get("PolyPhen", pd.Series(index=sub.index, dtype=object))) or join_unique(sub.get("Polyphen2_HDIV_pred", pd.Series(index=sub.index, dtype=object)))

        regulatory_bits = []
        if motif_names:
            regulatory_bits.append(f"motif={motif_names}")
        if tfs:
            regulatory_bits.append(f"TF={tfs}")
        if pd.notna(rep.get("MOTIF_SCORE_CHANGE")):
            regulatory_bits.append(f"motif_score_change={fmt_num(rep.get('MOTIF_SCORE_CHANGE'))}")

        protein_bits = []
        if "missense_variant" in cons_text:
            protein_bits.append("missense consequence")
        if sift:
            protein_bits.append(f"SIFT={sift}")
        if polyphen:
            protein_bits.append(f"PolyPhen={polyphen}")
        if not math.isnan(cadd_best):
            protein_bits.append(f"CADD_PHRED={fmt_num(cadd_best)}")

        splicing_bits = []
        if "splice" in cons_text.lower():
            splicing_bits.append("splice-related VEP consequence")
        if s(rep.get("HGVSc")):
            splicing_bits.append(f"HGVSc={s(rep.get('HGVSc'))}")

        pop_context = []
        if gnomad.notna().any():
            pop_context.append(f"gnomAD NFE AF={fmt_num(gnomad.max(skipna=True), 4)}")
        source = join_unique(sub.get("gnomAD_NFE_Source", pd.Series(index=sub.index, dtype=object)))
        if source:
            pop_context.append(source)
        eur_af = pd.to_numeric(sub.get("EUR_AF", pd.Series(index=sub.index, dtype=float)), errors="coerce")
        if eur_af.notna().any():
            pop_context.append(f"1000G EUR AF={fmt_num(eur_af.max(skipna=True), 4)}")
        if max_af.notna().any():
            pop_context.append(f"max population AF={fmt_num(max_af.max(skipna=True), 4)}")

        rows.append({
            "SNP_ID": snp_id,
            "rsID": extract_rsid(rep.get("Existing_variation")) or (snp_id if re.fullmatch(r"rs\d+", str(snp_id), flags=re.IGNORECASE) else ""),
            "Position": f"{s(rep.get('CHROM'))}:{int(rep.get('POS'))}" if pd.notna(rep.get("POS")) else s(rep.get("CHROM")),
            "Chromosome": s(rep.get("CHROM")),
            "POS": int(rep.get("POS")) if pd.notna(rep.get("POS")) else pd.NA,
            "REF": s(rep.get("REF")).upper(),
            "ALT": s(rep.get("ALT")).upper(),
            "HGVSc": s(rep.get("HGVSc")),
            "Gene": gene_text,
            "Consequence": s(rep.get("Consequence")) or cons_text,
            "Consequence_All": cons_text,
            "Impact": s(rep.get("IMPACT")),
            "Population_Frequency_Context": "; ".join(pop_context),
            "Published_PMIDs": "; ".join(pubmed_ids),
            "Published_DB_Notes": join_unique(sub.get("CLIN_SIG", pd.Series(index=sub.index, dtype=object))),
            "Regulatory_Annotation": "; ".join(regulatory_bits),
            "Protein_Annotation": "; ".join(protein_bits),
            "Splicing_Annotation": "; ".join(splicing_bits),
            "Coverage_Risk_Flag": bool(rep.get("Coverage_Risk_Flag", False)),
            "Coverage_Risk_Note": s(rep.get("Coverage_Risk_Note")),
        })
    return pd.DataFrame(rows)

def load_stage12(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    for sheet in ["Global", "Breast", "Endometrium"]:
        df = load_excel_sheet(path, sheet)
        if df.empty:
            continue
        df = df.copy()
        df["Comparison"] = sheet
        df["SNP_ID"] = df["SNP_ID"].astype(str).str.strip()
        frames.append(df)
    all_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    sig_df = all_df[all_df.get("FDR_Significant", pd.Series(index=all_df.index)).eq(True)].copy() if not all_df.empty else pd.DataFrame()
    return all_df, sig_df


def load_stage17(path: Path) -> pd.DataFrame:
    df = load_excel_sheet(path, "summary_significant")
    if df.empty:
        return df
    df = df.copy()
    df["Variant_ID"] = df["Variant_ID"].astype(str).str.strip()
    return df


def load_stage20(path: Path) -> pd.DataFrame:
    df = load_excel_sheet(path, "snp_isoform_assoc")
    if df.empty:
        return df
    df = df.copy()
    df["Variant_ID"] = df["Variant_ID"].astype(str).str.strip()
    return df[df.get("FDR_Sig", pd.Series(index=df.index)).eq(True)].copy()


def load_stage20b(path: Path) -> pd.DataFrame:
    df = load_excel_sheet(path, "snp_gene_assoc")
    if df.empty:
        return df
    df = df.copy()
    df["Variant_ID"] = df["Variant_ID"].astype(str).str.strip()
    return df[df.get("FDR_Sig", pd.Series(index=df.index)).eq(True)].copy()


def _row_get(row: pd.Series, *candidates: str, default: str = "") -> str:
    """Return the first non-null value found among candidate column names."""
    for col in candidates:
        val = row.get(col)
        if val is not None and not (isinstance(val, float) and math.isnan(val)):
            text = str(val).strip()
            if text and text.lower() not in {"nan", "none", "na"}:
                return text
    return default


def load_mapping_labels(path08: Path, path09: Path) -> pd.DataFrame:
    # Uses iterrows() + dict-style access so the function is robust to column
    # renames across pipeline versions (Variant_ID vs Variant_Key,
    # Consequence vs Consequence_Raw, etc.) without hard-coded attribute names.
    rows: list[dict[str, object]] = []
    stage08 = load_excel_sheet(path08, "labelled_variants")
    if not stage08.empty:
        for _, row in stage08.iterrows():
            snp_id = _row_get(row, "Variant_ID", "Variant_Key")
            cohort  = _row_get(row, "Cohort")
            tissue  = _row_get(row, "Tissue")
            cons    = _row_get(row, "Consequence", "Consequence_Raw")
            rows.append({
                "SNP_ID": snp_id,
                "Mapping_Source": "Stage08_Global_Labelled",
                "Mapping_Context": f"{cohort}/{tissue} {cons}",
            })
    for sheet in ["common_snp_points", "lower_freq_points"]:
        stage09 = load_excel_sheet(path09, sheet)
        if stage09.empty or "Labelled_On_Figure" not in stage09.columns:
            continue
        stage09 = stage09[stage09["Labelled_On_Figure"].eq(True)].copy()
        for _, row in stage09.iterrows():
            snp_id = _row_get(row, "Variant_ID", "Variant_Key")
            cohort  = _row_get(row, "Cohort")
            tissue  = _row_get(row, "Tissue")
            cons    = _row_get(row, "Consequence", "Consequence_Raw")
            rows.append({
                "SNP_ID": snp_id,
                "Mapping_Source": f"Stage09_{sheet}",
                "Mapping_Context": f"{cohort}/{tissue} {cons}",
            })
    return pd.DataFrame(rows).drop_duplicates() if rows else pd.DataFrame(columns=["SNP_ID", "Mapping_Source", "Mapping_Context"])


def load_stage22(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = load_excel_sheet(path, "Primary_Haplotypes")
    links = load_excel_sheet(path, "Haplotype_SNP_Comparison")
    if not links.empty:
        links = links.copy()
        links["SNP_ID"] = links.get("Linked_SNP_ID", pd.Series(index=links.index, dtype=object)).astype(str).str.strip()
        links = links[links["SNP_ID"].ne("") & links["SNP_ID"].ne("None")].copy()
    return primary, links


def load_stage24(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    snp = load_excel_sheet(path, "SNP_Comparison_All")
    hap = load_excel_sheet(path, "Haplotype_Comparison_All")
    if not snp.empty:
        snp = snp.copy()
        snp["SNP_ID"] = snp["SNP_ID"].astype(str).str.strip()
    return snp, hap


def build_priority_table(annotation_df: pd.DataFrame, stage12_sig: pd.DataFrame, stage17: pd.DataFrame, stage20: pd.DataFrame, stage20b: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    reasons: dict[str, list[str]] = {}

    def add(snp_id: str, reason: str) -> None:
        snp = s(snp_id)
        if not snp:
            return
        reasons.setdefault(snp, [])
        if reason not in reasons[snp]:
            reasons[snp].append(reason)

    for snp in stage12_sig.get("SNP_ID", pd.Series(dtype=object)):
        add(str(snp), "Stage12_FDR_significant")
    for snp in stage17.get("Variant_ID", pd.Series(dtype=object)):
        add(str(snp), "Stage17_significant")
    for snp in stage20.get("Variant_ID", pd.Series(dtype=object)):
        add(str(snp), "Stage20_isoform_signal")
    for snp in stage20b.get("Variant_ID", pd.Series(dtype=object)):
        add(str(snp), "Stage20b_gene_expression_signal")
    for snp in mapping.get("SNP_ID", pd.Series(dtype=object)):
        add(str(snp), "Mapping_labelled")

    available = set(annotation_df["SNP_ID"].astype(str))
    rows = []
    for snp_id, why in reasons.items():
        if snp_id not in available and not re.fullmatch(r"rs\d+", snp_id, flags=re.IGNORECASE):
            continue
        rows.append({
            "SNP_ID": snp_id,
            "Priority_Reasons": "; ".join(why),
            "Stage12_Flag": any(item.startswith("Stage12") for item in why),
            "Stage17_Flag": any(item.startswith("Stage17") for item in why),
            "Mapping_Flag": any(item.startswith("Mapping") for item in why),
            "RNA_Flag": any(item.startswith("Stage20") for item in why),
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["Priority_Rank"] = (
        out["Stage12_Flag"].astype(int) * 8
        + out["Stage17_Flag"].astype(int) * 4
        + out["RNA_Flag"].astype(int) * 2
        + out["Mapping_Flag"].astype(int)
    )
    return out.sort_values(["Priority_Rank", "SNP_ID"], ascending=[False, True], kind="stable").reset_index(drop=True)


def summarise_stage12_hits(snp_id: str, stage12_sig: pd.DataFrame) -> str:
    sub = stage12_sig[stage12_sig["SNP_ID"].eq(snp_id)].copy()
    if sub.empty:
        return "No FDR-significant stage-12 tumour-vs-control enrichment signal."
    parts = []
    sub = sub.sort_values(["Comparison", "FDR_P_Value", "P_Value"], kind="stable")
    for _, row in sub.iterrows():
        parts.append(
            f"{s(row.get('Comparison'))}: OR={fmt_num(row.get('Odds_Ratio'))}, FDR={fmt_num(row.get('FDR_P_Value'))}, tumour {fmt_num(row.get('Tumour_Freq_%'))}% vs control {fmt_num(row.get('Control_Freq_%'))}%"
        )
    return "; ".join(parts)


def summarise_stage17_hits(snp_id: str, stage17: pd.DataFrame) -> str:
    sub = stage17[stage17["Variant_ID"].eq(snp_id)].copy()
    if sub.empty:
        return "No stage-17 significant clinical or cancer-risk signal."
    parts = []
    sub = sub.sort_values(["Analysis_Type", "P_Unadj", "P_Value"], kind="stable")
    for row in sub.head(4).itertuples(index=False):
        analysis_type = s(getattr(row, "Analysis_Type", ""))
        cohort = s(getattr(row, "Cohort", ""))
        clin_label = s(getattr(row, "Clin_Label", ""))
        if analysis_type == "Tumour_vs_Control":
            parts.append(f"Tumour_vs_control {s(getattr(row, 'Analysis_Group', ''))}: FDR={fmt_num(getattr(row, 'FDR_P_Value', math.nan))}")
        elif analysis_type == "Cancer_Risk":
            parts.append(f"Cancer_risk {cohort or 'study'}: age-adjusted P={fmt_num(getattr(row, 'P_Adj_Age', math.nan))}, FDR={fmt_num(getattr(row, 'FDR_Adj_Age', math.nan))}")
        else:
            label = clin_label or analysis_type or cohort or "clinical"
            parts.append(f"{cohort}: {label} (P={fmt_num(getattr(row, 'P_Unadj', math.nan))}, FDR={fmt_num(getattr(row, 'FDR_Unadj', math.nan))})")
    return "; ".join(unique_keep_order(parts))


def summarise_rna_hits(snp_id: str, stage20: pd.DataFrame, stage20b: pd.DataFrame) -> str:
    parts = []
    sub20 = stage20[stage20["Variant_ID"].eq(snp_id)].copy()
    if not sub20.empty:
        for row in sub20.sort_values(["Cohort", "FDR", "Endpoint"], kind="stable").itertuples(index=False):
            parts.append(f"Isoform {row.Cohort}/{row.Endpoint}: slope={fmt_num(row.Slope)}, FDR={fmt_num(row.FDR)}")
    sub20b = stage20b[stage20b["Variant_ID"].eq(snp_id)].copy()
    if not sub20b.empty:
        for row in sub20b.sort_values(["Cohort", "FDR", "Gene"], kind="stable").itertuples(index=False):
            parts.append(f"Panel-gene {row.Cohort}/{row.Gene}: slope={fmt_num(row.Slope)}, FDR={fmt_num(row.FDR)}")
    return "; ".join(parts) if parts else "No FDR-significant project RNA / isoform hit was available."

def evidence_for_snp(curated: pd.DataFrame, snp_id: str) -> pd.DataFrame:
    return curated[curated["SNP_ID"].eq(snp_id)].copy() if not curated.empty else pd.DataFrame()


def her2_locus_evidence(curated: pd.DataFrame) -> pd.DataFrame:
    return curated[curated["Evidence_Kind"].eq("HER2_biology")].copy() if not curated.empty else pd.DataFrame()


def summarise_evidence_rows(rows: pd.DataFrame, kinds: set[str], fallback: str) -> str:
    if rows.empty:
        return fallback
    subset = rows[rows["Evidence_Kind"].isin(kinds)].copy()
    if subset.empty:
        return fallback
    parts = []
    for row in subset.itertuples(index=False):
        parts.append(f"{row.Evidence_Scope.replace('_', ' ')} [{str(row.Evidence_Strength).lower()}]: {row.Statement}")
    return " ".join(unique_keep_order(parts))


def summarise_other_function(annotation_row: pd.Series, snp_evidence: pd.DataFrame) -> str:
    bits = []
    if s(annotation_row.get("Protein_Annotation")):
        bits.append(s(annotation_row.get("Protein_Annotation")))
    if s(annotation_row.get("Regulatory_Annotation")):
        bits.append(s(annotation_row.get("Regulatory_Annotation")))
    if s(annotation_row.get("Published_DB_Notes")):
        bits.append(f"clinical_db={s(annotation_row.get('Published_DB_Notes'))}")
    extra = summarise_evidence_rows(
        snp_evidence,
        {"Other_functional_consequence", "General_association"},
        "",
    )
    if extra:
        bits.append(extra)
    if truthy(annotation_row.get("Coverage_Risk_Flag")):
        bits.append(s(annotation_row.get("Coverage_Risk_Note")))
    return "; ".join(unique_keep_order(bits)) or "No additional functional consequence beyond the current annotation layer was captured."


def summarise_published_relevance(annotation_row: pd.Series, snp_evidence: pd.DataFrame) -> str:
    if not snp_evidence.empty:
        parts = []
        for row in snp_evidence.head(3).itertuples(index=False):
            parts.append(f"{row.Evidence_Kind.replace('_', ' ')}: {row.Statement} ({row.Source_Label})")
        return " ".join(parts)
    pmids = s(annotation_row.get("Published_PMIDs"))
    if pmids:
        return f"VEP annotation lists PMID(s) {pmids}, but no curated SNP-specific literature summary has been added yet."
    return "No exact-SNP literature summary is currently loaded for this variant."


def summarise_her2(snp_id: str, annotation_row: pd.Series, stage17: pd.DataFrame, stage12_sig: pd.DataFrame, curated: pd.DataFrame) -> str:
    sub = stage17[(stage17["Variant_ID"].eq(snp_id)) & (stage17.get("Clin_Label", pd.Series(index=stage17.index, dtype=object)).isin(HER2_LABELS))].copy()
    if not sub.empty:
        parts = []
        for row in sub.sort_values(["P_Unadj", "FDR_Unadj"], kind="stable").itertuples(index=False):
            parts.append(f"Direct breast HER2 signal: {row.Clin_Label} (P={fmt_num(getattr(row, 'P_Unadj', math.nan))}, FDR={fmt_num(getattr(row, 'FDR_Unadj', math.nan))})")
        return "; ".join(unique_keep_order(parts))

    gene_set = genes_from_text(s(annotation_row.get("Gene")))
    breast_stage12 = stage12_sig[(stage12_sig["SNP_ID"].eq(snp_id)) & (stage12_sig["Comparison"].eq("Breast"))]
    locus_notes = her2_locus_evidence(curated)
    if not breast_stage12.empty and gene_set & LOCUS_GENES:
        extra = locus_notes.iloc[0]["Statement"] if not locus_notes.empty else "The variant sits in the ERBB2-adjacent 17q12-q21 locus."
        return f"HER2-compatible rather than HER2-specific: breast tumour association is present, but no direct HER2 subtype/copy result was significant in stage 17. {extra}"
    if gene_set & LOCUS_GENES:
        return "Not specifically HER2-linked in the project outputs, but biologically compatible with the wider ERBB2/GSDMB 17q12-q21 breast locus."
    return "Not specifically HER2-linked in the project outputs."


def summarise_haplotype_context(snp_id: str, links: pd.DataFrame, hap_external: pd.DataFrame) -> tuple[str, bool, str]:
    if links.empty:
        return "No link to a primary significant haplotype was recovered in stage 22.", False, ""
    sub = links[links["SNP_ID"].eq(snp_id)].copy()
    if sub.empty:
        return "No link to a primary significant haplotype was recovered in stage 22.", False, ""
    parts = []
    external_parts = []
    for row in sub.sort_values(["Comparison", "Region", "Haplotype_ID"], kind="stable").head(5).itertuples(index=False):
        parts.append(
            f"{row.Comparison} {row.Region}/{row.Haplotype_ID}: {row.Haplotype_Effect.lower()} haplotype; consistency={row.Consistency_Call}"
        )
        if not hap_external.empty:
            matched = hap_external[
                hap_external.get("Comparison", pd.Series(index=hap_external.index, dtype=object)).astype(str).eq(str(row.Comparison))
                & hap_external.get("Region", pd.Series(index=hap_external.index, dtype=object)).astype(str).eq(str(row.Region))
                & hap_external.get("Haplotype_ID", pd.Series(index=hap_external.index, dtype=object)).astype(str).eq(str(row.Haplotype_ID))
            ]
            if not matched.empty:
                ext = matched.iloc[0]
                external_parts.append(
                    f"{row.Comparison} {row.Region}/{row.Haplotype_ID}: study {fmt_num(ext.get('Study_Tumour_Frequency'))} vs 1000G {fmt_num(ext.get('Reference_1000G_Female_Frequency'))}; {s(ext.get('Validation_Status'))}"
                )
    return "; ".join(unique_keep_order(parts)), True, "; ".join(unique_keep_order(external_parts))


def summarise_external_support(snp_id: str, snp_external: pd.DataFrame, hap_external_summary: str) -> str:
    parts = []
    support_labels = []
    if not snp_external.empty:
        sub = snp_external[snp_external["SNP_ID"].eq(snp_id)].copy()
        for row in sub.sort_values(["Comparison", "FDR_P_Value"], kind="stable").itertuples(index=False):
            ref = safe_float(getattr(row, "Reference_1000G_Female_Carrier_Frequency", math.nan))
            study = safe_float(getattr(row, "Study_Tumour_Frequency", math.nan))
            delta = safe_float(getattr(row, "Study_vs_1000G_Delta", math.nan))
            phrase = f"{row.Comparison}: study {fmt_num(study)} vs 1000G female {fmt_num(ref)}; {s(getattr(row, 'Validation_Status', ''))}"
            parts.append(phrase)
            if not math.isnan(delta):
                support_labels.append("similar" if abs(delta) <= 0.10 else "divergent")
    if hap_external_summary:
        parts.append(f"Linked haplotypes: {hap_external_summary}")
    if not parts:
        return "No stage-24 SNP or linked-haplotype comparison against female-only 1000 Genomes was available."
    if support_labels and all(label == "similar" for label in support_labels):
        prefix = "Partial external support"
    elif "similar" in support_labels:
        prefix = "Mixed external support"
    elif support_labels:
        prefix = "External comparison suggests deviation from 1000G reference"
    else:
        prefix = "External comparison available"
    return prefix + ": " + " | ".join(parts)


def count_exact_support(snp_evidence: pd.DataFrame) -> int:
    if snp_evidence.empty:
        return 0
    return int(snp_evidence["Evidence_Scope"].astype(str).eq("Exact_SNP").sum())


def count_strong_support(snp_evidence: pd.DataFrame) -> int:
    if snp_evidence.empty:
        return 0
    return int(snp_evidence["Evidence_Strength"].astype(str).isin(["Strong", "Moderate"]).sum())


def classify_overall_interpretation(
    snp_id: str,
    annotation_row: pd.Series,
    stage12_sig: pd.DataFrame,
    stage17: pd.DataFrame,
    stage20: pd.DataFrame,
    stage20b: pd.DataFrame,
    in_key_haplotype: bool,
    snp_evidence: pd.DataFrame,
) -> tuple[str, int]:
    has_stage12 = not stage12_sig[stage12_sig["SNP_ID"].eq(snp_id)].empty
    has_stage17 = not stage17[stage17["Variant_ID"].eq(snp_id)].empty
    has_rna = not stage20[stage20["Variant_ID"].eq(snp_id)].empty or not stage20b[stage20b["Variant_ID"].eq(snp_id)].empty
    exact_support = count_exact_support(snp_evidence)
    strong_support = count_strong_support(snp_evidence)
    is_splice_or_missense = any(token in s(annotation_row.get("Consequence_All")).lower() for token in ["splice", "missense"])
    caution = truthy(annotation_row.get("Coverage_Risk_Flag"))

    score = 0
    score += 3 if has_stage12 else 0
    score += 2 if has_stage17 else 0
    score += 2 if has_rna else 0
    score += 2 if in_key_haplotype else 0
    score += strong_support
    score += 1 if is_splice_or_missense else 0
    score -= 2 if caution else 0

    if has_rna and (strong_support or is_splice_or_missense):
        label = "Supports observed project association"
        reason = "project RNA or isoform data converge with a functional mechanism"
    elif has_stage12 and strong_support:
        label = "Supports observed project association"
        reason = "project SNP enrichment sits on a literature-backed functional candidate"
    elif (has_stage12 or has_stage17 or in_key_haplotype) and (strong_support or is_splice_or_missense or exact_support):
        label = "Project association is novel but biologically plausible"
        reason = "project signal is supported mainly by annotation, locus regulation, or non-project functional evidence"
    elif has_stage12 or has_stage17 or in_key_haplotype:
        label = "Functional annotation is weak / indirect / inconsistent"
        reason = "the project signal is real, but the current functional layer remains indirect or under-powered"
    else:
        label = "Mapped or prioritised variant without a confirmed project association"
        reason = "the variant is retained for biological context, but not because it is a primary association hit"

    if caution:
        reason += "; interpret cautiously because the variant overlaps a flagged technical-risk amplicon"
    return f"{label}: {reason}.", score


def build_source_audit(annotation_df: pd.DataFrame, stage12_sig: pd.DataFrame, stage17: pd.DataFrame, stage20: pd.DataFrame, stage20b: pd.DataFrame, mapping: pd.DataFrame, stage22_links: pd.DataFrame, stage24_snp: pd.DataFrame, curated: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {"Source_Layer": "Annotated VEP backbone", "Available": True, "Rows_Used": int(annotation_df.shape[0]), "Notes": "Core SNP metadata, consequence, frequencies, PUBMED, motif, CADD, SIFT, PolyPhen."},
        {"Source_Layer": "Stage 12 SNP enrichment", "Available": not stage12_sig.empty, "Rows_Used": int(stage12_sig.shape[0]), "Notes": "Primary significant SNP set."},
        {"Source_Layer": "Stage 17 SNP clinical associations", "Available": not stage17.empty, "Rows_Used": int(stage17.shape[0]), "Notes": "Clinical significance and HER2-aware interpretation layer."},
        {"Source_Layer": "Stage 20 isoform associations", "Available": not stage20.empty, "Rows_Used": int(stage20.shape[0]), "Notes": "Project RNA isoform associations when available."},
        {"Source_Layer": "Stage 20b panel-gene expression", "Available": not stage20b.empty, "Rows_Used": int(stage20b.shape[0]), "Notes": "Project BAM-derived panel-gene associations when available."},
        {"Source_Layer": "Stage 08/09 mapping labels", "Available": not mapping.empty, "Rows_Used": int(mapping.shape[0]), "Notes": "Individually notable mapped SNPs from the landscape figures."},
        {"Source_Layer": "Stage 22 haplotype-first interpretation", "Available": not stage22_links.empty, "Rows_Used": int(stage22_links.shape[0]), "Notes": "Links SNPs to primary significant haplotypes."},
        {"Source_Layer": "Stage 24 external cross-validation", "Available": not stage24_snp.empty, "Rows_Used": int(stage24_snp.shape[0]), "Notes": "Female-only 1000 Genomes external support layer retained in the workflow."},
        {"Source_Layer": "Curated functional evidence TSV", "Available": not curated.empty, "Rows_Used": int(curated.shape[0]), "Notes": "Exact-SNP and locus-level literature notes with explicit citations and URLs."},
    ]
    return pd.DataFrame(rows)


def build_thesis_summary(final_df: pd.DataFrame, source_audit: pd.DataFrame) -> list[str]:
    prioritised = int(final_df.shape[0])
    strong = final_df.sort_values(["Biology_Score", "Significant_in_project", "SNP_ID"], ascending=[False, False, True], kind="stable").head(6)
    weak = final_df[final_df["Overall_biological_interpretation"].str.contains("weak / indirect / inconsistent", case=False, na=False)].head(6)
    lines = [
        f"Prioritised {prioritised} SNPs / polymorphisms from the significant project set plus the mapped standout variants.",
        "Functional layers used: " + "; ".join(source_audit[source_audit["Available"].eq(True)]["Source_Layer"].astype(str).tolist()) + ".",
        "HER2-related interpretation was added explicitly using the breast tumour clinical-significance layer and the ERBB2 / GSDMB 17q12-q21 locus context.",
        "External haplotype comparison against female-only 1000 Genomes remains included through the retained stage-24 layer.",
        "Most biologically meaningful SNPs:",
    ]
    for row in strong.itertuples(index=False):
        lines.append(f"- {row.SNP_ID}: {row.Overall_biological_interpretation}")
    if not weak.empty:
        lines.append("Findings that remain uncertain or weak:")
        for row in weak.itertuples(index=False):
            lines.append(f"- {row.SNP_ID}: {row.Overall_biological_interpretation}")
    return lines

def load_stage17b(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = load_excel_sheet(path, "summary_significant")
    catalogue = load_excel_sheet(path, "rare_variant_catalogue")
    if not summary.empty and "Exposure_ID" in summary.columns:
        summary = summary.copy()
        summary["Exposure_ID"] = summary["Exposure_ID"].astype(str).str.strip()
    if not catalogue.empty and "Variant_Key" in catalogue.columns:
        catalogue = catalogue.copy()
        catalogue["Variant_Key"] = catalogue["Variant_Key"].astype(str).str.strip()
    return summary, catalogue


def load_stage23_isoform_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    workbook = pd.ExcelFile(path)
    if "GSDMB Isoform Summary" in workbook.sheet_names:
        return pd.read_excel(path, sheet_name="GSDMB Isoform Summary")
    if "gsdmb_isoform_answer" in workbook.sheet_names:
        return pd.read_excel(path, sheet_name="gsdmb_isoform_answer")
    return pd.DataFrame()


def load_stage25_supervisor_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="	")


def _preferred_variant_label(rsid: str, gene: str, hgvsc: str, chrom: str, pos: object, ref: str, alt: str, fallback: str = "") -> str:
    if rsid:
        return rsid
    short_hgvsc = hgvsc.split(":")[-1] if hgvsc and ":" in hgvsc else hgvsc
    if gene and short_hgvsc:
        return f"{gene} {short_hgvsc}"
    genomic = ""
    if chrom and pd.notna(pos) and ref and alt:
        genomic = f"{chrom}:{int(pos)}:{ref}>{alt}"
    elif chrom and pd.notna(pos):
        genomic = f"{chrom}:{int(pos)}"
    if gene and genomic:
        return f"{gene} {genomic}"
    return genomic or fallback or rsid or gene or "Unknown variant"


def _lookup_annotation_row(annotation_df: pd.DataFrame, variant_id: str = "", rsid: str = "", chrom: str = "", pos: object = pd.NA, ref: str = "", alt: str = "") -> pd.Series:
    if annotation_df.empty:
        return pd.Series(dtype=object)
    if rsid:
        match = annotation_df[annotation_df.get("rsID", pd.Series(index=annotation_df.index, dtype=object)).astype(str).str.lower().eq(str(rsid).lower())]
        if not match.empty:
            return match.iloc[0]
    if variant_id:
        match = annotation_df[annotation_df["SNP_ID"].astype(str).str.lower().eq(str(variant_id).lower())]
        if not match.empty:
            return match.iloc[0]
    if chrom and pd.notna(pos):
        match = annotation_df[(annotation_df.get("Chromosome", pd.Series(index=annotation_df.index, dtype=object)).astype(str) == str(chrom)) & (pd.to_numeric(annotation_df.get("POS", pd.Series(index=annotation_df.index, dtype=float)), errors="coerce") == pd.to_numeric(pos, errors="coerce"))]
        if ref:
            match = match[match.get("REF", pd.Series(index=match.index, dtype=object)).astype(str).str.upper().eq(str(ref).upper())]
        if alt:
            match = match[match.get("ALT", pd.Series(index=match.index, dtype=object)).astype(str).str.upper().eq(str(alt).upper())]
        if not match.empty:
            return match.iloc[0]
    return pd.Series(dtype=object)


def _is_modifier_impact(value: object) -> bool:
    return bool(re.search(r"\bmodifier\b", s(value), flags=re.IGNORECASE))


def _is_alamut_shortlist_candidate(consequence: object, impact: object) -> bool:
    consequence_text = s(consequence).lower()
    return "intron" in consequence_text or "splice" in consequence_text or _is_modifier_impact(impact)


def _best_support_label(fdr_value: object) -> str:
    fdr_num = pd.to_numeric(fdr_value, errors="coerce")
    if pd.notna(fdr_num) and float(fdr_num) < FDR_THRESHOLD:
        return "FDR-supported"
    return "Nominal only"


def _stage17_best_p_fdr(row: pd.Series) -> tuple[float, float]:
    if s(row.get("Analysis_Type")) == "Tumour_vs_Control":
        return pd.to_numeric(row.get("P_Value"), errors="coerce"), pd.to_numeric(row.get("FDR_P_Value"), errors="coerce")
    if s(row.get("Analysis_Type")) == "Cancer_Risk":
        return (
            _min_numeric_across(row, ["P_Adj_Age", "P_Unadj", "P_Hom_vs_WT", "P_Het_vs_WT"]),
            _min_numeric_across(row, ["FDR_Adj_Age", "FDR_Unadj", "FDR_Hom_vs_WT", "FDR_Het_vs_WT"]),
        )
    return (
        _min_numeric_across(row, ["P_Unadj", "P_Adj_Age", "P_Adj_BMI", "P_Adj_AgeBMI"]),
        _min_numeric_across(row, ["FDR_Unadj", "FDR_Adj_Age", "FDR_Adj_BMI", "FDR_Adj_AgeBMI"]),
    )


def _stage17_context(row: pd.Series) -> str:
    analysis_type = s(row.get("Analysis_Type"))
    if analysis_type == "Tumour_vs_Control":
        return f"Stage 17 tumour vs control ({s(row.get('Analysis_Group'))})"
    if analysis_type == "Cancer_Risk":
        return f"Stage 17 cancer risk ({s(row.get('Cohort'))})"
    cohort = s(row.get("Cohort")) or s(row.get("Analysis_Group"))
    clin_label = s(row.get("Clin_Label")) or "clinical association"
    return f"Stage 17 clinical association ({cohort}: {clin_label})"


def build_intronic_alamut_tables(stage12_all: pd.DataFrame, stage17: pd.DataFrame, stage17b_summary: pd.DataFrame, stage17b_catalogue: pd.DataFrame, annotation_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []

    if not stage12_all.empty:
        stage12_intronic = stage12_all[stage12_all.get("Raw_Significant", pd.Series(index=stage12_all.index)).fillna(False)].copy()
        for _, row in stage12_intronic.iterrows():
            ann = _lookup_annotation_row(annotation_df, variant_id=s(row.get("SNP_ID")), rsid=s(row.get("rsID")))
            rsid = s(row.get("rsID")) or s(ann.get("rsID"))
            chrom = s(ann.get("Chromosome"))
            pos = ann.get("POS")
            ref = s(ann.get("REF")).upper()
            alt = s(ann.get("ALT")).upper()
            gene = _display_gene_name(row.get("Symbol"), ann.get("Gene"))
            consequence = s(row.get("Consequence")) or s(ann.get("Consequence"))
            impact = s(row.get("Impact")) or s(ann.get("Impact"))
            if not _is_alamut_shortlist_candidate(consequence, impact):
                continue
            hgvsc = s(ann.get("HGVSc"))
            splice_note = s(ann.get("Splicing_Annotation")) or hgvsc
            preferred = _preferred_variant_label(rsid, gene, hgvsc, chrom, pos, ref, alt, fallback=s(row.get("Display_Label")) or s(row.get("SNP_ID")))
            stable_key = f"{chrom}:{int(pos)}:{ref}>{alt}" if chrom and pd.notna(pos) and ref and alt else preferred
            rows.append({
                "Stable Variant Key": stable_key,
                "Variant Type": "Common SNP",
                "Preferred Label": preferred,
                "rsID": rsid,
                "Gene": gene,
                "Chromosome": chrom,
                "Position": int(pos) if pd.notna(pos) else pd.NA,
                "REF": ref,
                "ALT": alt,
                "Consequence": consequence,
                "Impact": impact,
                "Analysis Source": "Stage 12",
                "Cohort / Analysis Context": s(row.get("Analysis_Group")),
                "Analysis Context": f"Stage 12 tumour vs pooled healthy controls ({s(row.get('Analysis_Group'))})",
                "P-Value": pd.to_numeric(row.get("P_Value"), errors="coerce"),
                "FDR": pd.to_numeric(row.get("FDR_P_Value"), errors="coerce"),
                "Significance Support": _best_support_label(row.get("FDR_P_Value")),
                "Splice / HGVSc Annotation": splice_note,
                "Coverage Risk Flag": bool(row.get("Coverage_Risk_Flag", False)),
                "Coverage Risk Note": s(row.get("Coverage_Risk_Note")),
            })

    if not stage17.empty:
        stage17_intronic = stage17.copy()
        for _, row in stage17_intronic.iterrows():
            ann = _lookup_annotation_row(annotation_df, variant_id=s(row.get("Variant_ID")), rsid=s(row.get("Variant_ID")))
            rsid = s(ann.get("rsID")) or (s(row.get("Variant_ID")) if re.fullmatch(r"rs\d+", s(row.get("Variant_ID")), flags=re.IGNORECASE) else "")
            chrom = s(ann.get("Chromosome"))
            pos = ann.get("POS")
            ref = s(ann.get("REF")).upper()
            alt = s(ann.get("ALT")).upper()
            gene = _display_gene_name(row.get("Gene"), ann.get("Gene"))
            consequence = s(row.get("Consequence")) or s(ann.get("Consequence"))
            impact = s(row.get("IMPACT")) or s(ann.get("Impact"))
            if not _is_alamut_shortlist_candidate(consequence, impact):
                continue
            hgvsc = s(ann.get("HGVSc"))
            splice_note = s(ann.get("Splicing_Annotation")) or hgvsc
            best_p, best_fdr = _stage17_best_p_fdr(row)
            preferred = _preferred_variant_label(rsid, gene, hgvsc, chrom, pos, ref, alt, fallback=s(row.get("Variant_ID")))
            stable_key = f"{chrom}:{int(pos)}:{ref}>{alt}" if chrom and pd.notna(pos) and ref and alt else preferred
            rows.append({
                "Stable Variant Key": stable_key,
                "Variant Type": "Common SNP",
                "Preferred Label": preferred,
                "rsID": rsid,
                "Gene": gene,
                "Chromosome": chrom,
                "Position": int(pos) if pd.notna(pos) else pd.NA,
                "REF": ref,
                "ALT": alt,
                "Consequence": consequence,
                "Impact": impact,
                "Analysis Source": "Stage 17",
                "Cohort / Analysis Context": s(row.get("Cohort")) or s(row.get("Analysis_Group")),
                "Analysis Context": _stage17_context(row),
                "P-Value": best_p,
                "FDR": best_fdr,
                "Significance Support": _best_support_label(best_fdr),
                "Splice / HGVSc Annotation": splice_note,
                "Coverage Risk Flag": bool(row.get("Coverage_Risk_Flag", False)),
                "Coverage Risk Note": s(row.get("Coverage_Risk_Note")),
            })

    if not stage17b_summary.empty and not stage17b_catalogue.empty:
        rare_summary = stage17b_summary[stage17b_summary.get("Exposure_Type", pd.Series(index=stage17b_summary.index, dtype=object)).astype(str).eq("Rare_Variant")].copy()
        rare_catalogue = stage17b_catalogue[[c for c in ["Variant_Key", "Variant_Label", "Variant_ID", "Gene", "CHROM", "POS", "REF", "ALT", "Consequence", "IMPACT"] if c in stage17b_catalogue.columns]].drop_duplicates("Variant_Key")
        rare_summary = rare_summary.merge(rare_catalogue, left_on="Exposure_ID", right_on="Variant_Key", how="left", suffixes=("", "_catalogue"))
        for _, row in rare_summary.iterrows():
            rsid = s(row.get("Variant_ID")) if re.fullmatch(r"rs\d+", s(row.get("Variant_ID")), flags=re.IGNORECASE) else ""
            chrom = s(row.get("CHROM"))
            pos = pd.to_numeric(row.get("POS"), errors="coerce")
            ref = s(row.get("REF")).upper()
            alt = s(row.get("ALT")).upper()
            ann = _lookup_annotation_row(annotation_df, variant_id=rsid or s(row.get("Variant_Label")), rsid=rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
            gene = _display_gene_name(row.get("Gene"), ann.get("Gene"))
            consequence = s(row.get("Consequence")) or s(ann.get("Consequence"))
            impact = s(row.get("IMPACT")) or s(ann.get("Impact"))
            if not _is_alamut_shortlist_candidate(consequence, impact):
                continue
            hgvsc = s(ann.get("HGVSc"))
            splice_note = s(ann.get("Splicing_Annotation")) or hgvsc
            preferred = _preferred_variant_label(rsid, gene, hgvsc, chrom, pos, ref, alt, fallback=s(row.get("Variant_Label")) or s(row.get("Exposure_ID")))
            stable_key = f"{chrom}:{int(pos)}:{ref}>{alt}" if chrom and pd.notna(pos) and ref and alt else preferred
            rows.append({
                "Stable Variant Key": stable_key,
                "Variant Type": "Rare variant / mutation",
                "Preferred Label": preferred,
                "rsID": rsid,
                "Gene": gene,
                "Chromosome": chrom,
                "Position": int(pos) if pd.notna(pos) else pd.NA,
                "REF": ref,
                "ALT": alt,
                "Consequence": consequence,
                "Impact": impact,
                "Analysis Source": "Stage 17b",
                "Cohort / Analysis Context": s(row.get("Cohort")) or s(row.get("Analysis_Group")),
                "Analysis Context": f"Stage 17b {s(row.get('Analysis')).replace('_', ' ').strip()} ({s(row.get('Cohort')) or s(row.get('Analysis_Group'))})",
                "P-Value": pd.to_numeric(row.get("P_Primary"), errors="coerce"),
                "FDR": pd.to_numeric(row.get("FDR_Primary"), errors="coerce"),
                "Significance Support": _best_support_label(row.get("FDR_Primary")),
                "Splice / HGVSc Annotation": splice_note,
                "Coverage Risk Flag": bool(ann.get("Coverage_Risk_Flag", False)),
                "Coverage Risk Note": s(ann.get("Coverage_Risk_Note")),
            })

    if not rows:
        note = pd.DataFrame({"Notes": ["No significant intronic, splice-related, or other modifier-impact SNPs / mutations were available for Alamut shortlisting."]})
        return note.copy(), note.copy(), note

    analysis_context = pd.DataFrame(rows)
    support_rank = {"FDR-supported": 0, "Nominal only": 1}
    analysis_context["Support Rank"] = analysis_context["Significance Support"].map(support_rank).fillna(9)
    analysis_context = analysis_context.sort_values(["Support Rank", "P-Value", "Variant Type", "Preferred Label", "Analysis Context"], kind="stable").reset_index(drop=True)

    summary_rows = []
    for stable_key, sub in analysis_context.groupby("Stable Variant Key", sort=False):
        summary_rows.append({
            "Preferred Label": _first_nonempty_text(sub["Preferred Label"]),
            "rsID": _first_nonempty_text(sub["rsID"]),
            "Gene": join_unique(sub["Gene"]),
            "Chromosome": _first_nonempty_text(sub["Chromosome"]),
            "Position": pd.to_numeric(sub["Position"], errors="coerce").dropna().astype(int).iloc[0] if pd.to_numeric(sub["Position"], errors="coerce").notna().any() else pd.NA,
            "REF": _first_nonempty_text(sub["REF"]),
            "ALT": _first_nonempty_text(sub["ALT"]),
            "Consequence": join_unique(sub["Consequence"]),
            "Impact": join_unique(sub["Impact"]),
            "Variant Type": join_unique(sub["Variant Type"]),
            "Significant Contexts": " | ".join(dict.fromkeys(sub["Analysis Context"].astype(str).tolist())),
            "Best P-Value": pd.to_numeric(sub["P-Value"], errors="coerce").min(),
            "Best FDR": pd.to_numeric(sub["FDR"], errors="coerce").min(),
            "Best Support": "FDR-supported" if (sub["Significance Support"] == "FDR-supported").any() else "Nominal only",
            "Splice / HGVSc Annotation": join_unique(sub["Splice / HGVSc Annotation"]),
            "Coverage Risk Flag": bool(sub["Coverage Risk Flag"].fillna(False).any()),
            "Coverage Risk Note": join_unique(sub["Coverage Risk Note"]),
        })
    variant_summary = pd.DataFrame(summary_rows).sort_values(["Best FDR", "Best P-Value", "Preferred Label"], kind="stable").reset_index(drop=True)

    notes = pd.DataFrame({
        "Notes": [
            "The Variant Shortlist sheet collapses repeated analysis rows to one row per unique intronic, splice-related, or other modifier-impact variant.",
            "The Analysis Context sheet preserves the cohort / analysis source, P-value, FDR, and support level for each significant context.",
            "FDR-supported rows should be prioritised first for manual Alamut follow-up; nominal-only rows remain useful but need more caution.",
            "Where no rsID is available, the preferred label falls back to gene + HGVSc when available, otherwise gene + genomic position.",
        ]
    })
    return variant_summary, analysis_context.drop(columns=["Stable Variant Key", "Support Rank"]), notes


def build_corrected_outputs_table() -> pd.DataFrame:
    return pd.DataFrame([
        {"Output": "GSDMB_SNP_Enrichment_Results.xlsx", "Where To Look": "analysis_results/12_snp_enrichment", "What It Shows": "Common-SNP enrichment results with rsID-aware display labels and a concise significant-overview sheet.", "Main Update": "Volcano labels and shareable tables now prefer rsIDs, with a consistent fallback label when no rsID exists."},
        {"Output": "GSDMB_SNP_Enrichment_Volcano_RAW.png / GSDMB_SNP_Enrichment_Volcano_FDR.png", "Where To Look": "analysis_results/12_snp_enrichment", "What It Shows": "Stage-12 enrichment volcano plots for Breast, Endometrium, and pooled Global comparisons.", "Main Update": "Titles, axis labels, and legends now use consistent title case and rsID-aware labels."},
        {"Output": "GSDMB_SNP_Clinical_Association_Results.xlsx", "Where To Look": "analysis_results/17_snp_clinical_associations", "What It Shows": "Clinical, cancer-risk, genotype-dose, and tumour-vs-control SNP associations.", "Main Update": "A new Significant Overview sheet collapses the wide mixed-analysis summary into a supervisor-facing variant overview."},
        {"Output": "22_Haplotype_First_Interpretation.xlsx / 22_Haplotype_First_Genomic_Map.png", "Where To Look": "analysis_results/22_haplotype_first_interpretation", "What It Shows": "Primary haplotypes, linked significant SNPs, and the explicit LD-block map audit.", "Main Update": "Figure titles and wording were standardised, and all valid study LD blocks remain visible on the shared map track."},
        {"Output": "23_RNA_Genetic_Integration.xlsx", "Where To Look": "analysis_results/23_rna_integration", "What It Shows": "Integrated RNA associations for significant SNPs, mutations, and haplotypes, including the direct GSDMB isoform answer.", "Main Update": "A new GSDMB Isoform Summary sheet collapses repeated LD-equivalent exposure rows for easier supervisor review."},
        {"Output": "25_rs11078928_rs869402_haplotype_focus.xlsx / 25_Supervisor_Summary.tsv", "Where To Look": "analysis_results/25_rs11078928_rs869402_haplotype_focus", "What It Shows": "Focused rs11078928 / rs869402 genotype and mini-haplotype comparisons.", "Main Update": "A clean Supervisor Summary now replaces the previously stale manual TSV and keeps the pooled tumour-vs-pooled normal wording explicit."},
        {"Output": OUTPUT_INTRONIC_ALAMUT_XLSX + " / " + OUTPUT_INTRONIC_ALAMUT_TSV, "Where To Look": "analysis_results/26_snp_functional_interpretation", "What It Shows": "Shortlisted significant intronic, splice-related, and other modifier-impact common SNPs and rare variants for manual Alamut follow-up.", "Main Update": "Brings together stage-12, stage-17, and stage-17b significance with the best available splice / HGVSc annotation, while retaining modifier-impact variants even when they are not explicitly labelled intronic."},
    ])


def build_supervisor_summary_tables(intronic_variant_summary: pd.DataFrame, stage23_isoform_summary: pd.DataFrame, stage25_supervisor: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    corrected_outputs = build_corrected_outputs_table()

    finding_rows = []
    if not intronic_variant_summary.empty and "Preferred Label" in intronic_variant_summary.columns:
        top_intronic = intronic_variant_summary.head(6)
        for _, row in top_intronic.iterrows():
            finding_rows.append({
                "Theme": "Alamut shortlist",
                "Finding": f"{row['Preferred Label']} ({row['Gene']}) - {row['Best Support']}; best P={fmt_num(row['Best P-Value'])}; best FDR={fmt_num(row['Best FDR'])}; contexts: {row['Significant Contexts']}",
            })

    if not stage23_isoform_summary.empty:
        if "Summary Status" in stage23_isoform_summary.columns:
            work = stage23_isoform_summary.copy()
            subset = work[work["Summary Status"].astype(str).eq("FDR_significant")].copy()
            if subset.empty:
                subset = work.head(4)
            for _, row in subset.head(4).iterrows():
                finding_rows.append({
                    "Theme": "Isoform-expression integration",
                    "Finding": f"{row.get('Representative Exposure', row.get('Exposure_Label_Final', 'Exposure'))} in {row.get('Context', 'study')} - {row.get('RNA Label', row.get('RNA_Label', 'RNA readout'))}; P={fmt_num(row.get('P-Value', row.get('P_Value')))}; FDR={fmt_num(row.get('FDR'))}; {row.get('Interpretation', '')}",
                })
        else:
            subset = stage23_isoform_summary[stage23_isoform_summary.get("Summary_Status", pd.Series(index=stage23_isoform_summary.index, dtype=object)).astype(str).eq("FDR_significant")].copy()
            if subset.empty:
                subset = stage23_isoform_summary.head(4)
            for _, row in subset.head(4).iterrows():
                finding_rows.append({
                    "Theme": "Isoform-expression integration",
                    "Finding": f"{row.get('Exposure_Label_Final', 'Exposure')} in {row.get('Context', 'study')} - {row.get('RNA_Label', 'RNA readout')}; P={fmt_num(row.get('P_Value'))}; FDR={fmt_num(row.get('FDR'))}; {row.get('Interpretation', '')}",
                })

    if not stage25_supervisor.empty:
        pooled = stage25_supervisor[stage25_supervisor.get("Comparison", pd.Series(index=stage25_supervisor.index, dtype=object)).astype(str).eq("Pooled tumour vs pooled normal")].copy()
        if not pooled.empty:
            row = pooled.iloc[0]
            finding_rows.append({
                "Theme": "Focused rs11078928 / rs869402 comparison",
                "Finding": f"Pooled tumour vs pooled normal: {row.get('rs11078928 Summary', row.get('rs11078928 summary', ''))} {row.get('rs869402 Summary', row.get('rs869402 summary', ''))}",
            })
            finding_rows.append({
                "Theme": "Mini-haplotype comparison",
                "Finding": f"Pooled mini-haplotype result: {row.get('Mini-Haplotype Readout', row.get('Mini-haplotype readout', ''))}. {row.get('Mini-Haplotype Support', row.get('Mini-haplotype support', ''))}",
            })

    if not finding_rows:
        finding_rows.append({"Theme": "Summary", "Finding": "No refreshed finding rows were available when the supervisor summary workbook was built."})

    caveat_rows = [
        {"Caveat": "FDR-supported results should be prioritised over nominal-only rows in the Alamut shortlist."},
        {"Caveat": "Variants without rsIDs still need fallback labels based on gene + HGVSc or gene + genomic position."},
        {"Caveat": "Splice / nearest-exon style annotation depends on what is already present in the current VEP-derived annotation workbook; some variants still only have broad non-coding or intronic labels."},
    ]
    if not intronic_variant_summary.empty and "Coverage Risk Flag" in intronic_variant_summary.columns and intronic_variant_summary["Coverage Risk Flag"].fillna(False).any():
        caveat_rows.append({"Caveat": "At least one shortlisted Alamut variant overlaps a flagged technical-risk amplicon; review the coverage-risk note before manual interpretation."})

    return corrected_outputs, pd.DataFrame(finding_rows), pd.DataFrame(caveat_rows)


def main() -> None:
    annot = load_annotation_backbone(PATHS["annot"])
    curated = load_curated_evidence(PATHS["curated_evidence"])
    stage12_all, stage12_sig = load_stage12(PATHS["stage12"])
    stage17 = load_stage17(PATHS["stage17"])
    stage20 = load_stage20(PATHS["stage20"])
    stage20b = load_stage20b(PATHS["stage20b"])
    stage17b_summary, stage17b_catalogue = load_stage17b(PATHS["stage17b"])
    stage23_isoform_summary = load_stage23_isoform_summary(PATHS["stage23"])
    stage25_supervisor = load_stage25_supervisor_summary(PATHS["stage25"])
    mapping = load_mapping_labels(PATHS["stage08"], PATHS["stage09"])
    _, stage22_links = load_stage22(PATHS["stage22"])
    stage24_snp, stage24_hap = load_stage24(PATHS["stage24"])

    priority = build_priority_table(annot, stage12_sig, stage17, stage20, stage20b, mapping)
    if priority.empty:
        raise SystemExit("No prioritised SNPs were found from the available project outputs.")

    mapping_summary = mapping.groupby("SNP_ID", as_index=False).agg(
        Mapping_Source=("Mapping_Source", join_unique),
        Mapping_Context=("Mapping_Context", join_unique),
    ) if not mapping.empty else pd.DataFrame(columns=["SNP_ID", "Mapping_Source", "Mapping_Context"])

    stage12_summary = stage12_sig.groupby("SNP_ID", as_index=False).agg(
        Stage12_Comparisons=("Comparison", join_unique),
        Stage12_Best_FDR=("FDR_P_Value", "min"),
    ) if not stage12_sig.empty else pd.DataFrame(columns=["SNP_ID", "Stage12_Comparisons", "Stage12_Best_FDR"])

    stage17_summary = stage17.groupby("Variant_ID", as_index=False).agg(
        Stage17_Analyses=("Analysis_Type", join_unique),
        Stage17_Best_P=("P_Unadj", "min"),
    ).rename(columns={"Variant_ID": "SNP_ID"}) if not stage17.empty else pd.DataFrame(columns=["SNP_ID", "Stage17_Analyses", "Stage17_Best_P"])

    integrated_rows: list[dict[str, object]] = []
    for item in priority.itertuples(index=False):
        snp_id = item.SNP_ID
        annot_match = annot[annot["SNP_ID"].eq(snp_id)]
        if annot_match.empty:
            annotation_row = pd.Series({
                "SNP_ID": snp_id,
                "Position": "Unknown",
                "Gene": "Unknown",
                "Consequence": "Unknown",
                "Consequence_All": "",
                "Population_Frequency_Context": "Not available in annotation workbook",
                "Regulatory_Annotation": "",
                "Protein_Annotation": "",
                "Splicing_Annotation": "",
                "Published_PMIDs": "",
                "Published_DB_Notes": "",
                "Coverage_Risk_Flag": False,
                "Coverage_Risk_Note": "",
            })
        else:
            annotation_row = annot_match.iloc[0]

        snp_evidence = evidence_for_snp(curated, snp_id)
        gsdmb_expr = summarise_evidence_rows(
            snp_evidence,
            {"GSDMB_expression"},
            "No direct curated evidence for a GSDMB expression effect is currently loaded.",
        )
        other_expr = summarise_evidence_rows(
            snp_evidence,
            {"Other_gene_expression"},
            "No direct curated evidence for an effect on other regional genes is currently loaded.",
        )
        splicing = summarise_evidence_rows(
            snp_evidence,
            {"Splicing_or_isoform"},
            s(annotation_row.get("Splicing_Annotation")) or "No direct splicing or isoform evidence is currently loaded.",
        )
        stage12_text = summarise_stage12_hits(snp_id, stage12_sig)
        stage17_text = summarise_stage17_hits(snp_id, stage17)
        rna_text = summarise_rna_hits(snp_id, stage20, stage20b)
        hap_text, in_key_hap, hap_external_summary = summarise_haplotype_context(snp_id, stage22_links, stage24_hap)
        her2_text = summarise_her2(snp_id, annotation_row, stage17, stage12_sig, curated)
        external_text = summarise_external_support(snp_id, stage24_snp, hap_external_summary)
        published_text = summarise_published_relevance(annotation_row, snp_evidence)
        other_function = summarise_other_function(annotation_row, snp_evidence)
        overall_text, biology_score = classify_overall_interpretation(
            snp_id,
            annotation_row,
            stage12_sig,
            stage17,
            stage20,
            stage20b,
            in_key_hap,
            snp_evidence,
        )
        significant_in_project = any([
            not stage12_sig[stage12_sig["SNP_ID"].eq(snp_id)].empty,
            not stage17[stage17["Variant_ID"].eq(snp_id)].empty,
            not stage20[stage20["Variant_ID"].eq(snp_id)].empty,
            not stage20b[stage20b["Variant_ID"].eq(snp_id)].empty,
        ])

        integrated_rows.append({
            "SNP_ID": snp_id,
            "Priority_Reasons": item.Priority_Reasons,
            "Position": s(annotation_row.get("Position")),
            "Gene": s(annotation_row.get("Gene")),
            "Consequence": s(annotation_row.get("Consequence")),
            "Significant_in_project": bool(significant_in_project),
            "In_key_haplotype": bool(in_key_hap),
            "Evidence_for_GSDMB_expression_effect": gsdmb_expr,
            "Evidence_for_other_gene_expression_effect": other_expr,
            "Evidence_for_splicing_or_isoform_effect": splicing,
            "Other_functional_consequences": other_function,
            "Published_relevance_summary": published_text,
            "Project_SNP_result": stage12_text,
            "Project_clinical_result": stage17_text,
            "RNA_or_isoform_result": rna_text,
            "Haplotype_context": hap_text,
            "HER2_relevance": her2_text,
            "Supported_by_1000G_comparison": external_text,
            "Overall_biological_interpretation": overall_text,
            "Population_Frequency_Context": s(annotation_row.get("Population_Frequency_Context")),
            "Coverage_Risk_Flag": bool(annotation_row.get("Coverage_Risk_Flag", False)),
            "Coverage_Risk_Note": s(annotation_row.get("Coverage_Risk_Note")),
            "Published_PMIDs": s(annotation_row.get("Published_PMIDs")),
            "Biology_Score": int(biology_score),
        })

    final_df = pd.DataFrame(integrated_rows).sort_values(["Biology_Score", "Significant_in_project", "SNP_ID"], ascending=[False, False, True], kind="stable").reset_index(drop=True)

    source_audit = build_source_audit(annot, stage12_sig, stage17, stage20, stage20b, mapping, stage22_links, stage24_snp, curated)
    thesis_summary_lines = build_thesis_summary(final_df, source_audit)
    thesis_summary_df = pd.DataFrame({"Summary": thesis_summary_lines})

    integrated_table = final_df[[
        "SNP_ID",
        "Gene",
        "Consequence",
        "Significant_in_project",
        "Project_SNP_result",
        "Project_clinical_result",
        "RNA_or_isoform_result",
        "Haplotype_context",
        "HER2_relevance",
        "Supported_by_1000G_comparison",
        "Overall_biological_interpretation",
    ]].copy()

    literature_table = final_df[[
        "SNP_ID",
        "Published_relevance_summary",
        "Evidence_for_GSDMB_expression_effect",
        "Evidence_for_other_gene_expression_effect",
        "Evidence_for_splicing_or_isoform_effect",
        "Other_functional_consequences",
        "Published_PMIDs",
    ]].copy()

    prioritised_table = priority.merge(mapping_summary, on="SNP_ID", how="left").merge(stage12_summary, on="SNP_ID", how="left").merge(stage17_summary, on="SNP_ID", how="left")

    out_workbook = OUT_DIR / OUTPUT_WORKBOOK
    with pd.ExcelWriter(out_workbook, engine="openpyxl") as writer:
        pd.DataFrame({
            "Notes": [
                "This stage prioritises the significant project SNP set first, then adds individually notable mapped SNPs.",
                "Functional interpretation remains layered on top of the existing stage-12, stage-17, stage-20, stage-22, and stage-24 outputs rather than replacing them.",
                "External haplotype validation against female-only 1000 Genomes remains preserved through stage 24 and linked-haplotype summaries.",
                "Curated literature claims come only from docs/curated_snp_functional_evidence.tsv; when that file is silent, the script falls back to auditable annotation fields and reports that limitation explicitly.",
                f"Project FDR threshold used downstream: {FDR_THRESHOLD:.2f}",
            ]
        }).to_excel(writer, sheet_name="README", index=False)
        prioritised_table.to_excel(writer, sheet_name="Prioritised_SNPs", index=False)
        annot.merge(priority[["SNP_ID"]], on="SNP_ID", how="inner").to_excel(writer, sheet_name="Annotation_Backbone", index=False)
        curated.to_excel(writer, sheet_name="Curated_Evidence", index=False)
        final_df[[
            "SNP_ID",
            "Position",
            "Gene",
            "Consequence",
            "Significant_in_project",
            "In_key_haplotype",
            "Evidence_for_GSDMB_expression_effect",
            "Evidence_for_other_gene_expression_effect",
            "Evidence_for_splicing_or_isoform_effect",
            "Other_functional_consequences",
            "Published_relevance_summary",
            "HER2_relevance",
            "Supported_by_1000G_comparison",
            "Overall_biological_interpretation",
        ]].to_excel(writer, sheet_name="SNP_Interpretation_Table", index=False)
        literature_table.to_excel(writer, sheet_name="Literature_Summary", index=False)
        integrated_table.to_excel(writer, sheet_name="Integrated_Project_Table", index=False)
        source_audit.to_excel(writer, sheet_name="Source_Audit", index=False)
        thesis_summary_df.to_excel(writer, sheet_name="Thesis_Summary", index=False)

    final_df.to_csv(OUT_DIR / OUTPUT_TABLE_TSV, sep="\t", index=False)
    (OUT_DIR / OUTPUT_SUMMARY_TXT).write_text("\n".join(thesis_summary_lines) + "\n", encoding="utf-8")

    intronic_variant_summary, intronic_analysis_context, intronic_notes = build_intronic_alamut_tables(
        stage12_all, stage17, stage17b_summary, stage17b_catalogue, annot
    )
    alamut_workbook = OUT_DIR / OUTPUT_INTRONIC_ALAMUT_XLSX
    with pd.ExcelWriter(alamut_workbook, engine="openpyxl") as writer:
        intronic_notes.to_excel(writer, sheet_name="README", index=False)
        intronic_variant_summary.to_excel(writer, sheet_name="Variant Shortlist", index=False)
        intronic_analysis_context.to_excel(writer, sheet_name="Analysis Context", index=False)
    intronic_analysis_context.to_csv(OUT_DIR / OUTPUT_INTRONIC_ALAMUT_TSV, sep="	", index=False)

    corrected_outputs, main_findings, caveats = build_supervisor_summary_tables(
        intronic_variant_summary, stage23_isoform_summary, stage25_supervisor
    )
    supervisor_summary_workbook = OUT_DIR / OUTPUT_SUPERVISOR_SUMMARY_XLSX
    with pd.ExcelWriter(supervisor_summary_workbook, engine="openpyxl") as writer:
        corrected_outputs.to_excel(writer, sheet_name="Corrected Outputs", index=False)
        main_findings.to_excel(writer, sheet_name="Main Findings", index=False)
        caveats.to_excel(writer, sheet_name="Caveats", index=False)

    print(f"Wrote workbook:           {out_workbook}")
    print(f"Wrote table:              {OUT_DIR / OUTPUT_TABLE_TSV}")
    print(f"Wrote summary:            {OUT_DIR / OUTPUT_SUMMARY_TXT}")
    print(f"Wrote Alamut workbook:    {alamut_workbook}")
    print(f"Wrote Alamut TSV:         {OUT_DIR / OUTPUT_INTRONIC_ALAMUT_TSV}")
    print(f"Wrote supervisor summary: {supervisor_summary_workbook}")


if __name__ == "__main__":
    main()
