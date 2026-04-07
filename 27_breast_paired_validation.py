#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a breast tumour-normal paired validation workbook and audit tables.

This is a separate validation layer for the declared paired breast-cancer cases
from the clinical workbook sheet ``TANDA HER2 MAMA T + NT``. It does not alter
or replace the main unpaired cohort analyses.

Primary matching logic:
1. Recover clinically declared tumour-normal breast pairs.
2. Resolve tumour SNP codes from the nomenclature workbook sheet ``MT-T_N``
   using the clinical ``C?DIGO BB BLOQUE TUMORAL`` field.
3. Use the harmonised master only as a secondary consistency check.
4. Search breast tumour / normal BAM and VCF workflow folders.
5. Treat pairs as validated only when workbook linkage exists and both tumour
   and declared paired-normal BAMs are present.

Important design choice:
Declared ``...-NT`` paired normal IDs are not mapped to the generic ``MN_*``
normal cohort unless an explicit crosswalk exists in the workbook or
filesystem. This prevents the paired validation layer from being mixed with the
separate healthy-control breast-normal arm.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import unicodedata
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parent
DOCS_DIR = ROOT / "docs"


def _iter_docs_search_dirs() -> list[Path]:
    seen: set[str] = set()
    candidates: list[Path] = []

    env_docs = os.environ.get("TFM_DOCS_DIR", "").strip()
    if env_docs:
        env_path = Path(env_docs)
        candidates.extend([env_path / "clinical variables-snps", env_path])

    candidates.extend([DOCS_DIR / "source_workbooks", DOCS_DIR / "derived_workbooks", DOCS_DIR])

    users_root = Path("/mnt/c/Users")
    if users_root.exists():
        for user_dir in sorted(users_root.iterdir()):
            if not user_dir.is_dir():
                continue
            downloads = user_dir / "Downloads"
            if downloads.exists():
                candidates.append(downloads)
            for onedrive_dir in sorted(user_dir.glob("OneDrive*")):
                docs_dir = onedrive_dir / "Documents" / "TFM" / "Docs"
                candidates.extend([docs_dir / "clinical variables-snps", docs_dir])

    existing: list[Path] = []
    for candidate in candidates:
        key = str(candidate)
        if key in seen or not candidate.exists():
            continue
        seen.add(key)
        existing.append(candidate)
    return existing


def _discover_doc_file(filename: str, fallback_dir: Path) -> Path:
    for search_dir in _iter_docs_search_dirs():
        candidate = search_dir / filename
        if candidate.exists():
            return candidate
    return fallback_dir / filename


DEFAULT_CLINICAL_XLSX = _discover_doc_file(
    "MUESTRAS + VARIABLES CLINICAS PROY SNPs GM(4).xlsx",
    DOCS_DIR / "source_workbooks",
)
DEFAULT_NOMENCLATURE_XLSX = _discover_doc_file(
    "Muestras SNPs nomenclaturas, equivalencias y cuantificaciones.xlsx",
    DOCS_DIR / "source_workbooks",
)
DEFAULT_MASTER_XLSX = DOCS_DIR / "derived_workbooks/MASTER_SNP_plus_clinical_HARMONISED.xlsx"
DEFAULT_KEY_SNP_TSV = ROOT / "analysis_results/25_rs11078928_rs869402_haplotype_focus/25_Target_SNP_Metadata.tsv"
DEFAULT_OUTDIR = ROOT / "analysis_results/27_breast_paired_validation"
DEFAULT_OUT_XLSX = DOCS_DIR / "derived_workbooks/BREAST_PAIRED_TUMOUR_NORMAL_VALIDATION.xlsx"

TUMOUR_BAM_DIRS = [
    ROOT / "breast/tumour/dna",
    ROOT / "breast/tumour/dna_qc/pass_bams",
]
NORMAL_BAM_DIRS = [
    ROOT / "breast/normal/dna",
    ROOT / "breast/normal/dna_qc/pass_bams",
]
TUMOUR_DNA_CALLS_DIR = ROOT / "breast/tumour/dna_calls"
NORMAL_DNA_CALLS_DIR = ROOT / "breast/normal/dna_calls"


def _to_str(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _clean_id(value: object) -> str:
    text = _to_str(value).replace("\xa0", " ").strip()
    text = re.sub(r"\s+", "", text)
    return text.upper()


def _normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(col).strip() for col in out.columns]
    return out


def _slugify(text: object) -> str:
    value = unicodedata.normalize("NFKD", _to_str(text))
    value = value.encode("ascii", "ignore").decode("ascii")
    value = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    return value


def _col_by_slug(df: pd.DataFrame, slug: str) -> str:
    slug_map = {_slugify(col): str(col) for col in df.columns}
    if slug not in slug_map:
        raise KeyError(f"Column slug not found: {slug}")
    return slug_map[slug]


def _dedupe_paths(paths: Iterable[Path]) -> list[Path]:
    seen: set[str] = set()
    unique: list[Path] = []
    for path in paths:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        unique.append(path)
    return unique


def _discover_bams(roots: Iterable[Path]) -> list[Path]:
    bam_paths: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        bam_paths.extend(path for path in root.rglob("*.bam") if path.is_file())
    return _dedupe_paths(sorted(bam_paths))


def _prefer_path(paths: list[Path]) -> Path:
    def sort_key(path: Path) -> tuple[int, int, str]:
        text = str(path)
        pass_score = 0 if "pass_bams" in text else 1
        depth_score = len(path.parts)
        return (pass_score, depth_score, text)

    return sorted(paths, key=sort_key)[0]


def _match_files(paths: list[Path], tokens: list[str]) -> list[Path]:
    clean_tokens = [token.lower() for token in tokens if token]
    patterns = [re.compile(r"(?<![A-Za-z0-9])" + re.escape(token) + r"(?![A-Za-z0-9])") for token in clean_tokens]
    hits: list[Path] = []
    for path in paths:
        name = path.name.lower()
        if any(pattern.search(name) for pattern in patterns):
            hits.append(path)
    return hits


def _snp_tokens(snp_code: str) -> list[str]:
    code = _clean_id(snp_code)
    if not code:
        return []

    tokens = {
        code,
        code.replace("SNP_", ""),
        code.replace("-", "_"),
        code.replace("SNP_", "").replace("-", "_"),
        code.replace("_REP", ""),
        code.replace("SNP_", "").replace("_REP", ""),
    }

    mt_match = re.search(r"MT[-_]T[-_](\d+)", code)
    if mt_match:
        number = mt_match.group(1)
        tokens.update(
            {
                f"MT-T_{number}",
                f"MT_T_{number}",
                f"SNP_MT-T_{number}",
                f"SNP_MT_T_{number}",
                f"DNA_MT-T_{number}",
                f"DNA_MT_T_{number}",
                f"DNA_SNP_MT-T_{number}",
                f"DNA_SNP_MT_T_{number}",
            }
        )
    return sorted(tokens)


def _find_tumour_bam(row: pd.Series, bam_paths: list[Path]) -> Path | None:
    tokens = _snp_tokens(row.get("tumour_snp_code", ""))
    tokens.extend(
        [
            _clean_id(row.get("tumour_sample_id", "")),
            _clean_id(row.get("case_id", "")),
        ]
    )
    matches = _match_files(bam_paths, tokens)
    return _prefer_path(matches) if matches else None


def _find_normal_bam(row: pd.Series, bam_paths: list[Path]) -> tuple[Path | None, str]:
    search_tokens = [
        _clean_id(row.get("normal_sample_id", "")),
        _clean_id(row.get("normal_pathology_block_id", "")),
        _clean_id(row.get("case_id", "")) + "-NT" if _clean_id(row.get("case_id", "")) else "",
    ]
    matches = _match_files(bam_paths, search_tokens)
    if matches:
        return _prefer_path(matches), "paired_normal_id_or_block_match"
    return None, "no_explicit_paired_normal_bam_match"


def _find_sample_vcf(bam_path: str | None, dna_calls_dir: Path) -> str:
    if not bam_path:
        return ""
    bam_name = Path(bam_path).stem
    candidate_dir = dna_calls_dir / bam_name
    if candidate_dir.exists():
        for file_name in ["variants.norm.vcf.gz", "TSVC_variants.vcf.gz"]:
            candidate = candidate_dir / file_name
            if candidate.exists():
                return str(candidate)
    return ""


def load_candidate_pairs(clinical_xlsx: Path) -> pd.DataFrame:
    df = pd.read_excel(clinical_xlsx, sheet_name="TANDA HER2 MAMA T + NT")
    df = _normalise_columns(df)
    df = df.rename(
        columns={
            "CÓDIGO BB CASO": "case_id",
            "DX": "diagnosis",
            "COPIAS HER2": "her2_copies",
            "BLOQUE MAMA TUMORAL": "tumour_pathology_block_id",
            "CÓDIGO BB BLOQUE TUMORAL": "tumour_sample_id",
            "BLOQUE MAMA NORMAL PAREADO": "normal_pathology_block_id",
            "CÓDIGO BB BLOQUE NORMAL": "normal_sample_id",
            "REALIZADO": "tumour_material_note",
            "REALIZADO.1": "normal_material_note",
            "REPETICIÓN DE CORTES ": "repeat_cut_note",
        }
    )

    keep_cols = [
        "NHC",
        "case_id",
        "diagnosis",
        "her2_copies",
        "tumour_pathology_block_id",
        "tumour_sample_id",
        "tumour_material_note",
        "normal_pathology_block_id",
        "normal_sample_id",
        "normal_material_note",
        "AGO (GESTACIONES, PARTOS, CESÁREAS, ABORTOS)",
        "BMI",
        "Menarquia",
        "Menopausia",
        "repeat_cut_note",
    ]
    available_cols = [col for col in keep_cols if col in df.columns]
    pairs = df[available_cols].copy()

    for col in ["case_id", "tumour_sample_id", "normal_sample_id"]:
        if col in pairs.columns:
            pairs[col] = pairs[col].map(_clean_id)

    for col in ["tumour_pathology_block_id", "normal_pathology_block_id"]:
        if col in pairs.columns:
            pairs[col] = pairs[col].map(_to_str)

    pairs = pairs[
        (pairs["case_id"] != "")
        & (pairs["tumour_sample_id"] != "")
        & (pairs["normal_sample_id"] != "")
    ].copy()
    pairs = pairs.drop_duplicates(subset=["case_id"], keep="first").reset_index(drop=True)
    pairs["codigo_bb_bloque_tumoral_tumour"] = pairs["tumour_sample_id"]
    pairs["codigo_bb_bloque_tumoral_normal"] = pairs["normal_sample_id"]
    return pairs


def load_workbook_tumour_lookup(nomenclature_xlsx: Path) -> pd.DataFrame:
    df = pd.read_excel(nomenclature_xlsx, sheet_name="MT-T_N")
    df = _normalise_columns(df)
    lookup = df[
        [
            _col_by_slug(df, "codigo_bb_bloque_tumoral"),
            _col_by_slug(df, "codigo_snps"),
            _col_by_slug(df, "ac_nucleico"),
        ]
    ].copy()
    lookup.columns = ["tumour_sample_id", "tumour_snp_code_workbook", "tumour_nucleic_acid_workbook"]
    lookup["tumour_sample_id"] = lookup["tumour_sample_id"].map(_clean_id)
    lookup["tumour_snp_code_workbook"] = lookup["tumour_snp_code_workbook"].map(_clean_id)
    lookup["tumour_nucleic_acid_workbook"] = lookup["tumour_nucleic_acid_workbook"].map(lambda x: _to_str(x).upper())
    lookup = lookup[
        (lookup["tumour_sample_id"] != "")
        & (lookup["tumour_snp_code_workbook"] != "")
        & lookup["tumour_nucleic_acid_workbook"].eq("DNA")
    ].copy()
    lookup = lookup.drop_duplicates(subset=["tumour_sample_id"], keep="first")
    return lookup


def load_workbook_normal_lookup(nomenclature_xlsx: Path) -> pd.DataFrame:
    df = pd.read_excel(nomenclature_xlsx, sheet_name="MN")
    df = _normalise_columns(df)
    lookup = df[
        [
            _col_by_slug(df, "codigo_noray_bb"),
            _col_by_slug(df, "codigo_snp_mn_mama_normal"),
            _col_by_slug(df, "ac_nucleico"),
        ]
    ].copy()
    lookup.columns = ["normal_sample_id", "normal_snp_code_workbook", "normal_nucleic_acid_workbook"]
    lookup["normal_sample_id"] = lookup["normal_sample_id"].map(_clean_id)
    lookup["normal_snp_code_workbook"] = lookup["normal_snp_code_workbook"].map(_clean_id)
    lookup["normal_nucleic_acid_workbook"] = lookup["normal_nucleic_acid_workbook"].map(lambda x: _to_str(x).upper())
    lookup = lookup[
        (lookup["normal_sample_id"] != "")
        & (lookup["normal_snp_code_workbook"] != "")
        & lookup["normal_nucleic_acid_workbook"].eq("DNA")
    ].copy()
    lookup = lookup.drop_duplicates(subset=["normal_sample_id"], keep="first")
    return lookup


def load_master_tumour_lookup(master_xlsx: Path) -> pd.DataFrame:
    master = pd.read_excel(master_xlsx, sheet_name="harmonised_plus_canon")
    master = _normalise_columns(master)

    mask = (
        master["tissue"].astype(str).eq("Breast")
        & master["sheet"].astype(str).eq("MT-T_N")
        & master["sample_id"].notna()
    )
    breast = master.loc[mask].copy()
    breast["sample_id"] = breast["sample_id"].map(_clean_id)
    breast["case_id"] = breast["case_id"].map(_clean_id)
    breast["nucleic_acid"] = breast["nucleic_acid"].astype(str).str.upper()

    sort_rank = breast["nucleic_acid"].map({"DNA": 0, "RNA": 1}).fillna(2)
    breast = breast.assign(_sort_rank=sort_rank).sort_values(["sample_id", "case_id", "_sort_rank"])
    breast = breast.drop_duplicates(subset=["sample_id"], keep="first")
    return breast[
        [
            "case_id",
            "sample_id",
            "snp_code",
            "nucleic_acid",
            "block_type",
            "clin_her2__paired_normal_block_id",
            "clin_her2__BLOQUE MAMA NORMAL PAREADO",
        ]
    ].rename(
        columns={
            "sample_id": "tumour_sample_id_master",
            "snp_code": "tumour_snp_code_master",
            "nucleic_acid": "tumour_nucleic_acid_master",
            "block_type": "tumour_block_type_master",
            "clin_her2__paired_normal_block_id": "normal_sample_id_master",
            "clin_her2__BLOQUE MAMA NORMAL PAREADO": "normal_pathology_block_id_master",
        }
    )


def build_pair_manifest(
    candidate_pairs: pd.DataFrame,
    workbook_tumour_lookup: pd.DataFrame,
    workbook_normal_lookup: pd.DataFrame,
    master_tumour_lookup: pd.DataFrame,
    tumour_bams: list[Path],
    normal_bams: list[Path],
) -> pd.DataFrame:
    manifest = candidate_pairs.merge(workbook_tumour_lookup, how="left", on="tumour_sample_id")
    manifest = manifest.merge(workbook_normal_lookup, how="left", on="normal_sample_id")
    manifest = manifest.merge(
        master_tumour_lookup,
        how="left",
        left_on="tumour_sample_id",
        right_on="tumour_sample_id_master",
        suffixes=("", "_master"),
    )

    manifest["tumour_snp_code"] = manifest["tumour_snp_code_workbook"].fillna("")
    manifest["normal_snp_code"] = manifest["normal_snp_code_workbook"].fillna("")

    manifest["workbook_match_status"] = "unresolved_in_workbook"
    manifest.loc[manifest["tumour_snp_code"].astype(str).ne(""), "workbook_match_status"] = (
        "tumour_snp_resolved_from_MT-T_N"
    )
    manifest.loc[
        manifest["tumour_snp_code"].astype(str).ne("")
        & manifest["normal_snp_code"].astype(str).ne(""),
        "workbook_match_status",
    ] = "tumour_and_normal_snp_resolved"

    manifest["master_case_id_matches"] = manifest["case_id"] == manifest["case_id_master"].fillna("")
    manifest["master_snp_code_matches"] = (
        manifest["tumour_snp_code"].astype(str) == manifest["tumour_snp_code_master"].fillna("").astype(str)
    )
    manifest["master_consistency_status"] = "master_not_linked"
    manifest.loc[manifest["tumour_sample_id_master"].notna(), "master_consistency_status"] = (
        "master_linked_case_or_snp_mismatch"
    )
    manifest.loc[
        manifest["tumour_sample_id_master"].notna()
        & manifest["master_case_id_matches"]
        & manifest["master_snp_code_matches"],
        "master_consistency_status",
    ] = "master_linked_and_consistent"

    tumour_paths: list[str] = []
    normal_paths: list[str] = []
    normal_rules: list[str] = []
    for _, row in manifest.iterrows():
        tumour_bam = _find_tumour_bam(row, tumour_bams)
        normal_bam, normal_rule = _find_normal_bam(row, normal_bams)
        tumour_paths.append(str(tumour_bam) if tumour_bam else "")
        normal_paths.append(str(normal_bam) if normal_bam else "")
        normal_rules.append(normal_rule)

    manifest["tumour_bam_path"] = tumour_paths
    manifest["normal_bam_path"] = normal_paths
    manifest["tumour_bam_found"] = manifest["tumour_bam_path"].ne("")
    manifest["normal_bam_found"] = manifest["normal_bam_path"].ne("")
    manifest["tumour_vcf_path"] = manifest["tumour_bam_path"].map(
        lambda path: _find_sample_vcf(path, TUMOUR_DNA_CALLS_DIR)
    )
    manifest["normal_vcf_path"] = manifest["normal_bam_path"].map(
        lambda path: _find_sample_vcf(path, NORMAL_DNA_CALLS_DIR)
    )

    manifest["matching_rule_used"] = (
        "clinical paired row + MT-T_N tumour Codigo BB bloque tumoral -> CODIGO SNPs + explicit paired-normal ID/block BAM search"
    )
    manifest["normal_bam_matching_rule"] = normal_rules

    manifest["bam_match_status"] = "no_bam_match"
    manifest.loc[manifest["tumour_bam_found"] & ~manifest["normal_bam_found"], "bam_match_status"] = "tumour_bam_only"
    manifest.loc[manifest["tumour_bam_found"] & manifest["normal_bam_found"], "bam_match_status"] = "bam_confirmed_pair"

    manifest["tissue_pair_status"] = "candidate_pair_declared"
    manifest.loc[manifest["workbook_match_status"].eq("unresolved_in_workbook"), "tissue_pair_status"] = (
        "candidate_pair_without_workbook_tumour_resolution"
    )
    manifest.loc[
        manifest["workbook_match_status"].ne("unresolved_in_workbook")
        & manifest["tumour_bam_found"]
        & ~manifest["normal_bam_found"],
        "tissue_pair_status",
    ] = "workbook_resolved_tumour_with_tumour_bam_only"
    manifest.loc[
        manifest["workbook_match_status"].ne("unresolved_in_workbook")
        & manifest["tumour_bam_found"]
        & manifest["normal_bam_found"],
        "tissue_pair_status",
    ] = "bam_confirmed_pair"

    manifest["matching_confidence"] = "low_unresolved"
    manifest.loc[manifest["workbook_match_status"].ne("unresolved_in_workbook"), "matching_confidence"] = (
        "high_workbook_resolved_tumour"
    )
    manifest.loc[
        manifest["workbook_match_status"].ne("unresolved_in_workbook")
        & manifest["tumour_bam_found"],
        "matching_confidence",
    ] = "high_workbook_resolved_tumour_with_bam"
    manifest.loc[
        manifest["workbook_match_status"].ne("unresolved_in_workbook")
        & manifest["tumour_bam_found"]
        & manifest["normal_bam_found"],
        "matching_confidence",
    ] = "high_bam_confirmed_pair"

    notes: list[str] = []
    for _, row in manifest.iterrows():
        row_notes: list[str] = []
        if _to_str(row.get("tumour_snp_code")):
            row_notes.append("Tumour SNP code resolved from MT-T_N using the clinical tumour Codigo BB bloque tumoral field.")
        else:
            row_notes.append("Tumour SNP code could not be resolved from MT-T_N.")
        if not _to_str(row.get("normal_snp_code")):
            row_notes.append(
                "No explicit paired-normal SNP-code crosswalk was found for the declared normal sample in the current workbook structure."
            )
        if not bool(row.get("normal_bam_found")):
            row_notes.append(
                "No explicit paired-normal BAM/VCF was found from the declared ...-NT ID or paired normal block. Generic MN_* healthy-normal files were not treated as matched pairs."
            )
        if not bool(row.get("tumour_bam_found")):
            row_notes.append(
                "Tumour BAM could not be resolved in the current breast tumour workflow folders."
            )
        if _to_str(row.get("master_consistency_status")) == "master_linked_case_or_snp_mismatch":
            row_notes.append(
                "Master linkage exists but the case ID or SNP code did not fully match the workbook-resolved tumour mapping."
            )
        notes.append(" ".join(row_notes))
    manifest["notes"] = notes

    return manifest


def load_key_snp_metadata(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["rsID_clean", "POS_int", "Gene", "REF", "ALT"])
    df = pd.read_csv(path, sep="\t")
    keep = [col for col in ["rsID_clean", "POS_int", "Gene", "REF", "ALT"] if col in df.columns]
    return df[keep].copy()


def _parse_single_sample_vcf(vcf_path: str) -> pd.DataFrame:
    if not vcf_path:
        return pd.DataFrame(columns=["CHROM", "POS", "ID", "REF", "ALT", "GT"])

    records: list[dict[str, object]] = []
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            sample_field = fields[9]
            fmt_keys = fields[8].split(":")
            fmt_vals = sample_field.split(":")
            gt = ""
            if "GT" in fmt_keys:
                gt_index = fmt_keys.index("GT")
                if gt_index < len(fmt_vals):
                    gt = fmt_vals[gt_index]
            records.append(
                {
                    "CHROM": fields[0],
                    "POS": int(fields[1]),
                    "ID": fields[2],
                    "REF": fields[3],
                    "ALT": fields[4],
                    "GT": gt,
                }
            )
    return pd.DataFrame(records)


def build_concordance_tables(
    bam_confirmed_pairs: pd.DataFrame,
    key_snp_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    concordance_rows: list[dict[str, object]] = []
    key_detail_rows: list[dict[str, object]] = []

    for _, row in bam_confirmed_pairs.iterrows():
        tumour_vcf = _parse_single_sample_vcf(_to_str(row.get("tumour_vcf_path")))
        normal_vcf = _parse_single_sample_vcf(_to_str(row.get("normal_vcf_path")))
        if tumour_vcf.empty or normal_vcf.empty:
            concordance_rows.append(
                {
                    "case_id": row["case_id"],
                    "tumour_sample_id": row["tumour_sample_id"],
                    "normal_sample_id": row["normal_sample_id"],
                    "SNP_concordance_rate": pd.NA,
                    "key_SNP_discordances": "VCF missing for tumour or normal",
                    "haplotype_concordance": pd.NA,
                    "key_haplotype_discordances": "Not assessable",
                    "likely_explanation": "Missing VCF",
                    "validation_interpretation": "Pair had BAM support but no usable VCF comparison.",
                    "review_flag": True,
                }
            )
            continue

        tumour_key = tumour_vcf.assign(key=lambda df: df["CHROM"].astype(str) + ":" + df["POS"].astype(str))
        normal_key = normal_vcf.assign(key=lambda df: df["CHROM"].astype(str) + ":" + df["POS"].astype(str))

        merged = tumour_key[["key", "GT"]].rename(columns={"GT": "tumour_GT"}).merge(
            normal_key[["key", "GT"]].rename(columns={"GT": "normal_GT"}),
            on="key",
            how="inner",
        )
        usable = merged[
            merged["tumour_GT"].astype(str).ne("")
            & merged["normal_GT"].astype(str).ne("")
            & merged["tumour_GT"].astype(str).ne("./.")
            & merged["normal_GT"].astype(str).ne("./.")
        ].copy()

        if usable.empty:
            snp_concordance_rate = pd.NA
            discordance_note = "No overlapping callable variant genotypes in both VCFs"
            review_flag = True
        else:
            usable["match"] = usable["tumour_GT"] == usable["normal_GT"]
            snp_concordance_rate = float(usable["match"].mean())
            discordant_keys = usable.loc[~usable["match"], "key"].tolist()
            discordance_note = "; ".join(discordant_keys[:10]) if discordant_keys else ""
            review_flag = bool(discordant_keys)

        key_discordances: list[str] = []
        for _, key_row in key_snp_metadata.iterrows():
            chrom_key = f"chr17:{int(key_row['POS_int'])}"
            tumour_hit = tumour_key.loc[tumour_key["key"] == chrom_key]
            normal_hit = normal_key.loc[normal_key["key"] == chrom_key]
            tumour_gt = tumour_hit["GT"].iloc[0] if not tumour_hit.empty else ""
            normal_gt = normal_hit["GT"].iloc[0] if not normal_hit.empty else ""
            discordant = bool(tumour_gt and normal_gt and tumour_gt != normal_gt)
            if discordant:
                key_discordances.append(f"{key_row['rsID_clean']}:{tumour_gt}!={normal_gt}")
                review_flag = True
            key_detail_rows.append(
                {
                    "case_id": row["case_id"],
                    "tumour_sample_id": row["tumour_sample_id"],
                    "normal_sample_id": row["normal_sample_id"],
                    "rsID_clean": key_row["rsID_clean"],
                    "POS_int": key_row["POS_int"],
                    "tumour_GT": tumour_gt,
                    "normal_GT": normal_gt,
                    "discordant": discordant,
                }
            )

        concordance_rows.append(
            {
                "case_id": row["case_id"],
                "tumour_sample_id": row["tumour_sample_id"],
                "normal_sample_id": row["normal_sample_id"],
                "SNP_concordance_rate": snp_concordance_rate,
                "key_SNP_discordances": "; ".join(key_discordances) if key_discordances else discordance_note,
                "haplotype_concordance": pd.NA,
                "key_haplotype_discordances": (
                    "Not assessed from current files; no phased paired-normal haplotype layer was available."
                ),
                "likely_explanation": (
                    "Major discordances require manual review for somatic/LOH/coverage/sample issues."
                    if review_flag
                    else "No genotype discordance observed across overlapping callable sites."
                ),
                "validation_interpretation": (
                    "Concordant across overlapping callable sites."
                    if not review_flag
                    else "Discordance detected and needs manual review."
                ),
                "review_flag": review_flag,
            }
        )

    concordance_df = pd.DataFrame(
        concordance_rows,
        columns=[
            "case_id",
            "tumour_sample_id",
            "normal_sample_id",
            "SNP_concordance_rate",
            "key_SNP_discordances",
            "haplotype_concordance",
            "key_haplotype_discordances",
            "likely_explanation",
            "validation_interpretation",
            "review_flag",
        ],
    )
    key_detail_df = pd.DataFrame(key_detail_rows)
    return concordance_df, key_detail_df


def build_overview_table(manifest: pd.DataFrame, workbook_resolved: pd.DataFrame, confirmed: pd.DataFrame) -> pd.DataFrame:
    candidate_pairs = int(len(manifest))
    workbook_resolved_pairs = int(len(workbook_resolved))
    tumour_bam_count = int(manifest["tumour_bam_found"].sum())
    normal_bam_count = int(manifest["normal_bam_found"].sum())
    confirmed_pairs = int(len(confirmed))

    support_statement = (
        "Not assessable from matched tumour-normal BAM pairs in the current repo state."
        if confirmed_pairs == 0
        else "Assess with concordance outputs below."
    )

    rows = [
        {
            "metric": "candidate_pairs_identified",
            "value": candidate_pairs,
            "notes": "Clinically declared tumour-normal pairs recovered from TANDA HER2 MAMA T + NT.",
        },
        {
            "metric": "workbook_resolved_snp_code_pairs",
            "value": workbook_resolved_pairs,
            "notes": "Pairs whose tumour SNP code was resolved from MT-T_N using the tumour Codigo BB bloque tumoral field.",
        },
        {
            "metric": "tumour_bams_found",
            "value": tumour_bam_count,
            "notes": "Tumour BAMs resolved from breast/tumour DNA workflow folders using workbook-resolved SNP codes first.",
        },
        {
            "metric": "paired_normal_bams_found",
            "value": normal_bam_count,
            "notes": "Explicit paired-normal BAMs resolved from declared ...-NT IDs, normal block IDs, or an auditable normal SNP-code link.",
        },
        {
            "metric": "bam_confirmed_pairs",
            "value": confirmed_pairs,
            "notes": "Pairs with workbook tumour linkage plus both tumour and declared paired-normal BAM support.",
        },
        {
            "metric": "matching_logic",
            "value": "clinical paired row + MT-T_N tumour block-to-SNP mapping + explicit paired-normal BAM search",
            "notes": (
                "The harmonised master was used as a secondary consistency check only. "
                "Generic MN_* breast-normal BAMs were not treated as matched normals "
                "without an explicit crosswalk to the declared ...-NT samples."
            ),
        },
        {
            "metric": "paired_validation_supports_main_study",
            "value": support_statement,
            "notes": (
                "This validation layer is kept separate from the main unpaired cohort "
                "analyses and from the pooled healthy-control comparison."
            ),
        },
    ]
    return pd.DataFrame(rows)


def build_manual_review_table(manifest: pd.DataFrame) -> pd.DataFrame:
    review_rows: list[dict[str, object]] = []
    for _, row in manifest.iterrows():
        if row["workbook_match_status"] == "unresolved_in_workbook":
            review_rows.append(
                {
                    "case_id": row["case_id"],
                    "tumour_sample_id": row["tumour_sample_id"],
                    "normal_sample_id": row["normal_sample_id"],
                    "review_flag": "WORKBOOK_TUMOUR_SNP_UNRESOLVED",
                    "review_note": "Tumour SNP code could not be resolved from MT-T_N using the clinical tumour Codigo BB bloque tumoral field.",
                }
            )
        if not row["normal_bam_found"]:
            review_rows.append(
                {
                    "case_id": row["case_id"],
                    "tumour_sample_id": row["tumour_sample_id"],
                    "normal_sample_id": row["normal_sample_id"],
                    "review_flag": "DECLARED_PAIR_WITHOUT_MATCHED_NORMAL_BAM",
                    "review_note": (
                        "Clinical workbook declares a paired normal, but no explicit paired-normal BAM/VCF was found. "
                        "Review whether these samples exist outside the current repo or under a non-obvious lab-code crosswalk."
                    ),
                }
            )
        if not row["tumour_bam_found"]:
            review_rows.append(
                {
                    "case_id": row["case_id"],
                    "tumour_sample_id": row["tumour_sample_id"],
                    "normal_sample_id": row["normal_sample_id"],
                    "review_flag": "DECLARED_PAIR_WITHOUT_TUMOUR_BAM",
                    "review_note": "Declared paired case is present in the workbook, but the tumour BAM could not be resolved in the current breast tumour workflow folders.",
                }
            )
    if not review_rows:
        return pd.DataFrame(columns=["case_id", "tumour_sample_id", "normal_sample_id", "review_flag", "review_note"])
    return pd.DataFrame(review_rows)


def save_outputs(
    manifest: pd.DataFrame,
    workbook_resolved: pd.DataFrame,
    confirmed: pd.DataFrame,
    incomplete: pd.DataFrame,
    metadata: pd.DataFrame,
    overview: pd.DataFrame,
    concordance: pd.DataFrame,
    key_detail: pd.DataFrame,
    manual_review: pd.DataFrame,
    outdir: Path,
    out_xlsx: Path,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)

    manifest.to_csv(outdir / "candidate_pairs.tsv", sep="	", index=False)
    workbook_resolved.to_csv(outdir / "workbook_resolved_pairs.tsv", sep="	", index=False)
    confirmed.to_csv(outdir / "bam_confirmed_pairs.tsv", sep="	", index=False)
    incomplete.to_csv(outdir / "unmatched_or_incomplete_pairs.tsv", sep="	", index=False)
    metadata.to_csv(outdir / "pair_metadata_summary.tsv", sep="	", index=False)
    overview.to_csv(outdir / "paired_validation_overview.tsv", sep="	", index=False)
    concordance.to_csv(outdir / "paired_concordance_summary.tsv", sep="	", index=False)
    key_detail.to_csv(outdir / "paired_key_snp_details.tsv", sep="	", index=False)
    manual_review.to_csv(outdir / "manual_review_points.tsv", sep="	", index=False)

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        manifest.to_excel(writer, sheet_name="candidate_pairs", index=False)
        workbook_resolved.to_excel(writer, sheet_name="workbook_resolved_pairs", index=False)
        confirmed.to_excel(writer, sheet_name="bam_confirmed_pairs", index=False)
        incomplete.to_excel(writer, sheet_name="unmatched_incomplete", index=False)
        metadata.to_excel(writer, sheet_name="pair_metadata_summary", index=False)
        overview.to_excel(writer, sheet_name="overview", index=False)
        concordance.to_excel(writer, sheet_name="concordance_summary", index=False)
        key_detail.to_excel(writer, sheet_name="key_snp_details", index=False)
        manual_review.to_excel(writer, sheet_name="manual_review", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a breast paired tumour-normal validation workbook.")
    parser.add_argument("--clinical-xlsx", type=Path, default=DEFAULT_CLINICAL_XLSX)
    parser.add_argument("--nomenclature-xlsx", type=Path, default=DEFAULT_NOMENCLATURE_XLSX)
    parser.add_argument("--master-xlsx", type=Path, default=DEFAULT_MASTER_XLSX)
    parser.add_argument("--key-snp-tsv", type=Path, default=DEFAULT_KEY_SNP_TSV)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--out-xlsx", type=Path, default=DEFAULT_OUT_XLSX)
    args = parser.parse_args()

    candidate_pairs = load_candidate_pairs(args.clinical_xlsx)
    workbook_tumour_lookup = load_workbook_tumour_lookup(args.nomenclature_xlsx)
    workbook_normal_lookup = load_workbook_normal_lookup(args.nomenclature_xlsx)
    master_tumour_lookup = load_master_tumour_lookup(args.master_xlsx)
    tumour_bams = _discover_bams(TUMOUR_BAM_DIRS)
    normal_bams = _discover_bams(NORMAL_BAM_DIRS)

    manifest = build_pair_manifest(
        candidate_pairs=candidate_pairs,
        workbook_tumour_lookup=workbook_tumour_lookup,
        workbook_normal_lookup=workbook_normal_lookup,
        master_tumour_lookup=master_tumour_lookup,
        tumour_bams=tumour_bams,
        normal_bams=normal_bams,
    )
    workbook_resolved = manifest[manifest["workbook_match_status"].ne("unresolved_in_workbook")].copy()
    confirmed = manifest[
        manifest["workbook_match_status"].ne("unresolved_in_workbook")
        & manifest["tumour_bam_found"]
        & manifest["normal_bam_found"]
    ].copy()
    incomplete = manifest[
        (manifest["workbook_match_status"].eq("unresolved_in_workbook"))
        | ~(manifest["tumour_bam_found"] & manifest["normal_bam_found"])
    ].copy()

    metadata = manifest[
        [
            "case_id",
            "codigo_bb_bloque_tumoral_tumour",
            "codigo_bb_bloque_tumoral_normal",
            "tumour_sample_id",
            "normal_sample_id",
            "tumour_pathology_block_id",
            "normal_pathology_block_id",
            "tumour_snp_code",
            "normal_snp_code",
            "tumour_bam_found",
            "normal_bam_found",
            "workbook_match_status",
            "bam_match_status",
            "master_consistency_status",
            "matching_confidence",
            "matching_rule_used",
            "notes",
        ]
    ].copy()

    key_snp_metadata = load_key_snp_metadata(args.key_snp_tsv)
    concordance, key_detail = build_concordance_tables(confirmed, key_snp_metadata)
    overview = build_overview_table(manifest, workbook_resolved, confirmed)
    manual_review = build_manual_review_table(manifest)

    save_outputs(
        manifest=manifest,
        workbook_resolved=workbook_resolved,
        confirmed=confirmed,
        incomplete=incomplete,
        metadata=metadata,
        overview=overview,
        concordance=concordance,
        key_detail=key_detail,
        manual_review=manual_review,
        outdir=args.outdir,
        out_xlsx=args.out_xlsx,
    )

    print(f"Candidate pairs identified: {len(manifest)}")
    print(f"Workbook-resolved SNP-code pairs: {len(workbook_resolved)}")
    print(f"BAM-confirmed pairs: {len(confirmed)}")
    print(f"Workbook saved: {args.out_xlsx}")
    print(f"Analysis tables saved under: {args.outdir}")


if __name__ == "__main__":
    main()

