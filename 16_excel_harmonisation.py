#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
16_excel_harmonisation.py  (Option B: lossless merge + semantic harmonisation)  v4

OVERVIEW
--------
You have multiple cohorts/sheets with overlapping clinical concepts but different column names.
Goal:
  1) Keep a *lossless* raw merged master (all original columns preserved).
  2) Build an *analytical* harmonised layer by ADDING canonical columns (canon__*)
     that coalesce equivalent variables across cohorts + standardise values.

OUTPUT EXCEL (4 sheets)
-----------------------
1) raw_merged
   - Your merged master with all original columns preserved (lossless).

2) harmonised_plus_canon
   - raw_merged + canonical columns (canon__*) + provenance columns (canon__*__source)
   - plus derived columns (e.g., canon__er_status, canon__triple_negative_flag, etc.)

3) canon_audit
   - For each canon__* variable, how many rows are non-null (coverage)

4) harmonisation_map
   - Documentation: each canon__* variable and which original columns were used (priority order)

RUN (example)
---------------
python3 16_excel_harmonisation.py \
  --snp_xlsx "/path/to/Muestras SNPs nomenclaturas, equivalencias y cuantificaciones.xlsx" \
  --clinical_xlsx "/path/to/MUESTRAS + VARIABLES CLINICAS PROY SNPs GM(4).xlsx" \
  --au_xlsx "/path/to/FINAL_ECLAI_DB_Clinical_v12 febrero26 .xlsx" \
  --au_sheet "Clinical Data" \
  --out_xlsx "/path/to/MASTER_SNP_plus_clinical_HARMONISED.xlsx" \
  --manifests_dir "/path/to/manifests"

NOTE ON PERMISSION ERRORS
-------------------------
If writing to a synced cloud folder gives PermissionError, the file is usually
open in Excel or still syncing. Write locally first, then copy it back, or use
a new output filename.
"""

from __future__ import annotations

import argparse
import os
import re
import unicodedata
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from openpyxl.styles import PatternFill, Font
from openpyxl.utils import get_column_letter

from pipeline_validation import print_validation_summary, validate_file_exists, validate_nonempty


# =============================================================================
# DEFAULT CONFIG (override with CLI args)
# =============================================================================
# Public code keeps only repo-safe defaults. Your laptop can still provide
# machine-specific locations via CLI args or the TFM_DOCS_DIR environment
# variable (for example through pipeline_config.local.sh).

ROOT_DIR = Path(__file__).resolve().parent
DOCS_DIR = ROOT_DIR / "docs"


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


DEFAULT_SNP_XLSX = _discover_doc_file(
    "Muestras SNPs nomenclaturas, equivalencias y cuantificaciones.xlsx",
    DOCS_DIR / "source_workbooks",
)
DEFAULT_CLINICAL_XLSX = _discover_doc_file(
    "MUESTRAS + VARIABLES CLINICAS PROY SNPs GM(4).xlsx",
    DOCS_DIR / "source_workbooks",
)
DEFAULT_AU_XLSX = _discover_doc_file(
    "FINAL_ECLAI_DB_Clinical_v12 febrero26 .xlsx",
    DOCS_DIR / "source_workbooks",
)
DEFAULT_AU_SHEET = "Clinical Data"
DEFAULT_OUT_XLSX = ROOT_DIR / "MASTER_SNP_plus_clinical_HARMONISED.xlsx"

# Expected SNP workbook sheets we process
SNP_SHEETS = ["AT=AUs", "EN", "MT-T_N", "MN"]  # ignore OVSER


# =============================================================================
# BASIC HELPERS
# =============================================================================

def _to_str(x) -> str:
    """Safe string conversion for messy Excel cells."""
    return "" if pd.isna(x) else str(x).strip()


def _fix_mojibake_text(s: str) -> str:
    """Repair common UTF-8-as-Latin-1 mojibake without touching clean text."""
    if not isinstance(s, str):
        return s
    try:
        return s.encode("latin1").decode("utf-8")
    except UnicodeError:
        return s

def _strip_accents(s: str) -> str:
    """Remove accents after repairing any mojibake."""
    repaired = _fix_mojibake_text(s)
    return "".join(ch for ch in unicodedata.normalize("NFKD", repaired) if not unicodedata.combining(ch))



def standard_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Strip whitespace and repair header encoding drift."""
    df = df.copy()
    df.columns = [_fix_mojibake_text(str(c)).strip() for c in df.columns]
    return df


def normalize_snp_code(s: str) -> str:
    """Normalise SNP code strings by removing spaces."""
    return _to_str(s).replace(" ", "")


def parse_volume_to_ul(vol) -> Optional[float]:
    """Extract numeric part of a volume field like '20 uL' -> 20.0"""
    if pd.isna(vol):
        return None
    s = str(vol).strip().lower()
    m = re.search(r"(\d+(?:\.\d+)?)", s)
    return float(m.group(1)) if m else None


def _normalise_header_key(text: str) -> str:
    """Build an accent-insensitive lookup key for workbook headers."""
    return _strip_accents(str(text).strip()).upper()



def _pick_existing(df: pd.DataFrame, candidates: List[str]) -> str:
    """
    Choose the first existing column among candidates.
    Useful because some sheets have slightly different headers (accents/casing).
    """
    normalised = {_normalise_header_key(col): col for col in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        resolved = normalised.get(_normalise_header_key(c))
        if resolved is not None:
            return resolved
    raise KeyError(
        f"None of these columns were found: {candidates}. "
        f"Available: {df.columns.tolist()}"
    )



def _get_optional_series(df: pd.DataFrame, candidates: List[str]) -> pd.Series:
    """Return the first matching column, or an all-missing series if absent."""
    try:
        return df[_pick_existing(df, candidates)]
    except KeyError:
        return pd.Series([None] * len(df), index=df.index)


# =============================================================================
# CASE/BLOCK PARSING (supports both M####-T and C####-T)
# =============================================================================

def add_case_id_and_block_type(master: pd.DataFrame) -> pd.DataFrame:
    """
    Adds:
      - case_id: patient-level ID for case-based merges (e.g. HER2 tumour + paired normal)
      - block_type: 'T'/'NT' if sample_id looks like <Letter><digits>-T/NT

    Examples:
      M06152-T  -> case_id=M06152, block_type=T
      C14127-NT -> case_id=C14127, block_type=NT
    """
    out = master.copy()

    def parse_case_block(sid: str) -> tuple[str, Optional[str]]:
        sid = _to_str(sid)
        m = re.match(r"^([A-Z]\d+)-(NT|T)$", sid, flags=re.IGNORECASE)
        if m:
            return m.group(1).upper(), m.group(2).upper()
        return sid, None

    parsed = out["sample_id"].astype(str).apply(parse_case_block)
    out["case_id"] = parsed.apply(lambda x: x[0])
    out["block_type"] = parsed.apply(lambda x: x[1])
    return out


# =============================================================================
# AU ENDOMETRIAL ID NORMALISATION
# =============================================================================

def normalize_endo_id(x: str) -> str:
    """
    Canonical key to match AU endometrial tumour clinical (PATIENT_ID) with master sample_id.

    Input examples:
      BL039
      BL136 1:100
      EndoBL010
      MDA040
      EndoMDA_040
      EndoMDA_40

    Output:
      ENDOBL###    (3 digits)
      ENDOMDA_##   (no left padding; matches AU pattern like EndoMDA_40)
    """
    s = _to_str(x)
    if not s:
        return ""

    # Remove dilution/suffix e.g. "BL136 1:100" -> "BL136"
    s = s.split()[0].strip()
    s = s.replace(" ", "").upper()

    m = re.match(r"^ENDOBL(\d+)$", s)
    if m:
        return f"ENDOBL{int(m.group(1)):03d}"
    m = re.match(r"^BL(\d+)$", s)
    if m:
        return f"ENDOBL{int(m.group(1)):03d}"

    m = re.match(r"^ENDOMDA[_\-]?(\d+)$", s)
    if m:
        return f"ENDOMDA_{int(m.group(1))}"
    m = re.match(r"^MDA(\d+)$", s)
    if m:
        return f"ENDOMDA_{int(m.group(1))}"

    return s


# =============================================================================
# SNP SHEET PROCESSORS -> one standard schema
# =============================================================================

def process_AT(sheet_df: pd.DataFrame) -> pd.DataFrame:
    """
    AT=AUs sheet (tumour endometrial AU cohort)
    Standardises into the core schema.
    """
    df = standard_columns(sheet_df)
    out = pd.DataFrame({
        "snp_code": df["CODIGO SNPs"].map(normalize_snp_code),
        "sample_id": df["CODE"].map(_to_str),
        "nucleic_acid": df[_pick_existing(df, ["Ac. Nucleico"])].map(lambda x: _to_str(x).upper()),
        "concentration_ng_ul": pd.to_numeric(df["Qubit (ng/uL)"], errors="coerce"),
        "volume_ul": None,
        "extraction_flag": _get_optional_series(df, ["extraccion"]).map(_to_str),
        "pd_status": df.get("PD status", pd.Series([None] * len(df), index=df.index)).map(_to_str),
        "histology": df.get("Histology", pd.Series([None] * len(df), index=df.index)).map(_to_str),
        "tissue": "Endometrial",
        "tumour_normal": "Tumour",
        "sheet": "AT=AUs",
    })
    return out

def process_EN(sheet_df: pd.DataFrame) -> pd.DataFrame:
    """EN sheet (healthy endometrium)."""
    df = standard_columns(sheet_df)
    snp_col = "CODIGO SNP_EN (endometrio normal) FFPE"
    vol_col = "VOL (uL)"
    qubit_col = "Qubit (ng(uL)"  # header as in your file

    out = pd.DataFrame({
        "snp_code": df[snp_col].map(normalize_snp_code),
        "sample_id": df[_pick_existing(df, ["Codigo Noray"])].map(_to_str),
        "nucleic_acid": df[_pick_existing(df, ["Ac. Nucleico"])].map(lambda x: _to_str(x).upper()),
        "concentration_ng_ul": pd.to_numeric(df[qubit_col], errors="coerce"),
        "volume_ul": df[vol_col].map(parse_volume_to_ul),
        "extraction_flag": None,
        "pd_status": None,
        "histology": None,
        "tissue": "Endometrial",
        "tumour_normal": "Normal",
        "sheet": "EN",
    })
    return out

def process_MT_TN(sheet_df: pd.DataFrame) -> pd.DataFrame:
    """
    MT-T_N sheet: breast tumour + paired normal appear as two column blocks side-by-side.
    We split into "left" and "right" and stack.
    """
    df = standard_columns(sheet_df)

    left = pd.DataFrame({
        "sample_id": df[_pick_existing(df, ["CODIGO BB BLOQUE TUMORAL"])].map(_to_str),
        "snp_code": df[_pick_existing(df, ["CODIGO SNPs"])].map(normalize_snp_code),
        "nucleic_acid": df[_pick_existing(df, ["Ac. Nucleico"])].map(lambda x: _to_str(x).upper()),
        "volume_ul": df[_pick_existing(df, ["uL"])].map(parse_volume_to_ul),
        "concentration_ng_ul": pd.to_numeric(df[_pick_existing(df, ["EXTRACCION (ng/ul)"])], errors="coerce"),
    })

    right = pd.DataFrame({
        "sample_id": df[_pick_existing(df, ["CODIGO BB BLOQUE TUMORAL.1"])].map(_to_str),
        "snp_code": df[_pick_existing(df, ["CODIGO SNPs.1"])].map(normalize_snp_code),
        "nucleic_acid": df[_pick_existing(df, ["Ac. Nucleico.1"])].map(lambda x: _to_str(x).upper()),
        "volume_ul": df[_pick_existing(df, ["uL.1"])].map(parse_volume_to_ul),
        "concentration_ng_ul": pd.to_numeric(df[_pick_existing(df, ["EXTRACCION (ng/ul).1"])], errors="coerce"),
    })

    out = pd.concat([left, right], ignore_index=True).dropna(how="all")

    out["extraction_flag"] = None
    out["pd_status"] = None
    out["histology"] = None
    out["tissue"] = "Breast"
    out["tumour_normal"] = "Tumour"  # tumour/normal is handled later via case_id + block_type
    out["sheet"] = "MT-T_N"
    return out

def process_MN(sheet_df: pd.DataFrame) -> pd.DataFrame:
    """
    MN sheet (healthy breast): also appears as two column blocks side-by-side.
    We split into left/right and stack.
    """
    df = standard_columns(sheet_df)

    left = pd.DataFrame({
        "sample_id": df[_pick_existing(df, ["Codigo Noray-BB", "Codigo Noray"])].map(_to_str),
        "snp_code": df[_pick_existing(df, ["CODIGO SNP_MN (mama normal)"])].map(normalize_snp_code),
        "nucleic_acid": df[_pick_existing(df, ["Ac. Nucleico"])].map(lambda x: _to_str(x).upper()),
        "volume_ul": df[_pick_existing(df, ["VOL.", "VOL (uL)"])].map(parse_volume_to_ul),
        "concentration_ng_ul": pd.to_numeric(df[_pick_existing(df, ["Qubit (ng/uL)"])], errors="coerce"),
    })

    right = pd.DataFrame({
        "sample_id": df[_pick_existing(df, ["Codigo Noray-BB.1", "Codigo Noray.1"])].map(_to_str),
        "snp_code": df[_pick_existing(df, ["CODIGO SNP_MN (mama normal).1"])].map(normalize_snp_code),
        "nucleic_acid": df[_pick_existing(df, ["Ac. Nucleico.1"])].map(lambda x: _to_str(x).upper()),
        "volume_ul": df[_pick_existing(df, ["VOL..1", "VOL (uL).1"])].map(parse_volume_to_ul),
        "concentration_ng_ul": pd.to_numeric(df[_pick_existing(df, ["Qubit (ng/uL).1"])], errors="coerce"),
    })

    out = pd.concat([left, right], ignore_index=True).dropna(how="all")

    out["extraction_flag"] = None
    out["pd_status"] = None
    out["histology"] = None
    out["tissue"] = "Breast"
    out["tumour_normal"] = "Normal"
    out["sheet"] = "MN"
    return out

def build_master(snp_xlsx: Path) -> pd.DataFrame:
    """
    Reads the SNP workbook and constructs a single master table with a standard schema.
    """
    xls = pd.ExcelFile(snp_xlsx)
    frames: List[pd.DataFrame] = []

    if "AT=AUs" in xls.sheet_names:
        frames.append(process_AT(pd.read_excel(xls, sheet_name="AT=AUs")))
    if "EN" in xls.sheet_names:
        frames.append(process_EN(pd.read_excel(xls, sheet_name="EN")))
    if "MT-T_N" in xls.sheet_names:
        frames.append(process_MT_TN(pd.read_excel(xls, sheet_name="MT-T_N")))
    if "MN" in xls.sheet_names:
        frames.append(process_MN(pd.read_excel(xls, sheet_name="MN")))

    if not frames:
        raise ValueError(f"No expected SNP sheets found. Available: {xls.sheet_names}")

    frames = [frame.dropna(axis=1, how="all") for frame in frames if not frame.empty and not frame.dropna(how="all").empty]
    if not frames:
        raise ValueError("All processed SNP sheets were empty after dropping all-NA rows.")

    master = pd.concat(frames, ignore_index=True)

    # Filter: keep only DNA/RNA rows
    master["nucleic_acid"] = master["nucleic_acid"].astype(str).str.upper().replace({"": pd.NA, "NAN": pd.NA})
    master = master[master["nucleic_acid"].isin(["DNA", "RNA"])]

    # Filter: require SNP code
    master["snp_code"] = master["snp_code"].astype(str).replace({"": pd.NA, "nan": pd.NA})
    master = master[master["snp_code"].notna()].reset_index(drop=True)

    # Add case-level identifiers for HER2 case-based merges
    master = add_case_id_and_block_type(master)

    cols = [
        "snp_code", "nucleic_acid",
        "tissue", "tumour_normal",
        "sample_id", "case_id", "block_type",
        "concentration_ng_ul", "volume_ul",
        "extraction_flag", "pd_status", "histology",
        "sheet",
    ]
    return master[cols].reset_index(drop=True)


# =============================================================================
# CLINICAL LOADER + MERGER (lossless, prefixed)
# =============================================================================

def load_clinical_tables(clinical_xlsx: Path) -> Dict[str, pd.DataFrame]:
    """
    Loads and prefixes each clinical table so we never have column name collisions.

    Prefix convention:
      clin_mn__*      = healthy breast (TANDA MAMA SANA)
      clin_en__*      = healthy endometrium (TANDA ENDOMETRIO SANO)
      clin_her2__*    = HER2 case-level table (TANDA HER2 MAMA T + NT)
      clin_dcs__*     = DCs table (DCs MAMA HER2) merged into HER2 by NHC
    """
    xls = pd.ExcelFile(clinical_xlsx)
    tables: Dict[str, pd.DataFrame] = {}

    # Healthy breast (MN): join by sample_id
    if "TANDA MAMA SANA" in xls.sheet_names:
        df = standard_columns(pd.read_excel(xls, sheet_name="TANDA MAMA SANA"))
        sample_col = _pick_existing(df, ["Codigo Noray-BB"])
        df = df.rename(columns={sample_col: "sample_id"})
        df["sample_id"] = df["sample_id"].map(_to_str)
        tables["clinical_breast_healthy"] = df.rename(columns={c: f"clin_mn__{c}" for c in df.columns if c != "sample_id"})

    # Healthy endometrium (EN): join by sample_id
    if "TANDA ENDOMETRIO SANO" in xls.sheet_names:
        df = standard_columns(pd.read_excel(xls, sheet_name="TANDA ENDOMETRIO SANO"))
        sample_col = _pick_existing(df, ["Codigo Noray"])
        df = df.rename(columns={sample_col: "sample_id"})
        df["sample_id"] = df["sample_id"].map(_to_str)
        tables["clinical_endo_healthy"] = df.rename(columns={c: f"clin_en__{c}" for c in df.columns if c != "sample_id"})

    # HER2 case-level clinical: join by case_id
    her2 = None
    if "TANDA HER2 MAMA T + NT" in xls.sheet_names:
        her2 = standard_columns(pd.read_excel(xls, sheet_name="TANDA HER2 MAMA T + NT"))
        her2 = her2.rename(columns={
            _pick_existing(her2, ["CODIGO BB CASO"]): "case_id",
            _pick_existing(her2, ["CODIGO BB BLOQUE TUMORAL"]): "tumour_block_id",
            _pick_existing(her2, ["CODIGO BB BLOQUE NORMAL"]): "paired_normal_block_id",
        })
        her2["case_id"] = her2["case_id"].map(_to_str).str.upper()
        her2["NHC"] = her2["NHC"].map(_to_str)
        her2 = her2.rename(columns={c: f"clin_her2__{c}" for c in her2.columns if c not in {"case_id", "NHC"}})

    # Enrich HER2 with DCs via NHC
    if her2 is not None and "DCs MAMA HER2" in xls.sheet_names:
        dcs = standard_columns(pd.read_excel(xls, sheet_name="DCs MAMA HER2"))
        dcs["NHC"] = dcs["NHC"].map(_to_str)
        dcs = dcs.rename(columns={c: f"clin_dcs__{c}" for c in dcs.columns if c != "NHC"})
        her2 = her2.merge(dcs, how="left", on="NHC")

    if her2 is not None:
        tables["clinical_her2_case"] = her2

    return tables

def load_au_endo_clinical(au_xlsx: Path, sheet_name: str) -> Optional[pd.DataFrame]:
    """
    Loads the main endometrial clinical sheet and prefixes with clin_au_endo__.
    Join key is 'endo_id', derived from workbook PATIENT_ID (and from master sample_id).
    """
    if not au_xlsx.exists():
        print(f"AU endo clinical file not found at: {au_xlsx} (skipping).")
        return None

    au = pd.read_excel(au_xlsx, sheet_name=sheet_name)
    au = standard_columns(au)

    if "PATIENT_ID" not in au.columns:
        print(f"AU endo clinical: PATIENT_ID column not found in {sheet_name}; skipping merge.")
        return None

    au = au.copy()
    au["endo_id"] = au["PATIENT_ID"].astype(str).map(normalize_endo_id)
    au = au[au["endo_id"].astype(str).str.len() > 0].copy()

    # If duplicates exist, keep first (you can revise if you prefer latest/complete)
    if au["endo_id"].duplicated().any():
        print("[WARN] Duplicate endometrial PATIENT_ID values found after normalisation; keeping first row per endo_id.")
        au = au.drop_duplicates(subset=["endo_id"], keep="first").copy()

    au = au.rename(columns={c: f"clin_au_endo__{c}" for c in au.columns if c != "endo_id"})
    cols = ["endo_id"] + [c for c in au.columns if c != "endo_id"]
    return au[cols].copy()


def merge_master_with_clinical(
    master_df: pd.DataFrame,
    clinical_xlsx: Path,
    au_endo_xlsx: Path,
    au_sheet: str
) -> pd.DataFrame:
    """
    Lossless merge:
      - MN and EN join on sample_id
      - HER2 joins on case_id (patient-level)
      - AU endo joins on endo_id (normalised key)
    """
    tables = load_clinical_tables(clinical_xlsx)
    merged = master_df.copy()

    if "clinical_breast_healthy" in tables:
        merged = merged.merge(tables["clinical_breast_healthy"], how="left", on="sample_id")

    if "clinical_endo_healthy" in tables:
        merged = merged.merge(tables["clinical_endo_healthy"], how="left", on="sample_id")

    if "clinical_her2_case" in tables:
        merged = merged.merge(tables["clinical_her2_case"], how="left", on="case_id")

    au = load_au_endo_clinical(au_endo_xlsx, sheet_name=au_sheet)
    if au is not None:
        merged = merged.copy()
        merged["endo_id"] = merged["sample_id"].astype(str).map(normalize_endo_id)
        candidate_counts = au["endo_id"].value_counts().to_dict()
        merged["merge__au_endo_candidate_count"] = merged["endo_id"].map(candidate_counts).fillna(0).astype(int)
        merged = merged.merge(au, how="left", on="endo_id")
        merged["merge__au_endo_match_method"] = pd.NA
        merged["merge__au_endo_match_status"] = pd.NA
        endo_mask = merged["sheet"].eq("AT=AUs")
        patient_col = "clin_au_endo__PATIENT_ID"
        if patient_col in merged.columns:
            matched_mask = endo_mask & merged[patient_col].notna()
            ambiguous_mask = endo_mask & merged["merge__au_endo_candidate_count"].gt(1)
            merged.loc[matched_mask, "merge__au_endo_match_method"] = "sample_id->PATIENT_ID(normalized)"
            merged.loc[matched_mask, "merge__au_endo_match_status"] = "matched"
            merged.loc[endo_mask & ~matched_mask, "merge__au_endo_match_status"] = "unmatched"
            merged.loc[ambiguous_mask, "merge__au_endo_match_status"] = "ambiguous"

    return merged


def build_clinical_merge_audits(merged_df: pd.DataFrame, au_endo_xlsx: Path, au_sheet: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build explicit audit tables for the AU endometrial clinical merge."""
    endo_rows = merged_df.loc[merged_df["sheet"].eq("AT=AUs")].copy() if "sheet" in merged_df.columns else merged_df.iloc[0:0].copy()
    for col in ["merge__au_endo_match_status", "merge__au_endo_match_method", "merge__au_endo_candidate_count", "endo_id"]:
        if col not in endo_rows.columns:
            endo_rows[col] = pd.NA

    matched_ids = set(endo_rows.loc[endo_rows["merge__au_endo_match_status"].eq("matched"), "endo_id"].dropna().astype(str))
    unmatched_ids = set(endo_rows.loc[endo_rows["merge__au_endo_match_status"].eq("unmatched"), "endo_id"].dropna().astype(str))
    ambiguous_ids = set(endo_rows.loc[endo_rows["merge__au_endo_match_status"].eq("ambiguous"), "endo_id"].dropna().astype(str))

    au = load_au_endo_clinical(au_endo_xlsx, sheet_name=au_sheet)
    if au is None:
        clinical_only = pd.DataFrame(columns=["endo_id", "clin_au_endo__PATIENT_ID", "audit_note"])
        unique_au_ids = set()
    else:
        unique_au_ids = set(au["endo_id"].dropna().astype(str))
        clinical_only = au.loc[~au["endo_id"].astype(str).isin(set(endo_rows["endo_id"].dropna().astype(str)))].copy()
        if not clinical_only.empty:
            clinical_only["audit_note"] = "Present in AU clinical workbook but absent from the SNP master."

    summary = pd.DataFrame([
        {"Metric": "AT_AUs_rows_in_master", "Value": int(len(endo_rows)), "Note": "Rows from the SNP master sheet AT=AUs before harmonised coalescing."},
        {"Metric": "AT_AUs_unique_normalised_ids", "Value": int(endo_rows["endo_id"].dropna().astype(str).nunique()), "Note": "Unique normalised endometrial identifiers in the SNP master."},
        {"Metric": "AU_clinical_unique_normalised_ids", "Value": int(len(unique_au_ids)), "Note": "Unique normalised PATIENT_ID values in the ECLAI Clinical Data sheet."},
        {"Metric": "Matched_AT_AUs_rows", "Value": int(endo_rows["merge__au_endo_match_status"].eq("matched").sum()), "Note": "Rows matched by sample_id -> PATIENT_ID after normalisation."},
        {"Metric": "Matched_AT_AUs_unique_ids", "Value": int(len(matched_ids)), "Note": "Unique endometrial IDs matched to AU clinical data."},
        {"Metric": "Unmatched_AT_AUs_rows", "Value": int(endo_rows["merge__au_endo_match_status"].eq("unmatched").sum()), "Note": "Rows in the SNP master with no AU clinical match."},
        {"Metric": "Unmatched_AT_AUs_unique_ids", "Value": int(len(unmatched_ids)), "Note": "Unique endometrial IDs in the SNP master with no AU clinical match."},
        {"Metric": "Ambiguous_AT_AUs_unique_ids", "Value": int(len(ambiguous_ids)), "Note": "Unique IDs with more than one AU clinical candidate after normalisation."},
        {"Metric": "AU_clinical_only_unique_ids", "Value": int(len(unique_au_ids - set(endo_rows["endo_id"].dropna().astype(str)))), "Note": "Clinical IDs present in the AU workbook but absent from the SNP master."},
    ])

    match_detail_cols = [
        c for c in [
            "snp_code", "sample_id", "endo_id", "nucleic_acid", "tumour_normal", "tissue",
            "merge__au_endo_match_status", "merge__au_endo_match_method", "merge__au_endo_candidate_count",
            "clin_au_endo__PATIENT_ID", "clin_au_endo__STUDY", "clin_au_endo__COUNTRY OF SAMPLE ORIGIN",
            "clin_au_endo__HISTOLOGY_GROUP", "clin_au_endo__FIGO_STAGE", "clin_au_endo__PD_STATUS",
            "clin_au_endo__EXITUS", "clin_au_endo__OS", "clin_au_endo__PFS"
        ] if c in endo_rows.columns
    ]
    match_detail = endo_rows[match_detail_cols].sort_values(["merge__au_endo_match_status", "endo_id", "sample_id"], kind="stable") if match_detail_cols else endo_rows
    unmatched = match_detail.loc[match_detail["merge__au_endo_match_status"].isin(["unmatched", "ambiguous"])].copy() if not match_detail.empty else match_detail

    if clinical_only.empty:
        clinical_only = pd.DataFrame(columns=["endo_id", "clin_au_endo__PATIENT_ID", "audit_note"])

    return summary, match_detail, unmatched, clinical_only


# =============================================================================
# SEMANTIC NORMALISATION HELPERS (Option B)
# =============================================================================

def _norm_text(x) -> str:
    """Normalise text (remove NBSP, trim)."""
    s = _fix_mojibake_text(_to_str(x))
    return s.replace("\u00a0", " ").strip()


def has_meaningful_value(x) -> bool:
    """Treat blank strings and placeholder text as missing during coalescing."""
    if pd.isna(x):
        return False
    if isinstance(x, str) and _norm_text(x).upper() in {"", "NAN", "NONE", "NULL", "N/A"}:
        return False
    return True


def parse_yes_no(x) -> Optional[str]:
    """
    Standardise common yes/no encodings -> 'Yes'/'No'/None.
    Accepts: SI/S?/No, Yes/No, Y/N, 1/0, True/False, Pos/Neg.
    """
    if pd.isna(x):
        return None
    s = _strip_accents(_norm_text(x).lower())
    if s == "":
        return None

    # numeric
    if re.fullmatch(r"-?\d+(\.\d+)?", s):
        try:
            v = float(s)
            if v == 1:
                return "Yes"
            if v == 0:
                return "No"
        except Exception:
            pass

    # bool words
    if s in {"true", "t"}:
        return "Yes"
    if s in {"false", "f"}:
        return "No"

    # Spanish
    if s in {"si", "s"}:
        return "Yes"
    if s in {"no", "n"}:
        return "No"

    # English
    if s in {"yes", "y"}:
        return "Yes"

    # Clinical shorthand
    if s in {"positivo", "positiva", "pos", "+", "positive"}:
        return "Yes"
    if s in {"negativo", "negativa", "neg", "-", "negative"}:
        return "No"

    return None


def _normalise_treatment_text(x) -> str:
    """Normalise treatment-history text so simple keyword rules are stable."""
    if pd.isna(x):
        return ""
    s = _strip_accents(_norm_text(x).lower())
    s = re.sub(r"\s+", " ", s)
    return s


def derive_breast_treatment_exposure_flags(x) -> Dict[str, Optional[str]]:
    """
    Convert the free-text breast `Tto` field into coarse exposure flags.

    These flags are intentionally conservative and are meant for exploratory
    treatment-stratified analyses, not formal treatment-response modelling.
    """
    s = _normalise_treatment_text(x)
    if s in {"", "x", "na", "n/a", "none", "no consta"}:
        return {
            "canon__treatment_chemotherapy": None,
            "canon__treatment_anti_her2": None,
            "canon__treatment_endocrine": None,
            "canon__treatment_radiotherapy": None,
        }

    patterns = {
        "canon__treatment_chemotherapy": r"\bac\b|\bfec\b|\bcmf\b|taxol|paclitaxel|docetaxel|cbdca|carbo|quimio|chemo|\bqt\b",
        "canon__treatment_anti_her2": r"hercept|trastu|lapat|pertu|anti[\s-]?her[\s-]?2",
        "canon__treatment_endocrine": r"tamox|letroz|exemest|anastro|arimid|fulves|hormon|terapia hormonal|\bht\b",
        "canon__treatment_radiotherapy": r"\brt\b|radiot|radioter|rte|rdt",
    }
    return {
        key: ("Yes" if re.search(pattern, s) else "No")
        for key, pattern in patterns.items()
    }


def extract_first_number(x) -> Optional[float]:
    """
    Extract first numeric token from a messy cell.
    Examples:
      "12 aÃ±os" -> 12
      "13a" -> 13
      "~11" -> 11
      "11-12" -> 11
      ">=6" -> 6
      "6,2" -> 6.2
    """
    if pd.isna(x):
        return None
    s = _norm_text(x)
    if s == "":
        return None
    s = s.replace(",", ".")
    m = re.search(r"(\d+(?:\.\d+)?)", s)
    return float(m.group(1)) if m else None


def parse_bmi_value(x) -> Optional[float]:
    """
    Parse BMI while rejecting weight-only strings and free-text placeholders.
    """
    if pd.isna(x):
        return None
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        v = float(x)
        return round(v, 2) if 10.0 <= v <= 70.0 else None

    s = _norm_text(x).upper().replace(",", ".").replace("?", "-")
    if s in {"", "X", "SOBREPESO", "OBESIDAD", "SOBREPESO/OBESIDAD"}:
        return None

    m_weight = re.search(r"(\d+(?:\.\d+)?)\s*KG\b", s)
    m_height = re.search(r"(\d+(?:\.\d+)?)\s*CM\b", s)
    if m_weight and m_height:
        weight = float(m_weight.group(1))
        height_cm = float(m_height.group(1))
        if 25.0 <= weight <= 250.0 and 120.0 <= height_cm <= 230.0:
            bmi = weight / ((height_cm / 100.0) ** 2)
            return round(bmi, 2) if 10.0 <= bmi <= 70.0 else None

    if "KG" in s and "CM" not in s:
        return None

    m = re.search(r"(\d+(?:\.\d+)?)", s)
    if not m:
        return None
    v = float(m.group(1))
    return round(v, 2) if 10.0 <= v <= 70.0 else None


def parse_bmi_category(x) -> Optional[str]:
    """
    Parse a BMI value or a source label into a standard BMI category.

    This preserves meaningful source text such as "Sobrepeso" instead of
    discarding it as missing.
    """
    if pd.isna(x):
        return None

    if isinstance(x, str):
        s = _norm_text(x).upper().replace(",", ".").replace("?", "-").strip()
        if not s or s in {"X", "NA", "N/A", "NR", "NULL"}:
            return None
        if "SOBREPESO" in s or "OVERWEIGHT" in s:
            return "Overweight"
        if "OBESIDAD" in s or "OBESE" in s:
            return "Obese"
        if "BAJO PESO" in s or "UNDERWEIGHT" in s:
            return "Underweight"

    bmi = parse_bmi_value(x)
    if bmi is None:
        return None
    if bmi < 18.5:
        return "Underweight"
    if bmi < 25.0:
        return "Normal"
    if bmi < 30.0:
        return "Overweight"
    return "Obese"


def parse_weight_kg_value(x) -> Optional[float]:
    """Parse a weight value in kilograms from messy BMI/anthropometric text."""
    if pd.isna(x):
        return None
    s = _norm_text(x).upper().replace(",", ".").replace("?", "-")
    if not s or s in {"X", "SOBREPESO", "OBESIDAD", "SOBREPESO/OBESIDAD"}:
        return None
    m = re.search(r"(\d+(?:\.\d+)?)\s*KG\b", s)
    if not m:
        return None
    v = float(m.group(1))
    return round(v, 2) if 25.0 <= v <= 250.0 else None


def parse_height_cm_value(x) -> Optional[float]:
    """Parse a height value in centimetres from messy BMI/anthropometric text."""
    if pd.isna(x):
        return None
    s = _norm_text(x).upper().replace(",", ".").replace("?", "-")
    if not s or s in {"X", "SOBREPESO", "OBESIDAD", "SOBREPESO/OBESIDAD"}:
        return None
    m = re.search(r"(\d+(?:\.\d+)?)\s*CM\b", s)
    if not m:
        return None
    v = float(m.group(1))
    return round(v, 2) if 120.0 <= v <= 230.0 else None


def parse_excel_serial_date_percent(x) -> Optional[float]:
    """
    Recover integers that Excel displayed as 1899/1900 dates because the cell
    format was changed from numeric to date.
    """
    if pd.isna(x):
        return None

    if isinstance(x, pd.Timestamp):
        ts = x
    else:
        sx = _norm_text(x)
        if not (
            re.fullmatch(r"\d{4}-\d{2}-\d{2}(?:[ T]\d{2}:\d{2}:\d{2})?", sx)
            or re.fullmatch(r"\d{1,2}/\d{1,2}/\d{4}(?: \d{2}:\d{2}:\d{2})?", sx)
        ):
            return None
        ts = pd.to_datetime(sx, errors="coerce", dayfirst=False)
        if pd.isna(ts):
            return None

    if ts.year not in {1899, 1900}:
        return None

    serial = float((ts.normalize() - pd.Timestamp("1899-12-30")).days)
    return serial if 0.0 <= serial <= 100.0 else None


def parse_numeric_midpoint(x) -> Optional[float]:
    """Parse a single numeric value or numeric range midpoint from messy text."""
    if pd.isna(x):
        return None
    s = _norm_text(x).upper().replace(",", ".").replace("?", "-")
    if s == "":
        return None

    serial = parse_excel_serial_date_percent(x)
    if serial is not None:
        return serial

    m_range = re.search(r"(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)", s)
    if m_range:
        return (float(m_range.group(1)) + float(m_range.group(2))) / 2.0

    m_ineq = re.search(r"[<>]=?\s*(\d+(?:\.\d+)?)", s)
    if m_ineq:
        return float(m_ineq.group(1))

    m_num = re.search(r"(\d+(?:\.\d+)?)", s)
    if m_num:
        return float(m_num.group(1))
    return None


def parse_ki67_percent(x, source: Optional[str] = None) -> Optional[float]:
    """Normalise KI67 to percentage scale 0-100 while handling date-misformatted cells."""
    if pd.isna(x):
        return None
    s = _norm_text(x).upper()
    if s in {"", "X", "N/A", "ND", "NOT DONE", "PENDING"}:
        return None

    value = parse_numeric_midpoint(x)
    if value is None:
        return None

    source = _norm_text(source)
    if source == "clin_dcs__KI67":
        pct = value * 100.0 if value <= 1.0 else value
    elif source == "clin_au_endo__FFPE_KI67":
        pct = value
    else:
        pct = value * 100.0 if value < 1.0 else value

    return round(pct, 1) if 0.0 <= pct <= 100.0 else None


def parse_receptor_value(x) -> Tuple[Optional[float], Optional[str]]:
    """Parse mixed ER/PR inputs into a numeric percentage when available plus status."""
    if pd.isna(x):
        return None, None
    s = _norm_text(x).upper()
    if s in {"", "N/A", "ND", "NOT DONE"}:
        return None, None

    if any(token in s for token in ("NEGATIVE", "NEGATIVO", "NEG")) and "POS" not in s:
        return 0.0, "Negative"
    if any(token in s for token in ("POSITIVE", "POSITIVO", "FOCAL", "WEAK POSITIVE", "LOW POSITIVE")):
        return parse_numeric_midpoint(x), "Positive"

    value = parse_numeric_midpoint(x)
    if value is None:
        return None, None
    if 0.0 <= value <= 100.0:
        return value, ("Positive" if value >= 1.0 else "Negative")
    return None, None


def parse_p53_ihc_value(x) -> Tuple[Optional[float], Optional[str]]:
    """Parse p53 IHC values into a numeric score when possible plus Normal/Aberrant status."""
    if pd.isna(x):
        return None, None
    s = _norm_text(x).upper()
    if s in {"", "N/A", "ND", "NOT DONE"}:
        return None, None

    if any(token in s for token in ("WT", "WILD TYPE", "WILD-TYPE", "NORMAL")):
        return parse_numeric_midpoint(x), "Normal"

    aberrant_tokens = (
        "POSITIVE", "OVEREXPRESS", "OVER-EXPRESS", "HIGH", "ABERRANT",
        "MUT", "LOST", "NULL", "ABSENT"
    )
    if any(token in s for token in aberrant_tokens):
        return parse_numeric_midpoint(x), "Aberrant"

    value = parse_numeric_midpoint(x)
    if value is None:
        return None, None
    cutoff = 10.0 if value > 1.0 else 0.10
    return value, ("Aberrant" if value >= cutoff else "Normal")


def parse_grade(x) -> Optional[int]:
    """Extract grade as integer 1-4 from strings like 'G3', 'GRADO 2', '2'."""
    if pd.isna(x):
        return None
    s = _norm_text(x).upper()
    if s == "":
        return None
    m = re.search(r"\b([1-4])\b", s)
    if m:
        return int(m.group(1))
    m = re.search(r"\bG([1-4])\b", s)
    if m:
        return int(m.group(1))
    return None


def parse_figo_stage(x) -> Optional[str]:
    """
    Normalise FIGO stage into a consistent compact code:
      I, IA, IB, II, IIIA, IIIC1, IVB, etc.
    """
    if pd.isna(x):
        return None
    s = _norm_text(x).upper()
    if s == "":
        return None

    s = re.sub(r"FIGO", "", s).strip()
    s = s.replace("STAGE", "").replace("ESTADIO", "").strip()
    s = s.replace(" ", "").replace("-", "")

    arabic_map = {"1": "I", "2": "II", "3": "III", "4": "IV"}
    m = re.match(r"^([1-4])([A-C]?[0-9]?)$", s)
    if m:
        return arabic_map[m.group(1)] + m.group(2)

    m = re.match(r"^(I{1,3}|IV)([A-C]?[0-9]?)$", s)
    if m:
        return m.group(1) + m.group(2)

    m = re.match(r"^(I{1,3}|IV)(A|B|C|C1|C2)$", s)
    if m:
        return s

    return None


def parse_date(x) -> Optional[pd.Timestamp]:
    """Robust date parser (Excel serials + strings)."""
    if pd.isna(x):
        return None

    if isinstance(x, pd.Timestamp):
        return x

    if isinstance(x, (int, float)) and not isinstance(x, bool):
        try:
            if 1 <= float(x) <= 60000:
                dt = pd.to_datetime(float(x), unit="D", origin="1899-12-30", errors="coerce")
                return None if pd.isna(dt) else dt
        except Exception:
            return None

    s = _norm_text(x)
    if not s:
        return None

    iso_like = re.fullmatch(r"\d{4}-\d{1,2}-\d{1,2}(?:[ T]\d{1,2}:\d{2}(?::\d{2})?)?", s)
    if iso_like:
        dt = pd.to_datetime(s, errors="coerce")
        return None if pd.isna(dt) else dt

    slash_like = re.fullmatch(r"(\d{1,2})/(\d{1,2})/(\d{4})(?: (\d{1,2}:\d{2}(?::\d{2})?))?", s)
    if slash_like:
        first = int(slash_like.group(1))
        second = int(slash_like.group(2))
        time_part = slash_like.group(4) or ""
        if first > 12 and second <= 12:
            fmt = "%d/%m/%Y" + (" %H:%M:%S" if time_part.count(":") == 2 else " %H:%M" if time_part else "")
        elif second > 12 and first <= 12:
            fmt = "%m/%d/%Y" + (" %H:%M:%S" if time_part.count(":") == 2 else " %H:%M" if time_part else "")
        else:
            fmt = "%d/%m/%Y" + (" %H:%M:%S" if time_part.count(":") == 2 else " %H:%M" if time_part else "")
        dt = pd.to_datetime(s, format=fmt, errors="coerce")
        return None if pd.isna(dt) else dt

    dash_like = re.fullmatch(r"(\d{1,2})-(\d{1,2})-(\d{4})(?: (\d{1,2}:\d{2}(?::\d{2})?))?", s)
    if dash_like:
        time_part = dash_like.group(4) or ""
        fmt = "%d-%m-%Y" + (" %H:%M:%S" if time_part.count(":") == 2 else " %H:%M" if time_part else "")
        dt = pd.to_datetime(s, format=fmt, errors="coerce")
        return None if pd.isna(dt) else dt

    try:
        dt = pd.to_datetime(s, errors="coerce")
        return None if pd.isna(dt) else dt
    except Exception:
        return None


def coalesce_with_source(df: pd.DataFrame, cols: List[str]) -> Tuple[pd.Series, pd.Series]:
    """
    Coalesce: first non-null across candidate columns, returning both value and which source column was used.
    """
    existing: List[str] = []
    seen = set()
    for candidate in cols:
        try:
            resolved = _pick_existing(df, [candidate])
        except KeyError:
            continue
        if resolved not in seen:
            existing.append(resolved)
            seen.add(resolved)

    out = pd.Series([pd.NA] * len(df), index=df.index)
    src = pd.Series([pd.NA] * len(df), index=df.index)

    if not existing:
        return out, src

    for c in existing:
        mask = out.isna() & df[c].map(has_meaningful_value)
        out.loc[mask] = df.loc[mask, c]
        src.loc[mask] = c

    return out, src


# =============================================================================
# HARMONISATION MAP (your real headers)
# =============================================================================

def build_harmonisation_map() -> Dict[str, List[str]]:
    """
    Canonical variable -> list of candidate source columns (in priority order).
    The first non-null wins.

    v4 additions
    ------------
    - canon__bmi          : added clin_dcs__BMI (breast tumour; was missing)
    - canon__grade        : added clin_dcs__GRADO (breast tumour; was missing)
    - canon__ki67         : new â€” coalesces FFPE_KI67 (endo) + clin_dcs__KI67 (breast)
    - canon__p53_ihc      : new â€” coalesces TP53 IHC (endo) + p53 IHC (breast)
    - canon__date_diagnosis: new â€” diagnosis date (breast)
    - canon__date_last_fu  : new â€” last follow-up date (breast), used to derive OS
    - canon__date_recurrence: new â€” recurrence date (breast)
    - canon__recurrence    : new â€” recurrence binary flag (breast + endo)
    - canon__distant_mets  : new â€” distant metastasis flag (breast)
    - canon__lymph_nodes   : new â€” lymph node ratio string (breast) / N field (endo)
    - canon__myometrial_invasion: new â€” myometrial infiltration (endo)
    - canon__lvsi          : new â€” lymphovascular space invasion (endo)
    - canon__risk_group    : new â€” risk-of-recurrence group (endo)
    - canon__treatment     : new â€” treatment description (breast)
    - canon__vital_status  : new â€” verbose vital status (breast)
    - canon__ptnm          : already present but now also picks up clin_dcs__pTNM
    """
    return {
        # â”€â”€ Demographics â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        "canon__age": [
            "clin_au_endo__AGE",
            "clin_dcs__Edad dx",
            "clin_mn__Edad muestra",
            "clin_en__Edad muestra",
        ],
        "canon__bmi": [
            "clin_au_endo__BMI_CALCULATED",
            "clin_au_endo__BMI",
            "clin_dcs__BMI",        # â† v4: breast tumour BMI (was omitted before)
            "clin_her2__BMI",
            "clin_mn__BMI",
            "clin_en__BMI",
        ],
        "canon__weight_kg": [
            "clin_au_endo__BMI_CALCULATED",
            "clin_au_endo__BMI",
            "clin_dcs__BMI",
            "clin_her2__BMI",
            "clin_mn__BMI",
            "clin_en__BMI",
        ],
        "canon__height_cm": [
            "clin_au_endo__BMI_CALCULATED",
            "clin_au_endo__BMI",
            "clin_dcs__BMI",
            "clin_her2__BMI",
            "clin_mn__BMI",
            "clin_en__BMI",
        ],
        "canon__date_birth": [
            "clin_mn__Fecha nacimiento",
            "clin_en__Fecha nacimiento",
            "clin_dcs__Fecha nacimiento",
        ],
        "canon__menarche_age": [
            "clin_mn__Menarquia",
            "clin_en__Menarquia",
            "clin_her2__Menarquia",
        ],
        "canon__menopause_age_or_status": [
            "clin_mn__Menopausia",
            "clin_en__Menopausia",
            "clin_her2__Menopausia",
        ],
        "canon__pregnancy_history": [
            "clin_mn__AGO (GESTACIONES, PARTOS, CESÃREAS, ABORTOS)",
            "clin_en__AGO (GESTACIONES, PARTOS, CESÃREAS, ABORTOS)",
            "clin_her2__AGO (GESTACIONES, PARTOS, CESÃREAS, ABORTOS)",
        ],

        # â”€â”€ Diagnosis / histology / tumour â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        "canon__diagnosis": [
            "clin_dcs__Dx",             # â† v4: prefer DCS Dx (cleaner, more complete)
            "clin_her2__DX",
            "clin_au_endo__DIAGNOSIS",
        ],
        "canon__histology": [
            "clin_au_endo__HISTOLOGY_GROUP",
            "clin_au_endo__HISTOLOGY",
            "histology",
        ],
        "canon__grade": [
            "clin_au_endo__GRADE",
            "clin_dcs__GRADO",          # â† v4: breast tumour grade (was omitted before)
        ],
        "canon__figo_stage": [
            "clin_au_endo__FIGO_STAGE",
            "clin_au_endo__FIGO_STAGE_GROUP",
        ],
        "canon__ptnm": [
            "clin_dcs__pTNM",           # â† v4: breast pTNM
            "clin_au_endo__T",          #   (endo T/N/M kept in separate cols below)
        ],

        # â”€â”€ Pathological staging extras (endometrial) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        "canon__myometrial_invasion": [  # â† v4: new
            "clin_au_endo__MYOMETRIAL_INFILTRATION",
        ],
        "canon__lvsi": [                 # â† v4: new
            "clin_au_endo__LVSI",
        ],
        "canon__risk_group": [           # â† v4: new
            "clin_au_endo__RISK_OF_RECURRENCE",
        ],
        "canon__lymph_nodes": [          # â† v4: new â€” node ratio string (breast) or N (endo)
            "clin_dcs__Numero_ganglio",
            "clin_au_endo__N",
        ],

        # â”€â”€ Biomarkers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        "canon__ki67": [                 # â† v4: new unified KI67
            "clin_au_endo__FFPE_KI67",
            "clin_dcs__KI67",
        ],
        "canon__p53_ihc": [              # â† v4: new unified p53/TP53 IHC
            "clin_au_endo__FFPE_TP53_IHC_COMBINED",
            "clin_au_endo__FFPE_TP53_IHC_DYCOTOMIC",
            "clin_au_endo__FFPE_TP53_IHC",
            "clin_au_endo__FFPE_TP53_IHC POLAND RESULT",
            "clin_dcs__p53",
        ],

        # â”€â”€ ER / PR status â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        "canon__er_raw": [
            "clin_au_endo__FFPE_ESR1_RECEPTORS_GROUPED",
            "clin_au_endo__FFPE_ER1_RECEPTORS_RAW_VALUE",
            "clin_au_endo__FFPE_ER1_RECEPTORS_RAW_VALUE POLAND RESULTS",
        ],
        "canon__pr_raw": [
            "clin_au_endo__FFPE_PR_RECEPTORS_GROUPED",
            "clin_au_endo__FFPE_PR_RECEPTORS_RAW_VALUE",
            "clin_au_endo__FFPE_PR_RECEPTORS_RAW_VALUE POLAND RESULTS",
        ],
        "canon__er_status_breast": [
            "clin_dcs__RE",
        ],
        "canon__pr_status_breast": [
            "clin_dcs__RP",
        ],

        # â”€â”€ HER2 â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        "canon__her2_copies": [
            "clin_her2__COPIAS HER2",
            "clin_dcs__COPIAS HER2",
        ],
        "canon__her2_copies_note": [
            "clin_her2__COPIAS HER2",
            "clin_dcs__COPIAS HER2",
        ],

        # â”€â”€ Survival / outcomes â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        "canon__os_months": [
            "clin_au_endo__OS",
            # breast OS is DERIVED from dates in apply_semantic_transforms; not coalesced here
        ],
        "canon__pfs_months": [
            "clin_au_endo__PFS",
            # breast PFS not available in source
        ],
        "canon__pd_status": [
            "clin_au_endo__PD_STATUS",
        ],
        "canon__exitus": [
            "clin_au_endo__EXITUS",
            "clin_dcs__Exitus",         # â† v4: breast exitus
        ],
        "canon__vital_status": [         # â† v4: new â€” verbose status string
            "clin_dcs__STATUS",
        ],
        "canon__recurrence": [           # â† v4: new â€” any recurrence/progression
            "clin_dcs__Recaida/Progresión",
            "clin_au_endo__PD_STATUS",
        ],
        "canon__distant_mets": [         # â† v4: new â€” distant metastasis
            "clin_dcs__MTxDISTANCIA",
            "clin_au_endo__M",
            "clin_au_endo__PD_LOCATION",
        ],
        "canon__treatment": [            # â† v4: new â€” treatment description
            "clin_dcs__Tto",
            "clin_au_endo__PD_TREATMENT",
        ],

        # â”€â”€ Dates â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        "canon__date_diagnosis": [
            "clin_dcs__Fecha dx",
            "clin_au_endo__DATE_OF_SURGERY",
        ],
        "canon__date_last_fu": [         # â† v4: new â€” last follow-up (used to derive OS)
            "clin_dcs__Última fecha disponible",
            "clin_au_endo__LAST_UPDATED",
        ],
        "canon__date_recurrence": [      # â† v4: new
            "clin_dcs__Fecha Recidiva/Progresión",
            "clin_au_endo__PD_DATE",
        ],

        # â”€â”€ Molecular (endometrial) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        "canon__msi_status": [
            "clin_au_endo__MSI_STATUS_IHC",
            "clin_au_endo__UA_MSI_STATUS_NGS",
            "clin_au_endo__UA_MSI_STATUS_ddPCR",
        ],
        "canon__molecular_class": [
            "clin_au_endo__MOLECULAR CLASSIFICATION_according to IHC and/or NGS profile",
            "clin_au_endo__POLAND MOLECULAR CLASSIFICATION",
        ],

        # â”€â”€ Operational â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        "canon__performed_or_cuts": [
            "clin_mn__REALIZADO",
            "clin_en__CORTES",
            "clin_her2__REALIZADO",
            "clin_her2__REALIZADO.1",
            "clin_her2__REPETICI?N DE CORTES",
        ],
    }


# =============================================================================
# SEMANTIC TRANSFORMS (turn canon columns into comparable formats)
# =============================================================================

def apply_semantic_transforms(df: pd.DataFrame) -> pd.DataFrame:
    """
    After coalescing canon columns, standardise types/encodings and add derived variables.

    v4 additions
    ------------
    - canon__bmi         : strip 'X', 'Sobrepeso', weight-only entries before numeric coercion
    - canon__grade       : now also parses numeric values from DCS (already integers, safe)
    - canon__ki67        : new â€” unify endo (fraction 0-1) and breast (% or range string)
                           into a single canon__ki67_pct (percentage scale 0-100)
    - canon__p53_ihc     : new â€” unify endo (fraction) and breast (fraction) -> Positive/Negative
    - canon__exitus      : merged from two cohorts â€” normalise SI/NO -> Yes/No
    - canon__os_months   : v4 also DERIVES breast OS from canon__date_diagnosis +
                           canon__date_last_fu when canon__os_months is null
    - canon__recurrence_flag : new binary derived from canon__recurrence text
    - canon__distant_mets_flag: new binary derived from canon__distant_mets
    - canon__myometrial_invasion: new â€” normalise <50% / >50%
    - canon__lvsi        : new â€” normalise YES/NO
    - canon__risk_group  : new â€” normalise LOW/INTERMEDIATE/INTERMEDIATE-HIGH/HIGH
    - canon__vital_status: new â€” normalise verbose string -> Alive / Dead / Lost / Other
    - canon__performed_or_cuts__yesno : fix â€” was always null due to logic error; now fixed
    - MSI: 'UNSTABLE' (uppercase raw) now correctly normalised to 'Unstable'
    """
    out = df.copy()

    # â”€â”€ Strip whitespace from key categoricals â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    for _col in [
        "clin_au_endo__FIGO_STAGE",
        "clin_her2__DX",
        "clin_dcs__Dx",
        "canon__msi_status",
        "clin_au_endo__GRADE",
        "clin_au_endo__LVSI",
        "clin_au_endo__MYOMETRIAL_INFILTRATION",
        "clin_au_endo__HISTOLOGY_GROUP",
        "clin_au_endo__PD_STATUS",
        "clin_au_endo__EXITUS",
        "clin_au_endo__RISK_OF_RECURRENCE",
        "clin_dcs__GRADO",
        "clin_dcs__RE",
        "clin_dcs__RP",
        "clin_dcs__MTxDISTANCIA",
        "clin_dcs__Exitus",
        "clin_dcs__STATUS",
        "clin_dcs__Recaida/Progresión",
    ]:
        if _col in out.columns:
            out[_col] = out[_col].apply(lambda x: x.strip() if isinstance(x, str) else x)

    # â”€â”€ MSI status -> 'Stable' / 'Unstable' â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__msi_status" in out.columns:
        def _norm_msi(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if s in ("MSS", "STABLE", "MS-STABLE"): return "Stable"
            if s in ("MSI", "MSI-H", "UNSTABLE", "MSI-HIGH", "MSIH"): return "Unstable"
            return x
        out["canon__msi_status"] = out["canon__msi_status"].map(_norm_msi)

    # â”€â”€ Molecular class -> POLE / MMRd / NSMP / P53 â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__molecular_class" in out.columns:
        def _norm_molclass(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if "POLE" in s: return "POLE"
            if "MMR" in s or "MMRD" in s: return "MMRd"
            if "NSMP" in s or "NO SPECIFIC" in s: return "NSMP"
            if "P53" in s or "TP53" in s: return "P53"
            return x
        out["canon__molecular_class"] = out["canon__molecular_class"].map(_norm_molclass)

    # â”€â”€ Age â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__age" in out.columns:
        out["canon__age"] = pd.to_numeric(out["canon__age"], errors="coerce")

    # â”€â”€ BMI: strip non-numeric junk before coercion â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    # DCS has 'X' (unknown), 'Sobrepeso' (overweight text), '73 kg', '60 kg; 168 cm'
    if "canon__bmi" in out.columns:
        bmi_candidates = [
            "clin_au_endo__BMI_CALCULATED",
            "clin_au_endo__BMI",
            "clin_her2__BMI",
            "clin_dcs__BMI",
            "clin_mn__BMI",
            "clin_en__BMI",
        ]

        def _pick_bmi(row) -> Tuple[Optional[float], Optional[str]]:
            for col in bmi_candidates:
                if col not in row.index:
                    continue
                parsed = parse_bmi_value(row[col])
                if parsed is not None:
                    return parsed, col
            return None, None

        bmi_parsed = out.apply(_pick_bmi, axis=1, result_type="expand")
        out["canon__bmi"] = pd.to_numeric(bmi_parsed[0], errors="coerce")
        if "canon__bmi__source" in out.columns:
            out["canon__bmi__source"] = bmi_parsed[1].where(bmi_parsed[1].notna(), out["canon__bmi__source"])

        def _pick_weight_kg(row) -> Tuple[Optional[float], Optional[str]]:
            for col in bmi_candidates:
                if col not in row.index:
                    continue
                parsed = parse_weight_kg_value(row[col])
                if parsed is not None:
                    return parsed, col
            return None, None

        def _pick_height_cm(row) -> Tuple[Optional[float], Optional[str]]:
            for col in bmi_candidates:
                if col not in row.index:
                    continue
                parsed = parse_height_cm_value(row[col])
                if parsed is not None:
                    return parsed, col
            return None, None

        weight_parsed = out.apply(_pick_weight_kg, axis=1, result_type="expand")
        height_parsed = out.apply(_pick_height_cm, axis=1, result_type="expand")
        out["canon__weight_kg"] = pd.to_numeric(weight_parsed[0], errors="coerce")
        out["canon__height_cm"] = pd.to_numeric(height_parsed[0], errors="coerce")
        out["canon__weight_kg__source"] = weight_parsed[1]
        out["canon__height_cm__source"] = height_parsed[1]

        def _pick_bmi_category(row) -> Tuple[Optional[str], Optional[str]]:
            for col in bmi_candidates:
                if col not in row.index:
                    continue
                category = parse_bmi_category(row[col])
                if category is not None:
                    return category, col
            return None, None

        bmi_category = out.apply(_pick_bmi_category, axis=1, result_type="expand")
        out["canon__bmi_category"] = bmi_category[0]
        out["canon__bmi_category__source"] = bmi_category[1]

    # â”€â”€ Dates â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    for _dcol in ("canon__date_birth", "canon__date_diagnosis",
                  "canon__date_last_fu", "canon__date_recurrence"):
        if _dcol in out.columns:
            out[_dcol] = out[_dcol].apply(parse_date)

    # â”€â”€ Derive breast OS from diagnosis + last follow-up dates â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    # Only fills rows where canon__os_months is null (i.e. non-AU-endo rows)
    if "canon__date_diagnosis" in out.columns and "canon__date_last_fu" in out.columns:
        diag = out["canon__date_diagnosis"]
        last = out["canon__date_last_fu"]
        derived = ((last - diag).dt.days / 30.44).round(2)
        # Sanity gate: must be positive and < 600 months (~50 years)
        derived = derived.where((derived > 0) & (derived < 600), other=pd.NA)
        if "canon__os_months" not in out.columns:
            out["canon__os_months"] = pd.NA
        fill_mask = out["canon__os_months"].isna() & derived.notna()
        out.loc[fill_mask, "canon__os_months"] = derived[fill_mask]
        # label the source of OS for traceability
        if "canon__os_months__source" not in out.columns:
            out["canon__os_months__source"] = pd.NA
        out.loc[fill_mask, "canon__os_months__source"] = "derived_from_dates"

    # â”€â”€ Numeric coercion for survival endpoints â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    for _surv_col in ("canon__os_months", "canon__pfs_months"):
        if _surv_col in out.columns:
            out[_surv_col] = pd.to_numeric(out[_surv_col], errors="coerce")

    # â”€â”€ PD flag (endometrial) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__pd_status" in out.columns:
        out["canon__pd_flag"] = out["canon__pd_status"].map(
            lambda x: 1 if str(x).strip().upper() == "PD" else (0 if str(x).strip().upper() == "NO PD" else pd.NA)
            if pd.notna(x) else pd.NA
        )

    # â”€â”€ Exitus: normalise SI/NO + YES/NO -> Yes/No â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__exitus" in out.columns:
        def _norm_exitus(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if s in ("SI", "YES", "EXITUS", "1"): return "Yes"
            if s in ("NO", "0"): return "No"
            # Breast: 'SI: fallo hepÃ¡tico...' -> Yes
            if s.startswith("SI:"): return "Yes"
            return x
        out["canon__exitus"] = out["canon__exitus"].map(_norm_exitus)
        # Derived binary (1/0)
        out["canon__exitus_flag"] = out["canon__exitus"].map(
            lambda x: pd.NA if pd.isna(x) else (1 if x == "Yes" else (0 if x == "No" else pd.NA))
        )

    # â”€â”€ Vital status: normalise verbose string -> clean categories â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__vital_status" in out.columns:
        def _norm_status(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if any(t in s for t in ("SIN ENFERMEDAD", "VIVO", "VIVA")): return "Alive_NED"
            if "EXITUS" in s or "FALLEC" in s: return "Dead"
            if "PERDIDO" in s or "SEGUIMIENTO" in s: return "Lost_FU"
            if "PALIATIV" in s or "ECOG 4" in s or "RECIDIVA" in s: return "Alive_Disease"
            return "Other"
        out["canon__vital_status"] = out["canon__vital_status"].map(_norm_status)

    # â”€â”€ Recurrence flag â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__recurrence" in out.columns:
        def _norm_recurrence(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if s in ("NO", "NO PD"): return "No"
            if s in ("PD", "YES", "YES PD"): return "Yes"
            if len(s) > 3: return "Yes"   # any non-trivial text = recurrence described
            return x
        out["canon__recurrence"] = out["canon__recurrence"].map(_norm_recurrence)
        out["canon__recurrence_flag"] = out["canon__recurrence"].map(
            lambda x: pd.NA if pd.isna(x) else (1 if x == "Yes" else (0 if x == "No" else pd.NA))
        )

    # â”€â”€ Distant metastasis flag â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__distant_mets" in out.columns:
        def _norm_mets(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if s in ("SI", "M1"): return "Yes"
            if s.startswith("DISTANT"): return "Yes"
            if s in ("NO", "NO-LOCAL", "M0", "LOCAL"): return "No"
            return x
        out["canon__distant_mets"] = out["canon__distant_mets"].map(_norm_mets)
        out["canon__distant_mets_flag"] = out["canon__distant_mets"].map(
            lambda x: pd.NA if pd.isna(x) else (1 if x == "Yes" else (0 if x == "No" else pd.NA))
        )

    # â”€â”€ Myometrial invasion -> <50% / >50% / None â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__myometrial_invasion" in out.columns:
        def _norm_myo(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if "<50" in s: return "<50%"
            if ">50" in s: return ">50%"
            return x
        out["canon__myometrial_invasion"] = out["canon__myometrial_invasion"].map(_norm_myo)

    # â”€â”€ LVSI -> Yes / No â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__lvsi" in out.columns:
        out["canon__lvsi"] = out["canon__lvsi"].map(
            lambda x: "Yes" if str(x).strip().upper() == "YES"
            else ("No" if str(x).strip().upper() == "NO" else (pd.NA if pd.isna(x) else x))
        )

    # â”€â”€ Risk group: normalise to LOW / INTERMEDIATE / INTERMEDIATE-HIGH / HIGH â”€
    if "canon__risk_group" in out.columns:
        def _norm_risk(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if s == "LOW": return "Low"
            if s == "INTERMEDIATE": return "Intermediate"
            if s == "INTERMEDIATE-HIGH": return "Intermediate-High"
            if s == "HIGH": return "High"
            return x
        out["canon__risk_group"] = out["canon__risk_group"].map(_norm_risk)

    # â”€â”€ KI67: unify into a single percentage-scale column â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    # Endo source (FFPE_KI67): stored as percentage (0-100)
    # Breast DCS source (clin_dcs__KI67): stored as fraction (0.0-1.0) for most rows,
    #   but some are range strings like '20-25%', '>50%', '70-75%'
    if "canon__ki67" in out.columns:
        out["canon__ki67_raw"] = out["canon__ki67"]
        out["canon__ki67_pct"] = out.apply(
            lambda row: parse_ki67_percent(row["canon__ki67"], row.get("canon__ki67__source")),
            axis=1,
        )

    # â”€â”€ P53 IHC -> Positive / Negative (unified) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    # Endo: stored as fraction (FFPE_TP53_IHC); breast: fraction (clin_dcs__p53)
    # Convention: >= 0.1 (10%) = aberrant/positive p53
    if "canon__p53_ihc" in out.columns:
        out["canon__p53_ihc_raw"] = out["canon__p53_ihc"]
        p53_parsed = out["canon__p53_ihc"].apply(parse_p53_ihc_value)
        out["canon__p53_ihc_numeric"] = p53_parsed.map(lambda t: t[0])
        out["canon__p53_status"] = p53_parsed.map(lambda t: t[1] if t[1] is not None else pd.NA)

    # â”€â”€ Menarche â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__menarche_age" in out.columns:
        out["canon__menarche_age"] = out["canon__menarche_age"].apply(extract_first_number)

    # â”€â”€ Menopause â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    def _parse_menopause_age(x) -> Optional[float]:
        if pd.isna(x): return None
        s = _norm_text(x)
        nums = [float(n) for n in re.findall(r"\b(\d+(?:\.\d+)?)\b", s)]
        for n in nums:
            if 20 <= n <= 65:
                return n
        return None

    if "canon__menopause_age_or_status" in out.columns:
        s = out["canon__menopause_age_or_status"]
        out["canon__menopause_age"] = s.apply(_parse_menopause_age)
        out["canon__menopause_status"] = s.apply(lambda x: _norm_text(x) if _norm_text(x) else pd.NA)

    # â”€â”€ Grade: numeric coercion (covers both endo text 'G1' and breast integer) â”€
    if "canon__grade" in out.columns:
        out["canon__grade"] = out["canon__grade"].apply(parse_grade)

    # â”€â”€ FIGO stage â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__figo_stage" in out.columns:
        out["canon__figo_stage"] = out["canon__figo_stage"].apply(parse_figo_stage)

    # â”€â”€ AU endometrial ER/PR: numeric raw -> Positive/Negative (threshold >= 1%) â”€
    if "canon__er_raw" in out.columns:
        er_parsed = out["canon__er_raw"].apply(parse_receptor_value)
        out["canon__er_pct"] = er_parsed.map(lambda t: t[0])
        out["canon__er_status"] = er_parsed.map(lambda t: t[1] if t[1] is not None else pd.NA)
    if "canon__pr_raw" in out.columns:
        pr_parsed = out["canon__pr_raw"].apply(parse_receptor_value)
        out["canon__pr_pct"] = pr_parsed.map(lambda t: t[0])
        out["canon__pr_status"] = pr_parsed.map(lambda t: t[1] if t[1] is not None else pd.NA)
    for _col in ("canon__er_status", "canon__pr_status"):
        if _col not in out.columns:
            out[_col] = pd.NA

    # â”€â”€ SP breast ER/PR: POSITIVO/NEGATIVO -> Positive/Negative â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    _sp_er_map = {"POSITIVO": "Positive", "NEGATIVO": "Negative"}
    if "canon__er_status_breast" in out.columns:
        out["canon__er_status_breast"] = (
            out["canon__er_status_breast"].astype(str).str.strip().str.upper().map(
                lambda x: _sp_er_map.get(x, pd.NA)
            )
        )
        mask = out["canon__er_status"].isna() & out["canon__er_status_breast"].notna()
        out.loc[mask, "canon__er_status"] = out.loc[mask, "canon__er_status_breast"]

    if "canon__pr_status_breast" in out.columns:
        out["canon__pr_status_breast"] = (
            out["canon__pr_status_breast"].astype(str).str.strip().str.upper().map(
                lambda x: _sp_er_map.get(x, pd.NA)
            )
        )
        mask = out["canon__pr_status"].isna() & out["canon__pr_status_breast"].notna()
        out.loc[mask, "canon__pr_status"] = out.loc[mask, "canon__pr_status_breast"]

    # â”€â”€ HER2 copies â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__her2_copies" in out.columns:
        out["canon__her2_copies"] = out["canon__her2_copies"].apply(extract_first_number)

    # â”€â”€ HER2 note + triple negative flag â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if "canon__her2_copies_note" in out.columns:
        note = out["canon__her2_copies_note"].copy()
        note_num = note.apply(extract_first_number)
        out["canon__her2_copies_note"] = note.where(note_num.isna(), pd.NA)
        out["canon__triple_negative_flag"] = out["canon__her2_copies_note"].astype(str).str.upper().str.contains(
            "TRIPLE NEGATIVO", na=False
        ).map({True: "Yes", False: "No"})
        out.loc[out["canon__her2_copies_note"].isna(), "canon__triple_negative_flag"] = pd.NA

    # ?????? Breast treatment text -> exploratory exposure flags ??????????????????????????????????????????????????????
    if "canon__treatment" in out.columns:
        tx_flags = out["canon__treatment"].apply(derive_breast_treatment_exposure_flags).apply(pd.Series)
        for col in tx_flags.columns:
            out[col] = tx_flags[col]
            source_col = f"{col}__source"
            out[source_col] = pd.NA
            if "canon__treatment__source" in out.columns:
                valid_mask = out[col].notna()
                out.loc[valid_mask, source_col] = out.loc[valid_mask, "canon__treatment__source"]

    # â”€â”€ Performed/cuts: fix the always-null bug (was checking 'nan' string) â”€â”€â”€
    if "canon__performed_or_cuts" in out.columns:
        out["canon__performed_or_cuts__yesno"] = out["canon__performed_or_cuts"].apply(
            lambda x: "Yes" if pd.notna(x) and str(x).strip().lower() not in ("", "nan", "none") else pd.NA
        )

    return out


# =============================================================================
# HARMONISE PIPELINE
# =============================================================================

def harmonise(df_raw_merged: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Creates:
      - harm: raw + canon__* + canon__*__source + derived fields
      - audit: coverage summary for canon columns
      - mapping_doc: documents the harmonisation map
    """
    df = df_raw_merged.copy()

    # Drop Excel artefact columns that are just formatting leftovers
    unnamed_cols = [c for c in df.columns if str(c).startswith("Unnamed:")]
    df = df.drop(columns=unnamed_cols, errors="ignore")

    harm_map = build_harmonisation_map()

    harm = df.copy()
    map_rows = []

    # Create canonical columns + provenance
    for canon, sources in harm_map.items():
        series, src = coalesce_with_source(harm, sources)
        harm[canon] = series
        harm[f"{canon}__source"] = src
        map_rows.append({
            "canon_column": canon,
            "source_priority_list": " | ".join(sources),
            "n_sources_found_in_data": sum(1 for c in sources if c in harm.columns),
        })

    # Apply semantic transforms to canon columns (types, derived flags, etc.)
    harm = apply_semantic_transforms(harm)

    for canon in [
        "canon__treatment_chemotherapy",
        "canon__treatment_anti_her2",
        "canon__treatment_endocrine",
        "canon__treatment_radiotherapy",
    ]:
        if canon in harm.columns:
            map_rows.append({
                "canon_column": canon,
                "source_priority_list": "Derived from canon__treatment using keyword-based exploratory exposure flags",
                "n_sources_found_in_data": int("canon__treatment" in harm.columns),
            })

    if "canon__bmi_category" in harm.columns:
        map_rows.append({
            "canon_column": "canon__bmi_category",
            "source_priority_list": "Derived from BMI source text and numeric BMI values; preserves labels such as 'Sobrepeso' and standard BMI categories",
            "n_sources_found_in_data": int("canon__bmi" in harm.columns),
        })
    if "canon__weight_kg" in harm.columns:
        map_rows.append({
            "canon_column": "canon__weight_kg",
            "source_priority_list": "Derived from BMI source text when kg is present; preserves weight-only values such as '53.5 kg' and mixed strings like '60 kg; 168 cm'",
            "n_sources_found_in_data": int("canon__bmi" in harm.columns),
        })
    if "canon__height_cm" in harm.columns:
        map_rows.append({
            "canon_column": "canon__height_cm",
            "source_priority_list": "Derived from BMI source text when cm is present; preserves height information in mixed strings like '60 kg; 168 cm'",
            "n_sources_found_in_data": int("canon__bmi" in harm.columns),
        })

    # â”€â”€ Drop columns that are >99% empty â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    # Keeps the file clean by removing placeholder/unfilled columns (e.g. the
    # AU endometrial MUTATION TYPE detail columns which are all blank).
    # Canon columns and core identity columns are always protected regardless
    # of fill rate, so we never accidentally drop a harmonised output variable.
    PROTECTED_PREFIXES = ("canon__",)
    PROTECTED_COLS = {
        "snp_code", "nucleic_acid", "tissue", "tumour_normal",
        "sample_id", "case_id", "block_type", "sheet",
        "concentration_ng_ul", "volume_ul", "extraction_flag",
        "pd_status", "histology", "endo_id", "NHC", "is_sequenced",
    }
    n_rows = len(harm)
    threshold = 0.99  # drop if >= 99% null

    cols_to_drop = []
    for col in harm.columns:
        # Never drop protected columns
        if col in PROTECTED_COLS:
            continue
        if any(col.startswith(p) for p in PROTECTED_PREFIXES):
            continue
        null_frac = harm[col].isna().sum() / n_rows if n_rows > 0 else 1.0
        if null_frac >= threshold:
            cols_to_drop.append(col)

    if cols_to_drop:
        print(f"[INFO] Dropping {len(cols_to_drop)} columns that are â‰¥99% empty:")
        for c in cols_to_drop:
            print(f"  - {c}")
        harm = harm.drop(columns=cols_to_drop)
        print(f"[INFO] Columns after cleanup: {len(harm.columns)}")

    # Audit table: how many values for each canonical variable
    canon_cols = [c for c in harm.columns if c.startswith("canon__") and not c.endswith("__source")]
    audit = pd.DataFrame({
        "canon_column": canon_cols,
        "non_null_rows": [int(harm[c].notna().sum()) for c in canon_cols],
        "pct_filled": [round(100 * harm[c].notna().sum() / n_rows, 1) for c in canon_cols],
    }).sort_values("non_null_rows", ascending=False)

    mapping_doc = pd.DataFrame(map_rows).sort_values("canon_column")

    return harm, audit, mapping_doc


# =============================================================================
# MANIFEST LOADING + SEQUENCED SAMPLE TAGGING
# =============================================================================

DEFAULT_MANIFESTS_DIR = ROOT_DIR / "manifests"

def _extract_manifest_snp_code(sample_name: str) -> Optional[str]:
    """Parse sequencer BAM names into canonical snp_code values."""
    s = _to_str(sample_name)
    m = re.match(r"^(SNP_(?:AT|EN|MN|MT-T)_\d+)", s, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    m = re.match(r"^DNA_SNP_(EN|MN)_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_{m.group(1).upper()}_{m.group(2)}"
    m = re.match(r"^DNA_SNP_AT_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_AT_{m.group(1)}"
    m = re.match(r"^DNA_AT_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_AT_{m.group(1)}"
    m = re.match(r"^SNP_DNA_AT_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_AT_{m.group(1)}"
    m = re.match(r"^SNP_DNA_(MN|EN)_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_{m.group(1).upper()}_{m.group(2)}"
    m = re.match(r"^DNA_SNP_MT-T_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_MT-T_{m.group(1)}"
    m = re.match(r"^DNA_MT[-_]T_(\d+)_", s, re.IGNORECASE)
    if m:
        return f"SNP_MT-T_{m.group(1)}"
    return None


def load_sequenced_snp_codes(manifests_dir: Path) -> set:
    """
    Read all pass manifests and return canonical snp_code values for samples
    with passing BAMs. Replicate BAMs are ignored for highlighting purposes.
    """
    sequenced: set = set()

    if not manifests_dir.exists():
        print(f"[WARN] Manifests directory not found: {manifests_dir} - sequencing flags will be skipped.")
        return sequenced

    manifest_paths = sorted(manifests_dir.glob('*pass_manifest.txt'))
    if not manifest_paths:
        print(f"[WARN] No pass manifests found in: {manifests_dir}")
        return sequenced

    for fpath in manifest_paths:
        with open(fpath) as f:
            bams = [line.strip() for line in f if line.strip()]

        for bam in bams:
            bam_name = Path(bam).stem
            if 'repeticion' in bam_name.lower():
                continue
            code = _extract_manifest_snp_code(bam_name)
            if code is None:
                print(f"[WARN] Could not parse manifest sample name: {bam_name}")
                continue
            sequenced.add(code.upper())

    print(f"[INFO] Sequenced samples found in manifests: {len(sequenced)}")
    return sequenced


def tag_sequenced_samples(df: pd.DataFrame, sequenced_codes: set) -> pd.DataFrame:
    """
    Adds a boolean column 'is_sequenced' to the dataframe based on snp_code match
    against the set of passing BAM snp_codes.
    """
    def _check(code):
        if pd.isna(code):
            return False
        # Normalise: uppercase, remove spaces, strip _REP suffix
        c = str(code).upper().strip().replace(" ", "")
        c = re.sub(r"_REP(ETICION)?.*$", "", c, flags=re.IGNORECASE)
        return c in sequenced_codes

    df = df.copy()
    df["is_sequenced"] = df["snp_code"].apply(_check)
    return df


# =============================================================================
# EXCEL HIGHLIGHTING
# =============================================================================

GREEN_FILL  = PatternFill("solid", start_color="C8E6C9", end_color="C8E6C9")   # light green
HEADER_FILL = PatternFill("solid", start_color="1B5E20", end_color="1B5E20")   # dark green header
HEADER_FONT = Font(bold=True, color="FFFFFF")


def apply_sequenced_highlighting(ws, df: pd.DataFrame) -> None:
    """
    For the harmonised_plus_canon sheet:
      - Highlights the entire row green for sequenced samples (is_sequenced == True)
      - Adds a bold dark-green header to the is_sequenced column
      - Freezes the top row
    """
    if "is_sequenced" not in df.columns:
        return

    seq_col_idx = df.columns.get_loc("is_sequenced") + 1  # 1-based for openpyxl
    n_cols = len(df.columns)

    # Style the is_sequenced column header
    header_cell = ws.cell(row=1, column=seq_col_idx)
    header_cell.fill = HEADER_FILL
    header_cell.font = HEADER_FONT

    # Freeze top row
    ws.freeze_panes = "A2"

    # Highlight sequenced rows
    for row_idx, is_seq in enumerate(df["is_sequenced"], start=2):  # row 1 = header
        if is_seq:
            for col_idx in range(1, n_cols + 1):
                ws.cell(row=row_idx, column=col_idx).fill = GREEN_FILL


# =============================================================================
# MAIN
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build master SNP + clinical and create UPDATED Option-B harmonised clinical layer.")
    p.add_argument("--snp_xlsx", type=str, default=str(DEFAULT_SNP_XLSX), help="Path to SNP workbook.")
    p.add_argument("--clinical_xlsx", type=str, default=str(DEFAULT_CLINICAL_XLSX), help="Path to main clinical workbook.")
    p.add_argument("--au_xlsx", type=str, default=str(DEFAULT_AU_XLSX), help="Path to AU endometrial tumour clinical workbook.")
    p.add_argument("--au_sheet", type=str, default=DEFAULT_AU_SHEET, help="Sheet name in AU workbook.")
    p.add_argument("--out_xlsx", type=str, default=str(DEFAULT_OUT_XLSX), help="Output Excel path.")
    p.add_argument("--manifests_dir", type=str, default=str(DEFAULT_MANIFESTS_DIR),
                   help="Directory containing pass_manifest.txt files for sequencing status highlighting.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    snp_xlsx      = Path(args.snp_xlsx)
    clinical_xlsx = Path(args.clinical_xlsx)
    au_xlsx       = Path(args.au_xlsx)
    out_xlsx      = Path(args.out_xlsx)
    manifests_dir = Path(args.manifests_dir)

    # Input checks
    validate_file_exists(snp_xlsx, "Script 16b SNP workbook")
    validate_file_exists(clinical_xlsx, "Script 16b clinical workbook")

    if au_xlsx.exists():
        print(f"AU clinical workbook found: {au_xlsx} (sheet={args.au_sheet})")
    else:
        print(f"AU clinical workbook NOT found: {au_xlsx} (will skip AU merge)")

    # 1) Load sequenced sample codes from manifests
    sequenced_codes = load_sequenced_snp_codes(manifests_dir)

    # 2) Build SNP master
    master_df = build_master(snp_xlsx)
    validate_nonempty(master_df, "Script 16b master dataframe")
    print_validation_summary(master_df, None, "Script 16b master dataframe")
    print(f"MASTER: {master_df.shape}")

    # 3) Merge clinical (lossless)
    merged_df = merge_master_with_clinical(master_df, clinical_xlsx, au_xlsx, args.au_sheet)
    print(f"RAW MERGED: {merged_df.shape}")

    # 4) Harmonise (adds canon variables)
    harm_df, audit_df, map_df = harmonise(merged_df)
    clinical_merge_summary_df, clinical_match_detail_df, clinical_unmatched_df, clinical_only_df = build_clinical_merge_audits(merged_df, au_xlsx, args.au_sheet)
    validate_nonempty(harm_df, "Script 16b harmonised dataframe")
    print_validation_summary(harm_df, None, "Script 16b harmonised dataframe")
    print(f"HARMONISED: {harm_df.shape}")

    # 5) Tag sequenced samples in both sheets
    merged_df = tag_sequenced_samples(merged_df, sequenced_codes)
    harm_df   = tag_sequenced_samples(harm_df,   sequenced_codes)

    n_seq = harm_df["is_sequenced"].sum()
    print(f"Sequenced samples tagged: {n_seq} / {len(harm_df)} rows")

    # 6) Write output Excel with green highlighting on sequenced rows
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        merged_df.to_excel(xw, index=False, sheet_name="raw_merged")
        harm_df.to_excel(xw,   index=False, sheet_name="harmonised_plus_canon")
        audit_df.to_excel(xw,  index=False, sheet_name="canon_audit")
        map_df.to_excel(xw,    index=False, sheet_name="harmonisation_map")
        clinical_merge_summary_df.to_excel(xw, index=False, sheet_name="clinical_merge_audit")
        clinical_match_detail_df.to_excel(xw, index=False, sheet_name="au_endo_match_detail")
        clinical_unmatched_df.to_excel(xw, index=False, sheet_name="au_endo_unmatched")
        clinical_only_df.to_excel(xw, index=False, sheet_name="au_endo_clinical_only")

        # Apply green row highlighting to harmonised sheet
        ws_harm = xw.sheets["harmonised_plus_canon"]
        apply_sequenced_highlighting(ws_harm, harm_df)

        # Also highlight raw_merged sheet
        ws_raw = xw.sheets["raw_merged"]
        apply_sequenced_highlighting(ws_raw, merged_df)

    print(f"Saved: {out_xlsx}")
    print(f"Green rows = sequenced samples with passing BAMs ({n_seq} rows highlighted)")
    if not clinical_merge_summary_df.empty:
        unmatched_unique = clinical_merge_summary_df.loc[clinical_merge_summary_df["Metric"].eq("Unmatched_AT_AUs_unique_ids"), "Value"].iloc[0]
        ambiguous_unique = clinical_merge_summary_df.loc[clinical_merge_summary_df["Metric"].eq("Ambiguous_AT_AUs_unique_ids"), "Value"].iloc[0]
        print(f"AU endometrial clinical audit: unmatched unique IDs = {unmatched_unique}; ambiguous unique IDs = {ambiguous_unique}")


if __name__ == "__main__":
    main()
