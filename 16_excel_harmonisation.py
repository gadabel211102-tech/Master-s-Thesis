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

RUN (WSL example)
-----------------
python3 16_excel_harmonisation.py \
  --snp_xlsx "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/Muestras SNPs nomenclaturas, equivalencias y cuantificaciones.xlsx" \
  --clinical_xlsx "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/MUESTRAS + VARIABLES CLINICAS PROY SNPs GM(4).xlsx" \
  --au_xlsx "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/Clinical DATA AU_Endometrial cancer.xlsx" \
  --out_xlsx "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx" \
  --manifests_dir "/home/gadeaalonsoj/tfm/manifests"

NOTE ON PERMISSION ERRORS
-------------------------
If writing to OneDrive gives PermissionError, the file is usually open in Excel or syncing.
Write to ~/tfm first, then copy to OneDrive, or change --out_xlsx to a new filename.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from openpyxl.styles import PatternFill, Font
from openpyxl.utils import get_column_letter

from pipeline_validation import print_validation_summary, validate_file_exists, validate_nonempty


# =============================================================================
# DEFAULT CONFIG (override with CLI args)
# =============================================================================
# These defaults are convenient if you always run from WSL and keep files in OneDrive.
# You can override any path at runtime using the command-line arguments.

DEFAULT_SNP_XLSX = Path(
    "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/"
    "Muestras SNPs nomenclaturas, equivalencias y cuantificaciones.xlsx"
)
DEFAULT_CLINICAL_XLSX = Path(
    "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/"
    "MUESTRAS + VARIABLES CLINICAS PROY SNPs GM(4).xlsx"
)
DEFAULT_AU_XLSX = Path(
    "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/"
    "Clinical DATA AU_Endometrial cancer.xlsx"
)
DEFAULT_AU_SHEET = "Hoja1"
DEFAULT_OUT_XLSX = Path("/home/gadeaalonsoj/tfm/MASTER_SNP_plus_clinical_HARMONISED.xlsx")

# Expected SNP workbook sheets we process
SNP_SHEETS = ["AT=AUs", "EN", "MT-T_N", "MN"]  # ignore OVSER


# =============================================================================
# BASIC HELPERS
# =============================================================================

def _to_str(x) -> str:
    """Safe string conversion for messy Excel cells."""
    return "" if pd.isna(x) else str(x).strip()


def standard_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Strip whitespace from column names to reduce accidental mismatches."""
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
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


def _pick_existing(df: pd.DataFrame, candidates: List[str]) -> str:
    """
    Choose the first existing column among candidates.
    Useful because some sheets have slightly different headers (accents/casing).
    """
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        f"None of these columns were found: {candidates}. "
        f"Available: {df.columns.tolist()}"
    )


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
        "nucleic_acid": df["Ác. Nucleico"].map(lambda x: _to_str(x).upper()),
        "concentration_ng_ul": pd.to_numeric(df["Qubit (ng/uL)"], errors="coerce"),
        "volume_ul": None,
        "extraction_flag": df.get("extracción", pd.Series([None] * len(df))).map(_to_str),
        "pd_status": df.get("PD status", pd.Series([None] * len(df))).map(_to_str),
        "histology": df.get("Histology", pd.Series([None] * len(df))).map(_to_str),
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
        "sample_id": df["Código Noray"].map(_to_str),
        "nucleic_acid": df["Ác. Nucleico"].map(lambda x: _to_str(x).upper()),
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
        "sample_id": df[_pick_existing(df, ["CÓDIGO BB BLOQUE TUMORAL", "CODIGO BB BLOQUE TUMORAL"])].map(_to_str),
        "snp_code": df[_pick_existing(df, ["CODIGO SNPs"])].map(normalize_snp_code),
        "nucleic_acid": df[_pick_existing(df, ["Ác. Nucleico"])].map(lambda x: _to_str(x).upper()),
        "volume_ul": df[_pick_existing(df, ["uL"])].map(parse_volume_to_ul),
        "concentration_ng_ul": pd.to_numeric(df[_pick_existing(df, ["EXTRACCION (ng/ul)"])], errors="coerce"),
    })

    right = pd.DataFrame({
        "sample_id": df[_pick_existing(df, ["CÓDIGO BB BLOQUE TUMORAL.1", "CODIGO BB BLOQUE TUMORAL.1"])].map(_to_str),
        "snp_code": df[_pick_existing(df, ["CODIGO SNPs.1"])].map(normalize_snp_code),
        "nucleic_acid": df[_pick_existing(df, ["Ác. Nucleico.1"])].map(lambda x: _to_str(x).upper()),
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
        "sample_id": df[_pick_existing(df, ["Código Noray-BB", "Codigo Noray-BB", "Código Noray"])].map(_to_str),
        "snp_code": df[_pick_existing(df, ["CODIGO SNP_MN (mama normal)"])].map(normalize_snp_code),
        "nucleic_acid": df[_pick_existing(df, ["Ác. Nucleico"])].map(lambda x: _to_str(x).upper()),
        "volume_ul": df[_pick_existing(df, ["VOL.", "VOL (uL)"])].map(parse_volume_to_ul),
        "concentration_ng_ul": pd.to_numeric(df[_pick_existing(df, ["Qubit (ng/uL)"])], errors="coerce"),
    })

    right = pd.DataFrame({
        "sample_id": df[_pick_existing(df, ["Código Noray-BB.1", "Codigo Noray-BB.1", "Código Noray.1"])].map(_to_str),
        "snp_code": df[_pick_existing(df, ["CODIGO SNP_MN (mama normal).1"])].map(normalize_snp_code),
        "nucleic_acid": df[_pick_existing(df, ["Ác. Nucleico.1"])].map(lambda x: _to_str(x).upper()),
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
        df = df.rename(columns={"Código Noray-BB": "sample_id"})
        df["sample_id"] = df["sample_id"].map(_to_str)
        tables["clinical_breast_healthy"] = df.rename(columns={c: f"clin_mn__{c}" for c in df.columns if c != "sample_id"})

    # Healthy endometrium (EN): join by sample_id
    if "TANDA ENDOMETRIO SANO" in xls.sheet_names:
        df = standard_columns(pd.read_excel(xls, sheet_name="TANDA ENDOMETRIO SANO"))
        df = df.rename(columns={"Código Noray": "sample_id"})
        df["sample_id"] = df["sample_id"].map(_to_str)
        tables["clinical_endo_healthy"] = df.rename(columns={c: f"clin_en__{c}" for c in df.columns if c != "sample_id"})

    # HER2 case-level clinical: join by case_id
    her2 = None
    if "TANDA HER2 MAMA T + NT" in xls.sheet_names:
        her2 = standard_columns(pd.read_excel(xls, sheet_name="TANDA HER2 MAMA T + NT"))
        her2 = her2.rename(columns={
            "CÓDIGO BB CASO": "case_id",
            "NHC": "NHC",
            "CÓDIGO BB BLOQUE TUMORAL": "tumour_block_id",
            "CÓDIGO BB BLOQUE NORMAL": "paired_normal_block_id",
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
    Loads AU endometrial tumour clinical and prefixes with clin_au_endo__.
    Join key is 'endo_id', derived from AU PATIENT_ID (and from master sample_id).
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
        merged = merged.merge(au, how="left", on="endo_id")

    return merged


# =============================================================================
# SEMANTIC NORMALISATION HELPERS (Option B)
# =============================================================================

def _norm_text(x) -> str:
    """Normalise text (remove NBSP, trim)."""
    s = _to_str(x)
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
    Accepts: SI/Sí/No, Yes/No, Y/N, 1/0, True/False, Pos/Neg.
    """
    if pd.isna(x):
        return None
    s = _norm_text(x).lower()
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
    if s in {"si", "sí", "s"}:
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


def extract_first_number(x) -> Optional[float]:
    """
    Extract first numeric token from a messy cell.
    Examples:
      "12 años" -> 12
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
    try:
        dt = pd.to_datetime(x, errors="coerce", dayfirst=True)
        return None if pd.isna(dt) else dt
    except Exception:
        return None


def coalesce_with_source(df: pd.DataFrame, cols: List[str]) -> Tuple[pd.Series, pd.Series]:
    """
    Coalesce: first non-null across candidate columns, returning both value and which source column was used.
    """
    existing = [c for c in cols if c in df.columns]

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
    - canon__ki67         : new — coalesces FFPE_KI67 (endo) + clin_dcs__KI67 (breast)
    - canon__p53_ihc      : new — coalesces TP53 IHC (endo) + p53 IHC (breast)
    - canon__date_diagnosis: new — diagnosis date (breast)
    - canon__date_last_fu  : new — last follow-up date (breast), used to derive OS
    - canon__date_recurrence: new — recurrence date (breast)
    - canon__recurrence    : new — recurrence binary flag (breast + endo)
    - canon__distant_mets  : new — distant metastasis flag (breast)
    - canon__lymph_nodes   : new — lymph node ratio string (breast) / N field (endo)
    - canon__myometrial_invasion: new — myometrial infiltration (endo)
    - canon__lvsi          : new — lymphovascular space invasion (endo)
    - canon__risk_group    : new — risk-of-recurrence group (endo)
    - canon__treatment     : new — treatment description (breast)
    - canon__vital_status  : new — verbose vital status (breast)
    - canon__ptnm          : already present but now also picks up clin_dcs__pTNM
    """
    return {
        # ── Demographics ─────────────────────────────────────────────────────
        "canon__age": [
            "clin_au_endo__AGE",
            "clin_dcs__Edad dx",
            "clin_mn__Edad muestra",
            "clin_en__Edad muestra",
        ],
        "canon__bmi": [
            "clin_au_endo__BMI_CALCULATED",
            "clin_au_endo__BMI",
            "clin_dcs__BMI",        # ← v4: breast tumour BMI (was omitted before)
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
            "clin_mn__AGO (GESTACIONES, PARTOS, CESÁREAS, ABORTOS)",
            "clin_en__AGO (GESTACIONES, PARTOS, CESÁREAS, ABORTOS)",
            "clin_her2__AGO (GESTACIONES, PARTOS, CESÁREAS, ABORTOS)",
        ],

        # ── Diagnosis / histology / tumour ───────────────────────────────────
        "canon__diagnosis": [
            "clin_dcs__Dx",             # ← v4: prefer DCS Dx (cleaner, more complete)
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
            "clin_dcs__GRADO",          # ← v4: breast tumour grade (was omitted before)
        ],
        "canon__figo_stage": [
            "clin_au_endo__FIGO_STAGE",
            "clin_au_endo__FIGO_STAGE_GROUP",
        ],
        "canon__ptnm": [
            "clin_dcs__pTNM",           # ← v4: breast pTNM
            "clin_au_endo__T",          #   (endo T/N/M kept in separate cols below)
        ],

        # ── Pathological staging extras (endometrial) ────────────────────────
        "canon__myometrial_invasion": [  # ← v4: new
            "clin_au_endo__MYOMETRIAL_INFILTRATION",
        ],
        "canon__lvsi": [                 # ← v4: new
            "clin_au_endo__LVSI",
        ],
        "canon__risk_group": [           # ← v4: new
            "clin_au_endo__RISK_OF_RECURRENCE",
        ],
        "canon__lymph_nodes": [          # ← v4: new — node ratio string (breast) or N (endo)
            "clin_dcs__Numero_ganglio",
            "clin_au_endo__N",
        ],

        # ── Biomarkers ───────────────────────────────────────────────────────
        "canon__ki67": [                 # ← v4: new unified KI67
            "clin_au_endo__FFPE_KI67",
            "clin_dcs__KI67",
        ],
        "canon__p53_ihc": [              # ← v4: new unified p53/TP53 IHC
            "clin_au_endo__FFPE_TP53_IHC",
            "clin_au_endo__FFPE_TP53_IHC POLAND RESULT",
            "clin_dcs__p53",
        ],

        # ── ER / PR status ───────────────────────────────────────────────────
        "canon__er_raw": [
            "clin_au_endo__FFPE_ER1_RECEPTORS_RAW_VALUE",
            "clin_au_endo__FFPE_ER1_RECEPTORS_RAW_VALUE POLAND RESULTS",
        ],
        "canon__pr_raw": [
            "clin_au_endo__FFPE_PR_RECEPTORS_RAW_VALUE",
            "clin_au_endo__FFPE_PR_RECEPTORS_RAW_VALUE POLAND RESULTS",
        ],
        "canon__er_status_breast": [
            "clin_dcs__RE",
        ],
        "canon__pr_status_breast": [
            "clin_dcs__RP",
        ],

        # ── HER2 ─────────────────────────────────────────────────────────────
        "canon__her2_copies": [
            "clin_her2__COPIAS HER2",
            "clin_dcs__COPIAS HER2",
        ],
        "canon__her2_copies_note": [
            "clin_her2__COPIAS HER2",
            "clin_dcs__COPIAS HER2",
        ],

        # ── Survival / outcomes ──────────────────────────────────────────────
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
            "clin_dcs__Exitus",         # ← v4: breast exitus
        ],
        "canon__vital_status": [         # ← v4: new — verbose status string
            "clin_dcs__STATUS",
        ],
        "canon__recurrence": [           # ← v4: new — any recurrence/progression
            "clin_dcs__Recaida/Progresión",
        ],
        "canon__distant_mets": [         # ← v4: new — distant metastasis
            "clin_dcs__MTxDISTANCIA",
        ],
        "canon__treatment": [            # ← v4: new — treatment description
            "clin_dcs__Tto",
        ],

        # ── Dates ────────────────────────────────────────────────────────────
        "canon__date_diagnosis": [
            "clin_dcs__Fecha dx",
            "clin_au_endo__DATE_OF_SURGERY",
        ],
        "canon__date_last_fu": [         # ← v4: new — last follow-up (used to derive OS)
            "clin_dcs__Última fecha disponible",
            "clin_au_endo__LAST_UPDATED",
        ],
        "canon__date_recurrence": [      # ← v4: new
            "clin_dcs__Fecha Recidiva/Progresión",
            "clin_au_endo__PD_DATE",
        ],

        # ── Molecular (endometrial) ───────────────────────────────────────────
        "canon__msi_status": [
            "clin_au_endo__MSI_STATUS_IHC",
            "clin_au_endo__UA_MSI_STATUS_NGS",
            "clin_au_endo__UA_MSI_STATUS_ddPCR",
        ],
        "canon__molecular_class": [
            "clin_au_endo__MOLECULAR CLASSIFICATION_according to IHC and/or NGS profile",
            "clin_au_endo__POLAND MOLECULAR CLASSIFICATION",
        ],

        # ── Operational ──────────────────────────────────────────────────────
        "canon__performed_or_cuts": [
            "clin_mn__REALIZADO",
            "clin_en__CORTES",
            "clin_her2__REALIZADO",
            "clin_her2__REALIZADO.1",
            "clin_her2__REPETICION DE CORTES",
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
    - canon__ki67        : new — unify endo (fraction 0-1) and breast (% or range string)
                           into a single canon__ki67_pct (percentage scale 0-100)
    - canon__p53_ihc     : new — unify endo (fraction) and breast (fraction) -> Positive/Negative
    - canon__exitus      : merged from two cohorts — normalise SI/NO -> Yes/No
    - canon__os_months   : v4 also DERIVES breast OS from canon__date_diagnosis +
                           canon__date_last_fu when canon__os_months is null
    - canon__recurrence_flag : new binary derived from canon__recurrence text
    - canon__distant_mets_flag: new binary derived from canon__distant_mets
    - canon__myometrial_invasion: new — normalise <50% / >50%
    - canon__lvsi        : new — normalise YES/NO
    - canon__risk_group  : new — normalise LOW/INTERMEDIATE/INTERMEDIATE-HIGH/HIGH
    - canon__vital_status: new — normalise verbose string -> Alive / Dead / Lost / Other
    - canon__performed_or_cuts__yesno : fix — was always null due to logic error; now fixed
    - MSI: 'UNSTABLE' (uppercase raw) now correctly normalised to 'Unstable'
    """
    out = df.copy()

    # ── Strip whitespace from key categoricals ────────────────────────────────
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

    # ── MSI status -> 'Stable' / 'Unstable' ──────────────────────────────────
    if "canon__msi_status" in out.columns:
        def _norm_msi(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if s in ("MSS", "STABLE", "MS-STABLE"): return "Stable"
            if s in ("MSI", "MSI-H", "UNSTABLE", "MSI-HIGH", "MSIH"): return "Unstable"
            return x
        out["canon__msi_status"] = out["canon__msi_status"].map(_norm_msi)

    # ── Molecular class -> POLE / MMRd / NSMP / P53 ──────────────────────────
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

    # ── Age ───────────────────────────────────────────────────────────────────
    if "canon__age" in out.columns:
        out["canon__age"] = pd.to_numeric(out["canon__age"], errors="coerce")

    # ── BMI: strip non-numeric junk before coercion ───────────────────────────
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

    # ── Dates ─────────────────────────────────────────────────────────────────
    for _dcol in ("canon__date_birth", "canon__date_diagnosis",
                  "canon__date_last_fu", "canon__date_recurrence"):
        if _dcol in out.columns:
            out[_dcol] = out[_dcol].apply(parse_date)

    # ── Derive breast OS from diagnosis + last follow-up dates ────────────────
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

    # ── Numeric coercion for survival endpoints ───────────────────────────────
    for _surv_col in ("canon__os_months", "canon__pfs_months"):
        if _surv_col in out.columns:
            out[_surv_col] = pd.to_numeric(out[_surv_col], errors="coerce")

    # ── PD flag (endometrial) ─────────────────────────────────────────────────
    if "canon__pd_status" in out.columns:
        out["canon__pd_flag"] = out["canon__pd_status"].map(
            lambda x: 1 if str(x).strip().upper() == "PD" else (0 if str(x).strip().upper() == "NO PD" else pd.NA)
            if pd.notna(x) else pd.NA
        )

    # ── Exitus: normalise SI/NO + YES/NO -> Yes/No ────────────────────────────
    if "canon__exitus" in out.columns:
        def _norm_exitus(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if s in ("SI", "YES", "EXITUS", "1"): return "Yes"
            if s in ("NO", "0"): return "No"
            # Breast: 'SI: fallo hepático...' -> Yes
            if s.startswith("SI:"): return "Yes"
            return x
        out["canon__exitus"] = out["canon__exitus"].map(_norm_exitus)
        # Derived binary (1/0)
        out["canon__exitus_flag"] = out["canon__exitus"].map(
            lambda x: pd.NA if pd.isna(x) else (1 if x == "Yes" else (0 if x == "No" else pd.NA))
        )

    # ── Vital status: normalise verbose string -> clean categories ────────────
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

    # ── Recurrence flag ───────────────────────────────────────────────────────
    if "canon__recurrence" in out.columns:
        def _norm_recurrence(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if s == "NO": return "No"
            if s in ("PD", "YES"): return "Yes"
            if len(s) > 3: return "Yes"   # any non-trivial text = recurrence described
            return x
        out["canon__recurrence"] = out["canon__recurrence"].map(_norm_recurrence)
        out["canon__recurrence_flag"] = out["canon__recurrence"].map(
            lambda x: pd.NA if pd.isna(x) else (1 if x == "Yes" else (0 if x == "No" else pd.NA))
        )

    # ── Distant metastasis flag ───────────────────────────────────────────────
    if "canon__distant_mets" in out.columns:
        def _norm_mets(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if s == "SI": return "Yes"
            if s in ("NO", "NO-LOCAL"): return "No"
            return x
        out["canon__distant_mets"] = out["canon__distant_mets"].map(_norm_mets)
        out["canon__distant_mets_flag"] = out["canon__distant_mets"].map(
            lambda x: pd.NA if pd.isna(x) else (1 if x == "Yes" else (0 if x == "No" else pd.NA))
        )

    # ── Myometrial invasion -> <50% / >50% / None ────────────────────────────
    if "canon__myometrial_invasion" in out.columns:
        def _norm_myo(x):
            if pd.isna(x): return x
            s = str(x).strip().upper()
            if "<50" in s: return "<50%"
            if ">50" in s: return ">50%"
            return x
        out["canon__myometrial_invasion"] = out["canon__myometrial_invasion"].map(_norm_myo)

    # ── LVSI -> Yes / No ──────────────────────────────────────────────────────
    if "canon__lvsi" in out.columns:
        out["canon__lvsi"] = out["canon__lvsi"].map(
            lambda x: "Yes" if str(x).strip().upper() == "YES"
            else ("No" if str(x).strip().upper() == "NO" else (pd.NA if pd.isna(x) else x))
        )

    # ── Risk group: normalise to LOW / INTERMEDIATE / INTERMEDIATE-HIGH / HIGH ─
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

    # ── KI67: unify into a single percentage-scale column ────────────────────
    # Endo source (FFPE_KI67): stored as percentage (0-100)
    # Breast DCS source (clin_dcs__KI67): stored as fraction (0.0-1.0) for most rows,
    #   but some are range strings like '20-25%', '>50%', '70-75%'
    if "canon__ki67" in out.columns:
        out["canon__ki67_raw"] = out["canon__ki67"]
        out["canon__ki67_pct"] = out.apply(
            lambda row: parse_ki67_percent(row["canon__ki67"], row.get("canon__ki67__source")),
            axis=1,
        )

    # ── P53 IHC -> Positive / Negative (unified) ─────────────────────────────
    # Endo: stored as fraction (FFPE_TP53_IHC); breast: fraction (clin_dcs__p53)
    # Convention: >= 0.1 (10%) = aberrant/positive p53
    if "canon__p53_ihc" in out.columns:
        out["canon__p53_ihc_raw"] = out["canon__p53_ihc"]
        p53_parsed = out["canon__p53_ihc"].apply(parse_p53_ihc_value)
        out["canon__p53_ihc_numeric"] = p53_parsed.map(lambda t: t[0])
        out["canon__p53_status"] = p53_parsed.map(lambda t: t[1] if t[1] is not None else pd.NA)

    # ── Menarche ──────────────────────────────────────────────────────────────
    if "canon__menarche_age" in out.columns:
        out["canon__menarche_age"] = out["canon__menarche_age"].apply(extract_first_number)

    # ── Menopause ─────────────────────────────────────────────────────────────
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

    # ── Grade: numeric coercion (covers both endo text 'G1' and breast integer) ─
    if "canon__grade" in out.columns:
        out["canon__grade"] = out["canon__grade"].apply(parse_grade)

    # ── FIGO stage ────────────────────────────────────────────────────────────
    if "canon__figo_stage" in out.columns:
        out["canon__figo_stage"] = out["canon__figo_stage"].apply(parse_figo_stage)

    # ── AU endometrial ER/PR: numeric raw -> Positive/Negative (threshold >= 1%) ─
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

    # ── SP breast ER/PR: POSITIVO/NEGATIVO -> Positive/Negative ──────────────
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

    # ── HER2 copies ───────────────────────────────────────────────────────────
    if "canon__her2_copies" in out.columns:
        out["canon__her2_copies"] = out["canon__her2_copies"].apply(extract_first_number)

    # ── HER2 note + triple negative flag ─────────────────────────────────────
    if "canon__her2_copies_note" in out.columns:
        note = out["canon__her2_copies_note"].copy()
        note_num = note.apply(extract_first_number)
        out["canon__her2_copies_note"] = note.where(note_num.isna(), pd.NA)
        out["canon__triple_negative_flag"] = out["canon__her2_copies_note"].astype(str).str.upper().str.contains(
            "TRIPLE NEGATIVO", na=False
        ).map({True: "Yes", False: "No"})
        out.loc[out["canon__her2_copies_note"].isna(), "canon__triple_negative_flag"] = pd.NA

    # ── Performed/cuts: fix the always-null bug (was checking 'nan' string) ───
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

    # ── Drop columns that are >99% empty ──────────────────────────────────────
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
        print(f"[INFO] Dropping {len(cols_to_drop)} columns that are ≥99% empty:")
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

DEFAULT_MANIFESTS_DIR = Path("/home/gadeaalonsoj/tfm/manifests")

# Maps manifest filename pattern -> (snp_code regex, tissue, tumour_normal)
MANIFEST_PATTERNS = {
    "endometrium-tumour-pass_manifest.txt": (r"SNP_AT_(\d+)", "AT"),
    "endometrium-normal-pass_manifest.txt": (r"SNP_EN_(\d+)", "EN"),
    "breast-tumour-pass_manifest.txt":      (r"MT.T_(\d+)",   "MT"),
    "breast-normal-pass_manifest.txt":      (r"MN_(\d+)",     "MN"),
}


def load_sequenced_snp_codes(manifests_dir: Path) -> set:
    """
    Reads all pass_manifest.txt files from the manifests directory and returns
    a set of normalised snp_codes that have a passing BAM.

    The BAM filename encodes the SNP code (e.g. SNP_AT_10, DNA_SNP_MT-T_5, SNP_DNA_MN_3).
    We extract the numeric part and reconstruct the canonical snp_code as it
    appears in the master (SNP_AT_10, SNP_EN_3, SNP_MT-T_5, SNP_MN_3).
    """
    sequenced: set = set()

    if not manifests_dir.exists():
        print(f"[WARN] Manifests directory not found: {manifests_dir} — sequencing flags will be skipped.")
        return sequenced

    for fname, (num_pattern, prefix) in MANIFEST_PATTERNS.items():
        fpath = manifests_dir / fname
        if not fpath.exists():
            print(f"[WARN] Manifest not found: {fpath}")
            continue

        with open(fpath) as f:
            bams = [line.strip() for line in f if line.strip()]

        for bam in bams:
            bam_name = Path(bam).stem  # filename without .bam
            # Skip replicates — they have the same patient, don't double-count
            if "repeticion" in bam_name.lower():
                # Still add the base sample number
                pass
            m = re.search(num_pattern, bam_name, re.IGNORECASE)
            if m:
                num = int(m.group(1))
                # Reconstruct canonical snp_code as it appears in master
                if prefix == "AT":
                    code = f"SNP_AT_{num}"
                elif prefix == "EN":
                    code = f"SNP_EN_{num}"
                elif prefix == "MT":
                    code = f"SNP_MT-T_{num}"
                elif prefix == "MN":
                    code = f"SNP_MN_{num}"
                else:
                    code = bam_name
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

        # Apply green row highlighting to harmonised sheet
        ws_harm = xw.sheets["harmonised_plus_canon"]
        apply_sequenced_highlighting(ws_harm, harm_df)

        # Also highlight raw_merged sheet
        ws_raw = xw.sheets["raw_merged"]
        apply_sequenced_highlighting(ws_raw, merged_df)

    print(f"Saved: {out_xlsx}")
    print(f"Green rows = sequenced samples with passing BAMs ({n_seq} rows highlighted)")


if __name__ == "__main__":
    main()