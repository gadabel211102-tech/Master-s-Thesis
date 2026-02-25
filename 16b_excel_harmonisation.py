#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
16b_excel_harmonisation.py  (Option B: lossless merge + semantic harmonisation)

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
python3 16b_excel_harmonisation.py \
  --snp_xlsx "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/Muestras SNPs nomenclaturas, equivalencias y cuantificaciones.xlsx" \
  --clinical_xlsx "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/MUESTRAS + VARIABLES CLINICAS PROY SNPs GM(4).xlsx" \
  --au_xlsx "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/Clinical DATA AU_Endometrial cancer.xlsx" \
  --out_xlsx "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx"

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
DEFAULT_OUT_XLSX = Path(
    "/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/"
    "MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx"
)

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
        mask = out.isna() & df[c].notna()
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
    """
    return {
        # Demographics / reproductive history
        "canon__age": [
            "clin_au_endo__AGE",
        ],
        "canon__bmi": [
            "clin_au_endo__BMI_CALCULATED",
            "clin_au_endo__BMI",
            "clin_her2__BMI",
            "clin_mn__BMI",
            "clin_en__BMI",
        ],
        "canon__date_birth": [
            "clin_mn__Fecha nacimiento",
            "clin_en__Fecha nacimiento",
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

        # Diagnosis / histology / tumour
        "canon__diagnosis": [
            "clin_her2__DX",
            "clin_au_endo__DIAGNOSIS",
        ],
        "canon__histology": [
            "clin_au_endo__HISTOLOGY_GROUP",
            "clin_au_endo__HISTOLOGY",
            "histology",  # from SNP AT sheet (if present)
        ],
        "canon__grade": [
            "clin_au_endo__GRADE",
        ],
        "canon__figo_stage": [
            "clin_au_endo__FIGO_STAGE",
            "clin_au_endo__FIGO_STAGE_GROUP",
        ],

        # AU receptors (raw numeric) -> derived status later
        "canon__er_raw": [
            "clin_au_endo__FFPE_ER1_RECEPTORS_RAW_VALUE",
            "clin_au_endo__FFPE_ER1_RECEPTORS_RAW_VALUE POLAND RESULTS",
        ],
        "canon__pr_raw": [
            "clin_au_endo__FFPE_PR_RECEPTORS_RAW_VALUE",
            "clin_au_endo__FFPE_PR_RECEPTORS_RAW_VALUE POLAND RESULTS",
        ],

        # Breast HER2 copies:
        # - numeric values become canon__her2_copies
        # - non-numeric text becomes canon__her2_copies_note + derived triple-negative flag
        "canon__her2_copies": [
            "clin_her2__COPIAS HER2",
            "clin_dcs__COPIAS HER2",
        ],
        "canon__her2_copies_note": [
            "clin_her2__COPIAS HER2",
            "clin_dcs__COPIAS HER2",
        ],

        # Operational field: "REALIZADO" vs "CORTES"
        "canon__performed_or_cuts": [
            "clin_mn__REALIZADO",
            "clin_en__CORTES",
            "clin_her2__REALIZADO",
            "clin_her2__REALIZADO.1",
            "clin_her2__REPETICIÓN DE CORTES",
        ],
    }


# =============================================================================
# SEMANTIC TRANSFORMS (turn canon columns into comparable formats)
# =============================================================================

def apply_semantic_transforms(df: pd.DataFrame) -> pd.DataFrame:
    """
    After coalescing canon columns, standardise types/encodings and add derived variables.
    """
    out = df.copy()

    # --- Numeric fields ---
    if "canon__age" in out.columns:
        out["canon__age"] = pd.to_numeric(out["canon__age"], errors="coerce")

    if "canon__bmi" in out.columns:
        out["canon__bmi"] = pd.to_numeric(out["canon__bmi"], errors="coerce")

    # --- Dates ---
    if "canon__date_birth" in out.columns:
        out["canon__date_birth"] = out["canon__date_birth"].apply(parse_date)

    # --- Menarche: extract numeric age robustly ---
    if "canon__menarche_age" in out.columns:
        out["canon__menarche_age"] = out["canon__menarche_age"].apply(extract_first_number)

    # --- Menopause: keep original and add derived age + status ---
    if "canon__menopause_age_or_status" in out.columns:
        s = out["canon__menopause_age_or_status"]
        out["canon__menopause_age"] = s.apply(extract_first_number)
        out["canon__menopause_status"] = s.apply(lambda x: _norm_text(x) if _norm_text(x) else pd.NA)

    # --- Grade ---
    if "canon__grade" in out.columns:
        out["canon__grade"] = out["canon__grade"].apply(parse_grade)

    # --- FIGO stage ---
    if "canon__figo_stage" in out.columns:
        out["canon__figo_stage"] = out["canon__figo_stage"].apply(parse_figo_stage)

    # --- AU receptors: numeric + derived status ---
    # Threshold: >= 1 -> Positive (adjust if your project uses a different cutoff)
    if "canon__er_raw" in out.columns:
        out["canon__er_raw"] = pd.to_numeric(out["canon__er_raw"], errors="coerce")
        out["canon__er_status"] = out["canon__er_raw"].apply(
            lambda v: "Positive" if pd.notna(v) and v >= 1 else ("Negative" if pd.notna(v) else pd.NA)
        )

    if "canon__pr_raw" in out.columns:
        out["canon__pr_raw"] = pd.to_numeric(out["canon__pr_raw"], errors="coerce")
        out["canon__pr_status"] = out["canon__pr_raw"].apply(
            lambda v: "Positive" if pd.notna(v) and v >= 1 else ("Negative" if pd.notna(v) else pd.NA)
        )

    # --- HER2 copies: parse numeric robustly ---
    if "canon__her2_copies" in out.columns:
        out["canon__her2_copies"] = out["canon__her2_copies"].apply(extract_first_number)

    # --- HER2 note: keep only NON-numeric text + derive triple negative flag ---
    if "canon__her2_copies_note" in out.columns:
        note = out["canon__her2_copies_note"].copy()
        note_num = note.apply(extract_first_number)

        # If it contains a number, it's not a "note" -> blank it
        out["canon__her2_copies_note"] = note.where(note_num.isna(), pd.NA)

        # Derive triple negative flag (only meaningful when note exists)
        out["canon__triple_negative_flag"] = out["canon__her2_copies_note"].astype(str).str.upper().str.contains(
            "TRIPLE NEGATIVO", na=False
        ).map({True: "Yes", False: "No"})
        out.loc[out["canon__her2_copies_note"].isna(), "canon__triple_negative_flag"] = pd.NA

    # --- Operational yes/no derived ---
    if "canon__performed_or_cuts" in out.columns:
        out["canon__performed_or_cuts__yesno"] = out["canon__performed_or_cuts"].apply(parse_yes_no)

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

    # Audit table: how many values for each canonical variable
    canon_cols = [c for c in harm.columns if c.startswith("canon__") and not c.endswith("__source")]
    audit = pd.DataFrame({
        "canon_column": canon_cols,
        "non_null_rows": [int(harm[c].notna().sum()) for c in canon_cols],
    }).sort_values("non_null_rows", ascending=False)

    mapping_doc = pd.DataFrame(map_rows).sort_values("canon_column")

    return harm, audit, mapping_doc


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
    return p.parse_args()


def main() -> None:
    args = parse_args()
    snp_xlsx = Path(args.snp_xlsx)
    clinical_xlsx = Path(args.clinical_xlsx)
    au_xlsx = Path(args.au_xlsx)
    out_xlsx = Path(args.out_xlsx)

    # Input checks
    if not snp_xlsx.exists():
        raise FileNotFoundError(f"SNP workbook not found: {snp_xlsx}")
    if not clinical_xlsx.exists():
        raise FileNotFoundError(f"Clinical workbook not found: {clinical_xlsx}")

    if au_xlsx.exists():
        print(f"AU clinical workbook found: {au_xlsx} (sheet={args.au_sheet})")
    else:
        print(f"AU clinical workbook NOT found: {au_xlsx} (will skip AU merge)")

    # 1) Build SNP master
    master_df = build_master(snp_xlsx)
    print(f"MASTER: {master_df.shape}")

    # 2) Merge clinical (lossless)
    merged_df = merge_master_with_clinical(master_df, clinical_xlsx, au_xlsx, args.au_sheet)
    print(f"RAW MERGED: {merged_df.shape}")

    # 3) Harmonise (adds canon variables)
    harm_df, audit_df, map_df = harmonise(merged_df)
    print(f"HARMONISED: {harm_df.shape}")

    # 4) Write output Excel
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        merged_df.to_excel(xw, index=False, sheet_name="raw_merged")
        harm_df.to_excel(xw, index=False, sheet_name="harmonised_plus_canon")
        audit_df.to_excel(xw, index=False, sheet_name="canon_audit")
        map_df.to_excel(xw, index=False, sheet_name="harmonisation_map")

    print(f"Saved: {out_xlsx}")


if __name__ == "__main__":
    main()