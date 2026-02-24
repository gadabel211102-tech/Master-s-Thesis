"""
MASTER SNP BUILDER + CLINICAL JOINER (ONE MERGED FILE)
Updated to reflect: paired normal (NT) is treated like tumour for germline SNPs.
So for HER2 breast cases we join clinical data by CASE, not by T/NT block.

Inputs:
1) SNP workbook:
   "Muestras SNPs nomenclaturas, equivalencias y cuantificaciones.xlsx"
   Sheets used: AT=AUs, EN, MT-T_N, MN  (OVSER ignored)

2) Clinical workbook:
   "MUESTRAS + VARIABLES CLINICAS PROY SNPs GM(4).xlsx"
   Sheets:
     - TANDA MAMA SANA            (healthy breast; joins via NT250xxx sample code)
     - TANDA ENDOMETRIO SANO      (healthy endometrium; joins via NT2501xx / MO... code)
     - TANDA HER2 MAMA T + NT     (breast tumour + paired normal; joins via CASE ID)
     - DCs MAMA HER2              (detailed clinical; joins via NHC to HER2 sheet)

Output:
- One merged dataframe (and optional Excel) containing:
  master SNP rows (DNA/RNA) + all clinical variables that can be linked.

Key design:
- master has BOTH:
    sample_id  (e.g., M06152-T, NT250072)
    case_id    (e.g., M06152 for M06152-T / M06152-NT; else equals sample_id for NT codes)
- For HER2 cases, clinical join uses case_id (M06152), not block type.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, List

import pandas as pd


# =============================================================================
# CONFIG
# =============================================================================

SNP_XLSX = Path(r"/mnt/data/Muestras SNPs nomenclaturas, equivalencias y cuantificaciones.xlsx")
CLINICAL_XLSX = Path(r"/mnt/data/MUESTRAS + VARIABLES CLINICAS PROY SNPs GM(4).xlsx")

SNP_SHEETS = ["AT=AUs", "EN", "MT-T_N", "MN"]  # ignore OVSER


# =============================================================================
# HELPERS
# =============================================================================

def _to_str(x) -> str:
    return "" if pd.isna(x) else str(x).strip()


def standard_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df


def normalize_snp_code(s: str) -> str:
    return _to_str(s).replace(" ", "")


def parse_volume_to_ul(vol) -> Optional[float]:
    if pd.isna(vol):
        return None
    s = str(vol).strip().lower()
    m = re.search(r"(\d+(?:\.\d+)?)", s)
    return float(m.group(1)) if m else None


def _pick_existing(df: pd.DataFrame, candidates: List[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"None of these columns were found: {candidates}. Available: {df.columns.tolist()}")


def add_case_id_and_block_type(master: pd.DataFrame) -> pd.DataFrame:
    """
    Adds two columns to master:
      - case_id: germline case identifier used for HER2 breast tumour/paired-normal joins
      - block_type: 'T' / 'NT' if sample_id looks like M06152-T or M06152-NT, else None

    Rules:
    - If sample_id matches ^M\\d+-(T|NT)$ -> case_id = M\\d+ and block_type = T or NT
    - Otherwise (e.g., NT250072, MO221478, BL039) -> case_id = sample_id; block_type = None

    This supports your germline logic: T and NT share germline variants, so join by case_id.
    """
    out = master.copy()

    def parse_case_block(sid: str) -> tuple[str, Optional[str]]:
        sid = _to_str(sid)
        m = re.match(r"^(M\d+)-(NT|T)$", sid, flags=re.IGNORECASE)
        if m:
            case_id = m.group(1).upper()
            block_type = m.group(2).upper()
            return case_id, block_type
        return sid, None

    parsed = out["sample_id"].astype(str).apply(parse_case_block)
    out["case_id"] = parsed.apply(lambda x: x[0])
    out["block_type"] = parsed.apply(lambda x: x[1])
    return out


# =============================================================================
# SNP SHEET PROCESSORS
# =============================================================================

def process_AT(sheet_df: pd.DataFrame) -> pd.DataFrame:
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
    df = standard_columns(sheet_df)
    snp_col = "CODIGO SNP_EN (endometrio normal) FFPE"
    vol_col = "VOL (uL)"
    qubit_col = "Qubit (ng(uL)"  # as it appears in your file

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
    out["tumour_normal"] = "Tumour"  # NOTE: still labelled tumour, but germline joins use case_id
    out["sheet"] = "MT-T_N"
    return out


def process_MN(sheet_df: pd.DataFrame) -> pd.DataFrame:
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
    Build the SNP master (DNA+RNA long format) and add case_id + block_type.
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

    # Keep only DNA/RNA rows and non-empty SNP codes
    master["nucleic_acid"] = master["nucleic_acid"].astype(str).str.upper().replace({"": pd.NA, "NAN": pd.NA})
    master = master[master["nucleic_acid"].isin(["DNA", "RNA"])]
    master["snp_code"] = master["snp_code"].astype(str).replace({"": pd.NA, "nan": pd.NA})
    master = master[master["snp_code"].notna()].reset_index(drop=True)

    # Add case_id + block_type to support germline joins
    master = add_case_id_and_block_type(master)

    # Nice column order
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
# CLINICAL LOADER + MERGER (ONE MERGED FILE)
# =============================================================================

def load_clinical_tables(clinical_xlsx: Path) -> dict[str, pd.DataFrame]:
    """
    Load and lightly standardize the clinical workbook into 3 tables:

    1) clinical_breast_healthy (TANDA MAMA SANA) keyed by sample_id (NT250xxx)
    2) clinical_endo_healthy   (TANDA ENDOMETRIO SANO) keyed by sample_id (NT... / MO...)
    3) clinical_her2_case      (TANDA HER2 MAMA T + NT enriched with DCs) keyed by case_id (M06152)

    Returns dict of dataframes.
    """
    xls = pd.ExcelFile(clinical_xlsx)
    tables: dict[str, pd.DataFrame] = {}

    # ---- Healthy breast (NT250xxx) ----
    if "TANDA MAMA SANA" in xls.sheet_names:
        df = standard_columns(pd.read_excel(xls, sheet_name="TANDA MAMA SANA"))
        df = df.rename(columns={"Código Noray-BB": "sample_id"})
        df["sample_id"] = df["sample_id"].map(_to_str)
        # Prefix clinical columns to avoid collisions after merge
        keep = df.copy()
        rename_map = {c: f"clin_mn__{c}" for c in keep.columns if c != "sample_id"}
        keep = keep.rename(columns=rename_map)
        tables["clinical_breast_healthy"] = keep

    # ---- Healthy endometrium (NT.../MO...) ----
    if "TANDA ENDOMETRIO SANO" in xls.sheet_names:
        df = standard_columns(pd.read_excel(xls, sheet_name="TANDA ENDOMETRIO SANO"))
        df = df.rename(columns={"Código Noray": "sample_id"})
        df["sample_id"] = df["sample_id"].map(_to_str)
        keep = df.copy()
        rename_map = {c: f"clin_en__{c}" for c in keep.columns if c != "sample_id"}
        keep = keep.rename(columns=rename_map)
        tables["clinical_endo_healthy"] = keep

    # ---- HER2 breast case-level clinical ----
    # Base: TANDA HER2 MAMA T + NT (has NHC and CÓDIGO BB CASO)
    her2 = None
    if "TANDA HER2 MAMA T + NT" in xls.sheet_names:
        her2 = standard_columns(pd.read_excel(xls, sheet_name="TANDA HER2 MAMA T + NT"))
        # Standardize key columns
        her2 = her2.rename(columns={
            "CÓDIGO BB CASO": "case_id",
            "NHC": "NHC",
            "CÓDIGO BB BLOQUE TUMORAL": "tumour_block_id",
            "CÓDIGO BB BLOQUE NORMAL": "paired_normal_block_id",
        })
        her2["case_id"] = her2["case_id"].map(_to_str).str.upper()
        her2["NHC"] = her2["NHC"].map(_to_str)

        # Prefix non-key columns
        rename_map = {}
        for c in her2.columns:
            if c in {"case_id", "NHC"}:
                continue
            rename_map[c] = f"clin_her2__{c}"
        her2 = her2.rename(columns=rename_map)

    # Enrich: DCs MAMA HER2 (join to her2 via NHC, then keep case_id)
    if her2 is not None and "DCs MAMA HER2" in xls.sheet_names:
        dcs = standard_columns(pd.read_excel(xls, sheet_name="DCs MAMA HER2"))
        # standardize key column
        dcs["NHC"] = dcs["NHC"].map(_to_str)

        # Prefix DCs columns (except NHC)
        rename_map = {c: f"clin_dcs__{c}" for c in dcs.columns if c != "NHC"}
        dcs = dcs.rename(columns=rename_map)

        # Merge DCs into HER2 on NHC (left join keeps all HER2 cases)
        her2 = her2.merge(dcs, how="left", on="NHC")

    if her2 is not None:
        tables["clinical_her2_case"] = her2

    return tables


def merge_master_with_clinical(master_df: pd.DataFrame, clinical_xlsx: Path) -> pd.DataFrame:
    """
    Produce ONE merged file:
    - For NT-coded healthy breast: join by sample_id
    - For NT/MO-coded healthy endometrium: join by sample_id
    - For HER2 breast cases: join by case_id (germline logic), regardless of T vs NT blocks

    Returns merged dataframe.
    """
    tables = load_clinical_tables(clinical_xlsx)
    merged = master_df.copy()

    # Join healthy breast clinical (MN) by sample_id
    if "clinical_breast_healthy" in tables:
        merged = merged.merge(tables["clinical_breast_healthy"], how="left", on="sample_id")

    # Join healthy endometrium clinical (EN) by sample_id
    if "clinical_endo_healthy" in tables:
        merged = merged.merge(tables["clinical_endo_healthy"], how="left", on="sample_id")

    # Join HER2 clinical by case_id (important germline update)
    if "clinical_her2_case" in tables:
        merged = merged.merge(tables["clinical_her2_case"], how="left", on="case_id")

    return merged


# =============================================================================
# RUN (build + merge + save)
# =============================================================================

if __name__ == "__main__":
    # 1) Build SNP master
    master_df = build_master(SNP_XLSX)
    print("MASTER:", master_df.shape)
    print(master_df.head(10).to_string(index=False))

    # 2) Merge clinical into master (ONE merged file)
    merged_df = merge_master_with_clinical(master_df, CLINICAL_XLSX)
    print("\nMERGED:", merged_df.shape)
    print(merged_df.head(10).to_string(index=False))

    # 3) Save merged output (optional)
    out_path = SNP_XLSX.with_name("MASTER_SNP_plus_clinical__MERGED.xlsx")
    merged_df.to_excel(out_path, index=False)
    print("\nSaved:", out_path)