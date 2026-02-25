"""
MASTER SNP BUILDER + CLINICAL JOINER (ONE MERGED FILE)
+ AU ENDOMETRIAL TUMOUR CLINICAL JOIN (3rd Excel)

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

3) AU Endometrial tumour clinical workbook:
   "Clinical DATA AU_Endometrial cancer.xlsx"
   Sheet used: Hoja1
   Key used: PATIENT_ID (NOT "nombre en el chip"/SNP_AT_* because it is NOT stable across files)

Output:
- One merged dataframe (and Excel) containing:
  master SNP rows (DNA/RNA) + all clinical variables that can be linked.

Key design:
- master has BOTH:
    sample_id  (e.g., M06152-T, NT250072, BL039, EndoBL039, MDA040, EndoMDA_40)
    case_id    (e.g., M06152 for M06152-T / M06152-NT; else equals sample_id for NT/BL/MDA codes)
- For HER2 cases, clinical join uses case_id (M06152), not block type.
- For AU endometrial tumour clinical, join uses a normalised "endo_id" derived from PATIENT_ID and sample_id.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, List, Dict

import pandas as pd


# =============================================================================
# CONFIG
# =============================================================================

SNP_XLSX = Path("/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/Muestras SNPs nomenclaturas, equivalencias y cuantificaciones.xlsx")
CLINICAL_XLSX = Path("/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/MUESTRAS + VARIABLES CLINICAS PROY SNPs GM(4).xlsx")


# AU endometrial tumour clinical data (3rd Excel)
AU_ENDO_XLSX = Path("/mnt/c/Users/gadab/OneDrive - Uppsala universitet/Documents/TFM/Docs/Clinical DATA AU_Endometrial cancer.xlsx")
AU_ENDO_SHEET = "Hoja1"

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
    raise KeyError(
        f"None of these columns were found: {candidates}. "
        f"Available: {df.columns.tolist()}"
    )


def add_case_id_and_block_type(master: pd.DataFrame) -> pd.DataFrame:
    """
    Adds:
      - case_id: used for HER2 breast tumour/paired-normal joins (germline logic)
      - block_type: 'T' / 'NT' if sample_id looks like M06152-T or M06152-NT, else None
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


def normalize_endo_id(x: str) -> str:
    """
    Canonical key to match AU endometrial tumour clinical (PATIENT_ID)
    with SNP master sample_id values that can look like:
      BL039
      BL136 1:100
      EndoBL010
      MDA040
      EndoMDA_040
      EndoMDA_40
    Output always one of:
      ENDOBL###      (3 digits)
      ENDOMDA_##     (no zero padding; matches AU style EndoMDA_40)
    """
    s = _to_str(x)
    if not s:
        return ""

    # Remove dilution / trailing comments: "BL136 1:100" -> "BL136"
    s = s.split()[0].strip()

    # Upper + remove spaces
    s = s.replace(" ", "").upper()

    # BL
    m = re.match(r"^ENDOBL(\d+)$", s)
    if m:
        return f"ENDOBL{int(m.group(1)):03d}"
    m = re.match(r"^BL(\d+)$", s)
    if m:
        return f"ENDOBL{int(m.group(1)):03d}"

    # MDA (AU uses EndoMDA_40 style; keep underscore, no padding)
    m = re.match(r"^ENDOMDA[_\-]?(\d+)$", s)
    if m:
        return f"ENDOMDA_{int(m.group(1))}"
    m = re.match(r"^MDA(\d+)$", s)
    if m:
        return f"ENDOMDA_{int(m.group(1))}"

    return s


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
    out["tumour_normal"] = "Tumour"  # still labelled tumour, but germline joins use case_id
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
    """Build SNP master (DNA+RNA long format) and add case_id + block_type."""
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

    # Add case_id + block_type to support germline joins (HER2)
    master = add_case_id_and_block_type(master)

    # Column order
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

def load_clinical_tables(clinical_xlsx: Path) -> Dict[str, pd.DataFrame]:
    """
    Load and standardize the clinical workbook into tables:

    1) clinical_breast_healthy (TANDA MAMA SANA) keyed by sample_id (NT250xxx)
    2) clinical_endo_healthy   (TANDA ENDOMETRIO SANO) keyed by sample_id (NT... / MO...)
    3) clinical_her2_case      (TANDA HER2 MAMA T + NT enriched with DCs) keyed by case_id (M06152)
    """
    xls = pd.ExcelFile(clinical_xlsx)
    tables: Dict[str, pd.DataFrame] = {}

    # ---- Healthy breast (NT250xxx) ----
    if "TANDA MAMA SANA" in xls.sheet_names:
        df = standard_columns(pd.read_excel(xls, sheet_name="TANDA MAMA SANA"))
        df = df.rename(columns={"Código Noray-BB": "sample_id"})
        df["sample_id"] = df["sample_id"].map(_to_str)

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

        rename_map = {}
        for c in her2.columns:
            if c in {"case_id", "NHC"}:
                continue
            rename_map[c] = f"clin_her2__{c}"
        her2 = her2.rename(columns=rename_map)

    # Enrich HER2 with DCs
    if her2 is not None and "DCs MAMA HER2" in xls.sheet_names:
        dcs = standard_columns(pd.read_excel(xls, sheet_name="DCs MAMA HER2"))
        dcs["NHC"] = dcs["NHC"].map(_to_str)

        rename_map = {c: f"clin_dcs__{c}" for c in dcs.columns if c != "NHC"}
        dcs = dcs.rename(columns=rename_map)

        her2 = her2.merge(dcs, how="left", on="NHC")

    if her2 is not None:
        tables["clinical_her2_case"] = her2

    return tables


def load_au_endo_clinical(au_xlsx: Path, sheet_name: str = "Hoja1") -> Optional[pd.DataFrame]:
    """
    Load AU endometrial tumour clinical file and return a table keyed by endo_id.
    endo_id is derived from PATIENT_ID using normalize_endo_id().
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

    # Drop empty keys
    au = au[au["endo_id"].astype(str).str.len() > 0].copy()

    # Deduplicate if needed
    if au["endo_id"].duplicated().any():
        dups = au.loc[au["endo_id"].duplicated(), "endo_id"].unique().tolist()
        print(f"WARNING: AU endo clinical has duplicate PATIENT_ID keys after normalisation: {dups[:10]}{'...' if len(dups) > 10 else ''}")
        au = au.drop_duplicates(subset=["endo_id"], keep="first").copy()

    # Prefix columns (except endo_id)
    rename_map = {c: f"clin_au_endo__{c}" for c in au.columns if c != "endo_id"}
    au = au.rename(columns=rename_map)

    # Keep endo_id first
    cols = ["endo_id"] + [c for c in au.columns if c != "endo_id"]
    return au[cols].copy()


def merge_master_with_clinical(master_df: pd.DataFrame, clinical_xlsx: Path, au_endo_xlsx: Path) -> pd.DataFrame:
    """
    Produce ONE merged file:
    - Healthy breast: join by sample_id
    - Healthy endometrium: join by sample_id
    - HER2 breast: join by case_id (germline logic)
    - AU endometrial tumour clinical: join by normalised endo_id (PATIENT_ID <-> sample_id)
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

    # Join AU endometrial tumour clinical by normalised endo_id
    au = load_au_endo_clinical(au_endo_xlsx, sheet_name=AU_ENDO_SHEET)
    if au is not None:
        merged = merged.copy()
        merged["endo_id"] = merged["sample_id"].astype(str).map(normalize_endo_id)

        # Audit before merge
        au_keys = set(au["endo_id"].tolist())
        master_keys = set([k for k in merged["endo_id"].tolist() if k])
        matched = merged["endo_id"].isin(au_keys).sum()
        print("\nAU ENDO JOIN AUDIT")
        print("------------------")
        print(f"AU keys (unique):        {len(au_keys)}")
        print(f"Master endo_id (unique): {len(master_keys)}")
        print(f"Master rows matched:     {matched} / {len(merged)}")

        merged = merged.merge(au, how="left", on="endo_id")

    return merged


# =============================================================================
# RUN (build + merge + save)
# =============================================================================

if __name__ == "__main__":
    # 1) Build SNP master
    master_df = build_master(SNP_XLSX)
    print("MASTER:", master_df.shape)
    print(master_df.head(10).to_string(index=False))

    # 2) Merge all clinical into master (ONE merged file)
    merged_df = merge_master_with_clinical(master_df, CLINICAL_XLSX, AU_ENDO_XLSX)
    print("\nMERGED:", merged_df.shape)
    print(merged_df.head(10).to_string(index=False))

    # 3) Save merged output
    out_path = SNP_XLSX.with_name("MASTER_SNP_plus_clinical__MERGED.xlsx")
    merged_df.to_excel(out_path, index=False)
    print("\nSaved:", out_path)