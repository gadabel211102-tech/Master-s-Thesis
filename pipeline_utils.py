"""Shared utility helpers used across the TFM analysis pipeline.

The goal of this module is methodological consistency. Repeated tasks such as
loading configuration, harmonising labels, deriving pooled controls, extracting
rsIDs and calculating carrier percentages are defined once here so that every
script applies the same rules.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

_CONFIG_PATH = Path(__file__).with_name("pipeline_config.toml")


def load_pipeline_config(config_path: str | Path | None = None) -> dict[str, Any]:
    """Load the TOML configuration that defines paths, thresholds and grouping."""
    path = Path(config_path) if config_path else _CONFIG_PATH
    with path.open("rb") as handle:
        return tomllib.load(handle)


def get_paths(config_path: str | Path | None = None) -> dict[str, Path]:
    """Return configured filesystem locations as ``Path`` objects."""
    cfg = load_pipeline_config(config_path)
    return {key: Path(value) for key, value in cfg["paths"].items()}


def get_thresholds(config_path: str | Path | None = None) -> dict[str, Any]:
    """Return numerical and logical thresholds used by downstream analyses."""
    return load_pipeline_config(config_path)["thresholds"]


def get_grouping(config_path: str | Path | None = None) -> dict[str, Any]:
    """Return cohort-grouping settings, including pooled-control behaviour."""
    return load_pipeline_config(config_path)["grouping"]


def find_col(df: pd.DataFrame, target: str) -> str | None:
    """Locate a dataframe column by case-insensitive exact name matching."""
    target_norm = str(target).strip().upper()
    for col in df.columns:
        if str(col).strip().upper() == target_norm:
            return col
    return None


def standardize_tissue_labels(series: pd.Series) -> pd.Series:
    """Normalise tissue labels to the canonical tumour/healthy vocabulary."""
    return (
        series.astype(str)
        .str.strip()
        .replace({
            "Tumor": "Tumour",
            "Tumour": "Tumour",
            "Normal": "Healthy",
            "Healthy": "Healthy",
            "Control": "Healthy",
        })
    )


def standardize_cohort_labels(series: pd.Series) -> pd.Series:
    """Trim cohort labels without changing their substantive meaning."""
    return series.astype(str).str.strip()


def combine_gnomad_nfe(
    df: pd.DataFrame,
    exome_col: str = "gnomADe_NFE_AF",
    genome_col: str = "gnomADg_NFE_AF",
    out_col: str = "gnomAD_NFE_AF_combined",
    source_col: str = "gnomAD_NFE_Source",
) -> pd.DataFrame:
    """Create one harmonised NFE frequency field plus a provenance column."""
    df = df.copy()
    df[out_col] = df[exome_col].combine_first(df[genome_col])
    df[source_col] = np.where(
        df[exome_col].notna(), "Exome_NFE",
        np.where(df[genome_col].notna(), "Genome_NFE", "Missing")
    )
    return df


def extract_rsid(value: Any) -> str | None:
    """Extract the first rsID from a semi-structured annotation field."""
    text = str(value).strip()
    if not text or text.lower() == "nan" or text == "-":
        return None
    match = re.search(r"(rs\d+)", text, flags=re.IGNORECASE)
    return match.group(1).lower() if match else None


def build_variant_id_series(existing_series: pd.Series, symbol_series: pd.Series, hgvsp_series: pd.Series) -> pd.Series:
    """Build a stable variant identifier, preferring rsIDs over protein notation."""
    rsids = existing_series.astype(str).str.extract(r"(rs\d+)", expand=False)
    fallback = symbol_series.astype(str) + ":" + hgvsp_series.astype(str)
    return rsids.fillna(fallback)


def pooled_control_frame(
    df: pd.DataFrame,
    cohort_col: str,
    tissue_col: str,
    cohort_name: str | None = None,
    cohort_filter: str | None = None,
    tumour_label: str = "Tumour",
    control_label: str = "Healthy",
) -> pd.DataFrame:
    """Return a tumour subset combined with the pooled normal/control arm."""
    controls = df[df[tissue_col] == control_label].copy()
    if cohort_name is None or cohort_filter is None:
        tumours = df[df[tissue_col] == tumour_label].copy()
    else:
        tumours = df[
            (df[tissue_col] == tumour_label)
            & df[cohort_col].str.contains(cohort_filter, case=False, na=False)
        ].copy()
    return pd.concat([tumours, controls], ignore_index=True)


def compute_carrier_percentage(
    df: pd.DataFrame,
    sample_col: str,
    group_cols: list[str],
    denom_cols: list[str],
    count_name: str = "Carrier_Count",
    total_name: str = "Total_Samples",
    pct_name: str = "Carrier_Percentage",
) -> pd.DataFrame:
    """Summarise unique carrier counts and convert them to within-group percentages."""
    totals = df.groupby(denom_cols)[sample_col].nunique().reset_index(name=total_name)
    out = (
        df.groupby(group_cols)[sample_col]
        .nunique()
        .reset_index(name=count_name)
        .merge(totals, on=denom_cols, how="left")
    )
    out[pct_name] = out[count_name] / out[total_name] * 100
    return out


def ensure_directory(path: str | Path) -> Path:
    """Create a directory if needed and return it as a ``Path`` object."""
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def format_missing_columns(missing: Iterable[str]) -> str:
    """Format missing-column names for readable error messages."""
    return ", ".join(sorted(str(col) for col in missing))
