"""Validation helpers used to enforce basic data integrity across the pipeline.

These checks are intentionally lightweight: they focus on structural failures
that would undermine interpretation, such as missing files, missing columns,
empty post-filter tables, unexpected tissue labels, or percentage values that
fall outside the mathematically possible 0-100 range.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

from pipeline_utils import format_missing_columns


def validate_file_exists(path: str | Path, context: str) -> None:
    """Raise a clear error if an expected input file is absent."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{context}: file not found -> {p}")


def validate_required_columns(df: pd.DataFrame, required_cols: Iterable[str], context: str) -> None:
    """Ensure that every required analytical column is present in a dataframe."""
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"{context}: missing required columns -> {format_missing_columns(missing)}")


def validate_nonempty(df: pd.DataFrame, context: str) -> None:
    """Guard against silent continuation after an over-restrictive filter."""
    if df.empty:
        raise ValueError(f"{context}: dataframe is empty after filtering")


def validate_tissue_values(df: pd.DataFrame, tissue_col: str, context: str, allowed: tuple[str, ...] = ("Tumour", "Healthy")) -> None:
    """Check that tissue labels belong to the expected controlled vocabulary."""
    observed = set(df[tissue_col].dropna().astype(str).unique())
    unexpected = sorted(observed.difference(allowed))
    if unexpected:
        raise ValueError(f"{context}: unexpected tissue labels -> {', '.join(unexpected)}")


def validate_percentage_columns(df: pd.DataFrame, columns: Iterable[str], context: str) -> None:
    """Ensure that percentage-like outputs remain within their valid range."""
    for col in columns:
        if col not in df.columns:
            continue
        bad = df[col].dropna()
        if ((bad < 0) | (bad > 100)).any():
            raise ValueError(f"{context}: percentage column out of bounds -> {col}")


def print_validation_summary(df: pd.DataFrame, sample_col: str | None, context: str, group_cols: list[str] | None = None) -> None:
    """Print compact run-time metadata that helps document each analysis subset."""
    print(f"[VALIDATION] {context}")
    print(f"  rows: {len(df)}")
    if sample_col and sample_col in df.columns:
        print(f"  unique samples: {df[sample_col].nunique()}")
    if group_cols:
        existing = [col for col in group_cols if col in df.columns]
        if existing:
            print(f"  groups: {df[existing].drop_duplicates().shape[0]}")
