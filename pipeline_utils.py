"""Shared utility helpers used across the TFM analysis pipeline.

The goal of this module is methodological consistency. Repeated tasks such as
loading configuration, harmonising labels, deriving pooled controls, extracting
rsIDs and calculating carrier percentages are defined once here so that every
script applies the same rules.
"""

from __future__ import annotations

import os
import re
import tomllib
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

_CONFIG_PATH = Path(__file__).with_name("pipeline_config.toml")
_CONFIG_TEMPLATE_PATH = Path(__file__).with_name("pipeline_config.example.toml")


def resolve_config_path(config_path: str | Path | None = None) -> Path:
    """Resolve the config file used by the pipeline.

    Priority order:

    1. explicit ``config_path`` argument
    2. ``PIPELINE_CONFIG`` environment variable
    3. repo-local ``pipeline_config.toml``
    4. repo-local ``pipeline_config.example.toml`` as a documented fallback
    """
    if config_path:
        return Path(config_path)

    env_path = os.environ.get("PIPELINE_CONFIG", "").strip()
    if env_path:
        return Path(env_path)
    if _CONFIG_PATH.exists():
        return _CONFIG_PATH
    if _CONFIG_TEMPLATE_PATH.exists():
        return _CONFIG_TEMPLATE_PATH
    raise FileNotFoundError(
        "No pipeline configuration file was found. Expected pipeline_config.toml "
        "or pipeline_config.example.toml beside pipeline_utils.py."
    )


def load_pipeline_config(config_path: str | Path | None = None) -> dict[str, Any]:
    """Load the TOML configuration that defines paths, thresholds and grouping."""
    path = resolve_config_path(config_path)
    with path.open("rb") as handle:
        return tomllib.load(handle)


def get_paths(config_path: str | Path | None = None) -> dict[str, Path]:
    """Return configured filesystem locations as ``Path`` objects.

    Relative paths are resolved against the directory that contains the chosen
    configuration file. This keeps the public config template portable while
    preserving backwards compatibility with absolute local paths.
    """
    cfg = load_pipeline_config(config_path)
    config_dir = resolve_config_path(config_path).resolve().parent
    resolved: dict[str, Path] = {}
    for key, value in cfg["paths"].items():
        path = Path(value)
        resolved[key] = path if path.is_absolute() else (config_dir / path).resolve()
    return resolved


def get_thresholds(config_path: str | Path | None = None) -> dict[str, Any]:
    """Return numerical and logical thresholds used by downstream analyses."""
    return load_pipeline_config(config_path)["thresholds"]


def get_grouping(config_path: str | Path | None = None) -> dict[str, Any]:
    """Return cohort-grouping settings, including pooled-control behaviour."""
    return load_pipeline_config(config_path)["grouping"]


def normalise_chromosome_label(value: Any) -> str:
    """Normalise chromosome labels to a consistent ``chrN`` style."""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    if text.lower().startswith("chr"):
        suffix = text[3:]
    else:
        suffix = text
    suffix_upper = suffix.upper()
    if suffix_upper in {"M", "MT"}:
        return "chrM"
    if suffix_upper in {"X", "Y"}:
        return f"chr{suffix_upper}"
    return f"chr{suffix}"


def get_warning_amplicons(config_path: str | Path | None = None) -> pd.DataFrame:
    """Return the amplicon regions that should trigger cautious SNP interpretation."""
    cfg = load_pipeline_config(config_path)
    warning_rows = cfg.get("warning_amplicons", [])
    if not warning_rows:
        return pd.DataFrame(columns=["label", "chrom", "start", "end", "region"])

    warning_df = pd.DataFrame(warning_rows).copy()
    warning_df["chrom"] = warning_df["chrom"].map(normalise_chromosome_label)
    warning_df["start"] = pd.to_numeric(warning_df["start"], errors="coerce").astype("Int64")
    warning_df["end"] = pd.to_numeric(warning_df["end"], errors="coerce").astype("Int64")
    warning_df = warning_df.dropna(subset=["label", "chrom", "start", "end"]).copy()
    warning_df["start"] = warning_df["start"].astype(int)
    warning_df["end"] = warning_df["end"].astype(int)
    warning_df["region"] = (
        warning_df["chrom"] + ":" + warning_df["start"].astype(str) + "-" + warning_df["end"].astype(str)
    )
    return warning_df[["label", "chrom", "start", "end", "region"]]


def build_amplicon_warning_lookup(
    df: pd.DataFrame,
    variant_col: str = "Variant_ID",
    chrom_col: str = "CHROM",
    pos_col: str = "POS",
    config_path: str | Path | None = None,
) -> pd.DataFrame:
    """Build one warning row per variant for overlaps with flagged amplicons."""
    out_cols = [
        variant_col,
        "Coverage_Risk_Flag",
        "Coverage_Risk_Amplicon",
        "Coverage_Risk_Region",
        "Coverage_Risk_Note",
    ]
    warning_amplicons = get_warning_amplicons(config_path)
    if warning_amplicons.empty:
        return pd.DataFrame(columns=out_cols)
    if any(col not in df.columns for col in [variant_col, chrom_col, pos_col]):
        return pd.DataFrame(columns=out_cols)

    variant_df = df[[variant_col, chrom_col, pos_col]].dropna().copy()
    if variant_df.empty:
        return pd.DataFrame(columns=out_cols)

    variant_df[variant_col] = variant_df[variant_col].astype(str).str.strip()
    variant_df = variant_df[variant_df[variant_col] != ""].copy()
    variant_df[chrom_col] = variant_df[chrom_col].map(normalise_chromosome_label)
    variant_df[pos_col] = pd.to_numeric(variant_df[pos_col], errors="coerce")
    variant_df = variant_df.dropna(subset=[pos_col]).copy()
    if variant_df.empty:
        return pd.DataFrame(columns=out_cols)

    overlaps = variant_df.merge(warning_amplicons, left_on=chrom_col, right_on="chrom", how="inner")
    overlaps = overlaps[(overlaps[pos_col] >= overlaps["start"]) & (overlaps[pos_col] <= overlaps["end"])].copy()
    if overlaps.empty:
        return pd.DataFrame(columns=out_cols)

    def _unique_join(series: pd.Series) -> str:
        values = sorted({str(value).strip() for value in series if str(value).strip()})
        return "; ".join(values)

    lookup = (
        overlaps.groupby(variant_col, as_index=False)
        .agg({
            "label": _unique_join,
            "region": _unique_join,
        })
        .rename(columns={
            "label": "Coverage_Risk_Amplicon",
            "region": "Coverage_Risk_Region",
        })
    )
    lookup["Coverage_Risk_Flag"] = True
    lookup["Coverage_Risk_Note"] = (
        "WARNING: overlaps one of the 4 worst-performing amplicons from the technical audit; "
        "interpret cautiously."
    )
    return lookup[out_cols]


def attach_amplicon_warning_columns(
    df: pd.DataFrame,
    warning_lookup: pd.DataFrame,
    variant_col: str = "Variant_ID",
) -> pd.DataFrame:
    """Attach standard amplicon warning columns to a variant-level result table."""
    if variant_col not in df.columns:
        return df

    out = df.copy()
    warning_cols = [
        "Coverage_Risk_Flag",
        "Coverage_Risk_Amplicon",
        "Coverage_Risk_Region",
        "Coverage_Risk_Note",
    ]
    out = out.drop(columns=[col for col in warning_cols if col in out.columns], errors="ignore")

    if warning_lookup.empty:
        out["Coverage_Risk_Flag"] = False
        out["Coverage_Risk_Amplicon"] = ""
        out["Coverage_Risk_Region"] = ""
        out["Coverage_Risk_Note"] = ""
        return out

    out = out.merge(warning_lookup, on=variant_col, how="left")
    out["Coverage_Risk_Flag"] = out["Coverage_Risk_Flag"].eq(True)
    out["Coverage_Risk_Amplicon"] = out["Coverage_Risk_Amplicon"].fillna("")
    out["Coverage_Risk_Region"] = out["Coverage_Risk_Region"].fillna("")
    out["Coverage_Risk_Note"] = out["Coverage_Risk_Note"].fillna("")
    return out


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



def _phase_call_to_dosage(value: Any) -> float:
    """Convert one phased or unphased diploid genotype call to alt-allele dosage."""
    text = str(value).strip()
    if not text or text.lower() == "nan" or text in {".", "./.", ".|."}:
        return np.nan
    mapping = {
        "0|0": 0.0,
        "0/0": 0.0,
        "0|1": 1.0,
        "1|0": 1.0,
        "0/1": 1.0,
        "1/0": 1.0,
        "1|1": 2.0,
        "1/1": 2.0,
    }
    return mapping.get(text, np.nan)


def _define_ld_blocks_from_matrix(
    ld_matrix: pd.DataFrame,
    threshold: float = 0.80,
    min_snps: int = 2,
    max_snps: int = 12,
) -> list[list[int]]:
    """Define contiguous LD blocks using the same all-pairs coherence rule as stage 15."""
    n_snps = len(ld_matrix.index)
    if n_snps < min_snps:
        return []

    values = ld_matrix.to_numpy(dtype=float)

    def block_coherent(indices: list[int]) -> bool:
        if len(indices) < min_snps:
            return False
        sub = values[np.ix_(indices, indices)]
        lower = sub[np.tril_indices(len(indices), k=-1)]
        return len(lower) > 0 and np.all(np.isfinite(lower) & (lower >= threshold))

    blocks: list[list[int]] = []
    i = 0
    while i < n_snps:
        if i + min_snps > n_snps:
            break
        best_end: int | None = None
        max_j = min(n_snps - 1, i + max_snps - 1)
        for j in range(max_j, i + min_snps - 2, -1):
            indices = list(range(i, j + 1))
            if block_coherent(indices):
                best_end = j
                break
        if best_end is not None:
            blocks.append(list(range(i, best_end + 1)))
            i = best_end + 1
        else:
            i += 1
    return blocks


def load_common_snp_ld_reference(
    annotated_report: str | Path,
    phased_genotypes: str | Path,
    min_nfe_af: float,
    sample_metadata: str | Path | None = None,
    common_snp_whitelist: str | Path | None = None,
    threshold: float = 0.80,
    min_block_snps: int = 2,
    max_block_snps: int = 12,
) -> dict[str, Any]:
    """Compute pairwise LD and contiguous LD blocks for the established common-SNP backbone."""
    annotated_path = Path(annotated_report)
    phased_path = Path(phased_genotypes)
    metadata_path = Path(sample_metadata) if sample_metadata else phased_path.with_name('sample_metadata.tsv')

    annot = pd.read_excel(annotated_path, sheet_name='Biological_Annotations')
    sym_c = find_col(annot, 'SYMBOL')
    var_c = find_col(annot, 'Existing_variation')
    hgv_c = find_col(annot, 'HGVSp')
    exome_c = find_col(annot, 'gnomADe_NFE_AF')
    genome_c = find_col(annot, 'gnomADg_NFE_AF')
    pos_c = find_col(annot, 'POS')
    ref_c = find_col(annot, 'REF')
    alt_c = find_col(annot, 'ALT')
    required = [sym_c, var_c, hgv_c, exome_c, genome_c, pos_c, ref_c, alt_c]
    if any(col is None for col in required):
        raise ValueError('Annotated report is missing one or more columns required to build the LD backbone.')

    annot = combine_gnomad_nfe(annot, exome_c, genome_c)
    annot = annot[annot['gnomAD_NFE_AF_combined'] > float(min_nfe_af)].copy()
    if common_snp_whitelist:
        whitelist_path = Path(common_snp_whitelist)
        if whitelist_path.exists():
            whitelist_df = pd.read_excel(whitelist_path)
            variant_col = find_col(whitelist_df, 'Variant_ID')
            if variant_col is not None:
                whitelist_ids = set(whitelist_df[variant_col].dropna().astype(str).str.strip())
                annot['Variant_ID'] = build_variant_id_series(annot[var_c], annot[sym_c], annot[hgv_c])
                annot = annot[annot['Variant_ID'].astype(str).isin(whitelist_ids)].copy()
    annot['Variant_ID'] = build_variant_id_series(annot[var_c], annot[sym_c], annot[hgv_c])
    annot['POS_int'] = pd.to_numeric(annot[pos_c], errors='coerce').astype('Int64')
    annot['REF_clean'] = annot[ref_c].astype(str).str.strip().str.upper()
    annot['ALT_clean'] = annot[alt_c].astype(str).str.strip().str.upper()
    annot['Gene'] = annot[sym_c].astype(str).str.strip().replace({'': 'Unknown', 'nan': 'Unknown'})
    annot_map = (
        annot[['Variant_ID', 'POS_int', 'REF_clean', 'ALT_clean', 'Gene']]
        .dropna(subset=['Variant_ID', 'POS_int', 'REF_clean', 'ALT_clean'])
        .drop_duplicates(subset=['POS_int', 'REF_clean', 'ALT_clean'], keep='first')
        .sort_values(['POS_int', 'Variant_ID'])
        .reset_index(drop=True)
    )
    if annot_map.empty:
        raise ValueError('No common SNP backbone rows were available for LD computation.')

    phased = pd.read_csv(phased_path, sep='	', dtype=str)
    phased['POS_int'] = pd.to_numeric(phased['POS'], errors='coerce').astype('Int64')
    phased['REF_clean'] = phased['REF'].astype(str).str.strip().str.upper()
    phased['ALT_clean'] = phased['ALT'].astype(str).str.strip().str.upper()

    phased_backbone = (
        phased.merge(annot_map, on=['POS_int', 'REF_clean', 'ALT_clean'], how='inner')
        .sort_values(['POS_int', 'Variant_ID'])
        .drop_duplicates(subset=['Variant_ID'], keep='first')
        .reset_index(drop=True)
    )
    if phased_backbone.empty:
        raise ValueError('No common SNP backbone rows from the annotated report matched the phased genotype table.')

    sample_cols = list(phased.columns[5:])
    if metadata_path.exists():
        meta = pd.read_csv(metadata_path, sep='	', dtype=str)
        if 'Sample' in meta.columns:
            keep_samples = set(meta['Sample'].astype(str).str.strip())
            sample_cols = [col for col in sample_cols if col in keep_samples]
    if not sample_cols:
        raise ValueError('No phased study samples were available for LD computation.')

    dosage = phased_backbone[sample_cols].apply(lambda col: col.map(_phase_call_to_dosage))
    values = dosage.to_numpy(dtype=float)
    snp_ids = phased_backbone['Variant_ID'].astype(str).tolist()
    n_snps = len(snp_ids)
    ld = np.eye(n_snps, dtype=float)
    for i in range(n_snps):
        for j in range(i + 1, n_snps):
            xi = values[i, :]
            xj = values[j, :]
            keep = np.isfinite(xi) & np.isfinite(xj)
            if keep.sum() < 3:
                r2 = np.nan
            else:
                xi_keep = xi[keep]
                xj_keep = xj[keep]
                if np.nanstd(xi_keep) == 0 or np.nanstd(xj_keep) == 0:
                    r2 = np.nan
                else:
                    r = np.corrcoef(xi_keep, xj_keep)[0, 1]
                    r2 = np.nan if np.isnan(r) else float(r * r)
            ld[i, j] = r2
            ld[j, i] = r2

    ld_matrix = pd.DataFrame(ld, index=snp_ids, columns=snp_ids)
    blocks = _define_ld_blocks_from_matrix(
        ld_matrix,
        threshold=float(threshold),
        min_snps=min_block_snps,
        max_snps=max_block_snps,
    )

    threshold_label = f"r2_{int(round(float(threshold) * 100)):03d}"
    block_rows: list[dict[str, Any]] = []
    for block_number, block_indices in enumerate(blocks, start=1):
        block_name = f'Block{block_number}'
        for order, idx in enumerate(block_indices, start=1):
            block_rows.append({
                'LD_Threshold': float(threshold),
                'Threshold_Label': threshold_label,
                'Block': block_name,
                'SNP_Order': order,
                'rsID': phased_backbone.loc[idx, 'Variant_ID'],
                'POS': int(phased_backbone.loc[idx, 'POS_int']),
                'Gene': phased_backbone.loc[idx, 'Gene'],
            })
    block_membership = pd.DataFrame(block_rows)

    ld_long = (
        ld_matrix.stack(dropna=False)
        .rename('R2')
        .reset_index()
        .rename(columns={'level_0': 'SNP1', 'level_1': 'SNP2'})
    )

    missing_variants = sorted(set(annot_map['Variant_ID'].astype(str)) - set(phased_backbone['Variant_ID'].astype(str)))

    return {
        'snp_meta': phased_backbone[['Variant_ID', 'POS_int', 'Gene']].rename(columns={'POS_int': 'POS'}).copy(),
        'ld_matrix': ld_matrix,
        'ld_long': ld_long,
        'block_membership': block_membership,
        'missing_variants': missing_variants,
        'threshold': float(threshold),
        'threshold_label': threshold_label,
        'sample_count': len(sample_cols),
    }

def format_missing_columns(missing: Iterable[str]) -> str:
    """Format missing-column names for readable error messages."""
    return ", ".join(sorted(str(col) for col in missing))
