#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
17_snp_association.py  (v4 — SNP-only, genotypic model, all-cohort Cox, age+BMI adjusted)
=======================================================================
SNP–Clinical Variable Association Analysis

Joins the GSDMB variant report (GSDMB_Annotated_Report.xlsx) with the
harmonised clinical master (MASTER_SNP_plus_clinical__HARMONISED_B_v3.xlsx)
and tests each SNP for association with clinical variables.

SAMPLE FILTERING — SEQUENCED SAMPLES ONLY
-------------------------------------------
Analysis is restricted to samples confirmed as sequenced, derived from the
manifest sheets in the harmonised master.  Exclusion rules (verified from data):

  AT=AUs (AU endometrial):
    EXCLUDE extraction_flag == 'NO HAY'   — no biological material available
    EXCLUDE pd_status == 'NO HACER'       — explicitly flagged "do not process"
    → 114 / 118 DNA samples retained

  MT-T_N (breast tumour + paired normal):
    No exclusion flags present → all 71 DNA samples retained

  EN (endometrial healthy):
    No exclusion flags present → all 80 DNA samples retained

  MN (breast healthy):
    No exclusion flags present → all 100 DNA samples retained

The GSDMB variant report (inner join on snp_code) provides a second-level
filter: only snp_codes that produced variant calls will be present in that
file.  Together, these two filters ensure analysis is limited to samples
that have both been (a) approved for sequencing and (b) successfully run.

ANALYSES
--------
1. Tumour vs Control       — Fisher's exact test per cohort + globally
2. SNP × Clinical variable — Fisher's exact (binary/categorical) or
                              Mann–Whitney U (continuous), per tumour cohort
                              Age- and BMI-adjusted logistic regression for binary outcomes
3. Genotype-dose (WT/Het/Hom) — Overall trend + separate Het-vs-WT and
                              Hom-vs-WT contrasts for top SNPs, both cohorts
4. Survival (Cox PH)       ? Age- and BMI-adjusted Cox proportional hazards for
                              endometrial OS and PFS only in the current default
                              clinical set. KM curves are generated per genotype
                              class for endometrial survival endpoints.

VARIABLE COVERAGE (current default focus)
-----------------------------------------
BREAST (MT-T_N sheet, n=71 tumour samples)
  ? Recurrence        canon__recurrence_flag       (binary, n=65)
  ? Metastasis        canon__distant_mets_flag     (binary, n=65)
  ? Death             canon__exitus_flag           (binary, n=65)
  ? Age               canon__age                   (covariate)
  ~ BMI               canon__bmi                   (covariate; partial coverage)
  ? Overall survival  BREAST_OS_MONTHS_DERIVED     (currently empty in this master snapshot)

ENDOMETRIAL (AT=AUs sheet, n=114 sequenced DNA samples)
  ? Progression       canon__pd_flag               (binary, n=114)
  ? Death             canon__exitus_flag           (binary, n=106)
  ? Risk recurrence   canon__risk_group            (labelled groups, n=114)
  ? OS (months)       canon__os_months             (continuous, n=114)
  ? PFS (months)      canon__pfs_months            (continuous, n=114)
  ? Age               canon__age                   (covariate)
  ? BMI               canon__bmi                   (covariate)

OUTPUTS
-------
  17_SNP_Clinical_Association_Results.xlsx
    • tumour_vs_control
    • breast_clinical_assoc
    • endo_clinical_assoc
    • genotype_dose
    • survival_cox          (endometrial survival only in current default set)
    • summary_significant
    • sample_manifest

  17_SNP_Volcano_Plots_raw_p.png
  17_SNP_Volcano_Plots_FDR.png
  17_SNP_Heatmap_Breast_raw_p.png
  17_SNP_Heatmap_Breast_FDR.png
  17_SNP_Heatmap_Endometrial_raw_p.png
  17_SNP_Heatmap_Endometrial_FDR.png
  17_KM_Curves_Endometrial.pdf

RUN
---
python3 17_snp_association.py \\
  --gsdmb   /path/to/GSDMB_Annotated_Report.xlsx \\
  --master  /path/to/MASTER_SNP_plus_clinical_HARMONISED.xlsx \\
  --out_dir /path/to/output/
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import warnings
import textwrap
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.backends.backend_pdf as pdf_backend
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from scipy import stats
from scipy.stats import false_discovery_control, fisher_exact, mannwhitneyu
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning

from association_runtime import script17_defaults
from figure_style import COMPARATIVE_TAG, COHORT_COLORS, GENOTYPE_COLORS, IMPACT_COLORS, arm_color, cohort_color, tagged_title
from pipeline_utils import attach_amplicon_warning_columns, build_amplicon_warning_lookup
from pipeline_validation import print_validation_summary, validate_file_exists, validate_percentage_columns
from sample_identity_utils import attach_analysis_sample_ids

# lifelines optional — only needed for survival plots
try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test, multivariate_logrank_test
    _HAS_LIFELINES = True
except ImportError:
    _HAS_LIFELINES = False
    print("WARNING: lifelines not installed — survival plots will be skipped.")

# Suppress only known noisy third-party warnings; convergence warnings from
# statsmodels and lifelines are intentionally left visible so that failed or
# unreliable model fits are not silently swallowed.
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=".*tight_layout.*")
warnings.filterwarnings("ignore", message=".*More than 20 figures.*")

sns.set_theme(
    context="paper",
    style="whitegrid",
    font="DejaVu Sans",
    rc={
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.edgecolor": "#B9BEC9",
        "axes.labelcolor": "#2A2F36",
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.color": "#414854",
        "ytick.color": "#414854",
        "grid.color": "#D5DBE3",
        "grid.linestyle": ":",
        "grid.linewidth": 0.7,
        "legend.frameon": False,
    },
)

# Shared plotting / reporting helpers restored for thesis figures
_IMPACT_C = {
    "HIGH": "#D55E00",
    "MODERATE": "#E69F00",
    "LOW": "#009E73",
    "MODIFIER": "#7F7F7F",
    "UNKNOWN": "#BDBDBD",
}
_GENO_PLOT_C = {
    "WT": "#BDBDBD",
    "Het": "#E69F00",
    "Hom": "#D55E00",
}
_GENO_PLOT_HATCH = {
    "WT": "",
    "Het": "//",
    "Hom": "xx",
}
_GENO_PLOT_TEXT = {
    "WT": "#2D2D2D",
    "Het": "#202020",
    "Hom": "#FFFFFF",
}
_GENO_PLOT_MARKER = {
    "WT": "o",
    "Het": "s",
    "Hom": "D",
}
_COHORT_C = {
    "Breast": cohort_color("Breast"),
    "Endometrial": cohort_color("Endometrial"),
    "Endometrium": cohort_color("Endometrium"),
    "Global": "#355C7D",
}
_FIG_DPI = 450


def _blend_color(color: str, target: str = "#FFFFFF", weight: float = 0.35) -> str:
    base = np.array(to_rgb(color))
    target_rgb = np.array(to_rgb(target))
    mixed = np.clip((1 - weight) * base + weight * target_rgb, 0, 1)
    return to_hex(mixed)


def _cohort_family(label: str) -> Dict[str, str]:
    base = _COHORT_C.get(label, cohort_color(label, "#4C566A"))
    return {
        "base": base,
        "dark": _blend_color(base, target="#111827", weight=0.18),
        "mid": _blend_color(base, target="#FFFFFF", weight=0.20),
        "light": _blend_color(base, target="#FFFFFF", weight=0.62),
        "very_light": _blend_color(base, target="#FFFFFF", weight=0.84),
    }


def _cohort_heatmap_cmap(label: str, metric_name: Optional[str] = None) -> LinearSegmentedColormap:
    family = _cohort_family(label)
    if str(metric_name).lower() == "fdr":
        colors = [family["very_light"], family["light"], family["mid"], family["dark"]]
    else:
        colors = [family["very_light"], family["light"], family["base"], family["dark"]]
    return LinearSegmentedColormap.from_list(f"{label}_{metric_name or 'metric'}", colors)


def _style_ax(ax, grid=True):
    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#B9BEC9")
    ax.spines["bottom"].set_color("#B9BEC9")
    if grid:
        ax.yaxis.grid(True, linestyle=":", color="#D5DBE3", alpha=0.9, linewidth=0.7)
        ax.set_axisbelow(True)


def _sig_label(p):
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def _metric_prefix(metric_name: Optional[str] = None) -> str:
    return "q" if str(metric_name).lower() == "fdr" else "p"


def _format_probability(value: float, prefix: str = "p") -> str:
    if pd.isna(value):
        return ""
    value = float(value)
    if value < 1e-3:
        return f"{prefix}={value:.1e}"
    if value < 1e-2:
        return f"{prefix}={value:.3f}"
    if value < 0.1:
        return f"{prefix}={value:.2f}"
    return f"{prefix}={value:.2f}"


def _format_sig_annotation(value: float, metric_name: Optional[str] = None, show_ns: bool = False) -> str:
    if pd.isna(value):
        return ""
    sig = _sig_label(float(value))
    if sig == "ns" and not show_ns:
        return ""
    prefix = _metric_prefix(metric_name)
    return f"{sig}\n{_format_probability(float(value), prefix=prefix)}" if sig != "ns" else "ns"


def _save_figure(fig: plt.Figure, out_path: Path, pad_inches: float = 0.24):
    fig.canvas.draw_idle()
    fig.savefig(out_path, dpi=_FIG_DPI, bbox_inches="tight", pad_inches=pad_inches, facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def _safe_tight_layout(fig: plt.Figure, rect=None, pad: float = 1.0, h_pad: float = 1.0, w_pad: float = 1.0):
    try:
        fig.tight_layout(rect=rect, pad=pad, h_pad=h_pad, w_pad=w_pad)
    except Exception:
        pass


def _short_label(value, max_len=34):
    if pd.isna(value):
        return ""
    text = re.sub(r"\s+", " ", str(value)).strip()
    return text if len(text) <= max_len else text[: max_len - 1].rstrip() + "..."


def _wrap_label(value, width=20, max_lines=3):
    text = re.sub(r"\s+", " ", str(value)).strip()
    wrapped = textwrap.wrap(text, width=width, break_long_words=False, break_on_hyphens=False)
    if not wrapped:
        return text
    if len(wrapped) > max_lines:
        tail = " ".join(wrapped[max_lines - 1:])
        wrapped = wrapped[: max_lines - 1] + [_short_label(tail, max_len=width)]
    return "\n".join(wrapped)


def _set_wrapped_ticklabels(ax, axis="x", width=18, max_lines=3, rotation=0, ha="center", fontsize=9):
    labels = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    wrapped = [_wrap_label(label.get_text(), width=width, max_lines=max_lines) for label in labels]
    if axis == "x":
        ax.set_xticklabels(wrapped, rotation=rotation, ha=ha, fontsize=fontsize)
    else:
        ax.set_yticklabels(wrapped, rotation=rotation, ha=ha, fontsize=fontsize)


def _add_margin_labels(ax, points, side="right", max_labels=5, fontsize=7.5):
    if not points:
        return
    points = sorted(points, key=lambda item: item[1], reverse=True)[:max_labels]
    y_min, y_max = ax.get_ylim()
    pad = (y_max - y_min) * 0.06 if y_max != y_min else 0.5
    y_positions = np.linspace(y_max - pad, y_min + pad, len(points))
    x_text = 0.03 if side == "left" else 0.97
    ha = "left" if side == "left" else "right"
    for (x_val, y_val, label), y_text in zip(points, y_positions):
        ax.annotate(
            _short_label(label, 28),
            xy=(x_val, y_val),
            xycoords="data",
            xytext=(x_text, y_text),
            textcoords=ax.get_yaxis_transform(),
            ha=ha,
            va="center",
            fontsize=fontsize,
            clip_on=False,
            bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="#bdbdbd", alpha=0.94, linewidth=0.8),
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.8),
        )


def _neglog10(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    values = values.clip(lower=1e-300)
    return -np.log10(values)


def _extract_rsid_text(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    for token in re.split(r"[,&;\s]+", text):
        token = token.strip()
        if re.fullmatch(r"rs\d+", token, flags=re.IGNORECASE):
            return token
    return ""


def _extract_rsid(row: pd.Series) -> str:
    for key in ["rsID", "Existing_variation", "Variant_Label", "Variant_ID"]:
        rsid = _extract_rsid_text(row.get(key, ""))
        if rsid:
            return rsid
    return ""


def _preferred_variant_id(row: pd.Series) -> str:
    rsid = _extract_rsid(row)
    if rsid:
        return rsid
    for key in ["Variant_Label", "Variant_ID"]:
        text = str(row.get(key, "") or "").strip()
        if text and text.lower() not in {"nan", "none", "na", "<na>"}:
            return text
    return _extract_variant_id(row)


def _variant_label(row: pd.Series) -> str:
    return _variant_display_label(row, max_len=42, include_gene=True)


def _variant_display_label(row: pd.Series, max_len: int = 42, include_gene: bool = True) -> str:
    variant_id = _preferred_variant_id(row)
    gene_candidates = [
        str(row.get("SYMBOL", "") or "").strip(),
        str(row.get("Gene", "") or "").strip(),
    ]
    gene = ""
    for candidate in gene_candidates:
        if not candidate:
            continue
        if re.fullmatch(r"ENSG\d+(?:\.\d+)?", candidate, flags=re.IGNORECASE):
            continue
        gene = candidate
        break
    label = f"{gene} {variant_id}".strip() if include_gene and gene and gene != variant_id else variant_id
    return _short_label(label or variant_id or "Unlabelled variant", max_len=max_len)


def _format_ci95(lo: float, hi: float) -> object:
    if pd.isna(lo) or pd.isna(hi):
        return np.nan
    return f"[{float(lo):.3g}, {float(hi):.3g}]"


def _odds_ratio_ci95(a: int, b: int, c: int, d: int) -> Tuple[float, float, float, str]:
    a0, b0, c0, d0 = float(a), float(b), float(c), float(d)
    method = "Standard"
    if min(a0, b0, c0, d0) == 0:
        a0, b0, c0, d0 = [x + 0.5 for x in (a0, b0, c0, d0)]
        method = "Haldane-Anscombe"
    or_val = (a0 * d0) / (b0 * c0)
    se = np.sqrt((1.0 / a0) + (1.0 / b0) + (1.0 / c0) + (1.0 / d0))
    log_or = np.log(or_val)
    lo = float(np.exp(log_or - 1.96 * se))
    hi = float(np.exp(log_or + 1.96 * se))
    return float(or_val), lo, hi, method


def _plot_pivot_heatmap(
    df: pd.DataFrame,
    metric_col: str,
    out_path: Path,
    title: str,
    cbar_label: str,
    top_n: int = 12,
    significant_only: bool = False,
    annotate: bool = False,
    palette: str = "viridis",
    metric_name: Optional[str] = None,
):
    if df is None or df.empty or metric_col not in df.columns:
        _save_placeholder_plot(out_path, title, "No results available")
        return

    work = df.copy()
    work[metric_col] = pd.to_numeric(work[metric_col], errors="coerce")
    work = work[work[metric_col].notna()].copy()
    if work.empty:
        _save_placeholder_plot(out_path, title, "No results available")
        return

    if significant_only:
        sig_mask = pd.Series(False, index=work.index)
        if "Nominal_Sig_Unadj" in work.columns:
            sig_mask |= work["Nominal_Sig_Unadj"].fillna(False)
        if "FDR_Sig_Unadj" in work.columns:
            sig_mask |= work["FDR_Sig_Unadj"].fillna(False)
        work = work[sig_mask].copy()
        if work.empty:
            _save_placeholder_plot(out_path, title, "No significant rows")
            return

    work["Variant_Label"] = work.apply(_variant_label, axis=1)
    order = (
        work.groupby("Variant_Label")[metric_col]
        .min()
        .sort_values(ascending=True)
        .head(top_n)
        .index
        .tolist()
    )
    work = work[work["Variant_Label"].isin(order)].copy()
    work["Variant_Label"] = pd.Categorical(work["Variant_Label"], categories=order, ordered=True)
    work["_score"] = _neglog10(work[metric_col])
    column_key = "Clin_Label" if "Clin_Label" in work.columns else "Clinical_Var"

    label_order = list(dict.fromkeys(work.get("Clin_Label", pd.Series(dtype=object)).dropna().astype(str).tolist()))
    if not label_order and "Clinical_Var" in work.columns:
        label_order = list(dict.fromkeys(work["Clinical_Var"].dropna().astype(str).tolist()))

    pivot = work.pivot_table(
        index="Variant_Label",
        columns=column_key,
        values="_score",
        aggfunc="max",
    )
    metric_pivot = work.pivot_table(
        index="Variant_Label",
        columns=column_key,
        values=metric_col,
        aggfunc="min",
    )
    if label_order:
        pivot = pivot.reindex(columns=label_order)
        metric_pivot = metric_pivot.reindex(columns=label_order)
    pivot = pivot.sort_index(ascending=True)
    metric_pivot = metric_pivot.reindex(index=pivot.index, columns=pivot.columns)

    if pivot.empty:
        _save_placeholder_plot(out_path, title, "No plottable rows")
        return

    cohort_name = "Breast" if "Breast" in title else "Endometrial" if "Endometri" in title else "Global"
    cell_count = int(pivot.shape[0] * pivot.shape[1])
    show_annotations = bool(annotate) or (significant_only and cell_count <= 72)
    numeric_annotations = bool(annotate) and cell_count <= 56
    annotation_matrix = metric_pivot.copy().astype(object)
    for row_name in metric_pivot.index:
        for col_name in metric_pivot.columns:
            value = metric_pivot.loc[row_name, col_name]
            if pd.isna(value):
                annotation_matrix.loc[row_name, col_name] = ""
            elif numeric_annotations:
                annotation_matrix.loc[row_name, col_name] = _format_sig_annotation(value, metric_name=metric_name, show_ns=False)
            elif show_annotations:
                annotation_matrix.loc[row_name, col_name] = _sig_label(value) if float(value) < 0.05 else ""
            else:
                annotation_matrix.loc[row_name, col_name] = ""

    fig_w = max(9.6, 1.55 * max(3, len(pivot.columns)))
    fig_h = max(4.8, 0.52 * len(pivot.index) + 2.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = _cohort_heatmap_cmap(cohort_name, metric_name=metric_name) if palette == "viridis" else sns.color_palette(palette, as_cmap=True)
    hm = sns.heatmap(
        pivot,
        ax=ax,
        cmap=cmap,
        linewidths=0.8,
        linecolor="#F3F4F6",
        mask=pivot.isna(),
        annot=annotation_matrix if show_annotations else False,
        fmt="",
        cbar_kws={"label": cbar_label, "shrink": 0.88, "pad": 0.02, "aspect": 28},
        annot_kws={"size": 7.1, "fontweight": "bold", "ha": "center", "va": "center"},
    )
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(cbar_label, fontsize=9, labelpad=9)
    if metric_name:
        ax.set_title(_wrap_label(f"{title} ({metric_name})", width=64, max_lines=2), fontweight="bold", fontsize=11.5, loc="left", pad=12)
    else:
        ax.set_title(_wrap_label(title, width=64, max_lines=2), fontweight="bold", fontsize=11.5, loc="left", pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=0, labelsize=8, pad=4)
    ax.tick_params(axis="y", labelsize=8, pad=4)
    _set_wrapped_ticklabels(ax, axis="x", width=22, max_lines=3, rotation=0, ha="center", fontsize=8)
    _set_wrapped_ticklabels(ax, axis="y", width=28, max_lines=2, rotation=0, ha="right", fontsize=8)
    prefix = _metric_prefix(metric_name)
    if show_annotations:
        footnote = f"Cells show significance stars and {prefix}-values where space allows."
    else:
        footnote = f"Colour encodes -log10 {prefix}; annotations omitted to preserve readability."
    ax.text(0.0, -0.18, footnote, transform=ax.transAxes, ha="left", va="top", fontsize=8, color="#5B6472")
    _style_ax(ax, grid=False)
    _safe_tight_layout(fig, rect=[0.0, 0.03, 0.98, 0.95], pad=1.1)
    _save_figure(fig, out_path, pad_inches=0.20)


def _haldane_or(a: int, b: int, c: int, d: int) -> Tuple[float, str]:
    or_val, _, _, method = _odds_ratio_ci95(a, b, c, d)
    return or_val, method


def _apply_fdr(df: pd.DataFrame, p_col: str, out_col: str) -> pd.DataFrame:
    out = df.copy()
    pvals = pd.to_numeric(out.get(p_col), errors="coerce")
    mask = pvals.notna()
    out[out_col] = np.nan
    if mask.any():
        out.loc[mask, out_col] = false_discovery_control(pvals.loc[mask].values)
    return out


# ── DEFAULT PATHS ─────────────────────────────────────────────────────────────
DEFAULTS = script17_defaults()
DEFAULT_GSDMB  = DEFAULTS["gsdmb"]
DEFAULT_MASTER = DEFAULTS["master"]
DEFAULT_PHASED = DEFAULTS["phased"]
DEFAULT_OUT    = DEFAULTS["out_dir"]
DEFAULT_VARIANT_WHITELIST = DEFAULTS.get("variant_whitelist")

# Pass-BAM manifests - one file per cohort, listing QC-passed BAMs.
# These are the ground truth for which samples have actually been sequenced.
# Set to None to skip manifest filtering (falls back to extraction_flag logic).
DEFAULT_MANIFESTS = DEFAULTS.get("manifests", {
    "endometrium-tumour": Path(__file__).resolve().parent / "manifests/endometrium-tumour-pass_manifest.txt",
    "endometrium-normal": Path(__file__).resolve().parent / "manifests/endometrium-normal-pass_manifest.txt",
    "breast-tumour":      Path(__file__).resolve().parent / "manifests/breast-tumour-pass_manifest.txt",
    "breast-normal":      Path(__file__).resolve().parent / "manifests/breast-normal-pass_manifest.txt",
})

FDR_THRESHOLD  = float(DEFAULTS.get("fdr_threshold", 0.10))
MIN_CARRIERS   = int(DEFAULTS.get("min_carriers", 5))   # EPV-informed minimum; n=3 gives unreliable estimates
MIN_NFE_AF     = float(DEFAULTS.get("min_nfe_af", 0.01))
_VOLCANO_CLASS_ORDER = ["Not significant", "Nominal only", f"FDR < {FDR_THRESHOLD:g}"]
BMI_SOURCE_COLUMN = "canon__bmi"
BMI_CATEGORY_ORDER = ["Normal", "Overweight", "Obese"]
BMI_CATEGORY_COLORS = {
    "Normal": "#7FAE9B",
    "Overweight": "#C98A3D",
    "Obese": "#B56576",
}
BMI_COHORT_ORDER = ["Breast_Tumour", "Breast_Healthy", "Endometrial_Tumour", "Endometrial_Healthy"]

MIN_COMPARISON_CARRIERS = 3  # minimum cases/controls in each cell for stable comparison
MIN_EVENTS_FOR_LOGISTIC = 10  # EPV rule: ~10 events per predictor for stable logistic


# ── MANIFEST PARSING ──────────────────────────────────────────────────────────

def parse_manifests(manifest_paths: Dict[str, Path]) -> Tuple[set, set]:
    """
    Parse pass-BAM manifest files and return:
      sequenced_codes : set of primary snp_codes confirmed sequenced
      replicate_codes : set of formal _REP / Repeticion-derived codes observed in
                        the manifests. These are tracked for audit purposes only
                        and are not automatically excluded from analysis, because
                        the source workbooks indicate they can represent valid
                        resequenced replacement samples.

    Manifest format: one BAM path per line, e.g.:
      /path/to/pass_bams/DNA_SNP_MT-T_10_IonCode_0108_MamaHER2_LIB1y2_874.bam
      /path/to/pass_bams/DNA_MT_T_17_Repeticion_IonCode_0119.bam

    SNP code extraction mirrors the logic in _extract_snp_code().
    All four manifests are parsed and codes pooled across cohorts.
    """
    sequenced: set = set()
    replicates: set = set()

    for cohort_label, manifest_path in manifest_paths.items():
        if manifest_path is None or not Path(manifest_path).exists():
            print(f"      WARNING: manifest not found for {cohort_label}: {manifest_path}")
            continue

        with open(manifest_path) as fh:
            lines = [l.strip() for l in fh if l.strip()]

        for bam_path in lines:
            basename = os.path.basename(bam_path).replace(".bam", "")
            code = _extract_snp_code(basename)
            if code is None:
                print(f"      WARNING: could not parse snp_code from: {bam_path}")
                continue
            is_rep = bool(re.search(r"repeticion", basename, re.IGNORECASE))
            if is_rep:
                replicates.add(code)
            else:
                sequenced.add(code)

    print(f"      Manifest: {len(sequenced)} unique primary codes, "
          f"{len(replicates)} _REP/Repeticion resequenced codes")
    return sequenced, replicates


# ── SNP CODE EXTRACTION ───────────────────────────────────────────────────────

def _extract_snp_code(sample_name: str) -> Optional[str]:
    """
    Derive the SNP_XX_NN lab code from the sequencer sample name.
    Examples:
      SNP_AT_2__BL_039__IonCode_0118     ? SNP_AT_2
      DNA_SNP_EN_10_IonCode_0101_...     ? SNP_EN_10
      DNA_MT-T_37_IonCode_0121           ? SNP_MT-T_37
      DNA_MT_T_17_Repeticion_IonCode_0119 ? SNP_MT-T_17_REP
    """
    s = str(sample_name)
    m = re.match(r"^(SNP_MT-T_\d+_REP)", s, re.IGNORECASE)
    if m: return m.group(1).upper()
    m = re.match(r"^(SNP_(?:AT|EN|MN|MT-T)_\d+)", s, re.IGNORECASE)
    if m: return m.group(1).upper()
    m = re.match(r"^DNA_MT[-_]T_(\d+)_Repeticion_", s, re.IGNORECASE)
    if m: return f"SNP_MT-T_{m.group(1)}_REP"
    m = re.match(r"^DNA_SNP_(EN|MN)_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_{m.group(1).upper()}_{m.group(2)}"
    m = re.match(r"^DNA_SNP_AT_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_AT_{int(m.group(1))}"
    m = re.match(r"^DNA_AT_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_AT_{int(m.group(1))}"
    m = re.match(r"^DNA_SNP_(?:CK|RSB)_ECLAI_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_AT_{int(m.group(1))}"
    m = re.match(r"^SNP_DNA_AT_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_AT_{m.group(1)}"
    m = re.match(r"^SNP_DNA_(MN|EN)_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_{m.group(1).upper()}_{m.group(2)}"
    m = re.match(r"^DNA_SNP_MT-T_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_MT-T_{m.group(1)}"
    m = re.match(r"^DNA_MT[-_]T_(\d+)_", s, re.IGNORECASE)
    if m: return f"SNP_MT-T_{m.group(1)}"
    return None


# ── CLINICAL VARIABLE DERIVATIONS ─────────────────────────────────────────────

def _derive_breast_recurrence(val) -> Optional[int]:
    """1 if any recurrence/progression text, 0 if 'NO', else NaN."""
    if pd.isna(val): return None
    return 0 if str(val).strip().upper() == "NO" else 1


def _derive_breast_exitus(val) -> Optional[int]:
    """1 if SI/YES, 0 if NO."""
    if pd.isna(val): return None
    v = str(val).strip().upper()
    if v.startswith("SI") or v.startswith("YES"): return 1
    if v.startswith("NO"): return 0
    return None


def _derive_breast_metastasis(val) -> Optional[int]:
    """1 = distant metastasis (SI), 0 = NO or local only."""
    if pd.isna(val): return None
    v = str(val).strip().upper()
    if v == "SI": return 1
    if v in ("NO", "NO-LOCAL"): return 0
    return None


def _derive_her2_subtype(val) -> Optional[str]:
    """Collapse DX field to HER2+, TN, or Other."""
    if pd.isna(val): return None
    v = str(val).strip().upper()
    if "HER2+" in v or "HER2 +" in v: return "HER2+"
    if "TN" in v or "TRIPLE" in v:     return "TN"
    return "Other"


def _derive_figo_numeric(val) -> Optional[float]:
    """Convert FIGO stage string to numeric for ordinal trend tests."""
    if pd.isna(val): return None
    mapping = {
        "IA": 1.0, "IB": 1.5, "I": 1.5,
        "II": 2.0,
        "IIIA": 3.1, "IIIB": 3.2, "IIIC1": 3.3, "IIIC2": 3.4, "III": 3.0,
        "IVA": 4.1, "IVB": 4.2, "IV": 4.0,
    }
    return mapping.get(str(val).strip().upper(), None)


def _derive_endo_grade(val) -> Optional[float]:
    """G1→1, G2→2, G3→3."""
    if pd.isna(val): return None
    m = re.search(r"(\d)", str(val))
    return float(m.group(1)) if m else None


def _derive_risk_ordinal(val) -> Optional[float]:
    """LOW=1, INTERMEDIATE=2, INTERMEDIATE-HIGH=3, HIGH=4."""
    if pd.isna(val): return None
    v = str(val).strip().upper()
    mapping = {
        "LOW": 1.0, "INTERMEDIATE": 2.0,
        "INTERMEDIATE-HIGH": 3.0, "HIGH": 4.0,
    }
    return mapping.get(v, None)


def _derive_ki67_numeric(val) -> Optional[float]:
    """
    Coerce messy KI67 to a fraction in [0, 1].

    Handles inputs such as:
      ? Integer percentage:  40      ? 0.40
      ? String percentage:   '40%'   ? 0.40
      ? Decimal fraction:    0.40    ? 0.40  (already a fraction, left as-is)
      ? Range (percent):     '20-25%'? midpoint 22.5 ? 0.225
      ? Range (fraction):    '0.2-0.25' ? midpoint 0.225 (already a fraction)

    The heuristic: if the raw numeric value (or midpoint) is > 1 it is assumed
    to be a percentage and is divided by 100.  Values ? 1 are treated as
    already being a fraction and are returned unchanged.
    """
    if pd.isna(val):
        return None
    s = str(val).replace("%", "").strip()
    m = re.match(r"([\d.]+)\s*[-?]\s*([\d.]+)", s)
    if m:
        mid = (float(m.group(1)) + float(m.group(2))) / 2
        return mid / 100 if mid > 1 else mid
    try:
        v = float(s)
        return v / 100 if v > 1 else v
    except (ValueError, TypeError):
        return None


def _derive_breast_os(row) -> Optional[float]:
    """Fallback OS in months from diagnosis to last available date."""
    try:
        dx = pd.to_datetime(
            row.get("canon__date_diagnosis", row.get("clin_dcs__Fecha dx")),
            errors="coerce",
        )
        ult = pd.to_datetime(
            row.get("canon__date_last_fu", row.get("clin_dcs__Última fecha disponible")),
            errors="coerce",
        )
        if pd.isna(dx) or pd.isna(ult):
            return None
        return max((ult - dx).days / 30.44, 0)
    except Exception:
        return None



def _ki67_fraction_from_pct(val) -> Optional[float]:
    """Convert a 0-100 KI67 percentage into a 0-1 fraction."""
    if pd.isna(val):
        return None
    try:
        v = float(val)
    except (TypeError, ValueError):
        return None
    if 0 <= v <= 100:
        return round(v / 100.0, 4)
    return None


def _normalise_treatment_text(val) -> str:
    if pd.isna(val):
        return ""
    text = str(val).strip().lower()
    text = re.sub(r"\s+", " ", text)
    return text


def _binary_yes_no(series: pd.Series) -> pd.Series:
    values = pd.Series(index=series.index, dtype=float)
    numeric = pd.to_numeric(series, errors="coerce")
    values.loc[numeric.notna() & (numeric == 1)] = 1
    values.loc[numeric.notna() & (numeric == 0)] = 0

    text = series.astype(str).str.strip().str.upper()
    mapping = {
        "YES": 1, "Y": 1, "SI": 1, "SÍ": 1, "TRUE": 1, "POSITIVE": 1, "POS": 1, "1": 1, "1.0": 1,
        "NO": 0, "N": 0, "FALSE": 0, "NEGATIVE": 0, "NEG": 0, "0": 0, "0.0": 0,
    }
    text_mapped = text.map(mapping)
    values = values.combine_first(text_mapped.astype(float))
    return values


def _categorize_bmi_value(value) -> Optional[str]:
    bmi = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(bmi) or bmi < 20:
        return None
    if bmi < 25:
        return "Normal"
    if bmi < 30:
        return "Overweight"
    return "Obese"


def _derive_breast_treatment_exposure_bins(treatment_series: pd.Series) -> pd.DataFrame:
    cleaned = treatment_series.fillna("").astype(str).map(_normalise_treatment_text)
    meaningful = ~cleaned.isin({"", "x", "na", "n/a", "none", "no consta"})
    flags = {}
    patterns = {
        "BREAST_TX_CHEMOTHERAPY_BIN": r"\bac\b|\bfec\b|\bcmf\b|taxol|paclitaxel|docetaxel|cbdca|carbo|quimio|chemo|\bqt\b",
        "BREAST_TX_ANTI_HER2_BIN": r"hercept|trastu|lapat|pertu|anti[\s-]?her[\s-]?2",
        "BREAST_TX_ENDOCRINE_BIN": r"tamox|letroz|exemest|anastro|arimid|fulves|hormon|terapia hormonal|\bht\b",
        "BREAST_TX_RADIOTHERAPY_BIN": r"\brt\b|radiot|radioter|rte|rdt",
    }
    for col, pattern in patterns.items():
        out = pd.Series(np.nan, index=treatment_series.index, dtype=float)
        out.loc[meaningful] = cleaned.loc[meaningful].str.contains(pattern, regex=True).astype(int)
        flags[col] = out
    return pd.DataFrame(flags)


def _ensure_breast_treatment_exposure_columns(master: pd.DataFrame) -> pd.DataFrame:
    master = master.copy()
    canon_to_bin = {
        "canon__treatment_chemotherapy": "BREAST_TX_CHEMOTHERAPY_BIN",
        "canon__treatment_anti_her2": "BREAST_TX_ANTI_HER2_BIN",
        "canon__treatment_endocrine": "BREAST_TX_ENDOCRINE_BIN",
        "canon__treatment_radiotherapy": "BREAST_TX_RADIOTHERAPY_BIN",
    }
    for canon_col, bin_col in canon_to_bin.items():
        if canon_col in master.columns:
            master[bin_col] = _binary_yes_no(master[canon_col])
    missing = [bin_col for bin_col in canon_to_bin.values() if bin_col not in master.columns]
    if missing and "canon__treatment" in master.columns:
        derived = _derive_breast_treatment_exposure_bins(master["canon__treatment"])
        for bin_col in missing:
            master[bin_col] = derived[bin_col]
    return master


def _extract_variant_id(row: pd.Series) -> str:
    """
    Prefer rsIDs where available; otherwise fall back to a genomic label.
    This keeps downstream plots readable while remaining deterministic for
    novel or non-rs variants.
    """
    rsid = _extract_rsid(row)
    if rsid:
        return rsid

    existing = str(row.get("Existing_variation", "") or "")

    chrom = str(row.get("CHROM", "") or "").strip()
    pos = pd.to_numeric(row.get("POS"), errors="coerce")
    ref = str(row.get("REF", "") or "").strip()
    alt = str(row.get("ALT", "") or "").strip()
    if chrom and pd.notna(pos) and ref and alt:
        return f"{chrom}:{int(pos)}_{ref}>{alt}"
    if chrom and pd.notna(pos):
        return f"{chrom}:{int(pos)}"
    return existing.split("&")[0].strip() if existing else "UNKNOWN_VARIANT"


def _load_variant_whitelist_filters(whitelist_path: Path) -> Dict[str, set]:
    """
    Read a flexible Excel whitelist and return identifiers that can be used to
    filter the annotated variant table. The whitelist can contain rsIDs, a
    Variant_ID column, or genomic coordinates (CHROM/POS/REF/ALT).
    """
    wl = pd.read_excel(whitelist_path)
    wl.columns = [str(c).strip() for c in wl.columns]

    ids: set = set()
    rsids: set = set()
    coord_keys: set = set()
    positions: set = set()

    for col in wl.columns:
        series = wl[col].dropna().astype(str).str.strip()
        ids.update(v for v in series if v)
        rsids.update(v for v in series if re.fullmatch(r"rs\d+", v, flags=re.IGNORECASE))

    pos_col = next((c for c in wl.columns if re.fullmatch(r"POS|Position", c, flags=re.IGNORECASE)), None)
    ref_col = next((c for c in wl.columns if c.upper() == "REF"), None)
    alt_col = next((c for c in wl.columns if c.upper() == "ALT"), None)
    chrom_col = next((c for c in wl.columns if c.upper() == "CHROM"), None)

    if pos_col is not None:
        pos_vals = pd.to_numeric(wl[pos_col], errors="coerce").dropna().astype(int)
        positions.update(pos_vals.tolist())

    if all(c is not None for c in [chrom_col, pos_col, ref_col, alt_col]):
        tmp = wl[[chrom_col, pos_col, ref_col, alt_col]].copy()
        tmp[pos_col] = pd.to_numeric(tmp[pos_col], errors="coerce")
        tmp = tmp.dropna(subset=[pos_col, ref_col, alt_col])
        tmp[chrom_col] = tmp[chrom_col].astype(str).str.strip()
        tmp[ref_col] = tmp[ref_col].astype(str).str.strip().str.upper()
        tmp[alt_col] = tmp[alt_col].astype(str).str.strip().str.upper()
        coord_keys.update(
            f"{r[chrom_col]}|{int(r[pos_col])}|{r[ref_col]}|{r[alt_col]}"
            for _, r in tmp.iterrows()
        )

    return {
        "ids": {v for v in ids if v},
        "rsids": {v.lower() for v in rsids if v},
        "coord_keys": coord_keys,
        "positions": positions,
    }


CLINICAL_VARS_BREAST: Dict[str, Dict] = {
    # Lean default set: prognosis-first breast variables only.
    "canon__recurrence_flag":   {"type": "binary", "label": "Recurrence / progression", "note": "Harmonised recurrence/progression flag"},
    "canon__distant_mets_flag": {"type": "binary", "label": "Distant metastasis", "note": "Harmonised distant-metastasis flag"},
    "canon__exitus_flag":       {"type": "binary", "label": "Death", "note": "Harmonised overall death flag"},
    "canon__os_months":         {"type": "continuous", "label": "Overall survival (months)", "note": "Harmonised OS months"},
}

CLINICAL_VARS_ENDO: Dict[str, Dict] = {
    # Lean default set: prognosis-first endometrial variables with clear labels.
    "canon__pd_flag":     {"type": "binary",     "label": "Progression", "note": "Harmonised progressive-disease flag"},
    "canon__exitus_flag": {"type": "binary",     "label": "Death", "note": "Harmonised overall death flag"},
    "canon__os_months":   {"type": "continuous", "label": "Overall survival (months)", "note": "Harmonised OS months"},
    "canon__pfs_months":  {"type": "continuous", "label": "Progression-free survival (months)", "note": "Harmonised PFS months"},
    "canon__risk_group":  {"type": "nominal",    "label": "Risk of recurrence", "note": "Harmonised labelled risk-of-recurrence groups"},
}


BREAST_TREATMENT_EXPOSURES: Dict[str, Dict[str, str]] = {
    "BREAST_TX_CHEMOTHERAPY_BIN": {"label": "Chemotherapy exposure", "canon": "canon__treatment_chemotherapy"},
    "BREAST_TX_ANTI_HER2_BIN":    {"label": "Anti-HER2 exposure", "canon": "canon__treatment_anti_her2"},
    "BREAST_TX_ENDOCRINE_BIN":    {"label": "Endocrine therapy exposure", "canon": "canon__treatment_endocrine"},
    "BREAST_TX_RADIOTHERAPY_BIN": {"label": "Radiotherapy exposure", "canon": "canon__treatment_radiotherapy"},
}
MIN_TREATMENT_STRATUM_SAMPLES = 15


def load_clinical_master(master_path: Path,
                         manifest_paths: Optional[Dict[str, Path]] = None) -> pd.DataFrame:
    """
    Load the harmonised master, restrict to DNA samples, optionally apply the
    pass-BAM manifest filter, and derive the clinical columns required by the
    downstream association analyses.
    """
    print("  Loading harmonised clinical master ?")
    try:
        master = pd.read_excel(master_path, sheet_name="harmonised_plus_canon", engine="openpyxl")
    except Exception:
        print("  WARNING: 'harmonised_plus_canon' not found ? trying sheet 0")
        master = pd.read_excel(master_path, sheet_name=0, engine="openpyxl")

    master.columns = [str(c).strip() for c in master.columns]

    if "nucleic_acid" in master.columns:
        n_before = len(master)
        master = master[master["nucleic_acid"].astype(str).str.upper() == "DNA"].copy()
        print(f"  DNA rows retained: {len(master)} / {n_before}")

    if "snp_code" not in master.columns:
        for alt in ["SNP_code", "Sample", "SAMPLE", "sample_id"]:
            if alt in master.columns:
                master = master.rename(columns={alt: "snp_code"})
                break
    if "snp_code" not in master.columns:
        raise ValueError("Cannot find snp_code column in harmonised master workbook")

    master["snp_code"] = (
        master["snp_code"]
        .astype(str)
        .str.strip()
        .str.replace(r"\s+", "", regex=True)
    )

    replicate_codes: set = set()
    if manifest_paths:
        sequenced_codes, replicate_codes = parse_manifests(manifest_paths)
        allowed_codes = sequenced_codes | replicate_codes
        before = master["snp_code"].nunique()
        master = master[master["snp_code"].isin(allowed_codes)].copy()
        after = master["snp_code"].nunique()
        print(
            "  Manifest filter retained "
            f"{after} / {before} DNA samples "
            f"({len(sequenced_codes)} primary + {len(replicate_codes)} replicate codes)"
        )
    else:
        if "extraction_flag" in master.columns and "sheet" in master.columns:
            excl = (
                (master["sheet"] == "AT=AUs")
                & (
                    master["extraction_flag"].astype(str).str.upper().eq("NO HAY")
                    | master.get("pd_status", pd.Series("", index=master.index)).astype(str).str.upper().eq("NO HACER")
                )
            )
            dropped = int(excl.sum())
            if dropped:
                print(f"  Fallback filter removed {dropped} AU rows flagged NO HAY / NO HACER")
            master = master[~excl].copy()

    if "is_replicate" in master.columns:
        master["is_replicate"] = master["is_replicate"].fillna(False).astype(bool)
    else:
        master["is_replicate"] = False

    sheet_to_tissue = {
        "MT-T_N": "Tumour",
        "AT=AUs": "Tumour",
        "MN": "Healthy",
        "EN": "Healthy",
    }
    sheet_to_cohort = {
        "MT-T_N": "Breast_Tumour",
        "MN": "Breast_Healthy",
        "AT=AUs": "Endometrial_Tumour",
        "EN": "Endometrial_Healthy",
    }
    if "sheet" in master.columns:
        master["Tissue"] = master["sheet"].map(sheet_to_tissue).fillna(master.get("tissue"))
        master["Cohort"] = master["sheet"].map(sheet_to_cohort).fillna(master["sheet"])
    else:
        tissue_raw = master.get("tissue", pd.Series("", index=master.index)).astype(str).str.strip()
        master["Tissue"] = tissue_raw.replace({"Normal": "Healthy", "Tumour": "Tumour", "Tumor": "Tumour"})
        master["Cohort"] = master.get("cohort", "Unknown")

    if BMI_SOURCE_COLUMN not in master.columns:
        master[BMI_SOURCE_COLUMN] = np.nan
    master[BMI_SOURCE_COLUMN] = pd.to_numeric(master[BMI_SOURCE_COLUMN], errors="coerce")
    master["BMI_Category"] = master[BMI_SOURCE_COLUMN].apply(_categorize_bmi_value)
    master["BMI_Usable"] = master["BMI_Category"].notna()

    master = _ensure_breast_treatment_exposure_columns(master)

    if "canon__her2_copies" not in master.columns:
        master["canon__her2_copies"] = pd.to_numeric(
            master.get("clin_dcs__COPIAS HER2", master.get("clin_her2__COPIAS HER2", pd.Series(np.nan, index=master.index))),
            errors="coerce",
        )

    for src, dst, fn in [
        ("clin_dcs__Recaida/Progresión", "BREAST_RECURRENCE_DERIVED", _derive_breast_recurrence),
        ("clin_dcs__Exitus",             "BREAST_EXITUS_DERIVED",     _derive_breast_exitus),
        ("clin_dcs__MTxDISTANCIA",       "BREAST_METASTASIS_DERIVED", _derive_breast_metastasis),
        ("clin_her2__DX",                "BREAST_HER2_SUBTYPE",       _derive_her2_subtype),
    ]:
        if src in master.columns:
            master[dst] = master[src].apply(fn)

    recurrence_from_canon = pd.to_numeric(
        master.get("canon__recurrence_flag", pd.Series(pd.NA, index=master.index)),
        errors="coerce",
    )
    exitus_from_canon = pd.to_numeric(
        master.get("canon__exitus_flag", pd.Series(pd.NA, index=master.index)),
        errors="coerce",
    )
    mets_from_canon = master.get("canon__distant_mets", pd.Series(pd.NA, index=master.index)).map(
        lambda x: pd.NA if pd.isna(x) else (1 if str(x).strip().upper() == "YES" else (0 if str(x).strip().upper() == "NO" else pd.NA))
    )
    master["BREAST_RECURRENCE_DERIVED"] = pd.to_numeric(
        master.get("BREAST_RECURRENCE_DERIVED", pd.Series(np.nan, index=master.index)),
        errors="coerce",
    ).combine_first(recurrence_from_canon)
    master["BREAST_EXITUS_DERIVED"] = pd.to_numeric(
        master.get("BREAST_EXITUS_DERIVED", pd.Series(np.nan, index=master.index)),
        errors="coerce",
    ).combine_first(exitus_from_canon)
    master["BREAST_METASTASIS_DERIVED"] = pd.to_numeric(
        master.get("BREAST_METASTASIS_DERIVED", pd.Series(np.nan, index=master.index)),
        errors="coerce",
    ).combine_first(pd.to_numeric(mets_from_canon, errors="coerce"))

    if "canon__er_status" in master.columns:
        master["BREAST_ER_BIN"] = master["canon__er_status"].map({"Positive": 1, "Negative": 0})
        master["ENDO_ER_BIN"] = master["canon__er_status"].map({"Positive": 1, "Negative": 0})
    elif "clin_dcs__RE" in master.columns:
        master["BREAST_ER_BIN"] = (
            master["clin_dcs__RE"].astype(str).str.strip().str.upper().map({"POSITIVO": 1, "NEGATIVO": 0})
        )

    if "canon__pr_status" in master.columns:
        master["BREAST_PR_BIN"] = master["canon__pr_status"].map({"Positive": 1, "Negative": 0})
        master["ENDO_PR_BIN"] = master["canon__pr_status"].map({"Positive": 1, "Negative": 0})
    elif "clin_dcs__RP" in master.columns:
        master["BREAST_PR_BIN"] = (
            master["clin_dcs__RP"].astype(str).str.strip().str.upper().map({"POSITIVO": 1, "NEGATIVO": 0})
        )

    ki67_pct = master.get("canon__ki67_pct", pd.Series(np.nan, index=master.index))
    grade = master.get("canon__grade", pd.Series(np.nan, index=master.index))
    figo = master.get("canon__figo_stage", pd.Series(np.nan, index=master.index))
    risk = master.get("canon__risk_group", pd.Series(np.nan, index=master.index))
    lvsi = master.get("canon__lvsi", pd.Series(pd.NA, index=master.index))
    myoinv = master.get("canon__myometrial_invasion", pd.Series(pd.NA, index=master.index))
    pd_flag = master.get("canon__pd_flag", pd.Series(pd.NA, index=master.index))
    exitus_flag = master.get("canon__exitus_flag", pd.Series(pd.NA, index=master.index))
    breast_os_raw = master.get("canon__os_months", pd.Series(np.nan, index=master.index))

    master["BREAST_KI67_NUMERIC"] = ki67_pct.apply(_ki67_fraction_from_pct)
    master["BREAST_GRADE_NUMERIC"] = pd.to_numeric(grade, errors="coerce")
    master["BREAST_OS_MONTHS_DERIVED"] = pd.to_numeric(breast_os_raw, errors="coerce").combine_first(
        master.apply(_derive_breast_os, axis=1)
    )

    master["ENDO_FIGO_NUMERIC"] = figo.apply(_derive_figo_numeric)
    master["ENDO_GRADE_NUMERIC"] = pd.to_numeric(grade, errors="coerce")
    master["ENDO_RISK_ORDINAL"] = risk.apply(_derive_risk_ordinal)
    master["ENDO_LVSI_BIN"] = lvsi.map({"Yes": 1, "No": 0, "YES": 1, "NO": 0})
    master["ENDO_MYOINV_BIN"] = myoinv.map({"<50%": 0, ">50%": 1})
    master["ENDO_MSI_BIN"] = master.get("canon__msi_status", pd.Series(pd.NA, index=master.index)).map({"Unstable": 1, "Stable": 0})
    master["ENDO_NEEC_BIN"] = master.get("clin_au_endo__HISTOLOGY_GROUP", pd.Series(pd.NA, index=master.index)).map({"NEEC": 1, "EEC": 0})
    master["ENDO_PD_BIN"] = pd.to_numeric(pd_flag, errors="coerce")
    master["ENDO_EXITUS_BIN"] = pd.to_numeric(exitus_flag, errors="coerce")
    master["ENDO_EXITUS_DISEASE_BIN"] = (
        master.get("clin_au_endo__EXITUS_DISEASE", pd.Series(pd.NA, index=master.index))
        .astype(str).str.strip().str.upper().map({"YES": 1, "NO": 0})
    )
    master["ENDO_TP53_ABN_BIN"] = master.get("canon__p53_status", pd.Series(pd.NA, index=master.index)).map({"Aberrant": 1, "Normal": 0})
    if "clin_au_endo__M" in master.columns:
        master["ENDO_M_STAGE_BIN"] = master["clin_au_endo__M"].apply(
            lambda x: 0 if pd.notna(x) and str(x).strip().upper() == "M0" else (1 if pd.notna(x) else None)
        )

    master = master.drop_duplicates("snp_code").copy()
    print(f"  Clinical master ready: {len(master)} unique DNA samples")
    if "sheet" in master.columns:
        print(f"  Sheet breakdown: {master['sheet'].value_counts().to_dict()}")
    return master


def load_and_merge(gsdmb_path: Path,
                   master_path: Path,
                   manifest_paths: Optional[Dict[str, Path]] = None,
                   whitelist_path: Optional[Path] = None) -> pd.DataFrame:
    """
    Join the annotated GSDMB callset to the harmonised clinical master.

    The merge is intentionally a sample-level inner join on `snp_code`: only
    samples present in the clinical master and represented in the annotated
    callset are retained, matching the documented script-17 behaviour.
    """
    print("=== Loading GSDMB SNP calls and harmonised clinical data ===")
    master = load_clinical_master(master_path, manifest_paths=manifest_paths)

    print("  Loading annotated GSDMB variants ?")
    variants = pd.read_excel(gsdmb_path, sheet_name="Biological_Annotations")
    variants.columns = [str(c).strip() for c in variants.columns]

    if "Sample" not in variants.columns:
        raise ValueError("Annotated GSDMB workbook is missing the Sample column")

    variants["Sample"] = variants["Sample"].astype(str).str.strip()
    variants = attach_analysis_sample_ids(
        variants,
        raw_col="Sample",
        phased_path=DEFAULT_PHASED,
        extract_snp_code=_extract_snp_code,
    )
    failed_rows = int(variants["snp_code"].isna().sum())
    if failed_rows:
        examples = variants.loc[variants["snp_code"].isna(), "Sample"].drop_duplicates().head(8).tolist()
        print(f"  WARNING: could not derive snp_code for {failed_rows} rows; examples: {examples}")
    variants = variants.dropna(subset=["snp_code", "analysis_sample_id"]).copy()

    if "GT" in variants.columns:
        variants["GT"] = (
            variants["GT"].astype(str).str.strip()
            .replace({"0|0": "0/0", "0|1": "0/1", "1|0": "1/0", "1|1": "1/1"})
        )

    if "Cohort" in variants.columns:
        variants = variants.rename(columns={"Cohort": "Cohort_callset"})
    if "Tissue" in variants.columns:
        variants = variants.rename(columns={"Tissue": "Tissue_callset"})

    variants["rsID"] = variants.apply(_extract_rsid, axis=1).replace("", np.nan)
    variants["Variant_ID"] = variants.apply(_extract_variant_id, axis=1)
    variants["Variant_Label"] = variants["rsID"].fillna(variants["Variant_ID"])
    if all(c in variants.columns for c in ["CHROM", "POS", "REF", "ALT"]):
        variants["_coord_key"] = (
            variants["CHROM"].astype(str).str.strip() + "|"
            + pd.to_numeric(variants["POS"], errors="coerce").fillna(-1).astype(int).astype(str) + "|"
            + variants["REF"].astype(str).str.strip().str.upper() + "|"
            + variants["ALT"].astype(str).str.strip().str.upper()
        )
    else:
        variants["_coord_key"] = ""

    if whitelist_path is not None:
        filters = _load_variant_whitelist_filters(whitelist_path)
        rsid_mask = variants["Variant_ID"].astype(str).str.lower().isin(filters["rsids"])
        id_mask = variants["Variant_ID"].astype(str).isin(filters["ids"])
        coord_mask = variants["_coord_key"].isin(filters["coord_keys"])
        pos_mask = pd.to_numeric(variants.get("POS"), errors="coerce").fillna(-1).astype(int).isin(filters["positions"])
        keep_mask = id_mask | rsid_mask | coord_mask | pos_mask
        before = variants["Variant_ID"].nunique()
        variants = variants[keep_mask].copy()
        after = variants["Variant_ID"].nunique()
        print(f"  Variant whitelist retained {after} / {before} unique variants")

    # Sample identity is now tracked at analysis_sample_id level: raw runs are
    # collapsed only when their phased genotype columns are exactly identical.
    dedup_cols = [c for c in ["analysis_sample_id", "Variant_ID", "CHROM", "POS", "REF", "ALT", "GT"] if c in variants.columns]
    if dedup_cols:
        variants = variants.drop_duplicates(subset=dedup_cols).copy()

    merged = variants.merge(master, on="snp_code", how="inner")
    if "Sample" in merged.columns:
        merged["Sample_callset_raw"] = merged["Sample"]
    merged["Sample"] = merged["analysis_sample_id"]

    print(
        f"  Annotated analysis samples: {variants['analysis_sample_id'].nunique()} | "
        f"clinical samples: {master['snp_code'].nunique()} | "
        f"matched analysis samples: {merged['analysis_sample_id'].nunique()}"
    )
    print(f"  Matched rows: {len(merged)} | unique variants: {merged['Variant_ID'].nunique()}")
    if not merged.empty:
        print(f"  Tissue breakdown: {merged.drop_duplicates('Sample')['Tissue'].value_counts().to_dict()}")
        print(f"  Cohort breakdown: {merged.drop_duplicates('Sample')['Cohort'].value_counts().to_dict()}")

    return merged.drop(columns=["_coord_key"], errors="ignore")


def tumour_vs_control(merged: pd.DataFrame) -> pd.DataFrame:
    """
    Analysis 1: compare tumour vs pooled healthy controls for each variant.

    Cohort-specific comparisons use the relevant tumour arm (Breast or
    Endometrium) against the pooled healthy control arm; a Global comparison
    uses all tumours vs all healthy controls.
    """
    print("=== Analysis 1: Tumour vs Control ===")
    df = merged[~merged["is_replicate"]].copy()
    if df.empty:
        print("  No merged rows - skipping.\n")
        return pd.DataFrame()

    sample_manifest = df.drop_duplicates("Sample").copy()
    pooled_controls = sample_manifest[sample_manifest["Tissue"] == "Healthy"].copy()

    comparison_defs = {
        "Breast": {
            "mask": (sample_manifest["Tissue"] == "Tumour")
                    & sample_manifest["Cohort"].astype(str).str.contains("Breast", case=False, na=False),
        },
        "Endometrium": {
            "mask": (sample_manifest["Tissue"] == "Tumour")
                    & sample_manifest["Cohort"].astype(str).str.contains("Endometri", case=False, na=False),
        },
        "Global": {
            "mask": sample_manifest["Tissue"] == "Tumour",
        },
    }

    rows = []
    for analysis_group, cfg in comparison_defs.items():
        tumour_samples = sample_manifest.loc[cfg["mask"]].copy()
        if tumour_samples.empty or pooled_controls.empty:
            continue

        n_tumour = len(tumour_samples)
        n_control = len(pooled_controls)
        print(f"  {analysis_group}: tumour n={n_tumour}, pooled control n={n_control}")
        if n_tumour < MIN_CARRIERS or n_control < MIN_CARRIERS:
            continue

        tumour_df = df[df["Sample"].isin(tumour_samples["Sample"])]
        control_df = df[df["Sample"].isin(pooled_controls["Sample"])]
        all_variants = pd.concat(
            [tumour_df[["Variant_ID"]], control_df[["Variant_ID"]]],
            ignore_index=True,
        )["Variant_ID"].dropna().unique()

        for var_id in all_variants:
            tum_carriers = set(tumour_df.loc[tumour_df["Variant_ID"] == var_id, "Sample"].unique())
            ctl_carriers = set(control_df.loc[control_df["Variant_ID"] == var_id, "Sample"].unique())

            a = len(tum_carriers)
            b = n_tumour - a
            c = len(ctl_carriers)
            d = n_control - c

            if a + c < MIN_COMPARISON_CARRIERS:
                continue

            var_rows = df[df["Variant_ID"] == var_id]
            gene = var_rows["SYMBOL"].dropna().iloc[0] if "SYMBOL" in var_rows.columns and var_rows["SYMBOL"].notna().any() else ""
            consequence = var_rows["Consequence"].dropna().iloc[0] if "Consequence" in var_rows.columns and var_rows["Consequence"].notna().any() else np.nan
            impact = var_rows["IMPACT"].dropna().iloc[0] if "IMPACT" in var_rows.columns and var_rows["IMPACT"].notna().any() else np.nan

            _, p_val = fisher_exact([[a, b], [c, d]])
            or_val, or_lo, or_hi, or_method = _odds_ratio_ci95(a, b, c, d)

            rows.append({
                "Analysis_Group": analysis_group,
                "Variant_ID": var_id,
                "Gene": gene,
                "Consequence": consequence,
                "IMPACT": impact,
                "N_Tumour": n_tumour,
                "N_Control": n_control,
                "Carriers_Tumour": a,
                "Carriers_Control": c,
                "Freq_Tumour_%": round(a / n_tumour * 100, 2),
                "Freq_Control_%": round(c / n_control * 100, 2),
                "Odds_Ratio": round(or_val, 4),
                "Odds_Ratio_CI95": _format_ci95(or_lo, or_hi),
                "Odds_Ratio_CI95_Lower": or_lo,
                "Odds_Ratio_CI95_Upper": or_hi,
                "OR_Method": or_method,
                "P_Value": p_val,
            })

    res = pd.DataFrame(rows)
    if res.empty:
        print("  No results.\n")

        return res

    # FDR correction within each group
    parts = []
    for grp in res["Analysis_Group"].unique():
        sub = res[res["Analysis_Group"] == grp].copy()
        sub = _apply_fdr(sub, p_col="P_Value", out_col="FDR_P_Value")
        parts.append(sub)

    res = pd.concat(parts, ignore_index=True).sort_values(["Analysis_Group", "P_Value"])
    res["Nominal_Sig"] = res["P_Value"] < 0.05
    res["FDR_Sig"] = res["FDR_P_Value"] < FDR_THRESHOLD

    print(f"  {len(res)} tests across {res['Analysis_Group'].nunique()} groups")
    for grp in res["Analysis_Group"].unique():
        sub = res[res["Analysis_Group"] == grp]
        print(f"    {grp}: {sub['Nominal_Sig'].sum()} nominal | {sub['FDR_Sig'].sum()} FDR")
    print()
    return res


def annotate_tumour_vs_control_scheme(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    out['Association_Scheme'] = 'ALT carrier vs WT/non-carrier'
    out['Comparison_Definition'] = out['Analysis_Group'].astype(str) + ' tumour vs pooled normal'
    out['Genotype_Definition'] = 'Carrier = Het or Hom ALT; non-carrier = WT'
    return out


def annotate_clinical_scheme(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    out['Association_Scheme'] = 'SNP × clinical variable association'
    out['Comparison_Definition'] = out['Cohort'].astype(str) + ' tumour cohort'
    out['Genotype_Definition'] = 'WT / Het / Hom (forced-genotype source)'
    return out


def annotate_genotype_dose_scheme(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    out['Association_Scheme'] = 'Genotype dose and genotype-class association'
    out['Comparison_Definition'] = out['Cohort'].astype(str) + ' tumour cohort'
    out['Genotype_Definition'] = 'WT / Het / Hom genotype model'
    return out


def annotate_cancer_risk_scheme(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    out['Association_Scheme'] = 'Cancer-risk association against pooled normal controls'
    out['Comparison_Definition'] = out['Cohort'].astype(str) + ' tumour vs pooled normal'
    out['Genotype_Definition'] = 'Carrier model plus WT / Het / Hom contrasts'
    return out


def annotate_survival_scheme(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    out["Association_Scheme"] = "Genotype survival association with Cox PH and KM"
    out["Comparison_Definition"] = out["Cohort"].astype(str) + " tumour cohort survival endpoint"
    out["Genotype_Definition"] = "WT / Het / Hom genotype model"
    return out




def _sample_manifest(merged: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in ['Sample', 'snp_code', 'sample_id', 'Cohort', 'Tissue', 'sheet', 'tumour_normal', 'is_replicate', 'canon__age', BMI_SOURCE_COLUMN, 'BMI_Category', 'BMI_Usable', 'canon__recurrence_flag', 'canon__distant_mets_flag', 'canon__exitus_flag', 'canon__pd_flag', 'canon__os_months', 'canon__pfs_months', 'canon__risk_group'] if c in merged.columns]
    return merged.drop_duplicates('Sample')[cols].copy().reset_index(drop=True)


def _variant_map(merged: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in ['Variant_ID', 'rsID', 'Variant_Label', 'Gene', 'SYMBOL', 'CHROM', 'POS', 'REF', 'ALT', 'Existing_variation'] if c in merged.columns]
    out = merged.drop_duplicates('Variant_ID')[cols].copy()
    if 'SYMBOL' in out.columns and 'Gene' not in out.columns:
        out = out.rename(columns={'SYMBOL': 'Gene'})
    if 'rsID' not in out.columns:
        out['rsID'] = out.apply(_extract_rsid, axis=1)
    out['rsID'] = out['rsID'].replace('', np.nan)
    if 'Variant_Label' not in out.columns:
        out['Variant_Label'] = out.apply(_preferred_variant_id, axis=1)
    else:
        out['Variant_Label'] = out['Variant_Label'].astype(str).str.strip()
        out.loc[out['Variant_Label'].isin(['', 'nan', 'None', '<NA>']), 'Variant_Label'] = np.nan
        out['Variant_Label'] = out['Variant_Label'].fillna(out.apply(_preferred_variant_id, axis=1))
    out['Variant_Key'] = out['CHROM'].astype(str).str.strip() + ':' + pd.to_numeric(out['POS'], errors='coerce').fillna(-1).astype(int).astype(str) + ':' + out['REF'].astype(str).str.strip().str.upper() + ':' + out['ALT'].astype(str).str.strip().str.upper()
    return out


def _attach_variant_metadata(df: pd.DataFrame, merged: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty or 'Variant_ID' not in df.columns:
        return df
    lookup = _variant_map(merged)[['Variant_ID', 'rsID', 'Variant_Label', 'Gene']].drop_duplicates('Variant_ID')
    out = df.merge(lookup, on='Variant_ID', how='left', suffixes=('', '__lookup'))
    if 'Gene__lookup' in out.columns:
        if 'Gene' in out.columns:
            out['Gene'] = out['Gene'].fillna(out['Gene__lookup'])
        else:
            out = out.rename(columns={'Gene__lookup': 'Gene'})
        out = out.drop(columns=['Gene__lookup'], errors='ignore')
    if 'rsID' not in out.columns:
        out['rsID'] = out['Variant_ID']
    out['rsID'] = out['rsID'].replace('', np.nan)
    out['rsID'] = out['rsID'].fillna(out['Variant_ID'].where(out['Variant_ID'].astype(str).str.match(r'^rs\d+$', case=False), np.nan))
    out['Variant_Label'] = out.get('Variant_Label', pd.Series(np.nan, index=out.index))
    out['Variant_Label'] = out['Variant_Label'].fillna(out['rsID']).fillna(out['Variant_ID'])
    preferred = [c for c in ['Variant_ID', 'rsID', 'Variant_Label'] if c in out.columns]
    remaining = [c for c in out.columns if c not in preferred]
    return out[preferred + remaining]


def _cohort_display_name(cohort: str) -> str:
    mapping = {
        "Breast_Tumour": "Breast tumour",
        "Breast_Healthy": "Breast healthy",
        "Endometrial_Tumour": "Endometrial tumour",
        "Endometrial_Healthy": "Endometrial healthy",
    }
    return mapping.get(str(cohort), str(cohort).replace("_", " ").strip())


def _cohort_family_label(cohort: str) -> str:
    text = str(cohort)
    if "Breast" in text:
        return "Breast"
    if "Endometri" in text:
        return "Endometrial"
    return text


def _cohort_sort_key(cohort: str) -> int:
    try:
        return BMI_COHORT_ORDER.index(str(cohort))
    except ValueError:
        return len(BMI_COHORT_ORDER)


def build_bmi_category_summaries(merged: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    manifest = _sample_manifest(merged)
    if manifest.empty or BMI_SOURCE_COLUMN not in manifest.columns:
        return pd.DataFrame(), pd.DataFrame()

    manifest = manifest.copy()
    manifest[BMI_SOURCE_COLUMN] = pd.to_numeric(manifest[BMI_SOURCE_COLUMN], errors="coerce")
    if "BMI_Category" not in manifest.columns:
        manifest["BMI_Category"] = manifest[BMI_SOURCE_COLUMN].apply(_categorize_bmi_value)
    manifest["BMI_Usable"] = manifest["BMI_Category"].notna()
    manifest["Cohort_Display"] = manifest["Cohort"].map(_cohort_display_name)
    manifest["Cohort_Family"] = manifest["Cohort"].map(_cohort_family_label)

    availability_rows = []
    summary_rows = []
    for cohort in manifest["Cohort"].dropna().astype(str).unique():
        sub = manifest[manifest["Cohort"].astype(str) == cohort].copy()
        total_samples = int(len(sub))
        numeric_bmi = sub[BMI_SOURCE_COLUMN].notna()
        usable_mask = sub["BMI_Usable"].fillna(False)
        usable_total = int(usable_mask.sum())
        under20_total = int((numeric_bmi & sub[BMI_SOURCE_COLUMN].lt(20)).sum())
        missing_total = int(total_samples - int(numeric_bmi.sum()))
        availability_rows.append(
            {
                "Cohort": cohort,
                "Cohort_Display": _cohort_display_name(cohort),
                "Tissue": _first_nonempty_text(pd.Series(sub.get("Tissue", pd.Series(dtype=object)))),
                "BMI_Source_Column": BMI_SOURCE_COLUMN,
                "Total_Samples": total_samples,
                "Samples_With_Numeric_BMI": int(numeric_bmi.sum()),
                "Samples_With_Usable_BMI_Category": usable_total,
                "Samples_Below_20_Excluded": under20_total,
                "Samples_Missing_or_NonNumeric_BMI": missing_total,
            }
        )
        for category in BMI_CATEGORY_ORDER:
            count = int((sub["BMI_Category"] == category).sum())
            pct = round((count / usable_total) * 100, 2) if usable_total else np.nan
            summary_rows.append(
                {
                    "Cohort": cohort,
                    "Cohort_Display": _cohort_display_name(cohort),
                    "Tissue": _first_nonempty_text(pd.Series(sub.get("Tissue", pd.Series(dtype=object)))),
                    "BMI_Source_Column": BMI_SOURCE_COLUMN,
                    "Total_Samples_With_Usable_BMI": usable_total,
                    "Category": category,
                    "Count": count,
                    "Percentage": pct,
                }
            )

    availability = pd.DataFrame(availability_rows)
    availability["_sort_order"] = availability["Cohort"].map(_cohort_sort_key)
    availability = availability.sort_values(["_sort_order", "Cohort_Display"], kind="stable").drop(columns="_sort_order").reset_index(drop=True)
    summary = pd.DataFrame(summary_rows)
    summary["Category"] = pd.Categorical(summary["Category"], categories=BMI_CATEGORY_ORDER, ordered=True)
    summary["_sort_order"] = summary["Cohort"].map(_cohort_sort_key)
    summary = summary.sort_values(["_sort_order", "Category"], kind="stable").drop(columns="_sort_order").reset_index(drop=True)
    return summary, availability


def build_bmi_by_tissue_summary(bmi_summary: pd.DataFrame) -> pd.DataFrame:
    if bmi_summary.empty:
        return pd.DataFrame()

    work = bmi_summary.copy()
    work["Cohort_Family"] = work["Cohort"].map(_cohort_family_label)
    grouped = (
        work.groupby(["Cohort_Family", "Tissue", "Category"], dropna=False, observed=False)["Count"]
        .sum()
        .reset_index()
    )
    totals = (
        grouped.groupby(["Cohort_Family", "Tissue"], dropna=False)["Count"]
        .sum()
        .rename("Total_Samples_With_Usable_BMI")
        .reset_index()
    )
    grouped = grouped.merge(totals, on=["Cohort_Family", "Tissue"], how="left")
    grouped["Percentage"] = np.where(
        grouped["Total_Samples_With_Usable_BMI"] > 0,
        (grouped["Count"] / grouped["Total_Samples_With_Usable_BMI"] * 100).round(2),
        np.nan,
    )
    grouped["BMI_Source_Column"] = BMI_SOURCE_COLUMN
    grouped["Category"] = pd.Categorical(grouped["Category"], categories=BMI_CATEGORY_ORDER, ordered=True)
    return grouped.sort_values(["Cohort_Family", "Tissue", "Category"], kind="stable").reset_index(drop=True)


@lru_cache(maxsize=1)
def _load_forced_genotype_calls() -> pd.DataFrame:
    forced_path = Path('/home/gadeaalonsoj/tfm/analysis_results/05b_forced_genotypes/GSDMB_Forced_Genotypes_Union_Sites.xlsx')
    if not forced_path.exists():
        return pd.DataFrame()
    df = pd.read_excel(forced_path, sheet_name='Per_Sample_Genotypes', engine='openpyxl')
    df.columns = [str(c).strip() for c in df.columns]
    df['Variant_Key'] = df['CHROM'].astype(str).str.strip() + ':' + pd.to_numeric(df['POS'], errors='coerce').fillna(-1).astype(int).astype(str) + ':' + df['REF'].astype(str).str.strip().str.upper() + ':' + df['Union_ALT_List'].astype(str).str.strip().str.upper()
    df['snp_code'] = df['Sample'].map(_extract_snp_code)
    state_map = {'WT': 'WT', 'Het_ALT': 'Het', 'Hom_ALT': 'Hom'}
    df['Genotype_State'] = df['Genotype_Status'].map(state_map)
    df['Genotype_Dose'] = df['Genotype_State'].map({'WT': 0, 'Het': 1, 'Hom': 2})
    return df


def _fit_glm_binomial(df: pd.DataFrame, y_col: str, x_cols: List[str]) -> Tuple[float, float, float, float]:
    data = df[[y_col] + x_cols].dropna().copy()
    if data.empty or data[y_col].nunique() < 2:
        return np.nan, np.nan, np.nan, np.nan
    try:
        X = sm.add_constant(data[x_cols], has_constant='add')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PerfectSeparationWarning)
            res = sm.GLM(data[y_col].astype(float), X, family=sm.families.Binomial()).fit()
        coef = float(res.params[x_cols[0]])
        ci = res.conf_int().loc[x_cols[0]]
        coef = float(np.clip(coef, -20, 20))
        lo = float(np.clip(ci.iloc[0], -20, 20))
        hi = float(np.clip(ci.iloc[1], -20, 20))
        return float(np.exp(coef)), float(np.exp(lo)), float(np.exp(hi)), float(res.pvalues[x_cols[0]])
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


def _variant_sample_table(merged: pd.DataFrame, variant_id: str, sample_subset: Optional[pd.Series] = None) -> pd.DataFrame:
    forced = _load_forced_genotype_calls()
    if forced.empty:
        return pd.DataFrame()
    vmap = _variant_map(merged)
    vmap = vmap[vmap['Variant_ID'].astype(str) == str(variant_id)].copy()
    if vmap.empty:
        return pd.DataFrame()
    out = forced.merge(vmap[['Variant_ID', 'Variant_Key', 'Gene']], on='Variant_Key', how='inner')
    if sample_subset is not None:
        keep = set(pd.Series(sample_subset).astype(str))
        keep_codes = {k for k in (_extract_snp_code(v) for v in keep) if k}
        out = out[out['Sample'].astype(str).isin(keep) | out['snp_code'].astype(str).isin(keep) | out['snp_code'].astype(str).isin(keep_codes)].copy()
    return out


def clinical_associations(merged: pd.DataFrame):
    sample_manifest = _sample_manifest(merged)
    if sample_manifest.empty:
        return pd.DataFrame(), pd.DataFrame()
    results = []
    for cohort_label, var_map in [('Breast_Tumour', CLINICAL_VARS_BREAST), ('Endometrial_Tumour', CLINICAL_VARS_ENDO)]:
        cohort_samples = sample_manifest[sample_manifest['Cohort'].astype(str) == cohort_label].copy()
        if cohort_samples.empty:
            continue
        cohort_calls = merged[merged['Sample'].astype(str).isin(set(cohort_samples['Sample'].astype(str)))].copy()
        for var_id in pd.Series(cohort_calls['Variant_ID'].dropna().unique()).astype(str):
            call_rows = cohort_calls[cohort_calls['Variant_ID'].astype(str) == var_id].drop_duplicates('Sample').copy()
            if call_rows.empty:
                continue
            geno_rows = _variant_sample_table(merged, var_id, cohort_samples['Sample'])
            if geno_rows.empty:
                continue
            gene = _first_nonempty_text(geno_rows.get('Gene', pd.Series(dtype=object)))
            geno_counts = geno_rows['Genotype_State'].value_counts()
            n_wt, n_het, n_hom = int(geno_counts.get('WT', 0)), int(geno_counts.get('Het', 0)), int(geno_counts.get('Hom', 0))
            join_key = 'snp_code' if 'snp_code' in cohort_samples.columns and 'snp_code' in geno_rows.columns else 'Sample'
            sample_frame = cohort_samples.copy().merge(geno_rows[[join_key, 'Genotype_State', 'Genotype_Dose']], on=join_key, how='left')
            sample_frame['carrier'] = sample_frame['Genotype_State'].isin(['Het', 'Hom']).astype(int)
            for clin_var, meta in var_map.items():
                if clin_var not in sample_frame.columns:
                    continue
                y = _binary_yes_no(sample_frame[clin_var])
                cols = ['carrier', 'Genotype_Dose']
                if 'canon__age' in sample_frame.columns:
                    cols.append('canon__age')
                if 'canon__bmi' in sample_frame.columns:
                    cols.append('canon__bmi')
                valid = sample_frame[cols].copy()
                valid['y'] = y
                valid = valid.dropna(subset=['y']).copy()
                if valid.empty or valid['y'].nunique() < 2:
                    continue
                or_unadj, lo_unadj, hi_unadj, p_unadj = _fit_glm_binomial(valid, 'y', ['carrier'])
                or_age, lo_age, hi_age, p_age = _fit_glm_binomial(valid.assign(age=pd.to_numeric(valid.get('canon__age'), errors='coerce')), 'y', ['carrier', 'age']) if 'canon__age' in valid.columns else (np.nan, np.nan, np.nan, np.nan)
                or_bmi, lo_bmi, hi_bmi, p_bmi = _fit_glm_binomial(valid.assign(bmi=pd.to_numeric(valid.get('canon__bmi'), errors='coerce')), 'y', ['carrier', 'bmi']) if 'canon__bmi' in valid.columns else (np.nan, np.nan, np.nan, np.nan)
                or_ab, lo_ab, hi_ab, p_ab = _fit_glm_binomial(valid.assign(age=pd.to_numeric(valid.get('canon__age'), errors='coerce'), bmi=pd.to_numeric(valid.get('canon__bmi'), errors='coerce')), 'y', ['carrier', 'age', 'bmi']) if 'canon__age' in valid.columns and 'canon__bmi' in valid.columns else (np.nan, np.nan, np.nan, np.nan)
                add_or, add_lo, add_hi, add_p = _fit_glm_binomial(valid, 'y', ['Genotype_Dose'])
                results.append({'Cohort': cohort_label, 'Variant_ID': var_id, 'Gene': gene, 'Contrast': 'Carrier_vs_WT', 'N_WT': n_wt, 'N_Het': n_het, 'N_Hom': n_hom, 'Clinical_Var': clin_var, 'Clin_Label': meta['label'], 'Clin_Type': meta['type'], 'Note': meta.get('note', ''), 'Test_Unadj': 'GLM binomial', 'N_Carriers': int((valid['carrier'] == 1).sum()), 'N_NonCarriers': int((valid['carrier'] == 0).sum()), 'Events_Carriers': int(((valid['carrier'] == 1) & (valid['y'] == 1)).sum()), 'Events_NonCarriers': int(((valid['carrier'] == 0) & (valid['y'] == 1)).sum()), 'OR_Unadj': or_unadj, 'OR_Unadj_CI95': _format_ci95(lo_unadj, hi_unadj), 'OR_Unadj_Method': 'GLM binomial', 'P_Unadj': p_unadj, 'Test_Adj_Age': 'GLM binomial', 'OR_Adj_Age': or_age, 'OR_Adj_Age_CI95': _format_ci95(lo_age, hi_age), 'P_Adj_Age': p_age, 'N_Adj_Age': int(valid[['carrier', 'canon__age']].dropna().shape[0]) if 'canon__age' in valid.columns else np.nan, 'Adj_Age_Note': 'Adjusted for age' if 'canon__age' in valid.columns else 'Age unavailable', 'Test_Adj_BMI': 'GLM binomial', 'OR_Adj_BMI': or_bmi, 'OR_Adj_BMI_CI95': _format_ci95(lo_bmi, hi_bmi), 'P_Adj_BMI': p_bmi, 'N_Adj_BMI': int(valid[['carrier', 'canon__bmi']].dropna().shape[0]) if 'canon__bmi' in valid.columns else np.nan, 'Adj_BMI_Note': 'Adjusted for BMI' if 'canon__bmi' in valid.columns else 'BMI unavailable', 'Test_Adj_AgeBMI': 'GLM binomial', 'OR_Adj_AgeBMI': or_ab, 'OR_Adj_AgeBMI_CI95': _format_ci95(lo_ab, hi_ab), 'P_Adj_AgeBMI': p_ab, 'N_Adj_AgeBMI': int(valid[['carrier', 'canon__age', 'canon__bmi']].dropna().shape[0]) if 'canon__age' in valid.columns and 'canon__bmi' in valid.columns else np.nan, 'Adj_AgeBMI_Note': 'Adjusted for age + BMI' if 'canon__age' in valid.columns and 'canon__bmi' in valid.columns else 'Age/BMI unavailable', 'Test_Unadj_Additive': 'GLM binomial', 'N_Carriers_Additive': int((valid['Genotype_Dose'] > 0).sum()), 'N_NonCarriers_Additive': int((valid['Genotype_Dose'] == 0).sum()), 'P_Unadj_Additive': add_p, 'OR_Unadj_Additive': add_or, 'OR_Unadj_Additive_CI95': _format_ci95(add_lo, add_hi), 'OR_Adj_Age_Additive': np.nan, 'P_Adj_Age_Additive': np.nan, 'OR_Adj_BMI_Additive': np.nan, 'P_Adj_BMI_Additive': np.nan, 'OR_Adj_AgeBMI_Additive': np.nan, 'P_Adj_AgeBMI_Additive': np.nan})
    if not results:
        return pd.DataFrame(), pd.DataFrame()
    out = pd.DataFrame(results)
    out = _apply_fdr(out, 'P_Unadj', 'FDR_Unadj')
    out = _apply_fdr(out, 'P_Adj_Age', 'FDR_Adj_Age')
    out = _apply_fdr(out, 'P_Adj_BMI', 'FDR_Adj_BMI')
    out = _apply_fdr(out, 'P_Adj_AgeBMI', 'FDR_Adj_AgeBMI')
    out = _apply_fdr(out, 'P_Unadj_Additive', 'FDR_Unadj_Additive')
    out['Nominal_Sig_Unadj'] = out['P_Unadj'] < 0.05
    out['FDR_Sig_Unadj'] = out['FDR_Unadj'] < FDR_THRESHOLD
    out['Nominal_Sig_Adj_Age'] = out['P_Adj_Age'] < 0.05
    out['FDR_Sig_Adj_Age'] = out['FDR_Adj_Age'] < FDR_THRESHOLD
    out['Nominal_Sig_Adj_BMI'] = out['P_Adj_BMI'] < 0.05
    out['FDR_Sig_Adj_BMI'] = out['FDR_Adj_BMI'] < FDR_THRESHOLD
    out['Nominal_Sig_Adj_AgeBMI'] = out['P_Adj_AgeBMI'] < 0.05
    out['FDR_Sig_Adj_AgeBMI'] = out['FDR_Adj_AgeBMI'] < FDR_THRESHOLD
    out['Primary_Row_Model'] = out['Test_Unadj']
    out['Secondary_Additive_Model'] = out['Test_Unadj_Additive']
    breast = annotate_clinical_scheme(out[out['Cohort'] == 'Breast_Tumour'].copy())
    endo = annotate_clinical_scheme(out[out['Cohort'] == 'Endometrial_Tumour'].copy())
    return breast, endo


def exploratory_breast_treatment_analysis(merged: pd.DataFrame):
    sample_manifest = _sample_manifest(merged)
    breast = sample_manifest[sample_manifest['Cohort'].astype(str) == 'Breast_Tumour'].copy()
    if breast.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    rows = []
    for col, meta in BREAST_TREATMENT_EXPOSURES.items():
        if col in breast.columns:
            rows.append({'Treatment_Exposure_Var': col, 'Treatment_Exposure': meta['label'], 'N_Exposed': int(_binary_yes_no(breast[col]).sum()), 'N_With_OS': int(breast.get('canon__os_months', pd.Series(dtype=float)).notna().sum()) if 'canon__os_months' in breast.columns else np.nan, 'N_Recurrence_Events': int(_binary_yes_no(breast.get('canon__recurrence_flag', pd.Series(dtype=object))).sum()) if 'canon__recurrence_flag' in breast.columns else np.nan, 'N_Metastasis_Events': int(_binary_yes_no(breast.get('canon__distant_mets_flag', pd.Series(dtype=object))).sum()) if 'canon__distant_mets_flag' in breast.columns else np.nan, 'N_Death_Events': int(_binary_yes_no(breast.get('canon__exitus_flag', pd.Series(dtype=object))).sum()) if 'canon__exitus_flag' in breast.columns else np.nan})
    return pd.DataFrame(rows), pd.DataFrame(), pd.DataFrame()


def genotype_dose_analysis(merged: pd.DataFrame, variant_list: List[str], clinical_vars: Dict[str, Dict], cohort_label: str) -> pd.DataFrame:
    sample_manifest = _sample_manifest(merged)
    cohort_samples = sample_manifest[sample_manifest['Cohort'].astype(str) == cohort_label].copy()
    if cohort_samples.empty:
        return pd.DataFrame()
    rows = []
    for var_id in pd.Series(variant_list).dropna().astype(str).unique():
        geno = _variant_sample_table(merged, var_id, cohort_samples['Sample'])
        if geno.empty:
            continue
        gene = _first_nonempty_text(geno.get('Gene', pd.Series(dtype=object)))
        counts = geno['Genotype_State'].value_counts()
        n_wt, n_het, n_hom = int(counts.get('WT', 0)), int(counts.get('Het', 0)), int(counts.get('Hom', 0))
        join_key = 'snp_code' if 'snp_code' in cohort_samples.columns and 'snp_code' in geno.columns else 'Sample'
        frame = cohort_samples.copy().merge(geno[[join_key, 'Genotype_State', 'Genotype_Dose']], on=join_key, how='left')
        frame = frame[frame['Genotype_State'].isin(['WT', 'Het', 'Hom'])].copy()
        frame['Het'] = (frame['Genotype_State'] == 'Het').astype(int)
        frame['Hom'] = (frame['Genotype_State'] == 'Hom').astype(int)
        for clin_var, meta in clinical_vars.items():
            if clin_var not in frame.columns:
                continue
            y = _binary_yes_no(frame[clin_var])
            valid = frame.copy(); valid['y'] = y; valid = valid.dropna(subset=['y'])
            if valid.empty or valid['y'].nunique() < 2:
                continue
            _, _, _, p_tr = _fit_glm_binomial(valid, 'y', ['Genotype_Dose'])
            het = valid[valid['Genotype_State'].isin(['WT', 'Het'])].copy(); het['x'] = (het['Genotype_State'] == 'Het').astype(int)
            hom = valid[valid['Genotype_State'].isin(['WT', 'Hom'])].copy(); hom['x'] = (hom['Genotype_State'] == 'Hom').astype(int)
            _, _, _, p_het = _fit_glm_binomial(het, 'y', ['x'])
            _, _, _, p_hom = _fit_glm_binomial(hom, 'y', ['x'])
            rows.append({'Cohort': cohort_label, 'Variant_ID': var_id, 'Gene': gene, 'Clinical_Var': clin_var, 'Clin_Label': meta['label'], 'N_WT': n_wt, 'N_Het': n_het, 'N_Hom': n_hom, 'Test_Trend': 'GLM binomial', 'Stat_Trend': np.nan, 'P_Trend': p_tr, 'Test_Het_vs_WT': 'GLM binomial', 'P_Het_vs_WT': p_het, 'N_Het_used': int((het['Genotype_State'] == 'Het').sum()), 'Test_Hom_vs_WT': 'GLM binomial', 'P_Hom_vs_WT': p_hom, 'N_Hom_used': int((hom['Genotype_State'] == 'Hom').sum()), 'FDR_Trend': np.nan, 'Nominal_Sig': p_tr < 0.05, 'FDR_Sig': False, 'FDR_Het_vs_WT': np.nan, 'FDR_Hom_vs_WT': np.nan, 'P_Value': p_tr, 'FDR_P_Value': np.nan, 'Coverage_Risk_Flag': False, 'Coverage_Risk_Amplicon': '', 'Coverage_Risk_Region': '', 'Coverage_Risk_Note': '', 'Trend_Model': 'GLM binomial', 'Het_vs_WT_Model': 'GLM binomial', 'Hom_vs_WT_Model': 'GLM binomial', 'Genotype_Definition': 'WT / Het / Hom genotype model'})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = _apply_fdr(out, 'P_Trend', 'FDR_Trend')
    out = _apply_fdr(out, 'P_Het_vs_WT', 'FDR_Het_vs_WT')
    out = _apply_fdr(out, 'P_Hom_vs_WT', 'FDR_Hom_vs_WT')
    out['FDR_Sig'] = out['FDR_Trend'] < FDR_THRESHOLD
    return annotate_genotype_dose_scheme(out)


def survival_analysis(merged: pd.DataFrame):
    if not _HAS_LIFELINES:
        return pd.DataFrame(), []

    sample_manifest = _sample_manifest(merged)
    sample_manifest = sample_manifest[~sample_manifest["is_replicate"].fillna(False)].copy()
    if sample_manifest.empty:
        return pd.DataFrame(), []

    endpoint_defs = [
        {
            "Cohort": "Breast",
            "mask": (sample_manifest["Tissue"].astype(str) == "Tumour") & sample_manifest["Cohort"].astype(str).str.contains("Breast", case=False, na=False),
            "Endpoint": "Overall survival",
            "Duration_Col": "canon__os_months",
            "Event_Col": "canon__exitus_flag",
        },
        {
            "Cohort": "Endometrium",
            "mask": (sample_manifest["Tissue"].astype(str) == "Tumour") & sample_manifest["Cohort"].astype(str).str.contains("Endometri", case=False, na=False),
            "Endpoint": "Overall survival",
            "Duration_Col": "canon__os_months",
            "Event_Col": "canon__exitus_flag",
        },
        {
            "Cohort": "Endometrium",
            "mask": (sample_manifest["Tissue"].astype(str) == "Tumour") & sample_manifest["Cohort"].astype(str).str.contains("Endometri", case=False, na=False),
            "Endpoint": "Progression-free survival",
            "Duration_Col": "canon__pfs_months",
            "Event_Col": "canon__pd_flag",
        },
    ]
    model_defs = [
        {"Model_Label": "Genotype only", "Model_Slug": "genotype_only", "Plot_Mode": "genotype", "Covariates": ["Genotype_Dose"]},
        {"Model_Label": "Genotype, age adjusted", "Model_Slug": "genotype_age", "Plot_Mode": "genotype", "Covariates": ["Genotype_Dose", "canon__age"]},
        {"Model_Label": "Genotype, BMI adjusted", "Model_Slug": "genotype_bmi", "Plot_Mode": "genotype", "Covariates": ["Genotype_Dose", "canon__bmi"]},
        {"Model_Label": "Genotype, age + BMI adjusted", "Model_Slug": "genotype_age_bmi", "Plot_Mode": "genotype", "Covariates": ["Genotype_Dose", "canon__age", "canon__bmi"]},
        {"Model_Label": "Carrier vs non-carrier", "Model_Slug": "carrier_only", "Plot_Mode": "carrier", "Covariates": ["Carrier"]},
        {"Model_Label": "Carrier, age adjusted", "Model_Slug": "carrier_age", "Plot_Mode": "carrier", "Covariates": ["Carrier", "canon__age"]},
        {"Model_Label": "Carrier, BMI adjusted", "Model_Slug": "carrier_bmi", "Plot_Mode": "carrier", "Covariates": ["Carrier", "canon__bmi"]},
        {"Model_Label": "Carrier, age + BMI adjusted", "Model_Slug": "carrier_age_bmi", "Plot_Mode": "carrier", "Covariates": ["Carrier", "canon__age", "canon__bmi"]},
    ]

    variant_ids = _variant_map(merged)["Variant_ID"].dropna().astype(str).unique().tolist()
    rows = []
    km_store = {}

    for cfg in endpoint_defs:
        tumour_samples = sample_manifest.loc[cfg["mask"]].copy()
        if tumour_samples.empty:
            continue
        base_cols = ["Sample", "snp_code", cfg["Duration_Col"], cfg["Event_Col"]]
        if "canon__age" in tumour_samples.columns:
            base_cols.append("canon__age")
        if "canon__bmi" in tumour_samples.columns:
            base_cols.append("canon__bmi")
        base = tumour_samples[base_cols].copy()
        base = base.rename(columns={cfg["Duration_Col"]: "Duration", cfg["Event_Col"]: "Event"})
        base["Duration"] = pd.to_numeric(base["Duration"], errors="coerce")
        base["Event"] = pd.to_numeric(base["Event"], errors="coerce")
        base = base.dropna(subset=["Duration", "Event"]).copy()
        base["Event"] = base["Event"].astype(int)
        if base.empty:
            continue

        for var_id in variant_ids:
            geno = _variant_sample_table(merged, var_id, tumour_samples["Sample"])
            if geno.empty:
                continue
            join_key = "snp_code" if "snp_code" in base.columns and "snp_code" in geno.columns else "Sample"
            frame = base.merge(geno[[join_key, "Genotype_State", "Genotype_Dose"]], on=join_key, how="left")
            frame = frame[frame["Genotype_State"].isin(["WT", "Het", "Hom"])].copy()
            if frame.empty or frame["Genotype_State"].nunique() < 2:
                continue
            frame["Carrier"] = frame["Genotype_State"].isin(["Het", "Hom"]).astype(int)
            frame["Carrier_State"] = np.where(frame["Carrier"] == 1, "Carrier", "WT")

            try:
                lr = multivariate_logrank_test(frame["Duration"], frame["Genotype_State"], frame["Event"])
                p_logrank = float(lr.p_value)
            except Exception:
                p_logrank = np.nan

            gene = _first_nonempty_text(geno.get("Gene", pd.Series(dtype=object)))
            n_wt = int((frame["Genotype_State"] == "WT").sum())
            n_het = int((frame["Genotype_State"] == "Het").sum())
            n_hom = int((frame["Genotype_State"] == "Hom").sum())

            for model in model_defs:
                covariates = [col for col in model["Covariates"] if col in frame.columns]
                if len(covariates) != len(model["Covariates"]):
                    continue
                hr, lo, hi, p_cox, n_used, n_events, n_levels = _fit_cox_survival(frame, "Duration", "Event", covariates)
                if pd.isna(p_cox):
                    continue
                rows.append({
                    "Analysis_Type": "Survival",
                    "Cohort": cfg["Cohort"],
                    "Endpoint": cfg["Endpoint"],
                    "Model": model["Model_Label"],
                    "Model_Slug": model["Model_Slug"],
                    "Variant_ID": var_id,
                    "Gene": gene,
                    "Clin_Label": cfg["Endpoint"],
                    "Clin_Type": "survival",
                    "N": n_used,
                    "Events": n_events,
                    "N_WT": n_wt,
                    "N_Het": n_het,
                    "N_Hom": n_hom,
                    "P_Cox_Additive": p_cox,
                    "HR_Cox_Additive": hr,
                    "HR_Cox_Additive_CI95": f"[{lo:.3g}, {hi:.3g}]" if pd.notna(lo) else np.nan,
                    "Cox_Covariates": " + ".join(covariates),
                    "P_LogRank": p_logrank,
                    "Nominal_Sig_Cox_Additive": p_cox < 0.05,
                    "FDR_Sig_Cox_Additive": False,
                    "Nominal_Sig_LogRank": p_logrank < 0.05 if pd.notna(p_logrank) else False,
                    "FDR_Sig_LogRank": False,
                    "Association_Scheme": "Genotype survival association with Cox PH and KM",
                    "Comparison_Definition": f"{cfg['Cohort']} tumour cohort survival endpoint",
                    "Genotype_Definition": "WT / Het / Hom genotype model",
                    "Duration_Col": cfg["Duration_Col"],
                    "Event_Col": cfg["Event_Col"],
                })
                km_store[(cfg["Cohort"], cfg["Endpoint"], var_id, model["Model_Slug"])] = frame.copy()

    out = pd.DataFrame(rows)
    if out.empty:
        return out, []

    parts = []
    for (cohort, endpoint, model_slug), sub in out.groupby(["Cohort", "Endpoint", "Model_Slug"], dropna=False):
        sub = sub.copy()
        sub = _apply_fdr(sub, "P_Cox_Additive", "FDR_Cox_Additive")
        sub = _apply_fdr(sub, "P_LogRank", "FDR_LogRank")
        sub["Nominal_Sig_Cox_Additive"] = sub["P_Cox_Additive"] < 0.05
        sub["FDR_Sig_Cox_Additive"] = sub["FDR_Cox_Additive"] < FDR_THRESHOLD
        sub["Nominal_Sig_LogRank"] = sub["P_LogRank"] < 0.05 if "P_LogRank" in sub.columns else False
        sub["FDR_Sig_LogRank"] = sub["FDR_LogRank"] < FDR_THRESHOLD if "FDR_LogRank" in sub.columns else False
        parts.append(sub)
    out = pd.concat(parts, ignore_index=True).sort_values(["Cohort", "Endpoint", "Model_Slug", "P_Cox_Additive", "Variant_ID"], kind="stable")
    out = annotate_survival_scheme(out)

    km_pages = []
    for (cohort, endpoint, model_slug), sub in out.groupby(["Cohort", "Endpoint", "Model_Slug"], dropna=False):
        chosen = sub.sort_values(["FDR_Sig_Cox_Additive", "Nominal_Sig_Cox_Additive", "P_Cox_Additive"], ascending=[False, False, True]).head(3)
        for _, row in chosen.iterrows():
            frame = km_store.get((cohort, endpoint, row["Variant_ID"], model_slug))
            if frame is None or frame.empty:
                continue
            km_pages.append({
                "Cohort": cohort,
                "Endpoint": endpoint,
                "Model": row.get("Model", ""),
                "Model_Slug": model_slug,
                "Variant_ID": row["Variant_ID"],
                "Gene": row.get("Gene", ""),
                "Frame": frame,
                "P_Cox_Additive": row.get("P_Cox_Additive", np.nan),
                "FDR_Cox_Additive": row.get("FDR_Cox_Additive", np.nan),
                "P_LogRank": row.get("P_LogRank", np.nan),
                "HR_Cox_Additive": row.get("HR_Cox_Additive", np.nan),
            })

    return out, km_pages


def cancer_risk_analysis(merged: pd.DataFrame) -> pd.DataFrame:
    sample_manifest = _sample_manifest(merged)
    rows = []
    for cohort, cohort_match in [('Breast', 'Breast'), ('Endometrium', 'Endometri')]:
        cases = sample_manifest[(sample_manifest['Tissue'].astype(str) == 'Tumour') & (sample_manifest['Cohort'].astype(str).str.contains(cohort_match, case=False, na=False))].copy()
        ctrls = sample_manifest[sample_manifest['Tissue'].astype(str) == 'Healthy'].copy()
        if cases.empty or ctrls.empty:
            continue
        for var_id in _variant_map(merged)['Variant_ID'].astype(str).dropna().unique():
            case_geno = _variant_sample_table(merged, var_id, cases['Sample'])
            ctl_geno = _variant_sample_table(merged, var_id, ctrls['Sample'])
            if case_geno.empty or ctl_geno.empty:
                continue
            gene = _first_nonempty_text(case_geno.get('Gene', pd.Series(dtype=object)))
            n_cases, n_ctrl = int(cases.shape[0]), int(ctrls.shape[0])
            a = int(case_geno['Genotype_State'].isin(['Het', 'Hom']).sum())
            c = int(ctl_geno['Genotype_State'].isin(['Het', 'Hom']).sum())
            b, d = n_cases - a, n_ctrl - c
            geno = pd.concat([case_geno.assign(outcome=1), ctl_geno.assign(outcome=0)], ignore_index=True)
            geno = geno[geno['Genotype_State'].isin(['WT', 'Het', 'Hom'])].copy(); geno['carrier'] = geno['Genotype_State'].isin(['Het', 'Hom']).astype(int); geno['dose'] = geno['Genotype_Dose']
            or_unadj, lo_unadj, hi_unadj, p_unadj = _fit_glm_binomial(geno, 'outcome', ['carrier'])
            or_age, lo_age, hi_age, p_age = _fit_glm_binomial(geno.assign(age=pd.to_numeric(geno.get('canon__age'), errors='coerce')), 'outcome', ['carrier', 'age']) if 'canon__age' in geno.columns else (np.nan, np.nan, np.nan, np.nan)
            or_ab, lo_ab, hi_ab, p_ab = _fit_glm_binomial(geno.assign(age=pd.to_numeric(geno.get('canon__age'), errors='coerce'), bmi=pd.to_numeric(geno.get('canon__bmi'), errors='coerce')), 'outcome', ['carrier', 'age', 'bmi']) if 'canon__age' in geno.columns and 'canon__bmi' in geno.columns else (np.nan, np.nan, np.nan, np.nan)
            het = geno[geno['Genotype_State'].isin(['WT', 'Het'])].copy(); het['x'] = (het['Genotype_State'] == 'Het').astype(int)
            hom = geno[geno['Genotype_State'].isin(['WT', 'Hom'])].copy(); hom['x'] = (hom['Genotype_State'] == 'Hom').astype(int)
            or_het, lo_het, hi_het, p_het = _fit_glm_binomial(het, 'outcome', ['x'])
            or_hom, lo_hom, hi_hom, p_hom = _fit_glm_binomial(hom, 'outcome', ['x'])
            or_tr, lo_tr, hi_tr, p_tr = _fit_glm_binomial(geno, 'outcome', ['dose'])
            rows.append({'Cohort': cohort, 'Variant_ID': var_id, 'Gene': gene, 'N_Cases': n_cases, 'N_Controls': n_ctrl, 'Carriers_Cases': a, 'Carriers_Controls': c, 'Freq_Cases_%': round(a / n_cases * 100, 2), 'Freq_Controls_%': round(c / n_ctrl * 100, 2), 'OR_Unadj': or_unadj, 'OR_Unadj_CI95': _format_ci95(lo_unadj, hi_unadj), 'OR_Unadj_Method': 'GLM binomial', 'P_Unadj': p_unadj, 'OR_Adj_Age': or_age, 'OR_Adj_Age_CI95': _format_ci95(lo_age, hi_age), 'P_Adj_Age': p_age, 'N_Adj_Age': int(geno[['carrier', 'canon__age']].dropna().shape[0]) if 'canon__age' in geno.columns else np.nan, 'Note_Age': 'Adjusted for age' if 'canon__age' in geno.columns else 'Age unavailable', 'OR_Adj_AgeBMI': or_ab, 'OR_Adj_AgeBMI_CI95': _format_ci95(lo_ab, hi_ab), 'P_Adj_AgeBMI': p_ab, 'N_Adj_AgeBMI': int(geno[['carrier', 'canon__age', 'canon__bmi']].dropna().shape[0]) if 'canon__age' in geno.columns and 'canon__bmi' in geno.columns else np.nan, 'Note_AgeBMI': 'Adjusted for age + BMI' if 'canon__age' in geno.columns and 'canon__bmi' in geno.columns else 'Age/BMI unavailable', 'OR_Het_vs_WT': or_het, 'OR_Het_vs_WT_CI95': _format_ci95(lo_het, hi_het), 'P_Het_vs_WT': p_het, 'N_Het': int((geno['Genotype_State'] == 'Het').sum()), 'OR_Hom_vs_WT': or_hom, 'OR_Hom_vs_WT_CI95': _format_ci95(lo_hom, hi_hom), 'P_Hom_vs_WT': p_hom, 'N_Hom': int((geno['Genotype_State'] == 'Hom').sum()), 'FDR_Unadj': np.nan, 'FDR_Adj_Age': np.nan, 'FDR_Adj_AgeBMI': np.nan, 'FDR_Het_vs_WT': np.nan, 'FDR_Hom_vs_WT': np.nan, 'Nominal_Sig_Unadj': p_unadj < 0.05, 'FDR_Sig_Unadj': False, 'Nominal_Sig_Adj_Age': p_age < 0.05 if pd.notna(p_age) else False, 'FDR_Sig_Adj_Age': False, 'Coverage_Risk_Flag': False, 'Coverage_Risk_Amplicon': '', 'Coverage_Risk_Region': '', 'Coverage_Risk_Note': '', 'P_Trend': p_tr, 'Trend_OR': or_tr, 'Trend_OR_CI95': _format_ci95(lo_tr, hi_tr)})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = _apply_fdr(out, 'P_Unadj', 'FDR_Unadj')
    out = _apply_fdr(out, 'P_Adj_Age', 'FDR_Adj_Age')
    out = _apply_fdr(out, 'P_Adj_AgeBMI', 'FDR_Adj_AgeBMI')
    out = _apply_fdr(out, 'P_Het_vs_WT', 'FDR_Het_vs_WT')
    out = _apply_fdr(out, 'P_Hom_vs_WT', 'FDR_Hom_vs_WT')
    out['Nominal_Sig_Unadj'] = out['P_Unadj'] < 0.05
    out['FDR_Sig_Unadj'] = out['FDR_Unadj'] < FDR_THRESHOLD
    out['Nominal_Sig_Adj_Age'] = out['P_Adj_Age'] < 0.05 if 'P_Adj_Age' in out.columns else False
    out['FDR_Sig_Adj_Age'] = out['FDR_Adj_Age'] < FDR_THRESHOLD if 'FDR_Adj_Age' in out.columns else False
    return annotate_cancer_risk_scheme(out)


def _fit_cox_survival(df: pd.DataFrame, duration_col: str, event_col: str, covariates: List[str]):
    needed = [duration_col, event_col] + covariates
    data = df[needed].dropna().copy()
    if data.empty:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    data[duration_col] = pd.to_numeric(data[duration_col], errors="coerce")
    data[event_col] = pd.to_numeric(data[event_col], errors="coerce")
    data = data.dropna(subset=[duration_col, event_col]).copy()
    if data.empty:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    data[event_col] = data[event_col].astype(int)
    if data[event_col].sum() < 3 or data[duration_col].nunique() < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    for col in covariates:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    data = data.dropna(subset=covariates).copy()
    if data.empty or any(data[col].nunique() < 2 for col in covariates):
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    try:
        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(data, duration_col=duration_col, event_col=event_col, show_progress=False)
        key = covariates[0]
        coef = float(cph.params_.get(key, np.nan))
        ci = cph.confidence_intervals_.loc[key] if key in cph.confidence_intervals_.index else pd.Series([np.nan, np.nan])
        p = float(cph.summary.loc[key, "p"]) if key in cph.summary.index else np.nan
        hr = float(np.exp(np.clip(coef, -20, 20))) if pd.notna(coef) else np.nan
        lo = float(np.exp(np.clip(ci.iloc[0], -20, 20))) if pd.notna(ci.iloc[0]) else np.nan
        hi = float(np.exp(np.clip(ci.iloc[1], -20, 20))) if pd.notna(ci.iloc[1]) else np.nan
        return hr, lo, hi, p, int(len(data)), int(data[event_col].sum()), int(data[key].nunique())
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


def _save_placeholder_plot(out_path: Path, title: str, subtitle: str = ''):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.55, title, ha='center', va='center', fontsize=14, fontweight='bold')
    if subtitle:
        ax.text(0.5, 0.4, subtitle, ha='center', va='center', fontsize=10, color='#555555')
    ax.set_axis_off()
    _save_figure(fig, out_path, pad_inches=0.18)


def make_volcano_plots(tvh: pd.DataFrame, out_dir: Path):
    if tvh.empty:
        return
    for suffix, p_col, threshold, ylabel in [
        ("raw_p", "P_Value", 0.05, "-log10 raw p-value"),
        ("FDR", "FDR_P_Value", FDR_THRESHOLD, "-log10 FDR q-value"),
    ]:
        fig, axes = plt.subplots(1, 3, figsize=(17.8, 6.2), squeeze=False, gridspec_kw={"wspace": 0.18})
        legend_handles = []
        for ax, grp in zip(axes[0], ["Breast", "Endometrium", "Global"]):
            sub = tvh[tvh['Analysis_Group'] == grp].copy()
            if sub.empty:
                ax.set_visible(False)
                continue
            sub = sub[pd.to_numeric(sub[p_col], errors='coerce').notna()].copy()
            sub['x'] = np.log2(pd.to_numeric(sub['Odds_Ratio'], errors='coerce').replace(0, np.nan))
            sub['y'] = _neglog10(sub[p_col])
            sub['Display_Label'] = sub.apply(_variant_display_label, axis=1)
            sub["P_Value"] = pd.to_numeric(sub.get("P_Value"), errors="coerce")
            sub["FDR_P_Value"] = pd.to_numeric(sub.get("FDR_P_Value"), errors="coerce")
            sub["Sig_Class"] = "Not significant"
            sub.loc[sub["P_Value"] < 0.05, "Sig_Class"] = "Nominal only"
            sub.loc[sub["FDR_P_Value"] < FDR_THRESHOLD, "Sig_Class"] = f"FDR < {FDR_THRESHOLD:g}"
            sub["Sig_Class"] = pd.Categorical(sub["Sig_Class"], categories=_VOLCANO_CLASS_ORDER, ordered=True)
            family = _cohort_family(grp)
            class_colors = {
                "Not significant": "#C9CED6",
                "Nominal only": family["light"],
                f"FDR < {FDR_THRESHOLD:g}": family["dark"],
            }
            for sig_class in _VOLCANO_CLASS_ORDER:
                layer = sub[sub["Sig_Class"] == sig_class]
                if layer.empty:
                    continue
                sns.scatterplot(
                    data=layer,
                    x="x",
                    y="y",
                    s=52 if sig_class != "Not significant" else 34,
                    color=class_colors[sig_class],
                    edgecolor="white",
                    linewidth=0.45,
                    alpha=0.9 if sig_class != "Not significant" else 0.65,
                    ax=ax,
                    legend=False,
                )
            ax.axhline(_neglog10(pd.Series([threshold])).iloc[0], ls='--', lw=1.05, c=family["dark"], alpha=0.85)
            ax.axvline(0, ls=':', lw=1, c='#7B8290', alpha=0.85)
            label_pool = sub[sub["Sig_Class"] != "Not significant"].copy()
            if label_pool.empty:
                label_pool = sub.nsmallest(min(6, len(sub)), p_col).copy()
            label_pool = label_pool.sort_values(["FDR_P_Value", "P_Value", "y"], na_position="last").head(8)
            left_pts = [(r["x"], r["y"], r["Display_Label"]) for _, r in label_pool[label_pool["x"] < 0].iterrows()]
            right_pts = [(r["x"], r["y"], r["Display_Label"]) for _, r in label_pool[label_pool["x"] >= 0].iterrows()]
            _add_margin_labels(ax, left_pts, side="left", max_labels=4, fontsize=7.2)
            _add_margin_labels(ax, right_pts, side="right", max_labels=4, fontsize=7.2)
            x_abs = np.nanmax(np.abs(sub["x"])) if sub["x"].notna().any() else 1
            ax.set_xlim(-max(1.1, x_abs * 1.18), max(1.1, x_abs * 1.18))
            ax.set_ylim(0, max(1.6, float(sub["y"].max()) * 1.12))
            ax.text(
                0.02,
                0.98,
                f"n={len(sub)} tests\nnominal={(sub['P_Value'] < 0.05).sum()} | FDR={(sub['FDR_P_Value'] < FDR_THRESHOLD).sum()}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#D0D5DD", alpha=0.96),
            )
            ax.set_title(f'{grp} tumour vs pooled control', fontsize=10.5, fontweight='bold', loc="left", pad=10)
            ax.set_xlabel('log2 odds ratio (tumour vs control)')
            ax.set_ylabel(ylabel)
            _style_ax(ax, grid=False)
            if not legend_handles:
                legend_handles = [
                    mpatches.Patch(facecolor=class_colors["Not significant"], label="Not significant"),
                    mpatches.Patch(facecolor=class_colors["Nominal only"], label="Nominal p < 0.05"),
                    mpatches.Patch(facecolor=class_colors[f"FDR < {FDR_THRESHOLD:g}"], label=f"FDR q < {FDR_THRESHOLD:g}"),
                ]
        fig.legend(
            handles=legend_handles,
            loc="upper center",
            ncol=3,
            frameon=False,
            bbox_to_anchor=(0.5, 1.01),
            fontsize=9,
            title="Association support",
            title_fontsize=9,
        )
        fig.suptitle(f'Tumour vs control SNP association volcano plots ({suffix.replace("_", " ")})', fontsize=13, fontweight='bold', y=1.05)
        _safe_tight_layout(fig, rect=[0.01, 0.02, 0.99, 0.90], pad=1.0, w_pad=2.4)
        out = out_dir / f'17_SNP_Volcano_{suffix}.png'
        _save_figure(fig, out, pad_inches=0.22)


def make_heatmaps(breast_clin: pd.DataFrame, endo_clin: pd.DataFrame, out_dir: Path, significant_only: bool = False):
    suffix = "significant only" if significant_only else "all associations"
    _plot_pivot_heatmap(
        breast_clin,
        "P_Unadj",
        out_dir / "17_SNP_Heatmap_Breast_raw_p.png",
        f"Breast clinical association heatmap ({suffix})",
        cbar_label="-log10 raw p",
        top_n=12,
        significant_only=significant_only,
        annotate=False,
        metric_name="raw p",
    )
    _plot_pivot_heatmap(
        breast_clin,
        "FDR_Unadj",
        out_dir / "17_SNP_Heatmap_Breast_FDR.png",
        f"Breast clinical association heatmap ({suffix})",
        cbar_label="-log10 FDR",
        top_n=12,
        significant_only=significant_only,
        annotate=False,
        metric_name="FDR",
    )
    _plot_pivot_heatmap(
        endo_clin,
        "P_Unadj",
        out_dir / "17_SNP_Heatmap_Endometrial_raw_p.png",
        f"Endometrial clinical association heatmap ({suffix})",
        cbar_label="-log10 raw p",
        top_n=12,
        significant_only=significant_only,
        annotate=False,
        metric_name="raw p",
    )
    _plot_pivot_heatmap(
        endo_clin,
        "FDR_Unadj",
        out_dir / "17_SNP_Heatmap_Endometrial_FDR.png",
        f"Endometrial clinical association heatmap ({suffix})",
        cbar_label="-log10 FDR",
        top_n=12,
        significant_only=significant_only,
        annotate=False,
        metric_name="FDR",
    )


def make_main_text_heatmaps(breast_clin: pd.DataFrame, endo_clin: pd.DataFrame, out_dir: Path):
    _plot_pivot_heatmap(
        breast_clin,
        "FDR_Unadj",
        out_dir / "17_SNP_Heatmap_Breast_MainText.png",
        "Breast main-text clinical heatmap",
        cbar_label="-log10 FDR",
        top_n=8,
        significant_only=True,
        annotate=True,
        metric_name="FDR",
    )
    _plot_pivot_heatmap(
        endo_clin,
        "FDR_Unadj",
        out_dir / "17_SNP_Heatmap_Endometrial_MainText.png",
        "Endometrial main-text clinical heatmap",
        cbar_label="-log10 FDR",
        top_n=8,
        significant_only=True,
        annotate=True,
        metric_name="FDR",
    )


def make_forest_plots(breast_clin: pd.DataFrame, endo_clin: pd.DataFrame, out_dir: Path):
    for cohort, df in [('Breast', breast_clin), ('Endometrial', endo_clin)]:
        if df.empty:
            continue
        sub = df[pd.to_numeric(df.get('P_Unadj', np.nan), errors='coerce').notna()].copy()
        if sub.empty:
            _save_placeholder_plot(out_dir / f'17_Forest_{cohort}.png', f'{cohort} clinical forest plot')
            continue
        sub = sub.nsmallest(min(18, len(sub)), 'P_Unadj')
        fig, ax = plt.subplots(figsize=(10, max(5, 0.45 * len(sub) + 2)))
        y = np.arange(len(sub))[::-1]
        ors = pd.to_numeric(sub.get('OR_Adj_Age', sub.get('OR_Unadj', 1)), errors='coerce').fillna(1.0).to_numpy()
        if 'OR_Adj_Age_CI95' in sub.columns:
            ci_bounds = sub['OR_Adj_Age_CI95'].astype(str).str.extract(r'^\[\s*([^\],]+)\s*,\s*([^\]]+)\s*\]$')
            ci_lo = pd.to_numeric(ci_bounds[0], errors='coerce').to_numpy()
            ci_hi = pd.to_numeric(ci_bounds[1], errors='coerce').to_numpy()
            valid_ci = np.isfinite(ci_lo) & np.isfinite(ci_hi)
            ax.errorbar(
                ors[valid_ci],
                y[valid_ci],
                xerr=[np.clip(ors[valid_ci] - ci_lo[valid_ci], 0, None), np.clip(ci_hi[valid_ci] - ors[valid_ci], 0, None)],
                fmt='none',
                ecolor='#C9CED6',
                elinewidth=1.2,
                capsize=3,
                zorder=1,
            )
        ax.scatter(ors, y, s=55, color='#D55E00', edgecolor='white', linewidth=0.5, zorder=2)
        ax.axvline(1, ls='--', c='#666666')
        ax.set_yticks(y)
        ax.set_yticklabels([_wrap_label(f"{r['Clin_Label']} · {_preferred_variant_id(r)}", 38) for _, r in sub.iterrows()], fontsize=8)
        ax.set_xlabel('Odds ratio (95% CI)')
        ax.set_title(f'{cohort} clinical associations', fontweight='bold')
        _style_ax(ax, grid=False)
        out = out_dir / f'17_Forest_{cohort}.png'
        _safe_tight_layout(fig, pad=1.0)
        _save_figure(fig, out)


def make_distribution_plots(breast_clin: pd.DataFrame, endo_clin: pd.DataFrame, merged: pd.DataFrame, out_dir: Path):
    for cohort, df in [("Breast", breast_clin), ("Endometrial", endo_clin)]:
        out = out_dir / f"17_Distribution_{cohort}.png"
        if df.empty or "P_Unadj" not in df.columns:
            _save_placeholder_plot(out, f"{cohort} clinical p-value distribution", "No results available")
            continue
        vals = pd.to_numeric(df["P_Unadj"], errors="coerce").dropna()
        if vals.empty:
            _save_placeholder_plot(out, f"{cohort} clinical p-value distribution", "No results available")
            continue
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        sns.histplot(-np.log10(vals.clip(lower=1e-300)), bins=20, color="#4C72B0", edgecolor="white", ax=ax)
        ax.axvline(-np.log10(0.05), ls="--", lw=1.2, c="#D55E00", label="Nominal p=0.05")
        ax.axvline(-np.log10(FDR_THRESHOLD), ls=":", lw=1.2, c="#009E73", label=f"FDR={FDR_THRESHOLD:g}")
        ax.set_xlabel("-log10 raw p")
        ax.set_ylabel("Number of associations")
        ax.set_title(f"{cohort} clinical association p-value distribution", fontweight="bold")
        sig_nom = int((pd.to_numeric(df["P_Unadj"], errors="coerce") < 0.05).sum())
        sig_fdr = int((pd.to_numeric(df["FDR_Unadj"], errors="coerce") < FDR_THRESHOLD).sum()) if "FDR_Unadj" in df.columns else 0
        ax.text(
            0.98,
            0.98,
            f"Nominal: {sig_nom}\nFDR: {sig_fdr}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#cccccc", alpha=0.95),
        )
        ax.legend(frameon=False, loc="upper left")
        _style_ax(ax)
        _safe_tight_layout(fig, pad=1.0)
        _save_figure(fig, out)


def make_binary_bar_plots(breast_clin: pd.DataFrame, endo_clin: pd.DataFrame, merged: pd.DataFrame, out_dir: Path):
    for cohort, df in [("Breast", breast_clin), ("Endometrial", endo_clin)]:
        out = out_dir / f"17_BinaryBar_{cohort}.png"
        if df.empty or "Clin_Label" not in df.columns:
            _save_placeholder_plot(out, f"{cohort} clinical significance counts", "No results available")
            continue
        rows = []
        for label, sub in df.groupby("Clin_Label"):
            rows.append(
                {
                    "Clin_Label": label,
                    "Nominal": int(pd.Series(sub.get("Nominal_Sig_Unadj", pd.Series(dtype=bool))).fillna(False).sum()) if "Nominal_Sig_Unadj" in sub.columns else 0,
                    "FDR": int(pd.Series(sub.get("FDR_Sig_Unadj", pd.Series(dtype=bool))).fillna(False).sum()) if "FDR_Sig_Unadj" in sub.columns else 0,
                }
            )
        summary = pd.DataFrame(rows).set_index("Clin_Label")
        if summary.empty:
            _save_placeholder_plot(out, f"{cohort} clinical significance counts", "No results available")
            continue
        summary = summary.sort_values(["FDR", "Nominal"], ascending=True)
        labels = summary.index.tolist()
        y = np.arange(len(labels))
        fig, ax = plt.subplots(figsize=(9.2, max(4.5, 0.45 * len(labels) + 1.5)))
        ax.barh(y - 0.18, summary["Nominal"], height=0.34, color="#4C72B0", label="Nominal")
        ax.barh(y + 0.18, summary["FDR"], height=0.34, color="#DD8452", label=f"FDR<{FDR_THRESHOLD:g}")
        ax.set_yticks(y)
        ax.set_yticklabels([_wrap_label(label, 28) for label in labels], fontsize=9)
        ax.set_xlabel("Number of significant associations")
        ax.set_title(f"{cohort} clinical association counts by endpoint", fontweight="bold")
        ax.legend(frameon=False, loc="lower right")
        _style_ax(ax)
        _safe_tight_layout(fig, pad=1.0)
        _save_figure(fig, out)


def make_bmi_category_plot(bmi_summary: pd.DataFrame, bmi_availability: pd.DataFrame, out_dir: Path):
    out = out_dir / "17_BMI_Category_Composition_ByCohort.png"
    if bmi_summary.empty or bmi_availability.empty:
        _save_placeholder_plot(out, "BMI category composition by cohort", "No usable BMI values available")
        return

    pivot_pct = bmi_summary.pivot_table(
        index="Cohort_Display",
        columns="Category",
        values="Percentage",
        aggfunc="first",
    ).reindex(columns=BMI_CATEGORY_ORDER).fillna(0)
    pivot_count = bmi_summary.pivot_table(
        index="Cohort_Display",
        columns="Category",
        values="Count",
        aggfunc="first",
    ).reindex(columns=BMI_CATEGORY_ORDER).fillna(0)
    cohort_order = bmi_availability["Cohort_Display"].tolist()
    pivot_pct = pivot_pct.reindex(cohort_order)
    pivot_count = pivot_count.reindex(cohort_order)
    availability = bmi_availability.set_index("Cohort_Display").reindex(cohort_order)

    fig, ax = plt.subplots(figsize=(11.8, 6.8))
    x = np.arange(len(pivot_pct.index))
    bottoms = np.zeros(len(pivot_pct.index))
    for category in BMI_CATEGORY_ORDER:
        vals = pivot_pct[category].to_numpy(dtype=float)
        counts = pivot_count[category].to_numpy(dtype=float)
        bars = ax.bar(
            x,
            vals,
            bottom=bottoms,
            width=0.62,
            color=BMI_CATEGORY_COLORS[category],
            edgecolor="white",
            linewidth=1.0,
            label=category,
        )
        for bar, val, count, bottom in zip(bars, vals, counts, bottoms):
            if val >= 9:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bottom + val / 2,
                    f"{int(count)}\n{val:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=8.3,
                    color="white" if category == "Obese" else "#1F2933",
                    fontweight="bold",
                )
        bottoms += vals

    xlabels = []
    for cohort_display in pivot_pct.index:
        row = availability.loc[cohort_display]
        usable = int(row.get("Samples_With_Usable_BMI_Category", 0))
        xlabels.append(f"{cohort_display}\nBMI n={usable}")
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.set_ylabel("BMI category composition within cohort")
    ax.set_title("BMI category composition by cohort", loc="left", fontsize=12, fontweight="bold", pad=12)
    usable_bits = [f"{row['Cohort_Display']} n={int(row['Samples_With_Usable_BMI_Category'])}" for _, row in bmi_availability.iterrows()]
    excluded_under20 = int(bmi_availability["Samples_Below_20_Excluded"].sum())
    excluded_missing = int(bmi_availability["Samples_Missing_or_NonNumeric_BMI"].sum())
    ax.text(
        0.0,
        -0.22,
        "BMI source: canon__bmi. Percentages use only samples with usable BMI categories.\n"
        + "Usable BMI counts: "
        + " | ".join(usable_bits)
        + f" | Excluded overall: BMI<20 (n={excluded_under20}), missing/non-numeric BMI (n={excluded_missing})",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color="#5B6472",
    )
    legend = ax.legend(
        title="BMI category",
        loc="upper center",
        bbox_to_anchor=(0.5, 1.03),
        ncol=3,
        frameon=False,
        fontsize=9,
        title_fontsize=9,
    )
    _style_ax(ax)
    ax.xaxis.grid(False)
    _safe_tight_layout(fig, rect=[0, 0.09, 1, 0.93], pad=1.0)
    _save_figure(fig, out)


def make_km_plots(surv_res: pd.DataFrame, km_pages, out_dir: Path):
    if surv_res.empty or not km_pages or not _HAS_LIFELINES:
        return
    pages_by_model = {}
    for page in km_pages:
        pages_by_model.setdefault(page.get("Model_Slug", "genotype_only"), []).append(page)

    model_titles = {
        "genotype_only": "Genotype only",
        "genotype_age": "Genotype, age adjusted",
        "genotype_bmi": "Genotype, BMI adjusted",
        "genotype_age_bmi": "Genotype, age + BMI adjusted",
        "carrier_only": "Carrier vs non-carrier",
        "carrier_age": "Carrier, age adjusted",
        "carrier_bmi": "Carrier, BMI adjusted",
        "carrier_age_bmi": "Carrier, age + BMI adjusted",
    }
    plot_modes = {
        "genotype_only": "genotype",
        "genotype_age": "genotype",
        "genotype_bmi": "genotype",
        "genotype_age_bmi": "genotype",
        "carrier_only": "carrier",
        "carrier_age": "carrier",
        "carrier_bmi": "carrier",
        "carrier_age_bmi": "carrier",
    }
    carrier_colors = {"WT": "#BDBDBD", "Carrier": "#0072B2"}

    for model_slug, pages in pages_by_model.items():
        out_name = "17_KM_Curves_Survival.pdf" if model_slug == "genotype_only" else f"17_KM_Curves_Survival_{model_slug}.pdf"
        out_path = out_dir / out_name
        with pdf_backend.PdfPages(out_path) as pdf:
            for page in pages:
                frame = page["Frame"].copy()
                if frame.empty:
                    continue
                fig, ax = plt.subplots(figsize=(9.6, 6.8))
                plotted = False
                max_time = 0.0
                plot_mode = plot_modes.get(model_slug, "genotype")
                if plot_mode == "carrier":
                    state_col = "Carrier_State"
                    states = ["WT", "Carrier"]
                    colors = carrier_colors
                else:
                    state_col = "Genotype_State"
                    states = ["WT", "Het", "Hom"]
                    colors = _GENO_PLOT_C
                if state_col not in frame.columns or frame[state_col].nunique() < 2:
                    plt.close(fig)
                    continue
                for state in states:
                    sub = frame[frame[state_col] == state].copy()
                    if sub.empty:
                        continue
                    sub = sub.dropna(subset=["Duration", "Event"])
                    if sub.empty:
                        continue
                    plotted = True
                    max_time = max(max_time, float(sub["Duration"].max()))
                    label = f"{state} (n={len(sub)})"
                    kmf = KaplanMeierFitter()
                    kmf.fit(sub["Duration"], event_observed=sub["Event"], label=label)
                    kmf.plot_survival_function(
                        ax=ax,
                        ci_show=False,
                        color=colors.get(state, "#4C72B0"),
                        linewidth=2.4,
                    )
                if not plotted:
                    plt.close(fig)
                    continue
                try:
                    lr = multivariate_logrank_test(frame["Duration"], frame[state_col], frame["Event"])
                    logrank_p = float(lr.p_value)
                except Exception:
                    logrank_p = page.get("P_LogRank", np.nan)
                logrank_fdr = pd.to_numeric(page.get("FDR_LogRank", np.nan), errors="coerce")
                cox_p = pd.to_numeric(page.get("P_Cox_Additive", np.nan), errors="coerce")
                cox_fdr = pd.to_numeric(page.get("FDR_Cox_Additive", np.nan), errors="coerce")
                ax.set_xlim(0, max_time * 1.05 if max_time > 0 else None)
                ax.set_ylim(0, 1.05)
                ax.set_xlabel("Time (months)")
                ax.set_ylabel("Survival probability")
                variant_text = _variant_display_label(pd.Series(page), max_len=80, include_gene=True)
                model_label = model_titles.get(model_slug, page.get("Model", model_slug))
                title = f"{page['Cohort']} {page['Endpoint']} survival by {variant_text} ({model_label})".strip()
                ax.set_title(_wrap_label(title, 54, max_lines=2), fontsize=11.5, fontweight="bold", loc="left", pad=12)
                cox_bits = [f"Cox additive {_format_probability(cox_p, 'p')} {_sig_label(cox_p)}".strip()]
                if pd.notna(cox_fdr):
                    cox_bits.append(f"{_format_probability(cox_fdr, 'q')} ({'FDR sig' if cox_fdr < FDR_THRESHOLD else 'ns'})")
                lr_bits = [f"Log-rank {_format_probability(logrank_p, 'p')} {_sig_label(logrank_p)}".strip()]
                if pd.notna(logrank_fdr):
                    lr_bits.append(f"{_format_probability(logrank_fdr, 'q')} ({'FDR sig' if logrank_fdr < FDR_THRESHOLD else 'ns'})")
                state_counts = frame[state_col].value_counts()
                count_bits = " | ".join([f"{state} n={int(state_counts.get(state, 0))}" for state in states if int(state_counts.get(state, 0)) > 0])
                annot = "\n".join([" | ".join(cox_bits), " | ".join(lr_bits), count_bits])
                ax.text(
                    0.98,
                    0.03,
                    annot,
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=8.6,
                    bbox=dict(boxstyle="round,pad=0.34", fc="white", ec="#D0D5DD", alpha=0.97),
                )
                _style_ax(ax, grid=False)
                ax.yaxis.grid(True, linestyle=":", color="#D5DBE3", alpha=0.9, linewidth=0.7)
                ax.xaxis.grid(True, linestyle=":", color="#E4E7EC", alpha=0.75, linewidth=0.6)
                legend = ax.legend(title="Group", frameon=True, loc="upper right", fontsize=8.5, title_fontsize=8.5)
                legend.get_frame().set_edgecolor("#D0D5DD")
                legend.get_frame().set_facecolor("white")
                legend.get_frame().set_alpha(0.95)
                _safe_tight_layout(fig, rect=[0.0, 0.02, 1.0, 0.96], pad=1.0)
                pdf.savefig(fig, bbox_inches="tight", pad_inches=0.18)
                plt.close(fig)
        print(f"  Saved: {out_path}")


def make_summary_panel(tvh: pd.DataFrame, breast_clin: pd.DataFrame, endo_clin: pd.DataFrame, surv_res: pd.DataFrame, out_dir: Path):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    panels = [(tvh, 'Tumour vs control', 'Nominal_Sig', 'FDR_Sig'), (breast_clin, 'Breast clinical', 'Nominal_Sig_Unadj', 'FDR_Sig_Unadj'), (endo_clin, 'Endometrial clinical', 'Nominal_Sig_Unadj', 'FDR_Sig_Unadj'), (surv_res, 'Survival', 'P_Cox_Additive', 'FDR_Cox_Additive')]
    for ax, (df, title, nominal_col, fdr_col) in zip(axes.flat, panels):
        if df is None or df.empty:
            ax.axis('off')
            continue
        if nominal_col in df.columns and pd.api.types.is_bool_dtype(df[nominal_col]):
            n_nom = int(df[nominal_col].fillna(False).sum())
        else:
            n_nom = int((pd.to_numeric(df[nominal_col], errors='coerce') < 0.05).sum()) if nominal_col in df.columns else 0
        if fdr_col in df.columns and pd.api.types.is_bool_dtype(df[fdr_col]):
            n_fdr = int(df[fdr_col].fillna(False).sum())
        else:
            n_fdr = int((pd.to_numeric(df[fdr_col], errors='coerce') < FDR_THRESHOLD).sum()) if fdr_col in df.columns else 0
        ax.bar(['Nominal', 'FDR'], [n_nom, n_fdr], color=['#4C72B0', '#DD8452'])
        ax.set_title(title, fontweight='bold')
        _style_ax(ax)
    fig.suptitle('Script 17 Summary Panel', fontsize=13, fontweight='bold')
    _safe_tight_layout(fig, rect=[0, 0, 1, 0.96], pad=1.0)
    out = out_dir / 'GSDMB_SNP_Association_Summary_Panel.png'
    _save_figure(fig, out)


def make_genotype_distribution_plots(breast_dose: pd.DataFrame, endo_dose: pd.DataFrame, merged: pd.DataFrame, out_dir: Path):
    for cohort, df in [("Breast", breast_dose), ("Endometrial", endo_dose)]:
        out = out_dir / f"17_GenoDistribution_{cohort}.png"
        if df.empty:
            _save_placeholder_plot(out, f"{cohort} genotype contrast heatmap", "No results available")
            continue
        work = df.copy()
        contrast_cols = [c for c in ["P_Trend", "P_Het_vs_WT", "P_Hom_vs_WT"] if c in work.columns]
        if not contrast_cols:
            _save_placeholder_plot(out, f"{cohort} genotype contrast heatmap", "No results available")
            continue
        for c in contrast_cols:
            work[c] = pd.to_numeric(work[c], errors="coerce")
        work = work.dropna(subset=[contrast_cols[0]]).copy()
        if work.empty:
            _save_placeholder_plot(out, f"{cohort} genotype contrast heatmap", "No results available")
            continue
        work["Variant_Label"] = work.apply(_variant_label, axis=1)
        order = (
            work.groupby("Variant_Label")[contrast_cols[0]]
            .min()
            .sort_values(ascending=True)
            .head(8)
            .index
            .tolist()
        )
        work = work[work["Variant_Label"].isin(order)].copy()
        work["Variant_Label"] = pd.Categorical(work["Variant_Label"], categories=order, ordered=True)
        melted = []
        for col, label in [("P_Trend", "Trend"), ("P_Het_vs_WT", "Het vs WT"), ("P_Hom_vs_WT", "Hom vs WT")]:
            if col in work.columns:
                tmp = work[["Variant_Label", col]].copy()
                tmp["Contrast"] = label
                tmp["Score"] = _neglog10(tmp[col])
                melted.append(tmp[["Variant_Label", "Contrast", "Score"]])
        long_df = pd.concat(melted, ignore_index=True)
        pivot = long_df.pivot_table(index="Variant_Label", columns="Contrast", values="Score", aggfunc="max")
        pivot = pivot.reindex(columns=[c for c in ["Trend", "Het vs WT", "Hom vs WT"] if c in pivot.columns])
        fig, ax = plt.subplots(figsize=(8.8, max(4.3, 0.42 * len(pivot.index) + 1.6)))
        sns.heatmap(
            pivot,
            ax=ax,
            cmap=sns.color_palette("mako", as_cmap=True),
            linewidths=0.5,
            linecolor="#E6E6E6",
            annot=True,
            fmt=".1f",
            cbar_kws={"label": "-log10 p"},
            annot_kws={"size": 7},
        )
        ax.set_title(f"{cohort} genotype model contrasts", fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=20, labelsize=8)
        ax.tick_params(axis="y", labelsize=8)
        _style_ax(ax, grid=False)
        _safe_tight_layout(fig, pad=1.0)
        _save_figure(fig, out)


def make_genotype_binary_bar_plots(breast_dose: pd.DataFrame, endo_dose: pd.DataFrame, merged: pd.DataFrame, out_dir: Path):
    for cohort, df in [("Breast", breast_dose), ("Endometrial", endo_dose)]:
        out = out_dir / f"17_GenoBinaryBar_{cohort}.png"
        if df.empty:
            _save_placeholder_plot(out, f"{cohort} genotype significance counts", "No results available")
            continue
        counts = []
        for contrast, p_col, fdr_col in [
            ("Trend", "P_Trend", "FDR_Trend"),
            ("Het vs WT", "P_Het_vs_WT", "FDR_Het_vs_WT"),
            ("Hom vs WT", "P_Hom_vs_WT", "FDR_Hom_vs_WT"),
        ]:
            if p_col not in df.columns:
                continue
            pvals = pd.to_numeric(df[p_col], errors="coerce")
            fdrs = pd.to_numeric(df[fdr_col], errors="coerce") if fdr_col in df.columns else pd.Series(np.nan, index=df.index)
            counts.append(
                {
                    "Contrast": contrast,
                    "Nominal": int((pvals < 0.05).sum()),
                    "FDR": int((fdrs < FDR_THRESHOLD).sum()) if fdr_col in df.columns else 0,
                }
            )
        summary = pd.DataFrame(counts)
        if summary.empty:
            _save_placeholder_plot(out, f"{cohort} genotype significance counts", "No results available")
            continue
        fig, ax = plt.subplots(figsize=(8.2, 4.8))
        y = np.arange(len(summary))
        ax.barh(y - 0.18, summary["Nominal"], height=0.34, color="#4C72B0", label="Nominal")
        ax.barh(y + 0.18, summary["FDR"], height=0.34, color="#DD8452", label=f"FDR<{FDR_THRESHOLD:g}")
        ax.set_yticks(y)
        ax.set_yticklabels(summary["Contrast"], fontsize=9)
        ax.set_xlabel("Number of significant genotype-contrast tests")
        ax.set_title(f"{cohort} genotype model significance counts", fontweight="bold")
        ax.legend(frameon=False, loc="lower right")
        _style_ax(ax)
        _safe_tight_layout(fig, pad=1.0)
        _save_figure(fig, out)


def make_genotype_composition_plots(tvh: pd.DataFrame, merged: pd.DataFrame, out_dir: Path):
    sample_manifest = _sample_manifest(merged)
    legend_handles = [
        mpatches.Patch(facecolor=_GENO_PLOT_C[state], edgecolor='#444444', label=state)
        for state in ['WT', 'Het', 'Hom']
    ]
    for grp in ['Breast', 'Endometrium', 'Global']:
        sub = tvh[tvh['Analysis_Group'] == grp].copy()
        if sub.empty:
            continue
        top = sub.sort_values(['FDR_P_Value', 'P_Value'], na_position='last').head(6)
        if grp == 'Breast':
            target_samples = sample_manifest[sample_manifest['Cohort'].astype(str).str.contains('Breast', case=False, na=False)]
        elif grp == 'Endometrium':
            target_samples = sample_manifest[sample_manifest['Cohort'].astype(str).str.contains('Endometri', case=False, na=False)]
        else:
            target_samples = sample_manifest.copy()
        if target_samples.empty:
            continue
        pooled_controls = target_samples[target_samples['Tissue'].astype(str) == 'Healthy'].copy()
        tumour_samples = target_samples[target_samples['Tissue'].astype(str) == 'Tumour'].copy()
        fig, axes = plt.subplots(2, 3, figsize=(15, 8), squeeze=False)
        for ax, (_, row) in zip(axes.flat, top.iterrows()):
            geno = _variant_sample_table(merged, row['Variant_ID'], target_samples['Sample'])
            if geno.empty:
                ax.set_visible(False)
                continue
            frame = pd.concat([pooled_controls[['snp_code']].assign(Group='Pooled normal'), tumour_samples[['snp_code']].assign(Group='Tumour')], ignore_index=True).merge(geno[['snp_code', 'Genotype_State']], on='snp_code', how='left')
            counts = frame.groupby(['Group', 'Genotype_State']).size().unstack(fill_value=0).reindex(index=['Pooled normal', 'Tumour'], columns=['WT', 'Het', 'Hom'], fill_value=0)
            if counts.empty:
                ax.set_visible(False)
                continue
            perc = counts.div(counts.sum(axis=1), axis=0).fillna(0) * 100
            bottom = np.zeros(len(perc.index))
            for geno_state in ['WT', 'Het', 'Hom']:
                vals = perc[geno_state].to_numpy(dtype=float)
                ax.bar(np.arange(len(vals)), vals, bottom=bottom, color=_GENO_PLOT_C[geno_state], edgecolor='#444444', width=0.6, label=geno_state)
                bottom += vals
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f'Pooled normal\n(n={int(counts.sum(axis=1).iloc[0])})', f'Tumour\n(n={int(counts.sum(axis=1).iloc[1])})'], fontsize=8)
            ax.set_ylim(0, 100)
            ax.set_title(_wrap_label(_variant_display_label(row, max_len=24, include_gene=True), 24), fontsize=9)
            _style_ax(ax)
        fig.legend(handles=legend_handles, loc='upper center', ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02), title='Genotype state')
        fig.suptitle(f'{grp} genotype composition in tumour and pooled control samples', fontsize=12, fontweight='bold')
        _safe_tight_layout(fig, rect=[0, 0, 1, 0.95], pad=1.0)
        out = out_dir / f'17_GenotypeComposition_{grp}.png'
        _save_figure(fig, out)

def make_risk_forest_plot(risk_res: pd.DataFrame, out_dir: Path):
    """
    Forest plot of cancer risk ORs (age-adjusted logistic regression),
    one panel per cohort, sorted by p-value.
    Confidence intervals capped at OR_CAP for display.
    """
    if risk_res.empty:
        return

    OR_CAP = 15.0

    cohorts = risk_res["Cohort"].unique()
    fig, axes = plt.subplots(1, len(cohorts),
                             figsize=(9 * len(cohorts), max(5, len(risk_res) // len(cohorts) * 0.45 + 2)),
                             squeeze=False)

    for ax, cohort in zip(axes[0], cohorts):
        sub = risk_res[risk_res["Cohort"] == cohort].copy()
        sub = sub.dropna(subset=["OR_Adj_Age"]).sort_values("P_Adj_Age")
        if sub.empty:
            ax.set_visible(False)
            continue

        def _parse_ci(s):
            try:
                lo, hi = re.findall(r"[-0-9.eE+]+", str(s))[:2]
                return float(lo), float(hi)
            except Exception:
                return np.nan, np.nan

        sub[["CI_lo", "CI_hi"]] = sub["OR_Adj_Age_CI95"].apply(lambda s: pd.Series(_parse_ci(s)))
        sub["OR_plot"]    = sub["OR_Adj_Age"].clip(upper=OR_CAP)
        sub["CI_hi_plot"] = sub["CI_hi"].clip(upper=OR_CAP)
        sub["CI_lo_plot"] = sub["CI_lo"].clip(lower=0)
        sub["capped"]     = sub["CI_hi"] > OR_CAP

        col   = _COHORT_C.get(cohort, "#555555")
        ypos  = np.arange(len(sub))[::-1]
        sig_mask = sub["P_Adj_Age"] < 0.05
        point_c  = [col if s else "#999999" for s in sig_mask]

        ax.errorbar(sub["OR_plot"].values, ypos,
                    xerr=[np.clip((sub["OR_plot"] - sub["CI_lo_plot"]).values, 0, None),
                          np.clip((sub["CI_hi_plot"] - sub["OR_plot"]).values, 0, None)],
                    fmt="none", ecolor="#cccccc", elinewidth=1.2, capsize=3, zorder=2)

        for i, (_, row) in enumerate(sub.iterrows()):
            ax.scatter(row["OR_plot"], ypos[i], s=90, color=point_c[i],
                       edgecolors="white", linewidths=0.6, zorder=3)
            if row["capped"]:
                ax.annotate("→", xy=(OR_CAP, ypos[i]), fontsize=9,
                            color="#aaaaaa", va="center")

        ax.set_xlim(0, OR_CAP * 1.10)
        for i, (_, row) in enumerate(sub.iterrows()):
            p_str = f"p={row['P_Adj_Age']:.3f}{_sig_label(row['P_Adj_Age'])}"
            n_str = f"n={int(row['N_Adj_Age'])}" if pd.notna(row.get("N_Adj_Age")) else ""
            label_txt = f"{p_str}  {n_str}".strip()
            ax.text(1.02, ypos[i], label_txt, transform=ax.get_yaxis_transform(),
                    va="center", ha="left", fontsize=6.8, color="#444444", clip_on=False)

        ax.axvline(1, color="#555555", linestyle="--", linewidth=1, alpha=0.7)
        labels = sub.apply(lambda r: _variant_display_label(r, max_len=60, include_gene=True), axis=1)
        ax.set_yticks(ypos)
        ax.set_yticklabels([_wrap_label(v, width=30, max_lines=2) for v in labels.values], fontsize=8)
        ax.set_xlabel("Odds Ratio for cancer risk\n(age-adjusted, 95% CI)", fontsize=10)
        ax.set_title(
            f"Age-Adjusted Association Between GSDMB SNP Any-ALT Carrier Status and {cohort} Cancer Risk\n"
            f"(cases n={risk_res.loc[risk_res['Cohort']==cohort, 'N_Cases'].iloc[0]}, "
            f"controls n={risk_res.loc[risk_res['Cohort']==cohort, 'N_Controls'].iloc[0]})",
            fontsize=11,
            fontweight="bold",
            pad=10,
        )

        xlims = ax.get_xlim()
        ax.axvspan(1, xlims[1], alpha=0.04, color="#c62828")
        ax.axvspan(xlims[0], 1, alpha=0.04, color="#1976D2")
        ax.text(0.02, 0.01, "← protective", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=7.5, color="#1976D2", style="italic")
        ax.text(0.98, 0.01, "risk →", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7.5, color="#c62828", style="italic")
        _style_ax(ax, grid=False)
        ax.xaxis.grid(True, linestyle=":", color="#cccccc", alpha=0.6)
        ax.set_axisbelow(True)

    fig.suptitle(
        "Age-Adjusted Associations Between GSDMB SNPs and Cancer Risk",
        fontsize=13,
        fontweight="bold",
        y=1.01,
    )
    _safe_tight_layout(fig, rect=[0, 0, 0.9, 1], pad=1.0)
    out_path = out_dir / "17_Forest_CancerRisk.png"
    _save_figure(fig, out_path)


def make_risk_genotype_composition_plots(risk_res, merged, out_dir):
    """Plot case/control genotype composition for the cancer-risk analysis."""
    if risk_res.empty or "GT" not in merged.columns:
        return

    TOP_N = 6
    GT_MAP = {"0/0": "WT", "0/1": "Het", "1/0": "Het", "1/1": "Hom"}
    GENO_ORDER = ["WT", "Het", "Hom"]

    sample_manifest = merged[~merged["is_replicate"]].drop_duplicates("Sample").copy()
    pooled_controls = sample_manifest[sample_manifest["Tissue"] == "Healthy"].copy()

    comparison_defs = {
        "Breast": {
            "tumour_mask": (sample_manifest["Tissue"] == "Tumour") & sample_manifest["Cohort"].astype(str).str.contains("Breast", case=False, na=False),
            "tumour_label": "Breast tumour",
            "control_label": "Pooled control",
        },
        "Endometrial": {
            "tumour_mask": (sample_manifest["Tissue"] == "Tumour") & sample_manifest["Cohort"].astype(str).str.contains("Endometri", case=False, na=False),
            "tumour_label": "Endometrium tumour",
            "control_label": "Pooled control",
        },
    }

    for cohort, cfg in comparison_defs.items():
        sub = (
            risk_res[risk_res["Cohort"] == cohort]
            .sort_values(["FDR_Sig_Adj_Age", "Nominal_Sig_Adj_Age", "P_Adj_Age", "P_Unadj"],
                         ascending=[False, False, True, True])
            .drop_duplicates("Variant_ID")
            .copy()
        )
        if sub.empty:
            continue

        chosen = sub[sub.get("Nominal_Sig_Adj_Age", False) | sub.get("FDR_Sig_Adj_Age", False)].head(TOP_N).copy()
        if chosen.empty:
            chosen = sub.head(TOP_N).copy()
        if chosen.empty:
            continue

        tumour_samples = sample_manifest.loc[cfg["tumour_mask"]].copy()
        control_samples = pooled_controls.copy()
        if tumour_samples.empty or control_samples.empty:
            continue

        arm_df = pd.concat([
            control_samples.assign(_arm=cfg["control_label"]),
            tumour_samples.assign(_arm=cfg["tumour_label"]),
        ], ignore_index=True)[["Sample", "_arm"]]

        n_plots = len(chosen)
        ncols = min(3, n_plots)
        nrows = int(np.ceil(n_plots / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.2 * nrows), squeeze=False, facecolor="white")
        axes_flat = axes.flatten()
        for ax in axes_flat:
            ax.set_facecolor("white")

        for idx, (_, row) in enumerate(chosen.iterrows()):
            ax = axes_flat[idx]
            var_id = row["Variant_ID"]
            gt_lookup = (
                merged[(merged["Variant_ID"] == var_id) & (~merged["is_replicate"]) ]
                .drop_duplicates("Sample")
                .set_index("Sample")["GT"]
                .apply(lambda g: GT_MAP.get(str(g), "WT"))
            )

            plot_df = arm_df.copy()
            plot_df["Genotype"] = plot_df["Sample"].map(gt_lookup).fillna("WT")
            counts = (
                plot_df.groupby(["_arm", "Genotype"]).size().unstack(fill_value=0)
                .reindex(index=[cfg["control_label"], cfg["tumour_label"]], columns=GENO_ORDER, fill_value=0)
            )
            totals = counts.sum(axis=1)
            if (totals == 0).any():
                ax.set_visible(False)
                continue

            perc = counts.div(totals, axis=0) * 100
            xpos = np.arange(len(counts.index))
            bottoms = np.zeros(len(counts.index))
            for geno in GENO_ORDER:
                vals = perc[geno].to_numpy(dtype=float)
                ax.bar(xpos, vals, bottom=bottoms, color=_GENO_PLOT_C[geno], edgecolor="#4A4A4A", linewidth=0.9,
                       hatch=_GENO_PLOT_HATCH[geno], width=0.58, label=geno)
                for xi, val, bottom in zip(xpos, vals, bottoms):
                    if val >= 9:
                        ax.text(xi, bottom + val / 2, f"{geno}\n{val:.0f}%", ha="center", va="center",
                                fontsize=8.2, color=_GENO_PLOT_TEXT[geno], fontweight="bold")
                bottoms += vals

            ax.set_ylim(0, 100)
            ax.set_xticks(xpos)
            ax.set_xticklabels([f"{label}\n(n={int(totals.loc[label])})" for label in counts.index], fontsize=8.5)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
            ax.set_ylabel("Genotype composition (%)", fontsize=9)
            snp_label = _variant_display_label(row, max_len=72, include_gene=True)
            p_val = row.get("P_Adj_Age", np.nan)
            sig = "FDR" if bool(row.get("FDR_Sig_Adj_Age", False)) else ("p<0.05" if bool(row.get("Nominal_Sig_Adj_Age", False)) else "top hit")
            or_text = row.get("OR_Adj_Age", np.nan)
            trend_p = row.get("P_Adj_Trend", np.nan)
            extra = f" | trend p={trend_p:.3g}" if pd.notna(trend_p) else ""
            ax.set_title(f"Cancer Risk Genotype Composition for {snp_label}\nOR={or_text:.2f} | p={p_val:.3g} ({sig}){extra}",
                         fontsize=8.5, fontweight="bold")
            _style_ax(ax)

        for ax in axes_flat[n_plots:]:
            ax.set_visible(False)

        fig.legend([plt.Rectangle((0, 0), 1, 1, facecolor=_GENO_PLOT_C[g], edgecolor="#4A4A4A", hatch=_GENO_PLOT_HATCH[g]) for g in GENO_ORDER],
                   ["WT", "Het", "Hom"], title="Genotype state (WT / Het / Hom)",
                   loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
        fig.suptitle(tagged_title(f"Cancer Risk Genotype Composition in Tumour and Pooled Control Samples: {cohort}", COMPARATIVE_TAG),
                     fontsize=12, fontweight="bold", y=1.04)
        _safe_tight_layout(fig, rect=[0, 0, 1, 0.95], pad=1.0, h_pad=2.0, w_pad=1.4)
        out_path = out_dir / f"17_CancerRisk_GenotypeComposition_{cohort}.png"
        _save_figure(fig, out_path)


def _first_nonempty_text(series: pd.Series) -> str:
    for value in series:
        text = str(value).strip()
        if text and text.lower() not in {"nan", "none", "na", "<na>"}:
            return text
    return ""


def _min_numeric_across(row: pd.Series, columns: list[str]) -> float:
    values = [pd.to_numeric(row.get(col), errors="coerce") for col in columns if col in row.index]
    values = [value for value in values if pd.notna(value)]
    return float(min(values)) if values else np.nan


def build_significant_overview(summary: pd.DataFrame) -> pd.DataFrame:
    """Collapse the very wide mixed-analysis summary into a supervisor-facing overview."""
    if summary.empty or "Note" in summary.columns:
        return summary.copy()

    rows = []
    for variant_id, sub in summary.groupby("Variant_ID", dropna=False, sort=False):
        gene = _first_nonempty_text(sub.get("Gene", pd.Series(dtype=object)))
        rsid = _first_nonempty_text(sub.get("rsID", pd.Series(dtype=object)))
        variant_label = _first_nonempty_text(sub.get("Variant_Label", pd.Series(dtype=object))) or rsid or str(variant_id)
        analyses = list(dict.fromkeys(sub.get("Analysis_Type", pd.Series(dtype=object)).astype(str).tolist()))
        contexts = []
        best_p_values = []
        best_fdr_values = []
        coverage_flag = bool(sub.get("Coverage_Risk_Flag", pd.Series(dtype=bool)).fillna(False).any()) if "Coverage_Risk_Flag" in sub.columns else False
        coverage_note = _first_nonempty_text(sub.get("Coverage_Risk_Note", pd.Series(dtype=object)))

        for _, row in sub.iterrows():
            analysis_type = str(row.get("Analysis_Type", "")).strip()
            if analysis_type == "Tumour_vs_Control":
                best_p = pd.to_numeric(row.get("P_Value"), errors="coerce")
                best_fdr = pd.to_numeric(row.get("FDR_P_Value"), errors="coerce")
                context = (
                    f"{row.get('Analysis_Group', '')}: tumour {pd.to_numeric(row.get('Freq_Tumour_%'), errors='coerce'):.2f}% vs "
                    f"control {pd.to_numeric(row.get('Freq_Control_%'), errors='coerce'):.2f}% "
                    f"(P={best_p:.3g}; FDR={best_fdr:.3g})"
                )
            elif analysis_type == "Cancer_Risk":
                best_p = _min_numeric_across(row, ["P_Adj_Age", "P_Adj_Trend", "P_Unadj", "P_Hom_vs_WT", "P_Het_vs_WT"])
                best_fdr = _min_numeric_across(row, ["FDR_Adj_Age", "FDR_Adj_Trend", "FDR_Unadj", "FDR_Hom_vs_WT", "FDR_Het_vs_WT"])
                cohort = str(row.get("Cohort", "study")).strip() or "study"
                context = f"{cohort} cancer risk (best P={best_p:.3g}; best FDR={best_fdr:.3g})"
            else:
                best_p = _min_numeric_across(row, ["P_Unadj", "P_Adj_Age", "P_Adj_BMI", "P_Adj_AgeBMI"])
                best_fdr = _min_numeric_across(row, ["FDR_Unadj", "FDR_Adj_Age", "FDR_Adj_BMI", "FDR_Adj_AgeBMI"])
                cohort = str(row.get("Cohort", row.get("Analysis_Group", "study"))).strip() or "study"
                label = str(row.get("Clin_Label", "clinical association")).strip() or "clinical association"
                context = f"{cohort}: {label} (best P={best_p:.3g}; best FDR={best_fdr:.3g})"

            if pd.notna(best_p):
                best_p_values.append(float(best_p))
            if pd.notna(best_fdr):
                best_fdr_values.append(float(best_fdr))
            contexts.append(context)

        best_p = min(best_p_values) if best_p_values else np.nan
        best_fdr = min(best_fdr_values) if best_fdr_values else np.nan
        support = "FDR-supported" if pd.notna(best_fdr) and best_fdr < FDR_THRESHOLD else "Nominal only"
        rows.append({
            "Variant ID": variant_id,
            "rsID": rsid or np.nan,
            "Variant Label": variant_label,
            "Gene": gene,
            "Analyses": "; ".join(analyses),
            "Best P-Value": best_p,
            "Best FDR": best_fdr,
            "Best Support": support,
            "Significant Contexts": " | ".join(dict.fromkeys(contexts)),
            "Coverage Risk Flag": coverage_flag,
            "Coverage Risk Note": coverage_note,
        })

    return pd.DataFrame(rows).sort_values(["Best P-Value", "Variant ID"], kind="stable").reset_index(drop=True)


def parse_args():
    p = argparse.ArgumentParser(description="SNP association analysis for clinical and risk endpoints")
    p.add_argument("--gsdmb", default=str(DEFAULT_GSDMB))
    p.add_argument("--master", default=str(DEFAULT_MASTER))
    p.add_argument("--out_dir", default=str(DEFAULT_OUT))
    p.add_argument("--variant_whitelist", default=(str(DEFAULT_VARIANT_WHITELIST) if DEFAULT_VARIANT_WHITELIST else None),
                   help="Excel whitelist from script 11 defining the thesis SNP backbone")
    p.add_argument("--manifest_endo_tumour", default=str(DEFAULT_MANIFESTS["endometrium-tumour"]),
                   help="Pass-BAM manifest for AU endometrial tumour")
    p.add_argument("--manifest_endo_normal", default=str(DEFAULT_MANIFESTS["endometrium-normal"]),
                   help="Pass-BAM manifest for SP endometrial normal")
    p.add_argument("--manifest_breast_tumour", default=str(DEFAULT_MANIFESTS["breast-tumour"]),
                   help="Pass-BAM manifest for SP breast tumour")
    p.add_argument("--manifest_breast_normal", default=str(DEFAULT_MANIFESTS["breast-normal"]),
                   help="Pass-BAM manifest for SP breast normal")
    p.add_argument("--no_manifests", action="store_true",
                   help="Skip manifest filter (use extraction_flag fallback instead)")
    return p.parse_args()


def main():
    args    = parse_args()
    gsdmb   = Path(args.gsdmb)
    master  = Path(args.master)
    out_dir = Path(args.out_dir)
    validate_file_exists(gsdmb, "Script 17 GSDMB input")
    validate_file_exists(master, "Script 17 master input")
    variant_whitelist = None if args.variant_whitelist in {None, "", "None"} else Path(args.variant_whitelist)
    if variant_whitelist is not None:
        validate_file_exists(variant_whitelist, "Script 17 variant whitelist")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_xlsx = out_dir / "GSDMB_SNP_Clinical_Association_Results.xlsx"

    # Build manifest dict (None if --no_manifests flag set)
    if args.no_manifests:
        manifest_paths = None
    else:
        manifest_paths = {
            "endometrium-tumour": Path(args.manifest_endo_tumour),
            "endometrium-normal": Path(args.manifest_endo_normal),
            "breast-tumour":      Path(args.manifest_breast_tumour),
            "breast-normal":      Path(args.manifest_breast_normal),
        }

    merged = load_and_merge(gsdmb, master, manifest_paths=manifest_paths, whitelist_path=variant_whitelist)
    print_validation_summary(merged, "Sample", "Script 17 merged analysis input", ["Cohort", "Tissue"])
    warning_lookup = build_amplicon_warning_lookup(merged, variant_col="Variant_ID", chrom_col="CHROM", pos_col="POS")
    bmi_category_summary, bmi_data_availability = build_bmi_category_summaries(merged)
    bmi_by_tissue_summary = build_bmi_by_tissue_summary(bmi_category_summary)

    # Analysis 1
    tvh = tumour_vs_control(merged)
    if not tvh.empty:
        validate_percentage_columns(tvh, ["Freq_Tumour_%", "Freq_Control_%"], "Script 17 tumour vs control")
    tvh = attach_amplicon_warning_columns(tvh, warning_lookup)
    tvh = _attach_variant_metadata(tvh, merged)
    tvh = annotate_tumour_vs_control_scheme(tvh)

    # Analysis 2
    breast_clin, endo_clin = clinical_associations(merged)
    breast_clin = attach_amplicon_warning_columns(breast_clin, warning_lookup)
    endo_clin = attach_amplicon_warning_columns(endo_clin, warning_lookup)
    breast_clin = _attach_variant_metadata(breast_clin, merged)
    endo_clin = _attach_variant_metadata(endo_clin, merged)
    breast_clin = annotate_clinical_scheme(breast_clin)
    endo_clin = annotate_clinical_scheme(endo_clin)

    # Exploratory treatment-stratified breast analyses
    breast_tx_summary, breast_tx_clin, breast_tx_surv = exploratory_breast_treatment_analysis(merged)
    breast_tx_clin = attach_amplicon_warning_columns(breast_tx_clin, warning_lookup)
    breast_tx_surv = attach_amplicon_warning_columns(breast_tx_surv, warning_lookup)
    breast_tx_clin = _attach_variant_metadata(breast_tx_clin, merged)
    breast_tx_surv = _attach_variant_metadata(breast_tx_surv, merged)
    breast_tx_clin = annotate_clinical_scheme(breast_tx_clin)

    # Analysis 3 — genotype-dose
    # Run on ALL variants, not just nominal hits, to avoid selection bias.
    # FDR correction within genotype_dose_analysis handles the multiple-testing
    # burden; the FDR_Sig column in the output flags the variants that survive.
    all_breast_vars = breast_clin["Variant_ID"].unique().tolist() if not breast_clin.empty else []
    all_endo_vars   = endo_clin["Variant_ID"].unique().tolist()   if not endo_clin.empty   else []
    breast_dose = genotype_dose_analysis(merged, all_breast_vars, CLINICAL_VARS_BREAST, "Breast_Tumour")
    endo_dose   = genotype_dose_analysis(merged, all_endo_vars,   CLINICAL_VARS_ENDO,   "Endometrial_Tumour")
    breast_dose = attach_amplicon_warning_columns(breast_dose, warning_lookup)
    endo_dose = attach_amplicon_warning_columns(endo_dose, warning_lookup)
    breast_dose = _attach_variant_metadata(breast_dose, merged)
    endo_dose = _attach_variant_metadata(endo_dose, merged)
    breast_dose = annotate_genotype_dose_scheme(breast_dose)
    endo_dose = annotate_genotype_dose_scheme(endo_dose)

    # Analysis 4 — survival
    surv_res = pd.DataFrame()
    km_pages = []
    surv_out = survival_analysis(merged)
    if isinstance(surv_out, tuple):
        surv_res, km_pages = surv_out
    surv_res = attach_amplicon_warning_columns(surv_res, warning_lookup)
    surv_res = _attach_variant_metadata(surv_res, merged)
    surv_res_primary = surv_res[surv_res["Model_Slug"] == "genotype_only"].copy() if not surv_res.empty and "Model_Slug" in surv_res.columns else surv_res.copy()

    # Analysis 5 — cancer risk (case-control)
    risk_res = cancer_risk_analysis(merged)
    risk_res = attach_amplicon_warning_columns(risk_res, warning_lookup)
    risk_res = _attach_variant_metadata(risk_res, merged)
    risk_res = annotate_cancer_risk_scheme(risk_res)

    # Summary
    sig_parts = []
    if not tvh.empty:
        s = tvh[tvh["Nominal_Sig"]].copy()
        s["Analysis_Type"] = "Tumour_vs_Control"
        # Alias P_Value → P_Unadj so the shared sort key works across all result types
        s["P_Unadj"] = s["P_Value"]
        sig_parts.append(s[["Analysis_Type", "Analysis_Group", "Variant_ID", "rsID", "Variant_Label", "Gene",
                             "Association_Scheme", "Comparison_Definition", "Genotype_Definition",
                             "Coverage_Risk_Flag", "Coverage_Risk_Amplicon",
                             "Coverage_Risk_Region", "Coverage_Risk_Note",
                             "Consequence", "IMPACT", "Freq_Tumour_%", "Freq_Control_%",
                             "Odds_Ratio", "Odds_Ratio_CI95", "P_Value", "P_Unadj", "FDR_P_Value",
                             "Nominal_Sig", "FDR_Sig"]])
    for res, label in [
        (breast_clin, "Breast_Tumour"),
        (endo_clin, "Endometrial_Tumour"),
        (breast_tx_clin, "Breast_Treatment_Exploratory"),
    ]:
        if not res.empty:
            sig_mask = res["Nominal_Sig_Unadj"] | res["Nominal_Sig_Adj_Age"] | res["Nominal_Sig_Adj_BMI"] | res["Nominal_Sig_Adj_AgeBMI"]
            s = res[sig_mask].copy(); s["Analysis_Type"] = f"Clinical_{label}"
            keep_cols = ["Analysis_Type", "Treatment_Exposure", "Treatment_Exposure_Var", "Exposure_N",
                         "Cohort", "Variant_ID", "rsID", "Variant_Label", "Gene", "Primary_Row_Model", "Secondary_Additive_Model", "Genotype_Definition",
                         "Coverage_Risk_Flag", "Coverage_Risk_Amplicon",
                         "Coverage_Risk_Region", "Coverage_Risk_Note",
                         "Contrast", "N_WT", "N_Het", "N_Hom",
                         "Clin_Label", "Clin_Type",
                         "Test_Unadj", "OR_Unadj", "OR_Unadj_CI95", "P_Unadj", "FDR_Unadj",
                         "Nominal_Sig_Unadj", "FDR_Sig_Unadj",
                         "Test_Adj_Age", "OR_Adj_Age", "OR_Adj_Age_CI95", "P_Adj_Age", "FDR_Adj_Age",
                         "Nominal_Sig_Adj_Age", "FDR_Sig_Adj_Age", "Adj_Age_Note",
                         "Test_Adj_BMI", "OR_Adj_BMI", "OR_Adj_BMI_CI95", "P_Adj_BMI", "FDR_Adj_BMI",
                         "Nominal_Sig_Adj_BMI", "FDR_Sig_Adj_BMI", "Adj_BMI_Note",
                         "Test_Adj_AgeBMI", "OR_Adj_AgeBMI", "OR_Adj_AgeBMI_CI95", "P_Adj_AgeBMI", "FDR_Adj_AgeBMI",
                         "Nominal_Sig_Adj_AgeBMI", "FDR_Sig_Adj_AgeBMI", "Adj_AgeBMI_Note",
                         "P_Unadj_Additive", "OR_Unadj_Additive", "OR_Unadj_Additive_CI95",
                         "OR_Adj_Age_Additive", "P_Adj_Age_Additive",
                         "OR_Adj_BMI_Additive", "P_Adj_BMI_Additive",
                         "OR_Adj_AgeBMI_Additive", "P_Adj_AgeBMI_Additive"]
            s = s[[c for c in keep_cols if c in s.columns]]
            sig_parts.append(s)
    if not risk_res.empty:
        risk_sig = risk_res[risk_res["Nominal_Sig_Unadj"] | risk_res["Nominal_Sig_Adj_Age"] | risk_res.get("Nominal_Sig_Adj_Trend", False)].copy()
        risk_sig["Analysis_Type"] = "Cancer_Risk"
        sig_parts.append(risk_sig)
    if not surv_res.empty:
        surv_mask = pd.Series(False, index=surv_res.index)
        for col in ["Nominal_Sig_Cox_Additive", "FDR_Sig_Cox_Additive", "Nominal_Sig_LogRank", "FDR_Sig_LogRank"]:
            if col in surv_res.columns:
                surv_mask |= surv_res[col].fillna(False).astype(bool)
        surv_sig = surv_res[surv_mask].copy()
        if not surv_sig.empty:
            surv_sig["Analysis_Type"] = "Survival"
            surv_sig["P_Unadj"] = surv_sig.get("P_Cox_Additive", np.nan)
            surv_sig["FDR_Unadj"] = surv_sig.get("FDR_Cox_Additive", np.nan)
            surv_sig["Nominal_Sig"] = surv_sig.get("Nominal_Sig_Cox_Additive", False)
            surv_sig["FDR_Sig"] = surv_sig.get("FDR_Sig_Cox_Additive", False)
            sig_parts.append(surv_sig)
    if not breast_tx_surv.empty:
        tx_surv_sig = breast_tx_surv[(breast_tx_surv["P_Cox_Additive"] < 0.05) | (breast_tx_surv["P_Cox_Het"] < 0.05) | (breast_tx_surv["P_Cox_Hom"] < 0.05)].copy()
        if not tx_surv_sig.empty:
            tx_surv_sig["Analysis_Type"] = "Treatment_Stratified_Survival"
            sig_parts.append(tx_surv_sig)
    summary = (pd.concat(sig_parts, ignore_index=True).sort_values("P_Unadj")
               if sig_parts else pd.DataFrame({"Note": ["No nominally significant results."]}))
    summary = attach_amplicon_warning_columns(summary, warning_lookup)
    summary = _attach_variant_metadata(summary, merged)
    significant_overview = build_significant_overview(summary)

    workbook_readme = pd.DataFrame([
        {
            "Section": "Primary thesis outputs",
            "Details": "The main stage-17 sheets for thesis use are tumour_vs_control, breast_clinical_assoc, endo_clinical_assoc, genotype_dose, summary_significant, Significant Overview, and sample_manifest.",
        },
        {
            "Section": "Secondary outputs",
            "Details": "Treatment-stratified, survival, cancer-risk, and BMI outputs are retained for supplementary interpretation only.",
        },
        {
            "Section": "Design note",
            "Details": "This script follows the cleaned unpaired workflow. Historical source sheet labels do not imply matched tumour-normal modelling.",
        },
    ])

    workbook_index = pd.DataFrame([
        {"Order": 1, "Sheet": "README", "Priority": "Primary", "Purpose": "Workbook guide and scope note"},
        {"Order": 2, "Sheet": "At_A_Glance", "Priority": "Primary", "Purpose": "Compact row counts for the main and secondary result layers"},
        {"Order": 3, "Sheet": "summary_significant", "Priority": "Primary", "Purpose": "Combined significant SNP associations across the retained workflows"},
        {"Order": 4, "Sheet": "Significant Overview", "Priority": "Primary", "Purpose": "Supervisor-facing overview of the strongest SNP signals"},
        {"Order": 5, "Sheet": "tumour_vs_control", "Priority": "Primary", "Purpose": "Tumour-versus-pooled-control SNP associations"},
        {"Order": 6, "Sheet": "breast_clinical_assoc", "Priority": "Primary", "Purpose": "Breast tumour clinical SNP associations"},
        {"Order": 7, "Sheet": "endo_clinical_assoc", "Priority": "Primary", "Purpose": "Endometrial tumour clinical SNP associations"},
        {"Order": 8, "Sheet": "genotype_dose", "Priority": "Primary", "Purpose": "Genotype-dose follow-up for the main tumour cohorts"},
        {"Order": 9, "Sheet": "sample_manifest", "Priority": "Primary", "Purpose": "Manifest of sequenced samples entering the analysis"},
        {"Order": 10, "Sheet": "breast_tx_summary", "Priority": "Secondary", "Purpose": "Breast treatment-stratified sample summary"},
        {"Order": 11, "Sheet": "breast_tx_clinical", "Priority": "Secondary", "Purpose": "Breast treatment-stratified clinical SNP associations"},
        {"Order": 12, "Sheet": "breast_tx_survival", "Priority": "Secondary", "Purpose": "Breast treatment-stratified survival outputs"},
        {"Order": 13, "Sheet": "survival_cox", "Priority": "Secondary", "Purpose": "Secondary survival modelling outputs"},
        {"Order": 14, "Sheet": "cancer_risk", "Priority": "Secondary", "Purpose": "Secondary case-control risk summaries"},
        {"Order": 15, "Sheet": "bmi_category_summary", "Priority": "Secondary", "Purpose": "BMI category summary table"},
        {"Order": 16, "Sheet": "bmi_by_tissue_summary", "Priority": "Secondary", "Purpose": "BMI category summary by tissue family"},
        {"Order": 17, "Sheet": "bmi_data_availability", "Priority": "Secondary", "Purpose": "BMI availability and exclusions"},
    ])

    at_a_glance = pd.DataFrame([
        {"Section": "Tumour_vs_Control", "Rows": int(len(tvh)) if not tvh.empty else 0, "Nominal_or_Significant_Rows": int(tvh["Nominal_Sig"].sum()) if not tvh.empty and "Nominal_Sig" in tvh.columns else 0},
        {"Section": "Breast_Clinical", "Rows": int(len(breast_clin)) if not breast_clin.empty else 0, "Nominal_or_Significant_Rows": int(breast_clin["Nominal_Sig_Unadj"].sum()) if not breast_clin.empty and "Nominal_Sig_Unadj" in breast_clin.columns else 0},
        {"Section": "Endometrial_Clinical", "Rows": int(len(endo_clin)) if not endo_clin.empty else 0, "Nominal_or_Significant_Rows": int(endo_clin["Nominal_Sig_Unadj"].sum()) if not endo_clin.empty and "Nominal_Sig_Unadj" in endo_clin.columns else 0},
        {"Section": "Genotype_Dose", "Rows": int(len(pd.concat([breast_dose, endo_dose], ignore_index=True))) if (not breast_dose.empty or not endo_dose.empty) else 0, "Nominal_or_Significant_Rows": 0},
        {"Section": "Summary_Significant", "Rows": int(len(summary)) if not summary.empty else 0, "Nominal_or_Significant_Rows": int(len(summary)) if not summary.empty and "Note" not in summary.columns else 0},
        {"Section": "Secondary_Survival", "Rows": int(len(surv_res)) if not surv_res.empty else 0, "Nominal_or_Significant_Rows": int(surv_res["Nominal_Sig_Cox_Additive"].sum()) if not surv_res.empty and "Nominal_Sig_Cox_Additive" in surv_res.columns else 0},
        {"Section": "Secondary_Cancer_Risk", "Rows": int(len(risk_res)) if not risk_res.empty else 0, "Nominal_or_Significant_Rows": int(risk_res["Nominal_Sig_Unadj"].sum()) if not risk_res.empty and "Nominal_Sig_Unadj" in risk_res.columns else 0},
    ])

    # Sample manifest
    keep_clin = ["canon__recurrence_flag", "canon__distant_mets_flag",
                 "canon__exitus_flag",
                 "BREAST_TX_CHEMOTHERAPY_BIN", "BREAST_TX_ANTI_HER2_BIN",
                 "BREAST_TX_ENDOCRINE_BIN", "BREAST_TX_RADIOTHERAPY_BIN",
                 "canon__pd_flag", "canon__exitus_flag",
                 "canon__os_months", "canon__pfs_months",
                 "canon__risk_group"]
    manifest = merged.drop_duplicates("Sample")[[
        c for c in ["Sample", "snp_code", "sample_id", "Cohort", "Tissue", "sheet",
                    "tumour_normal", "is_replicate"] + keep_clin
        if c in merged.columns
    ]].sort_values(["Cohort", "Tissue", "Sample"]).reset_index(drop=True)

    # Write Excel
    print("=== Writing output Excel ===")
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        workbook_readme.to_excel(xw,                    sheet_name="README",                 index=False)
        workbook_index.to_excel(xw,                     sheet_name="Workbook_Index",         index=False)
        at_a_glance.to_excel(xw,                        sheet_name="At_A_Glance",            index=False)
        summary.to_excel(xw,                            sheet_name="summary_significant",    index=False)
        significant_overview.to_excel(xw,               sheet_name="Significant Overview",   index=False)
        if not tvh.empty:       tvh.to_excel(xw,         sheet_name="tumour_vs_control",     index=False)
        if not breast_clin.empty: breast_clin.to_excel(xw,sheet_name="breast_clinical_assoc", index=False)
        if not endo_clin.empty:   endo_clin.to_excel(xw,  sheet_name="endo_clinical_assoc",   index=False)
        if not breast_tx_summary.empty: breast_tx_summary.to_excel(xw, sheet_name="breast_tx_summary", index=False)
        if not breast_tx_clin.empty: breast_tx_clin.to_excel(xw, sheet_name="breast_tx_clinical", index=False)
        if not breast_tx_surv.empty: breast_tx_surv.to_excel(xw, sheet_name="breast_tx_survival", index=False)
        dose_all = pd.concat([breast_dose, endo_dose], ignore_index=True)
        if not dose_all.empty:  dose_all.to_excel(xw,    sheet_name="genotype_dose",          index=False)
        if not surv_res.empty:  surv_res.to_excel(xw,    sheet_name="survival_cox",           index=False)
        if not risk_res.empty:  risk_res.to_excel(xw,    sheet_name="cancer_risk",             index=False)
        if not bmi_category_summary.empty: bmi_category_summary.to_excel(xw, sheet_name="bmi_category_summary", index=False)
        if not bmi_by_tissue_summary.empty: bmi_by_tissue_summary.to_excel(xw, sheet_name="bmi_by_tissue_summary", index=False)
        if not bmi_data_availability.empty: bmi_data_availability.to_excel(xw, sheet_name="bmi_data_availability", index=False)
        manifest.to_excel(xw,                            sheet_name="sample_manifest",        index=False)
    print(f"  Saved: {out_xlsx}\n")

    # Plots
    print("=== Generating plots ===")
    supplementary_dir = out_dir / "supplementary_figures"
    supplementary_dir.mkdir(parents=True, exist_ok=True)
    significant_dir = supplementary_dir / "significant_only_figures"
    significant_dir.mkdir(parents=True, exist_ok=True)
    legacy_root_dir = supplementary_dir / "_legacy_root_outputs"

    def _archive_root_output(path: Path) -> None:
        if not path.exists():
            return
        legacy_root_dir.mkdir(parents=True, exist_ok=True)
        if path.is_dir():
            target = legacy_root_dir / path.name
            counter = 1
            while target.exists():
                target = legacy_root_dir / f"{path.name}__archived_{counter}"
                counter += 1
        else:
            target = legacy_root_dir / path.name
            counter = 1
            while target.exists():
                target = legacy_root_dir / f"{path.stem}__archived_{counter}{path.suffix}"
                counter += 1
        shutil.move(str(path), str(target))

    keep_root_figures = {
        "17_SNP_Heatmap_Breast_MainText.png",
        "17_SNP_Heatmap_Endometrial_MainText.png",
    }
    for stale in out_dir.glob("17_*.png"):
        if stale.name not in keep_root_figures:
            _archive_root_output(stale)
    for stale in out_dir.glob("17_*.pdf"):
        _archive_root_output(stale)
    for stale in out_dir.glob("GSDMB_SNP_*.pdf"):
        _archive_root_output(stale)
    _archive_root_output(out_dir / "GSDMB_SNP_Association_Summary_Panel.png")
    _archive_root_output(out_dir / "significant_only_figures")
    make_volcano_plots(tvh, supplementary_dir)
    make_heatmaps(breast_clin, endo_clin, supplementary_dir, significant_only=False)
    make_heatmaps(breast_clin, endo_clin, significant_dir, significant_only=True)
    make_main_text_heatmaps(breast_clin, endo_clin, out_dir)
    make_forest_plots(breast_clin, endo_clin, supplementary_dir)
    make_distribution_plots(breast_clin, endo_clin, merged, supplementary_dir)
    make_binary_bar_plots(breast_clin, endo_clin, merged, supplementary_dir)
    make_bmi_category_plot(bmi_category_summary, bmi_data_availability, supplementary_dir)
    make_km_plots(surv_res, km_pages, supplementary_dir)
    make_risk_forest_plot(risk_res, supplementary_dir)
    make_summary_panel(tvh, breast_clin, endo_clin, surv_res_primary, supplementary_dir)
    make_genotype_distribution_plots(breast_dose, endo_dose, merged, supplementary_dir)
    make_genotype_binary_bar_plots(breast_dose, endo_dose, merged, supplementary_dir)
    make_genotype_composition_plots(tvh, merged, supplementary_dir)
    make_risk_genotype_composition_plots(risk_res, merged, supplementary_dir)

    # Final summary
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"  Samples in manifest (sequenced): {merged[~merged['is_replicate']]['Sample'].nunique()}")
    print(f"  Variants tested  : {merged['Variant_ID'].nunique()}")
    if not tvh.empty:
        print(f"  Tumour vs control   — nominal: {tvh['Nominal_Sig'].sum()} | FDR<{FDR_THRESHOLD}: {tvh['FDR_Sig'].sum()}")
    if not breast_clin.empty:
        print(f"  Breast clinical     — nominal unadj: {breast_clin['Nominal_Sig_Unadj'].sum()} | FDR: {breast_clin['FDR_Sig_Unadj'].sum()} | nominal age-adj: {breast_clin['Nominal_Sig_Adj_Age'].sum()} | nominal BMI-adj: {breast_clin['Nominal_Sig_Adj_BMI'].sum()} | nominal age+BMI-adj: {breast_clin['Nominal_Sig_Adj_AgeBMI'].sum()}")
    if not endo_clin.empty:
        print(f"  Endometrial clinical— nominal unadj: {endo_clin['Nominal_Sig_Unadj'].sum()} | FDR: {endo_clin['FDR_Sig_Unadj'].sum()} | nominal age-adj: {endo_clin['Nominal_Sig_Adj_Age'].sum()} | nominal BMI-adj: {endo_clin['Nominal_Sig_Adj_BMI'].sum()} | nominal age+BMI-adj: {endo_clin['Nominal_Sig_Adj_AgeBMI'].sum()}")
    if not bmi_data_availability.empty:
        print(f"  BMI category source  — {BMI_SOURCE_COLUMN}")
        for _, row in bmi_data_availability.iterrows():
            print(
                f"  BMI ({row['Cohort_Display']}) — usable: {int(row['Samples_With_Usable_BMI_Category'])} / {int(row['Total_Samples'])} | "
                f"<20 excluded: {int(row['Samples_Below_20_Excluded'])} | missing/non-numeric: {int(row['Samples_Missing_or_NonNumeric_BMI'])}"
            )
    if not surv_res.empty:
        if "Model" in surv_res.columns:
            model_bits = [f"{m}: {int((surv_res['Model'] == m).sum())}" for m in surv_res["Model"].dropna().unique()]
            print(f"  Survival (Cox)      — {len(surv_res)} tests | " + " | ".join(model_bits))
        else:
            print(f"  Survival (Cox)      — {len(surv_res)} tests | "
                  f"cohorts: {', '.join(surv_res['Cohort'].unique())}")
    if not risk_res.empty:
        for cohort in risk_res["Cohort"].unique():
            sub = risk_res[risk_res["Cohort"] == cohort]
            print(f"  Cancer risk ({cohort}) — "
                  f"{sub['Nominal_Sig_Unadj'].sum()} nominal (unadj) | "
                  f"{sub['FDR_Sig_Unadj'].sum()} FDR (unadj) | "
                  f"{sub['Nominal_Sig_Adj_Age'].sum()} nominal (age-adj)")
    print(f"\n  Results : {out_xlsx}")
    print(f"  Plots   : {out_dir}")


if __name__ == "__main__":
    main()
