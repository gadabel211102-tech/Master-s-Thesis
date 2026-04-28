"""
Script 12: SNP Enrichment Testing
=================================

Purpose
-------
Test whether common SNP genotype classes differ between tumour samples and the
pooled normal control arm, both globally and within each tumour cohort.

Statistical approach
--------------------
For each SNP, the script constructs a 3x2 WT/Het/Hom contingency table and
applies a chi-square test of independence for the workbook summary. Carrier-
based volcano plots are retained for the ALT carrier vs non-carrier view,
while separate genotype volcanoes compare Het vs WT and Hom vs WT against the
pooled normal controls. False-discovery-rate correction is then applied within
each analysis group.

Control definition
------------------
Cohort-specific tumour analyses use the pooled control group agreed for this
project: breast normal plus endometrium normal.
"""

import pandas as pd
import textwrap
import numpy as np
from scipy.stats import chi2, chi2_contingency, fisher_exact, norm
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
import seaborn as sns
from statsmodels.stats.multitest import multipletests

from annotation_genotype_bridge import build_allele_genotype_long
from figure_style import COMPARATIVE_TAG, GENOTYPE_COLORS, IMPACT_COLORS, arm_color, cohort_color, tagged_title
from pipeline_utils import (
    attach_amplicon_warning_columns,
    build_amplicon_warning_lookup,
    common_nfe_variant_mask,
    ensure_directory,
    get_paths,
    get_thresholds,
    load_common_snp_ld_reference,
    pooled_control_frame,
)
from pipeline_validation import print_validation_summary, validate_file_exists, validate_nonempty, validate_percentage_columns
from variant_label_utils import build_variant_display_table as shared_build_variant_display_table

# --- 1. PATH CONFIGURATION ---
PATHS = get_paths()
THRESHOLDS = get_thresholds()
output_dir = ensure_directory(PATHS["snp_enrichment_dir"])
supplementary_output_dir = ensure_directory(output_dir / "supplementary_results")
input_file = str(PATHS["annotated_report"])
output_xlsx = str(output_dir / "GSDMB_SNP_Enrichment_Results.xlsx")

# Separate volcano outputs
output_plot_raw = str(output_dir / "GSDMB_SNP_Enrichment_Volcano_RAW.png")
output_plot_fdr = str(output_dir / "GSDMB_SNP_Enrichment_Volcano_FDR.png")
output_plot_het_raw = str(output_dir / "GSDMB_SNP_Enrichment_Volcano_HetVsWT_RAW.png")
output_plot_het_fdr = str(output_dir / "GSDMB_SNP_Enrichment_Volcano_HetVsWT_FDR.png")
output_plot_hom_raw = str(output_dir / "GSDMB_SNP_Enrichment_Volcano_HomVsWT_RAW.png")
output_plot_hom_fdr = str(output_dir / "GSDMB_SNP_Enrichment_Volcano_HomVsWT_FDR.png")
output_plot_genotype_heatmap = str(output_dir / "GSDMB_SNP_Genotype_Enrichment_Heatmap.png")
output_plot_genotype_composition_prefix = str(output_dir / "GSDMB_SNP_Genotype_Composition")
output_mutation_xlsx = str(supplementary_output_dir / "GSDMB_Mutation_Enrichment_Results.xlsx")
output_plot_mutation_burden = str(supplementary_output_dir / "GSDMB_Mutation_Burden_Forest.png")
output_plot_mutation_hotspots = str(supplementary_output_dir / "GSDMB_Mutation_Hotspots.png")
output_ld_summary_xlsx = str(output_dir / "GSDMB_SNP_Enrichment_LD_Block_Summary.xlsx")
output_ppt_summary_xlsx = str(supplementary_output_dir / "GSDMB_SNP_Enrichment_PPT_Summary.xlsx")
LD_R2_THRESHOLD = 0.80
RAW_P_THRESHOLD = 0.05
FDR_THRESHOLD = float(THRESHOLDS.get("fdr_threshold", 0.10))
TARGET_RSIDS = {"rs11078928", "rs869402"}
GENOTYPE_ORDER = ["WT", "Het_ALT", "Hom_ALT"]
GENOTYPE_LABELS = {"WT": "WT", "Het_ALT": "Het", "Hom_ALT": "Hom"}
MUTATION_GENOTYPE_MIN_ALT_COUNT = 5


def _genotype_count_map(series: pd.Series) -> dict[str, int]:
    counts = series.astype(str).value_counts()
    return {genotype: int(counts.get(genotype, 0)) for genotype in GENOTYPE_ORDER}


def _mean_alt_dose(counts: dict[str, int]) -> float | None:
    n_callable = sum(counts.values())
    if not n_callable:
        return None
    return (counts["Het_ALT"] + 2 * counts["Hom_ALT"]) / n_callable


def _pct(value: int, total: int) -> float | None:
    return (value / total) * 100 if total else None


def _short_label(value, max_len=28):
    text = " ".join(str(value).split())
    return text if len(text) <= max_len else text[: max_len - 1].rstrip() + "?"


def _first_nonempty(series):
    for value in series:
        text = " ".join(str(value).split()).strip()
        if text and text.lower() not in {"nan", "none", "na", "<na>"}:
            return text
    return ""


def _short_hgvsc(value):
    text = _first_nonempty([value])
    return text.split(":")[-1] if ":" in text else text


def build_variant_display_table(df, variant_col, existing_col, symbol_col, hgvsc_col, chrom_col, pos_col, ref_col, alt_col):
    out = shared_build_variant_display_table(
        df,
        variant_col=variant_col,
        existing_col=existing_col,
        symbol_cols=[symbol_col] if symbol_col else [],
        hgvs_col=hgvsc_col,
        chrom_col=chrom_col,
        pos_col=pos_col,
        ref_col=ref_col,
        alt_col=alt_col,
        output_id_col="Variant_ID",
    )
    if "HGVS_Short" in out.columns:
        out = out.rename(columns={"HGVS_Short": "HGVSc"})
    return out


def build_significant_overview_table(final_res):
    sig = final_res[final_res["Raw_Significant"].fillna(False)].copy()
    if sig.empty:
        return pd.DataFrame({"Note": ["No nominally significant stage-12 SNP genotype-association results were available."]})
    sig["Significance Support"] = np.where(sig["FDR_Significant"].fillna(False), "FDR-supported", "Nominal only")
    keep = [
        "Analysis_Group", "Display_Label", "rsID", "SNP_ID", "Symbol", "Consequence", "Impact",
        "Tumour_WT", "Tumour_Het_ALT", "Tumour_Hom_ALT",
        "Control_WT", "Control_Het_ALT", "Control_Hom_ALT",
        "Tumour_WT_Pct_Callable", "Tumour_Het_Pct_Callable", "Tumour_Hom_Pct_Callable",
        "Control_WT_Pct_Callable", "Control_Het_Pct_Callable", "Control_Hom_Pct_Callable",
        "Tumour_Callable_Samples", "Control_Callable_Samples",
        "Tumour_Mean_ALT_Dose", "Control_Mean_ALT_Dose", "Dose_Difference",
        "P_Value", "FDR_P_Value",
        "Significance Support", "Comparison_Definition", "Genotype_Definition",
        "Coverage_Risk_Flag", "Coverage_Risk_Note",
    ]
    keep = [col for col in keep if col in sig.columns]
    return sig[keep].sort_values(["Analysis_Group", "P_Value", "Display_Label"], kind="stable").reset_index(drop=True)


def _haldane_odds_ratio(a: int, b: int, c: int, d: int) -> tuple[float, str]:
    if min(a, b, c, d) == 0:
        a, b, c, d = [value + 0.5 for value in (a, b, c, d)]
        return float((a * d) / (b * c)), "Haldane-Anscombe"
    if b * c == 0:
        return float("inf"), "Standard"
    return float((a * d) / (b * c)), "Standard"


def _apply_grouped_fdr(df: pd.DataFrame, p_col: str, out_col: str, group_cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    out[out_col] = np.nan
    if out.empty or p_col not in out.columns:
        return out
    for _, idx in out.groupby(group_cols, dropna=False).groups.items():
        pvals = pd.to_numeric(out.loc[idx, p_col], errors="coerce")
        valid = pvals.notna()
        if valid.sum() == 0:
            continue
        if valid.sum() == 1:
            out.loc[pvals[valid].index, out_col] = pvals[valid].iloc[0]
            continue
        _, adj, _, _ = multipletests(pvals[valid], method="fdr_bh")
        out.loc[pvals[valid].index, out_col] = adj
    return out


def _comparison_definition(analysis_group: str) -> str:
    mapping = {
        "Breast": "Breast tumours versus pooled healthy controls (breast normal + endometrial normal)",
        "Endometrium": "Endometrial tumours versus pooled healthy controls (breast normal + endometrial normal)",
        "Global": "All tumours pooled versus pooled healthy controls (breast normal + endometrial normal)",
    }
    return mapping.get(analysis_group, "Tumour arm versus pooled healthy controls")


def _cochran_armitage_trend_test(tumour_counts: dict[str, int], control_counts: dict[str, int]) -> dict[str, object]:
    tumour = np.array([tumour_counts[g] for g in GENOTYPE_ORDER], dtype=float)
    control = np.array([control_counts[g] for g in GENOTYPE_ORDER], dtype=float)
    totals = tumour + control
    scores = np.array([0.0, 1.0, 2.0], dtype=float)
    total_n = float(totals.sum())
    tumour_n = float(tumour.sum())
    control_n = float(control.sum())
    observed_levels = int((totals > 0).sum())

    if total_n <= 1 or tumour_n == 0 or control_n == 0 or observed_levels < 2:
        return {
            "Genotype_Trend_Z": np.nan,
            "Genotype_Trend_P_Value": np.nan,
            "Genotype_Trend_Method": "Skipped",
            "Genotype_Trend_Note": "Genotype-dose trend not estimable because fewer than two genotype-dose levels were observed.",
        }

    mean_score = float(np.sum(totals * scores) / total_n)
    numerator = float(np.sum(scores * (tumour - (totals * tumour_n / total_n))))
    variance = float((tumour_n * control_n * np.sum(totals * (scores - mean_score) ** 2)) / (total_n * (total_n - 1)))
    if not np.isfinite(variance) or variance <= 0:
        return {
            "Genotype_Trend_Z": np.nan,
            "Genotype_Trend_P_Value": np.nan,
            "Genotype_Trend_Method": "Skipped",
            "Genotype_Trend_Note": "Genotype-dose trend variance was not estimable from the observed genotype distribution.",
        }

    z_score = numerator / np.sqrt(variance)
    p_value = float(2 * norm.sf(abs(z_score)))
    return {
        "Genotype_Trend_Z": float(z_score),
        "Genotype_Trend_P_Value": p_value,
        "Genotype_Trend_Method": "Cochran-Armitage trend test",
        "Genotype_Trend_Note": "Ordinal genotype-dose association using ALT-allele scores 0 / 1 / 2.",
    }


def _genotype_distribution_test(tumour_counts: dict[str, int], control_counts: dict[str, int]) -> dict[str, object]:
    tumour = np.array([tumour_counts[g] for g in GENOTYPE_ORDER], dtype=int)
    control = np.array([control_counts[g] for g in GENOTYPE_ORDER], dtype=int)
    table = np.vstack([tumour, control])
    observed_cols = table.sum(axis=0) > 0
    represented = [GENOTYPE_LABELS[genotype] for genotype, keep in zip(GENOTYPE_ORDER, observed_cols, strict=False) if keep]
    reduced = table[:, observed_cols]

    if reduced.shape[1] < 2:
        return {
            "Genotype_Distribution_Statistic": np.nan,
            "Genotype_Distribution_P_Value": np.nan,
            "Genotype_Distribution_DOF": np.nan,
            "Genotype_Distribution_Method": "Skipped",
            "Genotype_Distribution_Note": "Genotype-distribution association was skipped because fewer than two genotype classes were observed.",
            "Expected_Min_Cell": np.nan,
            "Sparse_Cell_Flag": False,
            "Genotype_Classes_Observed": ", ".join(represented),
        }

    if reduced.shape[1] == 2:
        odds_ratio, p_value = fisher_exact(reduced, alternative="two-sided")
        odds_ratio = float(odds_ratio) if np.isfinite(odds_ratio) else np.nan
        return {
            "Genotype_Distribution_Statistic": odds_ratio,
            "Genotype_Distribution_P_Value": float(p_value),
            "Genotype_Distribution_DOF": 1.0,
            "Genotype_Distribution_Method": "Fisher exact on reduced 2x2 genotype table",
            "Genotype_Distribution_Note": f"Only two genotype classes were observed overall ({', '.join(represented)}), so a reduced exact test was used.",
            "Expected_Min_Cell": np.nan,
            "Sparse_Cell_Flag": True,
            "Genotype_Classes_Observed": ", ".join(represented),
        }

    chi2_stat, p_value, dof, expected = chi2_contingency(reduced, correction=False)
    expected_min = float(np.min(expected)) if expected.size else np.nan
    sparse = bool(np.isfinite(expected_min) and expected_min < 5)
    if sparse:
        return {
            "Genotype_Distribution_Statistic": float(chi2_stat),
            "Genotype_Distribution_P_Value": np.nan,
            "Genotype_Distribution_DOF": float(dof),
            "Genotype_Distribution_Method": "Skipped sparse 2x3 genotype table",
            "Genotype_Distribution_Note": f"All three genotype classes were observed ({', '.join(represented)}), but the minimum expected cell count was {expected_min:.2f}; the full-table association was not exported as a primary p-value.",
            "Expected_Min_Cell": expected_min,
            "Sparse_Cell_Flag": True,
            "Genotype_Classes_Observed": ", ".join(represented),
        }

    return {
        "Genotype_Distribution_Statistic": float(chi2_stat),
        "Genotype_Distribution_P_Value": float(p_value),
        "Genotype_Distribution_DOF": float(dof),
        "Genotype_Distribution_Method": "Pearson chi-square on 2x3 genotype table",
        "Genotype_Distribution_Note": f"All three genotype classes were analysed jointly ({', '.join(represented)}).",
        "Expected_Min_Cell": expected_min,
        "Sparse_Cell_Flag": False,
        "Genotype_Classes_Observed": ", ".join(represented),
    }


def _binary_genotype_contrast(case_alt: int, case_ref: int, control_alt: int, control_ref: int, label: str) -> dict[str, object]:
    total_alt = int(case_alt + control_alt)
    total_ref = int(case_ref + control_ref)
    if total_alt == 0 or total_ref == 0:
        return {
            f"{label}_Odds_Ratio": np.nan,
            f"{label}_OR_Method": "",
            f"{label}_P_Value": np.nan,
            f"{label}_Method": f"{label.replace('_', ' ')} skipped",
            f"{label}_Note": "Contrast not estimable because the contrast or reference genotype was absent overall.",
        }

    odds_ratio, or_method = _haldane_odds_ratio(case_alt, case_ref, control_alt, control_ref)
    _, p_value = fisher_exact([[case_alt, case_ref], [control_alt, control_ref]], alternative="two-sided")
    return {
        f"{label}_Odds_Ratio": odds_ratio,
        f"{label}_OR_Method": or_method,
        f"{label}_P_Value": float(p_value),
        f"{label}_Method": f"{label.replace('_', ' ')} Fisher exact",
        f"{label}_Note": f"{label.replace('_', ' ')} contrasted against WT among callable samples.",
    }


def _best_available_p(row: pd.Series, columns: list[str]) -> float | None:
    values = pd.to_numeric(row[columns], errors="coerce")
    values = values[np.isfinite(values)]
    if values.empty:
        return None
    return float(values.min())


def _variant_group_sets(manifest: pd.DataFrame, analysis_group: str) -> tuple[set[str], set[str]]:
    pooled_controls = set(manifest.loc[manifest["Tissue"].astype(str).eq("Healthy"), "Sample"].astype(str))
    if analysis_group == "Global":
        tumour_samples = set(manifest.loc[manifest["Tissue"].astype(str).eq("Tumour"), "Sample"].astype(str))
    elif analysis_group == "Breast":
        tumour_samples = set(
            manifest.loc[
                manifest["Tissue"].astype(str).eq("Tumour")
                & manifest["Cohort"].astype(str).str.contains("Breast", case=False, na=False),
                "Sample",
            ].astype(str)
        )
    else:
        tumour_samples = set(
            manifest.loc[
                manifest["Tissue"].astype(str).eq("Tumour")
                & manifest["Cohort"].astype(str).str.contains("Endometri", case=False, na=False),
                "Sample",
            ].astype(str)
        )
    return tumour_samples, pooled_controls


def _save_publication_figure(fig: plt.Figure, path: str) -> None:
    fig.savefig(path, dpi=400, bbox_inches="tight", pad_inches=0.28, facecolor="white")
    plt.close(fig)


def compute_snp_genotype_enrichment(common_df: pd.DataFrame, display_lookup: pd.DataFrame, warning_lookup: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cohort_list = ["Global", "Breast", "Endometrium"]
    cohort_filters = {"Breast": "Breast", "Endometrium": "Endometri"}
    rows: list[dict[str, object]] = []

    for analysis_group in cohort_list:
        if analysis_group == "Global":
            cohort_df = pooled_control_frame(common_df, "Cohort", "Tissue")
        else:
            cohort_df = pooled_control_frame(common_df, "Cohort", "Tissue", analysis_group, cohort_filters[analysis_group])
        if cohort_df.empty:
            continue

        for variant_id, sub in cohort_df.groupby("Variant_Key", sort=False):
            tumour = sub[sub["Tissue"] == "Tumour"].copy()
            control = sub[sub["Tissue"] == "Healthy"].copy()
            tumour_callable = tumour[tumour["Callable_Allele_Genotype"]].copy()
            control_callable = control[control["Callable_Allele_Genotype"]].copy()
            tumour_callable_n = int(tumour_callable["Sample"].nunique())
            control_callable_n = int(control_callable["Sample"].nunique())
            if tumour_callable_n == 0 or control_callable_n == 0:
                continue

            tumour_counts = _genotype_count_map(tumour_callable["Variant_Genotype_Status"])
            control_counts = _genotype_count_map(control_callable["Variant_Genotype_Status"])
            distribution = _genotype_distribution_test(tumour_counts, control_counts)
            trend = _cochran_armitage_trend_test(tumour_counts, control_counts)
            het_contrast = _binary_genotype_contrast(
                tumour_counts["Het_ALT"], tumour_counts["WT"], control_counts["Het_ALT"], control_counts["WT"], "Het_vs_WT"
            )
            hom_contrast = _binary_genotype_contrast(
                tumour_counts["Hom_ALT"], tumour_counts["WT"], control_counts["Hom_ALT"], control_counts["WT"], "Hom_vs_WT"
            )

            tumour_mean_dose = _mean_alt_dose(tumour_counts)
            control_mean_dose = _mean_alt_dose(control_counts)
            dose_diff = None if tumour_mean_dose is None or control_mean_dose is None else tumour_mean_dose - control_mean_dose
            primary_p = distribution["Genotype_Distribution_P_Value"]
            primary_method = distribution["Genotype_Distribution_Method"]
            primary_note = distribution["Genotype_Distribution_Note"]
            if pd.isna(primary_p):
                primary_p = trend["Genotype_Trend_P_Value"]
                primary_method = trend["Genotype_Trend_Method"]
                primary_note = f"{distribution['Genotype_Distribution_Note']} Primary genotype ranking therefore falls back to the genotype-dose trend test."

            consequence_col = "Detailed_Consequence" if "Detailed_Consequence" in sub.columns else "Consequence"
            impact_col = "Impact_Severity" if "Impact_Severity" in sub.columns else "Impact"
            rows.append({
                "Analysis_Group": analysis_group,
                "SNP_ID": variant_id,
                "Symbol": sub["Gene_Symbol"].iloc[0],
                "Consequence": sub[consequence_col].iloc[0],
                "Impact": sub[impact_col].iloc[0],
                "gnomAD_NFE_AF": sub["gnomAD_NFE_AF_combined"].iloc[0],
                "gnomAD_NFE_Source": sub["gnomAD_NFE_Source"].iloc[0],
                "Tumour_WT": tumour_counts["WT"],
                "Tumour_Het_ALT": tumour_counts["Het_ALT"],
                "Tumour_Hom_ALT": tumour_counts["Hom_ALT"],
                "Healthy_WT": control_counts["WT"],
                "Healthy_Het_ALT": control_counts["Het_ALT"],
                "Healthy_Hom_ALT": control_counts["Hom_ALT"],
                "Tumour_Callable_Samples": tumour_callable_n,
                "Healthy_Callable_Samples": control_callable_n,
                "Tumour_WT_Pct": _pct(tumour_counts["WT"], tumour_callable_n),
                "Tumour_Het_Pct": _pct(tumour_counts["Het_ALT"], tumour_callable_n),
                "Tumour_Hom_Pct": _pct(tumour_counts["Hom_ALT"], tumour_callable_n),
                "Healthy_WT_Pct": _pct(control_counts["WT"], control_callable_n),
                "Healthy_Het_Pct": _pct(control_counts["Het_ALT"], control_callable_n),
                "Healthy_Hom_Pct": _pct(control_counts["Hom_ALT"], control_callable_n),
                "Tumour_Mean_ALT_Dose": tumour_mean_dose,
                "Healthy_Mean_ALT_Dose": control_mean_dose,
                "Dose_Difference": dose_diff,
                "Comparison_Definition": _comparison_definition(analysis_group),
                "Genotype_Definition": "WT = 0 ALT alleles; Het_ALT = 1; Hom_ALT = 2",
                "Association_Scheme": "Genotype distribution plus genotype-dose enrichment for common SNPs",
                "Primary_Genotype_P_Value": primary_p,
                "Primary_Genotype_Method": primary_method,
                "Primary_Genotype_Note": primary_note,
                **distribution,
                **trend,
                **het_contrast,
                **hom_contrast,
            })

    genotype_df = pd.DataFrame(rows)
    if genotype_df.empty:
        note = pd.DataFrame({"Note": ["No common-SNP genotype enrichment results were generated."]})
        return note, note

    for p_col, fdr_col in [
        ("Primary_Genotype_P_Value", "Primary_Genotype_FDR_P_Value"),
        ("Genotype_Distribution_P_Value", "Genotype_Distribution_FDR_P_Value"),
        ("Genotype_Trend_P_Value", "Genotype_Trend_FDR_P_Value"),
        ("Het_vs_WT_P_Value", "Het_vs_WT_FDR_P_Value"),
        ("Hom_vs_WT_P_Value", "Hom_vs_WT_FDR_P_Value"),
    ]:
        genotype_df = _apply_grouped_fdr(genotype_df, p_col, fdr_col, ["Analysis_Group"])

    genotype_df["Raw_P_Threshold"] = RAW_P_THRESHOLD
    genotype_df["FDR_Threshold"] = FDR_THRESHOLD
    genotype_df["Primary_Genotype_Significant"] = pd.to_numeric(genotype_df["Primary_Genotype_P_Value"], errors="coerce") < RAW_P_THRESHOLD
    genotype_df["Primary_Genotype_FDR_Significant"] = pd.to_numeric(genotype_df["Primary_Genotype_FDR_P_Value"], errors="coerce") < FDR_THRESHOLD
    genotype_df["FDR_Scope"] = "BH within analysis group; genotype distribution, trend, and contrast families are corrected separately, plus a primary genotype ranking correction."
    genotype_df = genotype_df.merge(display_lookup, left_on="SNP_ID", right_on="Variant_ID", how="left")
    genotype_df = genotype_df.drop(columns=["Variant_ID"], errors="ignore")
    genotype_df["Display_Label"] = genotype_df["Display_Label"].fillna(genotype_df["SNP_ID"])
    genotype_df["Genomic_Label"] = genotype_df["Genomic_Label"].fillna(genotype_df["SNP_ID"])
    genotype_df["rsID"] = genotype_df["rsID"].replace({"": pd.NA})
    genotype_df = attach_amplicon_warning_columns(
        genotype_df,
        warning_lookup.rename(columns={"Variant_Key": "SNP_ID"}),
        variant_col="SNP_ID",
    )
    genotype_df["Best_Genotype_P"] = genotype_df.apply(
        lambda row: _best_available_p(
            row,
            [
                "Primary_Genotype_P_Value",
                "Genotype_Distribution_P_Value",
                "Genotype_Trend_P_Value",
                "Het_vs_WT_P_Value",
                "Hom_vs_WT_P_Value",
            ],
        ),
        axis=1,
    )
    genotype_df["Best_Genotype_FDR"] = genotype_df.apply(
        lambda row: _best_available_p(
            row,
            [
                "Primary_Genotype_FDR_P_Value",
                "Genotype_Distribution_FDR_P_Value",
                "Genotype_Trend_FDR_P_Value",
                "Het_vs_WT_FDR_P_Value",
                "Hom_vs_WT_FDR_P_Value",
            ],
        ),
        axis=1,
    )

    count_cols = [
        "Analysis_Group", "Display_Label", "rsID", "SNP_ID", "Genomic_Label", "HGVSc", "Symbol", "Consequence", "Impact",
        "Tumour_WT", "Tumour_Het_ALT", "Tumour_Hom_ALT", "Healthy_WT", "Healthy_Het_ALT", "Healthy_Hom_ALT",
        "Tumour_Callable_Samples", "Healthy_Callable_Samples",
        "Tumour_WT_Pct", "Tumour_Het_Pct", "Tumour_Hom_Pct",
        "Healthy_WT_Pct", "Healthy_Het_Pct", "Healthy_Hom_Pct",
        "Tumour_Mean_ALT_Dose", "Healthy_Mean_ALT_Dose", "Dose_Difference",
        "Primary_Genotype_Method", "Primary_Genotype_P_Value", "Primary_Genotype_FDR_P_Value",
        "Genotype_Distribution_Method", "Genotype_Distribution_P_Value", "Genotype_Distribution_FDR_P_Value",
        "Genotype_Trend_Method", "Genotype_Trend_P_Value", "Genotype_Trend_FDR_P_Value",
        "Het_vs_WT_P_Value", "Het_vs_WT_FDR_P_Value", "Hom_vs_WT_P_Value", "Hom_vs_WT_FDR_P_Value",
        "Primary_Genotype_Note", "Coverage_Risk_Flag", "Coverage_Risk_Note",
    ]
    genotype_counts = genotype_df[[col for col in count_cols if col in genotype_df.columns]].copy()
    genotype_df = genotype_df.sort_values(["Analysis_Group", "Best_Genotype_P", "Display_Label"], na_position="last", kind="stable").reset_index(drop=True)
    genotype_counts = genotype_counts.sort_values(["Analysis_Group", "Primary_Genotype_P_Value", "Display_Label"], na_position="last", kind="stable").reset_index(drop=True)
    return genotype_df, genotype_counts


def plot_snp_genotype_heatmap(genotype_df: pd.DataFrame, output_path: str) -> None:
    group_order = ["Breast", "Endometrium", "Global"]
    metric_map = [
        ("Genotype_Distribution_P_Value", "Distribution"),
        ("Genotype_Trend_P_Value", "Dose trend"),
        ("Het_vs_WT_P_Value", "Het vs WT"),
        ("Hom_vs_WT_P_Value", "Hom vs WT"),
    ]
    fig, axes = plt.subplots(1, len(group_order), figsize=(17.5, 5.8), squeeze=False)
    any_data = False

    for ax, group in zip(axes.flat, group_order):
        sub = genotype_df[genotype_df["Analysis_Group"] == group].copy()
        if sub.empty:
            ax.axis("off")
            ax.text(0.5, 0.5, "No genotype results", ha="center", va="center", fontsize=10)
            continue
        sub["Rank_P"] = pd.to_numeric(sub["Best_Genotype_P"], errors="coerce")
        sub = sub.sort_values(["Rank_P", "Display_Label"], na_position="last").head(8).copy()
        if sub.empty:
            ax.axis("off")
            ax.text(0.5, 0.5, "No genotype results", ha="center", va="center", fontsize=10)
            continue
        any_data = True
        rows = []
        for _, row in sub.iterrows():
            for col, label in metric_map:
                p_val = pd.to_numeric(row.get(col), errors="coerce")
                rows.append({
                    "Display_Label": row["Display_Label"],
                    "Metric": label,
                    "Score": -np.log10(max(float(p_val), 1e-300)) if pd.notna(p_val) and p_val > 0 else np.nan,
                })
        heatmap_df = pd.DataFrame(rows).pivot(index="Display_Label", columns="Metric", values="Score")
        heatmap_df = heatmap_df.reindex(columns=[label for _, label in metric_map if label in heatmap_df.columns])
        sns.heatmap(
            heatmap_df,
            ax=ax,
            cmap=sns.light_palette(arm_color("Tumour"), as_cmap=True),
            linewidths=0.5,
            linecolor="#ECECEC",
            annot=True,
            fmt=".1f",
            annot_kws={"fontsize": 7},
            cbar=ax is axes.flat[-1],
            cbar_kws={"label": "-log10 p"} if ax is axes.flat[-1] else None,
        )
        ax.set_title(f"{group} genotype enrichment", fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=18, labelsize=8.5)
        ax.tick_params(axis="y", rotation=0, labelsize=8.5)

    fig.suptitle(tagged_title("Common-SNP genotype enrichment overview", COMPARATIVE_TAG), fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    if any_data:
        _save_publication_figure(fig, output_path)
    else:
        plt.close(fig)


def plot_snp_genotype_composition(common_df: pd.DataFrame, genotype_df: pd.DataFrame, output_prefix: str) -> list[str]:
    manifest = common_df[["Sample", "Cohort", "Tissue"]].drop_duplicates().copy()
    output_paths: list[str] = []
    for group in ["Breast", "Endometrium", "Global"]:
        sub = genotype_df[genotype_df["Analysis_Group"] == group].copy()
        sub["Rank_P"] = pd.to_numeric(sub["Best_Genotype_P"], errors="coerce")
        top = sub.sort_values(["Rank_P", "Display_Label"], na_position="last").head(6)
        if top.empty:
            continue

        if group == "Breast":
            tumour_samples = set(
                manifest.loc[
                    manifest["Tissue"].astype(str).eq("Tumour")
                    & manifest["Cohort"].astype(str).str.contains("Breast", case=False, na=False),
                    "Sample",
                ].astype(str)
            )
        elif group == "Endometrium":
            tumour_samples = set(
                manifest.loc[
                    manifest["Tissue"].astype(str).eq("Tumour")
                    & manifest["Cohort"].astype(str).str.contains("Endometri", case=False, na=False),
                    "Sample",
                ].astype(str)
            )
        else:
            tumour_samples = set(manifest.loc[manifest["Tissue"].astype(str).eq("Tumour"), "Sample"].astype(str))
        control_samples = set(manifest.loc[manifest["Tissue"].astype(str).eq("Healthy"), "Sample"].astype(str))

        fig, axes = plt.subplots(2, 3, figsize=(15.2, 8.8), squeeze=False)
        for ax, (_, row) in zip(axes.flat, top.iterrows()):
            var_df = common_df[common_df["Variant_Key"] == row["SNP_ID"]].copy()
            var_df = var_df[var_df["Callable_Allele_Genotype"]].copy()
            if var_df.empty:
                ax.axis("off")
                continue
            frame = pd.concat(
                [
                    var_df[var_df["Sample"].astype(str).isin(control_samples)][["Sample", "Variant_Genotype_Status"]].assign(Group="Pooled healthy"),
                    var_df[var_df["Sample"].astype(str).isin(tumour_samples)][["Sample", "Variant_Genotype_Status"]].assign(Group="Tumour"),
                ],
                ignore_index=True,
            )
            counts = (
                frame.groupby(["Group", "Variant_Genotype_Status"]).size().unstack(fill_value=0)
                .reindex(index=["Pooled healthy", "Tumour"], columns=GENOTYPE_ORDER, fill_value=0)
            )
            if counts.sum().sum() == 0:
                ax.axis("off")
                continue
            perc = counts.div(counts.sum(axis=1), axis=0).fillna(0.0) * 100
            bottoms = np.zeros(len(perc.index))
            for genotype in GENOTYPE_ORDER:
                label = GENOTYPE_LABELS[genotype]
                ax.bar(
                    np.arange(len(perc.index)),
                    perc[genotype].to_numpy(dtype=float),
                    bottom=bottoms,
                    width=0.62,
                    color=GENOTYPE_COLORS[label],
                    edgecolor="#4B4B4B",
                    linewidth=0.6,
                    label=label,
                )
                bottoms += perc[genotype].to_numpy(dtype=float)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(
                [
                    f"Pooled healthy\n(n={int(counts.sum(axis=1).iloc[0])})",
                    f"Tumour\n(n={int(counts.sum(axis=1).iloc[1])})",
                ],
                fontsize=8,
            )
            ax.set_ylim(0, 100)
            ax.set_ylabel("Genotype composition (%)", fontsize=8.5)
            ax.set_title(_short_label(row["Display_Label"], max_len=26), fontsize=9.5, fontweight="bold")
            ax.grid(True, axis="y", linestyle=":", alpha=0.35)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        for ax in axes.flat[len(top):]:
            ax.axis("off")
        handles, labels = axes.flat[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.01))
        fig.suptitle(f"{group} common-SNP genotype composition", fontsize=13, fontweight="bold", y=1.04)
        fig.tight_layout()
        out_path = f"{output_prefix}_{group}.png"
        _save_publication_figure(fig, out_path)
        output_paths.append(out_path)
    return output_paths


def _mutation_workflow_variant_filter(analysis_df: pd.DataFrame) -> pd.Series:
    return ~common_nfe_variant_mask(analysis_df, THRESHOLDS["min_nfe_af"])


def _mutation_burden_definitions(mutation_df: pd.DataFrame) -> list[dict[str, object]]:
    functional = mutation_df.get("Functional_Class", pd.Series("", index=mutation_df.index)).astype(str).str.lower()
    location = mutation_df.get("Location_Class", pd.Series("", index=mutation_df.index)).astype(str)
    impact = mutation_df.get("Impact_Severity", pd.Series("", index=mutation_df.index)).astype(str).str.upper()
    consequence = mutation_df.get("Detailed_Consequence", pd.Series("", index=mutation_df.index)).astype(str).str.lower()
    gene = mutation_df.get("Gene_Symbol", pd.Series("", index=mutation_df.index)).astype(str).str.upper()
    known_rare_mask = pd.to_numeric(mutation_df.get("gnomAD_NFE_AF_combined"), errors="coerce").lt(THRESHOLDS["min_nfe_af"]).fillna(False)

    definitions = [
        {
            "Burden_ID": "Any_Mutation",
            "Burden_Label": "Any mutation",
            "Mask": pd.Series(True, index=mutation_df.index),
            "Definition": "Any non-common variant in the mutation workflow (gnomAD NFE AF <= 1% or AF missing).",
        },
        {
            "Burden_ID": "Any_Rare_Mutation",
            "Burden_Label": "Any rare mutation",
            "Mask": known_rare_mask,
            "Definition": "Any mutation with observed combined gnomAD NFE AF below the common-variant threshold.",
        },
        {
            "Burden_ID": "HIGH_Impact_Mutation",
            "Burden_Label": "Any HIGH impact mutation",
            "Mask": impact.eq("HIGH"),
            "Definition": "Any mutation annotated as HIGH impact.",
        },
        {
            "Burden_ID": "MODERATE_or_HIGH_Impact_Mutation",
            "Burden_Label": "Any MODERATE/HIGH mutation",
            "Mask": impact.isin(["MODERATE", "HIGH"]),
            "Definition": "Any mutation annotated as MODERATE or HIGH impact.",
        },
        {
            "Burden_ID": "Coding_Mutation_Burden",
            "Burden_Label": "Coding mutation burden",
            "Mask": location.eq("Exonic/Coding") | functional.isin(["synonymous", "missense", "nonsense", "frameshift"]),
            "Definition": "Any mutation in coding / exonic consequence classes, including synonymous, missense, nonsense, and frameshift changes.",
        },
        {
            "Burden_ID": "Splice_Related_Mutation_Burden",
            "Burden_Label": "Splice-related mutation burden",
            "Mask": location.eq("Splicing") | functional.eq("splice_site") | consequence.str.contains("splice", na=False),
            "Definition": "Any splice-site, splice-region, or other splice-related mutation.",
        },
    ]

    gsdmb_mask = gene.eq("GSDMB")
    if gsdmb_mask.any() and not gsdmb_mask.all():
        definitions.append(
            {
                "Burden_ID": "GSDMB_Only_Mutation_Burden",
                "Burden_Label": "GSDMB-only mutation burden",
                "Mask": gsdmb_mask,
                "Definition": "Any mutation whose representative annotation maps to GSDMB specifically.",
            }
        )
    return definitions


def compute_mutation_enrichment_outputs(analysis_df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    mutation_df = analysis_df[_mutation_workflow_variant_filter(analysis_df)].copy()
    if mutation_df.empty:
        note = pd.DataFrame({"Note": ["No mutation-level enrichment results were generated because no non-common variants passed the input filters."]})
        return {
            "mutation_carrier_enrichment": note,
            "mutation_burden": note,
            "mutation_hotspots": note,
            "mutation_summary_significant": note,
        }

    display_lookup = build_variant_display_table(mutation_df, "Variant_Key", "rsID", "Gene_Symbol", "HGVSc", "CHROM", "POS", "REF", "ALT")
    warning_lookup = build_amplicon_warning_lookup(mutation_df, variant_col="Variant_Key", chrom_col="CHROM", pos_col="POS")
    manifest = mutation_df[["Sample", "Cohort", "Tissue"]].drop_duplicates().copy()
    manifest["Sample"] = manifest["Sample"].astype(str)
    cohort_order = ["Breast", "Endometrium", "Global"]
    cohort_filters = {"Breast": "Breast", "Endometrium": "Endometri"}

    carrier_rows: list[dict[str, object]] = []
    hotspot_rows: list[dict[str, object]] = []

    mutation_df["Position_ID"] = (
        mutation_df["CHROM"].astype(str).str.strip()
        + ":"
        + pd.to_numeric(mutation_df["POS"], errors="coerce").fillna(-1).astype(int).astype(str)
    )

    for analysis_group in cohort_order:
        if analysis_group == "Global":
            cohort_df = pooled_control_frame(mutation_df, "Cohort", "Tissue")
        else:
            cohort_df = pooled_control_frame(mutation_df, "Cohort", "Tissue", analysis_group, cohort_filters[analysis_group])
        if cohort_df.empty:
            continue

        for variant_id, sub in cohort_df.groupby("Variant_Key", sort=False):
            tumour = sub[sub["Tissue"] == "Tumour"].copy()
            control = sub[sub["Tissue"] == "Healthy"].copy()
            tumour_callable = tumour[tumour["Callable_Allele_Genotype"]].copy()
            control_callable = control[control["Callable_Allele_Genotype"]].copy()
            tumour_callable_n = int(tumour_callable["Sample"].nunique())
            control_callable_n = int(control_callable["Sample"].nunique())
            if tumour_callable_n == 0 or control_callable_n == 0:
                continue
            tumour_counts = _genotype_count_map(tumour_callable["Variant_Genotype_Status"])
            control_counts = _genotype_count_map(control_callable["Variant_Genotype_Status"])
            tumour_carriers = int(tumour_counts["Het_ALT"] + tumour_counts["Hom_ALT"])
            control_carriers = int(control_counts["Het_ALT"] + control_counts["Hom_ALT"])
            tumour_noncarriers = int(max(0, tumour_callable_n - tumour_carriers))
            control_noncarriers = int(max(0, control_callable_n - control_carriers))
            _, p_value = fisher_exact([[tumour_carriers, tumour_noncarriers], [control_carriers, control_noncarriers]], alternative="two-sided")
            odds_ratio, or_method = _haldane_odds_ratio(tumour_carriers, tumour_noncarriers, control_carriers, control_noncarriers)
            genotype_eligible = int(tumour_carriers + control_carriers) >= MUTATION_GENOTYPE_MIN_ALT_COUNT
            genotype_test = _genotype_distribution_test(tumour_counts, control_counts) if genotype_eligible else {}
            genotype_note = (
                genotype_test.get("Genotype_Distribution_Note", "")
                if genotype_eligible
                else f"Skipped genotype-level mutation analysis because only {tumour_carriers + control_carriers} ALT genotype observations were available across tumour and healthy groups."
            )
            row = {
                "Analysis_Group": analysis_group,
                "Variant_ID": variant_id,
                "Tumour_Callable_Samples": tumour_callable_n,
                "Healthy_Callable_Samples": control_callable_n,
                "Tumour_Carriers": tumour_carriers,
                "Healthy_Carriers": control_carriers,
                "Tumour_NonCarriers": tumour_noncarriers,
                "Healthy_NonCarriers": control_noncarriers,
                "Tumour_Carrier_Pct": _pct(tumour_carriers, tumour_callable_n),
                "Healthy_Carrier_Pct": _pct(control_carriers, control_callable_n),
                "Carrier_Frequency_Difference_Pct": None if tumour_callable_n == 0 or control_callable_n == 0 else _pct(tumour_carriers, tumour_callable_n) - _pct(control_carriers, control_callable_n),
                "Odds_Ratio": odds_ratio,
                "OR_Method": or_method,
                "P_Value": float(p_value),
                "Association_Scheme": "Per-mutation ALT carrier enrichment (tumour versus pooled healthy controls)",
                "Comparison_Definition": _comparison_definition(analysis_group),
                "Variant_Definition": "Mutation workflow = non-common variant (gnomAD NFE AF <= 1% or AF missing).",
                "Tumour_WT": tumour_counts["WT"],
                "Tumour_Het_ALT": tumour_counts["Het_ALT"],
                "Tumour_Hom_ALT": tumour_counts["Hom_ALT"],
                "Healthy_WT": control_counts["WT"],
                "Healthy_Het_ALT": control_counts["Het_ALT"],
                "Healthy_Hom_ALT": control_counts["Hom_ALT"],
                "Tumour_Mean_ALT_Dose": _mean_alt_dose(tumour_counts),
                "Healthy_Mean_ALT_Dose": _mean_alt_dose(control_counts),
                "Genotype_Analysis_Eligible": genotype_eligible,
                "Genotype_Analysis_Method": genotype_test.get("Genotype_Distribution_Method", "Skipped"),
                "Genotype_Analysis_P_Value": genotype_test.get("Genotype_Distribution_P_Value", np.nan),
                "Genotype_Analysis_Note": genotype_note,
            }
            carrier_rows.append(row)

        for position_id, sub in cohort_df.groupby("Position_ID", sort=False):
            tumour = sub[sub["Tissue"] == "Tumour"].copy()
            control = sub[sub["Tissue"] == "Healthy"].copy()
            tumour_callable = tumour[tumour["Callable_Allele_Genotype"]].copy()
            control_callable = control[control["Callable_Allele_Genotype"]].copy()
            tumour_callable_n = int(tumour_callable["Sample"].nunique())
            control_callable_n = int(control_callable["Sample"].nunique())
            if tumour_callable_n == 0 or control_callable_n == 0:
                continue
            tumour_carriers = int(tumour_callable.loc[tumour_callable["Variant_Is_Carrier"], "Sample"].nunique())
            control_carriers = int(control_callable.loc[control_callable["Variant_Is_Carrier"], "Sample"].nunique())
            if tumour_carriers + control_carriers == 0:
                continue
            tumour_noncarriers = int(max(0, tumour_callable_n - tumour_carriers))
            control_noncarriers = int(max(0, control_callable_n - control_carriers))
            _, p_value = fisher_exact([[tumour_carriers, tumour_noncarriers], [control_carriers, control_noncarriers]], alternative="two-sided")
            odds_ratio, or_method = _haldane_odds_ratio(tumour_carriers, tumour_noncarriers, control_carriers, control_noncarriers)
            hotspot_rows.append({
                "Analysis_Group": analysis_group,
                "Position_ID": position_id,
                "Gene_Symbol": _first_nonempty(sub["Gene_Symbol"]) if "Gene_Symbol" in sub.columns else "",
                "CHROM": _first_nonempty(sub["CHROM"]),
                "POS": pd.to_numeric(sub["POS"], errors="coerce").dropna().astype(int).iloc[0] if pd.to_numeric(sub["POS"], errors="coerce").notna().any() else np.nan,
                "Observed_Alleles": ", ".join(sorted(sub["ALT"].dropna().astype(str).str.upper().unique())),
                "Representative_HGVSc": _short_hgvsc(_first_nonempty(sub["HGVSc"])) if "HGVSc" in sub.columns else "",
                "Position_Consequence": _first_nonempty(sub["Detailed_Consequence"]) if "Detailed_Consequence" in sub.columns else "",
                "Highest_Impact": _first_nonempty(sub["Impact_Severity"]) if "Impact_Severity" in sub.columns else "",
                "Tumour_Callable_Samples": tumour_callable_n,
                "Healthy_Callable_Samples": control_callable_n,
                "Tumour_Carriers": tumour_carriers,
                "Healthy_Carriers": control_carriers,
                "Tumour_Carrier_Pct": _pct(tumour_carriers, tumour_callable_n),
                "Healthy_Carrier_Pct": _pct(control_carriers, control_callable_n),
                "Recurrence_Count": int(tumour_carriers + control_carriers),
                "Recurrence_Rate_Pct": _pct(int(tumour_carriers + control_carriers), tumour_callable_n + control_callable_n),
                "Odds_Ratio": odds_ratio,
                "OR_Method": or_method,
                "P_Value": float(p_value),
                "Hotspot_Definition": "Position-level recurrence across all ALT alleles observed at the same genomic coordinate.",
            })

    mutation_carrier = pd.DataFrame(carrier_rows)
    mutation_hotspots = pd.DataFrame(hotspot_rows)

    if mutation_carrier.empty:
        note = pd.DataFrame({"Note": ["No mutation-level carrier enrichment rows were generated after callable-sample filtering."]})
        return {
            "mutation_carrier_enrichment": note,
            "mutation_burden": note,
            "mutation_hotspots": note,
            "mutation_summary_significant": note,
        }

    mutation_carrier = _apply_grouped_fdr(mutation_carrier, "P_Value", "FDR_P_Value", ["Analysis_Group"])
    mutation_carrier = _apply_grouped_fdr(mutation_carrier, "Genotype_Analysis_P_Value", "Genotype_Analysis_FDR_P_Value", ["Analysis_Group"])
    mutation_carrier["Nominal_Significant"] = mutation_carrier["P_Value"] < RAW_P_THRESHOLD
    mutation_carrier["FDR_Significant"] = mutation_carrier["FDR_P_Value"] < FDR_THRESHOLD
    mutation_carrier["Genotype_Analysis_FDR_Significant"] = pd.to_numeric(mutation_carrier["Genotype_Analysis_FDR_P_Value"], errors="coerce") < FDR_THRESHOLD
    mutation_carrier["FDR_Scope"] = "BH within analysis group for per-mutation carrier enrichment; optional mutation-genotype p-values corrected separately."
    mutation_carrier = mutation_carrier.merge(display_lookup, left_on="Variant_ID", right_on="Variant_ID", how="left")
    mutation_carrier["Display_Label"] = mutation_carrier["Display_Label"].fillna(mutation_carrier["Variant_ID"])
    mutation_carrier["Genomic_Label"] = mutation_carrier["Genomic_Label"].fillna(mutation_carrier["Variant_ID"])
    mutation_carrier = attach_amplicon_warning_columns(
        mutation_carrier,
        warning_lookup.rename(columns={"Variant_Key": "Variant_ID"}),
        variant_col="Variant_ID",
    )

    burden_rows: list[dict[str, object]] = []
    carrier_samples = mutation_df.loc[mutation_df["Variant_Is_Carrier"], ["Sample", "Cohort", "Tissue"]].drop_duplicates().copy()
    carrier_samples["Sample"] = carrier_samples["Sample"].astype(str)
    burden_defs = _mutation_burden_definitions(mutation_df)
    for burden_def in burden_defs:
        burden_carrier_set = set(
            mutation_df.loc[burden_def["Mask"] & mutation_df["Variant_Is_Carrier"], "Sample"].astype(str)
        )
        if not burden_carrier_set:
            continue
        for analysis_group in cohort_order:
            tumour_samples, control_samples = _variant_group_sets(manifest, analysis_group)
            if len(tumour_samples) == 0 or len(control_samples) == 0:
                continue
            tumour_carriers = len(burden_carrier_set & tumour_samples)
            control_carriers = len(burden_carrier_set & control_samples)
            tumour_noncarriers = len(tumour_samples) - tumour_carriers
            control_noncarriers = len(control_samples) - control_carriers
            _, p_value = fisher_exact([[tumour_carriers, tumour_noncarriers], [control_carriers, control_noncarriers]], alternative="two-sided")
            odds_ratio, or_method = _haldane_odds_ratio(tumour_carriers, tumour_noncarriers, control_carriers, control_noncarriers)
            burden_rows.append({
                "Analysis_Group": analysis_group,
                "Burden_ID": burden_def["Burden_ID"],
                "Burden_Label": burden_def["Burden_Label"],
                "Definition": burden_def["Definition"],
                "Tumour_Samples": len(tumour_samples),
                "Healthy_Samples": len(control_samples),
                "Tumour_Carriers": tumour_carriers,
                "Healthy_Carriers": control_carriers,
                "Tumour_Carrier_Pct": _pct(tumour_carriers, len(tumour_samples)),
                "Healthy_Carrier_Pct": _pct(control_carriers, len(control_samples)),
                "Carrier_Frequency_Difference_Pct": _pct(tumour_carriers, len(tumour_samples)) - _pct(control_carriers, len(control_samples)),
                "Odds_Ratio": odds_ratio,
                "OR_Method": or_method,
                "P_Value": float(p_value),
                "Association_Scheme": "Mutation burden enrichment (carrier-based burden flag at the sample level)",
            })

    mutation_burden = pd.DataFrame(burden_rows)
    if not mutation_burden.empty:
        mutation_burden = _apply_grouped_fdr(mutation_burden, "P_Value", "FDR_P_Value", ["Analysis_Group"])
        mutation_burden["Nominal_Significant"] = mutation_burden["P_Value"] < RAW_P_THRESHOLD
        mutation_burden["FDR_Significant"] = mutation_burden["FDR_P_Value"] < FDR_THRESHOLD
        mutation_burden["FDR_Scope"] = "BH within analysis group for burden definitions."

    if not mutation_hotspots.empty:
        mutation_hotspots = _apply_grouped_fdr(mutation_hotspots, "P_Value", "FDR_P_Value", ["Analysis_Group"])
        mutation_hotspots["Nominal_Significant"] = mutation_hotspots["P_Value"] < RAW_P_THRESHOLD
        mutation_hotspots["FDR_Significant"] = mutation_hotspots["FDR_P_Value"] < FDR_THRESHOLD
        mutation_hotspots["FDR_Scope"] = "BH within analysis group for recurrent mutation-position hotspot summaries."
        mutation_hotspots = mutation_hotspots.sort_values(
            ["Analysis_Group", "Recurrence_Count", "P_Value", "Position_ID"],
            ascending=[True, False, True, True],
            kind="stable",
        ).reset_index(drop=True)

    summary_frames = []
    for workflow_label, df in [
        ("mutation_carrier_enrichment", mutation_carrier),
        ("mutation_burden", mutation_burden),
        ("mutation_hotspots", mutation_hotspots),
    ]:
        if df is None or df.empty:
            continue
        sig = df[df.get("Nominal_Significant", pd.Series(False, index=df.index)).fillna(False) | df.get("FDR_Significant", pd.Series(False, index=df.index)).fillna(False)].copy()
        if sig.empty:
            continue
        sig["Workflow"] = workflow_label
        summary_frames.append(sig)
    mutation_summary = (
        pd.concat(summary_frames, ignore_index=True).sort_values(["Workflow", "Analysis_Group", "P_Value"], kind="stable")
        if summary_frames else
        pd.DataFrame({"Note": ["No nominally significant mutation-enrichment, burden, or hotspot signals were detected."]})
    )

    return {
        "mutation_carrier_enrichment": mutation_carrier.sort_values(["Analysis_Group", "P_Value", "Display_Label"], kind="stable").reset_index(drop=True),
        "mutation_burden": mutation_burden.sort_values(["Analysis_Group", "P_Value", "Burden_ID"], kind="stable").reset_index(drop=True) if not mutation_burden.empty else pd.DataFrame(),
        "mutation_hotspots": mutation_hotspots,
        "mutation_summary_significant": mutation_summary,
    }


def plot_mutation_burden_forest(mutation_burden: pd.DataFrame, output_path: str) -> None:
    if mutation_burden is None or mutation_burden.empty:
        return
    group_order = ["Breast", "Endometrium", "Global"]
    burden_order = list(dict.fromkeys(mutation_burden["Burden_Label"].astype(str)))
    fig, axes = plt.subplots(1, len(group_order), figsize=(16.5, max(4.6, 0.55 * len(burden_order) + 2.0)), sharey=True, squeeze=False)
    for ax, group in zip(axes.flat, group_order):
        sub = mutation_burden[mutation_burden["Analysis_Group"] == group].copy()
        if sub.empty:
            ax.axis("off")
            ax.text(0.5, 0.5, "No burden results", ha="center", va="center")
            continue
        sub["Burden_Label"] = pd.Categorical(sub["Burden_Label"], categories=burden_order, ordered=True)
        sub = sub.sort_values("Burden_Label")
        x_vals = np.log2(pd.to_numeric(sub["Odds_Ratio"], errors="coerce").replace(0, np.nan))
        y_vals = np.arange(len(sub))
        point_sizes = np.where(sub["FDR_Significant"].fillna(False), 130, np.where(sub["Nominal_Significant"].fillna(False), 95, 60))
        ax.scatter(x_vals, y_vals, s=point_sizes, color=cohort_color(group), edgecolor="#2F2F2F", linewidth=0.7, alpha=0.88)
        ax.axvline(0, color="#6E6E6E", linestyle="--", linewidth=1.0)
        ax.set_title(group, fontweight="bold")
        ax.set_xlabel("log2(odds ratio)")
        ax.set_yticks(y_vals)
        ax.set_yticklabels(sub["Burden_Label"].astype(str), fontsize=8.5)
        ax.grid(True, axis="x", linestyle=":", alpha=0.35)
        for x_val, y_val, p_val in zip(x_vals, y_vals, sub["P_Value"], strict=False):
            if pd.notna(x_val) and pd.notna(p_val):
                ax.text(x_val, y_val + 0.13, f"p={p_val:.3g}", fontsize=7, ha="center", va="bottom")
    fig.suptitle("Mutation burden enrichment across tumour and pooled healthy groups", fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    _save_publication_figure(fig, output_path)


def plot_mutation_hotspots(mutation_hotspots: pd.DataFrame, output_path: str) -> None:
    if mutation_hotspots is None or mutation_hotspots.empty:
        return
    sub = mutation_hotspots[mutation_hotspots["Analysis_Group"] == "Global"].copy()
    if sub.empty:
        return
    sub = sub.sort_values(["Recurrence_Count", "P_Value"], ascending=[False, True], kind="stable").head(10)
    if sub.empty:
        return
    labels = [
        _short_label(f"{row['Gene_Symbol']} {row['Position_ID']}".strip(), max_len=30)
        for _, row in sub.iterrows()
    ]
    y = np.arange(len(sub))
    fig, ax = plt.subplots(figsize=(10.8, max(4.8, 0.55 * len(sub) + 1.8)))
    ax.barh(y + 0.18, sub["Tumour_Carrier_Pct"], height=0.34, color=arm_color("Tumour"), label="Tumour")
    ax.barh(y - 0.18, sub["Healthy_Carrier_Pct"], height=0.34, color=arm_color("Control"), label="Pooled healthy")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("Carrier frequency (%)")
    ax.set_title("Recurrent mutation-position hotspots in tumours and pooled healthy controls", fontweight="bold")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(True, axis="x", linestyle=":", alpha=0.35)
    for idx, (_, row) in enumerate(sub.iterrows()):
        ax.text(max(row["Tumour_Carrier_Pct"], row["Healthy_Carrier_Pct"]) + 0.6, idx, f"n={int(row['Recurrence_Count'])} | p={row['P_Value']:.3g}", fontsize=7.4, va="center")
    fig.tight_layout()
    _save_publication_figure(fig, output_path)


def _annotate_ranked_hits(ax, labelled, label_col, x_col, y_col, max_labels=5):
    if labelled.empty:
        return
    ranked = labelled.sort_values([y_col, x_col], ascending=[False, True]).head(max_labels)
    y_min, y_max = ax.get_ylim()
    pad = (y_max - y_min) * 0.08 if y_max != y_min else 0.5
    if len(ranked) == 1:
        y_positions = [min(y_max - pad, ranked.iloc[0][y_col] + pad * 0.25)]
    else:
        upper = max(y_max - pad, y_min + pad)
        lower = max(y_min + pad, y_min + (y_max - y_min) * 0.18)
        y_positions = np.linspace(upper, lower, len(ranked))

    text_transform = blended_transform_factory(ax.transAxes, ax.transData)
    x_frac = 0.96

    for (_, row), y_text in zip(ranked.iterrows(), y_positions):
        ax.annotate(
            _short_label(row[label_col], max_len=22),
            xy=(row[x_col], row[y_col]),
            xycoords="data",
            xytext=(x_frac, y_text),
            textcoords=text_transform,
            ha="right",
            va="center",
            fontsize=8.8,
            fontweight="bold",
            clip_on=False,
            bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="#9e9e9e", alpha=0.96, linewidth=0.8),
            arrowprops=dict(arrowstyle="-", color="#777777", lw=0.8, shrinkA=0, shrinkB=4),
            zorder=5,
        )



def _singleton_block_id(snp_id):
    return f"Singleton::{snp_id}"


def _initialise_ld_block_columns(res_df, ld_reference):
    block_df = ld_reference["block_membership"].copy()
    if not block_df.empty:
        block_df["LD_Block_ID"] = block_df["Block"].astype(str)
        block_map = dict(zip(block_df["rsID"], block_df["LD_Block_ID"]))
        block_sizes = block_df.groupby("LD_Block_ID")["rsID"].size().to_dict()
        block_members = (
            block_df.groupby("LD_Block_ID")["rsID"]
            .apply(lambda series: ", ".join(series.astype(str)))
            .to_dict()
        )
    else:
        block_map = {}
        block_sizes = {}
        block_members = {}

    out = res_df.copy()
    out["LD_Threshold"] = ld_reference["threshold"]
    out["Threshold_Label"] = ld_reference["threshold_label"]
    out["LD_Block_ID"] = out["SNP_ID"].map(block_map)
    out["LD_Block_Assignment"] = np.where(out["LD_Block_ID"].notna(), "LD_Block", "Singleton")
    out.loc[out["LD_Block_ID"].isna(), "LD_Block_ID"] = out.loc[out["LD_Block_ID"].isna(), "SNP_ID"].map(_singleton_block_id)
    out["LD_Block_Size"] = out["LD_Block_ID"].map(block_sizes).fillna(1).astype(int)
    out["LD_Block_Members"] = out["LD_Block_ID"].map(block_members).fillna(out["SNP_ID"])
    return out, block_df


def _apply_block_fdr(res_df, group_col, p_col, block_col, primary_fdr_col, legacy_fdr_col):
    final_parts = []
    for group_value, group_df in res_df.groupby(group_col, sort=False):
        group_df = group_df.copy()
        _, group_df[legacy_fdr_col], _, _ = multipletests(group_df[p_col], method="fdr_bh")
        block_parts = []
        for _, block_df in group_df.groupby(block_col, sort=False):
            block_df = block_df.copy()
            if len(block_df) == 1:
                block_df[primary_fdr_col] = block_df[p_col]
            else:
                _, block_df[primary_fdr_col], _, _ = multipletests(block_df[p_col], method="fdr_bh")
            block_df["LD_Block_Test_Count"] = int(len(block_df))
            block_parts.append(block_df)
        group_df = pd.concat(block_parts, ignore_index=True)
        threshold_label = group_df["Threshold_Label"].iloc[0] if "Threshold_Label" in group_df.columns else "LD"
        group_df["FDR_Scope"] = f"Within_{group_value}_{threshold_label}_blocks"
        final_parts.append(group_df)
    return pd.concat(final_parts, ignore_index=True)


def _refresh_block_member_labels(df):
    """Rebuild block-member text from rsIDs so exported summaries stay rsID-centric."""
    out = df.copy()
    if "LD_Block_ID" not in out.columns or "rsID" not in out.columns:
        return out
    member_map = (
        out.groupby("LD_Block_ID")["rsID"]
        .apply(lambda series: ", ".join(dict.fromkeys(series.astype(str))))
        .to_dict()
    )
    out["LD_Block_Members"] = out["LD_Block_ID"].map(member_map).fillna(out["rsID"].astype(str))
    return out


def _build_ld_summary_tables(final_res, ld_reference, block_membership_df):
    global_df = final_res[final_res["Analysis_Group"] == "Global"].copy()
    if global_df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    global_df["Legacy_Global_FDR_Significant"] = global_df["FDR_P_Value_Global_Legacy"] < FDR_THRESHOLD
    global_df["Block_FDR_Significant"] = global_df["FDR_P_Value"] < FDR_THRESHOLD

    block_context = (
        global_df.groupby("LD_Block_ID", as_index=False)
        .agg(
            Block_Test_Count=("LD_Block_Test_Count", "first"),
            Block_Min_P_Value=("P_Value", "min"),
            Block_Min_FDR_P_Value=("FDR_P_Value", "min"),
            Block_Min_Global_Legacy_FDR=("FDR_P_Value_Global_Legacy", "min"),
            Block_FDR_Significant_SNPs=("Block_FDR_Significant", "sum"),
            Legacy_Global_FDR_Significant_SNPs=("Legacy_Global_FDR_Significant", "sum"),
            Block_Tumour_Mean_Dose_Min=("Tumour_Mean_ALT_Dose", "min"),
            Block_Tumour_Mean_Dose_Max=("Tumour_Mean_ALT_Dose", "max"),
            Block_Control_Mean_Dose_Min=("Control_Mean_ALT_Dose", "min"),
            Block_Control_Mean_Dose_Max=("Control_Mean_ALT_Dose", "max"),
            Block_Dose_Difference_Min=("Dose_Difference", "min"),
            Block_Dose_Difference_Max=("Dose_Difference", "max"),
        )
    )

    pairwise = ld_reference["ld_long"].copy()
    sig_snps = sorted(global_df.loc[global_df["Legacy_Global_FDR_Significant"], "SNP_ID"].astype(str).unique())
    if not sig_snps:
        sig_snps = sorted(global_df.loc[global_df["Block_FDR_Significant"], "SNP_ID"].astype(str).unique())
    pairwise_sig = pairwise[
        pairwise["SNP1"].isin(sig_snps)
        & pairwise["SNP2"].isin(sig_snps)
        & (pairwise["SNP1"] < pairwise["SNP2"])
    ].copy()
    if not pairwise_sig.empty:
        block_lookup = dict(zip(global_df["SNP_ID"], global_df["LD_Block_ID"]))
        pairwise_sig["LD_Block_ID_1"] = pairwise_sig["SNP1"].map(block_lookup)
        pairwise_sig["LD_Block_ID_2"] = pairwise_sig["SNP2"].map(block_lookup)
        pairwise_sig["Same_LD_Block"] = pairwise_sig["LD_Block_ID_1"] == pairwise_sig["LD_Block_ID_2"]
        pairwise_sig["Pass_r2_0_80"] = pairwise_sig["R2"] >= LD_R2_THRESHOLD

    pairwise_strings = {}
    if not pairwise_sig.empty:
        for snp in sig_snps:
            rows = pairwise_sig[(pairwise_sig["SNP1"] == snp) | (pairwise_sig["SNP2"] == snp)].copy()
            if rows.empty:
                pairwise_strings[snp] = ""
                continue
            bits = []
            for _, row in rows.iterrows():
                other = row["SNP2"] if row["SNP1"] == snp else row["SNP1"]
                same = "same block" if row["Same_LD_Block"] else "different block"
                bits.append(f"{other} (r2={row['R2']:.3f}, {same})")
            pairwise_strings[snp] = "; ".join(bits)

    global_summary = global_df.merge(block_context, on="LD_Block_ID", how="left")
    global_summary["Pairwise_LD_to_Global_Significant_SNPs"] = global_summary["SNP_ID"].map(pairwise_strings).fillna("")
    ordered_cols = [
        "SNP_ID", "LD_Block_ID", "LD_Block_Assignment", "LD_Block_Size", "LD_Block_Members",
        "P_Value", "FDR_P_Value", "FDR_P_Value_Global_Legacy",
        "Raw_Significant", "FDR_Significant", "Legacy_Global_FDR_Significant",
        "Tumour_WT", "Tumour_Het_ALT", "Tumour_Hom_ALT",
        "Control_WT", "Control_Het_ALT", "Control_Hom_ALT",
        "Tumour_Callable_Samples", "Control_Callable_Samples",
        "Tumour_WT_Pct_Callable", "Tumour_Het_Pct_Callable", "Tumour_Hom_Pct_Callable",
        "Control_WT_Pct_Callable", "Control_Het_Pct_Callable", "Control_Hom_Pct_Callable",
        "Tumour_Mean_ALT_Dose", "Control_Mean_ALT_Dose", "Dose_Difference",
        "Block_Test_Count", "Block_Min_P_Value", "Block_Min_FDR_P_Value", "Block_Min_Global_Legacy_FDR",
        "Block_FDR_Significant_SNPs", "Legacy_Global_FDR_Significant_SNPs",
        "Block_Tumour_Mean_Dose_Min", "Block_Tumour_Mean_Dose_Max",
        "Block_Control_Mean_Dose_Min", "Block_Control_Mean_Dose_Max",
        "Block_Dose_Difference_Min", "Block_Dose_Difference_Max",
        "Pairwise_LD_to_Global_Significant_SNPs",
    ]
    global_summary = global_summary[ordered_cols + [c for c in global_summary.columns if c not in ordered_cols]]
    return global_summary, pairwise_sig, block_membership_df


def run_snp_association_analysis():
    """Run tumour-versus-control SNP genotype association tests and export the results."""
    print(">>> Starting Script 12: SNP genotype association with carrier volcano plots and genotype contrast volcanoes...")

    validate_file_exists(input_file, "Script 12 input")

    # --- 2. DATA LOADING & CLEANING ---
    analysis_df = build_allele_genotype_long(PATHS["annotated_report"], PATHS["forced_genotypes_dir"] / "GSDMB_Forced_Genotypes_Union_Sites.xlsx")
    analysis_df = analysis_df[analysis_df["Tissue"].isin(["Tumour", "Healthy"])].copy()
    long_df = analysis_df[common_nfe_variant_mask(analysis_df, THRESHOLDS["min_nfe_af"])].copy()
    validate_nonempty(long_df, "Script 12 common SNP subset")
    print_validation_summary(long_df, "Sample", "Script 12 filtered SNPs", ["Cohort", "Tissue"])

    warning_lookup = build_amplicon_warning_lookup(long_df, variant_col="Variant_Key", chrom_col="CHROM", pos_col="POS")
    display_lookup = build_variant_display_table(long_df, "Variant_Key", "rsID", "Gene_Symbol", "HGVSc", "CHROM", "POS", "REF", "ALT")

    cohort_list = ["Global", "Breast", "Endometrium"]
    cohort_filters = {"Breast": "Breast", "Endometrium": "Endometri"}
    all_results = []

    # --- 6. GENOTYPE ASSOCIATION TESTS WITH CHI-SQUARE FOR 3x2 TABLES ---
    for cohort_name in cohort_list:
        if cohort_name == "Global":
            cohort_df = pooled_control_frame(long_df, "Cohort", "Tissue")
        else:
            cohort_df = pooled_control_frame(long_df, "Cohort", "Tissue", cohort_name, cohort_filters[cohort_name])

        if cohort_df.empty:
            continue

        unique_snps = cohort_df["Variant_Key"].dropna().astype(str).unique()

        for var in unique_snps:
            var_data = cohort_df[cohort_df["Variant_Key"] == var].copy()
            normal = var_data[var_data["Tissue"] == "Healthy"]
            tumour = var_data[var_data["Tissue"] == "Tumour"]
            normal_callable = normal[normal["Callable_Allele_Genotype"]].copy()
            tumour_callable = tumour[tumour["Callable_Allele_Genotype"]].copy()

            normal_total = int(normal["Sample"].nunique())
            tumour_total = int(tumour["Sample"].nunique())
            normal_callable_n = int(normal_callable["Sample"].nunique())
            tumour_callable_n = int(tumour_callable["Sample"].nunique())

            if normal_callable_n == 0 or tumour_callable_n == 0:
                continue

            normal_counts = _genotype_count_map(normal_callable["Variant_Genotype_Status"])
            tumour_counts = _genotype_count_map(tumour_callable["Variant_Genotype_Status"])

            contingency = np.array([
                [normal_counts["WT"], normal_counts["Het_ALT"], normal_counts["Hom_ALT"]],
                [tumour_counts["WT"], tumour_counts["Het_ALT"], tumour_counts["Hom_ALT"]],
            ], dtype=float)
            contingency_test = contingency[contingency.sum(axis=1) > 0, :]
            calc_note = "Chi-square on a 3x2 WT/Het/Hom table"
            if contingency_test.shape[0] < 2 or contingency_test.shape[1] < 2:
                chi2 = 0.0
                p_value = 1.0
                dof = 0
                expected = contingency_test.astype(float)
                calc_note = "Degenerate genotype table; association not estimable"
            else:
                try:
                    chi2, p_value, dof, expected = chi2_contingency(contingency_test, correction=False)
                except ValueError:
                    contingency_test = contingency_test + 0.5
                    chi2, p_value, dof, expected = chi2_contingency(contingency_test, correction=False)
                    calc_note = "Chi-square on a 3x2 WT/Het/Hom table with 0.5 pseudocount fallback"
            expected_min = float(np.min(expected)) if expected.size else np.nan

            normal_mean_dose = _mean_alt_dose(normal_counts)
            tumour_mean_dose = _mean_alt_dose(tumour_counts)
            dose_diff = None if normal_mean_dose is None or tumour_mean_dose is None else tumour_mean_dose - normal_mean_dose
            if dose_diff is None:
                direction = "Not estimable"
            elif dose_diff < -1e-9:
                direction = "Lower in tumour"
            elif dose_diff > 1e-9:
                direction = "Higher in tumour"
            else:
                direction = "Same dose"

            all_results.append({
                "Analysis_Group": cohort_name,
                "SNP_ID": var,
                "Symbol": var_data["Gene_Symbol"].iloc[0],
                "Consequence": var_data["Detailed_Consequence"].iloc[0] if "Detailed_Consequence" in var_data.columns else var_data["Consequence"].iloc[0],
                "Impact": var_data["Impact_Severity"].iloc[0] if "Impact_Severity" in var_data.columns else var_data["Impact"].iloc[0],
                "gnomAD_NFE_AF": var_data["gnomAD_NFE_AF_combined"].iloc[0],
                "gnomAD_NFE_Source": var_data["gnomAD_NFE_Source"].iloc[0],
                "Tumour_WT": tumour_counts["WT"],
                "Tumour_Het_ALT": tumour_counts["Het_ALT"],
                "Tumour_Hom_ALT": tumour_counts["Hom_ALT"],
                "Control_WT": normal_counts["WT"],
                "Control_Het_ALT": normal_counts["Het_ALT"],
                "Control_Hom_ALT": normal_counts["Hom_ALT"],
                "Tumour_Callable_Samples": tumour_callable_n,
                "Control_Callable_Samples": normal_callable_n,
                "Tumour_Total": tumour_total,
                "Control_Total": normal_total,
                "Tumour_NonCallable": int(max(0, tumour_total - tumour_callable_n)),
                "Control_NonCallable": int(max(0, normal_total - normal_callable_n)),
                "Tumour_WT_Pct_Callable": _pct(tumour_counts["WT"], tumour_callable_n),
                "Tumour_Het_Pct_Callable": _pct(tumour_counts["Het_ALT"], tumour_callable_n),
                "Tumour_Hom_Pct_Callable": _pct(tumour_counts["Hom_ALT"], tumour_callable_n),
                "Control_WT_Pct_Callable": _pct(normal_counts["WT"], normal_callable_n),
                "Control_Het_Pct_Callable": _pct(normal_counts["Het_ALT"], normal_callable_n),
                "Control_Hom_Pct_Callable": _pct(normal_counts["Hom_ALT"], normal_callable_n),
                "Tumour_Mean_ALT_Dose": tumour_mean_dose,
                "Control_Mean_ALT_Dose": normal_mean_dose,
                "Dose_Difference": dose_diff,
                "Dose_Direction": direction,
                "Chi2_Statistic": float(chi2),
                "Chi2_Degrees_Freedom": int(dof),
                "Expected_Min_Cell": expected_min,
                "Sparse_Cell_Flag": bool(expected_min < 5 if not np.isnan(expected_min) else False),
                "P_Value": float(p_value),
                "Calculation_Method": f"{calc_note}; minimum expected cell {expected_min:.2f}",
                "Variant_Status": "shared",
            })

    if not all_results:
        print("No SNP genotype-association results generated.")
        return

    # --- 7. LD STRUCTURE + MULTIPLE TESTING CORRECTION ---
    res_df = pd.DataFrame(all_results)
    ld_reference = load_common_snp_ld_reference(
        annotated_report=PATHS["annotated_report"],
        phased_genotypes=PATHS["haplotype_phased"],
        min_nfe_af=THRESHOLDS["min_nfe_af"],
        common_snp_whitelist=PATHS["common_snps_dir"] / "GSDMB_Common_SNP_Frequency_Summary.xlsx",
        threshold=LD_R2_THRESHOLD,
    )
    res_df, ld_block_membership = _initialise_ld_block_columns(res_df, ld_reference)
    res_df = _refresh_block_member_labels(res_df)
    final_res = _apply_block_fdr(
        res_df,
        group_col="Analysis_Group",
        p_col="P_Value",
        block_col="LD_Block_ID",
        primary_fdr_col="FDR_P_Value",
        legacy_fdr_col="FDR_P_Value_Global_Legacy",
    ).sort_values(["Analysis_Group", "P_Value"])

    # --- 8. DERIVED COLUMNS FOR PLOTTING ---
    final_res["P_Value_Safe"] = final_res["P_Value"].replace(0, 1e-300)
    final_res["FDR_P_Value_Safe"] = final_res["FDR_P_Value"].replace(0, 1e-300)

    final_res["neglog10_raw_p"] = -np.log10(final_res["P_Value_Safe"])
    final_res["neglog10_fdr_p"] = -np.log10(final_res["FDR_P_Value_Safe"])
    final_res["Dose_Difference_Plot"] = final_res["Dose_Difference"].replace([np.inf, -np.inf], np.nan).fillna(0.0)

    # Separate significance flags
    final_res["Raw_Significant"] = final_res["P_Value"] < RAW_P_THRESHOLD
    final_res["FDR_Significant"] = final_res["FDR_P_Value"] < FDR_THRESHOLD
    final_res["FDR_Significant_Global_Legacy"] = final_res["FDR_P_Value_Global_Legacy"] < FDR_THRESHOLD
    final_res["Raw_P_Threshold"] = RAW_P_THRESHOLD
    final_res["FDR_Threshold"] = FDR_THRESHOLD
    final_res["Association_Scheme"] = "WT / Het / Hom genotype distribution versus pooled normal controls (3x2 chi-square)"
    final_res["Comparison_Definition"] = final_res["Analysis_Group"].map({
        "Breast": "Breast tumours versus pooled healthy controls (breast normal + endometrial normal)",
        "Endometrium": "Endometrial tumours versus pooled healthy controls (breast normal + endometrial normal)",
        "Global": "All tumours pooled versus pooled healthy controls (breast normal + endometrial normal)",
    }).fillna("Tumour arm versus pooled healthy controls")
    final_res["Genotype_Definition"] = "WT = 0 ALT alleles; Het_ALT = 1; Hom_ALT = 2"
    validate_percentage_columns(final_res, [
        "Tumour_WT_Pct_Callable", "Tumour_Het_Pct_Callable", "Tumour_Hom_Pct_Callable",
        "Control_WT_Pct_Callable", "Control_Het_Pct_Callable", "Control_Hom_Pct_Callable",
    ], "Script 12 results")
    final_res = attach_amplicon_warning_columns(
        final_res,
        warning_lookup.rename(columns={"Variant_Key": "SNP_ID"}),
        variant_col="SNP_ID",
    )
    final_res = final_res.merge(display_lookup, left_on="SNP_ID", right_on="Variant_ID", how="left")
    final_res = final_res.drop(columns=["Variant_ID"], errors="ignore")
    final_res["rsID"] = final_res["rsID"].replace({"": pd.NA})
    final_res["Display_Label"] = final_res["Display_Label"].fillna(final_res["SNP_ID"])
    final_res["Genomic_Label"] = final_res["Genomic_Label"].fillna(final_res["SNP_ID"])
    final_res = _refresh_block_member_labels(final_res)

    # Separate labels for each plot
    final_res["Label_RAW"] = ""
    final_res["Label_FDR"] = ""

    for group in final_res["Analysis_Group"].unique():
        group_mask = final_res["Analysis_Group"] == group

        raw_sig = final_res[group_mask & final_res["Raw_Significant"]]
        if not raw_sig.empty:
            top_raw = raw_sig.nsmallest(5, "P_Value").index
            final_res.loc[top_raw, "Label_RAW"] = final_res.loc[top_raw, "Display_Label"]

        fdr_sig = final_res[group_mask & final_res["FDR_Significant"]]
        if not fdr_sig.empty:
            top_fdr = fdr_sig.nsmallest(5, "FDR_P_Value").index
            final_res.loc[top_fdr, "Label_FDR"] = final_res.loc[top_fdr, "Display_Label"]

    # Order columns more cleanly for Excel
    preferred_cols = [
        "Analysis_Group", "Display_Label", "rsID", "SNP_ID", "Genomic_Label", "HGVSc", "Association_Scheme", "Comparison_Definition", "Genotype_Definition",
        "LD_Threshold", "Threshold_Label", "LD_Block_ID", "LD_Block_Assignment",
        "LD_Block_Size", "LD_Block_Test_Count", "LD_Block_Members", "FDR_Scope",
        "Coverage_Risk_Flag", "Coverage_Risk_Amplicon", "Coverage_Risk_Region", "Coverage_Risk_Note",
        "Symbol", "Consequence", "Impact", "gnomAD_NFE_AF", "gnomAD_NFE_Source",
        "Tumour_WT", "Tumour_Het_ALT", "Tumour_Hom_ALT", "Control_WT", "Control_Het_ALT", "Control_Hom_ALT",
        "Tumour_Callable_Samples", "Control_Callable_Samples", "Tumour_Total", "Control_Total",
        "Tumour_NonCallable", "Control_NonCallable",
        "Tumour_WT_Pct_Callable", "Tumour_Het_Pct_Callable", "Tumour_Hom_Pct_Callable",
        "Control_WT_Pct_Callable", "Control_Het_Pct_Callable", "Control_Hom_Pct_Callable",
        "Tumour_Mean_ALT_Dose", "Control_Mean_ALT_Dose", "Dose_Difference", "Dose_Direction",
        "Chi2_Statistic", "Chi2_Degrees_Freedom", "Expected_Min_Cell", "Sparse_Cell_Flag",
        "P_Value", "FDR_P_Value", "FDR_P_Value_Global_Legacy",
        "Raw_Significant", "FDR_Significant", "FDR_Significant_Global_Legacy",
        "Calculation_Method", "Variant_Status",
    ]
    other_cols = [c for c in final_res.columns if c not in preferred_cols]
    final_res = final_res[preferred_cols + other_cols]

    # --- 9. SAVE EXCEL RESULTS ---
    def _carrier_stats(row: pd.Series) -> pd.Series:
        tumour_carriers = int(row["Tumour_Het_ALT"] + row["Tumour_Hom_ALT"])
        control_carriers = int(row["Control_Het_ALT"] + row["Control_Hom_ALT"])
        tumour_noncarriers = int(max(0, row["Tumour_Callable_Samples"] - tumour_carriers))
        control_noncarriers = int(max(0, row["Control_Callable_Samples"] - control_carriers))
        _, p_value = fisher_exact(
            [[tumour_carriers, tumour_noncarriers], [control_carriers, control_noncarriers]],
            alternative="two-sided",
        )
        if 0 in {tumour_carriers, tumour_noncarriers, control_carriers, control_noncarriers}:
            odds_ratio = ((tumour_carriers + 0.5) * (control_noncarriers + 0.5)) / ((tumour_noncarriers + 0.5) * (control_carriers + 0.5))
        else:
            odds_ratio = (tumour_carriers * control_noncarriers) / (tumour_noncarriers * control_carriers) if (tumour_noncarriers * control_carriers) != 0 else np.inf
        tumour_freq = _pct(tumour_carriers, row["Tumour_Callable_Samples"])
        control_freq = _pct(control_carriers, row["Control_Callable_Samples"])
        return pd.Series({
            "Carrier_Tumour_Count": tumour_carriers,
            "Carrier_Control_Count": control_carriers,
            "Carrier_Tumour_NonCarrier": tumour_noncarriers,
            "Carrier_Control_NonCarrier": control_noncarriers,
            "Carrier_Tumour_Freq_pct": tumour_freq,
            "Carrier_Control_Freq_pct": control_freq,
            "Carrier_Frequency_Difference_pct": None if tumour_freq is None or control_freq is None else tumour_freq - control_freq,
            "Carrier_Odds_Ratio": odds_ratio,
            "Carrier_P_Value": float(p_value),
            "Carrier_Calculation_Method": "ALT carrier versus WT/non-carrier Fisher exact",
        })

    carrier_res = final_res.copy()
    carrier_extra = carrier_res.apply(_carrier_stats, axis=1)
    carrier_res = pd.concat([carrier_res, carrier_extra], axis=1)
    carrier_res = _refresh_block_member_labels(carrier_res)
    carrier_res["P_Value"] = carrier_res["Carrier_P_Value"]
    carrier_res["Raw_Significant"] = carrier_res["P_Value"] < RAW_P_THRESHOLD
    carrier_res = _apply_block_fdr(
        carrier_res,
        group_col="Analysis_Group",
        p_col="P_Value",
        block_col="LD_Block_ID",
        primary_fdr_col="FDR_P_Value",
        legacy_fdr_col="FDR_P_Value_Global_Legacy",
    ).sort_values(["Analysis_Group", "P_Value"])
    carrier_res["Carrier_FDR_P_Value"] = carrier_res["FDR_P_Value"]
    carrier_res["Carrier_FDR_P_Value_Global_Legacy"] = carrier_res["FDR_P_Value_Global_Legacy"]
    carrier_res["FDR_Significant"] = carrier_res["FDR_P_Value"] < FDR_THRESHOLD
    carrier_res["Raw_P_Threshold"] = RAW_P_THRESHOLD
    carrier_res["FDR_Threshold"] = FDR_THRESHOLD
    carrier_res["Association_Scheme"] = "Any ALT carrier versus WT / non-carrier (2x2 Fisher exact)"
    carrier_res["Comparison_Definition"] = carrier_res["Analysis_Group"].map({
        "Breast": "Breast tumours versus pooled healthy controls (breast normal + endometrial normal)",
        "Endometrium": "Endometrial tumours versus pooled healthy controls (breast normal + endometrial normal)",
        "Global": "All tumours pooled versus pooled healthy controls (breast normal + endometrial normal)",
    }).fillna("Tumour arm versus pooled healthy controls")
    carrier_res["Genotype_Definition"] = "Carrier = Het_ALT or Hom_ALT; non-carrier = WT"
    carrier_res["Carrier_Frequency_Difference_Plot"] = carrier_res["Carrier_Frequency_Difference_pct"].fillna(0.0)
    carrier_res["Carrier_P_Value_Safe"] = carrier_res["P_Value"].replace(0, 1e-300)
    carrier_res["Carrier_FDR_P_Value_Safe"] = carrier_res["FDR_P_Value"].replace(0, 1e-300)
    carrier_res["Carrier_neglog10_raw_p"] = -np.log10(carrier_res["Carrier_P_Value_Safe"])
    carrier_res["Carrier_neglog10_fdr_p"] = -np.log10(carrier_res["Carrier_FDR_P_Value_Safe"])
    carrier_res["Label_RAW"] = ""
    carrier_res["Label_FDR"] = ""
    for group in carrier_res["Analysis_Group"].unique():
        group_mask = carrier_res["Analysis_Group"] == group
        raw_sig = carrier_res[group_mask & carrier_res["Raw_Significant"]]
        if not raw_sig.empty:
            top_raw = raw_sig.nsmallest(5, "P_Value").index
            carrier_res.loc[top_raw, "Label_RAW"] = carrier_res.loc[top_raw, "Display_Label"]
        fdr_sig = carrier_res[group_mask & carrier_res["FDR_Significant"]]
        if not fdr_sig.empty:
            top_fdr = fdr_sig.nsmallest(5, "FDR_P_Value").index
            carrier_res.loc[top_fdr, "Label_FDR"] = carrier_res.loc[top_fdr, "Display_Label"]

    genotype_enrichment_df, genotype_counts_df = compute_snp_genotype_enrichment(long_df, display_lookup, warning_lookup)
    mutation_outputs = compute_mutation_enrichment_outputs(analysis_df)

    snp_summary_parts = []
    carrier_summary = carrier_res[
        carrier_res["Raw_Significant"].fillna(False) | carrier_res["FDR_Significant"].fillna(False)
    ].copy()
    if not carrier_summary.empty:
        carrier_summary["Workflow"] = "carrier_enrichment"
        carrier_summary["Primary_Result_P_Value"] = carrier_summary["P_Value"]
        carrier_summary["Primary_Result_FDR_P_Value"] = carrier_summary["FDR_P_Value"]
        carrier_summary["Primary_Result_Method"] = carrier_summary["Carrier_Calculation_Method"]
        snp_summary_parts.append(carrier_summary)
    if isinstance(genotype_enrichment_df, pd.DataFrame) and "Primary_Genotype_P_Value" in genotype_enrichment_df.columns:
        genotype_summary = genotype_enrichment_df[
            genotype_enrichment_df["Primary_Genotype_Significant"].fillna(False)
            | genotype_enrichment_df["Primary_Genotype_FDR_Significant"].fillna(False)
        ].copy()
        if not genotype_summary.empty:
            genotype_summary["Workflow"] = "genotype_enrichment"
            genotype_summary["Primary_Result_P_Value"] = genotype_summary["Primary_Genotype_P_Value"]
            genotype_summary["Primary_Result_FDR_P_Value"] = genotype_summary["Primary_Genotype_FDR_P_Value"]
            genotype_summary["Primary_Result_Method"] = genotype_summary["Primary_Genotype_Method"]
            snp_summary_parts.append(genotype_summary)
    summary_significant_df = (
        pd.concat(snp_summary_parts, ignore_index=True).sort_values(
            ["Workflow", "Analysis_Group", "Primary_Result_P_Value", "Display_Label"],
            na_position="last",
            kind="stable",
        )
        if snp_summary_parts
        else pd.DataFrame({"Note": ["No nominally significant common-SNP carrier or genotype-enrichment signals were detected."]})
    )

    significant_overview = build_significant_overview_table(carrier_res)
    global_ld_summary, global_pairwise_ld, ld_membership_export = _build_ld_summary_tables(carrier_res, ld_reference, ld_block_membership)
    target_carrier_summary = carrier_res[carrier_res["rsID"].astype(str).str.lower().isin(TARGET_RSIDS)].copy()

    main_workbook_readme = pd.DataFrame([
        {
            "Section": "Primary thesis result",
            "Details": "This workbook focuses on common-SNP tumour-versus-healthy enrichment. The headline sheets are summary_significant, Significant Overview, At_A_Glance, carrier_enrichment, and genotype_enrichment.",
        },
        {
            "Section": "Control definition",
            "Details": "Breast and endometrial tumour comparisons use pooled healthy controls: breast normal plus endometrial normal.",
        },
        {
            "Section": "Secondary outputs",
            "Details": "Mutation enrichment, hotspot summaries, and the PPT-style summary workbook are written to supplementary_results so they do not dilute the main stage-12 story.",
        },
        {
            "Section": "Genotype logic",
            "Details": "WT / Het_ALT / Hom_ALT states come from the forced-genotype backbone, which preserves callable denominators and distinguishes WT from missing calls.",
        },
    ])

    main_workbook_index = pd.DataFrame([
        {"Order": 1, "Sheet": "README", "Priority": "Primary", "Purpose": "Workbook guide and stage-12 scope note"},
        {"Order": 2, "Sheet": "At_A_Glance", "Priority": "Primary", "Purpose": "Compact count of nominal and FDR-supported findings by cohort and workflow"},
        {"Order": 3, "Sheet": "summary_significant", "Priority": "Primary", "Purpose": "Combined significant common-SNP rows from the carrier and genotype workflows"},
        {"Order": 4, "Sheet": "Significant Overview", "Priority": "Primary", "Purpose": "Supervisor-facing overview of the main common-SNP signals"},
        {"Order": 5, "Sheet": "genotype_enrichment", "Priority": "Primary", "Purpose": "Genotype-distribution, genotype-dose, and genotype-contrast common-SNP enrichment results"},
        {"Order": 6, "Sheet": "carrier_enrichment", "Priority": "Primary", "Purpose": "ALT-carrier enrichment results across all analysis groups"},
        {"Order": 7, "Sheet": "genotype_counts", "Priority": "Primary", "Purpose": "Review-ready tumour and pooled-control genotype counts and percentages"},
        {"Order": 8, "Sheet": "Target_Carrier_Summary", "Priority": "Secondary", "Purpose": "Focused carrier summary for the main target SNPs"},
        {"Order": 9, "Sheet": "LD_Block_Membership", "Priority": "Secondary", "Purpose": "LD block membership at r = 0.8 for the common-SNP backbone"},
        {"Order": 10, "Sheet": "Global_Significant_LD", "Priority": "Secondary", "Purpose": "Pairwise LD among significant common SNPs"},
        {"Order": 11, "Sheet": "Global_Block_Context", "Priority": "Secondary", "Purpose": "Compact LD-block context for significant SNPs"},
        {"Order": 12, "Sheet": "Breast", "Priority": "Secondary", "Purpose": "Per-cohort carrier-based breast enrichment rows"},
        {"Order": 13, "Sheet": "Endometrium", "Priority": "Secondary", "Purpose": "Per-cohort carrier-based endometrial enrichment rows"},
        {"Order": 14, "Sheet": "Global", "Priority": "Secondary", "Purpose": "Per-cohort carrier-based pooled tumour enrichment rows"},
    ])
    main_glance_rows = []
    for view_name, df, sig_col in [
        ("Carrier", carrier_res, "FDR_Significant"),
        ("Carrier", carrier_res, "Raw_Significant"),
        ("Genotype", genotype_enrichment_df if isinstance(genotype_enrichment_df, pd.DataFrame) else pd.DataFrame(), "Primary_Genotype_FDR_Significant"),
        ("Genotype", genotype_enrichment_df if isinstance(genotype_enrichment_df, pd.DataFrame) else pd.DataFrame(), "Primary_Genotype_Significant"),
    ]:
        sig_label = "FDR" if "FDR" in sig_col else "Nominal"
        for group_name in ["Breast", "Endometrium", "Global"]:
            if df.empty or "Analysis_Group" not in df.columns or sig_col not in df.columns:
                sub = pd.DataFrame()
            else:
                sub = df[df["Analysis_Group"].astype(str).eq(group_name) & df[sig_col].fillna(False)]
            main_glance_rows.append({
                "View": view_name,
                "Significance": sig_label,
                "Analysis_Group": group_name,
                "Significant_Rows": int(sub.shape[0]),
                "Distinct_rsIDs": int(sub["rsID"].nunique()) if "rsID" in sub.columns else 0,
            })
    main_glance = pd.DataFrame(main_glance_rows)

    def _join_nonempty(values):
        cleaned = []
        for value in values:
            text = " ".join(str(value).split()).strip()
            if not text or text.lower() in {"nan", "none", "na", "<na>"}:
                continue
            if text not in cleaned:
                cleaned.append(text)
        return ", ".join(cleaned)

    ld_block_summary = pd.DataFrame()
    if not global_ld_summary.empty:
        ld_block_summary = (
            global_ld_summary.groupby(["LD_Block_ID", "LD_Block_Assignment"], as_index=False)
            .agg(
                Block_Size=("LD_Block_Size", "first"),
                rsIDs=("rsID", _join_nonempty),
                LD_Block_Members=("LD_Block_Members", "first"),
                FDR_Significant_SNPs=("FDR_Significant", "sum"),
                Min_P_Value=("P_Value", "min"),
                Min_FDR_P_Value=("FDR_P_Value", "min"),
            )
            .sort_values(["FDR_Significant_SNPs", "Min_FDR_P_Value", "LD_Block_ID"], ascending=[False, True, True])
        )

    ld_workbook_index = pd.DataFrame([
        {"Order": 1, "Sheet": "Global_Block_Context", "Purpose": "Block-level genotype context for the significant LD blocks"},
        {"Order": 2, "Sheet": "Global_Significant_LD", "Purpose": "Pairwise LD among significant SNPs"},
        {"Order": 3, "Sheet": "LD_Block_Membership", "Purpose": "LD block membership at r = 0.8"},
        {"Order": 4, "Sheet": "LD_Block_Summary", "Purpose": "Compact summary of blocks, genes, and significant SNPs"},
    ])

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        main_workbook_readme.to_excel(writer, sheet_name="README", index=False)
        main_workbook_index.to_excel(writer, sheet_name="Workbook_Index", index=False)
        main_glance.to_excel(writer, sheet_name="At_A_Glance", index=False)
        summary_significant_df.to_excel(writer, sheet_name="summary_significant", index=False)
        significant_overview.to_excel(writer, sheet_name="Significant Overview", index=False)
        if isinstance(genotype_enrichment_df, pd.DataFrame):
            genotype_enrichment_df.to_excel(writer, sheet_name="genotype_enrichment", index=False)
        carrier_res.to_excel(writer, sheet_name="carrier_enrichment", index=False)
        if isinstance(genotype_counts_df, pd.DataFrame):
            genotype_counts_df.to_excel(writer, sheet_name="genotype_counts", index=False)
        if not target_carrier_summary.empty:
            target_carrier_summary.to_excel(writer, sheet_name="Target_Carrier_Summary", index=False)
        if not ld_membership_export.empty:
            ld_membership_export.to_excel(writer, sheet_name="LD_Block_Membership", index=False)
        if not global_pairwise_ld.empty:
            global_pairwise_ld.to_excel(writer, sheet_name="Global_Significant_LD", index=False)
        if not global_ld_summary.empty:
            global_ld_summary.to_excel(writer, sheet_name="Global_Block_Context", index=False)
        for group in cohort_list:
            if group in carrier_res["Analysis_Group"].unique():
                carrier_res[carrier_res["Analysis_Group"] == group].to_excel(
                    writer, sheet_name=group, index=False
                )

    with pd.ExcelWriter(output_ld_summary_xlsx, engine="openpyxl") as writer:
        ld_workbook_index.to_excel(writer, sheet_name="Workbook_Index", index=False)
        ld_block_summary.to_excel(writer, sheet_name="LD_Block_Summary", index=False)

    def _unique_join(series):
        values = []
        for value in series:
            text_value = " ".join(str(value).split()).strip()
            if not text_value or text_value.lower() in {"nan", "none", "na", "<na>"}:
                continue
            if text_value not in values:
                values.append(text_value)
        return ", ".join(values)

    def _build_ppt_summary_frame(source_df, significance_col):
        sig_rows = source_df[source_df[significance_col].fillna(False)].copy()
        if not sig_rows.empty:
            sig_rows["Comparison"] = sig_rows["Analysis_Group"]
            sig_rows["Gene"] = sig_rows["Symbol"].fillna("")
            sig_rows["Block_Type"] = sig_rows["LD_Block_Assignment"]
            sig_rows["Block_Size"] = sig_rows["LD_Block_Size"]
            row_level = sig_rows[[
                "Comparison", "Analysis_Group", "Association_Scheme", "Comparison_Definition", "Genotype_Definition",
                "Display_Label", "rsID", "SNP_ID", "Gene", "LD_Block_ID", "Block_Type", "Block_Size", "LD_Block_Members",
                "Impact", "P_Value", "FDR_P_Value", "Dose_Direction",
            ]].sort_values(["Comparison", "Block_Type", "LD_Block_ID", "P_Value", "rsID"], kind="stable")
        else:
            row_level = pd.DataFrame(columns=[
                "Comparison", "Analysis_Group", "Association_Scheme", "Comparison_Definition", "Genotype_Definition",
                "Display_Label", "rsID", "SNP_ID", "Gene", "LD_Block_ID", "Block_Type", "Block_Size", "LD_Block_Members",
                "Impact", "P_Value", "FDR_P_Value", "Dose_Direction",
            ])

        if not row_level.empty:
            block_rows = row_level[row_level["Block_Type"] == "LD_Block"].copy()
            block_summary = (
                block_rows.groupby(["LD_Block_ID", "Block_Type"], as_index=False)
                .agg(
                    Block_Size=("Block_Size", "first"),
                    Gene=("Gene", _unique_join),
                    Block_Members=("LD_Block_Members", "first"),
                    rsIDs=("rsID", _unique_join),
                    Comparisons=("Comparison", _unique_join),
                    FDR_Significant_SNPs=("rsID", "nunique"),
                    Min_P_Value=("P_Value", "min"),
                    Min_FDR_P_Value=("FDR_P_Value", "min"),
                )
            )
            comp_counts = (
                block_rows.pivot_table(
                    index="LD_Block_ID",
                    columns="Comparison",
                    values="rsID",
                    aggfunc="nunique",
                    fill_value=0,
                )
                .reset_index()
            )
            block_summary = block_summary.merge(comp_counts, on="LD_Block_ID", how="left")
            singleton_summary = row_level[row_level["Block_Type"] == "Singleton"].copy()
            comparison_counts = (
                row_level.groupby(["Comparison", "Block_Type"], as_index=False)
                .agg(FDR_Significant_SNPs=("rsID", "nunique"))
                .sort_values(["Comparison", "Block_Type"])
            )
        else:
            block_summary = pd.DataFrame(columns=[
                "LD_Block_ID", "Block_Type", "Block_Size", "Gene", "Block_Members", "rsIDs",
                "Comparisons", "FDR_Significant_SNPs", "Min_P_Value", "Min_FDR_P_Value",
            ])
            singleton_summary = row_level.copy()
            comparison_counts = pd.DataFrame(columns=["Comparison", "Block_Type", "FDR_Significant_SNPs"])

        return row_level, block_summary, singleton_summary, comparison_counts

    carrier_fdr_rows, carrier_fdr_block_summary, carrier_fdr_singleton_summary, carrier_fdr_comparison_counts = _build_ppt_summary_frame(
        carrier_res, "FDR_Significant"
    )
    carrier_nominal_rows, carrier_nominal_block_summary, carrier_nominal_singleton_summary, carrier_nominal_comparison_counts = _build_ppt_summary_frame(
        carrier_res, "Raw_Significant"
    )
    fdr_rows, fdr_block_summary, fdr_singleton_summary, fdr_comparison_counts = _build_ppt_summary_frame(
        final_res, "FDR_Significant"
    )
    nominal_rows, nominal_block_summary, nominal_singleton_summary, nominal_comparison_counts = _build_ppt_summary_frame(
        final_res, "Raw_Significant"
    )

    def _comparison_overview(source_df, significance_col, view_label):
        sig_rows = source_df[source_df[significance_col].fillna(False)].copy()
        if sig_rows.empty:
            return pd.DataFrame(columns=[
                "View", "Significance", "Analysis_Group", "Association_Scheme", "Genotype_Definition",
                "Significant_SNPs", "rsIDs", "Block_Ids",
            ])
        summary = (
            sig_rows.groupby(["Analysis_Group", "Association_Scheme", "Genotype_Definition"], as_index=False)
            .agg(
                Significant_SNPs=("rsID", "nunique"),
                rsIDs=("rsID", _unique_join),
                Block_Ids=("LD_Block_ID", _unique_join),
            )
            .sort_values(["Analysis_Group", "Significant_SNPs", "rsIDs"], ascending=[True, False, True])
        )
        summary.insert(0, "View", view_label)
        summary.insert(1, "Significance", "FDR" if significance_col == "FDR_Significant" else "Nominal")
        return summary

    comparison_overview = pd.concat([
        _comparison_overview(carrier_res, "FDR_Significant", "Carrier"),
        _comparison_overview(carrier_res, "Raw_Significant", "Carrier"),
    ], ignore_index=True)

    workbook_index = pd.DataFrame([
        {"Order": 1, "Sheet": "Comparison_Overview", "Purpose": "Carrier associations by cohort and significance"},
        {"Order": 2, "Sheet": "Carrier_FDR_Significant_Rows", "Purpose": "Row-level carrier FDR hits"},
        {"Order": 3, "Sheet": "Carrier_Nominal_Rows", "Purpose": "Row-level carrier nominal hits"},
        {"Order": 4, "Sheet": "Notes", "Purpose": "Definitions of carrier, pooled normal, and thresholds"},
    ])

    at_a_glance = (
        comparison_overview.groupby(["View", "Significance", "Analysis_Group"], as_index=False)
        .agg(Significant_SNP_Entries=("Significant_SNPs", "sum"))
        .sort_values(["View", "Significance", "Analysis_Group"], kind="stable")
    )

    with pd.ExcelWriter(output_ppt_summary_xlsx) as writer:
        workbook_index.to_excel(writer, sheet_name="Workbook_Index", index=False)
        at_a_glance.to_excel(writer, sheet_name="At_A_Glance", index=False)
        carrier_fdr_rows.to_excel(writer, sheet_name="Carrier_FDR_Significant_Rows", index=False)
        carrier_fdr_block_summary.to_excel(writer, sheet_name="Carrier_FDR_Block_Summary", index=False)
        carrier_fdr_singleton_summary.to_excel(writer, sheet_name="Carrier_FDR_Singleton_Summary", index=False)
        carrier_fdr_comparison_counts.to_excel(writer, sheet_name="Carrier_FDR_Comparison_Counts", index=False)
        carrier_nominal_rows.to_excel(writer, sheet_name="Carrier_Nominal_Rows", index=False)
        carrier_nominal_block_summary.to_excel(writer, sheet_name="Carrier_Nominal_Blocks", index=False)
        carrier_nominal_singleton_summary.to_excel(writer, sheet_name="Carrier_Nominal_Singletons", index=False)
        carrier_nominal_comparison_counts.to_excel(writer, sheet_name="Carrier_Nominal_Counts", index=False)
        comparison_overview.to_excel(writer, sheet_name="Comparison_Overview", index=False)
        notes = pd.DataFrame({
            "Item": [
                "Threshold",
                "Block LD threshold",
                "Control definition",
                "Carrier interpretation",
                "Nominal significance",
                "Carrier summary",
            ],
            "Value": [
                f"FDR < {FDR_THRESHOLD}",
                f"r >= {LD_R2_THRESHOLD}",
                "Pooled normal = breast normal + endometrial normal",
                "ALT carrier = Het_ALT or Hom_ALT; non-carrier = WT",
                f"Raw P < {RAW_P_THRESHOLD}",
                "Carrier sheets test ALT presence versus WT/non-carrier",
            ],
        })
        notes.to_excel(writer, sheet_name="Notes", index=False)

    # --- 10. CARRIER VOLCANO HELPERS ---
    def add_volcano_decorations(ax, group_name, data, y_col, label_col, threshold_p, y_label, x_col, x_label, zero_label="No frequency difference"):
        """Add threshold lines, shading and labels to one volcano axis."""
        y_thresh = -np.log10(threshold_p)

        ax.axhline(
            y_thresh,
            color=arm_color('Tumour'),
            linestyle="--",
            alpha=0.7,
            linewidth=2,
            label=f"Threshold = {threshold_p}"
        )

        ax.axvline(
            0,
            color=arm_color('Control'),
            linestyle="--",
            alpha=0.7,
            linewidth=2,
            label=zero_label
        )

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        ax.axvspan(xlim[0], 0, alpha=0.05, color=arm_color('Control'), label="Lower in tumour")
        ax.axvspan(0, xlim[1], alpha=0.05, color=arm_color('Tumour'), label="Higher in tumour")
        ax.axhspan(y_thresh, ylim[1], alpha=0.08, color=arm_color('Tumour'), label="Above Threshold")

        group_data = data[data["Analysis_Group"] == group_name]
        labelled = group_data[group_data[label_col] != ""]
        if not group_data.empty:
            ax.set_title(f"{group_name}\nCarrier frequency (ALT presence) vs pooled normals", fontweight='bold')

        _annotate_ranked_hits(ax, labelled, label_col, x_col, y_col, max_labels=6)

        ax.set_ylabel(y_label, fontweight="bold")
        ax.set_xlabel(x_label, fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.4)
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), loc="lower left", fontsize=7.5, frameon=True)

    def make_volcano_plot(data, y_col, label_col, threshold_p, y_label, title, output_path, x_col, x_label):
        group_order = ["Breast", "Endometrium", "Global"]
        plot_data = data[data["Analysis_Group"].isin(group_order)].copy()
        plot_data["Analysis_Group"] = pd.Categorical(
            plot_data["Analysis_Group"],
            categories=group_order,
            ordered=True
        )
        plot_data = plot_data.sort_values("Analysis_Group")

        g = sns.FacetGrid(
            plot_data,
            col="Analysis_Group",
            hue="Impact",
            col_order=group_order,
            palette=IMPACT_COLORS,
            height=6.4,
            aspect=1.52,
            despine=False,
            sharex=True,
            sharey=True
        )

        g.map_dataframe(
            sns.scatterplot,
            x=x_col,
            y=y_col,
            s=120,
            edgecolor="black",
            alpha=0.8
        )

        for ax, group_name in zip(g.axes.flat, group_order):
            add_volcano_decorations(
                ax=ax,
                group_name=group_name,
                data=plot_data,
                y_col=y_col,
                label_col=label_col,
                threshold_p=threshold_p,
                y_label=y_label,
                x_col=x_col,
                x_label=x_label,
            )

        g.add_legend(title="VEP Impact")
        if g._legend is not None:
            g._legend.set_bbox_to_anchor((0.90, 0.56), transform=g.fig.transFigure)
            g._legend._loc = 6
            g._legend.get_title().set_fontsize(11)
            for text in g._legend.texts:
                text.set_fontsize(10)

        g.fig.subplots_adjust(top=0.82, right=0.86, left=0.055, bottom=0.12, wspace=0.18)
        g.fig.suptitle(tagged_title(f"SNP Carrier Association Volcano: {title} (ALT carrier vs non-carrier, pooled normal)", COMPARATIVE_TAG), fontsize=16, fontweight="bold", y=0.97)

        g.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(g.fig)

    def _build_genotype_contrast_frame(source_df, alt_col, contrast_name, contrast_label):
        contrast_res = source_df.copy()
        tumour_alt_col = f"Tumour_{alt_col}"
        control_alt_col = f"Control_{alt_col}"
        contrast_res["Contrast_Name"] = contrast_name
        contrast_res["Contrast_Label"] = contrast_label
        contrast_res["Tumour_Contrast_Count"] = contrast_res[tumour_alt_col].astype(int)
        contrast_res["Control_Contrast_Count"] = contrast_res[control_alt_col].astype(int)
        contrast_res["Tumour_Reference_Count"] = contrast_res["Tumour_WT"].astype(int)
        contrast_res["Control_Reference_Count"] = contrast_res["Control_WT"].astype(int)
        contrast_res["Tumour_Contrast_Pct_Callable"] = pd.Series([
            _pct(alt, total) for alt, total in zip(contrast_res["Tumour_Contrast_Count"], contrast_res["Tumour_Callable_Samples"])
        ], index=contrast_res.index, dtype="float")
        contrast_res["Control_Contrast_Pct_Callable"] = pd.Series([
            _pct(alt, total) for alt, total in zip(contrast_res["Control_Contrast_Count"], contrast_res["Control_Callable_Samples"])
        ], index=contrast_res.index, dtype="float")
        contrast_res["Contrast_Frequency_Difference_pct"] = contrast_res["Tumour_Contrast_Pct_Callable"].astype(float) - contrast_res["Control_Contrast_Pct_Callable"].astype(float)

        p_values = []
        odds_ratios = []
        calc_methods = []
        for _, row in contrast_res.iterrows():
            table = np.array([
                [int(row["Tumour_Contrast_Count"]), int(row["Tumour_Reference_Count"])],
                [int(row["Control_Contrast_Count"]), int(row["Control_Reference_Count"])],
            ], dtype=float)
            odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
            if not np.isfinite(odds_ratio):
                odds_ratio = ((table[0, 0] + 0.5) * (table[1, 1] + 0.5)) / ((table[0, 1] + 0.5) * (table[1, 0] + 0.5))
            p_values.append(float(p_value))
            odds_ratios.append(float(odds_ratio))
            calc_methods.append(f"{contrast_label} Fisher exact")

        contrast_res["P_Value"] = p_values
        contrast_res["Odds_Ratio"] = odds_ratios
        contrast_res["Calculation_Method"] = calc_methods
        contrast_res["Contrast_Frequency_Difference_Plot"] = contrast_res["Contrast_Frequency_Difference_pct"].fillna(0.0)
        contrast_res["P_Value_Safe"] = contrast_res["P_Value"].replace(0, 1e-300)
        contrast_res["neglog10_raw_p"] = -np.log10(contrast_res["P_Value_Safe"])
        contrast_res["Raw_Significant"] = contrast_res["P_Value"] < RAW_P_THRESHOLD
        contrast_res = _apply_block_fdr(
            contrast_res,
            group_col="Analysis_Group",
            p_col="P_Value",
            block_col="LD_Block_ID",
            primary_fdr_col="FDR_P_Value",
            legacy_fdr_col="FDR_P_Value_Global_Legacy",
        ).sort_values(["Analysis_Group", "P_Value"])
        contrast_res["FDR_P_Value_Safe"] = contrast_res["FDR_P_Value"].replace(0, 1e-300)
        contrast_res["neglog10_fdr_p"] = -np.log10(contrast_res["FDR_P_Value_Safe"])
        contrast_res["FDR_Significant"] = contrast_res["FDR_P_Value"] < FDR_THRESHOLD
        contrast_res["Label_RAW"] = ""
        contrast_res["Label_FDR"] = ""
        contrast_res["Association_Scheme"] = f"{contrast_label} versus WT (2x2 Fisher exact)"
        contrast_res["Comparison_Definition"] = contrast_res["Analysis_Group"].map({
            "Breast": "Breast tumours versus pooled healthy controls (breast normal + endometrial normal)",
            "Endometrium": "Endometrial tumours versus pooled healthy controls (breast normal + endometrial normal)",
            "Global": "All tumours pooled versus pooled healthy controls (breast normal + endometrial normal)",
        }).fillna("Tumour arm versus pooled healthy controls")
        contrast_res["Genotype_Definition"] = f"{contrast_label} compared with WT"
        for group in contrast_res["Analysis_Group"].unique():
            group_mask = contrast_res["Analysis_Group"] == group
            raw_sig = contrast_res[group_mask & contrast_res["Raw_Significant"]]
            if not raw_sig.empty:
                top_raw = raw_sig.nsmallest(5, "P_Value").index
                contrast_res.loc[top_raw, "Label_RAW"] = contrast_res.loc[top_raw, "Display_Label"]
            fdr_sig = contrast_res[group_mask & contrast_res["FDR_Significant"]]
            if not fdr_sig.empty:
                top_fdr = fdr_sig.nsmallest(5, "FDR_P_Value").index
                contrast_res.loc[top_fdr, "Label_FDR"] = contrast_res.loc[top_fdr, "Display_Label"]
        return contrast_res

    def make_genotype_contrast_volcano_plot(data, y_col, label_col, threshold_p, y_label, title, output_path, x_col, x_label, contrast_label):
        group_order = ["Breast", "Endometrium", "Global"]
        plot_data = data[data["Analysis_Group"].isin(group_order)].copy()
        plot_data["Analysis_Group"] = pd.Categorical(
            plot_data["Analysis_Group"],
            categories=group_order,
            ordered=True
        )
        plot_data = plot_data.sort_values("Analysis_Group")

        g = sns.FacetGrid(
            plot_data,
            col="Analysis_Group",
            hue="Impact",
            col_order=group_order,
            palette=IMPACT_COLORS,
            height=6.4,
            aspect=1.52,
            despine=False,
            sharex=True,
            sharey=True
        )

        g.map_dataframe(
            sns.scatterplot,
            x=x_col,
            y=y_col,
            s=120,
            edgecolor="black",
            alpha=0.8
        )

        for ax, group_name in zip(g.axes.flat, group_order):
            group_data = plot_data[plot_data["Analysis_Group"] == group_name]
            labelled = group_data[group_data[label_col] != ""]
            if not group_data.empty:
                ax.set_title(f"{group_name}\n{contrast_label} vs pooled normals", fontweight='bold')
            _annotate_ranked_hits(ax, labelled, label_col, x_col, y_col, max_labels=6)
            y_thresh = -np.log10(threshold_p)
            ax.axhline(y_thresh, color=arm_color('Tumour'), linestyle="--", alpha=0.7, linewidth=2, label=f"Threshold = {threshold_p}")
            ax.axvline(0, color=arm_color('Control'), linestyle="--", alpha=0.7, linewidth=2, label="No genotype difference")
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            ax.axvspan(xlim[0], 0, alpha=0.05, color=arm_color('Control'), label="Lower in tumour")
            ax.axvspan(0, xlim[1], alpha=0.05, color=arm_color('Tumour'), label="Higher in tumour")
            ax.axhspan(y_thresh, ylim[1], alpha=0.08, color=arm_color('Tumour'), label="Above Threshold")
            ax.set_ylabel(y_label, fontweight="bold")
            ax.set_xlabel(x_label, fontweight="bold")
            ax.grid(True, linestyle=":", alpha=0.4)
            handles, labels = ax.get_legend_handles_labels()
            unique = dict(zip(labels, handles))
            ax.legend(unique.values(), unique.keys(), loc="lower left", fontsize=7.5, frameon=True)

        g.add_legend(title="VEP Impact")
        if g._legend is not None:
            g._legend.set_bbox_to_anchor((0.90, 0.56), transform=g.fig.transFigure)
            g._legend._loc = 6
            g._legend.get_title().set_fontsize(11)
            for text in g._legend.texts:
                text.set_fontsize(10)

        g.fig.subplots_adjust(top=0.82, right=0.86, left=0.055, bottom=0.12, wspace=0.18)
        g.fig.suptitle(tagged_title(f"SNP Genotype Association Volcano: {title} ({contrast_label}, pooled normal)", COMPARATIVE_TAG), fontsize=16, fontweight="bold", y=0.97)

        g.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(g.fig)

    # --- 11. CARRIER RAW P-VALUE VOLCANO ---
    make_volcano_plot(
        data=carrier_res,
        y_col="Carrier_neglog10_raw_p",
        label_col="Label_RAW",
        threshold_p=RAW_P_THRESHOLD,
        y_label="-log10(Carrier Raw P-value)",
        title="Raw P-values",
        output_path=output_plot_raw,
        x_col="Carrier_Frequency_Difference_Plot",
        x_label="Carrier frequency difference (tumour - pooled normal, pp)",
    )

    # --- 12. CARRIER FDR-ADJUSTED VOLCANO ---
    make_volcano_plot(
        data=carrier_res,
        y_col="Carrier_neglog10_fdr_p",
        label_col="Label_FDR",
        threshold_p=FDR_THRESHOLD,
        y_label="-log10(Carrier FDR-adjusted P-value)",
        title="FDR-adjusted P-values",
        output_path=output_plot_fdr,
        x_col="Carrier_Frequency_Difference_Plot",
        x_label="Carrier frequency difference (tumour - pooled normal, pp)",
    )

    het_res = _build_genotype_contrast_frame(final_res, "Het_ALT", "Het_vs_WT", "Het vs WT")
    hom_res = _build_genotype_contrast_frame(final_res, "Hom_ALT", "Hom_vs_WT", "Hom vs WT")

    # --- 13. HET VS WT VOLCANOES ---
    make_genotype_contrast_volcano_plot(
        data=het_res,
        y_col="neglog10_raw_p",
        label_col="Label_RAW",
        threshold_p=RAW_P_THRESHOLD,
        y_label="-log10(Het vs WT Raw P-value)",
        title="Het vs WT Raw P-values",
        output_path=output_plot_het_raw,
        x_col="Contrast_Frequency_Difference_Plot",
        x_label="Het frequency difference (tumour - pooled normal, pp)",
        contrast_label="Het vs WT",
    )

    make_genotype_contrast_volcano_plot(
        data=het_res,
        y_col="neglog10_fdr_p",
        label_col="Label_FDR",
        threshold_p=FDR_THRESHOLD,
        y_label="-log10(Het vs WT FDR-adjusted P-value)",
        title="Het vs WT FDR-adjusted P-values",
        output_path=output_plot_het_fdr,
        x_col="Contrast_Frequency_Difference_Plot",
        x_label="Het frequency difference (tumour - pooled normal, pp)",
        contrast_label="Het vs WT",
    )

    # --- 14. HOM VS WT VOLCANOES ---
    make_genotype_contrast_volcano_plot(
        data=hom_res,
        y_col="neglog10_raw_p",
        label_col="Label_RAW",
        threshold_p=RAW_P_THRESHOLD,
        y_label="-log10(Hom vs WT Raw P-value)",
        title="Hom vs WT Raw P-values",
        output_path=output_plot_hom_raw,
        x_col="Contrast_Frequency_Difference_Plot",
        x_label="Hom frequency difference (tumour - pooled normal, pp)",
        contrast_label="Hom vs WT",
    )

    make_genotype_contrast_volcano_plot(
        data=hom_res,
        y_col="neglog10_fdr_p",
        label_col="Label_FDR",
        threshold_p=FDR_THRESHOLD,
        y_label="-log10(Hom vs WT FDR-adjusted P-value)",
        title="Hom vs WT FDR-adjusted P-values",
        output_path=output_plot_hom_fdr,
        x_col="Contrast_Frequency_Difference_Plot",
        x_label="Hom frequency difference (tumour - pooled normal, pp)",
        contrast_label="Hom vs WT",
    )

    if isinstance(genotype_enrichment_df, pd.DataFrame) and "Primary_Genotype_P_Value" in genotype_enrichment_df.columns:
        plot_snp_genotype_heatmap(genotype_enrichment_df, output_plot_genotype_heatmap)
        genotype_composition_paths = plot_snp_genotype_composition(long_df, genotype_enrichment_df, output_plot_genotype_composition_prefix)
    else:
        genotype_composition_paths = []

    mutation_readme = pd.DataFrame({
        "Item": [
            "Workbook purpose",
            "Mutation definition",
            "Per-mutation test",
            "Burden test",
            "Hotspot summary",
            "Optional genotype treatment",
            "Multiple testing",
        ],
        "Value": [
            "Rare / mutation-focused enrichment output that complements the common-SNP workbook.",
            "Mutation workflow = non-common variant by gnomAD NFE (AF <= 1%) or AF missing, using the same forced-genotype backbone for denominators.",
            "Carrier enrichment uses two-sided Fisher exact tests on tumour versus pooled healthy controls.",
            "Burden enrichment uses sample-level burden flags with two-sided Fisher exact tests and Haldane-Anscombe OR correction when needed.",
            "Hotspots summarise recurrent mutation positions and test tumour versus healthy carrier enrichment at the position level.",
            f"Mutation genotype-distribution testing is only exported when at least {MUTATION_GENOTYPE_MIN_ALT_COUNT} ALT genotype observations are available across tumour and healthy groups.",
            "BH FDR is applied separately within each analysis group for per-mutation carrier tests, burden tests, and hotspot summaries.",
        ],
    })

    with pd.ExcelWriter(output_mutation_xlsx) as writer:
        mutation_readme.to_excel(writer, sheet_name="README", index=False)
        for sheet_name in [
            "mutation_carrier_enrichment",
            "mutation_burden",
            "mutation_hotspots",
            "mutation_summary_significant",
        ]:
            mutation_outputs[sheet_name].to_excel(writer, sheet_name=sheet_name, index=False)

    plot_mutation_burden_forest(mutation_outputs.get("mutation_burden", pd.DataFrame()), output_plot_mutation_burden)
    plot_mutation_hotspots(mutation_outputs.get("mutation_hotspots", pd.DataFrame()), output_plot_mutation_hotspots)

    print("\n>>> SUCCESS: All SNPs processed.")
    print(f">>> Results saved to: {output_xlsx}")
    print(f">>> Mutation enrichment workbook saved to: {output_mutation_xlsx}")
    print(f">>> LD-block summary saved to: {output_ld_summary_xlsx}")
    print(">>> Genotype logic used:")
    print("    - WT / Het_ALT / Hom_ALT from the forced-genotype workbook")
    print("    - pooled healthy controls = breast normal + endometrial normal")
    print(">>> Figures saved:")
    print(f"    - Carrier RAW: {output_plot_raw}")
    print(f"    - Carrier FDR: {output_plot_fdr}")
    print(f"    - Het vs WT RAW: {output_plot_het_raw}")
    print(f"    - Het vs WT FDR: {output_plot_het_fdr}")
    print(f"    - Hom vs WT RAW: {output_plot_hom_raw}")
    print(f"    - Hom vs WT FDR: {output_plot_hom_fdr}")
    print(f"    - SNP genotype heatmap: {output_plot_genotype_heatmap}")
    for path in genotype_composition_paths:
        print(f"    - SNP genotype composition: {path}")
    print(f"    - Mutation burden forest: {output_plot_mutation_burden}")
    print(f"    - Mutation hotspot summary: {output_plot_mutation_hotspots}")
    print(f"    - PPT summary workbook: {output_ppt_summary_xlsx}")


if __name__ == "__main__":
    run_snp_association_analysis()
