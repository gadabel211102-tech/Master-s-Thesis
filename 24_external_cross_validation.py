#!/usr/bin/env python3
"""Script 24: external SNP and haplotype cross-validation.

This stage is a secondary validation layer for thesis reporting. It does not
perform primary discovery. Instead, it reuses the established study backbone:

- common / established SNPs from stage 11 and the cohort-specific enrichment
  results from stage 12
- full-region and LD-block haplotypes from stage 15
- the female-filtered 1000 Genomes reference logic already used in stage 19

The active external comparison design is:

- study breast tumours versus female-only 1000 Genomes
- study endometrial tumours versus female-only 1000 Genomes

The earlier TCGA branch was removed from the active workflow because the
relevant TCGA germline genotype files were controlled-access and were not
available for this project.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
from textwrap import shorten

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns

from association_runtime import script24_defaults
from pipeline_utils import ensure_directory

sns.set_theme(style="whitegrid", context="talk")

DATASET_STYLE = {
    "Reference": {"label": "1000 Genomes female", "color": "#2166ac", "marker": "o"},
    "Study": {"label": "Study tumour", "color": "#b2182b", "marker": "s"},
}

COHORT_CONFIG = {
    "Breast": {"study_dataset": "Study_Breast_Tumour", "study_manifest_file": "24_Study_Breast_Tumour_Manifest.tsv"},
    "Endometrium": {"study_dataset": "Study_Endometrium_Tumour", "study_manifest_file": "24_Study_Endometrium_Tumour_Manifest.tsv"},
}


def load_stage19_module():
    module_path = Path(__file__).with_name("19_1000g_haplotype_comparison.py")
    spec = importlib.util.spec_from_file_location("stage19_helpers", module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load stage-19 helpers from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    defaults = script24_defaults()
    parser = argparse.ArgumentParser(
        description="Cross-validate breast and endometrial tumour SNP / haplotype signals against female-filtered 1000 Genomes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--vcf", required=True, help="Path to the phased 1000 Genomes VCF/BCF (indexed).")
    parser.add_argument("--panel", required=True, help="Path to the 1000 Genomes population panel.")
    parser.add_argument("--haplo-results", default=str(defaults["haplotype_results"]), help="Stage-15 haplotype workbook.")
    parser.add_argument("--snp-results", default=str(defaults["snp_results"]), help="Stage-12 SNP enrichment workbook.")
    parser.add_argument("--sample-metadata", default=str(defaults["study_sample_metadata"]), help="Stage-15 study sample metadata TSV.")
    parser.add_argument("--female-study-samples", default=str(defaults["study_female_samples"]), help="Explicit female-by-design study sample list.")
    parser.add_argument("--reference-population", default=str(defaults["reference_population"]), help="Primary 1000 Genomes population or super-population to use as the healthy reference.")
    parser.add_argument("--out-dir", default=str(defaults["out_dir"]), help="Directory where cross-validation outputs will be written.")
    parser.add_argument("--contig", default="17", help="Primary contig name to search in the external VCF.")
    return parser.parse_args()


def get_header_samples(vcf_path: Path) -> list[str]:
    with pysam.VariantFile(vcf_path) as vcf:
        return list(vcf.header.samples)


def load_study_sample_metadata(sample_metadata_path: Path, female_sample_file: Path | None) -> pd.DataFrame:
    metadata = pd.read_csv(sample_metadata_path, sep="\t")
    required = {"Sample", "Cohort", "Tissue"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Study sample metadata is missing columns: {sorted(missing)}")
    metadata = metadata.copy()
    metadata["Sample"] = metadata["Sample"].astype(str).str.strip()
    metadata["Cohort"] = metadata["Cohort"].astype(str).str.strip()
    metadata["Tissue"] = metadata["Tissue"].astype(str).str.strip()
    metadata["Female_By_Design_Filter_Used"] = True
    if female_sample_file and female_sample_file.exists():
        female_df = pd.read_csv(female_sample_file, sep="\t", header=None, names=["Sample"])
        female_samples = {str(value).strip() for value in female_df["Sample"] if str(value).strip()}
        metadata["Female_By_Design_Filter_Used"] = metadata["Sample"].isin(female_samples)
    return metadata


def select_study_manifest(metadata: pd.DataFrame, cohort: str, female_filter_active: bool) -> tuple[pd.DataFrame, dict[str, object]]:
    selected = metadata[metadata["Cohort"].eq(cohort) & metadata["Tissue"].eq("Tumour")].copy()
    if female_filter_active:
        selected = selected[selected["Female_By_Design_Filter_Used"]].copy()
        sex_status = "Filtered with explicit female-by-design study sample list"
    else:
        sex_status = "Female-by-design file not provided"
    selected["CrossValidation_Group"] = COHORT_CONFIG[cohort]["study_dataset"]
    summary = {
        "Dataset": COHORT_CONFIG[cohort]["study_dataset"],
        "Comparison": cohort,
        "Requested_Sex_Filter": "Female-by-design study subset",
        "Sex_Filter_Feasible": True,
        "Sex_Filter_Applied": female_filter_active,
        "Sex_Filter_Status": sex_status,
        "Age_Filter_Requested": False,
        "Age_Filter_Feasible": False,
        "Age_Filter_Applied": False,
        "Age_Filter_Status": "No age filter is used in the active 1000G-only validation step",
        "Selected_Samples": int(selected.shape[0]),
        "Notes": f"Study {cohort.lower()} tumours are compared against female-only 1000 Genomes only.",
    }
    return selected.sort_values("Sample"), summary


def load_significant_snps(snp_results_path: Path, comparison: str) -> pd.DataFrame:
    df = pd.read_excel(snp_results_path, sheet_name=comparison)
    if "FDR_Significant" not in df.columns:
        raise ValueError(f"Stage-12 {comparison} sheet is missing the FDR_Significant column.")
    df = df[df["FDR_Significant"].eq(True)].copy()
    df["SNP_ID"] = df["SNP_ID"].astype(str).str.strip()
    if "Gene" not in df.columns:
        df["Gene"] = df.get("Symbol", pd.Series(index=df.index)).astype(str).str.strip()
    df["Study_Tumour_Frequency"] = pd.to_numeric(df["Tumour_Freq_%"], errors="coerce") / 100.0
    df["Study_Pooled_Control_Frequency"] = pd.to_numeric(df["Control_Freq_%"], errors="coerce") / 100.0
    df["Comparison"] = comparison
    return df.sort_values(["FDR_P_Value", "P_Value", "SNP_ID"]).reset_index(drop=True)


def load_haplotype_catalogues(haplo_results_path: Path, stage19, comparison: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    snp_meta, ld_summary, _, _ = stage19.load_local_haplotype_data(haplo_results_path)
    workbook = pd.ExcelFile(haplo_results_path)

    full_freq = pd.read_excel(workbook, sheet_name="Full_Global_Freqs").copy()
    full_assoc = pd.read_excel(workbook, sheet_name="Full_Associations").copy()
    full_assoc = full_assoc[full_assoc["Comparison"].astype(str).eq(comparison)].copy()
    full_freq["Allele_String"] = full_freq["Allele_String"].map(lambda value: stage19.normalise_allele_string(value, snp_meta.shape[0]))
    full_catalogue = full_freq.merge(full_assoc[["Haplotype_ID", "Haplotype_String", "Haplotype_Def", "Freq_Tumour", "Freq_Healthy", "P_Fisher", "FDR", "Significant_Fisher", "Significant_FDR"]], on="Haplotype_ID", how="left")
    full_catalogue["Comparison"] = comparison
    full_catalogue["Haplotype_Set"] = "Full_Region"
    full_catalogue["Haplotype_Definition_Type"] = "Full_Region"
    full_catalogue["Region"] = "Full_Region"
    full_catalogue["Threshold_Label"] = pd.NA
    full_catalogue["Block"] = pd.NA
    full_catalogue["N_SNPs"] = snp_meta.shape[0]
    full_catalogue["SNPs"] = "; ".join(snp_meta["rsID_clean"].astype(str))
    full_catalogue["Defining_SNPs"] = full_catalogue["SNPs"]
    full_catalogue["Genes"] = "; ".join(sorted(set(snp_meta["Gene"].astype(str))))
    full_catalogue["Study_Tumour_Frequency"] = pd.to_numeric(full_catalogue["Freq_Tumour"], errors="coerce")
    full_catalogue["Study_Pooled_Control_Frequency"] = pd.to_numeric(full_catalogue["Freq_Healthy"], errors="coerce")

    ld_rows = []
    for row in ld_summary.itertuples(index=False):
        region = str(row.Region)
        n_snps = int(row.N_SNPs)
        freq_sheet = f"{region}_Freqs"
        assoc_sheet = f"{region}_Assoc"
        if freq_sheet not in workbook.sheet_names:
            continue
        freq_df = pd.read_excel(workbook, sheet_name=freq_sheet).copy()
        freq_df["Allele_String"] = freq_df["Allele_String"].map(lambda value: stage19.normalise_allele_string(value, n_snps))
        assoc_df = pd.read_excel(workbook, sheet_name=assoc_sheet).copy() if assoc_sheet in workbook.sheet_names else pd.DataFrame()
        if not assoc_df.empty:
            assoc_df = assoc_df[assoc_df["Comparison"].astype(str).eq(comparison)].copy()
        assoc_cols = ["Haplotype_ID", "Haplotype_String", "Haplotype_Def", "Freq_Tumour", "Freq_Healthy", "P_Fisher", "FDR", "Significant_Fisher", "Significant_FDR"]
        merged = freq_df.merge(assoc_df[assoc_cols] if not assoc_df.empty else pd.DataFrame(columns=assoc_cols), on="Haplotype_ID", how="left")
        merged["Comparison"] = comparison
        merged["Haplotype_Set"] = "LD_Block"
        merged["Haplotype_Definition_Type"] = "LD_Block"
        merged["Region"] = region
        merged["Threshold_Label"] = str(row.Threshold_Label)
        merged["Block"] = region.split("_")[-1]
        merged["N_SNPs"] = n_snps
        merged["SNPs"] = str(row.SNPs).replace(", ", "; ")
        merged["Defining_SNPs"] = merged["SNPs"]
        merged["Genes"] = str(row.Genes).replace(", ", "; ")
        merged["Study_Tumour_Frequency"] = pd.to_numeric(merged["Freq_Tumour"], errors="coerce")
        merged["Study_Pooled_Control_Frequency"] = pd.to_numeric(merged["Freq_Healthy"], errors="coerce")
        ld_rows.append(merged)
    ld_catalogue = pd.concat(ld_rows, ignore_index=True) if ld_rows else pd.DataFrame()
    return snp_meta, ld_summary, full_catalogue, ld_catalogue


def build_full_region_definition_table(full_catalogue: pd.DataFrame) -> pd.DataFrame:
    if full_catalogue.empty:
        return pd.DataFrame()
    columns = ["Comparison", "Haplotype_Definition_Type", "Haplotype_Set", "Region", "N_SNPs", "Defining_SNPs", "Genes", "Haplotype_ID", "Allele_String", "Haplotype_Def", "Variant_Content", "Study_Tumour_Frequency", "Study_Pooled_Control_Frequency", "P_Fisher", "FDR", "Significant_Fisher", "Significant_FDR"]
    return full_catalogue[[c for c in columns if c in full_catalogue.columns]].sort_values(["Comparison", "Haplotype_ID", "Allele_String"], kind="stable").reset_index(drop=True)


def build_ld_block_definition_table(ld_catalogue: pd.DataFrame) -> pd.DataFrame:
    if ld_catalogue.empty:
        return pd.DataFrame()
    columns = ["Comparison", "Haplotype_Definition_Type", "Haplotype_Set", "Region", "Block", "Threshold_Label", "N_SNPs", "Defining_SNPs", "Genes"]
    return ld_catalogue[[c for c in columns if c in ld_catalogue.columns]].drop_duplicates().sort_values(["Comparison", "Region", "Block"], kind="stable").reset_index(drop=True)


def build_ld_block_haplotype_table(ld_catalogue: pd.DataFrame) -> pd.DataFrame:
    if ld_catalogue.empty:
        return pd.DataFrame()
    columns = ["Comparison", "Haplotype_Definition_Type", "Haplotype_Set", "Region", "Block", "Threshold_Label", "N_SNPs", "Defining_SNPs", "Genes", "Haplotype_ID", "Allele_String", "Haplotype_Def", "Variant_Content", "Study_Tumour_Frequency", "Study_Pooled_Control_Frequency", "P_Fisher", "FDR", "Significant_Fisher", "Significant_FDR"]
    return ld_catalogue[[c for c in columns if c in ld_catalogue.columns]].sort_values(["Comparison", "Region", "Haplotype_ID", "Allele_String"], kind="stable").reset_index(drop=True)


def compute_dataset_snp_frequencies(a1: np.ndarray, a2: np.ndarray, matched_meta: pd.DataFrame, dataset_label: str) -> pd.DataFrame:
    rows = []
    for idx, row in enumerate(matched_meta.itertuples(index=False)):
        keep = ~(np.isnan(a1[:, idx]) | np.isnan(a2[:, idx]))
        dosage = a1[keep, idx] + a2[keep, idx]
        rows.append({"Dataset": dataset_label, "SNP_ID": row.rsID_clean, "Gene": row.Gene, "Callable_Samples": int(keep.sum()), "Carrier_Frequency": float(np.mean(dosage >= 1)) if dosage.size else np.nan, "ALT_Allele_Frequency": float(np.mean(dosage) / 2.0) if dosage.size else np.nan})
    return pd.DataFrame(rows)


def compute_region_reference_haplotypes(stage19, snp_meta: pd.DataFrame, ld_summary: pd.DataFrame, a1: np.ndarray, a2: np.ndarray, matched_meta: pd.DataFrame, dataset_label: str) -> tuple[pd.DataFrame, dict[str, str]]:
    matched_index = {str(rsid): idx for idx, rsid in enumerate(matched_meta["rsID_clean"])}
    notes: dict[str, str] = {}
    tables: list[pd.DataFrame] = []
    full_rsids = snp_meta.sort_values("POS_int")["rsID_clean"].astype(str).tolist()
    if all(rsid in matched_index for rsid in full_rsids):
        region_indices = [matched_index[rsid] for rsid in full_rsids]
        region_meta = matched_meta.iloc[region_indices].reset_index(drop=True)
        full_freq, callable_samples = stage19.haplotype_frequency_table(a1[:, region_indices], a2[:, region_indices], region_meta, population=dataset_label, region_name="Full_Region")
        if not full_freq.empty:
            full_freq["Callable_Samples"] = callable_samples
            tables.append(full_freq)
        else:
            notes["Full_Region"] = f"No callable {dataset_label} samples across the full region."
    else:
        notes["Full_Region"] = f"One or more full-region SNPs were missing in the {dataset_label} VCF."
    for row in ld_summary.itertuples(index=False):
        region = str(row.Region)
        region_rsids = [token.strip() for token in str(row.SNPs).split(",") if token.strip()]
        if not all(rsid in matched_index for rsid in region_rsids):
            notes[region] = f"One or more SNPs for {region} were missing in the {dataset_label} VCF."
            continue
        region_indices = [matched_index[rsid] for rsid in region_rsids]
        region_meta = matched_meta.iloc[region_indices].reset_index(drop=True)
        region_freq, callable_samples = stage19.haplotype_frequency_table(a1[:, region_indices], a2[:, region_indices], region_meta, population=dataset_label, region_name=region)
        if region_freq.empty:
            notes[region] = f"No callable {dataset_label} samples for {region}."
            continue
        region_freq["Callable_Samples"] = callable_samples
        tables.append(region_freq)
    return (pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(), notes)


def delta_direction(value: float | None) -> str:
    if value is None or pd.isna(value):
        return "Not assessed"
    if abs(float(value)) < 1e-9:
        return "Matched to reference"
    return "Higher in tumour" if float(value) > 0 else "Lower in tumour"


def classify_reference_alignment(study_delta: float | None) -> str:
    if study_delta is None or pd.isna(study_delta):
        return "Non-comparable: missing reference frequency"
    if abs(float(study_delta)) <= 0.10:
        return "Similar to 1000G reference"
    return "Divergent from 1000G reference"


def build_snp_comparison(comparison: str, study_snps: pd.DataFrame, reference_snps: pd.DataFrame) -> pd.DataFrame:
    ref = reference_snps.rename(columns={"Carrier_Frequency": "Reference_1000G_Female_Carrier_Frequency", "ALT_Allele_Frequency": "Reference_1000G_Female_ALT_Allele_Frequency", "Callable_Samples": "Reference_1000G_Female_Callable_Samples"})
    out = study_snps.merge(ref.drop(columns=["Dataset"]), on=["SNP_ID", "Gene"], how="left")
    out["Comparison"] = comparison
    out["Study_vs_1000G_Delta"] = out["Study_Tumour_Frequency"] - out["Reference_1000G_Female_Carrier_Frequency"]
    out["Study_vs_1000G_Direction"] = out["Study_vs_1000G_Delta"].map(delta_direction)
    out["Validation_Status"] = out["Study_vs_1000G_Delta"].map(classify_reference_alignment)
    out["CrossValidation_Notes"] = "Exact SNP identity matched across the study and female-only 1000 Genomes reference."
    return out


def build_haplotype_comparison(comparison: str, catalogue: pd.DataFrame, reference_haps: pd.DataFrame) -> pd.DataFrame:
    if catalogue.empty:
        return pd.DataFrame()
    selected = catalogue[catalogue["Significant_FDR"].eq(True) | catalogue["Significant_Fisher"].eq(True)].copy()
    if selected.empty:
        return pd.DataFrame()
    ref = reference_haps.rename(columns={"Frequency": "Reference_1000G_Female_Frequency", "Reference_Haplotype_ID": "Reference_1000G_Haplotype_ID", "Variant_Content": "Reference_1000G_Variant_Content", "Callable_Samples": "Reference_1000G_Callable_Samples"})
    out = selected.merge(ref[["Region", "Allele_String", "Reference_1000G_Female_Frequency", "Reference_1000G_Haplotype_ID", "Reference_1000G_Variant_Content", "Reference_1000G_Callable_Samples"]], on=["Region", "Allele_String"], how="left")
    out["Comparison"] = comparison
    out["Result_Table_Name"] = f"{comparison}_Haplotype_Comparison"
    out["Result_Definition_Label"] = out["Haplotype_Definition_Type"].fillna(out["Haplotype_Set"])
    out["Study_vs_1000G_Delta"] = out["Study_Tumour_Frequency"] - out["Reference_1000G_Female_Frequency"]
    out["Study_vs_1000G_Direction"] = out["Study_vs_1000G_Delta"].map(delta_direction)
    out["Validation_Status"] = out["Study_vs_1000G_Delta"].map(classify_reference_alignment)
    out["CrossValidation_Notes"] = "Exact haplotype matching uses the study-defined region plus the exact allele string against female-only 1000 Genomes."
    return out


def collect_status_counts(df: pd.DataFrame, signal_type: str, comparison: str) -> list[dict[str, object]]:
    if df.empty:
        return []
    return [{"Signal_Type": signal_type, "Comparison": comparison, "Validation_Status": status, "Count": int(count)} for status, count in df["Validation_Status"].value_counts(dropna=False).items()]


def compact_region_label(value: object) -> str:
    text = str(value).strip()
    if text == "Full_Region":
        return "Full_Region"
    if "Block" in text:
        return text.split("_")[-1]
    return text


def wrapped_label(value: object, width: int = 52) -> str:
    text = str(value).replace("_", " ").strip()
    return text if len(text) <= width else shorten(text, width=width, placeholder="...")


def save_figure(fig: plt.Figure, base_path: Path) -> list[Path]:
    outputs = [base_path.with_suffix(".png"), base_path.with_suffix(".pdf")]
    for output in outputs:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outputs


def draw_frequency_figure(df: pd.DataFrame, label_column: str, reference_column: str, study_column: str, title: str, subtitle: str, footnote: str, base_path: Path) -> list[Path]:
    if df.empty:
        return []
    plot_df = df.reset_index(drop=True).copy()
    height = max(4.5, 0.38 * len(plot_df) + 1.9)
    fig, ax = plt.subplots(figsize=(11.5, height))
    y_positions = np.arange(len(plot_df))
    for idx, row in plot_df.iterrows():
        values = [row.get(reference_column), row.get(study_column)]
        values = [float(value) for value in values if pd.notna(value)]
        if values:
            ax.hlines(y=idx, xmin=min(values), xmax=max(values), color="#bdbdbd", linewidth=2.0, zorder=1)
    for column, style_key in [(reference_column, "Reference"), (study_column, "Study")]:
        valid = plot_df[column].notna() if column in plot_df.columns else pd.Series(False, index=plot_df.index)
        if not valid.any():
            continue
        style = DATASET_STYLE[style_key]
        ax.scatter(plot_df.loc[valid, column], y_positions[valid.to_numpy()], s=80, marker=style["marker"], color=style["color"], edgecolor="white", linewidth=0.8, label=style["label"], zorder=3)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(plot_df[label_column].tolist(), fontsize=9)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("Frequency")
    ax.set_title(title, loc="left", fontsize=16, weight="bold")
    fig.text(0.125, 0.955, subtitle, ha="left", va="top", fontsize=10)
    if footnote:
        fig.text(0.99, 0.01, footnote, ha="right", va="bottom", fontsize=8, color="#555555")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", color="#dddddd", linewidth=0.8)
    ax.grid(axis="y", visible=False)
    fig.tight_layout(rect=(0.0, 0.05, 1.0, 0.93))
    return save_figure(fig, base_path)


def plot_snp_frequency_figure(comparison: str, snp_comparison: pd.DataFrame, figures_dir: Path) -> list[Path]:
    if snp_comparison.empty:
        return []
    plot_df = snp_comparison.copy()
    plot_df["Display_Label"] = plot_df.apply(lambda row: wrapped_label(f"{row['SNP_ID']} | {row.get('Gene', 'Unknown')}", 48), axis=1)
    plot_df = plot_df.sort_values(["Study_vs_1000G_Delta", "SNP_ID"], key=lambda s: np.abs(s) if s.name == "Study_vs_1000G_Delta" else s, ascending=[False, True]).reset_index(drop=True)
    return draw_frequency_figure(plot_df, "Display_Label", "Reference_1000G_Female_Carrier_Frequency", "Study_Tumour_Frequency", f"{comparison} SNP frequency cross-validation", "Study tumour versus female-only 1000 Genomes", "Exact SNP identities and frequencies are listed in the stage-24 SNP comparison tables.", figures_dir / f"24_{comparison}_SNP_Frequency_Comparison")


def plot_haplotype_frequency_figure(comparison: str, hap_comparison: pd.DataFrame, definition_type: str, figures_dir: Path) -> list[Path]:
    subset = hap_comparison[hap_comparison["Haplotype_Definition_Type"].astype(str).eq(definition_type)].copy()
    if subset.empty:
        return []
    subset["Display_Label"] = subset.apply(lambda row: wrapped_label(f"{compact_region_label(row['Region'])} | {row['Haplotype_ID']} | {row['Allele_String']}", 52), axis=1)
    subset = subset.sort_values(["Study_vs_1000G_Delta", "Region", "Haplotype_ID"], key=lambda s: np.abs(s) if s.name == "Study_vs_1000G_Delta" else s, ascending=[False, True, True]).reset_index(drop=True)
    pretty_type = "Full-region" if definition_type == "Full_Region" else "LD-block"
    if definition_type == "Full_Region":
        footnote = f"Exact SNP backbone: 24_{comparison}_Full_Region_Haplotype_Definitions.tsv"
        file_stub = f"24_{comparison}_Haplotype_Frequency_Full_Region"
    else:
        footnote = f"Exact LD-block SNPs: 24_{comparison}_LD_Block_Definitions.tsv; tested haplotypes: 24_{comparison}_LD_Block_Haplotypes_Tested.tsv"
        file_stub = f"24_{comparison}_Haplotype_Frequency_LD_Block"
    return draw_frequency_figure(subset, "Display_Label", "Reference_1000G_Female_Frequency", "Study_Tumour_Frequency", f"{comparison} {pretty_type} haplotype cross-validation", f"{pretty_type} haplotypes only; full-region and LD-block outputs are plotted separately by design", footnote, figures_dir / file_stub)


def plot_validation_status_summary(status_rows: list[dict[str, object]], figures_dir: Path) -> list[Path]:
    if not status_rows:
        return []
    df = pd.DataFrame(status_rows)
    if df.empty:
        return []
    order = ["Similar to 1000G reference", "Divergent from 1000G reference"]
    extras = sorted([value for value in df["Validation_Status"].dropna().unique() if value not in order])
    df["Validation_Status"] = pd.Categorical(df["Validation_Status"], categories=order + extras, ordered=True)
    fig, ax = plt.subplots(figsize=(11, 6))
    plot_df = df.copy()
    plot_df["Panel"] = plot_df["Comparison"] + " | " + plot_df["Signal_Type"]
    sns.barplot(data=plot_df, x="Count", y="Panel", hue="Validation_Status", ax=ax)
    ax.set_title("External cross-validation status summary", loc="left", fontsize=16, weight="bold")
    ax.set_xlabel("Count")
    ax.set_ylabel("")
    ax.legend(frameon=False, title="Validation status", loc="lower right")
    fig.tight_layout()
    return save_figure(fig, figures_dir / "24_Validation_Status_Summary")


def main() -> None:
    args = parse_args()
    stage19 = load_stage19_module()
    out_dir = ensure_directory(args.out_dir)
    figures_dir = ensure_directory(out_dir / "figures")

    female_sample_path = Path(args.female_study_samples) if args.female_study_samples else None
    female_filter_active = bool(female_sample_path and female_sample_path.exists())
    study_metadata = load_study_sample_metadata(Path(args.sample_metadata), female_sample_path)

    panel_df = stage19.read_population_panel(Path(args.panel))
    reference_samples, ref_meta = stage19.get_population_samples(panel_df, get_header_samples(Path(args.vcf)), str(args.reference_population).upper(), sex_filter="female")
    reference_manifest = panel_df[panel_df["Sample"].isin(reference_samples)].copy()
    reference_manifest["Dataset"] = f"1000G_{str(args.reference_population).upper()}"
    reference_manifest["Comparison"] = "Breast;Endometrium"
    reference_manifest["Selected_for_CrossValidation"] = True
    reference_manifest["Age_Filter_Feasible"] = False
    reference_manifest["Age_Filter_Status"] = "1000 Genomes panel lacks age metadata in the local repo inputs"
    reference_manifest.to_csv(out_dir / "24_1000G_Female_Control_Manifest.tsv", sep="\t", index=False)

    all_snp_meta, reference_ld_summary, _, _ = load_haplotype_catalogues(Path(args.haplo_results), stage19, "Breast")
    a1_ref, a2_ref, matched_meta_ref, variant_qc_ref = stage19.load_reference_backbone(Path(args.vcf), all_snp_meta, reference_samples, args.contig)
    reference_snp_freqs = compute_dataset_snp_frequencies(a1_ref, a2_ref, matched_meta_ref, dataset_label=f"1000G_{str(args.reference_population).upper()}_FEMALE")
    reference_haps, reference_notes = compute_region_reference_haplotypes(stage19, all_snp_meta, reference_ld_summary, a1_ref, a2_ref, matched_meta_ref, dataset_label=f"1000G_{str(args.reference_population).upper()}_FEMALE")

    dataset_audit_rows = [{"Dataset": f"1000G_{str(args.reference_population).upper()}_FEMALE", "Comparison": "Breast;Endometrium", "Requested_Sex_Filter": ref_meta["Requested_Sex_Filter"], "Sex_Filter_Feasible": ref_meta["Sex_Metadata_Available"], "Sex_Filter_Applied": ref_meta["Sex_Filter_Applied"], "Sex_Filter_Status": ref_meta["Sex_Filter_Status"], "Age_Filter_Requested": False, "Age_Filter_Feasible": False, "Age_Filter_Applied": False, "Age_Filter_Status": "The local 1000 Genomes panel does not contain age metadata", "Selected_Samples": len(reference_samples), "Notes": "Healthy external reference for both tumour-cohort cross-validation branches."}]
    limitation_rows = [
        {"Category": "1000G sex filtering", "Comparison": "Breast;Endometrium", "Status": "Implemented", "Details": f"Female-only filtering used the local panel metadata for {str(args.reference_population).upper()}."},
        {"Category": "1000G age filtering", "Comparison": "Breast;Endometrium", "Status": "Not feasible", "Details": "The local 1000 Genomes panel exposes sex but not age metadata."},
        {"Category": "TCGA branch", "Comparison": "Breast;Endometrium", "Status": "Removed", "Details": "The TCGA external-comparison branch was removed because the required germline genotype files were controlled-access and not available for this project."},
    ]
    status_rows: list[dict[str, object]] = []
    workbook_tables = {"Study_Breast_Manifest": pd.DataFrame(), "Study_Endometrium_Manifest": pd.DataFrame(), "Breast_SNP_Comparison": pd.DataFrame(), "Endometrium_SNP_Comparison": pd.DataFrame(), "Breast_Haplotype_Comparison": pd.DataFrame(), "Endometrium_Haplotype_Comparison": pd.DataFrame(), "Breast_Full_Region_Definitions": pd.DataFrame(), "Endometrium_Full_Region_Definitions": pd.DataFrame(), "Breast_LD_Block_Definitions": pd.DataFrame(), "Endometrium_LD_Block_Definitions": pd.DataFrame(), "Breast_LD_Block_Haplotypes_Tested": pd.DataFrame(), "Endometrium_LD_Block_Haplotypes_Tested": pd.DataFrame(), "Result_Definition_Guide": pd.DataFrame(), "Figure_Guide": pd.DataFrame()}
    combined_snp: list[pd.DataFrame] = []
    combined_hap: list[pd.DataFrame] = []
    combined_full_defs: list[pd.DataFrame] = []
    combined_ld_block_defs: list[pd.DataFrame] = []
    combined_ld_block_haps: list[pd.DataFrame] = []
    figure_rows: list[dict[str, object]] = []

    for comparison in ["Breast", "Endometrium"]:
        config = COHORT_CONFIG[comparison]
        study_manifest, study_summary = select_study_manifest(study_metadata, comparison, female_filter_active)
        workbook_tables[f"Study_{comparison}_Manifest"] = study_manifest
        study_manifest.to_csv(out_dir / config["study_manifest_file"], sep="\t", index=False)
        dataset_audit_rows.append(study_summary)

        study_snps = load_significant_snps(Path(args.snp_results), comparison)
        snp_meta, ld_summary, full_catalogue, ld_catalogue = load_haplotype_catalogues(Path(args.haplo_results), stage19, comparison)
        full_definition_table = build_full_region_definition_table(full_catalogue)
        ld_block_definition_table = build_ld_block_definition_table(ld_catalogue)
        ld_block_haplotype_table = build_ld_block_haplotype_table(ld_catalogue)
        workbook_tables[f"{comparison}_Full_Region_Definitions"] = full_definition_table
        workbook_tables[f"{comparison}_LD_Block_Definitions"] = ld_block_definition_table
        workbook_tables[f"{comparison}_LD_Block_Haplotypes_Tested"] = ld_block_haplotype_table
        full_definition_table.to_csv(out_dir / f"24_{comparison}_Full_Region_Haplotype_Definitions.tsv", sep="\t", index=False)
        ld_block_definition_table.to_csv(out_dir / f"24_{comparison}_LD_Block_Definitions.tsv", sep="\t", index=False)
        ld_block_haplotype_table.to_csv(out_dir / f"24_{comparison}_LD_Block_Haplotypes_Tested.tsv", sep="\t", index=False)
        ld_block_haplotype_table.to_csv(out_dir / f"24_{comparison}_LD_Block_Haplotype_Definitions.tsv", sep="\t", index=False)
        combined_full_defs.append(full_definition_table)
        combined_ld_block_defs.append(ld_block_definition_table)
        combined_ld_block_haps.append(ld_block_haplotype_table)

        snp_comparison = build_snp_comparison(comparison, study_snps, reference_snp_freqs)
        hap_comparison = pd.concat([build_haplotype_comparison(comparison, full_catalogue, reference_haps), build_haplotype_comparison(comparison, ld_catalogue, reference_haps)], ignore_index=True)
        snp_comparison.to_csv(out_dir / f"24_{comparison}_SNP_Frequency_Comparison.tsv", sep="\t", index=False)
        hap_comparison.to_csv(out_dir / f"24_{comparison}_Haplotype_Frequency_Comparison.tsv", sep="\t", index=False)
        workbook_tables[f"{comparison}_SNP_Comparison"] = snp_comparison
        workbook_tables[f"{comparison}_Haplotype_Comparison"] = hap_comparison
        combined_snp.append(snp_comparison)
        combined_hap.append(hap_comparison)
        status_rows.extend(collect_status_counts(snp_comparison, "SNP", comparison))
        status_rows.extend(collect_status_counts(hap_comparison, "Haplotype", comparison))

        for figure_path in plot_snp_frequency_figure(comparison, snp_comparison, figures_dir):
            figure_rows.append({"Figure_File": figure_path.name, "Comparison": comparison, "Signal_Type": "SNP", "Haplotype_Definition_Used": pd.NA, "Description": "Carrier-frequency comparison for significant study SNPs across the study tumour cohort and female-only 1000 Genomes."})
        for definition_type in ["Full_Region", "LD_Block"]:
            for figure_path in plot_haplotype_frequency_figure(comparison, hap_comparison, definition_type, figures_dir):
                figure_rows.append({"Figure_File": figure_path.name, "Comparison": comparison, "Signal_Type": "Haplotype", "Haplotype_Definition_Used": definition_type, "Description": "Frequency comparison for significant study haplotypes against female-only 1000 Genomes. Full-region and LD-block figures are separated explicitly by haplotype definition."})

    combined_snp_df = pd.concat(combined_snp, ignore_index=True) if combined_snp else pd.DataFrame()
    combined_hap_df = pd.concat(combined_hap, ignore_index=True) if combined_hap else pd.DataFrame()
    combined_full_defs_df = pd.concat(combined_full_defs, ignore_index=True) if combined_full_defs else pd.DataFrame()
    combined_ld_block_defs_df = pd.concat(combined_ld_block_defs, ignore_index=True) if combined_ld_block_defs else pd.DataFrame()
    combined_ld_block_haps_df = pd.concat(combined_ld_block_haps, ignore_index=True) if combined_ld_block_haps else pd.DataFrame()
    combined_snp_df.to_csv(out_dir / "24_SNP_Frequency_Comparison.tsv", sep="\t", index=False)
    combined_hap_df.to_csv(out_dir / "24_Haplotype_Frequency_Comparison.tsv", sep="\t", index=False)
    combined_full_defs_df.to_csv(out_dir / "24_Full_Region_Haplotype_Definitions.tsv", sep="\t", index=False)
    combined_ld_block_defs_df.to_csv(out_dir / "24_LD_Block_Definitions.tsv", sep="\t", index=False)
    combined_ld_block_haps_df.to_csv(out_dir / "24_LD_Block_Haplotypes_Tested.tsv", sep="\t", index=False)
    combined_ld_block_haps_df.to_csv(out_dir / "24_LD_Block_Haplotype_Definitions.tsv", sep="\t", index=False)

    for figure_path in plot_validation_status_summary(status_rows, figures_dir):
        figure_rows.append({"Figure_File": figure_path.name, "Comparison": "Breast;Endometrium", "Signal_Type": "Summary", "Haplotype_Definition_Used": "Mixed", "Description": "Validation-status summary across cohorts and signal types."})

    workbook_tables["Figure_Guide"] = pd.DataFrame(figure_rows)
    workbook_tables["Figure_Guide"].to_csv(out_dir / "24_Figure_Guide.tsv", sep="\t", index=False)
    workbook_tables["Result_Definition_Guide"] = pd.DataFrame([
        {"Output_Name": "Breast_Full_Defs", "Comparison": "Breast", "Signal_Type": "Haplotype", "Haplotype_Definition_Used": "Full_Region", "Description": "All full-region haplotypes and the exact SNP backbone used to define them."},
        {"Output_Name": "Endometrium_Full_Defs", "Comparison": "Endometrium", "Signal_Type": "Haplotype", "Haplotype_Definition_Used": "Full_Region", "Description": "All full-region haplotypes and the exact SNP backbone used to define them."},
        {"Output_Name": "Breast_LD_Defs", "Comparison": "Breast", "Signal_Type": "Haplotype", "Haplotype_Definition_Used": "LD_Block", "Description": "LD blocks, their defining SNPs, and the r2 threshold used to define them."},
        {"Output_Name": "Endometrium_LD_Defs", "Comparison": "Endometrium", "Signal_Type": "Haplotype", "Haplotype_Definition_Used": "LD_Block", "Description": "LD blocks, their defining SNPs, and the r2 threshold used to define them."},
        {"Output_Name": "Breast_LD_Haps", "Comparison": "Breast", "Signal_Type": "Haplotype", "Haplotype_Definition_Used": "LD_Block", "Description": "All LD-block haplotypes tested within each breast LD block."},
        {"Output_Name": "Endometrium_LD_Haps", "Comparison": "Endometrium", "Signal_Type": "Haplotype", "Haplotype_Definition_Used": "LD_Block", "Description": "All LD-block haplotypes tested within each endometrial LD block."},
        {"Output_Name": "Breast_Hap_Compare", "Comparison": "Breast", "Signal_Type": "Haplotype", "Haplotype_Definition_Used": "Both", "Description": "Cross-validation results for significant haplotypes against female-only 1000 Genomes. Use the Haplotype_Definition_Type and Defining_SNPs columns to distinguish Full_Region from LD_Block results."},
        {"Output_Name": "Endometrium_Hap_Compare", "Comparison": "Endometrium", "Signal_Type": "Haplotype", "Haplotype_Definition_Used": "Both", "Description": "Cross-validation results for significant haplotypes against female-only 1000 Genomes. Use the Haplotype_Definition_Type and Defining_SNPs columns to distinguish Full_Region from LD_Block results."},
        {"Output_Name": "Haplotype_Comparison_All", "Comparison": "Breast;Endometrium", "Signal_Type": "Haplotype", "Haplotype_Definition_Used": "Both", "Description": "Combined significant-haplotype comparison table across cohorts, with explicit definition labels carried per row."},
    ])

    readme = pd.DataFrame({"Notes": [
        "Purpose: secondary cross-validation of breast and endometrial tumour SNP and haplotype findings, not primary discovery.",
        f"Healthy external reference: female-only 1000 Genomes {str(args.reference_population).upper()} for both tumour cohorts.",
        "No TCGA branch is retained in the active stage-24 workflow because the relevant TCGA germline genotype files were controlled-access and were not available for this project.",
        "SNP backbone: established/common SNPs from stages 11-12.",
        "Haplotype definitions used in this validation step: both full-region haplotypes and LD-block haplotypes are analysed.",
        "Breast_Full_Defs and Endometrium_Full_Defs list all full-region haplotypes together with the exact SNP backbone used to define them.",
        "Breast_LD_Defs and Endometrium_LD_Defs list the exact SNP composition of each LD block at r2 = 0.80.",
        "Breast_LD_Haps and Endometrium_LD_Haps list the exact haplotypes tested within each LD block.",
        "Breast_Hap_Compare, Endometrium_Hap_Compare, and Haplotype_Comparison_All can contain both full-region and LD-block results; use Haplotype_Definition_Type and Defining_SNPs to distinguish them row by row.",
        "Stage-24 figures are written to the figures directory. Haplotype figures are separated into Full_Region and LD_Block panels by filename and in the Figure_Guide sheet.",
        "Validation-status labels: Similar to 1000G reference, Divergent from 1000G reference, or Non-comparable.",
    ]})

    with pd.ExcelWriter(out_dir / "24_External_Cross_Validation.xlsx", engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        pd.DataFrame(dataset_audit_rows).to_excel(writer, sheet_name="Dataset_Audit", index=False)
        reference_manifest.to_excel(writer, sheet_name="1000G_Female_Manifest", index=False)
        reference_snp_freqs.to_excel(writer, sheet_name="1000G_SNP_Frequencies", index=False)
        reference_haps.to_excel(writer, sheet_name="1000G_Haplotype_Freqs", index=False)
        variant_qc_ref.to_excel(writer, sheet_name="1000G_Variant_QC", index=False)
        workbook_tables["Study_Breast_Manifest"].to_excel(writer, sheet_name="Study_Breast_Manifest", index=False)
        workbook_tables["Study_Endometrium_Manifest"].to_excel(writer, sheet_name="Study_Endometrium_Manifest", index=False)
        workbook_tables["Breast_Full_Region_Definitions"].to_excel(writer, sheet_name="Breast_Full_Defs", index=False)
        workbook_tables["Endometrium_Full_Region_Definitions"].to_excel(writer, sheet_name="Endometrium_Full_Defs", index=False)
        workbook_tables["Breast_LD_Block_Definitions"].to_excel(writer, sheet_name="Breast_LD_Defs", index=False)
        workbook_tables["Endometrium_LD_Block_Definitions"].to_excel(writer, sheet_name="Endometrium_LD_Defs", index=False)
        workbook_tables["Breast_LD_Block_Haplotypes_Tested"].to_excel(writer, sheet_name="Breast_LD_Haps", index=False)
        workbook_tables["Endometrium_LD_Block_Haplotypes_Tested"].to_excel(writer, sheet_name="Endometrium_LD_Haps", index=False)
        workbook_tables["Breast_SNP_Comparison"].to_excel(writer, sheet_name="Breast_SNP_Compare", index=False)
        workbook_tables["Endometrium_SNP_Comparison"].to_excel(writer, sheet_name="Endometrium_SNP_Compare", index=False)
        workbook_tables["Breast_Haplotype_Comparison"].to_excel(writer, sheet_name="Breast_Hap_Compare", index=False)
        workbook_tables["Endometrium_Haplotype_Comparison"].to_excel(writer, sheet_name="Endometrium_Hap_Compare", index=False)
        combined_snp_df.to_excel(writer, sheet_name="SNP_Comparison_All", index=False)
        combined_hap_df.to_excel(writer, sheet_name="Haplotype_Comparison_All", index=False)
        workbook_tables["Result_Definition_Guide"].to_excel(writer, sheet_name="Result_Def_Guide", index=False)
        workbook_tables["Figure_Guide"].to_excel(writer, sheet_name="Figure_Guide", index=False)
        pd.DataFrame(status_rows).to_excel(writer, sheet_name="Validation_Summary", index=False)
        pd.DataFrame(limitation_rows).to_excel(writer, sheet_name="Limitations", index=False)
        if reference_notes:
            pd.DataFrame([{"Region": region, "Note": note} for region, note in sorted(reference_notes.items())]).to_excel(writer, sheet_name="1000G_Hap_Notes", index=False)

    print("External cross-validation complete")
    print(f"Output directory: {out_dir}")
    print(f"Workbook: {out_dir / '24_External_Cross_Validation.xlsx'}")
    print(f"1000 Genomes female reference population: {str(args.reference_population).upper()} ({len(reference_samples)} samples)")
    print(f"Figures directory: {figures_dir}")


if __name__ == "__main__":
    main()
