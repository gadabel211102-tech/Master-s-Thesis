"""
Script 11: Common SNP Identification and Frequency Benchmarking
===============================================================

Purpose
-------
Identify recurrent polymorphic variants in the annotated report and summarise
how frequently they are observed in each cohort/tissue group.

Definition of SNP used here
---------------------------
A variant is treated as an established SNP if its combined non-Finnish European
(NFE) allele frequency in gnomAD exceeds the configured threshold. Exome NFE is
used when present, with genome NFE as fallback.

Outputs
-------
- Excel summary table of per-group carrier frequencies.
- Gene-level and source-level summary figures.
- Benchmark scatter plots comparing study frequencies against gnomAD NFE AF.
"""

import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from pipeline_utils import build_variant_id_series, combine_gnomad_nfe, get_paths, get_thresholds, standardize_cohort_labels, standardize_tissue_labels
from pipeline_validation import print_validation_summary, validate_file_exists, validate_nonempty, validate_percentage_columns, validate_required_columns

# Shared configuration ensures that SNP reporting uses the same canonical input
# table and threshold values as the downstream association scripts.
PATHS = get_paths()
THRESHOLDS = get_thresholds()
input_file = str(PATHS["annotated_report"])
output_dir = str(PATHS["results_dir"])


def identify_snps_pipeline():
    """Generate a study-wide summary of common SNP carrier frequencies."""
    validate_file_exists(input_file, "Script 11 input")

    print(">>> Script 11: Identifying SNPs and generating frequency reports...")
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")

    required_cols = [
        "Existing_variation", "SYMBOL", "HGVSp", "Cohort", "Tissue",
        "Sample", "Consequence", "IMPACT", "gnomADe_NFE_AF", "gnomADg_NFE_AF"
    ]
    validate_required_columns(df, required_cols, "Script 11 Biological_Annotations")

    # Create a single harmonised gnomAD frequency field before thresholding so
    # that exome and genome sources are treated consistently.
    df = combine_gnomad_nfe(df)
    df_snps = df[df["gnomAD_NFE_AF_combined"] > THRESHOLDS["min_nfe_af"]].copy()
    validate_nonempty(df_snps, "Script 11 common SNP subset")

    # Prefer rsIDs when available; otherwise fall back to a gene plus protein
    # notation so every retained SNP remains traceable in the output tables.
    df_snps["Variant_ID"] = build_variant_id_series(
        df_snps["Existing_variation"], df_snps["SYMBOL"], df_snps["HGVSp"]
    )
    df_snps["Cohort"] = standardize_cohort_labels(df_snps["Cohort"])
    df_snps["Tissue"] = standardize_tissue_labels(df_snps["Tissue"])
    print_validation_summary(df_snps, "Sample", "Script 11 SNP subset", ["Cohort", "Tissue"])

    # Frequency is defined on unique carriers rather than raw rows so that one
    # sample cannot contribute multiple times through transcript-level repeats.
    total_samples_dict = df_snps.groupby(["Cohort", "Tissue"])["Sample"].nunique().to_dict()
    summary = df_snps.groupby([
        "Variant_ID", "Cohort", "Tissue", "SYMBOL",
        "Consequence", "IMPACT", "gnomAD_NFE_AF_combined", "gnomAD_NFE_Source"
    ]).agg({"Sample": "nunique"}).reset_index()
    summary.rename(columns={"Sample": "Carrier_Count"}, inplace=True)
    summary["Total_Group_Samples"] = summary.apply(
        lambda x: total_samples_dict.get((x["Cohort"], x["Tissue"]), 0), axis=1
    )
    summary["Frequency_%"] = (summary["Carrier_Count"] / summary["Total_Group_Samples"]) * 100
    validate_percentage_columns(summary, ["Frequency_%"], "Script 11 summary")

    # Pivot to a wide table because supervisors and collaborators often find the
    # group-by-group percentage layout easier to review in Excel.
    master_pivot = summary.pivot_table(
        index=[
            "Variant_ID", "SYMBOL", "Consequence", "IMPACT",
            "gnomAD_NFE_AF_combined", "gnomAD_NFE_Source"
        ],
        columns=["Cohort", "Tissue"],
        values="Frequency_%",
        aggfunc="first"
    )
    master_pivot.columns = [f"{col[0]}_{col[1]}_Frequency_%" for col in master_pivot.columns.values]
    master_pivot = master_pivot.fillna(0).reset_index()
    master_pivot.rename(columns={"gnomAD_NFE_AF_combined": "gnomAD_NFE_AF"}, inplace=True)

    cols = master_pivot.columns.tolist()
    preferred_order = [
        "Variant_ID", "SYMBOL", "Consequence", "IMPACT",
        "gnomAD_NFE_AF", "gnomAD_NFE_Source"
    ]
    freq_cols = [c for c in cols if c not in preferred_order]
    master_pivot = master_pivot[preferred_order + freq_cols]

    output_excel = os.path.join(output_dir, "11_Master_Unique_SNP_Summary.xlsx")
    master_pivot.to_excel(output_excel, index=False)

    sns.set_style("whitegrid")

    # Gene-level counts summarise how much of the common-variant repertoire is
    # contributed by each locus represented in the panel.
    gene_counts = (
        master_pivot.groupby("SYMBOL")["Variant_ID"]
        .nunique()
        .sort_values(ascending=False)
        .reset_index()
    )
    plt.figure(figsize=(10, 6))
    sns.barplot(data=gene_counts, x="SYMBOL", y="Variant_ID", hue="SYMBOL", palette="viridis", legend=False)
    plt.title("Total Unique SNPs Identified per Gene", fontsize=14, fontweight="bold")
    plt.ylabel("Unique SNP Count")
    plt.xlabel("SYMBOL")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "11_BarPlot_SNPs_Per_Gene.png"), dpi=300)
    plt.close()

    # Benchmark plots place study frequencies against the population reference
    # frequency to highlight variants enriched or depleted in each cohort.
    for cohort in summary["Cohort"].dropna().unique():
        plt.figure(figsize=(10, 8))
        df_c = summary[summary["Cohort"] == cohort].copy()
        sns.scatterplot(
            data=df_c,
            x="gnomAD_NFE_AF_combined",
            y="Frequency_%",
            hue="Tissue",
            style="gnomAD_NFE_Source",
            s=130,
            palette={"Tumour": "#e74c3c", "Healthy": "#3498db"},
            alpha=0.75
        )
        plt.plot([0, 1], [0, 100], "--", color="grey", alpha=0.4, label="Reference line")
        plt.title(f"SNP Frequency Benchmarking: {cohort} Cohort", fontsize=15, fontweight="bold")
        plt.xlabel("gnomAD NFE AF")
        plt.ylabel("Carrier Frequency in Study (%)")
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"11_SNP_Benchmark_{cohort}.png"), dpi=300)
        plt.close()

    # This final panel documents how often the retained SNP definition was driven
    # by exome versus genome reference frequencies.
    source_counts = master_pivot.groupby("gnomAD_NFE_Source")["Variant_ID"].count().reset_index()
    plt.figure(figsize=(8, 5))
    sns.barplot(data=source_counts, x="gnomAD_NFE_Source", y="Variant_ID", hue="gnomAD_NFE_Source", palette="mako", legend=False)
    plt.title("SNP Annotation Source", fontsize=13, fontweight="bold")
    plt.ylabel("Number of SNPs")
    plt.xlabel("Source")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "11_SNP_Source_Counts.png"), dpi=300)
    plt.close()

    print(f"    Results saved to: {output_excel}")


if __name__ == "__main__":
    identify_snps_pipeline()
