"""
Script 11: Common SNP Identification and Frequency Benchmarking
===============================================================

Purpose
-------
Identify recurrent polymorphic variants in the annotated report and summarise
how frequently they are observed in each cohort/tissue group.

Definition of SNP used here
---------------------------
A variant is treated as an established SNP if either the gnomAD exome NFE or
genome NFE allele frequency exceeds the configured threshold.

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

from pipeline_utils import (
    attach_amplicon_warning_columns,
    build_amplicon_warning_lookup,
    build_genomic_variant_id_series,
    combine_gnomad_nfe,
    common_nfe_variant_mask,
    ensure_directory,
    get_paths,
    get_thresholds,
    standardize_cohort_labels,
    standardize_tissue_labels,
)
from pipeline_validation import print_validation_summary, validate_file_exists, validate_nonempty, validate_percentage_columns, validate_required_columns

# Shared configuration ensures that SNP reporting uses the same canonical input
# table and threshold values as the downstream association scripts.
PATHS = get_paths()
THRESHOLDS = get_thresholds()
input_file = str(PATHS["annotated_report"])
output_dir = str(ensure_directory(PATHS["common_snps_dir"]))


def identify_snps_pipeline():
    """Generate a study-wide summary of common SNP carrier frequencies."""
    validate_file_exists(input_file, "Script 11 input")

    print(">>> Script 11: Identifying SNPs and generating frequency reports...")
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")

    required_cols = [
        "Variant_Key", "CHROM", "POS", "REF", "ALT", "Cohort", "Tissue",
        "Sample", "Consequence", "IMPACT", "GT", "gnomADe_NFE_AF", "gnomADg_NFE_AF"
    ]
    validate_required_columns(df, required_cols, "Script 11 Biological_Annotations")

    # Create a single harmonised gnomAD frequency field before thresholding so
    # that exome and genome sources are treated consistently.
    df = combine_gnomad_nfe(df)
    df_snps_all = df[common_nfe_variant_mask(df, THRESHOLDS["min_nfe_af"])].copy()
    validate_nonempty(df_snps_all, "Script 11 common SNP subset")

    gt_map = {"0/0": "WT", "0/1": "Het", "1/0": "Het", "1/1": "Hom"}
    df_snps_all["Genotype_Category"] = df_snps_all["GT"].astype(str).str.strip().map(gt_map)
    excluded_reference_rows = int((df_snps_all["Genotype_Category"] == "WT").sum())
    excluded_unknown_rows = int(df_snps_all["Genotype_Category"].isna().sum())
    df_snps_all["Cohort"] = standardize_cohort_labels(df_snps_all["Cohort"])
    df_snps_all["Tissue"] = standardize_tissue_labels(df_snps_all["Tissue"])
    df_snps = df_snps_all[df_snps_all["Genotype_Category"].isin(["Het", "Hom"])].copy()
    validate_nonempty(df_snps, "Script 11 common SNP carrier subset")

    # Use genomic variant identity so non-coding common SNPs are preserved in
    # the unique set instead of disappearing when no protein label is available.
    df_snps["Variant_ID"] = build_genomic_variant_id_series(
        variant_key_series=df_snps["Variant_Key"],
        chrom_series=df_snps["CHROM"],
        pos_series=df_snps["POS"],
        ref_series=df_snps["REF"],
        alt_series=df_snps["ALT"],
    )

    # Extract rsID from Existing_variation so that script 15 (haplotype analysis)
    # can match whitelist entries against rsID_clean from the annotated variants
    # workbook. Existing_variation may contain compound strings like
    # "rs11078928&COSV58780118" so we extract only the leading rs accession.
    if "Existing_variation" in df_snps.columns:
        df_snps["rsID"] = df_snps["Existing_variation"].astype(str).str.extract(r"(rs\d+)", expand=False)
    else:
        df_snps["rsID"] = pd.NA

    df_snps["SYMBOL"] = df_snps["SYMBOL"].fillna("Intergenic")
    df_snps["Consequence"] = df_snps["Consequence"].fillna("Unknown")
    df_snps["IMPACT"] = df_snps["IMPACT"].fillna("Unknown")
    warning_lookup = build_amplicon_warning_lookup(df_snps, variant_col="Variant_ID", chrom_col="CHROM", pos_col="POS")
    print_validation_summary(df_snps, "Sample", "Script 11 SNP carrier subset", ["Cohort", "Tissue"])
    print(f"    Excluded GT 0/0 rows from carrier counts: {excluded_reference_rows}")
    print(f"    Excluded rows with unusable genotype labels: {excluded_unknown_rows}")

    # Frequency is defined on unique carriers rather than raw rows so that one
    # sample cannot contribute multiple times through transcript-level repeats.
    # Denominators still come from the full common-SNP subset before removing WT rows.
    total_samples_dict = df_snps_all.groupby(["Cohort", "Tissue"])["Sample"].nunique().to_dict()

    # Include rsID in the groupby so it is carried through to the pivot table.
    summary = df_snps.groupby([
        "Variant_ID", "rsID", "Cohort", "Tissue", "SYMBOL",
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
            "Variant_ID", "rsID", "SYMBOL", "Consequence", "IMPACT",
            "gnomAD_NFE_AF_combined", "gnomAD_NFE_Source"
        ],
        columns=["Cohort", "Tissue"],
        values="Frequency_%",
        aggfunc="first"
    )
    master_pivot.columns = [f"{col[0]}_{col[1]}_Frequency_%" for col in master_pivot.columns.values]
    master_pivot = master_pivot.fillna(0).reset_index()
    master_pivot.rename(columns={"gnomAD_NFE_AF_combined": "gnomAD_NFE_AF"}, inplace=True)
    master_pivot = attach_amplicon_warning_columns(master_pivot, warning_lookup)

    cols = master_pivot.columns.tolist()
    # rsID is placed first so the whitelist workbook is easier to read
    # downstream and the genomic coordinate remains available as a fallback.
    preferred_order = [
        "rsID", "Variant_ID", "Coverage_Risk_Flag", "Coverage_Risk_Amplicon",
        "Coverage_Risk_Region", "Coverage_Risk_Note", "SYMBOL", "Consequence", "IMPACT",
        "gnomAD_NFE_AF", "gnomAD_NFE_Source"
    ]
    freq_cols = [c for c in cols if c not in preferred_order]
    master_pivot = master_pivot[preferred_order + freq_cols]

    n_with_rsid = master_pivot["rsID"].notna().sum()
    print(f"    rsID populated for {n_with_rsid} / {len(master_pivot)} SNPs in whitelist")

    output_excel = os.path.join(output_dir, "GSDMB_Common_SNP_Frequency_Summary.xlsx")
    group_sample_counts = (
        df_snps_all[["Cohort", "Tissue", "Sample"]]
        .drop_duplicates()
        .groupby(["Cohort", "Tissue"], as_index=False)
        .agg(Total_Group_Samples=("Sample", "nunique"))
        .sort_values(["Cohort", "Tissue"], kind="stable")
        .reset_index(drop=True)
    )
    summary = summary.sort_values(
        ["Variant_ID", "Cohort", "Tissue"],
        kind="stable",
    ).reset_index(drop=True)

    workbook_readme = pd.DataFrame(
        [
            {
                "Sheet": "Common_SNP_Summary",
                "Purpose": "Wide thesis-facing summary of common SNP carrier frequencies by cohort and tissue group.",
                "Notes": "This is the canonical whitelist sheet used by downstream scripts.",
            },
            {
                "Sheet": "README",
                "Purpose": "Workbook guide and SNP definition notes.",
                "Notes": "A common SNP is defined here by combined gnomAD NFE frequency thresholding.",
            },
            {
                "Sheet": "Carrier_Frequency_Long",
                "Purpose": "Long-form per-group carrier counts and carrier frequencies.",
                "Notes": "Useful for checking exact cohort and tissue frequencies without the wide pivot.",
            },
            {
                "Sheet": "Group_Sample_Counts",
                "Purpose": "Group denominators used for carrier-frequency calculations.",
                "Notes": "These totals come from the common-SNP-ready annotated input before removing WT rows.",
            },
        ]
    )

    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
        # Keep the canonical summary as the first sheet so older helper scripts
        # that read the workbook without an explicit sheet name remain compatible.
        master_pivot.to_excel(writer, sheet_name="Common_SNP_Summary", index=False)
        workbook_readme.to_excel(writer, sheet_name="README", index=False)
        summary.to_excel(writer, sheet_name="Carrier_Frequency_Long", index=False)
        group_sample_counts.to_excel(writer, sheet_name="Group_Sample_Counts", index=False)

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
    plt.title("Common SNPs per gene", fontsize=14, fontweight="bold")
    plt.ylabel("Unique SNP count")
    plt.xlabel("Gene")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Common_SNPs_Per_Gene.png"), dpi=300)
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
        plt.title(f"Study versus gnomAD NFE reference: {cohort}", fontsize=15, fontweight="bold")
        plt.xlabel("gnomAD NFE allele frequency")
        plt.ylabel("Study carrier frequency (%)")
        plt.xlim(0, 1)
        plt.ylim(0, 100)
        legend = plt.legend(frameon=False, title="Tissue / reference source")
        for text in legend.get_texts():
            text.set_text(
                text.get_text()
                .replace("gnomAD_NFE_Source", "Reference source")
                .replace("Exome_NFE", "gnomAD exome NFE")
                .replace("Genome_NFE", "gnomAD genome NFE")
                .replace("Tumour", "Tumour")
                .replace("Healthy", "Healthy")
            )
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"Common_SNP_Frequency_vs_gnomAD_{cohort}.png"), dpi=300)
        plt.close()

    # This final panel documents how often the retained SNP definition was driven
    # by exome versus genome reference frequencies.
    source_counts = master_pivot.groupby("gnomAD_NFE_Source")["Variant_ID"].count().reset_index()
    source_counts["Reference_Source_Display"] = source_counts["gnomAD_NFE_Source"].replace({
        "Exome_NFE": "gnomAD exome NFE",
        "Genome_NFE": "gnomAD genome NFE",
    })
    plt.figure(figsize=(8, 5))
    sns.barplot(
        data=source_counts,
        x="Reference_Source_Display",
        y="Variant_ID",
        hue="Reference_Source_Display",
        palette="mako",
        legend=False,
    )
    plt.title("gnomAD NFE reference source counts", fontsize=13, fontweight="bold")
    plt.ylabel("Number of SNPs")
    plt.xlabel("Reference allele-frequency source")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Common_SNP_Annotation_Source_Counts.png"), dpi=300)
    plt.close()

    print(f"    Results saved to: {output_excel}")


if __name__ == "__main__":
    identify_snps_pipeline()
