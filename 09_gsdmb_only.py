"""
Script 09: GSDMB-Focused Landscape Mapping
==========================================

Purpose
-------
Generate split genomic landscapes restricted to GSDMB variants so the study
gene can be reviewed without mixing established common SNPs and lower-frequency
variants in the same main figure.

Input
-----
- ``Biological_Annotations`` sheet from the consolidated annotated report.

Output
------
- ``GSDMB_Target_Gene_Common_SNP_Landscape.png``
- ``GSDMB_Target_Gene_Lower_Frequency_Variant_Landscape.png``
- ``GSDMB_Target_Gene_Common_SNP_Genotype_Landscape.png``
- ``GSDMB_Target_Gene_Lower_Frequency_Variant_Genotype_Landscape.png``
- ``GSDMB_Target_Gene_Variant_Landscape_Split.xlsx``
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from figure_style import IMPACT_COLORS
from pipeline_utils import (
    attach_amplicon_warning_columns,
    build_amplicon_warning_lookup,
    build_variant_id_series,
    combine_gnomad_nfe,
    ensure_directory,
    extract_rsid,
    find_col,
    get_paths,
    get_thresholds,
    standardize_cohort_labels,
    standardize_tissue_labels,
)
from pipeline_validation import print_validation_summary, validate_file_exists, validate_required_columns

PATHS = get_paths()
THRESHOLDS = get_thresholds()
input_file = str(PATHS["annotated_report"])
output_dir = ensure_directory(PATHS.get("gsdmb_landscape_dir", PATHS["results_dir"] / "09_gsdmb_landscape"))
output_common_image = output_dir / "GSDMB_Target_Gene_Common_SNP_Landscape.png"
output_rare_image = output_dir / "GSDMB_Target_Gene_Lower_Frequency_Variant_Landscape.png"
output_common_genotype_image = output_dir / "GSDMB_Target_Gene_Common_SNP_Genotype_Landscape.png"
output_rare_genotype_image = output_dir / "GSDMB_Target_Gene_Lower_Frequency_Variant_Genotype_Landscape.png"
output_excel = output_dir / "GSDMB_Target_Gene_Variant_Landscape_Split.xlsx"
LABEL_MODE = "high_only"
GT_MAP = {"0/0": "WT", "0/1": "Het", "1/0": "Het", "1/1": "Hom"}
GENOTYPE_COLORS = {"WT": "#6C757D", "Het": "#1F77B4", "Hom": "#D62728"}
GENOTYPE_MARKERS = {"WT": "o", "Het": "s", "Hom": "^"}
GENOTYPE_OFFSETS = {"WT": -15000, "Het": 0, "Hom": 15000}


def label_impacts():
    return ["HIGH"] if LABEL_MODE == "high_only" else ["HIGH", "MODERATE"]


def annotate_panel_variants(ax, panel_df, pos_c, impact_c, y_col):
    offsets = [(6, 6), (6, -10), (-10, 8), (-12, -10), (10, 14), (-10, 14)]
    placed = []

    priority = {impact: rank for rank, impact in enumerate(label_impacts())}
    panel_df = panel_df.assign(_label_rank=panel_df[impact_c].map(priority).fillna(99))
    panel_df = panel_df.sort_values(by=["_label_rank", y_col, pos_c], ascending=[True, False, True])

    for _, row in panel_df.iterrows():
        x_val = float(row[pos_c])
        y_val = float(row[y_col])
        too_close = any(abs(x_val - px) < 20000 and abs(y_val - py) < 5 for px, py in placed)
        if too_close:
            continue
        ax.annotate(
            row["Variant_Label"],
            xy=(x_val, y_val),
            xytext=offsets[len(placed) % len(offsets)],
            textcoords="offset points",
            fontsize=8,
            fontweight="bold" if row[impact_c] == "HIGH" else "normal",
            color="black",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.75),
            arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5, alpha=0.5),
        )
        placed.append((x_val, y_val))


def normalise_gt(value):
    text = str(value).strip()
    return GT_MAP.get(text)


def build_sample_totals(df, cohort_col, tissue_col, sample_col):
    return df.groupby([cohort_col, tissue_col])[sample_col].nunique().reset_index(name="Total_Samples")


def prepare_annotations(df, cols):
    df = combine_gnomad_nfe(df)
    df[cols["cohort"]] = standardize_cohort_labels(df[cols["cohort"]])
    df[cols["tissue"]] = standardize_tissue_labels(df[cols["tissue"]])
    df[cols["impact"]] = df[cols["impact"]].astype(str).str.strip().str.upper()
    df["Variant_ID"] = build_variant_id_series(df["Existing_variation"], df[cols["symbol"]], df["HGVSp"])
    df["Variant_Label"] = df["Existing_variation"].apply(extract_rsid).fillna(df["Variant_ID"])
    df["Genotype_Category"] = df[cols["gt"]].apply(normalise_gt)
    df["Is_Alt_Carrier"] = df["Genotype_Category"].isin(["Het", "Hom"])
    df["Has_Usable_Genotype"] = df["Genotype_Category"].isin(["WT", "Het", "Hom"])
    df["Is_Common_SNP_By_NFE"] = df["gnomAD_NFE_AF_combined"].fillna(-1).gt(THRESHOLDS["min_nfe_af"])
    df["Variant_Class"] = df["Is_Common_SNP_By_NFE"].map({True: "Established_SNP", False: "Lower_Frequency_Mutation"})
    df["NFE_SNP_Threshold"] = THRESHOLDS["min_nfe_af"]
    warning_lookup = build_amplicon_warning_lookup(df, variant_col="Variant_ID", chrom_col="CHROM", pos_col=cols["pos"])
    sample_totals = build_sample_totals(df, cols["cohort"], cols["tissue"], cols["sample"])
    return df, sample_totals, warning_lookup


def build_variant_summary(df, cols, sample_totals, warning_lookup):
    cohort_col = cols["cohort"]
    tissue_col = cols["tissue"]
    sample_col = cols["sample"]
    impact_col = cols["impact"]
    pos_col = cols["pos"]

    carrier_df = df[df["Is_Alt_Carrier"]].copy()
    group_cols = [
        "Variant_ID", "CHROM", pos_col, "REF", "ALT", cols["symbol"],
        "Consequence", impact_col, cohort_col, tissue_col,
        "gnomAD_NFE_AF_combined", "gnomAD_NFE_Source",
        "Is_Common_SNP_By_NFE", "Variant_Class", "NFE_SNP_Threshold",
    ]
    variant_counts = (
        carrier_df.groupby(group_cols, dropna=False)[sample_col]
        .nunique().reset_index(name="Carrier_Count")
    )
    variant_counts = variant_counts.merge(sample_totals, on=[cohort_col, tissue_col], how="left")
    variant_counts["Carrier_Percentage"] = (variant_counts["Carrier_Count"] / variant_counts["Total_Samples"]) * 100
    variant_labels = (
        carrier_df.groupby(group_cols, dropna=False)["Variant_Label"]
        .agg(lambda values: "; ".join(sorted({str(v).strip() for v in values if str(v).strip() and str(v).strip() != "None"})) or None)
        .reset_index()
    )
    variant_counts = variant_counts.merge(variant_labels, on=group_cols, how="left")
    variant_counts["Labelled_On_Figure"] = variant_counts[impact_col].isin(label_impacts()) & variant_counts["Variant_Label"].notna()
    variant_counts = attach_amplicon_warning_columns(variant_counts, warning_lookup)
    return variant_counts.sort_values([cohort_col, tissue_col, pos_col, "Variant_ID"], kind="stable")


def classify_pattern(row):
    if row["Unknown_Genotype_Carriers"] > 0:
        return "Carrier_GT_missing_or_ambiguous"
    if row["Hom_Carriers"] == 0 and row["Het_Carriers"] > 0:
        return "Het_only"
    if row["Het_Carriers"] == 0 and row["Hom_Carriers"] > 0:
        return "Hom_only"
    if row["Het_Carriers"] > row["Hom_Carriers"]:
        return "Mixed_Het_dominant"
    if row["Hom_Carriers"] > row["Het_Carriers"]:
        return "Mixed_Hom_dominant"
    if row["Het_Carriers"] == row["Hom_Carriers"] and row["Het_Carriers"] > 0:
        return "Mixed_Balanced"
    return "No_nonWT_carriers"


def build_genotype_summary(df, variant_summary, cols, warning_lookup):
    cohort_col = cols["cohort"]
    tissue_col = cols["tissue"]
    sample_col = cols["sample"]
    pos_col = cols["pos"]

    summary_cols = [
        "Variant_ID", "CHROM", pos_col, "REF", "ALT", cols["symbol"],
        "Consequence", cols["impact"], cohort_col, tissue_col,
        "gnomAD_NFE_AF_combined", "gnomAD_NFE_Source",
        "Is_Common_SNP_By_NFE", "Variant_Class", "NFE_SNP_Threshold",
        "Variant_Label", "Labelled_On_Figure",
        "Coverage_Risk_Flag", "Coverage_Risk_Amplicon", "Coverage_Risk_Region", "Coverage_Risk_Note",
    ]
    variant_level = variant_summary[summary_cols + ["Carrier_Count", "Total_Samples"]].copy()
    sample_level = df[df["Is_Alt_Carrier"] | df["Genotype_Category"].isna()].copy()
    sample_level = sample_level.drop_duplicates(["Variant_ID", cohort_col, tissue_col, sample_col, "Genotype_Category"])
    group_keys = ["Variant_ID", cohort_col, tissue_col]

    het_counts = sample_level[sample_level["Genotype_Category"] == "Het"].groupby(group_keys)[sample_col].nunique().reset_index(name="Het_Carriers")
    hom_counts = sample_level[sample_level["Genotype_Category"] == "Hom"].groupby(group_keys)[sample_col].nunique().reset_index(name="Hom_Carriers")
    usable_counts = sample_level[sample_level["Genotype_Category"].isin(["Het", "Hom"])] .groupby(group_keys)[sample_col].nunique().reset_index(name="Usable_Genotype_Carriers")

    pattern_df = variant_level.merge(het_counts, on=group_keys, how="left")
    pattern_df = pattern_df.merge(hom_counts, on=group_keys, how="left")
    pattern_df = pattern_df.merge(usable_counts, on=group_keys, how="left")
    for col in ["Het_Carriers", "Hom_Carriers", "Usable_Genotype_Carriers"]:
        pattern_df[col] = pattern_df[col].fillna(0).astype(int)
    pattern_df["Unknown_Genotype_Carriers"] = (pattern_df["Carrier_Count"] - pattern_df["Usable_Genotype_Carriers"]).clip(lower=0)
    pattern_df["WT_Samples"] = (pattern_df["Total_Samples"] - pattern_df["Carrier_Count"]).clip(lower=0).astype(int)
    pattern_df["Het_Percentage"] = (pattern_df["Het_Carriers"] / pattern_df["Total_Samples"]) * 100
    pattern_df["Hom_Percentage"] = (pattern_df["Hom_Carriers"] / pattern_df["Total_Samples"]) * 100
    pattern_df["WT_Percentage"] = (pattern_df["WT_Samples"] / pattern_df["Total_Samples"]) * 100
    pattern_df["Carrier_Genotype_Pattern"] = pattern_df.apply(classify_pattern, axis=1)

    genotype_rows = []
    for genotype_label, count_col, pct_col in [("WT", "WT_Samples", "WT_Percentage"), ("Het", "Het_Carriers", "Het_Percentage"), ("Hom", "Hom_Carriers", "Hom_Percentage")]:
        sub = pattern_df[summary_cols + ["Total_Samples", count_col, pct_col, "Carrier_Genotype_Pattern", "Unknown_Genotype_Carriers"]].copy()
        sub["Genotype_Category"] = genotype_label
        sub["Genotype_Count"] = sub[count_col].astype(int)
        sub["Genotype_Percentage"] = sub[pct_col]
        genotype_rows.append(sub.drop(columns=[count_col, pct_col]))

    genotype_summary = pd.concat(genotype_rows, ignore_index=True)
    genotype_summary = attach_amplicon_warning_columns(genotype_summary, warning_lookup)
    pattern_df = attach_amplicon_warning_columns(pattern_df, warning_lookup)
    genotype_summary = genotype_summary.sort_values([cohort_col, tissue_col, pos_col, "Variant_ID", "Genotype_Category"], kind="stable")
    pattern_df = pattern_df.sort_values([cohort_col, tissue_col, pos_col, "Variant_ID"], kind="stable")
    return genotype_summary, pattern_df


def plot_variant_landscape(plot_df, cols, sample_size_lookup, title, output_image):
    if plot_df.empty:
        print(f"SKIP: no data for {output_image.name}")
        return

    pos_c = cols["pos"]
    imp_c = cols["impact"]
    coh_c = cols["cohort"]
    tis_c = cols["tissue"]

    plot_df = plot_df.copy()
    plot_df[imp_c] = pd.Categorical(plot_df[imp_c], categories=["HIGH", "MODERATE", "MODIFIER", "LOW"], ordered=True)
    sns.set_style("whitegrid")
    marker_map = {"HIGH": "X", "MODERATE": "o", "MODIFIER": "s", "LOW": "v"}

    g = sns.relplot(
        data=plot_df,
        x=pos_c,
        y="Carrier_Percentage",
        hue=imp_c,
        style=imp_c,
        row=coh_c,
        col=tis_c,
        kind="scatter",
        markers=marker_map,
        palette=IMPACT_COLORS,
        s=150,
        alpha=0.8,
        edgecolor="black",
        height=5,
        aspect=1.4,
        facet_kws={"sharex": True, "sharey": True},
    )

    labelled = plot_df[plot_df["Labelled_On_Figure"]].copy()
    for (cohort_label, tissue_label), panel_df in labelled.groupby([coh_c, tis_c], sort=False):
        ax = g.axes_dict.get((cohort_label, tissue_label))
        if ax is not None:
            annotate_panel_variants(ax, panel_df, pos_c, imp_c, "Carrier_Percentage")

    g.set_axis_labels("Chr17 Position (Mb)", "Carrier Frequency (%)")
    for (cohort_label, tissue_label), ax in g.axes_dict.items():
        n_samples = sample_size_lookup.get((cohort_label, tissue_label), 0)
        ax.set_title(f"{cohort_label} | {tissue_label} (n={n_samples})", fontweight="bold")
    for ax in g.axes.flat:
        ax.set_ylim(0, 100)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x/1e6:.2f}"))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y:.0f}%"))

    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(title, fontsize=16, fontweight="bold")
    plt.savefig(output_image, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    print(f"SUCCESS: saved {output_image}")


def plot_genotype_landscape(plot_df, cols, sample_size_lookup, title, output_image):
    if plot_df.empty:
        print(f"SKIP: no genotype data for {output_image.name}")
        return

    pos_c = cols["pos"]
    imp_c = cols["impact"]
    coh_c = cols["cohort"]
    tis_c = cols["tissue"]

    plot_df = plot_df[plot_df["Genotype_Category"].isin(["WT", "Het", "Hom"])].copy()
    if plot_df.empty:
        print(f"SKIP: no WT/Het/Hom rows for {output_image.name}")
        return

    plot_df["Plot_POS"] = plot_df[pos_c] + plot_df["Genotype_Category"].map(GENOTYPE_OFFSETS).fillna(0)
    sns.set_style("whitegrid")
    g = sns.relplot(
        data=plot_df,
        x="Plot_POS",
        y="Genotype_Percentage",
        hue="Genotype_Category",
        style="Genotype_Category",
        row=coh_c,
        col=tis_c,
        kind="scatter",
        markers=GENOTYPE_MARKERS,
        palette=GENOTYPE_COLORS,
        s=115,
        alpha=0.82,
        edgecolor="black",
        height=5,
        aspect=1.4,
        facet_kws={"sharex": True, "sharey": True},
    )

    labelled = (
        plot_df[
            plot_df["Labelled_On_Figure"]
            & plot_df["Genotype_Category"].isin(["Het", "Hom"])
            & plot_df["Genotype_Count"].gt(0)
        ]
        .sort_values("Genotype_Percentage", ascending=False)
        .drop_duplicates(["Variant_ID", coh_c, tis_c])
    )
    for (cohort_label, tissue_label), panel_df in labelled.groupby([coh_c, tis_c], sort=False):
        ax = g.axes_dict.get((cohort_label, tissue_label))
        if ax is not None:
            panel_df = panel_df.assign(**{pos_c: panel_df["Plot_POS"]})
            annotate_panel_variants(ax, panel_df, pos_c, imp_c, "Genotype_Percentage")

    g.set_axis_labels("Chr17 Position (Mb)", "WT / Het / Hom frequency (%)")
    for (cohort_label, tissue_label), ax in g.axes_dict.items():
        n_samples = sample_size_lookup.get((cohort_label, tissue_label), 0)
        ax.set_title(f"{cohort_label} | {tissue_label} (n={n_samples})", fontweight="bold")
    for ax in g.axes.flat:
        ax.set_ylim(0, 100)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x/1e6:.2f}"))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y:.0f}%"))

    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(title, fontsize=16, fontweight="bold")
    plt.savefig(output_image, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    print(f"SUCCESS: saved {output_image}")


def export_workbook(common_variants, rare_variants, common_genotypes, rare_genotypes, common_patterns, rare_patterns, sample_totals, excluded_reference_rows, excluded_unknown_rows):
    summary_df = pd.DataFrame([
        {"Output": "SNP definition", "Value": f"gnomAD NFE AF > {THRESHOLDS['min_nfe_af']:.2f}", "Notes": "Shared definition reused from Script 11 and the SNP workflow."},
        {"Output": "Lower-frequency definition", "Value": f"gnomAD NFE AF <= {THRESHOLDS['min_nfe_af']:.2f} or missing", "Notes": "Includes variants outside the established SNP set."},
        {"Output": "Reference GT rows excluded from carrier maps", "Value": int(excluded_reference_rows), "Notes": "GT 0/0 rows are not counted as variant carriers."},
        {"Output": "Rows with unusable genotype excluded from carrier maps", "Value": int(excluded_unknown_rows), "Notes": "Rows outside WT/Het/Hom are not forced into genotype plots."},
    ])
    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="summary", index=False)
        sample_totals.to_excel(writer, sheet_name="sample_totals", index=False)
        common_variants.to_excel(writer, sheet_name="common_snp_points", index=False)
        rare_variants.to_excel(writer, sheet_name="lower_freq_points", index=False)
        common_genotypes.to_excel(writer, sheet_name="common_snp_genotypes", index=False)
        rare_genotypes.to_excel(writer, sheet_name="lower_freq_genotypes", index=False)
        common_patterns.to_excel(writer, sheet_name="common_snp_patterns", index=False)
        rare_patterns.to_excel(writer, sheet_name="lower_freq_patterns", index=False)


def generate_gsdmb_maps():
    print("--- Generating split GSDMB-only landscape maps ---")

    validate_file_exists(input_file, "Script 09 input")
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")

    pos_c = find_col(df, "Pos")
    sym_c = find_col(df, "Symbol")
    imp_c = find_col(df, "Impact")
    coh_c = find_col(df, "Cohort")
    tis_c = find_col(df, "Tissue")
    sam_c = find_col(df, "Sample")
    gt_c = find_col(df, "GT")

    required_cols = [pos_c, sym_c, imp_c, coh_c, tis_c, sam_c, gt_c, "Existing_variation", "HGVSp", "CHROM", "REF", "ALT", "Consequence", "gnomADe_NFE_AF", "gnomADg_NFE_AF"]
    validate_required_columns(df, required_cols, "Script 09 Biological_Annotations")

    df_gsdmb = df[df[sym_c].astype(str).str.contains("GSDMB", case=False, na=False)].copy()
    print_validation_summary(df_gsdmb, sam_c, "Script 09 GSDMB subset", [coh_c, tis_c])
    if df_gsdmb.empty:
        print("Warning: No GSDMB variants found in the dataset!")
        return

    cols = {"cohort": coh_c, "tissue": tis_c, "sample": sam_c, "impact": imp_c, "symbol": sym_c, "pos": pos_c, "gt": gt_c}
    df_gsdmb, sample_totals, warning_lookup = prepare_annotations(df_gsdmb, cols)
    sample_size_lookup = sample_totals.set_index([coh_c, tis_c])["Total_Samples"].to_dict()

    variant_summary = build_variant_summary(df_gsdmb, cols, sample_totals, warning_lookup)
    genotype_summary, pattern_summary = build_genotype_summary(df_gsdmb, variant_summary, cols, warning_lookup)

    common_variants = variant_summary[variant_summary["Is_Common_SNP_By_NFE"]].copy()
    rare_variants = variant_summary[~variant_summary["Is_Common_SNP_By_NFE"]].copy()
    common_genotypes = genotype_summary[genotype_summary["Is_Common_SNP_By_NFE"]].copy()
    rare_genotypes = genotype_summary[~genotype_summary["Is_Common_SNP_By_NFE"]].copy()
    common_patterns = pattern_summary[pattern_summary["Is_Common_SNP_By_NFE"]].copy()
    rare_patterns = pattern_summary[~pattern_summary["Is_Common_SNP_By_NFE"]].copy()

    plot_variant_landscape(common_variants, cols, sample_size_lookup, "GSDMB Landscape: Established SNPs Only (>1% gnomAD NFE AF)", output_common_image)
    plot_variant_landscape(rare_variants, cols, sample_size_lookup, "GSDMB Landscape: Lower-Frequency Mutations Only (<=1% or missing gnomAD NFE AF)", output_rare_image)
    plot_genotype_landscape(common_genotypes, cols, sample_size_lookup, "GSDMB Genotype Landscape: Established SNPs (WT / Het / Hom)", output_common_genotype_image)
    plot_genotype_landscape(rare_genotypes, cols, sample_size_lookup, "GSDMB Genotype Landscape: Lower-Frequency Mutations (WT / Het / Hom)", output_rare_genotype_image)

    excluded_reference_rows = int((df_gsdmb["Genotype_Category"] == "WT").sum())
    excluded_unknown_rows = int(df_gsdmb["Genotype_Category"].isna().sum())
    export_workbook(common_variants, rare_variants, common_genotypes, rare_genotypes, common_patterns, rare_patterns, sample_totals, excluded_reference_rows, excluded_unknown_rows)
    print(f"SUCCESS: split GSDMB workbook saved to: {output_excel}")


if __name__ == "__main__":
    generate_gsdmb_maps()
