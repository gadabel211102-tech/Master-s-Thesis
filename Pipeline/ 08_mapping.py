"""
Script 08: Global Variant Landscape Mapping
===========================================

Purpose
-------
Visualise all detected variants across the targeted Chr17 panel as a faceted
landscape plot. Each point represents one unique variant observed within a
specific cohort and tissue stratum.

Methodological rationale
------------------------
The y-axis reports carrier percentage rather than raw carrier count. This makes
between-panel comparison defensible when cohort sizes differ, because each point
is scaled to the number of unique samples available in the same cohort/tissue
subset.

Input
-----
- ``Biological_Annotations`` sheet from the consolidated annotated report.

Output
------
- ``GSDMB_Global_Variant_Landscape.png`` in the configured results directory.
- ``GSDMB_Global_Variant_Landscape.xlsx`` containing the plotted variants plus
  gnomAD NFE SNP context.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pipeline_utils import (
    attach_amplicon_warning_columns,
    build_amplicon_warning_lookup,
    build_variant_id_series,
    combine_gnomad_nfe,
    compute_carrier_percentage,
    ensure_directory,
    extract_rsid,
    find_col,
    get_paths,
    get_thresholds,
    standardize_cohort_labels,
    standardize_tissue_labels,
)
from figure_style import QUALITATIVE_COLORBLIND_SEQUENCE
from pipeline_validation import print_validation_summary, validate_file_exists, validate_required_columns

# Centralised path configuration keeps this plotting script aligned with the
# rest of the pipeline and reduces the risk of stale hard-coded locations.
PATHS = get_paths()
THRESHOLDS = get_thresholds()
input_file = str(PATHS["annotated_report"])
output_dir = ensure_directory(PATHS.get("variant_landscape_dir", PATHS["results_dir"] / "08_variant_landscape"))
output_image = str(output_dir / "GSDMB_Global_Variant_Landscape.png")
output_excel = str(output_dir / "GSDMB_Global_Variant_Landscape.xlsx")
LABEL_MODE = "high_only"  # valid options: "high_only", "both"


def label_impacts():
    """Return the impact classes that should receive text labels."""
    return ["HIGH"] if LABEL_MODE == "high_only" else ["HIGH", "MODERATE"]


def annotate_panel_variants(ax, panel_df, pos_c, impact_c):
    """Place a small set of readable labels without letting them pile up."""
    offsets = [(6, 6), (6, -10), (-10, 8), (-12, -10), (10, 14), (-10, 14)]
    placed = []

    priority = {impact: rank for rank, impact in enumerate(label_impacts())}
    panel_df = panel_df.assign(_label_rank=panel_df[impact_c].map(priority).fillna(99))
    panel_df = panel_df.sort_values(
        by=["_label_rank", "Carrier_Percentage", pos_c],
        ascending=[True, False, True]
    )

    for _, row in panel_df.iterrows():
        x_val = float(row[pos_c])
        y_val = float(row["Carrier_Percentage"])
        too_close = any(
            abs(x_val - placed_x) < 30000 and abs(y_val - placed_y) < 5
            for placed_x, placed_y in placed
        )
        if too_close:
            continue

        offset = offsets[len(placed) % len(offsets)]
        ax.annotate(
            row["Variant_Label"],
            xy=(x_val, y_val),
            xytext=offset,
            textcoords='offset points',
            fontsize=8,
            fontweight='bold' if row[impact_c] == "HIGH" else "normal",
            color='black',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.75),
            arrowprops=dict(arrowstyle='-', color='0.5', lw=0.5, alpha=0.5)
        )
        placed.append((x_val, y_val))


def build_landscape_summary(df, cols):
    """Build the point-level summary used both for plotting and workbook export."""
    cohort_col = cols["cohort"]
    tissue_col = cols["tissue"]
    sample_col = cols["sample"]
    impact_col = cols["impact"]
    symbol_col = cols["symbol"]
    position_col = cols["pos"]

    df = combine_gnomad_nfe(df)
    df[cohort_col] = standardize_cohort_labels(df[cohort_col])
    df[tissue_col] = standardize_tissue_labels(df[tissue_col])
    df["Variant_ID"] = build_variant_id_series(df["Existing_variation"], df[symbol_col], df["HGVSp"])
    df["Variant_Label"] = df["Existing_variation"].apply(extract_rsid)
    df["Is_Common_SNP_By_NFE"] = df["gnomAD_NFE_AF_combined"].fillna(-1).gt(THRESHOLDS["min_nfe_af"])
    df["NFE_SNP_Threshold"] = THRESHOLDS["min_nfe_af"]

    group_cols = [
        "Variant_ID", "CHROM", position_col, "REF", "ALT", symbol_col,
        "Consequence", impact_col, cohort_col, tissue_col,
        "gnomAD_NFE_AF_combined", "gnomAD_NFE_Source",
        "Is_Common_SNP_By_NFE", "NFE_SNP_Threshold",
    ]
    variant_counts = compute_carrier_percentage(
        df=df,
        sample_col=sample_col,
        group_cols=group_cols,
        denom_cols=[cohort_col, tissue_col],
        pct_name='Carrier_Percentage',
    )

    variant_labels = (
        df.groupby(group_cols, dropna=False)["Variant_Label"]
        .agg(lambda values: '; '.join(sorted({str(v).strip() for v in values if str(v).strip()})) or None)
        .reset_index()
    )
    variant_counts = variant_counts.merge(variant_labels, on=group_cols, how='left')
    variant_counts["Labelled_On_Figure"] = (
        variant_counts[impact_col].isin(label_impacts()) & variant_counts["Variant_Label"].notna()
    )
    variant_counts["SNP_Context"] = variant_counts["Is_Common_SNP_By_NFE"].map({
        True: "Common SNP by gnomAD NFE threshold",
        False: "Rare or unconfirmed SNP by gnomAD NFE threshold",
    })

    warning_lookup = build_amplicon_warning_lookup(df, variant_col="Variant_ID", chrom_col="CHROM", pos_col=position_col)
    variant_counts = attach_amplicon_warning_columns(variant_counts, warning_lookup)

    sample_sizes = df.groupby([cohort_col, tissue_col])[sample_col].nunique().to_dict()
    return df, variant_counts, sample_sizes


def export_landscape_workbook(variant_counts, cols):
    """Write the plotted point table plus labelled/common-SNP subsets to Excel."""
    impact_col = cols["impact"]
    export_df = variant_counts.copy()
    export_df = export_df.rename(columns={
        cols["cohort"]: "Cohort",
        cols["tissue"]: "Tissue",
        cols["pos"]: "POS",
        cols["symbol"]: "SYMBOL",
        impact_col: "IMPACT",
    })
    export_df = export_df.rename(columns={"gnomAD_NFE_AF_combined": "gnomAD_NFE_AF"})

    preferred_cols = [
        "Variant_ID", "Variant_Label", "Labelled_On_Figure",
        "Is_Common_SNP_By_NFE", "SNP_Context", "NFE_SNP_Threshold",
        "Coverage_Risk_Flag", "Coverage_Risk_Amplicon", "Coverage_Risk_Region", "Coverage_Risk_Note",
        "CHROM", "POS", "REF", "ALT", "SYMBOL", "Consequence", "IMPACT",
        "Cohort", "Tissue", "Carrier_Count", "Total_Samples", "Carrier_Percentage",
        "gnomAD_NFE_AF", "gnomAD_NFE_Source",
    ]
    remaining_cols = [col for col in export_df.columns if col not in preferred_cols]
    export_df = export_df[preferred_cols + remaining_cols]
    export_df = export_df.sort_values(["Cohort", "Tissue", "POS", "Variant_ID"], kind="stable")

    labelled_df = export_df[export_df["Labelled_On_Figure"]].copy()
    common_snp_df = export_df[export_df["Is_Common_SNP_By_NFE"]].copy()

    summary_df = pd.DataFrame([
        {
            "Output": "Plotted point rows",
            "Count": int(len(export_df)),
            "Notes": "One row per plotted variant point within a cohort/tissue stratum.",
        },
        {
            "Output": "Labelled figure rows",
            "Count": int(len(labelled_df)),
            "Notes": "Subset of points labelled directly on the landscape figure.",
        },
        {
            "Output": "Common SNP rows by NFE",
            "Count": int(len(common_snp_df)),
            "Notes": f"Uses the shared gnomAD NFE threshold of {THRESHOLDS['min_nfe_af']:.2f}.",
        },
    ])

    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="summary", index=False)
        export_df.to_excel(writer, sheet_name="landscape_points", index=False)
        labelled_df.to_excel(writer, sheet_name="labelled_variants", index=False)
        common_snp_df.to_excel(writer, sheet_name="common_snps_by_nfe", index=False)


def generate_clean_all_impact_map():
    """Generate the all-genes landscape plot used for descriptive comparison."""
    print("--- Generating clean landscape map across all impact categories ---")

    # Load the consolidated variant table and confirm that the structural fields
    # required for faceting, counting and labelling are available.
    validate_file_exists(input_file, "Script 08 input")
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')

    pos_c = find_col(df, 'Pos')
    sym_c = find_col(df, 'Symbol')
    imp_c = find_col(df, 'Impact')
    coh_c = find_col(df, 'Cohort')
    tis_c = find_col(df, 'Tissue')
    sam_c = find_col(df, 'Sample')

    required_cols = [
        pos_c, sym_c, imp_c, coh_c, tis_c, sam_c,
        'Existing_variation', 'HGVSp', 'CHROM', 'REF', 'ALT',
        'Consequence', 'gnomADe_NFE_AF', 'gnomADg_NFE_AF',
    ]
    validate_required_columns(df, required_cols, "Script 08 Biological_Annotations")
    print_validation_summary(df, sam_c, "Script 08 raw annotations", [coh_c, tis_c])

    # Harmonise the point table once so the figure and workbook always describe
    # exactly the same variant set and the same SNP definition.
    df[imp_c] = df[imp_c].astype(str).str.strip().str.upper()
    _, variant_counts, sample_sizes = build_landscape_summary(
        df,
        {
            "cohort": coh_c,
            "tissue": tis_c,
            "sample": sam_c,
            "impact": imp_c,
            "symbol": sym_c,
            "pos": pos_c,
        },
    )
    export_landscape_workbook(
        variant_counts,
        {
            "cohort": coh_c,
            "tissue": tis_c,
            "impact": imp_c,
            "symbol": sym_c,
            "pos": pos_c,
        },
    )

    # Plot all variants in the same genomic coordinate system, while using marker
    # shape to distinguish impact classes and colour to distinguish genes.
    sns.set_style("whitegrid")
    marker_map = {
        "HIGH": "X",
        "MODERATE": "o",
        "MODIFIER": "s",
        "LOW": "v"
    }

    variant_counts[imp_c] = pd.Categorical(
        variant_counts[imp_c],
        categories=["HIGH", "MODERATE", "MODIFIER", "LOW"],
        ordered=True
    )

    gene_count = variant_counts[sym_c].nunique()
    if gene_count <= len(QUALITATIVE_COLORBLIND_SEQUENCE):
        point_palette = QUALITATIVE_COLORBLIND_SEQUENCE
    else:
        point_palette = sns.color_palette("husl", n_colors=gene_count)

    g = sns.relplot(
        data=variant_counts,
        x=pos_c,
        y='Carrier_Percentage',
        hue=sym_c,
        style=imp_c,
        row=coh_c,
        col=tis_c,
        kind='scatter',
        markers=marker_map,
        s=120,
        alpha=0.7,
        edgecolor="black",
        palette=point_palette,
        height=5,
        aspect=1.6,
        facet_kws={'sharex': True, 'sharey': True}
    )

    # Only the higher-priority consequence classes receive labels, and those
    # labels are restricted to rsIDs to avoid long overlapping annotations.
    labelled = variant_counts[variant_counts['Labelled_On_Figure']].copy()
    for (cohort_label, tissue_label), panel_df in labelled.groupby([coh_c, tis_c], sort=False):
        ax = g.axes_dict.get((cohort_label, tissue_label))
        if ax is None:
            continue
        annotate_panel_variants(ax, panel_df, pos_c, imp_c)

    # Present coordinates in Mb and frequencies as percentages so the figure is
    # immediately publication-friendly without requiring post-processing.
    g.set_axis_labels("Chr17 Position (Mb)", "Carrier Frequency (%)")

    for (cohort_label, tissue_label), ax in g.axes_dict.items():
        n_samples = sample_sizes.get((cohort_label, tissue_label), 0)
        ax.set_title(f"{cohort_label} | {tissue_label} (n={n_samples})", fontweight='bold')

    for ax in g.axes.flat:
        ax.set_ylim(0, 100)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.2f}'))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f'{y:.0f}%'))

    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(
        'Global Variant Landscape',
        fontsize=16,
        fontweight='bold'
    )

    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"SUCCESS: Clean map with all impacts saved to: {output_image}")
    print(f"SUCCESS: Landscape workbook with SNP context saved to: {output_excel}")


if __name__ == "__main__":
    generate_clean_all_impact_map()
