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
- ``clean_landscape_all_impacts.png`` in the configured results directory.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pipeline_utils import compute_carrier_percentage, extract_rsid, find_col, get_paths
from pipeline_validation import print_validation_summary, validate_file_exists, validate_required_columns

# Centralised path configuration keeps this plotting script aligned with the
# rest of the pipeline and reduces the risk of stale hard-coded locations.
PATHS = get_paths()
input_file = str(PATHS["annotated_report"])
output_image = str(PATHS["results_dir"] / "clean_landscape_all_impacts.png")


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
    existing_c = find_col(df, 'Existing_variation')

    required_cols = [pos_c, sym_c, imp_c, coh_c, tis_c, sam_c]
    validate_required_columns(df, required_cols, "Script 08 Biological_Annotations")
    print_validation_summary(df, sam_c, "Script 08 raw annotations", [coh_c, tis_c])

    # Harmonise impact labels and extract concise rsID-based labels so that the
    # annotated figure remains interpretable in manuscript-style presentation.
    df[imp_c] = df[imp_c].astype(str).str.strip().str.upper()
    df['Variant_Label'] = df[existing_c].apply(extract_rsid) if existing_c else None

    # Collapse repeated transcript-level rows to sample-level carrier summaries.
    # The denominator is cohort+tissue specific, which prevents larger groups
    # from appearing artificially more variant-rich simply because they contain
    # more sequenced samples.
    group_cols = [pos_c, sym_c, imp_c, coh_c, tis_c]
    variant_counts = compute_carrier_percentage(
        df=df,
        sample_col=sam_c,
        group_cols=group_cols,
        denom_cols=[coh_c, tis_c],
        pct_name='Carrier_Percentage',
    )
    variant_labels = (
        df.groupby(group_cols)['Variant_Label']
        .agg(lambda values: '; '.join(pd.unique([v for v in values if v])) or None)
        .reset_index()
    )
    variant_counts = variant_counts.merge(variant_labels, on=group_cols, how='left')

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
        palette="husl",
        height=5,
        aspect=1.6,
        facet_kws={'sharex': True, 'sharey': False}
    )

    # Only the higher-priority consequence classes receive labels, and those
    # labels are restricted to rsIDs to avoid long overlapping annotations.
    label_offsets = {
        'HIGH': (6, 6),
        'MODERATE': (6, -10)
    }
    labelled = variant_counts[
        variant_counts[imp_c].isin(['HIGH', 'MODERATE']) & variant_counts['Variant_Label'].notna()
    ].copy()
    for _, row in labelled.iterrows():
        ax = g.axes_dict.get((row[coh_c], row[tis_c]))
        if ax is None:
            continue
        offset = label_offsets.get(row[imp_c], (6, 6))
        ax.annotate(
            row['Variant_Label'],
            xy=(row[pos_c], row['Carrier_Percentage']),
            xytext=offset,
            textcoords='offset points',
            fontsize=8,
            fontweight='bold' if row[imp_c] == 'HIGH' else 'normal',
            color='black',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.7)
        )

    # Present coordinates in Mb and frequencies as percentages so the figure is
    # immediately publication-friendly without requiring post-processing.
    g.set_axis_labels("Genomic Position on Chr17 (Mb)", "Carrier Percentage of Samples")
    g.set_titles("{row_name} | {col_name}", fontweight='bold')

    for ax in g.axes.flat:
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.2f}'))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f'{y:.0f}%'))

    plt.subplots_adjust(top=0.92)
    g.fig.suptitle('Landscape of Variants', fontsize=18, fontweight='bold')

    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"SUCCESS: Clean map with all impacts saved to: {output_image}")


if __name__ == "__main__":
    generate_clean_all_impact_map()
