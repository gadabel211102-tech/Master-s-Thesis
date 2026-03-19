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

from pipeline_utils import compute_carrier_percentage, ensure_directory, extract_rsid, find_col, get_paths
from figure_style import COMPARATIVE_TAG, QUALITATIVE_COLORBLIND_SEQUENCE
from pipeline_validation import print_validation_summary, validate_file_exists, validate_required_columns

# Centralised path configuration keeps this plotting script aligned with the
# rest of the pipeline and reduces the risk of stale hard-coded locations.
PATHS = get_paths()
input_file = str(PATHS["annotated_report"])
output_dir = ensure_directory(PATHS.get("variant_landscape_dir", PATHS["results_dir"] / "08_variant_landscape"))
output_image = str(output_dir / "GSDMB_Global_Variant_Landscape.png")
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
    sample_sizes = df.groupby([coh_c, tis_c])[sam_c].nunique().to_dict()

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
        palette=QUALITATIVE_COLORBLIND_SEQUENCE,
        height=5,
        aspect=1.6,
        facet_kws={'sharex': True, 'sharey': True}
    )

    # Only the higher-priority consequence classes receive labels, and those
    # labels are restricted to rsIDs to avoid long overlapping annotations.
    labelled = variant_counts[
        variant_counts[imp_c].isin(label_impacts()) & variant_counts['Variant_Label'].notna()
    ].copy()
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


if __name__ == "__main__":
    generate_clean_all_impact_map()
