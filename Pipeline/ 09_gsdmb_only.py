"""
Script 09: GSDMB-Focused Landscape Mapping
==========================================

Purpose
-------
Generate a faceted genomic landscape restricted to GSDMB variants so that the
study gene can be examined without the visual competition of the broader panel.

Methodological rationale
------------------------
Carrier percentage is plotted instead of raw counts. This ensures that apparent
burden differences between cohort/tissue panels reflect relative prevalence and
not unequal sample size.

Input
-----
- ``Biological_Annotations`` sheet from the consolidated annotated report.

Output
------
- ``GSDMB_exclusive_landscape.png`` in the configured results directory.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from figure_style import COMPARATIVE_TAG, IMPACT_COLORS
from pipeline_utils import compute_carrier_percentage, ensure_directory, extract_rsid, find_col, get_paths
from pipeline_validation import print_validation_summary, validate_file_exists, validate_required_columns

# Use the shared pipeline configuration so that this gene-focused figure reads
# from the same canonical inputs as the rest of the analytical workflow.
PATHS = get_paths()
input_file = str(PATHS["annotated_report"])
output_dir = ensure_directory(PATHS.get("gsdmb_landscape_dir", PATHS["results_dir"] / "09_gsdmb_landscape"))
output_image = str(output_dir / "GSDMB_Target_Gene_Variant_Landscape.png")
LABEL_MODE = "high_only"  # valid options: "high_only", "both"


def label_impacts():
    """Return the impact classes that should receive text labels."""
    return ["HIGH"] if LABEL_MODE == "high_only" else ["HIGH", "MODERATE"]


def annotate_panel_variants(ax, panel_df, pos_c, impact_c):
    """Place readable labels while skipping crowded overlaps in each panel."""
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
            abs(x_val - placed_x) < 20000 and abs(y_val - placed_y) < 5
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


def generate_gsdmb_map():
    """Build the GSDMB-only landscape plot used for targeted interpretation."""
    print("--- Generating GSDMB-exclusive mutation map ---")

    validate_file_exists(input_file, "Script 09 input")
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')

    pos_c = find_col(df, 'Pos')
    sym_c = find_col(df, 'Symbol')
    imp_c = find_col(df, 'Impact')
    coh_c = find_col(df, 'Cohort')
    tis_c = find_col(df, 'Tissue')
    sam_c = find_col(df, 'Sample')
    existing_c = find_col(df, 'Existing_variation')

    required_cols = [pos_c, sym_c, imp_c, coh_c, tis_c, sam_c]
    validate_required_columns(df, required_cols, "Script 09 Biological_Annotations")

    # Restrict the analysis to the study gene of interest before any summary
    # statistics are generated.
    df_gsdmb = df[df[sym_c].astype(str).str.contains('GSDMB', case=False, na=False)].copy()
    print_validation_summary(df_gsdmb, sam_c, "Script 09 GSDMB subset", [coh_c, tis_c])

    if df_gsdmb.empty:
        print("Warning: No GSDMB variants found in the dataset!")
        return

    # Standardise impact labels and keep only compact rsID-based annotations for
    # the higher-priority classes to preserve readability.
    df_gsdmb[imp_c] = df_gsdmb[imp_c].astype(str).str.strip().str.upper()
    df_gsdmb['Variant_Label'] = df_gsdmb[existing_c].apply(extract_rsid) if existing_c else None

    # Summarise unique carrier prevalence per variant within each cohort/tissue
    # panel, using unique samples as both numerator and denominator.
    group_cols = [pos_c, imp_c, coh_c, tis_c]
    variant_counts = compute_carrier_percentage(
        df=df_gsdmb,
        sample_col=sam_c,
        group_cols=group_cols,
        denom_cols=[coh_c, tis_c],
        pct_name='Carrier_Percentage',
    )
    variant_labels = (
        df_gsdmb.groupby(group_cols)['Variant_Label']
        .agg(lambda values: '; '.join(pd.unique([v for v in values if v])) or None)
        .reset_index()
    )
    variant_counts = variant_counts.merge(variant_labels, on=group_cols, how='left')

    sns.set_style("whitegrid")
    marker_map = {"HIGH": "X", "MODERATE": "o", "MODIFIER": "s", "LOW": "v"}

    variant_counts[imp_c] = pd.Categorical(
        variant_counts[imp_c],
        categories=["HIGH", "MODERATE", "MODIFIER", "LOW"],
        ordered=True
    )

    # Here colour and marker both encode impact class, because the gene is held
    # constant and the visual goal is to emphasise biological severity.
    g = sns.relplot(
        data=variant_counts,
        x=pos_c,
        y='Carrier_Percentage',
        hue=imp_c,
        style=imp_c,
        row=coh_c,
        col=tis_c,
        kind='scatter',
        markers=marker_map,
        palette=IMPACT_COLORS,
        s=150,
        alpha=0.8,
        edgecolor="black",
        height=5,
        aspect=1.4,
        facet_kws={'sharex': True, 'sharey': True}
    )

    labelled = variant_counts[
        variant_counts[imp_c].isin(label_impacts()) & variant_counts['Variant_Label'].notna()
    ].copy()
    for (cohort_label, tissue_label), panel_df in labelled.groupby([coh_c, tis_c], sort=False):
        ax = g.axes_dict.get((cohort_label, tissue_label))
        if ax is None:
            continue
        annotate_panel_variants(ax, panel_df, pos_c, imp_c)

    g.set_axis_labels("Chr17 Position (Mb)", "Carrier Frequency (%)")
    g.set_titles("{row_name} | {col_name}", fontweight='bold')

    for ax in g.axes.flat:
        ax.set_ylim(0, 100)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.2f}'))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f'{y:.0f}%'))

    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(
        'GSDMB Variant Landscape',
        fontsize=16,
        fontweight='bold'
    )

    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"SUCCESS: GSDMB-exclusive map saved to: {output_image}")


if __name__ == "__main__":
    generate_gsdmb_map()
