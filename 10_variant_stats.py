"""
Script 10: Descriptive Variant Composition Statistics
=====================================================

Purpose
-------
Summarise the distribution of variant consequences and impact classes across
study groups using descriptive count tables and stacked bar plots.

Interpretive note
-----------------
This script is descriptive rather than inferential. It reports raw event counts
and unique-variant counts to characterise composition, whereas sample-size-
normalised comparative analyses are handled in later scripts.

Input
-----
- ``GSDMB_Annotated_Report_Fixed.xlsx`` from the results directory.

Outputs
-------
- CSV summary tables for impact and consequence counts.
- Four stacked bar plots for total versus unique variants, each in linear and
  log-scaled form.
"""

import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from figure_style import DESCRIPTIVE_TAG, IMPACT_COLORS

# This script intentionally keeps its historical output directory so that the
# descriptive figures continue to land in the same folder used in earlier runs.
input_dir = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(input_dir, "GSDMB_Annotated_Report_Fixed.xlsx")

# The four exported figures provide the same descriptive content under linear
# and log scaling, and for total versus deduplicated variant counts.
outputs = {
    "total_log": "01_ALL_VARIANTS_TOTAL_LOG.png",
    "total_lin": "02_ALL_VARIANTS_TOTAL_LINEAR.png",
    "unique_log": "03_ALL_VARIANTS_UNIQUE_LOG.png",
    "unique_lin": "04_ALL_VARIANTS_UNIQUE_LINEAR.png",
    "total_pct": "05_ALL_VARIANTS_TOTAL_PERCENT.png",
    "unique_pct": "06_ALL_VARIANTS_UNIQUE_PERCENT.png",
}


def find_col(df, target):
    """Return a case-insensitive column match from the annotated table."""
    for col in df.columns:
        if col.strip().upper() == target.upper():
            return col
    return None


def generate_and_save_tables(df, suffix, output_dir):
    """Export descriptive impact and consequence count tables for one dataset."""
    imp_c = find_col(df, 'Impact')
    con_c = find_col(df, 'Consequence')

    print(f"\n>>> PROCESSING TABLES FOR: {suffix.upper()}")

    # Impact classes are ordered explicitly so that spreadsheets and figures use
    # the conventional VEP severity ranking rather than alphabetical order.
    impact_table = df.groupby(['Group', imp_c]).size().unstack(fill_value=0)
    for cat in ["HIGH", "MODERATE", "MODIFIER", "LOW"]:
        if cat not in impact_table.columns:
            impact_table[cat] = 0
    impact_table = impact_table[["HIGH", "MODERATE", "MODIFIER", "LOW"]]
    impact_table.to_csv(os.path.join(output_dir, f"table_impact_{suffix}.csv"))

    con_table = df.groupby(['Group', con_c]).size().unstack(fill_value=0)
    con_table.to_csv(os.path.join(output_dir, f"table_consequence_{suffix}.csv"))

    print(f"Successfully exported {suffix} tables.")


def create_stat_visuals(df_to_plot, title_suffix, filename, use_log=False, as_percentage=False):
    """Render paired stacked-bar summaries for consequence and impact profiles."""
    imp_c = find_col(df_to_plot, 'Impact')
    con_c = find_col(df_to_plot, 'Consequence')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 15))
    sns.set_style("whitegrid")

    con_pivot = df_to_plot.groupby(['Group', con_c]).size().unstack(fill_value=0)
    impact_pivot = df_to_plot.groupby(['Group', imp_c]).size().unstack(fill_value=0)

    all_cats = ["HIGH", "MODERATE", "MODIFIER", "LOW"]
    for c in all_cats:
        if c not in impact_pivot.columns:
            impact_pivot[c] = 0
    impact_pivot = impact_pivot[all_cats]

    if as_percentage:
        con_pivot = con_pivot.div(con_pivot.sum(axis=1), axis=0).fillna(0) * 100
        impact_pivot = impact_pivot.div(impact_pivot.sum(axis=1), axis=0).fillna(0) * 100

    con_pivot.plot(kind='bar', stacked=True, ax=ax1, colormap='tab20')
    ax1.set_title(f'Distribution of Functional Consequences: {title_suffix}', fontsize=16, fontweight='bold')
    if use_log and not as_percentage:
        ax1.set_yscale('log')
    ax1.set_ylabel('Within-group percentage (%)' if as_percentage else ('Count (Log Scale)' if use_log else 'Count'))
    ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small', title='Consequence')

    impact_colors = [IMPACT_COLORS[c] for c in all_cats]
    impact_pivot.plot(kind='bar', stacked=True, ax=ax2, color=impact_colors)
    ax2.set_title(f'Distribution of Predicted Biological Impact: {title_suffix}', fontsize=16, fontweight='bold')
    if use_log and not as_percentage:
        ax2.set_yscale('log')
    ax2.set_ylabel('Within-group percentage (%)' if as_percentage else ('Count (Log Scale)' if use_log else 'Count'))
    ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small', title='Impact')

    for ax, pivot in [(ax1, con_pivot), (ax2, impact_pivot)]:
        ax.set_xlabel('Study group')
        ax.tick_params(axis='x', rotation=0)
        for p in ax.patches:
            h = p.get_height()
            if h <= 0:
                continue
            y = p.get_y() + h / 2
            label = f'{h:.0f}%' if as_percentage else f'{int(round(h))}'
            if h >= (6 if as_percentage else 1):
                ax.text(p.get_x() + p.get_width() / 2, y, label, ha='center', va='center',
                        fontsize=7, fontweight='bold', color='white' if h > 10 else 'black')

    fig.suptitle(
        f"Summary of Variant Annotation Profiles: {'within-group composition' if as_percentage else 'raw variant burden'}",
        fontsize=15, fontweight='bold', y=0.98
    )
    fig.text(0.5, 0.01,
             'Descriptive summary only; inferential tumour-versus-control comparisons are handled in later scripts.',
             ha='center', fontsize=9, color='#666666', style='italic')
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(os.path.join(input_dir, filename), dpi=300, bbox_inches='tight')
    plt.close()


def main():
    """Load the annotated report, derive grouped summaries, and export plots."""
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')

    sym_c = find_col(df, 'Symbol')
    tis_c = find_col(df, 'Tissue')
    coh_c = find_col(df, 'Cohort')
    pos_c = find_col(df, 'Pos')
    hgv_c = find_col(df, 'HGVSp')

    # Group labels combine cohort and tissue so that descriptive summaries align
    # with the same biological partitions used across later scripts.
    df[tis_c] = df[tis_c].astype(str).str.strip()
    df[coh_c] = df[coh_c].astype(str).str.strip()
    df['Group'] = df[coh_c] + " - " + df[tis_c]

    # Total rows reflect the complete annotated call set, including repeated
    # observations of the same variant in different samples.
    generate_and_save_tables(df, "all_total", input_dir)
    create_stat_visuals(df, "All Variants", outputs['total_log'], use_log=True)
    create_stat_visuals(df, "All Variants", outputs['total_lin'], use_log=False)
    create_stat_visuals(df, "All Variants", outputs['total_pct'], as_percentage=True)

    # Unique rows collapse duplicate observations within each study group so the
    # descriptive focus shifts from burden to repertoire breadth.
    df_u = df.drop_duplicates(subset=[sym_c, pos_c, hgv_c, 'Group'])

    generate_and_save_tables(df_u, "all_unique", input_dir)
    create_stat_visuals(df_u, "Unique Variants", outputs['unique_log'], use_log=True)
    create_stat_visuals(df_u, "Unique Variants", outputs['unique_lin'], use_log=False)
    create_stat_visuals(df_u, "Unique Variants", outputs['unique_pct'], as_percentage=True)

    print(f"\n{'='*60}")
    print(f" SUCCESS: Results for ALL variants generated in:\n {input_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
