import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- Configuration ---
input_file = "/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_Annotated_Report_Fixed.xlsx"
output_image = "/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_exclusive_landscape.png"

def find_col(df, target):
    for col in df.columns:
        if col.upper() == target.upper(): return col
    return None

def generate_gsdmb_map():
    print("--- Generating GSDMB-Exclusive Mutation Map ---")
    
    if not os.path.exists(input_file):
        print(f"Error: File not found at {input_file}")
        return

    # 1. Load data
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')

    # 2. Identify columns
    pos_c = find_col(df, 'Pos')
    sym_c = find_col(df, 'Symbol')
    imp_c = find_col(df, 'Impact')
    coh_c = find_col(df, 'Cohort')
    tis_c = find_col(df, 'Tissue')

    # 3. Filter for GSDMB ONLY
    # Using .str.contains to catch any slight naming variations
    df_gsdmb = df[df[sym_c].astype(str).str.contains('GSDMB', case=False)].copy()
    
    if df_gsdmb.empty:
        print("Warning: No GSDMB variants found in the dataset!")
        return

    # 4. Standardise Impact for the legend
    df_gsdmb[imp_c] = df_gsdmb[imp_c].astype(str).str.strip().str.upper()

    # 5. Aggregate
    variant_counts = df_gsdmb.groupby([pos_c, imp_c, coh_c, tis_c]).size().reset_index(name='Frequency')

    # 6. Visualisation
    sns.set_style("whitegrid")
    marker_map = {"HIGH": "X", "MODERATE": "o", "MODIFIER": "s", "LOW": "v"}
    
    # Use a specific colour for GSDMB (e.g., a strong Blue or Red)
    g = sns.relplot(
        data=variant_counts,
        x=pos_c,
        y='Frequency',
        hue=imp_c,       # In an exclusive plot, colour by Impact is more useful
        style=imp_c,
        row=coh_c,
        col=tis_c,
        kind='scatter',
        markers=marker_map,
        palette="flare", # Warm palette to highlight high-impact areas
        s=150,
        alpha=0.8,
        edgecolor="black",
        height=5,
        aspect=1.4,
        facet_kws={'sharex': True, 'sharey': True}
    )

    # 7. Add "Exon" context if possible (approximate GSDMB boundaries)
    # GSDMB is roughly between 39,905,000 and 39,920,000 on Chr17
    
    # 8. Formatting
    g.set_axis_labels("Position on Chr17 (bp)", "Number of Samples")
    g.set_titles("{row_name} | {col_name}", fontweight='bold')
    
    # Format X-axis with comma separators for full base-pair precision
    for ax in g.axes.flat:
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
        plt.setp(ax.get_xticklabels(), rotation=45)

    plt.subplots_adjust(top=0.9)
    g.fig.suptitle('GSDMB Mutation Landscape: Comparison by Cohort and Tissue', fontsize=18, fontweight='bold')

    # 9. Save
    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"SUCCESS: GSDMB-exclusive map saved to: {output_image}")

if __name__ == "__main__":
    generate_gsdmb_map()