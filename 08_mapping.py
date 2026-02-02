import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- Configuration ---
input_file = "/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_Annotated_Report_Fixed.xlsx"
output_image = "/home/gadeaalonsoj/tfm/gsdmb_final_results/clean_landscape_all_impacts.png"

def find_col(df, target):
    """Case-insensitive column finder."""
    for col in df.columns:
        if col.upper() == target.upper(): return col
    return None

def generate_clean_all_impact_map():
    print("--- Generating Clean Map (All Impacts, No Annotations) ---")
    
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

    # 3. Clean the Impact column to ensure "HIGH" and "MODIFIER" are definitely included
    # We strip any whitespace and make it uppercase to avoid missing categories
    df[imp_c] = df[imp_c].astype(str).str.strip().str.upper()

    # 4. Aggregate data to count frequency per group
    group_cols = [pos_c, sym_c, imp_c, coh_c, tis_c]
    variant_counts = df.groupby(group_cols).size().reset_index(name='Frequency')

    # 5. Visualisation Setup
    sns.set_style("whitegrid")
    
    # Explicitly define markers for ALL impact types
    # This ensures HIGH and MODIFIER are distinct and visible
    marker_map = {
        "HIGH": "X",       # Bold Cross
        "MODERATE": "o",   # Circle
        "MODIFIER": "s",   # Square
        "LOW": "v"         # Downward Triangle
    }

    # Ensure the plot order respects the importance
    variant_counts[imp_c] = pd.Categorical(
        variant_counts[imp_c], 
        categories=["HIGH", "MODERATE", "MODIFIER", "LOW"], 
        ordered=True
    )

    # 6. Create the Plot
    g = sns.relplot(
        data=variant_counts,
        x=pos_c,
        y='Frequency',
        hue=sym_c,
        style=imp_c,
        row=coh_c,
        col=tis_c,
        kind='scatter',
        markers=marker_map,
        s=120,          # Size of the dots
        alpha=0.7,      # Slight transparency to handle overlaps
        edgecolor="black",
        palette="husl", # High-contrast colours for genes
        height=5,
        aspect=1.6,
        facet_kws={'sharex': True, 'sharey': False}
    )

    # 7. Formatting
    g.set_axis_labels("Genomic Position on Chr17 (Mb)", "Number of Samples")
    g.set_titles("{row_name} | {col_name}", fontweight='bold')
    
    # Format X-axis to show Mb (e.g., 39.91 Mb)
    for ax in g.axes.flat:
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.2f}'))

    plt.subplots_adjust(top=0.92)
    g.fig.suptitle('Landscape of Variants', fontsize=18, fontweight='bold')

    # 8. Save
    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"SUCCESS: Clean map with all impacts saved to: {output_image}")

if __name__ == "__main__":
    generate_clean_all_impact_map()