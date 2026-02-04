import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- Configuration ---
# Set input data path and the destination for the final PNG plot
input_file = "/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_Annotated_Report_Fixed.xlsx"
output_image = "/home/gadeaalonsoj/tfm/gsdmb_final_results/clean_landscape_all_impacts.png"

def find_col(df, target):
    """Case-insensitive column finder."""
    # Ensures the script doesn't crash if 'Pos' is renamed to 'POS' or 'pos'
    for col in df.columns:
        if col.upper() == target.upper(): return col
    return None

def generate_clean_all_impact_map():
    print("--- Generating Clean Map (All Impacts, No Annotations) ---")
    
    if not os.path.exists(input_file):
        print(f"Error: File not found at {input_file}")
        return

    # 1. Load the biological annotations sheet from the Excel report
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')

    # 2. Identify core columns needed for the plot
    pos_c = find_col(df, 'Pos')
    sym_c = find_col(df, 'Symbol')
    imp_c = find_col(df, 'Impact')
    coh_c = find_col(df, 'Cohort')
    tis_c = find_col(df, 'Tissue')

    # 3. Standardise the Impact column strings to prevent category duplication
    df[imp_c] = df[imp_c].astype(str).str.strip().str.upper()

    # 4. Count the number of samples carrying each specific variant
    group_cols = [pos_c, sym_c, imp_c, coh_c, tis_c]
    variant_counts = df.groupby(group_cols).size().reset_index(name='Frequency')

    # 5. Visualisation Setup
    sns.set_style("whitegrid")
    
    # Map VEP impact categories to specific geometric shapes for clarity
    marker_map = {
        "HIGH": "X",       # Bold Cross (Critical)
        "MODERATE": "o",   # Circle
        "MODIFIER": "s",   # Square
        "LOW": "v"          # Downward Triangle
    }

    # Set a logical order for the legend (High to Low importance)
    variant_counts[imp_c] = pd.Categorical(
        variant_counts[imp_c], 
        categories=["HIGH", "MODERATE", "MODIFIER", "LOW"], 
        ordered=True
    )

    # 6. Create the Plot (Relational Plot)
    # This creates a grid (Facets) based on Cohort (rows) and Tissue (columns)
    g = sns.relplot(
        data=variant_counts,
        x=pos_c,
        y='Frequency',
        hue=sym_c,      # Colour-code by Gene symbol
        style=imp_c,    # Marker shape by Impact level
        row=coh_c,      # Separate horizontal plots by Cohort
        col=tis_c,      # Separate vertical plots by Tissue
        kind='scatter',
        markers=marker_map,
        s=120,          # Marker size
        alpha=0.7,      # Transparency to show overlapping points
        edgecolor="black",
        palette="husl", # High-contrast colour palette
        height=5,
        aspect=1.6,
        facet_kws={'sharex': True, 'sharey': False}
    )

    # 7. Formatting and Labelling
    g.set_axis_labels("Genomic Position on Chr17 (Mb)", "Number of Samples")
    g.set_titles("{row_name} | {col_name}", fontweight='bold')
    
    # Transform raw genomic coordinates (base pairs) into Megabases (Mb) for the X-axis
    for ax in g.axes.flat:
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.2f}'))

    plt.subplots_adjust(top=0.92)
    g.fig.suptitle('Landscape of Variants', fontsize=18, fontweight='bold')

    # 8. Save the final high-resolution (300 DPI) image
    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"SUCCESS: Clean map with all impacts saved to: {output_image}")

if __name__ == "__main__":
    generate_clean_all_impact_map()
