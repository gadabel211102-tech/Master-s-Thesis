import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- Configuration ---
input_dir = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(input_dir, "GSDMB_Annotated_Report_Fixed.xlsx")

# Manifest of the 4 essential bar plots
outputs = {
    "total_log": "01_GSDMB_TOTAL_LOG.png",
    "total_lin": "02_GSDMB_TOTAL_LINEAR.png",
    "unique_log": "03_GSDMB_UNIQUE_LOG.png",
    "unique_lin": "04_GSDMB_UNIQUE_LINEAR.png"
}

def find_col(df, target):
    """Robust case-insensitive column finder."""
    for col in df.columns:
        if col.strip().upper() == target.upper(): return col
    return None

def generate_and_save_tables(df, suffix, output_dir):
    """Generates and exports CSV counts for Impact and Consequence."""
    imp_c = find_col(df, 'Impact')
    con_c = find_col(df, 'Consequence')
    
    print(f"\n>>> PROCESSING TABLES FOR: {suffix.upper()}")

    # 1. Impact Table (Categorised by the 4 standard levels)
    impact_table = df.groupby(['Group', imp_c]).size().unstack(fill_value=0)
    for cat in ["HIGH", "MODERATE", "MODIFIER", "LOW"]:
        if cat not in impact_table.columns: impact_table[cat] = 0
    impact_table = impact_table[["HIGH", "MODERATE", "MODIFIER", "LOW"]]
    impact_table.to_csv(os.path.join(output_dir, f"table_impact_{suffix}.csv"))

    # 2. Consequence Table
    con_table = df.groupby(['Group', con_c]).size().unstack(fill_value=0)
    con_table.to_csv(os.path.join(output_dir, f"table_consequence_{suffix}.csv"))
    
    print(f"Successfully exported {suffix} tables.")

def create_stat_visuals(df_to_plot, title_suffix, filename, use_log=False):
    """Creates a 2-panel vertical figure: Consequences (Top) and Impacts (Bottom)."""
    imp_c = find_col(df_to_plot, 'Impact')
    con_c = find_col(df_to_plot, 'Consequence')
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 15))
    sns.set_style("whitegrid")

    # Top Panel: Consequences
    con_pivot = df_to_plot.groupby(['Group', con_c]).size().unstack(fill_value=0)
    con_pivot.plot(kind='bar', stacked=True, ax=ax1, colormap='tab20')
    ax1.set_title(f'Mutation Consequences: {title_suffix}', fontsize=16, fontweight='bold')
    if use_log: ax1.set_yscale('log')
    ax1.set_ylabel('Count (Log Scale)' if use_log else 'Count')
    ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small')

    # Bottom Panel: Impacts
    impact_pivot = df_to_plot.groupby(['Group', imp_c]).size().unstack(fill_value=0)
    all_cats = ["HIGH", "MODERATE", "MODIFIER", "LOW"]
    for c in all_cats:
        if c not in impact_pivot.columns: impact_pivot[c] = 0
    impact_pivot = impact_pivot[all_cats]
    
    colors = ['#d9534f', '#5bc0de', '#f0ad4e', '#5cb85c'] # Standard TFM colors
    impact_pivot.plot(kind='bar', stacked=True, ax=ax2, color=colors)
    ax2.set_title(f'Biological Impact Level: {title_suffix}', fontsize=16, fontweight='bold')
    if use_log: ax2.set_yscale('log')
    ax2.set_ylabel('Count (Log Scale)' if use_log else 'Count')
    plt.xticks(rotation=0)

    # Add numeric labels to the Impact bars
    for p in ax2.patches:
        h = p.get_height()
        if h > 0:
            x, y = p.get_x() + p.get_width() / 2, p.get_y()
            y_pos = y * 1.15 if use_log and y > 0 else y + h/2
            ax2.text(x, y_pos, f'{int(h)}', ha='center', va='center', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(input_dir, filename), dpi=300)
    plt.close()

def main():
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    # Load Data
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')
    
    # Define columns
    sym_c = find_col(df, 'Symbol')
    tis_c = find_col(df, 'Tissue')
    coh_c = find_col(df, 'Cohort')
    pos_c = find_col(df, 'Pos')
    hgv_c = find_col(df, 'HGVSp')

    # Subset GSDMB and Sanitize
    df_gsdmb = df[df[sym_c].astype(str).str.contains('GSDMB', case=False, na=False)].copy()
    df_gsdmb[tis_c] = df_gsdmb[tis_c].astype(str).str.strip()
    df_gsdmb[coh_c] = df_gsdmb[coh_c].astype(str).str.strip()
    df_gsdmb['Group'] = df_gsdmb[coh_c] + " - " + df_gsdmb[tis_c]

    # --- 1. TOTALS ---
    generate_and_save_tables(df_gsdmb, "total", input_dir)
    create_stat_visuals(df_gsdmb, "All Records (Log)", outputs['total_log'], use_log=True)
    create_stat_visuals(df_gsdmb, "All Records (Linear)", outputs['total_lin'], use_log=False)

    # --- 2. UNIQUES (Removing duplicates within the same group) ---
    df_u = df_gsdmb.drop_duplicates(subset=[pos_c, hgv_c, 'Group'])
    generate_and_save_tables(df_u, "unique", input_dir)
    create_stat_visuals(df_u, "Unique Variants (Log)", outputs['unique_log'], use_log=True)
    create_stat_visuals(df_u, "Unique Variants (Linear)", outputs['unique_lin'], use_log=False)
    
    print(f"\n{'='*60}")
    print(f" SUCCESS: 4 Charts and 4 Tables generated in:\n {input_dir}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()