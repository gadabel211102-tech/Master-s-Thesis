import pandas as pd
import numpy as np
import scipy.stats as stats
import os
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

# --- PATH CONFIGURATION ---
base_path = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_xlsx = os.path.join(base_path, "18_Cancer_Specific_Association_Results.xlsx")
output_plot = os.path.join(base_path, "18_Cancer_Specific_Volcano_Plots.png")

def find_col(df, target):
    """Searches for columns regardless of capitalisation or hidden spaces."""
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None

def run_multi_cohort_association():
    print(">>> Starting Multi-Cohort Association Analysis (Global, Breast, Endometrium)...")
    
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        return

    # 1. LOAD DATA
    # We use the Biological_Annotations sheet as the source of truth for all variants
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')
    
    # 2. IDENTIFY AND STANDARDISE COLUMNS
    # Dynamically find column names to prevent script crashes if names shift slightly
    sym_c = find_col(df, 'SYMBOL')
    imp_c = find_col(df, 'IMPACT')
    con_c = find_col(df, 'CONSEQUENCE')
    tis_c = find_col(df, 'TISSUE')
    sam_c = find_col(df, 'SAMPLE')
    var_c = find_col(df, 'Existing_variation')
    hgv_c = find_col(df, 'HGVSp')
    coh_c = find_col(df, 'COHORT')

    # Standardise Tissue names to British spelling for the thesis
    df[tis_c] = df[tis_c].astype(str).str.strip().replace({'Tumor': 'Tumour', 'Normal': 'Healthy', 'Control': 'Healthy'})
    df[coh_c] = df[coh_c].astype(str).str.strip()
    
    # Create a unique Variant ID: Priority 1 = rsID; Priority 2 = Gene:Protein Change
    df['rsID'] = df[var_c].astype(str).str.extract(r'(rs\d+)')
    df['Variant_ID'] = df['rsID'].fillna(df[sym_c].astype(str) + ":" + df[hgv_c].astype(str))
    
    # NO-SKIP LOGIC: Generate a temporary ID for variants with zero metadata to ensure math remains accurate
    replacement_series = pd.Series("Unknown_SNP_" + df.index.astype(str), index=df.index)
    df['Variant_ID'] = df['Variant_ID'].replace('nan:nan', np.nan).fillna(replacement_series)

    # Define the three levels of analysis requested
    cohort_list = ['Global', 'Breast', 'Endometrium']
    all_results = []

    # 3. PROCESS EACH LEVEL
    for cohort_name in cohort_list:
        print(f"   [PROCESS] Analyzing {cohort_name}...")
        
        # Filter data for specific tissue or keep all for 'Global'
        if cohort_name == 'Global':
            cohort_df = df.copy()
        else:
            cohort_df = df[df[coh_c].str.contains(cohort_name, case=False, na=False)].copy()
            
        if cohort_df.empty:
            continue
            
        # Denominators: Total unique samples available in this specific cohort/tissue
        n_tumour = cohort_df[cohort_df[tis_c] == 'Tumour'][sam_c].nunique()
        n_healthy = cohort_df[cohort_df[tis_c] == 'Healthy'][sam_c].nunique()
        
        # Guard against division by zero errors
        if n_tumour == 0 or n_healthy == 0:
            print(f"      [SKIP] Missing groups for {cohort_name}")
            continue

        unique_variants = cohort_df['Variant_ID'].unique()
        
        for var in unique_variants:
            var_data = cohort_df[cohort_df['Variant_ID'] == var]
            
            # Count how many unique samples carry this specific variant
            count_tumour = var_data[var_data[tis_c] == 'Tumour'][sam_c].nunique()
            count_healthy = var_data[var_data[tis_c] == 'Healthy'][sam_c].nunique()
            
            # FISHER'S EXACT TEST: Compare mutation frequency in Cases vs Controls
            # Table Structure:
            # [[Samples WITH in Cases, Samples WITHOUT in Cases],
            #  [Samples WITH in Controls, Samples WITHOUT in Controls]]
            table = [
                [count_tumour, max(0, n_tumour - count_tumour)],
                [count_healthy, max(0, n_healthy - count_healthy)]
            ]
            
            # Odds Ratio (Risk magnitude) and P-Value (Statistical probability)
            odds_ratio, p_value = stats.fisher_exact(table)
            
            # Helper to pull the first available metadata (Impact, Symbol, etc.) for this SNP
            def get_metadata(col_name):
                if col_name in var_data.columns and not var_data[col_name].empty:
                    val = var_data[col_name].iloc[0]
                    return val if pd.notna(val) else "N/A"
                return "N/A"

            all_results.append({
                'Analysis_Group': cohort_name,
                'Variant_ID': var,
                'Symbol': get_metadata(sym_c),
                'Impact': get_metadata(imp_c),
                'Consequence': get_metadata(con_c),
                'Tumour_Count': count_tumour,
                'Healthy_Count': count_healthy,
                'Tumour_Freq_%': (count_tumour / n_tumour) * 100,
                'Healthy_Freq_%': (count_healthy / n_healthy) * 100,
                'Odds_Ratio': odds_ratio,
                'P_Value': p_value
            })

    # 4. MULTI-TESTING CORRECTION (FDR)
    # We apply Benjamini-Hochberg (FDR) correction WITHIN each subgroup to control false discoveries
    res_df = pd.DataFrame(all_results)
    final_dfs = []
    for group in res_df['Analysis_Group'].unique():
        sub = res_df[res_df['Analysis_Group'] == group].copy()
        # FDR_P_Value is the one you should cite as "Robust Significance"
        _, sub['FDR_P_Value'], _, _ = multipletests(sub['P_Value'], method='fdr_bh')
        final_dfs.append(sub)
    
    final_res = pd.concat(final_dfs).sort_values(['Analysis_Group', 'P_Value'])
    
    # Export to a Multi-Sheet Excel (One tab per cohort)
    with pd.ExcelWriter(output_xlsx) as writer:
        for group in cohort_list:
            if group in final_res['Analysis_Group'].unique():
                final_res[final_res['Analysis_Group'] == group].to_excel(writer, sheet_name=group, index=False)
    
    print(f">>> Full report saved to: {output_xlsx}")

    # 5. VISUALISATION: SIDE-BY-SIDE VOLCANO PLOTS
    # Convert p-values to -log10 for standard volcano visualisation
    final_res['-log10_p'] = -np.log10(final_res['P_Value'].replace(0, 1e-20))
    
    # Cap the Odds Ratio at 40 to prevent infinite/outlier values from ruining the plot scale
    max_or_limit = min(final_res['Odds_Ratio'].replace(np.inf, np.nan).max(), 40)
    
    # Create a FacetGrid (one plot for each Analysis Group)
    g = sns.FacetGrid(final_res, col="Analysis_Group", hue="Impact", palette="viridis", height=6, aspect=1)
    g.map(sns.scatterplot, "Odds_Ratio", "-log10_p", s=100, edgecolor='black', alpha=0.7)
    
    # Add significance thresholds and axis labels
    for ax in g.axes.flat:
        ax.axhline(-np.log10(0.05), color='red', linestyle='--', label='Nominal p=0.05')
        ax.set_xlim(-1, max_or_limit + 5)
        ax.set_xlabel("Odds Ratio (Risk Factor)")
        ax.set_ylabel("-log10(P-Value)")
        ax.grid(True, alpha=0.2)

    g.add_legend(title="VEP Impact")
    g.set_titles("{col_name} Association Analysis", fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    print(f">>> Multi-cohort Volcano Plots saved to: {output_plot}")

if __name__ == "__main__":
    run_multi_cohort_association()
