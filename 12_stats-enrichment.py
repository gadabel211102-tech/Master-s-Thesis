import pandas as pd
import numpy as np
import scipy.stats as stats
import os
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

# --- 1. PATH CONFIGURATION ---
base_path = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
# Standardising output names to reflect that this is now a SNP-specific analysis
output_xlsx = os.path.join(base_path, "12_SNP_Enrichment_Results.xlsx")
output_plot = os.path.join(base_path, "12_SNP_Volcano_Plots.png")

def find_col(df, target):
    """Searches for columns regardless of capitalisation or hidden spaces."""
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None

def run_snp_association_analysis():
    print(">>> Starting Script 12: SNP Statistical Enrichment (>1% AF Filter)...")
    
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        return

    # --- 2. DATA LOADING & CLEANING ---
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')
    
    # Dynamically identifying columns to handle potential formatting shifts
    sym_c = find_col(df, 'SYMBOL')
    imp_c = find_col(df, 'IMPACT')
    con_c = find_col(df, 'CONSEQUENCE')
    tis_c = find_col(df, 'TISSUE')
    sam_c = find_col(df, 'SAMPLE')
    var_c = find_col(df, 'Existing_variation')
    hgv_c = find_col(df, 'HGVSp')
    coh_c = find_col(df, 'COHORT')
    nfe_c = find_col(df, 'gnomADe_NFE_AF')

    # --- 3. THE STRICT SNP FILTER ---
    # We apply the >1% filter here to ensure Script 12 matches Script 11's data pool
    initial_count = len(df)
    df = df[df[nfe_c] > 0.01].copy()
    print(f"   [FILTER] Focus: SNPs (>1% AF). Kept {len(df)} of {initial_count} entries.")

    # Standardise to British spelling for the thesis
    df[tis_c] = df[tis_c].astype(str).str.strip().replace({'Tumor': 'Tumour', 'Normal': 'Healthy', 'Control': 'Healthy'})
    df[coh_c] = df[coh_c].astype(str).str.strip()
    
    # --- 4. ROBUST ID CREATION (Matches Script 11) ---
    df['rsID'] = df[var_c].astype(str).str.extract(r'(rs\d+)')
    df['Variant_ID'] = df['rsID'].fillna(df[sym_c].astype(str) + ":" + df[hgv_c].astype(str))
    
    # Fallback to ensure math accuracy even for variants with missing metadata
    replacement_series = pd.Series("SNP_" + df.index.astype(str), index=df.index)
    df['Variant_ID'] = df['Variant_ID'].replace('nan:nan', np.nan).fillna(replacement_series)

    # Define our three tiers of analysis
    cohort_list = ['Global', 'Breast', 'Endometrium']
    all_results = []

    # --- 5. FISHER'S EXACT TEST PROCESSING ---
    for cohort_name in cohort_list:
        print(f"   [STATS] Calculating enrichment for {cohort_name}...")
        
        # Filter for the specific cohort (Global includes everyone)
        if cohort_name == 'Global':
            cohort_df = df.copy()
        else:
            cohort_df = df[df[coh_c].str.contains(cohort_name, case=False, na=False)].copy()
            
        if cohort_df.empty: continue
            
        # Unique sample denominators for the cohort
        n_tumour = cohort_df[cohort_df[tis_c] == 'Tumour'][sam_c].nunique()
        n_healthy = cohort_df[cohort_df[tis_c] == 'Healthy'][sam_c].nunique()
        
        if n_tumour == 0 or n_healthy == 0: continue

        unique_snps = cohort_df['Variant_ID'].unique()
        
        for var in unique_snps:
            var_data = cohort_df[cohort_df['Variant_ID'] == var]
            
            # Count unique carriers per tissue
            count_tumour = var_data[var_data[tis_c] == 'Tumour'][sam_c].nunique()
            count_healthy = var_data[var_data[tis_c] == 'Healthy'][sam_c].nunique()
            
            # 2x2 Table: [[Carriers Case, Non-Carriers Case], [Carriers Ctrl, Non-Carriers Ctrl]]
            table = [
                [count_tumour, max(0, n_tumour - count_tumour)],
                [count_healthy, max(0, n_healthy - count_healthy)]
            ]
            
            # Fisher's Exact Test provides the Odds Ratio (magnitude) and p-value (probability)
            odds_ratio, p_value = stats.fisher_exact(table)
            
            all_results.append({
                'Analysis_Group': cohort_name,
                'SNP_ID': var,
                'Symbol': var_data[sym_c].iloc[0] if sym_c else "N/A",
                'Impact': var_data[imp_c].iloc[0] if imp_c else "N/A",
                'Consequence': var_data[con_c].iloc[0] if con_c else "N/A",
                'Tumour_Freq_%': (count_tumour / n_tumour) * 100,
                'Healthy_Freq_%': (count_healthy / n_healthy) * 100,
                'Odds_Ratio': odds_ratio,
                'P_Value': p_value
            })

    # --- 6. MULTI-TESTING CORRECTION ---
    res_df = pd.DataFrame(all_results)
    final_dfs = []
    # Correcting for False Discovery Rate (FDR) within each group independently
    for group in res_df['Analysis_Group'].unique():
        sub = res_df[res_df['Analysis_Group'] == group].copy()
        _, sub['FDR_P_Value'], _, _ = multipletests(sub['P_Value'], method='fdr_bh')
        final_dfs.append(sub)
    
    final_res = pd.concat(final_dfs).sort_values(['Analysis_Group', 'P_Value'])
    
    # Save to Excel with separate tabs for each cohort analysis
    with pd.ExcelWriter(output_xlsx) as writer:
        for group in cohort_list:
            if group in final_res['Analysis_Group'].unique():
                final_res[final_res['Analysis_Group'] == group].to_excel(writer, sheet_name=group, index=False)
    
    # --- 7. VISUALISATION: VOLCANO PLOTS ---
    final_res['-log10_p'] = -np.log10(final_res['P_Value'].replace(0, 1e-20))
    # Cap OR for plotting purposes to avoid outliers stretching the graph
    max_or_plot = min(final_res['Odds_Ratio'].replace(np.inf, np.nan).max(), 25)
    
    plt.figure(figsize=(15, 8))
    g = sns.FacetGrid(final_res, col="Analysis_Group", hue="Impact", palette="Set1", height=6, aspect=1)
    g.map(sns.scatterplot, "Odds_Ratio", "-log10_p", s=120, edgecolor='black', alpha=0.7)
    
    # Adding significance thresholds for visual aid
    for ax in g.axes.flat:
        ax.axhline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
        ax.set_xlim(-1, max_or_plot + 5)
        ax.set_xlabel("Odds Ratio (Risk Factor)")
        ax.set_ylabel("-log10(P-Value)")
        ax.grid(True, alpha=0.3)

    g.add_legend(title="VEP Impact")
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    
    print(f"\n>>> SUCCESS: SNP Enrichment complete.")
    print(f">>> Results: {output_xlsx}")
    print(f">>> Plot: {output_plot}")

if __name__ == "__main__":
    run_snp_association_analysis()
