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
output_xlsx = os.path.join(base_path, "12_SNP_Enrichment_Results.xlsx")
output_plot = os.path.join(base_path, "12_SNP_Volcano_Plots.png")

def find_col(df, target):
    """Searches for columns regardless of capitalisation or hidden spaces."""
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None

def run_snp_association_analysis():
    print(">>> Starting Script 12: SNP Statistical Enrichment with Haldane-Anscombe Correction...")
    
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        return

    # --- 2. DATA LOADING & CLEANING ---
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')
    
    # Identify key columns
    sym_c = find_col(df, 'SYMBOL')
    imp_c = find_col(df, 'IMPACT')
    con_c = find_col(df, 'CONSEQUENCE')
    tis_c = find_col(df, 'TISSUE')
    sam_c = find_col(df, 'SAMPLE')
    var_c = find_col(df, 'Existing_variation')
    hgv_c = find_col(df, 'HGVSp')
    coh_c = find_col(df, 'COHORT')
    nfe_c = find_col(df, 'gnomADe_NFE_AF')

    # --- 3. THE STRICT SNP FILTER (Matches Script 11) ---
    # Focus only on established SNPs (>1% AF)
    df = df[df[nfe_c] > 0.01].copy()
    
    # Standardising to British spelling
    df[tis_c] = df[tis_c].astype(str).str.strip().replace({'Tumor': 'Tumour', 'Normal': 'Healthy', 'Control': 'Healthy'})
    df[coh_c] = df[coh_c].astype(str).str.strip()
    
    # --- 4. ROBUST ID CREATION ---
    df['rsID'] = df[var_c].astype(str).str.extract(r'(rs\d+)')
    df['Variant_ID'] = df['rsID'].fillna(df[sym_c].astype(str) + ":" + df[hgv_c].astype(str))
    
    cohort_list = ['Global', 'Breast', 'Endometrium']
    all_results = []

    # --- 5. FISHER'S EXACT TEST WITH CORRECTION ---
    for cohort_name in cohort_list:
        if cohort_name == 'Global':
            cohort_df = df.copy()
        else:
            cohort_df = df[df[coh_c].str.contains(cohort_name, case=False, na=False)].copy()
            
        if cohort_df.empty: continue
            
        n_tumour = cohort_df[cohort_df[tis_c] == 'Tumour'][sam_c].nunique()
        n_healthy = cohort_df[cohort_df[tis_c] == 'Healthy'][sam_c].nunique()
        
        if n_tumour == 0 or n_healthy == 0: continue

        unique_snps = cohort_df['Variant_ID'].unique()
        
        for var in unique_snps:
            var_data = cohort_df[cohort_df['Variant_ID'] == var]
            
            count_tumour = var_data[var_data[tis_c] == 'Tumour'][sam_c].nunique()
            count_healthy = var_data[var_data[tis_c] == 'Healthy'][sam_c].nunique()
            
            # Contingency Table cells
            a, b = count_tumour, max(0, n_tumour - count_tumour)
            c, d = count_healthy, max(0, n_healthy - count_healthy)
            
            # P-value calculation (using original integer counts)
            _, p_value = stats.fisher_exact([[int(a), int(b)], [int(c), int(d)]])
            
            # Haldane-Anscombe Correction for Odds Ratio (OR)
            # If any cell is 0, add 0.5 to all cells to avoid infinite OR
            if a == 0 or b == 0 or c == 0 or d == 0:
                or_val = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                note = "Corrected (+0.5)"
            else:
                or_val = (a * d) / (b * c) if (b * c) != 0 else np.inf
                note = "Standard"
            
            all_results.append({
                'Analysis_Group': cohort_name,
                'SNP_ID': var,
                'Symbol': var_data[sym_c].iloc[0],
                'Impact': var_data[imp_c].iloc[0],
                'Tumour_Freq_%': (count_tumour / n_tumour) * 100,
                'Healthy_Freq_%': (count_healthy / n_healthy) * 100,
                'Odds_Ratio': or_val,
                'P_Value': p_value,
                'Calculation_Method': note
            })

    # --- 6. MULTI-TESTING CORRECTION ---
    res_df = pd.DataFrame(all_results)
    final_dfs = []
    for group in res_df['Analysis_Group'].unique():
        sub = res_df[res_df['Analysis_Group'] == group].copy()
        _, sub['FDR_P_Value'], _, _ = multipletests(sub['P_Value'], method='fdr_bh')
        final_dfs.append(sub)
    
    final_res = pd.concat(final_dfs).sort_values(['Analysis_Group', 'P_Value'])
    
    # Save results
    with pd.ExcelWriter(output_xlsx) as writer:
        for group in cohort_list:
            if group in final_res['Analysis_Group'].unique():
                final_res[final_res['Analysis_Group'] == group].to_excel(writer, sheet_name=group, index=False)
    
    # --- 7. VISUALISATION: FULL VOLCANO PLOTS ---
    final_res['-log10_p'] = -np.log10(final_res['P_Value'].replace(0, 1e-10))
    
    # We no longer drop infinite values, as the correction made them finite!
    plt.figure(figsize=(16, 9))
    g = sns.FacetGrid(final_res, col="Analysis_Group", hue="Impact", palette="Set1", height=6, aspect=1.2)
    g.map(sns.scatterplot, "Odds_Ratio", "-log10_p", s=150, edgecolor='black', alpha=0.8)
    
    for ax in g.axes.flat:
        ax.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.6, label='p=0.05')
        ax.set_xlabel("Odds Ratio (Risk Factor)")
        ax.set_ylabel("-log10(P-Value)")
        ax.grid(True, linestyle=':', alpha=0.6)

    g.add_legend(title="VEP Impact")
    plt.subplots_adjust(top=0.85)
    g.fig.suptitle("SNP Association Analysis: Tumour vs. Healthy Cohorts", fontsize=18, fontweight='bold')
    plt.savefig(output_plot, dpi=300)
    
    print(f"\n>>> SUCCESS: All {len(unique_snps)} SNPs processed.")
    print(f">>> Results saved to {output_xlsx}")

if __name__ == "__main__":
    run_snp_association_analysis()
