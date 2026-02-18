import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import chi2_contingency
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

def test_genotype_distribution(var_data, tis_c, sam_c):
    """
    Test if genotype distribution (AA/Aa/aa) differs between tumour and healthy.
    This is more powerful than just testing carrier vs non-carrier status.
    
    Args:
        var_data: DataFrame with variant data for one SNP
        tis_c: Tissue column name
        sam_c: Sample column name
        
    Returns:
        Dictionary with chi-square test results
    """
    # For SNPs, we need genotype info (homozygous ref, het, homozygous alt)
    # Since we don't have explicit GT column in the annotations,
    # we'll use a simplified model based on allele frequency expectations
    
    # Count samples in each group
    tumour_samples = var_data[var_data[tis_c] == 'Tumour'][sam_c].unique()
    healthy_samples = var_data[var_data[tis_c] == 'Healthy'][sam_c].unique()
    
    # For now, return None - this would require GT field from VCF
    # In practice, you'd extract genotypes during annotation
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
    
    # --- 7. VISUALISATION: ENHANCED VOLCANO PLOTS ---
    final_res['-log10_p'] = -np.log10(final_res['P_Value'].replace(0, 1e-10))
    final_res['log2_OR'] = np.log2(final_res['Odds_Ratio'].replace([0, np.inf], [0.001, 1000]))
    
    # Identify top significant SNPs for labeling
    final_res['Is_Significant'] = final_res['FDR_P_Value'] < 0.05
    final_res['Label'] = ''
    
    for group in final_res['Analysis_Group'].unique():
        group_mask = final_res['Analysis_Group'] == group
        sig_mask = (final_res['FDR_P_Value'] < 0.05) & group_mask
        
        if sig_mask.any():
            # Label top 5 most significant in each group
            top_indices = final_res[sig_mask].nsmallest(5, 'FDR_P_Value').index
            final_res.loc[top_indices, 'Label'] = final_res.loc[top_indices, 'SNP_ID']
    
    # Create enhanced plots
    plt.figure(figsize=(18, 10))
    g = sns.FacetGrid(final_res, col="Analysis_Group", hue="Impact", palette="Set1", 
                      height=6, aspect=1.3, despine=False)
    g.map(sns.scatterplot, "Odds_Ratio", "-log10_p", s=150, edgecolor='black', alpha=0.8)
    
    # Add reference lines and labels
    for ax, group_name in zip(g.axes.flat, final_res['Analysis_Group'].unique()):
        # Significance threshold
        ax.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.6, 
                  linewidth=2, label='p=0.05')
        
        # Odds ratio = 1 (no effect)
        ax.axvline(1, color='blue', linestyle='--', alpha=0.6, 
                  linewidth=2, label='OR=1 (no effect)')
        
        # Shade regions
        ax.axhspan(-np.log10(0.05), ax.get_ylim()[1], alpha=0.1, color='red', 
                  label='Significant (p<0.05)')
        ax.axvspan(0, 1, alpha=0.05, color='blue', label='Protective (OR<1)')
        ax.axvspan(1, ax.get_xlim()[1], alpha=0.05, color='red', label='Risk (OR>1)')
        
        # Add labels for top SNPs
        group_data = final_res[final_res['Analysis_Group'] == group_name]
        labeled = group_data[group_data['Label'] != '']
        
        for _, row in labeled.iterrows():
            ax.annotate(row['Label'], 
                       xy=(row['Odds_Ratio'], row['-log10_p']),
                       xytext=(10, 10), textcoords='offset points',
                       fontsize=8, fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0',
                                     color='black', lw=1))
        
        ax.set_xlabel("Odds Ratio (log scale recommended)", fontweight='bold')
        ax.set_ylabel("-log10(P-Value)", fontweight='bold')
        ax.set_xscale('log')
        ax.grid(True, linestyle=':', alpha=0.4)
        ax.legend(loc='best', fontsize=8)

    g.add_legend(title="VEP Impact", bbox_to_anchor=(1.02, 0.5), loc='center left')
    plt.subplots_adjust(top=0.90, right=0.95)
    g.fig.suptitle("SNP Association Analysis: Tumour vs. Healthy Cohorts\n(Enhanced Volcano Plots)", 
                   fontsize=18, fontweight='bold')
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Also create a version with log2(OR) on x-axis for better visualization
    output_plot_log2 = output_plot.replace('.png', '_log2OR.png')
    plt.figure(figsize=(18, 10))
    g2 = sns.FacetGrid(final_res, col="Analysis_Group", hue="Impact", palette="Set1", 
                       height=6, aspect=1.3)
    g2.map(sns.scatterplot, "log2_OR", "-log10_p", s=150, edgecolor='black', alpha=0.8)
    
    for ax in g2.axes.flat:
        ax.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.6, linewidth=2)
        ax.axvline(0, color='blue', linestyle='--', alpha=0.6, linewidth=2, label='OR=1')
        ax.set_xlabel("log2(Odds Ratio)", fontweight='bold')
        ax.set_ylabel("-log10(P-Value)", fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.4)
    
    g2.add_legend(title="VEP Impact")
    plt.subplots_adjust(top=0.90)
    g2.fig.suptitle("SNP Association Analysis: Log2 Odds Ratios", fontsize=18, fontweight='bold')
    plt.savefig(output_plot_log2, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n>>> SUCCESS: All {len(unique_snps)} SNPs processed.")
    print(f">>> Results saved to {output_xlsx}")
    print(f">>> Volcano plots saved:")
    print(f"    - {output_plot}")
    print(f"    - {output_plot_log2}")

if __name__ == "__main__":
    run_snp_association_analysis()