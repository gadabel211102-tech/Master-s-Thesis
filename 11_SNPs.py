import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from scipy.stats import chisquare

# --- 1. PATH CONFIGURATION ---
base_path = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_dir = base_path

def calculate_hwe_qc(row):
    """Secondary QC check: Calculates HWE p-value to flag potential technical artefacts."""
    n = row['Total_Group_Samples']
    if n < 10 or pd.isna(n): 
        return np.nan # Insufficient power for small cohorts
    
    obs_carriers = row['Carrier_Count']
    obs_non_carriers = n - obs_carriers
    
    # Estimate allele frequency (q)
    q = obs_carriers / n
    p = 1 - q
    
    # Expected counts under HWE (p^2 + 2pq + q^2 = 1)
    exp_carriers = (1 - (p**2)) * n 
    exp_non_carriers = (p**2) * n
    
    # Chi-square is unreliable if expected counts are very low
    if exp_carriers < 5 or exp_non_carriers < 5: 
        return 1.0
    
    _, p_val = chisquare([obs_carriers, obs_non_carriers], [exp_carriers, exp_non_carriers])
    return p_val

def identify_snps_pipeline():
    """
    Main pipeline: Identifies established SNPs, calculates frequencies, 
    and adds HWE status as a secondary QC flag with 0-filling for missing data.
    """
    if not os.path.exists(input_file):
        print(f"ERROR: File not found: {input_file}")
        return

    print(">>> Script 11: Identifying SNPs and Generating Frequency Reports...")
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')
    
    # --- 2. SNP IDENTIFICATION & CLEANING ---
    # Focus on variants with >1% population frequency
    df_snps = df[df['gnomADe_NFE_AF'] > 0.01].copy()
    
    # Creating robust Variant IDs
    df_snps['rsID'] = df_snps['Existing_variation'].astype(str).str.extract(r'(rs\d+)')
    df_snps['Variant_ID'] = df_snps['rsID'].fillna(df_snps['SYMBOL'].astype(str) + ":" + df_snps['HGVSp'].astype(str))
    
    # Standardising labels with British spelling
    df_snps['Cohort'] = df_snps['Cohort'].astype(str).str.strip()
    df_snps['Tissue'] = df_snps['Tissue'].astype(str).str.strip().replace({
        'Tumor': 'Tumour', 'Normal': 'Healthy', 'Control': 'Healthy'
    })

    # --- 3. FREQUENCY CALCULATIONS ---
    total_samples_dict = df_snps.groupby(['Cohort', 'Tissue'])['Sample'].nunique().to_dict()
    
    summary = df_snps.groupby(['Variant_ID', 'Cohort', 'Tissue', 'SYMBOL', 'IMPACT', 'gnomADe_NFE_AF']).agg({
        'Sample': 'nunique'
    }).reset_index()
    summary.rename(columns={'Sample': 'Carrier_Count'}, inplace=True)

    summary['Total_Group_Samples'] = summary.apply(
        lambda x: total_samples_dict.get((x['Cohort'], x['Tissue']), 0), axis=1
    )
    summary['Frequency_%'] = (summary['Carrier_Count'] / summary['Total_Group_Samples']) * 100

    # --- 4. SECONDARY QC: HWE FLAGGING ---
    summary['HWE_P'] = summary.apply(calculate_hwe_qc, axis=1)
    summary['QC_Flag'] = "PASS"
    # Only flag deviations in Healthy (control) tissue
    summary.loc[(summary['Tissue'] == 'Healthy') & (summary['HWE_P'] < 0.05), 'QC_Flag'] = "HWE_DEVIATION"

    # --- 5. EXCEL GENERATION (Pivot with 0-Fill) ---
    master_pivot = summary.pivot_table(
        index=['Variant_ID', 'SYMBOL', 'IMPACT', 'gnomADe_NFE_AF'],
        columns=['Cohort', 'Tissue'],
        values=['Frequency_%', 'QC_Flag'],
        aggfunc='first'
    )

    # Descriptive column names
    master_pivot.columns = [f"{col[1]}_{col[2]}_{col[0]}" for col in master_pivot.columns.values]
    
    # Fill frequency gaps with 0
    freq_cols = [c for c in master_pivot.columns if 'Frequency_%' in c]
    master_pivot[freq_cols] = master_pivot[freq_cols].fillna(0)
    
    # Ensure QC_Flag columns are readable
    qc_cols = [c for c in master_pivot.columns if 'QC_Flag' in c]
    master_pivot[qc_cols] = master_pivot[qc_cols].fillna("PASS")

    master_pivot.reset_index(inplace=True)
    master_pivot.rename(columns={'gnomADe_NFE_AF': 'gnomAD_NFE_AF'}, inplace=True)
    
    master_pivot.to_excel(os.path.join(output_dir, "11_Master_Unique_SNP_Summary.xlsx"), index=False)

    # --- 6. VISUALISATION ---
    # Bar plot: Total unique SNPs per gene
    gene_counts = master_pivot.groupby('SYMBOL')['Variant_ID'].nunique().sort_values(ascending=False).reset_index()
    plt.figure(figsize=(10, 6))
    sns.barplot(data=gene_counts, x='SYMBOL', y='Variant_ID', palette='viridis')
    plt.title("Total Unique SNPs Identified per Gene", fontsize=14, fontweight='bold')
    plt.ylabel("Unique SNP Count"), plt.xticks(rotation=45), plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "11_BarPlot_SNPs_Per_Gene.png"), dpi=300)
    plt.close()

    # Identity Plots: Population Benchmarking
    for cohort in summary['Cohort'].unique():
        plt.figure(figsize=(10, 8))
        df_c = summary[summary['Cohort'] == cohort]
        
        sns.scatterplot(data=df_c, x='gnomADe_NFE_AF', y='Frequency_%', 
                        hue='Tissue', style='QC_Flag', s=130, 
                        palette={'Tumour': '#e74c3c', 'Healthy': '#3498db'}, 
                        markers={'PASS': 'o', 'HWE_DEVIATION': 'X'}, alpha=0.7)
        
        plt.plot([0, 1], [0, 100], '--', color='grey', alpha=0.4, label="Global Baseline")
        plt.title(f"SNP Frequency Benchmarking: {cohort} Cohort", fontsize=15, fontweight='bold')
        plt.xlabel("gnomAD NFE Frequency"), plt.ylabel("Study Frequency (%)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"11_IdentityPlot_{cohort}.png"), dpi=300)
        plt.close()

    print(f"\n>>> PIPELINE SUCCESSFUL. Found {len(master_pivot)} unique SNPs.")

if __name__ == "__main__":
    identify_snps_pipeline()
