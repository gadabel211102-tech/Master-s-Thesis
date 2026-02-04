import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# --- PATH CONFIGURATION ---
base_path = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_dir = base_path

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def master_pipeline():
    if not os.path.exists(input_file):
        print(f"ERROR: File not found at: {input_file}")
        return

    print(">>> Starting Ancestry-Matched SNP Discovery Pipeline (gnomAD NFE)...")
    
    # 1. Load data
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')
    
    # 2. Data Cleaning & ID Creation
    df['rsID'] = df['Existing_variation'].astype(str).str.extract(r'(rs\d+)')
    df['Variant_ID'] = df['rsID'].fillna(df['SYMBOL'] + ":" + df['HGVSp'].astype(str))
    
    df['Cohort'] = df['Cohort'].astype(str).str.strip()
    df['Tissue'] = df['Tissue'].astype(str).str.strip().replace('Tumor', 'Tumour')

    # 3. Denominators
    total_samples = df.groupby(['Cohort', 'Tissue'])['Sample'].nunique().to_dict()

   # 4. Global SNP Aggregation 
    nfe_col = 'gnomADe_NFE_AF' 
    # Added dropna=False so variants without gnomAD data are NOT deleted
    snp_master = df.groupby(['Variant_ID', 'SYMBOL', 'Cohort', 'Tissue', 'IMPACT', 'Consequence', nfe_col], 
                            dropna=False).size().reset_index(name='Count')

    # 5. Study Frequency (%) Calculation
    def get_percentage(row):
        denom = total_samples.get((row['Cohort'], row['Tissue']), 1)
        return (row['Count'] / denom) * 100

    snp_master['Study_Frequency_%'] = snp_master.apply(get_percentage, axis=1)

   # 6. APPLY 1% FILTER (Based on Population, not Study) & Novel Variants
    # Condition 1: Common in European Population (gnomAD NFE > 0.01)
    # Condition 2: Novel/Rare (gnomAD NFE is empty/NaN)
    pop_filter = (snp_master[nfe_col] > 0.01) | (snp_master[nfe_col].isna())
    
    final_snps = snp_master[pop_filter].copy()
    
    # Create the percentage column for plotting
    final_snps['gnomAD_perc'] = final_snps[nfe_col] * 100
    # For Novel variants, set gnomAD_perc to 0 so they show up on the Y-axis in plots
    final_snps['gnomAD_perc'] = final_snps['gnomAD_perc'].fillna(0)

    # Debugging print to help you understand the counts
    n_common = (snp_master[nfe_col] > 0.01).sum()
    n_novel = (snp_master[nfe_col].isna()).sum()
    print(f"Total processed: {len(snp_master)}")
    print(f"Common SNPs kept (>1%): {n_common}")
    print(f"Novel/Rare kept (NaN): {n_novel}")

    # 7. Pivot for Comparative Table
    pivot_table = final_snps.pivot_table(
        index=['Variant_ID', 'SYMBOL', 'IMPACT', 'Consequence', nfe_col],
        columns=['Cohort', 'Tissue'],
        values='Study_Frequency_%',
        fill_value=0
    )
    
    pivot_table.columns = [f"{c[0]}_{c[1]}_%" for c in pivot_table.columns]
    pivot_table = pivot_table.reset_index()

    # 8. EXPORT TO EXCEL
    xlsx_out = os.path.join(output_dir, "13_MASTER_SNP_Catalogue_NFE_Corrected.xlsx")
    with pd.ExcelWriter(xlsx_out) as writer:
        pivot_table.to_excel(writer, sheet_name='All_Genes_SNPs', index=False)
        gsdmb_only = pivot_table[pivot_table['SYMBOL'] == 'GSDMB']
        gsdmb_only.to_excel(writer, sheet_name='GSDMB_Only_SNPs', index=False)
    
    print(f"Excel catalogue generated: {xlsx_out}")

    # --- PLOTTING SECTION ---
    sns.set_context("talk") 
    
    # Figure 14: Top Genes
    plt.figure(figsize=(12, 7))
    top_genes = final_snps['SYMBOL'].value_counts().head(15)
    sns.barplot(x=top_genes.values, y=top_genes.index, palette='viridis')
    plt.title("Top Genes by SNP Diversity (NFE Baseline)", fontsize=16, fontweight='bold')
    plt.xlabel("Number of Unique Variants Detected")
    plt.ylabel("Gene Symbol")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "14_Top_Genes_SNP_Density_NFE.png"), dpi=300)
    plt.close()

    # Figure 15: Frequency Comparison (NFE Baseline)
    cohorts = final_snps['Cohort'].unique()
    
    for cohort in cohorts:
        plt.figure(figsize=(12, 9))
        df_cohort = final_snps[final_snps['Cohort'] == cohort]
        
        sns.scatterplot(
            data=df_cohort, 
            x='gnomAD_perc', 
            y='Study_Frequency_%', 
            hue='Tissue', 
            style='IMPACT', 
            alpha=0.6,
            s=120,
            palette={'Tumour': '#e74c3c', 'Normal': '#3498db'},
            zorder=2
        )
        
        # HIGHLIGHT GSDMB
        gsdmb_cohort = df_cohort[df_cohort['SYMBOL'] == 'GSDMB']
        if not gsdmb_cohort.empty:
            plt.scatter(
                gsdmb_cohort['gnomAD_perc'], 
                gsdmb_cohort['Study_Frequency_%'], 
                facecolors='none', 
                edgecolors='black', 
                s=350, 
                linewidths=3, 
                label='GSDMB Target',
                zorder=10
            )
        
        plt.plot([0, 100], [0, 100], '--', color='grey', alpha=0.5, label="Identity Line")
        
        # UPDATED: Axis labels now specify European (NFE) frequency
        plt.title(f"Allele Frequency Comparison: {cohort} Cohort", fontsize=18, fontweight='bold')
        plt.xlabel("European Population Frequency (gnomAD NFE %)", fontsize=14)
        plt.ylabel("Observed Study Frequency (%)", fontsize=14)
        
        plt.xlim(-5, 105)
        plt.ylim(-5, 105)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"15_{cohort}_NFE_Comparison.png"), dpi=300)
        plt.close()

    print(f"\n>>> PIPELINE FINISHED: Ancestry-matched analysis complete.")

if __name__ == "__main__":
    master_pipeline()
