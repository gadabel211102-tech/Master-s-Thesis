import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# --- PATH CONFIGURATION ---
# Path updated for your environment
base_path = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_dir = base_path

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def master_pipeline():
    if not os.path.exists(input_file):
        print(f"ERROR: File not found at: {input_file}")
        return

    print(">>> Starting Global SNP Discovery Pipeline (1% Filter)...")
    
    # 1. Load data
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')
    
    # 2. Data Cleaning & ID Creation
    # Extract rsID from 'Existing_variation'
    df['rsID'] = df['Existing_variation'].astype(str).str.extract(r'(rs\d+)')
    # Create Unique Variant ID: prioritising rsID, then Gene:Protein change
    df['Variant_ID'] = df['rsID'].fillna(df['SYMBOL'] + ":" + df['HGVSp'].astype(str))
    
    # Use British spelling for cohorts and tissues
    df['Cohort'] = df['Cohort'].astype(str).str.strip()
    df['Tissue'] = df['Tissue'].astype(str).str.strip().replace('Tumor', 'Tumour')

    # 3. Denominators (Unique samples per group for percentage calculation)
    total_samples = df.groupby(['Cohort', 'Tissue'])['Sample'].nunique().to_dict()
    print(f"Unique samples per group: {total_samples}")

    # 4. Global SNP Aggregation (Grouping occurrences into unique variants)
    snp_master = df.groupby(['Variant_ID', 'SYMBOL', 'Cohort', 'Tissue', 'IMPACT', 'Consequence', 'gnomADe_AF']).size().reset_index(name='Count')

    # 5. Study Frequency (%) Calculation
    def get_percentage(row):
        denom = total_samples.get((row['Cohort'], row['Tissue']), 1)
        return (row['Count'] / denom) * 100

    snp_master['Study_Frequency_%'] = snp_master.apply(get_percentage, axis=1)

    # 6. APPLY 1% FILTER
    # This identifies variants common enough to be considered SNPs in your cohort
    final_snps = snp_master[snp_master['Study_Frequency_%'] > 1].copy()
    final_snps['gnomAD_perc'] = final_snps['gnomADe_AF'] * 100

    # 7. Pivot for Comparative Table (One row per variant, columns for each group)
    pivot_table = final_snps.pivot_table(
        index=['Variant_ID', 'SYMBOL', 'IMPACT', 'Consequence', 'gnomADe_AF'],
        columns=['Cohort', 'Tissue'],
        values='Study_Frequency_%',
        fill_value=0
    )
    
    # Formatting column headers (e.g., Breast_Tumour_%)
    pivot_table.columns = [f"{c[0]}_{c[1]}_%" for c in pivot_table.columns]
    pivot_table = pivot_table.reset_index()

    # 8. EXPORT TO EXCEL (Two Sheets)
    xlsx_out = os.path.join(output_dir, "13_MASTER_SNP_Catalogue_1percent.xlsx")
    with pd.ExcelWriter(xlsx_out) as writer:
        pivot_table.to_excel(writer, sheet_name='All_Genes_SNPs', index=False)
        # Separate sheet for GSDMB for your main analysis
        gsdmb_only = pivot_table[pivot_table['SYMBOL'] == 'GSDMB']
        gsdmb_only.to_excel(writer, sheet_name='GSDMB_Only_SNPs', index=False)
    
    print(f"Excel catalogue generated: {xlsx_out}")

    # --- PLOTTING SECTION ---
    sns.set_context("talk") 
    
    # Figure 14: Top Genes by SNP Count
    plt.figure(figsize=(12, 7))
    top_genes = final_snps['SYMBOL'].value_counts().head(15)
    sns.barplot(x=top_genes.values, y=top_genes.index, palette='viridis')
    plt.title("Top 15 Genes by SNP Diversity (>1% Frequency)", fontsize=16, fontweight='bold')
    plt.xlabel("Number of Unique Variants Detected")
    plt.ylabel("Gene Symbol")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "14_Top_Genes_SNP_Density.png"), dpi=300)
    plt.close()

    # Figure 15: Frequency Comparison - Separate Plots for Cohorts
    cohorts = final_snps['Cohort'].unique()
    
    for cohort in cohorts:
        plt.figure(figsize=(12, 9))
        df_cohort = final_snps[final_snps['Cohort'] == cohort]
        
        # Plot All SNPs: Colour by Tissue (Tumour vs Normal)
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
        
        # HIGHLIGHT GSDMB (Black target rings)
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
        
        # Identity Line (Diagonal reference)
        plt.plot([0, 100], [0, 100], '--', color='grey', alpha=0.5, label="Identity Line")
        
        plt.title(f"Allele Frequency Comparison: {cohort} Cohort", fontsize=18, fontweight='bold')
        plt.xlabel("Global Population Frequency (gnomAD %)", fontsize=14)
        plt.ylabel("Observed Study Frequency (%)", fontsize=14)
        
        plt.xlim(-5, 105)
        plt.ylim(-5, 105)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"15_{cohort}_Frequency_Comparison.png"), dpi=300)
        plt.close()

    print(f"\n>>> PIPELINE FINISHED")
    print(f"Total Unique SNPs Detected (>1%): {len(pivot_table)}")
    print(f"GSDMB SNPs Detected (>1%): {len(gsdmb_only)}")

if __name__ == "__main__":
    master_pipeline()