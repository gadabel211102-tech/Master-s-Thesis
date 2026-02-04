import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# --- PATH CONFIGURATION ---
base_path = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_dir = base_path

def master_pipeline():
    if not os.path.exists(input_file): return

    print(">>> Starting Ancestry-Matched SNP Discovery Pipeline (gnomAD NFE)...")
    
    # 1. Load data
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')
    
    # 2. Data Cleaning & ID Creation
    # Extract the rsID (e.g., rs12345) or build a custom ID if the variant is novel
    df['rsID'] = df['Existing_variation'].astype(str).str.extract(r'(rs\d+)')
    df['Variant_ID'] = df['rsID'].fillna(df['SYMBOL'] + ":" + df['HGVSp'].astype(str))
    
    df['Cohort'] = df['Cohort'].astype(str).str.strip()
    # Standardise British spelling for the report
    df['Tissue'] = df['Tissue'].astype(str).str.strip().replace('Tumor', 'Tumour')

    # 3. Denominators: Count unique samples per group to calculate percentages later
    total_samples = df.groupby(['Cohort', 'Tissue'])['Sample'].nunique().to_dict()

    # 4. Global SNP Aggregation 
    # nfe_col tracks the Allele Frequency in the European population baseline
    nfe_col = 'gnomADe_NFE_AF' 
    snp_master = df.groupby(['Variant_ID', 'SYMBOL', 'Cohort', 'Tissue', 'IMPACT', 'Consequence', nfe_col], 
                            dropna=False).size().reset_index(name='Count')

    # 5. Study Frequency (%) Calculation: Observed frequency in your dataset
    def get_percentage(row):
        denom = total_samples.get((row['Cohort'], row['Tissue']), 1)
        return (row['Count'] / denom) * 100

    snp_master['Study_Frequency_%'] = snp_master.apply(get_percentage, axis=1)

    # 6. APPLY 1% FILTER & IDENTIFY NOVEL VARIANTS
    # We keep variants if they are:
    # A) Common in Europeans (>1% frequency) OR B) Novel (NaN in gnomAD)
    pop_filter = (snp_master[nfe_col] > 0.01) | (snp_master[nfe_col].isna())
    final_snps = snp_master[pop_filter].copy()
    
    # For plotting, convert population frequency to percentage (0-100)
    final_snps['gnomAD_perc'] = final_snps[nfe_col].fillna(0) * 100

    # 7. Pivot for Comparative Table: Create columns for each Cohort_Tissue group
    pivot_table = final_snps.pivot_table(
        index=['Variant_ID', 'SYMBOL', 'IMPACT', 'Consequence', nfe_col],
        columns=['Cohort', 'Tissue'],
        values='Study_Frequency_%',
        fill_value=0
    )
    pivot_table.columns = [f"{c[0]}_{c[1]}_%" for c in pivot_table.columns]
    pivot_table = pivot_table.reset_index()

    # 8. EXPORT TO EXCEL: Generate the final SNP catalogue
    xlsx_out = os.path.join(output_dir, "13_MASTER_SNP_Catalogue_NFE_Corrected.xlsx")
    with pd.ExcelWriter(xlsx_out) as writer:
        pivot_table.to_excel(writer, sheet_name='All_Genes_SNPs', index=False)
        # Create a specific sheet for your target gene: GSDMB
        gsdmb_only = pivot_table[pivot_table['SYMBOL'] == 'GSDMB']
        gsdmb_only.to_excel(writer, sheet_name='GSDMB_Only_SNPs', index=False)

    # --- PLOTTING SECTION ---
    sns.set_context("talk") 
    
    # Figure 14: Bar plot of top 15 genes with the most unique variants detected
    plt.figure(figsize=(12, 7))
    top_genes = final_snps['SYMBOL'].value_counts().head(15)
    sns.barplot(x=top_genes.values, y=top_genes.index, palette='viridis')
    plt.title("Top Genes by SNP Diversity (NFE Baseline)", fontsize=16, fontweight='bold')
    plt.savefig(os.path.join(output_dir, "14_Top_Genes_SNP_Density_NFE.png"), dpi=300)

    # Figure 15: Frequency Comparison (Study vs. Population)
    # This plot helps identify variants that are "enriched" in your study
    for cohort in final_snps['Cohort'].unique():
        plt.figure(figsize=(12, 9))
        df_cohort = final_snps[final_snps['Cohort'] == cohort]
        
        # Scatter plot: X = European pop freq, Y = Your study freq
        sns.scatterplot(
            data=df_cohort, 
            x='gnomAD_perc', 
            y='Study_Frequency_%', 
            hue='Tissue', 
            style='IMPACT', 
            palette={'Tumour': '#e74c3c', 'Normal': '#3498db'},
            s=120, alpha=0.6, zorder=2
        )
        
        # SPECIAL HIGHLIGHT: Draw a thick circle around GSDMB variants
        gsdmb_cohort = df_cohort[df_cohort['SYMBOL'] == 'GSDMB']
        if not gsdmb_cohort.empty:
            plt.scatter(
                gsdmb_cohort['gnomAD_perc'], 
                gsdmb_cohort['Study_Frequency_%'], 
                facecolors='none', edgecolors='black', 
                s=350, linewidths=3, label='GSDMB Target', zorder=10
            )
        
        # Identity Line: Variants on this line have identical study/pop frequencies
        plt.plot([0, 100], [0, 100], '--', color='grey', alpha=0.5, label="Identity Line")
        
        plt.title(f"Allele Frequency Comparison: {cohort} Cohort", fontsize=18, fontweight='bold')
        plt.xlabel("European Population Frequency (gnomAD NFE %)")
        plt.ylabel("Observed Study Frequency (%)")
        plt.xlim(-5, 105); plt.ylim(-5, 105)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig(os.path.join(output_dir, f"15_{cohort}_NFE_Comparison.png"), dpi=300)
        plt.close()

    print(f"\n>>> PIPELINE FINISHED: Ancestry-matched analysis complete.")
