import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# --- 1. PATH CONFIGURATION ---
base_path = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_dir = base_path


def identify_snps_pipeline():
    """
    Main pipeline: Identifies established SNPs and calculates frequencies.
    Now includes Consequence information for each SNP.
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
    
    # Updated groupby to include Consequence
    summary = df_snps.groupby([
        'Variant_ID', 'Cohort', 'Tissue', 'SYMBOL', 
        'Consequence', 'IMPACT', 'gnomADe_NFE_AF'
    ]).agg({
        'Sample': 'nunique'
    }).reset_index()
    summary.rename(columns={'Sample': 'Carrier_Count'}, inplace=True)

    summary['Total_Group_Samples'] = summary.apply(
        lambda x: total_samples_dict.get((x['Cohort'], x['Tissue']), 0), axis=1
    )
    summary['Frequency_%'] = (summary['Carrier_Count'] / summary['Total_Group_Samples']) * 100

    # --- 5. EXCEL GENERATION (Pivot with 0-Fill) ---
    master_pivot = summary.pivot_table(
        index=['Variant_ID', 'SYMBOL', 'Consequence', 'IMPACT', 'gnomADe_NFE_AF'],
        columns=['Cohort', 'Tissue'],
        values='Frequency_%',
        aggfunc='first'
    )

    # Descriptive column names
    master_pivot.columns = [f"{col[0]}_{col[1]}_Frequency_%" for col in master_pivot.columns.values]
    
    # Fill frequency gaps with 0
    master_pivot = master_pivot.fillna(0)

    master_pivot.reset_index(inplace=True)
    master_pivot.rename(columns={'gnomADe_NFE_AF': 'gnomAD_NFE_AF'}, inplace=True)
    
    # Reorder columns to place Consequence after SYMBOL for better readability
    cols = master_pivot.columns.tolist()
    # Move Consequence to position after SYMBOL if it exists
    if 'SYMBOL' in cols and 'Consequence' in cols:
        symbol_idx = cols.index('SYMBOL')
        cols.remove('Consequence')
        cols.insert(symbol_idx + 1, 'Consequence')
        master_pivot = master_pivot[cols]
    
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
                        hue='Tissue', s=130, 
                        palette={'Tumour': '#e74c3c', 'Healthy': '#3498db'}, 
                        alpha=0.7)
        
        plt.plot([0, 1], [0, 100], '--', color='grey', alpha=0.4, label="Global Baseline")
        plt.title(f"SNP Frequency Benchmarking: {cohort} Cohort", fontsize=15, fontweight='bold')
        plt.xlabel("gnomAD NFE Frequency"), plt.ylabel("Study Frequency (%)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"11_IdentityPlot_{cohort}.png"), dpi=300)
        plt.close()

    # New visualisation: SNP distribution by Consequence type
    plt.figure(figsize=(12, 6))
    consequence_counts = master_pivot['Consequence'].value_counts().head(10)
    sns.barplot(x=consequence_counts.values, y=consequence_counts.index, palette='mako')
    plt.title("Top 10 Most Common SNP Consequences", fontsize=14, fontweight='bold')
    plt.xlabel("Number of Unique SNPs")
    plt.ylabel("Consequence Type")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "11_SNP_Consequences_Distribution.png"), dpi=300)
    plt.close()

    print(f"\n>>> PIPELINE SUCCESSFUL. Found {len(master_pivot)} unique SNPs.")
    print(f"    Consequence types identified: {master_pivot['Consequence'].nunique()}")

if __name__ == "__main__":
    identify_snps_pipeline()