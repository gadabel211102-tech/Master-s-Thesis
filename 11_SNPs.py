import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# --- 1. PATH CONFIGURATION ---
# Setting the directory for input data and where all results/plots will be saved
base_path = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_dir = base_path

def identify_snps_pipeline():
    """
    Main pipeline to identify established SNPs, calculate their frequencies 
    across cohorts/tissues, and generate visual reports.
    """
    if not os.path.exists(input_file):
        print(f"ERROR: File not found: {input_file}")
        return

    print(">>> Script 11: Executing Full SNP Annotation and Visualisation Pipeline...")
    
    # --- 2. DATA LOADING ---
    # Reading the VEP-annotated data from the 'Biological_Annotations' sheet
    df = pd.read_excel(input_file, sheet_name='Biological_Annotations')
    
    # --- 3. APPLY SNP DEFINITION FILTER ---
    # A 'SNP' is defined as a variant with >1% frequency in the general population.
    # We filter using the gnomAD Non-Finnish European (NFE) Allele Frequency.
    df_snps = df[df['gnomADe_NFE_AF'] > 0.01].copy()
    print(f"   [FILTER] Retained variants with >1% population frequency.")

    # --- 4. ROBUST ID CREATION (The 'No-Skip' Logic) ---
    # We create a unique 'Variant_ID'. 
    # Priority 1: Use the rsID if available.
    # Priority 2: If no rsID, combine Gene Symbol and Protein change.
    # Priority 3: Final fallback to a position-based ID to ensure no data is lost.
    df_snps['rsID'] = df_snps['Existing_variation'].astype(str).str.extract(r'(rs\d+)')
    df_snps['Variant_ID'] = df_snps['rsID'].fillna(df_snps['SYMBOL'].astype(str) + ":" + df_snps['HGVSp'].astype(str))
    
    # Formatting columns for consistency and British spelling
    df_snps['Cohort'] = df_snps['Cohort'].astype(str).str.strip()
    df_snps['Tissue'] = df_snps['Tissue'].astype(str).str.strip().replace({
        'Tumor': 'Tumour', 
        'Normal': 'Healthy', 
        'Control': 'Healthy'
    })

    # --- 5. FREQUENCY CALCULATIONS ---
    # Determining denominators: the total number of unique samples in each group
    total_samples = df_snps.groupby(['Cohort', 'Tissue'])['Sample'].nunique().to_dict()
    
    # Aggregating detections: how many people carry each SNP in each group?
    summary = df_snps.groupby(['Variant_ID', 'Cohort', 'Tissue', 'SYMBOL', 'IMPACT', 'gnomADe_NFE_AF']).agg({
        'Sample': 'nunique'
    }).reset_index()
    summary.rename(columns={'Sample': 'Carrier_Count'}, inplace=True)

    # Calculating percentage frequency (%) for each SNP per group
    summary['Frequency_%'] = summary.apply(
        lambda x: (x['Carrier_Count'] / total_samples.get((x['Cohort'], x['Tissue']), 1)) * 100, axis=1
    )

    # --- 6. MASTER SUMMARY EXCEL (The Pivot Table) ---
    # Creating a wide-format table showing each SNP on one row with across-group comparisons
    master_pivot = summary.pivot_table(
        index=['Variant_ID', 'SYMBOL', 'IMPACT', 'gnomADe_NFE_AF'],
        columns=['Cohort', 'Tissue'],
        values='Frequency_%'
    ).reset_index()

    # Flattening multi-level headers for Excel compatibility
    master_pivot.columns = [
        f"{col[0]}_{col[1]}".strip('_') if isinstance(col, tuple) else col 
        for col in master_pivot.columns.values
    ]
    
    master_pivot.rename(columns={'gnomADe_NFE_AF': 'gnomAD_NFE_Frequency'}, inplace=True)
    master_pivot.to_excel(os.path.join(output_dir, "11_Master_Unique_SNP_Summary.xlsx"), index=False)
    print("   [OUTPUT] Master Summary Excel generated.")

    # --- 7. BAR PLOT: TOTAL UNIQUE SNPS PER GENE ---
    # Visualising the variety of SNPs found in each gene targeted by the panel
    gene_counts = master_pivot.groupby('SYMBOL')['Variant_ID'].nunique().sort_values(ascending=False).reset_index()
    gene_counts.columns = ['Gene', 'Unique_SNPs']

    plt.figure(figsize=(10, 6))
    sns.barplot(data=gene_counts, x='Gene', y='Unique_SNPs', palette='viridis')
    plt.title("Total Unique SNPs Identified per Gene (>1% Frequency)", fontsize=15, fontweight='bold')
    plt.ylabel("Number of Unique SNPs", fontsize=12)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "11_BarPlot_Total_SNPs_Per_Gene.png"), dpi=300)
    plt.close()

    # --- 8. BAR PLOT: SPLIT BY COHORT AND TISSUE (The requested modification) ---
    # Counting unique SNPs detected in each subgroup (Breast Tumour, Breast Healthy, etc.)
    subgroup_counts = summary.groupby(['Cohort', 'Tissue', 'SYMBOL'])['Variant_ID'].nunique().reset_index()
    subgroup_counts.columns = ['Cohort', 'Tissue', 'Gene', 'Unique_SNPs']
    
    # Creating a helper column for plotting the nested groups
    subgroup_counts['Group'] = subgroup_counts['Cohort'] + " (" + subgroup_counts['Tissue'] + ")"

    plt.figure(figsize=(14, 8))
    # Using 'hue' to show both Cohort and Tissue health status side-by-side
    sns.barplot(data=subgroup_counts, x='Gene', y='Unique_SNPs', hue='Group', palette='RdBu_r')
    plt.title("Unique SNPs Identified per Gene split by Cohort and Tissue", fontsize=15, fontweight='bold')
    plt.ylabel("Number of Unique SNPs Detected", fontsize=12)
    plt.xlabel("Gene Symbol", fontsize=12)
    plt.xticks(rotation=45)
    plt.legend(title='Subgroup', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "11_BarPlot_SNPs_By_Tissue_Subgroup.png"), dpi=300)
    plt.close()
    print("   [OUTPUT] Tissue-stratified Bar Plot generated.")

    # --- 9. IDENTITY PLOTS (Population Comparison) ---
    # Creating scatter plots to check if study frequencies match population baselines
    sns.set_style("whitegrid")
    for cohort in summary['Cohort'].unique():
        plt.figure(figsize=(11, 9))
        df_c = summary[summary['Cohort'] == cohort]
        
        # Colouring by Tissue (Tumour vs Healthy)
        sns.scatterplot(data=df_c, x='gnomADe_NFE_AF', y='Frequency_%', 
                        hue='Tissue', style='IMPACT', s=150, 
                        palette={'Tumour': '#e74c3c', 'Healthy': '#3498db'}, alpha=0.7)
        
        # Circling GSDMB specific targets
        gsdmb = df_c[df_c['SYMBOL'] == 'GSDMB']
        if not gsdmb.empty:
            plt.scatter(gsdmb['gnomADe_NFE_AF'], gsdmb['Frequency_%'], 
                        facecolors='none', edgecolors='black', s=450, linewidths=3, label='GSDMB Target')

        # The Identity line (Expected frequency vs Observed)
        plt.plot([0, 1], [0, 100], '--', color='grey', alpha=0.5, label="Population Baseline")
        
        plt.title(f"Allele Frequency Comparison: {cohort} Cohort", fontsize=16, fontweight='bold')
        plt.xlabel("European Population Frequency (gnomAD NFE)", fontsize=12)
        plt.ylabel("Observed Study Frequency (%)", fontsize=12)
        plt.xlim(-0.05, 1.05)
        plt.ylim(-5, 105)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"11_IdentityPlot_{cohort}.png"), dpi=300)
        plt.close()

    print(f"\n>>> PIPELINE SUCCESSFUL.")
    print(f">>> Found {len(master_pivot)} unique SNPs across all groups.")
    print(f">>> Outputs saved in: {output_dir}")

if __name__ == "__main__":
    identify_snps_pipeline()
