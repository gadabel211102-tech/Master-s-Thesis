import pandas as pd
import os
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- PATH CONFIGURATION ---
base_path = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
# Filtered SNP list from previous ancestry-matching step
input_file = os.path.join(base_path, "13_MASTER_SNP_Catalogue_NFE_Corrected.xlsx")
# Raw data used to calculate the denominators (total sample sizes)
raw_data_for_denoms = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_dir = base_path

def run_enrichment_analysis():
    print(">>> Starting Statistical Enrichment (Population SNPs >1% and Novel)...")
    
    # 1. Load the "Shopping List" of variants to test
    if not os.path.exists(input_file):
        print(f"ERROR: Filtered file not found at {input_file}.")
        return
    df_snps = pd.read_excel(input_file, sheet_name='All_Genes_SNPs')
    
    # 2. Load Raw Data for frequency counting
    df_raw = pd.read_excel(raw_data_for_denoms, sheet_name='Biological_Annotations')
    
    # Clean and standardise the raw data to match the SNP catalogue
    df_raw['rsID'] = df_raw['Existing_variation'].astype(str).str.extract(r'(rs\d+)')
    df_raw['Variant_ID'] = df_raw['rsID'].fillna(df_raw['SYMBOL'] + ":" + df_raw['HGVSp'].astype(str))
    df_raw['Cohort'] = df_raw['Cohort'].astype(str).str.strip()
    df_raw['Tissue'] = df_raw['Tissue'].astype(str).str.strip().replace('Tumor', 'Tumour')

    # 3. Denominators: Calculate total unique samples per Cohort and Tissue
    # This is critical for the "No Variant" cells in the 2x2 contingency table
    sample_counts = df_raw.groupby(['Cohort', 'Tissue'])['Sample'].nunique().to_dict()

    results = []

    # 4. Statistical Testing (Fisher's Exact Test)
    # We iterate through each cohort (e.g., Breast, Endometrium)
    for cohort in df_raw['Cohort'].unique():
        n_tumour_total = sample_counts.get((cohort, 'Tumour'), 0)
        n_normal_total = sample_counts.get((cohort, 'Normal'), 0)
        
        # We only perform the test if we have both groups to compare
        if n_tumour_total == 0 or n_normal_total == 0:
            continue

        print(f"Analysing {cohort} ({n_tumour_total} Tumour vs {n_normal_total} Normal)...")

        for _, row in df_snps.iterrows():
            var = row['Variant_ID']
            symbol = row['SYMBOL']
            
            # Count how many unique samples in each group carry this specific variant
            count_t = len(df_raw[(df_raw['Variant_ID'] == var) & 
                                 (df_raw['Cohort'] == cohort) & 
                                 (df_raw['Tissue'] == 'Tumour')]['Sample'].unique())
            
            count_n = len(df_raw[(df_raw['Variant_ID'] == var) & 
                                 (df_raw['Cohort'] == cohort) & 
                                 (df_raw['Tissue'] == 'Normal')]['Sample'].unique())

            # Build the 2x2 Contingency Table
            # [ [Variant_In_Tumour, No_Variant_In_Tumour], 
            #   [Variant_In_Normal, No_Variant_In_Normal] ]
            table = [[count_t, n_tumour_total - count_t],
                     [count_n, n_normal_total - count_n]]
            
            # Run the test: odds_ratio > 1 suggests enrichment in Tumour
            odds_ratio, p_value = stats.fisher_exact(table)
            
            results.append({
                'Cohort': cohort,
                'Variant_ID': var,
                'SYMBOL': symbol,
                'Tumour_Count': count_t,
                'Normal_Count': count_n,
                'Tumour_Freq_%': (count_t / n_tumour_total) * 100,
                'Normal_Freq_%': (count_n / n_normal_total) * 100,
                'Odds_Ratio': odds_ratio,
                'P_Value': p_value,
                '-log10_p': -np.log10(max(p_value, 1e-50)) if p_value > 0 else 50
            })

    # 5. Export Results
    stats_df = pd.DataFrame(results)
    stats_df = stats_df.sort_values(['P_Value', 'Odds_Ratio'], ascending=[True, False])
    out_xlsx = os.path.join(output_dir, "16_Statistical_Enrichment_Results.xlsx")
    stats_df.to_excel(out_xlsx, index=False)
    
    # 6. Volcano Plot: Visualise significance vs. magnitude of effect
    
    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=stats_df, x='Odds_Ratio', y='-log10_p', hue='SYMBOL', alpha=0.7, s=100)
    
    # Add a horizontal line at p=0.05
    plt.axhline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05 (Sig.)')
    plt.axvline(1, color='black', linestyle='-', alpha=0.3) # Neutral Odds Ratio line
    
    plt.xscale('log') # Odds ratios viewed on a log scale
    plt.title("Statistical Enrichment: Population SNPs & Novel Variants", fontsize=14, fontweight='bold')
    plt.xlabel("Odds Ratio (Log Scale - >1 is Tumour Enriched)")
    plt.ylabel("-log10(p-value)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "17_Volcano_Plot_Enrichment.png"), dpi=300)
    
    print(f"\n>>> Results saved to: {out_xlsx}")

if __name__ == "__main__":
    run_enrichment_analysis()
