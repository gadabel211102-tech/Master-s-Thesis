import pandas as pd
import os
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- PATH CONFIGURATION ---
base_path = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
# This is the filtered list from Script 11 (the "Shopping List")
input_file = os.path.join(base_path, "13_MASTER_SNP_Catalogue_NFE_Corrected.xlsx")
# This is the raw data used for counting patient occurrences
raw_data_for_denoms = os.path.join(base_path, "GSDMB_Annotated_Report_Fixed.xlsx")
output_dir = base_path

def run_enrichment_analysis():
    print(">>> Starting Statistical Enrichment (Population SNPs >1% and Novel)...")
    
    # 1. Load Filtered SNP Data from Script 11
    # We use this to know WHICH variants to look for
    if not os.path.exists(input_file):
        print(f"ERROR: Filtered file not found at {input_file}. Run Script 11 first.")
        return
    df_snps = pd.read_excel(input_file, sheet_name='All_Genes_SNPs')
    
    # 2. Load Raw Data
    df_raw = pd.read_excel(raw_data_for_denoms, sheet_name='Biological_Annotations')
    
    # --- FIX: Create the Variant_ID in the raw data so the script can count them ---
    df_raw['rsID'] = df_raw['Existing_variation'].astype(str).str.extract(r'(rs\d+)')
    df_raw['Variant_ID'] = df_raw['rsID'].fillna(df_raw['SYMBOL'] + ":" + df_raw['HGVSp'].astype(str))
    
    # Standardise cohort and tissue names
    df_raw['Cohort'] = df_raw['Cohort'].astype(str).str.strip()
    df_raw['Tissue'] = df_raw['Tissue'].astype(str).str.strip().replace('Tumor', 'Tumour')
    # -------------------------------------------------------------------------------

    # 3. Get Denominators (Total unique samples per group)
    sample_counts = df_raw.groupby(['Cohort', 'Tissue'])['Sample'].nunique().to_dict()

    results = []

    # 4. Statistical Testing (Fisher's Exact)
    # We iterate through the unique cohorts in the study
    for cohort in df_raw['Cohort'].unique():
        n_tumour_total = sample_counts.get((cohort, 'Tumour'), 0)
        n_normal_total = sample_counts.get((cohort, 'Normal'), 0)
        
        # Skip if we don't have both groups for comparison
        if n_tumour_total == 0 or n_normal_total == 0:
            print(f"Skipping {cohort}: Missing Tumour or Normal group.")
            continue

        print(f"Analysing {cohort} ({n_tumour_total} Tumour vs {n_normal_total} Normal)...")

        for _, row in df_snps.iterrows():
            var = row['Variant_ID']
            symbol = row['SYMBOL']
            
            # Count unique samples carrying this variant in each tissue
            count_t = len(df_raw[(df_raw['Variant_ID'] == var) & 
                                 (df_raw['Cohort'] == cohort) & 
                                 (df_raw['Tissue'] == 'Tumour')]['Sample'].unique())
            
            count_n = len(df_raw[(df_raw['Variant_ID'] == var) & 
                                 (df_raw['Cohort'] == cohort) & 
                                 (df_raw['Tissue'] == 'Normal')]['Sample'].unique())

            # Build the 2x2 Contingency Table for Fisher's Exact Test
            # [ [Has Variant in Tumour, No Variant in Tumour], 
            #   [Has Variant in Normal, No Variant in Normal] ]
            table = [[count_t, n_tumour_total - count_t],
                     [count_n, n_normal_total - count_n]]
            
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
                '-log10_p': -np.log10(p_value) if p_value > 0 else 0
            })

    stats_df = pd.DataFrame(results)
    
    # 5. Export Results
    stats_df = stats_df.sort_values(['P_Value', 'Odds_Ratio'], ascending=[True, False])
    out_xlsx = os.path.join(output_dir, "16_Statistical_Enrichment_Results.xlsx")
    stats_df.to_excel(out_xlsx, index=False)
    
    # 6. Volcano Plot (Filtered View)
    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=stats_df, x='Odds_Ratio', y='-log10_p', hue='SYMBOL', alpha=0.7, s=100)
    
    plt.axhline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05 (Sig.)')
    plt.axvline(1, color='black', linestyle='-', alpha=0.3) # Neutral line
    
    plt.xscale('log')
    plt.title("Statistical Enrichment: Population SNPs & Novel Variants", fontsize=14, fontweight='bold')
    plt.xlabel("Odds Ratio (Log Scale - >1 is Tumour Enriched)")
    plt.ylabel("-log10(p-value)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, which="both", ls="-", alpha=0.1)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "17_Volcano_Plot_Enrichment.png"), dpi=300)
    
    print(f"\n>>> Analysis Finished.")
    print(f"Results saved to: {out_xlsx}")

if __name__ == "__main__":
    run_enrichment_analysis()