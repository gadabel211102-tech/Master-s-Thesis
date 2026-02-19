#!/usr/bin/env python3
"""
Script 13: Permutation Testing for SNP Associations
===================================================
Addresses small sample size limitations using non-parametric permutation tests.
More robust than Fisher's exact test when sample sizes are imbalanced.

Date: 2026-02-12
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import permutation_test
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# --- CONFIGURATION ---
BASE_PATH = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
INPUT_FILE = os.path.join(BASE_PATH, "GSDMB_Annotated_Report_Fixed.xlsx")
OUTPUT_XLSX = os.path.join(BASE_PATH, "13_Permutation_Test_Results.xlsx")
OUTPUT_PLOT = os.path.join(BASE_PATH, "13_Permutation_Comparison.png")

# Permutation parameters
N_PERMUTATIONS = 10000  # Recommended: 10,000 for publication
RANDOM_SEED = 42


def find_col(df, target):
    """Case-insensitive column finder."""
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None


def permutation_test_snp(tumour_carriers, tumour_total, 
                         healthy_carriers, healthy_total,
                         n_permutations=10000, random_state=42):
    """
    Perform permutation test for SNP association.
    
    Null hypothesis: The SNP carrier status is independent of tissue type.
    
    Parameters:
    -----------
    tumour_carriers : int
        Number of tumour samples carrying the variant
    tumour_total : int
        Total number of tumour samples
    healthy_carriers : int
        Number of healthy samples carrying the variant
    healthy_total : int
        Total number of healthy samples
    n_permutations : int
        Number of permutations to perform
    random_state : int
        Random seed for reproducibility
        
    Returns:
    --------
    dict : Test statistics and p-value
    """
    np.random.seed(random_state)
    
    # Create binary arrays (1 = carrier, 0 = non-carrier)
    tumour_data = np.array([1] * tumour_carriers + [0] * (tumour_total - tumour_carriers))
    healthy_data = np.array([1] * healthy_carriers + [0] * (healthy_total - healthy_carriers))
    
    # Observed difference in proportions
    obs_tumour_prop = tumour_carriers / tumour_total if tumour_total > 0 else 0
    obs_healthy_prop = healthy_carriers / healthy_total if healthy_total > 0 else 0
    observed_diff = obs_tumour_prop - obs_healthy_prop
    
    # Combine all data
    all_data = np.concatenate([tumour_data, healthy_data])
    n_tumour = len(tumour_data)
    
    # Permutation test
    perm_diffs = []
    for _ in range(n_permutations):
        # Shuffle labels
        shuffled = np.random.permutation(all_data)
        
        # Split into two groups
        perm_tumour = shuffled[:n_tumour]
        perm_healthy = shuffled[n_tumour:]
        
        # Calculate difference
        perm_tumour_prop = np.mean(perm_tumour)
        perm_healthy_prop = np.mean(perm_healthy)
        perm_diffs.append(perm_tumour_prop - perm_healthy_prop)
    
    perm_diffs = np.array(perm_diffs)
    
    # Two-tailed p-value
    p_value = np.mean(np.abs(perm_diffs) >= np.abs(observed_diff))
    
    # Calculate effect size (Cohen's h for proportions)
    h = 2 * (np.arcsin(np.sqrt(obs_tumour_prop)) - 
             np.arcsin(np.sqrt(obs_healthy_prop)))
    
    return {
        'observed_diff': observed_diff,
        'p_value': p_value,
        'effect_size_h': h,
        'perm_diffs': perm_diffs,
        'tumour_prop': obs_tumour_prop,
        'healthy_prop': obs_healthy_prop
    }


def run_permutation_analysis():
    """Main pipeline for permutation testing."""
    print("="*70)
    print("SCRIPT 13: PERMUTATION TESTING FOR SNP ASSOCIATIONS")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  - Input: {INPUT_FILE}")
    print(f"  - Permutations: {N_PERMUTATIONS:,}")
    print(f"  - Random seed: {RANDOM_SEED}")
    print(f"  - SNP threshold: >1% gnomAD NFE AF")
    print("\n" + "="*70 + "\n")
    
    # Load data
    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: File not found: {INPUT_FILE}")
        return
    
    df = pd.read_excel(INPUT_FILE, sheet_name='Biological_Annotations')
    
    # Identify columns
    sym_c = find_col(df, 'SYMBOL')
    tis_c = find_col(df, 'TISSUE')
    sam_c = find_col(df, 'SAMPLE')
    var_c = find_col(df, 'Existing_variation')
    hgv_c = find_col(df, 'HGVSp')
    coh_c = find_col(df, 'COHORT')
    nfe_c = find_col(df, 'gnomADe_NFE_AF')
    imp_c = find_col(df, 'IMPACT')
    con_c = find_col(df, 'CONSEQUENCE')
    
    # Filter for SNPs (>1% AF)
    df = df[df[nfe_c] > 0.01].copy()
    
    # Standardize tissue labels
    df[tis_c] = df[tis_c].astype(str).str.strip().replace({
        'Tumor': 'Tumour', 'Normal': 'Healthy', 'Control': 'Healthy'
    })
    df[coh_c] = df[coh_c].astype(str).str.strip()
    
    # Create variant IDs
    df['rsID'] = df[var_c].astype(str).str.extract(r'(rs\d+)')
    df['Variant_ID'] = df['rsID'].fillna(
        df[sym_c].astype(str) + ":" + df[hgv_c].astype(str)
    )
    
    # Analysis groups
    cohort_list = ['Global', 'Breast', 'Endometrium']
    all_results = []
    
    for cohort_name in cohort_list:
        print(f"\n{'─'*70}")
        print(f"Processing: {cohort_name} Cohort")
        print(f"{'─'*70}")
        
        if cohort_name == 'Global':
            cohort_df = df.copy()
        else:
            cohort_df = df[df[coh_c].str.contains(cohort_name, case=False, na=False)].copy()
        
        if cohort_df.empty:
            print(f"  ⚠ No data found for {cohort_name}")
            continue
        
        # Count total samples
        n_tumour = cohort_df[cohort_df[tis_c] == 'Tumour'][sam_c].nunique()
        n_healthy = cohort_df[cohort_df[tis_c] == 'Healthy'][sam_c].nunique()
        
        if n_tumour == 0 or n_healthy == 0:
            print(f"  ⚠ Insufficient data: Tumour={n_tumour}, Healthy={n_healthy}")
            continue
        
        print(f"  Total samples: Tumour={n_tumour}, Healthy={n_healthy}")
        
        unique_snps = cohort_df['Variant_ID'].unique()
        print(f"  Unique SNPs to test: {len(unique_snps)}")
        
        # Progress bar
        for var in tqdm(unique_snps, desc=f"  Testing {cohort_name} SNPs"):
            var_data = cohort_df[cohort_df['Variant_ID'] == var]
            
            # Count carriers
            count_tumour = var_data[var_data[tis_c] == 'Tumour'][sam_c].nunique()
            count_healthy = var_data[var_data[tis_c] == 'Healthy'][sam_c].nunique()
            
            # Perform permutation test
            perm_result = permutation_test_snp(
                tumour_carriers=count_tumour,
                tumour_total=n_tumour,
                healthy_carriers=count_healthy,
                healthy_total=n_healthy,
                n_permutations=N_PERMUTATIONS,
                random_state=RANDOM_SEED
            )
            
            # Also run Fisher's exact for comparison
            from scipy.stats import fisher_exact
            a, b = count_tumour, n_tumour - count_tumour
            c, d = count_healthy, n_healthy - count_healthy
            fisher_or, fisher_p = fisher_exact([[a, b], [c, d]])
            
            # Store results
            all_results.append({
                'Cohort': cohort_name,
                'Variant_ID': var,
                'Symbol': var_data[sym_c].iloc[0],
                'Impact': var_data[imp_c].iloc[0],
                'Consequence': var_data[con_c].iloc[0],
                'gnomAD_NFE_AF': var_data[nfe_c].iloc[0],
                'Tumour_Carriers': count_tumour,
                'Tumour_Total': n_tumour,
                'Tumour_Freq_%': perm_result['tumour_prop'] * 100,
                'Healthy_Carriers': count_healthy,
                'Healthy_Total': n_healthy,
                'Healthy_Freq_%': perm_result['healthy_prop'] * 100,
                'Freq_Difference_%': perm_result['observed_diff'] * 100,
                'Permutation_P_Value': perm_result['p_value'],
                'Effect_Size_h': perm_result['effect_size_h'],
                'Fisher_P_Value': fisher_p,
                'Fisher_OR': fisher_or,
                'N_Permutations': N_PERMUTATIONS
            })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Multiple testing correction (FDR)
    from statsmodels.stats.multitest import multipletests
    
    corrected_results = []
    for cohort in results_df['Cohort'].unique():
        cohort_data = results_df[results_df['Cohort'] == cohort].copy()
        
        # FDR correction for permutation p-values
        _, cohort_data['Permutation_FDR'], _, _ = multipletests(
            cohort_data['Permutation_P_Value'], method='fdr_bh'
        )
        
        # FDR correction for Fisher p-values
        _, cohort_data['Fisher_FDR'], _, _ = multipletests(
            cohort_data['Fisher_P_Value'], method='fdr_bh'
        )
        
        corrected_results.append(cohort_data)
    
    final_results = pd.concat(corrected_results).sort_values(
        ['Cohort', 'Permutation_P_Value']
    )
    
    # Classify significance
    final_results['Permutation_Significant'] = final_results['Permutation_P_Value'] < 0.05
    final_results['Fisher_Significant'] = final_results['Fisher_P_Value'] < 0.05
    final_results['Agreement'] = (
        final_results['Permutation_Significant'] == final_results['Fisher_Significant']
    )
    
    # Save results
    print(f"\n{'='*70}")
    print("SAVING RESULTS")
    print(f"{'='*70}\n")
    
    with pd.ExcelWriter(OUTPUT_XLSX, engine='openpyxl') as writer:
        # Summary sheet
        summary = final_results.groupby('Cohort').agg({
            'Variant_ID': 'count',
            'Permutation_Significant': 'sum',
            'Fisher_Significant': 'sum',
            'Agreement': 'sum'
        }).reset_index()
        summary.columns = ['Cohort', 'Total_SNPs', 'Perm_Significant', 
                          'Fisher_Significant', 'Test_Agreement']
        summary.to_excel(writer, sheet_name='Summary', index=False)
        
        # Individual cohort sheets
        for cohort in cohort_list:
            cohort_data = final_results[final_results['Cohort'] == cohort]
            if not cohort_data.empty:
                cohort_data.to_excel(writer, sheet_name=cohort, index=False)
    
    print(f"✓ Results saved to: {OUTPUT_XLSX}\n")
    
    # Generate comparison visualization
    generate_comparison_plots(final_results)
    
    # Print summary statistics
    print(f"{'='*70}")
    print("SUMMARY STATISTICS")
    print(f"{'='*70}\n")
    
    for cohort in cohort_list:
        cohort_data = final_results[final_results['Cohort'] == cohort]
        if cohort_data.empty:
            continue
        
        print(f"{cohort} Cohort:")
        print(f"  Total SNPs tested: {len(cohort_data)}")
        print(f"  Permutation sig. (FDR<0.05): {cohort_data['Permutation_Significant'].sum()}")
        print(f"  Fisher sig. (FDR<0.05): {cohort_data['Fisher_Significant'].sum()}")
        print(f"  Test agreement: {cohort_data['Agreement'].sum()}/{len(cohort_data)} "
              f"({cohort_data['Agreement'].mean()*100:.1f}%)")
        
        # Top 5 most significant by permutation
        top5 = cohort_data.nsmallest(5, 'Permutation_P_Value')
        if not top5.empty:
            print(f"\n  Top 5 SNPs by permutation test:")
            for idx, row in top5.iterrows():
                print(f"    {row['Variant_ID']} (p={row['Permutation_P_Value']:.4e}, "
                      f"FDR={row['Permutation_FDR']:.4e})")
        print()
    
    print(f"{'='*70}")
    print("✓ Permutation analysis complete!")
    print(f"{'='*70}\n")


def generate_comparison_plots(results_df):
    """Generate visualization comparing permutation vs Fisher's exact test."""
    print("Generating comparison plots...")
    
    # Setup
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    
    cohorts = results_df['Cohort'].unique()
    colors = {'Global': '#3498db', 'Breast': '#e74c3c', 'Endometrium': '#2ecc71'}
    
    # Plot 1: P-value comparison scatter
    ax1 = axes[0, 0]
    for cohort in cohorts:
        data = results_df[results_df['Cohort'] == cohort]
        ax1.scatter(
            -np.log10(data['Fisher_P_Value']),
            -np.log10(data['Permutation_P_Value']),
            label=cohort, alpha=0.6, s=80,
            color=colors.get(cohort, '#95a5a6')
        )
    
    # Add diagonal line and significance thresholds
    max_val = max(ax1.get_xlim()[1], ax1.get_ylim()[1])
    ax1.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='x=y')
    sig_line = -np.log10(0.05)
    ax1.axhline(sig_line, color='red', linestyle=':', alpha=0.5, label='p=0.05')
    ax1.axvline(sig_line, color='red', linestyle=':', alpha=0.5)
    
    ax1.set_xlabel("Fisher's Exact Test (-log10 p)", fontweight='bold')
    ax1.set_ylabel("Permutation Test (-log10 p)", fontweight='bold')
    ax1.set_title("P-Value Comparison: Permutation vs Fisher", fontweight='bold', fontsize=12)
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Effect size vs Permutation P-value
    ax2 = axes[0, 1]
    for cohort in cohorts:
        data = results_df[results_df['Cohort'] == cohort]
        ax2.scatter(
            data['Effect_Size_h'],
            -np.log10(data['Permutation_P_Value']),
            label=cohort, alpha=0.6, s=80,
            color=colors.get(cohort, '#95a5a6')
        )
    
    ax2.axhline(sig_line, color='red', linestyle=':', alpha=0.5, label='p=0.05')
    ax2.set_xlabel("Effect Size (Cohen's h)", fontweight='bold')
    ax2.set_ylabel("Permutation Test (-log10 p)", fontweight='bold')
    ax2.set_title("Effect Size vs Statistical Significance", fontweight='bold', fontsize=12)
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Frequency difference distribution
    ax3 = axes[1, 0]
    for cohort in cohorts:
        data = results_df[results_df['Cohort'] == cohort]
        ax3.hist(
            data['Freq_Difference_%'],
            bins=30, alpha=0.5, label=cohort,
            color=colors.get(cohort, '#95a5a6'),
            edgecolor='black'
        )
    
    ax3.axvline(0, color='black', linestyle='--', alpha=0.5, label='No difference')
    ax3.set_xlabel("Frequency Difference (Tumour - Healthy %)", fontweight='bold')
    ax3.set_ylabel("Number of SNPs", fontweight='bold')
    ax3.set_title("Distribution of Frequency Differences", fontweight='bold', fontsize=12)
    ax3.legend(loc='best')
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Agreement between tests
    ax4 = axes[1, 1]
    agreement_data = []
    for cohort in cohorts:
        data = results_df[results_df['Cohort'] == cohort]
        both_sig = ((data['Permutation_Significant']) & (data['Fisher_Significant'])).sum()
        perm_only = ((data['Permutation_Significant']) & (~data['Fisher_Significant'])).sum()
        fisher_only = ((~data['Permutation_Significant']) & (data['Fisher_Significant'])).sum()
        neither = ((~data['Permutation_Significant']) & (~data['Fisher_Significant'])).sum()
        
        agreement_data.append({
            'Cohort': cohort,
            'Both Significant': both_sig,
            'Perm Only': perm_only,
            'Fisher Only': fisher_only,
            'Neither': neither
        })
    
    agreement_df = pd.DataFrame(agreement_data).set_index('Cohort')
    agreement_df.plot(kind='bar', stacked=True, ax=ax4, 
                     color=['#2ecc71', '#f39c12', '#e74c3c', '#95a5a6'])
    ax4.set_xlabel("Cohort", fontweight='bold')
    ax4.set_ylabel("Number of SNPs", fontweight='bold')
    ax4.set_title("Test Agreement Analysis", fontweight='bold', fontsize=12)
    ax4.legend(title="Significance", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.setp(ax4.xaxis.get_majorticklabels(), rotation=0)
    ax4.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle(
        f'Permutation Testing vs Fisher\'s Exact Test Comparison\n'
        f'({N_PERMUTATIONS:,} permutations per SNP)',
        fontsize=16, fontweight='bold', y=0.995
    )
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Comparison plots saved to: {OUTPUT_PLOT}\n")


if __name__ == "__main__":
    run_permutation_analysis()
