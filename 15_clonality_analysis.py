#!/usr/bin/env python3
"""
Script 16: Clonality Analysis - Tumor Heterogeneity Assessment
===============================================================
Analyzes allele frequency distributions to identify:
- Clonal mutations (present in all tumor cells)
- Subclonal mutations (present in subset of cells)
- Tumor purity and heterogeneity metrics

Author: Enhanced Genomics Analysis Pipeline
Date: 2026-02-12
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings('ignore')

# --- CONFIGURATION ---
BASE_PATH = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
INPUT_FILE = os.path.join(BASE_PATH, "GSDMB_Annotated_Report_Fixed.xlsx")
OUTPUT_XLSX = os.path.join(BASE_PATH, "16_Clonality_Analysis_Results.xlsx")
OUTPUT_PLOT = os.path.join(BASE_PATH, "16_Clonality_Plots.png")

# Clonality thresholds
CLONAL_THRESHOLD = 0.35  # VAF > 0.35 typically indicates clonal
SUBCLONAL_LOW = 0.10     # VAF < 0.10 indicates rare subclone
SUBCLONAL_MID = 0.25     # 0.10 < VAF < 0.25 is subclonal

# Minimum read depth for reliable VAF estimation
MIN_DEPTH_FOR_VAF = 30


def find_col(df, target):
    """Case-insensitive column finder."""
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None


def extract_vaf_from_annotations(df):
    """
    Extract real VAF and depth from the AF and DP columns produced by script 07.
    Both come directly from VCF FORMAT fields (AO/RO-derived AF and DP).
    """
    # Locate AF and DP columns case-insensitively
    af_col = next((c for c in df.columns if c.upper() == 'AF'), None)
    dp_col = next((c for c in df.columns if c.upper() == 'DP'), None)

    if af_col is None:
        raise ValueError(
            "No AF column found in the annotated Excel. "
            "Re-run script 07 to ensure AF is extracted from the VCFs."
        )
    if dp_col is None:
        raise ValueError(
            "No DP column found in the annotated Excel. "
            "Re-run script 07 to ensure DP is extracted from the VCFs."
        )

    # Coerce to numeric — VEP occasionally writes '.' for missing values
    df['VAF'] = pd.to_numeric(df[af_col], errors='coerce')
    df['Read_Depth'] = pd.to_numeric(df[dp_col], errors='coerce')

    # Drop rows with no real VAF (e.g. reference-only sites)
    before = len(df)
    df = df.dropna(subset=['VAF', 'Read_Depth'])
    dropped = before - len(df)
    if dropped > 0:
        print(f"  ℹ Dropped {dropped} rows with missing AF or DP values")

    # Clamp VAF to [0.001, 1.0] to protect downstream log transforms
    df['VAF'] = df['VAF'].clip(0.001, 1.0)

    print(f"✓ Real VAF extracted from '{af_col}' column")
    print(f"  VAF range: {df['VAF'].min():.3f} – {df['VAF'].max():.3f}")
    print(f"  Median depth: {df['Read_Depth'].median():.0f}x\n")

    return df


def classify_clonality(vaf, tissue='Tumour'):
    """
    Classify variant as clonal, subclonal, or germline based on VAF.
    
    Parameters:
        vaf: Variant allele frequency (0-1)
        tissue: 'Tumour' or 'Healthy'
    
    Returns:
        String classification
    """
    if tissue != 'Tumour':
        # Germline variants
        if vaf >= 0.9:
            return 'Germline_Homozygous'
        elif 0.4 <= vaf <= 0.6:
            return 'Germline_Heterozygous'
        else:
            return 'Germline_Unknown'
    
    # Somatic variants in tumor
    if vaf >= CLONAL_THRESHOLD:
        return 'Clonal'
    elif vaf >= SUBCLONAL_MID:
        return 'Subclonal_Major'
    elif vaf >= SUBCLONAL_LOW:
        return 'Subclonal_Minor'
    else:
        return 'Rare_Subclone'


def calculate_heterogeneity_metrics(sample_vafs):
    """
    Calculate tumor heterogeneity metrics for a single sample.
    
    Metrics:
    - Shannon entropy: Higher = more heterogeneous
    - VAF variance: Higher = more diverse subclones
    - Clonal fraction: % of variants that are clonal
    """
    if len(sample_vafs) == 0:
        return {}
    
    # Shannon entropy of VAF distribution
    hist, _ = np.histogram(sample_vafs, bins=10, range=(0, 1))
    hist = hist / hist.sum()  # Normalize
    hist = hist[hist > 0]  # Remove zeros
    entropy = -np.sum(hist * np.log2(hist))
    
    # VAF statistics
    vaf_mean = np.mean(sample_vafs)
    vaf_std = np.std(sample_vafs)
    vaf_variance = np.var(sample_vafs)
    
    # Clonal fraction
    clonal_count = np.sum(sample_vafs >= CLONAL_THRESHOLD)
    clonal_fraction = clonal_count / len(sample_vafs)
    
    # Subclonal diversity
    subclonal_vafs = sample_vafs[sample_vafs < CLONAL_THRESHOLD]
    subclonal_diversity = len(np.unique(np.round(subclonal_vafs, 2)))
    
    return {
        'Shannon_Entropy': entropy,
        'Mean_VAF': vaf_mean,
        'VAF_StdDev': vaf_std,
        'VAF_Variance': vaf_variance,
        'Clonal_Fraction': clonal_fraction,
        'N_Clonal_Variants': clonal_count,
        'N_Subclonal_Variants': len(sample_vafs) - clonal_count,
        'Subclonal_Diversity': subclonal_diversity,
        'Total_Variants': len(sample_vafs)
    }


def cluster_vaf_distribution(vafs, n_clusters=3):
    """
    Use k-means to identify distinct subclonal populations.
    
    Returns:
        Cluster assignments and centers
    """
    if len(vafs) < n_clusters:
        return None, None
    
    vafs_2d = vafs.reshape(-1, 1)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(vafs_2d)
    centers = kmeans.cluster_centers_.flatten()
    
    return labels, sorted(centers)


def run_clonality_analysis():
    """Main pipeline for clonality analysis."""
    print("="*70)
    print("SCRIPT 16: TUMOR CLONALITY & HETEROGENEITY ANALYSIS")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  - Input: {INPUT_FILE}")
    print(f"  - Clonal threshold: VAF ≥ {CLONAL_THRESHOLD}")
    print(f"  - Subclonal ranges: {SUBCLONAL_LOW}-{SUBCLONAL_MID} (minor), "
          f"{SUBCLONAL_MID}-{CLONAL_THRESHOLD} (major)")
    print(f"  - Minimum depth: {MIN_DEPTH_FOR_VAF}x")
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
    coh_c = find_col(df, 'COHORT')
    imp_c = find_col(df, 'IMPACT')
    pos_c = find_col(df, 'Pos')
    
    # Standardize tissue labels
    df[tis_c] = df[tis_c].astype(str).str.strip().replace({
        'Tumor': 'Tumour', 'Normal': 'Healthy', 'Control': 'Healthy'
    })
    df[coh_c] = df[coh_c].astype(str).str.strip()
    
    # Extract/simulate VAF
    print("STEP 1: Extracting Variant Allele Frequencies")
    print("─"*70)
    df = extract_vaf_from_annotations(df)
    print(f"✓ VAF extracted for {len(df)} variants\n")
    
    # Filter high-quality variants
    df_filtered = df[df['Read_Depth'] >= MIN_DEPTH_FOR_VAF].copy()
    print(f"✓ Filtered to {len(df_filtered)} variants with ≥{MIN_DEPTH_FOR_VAF}x depth\n")
    
    # Classify clonality
    print("STEP 2: Classifying Variant Clonality")
    print("─"*70)
    df_filtered['Clonality'] = df_filtered.apply(
        lambda x: classify_clonality(x['VAF'], x[tis_c]), axis=1
    )
    
    clonality_counts = df_filtered[df_filtered[tis_c] == 'Tumour']['Clonality'].value_counts()
    print("Tumor variant classification:")
    for cat, count in clonality_counts.items():
        print(f"  {cat}: {count}")
    print()
    
    # Calculate per-sample heterogeneity
    print("STEP 3: Calculating Tumor Heterogeneity Metrics")
    print("─"*70)
    
    heterogeneity_results = []
    tumor_samples = df_filtered[df_filtered[tis_c] == 'Tumour'][sam_c].unique()
    
    for sample in tumor_samples:
        sample_data = df_filtered[
            (df_filtered[sam_c] == sample) & 
            (df_filtered[tis_c] == 'Tumour')
        ]
        
        if len(sample_data) == 0:
            continue
        
        cohort = sample_data[coh_c].iloc[0]
        vafs = sample_data['VAF'].values
        
        # Calculate metrics
        metrics = calculate_heterogeneity_metrics(vafs)
        
        # Add sample info
        metrics['Sample'] = sample
        metrics['Cohort'] = cohort
        
        # Cluster analysis
        if len(vafs) >= 10:  # Need enough variants for clustering
            labels, centers = cluster_vaf_distribution(vafs, n_clusters=3)
            if centers is not None:
                metrics['N_Subclones_Detected'] = len(centers)
                metrics['Subclone_VAFs'] = ', '.join([f"{c:.2f}" for c in centers])
        
        heterogeneity_results.append(metrics)
    
    heterogeneity_df = pd.DataFrame(heterogeneity_results)
    print(f"✓ Analyzed {len(heterogeneity_df)} tumor samples\n")
    
    # Summary statistics
    print("Heterogeneity Summary:")
    print(f"  Mean Shannon Entropy: {heterogeneity_df['Shannon_Entropy'].mean():.3f}")
    print(f"  Mean Clonal Fraction: {heterogeneity_df['Clonal_Fraction'].mean():.2%}")
    print(f"  Mean VAF: {heterogeneity_df['Mean_VAF'].mean():.3f}\n")
    
    # Compare cohorts
    print("STEP 4: Comparing Heterogeneity Between Cancer Types")
    print("─"*70)
    
    cohort_comparison = heterogeneity_df.groupby('Cohort').agg({
        'Shannon_Entropy': ['mean', 'std'],
        'Clonal_Fraction': ['mean', 'std'],
        'VAF_Variance': ['mean', 'std'],
        'Total_Variants': 'mean'
    }).round(3)
    
    print(cohort_comparison.to_string())
    print()
    
    # Statistical tests
    if len(heterogeneity_df['Cohort'].unique()) >= 2:
        cohorts = heterogeneity_df['Cohort'].unique()
        if len(cohorts) == 2:
            c1_data = heterogeneity_df[heterogeneity_df['Cohort'] == cohorts[0]]['Shannon_Entropy']
            c2_data = heterogeneity_df[heterogeneity_df['Cohort'] == cohorts[1]]['Shannon_Entropy']
            
            t_stat, p_val = stats.ttest_ind(c1_data, c2_data)
            print(f"T-test for Shannon Entropy difference:")
            print(f"  {cohorts[0]} vs {cohorts[1]}: t={t_stat:.3f}, p={p_val:.4f}\n")
    
    # Save results
    print("STEP 5: Saving Results")
    print("─"*70)
    
    with pd.ExcelWriter(OUTPUT_XLSX, engine='openpyxl') as writer:
        # Sample heterogeneity metrics
        heterogeneity_df.to_excel(writer, sheet_name='Sample_Heterogeneity', index=False)
        
        # All variants with clonality
        df_filtered.to_excel(writer, sheet_name='All_Variants_With_Clonality', index=False)
        
        # Clonal vs subclonal by gene
        clonal_by_gene = df_filtered[df_filtered[tis_c] == 'Tumour'].groupby(
            [sym_c, 'Clonality']
        ).size().unstack(fill_value=0)
        clonal_by_gene.to_excel(writer, sheet_name='Clonality_By_Gene')
        
        # Cohort summary
        cohort_comparison.to_excel(writer, sheet_name='Cohort_Comparison')
    
    print(f"✓ Results saved to: {OUTPUT_XLSX}\n")
    
    # Generate visualizations
    generate_clonality_plots(df_filtered, heterogeneity_df, sym_c, tis_c, sam_c, coh_c, imp_c)
    
    print("="*70)
    print("✓ Clonality analysis complete!")
    print("="*70 + "\n")


def generate_clonality_plots(df, heterogeneity_df, sym_c, tis_c, sam_c, coh_c, imp_c='IMPACT'):
    """Generate comprehensive clonality visualizations."""
    print("Generating clonality plots...")
    
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3)
    
    tumor_data = df[df[tis_c] == 'Tumour']
    healthy_data = df[df[tis_c] == 'Healthy']
    
    # Plot 1: VAF Distribution - Tumor vs Healthy
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.hist(tumor_data['VAF'], bins=50, alpha=0.6, label='Tumor', 
            color='#e74c3c', edgecolor='black')
    ax1.hist(healthy_data['VAF'], bins=50, alpha=0.6, label='Healthy', 
            color='#3498db', edgecolor='black')
    
    # Add threshold lines
    ax1.axvline(CLONAL_THRESHOLD, color='red', linestyle='--', 
               linewidth=2, label=f'Clonal threshold ({CLONAL_THRESHOLD})')
    ax1.axvline(0.5, color='green', linestyle=':', 
               linewidth=2, label='Expected germline (0.5)')
    
    ax1.set_xlabel("Variant Allele Frequency (VAF)", fontweight='bold')
    ax1.set_ylabel("Number of Variants", fontweight='bold')
    ax1.set_title("VAF Distribution: Tumor vs Healthy Tissue", 
                 fontweight='bold', fontsize=12)
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Clonality Classification
    ax2 = fig.add_subplot(gs[0, 2])
    clonality_counts = tumor_data['Clonality'].value_counts()
    colors_clon = {
        'Clonal': '#e74c3c',
        'Subclonal_Major': '#f39c12',
        'Subclonal_Minor': '#3498db',
        'Rare_Subclone': '#95a5a6'
    }
    colors_list = [colors_clon.get(cat, '#95a5a6') for cat in clonality_counts.index]
    
    wedges, texts, autotexts = ax2.pie(
        clonality_counts.values,
        labels=clonality_counts.index,
        autopct='%1.1f%%',
        colors=colors_list,
        startangle=90
    )
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    ax2.set_title("Tumor Clonality\nClassification", fontweight='bold', fontsize=11)
    
    # Plot 3: VAF by Impact Level
    ax3 = fig.add_subplot(gs[1, 0])
    impact_order = ['HIGH', 'MODERATE', 'MODIFIER', 'LOW']
    impact_colors = {'HIGH': '#d9534f', 'MODERATE': '#5bc0de', 
                    'MODIFIER': '#f0ad4e', 'LOW': '#5cb85c'}
    
    vaf_by_impact = []
    impact_labels = []
    for impact in impact_order:
        data = tumor_data[tumor_data[imp_c] == impact]['VAF']
        if len(data) > 0:
            vaf_by_impact.append(data)
            impact_labels.append(impact)
    
    bp = ax3.boxplot(vaf_by_impact, labels=impact_labels, patch_artist=True)
    for patch, label in zip(bp['boxes'], impact_labels):
        patch.set_facecolor(impact_colors.get(label, '#95a5a6'))
        patch.set_alpha(0.7)
    
    ax3.axhline(CLONAL_THRESHOLD, color='red', linestyle='--', alpha=0.5)
    ax3.set_ylabel("VAF", fontweight='bold')
    ax3.set_title("VAF Distribution by Impact Level", fontweight='bold', fontsize=11)
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Sample Heterogeneity (Shannon Entropy)
    ax4 = fig.add_subplot(gs[1, 1])
    cohorts = heterogeneity_df['Cohort'].unique()
    hetero_data = [heterogeneity_df[heterogeneity_df['Cohort'] == c]['Shannon_Entropy']
                  for c in cohorts]
    
    bp2 = ax4.boxplot(hetero_data, labels=cohorts, patch_artist=True)
    colors_cohort = {'Breast': '#e74c3c', 'Endometrium': '#2ecc71', 'Global': '#3498db'}
    for patch, cohort in zip(bp2['boxes'], cohorts):
        patch.set_facecolor(colors_cohort.get(cohort, '#95a5a6'))
        patch.set_alpha(0.7)
    
    ax4.set_ylabel("Shannon Entropy", fontweight='bold')
    ax4.set_title("Tumor Heterogeneity by Cohort", fontweight='bold', fontsize=11)
    ax4.grid(True, alpha=0.3, axis='y')
    
    # Plot 5: Clonal Fraction vs Total Variants
    ax5 = fig.add_subplot(gs[1, 2])
    for cohort in cohorts:
        cohort_data = heterogeneity_df[heterogeneity_df['Cohort'] == cohort]
        ax5.scatter(
            cohort_data['Total_Variants'],
            cohort_data['Clonal_Fraction'],
            label=cohort, alpha=0.6, s=100,
            color=colors_cohort.get(cohort, '#95a5a6'),
            edgecolor='black'
        )
    
    ax5.set_xlabel("Total Variants", fontweight='bold')
    ax5.set_ylabel("Clonal Fraction", fontweight='bold')
    ax5.set_title("Clonal Burden Analysis", fontweight='bold', fontsize=11)
    ax5.legend(loc='best')
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: VAF Heatmap by Sample
    ax6 = fig.add_subplot(gs[2, :])
    
    # Create matrix of top genes × samples
    top_genes = tumor_data.groupby(sym_c).size().nlargest(15).index
    samples_subset = tumor_data[sam_c].unique()[:20]  # Limit to 20 samples
    
    vaf_matrix = []
    for gene in top_genes:
        gene_vafs = []
        for sample in samples_subset:
            sample_gene_data = tumor_data[
                (tumor_data[sym_c] == gene) & 
                (tumor_data[sam_c] == sample)
            ]
            if len(sample_gene_data) > 0:
                gene_vafs.append(sample_gene_data['VAF'].mean())
            else:
                gene_vafs.append(0)
        vaf_matrix.append(gene_vafs)
    
    vaf_matrix = np.array(vaf_matrix)
    
    im = ax6.imshow(vaf_matrix, aspect='auto', cmap='YlOrRd', vmin=0, vmax=1)
    ax6.set_yticks(range(len(top_genes)))
    ax6.set_yticklabels(top_genes)
    ax6.set_xticks(range(min(10, len(samples_subset))))
    ax6.set_xticklabels([f"S{i+1}" for i in range(min(10, len(samples_subset)))], 
                        rotation=45)
    ax6.set_xlabel("Samples", fontweight='bold')
    ax6.set_ylabel("Genes", fontweight='bold')
    ax6.set_title("VAF Heatmap: Top Genes × Samples", fontweight='bold', fontsize=12)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax6, fraction=0.046, pad=0.04)
    cbar.set_label('VAF', rotation=270, labelpad=20, fontweight='bold')
    
    plt.suptitle('Tumor Clonality & Heterogeneity Analysis', 
                fontsize=16, fontweight='bold', y=0.995)
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Clonality plots saved to: {OUTPUT_PLOT}\n")


if __name__ == "__main__":
    run_clonality_analysis()