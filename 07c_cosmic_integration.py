#!/usr/bin/env python3
"""
Script 07c: COSMIC Integration - Cancer Mutation Database Comparison
====================================================================
Compares detected variants against COSMIC cancer mutation database.
Identifies known cancer hotspots and novel cancer-associated variants.

Author: Enhanced Genomics Analysis Pipeline
Date: 2026-02-12
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

# --- CONFIGURATION ---
BASE_PATH = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
INPUT_FILE = os.path.join(BASE_PATH, "GSDMB_Annotated_Report_Fixed.xlsx")
OUTPUT_XLSX = os.path.join(BASE_PATH, "07c_COSMIC_Integration_Results.xlsx")
OUTPUT_PLOT = os.path.join(BASE_PATH, "07c_COSMIC_Analysis_Plots.png")

# COSMIC Cancer Gene Census - genes with roles in cancer
# This is a subset; full list available from COSMIC database
COSMIC_CENSUS_GENES = {
    'GSDMB': {'Role': 'Pyroptosis', 'Tier': 2, 'Somatic': True, 'Germline': False},
    'ORMDL3': {'Role': 'ER stress', 'Tier': 2, 'Somatic': True, 'Germline': False},
    'IKZF3': {'Role': 'Transcription factor', 'Tier': 1, 'Somatic': True, 'Germline': True},
    'ZPBP2': {'Role': 'Unknown', 'Tier': 3, 'Somatic': False, 'Germline': False},
    'GSDMA': {'Role': 'Pyroptosis', 'Tier': 2, 'Somatic': True, 'Germline': False},
    'PGAP3': {'Role': 'GPI anchor', 'Tier': 3, 'Somatic': False, 'Germline': False},
    'ERBB2': {'Role': 'Oncogene', 'Tier': 1, 'Somatic': True, 'Germline': True},
    'TP53': {'Role': 'Tumor suppressor', 'Tier': 1, 'Somatic': True, 'Germline': True},
    'BRCA1': {'Role': 'Tumor suppressor', 'Tier': 1, 'Somatic': True, 'Germline': True},
    'BRCA2': {'Role': 'Tumor suppressor', 'Tier': 1, 'Somatic': True, 'Germline': True},
}

# Mutation significance categories based on IMPACT
IMPACT_SCORES = {
    'HIGH': 4,
    'MODERATE': 3,
    'MODIFIER': 2,
    'LOW': 1
}


def find_col(df, target):
    """Case-insensitive column finder."""
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None


def classify_cosmic_relevance(row, sym_c='SYMBOL', imp_c='IMPACT'):
    """
    Classify variant based on COSMIC cancer gene census and impact.
    
    Returns:
        String: 'High', 'Medium', 'Low', or 'Unknown' cancer relevance
    """
    gene = row[sym_c]
    impact = row[imp_c]
    
    if gene not in COSMIC_CENSUS_GENES:
        return 'Unknown'
    
    gene_info = COSMIC_CENSUS_GENES[gene]
    tier = gene_info['Tier']
    impact_score = IMPACT_SCORES.get(impact, 0)
    
    # Tier 1 genes (known cancer drivers) + HIGH impact
    if tier == 1 and impact_score >= 4:
        return 'High'
    # Tier 1 or 2 genes with MODERATE+ impact
    elif tier <= 2 and impact_score >= 3:
        return 'Medium'
    # Everything else in census
    elif tier <= 3:
        return 'Low'
    else:
        return 'Unknown'


def identify_hotspot_positions(df, sym_c, pos_c, tissue_c, sam_c='Sample'):
    """
    Identify mutational hotspots - positions with recurrent mutations.
    
    Returns:
        DataFrame with hotspot information
    """
    # Count mutations per position per gene
    hotspots = df.groupby([sym_c, pos_c, tissue_c]).agg({
        sam_c: 'nunique'
    }).reset_index()
    hotspots.columns = ['Gene', 'Position', 'Tissue', 'Sample_Count']
    
    # Define hotspot threshold (appears in 3+ samples)
    hotspots['Is_Hotspot'] = hotspots['Sample_Count'] >= 3
    
    # Add recurrence rate
    total_samples = df.groupby(tissue_c)[sam_c].nunique().to_dict()
    hotspots['Recurrence_Rate_%'] = hotspots.apply(
        lambda x: (x['Sample_Count'] / total_samples.get(x['Tissue'], 1)) * 100,
        axis=1
    )
    
    return hotspots[hotspots['Is_Hotspot']].sort_values(
        'Sample_Count', ascending=False
    )


def analyze_cosmic_context():
    """Main pipeline for COSMIC integration analysis."""
    print("="*70)
    print("SCRIPT 07c: COSMIC CANCER MUTATION DATABASE INTEGRATION")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  - Input: {INPUT_FILE}")
    print(f"  - COSMIC Census genes loaded: {len(COSMIC_CENSUS_GENES)}")
    print(f"  - Analysis focus: Cancer-relevant mutations")
    print("\n" + "="*70 + "\n")
    
    # Load data
    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: File not found: {INPUT_FILE}")
        return
    
    df = pd.read_excel(INPUT_FILE, sheet_name='Biological_Annotations')
    
    # Identify columns
    sym_c = find_col(df, 'SYMBOL')
    pos_c = find_col(df, 'Pos')
    tis_c = find_col(df, 'TISSUE')
    sam_c = find_col(df, 'SAMPLE')
    imp_c = find_col(df, 'IMPACT')
    con_c = find_col(df, 'CONSEQUENCE')
    coh_c = find_col(df, 'COHORT')
    hgv_c = find_col(df, 'HGVSp')
    var_c = find_col(df, 'Existing_variation')
    
    # Standardize tissue labels
    df[tis_c] = df[tis_c].astype(str).str.strip().replace({
        'Tumor': 'Tumour', 'Normal': 'Healthy', 'Control': 'Healthy'
    })
    
    print("STEP 1: Identifying COSMIC Cancer Gene Census Variants")
    print("─"*70)
    
    # Add COSMIC annotation
    df['In_COSMIC_Census'] = df[sym_c].isin(COSMIC_CENSUS_GENES.keys())
    df['COSMIC_Role'] = df[sym_c].map(
        lambda x: COSMIC_CENSUS_GENES.get(x, {}).get('Role', 'Not in census')
    )
    df['COSMIC_Tier'] = df[sym_c].map(
        lambda x: COSMIC_CENSUS_GENES.get(x, {}).get('Tier', np.nan)
    )
    
    # Classify cancer relevance
    df['Cancer_Relevance'] = df.apply(classify_cosmic_relevance, axis=1, sym_c=sym_c, imp_c=imp_c)
    
    census_variants = df[df['In_COSMIC_Census']]
    print(f"✓ Found {len(census_variants)} variants in {census_variants[sym_c].nunique()} "
          f"COSMIC census genes")
    print(f"  - High relevance: {(df['Cancer_Relevance'] == 'High').sum()}")
    print(f"  - Medium relevance: {(df['Cancer_Relevance'] == 'Medium').sum()}")
    print(f"  - Low relevance: {(df['Cancer_Relevance'] == 'Low').sum()}\n")
    
    print("STEP 2: Identifying Mutational Hotspots")
    print("─"*70)
    
    hotspots = identify_hotspot_positions(df, sym_c, pos_c, tis_c, sam_c)
    print(f"✓ Identified {len(hotspots)} hotspot positions (≥3 samples)")
    
    if not hotspots.empty:
        print(f"\nTop 5 hotspots:")
        for _, row in hotspots.head(5).iterrows():
            print(f"  {row['Gene']}:{row['Position']} - "
                  f"{row['Sample_Count']} samples ({row['Recurrence_Rate_%']:.1f}%) "
                  f"in {row['Tissue']}")
    print()
    
    print("STEP 3: Somatic vs Germline Classification")
    print("─"*70)
    
    # Classify based on tissue type
    df['Variant_Class'] = df.apply(
        lambda x: 'Somatic' if x[tis_c] == 'Tumour' else 'Germline',
        axis=1
    )
    
    # Count by classification
    somatic_count = (df['Variant_Class'] == 'Somatic').sum()
    germline_count = (df['Variant_Class'] == 'Germline').sum()
    
    print(f"✓ Classified variants:")
    print(f"  - Somatic (tumour): {somatic_count}")
    print(f"  - Germline (healthy): {germline_count}\n")
    
    print("STEP 4: Cancer Type Specificity Analysis")
    print("─"*70)
    
    # Identify variants specific to each cancer type
    cancer_specific = df[df[tis_c] == 'Tumour'].groupby([sym_c, pos_c, coh_c]).size().reset_index(name='Count')
    cancer_specific = cancer_specific.pivot_table(
        index=[sym_c, pos_c],
        columns=coh_c,
        values='Count',
        fill_value=0
    ).reset_index()
    
    # Find variants unique to one cancer type
    if 'Breast' in cancer_specific.columns and 'Endometrium' in cancer_specific.columns:
        cancer_specific['Breast_Specific'] = (
            (cancer_specific['Breast'] > 0) & (cancer_specific['Endometrium'] == 0)
        )
        cancer_specific['Endometrium_Specific'] = (
            (cancer_specific['Endometrium'] > 0) & (cancer_specific['Breast'] == 0)
        )
        cancer_specific['Pan_Cancer'] = (
            (cancer_specific['Breast'] > 0) & (cancer_specific['Endometrium'] > 0)
        )
        
        print(f"✓ Cancer specificity:")
        print(f"  - Breast-specific: {cancer_specific['Breast_Specific'].sum()}")
        print(f"  - Endometrium-specific: {cancer_specific['Endometrium_Specific'].sum()}")
        print(f"  - Pan-cancer: {cancer_specific['Pan_Cancer'].sum()}\n")
    
    print("STEP 5: Generating Summary Reports")
    print("─"*70)
    
    # Create comprehensive summary
    summary_by_gene = df[df['In_COSMIC_Census']].groupby(sym_c).agg({
        sam_c: 'nunique',
        var_c: 'count',
        'Cancer_Relevance': lambda x: (x == 'High').sum()
    }).reset_index()
    summary_by_gene.columns = ['Gene', 'Unique_Samples', 'Total_Variants', 'High_Relevance_Variants']
    
    # Add COSMIC info
    summary_by_gene['COSMIC_Role'] = summary_by_gene['Gene'].map(
        lambda x: COSMIC_CENSUS_GENES.get(x, {}).get('Role', '')
    )
    summary_by_gene['COSMIC_Tier'] = summary_by_gene['Gene'].map(
        lambda x: COSMIC_CENSUS_GENES.get(x, {}).get('Tier', np.nan)
    )
    
    summary_by_gene = summary_by_gene.sort_values('Total_Variants', ascending=False)
    
    # Save results
    with pd.ExcelWriter(OUTPUT_XLSX, engine='openpyxl') as writer:
        # Summary sheet
        summary_by_gene.to_excel(writer, sheet_name='COSMIC_Gene_Summary', index=False)
        
        # All COSMIC variants
        cosmic_vars = df[df['In_COSMIC_Census']].copy()
        cosmic_vars.to_excel(writer, sheet_name='All_COSMIC_Variants', index=False)
        
        # High cancer relevance only
        high_relevance = df[df['Cancer_Relevance'] == 'High'].copy()
        if not high_relevance.empty:
            high_relevance.to_excel(writer, sheet_name='High_Cancer_Relevance', index=False)
        
        # Hotspots
        if not hotspots.empty:
            hotspots.to_excel(writer, sheet_name='Mutational_Hotspots', index=False)
        
        # Cancer-specific variants
        if 'Breast' in cancer_specific.columns:
            cancer_specific.to_excel(writer, sheet_name='Cancer_Specificity', index=False)
    
    print(f"✓ Results saved to: {OUTPUT_XLSX}\n")
    
    # Generate visualizations
    generate_cosmic_plots(df, summary_by_gene, hotspots, cancer_specific,
                         sym_c, tis_c, coh_c, imp_c)
    
    # Print final summary
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\nTotal variants analyzed: {len(df)}")
    print(f"Variants in COSMIC census genes: {len(census_variants)} "
          f"({len(census_variants)/len(df)*100:.1f}%)")
    print(f"\nCancer relevance breakdown:")
    for level in ['High', 'Medium', 'Low', 'Unknown']:
        count = (df['Cancer_Relevance'] == level).sum()
        print(f"  {level}: {count} ({count/len(df)*100:.1f}%)")
    
    print(f"\nTop 5 most mutated COSMIC genes:")
    for _, row in summary_by_gene.head(5).iterrows():
        print(f"  {row['Gene']} ({row['COSMIC_Role']}): "
              f"{row['Total_Variants']} variants in {row['Unique_Samples']} samples")
    
    print("\n" + "="*70)
    print("✓ COSMIC integration analysis complete!")
    print("="*70 + "\n")


def generate_cosmic_plots(df, summary_by_gene, hotspots, cancer_specific,
                         sym_c, tis_c, coh_c, imp_c='IMPACT'):
    """Generate comprehensive COSMIC analysis visualizations."""
    print("Generating COSMIC analysis plots...")
    
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Plot 1: COSMIC Census genes - mutation burden
    ax1 = fig.add_subplot(gs[0, :2])
    top_genes = summary_by_gene.head(10)
    colors_tier = {1: '#e74c3c', 2: '#f39c12', 3: '#95a5a6'}
    bar_colors = [colors_tier.get(tier, '#95a5a6') for tier in top_genes['COSMIC_Tier']]
    
    bars = ax1.barh(range(len(top_genes)), top_genes['Total_Variants'], color=bar_colors)
    ax1.set_yticks(range(len(top_genes)))
    ax1.set_yticklabels([f"{gene} ({role})" 
                         for gene, role in zip(top_genes['Gene'], top_genes['COSMIC_Role'])])
    ax1.set_xlabel("Total Variants", fontweight='bold')
    ax1.set_title("Top 10 COSMIC Census Genes by Mutation Burden", 
                 fontweight='bold', fontsize=12)
    ax1.grid(True, alpha=0.3, axis='x')
    
    # Add tier legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e74c3c', label='Tier 1 (Strong evidence)'),
        Patch(facecolor='#f39c12', label='Tier 2 (Moderate evidence)'),
        Patch(facecolor='#95a5a6', label='Tier 3 (Limited evidence)')
    ]
    ax1.legend(handles=legend_elements, loc='lower right', fontsize=9)
    
    # Plot 2: Cancer Relevance Distribution
    ax2 = fig.add_subplot(gs[0, 2])
    relevance_counts = df['Cancer_Relevance'].value_counts()
    colors_rel = {'High': '#e74c3c', 'Medium': '#f39c12', 
                 'Low': '#3498db', 'Unknown': '#95a5a6'}
    colors_list = [colors_rel.get(cat, '#95a5a6') for cat in relevance_counts.index]
    
    ax2.pie(relevance_counts.values, labels=relevance_counts.index, 
           autopct='%1.1f%%', colors=colors_list, startangle=90)
    ax2.set_title("Cancer Relevance\nClassification", fontweight='bold', fontsize=11)
    
    # Plot 3: Somatic vs Germline by Impact
    ax3 = fig.add_subplot(gs[1, 0])
    impact_class = df.groupby(['Variant_Class', imp_c]).size().unstack(fill_value=0)
    impact_class = impact_class[['HIGH', 'MODERATE', 'MODIFIER', 'LOW']]
    impact_class.plot(kind='bar', stacked=True, ax=ax3,
                     color=['#d9534f', '#5bc0de', '#f0ad4e', '#5cb85c'])
    ax3.set_xlabel("Variant Classification", fontweight='bold')
    ax3.set_ylabel("Variant Count", fontweight='bold')
    ax3.set_title("Somatic vs Germline by Impact", fontweight='bold', fontsize=11)
    plt.setp(ax3.xaxis.get_majorticklabels(), rotation=0)
    ax3.legend(title="Impact", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    
    # Plot 4: Hotspot Recurrence
    ax4 = fig.add_subplot(gs[1, 1])
    if not hotspots.empty:
        top_hotspots = hotspots.nsmallest(10, 'Recurrence_Rate_%')
        ax4.barh(range(len(top_hotspots)), top_hotspots['Recurrence_Rate_%'],
                color='#e74c3c', alpha=0.7)
        ax4.set_yticks(range(len(top_hotspots)))
        ax4.set_yticklabels([f"{row['Gene']}:{row['Position']}" 
                            for _, row in top_hotspots.iterrows()], fontsize=8)
        ax4.set_xlabel("Recurrence Rate (%)", fontweight='bold')
        ax4.set_title("Top Mutational Hotspots", fontweight='bold', fontsize=11)
        ax4.grid(True, alpha=0.3, axis='x')
    else:
        ax4.text(0.5, 0.5, 'No hotspots identified', 
                ha='center', va='center', fontsize=12)
        ax4.axis('off')
    
    # Plot 5: Cancer Type Specificity
    ax5 = fig.add_subplot(gs[1, 2])
    if 'Breast_Specific' in cancer_specific.columns:
        specificity_counts = {
            'Breast-specific': cancer_specific['Breast_Specific'].sum(),
            'Endometrium-specific': cancer_specific['Endometrium_Specific'].sum(),
            'Pan-cancer': cancer_specific['Pan_Cancer'].sum()
        }
        colors_spec = ['#e74c3c', '#2ecc71', '#3498db']
        ax5.bar(specificity_counts.keys(), specificity_counts.values(), 
               color=colors_spec, alpha=0.7, edgecolor='black')
        ax5.set_ylabel("Number of Variants", fontweight='bold')
        ax5.set_title("Cancer Type Specificity", fontweight='bold', fontsize=11)
        plt.setp(ax5.xaxis.get_majorticklabels(), rotation=15)
        ax5.grid(True, alpha=0.3, axis='y')
    
    # Plot 6: COSMIC genes across cohorts
    ax6 = fig.add_subplot(gs[2, :])
    cosmic_by_cohort = df[df['In_COSMIC_Census']].groupby([sym_c, coh_c]).size().unstack(fill_value=0)
    
    if len(cosmic_by_cohort) > 0:
        # Show top 15 genes
        top_cosmic = cosmic_by_cohort.sum(axis=1).nlargest(15)
        plot_data = cosmic_by_cohort.loc[top_cosmic.index]
        
        plot_data.plot(kind='barh', stacked=True, ax=ax6, 
                      color=['#e74c3c', '#2ecc71', '#3498db'])
        ax6.set_xlabel("Variant Count", fontweight='bold')
        ax6.set_ylabel("Gene", fontweight='bold')
        ax6.set_title("COSMIC Census Genes Across Cancer Types", 
                     fontweight='bold', fontsize=12)
        ax6.legend(title="Cohort", bbox_to_anchor=(1.02, 1), loc='upper left')
        ax6.grid(True, alpha=0.3, axis='x')
    
    plt.suptitle('COSMIC Cancer Gene Census Integration Analysis', 
                fontsize=16, fontweight='bold', y=0.995)
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ COSMIC plots saved to: {OUTPUT_PLOT}\n")


if __name__ == "__main__":
    analyze_cosmic_context()
