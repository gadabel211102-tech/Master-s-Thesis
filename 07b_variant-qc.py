#!/usr/bin/env python3
"""
Variant-Level Quality Control Script (Enhanced)
================================================
Performs comprehensive QC checks on called variants including:
- Ti/Tv ratio (transition/transversion)
- Allelic balance for heterozygous calls
- Batch effects detection
- Sequencing run date clustering
- Color-coded threshold visualization for amplicon QC

Expected to run AFTER variant calling (script 05) and BEFORE annotation (script 06)

"""

import pandas as pd
import numpy as np
import glob
import os
import sys
import argparse
import logging
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle

# ============================================================================
# QC THRESHOLDS CONFIGURATION
# ============================================================================

QC_THRESHOLDS = {
    'titv': {
        'optimal_min': 2.0,
        'optimal_max': 2.1,
        'acceptable_min': 1.5,
        'acceptable_max': 3.0
    },
    'allelic_balance': {
        'lower': 0.25,
        'upper': 0.75,
        'optimal_lower': 0.4,
        'optimal_upper': 0.6
    },
    'batch_cv': {
        'good': 15,
        'acceptable': 30,
        'poor': 50
    },
    'depth': {
        'min_acceptable': 20,
        'min_good': 50,
        'optimal': 100
    },
    'quality': {
        'min_acceptable': 20,
        'min_good': 30,
        'optimal': 40
    }
}

# ============================================================================
# LOGGING SETUP
# ============================================================================

def setup_logging(verbose: bool = False) -> logging.Logger:
    """Configure logging with appropriate level."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%H:%M:%S'
    )
    return logging.getLogger(__name__)

# ============================================================================
# TI/TV RATIO CALCULATION
# ============================================================================

def calculate_titv_ratio(vcf_df: pd.DataFrame, logger: logging.Logger) -> dict:
    """
    Calculate Transition/Transversion ratio from variants.
    Expected value for WES/targeted panels: ~2.0-2.1
    Lower values suggest poor quality calls.
    
    Args:
        vcf_df: DataFrame with CHROM, POS, REF, ALT columns
        logger: Logger instance
        
    Returns:
        Dictionary with Ti/Tv metrics
    """
    transitions = {'A>G', 'G>A', 'C>T', 'T>C'}
    transversions = {'A>C', 'C>A', 'A>T', 'T>A', 'G>C', 'C>G', 'G>T', 'T>G'}
    
    ti_count = 0
    tv_count = 0
    invalid = 0
    
    # Case-insensitive column lookup (columns may be REF/ref/Ref etc.)
    ref_col = next((c for c in vcf_df.columns if c.upper() == 'REF'), None)
    alt_col = next((c for c in vcf_df.columns if c.upper() == 'ALT'), None)
    if ref_col is None or alt_col is None:
        logger.warning("âš ï¸  REF or ALT column not found - cannot calculate Ti/Tv")
        return {'ti_count': 0, 'tv_count': 0, 'ratio': float('nan'),
                'qc_status': 'FAIL_NO_REF_ALT', 'qc_color': '#e74c3c'}

    for _, row in vcf_df.iterrows():
        ref = str(row.get(ref_col, '')).upper().strip()
        alt = str(row.get(alt_col, '')).upper().strip()
        
        # Skip if not SNV (single nucleotide variant)
        if len(ref) != 1 or len(alt) != 1:
            continue
        
        # Skip if not standard bases
        if ref not in 'ACGT' or alt not in 'ACGT':
            invalid += 1
            continue
        
        mutation = f"{ref}>{alt}"
        
        if mutation in transitions:
            ti_count += 1
        elif mutation in transversions:
            tv_count += 1
        else:
            invalid += 1
    
    if tv_count == 0:
        logger.warning("âš ï¸  No transversions found - Ti/Tv calculation may be unreliable")
        return {
            'ti_count': ti_count,
            'tv_count': tv_count,
            'ratio': np.nan,
            'qc_status': 'FAIL_NO_TV',
            'qc_color': '#e74c3c'
        }
    
    ratio = ti_count / tv_count
    
    # QC interpretation with color coding
    thresholds = QC_THRESHOLDS['titv']
    if ratio < thresholds['acceptable_min']:
        status = 'FAIL_LOW_TITV'
        color = '#e74c3c'  # Red
        logger.warning(f"âš ï¸  Low Ti/Tv ratio ({ratio:.3f}) - expected ~2.0-2.1")
        logger.warning(f"     Possible causes: poor sequencing quality, PCR errors, or contamination")
    elif ratio > thresholds['acceptable_max']:
        status = 'WARN_HIGH_TITV'
        color = '#f39c12'  # Orange
        logger.warning(f"âš ï¸  High Ti/Tv ratio ({ratio:.3f}) - possible over-filtering or bias")
    elif thresholds['optimal_min'] <= ratio <= thresholds['optimal_max']:
        status = 'OPTIMAL'
        color = '#2ecc71'  # Green
        logger.info(f"âœ… Ti/Tv Ratio: {ratio:.3f} (optimal range)")
    else:
        status = 'ACCEPTABLE'
        color = '#3498db'  # Blue
        logger.info(f"âœ… Ti/Tv Ratio: {ratio:.3f} (acceptable range)")
    
    logger.info(f"   Transitions: {ti_count}, Transversions: {tv_count}")
    
    return {
        'ti_count': ti_count,
        'tv_count': tv_count,
        'ratio': ratio,
        'qc_status': status,
        'qc_color': color
    }

# ============================================================================
# ALLELIC BALANCE CHECK
# ============================================================================

def check_allelic_balance(vcf_df: pd.DataFrame, logger: logging.Logger, 
                          ab_lower: float = None, ab_upper: float = None,
                          sample_col: str = None) -> pd.DataFrame:
    """
    Flag heterozygous variants with abnormal allelic balance.
    Expected: ~50% for true heterozygotes (0.25-0.75 acceptable range)
    Deviations suggest: contamination, LOH, or copy number changes
    
    Args:
        vcf_df: DataFrame with genotype info
        logger: Logger instance
        ab_lower: Lower threshold for allelic fraction (default from QC_THRESHOLDS)
        ab_upper: Upper threshold for allelic fraction (default from QC_THRESHOLDS)
        sample_col: Sample column name if multi-sample VCF
        
    Returns:
        DataFrame with AB_Flag and AF_calc columns added
    """
    # Use configured thresholds if not provided
    if ab_lower is None:
        ab_lower = QC_THRESHOLDS['allelic_balance']['lower']
    if ab_upper is None:
        ab_upper = QC_THRESHOLDS['allelic_balance']['upper']
    
    logger.info("Checking allelic balance for heterozygous calls...")
    logger.info(f"   Thresholds: {ab_lower:.2f} - {ab_upper:.2f}")
    
    vcf_df['AB_Flag'] = 'PASS'
    vcf_df['AF_calc'] = np.nan
    vcf_df['AB_Quality'] = 'UNKNOWN'  # New column for color coding
    flagged_count = 0
    het_count = 0
    
    # Determine which columns are available (case-insensitive)
    col_upper = {c.upper(): c for c in vcf_df.columns}
    has_ao_ro = 'AO' in col_upper and 'RO' in col_upper
    has_af    = 'AF' in col_upper
    has_dp    = 'DP' in col_upper
    # Remap to actual column names so .get() works below
    if has_ao_ro:
        vcf_df = vcf_df.rename(columns={col_upper['AO']: 'AO', col_upper['RO']: 'RO'})
    if has_af:
        vcf_df = vcf_df.rename(columns={col_upper['AF']: 'AF'})
    if 'GT' in col_upper:
        vcf_df = vcf_df.rename(columns={col_upper['GT']: 'GT'})
    if has_dp:
        vcf_df = vcf_df.rename(columns={col_upper['DP']: 'DP'})
    
    if not has_ao_ro and not has_af:
        logger.warning("âš ï¸  No AO/RO or AF columns found - skipping allelic balance check")
        logger.warning("     This check requires variant caller output with allele depths")
        return vcf_df
    
    optimal_lower = QC_THRESHOLDS['allelic_balance']['optimal_lower']
    optimal_upper = QC_THRESHOLDS['allelic_balance']['optimal_upper']
    
    # Calculate AF from AO/RO if available
    if has_ao_ro:
        for idx, row in vcf_df.iterrows():
            # Get genotype - handle both phased (|) and unphased (/) separators
            gt = str(row.get('GT', '.')).replace('|', '/').split('/')[0] if '/' in str(row.get('GT', '.')).replace('|', '/') else str(row.get('GT', '.'))
            
            # Only check heterozygotes (0/1, 1/0, etc.)
            if '/' not in str(row.get('GT', '.')).replace('|', '/'):
                continue
            
            gt_parts = str(row.get('GT', '.')).replace('|', '/').split('/')
            if len(gt_parts) != 2:
                continue
                
            if gt_parts[0] == gt_parts[1]:  # Homozygous
                continue
            
            het_count += 1
            
            try:
                ao = float(row.get('AO', 0))
                ro = float(row.get('RO', 0))
                total = ao + ro
                
                if total > 0:
                    af = ao / total
                    vcf_df.at[idx, 'AF_calc'] = af
                    
                    # Assign quality level based on thresholds
                    if af < ab_lower or af > ab_upper:
                        vcf_df.at[idx, 'AB_Flag'] = f'IMBALANCED_AF={af:.3f}'
                        vcf_df.at[idx, 'AB_Quality'] = 'POOR'
                        flagged_count += 1
                    elif optimal_lower <= af <= optimal_upper:
                        vcf_df.at[idx, 'AB_Quality'] = 'OPTIMAL'
                    else:
                        vcf_df.at[idx, 'AB_Quality'] = 'ACCEPTABLE'
            except (ValueError, TypeError):
                continue
                
    elif has_af:
        for idx, row in vcf_df.iterrows():
            gt = str(row.get('GT', '.')).replace('|', '/')
            
            if '/' not in gt:
                continue
                
            gt_parts = gt.split('/')
            if len(gt_parts) != 2 or gt_parts[0] == gt_parts[1]:
                continue
            
            het_count += 1
            
            try:
                af = float(row.get('AF', 0.5))
                vcf_df.at[idx, 'AF_calc'] = af
                
                if af < ab_lower or af > ab_upper:
                    vcf_df.at[idx, 'AB_Flag'] = f'IMBALANCED_AF={af:.3f}'
                    vcf_df.at[idx, 'AB_Quality'] = 'POOR'
                    flagged_count += 1
                elif optimal_lower <= af <= optimal_upper:
                    vcf_df.at[idx, 'AB_Quality'] = 'OPTIMAL'
                else:
                    vcf_df.at[idx, 'AB_Quality'] = 'ACCEPTABLE'
            except (ValueError, TypeError):
                continue
    
    if het_count == 0:
        logger.warning("âš ï¸  No heterozygous calls found")
    else:
        pct_flagged = (flagged_count / het_count) * 100
        optimal_count = vcf_df[vcf_df['AB_Quality'] == 'OPTIMAL'].shape[0]
        acceptable_count = vcf_df[vcf_df['AB_Quality'] == 'ACCEPTABLE'].shape[0]
        
        logger.info(f"   Heterozygous variants checked: {het_count}")
        logger.info(f"   Optimal allelic balance: {optimal_count} ({optimal_count/het_count*100:.1f}%)")
        logger.info(f"   Acceptable: {acceptable_count} ({acceptable_count/het_count*100:.1f}%)")
        logger.info(f"   Flagged for allelic imbalance: {flagged_count} ({pct_flagged:.1f}%)")
        
        if pct_flagged > 20:
            logger.warning(f"âš ï¸  High percentage of allelic imbalance ({pct_flagged:.1f}%)")
            logger.warning(f"     Possible causes: LOH, contamination, CNV, or subclonal mutations")
    
    return vcf_df

# ============================================================================
# BATCH EFFECTS DETECTION
# ============================================================================

def detect_batch_effects(variant_summary: pd.DataFrame, logger: logging.Logger,
                        cohort_col: str = 'Cohort', 
                        tissue_col: str = 'Tissue',
                        sample_col: str = 'Sample') -> pd.DataFrame:
    """
    Detect potential batch effects by analyzing variant count distributions
    across cohorts and tissue types.
    
    Args:
        variant_summary: DataFrame with cohort/tissue/sample info
        logger: Logger instance
        cohort_col: Column name for cohort
        tissue_col: Column name for tissue type
        sample_col: Column name for sample ID
        
    Returns:
        DataFrame with batch effect statistics
    """
    logger.info("Analyzing potential batch effects...")
    
    # Count variants per sample
    sample_counts = variant_summary.groupby([cohort_col, tissue_col, sample_col]).size().reset_index(name='Variant_Count')
    
    # Calculate statistics per group
    batch_stats = []
    thresholds = QC_THRESHOLDS['batch_cv']
    
    for (cohort, tissue), group in sample_counts.groupby([cohort_col, tissue_col]):
        mean_count = group['Variant_Count'].mean()
        median_count = group['Variant_Count'].median()
        std_count = group['Variant_Count'].std()
        cv = (std_count / mean_count * 100) if mean_count > 0 else 0
        
        # Assign quality level based on CV thresholds
        if cv < thresholds['good']:
            qc_flag = 'EXCELLENT'
            qc_color = '#2ecc71'  # Green
        elif cv < thresholds['acceptable']:
            qc_flag = 'ACCEPTABLE'
            qc_color = '#3498db'  # Blue
        elif cv < thresholds['poor']:
            qc_flag = 'HIGH_VARIATION'
            qc_color = '#f39c12'  # Orange
        else:
            qc_flag = 'VERY_HIGH_VARIATION'
            qc_color = '#e74c3c'  # Red
        
        batch_stats.append({
            'Group': f"{cohort}_{tissue}",
            'Cohort': cohort,
            'Tissue': tissue,
            'N_Samples': len(group),
            'Mean_Variants': mean_count,
            'Median_Variants': median_count,
            'Std_Variants': std_count,
            'CV_%': cv,
            'QC_Flag': qc_flag,
            'QC_Color': qc_color
        })
    
    batch_df = pd.DataFrame(batch_stats)
    
    # Report findings
    logger.info("\n" + "="*60)
    logger.info("BATCH EFFECT ANALYSIS")
    logger.info("="*60)
    print(batch_df[['Group', 'N_Samples', 'Mean_Variants', 'CV_%', 'QC_Flag']].to_string(index=False))
    
    high_cv_groups = batch_df[batch_df['CV_%'] > thresholds['acceptable']]
    if not high_cv_groups.empty:
        logger.warning(f"\nâš ï¸  {len(high_cv_groups)} groups with high variation (CV > {thresholds['acceptable']}%)")
        logger.warning("   This may indicate batch effects or sample quality issues")
    else:
        logger.info(f"\nâœ… No significant batch effects detected (all CV < {thresholds['acceptable']}%)")
    
    return batch_df

# ============================================================================
# VISUALIZATION WITH THRESHOLD HIGHLIGHTING
# ============================================================================

def plot_qc_summary(titv_results: dict, ab_df: pd.DataFrame, batch_df: pd.DataFrame,
                    output_dir: str, logger: logging.Logger):
    """
    Create summary visualizations of QC metrics with color-coded thresholds.
    
    Args:
        titv_results: Ti/Tv ratio results
        ab_df: Allelic balance DataFrame
        batch_df: Batch effect statistics
        output_dir: Output directory
        logger: Logger instance
    """
    logger.info("Generating QC summary plots with threshold highlighting...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Ti/Tv Ratio with Threshold Zones
    ax1 = axes[0, 0]
    if not np.isnan(titv_results['ratio']):
        categories = ['Transitions', 'Transversions']
        counts = [titv_results['ti_count'], titv_results['tv_count']]
        colors = ['#3498db', '#e74c3c']
        
        ax1.bar(categories, counts, color=colors, edgecolor='black', linewidth=1.5)
        
        # Color-coded title based on QC status
        title_color = titv_results['qc_color']
        ax1.set_title(f"Ti/Tv Ratio: {titv_results['ratio']:.3f}\nStatus: {titv_results['qc_status']}", 
                     fontweight='bold', fontsize=12, color=title_color)
        ax1.set_ylabel('Count')
        
        # Add threshold reference zones
        thresholds = QC_THRESHOLDS['titv']
        ax1.text(0.5, 0.95, 
                f'Optimal: {thresholds["optimal_min"]:.1f}-{thresholds["optimal_max"]:.1f} | '
                f'Acceptable: {thresholds["acceptable_min"]:.1f}-{thresholds["acceptable_max"]:.1f}', 
                transform=ax1.transAxes, ha='center', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                fontsize=9)
    else:
        ax1.text(0.5, 0.5, 'Ti/Tv Calculation Failed', 
                transform=ax1.transAxes, ha='center', va='center',
                fontsize=14, color='red')
    
    # 2. Allelic Balance Distribution with Threshold Zones
    ax2 = axes[0, 1]
    if 'AF_calc' in ab_df.columns and ab_df['AF_calc'].notna().any():
        het_af = ab_df[ab_df['AF_calc'].notna()]['AF_calc']
        
        # Create histogram with color-coded regions
        thresholds_ab = QC_THRESHOLDS['allelic_balance']
        
        # Add colored background zones
        ax2.axvspan(0, thresholds_ab['lower'], alpha=0.2, color='#e74c3c', label='Poor')
        ax2.axvspan(thresholds_ab['lower'], thresholds_ab['optimal_lower'], 
                   alpha=0.2, color='#3498db', label='Acceptable')
        ax2.axvspan(thresholds_ab['optimal_lower'], thresholds_ab['optimal_upper'], 
                   alpha=0.2, color='#2ecc71', label='Optimal')
        ax2.axvspan(thresholds_ab['optimal_upper'], thresholds_ab['upper'], 
                   alpha=0.2, color='#3498db')
        ax2.axvspan(thresholds_ab['upper'], 1.0, alpha=0.2, color='#e74c3c')
        
        ax2.hist(het_af, bins=30, color='#9b59b6', edgecolor='black', alpha=0.7, zorder=3)
        ax2.axvline(0.5, color='darkgreen', linestyle='--', linewidth=2, 
                   label='Expected (0.5)', zorder=4)
        
        ax2.set_xlabel('Allelic Fraction')
        ax2.set_ylabel('Count')
        ax2.set_title('Allelic Balance Distribution\n(Heterozygous Variants)', 
                     fontweight='bold', fontsize=12)
        ax2.legend(loc='upper right', fontsize=8)
        ax2.set_xlim(0, 1)
    else:
        ax2.text(0.5, 0.5, 'No Allelic Balance Data', 
                transform=ax2.transAxes, ha='center', va='center',
                fontsize=14)
    
    # 3. Variant Counts per Group with Color-Coded Bars
    ax3 = axes[1, 0]
    if not batch_df.empty:
        x_pos = np.arange(len(batch_df))
        
        # Use the QC_Color column from batch_df
        colors_batch = batch_df['QC_Color'].tolist()
        
        ax3.bar(x_pos, batch_df['Mean_Variants'], yerr=batch_df['Std_Variants'],
               color=colors_batch, edgecolor='black', linewidth=1.5, alpha=0.7,
               capsize=5)
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(batch_df['Group'], rotation=45, ha='right')
        ax3.set_ylabel('Mean Variant Count')
        ax3.set_title('Variant Distribution Across Groups\n(Error bars = Std Dev)', 
                     fontweight='bold', fontsize=12)
        ax3.grid(axis='y', alpha=0.3)
        
        # Add legend for color coding
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#2ecc71', label='Excellent (CV<15%)'),
            Patch(facecolor='#3498db', label='Acceptable (CV<30%)'),
            Patch(facecolor='#f39c12', label='High Var (CV<50%)'),
            Patch(facecolor='#e74c3c', label='Very High (CVâ‰¥50%)')
        ]
        ax3.legend(handles=legend_elements, loc='upper right', fontsize=8)
    else:
        ax3.text(0.5, 0.5, 'No Batch Data', 
                transform=ax3.transAxes, ha='center', va='center',
                fontsize=14)
    
    # 4. Coefficient of Variation with Threshold Lines
    ax4 = axes[1, 1]
    if not batch_df.empty:
        thresholds_cv = QC_THRESHOLDS['batch_cv']
        
        # Color bars based on CV thresholds
        colors_cv = batch_df['QC_Color'].tolist()
        
        ax4.barh(batch_df['Group'], batch_df['CV_%'], color=colors_cv, 
                edgecolor='black', linewidth=1.5, alpha=0.7)
        
        # Add threshold lines
        ax4.axvline(thresholds_cv['good'], color='#2ecc71', linestyle='--', 
                   linewidth=2, label=f'Good (<{thresholds_cv["good"]}%)')
        ax4.axvline(thresholds_cv['acceptable'], color='#f39c12', linestyle='--', 
                   linewidth=2, label=f'Acceptable (<{thresholds_cv["acceptable"]}%)')
        ax4.axvline(thresholds_cv['poor'], color='#e74c3c', linestyle='--', 
                   linewidth=2, label=f'Poor (<{thresholds_cv["poor"]}%)')
        
        ax4.set_xlabel('Coefficient of Variation (%)')
        ax4.set_title('Sample Variation Within Groups\n(CV% - lower is better)', 
                     fontweight='bold', fontsize=12)
        ax4.legend(loc='lower right', fontsize=8)
        ax4.grid(axis='x', alpha=0.3)
    else:
        ax4.text(0.5, 0.5, 'No CV Data', 
                transform=ax4.transAxes, ha='center', va='center',
                fontsize=14)
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, "Variant_QC_Summary_Enhanced.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"âœ… Enhanced QC summary plot saved: {output_file}")

# ============================================================================
# AMPLICON-SPECIFIC QC (NEW FUNCTION)
# ============================================================================

def _assign_amplicons_from_bed(df: pd.DataFrame, bed_path: str,
                               logger: logging.Logger) -> pd.DataFrame:
    """
    Assign each variant to an amplicon by positional overlap with the BED file.

    BED columns expected: chr, start (0-based), end (0-based exclusive), annotation/id.
    Each variant gains an 'Amplicon_ID' column. Variants not overlapping any amplicon
    receive 'No_Amplicon'.
    """
    try:
        bed = pd.read_csv(bed_path, sep='\t', header=None,
                          usecols=[0, 1, 2, 3],
                          names=['bed_chr', 'bed_start', 'bed_end', 'bed_id'])
        bed['bed_chr'] = bed['bed_chr'].astype(str).apply(
            lambda c: c if c.startswith('chr') else f'chr{c}'
        )
        logger.info(f"Loaded {len(bed)} amplicons from BED: {bed_path}")
    except Exception as e:
        logger.warning(f"Could not load BED file ({e}) - falling back to SYMBOL grouping")
        return df

    # VCF POS is 1-based; BED is 0-based half-open.
    # Overlap: bed_start < POS <= bed_end
    amplicon_ids = []
    for _, var in df.iterrows():
        chrom = str(var.get('CHROM', ''))
        pos   = var.get('POS', None)
        if pd.isna(pos):
            amplicon_ids.append('No_Amplicon')
            continue
        pos = int(pos)
        match = bed[
            (bed['bed_chr'] == chrom) &
            (bed['bed_start'] < pos) &
            (bed['bed_end']   >= pos)
        ]
        if len(match) == 1:
            amplicon_ids.append(str(match.iloc[0]['bed_id']))
        elif len(match) > 1:
            # Overlapping amplicons - use the smallest (most specific) one
            match = match.copy()
            match['span'] = match['bed_end'] - match['bed_start']
            amplicon_ids.append(str(match.loc[match['span'].idxmin(), 'bed_id']))
        else:
            amplicon_ids.append('No_Amplicon')

    df = df.copy()
    df['Amplicon_ID'] = amplicon_ids
    assigned = (df['Amplicon_ID'] != 'No_Amplicon').sum()
    logger.info(f"Assigned {assigned}/{len(df)} variant rows to amplicons via BED overlap")
    return df


def analyze_amplicon_metrics(df: pd.DataFrame, logger: logging.Logger,
                             output_dir: str,
                             bed_path: str = None) -> pd.DataFrame:
    """
    Analyze amplicon-specific metrics with threshold-based color coding.

    If a BED file path is provided, variants are assigned to individual amplicons
    by positional overlap (true amplicon-level resolution - one bar per amplicon).
    Otherwise falls back to grouping by SYMBOL (gene-level, as before).

    Args:
        df:         DataFrame with variant annotations (needs CHROM, POS, DP columns)
        logger:     Logger instance
        output_dir: Output directory for plots
        bed_path:   Optional path to targets BED file for amplicon-level grouping

    Returns:
        DataFrame with per-amplicon QC statistics
    """
    logger.info("Analyzing amplicon-specific QC metrics...")

    # Find DP column case-insensitively
    dp_col = next((c for c in df.columns if c.upper() == 'DP'), None)
    if dp_col is None:
        logger.warning("No depth (DP) column found - skipping amplicon QC")
        return pd.DataFrame()
    if dp_col != 'DP':
        df = df.rename(columns={dp_col: 'DP'})

    # ------------------------------------------------------------------
    # Determine grouping column
    # ------------------------------------------------------------------
    if bed_path is not None:
        df = _assign_amplicons_from_bed(df, bed_path, logger)
        if 'Amplicon_ID' in df.columns and (df['Amplicon_ID'] != 'No_Amplicon').any():
            amplicon_col = 'Amplicon_ID'
            logger.info("Grouping by individual amplicon ID (BED positional overlap)")
        else:
            logger.warning("BED overlap produced no assignments - falling back to SYMBOL")
            amplicon_col = 'SYMBOL' if 'SYMBOL' in df.columns else None
    else:
        amplicon_col = None
        for col in ['Amplicon', 'Region', 'Gene', 'SYMBOL']:
            if col in df.columns:
                amplicon_col = col
                break

    if amplicon_col is None:
        logger.warning("No amplicon/region identifier found - using global stats")
        amplicon_col = 'Global'
        df = df.copy()
        df['Global'] = 'All_Variants'

    if bed_path is None:
        logger.warning(
            f"No --bed file provided: grouping by '{amplicon_col}' (gene-level). "
            "Pass --bed targets.sorted.bed for true per-amplicon resolution."
        )
    
    # Calculate statistics per amplicon
    amplicon_stats = []
    thresholds = QC_THRESHOLDS['depth']
    
    for amplicon, group in df.groupby(amplicon_col):
        depths = group['DP'].dropna()
        
        if len(depths) == 0:
            continue
        
        mean_depth = depths.mean()
        median_depth = depths.median()
        min_depth = depths.min()
        max_depth = depths.max()
        std_depth = depths.std()
        
        # Assign quality level based on mean depth
        if mean_depth >= thresholds['optimal']:
            qc_flag = 'EXCELLENT'
            qc_color = '#2ecc71'  # Green
        elif mean_depth >= thresholds['min_good']:
            qc_flag = 'GOOD'
            qc_color = '#3498db'  # Blue
        elif mean_depth >= thresholds['min_acceptable']:
            qc_flag = 'ACCEPTABLE'
            qc_color = '#f39c12'  # Orange
        else:
            qc_flag = 'POOR'
            qc_color = '#e74c3c'  # Red
        
        amplicon_stats.append({
            'Amplicon': amplicon,
            'N_Variants': len(group),
            'Mean_Depth': mean_depth,
            'Median_Depth': median_depth,
            'Min_Depth': min_depth,
            'Max_Depth': max_depth,
            'Std_Depth': std_depth,
            'QC_Flag': qc_flag,
            'QC_Color': qc_color
        })
    
    amplicon_df = pd.DataFrame(amplicon_stats)
    
    if amplicon_df.empty:
        return amplicon_df
    
    # Plot amplicon metrics with color coding
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # 1. Mean Depth per Amplicon
    ax1 = axes[0]
    x_pos = np.arange(len(amplicon_df))
    colors = amplicon_df['QC_Color'].tolist()
    
    ax1.bar(x_pos, amplicon_df['Mean_Depth'], color=colors, 
           edgecolor='black', linewidth=1.5, alpha=0.7)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(amplicon_df['Amplicon'], rotation=45, ha='right')
    ax1.set_ylabel('Mean Depth')
    ax1.set_title('Mean Coverage Depth per Amplicon/Region', 
                 fontweight='bold', fontsize=12)
    
    # Add threshold lines
    ax1.axhline(thresholds['min_acceptable'], color='#e74c3c', linestyle='--', 
               linewidth=2, label=f'Min Acceptable ({thresholds["min_acceptable"]}x)')
    ax1.axhline(thresholds['min_good'], color='#3498db', linestyle='--', 
               linewidth=2, label=f'Good ({thresholds["min_good"]}x)')
    ax1.axhline(thresholds['optimal'], color='#2ecc71', linestyle='--', 
               linewidth=2, label=f'Optimal ({thresholds["optimal"]}x)')
    
    ax1.legend(loc='upper right')
    ax1.grid(axis='y', alpha=0.3)
    
    # 2. Depth Distribution with Threshold Zones
    ax2 = axes[1]
    
    # Add colored background zones
    ax2.axhspan(0, thresholds['min_acceptable'], alpha=0.2, color='#e74c3c', label='Poor')
    ax2.axhspan(thresholds['min_acceptable'], thresholds['min_good'], 
               alpha=0.2, color='#f39c12', label='Acceptable')
    ax2.axhspan(thresholds['min_good'], thresholds['optimal'], 
               alpha=0.2, color='#3498db', label='Good')
    ax2.axhspan(thresholds['optimal'], df['DP'].max() if 'DP' in df.columns else 200, 
               alpha=0.2, color='#2ecc71', label='Excellent')
    
    if 'DP' in df.columns:
        ax2.hist(df['DP'].dropna(), bins=50, color='#9b59b6', 
                edgecolor='black', alpha=0.7, zorder=3)
    
    ax2.set_xlabel('Depth (DP)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Overall Depth Distribution\n(All Variants)', 
                 fontweight='bold', fontsize=12)
    ax2.legend(loc='upper right')
    ax2.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, "Amplicon_QC_Metrics.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"âœ… Amplicon QC plot saved: {output_file}")
    
    # Print summary
    logger.info("\n" + "="*60)
    logger.info("AMPLICON QC SUMMARY")
    logger.info("="*60)
    print(amplicon_df[['Amplicon', 'N_Variants', 'Mean_Depth', 'QC_Flag']].to_string(index=False))
    
    poor_amplicons = amplicon_df[amplicon_df['Mean_Depth'] < thresholds['min_acceptable']]
    if not poor_amplicons.empty:
        logger.warning(f"\nâš ï¸  {len(poor_amplicons)} amplicons with poor coverage (<{thresholds['min_acceptable']}x)")
        logger.warning("   Consider re-sequencing or excluding these regions")
    else:
        logger.info(f"\nâœ… All amplicons meet minimum coverage threshold (â‰¥{thresholds['min_acceptable']}x)")
    
    return amplicon_df

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive variant-level QC analysis with threshold highlighting',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLE USAGE:
  # Process annotated variants
  python 07b_variant-qc.py \\
    --input /path/to/GSDMB_Annotated_Report_Fixed.xlsx \\
    --bed /path/to/targets.sorted.bed \\
    --output /path/to/qc_results \\
    --verbose

  # With custom thresholds
  python 05b_variant_qc_enhanced.py \\
    --input annotations.xlsx \\
    --ab-lower 0.2 \\
    --ab-upper 0.8 \\
    --min-depth 30
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Path to annotated variants Excel file (from script 07)'
    )

    parser.add_argument(
        '--bed',
        default=None,
        help='Path to targets BED file for true amplicon-level grouping by positional '
             'overlap. If not provided, falls back to grouping by SYMBOL.'
    )

    parser.add_argument(
        '--output', '-o',
        default='.',
        help='Output directory for QC reports'
    )
    
    parser.add_argument(
        '--ab-lower',
        type=float,
        default=None,
        help=f'Lower threshold for allelic balance (default: {QC_THRESHOLDS["allelic_balance"]["lower"]})'
    )
    
    parser.add_argument(
        '--ab-upper',
        type=float,
        default=None,
        help=f'Upper threshold for allelic balance (default: {QC_THRESHOLDS["allelic_balance"]["upper"]})'
    )
    
    parser.add_argument(
        '--min-depth',
        type=int,
        default=None,
        help=f'Minimum acceptable depth (default: {QC_THRESHOLDS["depth"]["min_acceptable"]})'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose logging'
    )
    
    args = parser.parse_args()
    logger = setup_logging(args.verbose)
    
    # Update thresholds if provided
    if args.min_depth is not None:
        QC_THRESHOLDS['depth']['min_acceptable'] = args.min_depth
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    logger.info("="*70)
    logger.info("VARIANT-LEVEL QUALITY CONTROL ANALYSIS (ENHANCED)")
    logger.info("="*70)
    logger.info(f"Input file: {args.input}")
    logger.info(f"Output dir: {args.output}")
    logger.info("="*70)
    
    # Load annotated variants
    logger.info("\nLoading annotated variants...")
    try:
        df = pd.read_excel(args.input, sheet_name='Biological_Annotations')
        logger.info(f"âœ… Loaded {len(df)} variants")
    except Exception as e:
        logger.error(f"Failed to load input file: {e}")
        sys.exit(1)
    
    # 1. Calculate Ti/Tv Ratio
    logger.info("\n" + "="*70)
    logger.info("1. CALCULATING Ti/Tv RATIO")
    logger.info("="*70)
    titv_results = calculate_titv_ratio(df, logger)
    
    # 2. Check Allelic Balance
    logger.info("\n" + "="*70)
    logger.info("2. CHECKING ALLELIC BALANCE")
    logger.info("="*70)
    df = check_allelic_balance(df, logger, args.ab_lower, args.ab_upper)
    
    # 3. Detect Batch Effects
    logger.info("\n" + "="*70)
    logger.info("3. DETECTING BATCH EFFECTS")
    logger.info("="*70)
    batch_df = detect_batch_effects(df, logger)
    
    # 4. Analyze Amplicon Metrics (NEW)
    logger.info("\n" + "="*70)
    logger.info("4. ANALYZING AMPLICON METRICS")
    logger.info("="*70)
    amplicon_df = analyze_amplicon_metrics(df, logger, args.output, bed_path=args.bed)
    
    # 5. Generate Visualizations
    logger.info("\n" + "="*70)
    logger.info("5. GENERATING QC VISUALIZATIONS")
    logger.info("="*70)
    plot_qc_summary(titv_results, df, batch_df, args.output, logger)
    
    # 6. Save Results
    logger.info("\n" + "="*70)
    logger.info("6. SAVING QC RESULTS")
    logger.info("="*70)
    
    # Save Ti/Tv results
    titv_df = pd.DataFrame([titv_results])
    titv_file = os.path.join(args.output, "TiTv_Ratio_Report.csv")
    titv_df.to_csv(titv_file, index=False)
    logger.info(f"âœ… Ti/Tv report saved: {titv_file}")
    
    # Save batch analysis
    batch_file = os.path.join(args.output, "Batch_Effects_Report.csv")
    batch_df.to_csv(batch_file, index=False)
    logger.info(f"âœ… Batch effects report saved: {batch_file}")
    
    # Save amplicon analysis
    if not amplicon_df.empty:
        amplicon_file = os.path.join(args.output, "Amplicon_QC_Report.csv")
        amplicon_df.to_csv(amplicon_file, index=False)
        logger.info(f"âœ… Amplicon QC report saved: {amplicon_file}")
    
    # Save updated annotations with AB flags
    output_excel = os.path.join(args.output, "GSDMB_Annotated_Report_QC.xlsx")
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Biological_Annotations', index=False)
    logger.info(f"âœ… Updated annotations saved: {output_excel}")
    
    # Final Summary
    logger.info("\n" + "="*70)
    logger.info("QC SUMMARY")
    logger.info("="*70)
    logger.info(f"Ti/Tv Ratio: {titv_results['ratio']:.3f} ({titv_results['qc_status']})")
    
    if 'AB_Flag' in df.columns:
        flagged = df[df['AB_Flag'] != 'PASS'].shape[0]
        optimal = df[df['AB_Quality'] == 'OPTIMAL'].shape[0]
        logger.info(f"Allelic Balance: {optimal} optimal, {flagged} flagged")
    
    high_cv = batch_df[batch_df['CV_%'] > QC_THRESHOLDS['batch_cv']['acceptable']].shape[0]
    logger.info(f"Batch Effects: {high_cv} groups with high variation")
    
    if not amplicon_df.empty:
        poor_amp = amplicon_df[amplicon_df['Mean_Depth'] < QC_THRESHOLDS['depth']['min_acceptable']].shape[0]
        logger.info(f"Amplicon QC: {poor_amp} regions below minimum coverage")
    
    logger.info("\n" + "="*70)
    logger.info("âœ… VARIANT QC COMPLETE")
    logger.info("="*70)

if __name__ == "__main__":
    main()