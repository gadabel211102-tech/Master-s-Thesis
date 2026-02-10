#!/usr/bin/env python3
"""
Unified Technical Audit for Ion Torrent Amplicon Sequencing
============================================================
Analyzes sequence complexity, coverage patterns, and technical failures
across multiple cohorts to identify problematic amplicons.

Author: [Your name]
Date: February 2026
"""

import pandas as pd
import numpy as np
import glob
import os
import re
import argparse
import logging
from pathlib import Path
from typing import Tuple, List, Dict, Optional
from Bio import SeqIO
from Bio.Seq import Seq

# ============================================================================
# CONFIGURATION
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
# SEQUENCE ANALYSIS
# ============================================================================

def analyse_sequence_complexity(seq: str) -> Tuple[float, float, float, float, int, int, int]:
    """
    Calculates technical failure metrics for Ion Torrent sequencing.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Tuple of (GC%, GC_skew, Max_Window_GC, Min_Window_GC, 
                  Max_Homopolymer, Max_SSR_Run, Hairpin_Risk)
    """
    if not seq: 
        return tuple([np.nan] * 7)
    
    seq = str(seq).upper()
    L = len(seq)
    
    if L == 0:
        return tuple([np.nan] * 7)
    
    # GC Metrics
    g, c = seq.count('G'), seq.count('C')
    gc_pct = (g + c) / L * 100
    gc_skew = (g - c) / (g + c) if (g + c) > 0 else 0
    
    # Sliding Window GC (20bp windows for local GC extremes)
    window_size = min(20, L)  # Handle short sequences
    if L >= window_size:
        windows = [seq[i:i+window_size] for i in range(L - window_size + 1)]
        window_gcs = [(w.count('G') + w.count('C')) / window_size * 100 for w in windows]
        max_window_gc = max(window_gcs)
        min_window_gc = min(window_gcs)
    else:
        max_window_gc = min_window_gc = gc_pct
    
    # Homopolymers (Ion Torrent struggles with long runs)
    homopolymer_runs = re.findall(r'([A-Z])\1*', seq)
    max_homo = max([len(run) for run in homopolymer_runs]) if homopolymer_runs else 0
    
    # SSRs (Simple Sequence Repeats - di/tri-nucleotide repeats)
    multi_repeats = re.findall(r'([A-Z]{2,3})\1+', seq)
    max_ssr_run = max([len(r) for r in multi_repeats]) if multi_repeats else 0

    # Hairpin Risk (self-complementary 6-mers can cause secondary structure)
    hairpin_risk = 0
    if L >= 12:  # Need at least 6bp stem + 6bp gap
        for i in range(L - 11):
            stem = seq[i:i+6]
            downstream = seq[i+6:]
            if str(Seq(stem).reverse_complement()) in downstream:
                hairpin_risk = 1
                break

    return (round(gc_pct, 2), round(gc_skew, 3), round(max_window_gc, 1), 
            round(min_window_gc, 1), max_homo, max_ssr_run, hairpin_risk)

# ============================================================================
# DATA LOADING
# ============================================================================

def load_bed_file(bed_path: str, logger: logging.Logger) -> pd.DataFrame:
    """
    Load and process BED file with amplicon coordinates.
    
    Args:
        bed_path: Path to BED file
        logger: Logger instance
        
    Returns:
        DataFrame with parsed BED data
    """
    logger.info(f"Loading BED file: {bed_path}")
    
    bed = pd.read_csv(
        bed_path, 
        sep='\t', 
        skiprows=1, 
        header=None, 
        names=['chr', 'start', 'end', 'numeric_id', 'gene_snp', 'extra']
    )
    
    # Create unique match ID for joining data
    bed['match_id'] = bed['chr'] + ":" + bed['start'].astype(str) + "-" + bed['end'].astype(str)
    bed['length'] = bed['end'] - bed['start']
    
    logger.info(f"Loaded {len(bed)} amplicons")
    return bed

def extract_sequence_features(bed: pd.DataFrame, fasta_path: str, logger: logging.Logger) -> pd.DataFrame:
    """
    Extract sequence complexity features from reference genome.
    
    Args:
        bed: DataFrame with amplicon coordinates
        fasta_path: Path to reference FASTA
        logger: Logger instance
        
    Returns:
        DataFrame with added sequence complexity columns
    """
    if not os.path.exists(fasta_path):
        logger.warning(f"FASTA not found: {fasta_path} - skipping sequence analysis")
        return bed
    
    logger.info("Extracting sequence complexity metrics from FASTA...")
    
    # Load entire genome into memory (one-time cost)
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    
    results = []
    failed_extractions = 0
    
    for idx, row in bed.iterrows():
        try:
            sequence = str(genome[row['chr']].seq[row['start']:row['end']])
            metrics = analyse_sequence_complexity(sequence)
            
            results.append({
                'match_id': row['match_id'], 
                'GC_pct': metrics[0], 
                'GC_skew': metrics[1],
                'Max_Window_GC': metrics[2], 
                'Min_Window_GC': metrics[3],
                'Max_Homopolymer': metrics[4], 
                'Max_SSR_Run': metrics[5], 
                'Hairpin_Risk': metrics[6]
            })
        except KeyError:
            # Chromosome not in reference
            logger.warning(f"Chromosome {row['chr']} not found in reference")
            results.append({'match_id': row['match_id']})
            failed_extractions += 1
        except Exception as e:
            # Other extraction errors
            logger.debug(f"Failed to extract sequence for {row['match_id']}: {e}")
            results.append({'match_id': row['match_id']})
            failed_extractions += 1
    
    if failed_extractions > 0:
        logger.warning(f"Failed to extract {failed_extractions}/{len(bed)} sequences")
    
    bed = bed.merge(pd.DataFrame(results), on='match_id', how='left')
    logger.info(f"Added sequence complexity metrics")
    
    return bed

# ============================================================================
# COVERAGE ANALYSIS
# ============================================================================

def process_coverage_data(cov_pattern: str, logger: logging.Logger) -> Tuple[pd.DataFrame, List[Dict]]:
    """
    Process per-sample coverage files and compute statistics.
    
    Args:
        cov_pattern: Glob pattern for coverage files
        logger: Logger instance
        
    Returns:
        Tuple of (master_coverage_df, sample_stats_list)
    """
    all_cov_files = glob.glob(cov_pattern, recursive=True)
    
    if not all_cov_files:
        logger.warning(f"No coverage files found matching: {cov_pattern}")
        return pd.DataFrame(), []
    
    logger.info(f"Processing {len(all_cov_files)} coverage files...")
    
    cov_rows = []
    sample_stats = []
    failed_files = 0
    
    for idx, f in enumerate(all_cov_files, 1):
        sample_name = os.path.basename(f).replace('.coverage.tsv.gz', '')
        cohort = os.path.basename(os.path.dirname(f))
        
        if idx % 50 == 0:
            logger.debug(f"  Processed {idx}/{len(all_cov_files)} samples...")
        
        try:
            df = pd.read_csv(f, sep='\t', compression='gzip')
            
            # Validate expected columns
            if 'depth' not in df.columns or 'id' not in df.columns:
                logger.warning(f"Missing required columns in {f}")
                failed_files += 1
                continue
            
            sample_stats.append({
                'Sample': sample_name, 
                'Cohort': cohort, 
                'Avg_Depth': df['depth'].mean(),
                'Median_Depth': df['depth'].median(),
                'Min_Depth': df['depth'].min(),
                'Max_Depth': df['depth'].max(),
                'Zero_Coverage_Amplicons': (df['depth'] == 0).sum()
            })
            
            df['is_failed'] = (df['depth'] == 0).astype(int)
            df['cohort'] = cohort
            cov_rows.append(df[['id', 'is_failed', 'depth', 'cohort']])
            
        except Exception as e:
            logger.warning(f"Failed to process {f}: {e}")
            failed_files += 1
            continue
    
    if failed_files > 0:
        logger.warning(f"Failed to process {failed_files}/{len(all_cov_files)} files")
    
    if not cov_rows:
        logger.error("No coverage data successfully loaded!")
        return pd.DataFrame(), sample_stats
    
    master_cov = pd.concat(cov_rows, ignore_index=True)
    logger.info(f"Combined coverage data: {len(master_cov):,} records from {len(cov_rows)} samples")
    
    return master_cov, sample_stats

def add_coverage_metrics(bed: pd.DataFrame, master_cov: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Add global and per-cohort coverage metrics to BED dataframe.
    
    Args:
        bed: DataFrame with amplicon data
        master_cov: Combined coverage data
        logger: Logger instance
        
    Returns:
        DataFrame with added coverage metrics
    """
    if master_cov.empty:
        logger.warning("No coverage data to add")
        return bed
    
    logger.info("Computing coverage statistics...")
    
    # Global statistics across all samples
    global_fail = master_cov.groupby('id')['is_failed'].mean() * 100
    global_depth = master_cov.groupby('id')['depth'].mean()
    global_median_depth = master_cov.groupby('id')['depth'].median()
    
    # Per-cohort statistics
    cohort_fail = (master_cov.groupby(['id', 'cohort'])['is_failed'].mean() * 100).unstack().add_suffix('_fail_%')
    cohort_depth = master_cov.groupby(['id', 'cohort'])['depth'].mean().unstack().add_suffix('_avg_depth')
    
    # Add to BED
    bed['Global_Failure_%'] = bed['match_id'].map(global_fail)
    bed['Global_Mean_Depth'] = bed['match_id'].map(global_depth)
    bed['Global_Median_Depth'] = bed['match_id'].map(global_median_depth)
    
    bed = bed.merge(cohort_fail, left_on='match_id', right_index=True, how='left')
    bed = bed.merge(cohort_depth, left_on='match_id', right_index=True, how='left')
    
    logger.info(f"Added global and {len(cohort_fail.columns)} per-cohort metrics")
    
    return bed

# ============================================================================
# ZERO-COVERAGE BASE ANALYSIS
# ============================================================================

def process_zero_coverage_reports(zero_pattern: str, logger: logging.Logger) -> Optional[pd.DataFrame]:
    """
    Merge zero-base coverage reports from multiple cohorts.
    
    Args:
        zero_pattern: Glob pattern for zero-coverage TSV files
        logger: Logger instance
        
    Returns:
        DataFrame with per-cohort empty base counts, or None if no files found
    """
    all_zero_files = glob.glob(zero_pattern, recursive=True)
    
    if not all_zero_files:
        logger.warning(f"No zero-coverage files found matching: {zero_pattern}")
        return None
    
    logger.info(f"Merging {len(all_zero_files)} zero-coverage reports...")
    
    merged_zero = None
    failed_merges = 0
    
    for z_file in all_zero_files:
        cohort_id = os.path.basename(z_file).replace('zero_cov_', '').replace('.tsv', '')
        
        try:
            z_df = pd.read_csv(z_file, sep='\t')
            
            # Validate structure (first 4 columns should be chr, start, end, annotation)
            if len(z_df.columns) < 5:
                logger.warning(f"Insufficient columns in {z_file}")
                failed_merges += 1
                continue
            
            # Average zero bases across all samples (columns 4+)
            sample_cols = z_df.columns[4:]
            z_df[f'{cohort_id}_empty_bases'] = z_df[sample_cols].mean(axis=1)
            
            # Create match ID
            z_df['match_id'] = z_df['chr'] + ":" + z_df['start'].astype(str) + "-" + z_df['end'].astype(str)
            
            temp = z_df[['match_id', f'{cohort_id}_empty_bases']]
            
            if merged_zero is None:
                merged_zero = temp
            else:
                merged_zero = pd.merge(merged_zero, temp, on='match_id', how='outer')
                
        except Exception as e:
            logger.warning(f"Failed to process {z_file}: {e}")
            failed_merges += 1
            continue
    
    if failed_merges > 0:
        logger.warning(f"Failed to merge {failed_merges}/{len(all_zero_files)} zero-coverage files")
    
    return merged_zero

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def run_unified_audit(
    bed_path: str,
    fasta_path: str,
    cov_pattern: str,
    zero_pattern: str,
    output_dir: str = ".",
    verbose: bool = False
) -> None:
    """
    Execute the unified technical audit pipeline.
    
    Args:
        bed_path: Path to BED file with amplicon coordinates
        fasta_path: Path to reference genome FASTA
        cov_pattern: Glob pattern for per-sample coverage files
        zero_pattern: Glob pattern for zero-coverage reports
        output_dir: Directory for output files
        verbose: Enable debug logging
    """
    logger = setup_logging(verbose)
    logger.info("=" * 70)
    logger.info("UNIFIED TECHNICAL AUDIT - Ion Torrent Amplicon QC")
    logger.info("=" * 70)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Load BED file
    bed = load_bed_file(bed_path, logger)
    
    # 2. Extract sequence features
    bed = extract_sequence_features(bed, fasta_path, logger)
    
    # 3. Process coverage data
    master_cov, sample_stats = process_coverage_data(cov_pattern, logger)
    
    if sample_stats:
        sample_df = pd.DataFrame(sample_stats)
        sample_output = os.path.join(output_dir, "Sample_Average_Coverage.csv")
        sample_df.to_csv(sample_output, index=False)
        logger.info(f"Saved sample statistics: {sample_output}")
    
    if not master_cov.empty:
        bed = add_coverage_metrics(bed, master_cov, logger)
    
    # 4. Process zero-coverage base counts
    merged_zero = process_zero_coverage_reports(zero_pattern, logger)
    
    if merged_zero is not None:
        bed = bed.merge(merged_zero, on='match_id', how='left').fillna(0)
        logger.info(f"Added zero-base metrics for {len(merged_zero.columns)-1} cohorts")
    
    # 5. Export final audit
    output_file = os.path.join(output_dir, "TFM_GSDMB_Final_Audit.csv")
    bed.sort_values(by='Global_Failure_%', ascending=False, na_position='last').to_csv(output_file, index=False)
    
    logger.info("=" * 70)
    logger.info("ANALYSIS SUMMARY")
    logger.info("=" * 70)
    
    # Print top 10 worst amplicons
    if 'Global_Failure_%' in bed.columns:
        logger.info("\nTOP 10 WORST AMPLICONS (by failure rate):")
        worst_10 = bed.sort_values(
            by=['Global_Failure_%', 'Global_Mean_Depth'], 
            ascending=[False, True],
            na_position='last'
        ).head(10)
        
        display_cols = ['numeric_id', 'gene_snp', 'length', 'Global_Failure_%', 'Global_Mean_Depth']
        print(worst_10[display_cols].to_string(index=False))
    
    # Summary statistics
    logger.info("\nOVERALL STATISTICS:")
    logger.info(f"  Total amplicons: {len(bed)}")
    
    if 'Global_Failure_%' in bed.columns:
        logger.info(f"  Mean failure rate: {bed['Global_Failure_%'].mean():.2f}%")
        logger.info(f"  Amplicons with >50% failure: {(bed['Global_Failure_%'] > 50).sum()}")
        logger.info(f"  Amplicons with 100% failure: {(bed['Global_Failure_%'] == 100).sum()}")
    
    if 'Global_Mean_Depth' in bed.columns:
        logger.info(f"  Mean coverage depth: {bed['Global_Mean_Depth'].mean():.1f}x")
        logger.info(f"  Amplicons with <100x coverage: {(bed['Global_Mean_Depth'] < 100).sum()}")
    
    logger.info(f"\nOUTPUT FILES:")
    logger.info(f"  Main audit:      {output_file}")
    if sample_stats:
        logger.info(f"  Sample stats:    {sample_output}")
    
    logger.info("\n✅ ANALYSIS COMPLETE")

# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

def main():
    """Command-line interface for the audit script."""
    parser = argparse.ArgumentParser(
        description='Unified technical audit for Ion Torrent amplicon sequencing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with hardcoded paths
  python 13_failure.py
  
  # Custom paths
  python 13_failure.py \\
    --bed /path/to/targets.bed \\
    --fasta /path/to/genome.fa \\
    --coverage "/data/**/coverage/*.tsv.gz" \\
    --zero "/data/**/zero_cov_*.tsv" \\
    --output results/
  
  # Verbose output
  python 13_failure.py -v
        """
    )
    
    parser.add_argument(
        '--bed', 
        default="/home/gadeaalonsoj/tfm/IAD255368_167_Submitted.bed",
        help='Path to BED file with amplicon coordinates'
    )
    
    parser.add_argument(
        '--fasta',
        default="/home/gadeaalonsoj/ref_alt/hg38_alt.fa",
        help='Path to reference genome FASTA'
    )
    
    parser.add_argument(
        '--coverage',
        default="/home/gadeaalonsoj/tfm/gsdmb_final_results/individual_sample_coverage/**/*.coverage.tsv.gz",
        help='Glob pattern for per-sample coverage files (quote the pattern!)'
    )
    
    parser.add_argument(
        '--zero',
        default="/home/gadeaalonsoj/tfm/**/zero_cov_fast/zero_cov_*.tsv",
        help='Glob pattern for zero-coverage reports (quote the pattern!)'
    )
    
    parser.add_argument(
        '--output', '-o',
        default=".",
        help='Output directory for results (default: current directory)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose debug logging'
    )
    
    args = parser.parse_args()
    
    run_unified_audit(
        bed_path=args.bed,
        fasta_path=args.fasta,
        cov_pattern=args.coverage,
        zero_pattern=args.zero,
        output_dir=args.output,
        verbose=args.verbose
    )

if __name__ == "__main__":
    main()