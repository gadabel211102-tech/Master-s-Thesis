#!/usr/bin/env python3
"""
Unified Amplicon Technical Audit
=================================
Combines comprehensive technical analysis with cohort-specific reporting.

Two modes:
1. COMPREHENSIVE: Full technical audit with sequence complexity
2. QUICK: Excel report with worst 10 amplicons per cohort


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
# SHARED UTILITIES
# ============================================================================

def normalise_chr(c: str) -> str:
    """
    Normalize chromosome names: '17' -> 'chr17', keep existing 'chr' prefix.
    
    Args:
        c: Chromosome name
        
    Returns:
        Normalized chromosome name with 'chr' prefix
    """
    c = str(c).strip()
    if c.startswith("chr"):
        return c
    if c in {"X", "Y", "M", "MT"}:
        return "chrM" if c in {"M", "MT"} else f"chr{c}"
    try:
        return f"chr{int(c)}"
    except Exception:
        return c if c.startswith("chr") else f"chr{c}"

def clean_sheet_name(name: str, max_length: int = 31) -> str:
    """
    Clean sheet name for Excel compatibility.
    
    Args:
        name: Original name
        max_length: Maximum length (Excel limit is 31)
        
    Returns:
        Cleaned name
    """
    # Remove invalid Excel characters
    cleaned = re.sub(r'[\\/*?:\[\]]', '', name)
    return cleaned[:max_length]

def load_bed_annotations(bed_path: str, logger: logging.Logger) -> pd.DataFrame:
    """
    Load BED file and extract annotations.
    
    Columns expected:
    0: chr, 1: start, 2: end, 3: numeric_id, 4: rs_id/gene_snp, 5: gene (optional)
    
    Args:
        bed_path: Path to BED file
        logger: Logger instance
        
    Returns:
        DataFrame with chr, start, end, annotation, match_id
    """
    logger.info(f"Loading BED annotations: {bed_path}")
    
    try:
        # Try reading with flexible column detection
        bed_df = pd.read_csv(
            bed_path, 
            sep='\t', 
            comment='t',
            skiprows=1,  # Skip track line if present
            header=None
        )
        
        # Map columns based on what we have
        if len(bed_df.columns) >= 6:
            bed_df.columns = ['chr', 'start', 'end', 'numeric_id', 'gene_snp', 'extra'] + \
                            [f'col{i}' for i in range(6, len(bed_df.columns))]
        elif len(bed_df.columns) >= 5:
            bed_df.columns = ['chr', 'start', 'end', 'numeric_id', 'gene_snp'] + \
                            [f'col{i}' for i in range(5, len(bed_df.columns))]
        else:
            bed_df.columns = ['chr', 'start', 'end', 'numeric_id'] + \
                            [f'col{i}' for i in range(4, len(bed_df.columns))]
        
        # Normalize chromosome names
        bed_df['chr'] = bed_df['chr'].astype(str).apply(normalise_chr)
        
        # Create annotation - prefer rsIDs, filter out non-informative values
        if 'gene_snp' in bed_df.columns:
            # Start with gene_snp column
            bed_df['annotation'] = bed_df['gene_snp'].astype(str)
            
            # Replace non-informative values with NaN
            non_informative = ['.', 'in', 'ex', 'intron', 'exon', 'intergenic', 
                             'intronic', 'exonic', 'nan', 'none', '']
            bed_df['annotation'] = bed_df['annotation'].replace(
                non_informative, np.nan
            )
            
            # Filter to only keep values starting with 'rs' (rsIDs)
            # If not an rsID, set to NaN
            bed_df['annotation'] = bed_df['annotation'].apply(
                lambda x: x if pd.notna(x) and str(x).lower().startswith('rs') else np.nan
            )
            
            # Fall back to 'extra' column if available
            if 'extra' in bed_df.columns:
                bed_df['annotation'] = bed_df['annotation'].fillna(bed_df['extra'])
            
            # Final fallback: create descriptive ID from coordinates
            bed_df['annotation'] = bed_df['annotation'].fillna(
                bed_df['chr'] + ':' + bed_df['start'].astype(str) + '-' + bed_df['end'].astype(str)
            )
        else:
            # No gene_snp column - use coordinates
            bed_df['annotation'] = (bed_df['chr'] + ':' + 
                                   bed_df['start'].astype(str) + '-' + 
                                   bed_df['end'].astype(str))
        
        # Create match_id for joining
        bed_df['match_id'] = (bed_df['chr'] + ":" + 
                             bed_df['start'].astype(str) + "-" + 
                             bed_df['end'].astype(str))
        
        # Calculate length
        bed_df['length'] = bed_df['end'] - bed_df['start']
        
        logger.info(f"Loaded {len(bed_df)} amplicons from BED")
        
        return bed_df[['chr', 'start', 'end', 'annotation', 'match_id', 'length', 'numeric_id']]
        
    except Exception as e:
        logger.error(f"Failed to load BED file: {e}")
        raise

# ============================================================================
# SEQUENCE COMPLEXITY ANALYSIS
# ============================================================================

def analyse_sequence_complexity(seq: str) -> Tuple[float, float, float, float, int, int, int]:
    """
    Calculate technical failure metrics for Ion Torrent sequencing.
    
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
    
    # Sliding Window GC
    window_size = min(20, L)
    if L >= window_size:
        windows = [seq[i:i+window_size] for i in range(L - window_size + 1)]
        window_gcs = [(w.count('G') + w.count('C')) / window_size * 100 for w in windows]
        max_window_gc = max(window_gcs)
        min_window_gc = min(window_gcs)
    else:
        max_window_gc = min_window_gc = gc_pct
    
    # Homopolymers
    homopolymer_runs = re.findall(r'([A-Z])\1*', seq)
    max_homo = max([len(run) for run in homopolymer_runs]) if homopolymer_runs else 0
    
    # SSRs
    multi_repeats = re.findall(r'([A-Z]{2,3})\1+', seq)
    max_ssr_run = max([len(r) for r in multi_repeats]) if multi_repeats else 0

    # Hairpin Risk
    hairpin_risk = 0
    if L >= 12:
        for i in range(L - 11):
            stem = seq[i:i+6]
            downstream = seq[i+6:]
            if str(Seq(stem).reverse_complement()) in downstream:
                hairpin_risk = 1
                break

    return (round(gc_pct, 2), round(gc_skew, 3), round(max_window_gc, 1), 
            round(min_window_gc, 1), max_homo, max_ssr_run, hairpin_risk)

def extract_sequence_features(bed: pd.DataFrame, fasta_path: str, logger: logging.Logger) -> pd.DataFrame:
    """Extract sequence complexity features from reference genome."""
    if not os.path.exists(fasta_path):
        logger.warning(f"FASTA not found: {fasta_path} - skipping sequence analysis")
        return bed
    
    logger.info("Extracting sequence complexity metrics...")
    
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    results = []
    failed = 0
    
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
        except Exception as e:
            logger.debug(f"Failed sequence extraction for {row['match_id']}: {e}")
            results.append({'match_id': row['match_id']})
            failed += 1
    
    if failed > 0:
        logger.warning(f"Failed to extract {failed}/{len(bed)} sequences")
    
    return bed.merge(pd.DataFrame(results), on='match_id', how='left')

# ============================================================================
# COVERAGE DATA PROCESSING
# ============================================================================

def process_individual_coverage(cov_pattern: str, logger: logging.Logger) -> Tuple[pd.DataFrame, List[Dict]]:
    """
    Process per-sample coverage files.
    
    Args:
        cov_pattern: Glob pattern for coverage files
        logger: Logger instance
        
    Returns:
        Tuple of (master_coverage_df, sample_stats_list)
    """
    all_cov_files = glob.glob(cov_pattern, recursive=True)
    
    if not all_cov_files:
        logger.warning(f"No coverage files found: {cov_pattern}")
        return pd.DataFrame(), []
    
    logger.info(f"Processing {len(all_cov_files)} individual coverage files...")
    
    cov_rows = []
    sample_stats = []
    failed = 0
    
    for idx, f in enumerate(all_cov_files, 1):
        sample_name = os.path.basename(f).replace('.coverage.tsv.gz', '')
        cohort = os.path.basename(os.path.dirname(f))
        
        if idx % 50 == 0:
            logger.debug(f"  Processed {idx}/{len(all_cov_files)}...")
        
        try:
            df = pd.read_csv(f, sep='\t', compression='gzip')
            
            if 'depth' not in df.columns or 'id' not in df.columns:
                logger.warning(f"Missing columns in {f}")
                failed += 1
                continue
            
            sample_stats.append({
                'Sample': sample_name,
                'Cohort': cohort,
                'Avg_Depth': df['depth'].mean(),
                'Median_Depth': df['depth'].median(),
                'Zero_Coverage_Amplicons': (df['depth'] == 0).sum()
            })
            
            df['is_failed'] = (df['depth'] == 0).astype(int)
            df['cohort'] = cohort
            cov_rows.append(df[['id', 'is_failed', 'depth', 'cohort']])
            
        except Exception as e:
            logger.warning(f"Failed to process {f}: {e}")
            failed += 1
    
    if failed:
        logger.warning(f"Failed: {failed}/{len(all_cov_files)} files")
    
    if cov_rows:
        return pd.concat(cov_rows, ignore_index=True), sample_stats
    return pd.DataFrame(), sample_stats

def process_panel_coverage(panel_pattern: str, logger: logging.Logger) -> List[Tuple[str, pd.DataFrame]]:
    """
    Process GSDMB_*_FullPanel.csv files for cohort-specific analysis.
    
    Args:
        panel_pattern: Glob pattern for panel files
        logger: Logger instance
        
    Returns:
        List of (cohort_name, dataframe) tuples
    """
    panel_files = glob.glob(panel_pattern, recursive=True)
    
    if not panel_files:
        logger.warning(f"No panel files found: {panel_pattern}")
        return []
    
    logger.info(f"Loading {len(panel_files)} panel coverage files...")
    
    cohort_data = []
    
    for f in panel_files:
        cohort_name = (
            os.path.basename(f)
            .replace("GSDMB_", "")
            .replace("_FullPanel.csv", "")
        )
        
        try:
            df = pd.read_csv(f)
            
            # Normalize chromosomes if present
            if 'chr' in df.columns:
                df['chr'] = df['chr'].astype(str).apply(normalise_chr)
            
            # Ensure numeric types for coordinates
            for col in ['start', 'end']:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
            
            cohort_data.append((cohort_name, df))
            logger.debug(f"  Loaded {cohort_name}: {len(df)} records")
            
        except Exception as e:
            logger.warning(f"Failed to load {f}: {e}")
    
    return cohort_data

# ============================================================================
# ZERO COVERAGE ANALYSIS
# ============================================================================

def process_zero_coverage(zero_pattern: str, logger: logging.Logger) -> Optional[pd.DataFrame]:
    """
    Merge zero-coverage reports from multiple cohorts.

    Handles the updated output format from 02b_more-qc.sh, which produces
    four columns per sample:
        <sample>_zero_bases  - number of bases with 0x coverage
        <sample>_pct_zero    - % of amplicon with 0x coverage
        <sample>_low_bases   - number of bases with <200x coverage
        <sample>_pct_low     - % of amplicon with <200x coverage

    For each cohort file, this function computes per-amplicon summary metrics
    averaged across all samples:
        <cohort>_mean_zero_bases  - mean bases with 0x coverage across samples
        <cohort>_mean_pct_zero    - mean % of amplicon with 0x coverage
        <cohort>_mean_low_bases   - mean bases with <200x coverage
        <cohort>_mean_pct_low     - mean % of amplicon with <200x coverage
        <cohort>_n_samples_zero   - number of samples with ANY 0x bases
        <cohort>_n_samples_low    - number of samples with ANY <200x bases
    """
    all_zero_files = glob.glob(zero_pattern, recursive=True)

    if not all_zero_files:
        logger.warning(f"No zero-coverage files found: {zero_pattern}")
        return None

    logger.info(f"Merging {len(all_zero_files)} zero-coverage reports...")

    merged_zero = None
    failed = 0

    # Fixed metadata columns present in every file from 02b
    FIXED_COLS = {'chr', 'start', 'end', 'annotation', 'amplicon_length_bp'}

    for z_file in all_zero_files:
        cohort_id = os.path.basename(z_file).replace('zero_cov_', '').replace('.tsv', '')

        try:
            z_df = pd.read_csv(z_file, sep='\t')

            # Require at least the 5 fixed columns plus one set of sample metrics
            if len(z_df.columns) < 9:
                logger.warning(f"Insufficient columns in {z_file} - may be old format")
                failed += 1
                continue

            # Normalize chromosomes and build match_id for joining
            z_df['chr'] = z_df['chr'].astype(str).apply(normalise_chr)
            z_df['match_id'] = (z_df['chr'] + ":" +
                                z_df['start'].astype(str) + "-" +
                                z_df['end'].astype(str))

            # Identify the four groups of per-sample columns by suffix
            all_cols = set(z_df.columns)
            data_cols = all_cols - FIXED_COLS - {'match_id'}

            zero_base_cols = [c for c in data_cols if c.endswith('_zero_bases')]
            pct_zero_cols  = [c for c in data_cols if c.endswith('_pct_zero')]
            low_base_cols  = [c for c in data_cols if c.endswith('_low_bases')]
            pct_low_cols   = [c for c in data_cols if c.endswith('_pct_low')]

            if not zero_base_cols:
                logger.warning(f"No '_zero_bases' columns found in {z_file} - skipping")
                failed += 1
                continue

            logger.debug(f"  {cohort_id}: {len(zero_base_cols)} samples detected")

            # Compute per-amplicon summary metrics averaged across samples
            summary = z_df[['match_id']].copy()
            summary[f'{cohort_id}_mean_zero_bases'] = z_df[zero_base_cols].mean(axis=1).round(2)
            summary[f'{cohort_id}_mean_pct_zero']   = z_df[pct_zero_cols].mean(axis=1).round(2)
            summary[f'{cohort_id}_mean_low_bases']  = z_df[low_base_cols].mean(axis=1).round(2)
            summary[f'{cohort_id}_mean_pct_low']    = z_df[pct_low_cols].mean(axis=1).round(2)

            # Count how many samples had at least one affected base
            summary[f'{cohort_id}_n_samples_zero'] = (z_df[zero_base_cols] > 0).sum(axis=1)
            summary[f'{cohort_id}_n_samples_low']  = (z_df[low_base_cols] > 0).sum(axis=1)

            merged_zero = summary if merged_zero is None else pd.merge(
                merged_zero, summary, on='match_id', how='outer'
            )

        except Exception as e:
            logger.warning(f"Failed to process {z_file}: {e}")
            failed += 1

    if failed:
        logger.warning(f"Failed: {failed}/{len(all_zero_files)} files")

    return merged_zero

# ============================================================================
# COMPREHENSIVE AUDIT MODE
# ============================================================================

def run_comprehensive_audit(
    bed_path: str,
    fasta_path: str,
    cov_pattern: str,
    zero_pattern: str,
    output_dir: str,
    logger: logging.Logger
) -> pd.DataFrame:
    """
    Run comprehensive technical audit with sequence complexity.
    
    Returns the main audit DataFrame for potential reuse.
    """
    logger.info("=" * 70)
    logger.info("COMPREHENSIVE TECHNICAL AUDIT")
    logger.info("=" * 70)
    
    # Load BED
    bed = load_bed_annotations(bed_path, logger)
    
    # Sequence features
    bed = extract_sequence_features(bed, fasta_path, logger)
    
    # Individual coverage
    master_cov, sample_stats = process_individual_coverage(cov_pattern, logger)
    
    if sample_stats:
        sample_df = pd.DataFrame(sample_stats)
        sample_output = os.path.join(output_dir, "Sample_Average_Coverage.csv")
        sample_df.to_csv(sample_output, index=False)
        logger.info(f"Saved sample stats: {sample_output}")
    
    # Add coverage metrics
    if not master_cov.empty:
        logger.info("Computing coverage statistics...")
        
        global_fail = master_cov.groupby('id')['is_failed'].mean() * 100
        global_depth = master_cov.groupby('id')['depth'].mean()
        global_median = master_cov.groupby('id')['depth'].median()
        
        cohort_fail = (master_cov.groupby(['id', 'cohort'])['is_failed']
                      .mean() * 100).unstack().add_suffix('_fail_%')
        cohort_depth = (master_cov.groupby(['id', 'cohort'])['depth']
                       .mean()).unstack().add_suffix('_avg_depth')
        
        bed['Global_Failure_%'] = bed['match_id'].map(global_fail)
        bed['Global_Mean_Depth'] = bed['match_id'].map(global_depth)
        bed['Global_Median_Depth'] = bed['match_id'].map(global_median)
        
        bed = bed.merge(cohort_fail, left_on='match_id', right_index=True, how='left')
        bed = bed.merge(cohort_depth, left_on='match_id', right_index=True, how='left')
    
    # Zero coverage
    merged_zero = process_zero_coverage(zero_pattern, logger)
    if merged_zero is not None:
        bed = bed.merge(merged_zero, on='match_id', how='left').fillna(0)
        zero_cols_added = [c for c in merged_zero.columns if c != 'match_id']
        logger.info(f"Added zero/low-coverage metrics ({len(zero_cols_added)} columns): "
                    f"mean_zero_bases, mean_pct_zero, mean_low_bases, mean_pct_low, "
                    f"n_samples_zero, n_samples_low per cohort")
    
    # Export
    output_file = os.path.join(output_dir, "GSDMB_Technical_Audit.csv")
    bed.sort_values(by='Global_Failure_%', ascending=False, na_position='last').to_csv(
        output_file, index=False
    )
    
    logger.info(f"\nâœ… Comprehensive audit saved: {output_file}")
    
    # Summary
    if 'Global_Failure_%' in bed.columns:
        logger.info("\nTOP 10 WORST AMPLICONS:")
        worst_10 = bed.sort_values(
            by=['Global_Failure_%', 'Global_Mean_Depth'],
            ascending=[False, True],
            na_position='last'
        ).head(10)
        
        print(worst_10[['numeric_id', 'annotation', 'length', 
                       'Global_Failure_%', 'Global_Mean_Depth']].to_string(index=False))
    
    return bed

# ============================================================================
# QUICK EXCEL REPORT MODE
# ============================================================================

def run_quick_excel_report(
    bed_path: str,
    panel_pattern: str,
    output_dir: str,
    logger: logging.Logger
) -> None:
    """
    Generate Excel report with worst 10 amplicons per cohort.
    """
    logger.info("=" * 70)
    logger.info("QUICK EXCEL REPORT - WORST 10 PER COHORT")
    logger.info("=" * 70)
    
    # Load annotations
    annotations = load_bed_annotations(bed_path, logger)
    
    # Load panel data
    cohort_data = process_panel_coverage(panel_pattern, logger)
    
    if not cohort_data:
        logger.error("No panel data found!")
        return
    
    # Create Excel file
    output_file = os.path.join(output_dir, "GSDMB_Worst_Amplicons_PerCohort.xlsx")
    
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        used_sheets = set()
        
        for cohort_name, df in cohort_data:
            # Merge with annotations
            df_merged = pd.merge(
                df,
                annotations[['chr', 'start', 'end', 'annotation']],
                on=['chr', 'start', 'end'],
                how='left',
                suffixes=('_data', '_bed')
            )
            
            # Use BED annotation if available
            if 'annotation_bed' in df_merged.columns:
                df_merged['annotation'] = df_merged['annotation_bed'].fillna(
                    df_merged.get('annotation_data', '')
                )
            
            # Identify numeric sample columns
            metadata_cols = {'chr', 'start', 'end', 'id', 'annotation', 
                           'annotation_data', 'annotation_bed', 'match_id'}
            sample_cols = [c for c in df_merged.columns 
                          if c not in metadata_cols and 
                          pd.api.types.is_numeric_dtype(df_merged[c])]
            
            if not sample_cols:
                logger.warning(f"No numeric columns for {cohort_name}")
                continue
            
            # Calculate median and find worst 10
            df_merged['median_depth'] = df_merged[sample_cols].median(axis=1)
            worst_10 = df_merged.sort_values('median_depth', ascending=True).head(10)
            
            # Create unique sheet name
            base_sheet = clean_sheet_name(cohort_name)
            sheet_name = base_sheet
            i = 1
            while sheet_name in used_sheets:
                sheet_name = clean_sheet_name(f"{base_sheet}_{i}")
                i += 1
            used_sheets.add(sheet_name)
            
            # Export
            output_cols = ['chr', 'start', 'end', 'annotation', 'median_depth']
            worst_10[output_cols].to_excel(writer, sheet_name=sheet_name, index=False)
            
            logger.info(f"  Added sheet: {sheet_name}")
    
    logger.info(f"\nâœ… Excel report saved: {output_file}")

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description='Unified amplicon technical audit with two modes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
MODES:
  comprehensive   Full audit with sequence complexity (run after QC)
  quick          Excel report with worst 10 per cohort (quick review)
  both           Run both modes (recommended)

EXAMPLES:
  # Run both modes with defaults
  python unified_audit.py both
  
  # Just comprehensive audit
  python unified_audit.py comprehensive -v
  
  # Just Excel report
  python unified_audit.py quick
  
  # Custom paths
  python unified_audit.py both \\
    --bed /path/to/targets.bed \\
    --fasta /path/to/genome.fa \\
    --coverage "/data/**/coverage/*.tsv.gz" \\
    --panels "/data/GSDMB_*_FullPanel.csv"
        """
    )
    
    parser.add_argument(
        'mode',
        choices=['comprehensive', 'quick', 'both'],
        help='Analysis mode to run'
    )
    
    parser.add_argument(
        '--bed',
        default=str(Path(__file__).resolve().parent / "dna_bed/IAD255368_167_Submitted.bed"),
        help='Path to BED file'
    )
    
    parser.add_argument(
        '--fasta',
        default=str(Path(__file__).resolve().parent / "ref_alt/hg38_alt.fa"),
        help='Path to reference FASTA (for comprehensive mode)'
    )
    
    parser.add_argument(
        '--coverage',
        default=str(Path(__file__).resolve().parent / "analysis_results/03_qc_visualisation/sample_level_coverage_tables/**/*.coverage.tsv.gz"),
        help='Glob pattern for individual coverage files (for comprehensive mode)'
    )
    
    parser.add_argument(
        '--zero',
        default=str(Path(__file__).resolve().parent / "**/zero_coverage/zero_cov_*.tsv"),
        help='Glob pattern for zero-coverage reports (for comprehensive mode)'
    )
    
    parser.add_argument(
        '--panels',
        default=str(Path(__file__).resolve().parent / "analysis_results/03_qc_visualisation/GSDMB_*_FullPanel.csv"),
        help='Glob pattern for panel files (for quick mode)'
    )
    
    parser.add_argument(
        '--output', '-o',
        default=str(Path(__file__).resolve().parent / "analysis_results/04_technical_audit"),
        help='Output directory'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose logging'
    )
    
    args = parser.parse_args()
    logger = setup_logging(args.verbose)
    
    os.makedirs(args.output, exist_ok=True)
    
    if args.mode in ['comprehensive', 'both']:
        run_comprehensive_audit(
            bed_path=args.bed,
            fasta_path=args.fasta,
            cov_pattern=args.coverage,
            zero_pattern=args.zero,
            output_dir=args.output,
            logger=logger
        )
    
    if args.mode in ['quick', 'both']:
        run_quick_excel_report(
            bed_path=args.bed,
            panel_pattern=args.panels,
            output_dir=args.output,
            logger=logger
        )
    
    logger.info("\n" + "=" * 70)
    logger.info("ALL ANALYSES COMPLETE")
    logger.info("=" * 70)

if __name__ == "__main__":
    main()