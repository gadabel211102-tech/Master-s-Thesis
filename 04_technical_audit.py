#!/usr/bin/env python3
"""
Amplicon technical audit combining sequence complexity metrics with coverage
failure statistics. Operates in two modes:

  comprehensive  — extracts sequence features from the reference genome and
                   merges them with per-sample coverage and zero-coverage data;
                   produces a ranked CSV of all amplicons
  quick          — generates an Excel workbook listing the 10 worst-performing
                   amplicons per cohort based on median depth; useful for a
                   rapid overview without needing the reference FASTA

Usage:
  python 04_technical_audit.py both           # recommended; runs both modes
  python 04_technical_audit.py comprehensive -v
  python 04_technical_audit.py quick
  python 04_technical_audit.py both --bed /path/to/targets.bed --fasta /path/to/hg38.fa
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


# =============================================================================
# Logging
# =============================================================================

def setup_logging(verbose: bool = False) -> logging.Logger:
    """Configure timestamped logging; DEBUG level when verbose is True."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%H:%M:%S'
    )
    return logging.getLogger(__name__)


# =============================================================================
# Shared utilities
# =============================================================================

def normalise_chr(c: str) -> str:
    """Ensure chromosome names carry the 'chr' prefix (e.g. '17' -> 'chr17')."""
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
    """Strip characters that are illegal in Excel sheet names and enforce the 31-character limit."""
    return re.sub(r'[\\/*?:\[\]]', '', name)[:max_length]


def load_bed_annotations(bed_path: str, logger: logging.Logger) -> pd.DataFrame:
    """
    Load the target BED file and return one row per amplicon with a
    human-readable annotation label and a coordinate-based match_id for joining.

    Annotation priority: rsID (col 5, must start with 'rs') > extra column (col 6)
    > genomic coordinates as fallback. Generic labels such as 'intron' or '.'
    are discarded as uninformative.
    """
    logger.info(f"Loading BED annotations: {bed_path}")
    try:
        bed_df = pd.read_csv(bed_path, sep='\t', comment='t',
                             skiprows=1, header=None)

        # Assign column names flexibly based on how many columns are present
        if len(bed_df.columns) >= 6:
            bed_df.columns = (['chr', 'start', 'end', 'numeric_id', 'gene_snp', 'extra'] +
                              [f'col{i}' for i in range(6, len(bed_df.columns))])
        elif len(bed_df.columns) >= 5:
            bed_df.columns = (['chr', 'start', 'end', 'numeric_id', 'gene_snp'] +
                              [f'col{i}' for i in range(5, len(bed_df.columns))])
        else:
            bed_df.columns = (['chr', 'start', 'end', 'numeric_id'] +
                              [f'col{i}' for i in range(4, len(bed_df.columns))])

        bed_df['chr'] = bed_df['chr'].astype(str).apply(normalise_chr)

        if 'gene_snp' in bed_df.columns:
            bed_df['annotation'] = bed_df['gene_snp'].astype(str)

            # Discard generic positional labels that provide no variant-level information
            non_informative = ['.', 'in', 'ex', 'intron', 'exon', 'intergenic',
                               'intronic', 'exonic', 'nan', 'none', '']
            bed_df['annotation'] = bed_df['annotation'].replace(non_informative, np.nan)

            # Retain only values that are rsIDs; gene names without an rsID are
            # less specific and may be ambiguous across multiple amplicons
            bed_df['annotation'] = bed_df['annotation'].apply(
                lambda x: x if pd.notna(x) and str(x).lower().startswith('rs') else np.nan
            )

            # Use the extra column (col 6) if rsID is absent
            if 'extra' in bed_df.columns:
                bed_df['annotation'] = bed_df['annotation'].fillna(bed_df['extra'])

            # Final fallback: genomic coordinates as a unique identifier
            bed_df['annotation'] = bed_df['annotation'].fillna(
                bed_df['chr'] + ':' + bed_df['start'].astype(str) + '-' + bed_df['end'].astype(str)
            )
        else:
            bed_df['annotation'] = (bed_df['chr'] + ':' +
                                    bed_df['start'].astype(str) + '-' + bed_df['end'].astype(str))

        # match_id is used as the join key across all data sources
        bed_df['match_id'] = (bed_df['chr'] + ":" +
                              bed_df['start'].astype(str) + "-" + bed_df['end'].astype(str))
        bed_df['length'] = bed_df['end'] - bed_df['start']

        logger.info(f"Loaded {len(bed_df)} amplicons from BED")
        return bed_df[['chr', 'start', 'end', 'annotation', 'match_id', 'length', 'numeric_id']]

    except Exception as e:
        logger.error(f"Failed to load BED file: {e}")
        raise


# =============================================================================
# Sequence complexity analysis
# =============================================================================

def analyse_sequence_complexity(seq: str) -> Tuple[float, float, float, float, int, int, int]:
    """
    Compute sequence features that predict technical failure on Ion Torrent platforms.

    Ion Torrent sequencing is particularly sensitive to:
      - Extreme GC content (poor polymerase processivity)
      - Long homopolymer runs (signal overflow causes indel errors)
      - Short tandem repeats / SSRs (slippage during amplification)
      - Hairpin-forming sequences (secondary structure blocks sequencing)

    Returns a tuple of:
      GC_pct         — overall GC content (%)
      GC_skew        — (G-C)/(G+C); asymmetry between strands
      Max_Window_GC  — highest GC% in any 20 bp sliding window
      Min_Window_GC  — lowest GC% in any 20 bp sliding window
      Max_Homopolymer — length of the longest single-base run
      Max_SSR_Run    — length of the longest di/trinucleotide repeat
      Hairpin_Risk   — 1 if a 6 bp stem with complementary downstream sequence
                       is found, 0 otherwise
    """
    if not seq:
        return tuple([np.nan] * 7)
    seq = str(seq).upper()
    L   = len(seq)
    if L == 0:
        return tuple([np.nan] * 7)

    g, c   = seq.count('G'), seq.count('C')
    gc_pct  = (g + c) / L * 100
    gc_skew = (g - c) / (g + c) if (g + c) > 0 else 0

    # Sliding window GC: detects local extremes that the overall GC% would miss
    window_size = min(20, L)
    if L >= window_size:
        windows      = [seq[i:i+window_size] for i in range(L - window_size + 1)]
        window_gcs   = [(w.count('G') + w.count('C')) / window_size * 100 for w in windows]
        max_window_gc = max(window_gcs)
        min_window_gc = min(window_gcs)
    else:
        max_window_gc = min_window_gc = gc_pct

    # Homopolymer detection: runs of a single nucleotide (e.g. AAAAAAA)
    homopolymer_runs = re.findall(r'([A-Z])\1*', seq)
    max_homo = max((len(r) for r in homopolymer_runs), default=0)

    # SSR detection: di- or trinucleotide repeats (e.g. ACACAC, ATGATGATG)
    multi_repeats = re.findall(r'([A-Z]{2,3})\1+', seq)
    max_ssr_run   = max((len(r) for r in multi_repeats), default=0)

    # Hairpin risk: search for a 6 bp stem whose reverse complement appears
    # within the same sequence downstream — a simple but sensitive heuristic
    hairpin_risk = 0
    if L >= 12:
        for i in range(L - 11):
            stem       = seq[i:i+6]
            downstream = seq[i+6:]
            if str(Seq(stem).reverse_complement()) in downstream:
                hairpin_risk = 1
                break

    return (round(gc_pct, 2), round(gc_skew, 3), round(max_window_gc, 1),
            round(min_window_gc, 1), max_homo, max_ssr_run, hairpin_risk)


def extract_sequence_features(bed: pd.DataFrame, fasta_path: str,
                               logger: logging.Logger) -> pd.DataFrame:
    """
    Extract the amplicon sequence for each BED region from the reference FASTA
    and compute complexity metrics. Rows that fail (e.g. chromosome not in FASTA)
    are retained with NaN metrics rather than dropped.
    """
    if not os.path.exists(fasta_path):
        logger.warning(f"FASTA not found: {fasta_path} — skipping sequence analysis")
        return bed

    logger.info("Extracting sequence complexity metrics...")
    # Load the entire genome into memory as a dictionary for fast random access
    genome  = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    results = []
    failed  = 0

    for _, row in bed.iterrows():
        try:
            sequence = str(genome[row['chr']].seq[row['start']:row['end']])
            metrics  = analyse_sequence_complexity(sequence)
            results.append({
                'match_id':        row['match_id'],
                'GC_pct':          metrics[0],
                'GC_skew':         metrics[1],
                'Max_Window_GC':   metrics[2],
                'Min_Window_GC':   metrics[3],
                'Max_Homopolymer': metrics[4],
                'Max_SSR_Run':     metrics[5],
                'Hairpin_Risk':    metrics[6]
            })
        except Exception as e:
            logger.debug(f"Sequence extraction failed for {row['match_id']}: {e}")
            results.append({'match_id': row['match_id']})   # NaN metrics retained
            failed += 1

    if failed:
        logger.warning(f"Sequence extraction failed for {failed}/{len(bed)} amplicons")

    return bed.merge(pd.DataFrame(results), on='match_id', how='left')


# =============================================================================
# Coverage data processing
# =============================================================================

def process_individual_coverage(cov_pattern: str,
                                 logger: logging.Logger) -> Tuple[pd.DataFrame, List[Dict]]:
    """
    Load all per-sample coverage TSV files matching the glob pattern and
    concatenate them into a master DataFrame.

    Each file contains one row per amplicon with a 'depth' column. A binary
    'is_failed' flag (1 = zero depth) is added here so that failure rates can
    be computed later as simple column means.

    Also returns a list of per-sample summary statistics for export.
    """
    all_cov_files = glob.glob(cov_pattern, recursive=True)
    if not all_cov_files:
        logger.warning(f"No coverage files found matching: {cov_pattern}")
        return pd.DataFrame(), []

    logger.info(f"Processing {len(all_cov_files)} individual coverage files...")
    cov_rows, sample_stats = [], []
    failed = 0

    for idx, f in enumerate(all_cov_files, 1):
        sample_name = os.path.basename(f).replace('.coverage.tsv.gz', '')
        cohort      = os.path.basename(os.path.dirname(f))

        if idx % 50 == 0:
            logger.debug(f"  Processed {idx}/{len(all_cov_files)}...")

        try:
            df = pd.read_csv(f, sep='\t', compression='gzip')
            if 'depth' not in df.columns or 'id' not in df.columns:
                logger.warning(f"Missing required columns in {f}"); failed += 1; continue

            sample_stats.append({
                'Sample':                  sample_name,
                'Cohort':                  cohort,
                'Avg_Depth':               df['depth'].mean(),
                'Median_Depth':            df['depth'].median(),
                'Zero_Coverage_Amplicons': (df['depth'] == 0).sum()
            })

            df['is_failed'] = (df['depth'] == 0).astype(int)   # 1 if amplicon has zero depth
            df['cohort']    = cohort
            cov_rows.append(df[['id', 'is_failed', 'depth', 'cohort']])

        except Exception as e:
            logger.warning(f"Failed to process {f}: {e}"); failed += 1

    if failed:
        logger.warning(f"{failed}/{len(all_cov_files)} files could not be processed")

    return (pd.concat(cov_rows, ignore_index=True) if cov_rows else pd.DataFrame()), sample_stats


def process_panel_coverage(panel_pattern: str,
                            logger: logging.Logger) -> List[Tuple[str, pd.DataFrame]]:
    """
    Load the per-cohort FullPanel CSV files produced by script 03 and return
    them as a list of (cohort_name, DataFrame) tuples for use in the quick
    Excel report mode.
    """
    panel_files = glob.glob(panel_pattern, recursive=True)
    if not panel_files:
        logger.warning(f"No panel files found matching: {panel_pattern}")
        return []

    logger.info(f"Loading {len(panel_files)} panel coverage files...")
    cohort_data = []

    for f in panel_files:
        cohort_name = (os.path.basename(f)
                       .replace("GSDMB_", "")
                       .replace("_FullPanel.csv", ""))
        try:
            df = pd.read_csv(f)
            if 'chr' in df.columns:
                df['chr'] = df['chr'].astype(str).apply(normalise_chr)
            for col in ['start', 'end']:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
            cohort_data.append((cohort_name, df))
            logger.debug(f"  {cohort_name}: {len(df)} rows")
        except Exception as e:
            logger.warning(f"Failed to load {f}: {e}")

    return cohort_data


# =============================================================================
# Zero coverage processing
# =============================================================================

def process_zero_coverage(zero_pattern: str,
                           logger: logging.Logger) -> Optional[pd.DataFrame]:
    """
    Merge zero-coverage TSV files from 02b_more-qc.sh across all cohorts.

    Each file contains four groups of per-sample columns:
      <sample>_zero_bases / _pct_zero  — bases and percentage at 0x depth
      <sample>_low_bases  / _pct_low   — bases and percentage below 200x depth

    For each cohort, these are reduced to six per-amplicon summary columns:
      <cohort>_mean_zero_bases / _mean_pct_zero   — average across samples
      <cohort>_mean_low_bases  / _mean_pct_low    — average across samples
      <cohort>_n_samples_zero  — number of samples with any zero-coverage bases
      <cohort>_n_samples_low   — number of samples with any sub-200x bases

    All cohort summaries are outer-joined on match_id so that amplicons present
    in only a subset of cohorts are still represented (with NaN for absent cohorts).
    """
    all_zero_files = glob.glob(zero_pattern, recursive=True)
    if not all_zero_files:
        logger.warning(f"No zero-coverage files found matching: {zero_pattern}")
        return None

    logger.info(f"Merging {len(all_zero_files)} zero-coverage reports...")
    merged_zero = None
    failed      = 0
    FIXED_COLS  = {'chr', 'start', 'end', 'annotation', 'amplicon_length_bp'}

    for z_file in all_zero_files:
        cohort_id = os.path.basename(z_file).replace('zero_cov_', '').replace('.tsv', '')
        try:
            z_df = pd.read_csv(z_file, sep='\t')

            # Minimum of 9 columns = 5 fixed + at least one sample's 4 metrics
            if len(z_df.columns) < 9:
                logger.warning(f"Too few columns in {z_file} — may be old format"); failed += 1; continue

            z_df['chr']      = z_df['chr'].astype(str).apply(normalise_chr)
            z_df['match_id'] = (z_df['chr'] + ":" +
                                z_df['start'].astype(str) + "-" + z_df['end'].astype(str))

            # Identify per-sample column groups by suffix
            data_cols      = set(z_df.columns) - FIXED_COLS - {'match_id'}
            zero_base_cols = [c for c in data_cols if c.endswith('_zero_bases')]
            pct_zero_cols  = [c for c in data_cols if c.endswith('_pct_zero')]
            low_base_cols  = [c for c in data_cols if c.endswith('_low_bases')]
            pct_low_cols   = [c for c in data_cols if c.endswith('_pct_low')]

            if not zero_base_cols:
                logger.warning(f"No '_zero_bases' columns in {z_file} — skipping"); failed += 1; continue

            logger.debug(f"  {cohort_id}: {len(zero_base_cols)} samples")

            summary = z_df[['match_id']].copy()
            summary[f'{cohort_id}_mean_zero_bases'] = z_df[zero_base_cols].mean(axis=1).round(2)
            summary[f'{cohort_id}_mean_pct_zero']   = z_df[pct_zero_cols].mean(axis=1).round(2)
            summary[f'{cohort_id}_mean_low_bases']  = z_df[low_base_cols].mean(axis=1).round(2)
            summary[f'{cohort_id}_mean_pct_low']    = z_df[pct_low_cols].mean(axis=1).round(2)
            # Count rather than average: tells us whether failures are isolated or widespread
            summary[f'{cohort_id}_n_samples_zero']  = (z_df[zero_base_cols] > 0).sum(axis=1)
            summary[f'{cohort_id}_n_samples_low']   = (z_df[low_base_cols]  > 0).sum(axis=1)

            merged_zero = summary if merged_zero is None else pd.merge(
                merged_zero, summary, on='match_id', how='outer'
            )

        except Exception as e:
            logger.warning(f"Failed to process {z_file}: {e}"); failed += 1

    if failed:
        logger.warning(f"{failed}/{len(all_zero_files)} zero-coverage files could not be processed")

    return merged_zero


# =============================================================================
# Comprehensive audit mode
# =============================================================================

def run_comprehensive_audit(
    bed_path: str,
    fasta_path: str,
    cov_pattern: str,
    zero_pattern: str,
    output_dir: str,
    logger: logging.Logger
) -> pd.DataFrame:
    """
    Build a per-amplicon audit table combining sequence complexity features,
    per-cohort coverage failure rates, and zero/low-coverage summaries.
    The final CSV is sorted by Global_Failure_% (descending) so the most
    problematic amplicons appear at the top.
    """
    logger.info("=" * 70)
    logger.info("COMPREHENSIVE TECHNICAL AUDIT")
    logger.info("=" * 70)

    bed = load_bed_annotations(bed_path, logger)
    bed = extract_sequence_features(bed, fasta_path, logger)

    master_cov, sample_stats = process_individual_coverage(cov_pattern, logger)

    if sample_stats:
        sample_df = pd.DataFrame(sample_stats)
        sample_output = os.path.join(output_dir, "Sample_Average_Coverage.csv")
        sample_df.to_csv(sample_output, index=False)
        logger.info(f"Sample statistics saved: {sample_output}")

    if not master_cov.empty:
        logger.info("Computing per-amplicon coverage statistics...")

        # Global failure rate: fraction of samples in which an amplicon has zero depth
        global_fail   = master_cov.groupby('id')['is_failed'].mean() * 100
        global_depth  = master_cov.groupby('id')['depth'].mean()
        global_median = master_cov.groupby('id')['depth'].median()

        # Per-cohort failure rate and mean depth — unstack pivots cohort into columns
        cohort_fail  = (master_cov.groupby(['id', 'cohort'])['is_failed']
                        .mean() * 100).unstack().add_suffix('_fail_%')
        cohort_depth = (master_cov.groupby(['id', 'cohort'])['depth']
                        .mean()).unstack().add_suffix('_avg_depth')

        bed['Global_Failure_%']   = bed['match_id'].map(global_fail)
        bed['Global_Mean_Depth']  = bed['match_id'].map(global_depth)
        bed['Global_Median_Depth'] = bed['match_id'].map(global_median)

        bed = bed.merge(cohort_fail,  left_on='match_id', right_index=True, how='left')
        bed = bed.merge(cohort_depth, left_on='match_id', right_index=True, how='left')

    # Merge zero-coverage summaries from 02b (NaN-filled where absent)
    merged_zero = process_zero_coverage(zero_pattern, logger)
    if merged_zero is not None:
        bed = bed.merge(merged_zero, on='match_id', how='left').fillna(0)
        zero_cols_added = [c for c in merged_zero.columns if c != 'match_id']
        logger.info(f"Zero/low-coverage metrics added: {len(zero_cols_added)} columns")

    output_file = os.path.join(output_dir, "TFM_GSDMB_Final_Audit.csv")
    bed.sort_values(by='Global_Failure_%', ascending=False, na_position='last').to_csv(
        output_file, index=False
    )
    logger.info(f"Comprehensive audit saved: {output_file}")

    if 'Global_Failure_%' in bed.columns:
        logger.info("\nTop 10 worst amplicons:")
        worst_10 = bed.sort_values(
            by=['Global_Failure_%', 'Global_Mean_Depth'],
            ascending=[False, True],
            na_position='last'
        ).head(10)
        print(worst_10[['numeric_id', 'annotation', 'length',
                         'Global_Failure_%', 'Global_Mean_Depth']].to_string(index=False))

    return bed


# =============================================================================
# Quick Excel report mode
# =============================================================================

def run_quick_excel_report(
    bed_path: str,
    panel_pattern: str,
    output_dir: str,
    logger: logging.Logger
) -> None:
    """
    Generate a one-sheet-per-cohort Excel workbook listing the 10 amplicons
    with the lowest median depth. Intended as a rapid review tool that does
    not require the reference FASTA.
    """
    logger.info("=" * 70)
    logger.info("QUICK EXCEL REPORT — WORST 10 PER COHORT")
    logger.info("=" * 70)

    annotations = load_bed_annotations(bed_path, logger)
    cohort_data = process_panel_coverage(panel_pattern, logger)

    if not cohort_data:
        logger.error("No panel data found — cannot generate Excel report"); return

    output_file = os.path.join(output_dir, "GSDMB_Worst_Amplicons_PerCohort.xlsx")

    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        used_sheets = set()

        for cohort_name, df in cohort_data:
            df_merged = pd.merge(
                df,
                annotations[['chr', 'start', 'end', 'annotation']],
                on=['chr', 'start', 'end'],
                how='left',
                suffixes=('_data', '_bed')
            )

            # Prefer the BED annotation over any annotation column already in the panel CSV
            if 'annotation_bed' in df_merged.columns:
                df_merged['annotation'] = df_merged['annotation_bed'].fillna(
                    df_merged.get('annotation_data', '')
                )

            metadata_cols = {'chr', 'start', 'end', 'id', 'annotation',
                             'annotation_data', 'annotation_bed', 'match_id'}
            sample_cols   = [c for c in df_merged.columns
                             if c not in metadata_cols and
                             pd.api.types.is_numeric_dtype(df_merged[c])]

            if not sample_cols:
                logger.warning(f"No numeric depth columns for {cohort_name}"); continue

            # Rank by median depth across samples — median is preferred over mean
            # here because it is robust to the occasional extreme-depth outlier amplicon
            df_merged['median_depth'] = df_merged[sample_cols].median(axis=1)
            worst_10 = df_merged.sort_values('median_depth', ascending=True).head(10)

            # Deduplicate sheet names to satisfy Excel's uniqueness requirement
            base_sheet, i = clean_sheet_name(cohort_name), 1
            sheet_name    = base_sheet
            while sheet_name in used_sheets:
                sheet_name = clean_sheet_name(f"{base_sheet}_{i}"); i += 1
            used_sheets.add(sheet_name)

            worst_10[['chr', 'start', 'end', 'annotation', 'median_depth']].to_excel(
                writer, sheet_name=sheet_name, index=False
            )
            logger.info(f"  Sheet added: {sheet_name}")

    logger.info(f"Excel report saved: {output_file}")


# =============================================================================
# Entry point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Amplicon technical audit — comprehensive sequence analysis '
                    'and/or quick Excel failure summary.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:
  comprehensive   Full audit with sequence complexity (requires reference FASTA)
  quick           Excel report — worst 10 amplicons per cohort (no FASTA needed)
  both            Run both modes (recommended)

Examples:
  python 04_technical_audit.py both
  python 04_technical_audit.py comprehensive -v
  python 04_technical_audit.py both \\
    --bed /path/to/targets.bed \\
    --fasta /path/to/hg38.fa \\
    --coverage "/data/**/coverage/*.tsv.gz" \\
    --panels "/data/GSDMB_*_FullPanel.csv"
        """
    )
    parser.add_argument('mode', choices=['comprehensive', 'quick', 'both'],
                        help='Analysis mode')
    parser.add_argument('--bed',
                        default="/home/gadeaalonsoj/tfm/IAD255368_167_Submitted.bed",
                        help='Path to target BED file')
    parser.add_argument('--fasta',
                        default="/home/gadeaalonsoj/tfm/ref_alt/hg38_alt.fa",
                        help='Reference genome FASTA (comprehensive mode only)')
    parser.add_argument('--coverage',
                        default="/home/gadeaalonsoj/tfm/gsdmb_final_results/individual_sample_coverage/**/*.coverage.tsv.gz",
                        help='Glob pattern for per-sample coverage files')
    parser.add_argument('--zero',
                        default="/home/gadeaalonsoj/tfm/**/zero_cov_fast/zero_cov_*.tsv",
                        help='Glob pattern for zero-coverage reports from 02b')
    parser.add_argument('--panels',
                        default="/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_*_FullPanel.csv",
                        help='Glob pattern for FullPanel CSV files (quick mode only)')
    parser.add_argument('--output', '-o', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable DEBUG-level logging')

    args   = parser.parse_args()
    logger = setup_logging(args.verbose)
    os.makedirs(args.output, exist_ok=True)

    if args.mode in ['comprehensive', 'both']:
        run_comprehensive_audit(
            bed_path    = args.bed,
            fasta_path  = args.fasta,
            cov_pattern = args.coverage,
            zero_pattern= args.zero,
            output_dir  = args.output,
            logger      = logger
        )

    if args.mode in ['quick', 'both']:
        run_quick_excel_report(
            bed_path      = args.bed,
            panel_pattern = args.panels,
            output_dir    = args.output,
            logger        = logger
        )

    logger.info("\n" + "=" * 70)
    logger.info("All analyses complete.")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
