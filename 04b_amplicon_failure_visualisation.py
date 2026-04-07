#!/usr/bin/env python3
"""
Visualise flagged DNA amplicon failures using cohort coverage tables and BAM depth.

For each flagged locus this script:
1. Identifies the worst cohort from the technical audit.
2. Selects one failing DNA sample (mean amplicon depth == 0) and one passing DNA sample.
3. Plots the cohort-wide amplicon-depth distribution.
4. Extracts per-base depth from the corresponding BAMs with samtools depth.
5. Saves one figure per locus plus a summary table and an IGV batch file.
"""

from __future__ import annotations

import argparse
import gzip
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
AUDIT_CSV = ROOT / 'analysis_results/04_technical_audit/GSDMB_Technical_Audit.csv'
SAMPLE_TABLES = ROOT / 'analysis_results/03_qc_visualisation/sample_level_coverage_tables'
OUTPUT_DIR = ROOT / 'analysis_results/04b_amplicon_failure_visualisation'
REFERENCE = ROOT / 'ref/hg38.fa'

FLAGGED_LABELS = [
    'rs921650',
    'rs12936409',
    'rs2305479',
    'chr17:39907450-39907625',
    'GSDMB@39905415',
    'GSDMB@39908947',
]

COHORT_LABELS = {
    'Breast control': 'Breast_Control',
    'Breast ctrl': 'Breast_Control',
    'Breast tumour': 'Breast_Tumour',
    'Endom ctrl': 'Endometrium_Control',
    'Endom tum': 'Endometrium_Tumour',
    'Endometrium ctrl': 'Endometrium_Control',
    'Endometrium tumour': 'Endometrium_Tumour',
}


@dataclass
class SamplePick:
    sample_key: str
    cohort_depth: float
    coverage_table: Path
    bam_path: Path


def normalise_name(text: str) -> str:
    text = Path(text).name.lower()
    for suffix in ('.coverage.tsv.gz', '.bam', '.bai'):
        if text.endswith(suffix):
            text = text[: -len(suffix)]
    text = text.replace('(', '_').replace(')', '_')
    text = re.sub(r'[^a-z0-9]+', '_', text)
    text = re.sub(r'_+', '_', text).strip('_')
    return text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', default=str(OUTPUT_DIR), help='Directory for figures and summary outputs')
    parser.add_argument('--padding', type=int, default=80, help='Flanking bases to include around each amplicon in BAM depth plots')
    return parser.parse_args()


def resolve_flagged_label(row: pd.Series) -> str:
    if row['annotation'] == 'GSDMB':
        return f"GSDMB@{int(row['start'])}"
    return row['annotation'] if str(row['annotation']).startswith('rs') else row['match_id']


def load_audit() -> pd.DataFrame:
    audit = pd.read_csv(AUDIT_CSV)
    audit['Flagged label'] = audit.apply(resolve_flagged_label, axis=1)
    return audit


def worst_cohort_from_row(row: pd.Series) -> str:
    fail_cols = [
        'Breast_Control_fail_%',
        'Breast_Tumour_fail_%',
        'Endometrium_Control_fail_%',
        'Endometrium_Tumour_fail_%',
    ]
    available = {col.replace('_fail_%', ''): float(row[col]) for col in fail_cols if col in row and pd.notna(row[col])}
    if not available:
        raise KeyError('No cohort-specific failure columns found in audit row')
    return max(available, key=available.get)


def build_bam_index() -> dict[str, Path]:
    bam_index: dict[str, Path] = {}
    for pattern in ('breast/*/dna/*.bam', 'endometrium/*/dna/*.bam'):
        for bam_path in ROOT.glob(pattern):
            bam_index[normalise_name(str(bam_path))] = bam_path
    return bam_index


def load_cohort_zero_counts(cohort: str) -> dict[str, int]:
    audit_path = ROOT / f'analysis_results/03_qc_visualisation/CoverageAudit_{cohort}.csv'
    if not audit_path.exists():
        return {}
    audit = pd.read_csv(audit_path)
    return {normalise_name(sample): int(zeros) for sample, zeros in zip(audit['sample'], audit['zeros'])}


def read_amplicon_depths(match_id: str, cohort: str, bam_index: dict[str, Path]) -> list[SamplePick]:
    rows: list[SamplePick] = []
    cohort_dir = SAMPLE_TABLES / cohort
    for coverage_file in sorted(cohort_dir.glob('*.coverage.tsv.gz')):
        depth = None
        with gzip.open(coverage_file, 'rt') as handle:
            next(handle)
            for line in handle:
                fields = line.rstrip('\n').split('\t')
                if fields[0] == match_id:
                    depth = float(fields[4])
                    break
        if depth is None:
            continue
        sample_key = normalise_name(str(coverage_file))
        bam_path = bam_index.get(sample_key)
        if bam_path is None:
            continue
        rows.append(SamplePick(sample_key=sample_key, cohort_depth=depth, coverage_table=coverage_file, bam_path=bam_path))
    return rows


def choose_examples(rows: list[SamplePick], zero_counts: dict[str, int]) -> tuple[SamplePick | None, SamplePick | None]:
    failing = [row for row in rows if row.cohort_depth == 0]
    failing.sort(key=lambda row: (zero_counts.get(row.sample_key, 9999), row.sample_key))
    fail = failing[0] if failing else None

    passing = [row for row in rows if row.cohort_depth > 0]
    passing.sort(key=lambda row: (zero_counts.get(row.sample_key, 9999), abs(row.cohort_depth - np.median([p.cohort_depth for p in passing]) if passing else 0)))
    good = passing[0] if passing else None
    return fail, good


def samtools_depth(bam_path: Path, region: str) -> pd.DataFrame:
    cmd = ['samtools', 'depth', '-aa', '-r', region, str(bam_path)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    records = []
    for line in result.stdout.splitlines():
        chrom, pos, depth = line.split('\t')
        records.append((chrom, int(pos), int(depth)))
    return pd.DataFrame(records, columns=['chrom', 'pos', 'depth'])


def per_base_metrics(depth_df: pd.DataFrame, start: int, end: int) -> dict[str, float]:
    target = depth_df[(depth_df['pos'] >= start) & (depth_df['pos'] <= end)]
    flank = depth_df[(depth_df['pos'] < start) | (depth_df['pos'] > end)]
    return {
        'target_mean_depth': float(target['depth'].mean()) if not target.empty else np.nan,
        'target_zero_frac': float((target['depth'] == 0).mean()) if not target.empty else np.nan,
        'flank_mean_depth': float(flank['depth'].mean()) if not flank.empty else np.nan,
        'flank_zero_frac': float((flank['depth'] == 0).mean()) if not flank.empty else np.nan,
    }


def infer_pattern(fail_metrics: dict[str, float], pass_metrics: dict[str, float]) -> str:
    fail_target = fail_metrics.get('target_mean_depth', np.nan)
    fail_flank = fail_metrics.get('flank_mean_depth', np.nan)
    pass_target = pass_metrics.get('target_mean_depth', np.nan)
    if np.isnan(fail_target) or np.isnan(pass_target):
        return 'insufficient depth data'
    if fail_target == 0 and fail_flank > 50:
        return 'focal dropout within target while nearby sequence still has reads'
    if fail_target == 0 and fail_flank <= 5:
        return 'broad local absence of reads in the failed sample'
    if fail_target < max(pass_target * 0.1, 20):
        return 'strong sample-specific depth collapse across this amplicon'
    return 'partial depth reduction rather than a full dropout'


def safe_slug(text: str) -> str:
    text = re.sub(r'[^A-Za-z0-9]+', '_', text).strip('_')
    return text or 'locus'


def rank_plot_order(rows: Iterable[SamplePick]) -> pd.DataFrame:
    df = pd.DataFrame([{'sample_key': row.sample_key, 'cohort_depth': row.cohort_depth} for row in rows]).sort_values('cohort_depth', ascending=True)
    df['rank'] = np.arange(1, len(df) + 1)
    return df


def make_figure(row: pd.Series, cohort: str, all_rows: list[SamplePick], fail_pick: SamplePick | None, pass_pick: SamplePick | None, fail_depth: pd.DataFrame | None, pass_depth: pd.DataFrame | None, fail_metrics: dict[str, float] | None, pass_metrics: dict[str, float] | None, output_path: Path, padding: int) -> None:
    start = int(row['start'])
    end = int(row['end'])
    region_start = start - padding
    region_end = end + padding
    rank_df = rank_plot_order(all_rows)

    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.5], hspace=0.28)

    ax1 = fig.add_subplot(gs[0])
    ax1.scatter(rank_df['rank'], rank_df['cohort_depth'], color='#9aa4b2', s=28)
    ax1.set_yscale('symlog', linthresh=1)
    ax1.set_ylabel('Mean amplicon depth')
    ax1.set_xlabel(f"{cohort.replace('_', ' ')} samples ranked by amplicon depth")
    ax1.set_title(f"{row['Flagged label']} | {cohort.replace('_', ' ')} | {row['Global_Failure_%']:.1f}% global fail", fontsize=14, weight='bold')
    if fail_pick:
        fx = int(rank_df.loc[rank_df['sample_key'] == fail_pick.sample_key, 'rank'].iloc[0])
        ax1.scatter([fx], [fail_pick.cohort_depth], color='#c0392b', s=65, label='Fail example')
        ax1.annotate(fail_pick.sample_key, (fx, max(fail_pick.cohort_depth, 0.15)), xytext=(8, 8), textcoords='offset points', fontsize=8, color='#c0392b')
    if pass_pick:
        px = int(rank_df.loc[rank_df['sample_key'] == pass_pick.sample_key, 'rank'].iloc[0])
        ax1.scatter([px], [pass_pick.cohort_depth], color='#1f77b4', s=65, label='Pass example')
        ax1.annotate(pass_pick.sample_key, (px, pass_pick.cohort_depth), xytext=(8, -12), textcoords='offset points', fontsize=8, color='#1f77b4')
    ax1.legend(loc='upper left', frameon=False)
    ax1.grid(alpha=0.2)

    ax2 = fig.add_subplot(gs[1])
    if fail_depth is not None and not fail_depth.empty:
        ax2.plot(fail_depth['pos'], fail_depth['depth'], color='#c0392b', lw=1.6, label='Fail sample')
    if pass_depth is not None and not pass_depth.empty:
        ax2.plot(pass_depth['pos'], pass_depth['depth'], color='#1f77b4', lw=1.6, label='Pass sample')
    ax2.axvspan(start, end, color='#f1c40f', alpha=0.25, label='Amplicon')
    ax2.set_xlim(region_start, region_end)
    ax2.set_ylabel('Per-base depth')
    ax2.set_xlabel(f'chr17 position ({padding} bp flank each side)')
    ax2.grid(alpha=0.2)
    ax2.legend(loc='upper right', frameon=False)

    info_lines = [
        f"Target: chr17:{start}-{end} ({int(row['length'])} bp)",
        f"Worst cohort: {cohort.replace('_', ' ')}",
        f"Hairpin risk: {int(row['Hairpin_Risk'])} | GC: {row['GC_pct']}%",
    ]
    if fail_metrics and pass_metrics:
        info_lines.append('Fail target mean {:.1f}x, flank {:.1f}x | Pass target mean {:.1f}x'.format(fail_metrics['target_mean_depth'], fail_metrics['flank_mean_depth'], pass_metrics['target_mean_depth']))
    fig.text(0.01, 0.01, '\n'.join(info_lines), fontsize=9, va='bottom')
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def build_igv_batch(records: list[dict[str, object]], output_dir: Path, padding: int) -> None:
    lines = ['new', f'genome {REFERENCE}', 'snapshotDirectory .']
    for record in records:
        fail_bam = record.get('fail_bam')
        pass_bam = record.get('pass_bam')
        if not fail_bam or not pass_bam:
            continue
        start = int(record['start']) - padding
        end = int(record['end']) + padding
        region = f'chr17:{start}-{end}'
        lines.append(f'load {fail_bam}')
        lines.append(f'load {pass_bam}')
        lines.append(f'goto {region}')
        lines.append('sort base')
        lines.append(f"snapshot {record['slug']}.png")
        lines.append('new')
        lines.append(f'genome {REFERENCE}')
    (output_dir / 'igv_batch_flagged_amplicons.txt').write_text('\n'.join(lines) + '\n')


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    audit = load_audit()
    audit = audit[audit['Flagged label'].isin(FLAGGED_LABELS)].copy()
    bam_index = build_bam_index()
    summary_records: list[dict[str, object]] = []

    for _, row in audit.iterrows():
        cohort = worst_cohort_from_row(row)
        all_rows = read_amplicon_depths(row['match_id'], cohort, bam_index)
        zero_counts = load_cohort_zero_counts(cohort)
        fail_pick, pass_pick = choose_examples(all_rows, zero_counts)

        fail_depth = pass_depth = None
        fail_metrics = pass_metrics = None
        region = f"chr17:{int(row['start']) - args.padding}-{int(row['end']) + args.padding}"

        if fail_pick:
            fail_depth = samtools_depth(fail_pick.bam_path, region)
            fail_metrics = per_base_metrics(fail_depth, int(row['start']), int(row['end']))
        if pass_pick:
            pass_depth = samtools_depth(pass_pick.bam_path, region)
            pass_metrics = per_base_metrics(pass_depth, int(row['start']), int(row['end']))

        slug = safe_slug(str(row['Flagged label']))
        make_figure(row=row, cohort=cohort, all_rows=all_rows, fail_pick=fail_pick, pass_pick=pass_pick, fail_depth=fail_depth, pass_depth=pass_depth, fail_metrics=fail_metrics, pass_metrics=pass_metrics, output_path=output_dir / f'{slug}.png', padding=args.padding)

        summary_records.append({
            'flagged_label': row['Flagged label'],
            'match_id': row['match_id'],
            'start': int(row['start']),
            'end': int(row['end']),
            'global_fail_pct': row['Global_Failure_%'],
            'global_mean_depth': row['Global_Mean_Depth'],
            'worst_cohort': cohort,
            'hairpin_risk': int(row['Hairpin_Risk']),
            'gc_pct': row['GC_pct'],
            'fail_sample': fail_pick.sample_key if fail_pick else '',
            'fail_depth_table': fail_pick.cohort_depth if fail_pick else np.nan,
            'fail_zero_amplicons_total': zero_counts.get(fail_pick.sample_key, np.nan) if fail_pick else np.nan,
            'fail_bam': str(fail_pick.bam_path) if fail_pick else '',
            'pass_sample': pass_pick.sample_key if pass_pick else '',
            'pass_depth_table': pass_pick.cohort_depth if pass_pick else np.nan,
            'pass_zero_amplicons_total': zero_counts.get(pass_pick.sample_key, np.nan) if pass_pick else np.nan,
            'pass_bam': str(pass_pick.bam_path) if pass_pick else '',
            'fail_target_mean_depth': fail_metrics['target_mean_depth'] if fail_metrics else np.nan,
            'fail_target_zero_frac': fail_metrics['target_zero_frac'] if fail_metrics else np.nan,
            'fail_flank_mean_depth': fail_metrics['flank_mean_depth'] if fail_metrics else np.nan,
            'pass_target_mean_depth': pass_metrics['target_mean_depth'] if pass_metrics else np.nan,
            'pass_target_zero_frac': pass_metrics['target_zero_frac'] if pass_metrics else np.nan,
            'pass_flank_mean_depth': pass_metrics['flank_mean_depth'] if pass_metrics else np.nan,
            'suggested_pattern': infer_pattern(fail_metrics or {}, pass_metrics or {}),
            'slug': slug,
        })

    summary_df = pd.DataFrame(summary_records).sort_values(['global_fail_pct', 'global_mean_depth'], ascending=[False, True])
    summary_df.to_csv(output_dir / 'flagged_amplicon_visualisation_summary.tsv', sep='\t', index=False)
    build_igv_batch(summary_records, output_dir, args.padding)


if __name__ == '__main__':
    main()
