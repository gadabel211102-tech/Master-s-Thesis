#!/usr/bin/env python3
from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path('/home/gadeaalonsoj/tfm')
SUMMARY = ROOT / 'analysis_results/04b_amplicon_failure_visualisation/flagged_amplicon_visualisation_summary.tsv'
OUT = ROOT / 'analysis_results/04b_amplicon_failure_visualisation/Worst_Flagged_Amplicons_Overview.png'
PADDING = 80


def samtools_depth(bam_path: str, region: str) -> pd.DataFrame:
    cmd = ['samtools', 'depth', '-aa', '-r', region, bam_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    rows = []
    for line in result.stdout.splitlines():
        chrom, pos, depth = line.split('\t')
        rows.append((chrom, int(pos), int(depth)))
    return pd.DataFrame(rows, columns=['chrom', 'pos', 'depth'])


def pretty_label(label: str) -> str:
    return label.replace('chr17:', 'chr17 ').replace('-', ' to ')


def main() -> None:
    df = pd.read_csv(SUMMARY, sep='\t').sort_values(['global_fail_pct', 'global_mean_depth'], ascending=[False, True])
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    axes = axes.flatten()

    for ax, (_, row) in zip(axes, df.iterrows()):
        start = int(row['start'])
        end = int(row['end'])
        region = f'chr17:{start-PADDING}-{end+PADDING}'
        fail_depth = samtools_depth(row['fail_bam'], region)
        pass_depth = samtools_depth(row['pass_bam'], region)

        ax.plot(fail_depth['pos'], fail_depth['depth'], color='#c0392b', lw=1.6, label='Fail')
        ax.plot(pass_depth['pos'], pass_depth['depth'], color='#1f77b4', lw=1.6, label='Pass')
        ax.axvspan(start, end, color='#f1c40f', alpha=0.25)
        ax.set_xlim(start-PADDING, end+PADDING)
        ax.set_title(f"{row['flagged_label']} | {row['global_fail_pct']:.1f}% fail", fontsize=12, weight='bold')
        ax.text(0.02, 0.96, f"{row['worst_cohort'].replace('_', ' ')}\nfailed sample zero amplicons: {int(row['fail_zero_amplicons_total'])}", transform=ax.transAxes, va='top', fontsize=9, bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#d0d7de', alpha=0.9))
        ax.grid(alpha=0.18)
        ax.set_xlabel('chr17 position')
        ax.set_ylabel('Per-base depth')

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.01))
    fig.suptitle('Six worst flagged amplicons: representative failed vs passed DNA BAM depth', fontsize=18, weight='bold')
    fig.savefig(OUT, dpi=220, bbox_inches='tight')


if __name__ == '__main__':
    main()
