import argparse
import os

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from collections import Counter
from pathlib import Path
import textwrap

ROOT = Path(__file__).resolve().parent
DOCS_DIR = ROOT / "docs"


def _iter_docs_search_dirs() -> list[Path]:
    seen: set[str] = set()
    candidates: list[Path] = []

    env_docs = os.environ.get("TFM_DOCS_DIR", "").strip()
    if env_docs:
        env_path = Path(env_docs)
        candidates.extend([env_path / "clinical variables-snps", env_path])

    candidates.extend([DOCS_DIR / "source_workbooks", DOCS_DIR / "derived_workbooks", DOCS_DIR])

    users_root = Path("/mnt/c/Users")
    if users_root.exists():
        for user_dir in sorted(users_root.iterdir()):
            if not user_dir.is_dir():
                continue
            downloads = user_dir / "Downloads"
            if downloads.exists():
                candidates.append(downloads)
            for onedrive_dir in sorted(user_dir.glob("OneDrive*")):
                docs_dir = onedrive_dir / "Documents" / "TFM" / "Docs"
                candidates.extend([docs_dir / "clinical variables-snps", docs_dir])

    existing: list[Path] = []
    for candidate in candidates:
        key = str(candidate)
        if key in seen or not candidate.exists():
            continue
        seen.add(key)
        existing.append(candidate)
    return existing


def _discover_doc_file(filename: str, fallback_dir: Path) -> Path:
    for search_dir in _iter_docs_search_dirs():
        candidate = search_dir / filename
        if candidate.exists():
            return candidate
    return fallback_dir / filename


DEFAULT_PHASED_PATH = ROOT / 'analysis_results/15_haplotype_phasing/phased_genotypes.tsv'
DEFAULT_META_PATH = ROOT / 'analysis_results/15_haplotype_phasing/sample_metadata.tsv'
DEFAULT_COLLAB_PATH = _discover_doc_file('SNPS PROYECTO MAMA_ENDOMETRIO(8).xlsx', DOCS_DIR / 'source_workbooks')
DEFAULT_OUTDIR = ROOT / 'analysis_results/19_1000g_haplotype_comparison/core_haplotype_comparison'

SNP_SETS = {
    '4SNP': [
        ('rs2290400', 39909987),
        ('rs1008723', 39910014),
        ('rs869402', 39911790),
        ('rs870829', 39912129),
    ],
    '5SNP': [
        ('rs2290400', 39909987),
        ('rs1008723', 39910014),
        ('rs869402', 39911790),
        ('rs870829', 39912129),
        ('rs7216389', 39913696),
    ],
}
GROUPS = {
    'Study controls': lambda m: m['Tissue'].eq('Healthy'),
    'Study tumours': lambda m: m['Tissue'].eq('Tumour'),
    'Breast tumour': lambda m: (m['Cohort'].eq('Breast') & m['Tissue'].eq('Tumour')),
    'Endometrial tumour': lambda m: (m['Cohort'].eq('Endometrium') & m['Tissue'].eq('Tumour')),
    'Study all': lambda m: pd.Series(True, index=m.index),
}

COLORS = {
    'Collaborator 1000G EUR': '#4C78A8',
    'Study controls': '#72B7B2',
    'Study tumours': '#E45756',
}


def parse_args():
    parser = argparse.ArgumentParser(
        description='Build compact core-haplotype comparisons between the study cohort and the collaborator workbook.'
    )
    parser.add_argument('--phased-path', default=str(DEFAULT_PHASED_PATH), help='Phased study genotype table from stage 15.')
    parser.add_argument('--meta-path', default=str(DEFAULT_META_PATH), help='Study sample metadata table from stage 15.')
    parser.add_argument('--collab-xlsx', default=str(DEFAULT_COLLAB_PATH), help='Collaborator workbook with the core_haplotype sheet.')
    parser.add_argument('--outdir', default=str(DEFAULT_OUTDIR), help='Output directory for stage-21 summaries and figures.')
    return parser.parse_args()


def parse_collab_core(path: Path):
    raw = pd.read_excel(path, sheet_name='core_haplotype')
    rows = []
    current = None
    for _, row in raw.iterrows():
        snps = row.get('SNPs')
        hap = row.get('hapl_eu')
        count = row.get('Count')
        freq = row.get('Freq')
        if isinstance(snps, str) and 'rs2290400' in snps:
            snp_list = [s.strip() for s in snps.split('\n') if s and str(s).startswith('rs')]
            current = '4SNP' if len(snp_list) == 4 else '5SNP' if len(snp_list) == 5 else None
        if current is None:
            continue
        if pd.isna(hap) or str(hap).strip() in {'', 'hapl_eu'}:
            continue
        try:
            count_val = int(float(count))
            freq_val = float(freq)
        except Exception:
            continue
        rows.append({
            'Set': current,
            'Dataset': 'Collaborator 1000G EUR',
            'Haplotype': str(hap).strip(),
            'Count': count_val,
            'Frequency': freq_val,
        })
    df = pd.DataFrame(rows)

    # rank() returns NaN for rows where Frequency is NaN.
    # Use pandas nullable integer 'Int64' (capital I) so NaN ranks are
    # preserved as pd.NA rather than crashing on .astype(int).
    n_nan = df['Frequency'].isna().sum()
    if n_nan > 0:
        print(f"  WARNING: {n_nan} collaborator rows have NaN Frequency — their ranks will be pd.NA")
    df['Rank'] = (
        df.groupby('Set')['Frequency']
        .rank(method='first', ascending=False)
        .astype('Int64')
    )
    return df.sort_values(['Set', 'Rank', 'Haplotype'])


def compute_study_core(phased_path: Path, meta_path: Path):
    phased = pd.read_csv(phased_path, sep='\t', dtype=str)
    phased['POS'] = phased['POS'].astype(int)
    meta = pd.read_csv(meta_path, sep='\t', dtype=str).set_index('Sample')
    sample_cols = list(phased.columns[5:])
    rows = []
    summary_rows = []
    for set_name, pairs in SNP_SETS.items():
        ordered_pos = [pos for _, pos in pairs]
        sub = phased[phased['POS'].isin(ordered_pos)].copy().set_index('POS').loc[ordered_pos].reset_index()
        for dataset, selector in GROUPS.items():
            samples = meta[selector(meta)].index.intersection(sample_cols)
            counts = Counter()
            total_chrom = 0
            complete_samples = 0
            for sample in samples:
                vals = sub[sample].tolist()
                if any(pd.isna(v) or v in {'.|.', './.', '.', 'nan'} for v in vals):
                    continue
                a = ''.join(v.split('|')[0] for v in vals)
                b = ''.join(v.split('|')[1] for v in vals)
                counts[a] += 1
                counts[b] += 1
                total_chrom += 2
                complete_samples += 1
            freqs = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
            for rank, (hap, cnt) in enumerate(freqs, start=1):
                rows.append({
                    'Set': set_name,
                    'Dataset': dataset,
                    'Haplotype': hap,
                    'Count': cnt,
                    'Frequency': cnt / total_chrom if total_chrom else 0.0,
                    'Rank': rank,
                })
            top5_mass = sum(cnt for _, cnt in freqs[:5]) / total_chrom if total_chrom else 0.0
            summary_rows.append({
                'Set': set_name,
                'Dataset': dataset,
                'Complete_Samples': complete_samples,
                'Chromosomes': total_chrom,
                'Observed_Haplotypes': len(freqs),
                'Top_Haplotype': freqs[0][0] if freqs else None,
                'Top_Haplotype_Freq': freqs[0][1] / total_chrom if freqs and total_chrom else 0.0,
                'Top5_Frequency_Mass': top5_mass,
            })
    return pd.DataFrame(rows), pd.DataFrame(summary_rows)


def describe_haplotype(set_name: str, haplotype: str) -> str:
    hap = str(haplotype).strip()
    pairs = SNP_SETS.get(set_name, [])
    if not pairs or len(hap) != len(pairs) or any(bit not in {'0', '1'} for bit in hap):
        return hap
    alt_idx = [str(i + 1) for i, bit in enumerate(hap) if bit == '1']
    if not alt_idx:
        return 'Ref'
    if len(alt_idx) == len(pairs):
        return 'ALT at all\nSNPs'
    joined = '+'.join(alt_idx)
    return f'Alt {joined}'


def snp_key_text(set_name: str) -> str:
    pairs = SNP_SETS[set_name]
    return '  '.join(f'SNP{i}={rsid}' for i, (rsid, _) in enumerate(pairs, start=1))


def plot_set(set_name: str, combined: pd.DataFrame, outpath: Path):
    plot_df = combined[(combined['Set'] == set_name) & (combined['Dataset'].isin(['Collaborator 1000G EUR', 'Study controls', 'Study tumours']))].copy()
    fig, axes = plt.subplots(1, 3, figsize=(18, 7.0), sharey=True)
    for ax, dataset in zip(axes, ['Collaborator 1000G EUR', 'Study controls', 'Study tumours']):
        sub = plot_df[plot_df['Dataset'] == dataset].sort_values('Frequency', ascending=False).head(8).copy()
        sub['Display_Label'] = sub['Haplotype'].map(lambda hap: describe_haplotype(set_name, hap))
        bars = ax.bar(range(len(sub)), sub['Frequency'], color=COLORS[dataset], alpha=0.9)
        ax.set_title(dataset, fontsize=13, weight='bold')
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(sub['Display_Label'], rotation=0, fontsize=9)
        ax.tick_params(axis='x', pad=10)
        ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
        ax.grid(axis='y', linestyle=':', alpha=0.4)
        ax.set_axisbelow(True)
        for bar, freq in zip(bars, sub['Frequency']):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008, f'{freq:.1%}', ha='center', va='bottom', fontsize=9)
    axes[0].set_ylabel('Haplotype frequency')
    fig.suptitle(f'Core {set_name} haplotypes: collaborator EUR vs study frequencies', fontsize=18, weight='bold', y=0.98)
    fig.text(0.5, 0.045, 'Ref = reference haplotype. Alt 1+4 means the alternate allele is present at SNP1 and SNP4.', ha='center', fontsize=10)
    fig.text(0.5, 0.015, snp_key_text(set_name), ha='center', fontsize=9)
    fig.tight_layout(rect=[0, 0.14, 1, 0.93])
    fig.savefig(outpath, dpi=220, bbox_inches='tight')
    plt.close(fig)


def build_rank_table(combined: pd.DataFrame, summary: pd.DataFrame, outpath: Path):
    rows = []
    for set_name in ['4SNP', '5SNP']:
        collab = combined[(combined['Set'] == set_name) & (combined['Dataset'] == 'Collaborator 1000G EUR')].sort_values('Frequency', ascending=False).head(6)
        controls = combined[(combined['Set'] == set_name) & (combined['Dataset'] == 'Study controls')].sort_values('Frequency', ascending=False).head(6)
        tumours = combined[(combined['Set'] == set_name) & (combined['Dataset'] == 'Study tumours')].sort_values('Frequency', ascending=False).head(6)
        rows.append({
            'Set': set_name,
            'Collaborator_haplotypes': int(summary[(summary['Set'] == set_name) & (summary['Dataset'] == 'Collaborator 1000G EUR')]['Observed_Haplotypes'].iloc[0]),
            'Study_all_haplotypes': int(summary[(summary['Set'] == set_name) & (summary['Dataset'] == 'Study all')]['Observed_Haplotypes'].iloc[0]),
            'Controls_top': '; '.join(f"{r.Haplotype} ({r.Frequency:.1%})" for r in controls.itertuples()),
            'Tumours_top': '; '.join(f"{r.Haplotype} ({r.Frequency:.1%})" for r in tumours.itertuples()),
            'Collaborator_top': '; '.join(f"{r.Haplotype} ({r.Frequency:.1%})" for r in collab.itertuples()),
            'Study_all_top5_mass': float(summary[(summary['Set'] == set_name) & (summary['Dataset'] == 'Study all')]['Top5_Frequency_Mass'].iloc[0]),
        })
    pd.DataFrame(rows).to_csv(outpath, sep='\t', index=False)


def main():
    args = parse_args()
    phased_path = Path(args.phased_path)
    meta_path = Path(args.meta_path)
    collab_path = Path(args.collab_xlsx)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not phased_path.exists():
        raise FileNotFoundError(f'Phased genotype table not found: {phased_path}')
    if not meta_path.exists():
        raise FileNotFoundError(f'Sample metadata table not found: {meta_path}')
    if not collab_path.exists():
        raise FileNotFoundError(f'Collaborator workbook not found: {collab_path}')

    collab = parse_collab_core(collab_path)
    study, study_summary = compute_study_core(phased_path, meta_path)
    collab_summary = collab.groupby(['Set', 'Dataset'], as_index=False).agg(
        Complete_Samples=('Count', lambda s: pd.NA),
        Chromosomes=('Count', 'sum'),
        Observed_Haplotypes=('Haplotype', 'nunique'),
        Top_Haplotype=('Haplotype', lambda s: collab.loc[s.index].sort_values('Frequency', ascending=False).iloc[0]['Haplotype']),
        Top_Haplotype_Freq=('Frequency', 'max'),
        Top5_Frequency_Mass=('Frequency', lambda s: s.sort_values(ascending=False).head(5).sum()),
    )
    combined = pd.concat([collab, study], ignore_index=True)
    summary = pd.concat([collab_summary, study_summary], ignore_index=True)
    combined.to_csv(outdir / 'core_haplotype_comparison_all.tsv', sep='\t', index=False)
    summary.to_csv(outdir / 'core_haplotype_comparison_summary.tsv', sep='\t', index=False)
    build_rank_table(combined, summary, outdir / 'core_haplotype_rank_summary.tsv')
    for set_name in SNP_SETS:
        plot_set(set_name, combined, outdir / f'Core_Haplotype_{set_name}_Collaborator_vs_Study.png')
    print('Wrote', outdir)

if __name__ == '__main__':
    main()