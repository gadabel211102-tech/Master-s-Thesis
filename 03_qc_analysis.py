
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# ==============================================================================
# 1. VALIDATED QC THRESHOLDS & CONFIGURATION
# ==============================================================================

QC_LIMITS = {
    'total_reads': 120000,
    'mapped_pct': 90.0,
    'on_target_pct': 75.0,
    'mean_cov': 200.0,
    'uniformity_100x': 95.0,
    'worst_amplicon_floor': 50.0
}

PATHS = {
    "out_dir": "./gsdmb_final_results",
    "samples_dir": "./gsdmb_final_results/individual_sample_coverage",
    "plots_dir": "./gsdmb_final_results/plots",
    "raw_root": "./"
}

COHORT_MAP = {
    'Breast Tumour': 'breast/tumour',
    'Endometrium Tumour': 'endometrium/tumour',
    'Endometrium Control': 'endometrium/normal'
}

PALETTE = {
    'Breast Tumour': '#3498db',
    'Endometrium Tumour': '#9b59b6',
    'Endometrium Control': '#2ecc71'
}

STATUS_PALETTE = {'Pass': '#2ecc71', 'Fail': '#e74c3c'}

# Fail loudly if a sample looks empty
STRICT_MOSDEPTH = True
MIN_NUMERIC_FRACTION = 0.50   # warn if <50% of rows are numeric for a sample

# ==============================================================================
# Utilities
# ==============================================================================

def normalise_sample(s: str) -> str:
    s = str(s).strip()
    s = s.replace(" ", "_")
    s = s.replace(".regions", "").replace(".bed", "").replace(".gz", "")
    return s

def normalise_chr(c: str) -> str:
    """Normalise chromosome labels to 'chr*'. Assumes hg-style naming."""
    c = str(c)
    if c.startswith("chr"):
        return c
    # Common numerics or single letters
    if c in {"X", "Y", "M", "MT"}:
        return "chrM" if c in {"M", "MT"} else f"chr{c}"
    # Numeric contigs
    try:
        _n = int(c)
        return f"chr{_n}"
    except ValueError:
        # Leave as-is for alt contigs, but ensure it starts with 'chr'
        return c if c.startswith("chr") else f"chr{c}"

def pick_depth_column(df: pd.DataFrame) -> int:
    """
    Determine which column is the depth column.
    Heuristic:
      - If 4 columns, the 4th is depth (mosdepth regions standard).
      - If >4, try last column first; fall back to any column that is mostly numeric.
      - If 3 columns, there is no depth -> return -1.
    """
    ncols = df.shape[1]
    if ncols == 4:
        return 3
    if ncols < 4:
        return -1
    # Try last col
    last = df.columns[-1]
    frac_numeric = pd.to_numeric(df[last], errors='coerce').notna().mean()
    if frac_numeric > 0.8:
        return ncols - 1
    # Otherwise, search all columns after the first three
    candidates = range(3, ncols)
    best_col = -1
    best_frac = -1.0
    for i in candidates:
        frac = pd.to_numeric(df.iloc[:, i], errors='coerce').notna().mean()
        if frac > best_frac:
            best_frac = frac
            best_col = i
    return best_col if best_frac >= 0.5 else -1

def load_one_regions(file_path: str, sample_name: str) -> pd.DataFrame:
    """
    Load a mosdepth .regions.bed.gz robustly and return a dataframe:
        index = id (chr:start-end)
        columns = ['chr', 'start', 'end', sample_name]  (depth numeric)
    Raises a clear error if depth looks invalid and STRICT_MOSDEPTH is True.
    """
    raw = pd.read_csv(file_path, sep='\t', header=None, compression='gzip', dtype=str)
    if raw.empty:
        msg = f"[ERROR] {file_path} is empty."
        if STRICT_MOSDEPTH: raise ValueError(msg)
        print(msg)
        return pd.DataFrame(columns=['chr', 'start', 'end', sample_name]).set_index([])

    # Ensure we have at least 3 coordinate columns
    if raw.shape[1] < 3:
        raise ValueError(f"[ERROR] {file_path} has <3 columns; cannot parse coordinates.")

    # Coerce coord columns
    raw.iloc[:, 0] = raw.iloc[:, 0].map(normalise_chr)
    # Start/end as integers where possible
    raw.iloc[:, 1] = pd.to_numeric(raw.iloc[:, 1], errors='coerce')
    raw.iloc[:, 2] = pd.to_numeric(raw.iloc[:, 2], errors='coerce')

    # Find depth column
    depth_col = pick_depth_column(raw)
    if depth_col == -1:
        msg = (f"[ERROR] {file_path} has no usable depth column; "
               f"detected {raw.shape[1]} columns. First rows:\n{raw.head(3)}")
        if STRICT_MOSDEPTH: raise ValueError(msg)
        print(msg)
        depth = pd.Series(np.nan, index=raw.index)
    else:
        depth = pd.to_numeric(raw.iloc[:, depth_col], errors='coerce')

    df = pd.DataFrame({
        'chr': raw.iloc[:, 0].astype(str),
        'start': raw.iloc[:, 1].astype('Int64'),
        'end': raw.iloc[:, 2].astype('Int64'),
        sample_name: depth
    })

    # Build robust ID
    df['id'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    df = df.drop_duplicates('id').set_index('id')

    # Quick per-sample audit
    frac_numeric = df[sample_name].notna().mean()
    if frac_numeric < MIN_NUMERIC_FRACTION:
        msg = (f"[WARN] {sample_name}: only {frac_numeric:.1%} of amplicons have numeric depth "
               f"(file {os.path.basename(file_path)}). "
               f"This will appear blank or nearly blank in heatmaps.")
        print(msg)
        if STRICT_MOSDEPTH:
            # Make it a hard error so you catch it immediately
            raise ValueError(msg)

    return df

def join_by_id(sample_frames: list) -> pd.DataFrame:
    """
    Outer-join all per-sample frames on 'id' and keep coord columns from the first.
    """
    base = None
    for one in sample_frames:
        # 'one' is indexed by 'id'; last col is the sample depth
        if base is None:
            base = one.copy()
        else:
            # Identify the sample column to join in
            sample_col = [c for c in one.columns if c not in ['chr', 'start', 'end']]
            base = base.join(one[sample_col], how='outer')
    return base

# ==============================================================================
# 2. PIPELINE
# ==============================================================================

def run_validated_pipeline():

    # Create directories
    for p in PATHS.values():
        os.makedirs(p, exist_ok=True)

    print("[*] Phase 1: Consolidating Panel Data and Harmonising Coordinates...")

    master_qc_list = []
    cohort_dfs = {}

    # Load QC + coverage per cohort
    for cohort, subpath in COHORT_MAP.items():
        full_subpath = os.path.join(PATHS['raw_root'], subpath)
        if not os.path.exists(full_subpath):
            continue

        # QC summaries
        qc_file = os.path.join(full_subpath, "dna_qc/qc_summary.tsv")
        if os.path.exists(qc_file):
            temp_df = pd.read_csv(qc_file, sep="\t").assign(cohort=cohort)
            temp_df["status"] = temp_df["status"].str.capitalize()
            temp_df["sample_norm"] = temp_df["sample"].apply(normalise_sample)
            master_qc_list.append(temp_df)

        # Coverage files
        bed_pattern = os.path.join(full_subpath, "**/*.regions.bed.gz")
        files = [
            f for f in glob.glob(bed_pattern, recursive=True)
            if f.endswith(".regions.bed.gz") and ".bai" not in f.lower()
        ]
        if not files:
            continue

        sample_frames = []
        debug_rows = []
        for f in files:
            s_name = normalise_sample(os.path.basename(f).replace(".regions.bed.gz", ""))
            try:
                one = load_one_regions(f, s_name)
                sample_frames.append(one)
                total = one.shape[0]
                numeric = one[s_name].notna().sum()
                zeros = (one[s_name] == 0).sum() if numeric > 0 else 0
                debug_rows.append({
                    'sample': s_name,
                    'file': os.path.basename(f),
                    'amplicons_total': int(total),
                    'numeric_depth': int(numeric),
                    'zeros': int(zeros),
                    'numeric_fraction': float(numeric / total) if total else 0.0
                })
            except Exception as e:
                print(f"[ERROR] Failed to load {f}: {e}")
                # Create a fully-NaN column for visibility, but mark it
                dummy = pd.DataFrame(columns=['chr','start','end', s_name])
                sample_frames.append(dummy.set_index([]))

        if not sample_frames:
            print(f"[!] No usable coverage frames for {cohort}")
            continue

        combined = join_by_id(sample_frames).reset_index()
        cohort_dfs[cohort] = combined

        # Persist cohort CSV and debug audit
        csv_name = f"GSDMB_{cohort.replace(' ', '_')}_Full_Panel.csv"
        combined.to_csv(os.path.join(PATHS['out_dir'], csv_name), index=False)

        dbg = pd.DataFrame(debug_rows)
        dbg_name = f"DEBUG_cohort_coverage_audit_{cohort.replace(' ', '_')}.csv"
        dbg.to_csv(os.path.join(PATHS['out_dir'], dbg_name), index=False)

    if not cohort_dfs:
        print("[!] ERROR: No data found in directories.")
        return

    full_qc_df = pd.concat(master_qc_list) if master_qc_list else pd.DataFrame()

    # Global matrix: union of amplicons; keep NaN (missing) so they’re masked in heatmaps
    all_cov = pd.concat(
        [
            df.set_index("id").drop(columns=['chr', 'start', 'end'], errors='ignore')
            for df in cohort_dfs.values()
        ],
        axis=1
    )

    failing_samples = full_qc_df[full_qc_df["status"] == "Fail"]["sample_norm"].tolist()

    # ==============================================================================
    # 3. PLOTTING (unchanged, but with masks where appropriate)
    # ==============================================================================

    print("[*] Phase 2: Generating Harmonised TFM Visualisations (16 Figures)...")

    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    # 01. PCA
    if not full_qc_df.empty:
        feats = ['mean_cov', 'on_target_pct', 'total_reads', 'mapped_pct']
        x_scaled = StandardScaler().fit_transform(full_qc_df[feats].fillna(0))
        pca_res = PCA(n_components=2).fit_transform(x_scaled)

        full_qc_df["PC1"], full_qc_df["PC2"] = pca_res[:, 0], pca_res[:, 1]

        plt.figure(figsize=(8, 6))
        ax = sns.scatterplot(
            data=full_qc_df,
            x="PC1", y="PC2",
            hue="cohort", style="status",
            palette=PALETTE, s=100
        )
        plt.xlabel("Principal Component 1")
        plt.ylabel("Principal Component 2")
        plt.title("Technical Map (PCA)")
        ax.legend(title="Cohort and Status", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/01_pca.png")

    # 02–04. Cohort heatmaps (mask NaN)
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in ['chr', 'start', 'end', 'id']]
        dmat = df.set_index("id")[s_cols]
        logmat = np.log10(dmat + 1)
        mask = dmat.isna()

        plt.figure(figsize=(12, 8))
        sns.heatmap(logmat, cmap="magma", mask=mask, xticklabels=True, vmin=0, vmax=4)
        plt.xlabel("Sample Identifier")
        plt.ylabel("Amplicon Genomic Position")
        plt.title(f"Heatmap – {cohort} Depth Intensity")

        ax = plt.gca()
        plt.xticks(rotation=45, ha='right', fontsize=8)
        for label in ax.get_xticklabels():
            if normalise_sample(label.get_text()) in failing_samples:
                label.set_color('red'); label.set_weight('bold')

        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/02_heatmap_{cohort.replace(' ', '_')}.png")

    # 05. Global heatmap (mask NaN)
    plt.figure(figsize=(16, 10))
    sns.heatmap(np.log10(all_cov + 1), cmap="magma", mask=all_cov.isna(),
                xticklabels=True, vmin=0, vmax=4)
    plt.xlabel("Global Sample Set")
    plt.ylabel("Amplicon Genomic Position")
    plt.title("Global Overall Coverage Intensity")

    ax = plt.gca()
    plt.xticks(rotation=45, ha='right', fontsize=6)
    for label in ax.get_xticklabels():
        if normalise_sample(label.get_text()) in failing_samples:
            label.set_color('red'); label.set_weight('bold')
    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/03_heatmap_global_all.png")

    # 06–16 unchanged except NaN-safe logic
    if not full_qc_df.empty:
        plt.figure(figsize=(10, 6))
        sns.barplot(data=full_qc_df, x='cohort', y='mapped_pct', hue='cohort',
                    palette=PALETTE, errorbar='sd', legend=False)
        plt.axhline(QC_LIMITS['mapped_pct'], color='red', ls='--')
        plt.xlabel("Patient Cohort"); plt.ylabel("Mapped Reads (%)")
        plt.title("Mean Mapping Rate per Cohort")
        plt.savefig(f"{PATHS['plots_dir']}/06_mapped_reads_bar.png")

    plt.figure(figsize=(15, 6))
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in ['chr', 'start', 'end', 'id']]
        med = df[s_cols].median(axis=1)
        plt.plot(range(len(df)), med, label=cohort, color=PALETTE[cohort], lw=2)
    plt.axhline(QC_LIMITS['worst_amplicon_floor'], color='red', ls='--')
    plt.axhline(QC_LIMITS['mean_cov'], color='gold', ls=':')
    plt.yscale('log'); plt.xlabel("Amplicon Position (Index)"); plt.ylabel("Median Depth (X)")
    plt.title("Panel Landscape vs. Validated Thresholds")
    plt.legend(title="Cohort", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout(); plt.savefig(f"{PATHS['plots_dir']}/07_panel_landscape.png")

    if not full_qc_df.empty:
        plt.figure(figsize=(9, 6))
        sns.violinplot(data=full_qc_df, x='cohort', y='on_target_pct',
                       hue='cohort', palette=PALETTE, inner="quart", legend=False)
        plt.axhline(QC_LIMITS['on_target_pct'], color='darkred', ls='--')
        plt.xlabel("Patient Cohort"); plt.ylabel("On-Target Reads (%)")
        plt.title("Sequencing Specificity Distribution")
        plt.savefig(f"{PATHS['plots_dir']}/08_specificity_violin.png")

    # 09. Gaps audit (ignore NaN)
    gap_stats = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in ['chr', 'start', 'end', 'id']]
        for s in s_cols:
            vals = df[s]
            gaps = (vals < QC_LIMITS['worst_amplicon_floor']) & (~vals.isna())
            gap_stats.append({'cohort': cohort, 'sample': s, 'gaps': int(gaps.sum())})
    plt.figure(figsize=(9, 6))
    sns.boxplot(data=pd.DataFrame(gap_stats), x='cohort', y='gaps',
                hue='cohort', palette=PALETTE, legend=False)
    plt.xlabel("Patient Cohort"); plt.ylabel("Failed Amplicons (<50x)")
    plt.title("Coverage Gaps")
    plt.savefig(f"{PATHS['plots_dir']}/09_gaps_audit.png")

    if not full_qc_df.empty:
        plt.figure(figsize=(10, 7))
        sns.scatterplot(data=full_qc_df, x='mean_cov', y='on_target_pct',
                        hue='status', style='cohort', palette=STATUS_PALETTE, s=120)
        plt.axvline(QC_LIMITS['mean_cov'], color='gold', ls=':')
        plt.axhline(QC_LIMITS['on_target_pct'], color='red', ls='--')
        plt.xlabel("Mean Coverage Depth (X)"); plt.ylabel("On-Target Reads (%)")
        plt.title("Mean Coverage vs. Specificity")
        plt.legend(title="Status", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout(); plt.savefig(f"{PATHS['plots_dir']}/10_correlation.png")

    if not full_qc_df.empty:
        plt.figure(figsize=(10, 7))
        full_qc_df.groupby(['cohort', 'status']).size().unstack().fillna(0).plot(
            kind='bar', stacked=True, color=STATUS_PALETTE, ax=plt.gca()
        )
        plt.xlabel("Patient Cohort"); plt.ylabel("Sample Count")
        plt.title("Retention Rate")
        plt.legend(title="Status"); plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/11_yield_bar.png")

    plt.figure(figsize=(15, 5))
    cv = all_cov.std(axis=1) / all_cov.mean(axis=1)
    plt.bar(range(len(cv)), cv, color='dimgray', alpha=0.7)
    plt.xlabel("Amplicon Position (Index)"); plt.ylabel("Coefficient of Variation (CV)")
    plt.title("Systemic Amplicon Instability (CV)")
    plt.savefig(f"{PATHS['plots_dir']}/12_variance.png")

    if not full_qc_df.empty:
        plt.figure(figsize=(8, 6))
        sns.stripplot(data=full_qc_df, x='cohort', y='mapped_pct',
                      hue='status', palette=STATUS_PALETTE,
                      size=8, alpha=0.7, jitter=True)
        plt.axhline(QC_LIMITS['mapped_pct'], color='red', ls='--')
        plt.xlabel("Patient Cohort"); plt.ylabel("Mapped Reads (%)")
        plt.title("Mapping Efficiency Audit")
        plt.legend(title="Status", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout(); plt.savefig(f"{PATHS['plots_dir']}/13_mapping_audit.png")

    if not full_qc_df.empty:
        plt.figure(figsize=(10, 6))
        sns.stripplot(data=full_qc_df, x='cohort', y='on_target_pct',
                      hue='status', palette=STATUS_PALETTE,
                      jitter=True, size=7)
        plt.axhline(QC_LIMITS['on_target_pct'], color='firebrick', ls='--')
        plt.xlabel("Patient Cohort"); plt.ylabel("On-Target Reads (%)")
        plt.title("Sequencing Specificity Jitter")
        plt.legend(title="Status", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout(); plt.savefig(f"{PATHS['plots_dir']}/14_specificity_jitter.png")

    worst_data = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in ['chr', 'start', 'end', 'id']]
        for s in s_cols:
            if s in df:
                md = df[s].min()
                if pd.notna(md):
                    worst_data.append({'cohort': cohort, 'min_depth': md})
    worst_df = pd.DataFrame(worst_data)
    plt.figure(figsize=(10, 6))
    if not worst_df.empty:
        sns.swarmplot(data=worst_df, x='cohort', y='min_depth',
                      hue='cohort', palette=PALETTE, size=6)
    plt.axhline(QC_LIMITS['worst_amplicon_floor'], color='red', ls='--')
    plt.yscale('symlog', linthresh=10); plt.xlabel("Patient Cohort")
    plt.ylabel("Lowest Amplicon Depth (X)")
    plt.title("Lowest Coverage per Sample")
    h, l = plt.gca().get_legend_handles_labels()
    if h: plt.legend(h, l, title="Cohort", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout(); plt.savefig(f"{PATHS['plots_dir']}/15_worst_amplicon_depth.png")

    fail_mask = (all_cov < QC_LIMITS['worst_amplicon_floor']) & (~all_cov.isna())
    fail_counts = fail_mask.sum(axis=1).reset_index()
    fail_counts.columns = ['id', 'fail_count']
    if not fail_counts.empty:
        plt.figure(figsize=(12, 10))
        top_fails = fail_counts.sort_values('fail_count', ascending=False).head(30)
        top_fails = top_fails[top_fails['fail_count'] > 0]
        if not top_fails.empty:
            sns.barplot(data=top_fails, x='fail_count', y='id', palette="Reds_r")
            plt.title("Systematic Amplicon Failure Audit (<50x Depth)")
            plt.xlabel("Number of Samples Failed"); plt.ylabel("Amplicon Coordinate (ID)")
            plt.tight_layout(); plt.savefig(f"{PATHS['plots_dir']}/16_systematic_amplicon_failures.png")

    plt.close('all')
    print("\n[SUCCESS] Harmonised CSVs, diagnostics, and 16 visualisations finalised.")

# ==============================================================================
# MAIN
# ==============================================================================
if __name__ == "__main__":
    run_validated_pipeline()
