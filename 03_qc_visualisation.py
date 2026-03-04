#!/usr/bin/env python3
"""
Generates 21 QC visualisation plots from mosdepth coverage outputs and
qc_summary.tsv files produced by scripts 02 and 02b. Covers metric
distributions, PCA, coverage heatmaps, panel landscapes, amplicon failure
analysis, and sample retention summaries.

Configuration: edit PATHS and COHORT_MAP below before running.
"""

import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.ticker as mticker


# --- QC thresholds (must match values used in script 02) ---------------------
QC_LIMITS = {
    'total_reads':           120000,
    'mapped_pct':            90.0,
    'on_target_pct':         75.0,
    'mean_cov':              200.0,
    'uniformity_100x':       95.0,
    'worst_amplicon_floor':  50.0   # minimum acceptable depth in any single amplicon
}

# --- File paths --------------------------------------------------------------
PATHS = {
    "out_dir":     "./gsdmb_final_results",
    "samples_dir": "./gsdmb_final_results/individual_sample_coverage",
    "plots_dir":   "./gsdmb_final_results/plots",
    "raw_root":    "./",
    # Optional: path to a single BED file for rsID annotations.
    # If None, the script looks for dna_qc/targets.sorted.bed within each cohort directory.
    "bed_file":    None
}

# --- Cohort labels mapped to relative directory paths -----------------------
COHORT_MAP = {
    'Breast Tumour':        'breast/tumour',
    'Breast Control':       'breast/normal',
    'Endometrium Tumour':   'endometrium/tumour',
    'Endometrium Control':  'endometrium/normal',
}

# --- Colour palettes ---------------------------------------------------------
PALETTE = {
    'Breast Tumour':       '#3498db',
    'Breast Control':      '#d1d946',
    'Endometrium Tumour':  '#9b59b6',
    'Endometrium Control': '#2ecc71'
}

STATUS_PALETTE = {'Pass': '#2ecc71', 'Fail': '#e74c3c'}

# When True, raises an error if a mosdepth file is empty or mostly non-numeric;
# set to False only during development to tolerate incomplete outputs
STRICT_MOSDEPTH        = True
MIN_NUMERIC_FRACTION   = 0.90   # minimum fraction of depth values that must be numeric


# =============================================================================
# Utility functions
# =============================================================================

def normalise_sample(s: str) -> str:
    """Standardise sample names to lowercase with underscores for consistent matching."""
    s = str(s).strip().lower()
    return (s.replace(" ", "_")
             .replace(".regions", "")
             .replace(".bed", "")
             .replace(".gz", ""))


def normalise_chr(c: str) -> str:
    """Ensure all chromosome names carry the 'chr' prefix (e.g. '17' -> 'chr17')."""
    c = str(c)
    if c.startswith("chr"):
        return c
    if c in {"X", "Y", "M", "MT"}:
        return "chrM" if c in {"M", "MT"} else f"chr{c}"
    try:
        return f"chr{int(c)}"
    except Exception:
        return c if c.startswith("chr") else f"chr{c}"


def pick_depth_column(df: pd.DataFrame) -> int:
    """
    Identify the column index containing coverage depth in a mosdepth BED file.
    Standard 4-column output has depth in column 3; extended formats (produced
    when a BED annotation file is passed to mosdepth) may shift this column.
    The function selects the rightmost column that is at least 80% numeric,
    falling back to a column-by-column search if that heuristic fails.
    """
    ncols = df.shape[1]
    if ncols == 4:
        return 3
    if ncols < 4:
        return -1
    # Check the last column first — mosdepth appends depth at the end
    frac = pd.to_numeric(df.iloc[:, -1], errors='coerce').notna().mean()
    if frac > 0.8:
        return ncols - 1
    # Fall back to scanning columns 3 onwards
    best, best_frac = -1, -1
    for i in range(3, ncols):
        f = pd.to_numeric(df.iloc[:, i], errors='coerce').notna().mean()
        if f > best_frac:
            best, best_frac = i, f
    return best if best_frac >= 0.5 else -1


def create_display_label(row: pd.Series) -> str:
    """
    Return the annotation label for an amplicon row, falling back to genomic
    coordinates if the annotation field is absent or uninformative.
    Annotation priority (set during loading): rsID > gene name > amplicon ID.
    """
    annot = row.get("annotation", "")
    if pd.notna(annot) and str(annot).strip() not in ["", "."]:
        return str(annot).strip()
    return f"{row['chr']}:{row['start']}-{row['end']}"


# =============================================================================
# Data loading
# =============================================================================

def load_annotation_for_cohort(root_path: str, bed_file_override: str = None) -> pd.DataFrame:
    """
    Load the target BED file and extract the best available annotation per amplicon.
    Priority: rsID (col 5) > gene name (col 6) > amplicon ID (col 4).
    A global BED file path takes precedence over the per-cohort default location.
    """
    if bed_file_override and os.path.exists(bed_file_override):
        annot_file = bed_file_override
    else:
        annot_file = os.path.join(root_path, "dna_qc", "targets.sorted.bed")
        if not os.path.exists(annot_file):
            raise FileNotFoundError(f"Annotation file not found: {annot_file}")

    with open(annot_file, 'r') as f:
        first_line = f.readline()
    skip_rows = 1 if first_line.startswith('track') else 0   # skip UCSC track header if present

    ann = pd.read_csv(annot_file, sep="\t", header=None, skiprows=skip_rows)

    if ann.shape[1] >= 6:
        ann.columns = (["chr", "start", "end", "amplicon_id", "rsid", "gene"] +
                       [f"extra_{i}" for i in range(ann.shape[1] - 6)])

        def make_annotation(row):
            rsid   = row["rsid"]
            gene   = row["gene"]
            amp_id = row["amplicon_id"]
            if pd.notna(rsid) and str(rsid).strip() not in ["", "."]:
                return str(rsid).strip()
            if pd.notna(gene) and str(gene).strip() not in ["", "."]:
                return str(gene).strip()
            return str(amp_id)

        ann["annotation"] = ann.apply(make_annotation, axis=1)

    elif ann.shape[1] >= 5:
        ann.columns = (["chr", "start", "end", "amplicon_id", "annotation"] +
                       [f"extra_{i}" for i in range(ann.shape[1] - 5)])
    else:
        ann.columns = ["chr", "start", "end", "annotation"]

    ann["chr"]      = ann["chr"].astype(str)
    ann["annot_id"] = ann["annotation"]

    return ann[["chr", "start", "end", "annotation", "annot_id"]].drop_duplicates(["chr", "start", "end"])


def load_coverage_for_file(path: str, sample_name: str) -> pd.DataFrame:
    """
    Load one mosdepth .regions.bed.gz file, normalise chromosome names,
    identify the depth column, and return a DataFrame indexed by region ID.
    Raises ValueError if the file is empty or fails the numeric fraction check
    (only when STRICT_MOSDEPTH is True).
    """
    raw = pd.read_csv(path, sep="\t", header=None, compression="gzip", dtype=str)
    if raw.empty:
        if STRICT_MOSDEPTH:
            raise ValueError(f"{path} is empty")
        return pd.DataFrame(columns=["chr", "start", "end", sample_name])

    raw.iloc[:, 0] = raw.iloc[:, 0].map(normalise_chr)
    raw.iloc[:, 1] = pd.to_numeric(raw.iloc[:, 1], errors='coerce')
    raw.iloc[:, 2] = pd.to_numeric(raw.iloc[:, 2], errors='coerce')

    depth_col = pick_depth_column(raw)
    depth = (pd.to_numeric(raw.iloc[:, depth_col], errors='coerce')
             if depth_col != -1 else pd.Series(np.nan, index=raw.index))

    df = pd.DataFrame({
        "chr":       raw.iloc[:, 0].astype(str),
        "start":     raw.iloc[:, 1].astype("Int64"),
        "end":       raw.iloc[:, 2].astype("Int64"),
        sample_name: depth
    })
    # A string key (chr:start-end) allows region-level joining across samples
    df["id"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    df = df.drop_duplicates("id").set_index("id")

    frac = df[sample_name].notna().mean()
    if frac < MIN_NUMERIC_FRACTION and STRICT_MOSDEPTH:
        raise ValueError(f"{sample_name}: only {frac:.1%} depth values are numeric")
    return df


def load_coverage_for_cohort(cohort_name: str, root_path: str) -> pd.DataFrame:
    """
    Discover all mosdepth .regions.bed.gz files within a cohort directory,
    load each one, and outer-join them on region ID to produce a wide matrix
    (rows = amplicons, columns = samples).
    """
    pattern = os.path.join(root_path, "**", "*.regions.bed.gz")
    files = [f for f in glob.glob(pattern, recursive=True)
             if f.endswith(".regions.bed.gz") and ".bai" not in f.lower()]
    if not files:
        return pd.DataFrame()

    frames = []
    for f in sorted(files):
        name = normalise_sample(os.path.basename(f).replace(".regions.bed.gz", ""))
        try:
            frames.append(load_coverage_for_file(f, name))
        except Exception as e:
            raise RuntimeError(f"[MOSDEPTH ERROR] {cohort_name} — {name}: {e}")

    # Outer join so that amplicons missing in one sample are retained as NaN
    combined = frames[0].copy()
    for frame in frames[1:]:
        cols = [c for c in frame.columns if c not in ["chr", "start", "end"]]
        combined = combined.join(frame[cols], how="outer")

    return combined.reset_index()


def load_qc_for_cohort(cohort_name: str, root_path: str) -> pd.DataFrame:
    """Load qc_summary.tsv from script 02 for a single cohort."""
    qc_file = os.path.join(root_path, "dna_qc", "qc_summary.tsv")
    if not os.path.exists(qc_file):
        return pd.DataFrame()
    df = pd.read_csv(qc_file, sep="\t")
    df["cohort"]      = cohort_name
    df["status"]      = df["status"].str.capitalize()   # normalise PASS/FAIL to Pass/Fail
    df["sample_norm"] = df["sample"].apply(normalise_sample)
    return df


# =============================================================================
# Pipeline assembly
# =============================================================================

def build_pipeline():
    """
    Load QC summaries, per-sample coverage, and BED annotations for all cohorts.
    Returns:
        full_qc      — concatenated QC metrics across all cohorts
        cohort_dfs   — dict mapping cohort name to annotated coverage matrix
        annotation_df — deduplicated annotation table
        all_cov      — global coverage matrix (all cohorts combined)
        failing_samples — list of normalised sample names that failed script 02 QC
    """
    for key, p in PATHS.items():
        if key != "bed_file" and p and isinstance(p, str):
            os.makedirs(p, exist_ok=True)

    qc_frames, annotation_frames = [], []
    cohort_dfs = {}

    # Load a single global BED file if specified; otherwise each cohort loads its own
    global_bed = None
    if PATHS.get("bed_file") and os.path.exists(PATHS["bed_file"]):
        print(f"[*] Loading global BED: {PATHS['bed_file']}")
        global_bed = load_annotation_for_cohort("", bed_file_override=PATHS["bed_file"])
        annotation_frames.append(global_bed)

    print("[*] Loading data...")
    for cohort, relpath in COHORT_MAP.items():
        root = os.path.join(PATHS["raw_root"], relpath)
        if not os.path.exists(root):
            print(f"[WARN] Path not found, skipping: {root}")
            continue

        qc = load_qc_for_cohort(cohort, root)
        if not qc.empty:
            qc_frames.append(qc)

        cov = load_coverage_for_cohort(cohort, root)
        if cov.empty:
            print(f"[WARN] No coverage data found for {cohort}")
            continue

        ann = global_bed.copy() if global_bed is not None else load_annotation_for_cohort(root)
        if global_bed is None:
            annotation_frames.append(ann)

        # Merge annotation onto coverage matrix by genomic coordinates
        cov = cov.merge(ann[["chr", "start", "end", "annotation", "annot_id"]],
                        on=["chr", "start", "end"], how="left")
        cohort_dfs[cohort] = cov

        # Save full per-cohort coverage matrix for downstream use
        cov.to_csv(os.path.join(PATHS["out_dir"],
                                f"GSDMB_{cohort.replace(' ', '_')}_FullPanel.csv"),
                   index=False)

        # Write a per-cohort audit table recording how many depth values are numeric
        # (a sanity check that mosdepth outputs were parsed correctly)
        meta_cols = {"chr", "start", "end", "id", "annotation", "annot_id"}
        s_cols = [c for c in cov.columns if c not in meta_cols]
        audit_rows = []
        for s in s_cols:
            depth = cov[s]
            tot   = len(depth)
            audit_rows.append({
                "sample":           s,
                "amplicons_total":  tot,
                "numeric_depth":    depth.notna().sum(),
                "zeros":            (depth == 0).sum(),
                "numeric_fraction": depth.notna().sum() / tot if tot else 0
            })
        pd.DataFrame(audit_rows).to_csv(
            os.path.join(PATHS["out_dir"], f"CoverageAudit_{cohort.replace(' ', '_')}.csv"),
            index=False)

    full_qc = pd.concat(qc_frames, ignore_index=True) if qc_frames else pd.DataFrame()

    annotation_df = (
        pd.concat(annotation_frames, ignore_index=True)
          .drop_duplicates(["chr", "start", "end"])
        if annotation_frames
        else pd.DataFrame(columns=["chr", "start", "end", "annotation", "annot_id"])
    )

    # Build a global coverage matrix by concatenating all cohort matrices column-wise
    all_cov_frames = []
    for df in cohort_dfs.values():
        tmp = df.set_index("id").drop(
            columns=["chr", "start", "end", "annotation", "annot_id"], errors="ignore")
        all_cov_frames.append(tmp)
    all_cov = pd.concat(all_cov_frames, axis=1) if all_cov_frames else pd.DataFrame()

    failing_samples = (full_qc[full_qc["status"] == "Fail"]["sample_norm"].tolist()
                       if not full_qc.empty else [])

    print("[*] Data loading complete.")
    return full_qc, cohort_dfs, annotation_df, all_cov, failing_samples

# =============================================================================
# Plotting functions
# =============================================================================

def plot_qc_distributions(full_qc_df: pd.DataFrame, plots_dir: str):
    """
    Plots 01–05: Per-metric distributions as paired violin + box plots.
    Left panel shows distribution by cohort; right panel compares PASS vs FAIL.
    One figure is saved per metric.
    """
    if full_qc_df.empty:
        print("[WARN] No QC data for distribution plots"); return

    print("[*] Generating QC distribution plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    metrics = [
        ('total_reads',    'Total Reads',          QC_LIMITS['total_reads'],    '01'),
        ('mapped_pct',     'Mapped %',              QC_LIMITS['mapped_pct'],     '02'),
        ('on_target_pct',  'On Target %',           QC_LIMITS['on_target_pct'], '03'),
        ('mean_cov',       'Mean Coverage',         QC_LIMITS['mean_cov'],       '04'),
        ('uniformity_100x','Uniformity >100× %',    QC_LIMITS['uniformity_100x'],'05'),
    ]

    for metric, label, threshold, plot_num in metrics:
        if metric not in full_qc_df.columns:
            print(f"[WARN] Metric '{metric}' not in QC table — skipping"); continue

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Left: distribution per cohort — violin shows full density, useful for
        # detecting bimodality or cohort-level batch effects
        sns.violinplot(data=full_qc_df, x='cohort', y=metric,
                       hue='cohort', palette=PALETTE, ax=axes[0], legend=False)
        axes[0].axhline(threshold, color='red', linestyle='--', linewidth=2,
                        label=f'Threshold: {threshold}')
        axes[0].set_title(f'{label} Distribution by Cohort', fontweight='bold')
        axes[0].set_xlabel('Cohort', fontweight='bold')
        axes[0].set_ylabel(label, fontweight='bold')
        axes[0].legend()
        axes[0].tick_params(axis='x', rotation=45)

        # Right: box plot stratified by QC status — confirms that the threshold
        # cleanly separates PASS and FAIL distributions
        if 'status' in full_qc_df.columns:
            sns.boxplot(data=full_qc_df, x='status', y=metric,
                        hue='status', palette=STATUS_PALETTE, ax=axes[1], legend=False)
            axes[1].axhline(threshold, color='red', linestyle='--', linewidth=2,
                            label=f'Threshold: {threshold}')
            axes[1].set_title(f'{label} by QC Status', fontweight='bold')
            axes[1].set_xlabel('QC Status', fontweight='bold')
            axes[1].set_ylabel(label, fontweight='bold')
            axes[1].legend()

        plt.tight_layout()
        plt.savefig(f"{plots_dir}/{plot_num}_qc_{metric}_distribution.png", dpi=300)
        plt.close()


def plot_pca_analysis(full_qc_df: pd.DataFrame, plots_dir: str):
    """
    Plots 06–07: PCA of QC metrics.
    PCA reduces five correlated QC dimensions to two principal components,
    making it easier to spot outlier samples and cohort-level clustering.
    Plot 06 includes all samples; plot 07 is restricted to PASS samples to
    reveal within-cohort structure without the influence of failing outliers.
    """
    if full_qc_df.empty:
        print("[WARN] No QC data for PCA"); return

    print("[*] Generating PCA plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    full_qc_df = full_qc_df.copy()
    # Invert uniformity so that all features increase with worsening QC,
    # which makes PCA loadings easier to interpret
    full_qc_df["uniformity_gap"] = 100.0 - full_qc_df["uniformity_100x"]
    feats = ['mean_cov', 'on_target_pct', 'total_reads', 'mapped_pct', 'uniformity_gap']

    # Plot 06: All samples
    X   = StandardScaler().fit_transform(full_qc_df[feats].fillna(0))
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)
    var1, var2 = pca.explained_variance_ratio_[0] * 100, pca.explained_variance_ratio_[1] * 100

    full_qc_df["PC1"], full_qc_df["PC2"] = pcs[:, 0], pcs[:, 1]

    plt.figure(figsize=(10, 7))
    ax = sns.scatterplot(data=full_qc_df, x="PC1", y="PC2",
                         hue="cohort", style="status", palette=PALETTE, s=100)
    ax.set_xlabel(f"PC1 ({var1:.1f}% variance)", fontweight='bold')
    ax.set_ylabel(f"PC2 ({var2:.1f}% variance)", fontweight='bold')
    plt.title("PCA of QC Metrics — All Samples", fontweight='bold', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/06_pca_all_samples.png", dpi=300)
    plt.close()

    # Plot 07: PASS samples only — requires at least 3 points for meaningful PCA
    pass_only = full_qc_df[full_qc_df["status"] == "Pass"].copy()
    if len(pass_only) >= 3:
        X   = StandardScaler().fit_transform(pass_only[feats].fillna(0))
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X)
        var1, var2 = pca.explained_variance_ratio_[0] * 100, pca.explained_variance_ratio_[1] * 100

        pass_only["PC1"], pass_only["PC2"] = pcs[:, 0], pcs[:, 1]

        plt.figure(figsize=(10, 7))
        ax = sns.scatterplot(data=pass_only, x="PC1", y="PC2",
                             hue="cohort", palette=PALETTE, s=100)
        ax.set_xlabel(f"PC1 ({var1:.1f}% variance)", fontweight='bold')
        ax.set_ylabel(f"PC2 ({var2:.1f}% variance)", fontweight='bold')
        plt.title("PCA of QC Metrics — Pass Samples Only", fontweight='bold', fontsize=14)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/07_pca_pass_samples.png", dpi=300)
        plt.close()


def plot_coverage_heatmaps(cohort_dfs: dict, all_cov: pd.DataFrame,
                           failing_samples: list, plots_dir: str):
    """
    Plots 08–09: Coverage heatmaps (amplicons × samples).
    Depth is log10-transformed before plotting to compress the wide dynamic
    range typical of amplicon sequencing (some amplicons may differ by >100×).
    Failing samples are highlighted in red on the x-axis.
    Plot 08: one heatmap per cohort. Plot 09: global heatmap across all cohorts.
    """
    print("[*] Generating coverage heatmaps...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    meta = {"chr", "start", "end", "id", "annotation", "annot_id"}
    clean_failing = {normalise_sample(s) for s in failing_samples}

    # Plot 08: per-cohort heatmaps
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]

        df_plot = df.copy()
        df_plot["display_label"] = df_plot.apply(create_display_label, axis=1)
        dmat    = df_plot.set_index("display_label")[s_cols]
        plotmat = np.log10(dmat.fillna(0) + 1)   # +1 avoids log(0); grey cells = missing data

        cmap = sns.color_palette("magma", as_cmap=True)
        cmap.set_bad(color="#B0B0B0")

        # Scale figure height proportionally to amplicon count to keep rows legible
        fig_height = max(8, min(20, len(dmat) * 0.15))
        plt.figure(figsize=(20, fig_height))
        ax = sns.heatmap(plotmat, cmap=cmap, mask=dmat.isna(), vmin=0, vmax=4,
                         cbar_kws={'label': 'log10(Coverage + 1)'})

        ax.xaxis.set_major_locator(mticker.FixedLocator(np.arange(len(s_cols)) + 0.5))
        ax.set_xticklabels(s_cols, rotation=90, ha='center', fontsize=6)
        ax.tick_params(axis='x', which='both', length=0, pad=1)

        for x in range(len(s_cols) + 1):
            ax.axvline(x, color="white", lw=0.3, alpha=0.6)

        # Limit y-axis labels when there are many amplicons to avoid overplotting
        n_yticks = min(len(dmat), 50)
        if len(dmat) > n_yticks:
            step = len(dmat) // n_yticks
            positions = list(range(0, len(dmat), step))
            ax.set_yticks([i + 0.5 for i in positions])
            ax.set_yticklabels([dmat.index[i] for i in positions], fontsize=7)
        else:
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)

        # Colour-code failing sample labels on the x-axis
        for label in ax.get_xticklabels():
            if normalise_sample(label.get_text()) in clean_failing:
                label.set_color("red")
                label.set_weight("bold")

        plt.xlabel("Sample ID", fontsize=10, fontweight='bold')
        plt.ylabel("Amplicon (rsID / Annotation)", fontsize=10, fontweight='bold')
        plt.title(f"Coverage Heatmap: {cohort}", fontweight='bold', fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/08_heatmap_{cohort.replace(' ', '_')}.png", dpi=300)
        plt.close()

    # Plot 09: global heatmap combining all cohorts
    if not all_cov.empty:
        plt.figure(figsize=(30, 12))
        cmap = sns.color_palette("magma", as_cmap=True)
        cmap.set_bad(color="#B0B0B0")
        sns.heatmap(np.log10(all_cov.fillna(0) + 1),
                    cmap=cmap, mask=all_cov.isna(), vmin=0, vmax=4,
                    cbar_kws={'label': 'log10(Coverage + 1)'})
        ax = plt.gca()
        plt.xticks(rotation=90, fontsize=5)
        for label in ax.get_xticklabels():
            if normalise_sample(label.get_text()) in clean_failing:
                label.set_color("red")
                label.set_weight("bold")
        plt.xlabel("Sample ID", fontsize=10, fontweight='bold')
        plt.ylabel("Genomic Regions", fontsize=10, fontweight='bold')
        plt.title("Global Coverage Heatmap — All Cohorts", fontweight='bold', fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/09_heatmap_global.png", dpi=300)
        plt.close()


def plot_panel_landscapes(cohort_dfs: dict, plots_dir: str):
    """
    Plots 10–12: Median coverage per amplicon across three orderings.
    10: amplicon index order (as they appear in the BED file)
    11: genomic coordinate order (chr1 -> chrX, sorted by position)
    12: grouped by annotated gene/region
    Log scale is used on the y-axis because coverage can span several orders
    of magnitude across amplicons, and log scale preserves low-coverage detail.
    """
    print("[*] Generating panel landscape plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    meta       = {"chr", "start", "end", "id", "annotation", "annot_id"}
    mean_line  = QC_LIMITS['mean_cov']
    floor_line = QC_LIMITS['worst_amplicon_floor']

    # Plot 10: index order
    plt.figure(figsize=(15, 6))
    ax = plt.gca()
    cohort_handles = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        med    = df[s_cols].median(axis=1)
        line, = ax.plot(range(len(med)), med, lw=2, label=cohort, color=PALETTE[cohort])
        cohort_handles.append(line)

    h_mean  = ax.axhline(mean_line,  color="gold", ls=":",  lw=2,
                         label=f"Mean Coverage Target >= {mean_line}x")
    h_floor = ax.axhline(floor_line, color="red",  ls="--", lw=2,
                         label=f"Minimum Floor {floor_line}x")
    ax.set_yscale("log")
    ax.set_xlabel("Amplicon Index (Order in File)", fontweight='bold')
    ax.set_ylabel("Median Coverage Depth (x)", fontweight='bold')
    ax.set_title("Panel Coverage Landscape — Index Order", fontweight='bold', fontsize=14)
    ax.grid(alpha=0.3)
    handles = cohort_handles + [h_mean, h_floor]
    ax.legend(handles, [h.get_label() for h in handles], bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/10_landscape_index_order.png", dpi=300)
    plt.close()

    # Helper: numeric sort key for chromosome names
    def _chr_order(c):
        c = str(c).replace("chr", "")
        mapping = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        try:
            return int(c)
        except Exception:
            return mapping.get(c, 1000)

    # Plot 11: genomic coordinate order
    coord_frames = []
    for df in cohort_dfs.values():
        tmp = df[["chr", "start", "end", "id", "annotation"]].drop_duplicates()
        tmp = tmp.copy()
        tmp["display_label"] = tmp.apply(create_display_label, axis=1)
        coord_frames.append(tmp)

    kdf = pd.concat(coord_frames).drop_duplicates(subset=["chr", "start", "end"])
    kdf["key"] = kdf["chr"].apply(_chr_order)
    kdf = kdf.sort_values(["key", "start", "end"])
    ordered_ids    = kdf["id"].tolist()
    ordered_labels = kdf["display_label"].tolist()

    plt.figure(figsize=(18, 6))
    ax = plt.gca()
    cohort_handles = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        med    = df[s_cols].median(axis=1)
        t = (pd.DataFrame({"id": df["id"], "med": med})
               .drop_duplicates("id").set_index("id").reindex(ordered_ids))
        line, = ax.plot(np.arange(len(ordered_ids)), t["med"].values,
                        lw=2, label=cohort, color=PALETTE[cohort])
        cohort_handles.append(line)

    h_mean  = ax.axhline(mean_line,  color="gold", ls=":",  lw=2,
                         label=f"Mean Coverage Target >= {mean_line}x")
    h_floor = ax.axhline(floor_line, color="red",  ls="--", lw=2,
                         label=f"Minimum Floor {floor_line}x")

    if ordered_labels:
        n    = len(ordered_labels)
        step = max(1, n // 40)   # sample every nth label to avoid x-axis crowding
        sel  = list(range(0, n, step))
        ax.set_xticks(sel)
        ax.set_xticklabels([ordered_labels[i] for i in sel], rotation=45, ha='right', fontsize=7)

    ax.set_yscale("log")
    ax.set_xlabel("Genomic Position (rsID / Annotation)", fontweight='bold')
    ax.set_ylabel("Median Coverage Depth (x)", fontweight='bold')
    ax.set_title("Panel Coverage Landscape — Genomic Order", fontweight='bold', fontsize=14)
    ax.grid(alpha=0.3)
    handles = cohort_handles + [h_mean, h_floor]
    ax.legend(handles, [h.get_label() for h in handles], bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/11_landscape_genomic_order.png", dpi=300)
    plt.close()

    # Plot 12: grouped by gene / annotation label
    def parse_gene_or_annotation(a):
        if pd.isna(a):
            return None
        t = str(a)
        return t.split("\n")[1] if "\n" in t else t.strip()

    gene_frames = []
    for cohort, df in cohort_dfs.items():
        if "annotation" not in df.columns:
            continue
        tmp = df[["annotation", "chr", "start", "end"]].dropna(subset=["annotation"]).copy()
        if tmp.empty:
            continue
        tmp["gene"] = tmp["annotation"].apply(parse_gene_or_annotation)
        tmp = tmp.dropna(subset=["gene"])
        grp = tmp.groupby("gene", as_index=False).agg({"chr": "first", "start": "min", "end": "max"})
        grp["key"] = grp["chr"].apply(_chr_order)
        gene_frames.append(grp)

    if gene_frames:
        gdf = (pd.concat(gene_frames, ignore_index=True)
                 .drop_duplicates("gene")
                 .sort_values(["key", "start", "end"]))
        gene_order = gdf["gene"].tolist()

        plt.figure(figsize=(18, 6))
        ax = plt.gca()
        cohort_handles = []
        for cohort, df in cohort_dfs.items():
            s_cols = [c for c in df.columns if c not in meta]
            med_amp = df[s_cols].median(axis=1)
            tmp = pd.DataFrame({"annotation": df["annotation"], "med": med_amp}).copy()
            tmp["gene"] = tmp["annotation"].apply(parse_gene_or_annotation)
            tmp = tmp.dropna(subset=["gene"])
            gmed = tmp.groupby("gene")["med"].median().reindex(gene_order)
            line, = ax.plot(np.arange(len(gene_order)), gmed.values,
                            lw=2, label=cohort, color=PALETTE[cohort])
            cohort_handles.append(line)

        h_mean  = ax.axhline(mean_line,  color="gold", ls=":",  lw=2,
                             label=f"Mean Coverage Target >= {mean_line}x")
        h_floor = ax.axhline(floor_line, color="red",  ls="--", lw=2,
                             label=f"Minimum Floor {floor_line}x")

        n    = len(gene_order)
        step = max(1, n // 40)
        sel  = list(range(0, n, step))
        ax.set_xticks(sel)
        ax.set_xticklabels([gene_order[i] for i in sel], rotation=45, ha='right', fontsize=7)
        ax.set_yscale("log")
        ax.set_xlabel("Annotated Region / Gene (Genomic Order)", fontweight='bold')
        ax.set_ylabel("Median Coverage Depth (x)", fontweight='bold')
        ax.set_title("Panel Coverage Landscape — By Annotation", fontweight='bold', fontsize=14)
        ax.grid(alpha=0.3)
        handles = cohort_handles + [h_mean, h_floor]
        ax.legend(handles, [h.get_label() for h in handles],
                  bbox_to_anchor=(1.02, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/12_landscape_by_annotation.png", dpi=300)
        plt.close()
    else:
        print("[INFO] No gene annotations found — skipping plot 12")


def plot_coverage_gaps(cohort_dfs: dict, failing_samples: list, plots_dir: str):
    """
    Plot 13: Box plot showing the number of amplicons below the minimum depth
    floor per sample, grouped by cohort. Outlier points corresponding to
    failing samples are highlighted in red to contextualise the gap counts
    within the overall QC framework.
    """
    print("[*] Generating coverage gaps plot...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    meta = {"chr", "start", "end", "id", "annotation", "annot_id"}
    gap_stats = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        for s in s_cols:
            vals = df[s]
            gaps = int(((vals < QC_LIMITS['worst_amplicon_floor']) & (~vals.isna())).sum())
            gap_stats.append({"cohort": cohort, "sample": s, "gaps": gaps})

    gap_df = pd.DataFrame(gap_stats)
    plt.figure(figsize=(10, 7))
    ax = sns.boxplot(data=gap_df, x="cohort", y="gaps",
                     hue="cohort", palette=PALETTE, legend=False)

    # Identify outlier markers rendered by seaborn and recolour those belonging
    # to failing samples; seaborn draws 6 structural lines per box before outliers
    cohort_levels   = list(gap_df["cohort"].unique())
    n_cohorts       = len(cohort_levels)
    structural_lines = n_cohorts * 6

    for line in ax.lines[structural_lines:]:
        if line.get_marker() not in ["o", "s", "D", "^", "v"]:
            continue
        xs, ys = line.get_xdata(), line.get_ydata()
        if len(xs) != 1:
            continue
        x, y         = float(xs[0]), float(ys[0])
        nearest_x    = int(round(x))
        if nearest_x < 0 or nearest_x >= n_cohorts:
            continue
        this_cohort  = cohort_levels[nearest_x]
        matches      = gap_df[(gap_df["cohort"] == this_cohort) & (gap_df["gaps"] == int(y))]
        is_fail      = any(normalise_sample(r["sample"]) in failing_samples
                           for _, r in matches.iterrows())
        if is_fail:
            line.set_color("red")
            line.set_markeredgecolor("black")
            line.set_markersize(7)

    plt.title(f"Coverage Gaps (<{QC_LIMITS['worst_amplicon_floor']}x) per Sample",
              fontweight='bold', fontsize=14)
    plt.xlabel("Cohort", fontweight='bold')
    plt.ylabel(f"Number of Amplicons <{QC_LIMITS['worst_amplicon_floor']}x", fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/13_coverage_gaps.png", dpi=300)
    plt.close()


def plot_retention_and_failures(full_qc_df: pd.DataFrame, plots_dir: str):
    """
    Plots 14–15: Sample retention and failure reason breakdown.
    14: stacked bar chart of PASS/FAIL counts per cohort.
    15: stacked bar chart disaggregating failure reasons per cohort, useful
        for identifying whether failures are systematic (e.g. consistently
        low on-target rate) or random.
    """
    if full_qc_df.empty:
        print("[WARN] No QC data for retention plots"); return

    print("[*] Generating retention and failure plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    # Plot 14: retention rate
    plt.figure(figsize=(10, 7))
    full_qc_df.groupby(["cohort", "status"]).size().unstack().fillna(0).plot(
        kind="bar", stacked=True, color=STATUS_PALETTE, ax=plt.gca()
    )
    plt.title("Sample Retention Rate by Cohort", fontweight='bold', fontsize=14)
    plt.xlabel("Cohort", fontweight='bold')
    plt.ylabel("Number of Samples", fontweight='bold')
    plt.legend(title="QC Status", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/14_retention_rate.png", dpi=300)
    plt.close()

    # Plot 15: failure reasons — each failing sample may contribute more than one reason
    df_rr = full_qc_df.copy()
    df_rr["fail_reasons"] = df_rr["fail_reasons"].fillna("")
    df_rr["fail_list"] = df_rr["fail_reasons"].apply(
        lambda s: [x.strip() for x in str(s).split(",") if x.strip()]
    )
    failed = df_rr[df_rr["status"] == "Fail"].explode("fail_list")
    failed["fail_list"] = failed["fail_list"].replace("", np.nan).fillna("unspecified")

    if not failed.empty:
        reason_map = {
            "low_reads":       "Low Read Count",
            "low_uniformity":  "Low Uniformity (<95% @ 100x)",
            "low_on_target":   "Low On-Target Rate",
            "low_mapping":     "Low Mapping %",
            "low_mean_cov":    "Low Mean Coverage",
            "unspecified":     "Unspecified Error"
        }
        failed["fail_list"] = failed["fail_list"].map(lambda x: reason_map.get(x, x))
        counts = (failed.groupby(["cohort", "fail_list"])
                        .size().reset_index(name="n")
                        .pivot(index="cohort", columns="fail_list", values="n")
                        .fillna(0))
        # Sort columns by total frequency so the most common failure is plotted first
        counts = counts[counts.sum(axis=0).sort_values(ascending=False).index]

        plt.figure(figsize=(12, 7))
        counts.plot(kind="bar", stacked=True, ax=plt.gca(), cmap="tab20")
        plt.legend(title="Failure Reason", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title("QC Failure Reasons by Cohort", fontweight='bold', fontsize=14)
        plt.xlabel("Cohort", fontsize=10, fontweight='bold')
        plt.ylabel("Number of Failed Samples", fontsize=10, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/15_failure_reasons.png", dpi=300)
        plt.close()


def plot_qc_audits(full_qc_df: pd.DataFrame, plots_dir: str):
    """
    Plots 16–17: Strip plots for mapping efficiency and on-target specificity.
    Each point is one sample; colour indicates PASS/FAIL status. Strip plots
    are preferred over box plots here because sample counts per cohort are
    small enough that showing individual points is more informative.
    """
    if full_qc_df.empty:
        print("[WARN] No QC data for audit plots"); return

    print("[*] Generating QC audit plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    # Plot 16: mapping efficiency
    plt.figure(figsize=(10, 7))
    sns.stripplot(data=full_qc_df, x="cohort", y="mapped_pct",
                  hue="status", palette=STATUS_PALETTE, jitter=True, size=6)
    plt.axhline(QC_LIMITS["mapped_pct"], color="red", ls="--", lw=2,
                label=f"Threshold: {QC_LIMITS['mapped_pct']}%")
    plt.title("Mapping Efficiency Audit", fontweight='bold', fontsize=14)
    plt.xlabel("Cohort", fontsize=10, fontweight='bold')
    plt.ylabel("Mapped Reads (%)", fontsize=10, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/16_mapping_efficiency.png", dpi=300)
    plt.close()

    # Plot 17: on-target specificity
    plt.figure(figsize=(10, 7))
    sns.stripplot(data=full_qc_df, x="cohort", y="on_target_pct",
                  hue="status", palette=STATUS_PALETTE, jitter=True, size=6)
    plt.axhline(QC_LIMITS["on_target_pct"], color="firebrick", ls="--", lw=2,
                label=f"Threshold: {QC_LIMITS['on_target_pct']}%")
    plt.title("On-Target Specificity Audit", fontweight='bold', fontsize=14)
    plt.xlabel("Cohort", fontsize=10, fontweight='bold')
    plt.ylabel("On-Target Reads (%)", fontsize=10, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/17_specificity_audit.png", dpi=300)
    plt.close()


def plot_worst_amplicon_analysis(cohort_dfs: dict, plots_dir: str):
    """
    Plot 18: Swarm plot of the minimum depth observed in any amplicon per sample.
    The worst amplicon sets a hard lower bound on variant calling reliability;
    samples whose minimum amplicon depth falls below the floor threshold may
    have missed variants in that region. Symmetric log scale is used because
    values cluster near zero (true zeros are common) but can reach thousands.
    """
    print("[*] Generating worst amplicon plot...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    meta  = {"chr", "start", "end", "id", "annotation", "annot_id"}
    worst = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        for s in s_cols:
            md = df[s].min()
            if pd.notna(md):
                worst.append({"cohort": cohort, "min_depth": md})

    worst_df = pd.DataFrame(worst)
    plt.figure(figsize=(10, 7))
    sns.swarmplot(data=worst_df, x="cohort", y="min_depth",
                  hue="cohort", palette=PALETTE, size=6)
    plt.axhline(QC_LIMITS["worst_amplicon_floor"], color="red", ls="--", lw=2,
                label=f"Minimum Threshold: {QC_LIMITS['worst_amplicon_floor']}x")
    # symlog scale linearises the region near zero (linthresh) while compressing
    # the upper range, accommodating both zero-depth and high-depth samples
    plt.yscale("symlog", linthresh=10)
    plt.title("Worst Amplicon Coverage per Sample", fontweight='bold', fontsize=14)
    plt.xlabel("Cohort", fontsize=10, fontweight='bold')
    plt.ylabel("Minimum Amplicon Depth per Sample (x)", fontsize=10, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/18_worst_amplicon_per_sample.png", dpi=300)
    plt.close()


def plot_systemic_amplicon_failures(cohort_dfs: dict, plots_dir: str):
    """
    Plots 19–20: Amplicons that fail the depth floor in multiple samples.
    Systemic failures (same amplicon failing across many samples) indicate
    a panel design problem rather than a sample quality issue.
    Plot 19 uses rsID / annotation labels; plot 20 uses numeric amplicon IDs,
    providing a cross-reference for panel redesign decisions.
    """
    if not cohort_dfs:
        print("[WARN] No cohort data for failure plots"); return

    print("[*] Generating systemic amplicon failure plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    meta    = {"chr", "start", "end", "id", "annotation", "annot_id"}
    floor   = QC_LIMITS['worst_amplicon_floor']
    cohorts = list(cohort_dfs.keys())
    nrows   = len(cohorts)
    heights = [6] + [4.0] * (nrows - 1) if nrows > 1 else [6]

    def _failure_bar_plot(label_col, filename, ylabel):
        """Shared helper for plots 19 and 20 — differs only in the label column used."""
        fig, axes = plt.subplots(nrows=nrows, ncols=1,
                                 figsize=(16, sum(heights)),
                                 gridspec_kw={'height_ratios': heights},
                                 sharex=True)
        if nrows == 1:
            axes = [axes]
        max_x = 0
        for ax, cohort in zip(axes, cohorts):
            df     = cohort_dfs[cohort]
            s_cols = [c for c in df.columns if c not in meta]
            if not s_cols:
                ax.set_axis_off(); ax.set_title(f"{cohort} — No Data"); continue

            df_copy = df.copy()
            if label_col == "display_label":
                df_copy["display_label"] = df_copy.apply(create_display_label, axis=1)

            id_col    = label_col if label_col in df_copy.columns else "id"
            fail_mask = (df[s_cols] < floor) & (~df[s_cols].isna())
            tmp = pd.DataFrame({
                "label":      df_copy[id_col].values,
                "fail_count": fail_mask.sum(axis=1).values
            }).dropna(subset=["label"])

            agg = (tmp.groupby("label", as_index=False)["fail_count"]
                      .sum().sort_values("fail_count", ascending=False))
            agg = agg[agg["fail_count"] > 0].head(30)

            if agg.empty:
                ax.text(0.5, 0.5, f"No amplicons below {floor}x", ha="center", va="center")
                ax.set_axis_off()
                ax.set_title(f"{cohort} — Systemic Amplicon Failures (<{floor}x)", fontweight='bold')
                continue

            colours = sns.color_palette("Reds", n_colors=len(agg))
            plot_df = agg.iloc[::-1]   # reverse so highest bar is at the top
            ax.barh(plot_df["label"], plot_df["fail_count"],
                    color=colours[::-1], edgecolor="none")
            ax.set_title(f"{cohort} — Systemic Amplicon Failures (<{floor}x)", fontweight='bold')
            ax.set_ylabel(ylabel, fontsize=10, fontweight='bold')
            ax.grid(axis="x", linestyle=":", alpha=0.4)
            max_x = max(max_x, int(plot_df["fail_count"].max()))

        axes[-1].set_xlabel("Number of Failed Samples", fontsize=10, fontweight='bold')
        for ax in axes:
            ax.set_xlim(0, max_x + 1)
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/{filename}", dpi=300)
        plt.close()

    _failure_bar_plot("display_label", "19_systemic_failures_by_rsid.png",    "rsID / Annotation")
    _failure_bar_plot("annot_id",      "20_systemic_failures_by_amplicon_id.png", "Amplicon ID")


def plot_coverage_summary(all_cov: pd.DataFrame, plots_dir: str):
    """
    Plot 21: Histogram with KDE of coverage depth across all amplicons and samples.
    Provides a single-panel summary of the global depth distribution, with
    reference lines marking the mean coverage target and the minimum floor.
    """
    if all_cov.empty:
        print("[WARN] No coverage data for summary plot"); return

    print("[*] Generating coverage summary plot...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    all_depths = []
    for col in all_cov.columns:
        all_depths.extend(all_cov[col].dropna().tolist())

    if not all_depths:
        print("[WARN] No valid coverage values for summary"); return

    plt.figure(figsize=(12, 6))
    sns.histplot(all_depths, bins=50, kde=True, color='steelblue')
    plt.axvline(QC_LIMITS['mean_cov'], color='red', linestyle='--', lw=2,
                label=f"Target Mean: {QC_LIMITS['mean_cov']}x")
    plt.axvline(QC_LIMITS['worst_amplicon_floor'], color='orange', linestyle='--', lw=2,
                label=f"Minimum Floor: {QC_LIMITS['worst_amplicon_floor']}x")
    plt.xlabel('Coverage Depth (x)', fontweight='bold')
    plt.ylabel('Frequency', fontweight='bold')
    plt.title('Global Coverage Distribution — All Samples', fontweight='bold', fontsize=14)
    plt.legend()
    plt.xlim(0, min(1000, max(all_depths)))   # cap x-axis at 1000x to avoid extreme outliers
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/21_coverage_summary.png", dpi=300)
    plt.close()


# =============================================================================
# Entry point
# =============================================================================

def run_complete_pipeline():
    """Load all data and execute all 21 visualisation functions in sequence."""
    print("\n" + "=" * 70)
    print("COMPLETE QC VISUALISATION PIPELINE")
    print("=" * 70 + "\n")

    full_qc_df, cohort_dfs, annotation_df, all_cov, failing_samples = build_pipeline()

    # Export per-sample coverage tables (one compressed TSV per sample per cohort)
    print("\n[*] Exporting per-sample coverage files...")
    meta = {"chr", "start", "end", "id", "annotation", "annot_id"}
    for cohort, df in cohort_dfs.items():
        s_cols    = [c for c in df.columns if c not in meta]
        coords    = df[["id", "chr", "start", "end"]].drop_duplicates().set_index("id")
        cohort_dir = os.path.join(PATHS["samples_dir"], cohort.replace(" ", "_"))
        os.makedirs(cohort_dir, exist_ok=True)
        for s in s_cols:
            out = coords.join(df.set_index("id")[s].rename("depth"))
            out.to_csv(os.path.join(cohort_dir, f"{s}.coverage.tsv.gz"),
                       sep="\t", compression="gzip")

    print("\n" + "=" * 70)
    print("GENERATING ALL VISUALISATIONS")
    print("=" * 70 + "\n")

    plot_qc_distributions(full_qc_df, PATHS["plots_dir"])                          # 01–05
    plot_pca_analysis(full_qc_df, PATHS["plots_dir"])                              # 06–07
    plot_coverage_heatmaps(cohort_dfs, all_cov, failing_samples, PATHS["plots_dir"])  # 08–09
    plot_panel_landscapes(cohort_dfs, PATHS["plots_dir"])                          # 10–12
    plot_coverage_gaps(cohort_dfs, failing_samples, PATHS["plots_dir"])            # 13
    plot_retention_and_failures(full_qc_df, PATHS["plots_dir"])                    # 14–15
    plot_qc_audits(full_qc_df, PATHS["plots_dir"])                                 # 16–17
    plot_worst_amplicon_analysis(cohort_dfs, PATHS["plots_dir"])                   # 18
    plot_systemic_amplicon_failures(cohort_dfs, PATHS["plots_dir"])                # 19–20
    plot_coverage_summary(all_cov, PATHS["plots_dir"])                             # 21

    print("\n" + "=" * 70)
    print("Visualisation pipeline complete.")
    print(f"  Plots:              {PATHS['plots_dir']}")
    print(f"  Per-sample coverage:{PATHS['samples_dir']}")
    print(f"  Full panel CSVs:    {PATHS['out_dir']}")
    print(f"  Total plots:        21")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    run_complete_pipeline()
