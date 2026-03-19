#!/usr/bin/env python3
"""
Complete QC Visualization Pipeline for Ion Torrent Amplicon Sequencing
=======================================================================
Comprehensive quality control visualizations from mosdepth coverage
and QC summary data across multiple cohorts.

This script combines all visualizations from previous versions:
- QC metric distributions (violin + box plots)
- PCA analysis (all samples and pass-only)
- Coverage heatmaps (per-cohort and global)
- Panel landscapes (index, genomic, gene-based)
- Amplicon failure analysis
- Sample retention and QC audit plots

CONFIGURATION:
--------------
To use a specific BED file for rsID annotations, edit the PATHS dictionary below:
    PATHS = {
        ...
        "bed_file": "/path/to/your/IAD255368_167_Submitted.bed"  # <-- Set this!
    }

Or set to None to use the default location (dna_qc/targets.clean.bed):
    "bed_file": None

Date: February 2026
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

# ==============================================================================
# 1. CONFIGURATION
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
    "out_dir": "./analysis_results/03_qc_visualisation",
    "samples_dir": "./analysis_results/03_qc_visualisation/sample_level_coverage_tables",
    "plots_dir": "./analysis_results/03_qc_visualisation/figures",
    "raw_root": "./",
    # ========================================================================
    # SET YOUR BED FILE PATH HERE for rsID annotations!
    # Example: "bed_file": "/path/to/IAD255368_167_Submitted.bed"
    # Or leave as None to use default: cohort_dir/dna_qc/targets.sorted.bed
    # ========================================================================
    "bed_file": None  
}

COHORT_MAP = {
    'Breast Tumour': 'breast/tumour',
    'Breast Control': 'breast/normal',
    'Endometrium Tumour': 'endometrium/tumour',
    'Endometrium Control': 'endometrium/normal',
}

PALETTE = {
    'Breast Tumour': '#3498db',
    'Breast Control': "#d1d946",
    'Endometrium Tumour': '#9b59b6',
    'Endometrium Control': '#2ecc71'
}

STATUS_PALETTE = {
    'Pass': '#2ecc71',
    'Fail': '#e74c3c'
}

# With clean mosdepth outputs, keep strict
STRICT_MOSDEPTH = True
MIN_NUMERIC_FRACTION = 0.90


# ==============================================================================
# 2. UTILITIES
# ==============================================================================
def normalise_sample(s: str) -> str:
    """Standardise sample names."""
    s = str(s).strip().lower()
    return (
        s.replace(" ", "_")
         .replace(".regions", "")
         .replace(".bed", "")
         .replace(".gz", "")
    )


def normalise_chr(c: str) -> str:
    """Ensure chromosome naming is consistent (chr17 etc.)."""
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
    """Choose appropriate depth column if BED+mosdepth output has >4 cols."""
    ncols = df.shape[1]
    if ncols == 4:
        return 3
    if ncols < 4:
        return -1
    last = df.columns[-1]
    frac = pd.to_numeric(df[last], errors='coerce').notna().mean()
    if frac > 0.8:
        return ncols - 1
    best, best_frac = -1, -1
    for i in range(3, ncols):
        f = pd.to_numeric(df.iloc[:, i], errors='coerce').notna().mean()
        if f > best_frac:
            best = i
            best_frac = f
    return best if best_frac >= 0.5 else -1


def create_display_label(row: pd.Series) -> str:
    """
    Create a clean display label from annotation.
    Since load_annotation_for_cohort() already prioritizes rsID > gene > coordinate,
    we just need to return the annotation or fall back to coordinate.
    """
    annot = row.get("annotation", "")
    
    # If we have a valid annotation, use it
    if pd.notna(annot) and str(annot).strip() not in ["", "."]:
        return str(annot).strip()
    
    # Fallback to coordinate
    return f"{row['chr']}:{row['start']}-{row['end']}"


# ==============================================================================
# 3. LOADING: annotation, per-file coverage, per-cohort coverage, QC
# ==============================================================================
def load_annotation_for_cohort(root_path: str, bed_file_override: str = None) -> pd.DataFrame:
    """
    Load per-cohort annotation BED and extract rsIDs properly.
    
    Args:
        root_path: Path to cohort directory (used if bed_file_override is None)
        bed_file_override: Direct path to BED file (takes priority if provided)
    """
    if bed_file_override and os.path.exists(bed_file_override):
        annot_file = bed_file_override
        print(f"[*] Using specified BED file: {annot_file}")
    else:
        annot_file = os.path.join(root_path, "dna_qc", "targets.sorted.bed")
        if not os.path.exists(annot_file):
            raise FileNotFoundError(f"Annotation file not found: {annot_file}")

    # Read BED file, skip track line if present
    with open(annot_file, 'r') as f:
        first_line = f.readline()
    
    skip_rows = 1 if first_line.startswith('track') else 0
    
    # Read BED - format: chr, start, end, amplicon_id, rsid, gene
    ann = pd.read_csv(annot_file, sep="\t", header=None, skiprows=skip_rows)
    
    # Handle different BED formats
    if ann.shape[1] >= 6:
        # Full format: chr, start, end, amplicon_id, rsid, gene
        ann.columns = ["chr", "start", "end", "amplicon_id", "rsid", "gene"] + \
                      [f"extra_{i}" for i in range(ann.shape[1] - 6)]
        
        # Create annotation field: prioritize rsID, then gene, then amplicon_id
        def make_annotation(row):
            rsid = row["rsid"]
            gene = row["gene"]
            amp_id = row["amplicon_id"]
            
            # Priority 1: rsID (if not empty or just ".")
            if pd.notna(rsid) and str(rsid).strip() not in ["", "."]:
                return str(rsid).strip()
            
            # Priority 2: Gene name (if not empty or just ".")
            if pd.notna(gene) and str(gene).strip() not in ["", "."]:
                return str(gene).strip()
            
            # Priority 3: Amplicon ID as fallback
            return str(amp_id)
        
        ann["annotation"] = ann.apply(make_annotation, axis=1)
        
    elif ann.shape[1] >= 5:
        # 5 columns: chr, start, end, amplicon_id, annotation
        ann.columns = ["chr", "start", "end", "amplicon_id", "annotation"] + \
                      [f"extra_{i}" for i in range(ann.shape[1] - 5)]
    else:
        # Basic 4-column BED
        ann.columns = ["chr", "start", "end", "annotation"]
    
    ann["chr"] = ann["chr"].astype(str)
    ann["annot_id"] = ann["annotation"]
    
    return ann[["chr", "start", "end", "annotation", "annot_id"]].drop_duplicates(["chr", "start", "end"])


def load_coverage_for_file(path: str, sample_name: str) -> pd.DataFrame:
    """Load one mosdepth .regions.bed.gz file and normalise."""
    raw = pd.read_csv(path, sep="\t", header=None, compression="gzip", dtype=str)
    if raw.empty:
        if STRICT_MOSDEPTH:
            raise ValueError(f"{path} is empty")
        return pd.DataFrame(columns=["chr", "start", "end", sample_name]).set_index([])

    raw.iloc[:, 0] = raw.iloc[:, 0].map(normalise_chr)
    raw.iloc[:, 1] = pd.to_numeric(raw.iloc[:, 1], errors='coerce')
    raw.iloc[:, 2] = pd.to_numeric(raw.iloc[:, 2], errors='coerce')

    depth_col = pick_depth_column(raw)
    depth = (
        pd.to_numeric(raw.iloc[:, depth_col], errors='coerce')
        if depth_col != -1 else pd.Series(np.nan, index=raw.index)
    )

    df = pd.DataFrame({
        "chr": raw.iloc[:, 0].astype(str),
        "start": raw.iloc[:, 1].astype("Int64"),
        "end": raw.iloc[:, 2].astype("Int64"),
        sample_name: depth
    })
    df["id"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    df = df.drop_duplicates("id").set_index("id")

    frac = df[sample_name].notna().mean()
    if frac < MIN_NUMERIC_FRACTION and STRICT_MOSDEPTH:
        raise ValueError(f"{sample_name}: only {frac:.1%} numeric")
    return df


def load_coverage_for_cohort(cohort_name: str, root_path: str) -> pd.DataFrame:
    """Load all mosdepth outputs for a cohort and join by ID."""
    pattern = os.path.join(root_path, "**", "*.regions.bed.gz")
    files = [
        f for f in glob.glob(pattern, recursive=True)
        if f.endswith(".regions.bed.gz") and ".bai" not in f.lower()
    ]
    if not files:
        return pd.DataFrame()

    frames = []
    for f in sorted(files):
        name = normalise_sample(os.path.basename(f).replace(".regions.bed.gz", ""))
        try:
            frames.append(load_coverage_for_file(f, name))
        except Exception as e:
            raise RuntimeError(f"[MOSDEPTH ERROR] {cohort_name} {name}: {e}")

    combined = frames[0].copy()
    for frame in frames[1:]:
        cols = [c for c in frame.columns if c not in ["chr", "start", "end"]]
        combined = combined.join(frame[cols], how="outer")

    return combined.reset_index()


def load_qc_for_cohort(cohort_name: str, root_path: str) -> pd.DataFrame:
    """Load QC summary for one cohort."""
    qc_file = os.path.join(root_path, "dna_qc", "qc_summary.tsv")
    if not os.path.exists(qc_file):
        return pd.DataFrame()
    df = pd.read_csv(qc_file, sep="\t")
    df["cohort"] = cohort_name
    df["status"] = df["status"].str.capitalize()
    df["sample_norm"] = df["sample"].apply(normalise_sample)
    return df


# ==============================================================================
# 4. Pipeline Assembly
# ==============================================================================
def build_pipeline():
    """Load everything: QC, coverage, annotation; produce merged tables."""
    
    # DEBUG: Check BED file configuration
    print("="*80)
    print("CONFIGURATION CHECK:")
    print("="*80)
    print(f"BED file in PATHS: {PATHS.get('bed_file')}")
    if PATHS.get('bed_file'):
        print(f"BED file exists: {os.path.exists(PATHS['bed_file'])}")
        if os.path.exists(PATHS['bed_file']):
            print(f"BED file path (absolute): {os.path.abspath(PATHS['bed_file'])}")
    print("="*80)
    print()
    
    # Create output directories (skip bed_file as it's a file path, not a directory)
    for key, p in PATHS.items():
        if key != "bed_file" and p and isinstance(p, str):
            os.makedirs(p, exist_ok=True)

    qc_frames = []
    cohort_dfs = {}
    annotation_frames = []
    
    # If a BED file is specified, load it once and use for all cohorts
    global_bed = None
    if PATHS.get("bed_file") and os.path.exists(PATHS["bed_file"]):
        print(f"[*] Loading global BED file: {PATHS['bed_file']}")
        global_bed = load_annotation_for_cohort("", bed_file_override=PATHS["bed_file"])
        annotation_frames.append(global_bed)
        
        # DEBUG: Show what annotations were loaded
        print(f"[DEBUG] Loaded {len(global_bed)} annotations from BED file")
        print(f"[DEBUG] First 10 annotations:")
        print(global_bed[["chr", "start", "end", "annotation"]].head(10))
        print(f"[DEBUG] Annotation types:")
        rsid_count = global_bed["annotation"].str.startswith("rs", na=False).sum()
        gene_count = (global_bed["annotation"].str.match(r'^[A-Z]+$', na=False)).sum()
        other_count = len(global_bed) - rsid_count - gene_count
        print(f"  rsIDs: {rsid_count} ({rsid_count/len(global_bed)*100:.1f}%)")
        print(f"  Gene names: {gene_count} ({gene_count/len(global_bed)*100:.1f}%)")
        print(f"  Other: {other_count} ({other_count/len(global_bed)*100:.1f}%)")
        print()

    print("[*] Loading data...")
    for cohort, relpath in COHORT_MAP.items():
        root = os.path.join(PATHS["raw_root"], relpath)
        if not os.path.exists(root):
            print(f"[WARN] Missing path {root}")
            continue

        qc = load_qc_for_cohort(cohort, root)
        if not qc.empty:
            qc_frames.append(qc)

        cov = load_coverage_for_cohort(cohort, root)
        if cov.empty:
            print(f"[WARN] No coverage for {cohort}")
            continue

        # Use global BED if specified, otherwise load per-cohort
        if global_bed is not None:
            ann = global_bed.copy()
        else:
            ann = load_annotation_for_cohort(root)
            annotation_frames.append(ann)

        cov = cov.merge(
            ann[["chr", "start", "end", "annotation", "annot_id"]],
            on=["chr", "start", "end"],
            how="left"
        )
        
        # DEBUG: Check what annotations ended up in the coverage data
        print(f"[DEBUG] {cohort}: After merge, coverage data has {len(cov)} rows")
        print(f"[DEBUG] {cohort}: Sample annotations in merged data:")
        print(cov[["chr", "start", "end", "annotation"]].head(10))
        print()
        
        cohort_dfs[cohort] = cov

        out_full = os.path.join(
            PATHS["out_dir"], f"GSDMB_{cohort.replace(' ', '_')}_FullPanel.csv"
        )
        cov.to_csv(out_full, index=False)

        s_cols = [c for c in cov.columns
                  if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]
        rows = []
        for s in s_cols:
            depth = cov[s]
            tot = len(depth)
            num = depth.notna().sum()
            zeros = (depth == 0).sum()
            rows.append({
                "sample": s,
                "amplicons_total": tot,
                "numeric_depth": num,
                "zeros": zeros,
                "numeric_fraction": num / tot if tot else 0
            })
        audit = pd.DataFrame(rows)
        audit.to_csv(
            os.path.join(PATHS["out_dir"], f"CoverageAudit_{cohort.replace(' ', '_')}.csv"),
            index=False
        )

    full_qc = pd.concat(qc_frames, ignore_index=True) if qc_frames else pd.DataFrame()
    annotation_df = (
        pd.concat(annotation_frames, ignore_index=True)
          .drop_duplicates(["chr", "start", "end"])
        if annotation_frames else
        pd.DataFrame(columns=["chr", "start", "end", "annotation", "annot_id"])
    )

    print("[*] Building global coverage...")
    all_cov_frames = []
    for cohort, df in cohort_dfs.items():
        tmp = df.set_index("id").drop(
            columns=["chr", "start", "end", "annotation", "annot_id"],
            errors="ignore"
        )
        all_cov_frames.append(tmp)
    all_cov = pd.concat(all_cov_frames, axis=1) if all_cov_frames else pd.DataFrame()

    failing_samples = (
        full_qc[full_qc["status"] == "Fail"]["sample_norm"].tolist()
        if not full_qc.empty else []
    )

    print("[*] Data loading complete.")
    return full_qc, cohort_dfs, annotation_df, all_cov, failing_samples


# ==============================================================================
# 5. PLOTTING FUNCTIONS
# ==============================================================================
def plot_qc_distributions(full_qc_df: pd.DataFrame, plots_dir: str):
    """
    Plot 01-05: QC metric distributions (violin + box plots)
    One plot per metric showing cohort distribution and pass/fail comparison.
    """
    if full_qc_df.empty:
        print("[WARN] No QC data for distribution plots")
        return
    
    print("[*] Generating QC distribution plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300
    
    metrics = [
        ('total_reads', 'Total Reads', QC_LIMITS['total_reads'], '01'),
        ('mapped_pct', 'Mapped %', QC_LIMITS['mapped_pct'], '02'),
        ('on_target_pct', 'On Target %', QC_LIMITS['on_target_pct'], '03'),
        ('mean_cov', 'Mean Coverage', QC_LIMITS['mean_cov'], '04'),
        ('uniformity_100x', 'Uniformity >100x %', QC_LIMITS['uniformity_100x'], '05')
    ]
    
    for metric, label, threshold, plot_num in metrics:
        if metric not in full_qc_df.columns:
            print(f"[WARN] Metric {metric} not found")
            continue
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Left: Violin plot by cohort
        sns.violinplot(
            data=full_qc_df,
            x='cohort',
            y=metric,
            hue='cohort',
            palette=PALETTE,
            ax=axes[0],
            legend=False
        )
        axes[0].axhline(threshold, color='red', linestyle='--', 
                       label=f'Threshold: {threshold}', linewidth=2)
        axes[0].set_title(f'{label} by Cohort', fontweight='bold')
        axes[0].set_xlabel('Cohort', fontweight='bold')
        axes[0].set_ylabel(label, fontweight='bold')
        axes[0].legend()
        axes[0].tick_params(axis='x', rotation=45)
        
        # Right: Box plot by status
        if 'status' in full_qc_df.columns:
            sns.boxplot(
                data=full_qc_df,
                x='status',
                y=metric,
                hue='status',
                palette=STATUS_PALETTE,
                ax=axes[1],
                legend=False
            )
            axes[1].axhline(threshold, color='red', linestyle='--', 
                          label=f'Threshold: {threshold}', linewidth=2)
            axes[1].set_title(f'{label} by QC Status', fontweight='bold')
            axes[1].set_xlabel('QC Status', fontweight='bold')
            axes[1].set_ylabel(label, fontweight='bold')
            axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/{plot_num}_qc_{metric}_distribution.png", dpi=300)
        plt.close()


def plot_pca_analysis(full_qc_df: pd.DataFrame, plots_dir: str):
    """
    Plot 06-07: PCA analysis of QC metrics
    06: All samples (colored by cohort and styled by status)
    07: Pass-only samples (colored by cohort)
    """
    if full_qc_df.empty:
        print("[WARN] No QC data for PCA")
        return
    
    print("[*] Generating PCA plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300
    
    full_qc_df = full_qc_df.copy()
    full_qc_df["uniformity_gap"] = 100.0 - full_qc_df["uniformity_100x"]
    feats = ['mean_cov', 'on_target_pct', 'total_reads', 'mapped_pct', 'uniformity_gap']
    
    # Plot 06: All samples PCA
    X = StandardScaler().fit_transform(full_qc_df[feats].fillna(0))
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)
    var1 = pca.explained_variance_ratio_[0] * 100
    var2 = pca.explained_variance_ratio_[1] * 100
    
    full_qc_df["PC1"], full_qc_df["PC2"] = pcs[:, 0], pcs[:, 1]
    
    plt.figure(figsize=(10, 7))
    ax = sns.scatterplot(
        data=full_qc_df, x="PC1", y="PC2",
        hue="cohort", style="status",
        palette=PALETTE, s=100
    )
    ax.set_xlabel(f"PC1 ({var1:.1f}% variance)", fontweight='bold')
    ax.set_ylabel(f"PC2 ({var2:.1f}% variance)", fontweight='bold')
    plt.title("QC PCA: All Samples", fontweight='bold', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/06_pca_all_samples.png", dpi=300)
    plt.close()
    
    # Plot 07: Pass-only PCA
    pass_only = full_qc_df[full_qc_df["status"] == "Pass"].copy()
    if len(pass_only) >= 3:
        X = StandardScaler().fit_transform(pass_only[feats].fillna(0))
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X)
        var1 = pca.explained_variance_ratio_[0] * 100
        var2 = pca.explained_variance_ratio_[1] * 100
        
        pass_only["PC1"], pass_only["PC2"] = pcs[:, 0], pcs[:, 1]
        
        plt.figure(figsize=(10, 7))
        ax = sns.scatterplot(
            data=pass_only, x="PC1", y="PC2",
            hue="cohort", palette=PALETTE, s=100
        )
        ax.set_xlabel(f"PC1 ({var1:.1f}% variance)", fontweight='bold')
        ax.set_ylabel(f"PC2 ({var2:.1f}% variance)", fontweight='bold')
        plt.title("QC PCA: Passing Samples", fontweight='bold', fontsize=14)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/07_pca_pass_samples.png", dpi=300)
        plt.close()


def plot_coverage_heatmaps(cohort_dfs: dict, all_cov: pd.DataFrame, 
                          failing_samples: list, plots_dir: str):
    """
    Plot 08-09: Coverage heatmaps
    08: Per-cohort heatmaps (amplicons × samples)
    09: Global heatmap (all cohorts combined)
    """
    print("[*] Generating coverage heatmaps...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300
    
    meta = {"chr", "start", "end", "id", "annotation", "annot_id"}
    clean_failing_samples = {normalise_sample(s) for s in failing_samples}
    
    # Plot 08: Per-cohort heatmaps
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        
        # Create display labels for amplicons
        df_plot = df.copy()
        df_plot["display_label"] = df_plot.apply(create_display_label, axis=1)
        
        dmat = df_plot.set_index("display_label")[s_cols]
        
        plotmat = np.log10(dmat.fillna(0) + 1)
        cmap = sns.color_palette("magma", as_cmap=True)
        cmap.set_bad(color="#B0B0B0")
        
        # Adjust figure size based on number of amplicons
        n_amplicons = len(dmat)
        fig_height = max(8, min(20, n_amplicons * 0.15))
        
        plt.figure(figsize=(20, fig_height))
        ax = sns.heatmap(plotmat, cmap=cmap, mask=dmat.isna(), vmin=0, vmax=4,
                        cbar_kws={'label': 'log10(Coverage + 1)'})
        
        # Force tick at every column (samples on x-axis)
        ax.xaxis.set_major_locator(mticker.FixedLocator(np.arange(len(s_cols)) + 0.5))
        ax.set_xticklabels(s_cols, rotation=90, ha='center', fontsize=6)
        ax.tick_params(axis='x', which='both', length=0, pad=1)
        
        # Add vertical separators for samples
        for x in range(len(s_cols) + 1):
            ax.axvline(x, color="white", lw=0.3, alpha=0.6)
        
        # Y-axis: show rsID labels (limit to reasonable number if too many)
        n_yticks = min(len(dmat), 50)  # Limit y-axis labels
        if len(dmat) > n_yticks:
            step = len(dmat) // n_yticks
            ytick_positions = list(range(0, len(dmat), step))
            ytick_labels = [dmat.index[i] for i in ytick_positions]
            ax.set_yticks([i + 0.5 for i in ytick_positions])
            ax.set_yticklabels(ytick_labels, fontsize=7)
        else:
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
        
        # Highlight failing samples on x-axis
        for label in ax.get_xticklabels():
            sample_on_plot = normalise_sample(label.get_text())
            if sample_on_plot in clean_failing_samples:
                label.set_color("red")
                label.set_weight("bold")
        
        plt.xlabel("Sample ID", fontsize=10, fontweight='bold')
        plt.ylabel("Amplicon (rsID / Annotation)", fontsize=10, fontweight='bold')
        plt.title(f"{cohort} Coverage Heatmap", fontweight='bold', fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/08_heatmap_{cohort.replace(' ', '_')}.png", dpi=300)
        plt.close()
    
    # Plot 09: Global heatmap
    if not all_cov.empty:
        plt.figure(figsize=(30, 12))
        cmap = sns.color_palette("magma", as_cmap=True)
        cmap.set_bad(color="#B0B0B0")
        
        sns.heatmap(np.log10(all_cov.fillna(0) + 1),
                   cmap=cmap, mask=all_cov.isna(),
                   vmin=0, vmax=4,
                   cbar_kws={'label': 'log10(Coverage + 1)'})
        
        ax = plt.gca()
        plt.xticks(rotation=90, fontsize=5)
        for label in ax.get_xticklabels():
            sample_on_plot = normalise_sample(label.get_text())
            if sample_on_plot in clean_failing_samples:
                label.set_color("red")
                label.set_weight("bold")
        
        plt.xlabel("Sample ID", fontsize=10, fontweight='bold')
        plt.ylabel("Genomic Regions (All Amplicons)", fontsize=10, fontweight='bold')
        plt.title("Global Coverage Heatmap", fontweight='bold', fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/09_heatmap_global.png", dpi=300)
        plt.close()


def plot_panel_landscapes(cohort_dfs: dict, plots_dir: str):
    """
    Plot 10-12: Panel coverage landscapes
    10: By amplicon index order
    11: By genomic coordinates
    12: By annotated gene/region
    """
    print("[*] Generating panel landscape plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300
    
    meta = {"chr", "start", "end", "id", "annotation", "annot_id"}
    mean_line = QC_LIMITS['mean_cov']
    floor_line = QC_LIMITS['worst_amplicon_floor']
    
    # Plot 10: Index order landscape
    plt.figure(figsize=(15, 6))
    ax = plt.gca()
    cohort_handles = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        med = df[s_cols].median(axis=1)
        line, = ax.plot(range(len(med)), med, lw=2, label=cohort, color=PALETTE[cohort])
        cohort_handles.append(line)
    
    h_mean = ax.axhline(mean_line, color="gold", ls=":", lw=2, 
                       label=f"Mean Coverage Target ≥ {mean_line}×")
    h_floor = ax.axhline(floor_line, color="red", ls="--", lw=2, 
                        label=f"Minimum Floor {floor_line}×")
    
    ax.set_yscale("log")
    ax.set_xlabel("Amplicon Index (Order in File)", fontweight='bold')
    ax.set_ylabel("Median Coverage Depth (×)", fontweight='bold')
    ax.set_title("Coverage Landscape: Input Order", fontweight='bold', fontsize=14)
    ax.grid(alpha=0.3)
    
    handles = cohort_handles + [h_mean, h_floor]
    labels = [h.get_label() for h in handles]
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/10_landscape_index_order.png", dpi=300)
    plt.close()
    
    # Plot 11: Genomic coordinate order landscape
    def _chr_order(c):
        c = str(c).replace("chr", "")
        mapping = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        try:
            return int(c)
        except Exception:
            return mapping.get(c, 1000)
    
    coord_frames = []
    for df in cohort_dfs.values():
        tmp = df[["chr", "start", "end", "id", "annotation"]].drop_duplicates()
        tmp["display_label"] = tmp.apply(create_display_label, axis=1)
        coord_frames.append(tmp)
    
    kdf = pd.concat(coord_frames).drop_duplicates(subset=["chr", "start", "end"])
    kdf["key"] = kdf["chr"].apply(_chr_order)
    kdf = kdf.sort_values(["key", "start", "end"])
    ordered_ids = kdf["id"].tolist()
    ordered_labels = kdf["display_label"].tolist()
    
    plt.figure(figsize=(18, 6))
    ax = plt.gca()
    cohort_handles = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        med = df[s_cols].median(axis=1)
        t = (pd.DataFrame({"id": df["id"], "med": med})
              .drop_duplicates("id").set_index("id").reindex(ordered_ids))
        xvals = np.arange(len(ordered_ids))
        line, = ax.plot(xvals, t["med"].values, lw=2, label=cohort, color=PALETTE[cohort])
        cohort_handles.append(line)
    
    h_mean = ax.axhline(mean_line, color="gold", ls=":", lw=2, 
                       label=f"Mean Coverage Target ≥ {mean_line}×")
    h_floor = ax.axhline(floor_line, color="red", ls="--", lw=2, 
                        label=f"Minimum Floor {floor_line}×")
    
    if ordered_labels:
        n = len(ordered_labels)
        step = max(1, n // 40)
        sel = list(range(0, n, step))
        ax.set_xticks(sel)
        ax.set_xticklabels([ordered_labels[i] for i in sel], rotation=45, ha='right', fontsize=7)
    
    ax.set_yscale("log")
    ax.set_xlabel("Genomic Position (rsID / Annotation)", fontweight='bold')
    ax.set_ylabel("Median Coverage Depth (×)", fontweight='bold')
    ax.set_title("Coverage Landscape: Genomic Order", fontweight='bold', fontsize=14)
    ax.grid(alpha=0.3)
    
    handles = cohort_handles + [h_mean, h_floor]
    labels = [h.get_label() for h in handles]
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/11_landscape_genomic_order.png", dpi=300)
    plt.close()
    
    # Plot 12: Gene/annotation landscape
    def parse_gene_or_annotation(a):
        if pd.isna(a):
            return None
        t = str(a)
        if "\n" in t:
            return t.split("\n")[1]
        return t.strip()
    
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
        gdf = pd.concat(gene_frames, ignore_index=True).drop_duplicates("gene")
        gdf = gdf.sort_values(["key", "start", "end"])
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
            
            xvals = np.arange(len(gene_order))
            line, = ax.plot(xvals, gmed.values, lw=2, label=cohort, color=PALETTE[cohort])
            cohort_handles.append(line)
        
        h_mean = ax.axhline(mean_line, color="gold", ls=":", lw=2,
                           label=f"Mean Coverage Target ≥ {mean_line}×")
        h_floor = ax.axhline(floor_line, color="red", ls="--", lw=2,
                            label=f"Minimum Floor {floor_line}×")
        
        n = len(gene_order)
        step = max(1, n // 40)
        sel = list(range(0, n, step))
        ax.set_xticks(sel)
        ax.set_xticklabels([gene_order[i] for i in sel], rotation=45, ha='right', fontsize=7)
        
        ax.set_yscale("log")
        ax.set_xlabel("Annotated Region / Gene (Genomic Order)", fontweight='bold')
        ax.set_ylabel("Median Coverage Depth (×)", fontweight='bold')
        ax.set_title("Coverage Landscape: By Annotation", fontweight='bold', fontsize=14)
        ax.grid(alpha=0.3)
        
        handles = cohort_handles + [h_mean, h_floor]
        labels = [h.get_label() for h in handles]
        ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left')
        
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/12_landscape_by_annotation.png", dpi=300)
        plt.close()
    else:
        print("[INFO] Skipping gene landscape plot: no gene annotations found")


def plot_coverage_gaps(cohort_dfs: dict, failing_samples: list, plots_dir: str):
    """
    Plot 13: Coverage gaps analysis
    Box plot showing number of amplicons below threshold per sample, by cohort.
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
            gaps = (vals < QC_LIMITS['worst_amplicon_floor']) & (~vals.isna())
            gap_stats.append({"cohort": cohort, "sample": s, "gaps": int(gaps.sum())})
    
    gap_df = pd.DataFrame(gap_stats)
    
    plt.figure(figsize=(10, 7))
    ax = sns.boxplot(
        data=gap_df, x="cohort", y="gaps",
        hue="cohort", palette=PALETTE, legend=False
    )
    
    # Highlight failing sample outliers
    cohort_levels = list(gap_df["cohort"].unique())
    n_cohorts = len(cohort_levels)
    structural_lines = n_cohorts * 6
    
    for line in ax.lines[structural_lines:]:
        if line.get_marker() not in ["o", "s", "D", "^", "v"]:
            continue
        xs = line.get_xdata()
        ys = line.get_ydata()
        if len(xs) != 1:
            continue
        x, y = float(xs[0]), float(ys[0])
        nearest_x = int(round(x))
        if nearest_x < 0 or nearest_x >= n_cohorts:
            continue
        this_cohort = cohort_levels[nearest_x]
        matches = gap_df[(gap_df["cohort"] == this_cohort) & (gap_df["gaps"] == int(y))]
        is_fail = any(normalise_sample(row["sample"]) in failing_samples
                     for _, row in matches.iterrows())
        if is_fail:
            line.set_color("red")
            line.set_markeredgecolor("black")
            line.set_markersize(7)
    
    plt.title(f"Coverage Gaps per Sample (<{QC_LIMITS['worst_amplicon_floor']}x)", 
             fontweight='bold', fontsize=14)
    plt.xlabel("Cohort", fontweight='bold')
    plt.ylabel(f"Number of Amplicons <{QC_LIMITS['worst_amplicon_floor']}×", fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/13_coverage_gaps.png", dpi=300)
    plt.close()


def plot_retention_and_failures(full_qc_df: pd.DataFrame, plots_dir: str):
    """
    Plot 14-15: Sample retention and failure analysis
    14: Retention rate (pass/fail by cohort)
    15: Detailed failure reasons by cohort
    """
    if full_qc_df.empty:
        print("[WARN] No QC data for retention plots")
        return
    
    print("[*] Generating retention and failure plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300
    
    # Plot 14: Retention rate
    plt.figure(figsize=(10, 7))
    full_qc_df.groupby(["cohort", "status"]).size().unstack().fillna(0).plot(
        kind="bar", stacked=True, color=STATUS_PALETTE, ax=plt.gca()
    )
    plt.title("QC Retention by Cohort", fontweight='bold', fontsize=14)
    plt.xlabel("Cohort", fontweight='bold')
    plt.ylabel("Number of Samples", fontweight='bold')
    plt.legend(title="QC Status", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/14_retention_rate.png", dpi=300)
    plt.close()
    
    # Plot 15: Failure reasons
    df_rr = full_qc_df.copy()
    df_rr["fail_reasons"] = df_rr["fail_reasons"].fillna("")
    df_rr["fail_list"] = df_rr["fail_reasons"].apply(
        lambda s: [x.strip() for x in str(s).split(",") if x.strip()]
    )
    failed = df_rr[df_rr["status"] == "Fail"].explode("fail_list")
    failed["fail_list"] = failed["fail_list"].replace("", np.nan).fillna("unspecified")
    
    if not failed.empty:
        reason_map = {
            "low_reads": "Low Read Count",
            "low_uniformity": "Low Uniformity (<95% @ 100×)",
            "low_on_target": "Low On‑Target Rate",
            "low_mapping": "Low Mapping %",
            "low_mean_cov": "Low Mean Coverage",
            "unspecified": "Unspecified Error"
        }
        failed["fail_list"] = failed["fail_list"].map(lambda x: reason_map.get(x, x))
        
        counts = (failed.groupby(["cohort", "fail_list"])
                        .size().reset_index(name="n")
                        .pivot(index="cohort", columns="fail_list", values="n")
                        .fillna(0))
        counts = counts[counts.sum(axis=0).sort_values(ascending=False).index]
        
        plt.figure(figsize=(12, 7))
        counts.plot(kind="bar", stacked=True, ax=plt.gca(), cmap="tab20")
        plt.legend(title="Failure Reason", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title("QC Failure Reasons", fontweight='bold', fontsize=14)
        plt.xlabel("Cohort", fontsize=10, fontweight='bold')
        plt.ylabel("Number of Failed Samples", fontsize=10, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f"{plots_dir}/15_failure_reasons.png", dpi=300)
        plt.close()


def plot_qc_audits(full_qc_df: pd.DataFrame, plots_dir: str):
    """
    Plot 16-17: QC metric audits
    16: Mapping efficiency audit
    17: Specificity (on-target) audit
    """
    if full_qc_df.empty:
        print("[WARN] No QC data for audit plots")
        return
    
    print("[*] Generating QC audit plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300
    
    # Plot 16: Mapping efficiency
    plt.figure(figsize=(10, 7))
    sns.stripplot(
        data=full_qc_df, x="cohort", y="mapped_pct",
        hue="status", palette=STATUS_PALETTE, jitter=True, size=6
    )
    plt.axhline(QC_LIMITS["mapped_pct"], color="red", ls="--", lw=2,
               label=f"Threshold: {QC_LIMITS['mapped_pct']}%")
    plt.title("Mapping Efficiency by Cohort", fontweight='bold', fontsize=14)
    plt.xlabel("Cohort", fontsize=10, fontweight='bold')
    plt.ylabel("Mapped Reads (%)", fontsize=10, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/16_mapping_efficiency.png", dpi=300)
    plt.close()
    
    # Plot 17: Specificity (on-target)
    plt.figure(figsize=(10, 7))
    sns.stripplot(
        data=full_qc_df, x="cohort", y="on_target_pct",
        hue="status", palette=STATUS_PALETTE, jitter=True, size=6
    )
    plt.axhline(QC_LIMITS["on_target_pct"], color="firebrick", ls="--", lw=2,
               label=f"Threshold: {QC_LIMITS['on_target_pct']}%")
    plt.title("On-Target Specificity by Cohort", fontweight='bold', fontsize=14)
    plt.xlabel("Cohort", fontsize=10, fontweight='bold')
    plt.ylabel("On-Target Reads (%)", fontsize=10, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/17_specificity_audit.png", dpi=300)
    plt.close()


def plot_worst_amplicon_analysis(cohort_dfs: dict, plots_dir: str):
    """
    Plot 18: Worst amplicon coverage per sample
    Shows the minimum depth observed in any amplicon for each sample.
    """
    print("[*] Generating worst amplicon plot...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300
    
    meta = {"chr", "start", "end", "id", "annotation", "annot_id"}
    worst = []
    
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        for s in s_cols:
            md = df[s].min()
            if pd.notna(md):
                worst.append({"cohort": cohort, "min_depth": md})
    
    worst_df = pd.DataFrame(worst)
    plt.figure(figsize=(10, 7))
    sns.swarmplot(
        data=worst_df, x="cohort", y="min_depth",
        hue="cohort", palette=PALETTE, size=6
    )
    plt.axhline(QC_LIMITS["worst_amplicon_floor"], color="red", ls="--", lw=2,
               label=f"Minimum Threshold: {QC_LIMITS['worst_amplicon_floor']}×")
    plt.yscale("symlog", linthresh=10)
    plt.title("Lowest Amplicon Coverage per Sample", fontweight='bold', fontsize=14)
    plt.xlabel("Cohort", fontsize=10, fontweight='bold')
    plt.ylabel("Minimum Amplicon Depth per Sample (×)", fontsize=10, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/18_worst_amplicon_per_sample.png", dpi=300)
    plt.close()


def plot_systemic_amplicon_failures(cohort_dfs: dict, plots_dir: str):
    """
    Plot 19-20: Systemic amplicon failure analysis
    19: Failures by rsID/annotation
    20: Failures by amplicon ID (numeric)
    """
    if not cohort_dfs:
        print("[WARN] No cohort data for failure plots")
        return
    
    print("[*] Generating systemic amplicon failure plots...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300
    
    meta = {"chr", "start", "end", "id", "annotation", "annot_id"}
    floor = QC_LIMITS['worst_amplicon_floor']
    cohorts = list(cohort_dfs.keys())
    
    # Plot 19: By rsID/annotation (primary view)
    nrows = len(cohorts)
    heights = [6] + [4.0] * (nrows - 1) if nrows > 1 else [6]
    
    fig, axes = plt.subplots(
        nrows=nrows, ncols=1,
        figsize=(16, sum(heights)),
        gridspec_kw={'height_ratios': heights},
        sharex=True
    )
    if nrows == 1:
        axes = [axes]
    
    max_x = 0
    for ax, cohort in zip(axes, cohorts):
        df = cohort_dfs[cohort]
        s_cols = [c for c in df.columns if c not in meta]
        
        if len(s_cols) == 0:
            ax.set_axis_off()
            ax.set_title(f"{cohort}: No Coverage Data Available")
            continue
        
        # Create display labels
        df_with_labels = df.copy()
        df_with_labels["display_label"] = df_with_labels.apply(create_display_label, axis=1)
        
        fail_mask = (df[s_cols] < floor) & (~df[s_cols].isna())
        tmp = pd.DataFrame({
            "label": df_with_labels["display_label"].values,
            "fail_count": fail_mask.sum(axis=1).values
        }).dropna(subset=["label"])
        
        agg = (tmp.groupby("label", as_index=False)["fail_count"]
                  .sum().sort_values("fail_count", ascending=False))
        agg = agg[agg["fail_count"] > 0].head(30)
        
        if agg.empty:
            ax.text(0.5, 0.5, f"No amplicons below {floor}×", ha="center", va="center")
            ax.set_axis_off()
            ax.set_title(f"{cohort}: Recurrent Failures < {floor}x", fontweight='bold')
            continue
        
        colours = sns.color_palette("Reds", n_colors=len(agg))
        plot_df = agg.iloc[::-1]
        ax.barh(plot_df["label"], plot_df["fail_count"],
               color=colours[::-1], edgecolor="none")
        
        ax.set_title(f"{cohort}: Recurrent Failures < {floor}x", fontweight='bold')
        ax.set_ylabel("rsID / Annotation", fontsize=10, fontweight='bold')
        ax.grid(axis="x", linestyle=":", alpha=0.4)
        max_x = max(max_x, int(plot_df["fail_count"].max()))
    
    axes[-1].set_xlabel("Number of Failed Samples", fontsize=10, fontweight='bold')
    for ax in axes:
        ax.set_xlim(0, max_x + 1)
    
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/19_systemic_failures_by_rsid.png", dpi=300)
    plt.close()
    
    # Plot 20: By amplicon ID (alternative view for tracking)
    fig, axes = plt.subplots(
        nrows=nrows, ncols=1,
        figsize=(16, sum(heights)),
        gridspec_kw={'height_ratios': heights},
        sharex=True
    )
    if nrows == 1:
        axes = [axes]
    
    max_x = 0
    for ax, cohort in zip(axes, cohorts):
        df = cohort_dfs[cohort]
        s_cols = [c for c in df.columns if c not in meta]
        
        if len(s_cols) == 0:
            ax.set_axis_off()
            continue
        
        fail_mask = (df[s_cols] < floor) & (~df[s_cols].isna())
        
        # Use annot_id (amplicon ID) if available, otherwise use id
        id_col = "annot_id" if "annot_id" in df.columns else "id"
        
        tmp = pd.DataFrame({
            "amplicon_id": df[id_col].values,
            "fail_count": fail_mask.sum(axis=1).values
        }).dropna(subset=["amplicon_id"])
        
        agg = (tmp.groupby("amplicon_id", as_index=False)["fail_count"]
                  .sum().sort_values("fail_count", ascending=False))
        agg = agg[agg["fail_count"] > 0].head(30)
        
        if agg.empty:
            ax.text(0.5, 0.5, f"No amplicons below {floor}×", ha="center", va="center")
            ax.set_axis_off()
            ax.set_title(f"{cohort}: Recurrent Failures < {floor}x", fontweight='bold')
            continue
        
        colours = sns.color_palette("Reds", n_colors=len(agg))
        plot_df = agg.iloc[::-1]
        ax.barh(plot_df["amplicon_id"], plot_df["fail_count"],
               color=colours[::-1], edgecolor="none")
        
        ax.set_title(f"{cohort}: Recurrent Failures < {floor}x", fontweight='bold')
        ax.set_ylabel("Amplicon ID", fontsize=10, fontweight='bold')
        ax.grid(axis="x", linestyle=":", alpha=0.4)
        max_x = max(max_x, int(plot_df["fail_count"].max()))
    
    axes[-1].set_xlabel("Number of Failed Samples", fontsize=10, fontweight='bold')
    for ax in axes:
        ax.set_xlim(0, max_x + 1)
    
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/20_systemic_failures_by_amplicon_id.png", dpi=300)
    plt.close()


def plot_coverage_summary(all_cov: pd.DataFrame, plots_dir: str):
    """
    Plot 21: Global coverage distribution summary
    Shows overall distribution and comparison by cohort if available.
    """
    if all_cov.empty:
        print("[WARN] No coverage data for summary plot")
        return
    
    print("[*] Generating coverage summary plot...")
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300
    
    # Flatten coverage data for plotting
    all_depths = []
    for col in all_cov.columns:
        depths = all_cov[col].dropna()
        all_depths.extend(depths.tolist())
    
    if not all_depths:
        print("[WARN] No valid coverage depths for summary")
        return
    
    plt.figure(figsize=(12, 6))
    
    # Histogram with KDE
    sns.histplot(all_depths, bins=50, kde=True, color='steelblue')
    plt.axvline(QC_LIMITS['mean_cov'], color='red', linestyle='--', lw=2,
               label=f"Target Mean: {QC_LIMITS['mean_cov']}×")
    plt.axvline(QC_LIMITS['worst_amplicon_floor'], color='orange', linestyle='--', lw=2,
               label=f"Minimum Floor: {QC_LIMITS['worst_amplicon_floor']}×")
    
    plt.xlabel('Coverage Depth (×)', fontweight='bold')
    plt.ylabel('Frequency', fontweight='bold')
    plt.title('Coverage Depth Distribution', fontweight='bold', fontsize=14)
    plt.legend()
    plt.xlim(0, min(1000, max(all_depths)))
    plt.tight_layout()
    plt.savefig(f"{plots_dir}/21_coverage_summary.png", dpi=300)
    plt.close()


# ==============================================================================
# 6. MAIN EXECUTION
# ==============================================================================
def run_complete_pipeline():
    """
    Execute complete visualization pipeline with all plots.
    
    Generates 21 comprehensive plots:
    01-05: QC metric distributions
    06-07: PCA analysis
    08-09: Coverage heatmaps
    10-12: Panel landscapes
    13: Coverage gaps
    14-15: Retention and failures
    16-17: QC audits
    18: Worst amplicon
    19-20: Systemic failures
    21: Coverage summary
    """
    print("\n" + "="*70)
    print("COMPLETE QC VISUALIZATION PIPELINE")
    print("="*70 + "\n")
    
    # Load all data
    full_qc_df, cohort_dfs, annotation_df, all_cov, failing_samples = build_pipeline()
    
    # Export per-sample coverage tables
    print("\n[*] Exporting per-sample coverage files...")
    meta = {"chr", "start", "end", "id", "annotation", "annot_id"}
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        coords = df[["id", "chr", "start", "end"]].drop_duplicates().set_index("id")
        cohort_dir = os.path.join(PATHS["samples_dir"], cohort.replace(" ", "_"))
        os.makedirs(cohort_dir, exist_ok=True)
        
        for s in s_cols:
            out = coords.join(df.set_index("id")[s].rename("depth"))
            out_path = os.path.join(cohort_dir, f"{s}.coverage.tsv.gz")
            out.to_csv(out_path, sep="\t", compression="gzip")
    
    # Generate ALL visualizations
    print("\n" + "="*70)
    print("GENERATING ALL VISUALIZATIONS")
    print("="*70 + "\n")
    
    plot_qc_distributions(full_qc_df, PATHS["plots_dir"])              # 01-05
    plot_pca_analysis(full_qc_df, PATHS["plots_dir"])                 # 06-07
    plot_coverage_heatmaps(cohort_dfs, all_cov, failing_samples, 
                          PATHS["plots_dir"])                          # 08-09
    plot_panel_landscapes(cohort_dfs, PATHS["plots_dir"])             # 10-12
    plot_coverage_gaps(cohort_dfs, failing_samples, PATHS["plots_dir"]) # 13
    plot_retention_and_failures(full_qc_df, PATHS["plots_dir"])       # 14-15
    plot_qc_audits(full_qc_df, PATHS["plots_dir"])                    # 16-17
    plot_worst_amplicon_analysis(cohort_dfs, PATHS["plots_dir"])      # 18
    plot_systemic_amplicon_failures(cohort_dfs, PATHS["plots_dir"])   # 19-20
    plot_coverage_summary(all_cov, PATHS["plots_dir"])                # 21
    
    print("\n" + "="*70)
    print("✅ COMPLETE VISUALISATION PIPELINE FINISHED")
    print("="*70)
    print(f"\nOutputs:")
    print(f"  📊 Plots: {PATHS['plots_dir']}")
    print(f"  📁 Per-sample coverage: {PATHS['samples_dir']}")
    print(f"  📄 Full panel CSVs: {PATHS['out_dir']}")
    print(f"\n  Total plots generated: 21")
    print("="*70 + "\n")


if __name__ == "__main__":
    run_complete_pipeline()
