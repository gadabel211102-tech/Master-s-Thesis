#!/usr/bin/env python3
import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

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

STATUS_PALETTE = {
    'Pass': '#2ecc71',
    'Fail': '#e74c3c'
}

STRICT_MOSDEPTH = True
MIN_NUMERIC_FRACTION = 0.50


# ==============================================================================
# 2. UTILITIES
# ==============================================================================
def normalise_sample(s: str) -> str:
    s = str(s).strip()
    return (
        s.replace(" ", "_")
        .replace(".regions", "")
        .replace(".bed", "")
        .replace(".gz", "")
    )


def normalise_chr(c: str) -> str:
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
            best, best_frac = i, f
    return best if best_frac >= 0.5 else -1


# ==============================================================================
# 3. LOADING
# ==============================================================================
def load_annotation_for_cohort(root_path: str) -> pd.DataFrame:
    """
    YOUR BED FORMAT:
    track ...
    chr start end annotationID (e.g. 557315\nrs294150)
    So we load 4 columns only.
    """
    annot_file = os.path.join(root_path, "dna_qc", "targets.annotated.bed")
    if not os.path.exists(annot_file):
        raise FileNotFoundError(f"Annotation file not found: {annot_file}")
    ann = pd.read_csv(
        annot_file,
        sep="\t",
        header=None,
        comment='t',  # skip "track" line
        names=["chr", "start", "end", "annotation"]
    )
    ann["chr"] = ann["chr"].astype(str)
    ann["annot_id"] = ann["annotation"]
    ann = ann.drop_duplicates(["chr", "start", "end"])
    return ann


def load_coverage_for_file(path: str, sample_name: str) -> pd.DataFrame:
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
    pattern = os.path.join(root_path, "**", "*.regions.bed.gz")
    files = [
        f for f in glob.glob(pattern, recursive=True)
        if f.endswith(".regions.bed.gz") and ".bai" not in f.lower()
    ]
    if not files:
        return pd.DataFrame()

    frames = []
    for f in files:
        name = normalise_sample(os.path.basename(f).replace(".regions.bed.gz", ""))
        try:
            frames.append(load_coverage_for_file(f, name))
        except Exception as e:
            print(f"[ERROR] {cohort_name}: {f}: {e}")
            dummy = pd.DataFrame(columns=["chr", "start", "end", name])
            frames.append(dummy.set_index([]))

    combined = frames[0].copy()
    for frame in frames[1:]:
        cols = [c for c in frame.columns if c not in ["chr", "start", "end"]]
        combined = combined.join(frame[cols], how="outer")
    return combined.reset_index()


def load_qc_for_cohort(cohort_name: str, root_path: str) -> pd.DataFrame:
    qc_file = os.path.join(root_path, "dna_qc", "qc_summary.tsv")
    if not os.path.exists(qc_file):
        return pd.DataFrame()
    df = pd.read_csv(qc_file, sep="\t")
    df["cohort"] = cohort_name
    df["status"] = df["status"].str.capitalize()
    df["sample_norm"] = df["sample"].apply(normalise_sample)
    return df


# ==============================================================================
# 4. PIPELINE ASSEMBLY
# ==============================================================================
def build_pipeline():
    for p in PATHS.values():
        os.makedirs(p, exist_ok=True)

    qc_frames = []
    cohort_dfs = {}
    annotation_frames = []

    print("[*] Loading...")
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

        ann = load_annotation_for_cohort(root)
        annotation_frames.append(ann)

        # *** FIXED MERGE ***
        cov = cov.merge(
            ann[["chr", "start", "end", "annotation", "annot_id"]],
            on=["chr", "start", "end"],
            how="left"
        )
        cohort_dfs[cohort] = cov

        out_full = os.path.join(
            PATHS["out_dir"], f"GSDMB_{cohort.replace(' ', '_')}_FullPanel.csv"
        )
        cov.to_csv(out_full, index=False)

        # Audit
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

    full_qc = (pd.concat(qc_frames, ignore_index=True) if qc_frames else pd.DataFrame())
    annotation_df = (
        pd.concat(annotation_frames, ignore_index=True)
        .drop_duplicates(["chr", "start", "end"])
        if annotation_frames
        else pd.DataFrame(columns=["chr", "start", "end", "annotation", "annot_id"])
    )

    print("[*] Building global coverage...")
    all_cov_frames = []
    for cohort, df in cohort_dfs.items():
        tmp = df.set_index("id").drop(columns=["chr", "start", "end", "annotation", "annot_id"],
                                      errors="ignore")
        all_cov_frames.append(tmp)
    all_cov = pd.concat(all_cov_frames, axis=1) if all_cov_frames else pd.DataFrame()

    failing_samples = (
        full_qc[full_qc["status"] == "Fail"]["sample_norm"].tolist()
        if not full_qc.empty else []
    )

    print("[*] Done loading.")
    return full_qc, cohort_dfs, annotation_df, all_cov, failing_samples


# ==============================================================================
# 5. PLOTTING
# ==============================================================================
def plot_all(full_qc_df, cohort_dfs, annotation_df, all_cov, failing_samples):
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    # ======================================================================
    # 01. PCA
    # ======================================================================
    if not full_qc_df.empty:
        feats = ['mean_cov', 'on_target_pct', 'total_reads', 'mapped_pct']
        X = StandardScaler().fit_transform(full_qc_df[feats].fillna(0))

        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X)
        var1 = pca.explained_variance_ratio_[0] * 100
        var2 = pca.explained_variance_ratio_[1] * 100

        full_qc_df["PC1"], full_qc_df["PC2"] = pcs[:, 0], pcs[:, 1]

        plt.figure(figsize=(8, 6))
        ax = sns.scatterplot(
            data=full_qc_df, x="PC1", y="PC2",
            hue="cohort", style="status",
            palette=PALETTE, s=100
        )
        ax.set_xlabel(f"PC1 ({var1:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({var2:.1f}% variance)")
        plt.title("PCA of QC Metrics")
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/01_pca.png")
        plt.close()

        # 01b. Pass-only PCA
        pass_only = full_qc_df[full_qc_df["status"] == "Pass"].copy()
        if len(pass_only) >= 3:
            X = StandardScaler().fit_transform(pass_only[feats].fillna(0))

            pca = PCA(n_components=2)
            pcs = pca.fit_transform(X)
            var1 = pca.explained_variance_ratio_[0] * 100
            var2 = pca.explained_variance_ratio_[1] * 100

            pass_only["PC1"], pass_only["PC2"] = pcs[:, 0], pcs[:, 1]

            plt.figure(figsize=(8, 6))
            ax = sns.scatterplot(
                data=pass_only, x="PC1", y="PC2",
                hue="cohort", palette=PALETTE, s=100
            )
            ax.set_xlabel(f"PC1 ({var1:.1f}% variance)")
            ax.set_ylabel(f"PC2 ({var2:.1f}% variance)")
            plt.title("PCA (Pass-Only Samples)")
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(f"{PATHS['plots_dir']}/01b_pca_pass_only.png")
            plt.close()

   # 1. First, ensure your failure list is "cleaned" once at the start
# This ensures it matches the format the labels will use
clean_failing_samples = {normalise_sample(s) for s in failing_samples}

# ======================================================================
# 02. Cohort heatmaps
# ======================================================================
for cohort, df in cohort_dfs.items():
    s_cols = [c for c in df.columns
              if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]
    dmat = df.set_index("id")[s_cols]
    logmat = np.log10(dmat + 1)
    mask = dmat.isna()

    plt.figure(figsize=(12, 8))
    sns.heatmap(logmat, cmap="magma", mask=mask, vmin=0, vmax=4)
    ax = plt.gca()
    plt.xticks(rotation=45, ha='right', fontsize=8)
    
    # HIGHLIGHT LOGIC
    for label in ax.get_xticklabels():
        # Normalise the label text before checking the list
        sample_on_plot = normalise_sample(label.get_text())
        if sample_on_plot in clean_failing_samples:
            label.set_color("red")
            label.set_weight("bold")
        else:
            # UNCOMMENT THE LINE BELOW TO DEBUG IN THE TERMINAL:
            # print(f"DEBUG: '{sample_on_plot}' not found in {clean_failing_samples}")
            pass

    plt.xlabel("Sample ID", fontsize=10, fontweight='bold')
    plt.ylabel("Genomic Coordinates", fontsize=10, fontweight='bold') # Fixed the extra bracket here
    plt.title(f"Coverage Heatmap: {cohort}")
    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/02_heatmap_{cohort.replace(' ', '_')}.png")
    plt.close()

	# 1. First, ensure your failure list is "cleaned" once at the start
# This ensures it matches the format the labels will use
clean_failing_samples = {normalise_sample(s) for s in failing_samples}

# ======================================================================
# 02. Cohort heatmaps
# ======================================================================
for cohort, df in cohort_dfs.items():
    s_cols = [c for c in df.columns
              if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]
    dmat = df.set_index("id")[s_cols]
    logmat = np.log10(dmat + 1)
    mask = dmat.isna()

    plt.figure(figsize=(12, 8))
    sns.heatmap(logmat, cmap="magma", mask=mask, vmin=0, vmax=4)
    ax = plt.gca()
    plt.xticks(rotation=45, ha='right', fontsize=8)
    
    # HIGHLIGHT LOGIC
    for label in ax.get_xticklabels():
        # Normalise the label text before checking the list
        sample_on_plot = normalise_sample(label.get_text())
        if sample_on_plot in clean_failing_samples:
            label.set_color("red")
            label.set_weight("bold")
        else:
            # UNCOMMENT THE LINE BELOW TO DEBUG IN THE TERMINAL:
            # print(f"DEBUG: '{sample_on_plot}' not found in {clean_failing_samples}")
            pass

    plt.xlabel("Sample ID", fontsize=10, fontweight='bold')
    plt.ylabel("Genomic Coordinates", fontsize=10, fontweight='bold') # Fixed the extra bracket here
    plt.title(f"Coverage Heatmap: {cohort}")
    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/02_heatmap_{cohort.replace(' ', '_')}.png")
    plt.close()

# ======================================================================
# 03. Global Heatmap
# ======================================================================
if not all_cov.empty:
    plt.figure(figsize=(16, 10))
    sns.heatmap(np.log10(all_cov + 1), cmap="magma",
                mask=all_cov.isna(), vmin=0, vmax=4)
    ax = plt.gca()
    plt.xticks(rotation=45, ha='right', fontsize=6)
    
    for label in ax.get_xticklabels():
        sample_on_plot = normalise_sample(label.get_text())
        if sample_on_plot in clean_failing_samples:
            label.set_color("red")
            label.set_weight("bold")

    plt.xlabel("Sample ID", fontsize=10, fontweight='bold')
    plt.ylabel("Genomic Coordinates", fontsize=10, fontweight='bold') # Fixed the extra bracket here
    plt.title("Global Coverage Heatmap")
    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/03_heatmap_global.png")
    plt.close()
    # ======================================================================
    # 05. Panel Landscape (index order)
    # ======================================================================
    plt.figure(figsize=(15, 6))
    ax = plt.gca()
    cohort_handles = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns
                  if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]
        med = df[s_cols].median(axis=1)
        line, = ax.plot(
            range(len(med)),
            med,
            lw=2,
            label=cohort,
            color=PALETTE[cohort]
        )
        cohort_handles.append(line)

    mean_line = QC_LIMITS['mean_cov']
    floor_line = QC_LIMITS['worst_amplicon_floor']
    h_mean = ax.axhline(mean_line, color="gold", ls=":", lw=2, label=f"Mean ≥ {mean_line}×")
    h_floor = ax.axhline(floor_line, color="red", ls="--", lw=2, label=f"Floor {floor_line}×")

    ax.set_yscale("log")
    ax.set_xlabel("Amplicon Index")
    ax.set_ylabel("Median Depth (×)")
    ax.set_title("Panel Landscape vs Thresholds")

    # Fixed legend
    handles = cohort_handles + [h_mean, h_floor]
    labels = [h.get_label() for h in handles]
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/05_landscape_index.png")
    plt.close()

    # ======================================================================
    # 06. Panel landscape ordered by genomic coordinate
    # ======================================================================
    def _chr_order(c):
        c = str(c).replace("chr", "")
        mapping = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        try:
            return int(c)
        except Exception:
            return mapping.get(c, 1000)

    coord_frames = []
    for cohort, df in cohort_dfs.items():
        coord_frames.append(df[["chr", "start", "end", "id"]].drop_duplicates())
    kdf = pd.concat(coord_frames).drop_duplicates()
    kdf["key"] = kdf["chr"].apply(_chr_order)
    kdf = kdf.sort_values(["key", "start", "end"])
    ordered_ids = kdf["id"].tolist()

    plt.figure(figsize=(18, 6))
    ax = plt.gca()
    cohort_handles = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns
                  if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]
        med = df[s_cols].median(axis=1)
        t = (
            pd.DataFrame({"id": df["id"], "med": med})
            .drop_duplicates("id")
            .set_index("id")
            .reindex(ordered_ids)
        )
        xvals = np.arange(len(ordered_ids))
        line, = ax.plot(
            xvals,
            t["med"].values,
            lw=2,
            label=cohort,
            color=PALETTE[cohort]
        )
        cohort_handles.append(line)

    h_mean = ax.axhline(mean_line, color="gold", ls=":", lw=2, label=f"Mean ≥ {mean_line}×")
    h_floor = ax.axhline(floor_line, color="red", ls="--", lw=2, label=f"Floor {floor_line}×")

    if ordered_ids:
        n = len(ordered_ids)
        step = max(1, n // 40)
        sel = list(range(0, n, step))
        ax.set_xticks(sel)
        ax.set_xticklabels(
            [ordered_ids[i] for i in sel],
            rotation=45,
            ha='right',
            fontsize=7
        )

    ax.set_yscale("log")
    ax.set_xlabel("Amplicon ID (Genomic Coordinates)")
    ax.set_ylabel("Median Depth (×)")
    ax.set_title("Panel Landscape (Genomic Coordinates)")

    handles = cohort_handles + [h_mean, h_floor]
    labels = [h.get_label() for h in handles]
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/06_landscape_genomic_ids.png")
    plt.close()

    # ======================================================================
    # 07. Panel landscape by gene
    # ======================================================================
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
        grp = (
            tmp.groupby("gene", as_index=False)
               .agg({"chr": "first", "start": "min", "end": "max"})
        )
        grp["key"] = grp["chr"].apply(_chr_order)
        gene_frames.append(grp)

    if not gene_frames:
        print("[INFO] 07 (genes) skipped: no gene annotations found.")
    else:
        gdf = pd.concat(gene_frames, ignore_index=True).drop_duplicates("gene")
        gdf = gdf.sort_values(["key", "start", "end"])
        gene_order = gdf["gene"].tolist()

        plt.figure(figsize=(18, 6))
        ax = plt.gca()
        cohort_handles = []
        for cohort, df in cohort_dfs.items():
            s_cols = [c for c in df.columns
                      if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]
            med_amp = df[s_cols].median(axis=1)
            tmp = pd.DataFrame({"annotation": df["annotation"], "med": med_amp}).copy()
            tmp["gene"] = tmp["annotation"].apply(parse_gene_or_annotation)
            tmp = tmp.dropna(subset=["gene"])
            gmed = tmp.groupby("gene")["med"].median().reindex(gene_order)
            xvals = np.arange(len(gene_order))
            line, = ax.plot(
                xvals,
                gmed.values,
                lw=2,
                label=cohort,
                color=PALETTE[cohort]
            )
            cohort_handles.append(line)

        h_mean = ax.axhline(mean_line, color="gold", ls=":", lw=2)
        h_floor = ax.axhline(floor_line, color="red", ls="--", lw=2)

        n = len(gene_order)
        step = max(1, n // 40)
        sel = list(range(0, n, step))
        ax.set_xticks(sel)
        ax.set_xticklabels(
            [gene_order[i] for i in sel],
            rotation=45,
            ha='right',
            fontsize=7
        )

        ax.set_yscale("log")
        ax.set_xlabel("Annotated Region (Genomic Order)")
        ax.set_ylabel("Median Depth (×)")
        ax.set_title("Panel Landscape (Region)")

        handles = cohort_handles + [h_mean, h_floor]
        labels = [h.get_label() for h in handles]
        ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left')

        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/07_landscape_genes.png")
        plt.close()

	
    # ======================================================================
    # 08. Coverage Gaps (<50×)
    # ======================================================================
    gap_stats = []
    for cohort, df in cohort_dfs.items():
        s_cols = [
            c for c in df.columns
            if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]
        ]
        for s in s_cols:
            vals = df[s]
            gaps = (vals < QC_LIMITS['worst_amplicon_floor']) & (~vals.isna())
            gap_stats.append({
                "cohort": cohort,
                "sample": s,
                "gaps": int(gaps.sum())
            })

    gap_df = pd.DataFrame(gap_stats)

    plt.figure(figsize=(10, 7))
    ax = sns.boxplot(
        data=gap_df,
        x="cohort",
        y="gaps",
        hue="cohort",
        palette=PALETTE,
        legend=False
    )

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
        cohort = cohort_levels[nearest_x]
        matches = gap_df[(gap_df["cohort"] == cohort) & (gap_df["gaps"] == int(y))]
        is_fail = any(
            normalise_sample(row["sample"]) in failing_samples
            for _, row in matches.iterrows()
        )
        if is_fail:
            line.set_color("red")
            line.set_markeredgecolor("black")
            line.set_markersize(7)

    plt.title("Coverage Gaps (<50×)")
    plt.xlabel("Cohort")
    plt.ylabel("Amplicons <50×")
    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/08_coverage_gaps.png")
    plt.close()


    # ======================================================================
    # 09. Retention Rate
    # ======================================================================
    if not full_qc_df.empty:
        plt.figure(figsize=(10, 7))
        ax = full_qc_df.groupby(["cohort", "status"]).size().unstack().fillna(0).plot(
            kind="bar",
            stacked=True,
            color=STATUS_PALETTE,
            ax=plt.gca()
        )
        plt.title("Retention Rate")
        plt.xlabel("Cohort")
        plt.ylabel("Sample Count")
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/09_retention.png")
        plt.close()

	# ======================================================================
    # 09b. Failure reasons by cohort
    # ======================================================================
    if not full_qc_df.empty:
        df_rr = full_qc_df.copy()
        df_rr["fail_reasons"] = df_rr["fail_reasons"].fillna("")
        df_rr["fail_list"] = df_rr["fail_reasons"].apply(
            lambda s: [x.strip() for x in str(s).split(",") if x.strip()]
        )
        failed = df_rr[df_rr["status"] == "Fail"].explode("fail_list")
        failed["fail_list"] = failed["fail_list"].replace("", np.nan).fillna("unspecified")

        if not failed.empty:
            # --- DEFINE PROPER NAMES HERE ---
            reason_map = {
                "low_reads": "Low Read Count",
                "low_uniformity": "Low Uniformity (<95% @ 100x)",
                "low_on_target": "Low On-Target Rate",
                "unspecified": "Unspecified Error"
            }
            # Apply the mapping to the column
            failed["fail_list"] = failed["fail_list"].map(lambda x: reason_map.get(x, x))
            # --------------------------------

            counts = (
                failed.groupby(["cohort", "fail_list"])
                      .size()
                      .reset_index(name="n")
                      .pivot(index="cohort", columns="fail_list", values="n")
                      .fillna(0)
            )
            counts = counts[counts.sum(axis=0).sort_values(ascending=False).index]

            plt.figure(figsize=(12, 7))
            counts.plot(kind="bar", stacked=True, ax=plt.gca(), cmap="tab20")
            
            # Additional Legend Formatting
            plt.legend(title="Reason for Failure", bbox_to_anchor=(1.05, 1), loc='upper left')
            
            plt.title("Failure Reasons by Cohort")
            plt.xlabel("Cohort", fontsize=10, fontweight='bold')
            plt.ylabel("Failed Samples", fontsize=10, fontweight='bold')
            plt.tight_layout()
            plt.savefig(f"{PATHS['plots_dir']}/09b_fail_reasons.png")
            plt.close()

    # ======================================================================
    # 10. Mapping efficiency audit
    # ======================================================================
    if not full_qc_df.empty:
        plt.figure(figsize=(8, 6))
        sns.stripplot(
            data=full_qc_df,
            x="cohort",
            y="mapped_pct",
            hue="status",
            palette=STATUS_PALETTE,
            jitter=True
        )
        plt.axhline(QC_LIMITS["mapped_pct"], color="red", ls="--")
        plt.title("Mapping Efficiency Audit")
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/10_mapping_audit.png")
        plt.close()

    # ======================================================================
    # 11. Specificity jitter
    # ======================================================================
    if not full_qc_df.empty:
        plt.figure(figsize=(8, 6))
        sns.stripplot(
            data=full_qc_df,
            x="cohort",
            y="on_target_pct",
            hue="status",
            palette=STATUS_PALETTE,
            jitter=True
        )
        plt.axhline(QC_LIMITS["on_target_pct"], color="firebrick", ls="--")
        plt.title("Specificity Jitter")
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/11_specificity_jitter.png")
        plt.close()

    # ======================================================================
    # 12. Lowest coverage per sample
    # ======================================================================
    worst = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns
                  if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]
        for s in s_cols:
            md = df[s].min()
            if pd.notna(md):
                worst.append({"cohort": cohort, "min_depth": md})

    worst_df = pd.DataFrame(worst)
    plt.figure(figsize=(10, 6))
    sns.swarmplot(
        data=worst_df,
        x="cohort",
        y="min_depth",
        hue="cohort",
        palette=PALETTE,
        size=6
    )
    plt.axhline(QC_LIMITS["worst_amplicon_floor"], color="red", ls="--")
    plt.yscale("symlog", linthresh=10)
    plt.title("Lowest Coverage per Sample")
    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/12_worst_amplicon.png")
    plt.close()

	
    # ======================================================================
    # 13. Systemic Amplicon Failures — by cohort
    # ======================================================================
    if cohort_dfs:
        cohorts = list(cohort_dfs.keys())
        nrows = len(cohorts)
        height_per = 4.0
        
        # Patch 3: Enlarged first cohort
        heights = [6] + [height_per] * (nrows - 1)
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
            s_cols = [c for c in df.columns
                      if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]

            if len(s_cols) == 0:
                ax.set_axis_off()
                ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)")
                continue

            floor = QC_LIMITS['worst_amplicon_floor']
            fail_mask = (df[s_cols] < floor) & (~df[s_cols].isna())

            tmp = pd.DataFrame({
                "annotation": df["annotation"].values,
                "fail_count": fail_mask.sum(axis=1).values
            }).dropna(subset=["annotation"])

            agg = (
                tmp.groupby("annotation", as_index=False)["fail_count"]
                   .sum()
                   .sort_values("fail_count", ascending=False)
            )

            agg = agg[agg["fail_count"] > 0]
            top_n = 30
            agg = agg.head(top_n)

            if agg.empty:
                ax.text(0.5, 0.5, "No amplicons below 50×", ha="center", va="center")
                ax.set_axis_off()
                ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)")
                continue

            colours = sns.color_palette("Reds", n_colors=len(agg))
            plot_df = agg.iloc[::-1]
            ax.barh(
                plot_df["annotation"],
                plot_df["fail_count"],
                color=colours[::-1],
                edgecolor="none"
            )

            ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)")
            ax.set_ylabel("Amplicon")
            ax.grid(axis="x", linestyle=":", alpha=0.4)
            max_x = max(max_x, int(plot_df["fail_count"].max()))

        axes[-1].set_xlabel("Failed Samples")
        for ax in axes:
            ax.set_xlim(0, max_x + 1)

        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/13_failures_global.png")
        plt.close()

    return
# ======================================================================
    # 13b. Systemic Amplicon Failures — by cohort (Genomic Coordinates)
    # ======================================================================
    if cohort_dfs:
        cohorts = list(cohort_dfs.keys())
        nrows = len(cohorts)
        height_per = 4.0
        
        # Patch 3: Enlarged first cohort
        heights = [6] + [height_per] * (nrows - 1)
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
            s_cols = [c for c in df.columns
                      if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]

            if len(s_cols) == 0:
                ax.set_axis_off()
                ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)")
                continue

            floor = QC_LIMITS['worst_amplicon_floor']
            fail_mask = (df[s_cols] < floor) & (~df[s_cols].isna())

            # --- CHANGE: Use 'id' (coordinates) instead of 'annotation' ---
            tmp = pd.DataFrame({
                "coord": df["id"].values, # Changed from annotation to id
                "fail_count": fail_mask.sum(axis=1).values
            }).dropna(subset=["coord"])

            agg = (
                tmp.groupby("coord", as_index=False)["fail_count"]
                   .sum()
                   .sort_values("fail_count", ascending=False)
            )
            # -------------------------------------------------------------

            agg = agg[agg["fail_count"] > 0]
            top_n = 30
            agg = agg.head(top_n)

            if agg.empty:
                ax.text(0.5, 0.5, "No amplicons below 50×", ha="center", va="center")
                ax.set_axis_off()
                ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)")
                continue

            colours = sns.color_palette("Reds", n_colors=len(agg))
            plot_df = agg.iloc[::-1]
            ax.barh(
                plot_df["coord"], # Using coordinates here
                plot_df["fail_count"],
                color=colours[::-1],
                edgecolor="none"
            )

            ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)")
            ax.set_ylabel("Genomic Coordinate (Chr17)", fontsize=10, fontweight='bold')
            ax.grid(axis="x", linestyle=":", alpha=0.4)
            max_x = max(max_x, int(plot_df["fail_count"].max()))

        axes[-1].set_xlabel("Number of Failed Samples", fontsize=10, fontweight='bold')
        for ax in axes:
            ax.set_xlim(0, max_x + 1)

        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/13_failures_coordinates.png") # Saved with a new name
        plt.close()

# ==============================================================================
# 6. MAIN EXECUTION
# ==============================================================================
def run_pipeline():
    full_qc_df, cohort_dfs, annotation_df, all_cov, failing_samples = build_pipeline()

    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns
                  if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]
        coords = df[["id", "chr", "start", "end"]].drop_duplicates().set_index("id")
        cohort_dir = os.path.join(PATHS["samples_dir"], cohort.replace(" ", "_"))
        os.makedirs(cohort_dir, exist_ok=True)
        for s in s_cols:
            out = coords.join(df.set_index("id")[s].rename("depth"))
            out_path = os.path.join(cohort_dir, f"{s}.coverage.tsv.gz")
            out.to_csv(out_path, sep="\t", compression="gzip")

    plot_all(full_qc_df, cohort_dfs, annotation_df, all_cov, failing_samples)
    print("\n[SUCCESS] Harmonised CSVs, diagnostics, and visualisations completed.")


if __name__ == "__main__":
    run_pipeline()
