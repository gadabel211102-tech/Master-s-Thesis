#!/usr/bin/env python3  # Tells the computer to run this script using Python 3
import os               # Import library for interacting with the operating system (folders/paths)
import glob             # Import library for finding files using patterns (like *.bam)
import pandas as pd     # Data analysis library (used for tables/spreadsheets)
import numpy as np      # Math library for working with numbers and arrays
import seaborn as sns # type: ignore  # Library for making attractive statistical plots
import matplotlib.pyplot as plt      # The core library for creating visualisations/graphs
from sklearn.preprocessing import StandardScaler # type: ignore # Tool to scale data for fair comparison
from sklearn.decomposition import PCA # type: ignore        # Tool for Principal Component Analysis (simplifies data)
import matplotlib.ticker as mticker  # Allows for fine control over graph axis labels

# ==============================================================================
# 1. CONFIGURATION (Defining the "Gold Standards")
# ==============================================================================
QC_LIMITS = {           # A dictionary (list of rules) for Quality Control
    'total_reads': 120000,      # Minimum number of total reads required
    'mapped_pct': 90.0,         # Minimum percentage of reads that must map to the genome
    'on_target_pct': 75.0,      # Minimum percentage of reads on the specific DNA targets
    'mean_cov': 200.0,          # Minimum average "depth" (how many times a base is read)
    'uniformity_100x': 95.0,    # Target for how evenly the coverage is spread
    'worst_amplicon_floor': 50.0 # The absolute minimum coverage allowed for the worst region
}

PATHS = {               # A list of folder paths the script will use/create
    "out_dir": "./gsdmb_final_results",
    "samples_dir": "./gsdmb_final_results/individual_sample_coverage",
    "plots_dir": "./gsdmb_final_results/plots",
    "raw_root": "./"
}

COHORT_MAP = {          # Maps "Human-friendly" names to the actual folder structure
    'Breast Tumour': 'breast/tumour',
    'Breast Control': 'breast/normal',
    'Endometrium Tumour': 'endometrium/tumour',
    'Endometrium Control': 'endometrium/normal',
}

PALETTE = {             # Defines specific colours for the graphs (HEX codes)
    'Breast Tumour': '#3498db',
    'Breast Control': "#d1d946",
    'Endometrium Tumour': '#9b59b6',
    'Endometrium Control': '#2ecc71'
}

STATUS_PALETTE = {      # Colours for "Pass" (Green) and "Fail" (Red) results
    'Pass': '#2ecc71',
    'Fail': '#e74c3c'
}

# Configuration flags for handling 'mosdepth' (a tool that measures DNA sequencing depth)
STRICT_MOSDEPTH = True           # If True, the script is less tolerant of messy data
MIN_NUMERIC_FRACTION = 0.90      # At least 90% of a column must be numbers to be valid


# ==============================================================================
# 2. UTILITIES (Small helper functions to clean up data)
# ==============================================================================
def normalise_sample(s: str) -> str:  # Function to make sample names consistent
    """Standardise sample names."""
    s = str(s).strip().lower()        # Remove spaces and make it lowercase
    return (
        s.replace(" ", "_")           # Replace any remaining spaces with underscores
         .replace(".regions", "")     # Remove file extensions from the name
         .replace(".bed", "")         
         .replace(".gz", "")          
    )                                 # Returns the "clean" name


def normalise_chr(c: str) -> str:     # Function to ensure chromosome names match (e.g., '17' vs 'chr17')
    """Ensure chromosome naming is consistent (chr17 etc.)."""
    c = str(c)                        # Convert input to a string
    if c.startswith("chr"):           # If it already starts with 'chr', leave it alone
        return c
    if c in {"X", "Y", "M", "MT"}:    # If it's a special chromosome...
        return "chrM" if c in {"M", "MT"} else f"chr{c}" # ...format it as chrM, chrX, or chrY
    try:
        return f"chr{int(c)}"         # If it's a number (like 17), turn it into 'chr17'
    except Exception:                 # If something weird happens...
        return c if c.startswith("chr") else f"chr{c}" # ...try one last time to add 'chr'


def pick_depth_column(df: pd.DataFrame) -> int: # Function to find which column in a table has the 'depth' data
    """Choose appropriate depth column if BED+mosdepth output has >4 cols."""
    ncols = df.shape[1]               # Get the total number of columns in the table
    if ncols == 4:                    # If there are exactly 4, the 4th (index 3) is usually the depth
        return 3
    if ncols < 4:                     # If fewer than 4, it's not a valid file
        return -1
    last = df.columns[-1]             # Look at the very last column
    # Check if the last column is actually made of numbers (not text)
    frac = pd.to_numeric(df[last], errors='coerce').notna().mean()
    if frac > 0.8:                    # If more than 80% is numeric, this is our depth column
        return ncols - 1
    best, best_frac = -1, -1          # Otherwise, search through columns to find the most "numeric" one
    for i in range(3, ncols):         # Check from the 4th column onwards
        f = pd.to_numeric(df.iloc[:, i], errors='coerce').notna().mean()
        if f > best_frac:             # Keep track of which column has the highest percentage of numbers
            best = i
            best_frac = f
    return best if best_frac >= 0.5 else -1 # Return the best column found, or -1 if none look like depth
# ==============================================================================
# 3. LOADING: annotation, per-file coverage, per-cohort coverage, QC
# ==============================================================================
def load_annotation_for_cohort(root_path: str) -> pd.DataFrame:
    """Load per-cohort annotation BED."""
    # Build the path to the BED file which contains the DNA target regions
    annot_file = os.path.join(root_path, "dna_qc", "targets.clean.bed")
    if not os.path.exists(annot_file): # If file is missing, stop and report an error
        raise FileNotFoundError(f"Annotation file not found: {annot_file}")

    # Read the file into a table (Pandas DataFrame)
    ann = pd.read_csv(
        annot_file,
        sep="\t",           # Data is separated by Tabs
        header=None,        # File has no title row
        names=["chr", "start", "end", "annotation"] # Assign these column names
    )
    ann["chr"] = ann["chr"].astype(str) # Ensure chromosome names are treated as text
    ann["annot_id"] = ann["annotation"] # Create a copy of the annotation column for identification
    return ann.drop_duplicates(["chr", "start", "end"]) # Remove any accidental duplicate rows


def load_coverage_for_file(path: str, sample_name: str) -> pd.DataFrame:
    """Load one mosdepth .regions.bed.gz file and normalise."""
    # Read a compressed (.gz) file directly into a table
    raw = pd.read_csv(path, sep="\t", header=None, compression="gzip", dtype=str)
    if raw.empty: # Handle empty files
        if STRICT_MOSDEPTH: # If we are in "strict mode," crash on empty files
            raise ValueError(f"{path} is empty")
        return pd.DataFrame(columns=["chr", "start", "end", sample_name]).set_index([])

    raw.iloc[:, 0] = raw.iloc[:, 0].map(normalise_chr) # Clean chromosome names (e.g., add 'chr')
    raw.iloc[:, 1] = pd.to_numeric(raw.iloc[:, 1], errors='coerce') # Convert start position to number
    raw.iloc[:, 2] = pd.to_numeric(raw.iloc[:, 2], errors='coerce') # Convert end position to number

    depth_col = pick_depth_column(raw) # Find which column holds the coverage depth
    depth = (
        pd.to_numeric(raw.iloc[:, depth_col], errors='coerce') # Convert depth to numbers
        if depth_col != -1 else pd.Series(np.nan, index=raw.index) # If missing, use NaN (empty)
    )

    # Build a clean table for this specific sample
    df = pd.DataFrame({
        "chr": raw.iloc[:, 0].astype(str),
        "start": raw.iloc[:, 1].astype("Int64"),
        "end": raw.iloc[:, 2].astype("Int64"),
        sample_name: depth
    })
    # Create a unique ID string like "chr17:100-200" to identify each region
    df["id"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str) # type: ignore
    df = df.drop_duplicates("id").set_index("id") # Use the ID as the table index (key)

    frac = df[sample_name].notna().mean() # Calculate what percentage of data is actually numeric
    if frac < MIN_NUMERIC_FRACTION and STRICT_MOSDEPTH: # If data is too "empty," report error
        raise ValueError(f"{sample_name}: only {frac:.1%} numeric")
    return df


def load_coverage_for_cohort(cohort_name: str, root_path: str) -> pd.DataFrame:
    """Load all mosdepth outputs for a cohort and join by ID."""
    # Search for all files ending in .regions.bed.gz in all subfolders (**)
    pattern = os.path.join(root_path, "**", "*.regions.bed.gz")
    files = [
        f for f in glob.glob(pattern, recursive=True)
        if f.endswith(".regions.bed.gz") and ".bai" not in f.lower() # Ignore index files
    ]
    if not files: # If no files found, return an empty table
        return pd.DataFrame()

    frames = [] # Initialise list to store tables for each sample
    for f in sorted(files): # Loop through found files alphabetically
        name = normalise_sample(os.path.basename(f).replace(".regions.bed.gz", "")) # Clean sample name
        try:
            frames.append(load_coverage_for_file(f, name)) # Load and add to our list
        except Exception as e: # Catch errors to explain exactly which file failed
            raise RuntimeError(f"[MOSDEPTH ERROR] {cohort_name} {name}: {e}")

    combined = frames[0].copy() # Start with the first sample table
    for frame in frames[1:]: # Join all other sample tables to the first one
        cols = [c for c in frame.columns if c not in ["chr", "start", "end"]] # Only get depth column
        combined = combined.join(frame[cols], how="outer") # Combine based on matching region IDs

    return combined.reset_index() # Return the final multi-sample table


def load_qc_for_cohort(cohort_name: str, root_path: str) -> pd.DataFrame:
    """Load QC summary for one cohort."""
    qc_file = os.path.join(root_path, "dna_qc", "qc_summary.tsv") # Locate the QC file
    if not os.path.exists(qc_file):
        return pd.DataFrame()
    df = pd.read_csv(qc_file, sep="\t") # Read the tab-separated QC table
    df["cohort"] = cohort_name # Add a column to identify which group this belongs to
    df["status"] = df["status"].str.capitalize() # Standardise status (e.g., 'Pass', 'Fail')
    df["sample_norm"] = df["sample"].apply(normalise_sample) # Standardise sample names
    return df


# ==============================================================================
# 4. Pipeline Assembly
# ==============================================================================
def build_pipeline():
    """Load everything: QC, coverage, annotation; produce merged tables."""
    for p in PATHS.values(): # Create all necessary output folders
        os.makedirs(p, exist_ok=True)

    qc_frames = [] # List for all QC data
    cohort_dfs = {} # Dictionary to store coverage tables per cohort
    annotation_frames = [] # List for all target annotations

    print("[*] Loading...")
    for cohort, relpath in COHORT_MAP.items(): # Loop through the groups (Breast, Endometrium, etc.)
        root = os.path.join(PATHS["raw_root"], relpath)
        if not os.path.exists(root):
            print(f"[WARN] Missing path {root}")
            continue

        qc = load_qc_for_cohort(cohort, root) # Load QC data
        if not qc.empty:
            qc_frames.append(qc)

        cov = load_coverage_for_cohort(cohort, root) # Load Coverage data
        if cov.empty:
            print(f"[WARN] No coverage for {cohort}")
            continue

        ann = load_annotation_for_cohort(root) # Load BED annotations
        annotation_frames.append(ann)

        # Merge the coverage data with the annotation info (like chromosome and gene name)
        cov = cov.merge(
            ann[["chr", "start", "end", "annotation", "annot_id"]],
            on=["chr", "start", "end"],
            how="left"
        )
        cohort_dfs[cohort] = cov # Store the merged table

        # Save this cohort's full data to a CSV file
        out_full = os.path.join(
            PATHS["out_dir"], f"GSDMB_{cohort.replace(' ', '_')}_FullPanel.csv"
        )
        cov.to_csv(out_full, index=False)

        # Audit Section: Calculate statistics for each sample
        s_cols = [c for c in cov.columns
                  if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]
        rows = []
        for s in s_cols: # For every sample...
            depth = cov[s]
            tot = len(depth)
            num = depth.notna().sum() # How many regions have actual numbers
            zeros = (depth == 0).sum() # How many regions have 0 coverage
            rows.append({
                "sample": s,
                "amplicons_total": tot,
                "numeric_depth": num,
                "zeros": zeros,
                "numeric_fraction": num / tot if tot else 0
            })
        audit = pd.DataFrame(rows) # Save these stats to a separate "Audit" file
        audit.to_csv(
            os.path.join(PATHS["out_dir"], f"CoverageAudit_{cohort.replace(' ', '_')}.csv"),
            index=False
        )

    # Combine all QC data from all cohorts into one master table
    full_qc = pd.concat(qc_frames, ignore_index=True) if qc_frames else pd.DataFrame()
    # Combine all annotation data into one master list
    annotation_df = (
        pd.concat(annotation_frames, ignore_index=True)
          .drop_duplicates(["chr", "start", "end"])
        if annotation_frames else
        pd.DataFrame(columns=["chr", "start", "end", "annotation", "annot_id"])
    )

    print("[*] Building global coverage...")
    all_cov_frames = []
    for cohort, df in cohort_dfs.items(): # Prepare data to merge all cohorts together
        tmp = df.set_index("id").drop(
            columns=["chr", "start", "end", "annotation", "annot_id"],
            errors="ignore"
        )
        all_cov_frames.append(tmp)
    # This creates one giant table with EVERY sample from EVERY cohort
    all_cov = pd.concat(all_cov_frames, axis=1) if all_cov_frames else pd.DataFrame()

    # Create a list of samples that failed the QC check
    failing_samples = (
        full_qc[full_qc["status"] == "Fail"]["sample_norm"].tolist()
        if not full_qc.empty else []
    )

    print("[*] Done loading.")
    return full_qc, cohort_dfs, annotation_df, all_cov, failing_samples # Return everything to the main script

# ==============================================================================
# 5. PLOTTING: Visualising results using graphs and heatmaps
# ==============================================================================
def plot_all(full_qc_df, cohort_dfs, annotation_df, all_cov, failing_samples):
    # Set the general look for the graphs (font size and high resolution)
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams["figure.dpi"] = 300

    # ==========================================================================
    # 01. PCA (Principal Component Analysis)
    # ==========================================================================
    if not full_qc_df.empty:
        # Define which QC metrics we want to compare
        feats = ['mean_cov', 'on_target_pct', 'total_reads', 'mapped_pct']
        
        # Standardise the data (make them all comparable, even with different units)
        X = StandardScaler().fit_transform(full_qc_df[feats].fillna(0))
        
        # Collapse the 4 metrics into 2 "Principal Components" (PC1 and PC2)
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X)
        
        # Calculate how much "truth" about the data is captured in each PC
        var1 = pca.explained_variance_ratio_[0] * 100
        var2 = pca.explained_variance_ratio_[1] * 100

        # Save the new coordinates back into our main table
        full_qc_df["PC1"], full_qc_df["PC2"] = pcs[:, 0], pcs[:, 1]

        # Create the plot
        plt.figure(figsize=(8, 6))
        ax = sns.scatterplot(
            data=full_qc_df, x="PC1", y="PC2",
            hue="cohort", style="status", # Colour by group, shape by Pass/Fail
            palette=PALETTE, s=100
        )
        ax.set_xlabel(f"PC1 ({var1:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({var2:.1f}% variance)")
        plt.title("PCA of QC Metrics", fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left') # Put legend outside
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/01_pca.png")
        plt.close() # Close plot to free up computer memory

        # 01b. Pass-only PCA: Repeat the steps above but ONLY for samples that passed QC
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
            plt.title("PCA (Pass-Only Samples)", fontweight='bold')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(f"{PATHS['plots_dir']}/01b_pca_pass_only.png")
            plt.close()

    # Create a set of "clean" names for samples that failed (to colour them red later)
    clean_failing_samples = {normalise_sample(s) for s in failing_samples}

    # ==========================================================================
    # 02. Cohort heatmaps
    # ==========================================================================
    meta = {"chr", "start", "end", "id", "annotation", "annot_id"} # Non-data columns

    for cohort, df in cohort_dfs.items():
        # Get only the columns that represent actual samples
        s_cols = [c for c in df.columns if c not in meta]
        dmat = df.set_index("id")[s_cols] # dmat stands for 'data matrix'

        # Log10 transformation: makes differences between depth (e.g., 10 vs 1000) easier to see
        # We add +1 so that 0 depth doesn't break the math (log of 0 is impossible)
        plotmat = np.log10(dmat.fillna(0) + 1)
        cmap = sns.color_palette("magma", as_cmap=True) # Use a purple-to-yellow colour scheme
        cmap.set_bad(color="#B0B0B0") # Colour missing data grey

        plt.figure(figsize=(20, 8)) # Wide canvas for many samples
        ax = sns.heatmap(plotmat, cmap=cmap, mask=dmat.isna(), vmin=0, vmax=4)

        # Force labels for every single column so none are skipped
        ax.xaxis.set_major_locator(mticker.FixedLocator(np.arange(len(s_cols)) + 0.5))
        ax.set_xticklabels(s_cols, rotation=90, ha='center', fontsize=6)
        ax.tick_params(axis='x', which='both', length=0, pad=1)

        # Draw thin white vertical lines to separate columns visually
        for x in range(len(s_cols) + 1):
            ax.axvline(x, color="white", lw=0.3, alpha=0.6)

        # Logic to turn the label RED if the sample failed QC
        for label in ax.get_xticklabels():
            sample_on_plot = normalise_sample(label.get_text())
            if sample_on_plot in clean_failing_samples:
                label.set_color("red")
                label.set_weight("bold")

        plt.xlabel("Sample ID", fontsize=10, fontweight='bold')
        plt.ylabel("Genomic Coordinates", fontsize=10, fontweight='bold')
        plt.title(f"Coverage Heatmap: {cohort}", fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/02_heatmap_{cohort.replace(' ', '_')}.png", dpi=300)
        plt.close()

    # ==========================================================================
    # 03. Global Heatmap: Same as above but with EVERY cohort together
    # ==========================================================================
    if not all_cov.empty:
        plt.figure(figsize=(30, 12)) # Even wider for the global view
        cmap = sns.color_palette("magma", as_cmap=True)
        cmap.set_bad(color="#B0B0B0")

        sns.heatmap(np.log10(all_cov.fillna(0) + 1),
                    cmap=cmap, mask=all_cov.isna(),
                    vmin=0, vmax=4)

        ax = plt.gca() # 'Get Current Axis' to modify labels
        plt.xticks(rotation=90, fontsize=5)
        for label in ax.get_xticklabels():
            sample_on_plot = normalise_sample(label.get_text())
            if sample_on_plot in clean_failing_samples:
                label.set_color("red")
                label.set_weight("bold")

        plt.xlabel("Sample ID", fontsize=10, fontweight='bold')
        plt.ylabel("Genomic Coordinates", fontsize=10, fontweight='bold')
        plt.title("Global Coverage Heatmap", fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/03_heatmap_global.png")
        plt.close()

    # ==========================================================================
    # 05. Panel Landscape: Shows average coverage across the whole target panel
    # ==========================================================================
    plt.figure(figsize=(15, 6))
    ax = plt.gca()
    cohort_handles = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        med = df[s_cols].median(axis=1) # Calculate the median coverage per region
        # Plot the median as a line
        line, = ax.plot(range(len(med)), med, lw=2, label=cohort, color=PALETTE[cohort])
        cohort_handles.append(line)

    # Draw horizontal goal-lines using the QC_LIMITS we defined at the very start
    mean_line = QC_LIMITS['mean_cov']
    floor_line = QC_LIMITS['worst_amplicon_floor']
    h_mean = ax.axhline(mean_line, color="gold", ls=":", lw=2, label=f"Mean ≥ {mean_line}×")
    h_floor = ax.axhline(floor_line, color="red", ls="--", lw=2, label=f"Floor {floor_line}×")

    ax.set_yscale("log") # Use log scale because coverage can jump from 10 to 10,000
    ax.set_xlabel("Amplicon Index", fontweight='bold')
    ax.set_ylabel("Median Depth (×)", fontweight='bold')
    ax.set_title("Panel Landscape vs Thresholds", fontweight='bold')

    # Combine all lines into one legend
    handles = cohort_handles + [h_mean, h_floor]
    labels = [h.get_label() for h in handles]
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/05_landscape_index.png")
    plt.close()
# ==========================================================================
    # 06. Panel landscape ordered by genomic coordinate
    # ==========================================================================
    def _chr_order(c): # Helper function to sort chromosomes correctly (1, 2, ... 22, X, Y, M)
        c = str(c).replace("chr", "") # Remove "chr" prefix if it exists
        mapping = {"X": 23, "Y": 24, "M": 25, "MT": 25} # Assign numbers to non-numeric chromosomes
        try:
            return int(c) # Try to return as a number (e.g., "17" becomes 17)
        except Exception:
            return mapping.get(c, 1000) # If it's X or Y, use the map; otherwise, put at the end (1000)

    # Collect all unique target regions from all cohorts
    coord_frames = [df[["chr", "start", "end", "id"]].drop_duplicates()
                    for df in cohort_dfs.values()]
    kdf = pd.concat(coord_frames).drop_duplicates() # Combine into one master list
    kdf["key"] = kdf["chr"].apply(_chr_order) # Apply our sorting helper
    # Sort the list by chromosome number, then start position, then end position
    kdf = kdf.sort_values(["key", "start", "end"])
    ordered_ids = kdf["id"].tolist() # Create a list of IDs in the "correct" genomic order

    plt.figure(figsize=(18, 6)) # Create a wide plot
    ax = plt.gca()
    cohort_handles = []
    for cohort, df in cohort_dfs.items(): # Loop through each cohort
        s_cols = [c for c in df.columns if c not in meta]
        med = df[s_cols].median(axis=1) # Calculate median coverage for this cohort
        # Align this cohort's data to the master "genomic order" list
        t = (pd.DataFrame({"id": df["id"], "med": med})
               .drop_duplicates("id").set_index("id").reindex(ordered_ids))
        xvals = np.arange(len(ordered_ids)) # Create x-axis points
        # Plot the median coverage as a line
        line, = ax.plot(xvals, t["med"].values, lw=2, label=cohort, color=PALETTE[cohort])
        cohort_handles.append(line)

    # Add the horizontal goal-lines (Mean and Floor thresholds)
    h_mean = ax.axhline(mean_line, color="gold", ls=":", lw=2, label=f"Mean ≥ {mean_line}×")
    h_floor = ax.axhline(floor_line, color="red", ls="--", lw=2, label=f"Floor {floor_line}×")

    if ordered_ids: # Logic to decide which labels to show on the x-axis (so it doesn't get crowded)
        n = len(ordered_ids)
        step = max(1, n // 40) # Only show about 40 labels total
        sel = list(range(0, n, step))
        ax.set_xticks(sel)
        ax.set_xticklabels([ordered_ids[i] for i in sel], rotation=45, ha='right', fontsize=7)

    ax.set_yscale("log") # Use logarithmic scale for depth
    ax.set_xlabel("Amplicon ID", fontweight='bold')
    ax.set_ylabel("Median Depth (×)", fontweight='bold')
    ax.set_title("Panel Landscape (Genomic Coordinates)", fontweight='bold')

    # Rebuild the legend with the cohort lines and the threshold lines
    handles = cohort_handles + [h_mean, h_floor]
    labels = [h.get_label() for h in handles]
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/06_landscape_genomic_ids.png")
    plt.close()

    # ==========================================================================
    # 07. Panel landscape by gene (Grouping amplicons by their gene names)
    # ==========================================================================
    def parse_gene_or_annotation(a): # Helper to extract gene name from annotation text
        if pd.isna(a):
            return None
        t = str(a)
        if "\n" in t: # If the text has a newline, the gene name is usually on the second line
            return t.split("\n")[1]
        return t.strip()

    gene_frames = []
    for cohort, df in cohort_dfs.items():
        if "annotation" not in df.columns: # Skip if there are no annotations
            continue
        tmp = df[["annotation", "chr", "start", "end"]].dropna(subset=["annotation"]).copy()
        if tmp.empty:
            continue
        tmp["gene"] = tmp["annotation"].apply(parse_gene_or_annotation) # Get gene names
        tmp = tmp.dropna(subset=["gene"])
        # Group amplicons together if they belong to the same gene
        grp = tmp.groupby("gene", as_index=False).agg({"chr": "first", "start": "min", "end": "max"})
        grp["key"] = grp["chr"].apply(_chr_order) # Sort genes by chromosome position
        gene_frames.append(grp)

    if gene_frames: # If we successfully found gene names...
        gdf = pd.concat(gene_frames, ignore_index=True).drop_duplicates("gene")
        gdf = gdf.sort_values(["key", "start", "end"]) # Sort genes genomically
        gene_order = gdf["gene"].tolist()

        plt.figure(figsize=(18, 6))
        ax = plt.gca()
        cohort_handles = []
        for cohort, df in cohort_dfs.items():
            s_cols = [c for c in df.columns if c not in meta]
            med_amp = df[s_cols].median(axis=1) # Median per amplicon
            tmp = pd.DataFrame({"annotation": df["annotation"], "med": med_amp}).copy()
            tmp["gene"] = tmp["annotation"].apply(parse_gene_or_annotation)
            tmp = tmp.dropna(subset=["gene"])
            # Calculate the median coverage for the WHOLE gene (all its amplicons combined)
            gmed = tmp.groupby("gene")["med"].median().reindex(gene_order)

            xvals = np.arange(len(gene_order))
            line, = ax.plot(xvals, gmed.values, lw=2, label=cohort, color=PALETTE[cohort])
            cohort_handles.append(line)

        # Standard threshold lines
        h_mean = ax.axhline(mean_line, color="gold", ls=":", lw=2)
        h_floor = ax.axhline(floor_line, color="red", ls="--", lw=2)

        # X-axis label logic (limit to 40 labels)
        n = len(gene_order)
        step = max(1, n // 40)
        sel = list(range(0, n, step))
        ax.set_xticks(sel)
        ax.set_xticklabels([gene_order[i] for i in sel], rotation=45, ha='right', fontsize=7)

        ax.set_yscale("log")
        ax.set_xlabel("Annotated Region (Genomic Order)", fontweight='bold')
        ax.set_ylabel("Median Depth (×)", fontweight='bold')
        ax.set_title("Panel Landscape (Region)", fontweight='bold')

        handles = cohort_handles + [h_mean, h_floor]
        labels = [h.get_label() for h in handles]
        ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left')

        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/07_landscape_genes.png")
        plt.close()
    else:
        print("[INFO] 07 (genes) skipped: no gene annotations found.")
# ==========================================================================
    # 08. Coverage Gaps (<50×)
    # ==========================================================================
    gap_stats = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        for s in s_cols:
            vals = df[s]
            # Identify "gaps": regions where coverage is below the floor (50x) but not empty (NaN)
            gaps = (vals < QC_LIMITS['worst_amplicon_floor']) & (~vals.isna())
            gap_stats.append({"cohort": cohort, "sample": s, "gaps": int(gaps.sum())})

    gap_df = pd.DataFrame(gap_stats) # Create a table showing gap counts per sample

    plt.figure(figsize=(10, 7))
    # Create a boxplot to show the distribution of gaps across different cohorts
    ax = sns.boxplot(
        data=gap_df, x="cohort", y="gaps",
        hue="cohort", palette=PALETTE, legend=False
    )

    # ADVANCED LOGIC: Highlight outlier samples that failed QC in RED
    cohort_levels = list(gap_df["cohort"].unique())
    n_cohorts = len(cohort_levels)
    structural_lines = n_cohorts * 6  # Skips the lines used to draw the actual boxes

    # Loop through the "fliers" (the individual dots representing outliers on the boxplot)
    for line in ax.lines[structural_lines:]:
        if line.get_marker() not in ["o", "s", "D", "^", "v"]: # Only look at actual dots
            continue
        xs = line.get_xdata()
        ys = line.get_ydata()
        if len(xs) != 1:
            continue
        x, y = float(xs[0]), float(ys[0])
        nearest_x = int(round(x)) # Find which cohort column the dot is in
        
        if nearest_x < 0 or nearest_x >= n_cohorts:
            continue
            
        this_cohort = cohort_levels[nearest_x]
        # Match the dot back to the specific sample name
        matches = gap_df[(gap_df["cohort"] == this_cohort) & (gap_df["gaps"] == int(y))]
        # Check if any sample matching this gap count is in the 'failing_samples' list
        is_fail = any(normalise_sample(row["sample"]) in failing_samples
                      for _, row in matches.iterrows())
        
        if is_fail: # If it's a failed sample, turn the dot RED and make it larger
            line.set_color("red")
            line.set_markeredgecolor("black")
            line.set_markersize(7)

    plt.title("Coverage Gaps (<50×)", fontweight='bold')
    plt.xlabel("Cohort", fontweight='bold')
    plt.ylabel("Amplicons <50×", fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/08_coverage_gaps.png")
    plt.close()

    # ==========================================================================
    # 09. Retention Rate: Comparing Pass vs Fail counts
    # ==========================================================================
    if not full_qc_df.empty:
        plt.figure(figsize=(10, 7))
        # Group samples by cohort and status, then create a stacked bar chart
        full_qc_df.groupby(["cohort", "status"]).size().unstack().fillna(0).plot(
            kind="bar", stacked=True, color=STATUS_PALETTE, ax=plt.gca()
        )
        plt.title("Retention Rate", fontweight='bold')
        plt.xlabel("Cohort", fontweight='bold')
        plt.ylabel("Sample Count", fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/09_retention.png")
        plt.close()

    # ==========================================================================
    # 09b. Failure reasons by cohort: Diagnosing WHY samples failed
    # ==========================================================================
    if not full_qc_df.empty:
        df_rr = full_qc_df.copy()
        df_rr["fail_reasons"] = df_rr["fail_reasons"].fillna("")
        # Turn comma-separated strings (e.g., "low_reads, low_mapping") into lists
        df_rr["fail_list"] = df_rr["fail_reasons"].apply(
            lambda s: [x.strip() for x in str(s).split(",") if x.strip()]
        )
        # 'Explode' the list: if a sample has 3 reasons, it gets 3 rows in the table
        failed = df_rr[df_rr["status"] == "Fail"].explode("fail_list")
        failed["fail_list"] = failed["fail_list"].replace("", np.nan).fillna("unspecified")

        if not failed.empty:
            # Map code-like reasons to human-readable labels for the graph legend
            reason_map = {
                "low_reads": "Low Read Count",
                "low_uniformity": "Low Uniformity (<95% @ 100×)",
                "low_on_target": "Low On‑Target Rate",
                "low_mapping": "Low Mapping %",
                "low_mean_cov": "Low Mean Coverage",
                "unspecified": "Unspecified Error"
            }
            failed["fail_list"] = failed["fail_list"].map(lambda x: reason_map.get(x, x))

            # Count occurrences of each reason per cohort and pivot for plotting
            counts = (failed.groupby(["cohort", "fail_list"])
                             .size().reset_index(name="n")
                             .pivot(index="cohort", columns="fail_list", values="n")
                             .fillna(0))
            # Sort the legend so the most common reasons appear at the top
            counts = counts[counts.sum(axis=0).sort_values(ascending=False).index]

            plt.figure(figsize=(12, 7))
            counts.plot(kind="bar", stacked=True, ax=plt.gca(), cmap="tab20")
            plt.legend(title="Reason for Failure", bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.title("Failure Reasons by Cohort", fontweight='bold')
            plt.xlabel("Cohort", fontsize=10, fontweight='bold')
            plt.ylabel("Failed Samples", fontsize=10, fontweight='bold')
            plt.tight_layout()
            plt.savefig(f"{PATHS['plots_dir']}/09b_fail_reasons.png")
            plt.close()
   # ==========================================================================
    # 10. Mapping efficiency audit: Visualising how well reads match the genome
    # ==========================================================================
    if not full_qc_df.empty:
        plt.figure(figsize=(8, 6))
        # Create a 'stripplot' (a scatter plot for categories)
        sns.stripplot(
            data=full_qc_df, x="cohort", y="mapped_pct",
            hue="status", palette=STATUS_PALETTE, jitter=True # status determines colour (Pass/Fail)
        )
        # Draw a red dashed line at the 90% mapping threshold defined in config
        plt.axhline(QC_LIMITS["mapped_pct"], color="red", ls="--")
        plt.title("Mapping Efficiency Audit", fontweight='bold')
        plt.xlabel("Cohort", fontsize=10, fontweight='bold')
        plt.ylabel("% mapped", fontsize=10, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/10_mapping_audit.png")
        plt.close()

    # ==========================================================================
    # 11. Specificity jitter: Visualising "on-target" accuracy
    # ==========================================================================
    if not full_qc_df.empty:
        plt.figure(figsize=(8, 6))
        sns.stripplot(
            data=full_qc_df, x="cohort", y="on_target_pct",
            hue="status", palette=STATUS_PALETTE, jitter=True
        )
        # Threshold line for percentage of reads that actually hit the intended targets
        plt.axhline(QC_LIMITS["on_target_pct"], color="firebrick", ls="--")
        plt.title("Specificity Jitter", fontweight='bold')
        plt.xlabel("Cohort", fontsize=10, fontweight='bold')
        plt.ylabel("On target %", fontsize=10, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/11_specificity_jitter.png")
        plt.close()

    # ==========================================================================
    # 12. Lowest coverage per sample: Finding the "weakest link" in each sample
    # ==========================================================================
    worst = []
    for cohort, df in cohort_dfs.items():
        s_cols = [c for c in df.columns if c not in meta]
        for s in s_cols:
            md = df[s].min() # Find the single lowest coverage value for this sample
            if pd.notna(md):
                worst.append({"cohort": cohort, "min_depth": md})

    worst_df = pd.DataFrame(worst)
    plt.figure(figsize=(10, 6))
    # Swarmplot prevents dots from overlapping, making the density easier to see
    sns.swarmplot(
        data=worst_df, x="cohort", y="min_depth",
        hue="cohort", palette=PALETTE, size=6
    )
    plt.axhline(QC_LIMITS["worst_amplicon_floor"], color="red", ls="--")
    # 'symlog' scale: handles values near zero better than a standard log scale
    plt.yscale("symlog", linthresh=10)
    plt.title("Lowest Coverage per Sample", fontweight='bold')
    plt.xlabel("Cohort", fontsize=10, fontweight='bold')
    plt.ylabel("Minimum amplicon depth per sample", fontsize=10, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{PATHS['plots_dir']}/12_worst_amplicon.png")
    plt.close()


    # ==========================================================================
    # 13. Systemic Amplicon Failures: Identifying "troublesome" DNA regions
    # ==========================================================================
    if cohort_dfs:
        cohorts = list(cohort_dfs.keys())
        nrows = len(cohorts)
        # Make the first subplot slightly taller to accommodate more legend info
        heights = [6] + [4.0] * (nrows - 1)

        # Create a vertical stack of plots (one for each cohort)
        fig, axes = plt.subplots(
            nrows=nrows, ncols=1,
            figsize=(16, sum(heights)),
            gridspec_kw={'height_ratios': heights},
            sharex=True # All subplots share the same horizontal (count) axis
        )
        if nrows == 1:
            axes = [axes]

        max_x = 0 # Track the highest failure count to standardise the x-axis later
        for ax, cohort in zip(axes, cohorts):
            df = cohort_dfs[cohort]
            s_cols = [c for c in df.columns if c not in meta]
            if len(s_cols) == 0:
                ax.set_axis_off()
                ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)")
                continue

            floor = QC_LIMITS['worst_amplicon_floor']
            # Create a "True/False" table identifying every instance of coverage < 50x
            fail_mask = (df[s_cols] < floor) & (~df[s_cols].isna())

            # Map coordinates to the number of samples where that coordinate failed
            tmp = pd.DataFrame({
                "coord": df["id"].values,
                "fail_count": fail_mask.sum(axis=1).values
            }).dropna(subset=["coord"])

            # Sort by highest failure count and take only the Top 30 "worst" regions
            agg = (tmp.groupby("coord", as_index=False)["fail_count"]
                         .sum().sort_values("fail_count", ascending=False))
            agg = agg[agg["fail_count"] > 0].head(30)

            if agg.empty:
                ax.text(0.5, 0.5, "No amplicons below 50×", ha="center", va="center")
                ax.set_axis_off()
                ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)", fontweight='bold')
                continue

            # Create a horizontal bar chart
            colours = sns.color_palette("Reds", n_colors=len(agg))
            plot_df = agg.iloc[::-1] # Reverse order so the worst is at the top
            ax.barh(plot_df["coord"], plot_df["fail_count"],
                    color=colours[::-1], edgecolor="none")

            ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)", fontweight='bold')
            ax.set_ylabel("Genomic Coordinate (Chr17)", fontsize=10, fontweight='bold',)
            ax.grid(axis="x", linestyle=":", alpha=0.4)
            # Update the max failure count seen so far
            max_x = max(max_x, int(plot_df["fail_count"].max()))

        # Label the bottom-most x-axis and ensure all plots use the same scale
        axes[-1].set_xlabel("Number of Failed Samples", fontsize=10, fontweight='bold')
        for ax in axes:
            ax.set_xlim(0, max_x + 1) # type: ignore

        plt.tight_layout()
        plt.savefig(f"{PATHS['plots_dir']}/13_failures_coordinates.png")
        plt.close()
# ==============================================================================
# 13b. Systemic Amplicon Failures — by cohort (Annotation → extract rsID)
# ==============================================================================

    # Helper function to find the "rsID" inside a messy text string
    def extract_rs(a):
        if pd.isna(a):
            return None
        a = str(a)

        # "Tokenise" the text: replace common separators (;, | \n) with spaces
        # so we can split the text into individual words
        tokens = (
            a.replace(";", " ")
            .replace(",", " ")
            .replace("|", " ")
            .replace("\n", " ")
            .split()
        )

        # Loop through the words and return the first one that starts with "rs"
        # rsIDs look like "rs12345"
        for tok in tokens:
            if tok.lower().startswith("rs"):
                return tok

        # If no rsID is found, just return the original text as a backup
        return a

    if cohort_dfs:
        cohorts = list(cohort_dfs.keys())
        nrows = len(cohorts)

        # Set up a vertical stack of plots (standardized layout)
        heights = [6] + [4.0] * (nrows - 1)
        fig, axes = plt.subplots(
            nrows=nrows, ncols=1,
            figsize=(16, sum(heights)),
            gridspec_kw={'height_ratios': heights},
            sharex=True
        )

        if nrows == 1:
            axes = [axes]

        max_x = 0 # To keep x-axis scale consistent across all cohorts

        for ax, cohort in zip(axes, cohorts):
            df = cohort_dfs[cohort]
            s_cols = [c for c in df.columns if c not in meta]

            if len(s_cols) == 0:
                ax.set_axis_off()
                ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×, Annotation)")
                continue

            floor = QC_LIMITS['worst_amplicon_floor'] # Target is 50x

            # Identify everywhere coverage is below 50 but not missing (NaN)
            fail_mask = (df[s_cols] < floor) & (~df[s_cols].isna())

            # Create a table that links the rsID to how many times it failed
            tmp = pd.DataFrame({
                "annot": df["annotation"].apply(extract_rs).values,
                "fail_count": fail_mask.sum(axis=1).values
            }).dropna(subset=["annot"])

            # Group by rsID and sum up the failures
            agg = (
                tmp.groupby("annot", as_index=False)["fail_count"]
                .sum()
                .sort_values("fail_count", ascending=False)
            )

            # Keep only the top 30 most frequently failing SNPs
            agg = agg[agg["fail_count"] > 0].head(30)

            if agg.empty:
                ax.text(0.5, 0.5, "No amplicons below 50×",
                        ha="center", va="center")
                ax.set_axis_off()
                ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)", fontweight='bold')
                continue

            # Plotting the horizontal bar chart
            colours = sns.color_palette("Reds", n_colors=len(agg))
            plot_df = agg.iloc[::-1] # Reverse for "worst at top" look

            ax.barh(
                plot_df["annot"],
                plot_df["fail_count"],
                color=colours[::-1],
                edgecolor="none"
            )

            ax.set_title(f"{cohort} — Systemic Amplicon Failure (<50×)", fontweight='bold')
            ax.set_ylabel("rsID", fontsize=10, fontweight='bold')
            ax.grid(axis="x", linestyle=":", alpha=0.4)

            # Update the maximum x-value so we can sync the scale later
            max_x = max(max_x, int(plot_df["fail_count"].max()))

        # Set labels and ensure all charts use the same X-axis limit
        axes[-1].set_xlabel("Number of Failed Samples", fontsize=10, fontweight='bold')

        for ax in axes:
            ax.set_xlim(0, max_x + 1) # type: ignore

        plt.tight_layout()
        # Save as a separate file so you can compare rsID failures vs Coordinate failures
        plt.savefig(f"{PATHS['plots_dir']}/13b_failures_annotation.png")
        plt.close()

# ==============================================================================
# 6. MAIN EXECUTION: The final logic that runs the whole show
# ==============================================================================
def run_pipeline():
    # 1. LOAD EVERYTHING: 
    # This runs the 'build_pipeline' function to gather QC, coverage, and annotations
    full_qc_df, cohort_dfs, annotation_df, all_cov, failing_samples = build_pipeline()

    # 2. EXPORT PER-SAMPLE TABLES:
    # Instead of one giant file, this creates individual files for every single patient/sample
    for cohort, df in cohort_dfs.items():
        # Identify which columns are samples and which are metadata
        s_cols = [c for c in df.columns
                  if c not in ["chr", "start", "end", "id", "annotation", "annot_id"]]

        # Extract the genomic locations
        coords = df[["id", "chr", "start", "end"]].drop_duplicates().set_index("id")
        
        # Create a specific folder for this cohort (e.g., 'Breast_Tumour')
        cohort_dir = os.path.join(PATHS["samples_dir"], cohort.replace(" ", "_"))
        os.makedirs(cohort_dir, exist_ok=True)

        for s in s_cols:
            # Join the coordinates with the coverage depth for just THIS sample
            out = coords.join(df.set_index("id")[s].rename("depth"))
            
            # Save as a compressed (.gz) Tab-Separated Value file to save disk space
            out_path = os.path.join(cohort_dir, f"{s}.coverage.tsv.gz")
            out.to_csv(out_path, sep="\t", compression="gzip") 

    # 3. GENERATE VISUALISATIONS:
    # This triggers the massive plotting function that creates all 13+ types of graphs
    plot_all(full_qc_df, cohort_dfs, annotation_df, all_cov, failing_samples)

    print("\n[SUCCESS] Harmonised CSVs, diagnostics, and visualisations completed.")


# This 'if' statement ensures the script only runs if you open THIS file directly.
# It prevents the script from accidentally starting if it's imported by another program.
if __name__ == "__main__":
    run_pipeline()
