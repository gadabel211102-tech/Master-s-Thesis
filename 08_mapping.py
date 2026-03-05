#!/usr/bin/env python3
"""
================================================================================
Script 08 — Genomic Variant Landscape Map
================================================================================
Pipeline    : Step 08 of 17 — runs AFTER variant QC (07b) and COSMIC
              integration (07c), and BEFORE GSDMB-specific analysis (09)

--------------------------------------------------------------------------------
PURPOSE
--------------------------------------------------------------------------------
This script produces a genomic landscape map — a scatter plot showing where
variants occur along chromosome 17 within the GSDMB locus, how frequently
each variant position is observed across samples, and how variants are
distributed by gene, functional impact, cancer cohort, and tissue type.

The output is a single high-resolution figure arranged as a grid of panels:
  - Rows   = cancer cohort (e.g. Breast, Endometrium)
  - Columns = tissue type (Tumour, Healthy)

Within each panel, each dot represents a unique variant (grouped by position,
gene, impact category, cohort, and tissue). The dot's properties encode
biological information:
  - X position  → genomic coordinate on chromosome 17 (in megabases, Mb)
  - Y position  → how many samples carry that variant at that position
  - Colour      → which gene the variant falls in (GSDMB, ERBB2, GSDMA, etc.)
  - Shape       → predicted functional impact (HIGH, MODERATE, MODIFIER, LOW)

This kind of visualisation is sometimes called a "mutation landscape" or
"lollipop plot extension" and is standard in cancer genomics publications
for presenting the spatial distribution of variants across a genomic region.

--------------------------------------------------------------------------------
WHY GENOMIC POSITION IS PLOTTED IN MEGABASES (Mb)
--------------------------------------------------------------------------------
Raw genomic coordinates on chromosome 17 are large integers (e.g. 39,910,842).
Dividing by 1,000,000 converts these to megabases (e.g. 39.91 Mb), which are
the conventional unit used in genomics figures for readability. The GSDMB
locus spans approximately 39.9–40.0 Mb on chromosome 17q21.32.

--------------------------------------------------------------------------------
UNDERSTANDING THE IMPACT CATEGORIES (VEP)
--------------------------------------------------------------------------------
Each variant is annotated with a predicted functional impact by VEP (Variant
Effect Predictor). The four categories, from most to least severe, are:

  HIGH     (✗ marker) — likely causes loss of protein function.
                        Examples: frameshift insertion/deletion, premature
                        stop codon, splice site disruption.

  MODERATE (● marker) — likely changes protein function but does not
                        completely abolish it.
                        Examples: missense substitution (amino acid change),
                        in-frame insertion or deletion.

  MODIFIER (■ marker) — falls in a non-coding region or is otherwise unlikely
                        to directly affect the protein product.
                        Examples: intronic variant, upstream/downstream gene
                        variant, 3′ or 5′ UTR variant.

  LOW      (▼ marker) — predicted to have little or no functional consequence.
                        Examples: synonymous substitution (same amino acid),
                        stop codon retained variant.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        Required columns: Pos, Symbol, Impact, Cohort, Tissue

--------------------------------------------------------------------------------
OUTPUT
--------------------------------------------------------------------------------
  clean_landscape_all_impacts.png
    Multi-panel scatter figure (300 dpi), one panel per cohort × tissue
    combination. Suitable for direct inclusion in a thesis or publication.

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python 08_mapping.py

  (Paths are configured in the CONFIGURATION block below.)

--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------
  Python ≥ 3.8
  pandas, matplotlib, seaborn, openpyxl

  Install: pip install pandas matplotlib seaborn openpyxl
================================================================================
"""

# ──────────────────────────────────────────────────────────────────────────────
# STANDARD LIBRARY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import os  # File existence check and path construction

# ──────────────────────────────────────────────────────────────────────────────
# THIRD-PARTY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt  # Figure creation and axis formatting
import pandas as pd              # Data loading and aggregation
import seaborn as sns            # High-level plotting (relplot grid)


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
# Update BASE_PATH to match your working directory before running.

input_file   = "/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_Annotated_Report_Fixed.xlsx"
output_image = "/home/gadeaalonsoj/tfm/gsdmb_final_results/clean_landscape_all_impacts.png"


# ══════════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def find_col(df: pd.DataFrame, target: str):
    """
    Case-insensitive column name lookup.

    Excel files exported from different tools may capitalise column names
    inconsistently (e.g. "SYMBOL" vs "Symbol" vs "symbol"). This function
    normalises the comparison so the script works regardless of capitalisation.

    Parameters
    ----------
    df     : pd.DataFrame — the DataFrame to search
    target : str          — the desired column name (any case)

    Returns
    -------
    str or None
        The actual column name as it appears in the DataFrame, or None if the
        column is not found.
    """
    for col in df.columns:
        if col.upper() == target.upper():
            return col
    return None


# ══════════════════════════════════════════════════════════════════════════════
# MAIN FUNCTION: GENOMIC LANDSCAPE MAP
# ══════════════════════════════════════════════════════════════════════════════

def generate_clean_all_impact_map():
    """
    Generate and save the genomic variant landscape figure.

    Execution steps
    ---------------
    1. Load the annotated variant table from Excel.
    2. Identify required columns (case-insensitively).
    3. Normalise the IMPACT column (strip whitespace, uppercase) to ensure
       all four impact categories are detected correctly.
    4. Aggregate: count how many times each unique (position, gene, impact,
       cohort, tissue) combination appears across samples. This "frequency"
       becomes the Y-axis value, representing how many samples share that
       variant.
    5. Build the figure using seaborn's relplot, which automatically creates
       a grid of subplots (facets) — one per cohort × tissue combination.
    6. Format axes: convert raw genomic coordinates to megabases (Mb).
    7. Save as a 300 dpi PNG.

    Design decisions
    ----------------
    - sharex=True ensures all panels share the same genomic X-axis range,
      making it easy to visually compare positions across cohorts and tissues.
    - sharey=False allows each panel to have its own Y-axis scale, because
      variant frequencies may differ substantially between tissue types.
    - The "husl" colour palette is used for gene colours because it generates
      maximally distinct, perceptually uniform hues — important when multiple
      closely located genes need to be visually separated.
    - Marker shapes are explicitly mapped to impact categories rather than
      using seaborn defaults, ensuring consistent and intuitive encoding:
        HIGH (×) draws the eye immediately; MODIFIER (■) is visually distinct
        from the primary MODERATE (●) category.
    """
    print("=" * 65)
    print("SCRIPT 08 — Genomic Variant Landscape Map")
    print("=" * 65)
    print(f"Input : {input_file}")
    print(f"Output: {output_image}\n")

    # ── Guard: check input file exists before attempting to load ─────────────
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        print("       Please update 'input_file' in the CONFIGURATION block.")
        return

    # ── Step 1: Load annotated variant table ─────────────────────────────────
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")
    print(f"Loaded {len(df)} rows from sheet 'Biological_Annotations'")

    # ── Step 2: Identify required columns (case-insensitive) ─────────────────
    pos_c = find_col(df, "Pos")      # Genomic position (integer, e.g. 39910842)
    sym_c = find_col(df, "Symbol")   # Gene symbol (e.g. "GSDMB", "ERBB2")
    imp_c = find_col(df, "Impact")   # VEP impact category (HIGH/MODERATE/MODIFIER/LOW)
    coh_c = find_col(df, "Cohort")   # Cancer cohort (e.g. "Breast", "Endometrium")
    tis_c = find_col(df, "Tissue")   # Tissue type ("Tumour" or "Healthy")

    # ── Step 3: Normalise the IMPACT column ───────────────────────────────────
    # Strip leading/trailing whitespace and convert to uppercase.
    # This prevents subtle mismatches where, for example, " High" (with a
    # leading space) would be treated as a separate category from "HIGH",
    # causing some impact types to appear missing from the plot.
    df[imp_c] = df[imp_c].astype(str).str.strip().str.upper()

    # ── Step 4: Aggregate variant frequency ───────────────────────────────────
    # Group by all five identifying columns and count the number of rows in
    # each group. Each row in the annotated table corresponds to one variant
    # in one sample (after VEP expansion), so the count gives the number of
    # samples sharing that (position, gene, impact, cohort, tissue) combination.
    group_cols     = [pos_c, sym_c, imp_c, coh_c, tis_c]
    variant_counts = df.groupby(group_cols).size().reset_index(name="Frequency")
    print(f"Unique (position × gene × impact × cohort × tissue) groups: {len(variant_counts)}")

    # ── Step 5: Define visual encoding ───────────────────────────────────────
    # Map each VEP impact category to a distinct marker shape.
    # Shapes were chosen to be easily distinguishable in both colour and
    # black-and-white printing:
    #   X (cross)            → HIGH impact   — visually prominent, signals importance
    #   ● (circle)           → MODERATE      — standard default marker
    #   ■ (square)           → MODIFIER      — clearly distinct from circle
    #   ▼ (triangle down)    → LOW           — indicates lower importance
    marker_map = {
        "HIGH":     "X",   # Bold cross
        "MODERATE": "o",   # Circle
        "MODIFIER": "s",   # Square
        "LOW":      "v",   # Downward-pointing triangle
    }

    # Convert the impact column to an ordered categorical type.
    # This ensures seaborn renders the legend entries in the intended order
    # (HIGH first, LOW last) rather than alphabetically.
    variant_counts[imp_c] = pd.Categorical(
        variant_counts[imp_c],
        categories=["HIGH", "MODERATE", "MODIFIER", "LOW"],
        ordered=True,
    )

    # ── Step 6: Build the faceted scatter plot ────────────────────────────────
    # seaborn's relplot creates a grid of axes automatically:
    #   row=coh_c  → one row of panels per cancer cohort
    #   col=tis_c  → one column of panels per tissue type
    # Within each panel, dots are coloured by gene (hue) and shaped by impact (style).
    sns.set_style("whitegrid")

    g = sns.relplot(
        data     = variant_counts,
        x        = pos_c,          # Genomic position on X-axis
        y        = "Frequency",    # Sample count on Y-axis
        hue      = sym_c,          # Colour encodes gene identity
        style    = imp_c,          # Shape encodes VEP impact category
        row      = coh_c,          # Grid rows = cancer cohorts
        col      = tis_c,          # Grid columns = tissue types
        kind     = "scatter",
        markers  = marker_map,     # Use our explicitly defined shape mapping
        s        = 120,            # Dot size (points²) — large enough to read
        alpha    = 0.7,            # Slight transparency to reveal overlapping dots
        edgecolor= "black",        # Thin black outline improves dot visibility
        palette  = "husl",         # Maximally distinct, perceptually uniform colours
        height   = 5,              # Height of each individual panel (inches)
        aspect   = 1.6,            # Width = height × aspect (landscape proportions)
        facet_kws= {
            "sharex": True,    # All panels share the same genomic X range
            "sharey": False,   # Each panel has its own Y scale
        },
    )

    # ── Step 7: Axis formatting ───────────────────────────────────────────────
    # Set axis labels shared across all panels
    g.set_axis_labels(
        "Genomic Position on Chr17 (Mb)",
        "Number of Samples",
    )
    # Panel titles: show cohort and tissue type for each facet
    g.set_titles("{row_name} | {col_name}", fontweight="bold")

    # Convert raw genomic coordinates (integers) to megabases for readability.
    # Lambda function: divides each tick value by 1,000,000 and formats to
    # two decimal places (e.g. 39910842 → "39.91 Mb").
    for ax in g.axes.flat:
        ax.xaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: f"{x / 1e6:.2f}")
        )

    # Add overall figure title
    plt.subplots_adjust(top=0.92)
    g.fig.suptitle(
        "Landscape of Variants Across the GSDMB Locus (Chr17q21.32)",
        fontsize=16,
        fontweight="bold",
    )

    # ── Step 8: Save figure ───────────────────────────────────────────────────
    plt.savefig(output_image, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"\nFigure saved: {output_image}")
    print("=" * 65)
    print("Script 08 complete.")
    print("=" * 65)


# ── Script entry point ─────────────────────────────────────────────────────────
# Ensures generate_clean_all_impact_map() is only called when this script is
# run directly, not when imported as a module by another script.
if __name__ == "__main__":
    generate_clean_all_impact_map()
