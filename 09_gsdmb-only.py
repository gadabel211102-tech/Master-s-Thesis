#!/usr/bin/env python3
"""
================================================================================
Script 09 — GSDMB-Exclusive Mutation Landscape Map
================================================================================
Pipeline    : Step 09 of 17 — runs AFTER the full-locus landscape map (08)
              and BEFORE variant statistics (10)

--------------------------------------------------------------------------------
PURPOSE
--------------------------------------------------------------------------------
This script produces a focused mutation landscape map showing only variants
that fall within the GSDMB gene itself — as opposed to script 08, which plots
all genes across the entire GSDMB locus (GSDMB, ORMDL3, GSDMA, ERBB2, etc.).

By restricting to GSDMB alone, this figure allows a more detailed examination
of the internal variant distribution within the gene: which positions are most
frequently mutated, how impact categories are distributed along the gene body,
and whether mutation patterns differ between cancer cohorts and tissue types.

--------------------------------------------------------------------------------
DIFFERENCE FROM SCRIPT 08
--------------------------------------------------------------------------------
Script 08 (full locus):
  - All genes in the region plotted together
  - Colour encodes gene identity (each gene a different colour)
  - X-axis spans the full locus (~39.9–40.0 Mb), displayed in megabases
  - sharey=False: each panel has its own Y-axis scale

Script 09 (GSDMB only):
  - Only GSDMB variants plotted
  - Colour encodes functional impact category (more informative when a single
    gene is shown, since gene identity is no longer variable)
  - X-axis shows full base-pair precision with comma-formatted coordinates
    (e.g. 39,910,842) to allow identification of individual variant positions
  - sharey=True: all panels share the same Y-axis scale, enabling direct
    comparison of variant frequencies across cohorts and tissue types

--------------------------------------------------------------------------------
WHY BASE-PAIR PRECISION INSTEAD OF MEGABASES
--------------------------------------------------------------------------------
In script 08, the wide genomic range (~100 kb) made megabase (Mb) units
necessary for readability. Here, because we are zoomed into a single gene
spanning only ~15 kb, full base-pair coordinates are used instead. This allows
individual variant positions to be identified and cross-referenced with the
GSDMB transcript structure or external databases (e.g. Ensembl, gnomAD).
Coordinates are formatted with comma separators (e.g. 39,910,842) to remain
readable at full precision.

--------------------------------------------------------------------------------
UNDERSTANDING THE COLOUR / SHAPE ENCODING
--------------------------------------------------------------------------------
Since only one gene (GSDMB) is shown, colour is reassigned from gene identity
to VEP impact category. This makes the impact distribution the primary visual
signal within each panel. Both colour (hue) and shape (style) encode the same
variable (impact), which provides redundant encoding — a visualisation best
practice that makes the figure accessible to colour-blind readers.

  HIGH     (warm red,  ✗) — likely loss of protein function
  MODERATE (amber,     ●) — predicted change in protein function
  MODIFIER (light,     ■) — non-coding or negligible protein effect
  LOW      (pale,      ▼) — synonymous or otherwise benign

The "flare" palette (warm tones from pale yellow to deep red) naturally maps
to the severity ordering, reinforcing the visual hierarchy from LOW to HIGH.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        Required columns: Pos, Symbol, Impact, Cohort, Tissue

--------------------------------------------------------------------------------
OUTPUT
--------------------------------------------------------------------------------
  GSDMB_exclusive_landscape.png
    Multi-panel scatter figure (300 dpi), one panel per cohort × tissue
    combination, restricted to GSDMB variants only.

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python 09_gsdmb-only.py

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

import os  # File existence check

# ──────────────────────────────────────────────────────────────────────────────
# THIRD-PARTY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt  # Figure creation and axis formatting
import pandas as pd              # Data loading and aggregation
import seaborn as sns            # High-level faceted scatter plot


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
# Update these paths to match your working directory before running.

input_file   = "/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_Annotated_Report_Fixed.xlsx"
output_image = "/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_exclusive_landscape.png"


# ══════════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def find_col(df: pd.DataFrame, target: str):
    """
    Case-insensitive column name lookup.

    Excel files exported from different tools may capitalise column names
    inconsistently. This function normalises the comparison so the script
    works regardless of capitalisation.

    Parameters
    ----------
    df     : pd.DataFrame — the DataFrame to search
    target : str          — the desired column name (any case)

    Returns
    -------
    str or None
        The actual column name as it appears in the DataFrame, or None if
        the column is not found.
    """
    for col in df.columns:
        if col.upper() == target.upper():
            return col
    return None


# ══════════════════════════════════════════════════════════════════════════════
# MAIN FUNCTION: GSDMB-EXCLUSIVE LANDSCAPE MAP
# ══════════════════════════════════════════════════════════════════════════════

def generate_gsdmb_map():
    """
    Generate and save the GSDMB-exclusive mutation landscape figure.

    Execution steps
    ---------------
    1. Load the annotated variant table from Excel.
    2. Identify required columns (case-insensitively).
    3. Filter the table to retain only rows where the gene symbol contains
       "GSDMB". The .str.contains() approach is used (rather than exact
       equality) to catch any minor naming variations in the source file
       (e.g. "GSDMB" vs "gsdmb" vs "GSDMB ").
    4. Normalise the IMPACT column (strip whitespace, uppercase).
    5. Aggregate: count how many samples share each (position, impact,
       cohort, tissue) combination — this becomes the Y-axis frequency.
    6. Build the faceted scatter plot using seaborn relplot.
    7. Format the X-axis using full base-pair coordinates with comma
       separators, and rotate tick labels for readability.
    8. Save as a 300 dpi PNG.
    """
    print("=" * 65)
    print("SCRIPT 09 — GSDMB-Exclusive Mutation Landscape Map")
    print("=" * 65)
    print(f"Input : {input_file}")
    print(f"Output: {output_image}\n")

    # ── Guard: check input file exists ───────────────────────────────────────
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        print("       Please update 'input_file' in the CONFIGURATION block.")
        return

    # ── Step 1: Load annotated variant table ─────────────────────────────────
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")
    print(f"Loaded {len(df)} rows from sheet 'Biological_Annotations'")

    # ── Step 2: Identify required columns (case-insensitive) ─────────────────
    pos_c = find_col(df, "Pos")      # Genomic position (integer base-pair coordinate)
    sym_c = find_col(df, "Symbol")   # Gene symbol
    imp_c = find_col(df, "Impact")   # VEP functional impact category
    coh_c = find_col(df, "Cohort")   # Cancer cohort (Breast / Endometrium)
    tis_c = find_col(df, "Tissue")   # Tissue type (Tumour / Healthy)

    # ── Step 3: Filter to GSDMB variants only ────────────────────────────────
    # str.contains('GSDMB', case=False) is used rather than exact equality to
    # handle minor variations in gene symbol formatting in the source file.
    df_gsdmb = df[df[sym_c].astype(str).str.contains("GSDMB", case=False)].copy()

    if df_gsdmb.empty:
        print("WARNING: No GSDMB variants found in the dataset.")
        print("         Check that the Symbol column contains 'GSDMB' entries.")
        return

    print(f"GSDMB variants retained: {len(df_gsdmb)} rows "
          f"(from {df_gsdmb[pos_c].nunique()} unique positions)")

    # ── Step 4: Normalise the IMPACT column ───────────────────────────────────
    # Strip whitespace and convert to uppercase to prevent category mismatches
    # caused by formatting inconsistencies in the source Excel file.
    df_gsdmb[imp_c] = df_gsdmb[imp_c].astype(str).str.strip().str.upper()

    # ── Step 5: Aggregate variant frequency ───────────────────────────────────
    # For each unique (position, impact, cohort, tissue) combination, count
    # how many rows exist. Each row corresponds to one variant in one sample
    # (after VEP transcript expansion and deduplication), so this count
    # reflects how many samples share that variant at that position.
    variant_counts = (
        df_gsdmb
        .groupby([pos_c, imp_c, coh_c, tis_c])
        .size()
        .reset_index(name="Frequency")
    )
    print(f"Unique (position × impact × cohort × tissue) groups: {len(variant_counts)}")

    # ── Step 6: Define visual encoding ───────────────────────────────────────
    # Shape marker mapping — same convention as script 08 for consistency
    # across the figure set in the thesis:
    marker_map = {
        "HIGH":     "X",   # Bold cross — visually prominent, signals high severity
        "MODERATE": "o",   # Circle — standard default marker
        "MODIFIER": "s",   # Square — distinct from circle
        "LOW":      "v",   # Downward triangle — visually de-emphasised
    }

    # ── Step 7: Build the faceted scatter plot ────────────────────────────────
    # Both hue (colour) and style (shape) encode the same variable (impact).
    # This "redundant encoding" ensures the plot is interpretable even when
    # printed in greyscale or viewed by colour-blind readers.
    sns.set_style("whitegrid")

    g = sns.relplot(
        data      = variant_counts,
        x         = pos_c,         # Genomic position (base pairs) on X-axis
        y         = "Frequency",   # Sample count on Y-axis
        hue       = imp_c,         # Colour encodes impact category
        style     = imp_c,         # Shape also encodes impact (redundant encoding)
        row       = coh_c,         # Grid rows = cancer cohorts
        col       = tis_c,         # Grid columns = tissue types
        kind      = "scatter",
        markers   = marker_map,
        palette   = "flare",       # Warm gradient palette: pale (LOW) → deep red (HIGH)
        s         = 150,           # Dot size — slightly larger than script 08 since
                                   # fewer points are plotted (GSDMB only)
        alpha     = 0.8,           # Slight transparency for overlapping points
        edgecolor = "black",       # Thin outline improves dot visibility on pale colours
        height    = 5,             # Height per panel (inches)
        aspect    = 1.4,           # Width = height × aspect
        facet_kws = {
            "sharex": True,    # All panels share the same X range (same gene)
            "sharey": True,    # All panels share the same Y scale — enables direct
                               # comparison of frequencies across cohorts/tissues
        },
    )

    # ── Step 8: Axis formatting ───────────────────────────────────────────────
    g.set_axis_labels("Genomic Position on Chr17 (bp)", "Number of Samples")
    g.set_titles("{row_name} | {col_name}", fontweight="bold")

    # Format X-axis with full base-pair precision and comma separators.
    # Example: 39910842 → "39,910,842"
    # Tick labels are rotated 45° to prevent overlap at this level of precision.
    for ax in g.axes.flat:
        ax.xaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: format(int(x), ","))
        )
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(
        "GSDMB Mutation Landscape: Variant Distribution by Cohort and Tissue Type",
        fontsize=16,
        fontweight="bold",
    )

    # ── Step 9: Save figure ───────────────────────────────────────────────────
    plt.savefig(output_image, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"\nFigure saved: {output_image}")
    print("=" * 65)
    print("Script 09 complete.")
    print("=" * 65)


# ── Script entry point ─────────────────────────────────────────────────────────
# Ensures generate_gsdmb_map() is only called when this script is run
# directly, not when imported as a module by another script.
if __name__ == "__main__":
    generate_gsdmb_map()
