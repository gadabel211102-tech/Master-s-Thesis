#!/usr/bin/env python3
"""
================================================================================
Script 10 — Variant Statistics: Impact and Consequence Distributions
================================================================================
Pipeline    : Step 10 of 17 — runs AFTER the GSDMB-exclusive landscape map (09)
              and BEFORE SNP analysis (11)

--------------------------------------------------------------------------------
PURPOSE
--------------------------------------------------------------------------------
This script generates descriptive statistics and bar chart figures summarising
the functional impact and molecular consequence of all variants detected across
the GSDMB locus, broken down by cohort and tissue type.

Two parallel analyses are performed:

  1. Total variants
     All rows in the annotated table are included. Because VEP expands each
     variant across all transcripts it overlaps, the same variant may appear
     multiple times (once per transcript). This "total" count therefore
     reflects the full annotation burden — useful for understanding how many
     gene-transcript combinations are affected.

  2. Unique variants
     Rows are deduplicated so that each unique combination of gene symbol,
     genomic position, protein-level change (HGVSp), and sample group is
     counted only once. This gives a less inflated view of variant diversity,
     closer to the true number of distinct mutational events.

For each of the two analyses, four output files are produced:
  - A log-scale bar chart figure  (useful when counts span several orders of
    magnitude, e.g. MODIFIER variants vastly outnumber HIGH impact variants)
  - A linear-scale bar chart figure (easier to read absolute differences)
  - A CSV table of impact counts per group
  - A CSV table of consequence counts per group

--------------------------------------------------------------------------------
GROUPING STRATEGY
--------------------------------------------------------------------------------
Variants are grouped by the combination of cohort and tissue type, producing
labels such as "Breast - Tumour", "Breast - Healthy", "Endometrium - Tumour",
"Endometrium - Healthy". This grouping allows simultaneous comparison across
cancer types (cohort) and between tumour and matched healthy tissue within
each cancer type.

--------------------------------------------------------------------------------
UNDERSTANDING THE TWO FIGURE PANELS
--------------------------------------------------------------------------------
Each figure contains two stacked bar chart panels:

  Top panel — Functional Consequences
    Shows the distribution of specific molecular consequence types as
    annotated by VEP (e.g. missense_variant, synonymous_variant, intron_variant,
    3_prime_UTR_variant, splice_region_variant, etc.). There are many possible
    consequence types, so a tab20 colour palette is used for maximum
    distinction. Numeric labels are printed to the right of each bar segment
    for precise reading.

  Bottom panel — Biological Impact Level
    Shows the four high-level VEP impact categories (HIGH, MODERATE, MODIFIER,
    LOW) aggregated from the detailed consequences above. Uses a fixed colour
    scheme consistent across all scripts in this pipeline:
      HIGH     → red     (#d9534f)
      MODERATE → blue    (#5bc0de)
      MODIFIER → amber   (#f0ad4e)
      LOW      → green   (#5cb85c)

--------------------------------------------------------------------------------
LOG SCALE vs LINEAR SCALE
--------------------------------------------------------------------------------
Both scale versions are saved because they reveal different aspects of the data:

  Log scale:
    Compresses the range, making rare categories (e.g. HIGH impact) visible
    even when MODIFIER variants dominate numerically. Useful for checking
    whether any HIGH impact variants are present at all.

  Linear scale:
    Preserves the true proportional differences between categories.
    More intuitive for general audiences and thesis figures where the
    absolute dominance of non-coding / MODIFIER variants is the key message.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        Required columns: Symbol, Pos, HGVSp, Tissue, Cohort, Impact,
                          Consequence

--------------------------------------------------------------------------------
OUTPUT FILES (saved to input_dir)
--------------------------------------------------------------------------------
  Figures (PNG, 300 dpi):
    01_ALL_VARIANTS_TOTAL_LOG.png     — total counts, log scale
    02_ALL_VARIANTS_TOTAL_LINEAR.png  — total counts, linear scale
    03_ALL_VARIANTS_UNIQUE_LOG.png    — unique counts, log scale
    04_ALL_VARIANTS_UNIQUE_LINEAR.png — unique counts, linear scale

  Tables (CSV):
    table_impact_all_total.csv        — impact counts, total variants
    table_consequence_all_total.csv   — consequence counts, total variants
    table_impact_all_unique.csv       — impact counts, unique variants
    table_consequence_all_unique.csv  — consequence counts, unique variants

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python 10_variant-stats.py

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

import os  # File path construction and existence checks

# ──────────────────────────────────────────────────────────────────────────────
# THIRD-PARTY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt  # Figure creation, axis formatting, text labels
import pandas as pd              # Data loading, grouping, pivot tables
import seaborn as sns            # Plot styling


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
# Update input_dir to match your working directory before running.

input_dir  = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
input_file = os.path.join(input_dir, "GSDMB_Annotated_Report_Fixed.xlsx")

# Output filenames for the four bar chart figures.
# Numbers prefix the filenames so they sort in logical order in the file browser.
outputs = {
    "total_log":  "01_ALL_VARIANTS_TOTAL_LOG.png",
    "total_lin":  "02_ALL_VARIANTS_TOTAL_LINEAR.png",
    "unique_log": "03_ALL_VARIANTS_UNIQUE_LOG.png",
    "unique_lin": "04_ALL_VARIANTS_UNIQUE_LINEAR.png",
}

# Fixed colour scheme for the four VEP impact categories — consistent across
# all figures in this pipeline.
IMPACT_COLOURS = {
    "HIGH":     "#d9534f",  # Red     — loss of function / high severity
    "MODERATE": "#5bc0de",  # Blue    — functional change
    "MODIFIER": "#f0ad4e",  # Amber   — non-coding / low direct effect
    "LOW":      "#5cb85c",  # Green   — benign / synonymous
}

# Ordered list of impact categories (most to least severe).
# Used to ensure consistent column ordering in pivot tables and figures.
IMPACT_ORDER = ["HIGH", "MODERATE", "MODIFIER", "LOW"]


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
        if col.strip().upper() == target.upper():
            return col
    return None


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION: EXPORT SUMMARY TABLES
# ══════════════════════════════════════════════════════════════════════════════

def generate_and_save_tables(df: pd.DataFrame, suffix: str, output_dir: str):
    """
    Export impact and consequence count tables as CSV files.

    Two tables are produced per call:

      Impact table
        Rows = cohort-tissue groups (e.g. "Breast - Tumour").
        Columns = the four VEP impact categories (HIGH, MODERATE, MODIFIER, LOW).
        Values = number of variant rows in each group × impact combination.
        Columns are fixed in severity order (HIGH first) for consistent
        presentation regardless of which categories happen to be present.

      Consequence table
        Rows = cohort-tissue groups.
        Columns = all VEP consequence types present in the data (e.g.
                  missense_variant, intron_variant, 3_prime_UTR_variant, etc.).
        Values = number of variant rows per group × consequence combination.

    Parameters
    ----------
    df         : pd.DataFrame — variant table (must contain 'Group' column,
                                added in main() before calling this function)
    suffix     : str          — label appended to output filenames to distinguish
                                "all_total" from "all_unique" outputs
    output_dir : str          — directory where CSV files will be saved
    """
    imp_c = find_col(df, "Impact")
    con_c = find_col(df, "Consequence")

    print(f"\n  Exporting tables for: {suffix.upper()}")

    # ── Impact table ──────────────────────────────────────────────────────────
    impact_table = df.groupby(["Group", imp_c]).size().unstack(fill_value=0)

    # Ensure all four impact categories are present as columns, even if no
    # variants of that type were found (adds zero-filled columns for completeness)
    for cat in IMPACT_ORDER:
        if cat not in impact_table.columns:
            impact_table[cat] = 0
    impact_table = impact_table[IMPACT_ORDER]  # Enforce severity ordering

    impact_table.to_csv(os.path.join(output_dir, f"table_impact_{suffix}.csv"))
    print(f"    Saved: table_impact_{suffix}.csv")

    # ── Consequence table ─────────────────────────────────────────────────────
    con_table = df.groupby(["Group", con_c]).size().unstack(fill_value=0)
    con_table.to_csv(os.path.join(output_dir, f"table_consequence_{suffix}.csv"))
    print(f"    Saved: table_consequence_{suffix}.csv")


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION: GENERATE BAR CHART FIGURE
# ══════════════════════════════════════════════════════════════════════════════

def create_stat_visuals(
    df_to_plot: pd.DataFrame,
    title_suffix: str,
    filename: str,
    use_log: bool = False,
):
    """
    Create and save a two-panel stacked bar chart figure.

    The figure has a vertical layout with two panels:
      - Top panel:    Functional Consequences (detailed VEP consequence types)
      - Bottom panel: Biological Impact Level (HIGH / MODERATE / MODIFIER / LOW)

    Numeric count labels are printed to the right of each bar segment,
    providing precise values without requiring the reader to measure against
    the axis grid.

    Parameters
    ----------
    df_to_plot   : pd.DataFrame — variant table to visualise (must contain
                                  'Group', Impact, and Consequence columns)
    title_suffix : str          — appended to each panel title (e.g. "All Variants"
                                  or "Unique Variants")
    filename     : str          — output PNG filename (saved to input_dir)
    use_log      : bool         — if True, Y-axis is displayed on a log scale;
                                  if False, linear scale is used
    """
    imp_c = find_col(df_to_plot, "Impact")
    con_c = find_col(df_to_plot, "Consequence")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 15))
    sns.set_style("whitegrid")

    # ── Top panel: Functional Consequences ───────────────────────────────────
    # Pivot: rows = groups, columns = consequence types, values = counts
    con_pivot = df_to_plot.groupby(["Group", con_c]).size().unstack(fill_value=0)

    # tab20 provides 20 distinct colours — suitable for the many VEP consequence
    # types (missense_variant, intron_variant, 3_prime_UTR_variant, etc.)
    con_pivot.plot(kind="bar", stacked=True, ax=ax1, colormap="tab20")

    ax1.set_title(f"Functional Consequences: {title_suffix}",
                  fontsize=16, fontweight="bold")
    ax1.set_ylabel("Count (Log Scale)" if use_log else "Count")
    ax1.set_xlabel("")

    if use_log:
        ax1.set_yscale("log")

    # Rotate x-axis group labels for readability
    plt.setp(ax1.get_xticklabels(), rotation=0)

    # Move legend outside the plot area to avoid overlapping bars
    ax1.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize="small",
               title="Consequence Type")

    # Add count labels to the right of each bar segment.
    # For each coloured segment (patch), compute its centre position and
    # print the count value. Only segments with count > 0 are labelled
    # to avoid cluttering the figure with zeros.
    for patch in ax1.patches:
        height = patch.get_height()
        if height > 0:
            x = patch.get_x() + patch.get_width() + 0.02  # Slightly right of bar edge
            y = patch.get_y() + height / 2                 # Vertically centred on segment
            ax1.text(x, y, f"{int(height)}", ha="left", va="center",
                     fontsize=7, fontweight="bold", clip_on=False)

    # ── Bottom panel: Biological Impact Level ────────────────────────────────
    # Pivot: rows = groups, columns = impact categories, values = counts
    impact_pivot = df_to_plot.groupby(["Group", imp_c]).size().unstack(fill_value=0)

    # Ensure all four categories exist and enforce severity ordering
    for cat in IMPACT_ORDER:
        if cat not in impact_pivot.columns:
            impact_pivot[cat] = 0
    impact_pivot = impact_pivot[IMPACT_ORDER]

    # Use the fixed impact colour scheme defined in IMPACT_COLOURS
    bar_colors = [IMPACT_COLOURS[cat] for cat in IMPACT_ORDER]
    impact_pivot.plot(kind="bar", stacked=True, ax=ax2, color=bar_colors)

    ax2.set_title(f"Biological Impact Level: {title_suffix}",
                  fontsize=16, fontweight="bold")
    ax2.set_ylabel("Count (Log Scale)" if use_log else "Count")
    ax2.set_xlabel("")

    if use_log:
        ax2.set_yscale("log")

    plt.setp(ax2.get_xticklabels(), rotation=0)
    ax2.legend(title="Impact Category", loc="upper right")

    # Add count labels to the right of each impact bar segment
    for patch in ax2.patches:
        height = patch.get_height()
        if height > 0:
            x = patch.get_x() + patch.get_width() + 0.02
            y = patch.get_y() + height / 2
            ax2.text(x, y, f"{int(height)}", ha="left", va="center",
                     fontsize=8, fontweight="bold", clip_on=False)

    plt.tight_layout()
    plt.savefig(os.path.join(input_dir, filename), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"    Saved: {filename}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """
    Orchestrate the full variant statistics pipeline.

    Execution order
    ---------------
    1. Load the annotated variant table from Excel.
    2. Identify required columns (case-insensitively).
    3. Construct a 'Group' label combining cohort and tissue type for each row
       (e.g. "Breast - Tumour"), used as the X-axis grouping in all figures.
    4. Analysis A — Total variants:
         a. Export impact and consequence CSV tables.
         b. Save log-scale bar chart figure.
         c. Save linear-scale bar chart figure.
    5. Analysis B — Unique variants:
         a. Deduplicate to one row per (gene, position, protein change, group).
         b. Export impact and consequence CSV tables.
         c. Save log-scale bar chart figure.
         d. Save linear-scale bar chart figure.
    6. Print a completion summary.
    """
    print("=" * 65)
    print("SCRIPT 10 — Variant Statistics: Impact and Consequence Distributions")
    print("=" * 65)
    print(f"Input : {input_file}")
    print(f"Output: {input_dir}\n")

    # ── Guard: check input file exists ───────────────────────────────────────
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        print("       Please update 'input_dir' in the CONFIGURATION block.")
        return

    # ── Step 1: Load annotated variant table ─────────────────────────────────
    df = pd.read_excel(input_file, sheet_name="Biological_Annotations")
    print(f"Loaded {len(df)} rows from sheet 'Biological_Annotations'")

    # ── Step 2: Identify required columns ────────────────────────────────────
    sym_c = find_col(df, "Symbol")      # Gene symbol
    tis_c = find_col(df, "Tissue")      # Tissue type (Tumour / Healthy)
    coh_c = find_col(df, "Cohort")      # Cancer cohort (Breast / Endometrium)
    pos_c = find_col(df, "Pos")         # Genomic position
    hgv_c = find_col(df, "HGVSp")      # Protein-level variant notation
                                        # (e.g. p.Arg175His) — used for
                                        # deduplication in the unique analysis

    # ── Step 3: Construct the grouping label ─────────────────────────────────
    # Combine cohort and tissue into a single readable label used as the X-axis
    # category in all bar charts. Strip whitespace first to avoid inconsistent
    # labels caused by trailing spaces in the source data.
    df[tis_c] = df[tis_c].astype(str).str.strip()
    df[coh_c] = df[coh_c].astype(str).str.strip()
    df["Group"] = df[coh_c] + " - " + df[tis_c]

    unique_groups = sorted(df["Group"].unique())
    print(f"Groups identified: {unique_groups}\n")

    # ══════════════════════════════════════════════════════════════════════════
    # ANALYSIS A: TOTAL VARIANTS
    # (All rows including transcript-level VEP expansions)
    # ══════════════════════════════════════════════════════════════════════════
    print("── Analysis A: Total Variants (all transcript rows) ────────────")

    generate_and_save_tables(df, "all_total", input_dir)
    create_stat_visuals(df, "All Variants", outputs["total_log"],  use_log=True)
    create_stat_visuals(df, "All Variants", outputs["total_lin"],  use_log=False)

    # ══════════════════════════════════════════════════════════════════════════
    # ANALYSIS B: UNIQUE VARIANTS
    # (One row per distinct gene × position × protein change × group)
    # ══════════════════════════════════════════════════════════════════════════
    print("\n── Analysis B: Unique Variants (deduplicated) ──────────────────")

    # Deduplicate by (Symbol, Pos, HGVSp, Group).
    # Including Symbol ensures that two different genes overlapping at a similar
    # position (possible in dense loci) are treated as separate variants.
    # Including HGVSp ensures that two different amino acid changes at the same
    # position (e.g. p.Arg175His vs p.Arg175Cys) are treated as distinct events.
    df_unique = df.drop_duplicates(subset=[sym_c, pos_c, hgv_c, "Group"])
    print(f"  Unique variants after deduplication: {len(df_unique)} rows "
          f"(from {len(df)} total)")

    generate_and_save_tables(df_unique, "all_unique", input_dir)
    create_stat_visuals(df_unique, "Unique Variants", outputs["unique_log"], use_log=True)
    create_stat_visuals(df_unique, "Unique Variants", outputs["unique_lin"], use_log=False)

    # ── Completion summary ────────────────────────────────────────────────────
    print(f"\n{'=' * 65}")
    print("Script 10 complete. Output files:")
    for key, fname in outputs.items():
        print(f"  {fname}")
    print(f"  table_impact_all_total.csv")
    print(f"  table_consequence_all_total.csv")
    print(f"  table_impact_all_unique.csv")
    print(f"  table_consequence_all_unique.csv")
    print(f"  All saved to: {input_dir}")
    print(f"{'=' * 65}")


# ── Script entry point ─────────────────────────────────────────────────────────
# Ensures main() is only called when this script is run directly, not when
# imported as a module by another script.
if __name__ == "__main__":
    main()
