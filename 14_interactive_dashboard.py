#!/usr/bin/env python3
"""
================================================================================
Script 14 — Interactive Visualisation Dashboard
================================================================================
Pipeline    : Step 14 of 17 — runs AFTER permutation testing (13) and BEFORE
              haplotype analysis (15)

--------------------------------------------------------------------------------
PURPOSE
--------------------------------------------------------------------------------
This script synthesises the results of the entire analysis pipeline into a
single self-contained, interactive HTML dashboard. Unlike the static PNG
figures produced by earlier scripts, every plot in this dashboard is
interactive: the user can zoom, pan, hover over individual data points to
see their annotations, click legend entries to isolate specific groups, and
rotate the 3D plot.

The dashboard is generated using Plotly — a Python library that converts data
into interactive web-based figures and embeds them directly into an HTML file.
The output is a single .html file that can be opened in any modern web browser
without requiring Python, Plotly, or any other software to be installed.

--------------------------------------------------------------------------------
PLOTLY vs MATPLOTLIB / SEABORN
--------------------------------------------------------------------------------
All earlier scripts in this pipeline use Matplotlib or Seaborn to produce
static PNG figures suitable for thesis chapters and publications. This script
uses Plotly instead, which is better suited for exploratory data inspection:

  Matplotlib/Seaborn:
    - Static PNG output
    - Publication-quality formatting
    - Best for thesis figures, papers, presentations

  Plotly (this script):
    - Interactive HTML output
    - Hover tooltips showing variant-level details (rsID, HGVSp, consequence)
    - Zoom, pan, filter by legend, rotate 3D
    - Best for exploratory data inspection and sharing with collaborators

The two approaches are complementary. The static figures provide polished
representations of key results; the dashboard provides an interface for
detailed, ad-hoc exploration of the full dataset.

--------------------------------------------------------------------------------
HOW THE HTML OUTPUT WORKS
--------------------------------------------------------------------------------
Plotly converts each figure into a self-contained block of JavaScript code.
This script assembles all figures together with custom CSS styling and embeds
them into a single HTML document. The Plotly JavaScript library is loaded once
from a CDN (Content Delivery Network) for the first plot
(include_plotlyjs='cdn') and reused for all subsequent plots
(include_plotlyjs=False), keeping the file size manageable.

The resulting HTML file is entirely self-contained: no internet connection is
needed after the initial Plotly CDN load (or if opened on a machine with the
.js file cached). The file can be shared with collaborators by email or USB
and opened directly in a browser.

--------------------------------------------------------------------------------
DASHBOARD STRUCTURE — 10 INTERACTIVE PLOTS
--------------------------------------------------------------------------------
Plots 1–7 are generated from the main annotated variant table.
Plots 8–10 are generated from the SNP summary table produced by script 11
and are only included if that file exists.

  Plot 1 — Genomic Landscape Explorer
    Scatter: X = chromosomal position, Y = gnomAD NFE allele frequency.
    Colour = gene, shape = impact. Faceted by cohort (row) and tissue (column).
    Hover shows HGVSp (protein change), consequence type, and sample ID.
    This is the interactive version of the static landscape from scripts 08–09.

  Plot 2 — Variant Hierarchy Sunburst
    Hierarchical ring chart: outer ring = gene → middle ring = impact
    → inner ring = specific consequence type. Segment size = variant count.
    Click any ring segment to zoom into that subtree.
    Useful for understanding the nested composition of variant types per gene.

  Plot 3 — Variant Burden Distribution
    Box plots with individual dots (strip chart overlay): one box per
    cohort × tissue group, Y = number of variants per sample.
    Reveals whether individual samples are outliers in their group.
    Hover shows the Sample ID for each dot.

  Plot 4 — Impact vs Consequence Heatmap
    Grid heatmap: rows = VEP impact categories (HIGH → LOW), columns =
    specific consequence types. Cell colour and number = variant count.
    Shows which consequence types drive each impact level.

  Plot 5 — Gene Mutation Frequency by Cohort
    Grouped bar chart: top 10 genes by variant count, bars split by cohort.
    Hover shows exact sample counts per gene per cohort.

  Plot 6 — 3D Variant Explorer
    Three-dimensional scatter: X = chromosomal position, Y = gnomAD NFE AF,
    Z = sequencing read depth (DP). Colour = impact level.
    Requires a DP column in the input data; renders a placeholder if absent.
    Drag to rotate; scroll to zoom. Shows technical quality (depth) vs.
    population frequency vs. genomic position simultaneously.

  Plot 7 — Interactive Data Table (top 100 variants)
    Scrollable, formatted table of the first 100 rows showing key columns.
    Allows direct reading of individual variant details.

  Plot 8 — SNP Frequency Benchmarking (script 11 data)
    Interactive version of the identity plots from script 11.
    X = gnomAD NFE population frequency, Y = observed study frequency (%).
    Diagonal = line of identity (perfect agreement with population).
    Points above diagonal = enriched in our cohort vs. general population.
    Colour = tissue, shape = cohort. Hover shows rsID, gene, consequence.

  Plot 9 — SNP Carrier Frequency Heatmap (top 40 SNPs, script 11 data)
    Heatmap: rows = top 40 SNPs by mean frequency, columns = cohort × tissue
    groups. Cell colour and label = carrier frequency (%).
    Immediately shows which SNPs are pan-cohort vs. group-specific.

  Plot 10 — SNP Consequence Breakdown (script 11 data)
    Stacked horizontal bar chart: consequence types on Y-axis, count on X-axis,
    bars split by gene. Click legend entries to isolate individual genes.

--------------------------------------------------------------------------------
SNP PANELS: CONDITIONAL INCLUSION
--------------------------------------------------------------------------------
Plots 8–10 depend on the file "11_Master_Unique_SNP_Summary.xlsx" produced
by script 11. If this file does not exist when script 14 is run (e.g. if
script 11 has not been executed yet), the SNP panels are silently omitted
and the dashboard is generated with 7 plots instead of 10. A warning is
printed to the console indicating that script 11 should be run first.

--------------------------------------------------------------------------------
3D PLOT: READ DEPTH (DP) COLUMN
--------------------------------------------------------------------------------
The Z-axis of plot 6 uses the sequencing read depth (DP) at each variant
site. DP is a quality metric from the VCF FORMAT field representing how many
sequencing reads covered that genomic position. Higher depth = more confident
variant call. If no DP column is found in the input data, a placeholder empty
figure is generated with an explanatory title.

To prevent extreme outliers (very high-coverage positions) from compressing
the rest of the Z-axis, DP values are capped at the 99th percentile before
plotting. To maintain rendering performance in the browser, a random sample
of up to 500 points is used for the 3D plot.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
  GSDMB_Annotated_Report_Fixed.xlsx
    └── Sheet: "Biological_Annotations"
        Required columns: SYMBOL, Pos, IMPACT, CONSEQUENCE, COHORT, TISSUE,
                          HGVSp, Sample, gnomADe_NFE_AF
        Optional column : DP (read depth — for 3D plot)

  11_Master_Unique_SNP_Summary.xlsx  (optional — produced by script 11)
    Required if SNP panels (plots 8–10) are desired.

--------------------------------------------------------------------------------
OUTPUT
--------------------------------------------------------------------------------
  14_Interactive_Dashboard.html
    Single self-contained HTML file. Open in any modern web browser
    (Chrome, Firefox, Edge, Safari). No additional software required.
    File size is typically 5–20 MB depending on dataset size.

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python 14_interactive_dashboard.py

  Then open 14_Interactive_Dashboard.html in a web browser.
  (Paths are configured in the CONFIGURATION block below.)

--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------
  Python ≥ 3.8
  pandas, numpy, plotly, openpyxl

  Install: pip install pandas numpy plotly openpyxl
================================================================================
"""

# ──────────────────────────────────────────────────────────────────────────────
# STANDARD LIBRARY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import os       # File path construction and existence checks
import warnings # Suppress non-critical library warnings

# ──────────────────────────────────────────────────────────────────────────────
# THIRD-PARTY IMPORTS
# ──────────────────────────────────────────────────────────────────────────────

import numpy as np                          # Numerical operations
import pandas as pd                         # Data loading and manipulation
import plotly.express as px                 # High-level Plotly figure constructors
import plotly.graph_objects as go           # Low-level Plotly figure objects
from plotly.subplots import make_subplots   # Multi-panel figure layouts

warnings.filterwarnings("ignore")


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

BASE_PATH   = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
INPUT_FILE  = os.path.join(BASE_PATH, "GSDMB_Annotated_Report_Fixed.xlsx")
OUTPUT_HTML = os.path.join(BASE_PATH, "14_Interactive_Dashboard.html")


# ══════════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def find_col(df: pd.DataFrame, target: str):
    """
    Case-insensitive column name lookup, also stripping leading/trailing spaces.

    Parameters
    ----------
    df     : pd.DataFrame — the DataFrame to search
    target : str          — the desired column name (any case)

    Returns
    -------
    str or None
        The actual column name as it appears in the DataFrame, or None.
    """
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — GENOMIC LANDSCAPE EXPLORER
# ══════════════════════════════════════════════════════════════════════════════

def create_interactive_scatter(df, sym_c, pos_c, imp_c, coh_c, tis_c):
    """
    Interactive scatter plot of all variants along chromosome 17.

    Axes and encodings
    ------------------
    X-axis   : Genomic position on chromosome 17
    Y-axis   : gnomAD NFE population allele frequency (log scale)
               Log scale is used because frequencies span several orders of
               magnitude — rare variants (AF ~ 0.001) and common SNPs
               (AF ~ 0.5) would not be distinguishable on a linear scale.
    Colour   : Gene (GSDMB, ORMDL3, ERBB2, etc.)
    Shape    : VEP impact category (HIGH / MODERATE / MODIFIER / LOW)
    Facets   : Row = cancer cohort; Column = tissue type

    Hover tooltip shows the protein-level change (HGVSp), specific consequence
    type (e.g. missense_variant), and the sample ID carrying that variant.

    This is the interactive counterpart to the static landscape figures from
    scripts 08 (all genes) and 09 (GSDMB only).

    Parameters
    ----------
    df    : pd.DataFrame — full annotated variant table
    sym_c : str          — gene symbol column name
    pos_c : str          — genomic position column name
    imp_c : str          — VEP impact column name
    coh_c : str          — cohort column name
    tis_c : str          — tissue type column name

    Returns
    -------
    plotly.graph_objects.Figure
    """
    fig = px.scatter(
        df,
        x           = pos_c,
        y           = "gnomADe_NFE_AF",
        color       = sym_c,
        symbol      = imp_c,
        facet_row   = coh_c,
        facet_col   = tis_c,
        hover_data  = ["HGVSp", "Consequence", "Sample"],
        title       = "Interactive Genomic Landscape",
        labels      = {
            pos_c:            "Genomic Position (Chr17)",
            "gnomADe_NFE_AF": "Population Frequency (gnomAD NFE)",
            sym_c:            "Gene",
            imp_c:            "Impact",
        },
        height                  = 800,
        color_discrete_sequence = px.colors.qualitative.Set3,
    )

    fig.update_layout(
        hovermode      = "closest",
        font           = dict(size=11),
        title_font_size= 16,
    )

    fig.update_xaxes(title_text="Position on Chr17 (bp)")
    # Log scale on Y-axis so rare and common variants are both visible
    fig.update_yaxes(type="log", title_text="Allele Frequency (log scale)")

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — VARIANT HIERARCHY SUNBURST
# ══════════════════════════════════════════════════════════════════════════════

def create_variant_sunburst(df, sym_c, imp_c, con_c):
    """
    Hierarchical sunburst chart showing the composition of variants by gene,
    impact, and consequence type.

    What is a sunburst chart?
    -------------------------
    A sunburst is a multi-level pie chart where each concentric ring represents
    one level of a hierarchy. Here the three levels are:
      Outer ring : Gene (e.g. GSDMB, ERBB2)
      Middle ring : Impact category within that gene (HIGH, MODERATE, etc.)
      Inner ring  : Specific consequence type within that impact level
                    (e.g. missense_variant, intron_variant)

    Segment area is proportional to variant count. Clicking a segment zooms
    into that subtree. This makes it easy to see, for example, what proportion
    of GSDMB variants are HIGH impact, and what specific consequence types
    drive that.

    Parameters
    ----------
    df    : pd.DataFrame — full annotated variant table
    sym_c : str          — gene symbol column name
    imp_c : str          — VEP impact column name
    con_c : str          — VEP consequence column name

    Returns
    -------
    plotly.graph_objects.Figure
    """
    # Aggregate variant counts for each (gene, impact, consequence) combination
    sunburst_data = (
        df.groupby([sym_c, imp_c, con_c])
        .size()
        .reset_index(name="Count")
    )

    fig = px.sunburst(
        sunburst_data,
        path                  = [sym_c, imp_c, con_c],  # Hierarchy order
        values                = "Count",
        title                 = "Variant Hierarchy: Gene → Impact → Consequence",
        color                 = "Count",
        color_continuous_scale= "RdYlBu_r",   # Red = high count, blue = low count
        height                = 700,
    )

    # Show both the label and the percentage of the parent segment in each slice
    fig.update_traces(textinfo="label+percent parent")
    fig.update_layout(font=dict(size=12))

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 3 — VARIANT BURDEN DISTRIBUTION (BOX PLOTS)
# ══════════════════════════════════════════════════════════════════════════════

def create_cohort_comparison_box(df, coh_c, tis_c):
    """
    Box plots showing the distribution of per-sample variant counts across
    cohorts and tissue types.

    What does "variant burden" mean?
    ---------------------------------
    Each sample in the dataset was sequenced and had variants called. The
    total number of variants detected in a sample is its "variant burden".
    Differences in burden between tumour and healthy samples, or between
    cancer types, may reflect biological differences in mutation load or
    technical differences in sequencing depth.

    Individual sample dots are overlaid on the boxes (points='all') so that
    outlier samples — those with unusually high or low variant counts — can
    be identified. Hovering over a dot shows the Sample ID.

    Parameters
    ----------
    df    : pd.DataFrame — full annotated variant table (one row per variant)
    coh_c : str          — cohort column name
    tis_c : str          — tissue type column name

    Returns
    -------
    plotly.graph_objects.Figure
    """
    # Count how many variant rows belong to each (sample, cohort, tissue) group.
    # This gives the per-sample variant count (variant burden).
    sample_counts = (
        df.groupby(["Sample", coh_c, tis_c])
        .size()
        .reset_index(name="Variant_Count")
    )

    fig = px.box(
        sample_counts,
        x              = coh_c,
        y              = "Variant_Count",
        color          = tis_c,
        title          = "Variant Burden Distribution by Cohort and Tissue Type",
        labels         = {
            coh_c:           "Cancer Type",
            "Variant_Count": "Number of Variants per Sample",
            tis_c:           "Tissue Type",
        },
        points             = "all",    # Show individual sample dots overlaid on box
        hover_data         = ["Sample"],
        height             = 500,
        color_discrete_map = {"Tumour": "#e74c3c", "Healthy": "#3498db"},
    )

    fig.update_layout(
        boxmode = "group",   # Place tumour and healthy boxes side by side per cohort
        font    = dict(size=12),
    )

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 4 — IMPACT vs CONSEQUENCE HEATMAP
# ══════════════════════════════════════════════════════════════════════════════

def create_impact_consequence_heatmap(df, imp_c, con_c):
    """
    Interactive heatmap cross-tabulating VEP impact categories against specific
    consequence types.

    Reading the heatmap
    -------------------
    Each cell shows how many variants have a given impact level AND a given
    consequence type. For example, the cell at (HIGH, splice_donor_variant)
    shows how many variants are both HIGH impact and a splice donor disruption.

    Rows are ordered from most to least severe impact (HIGH → LOW). Columns
    are the specific VEP consequence strings (e.g. missense_variant,
    intron_variant, 3_prime_UTR_variant). The colour scale (yellow → red)
    reflects variant count: darker = more variants.

    Hovering over any cell shows the exact count.

    Parameters
    ----------
    df    : pd.DataFrame — full annotated variant table
    imp_c : str          — VEP impact column name
    con_c : str          — VEP consequence column name

    Returns
    -------
    plotly.graph_objects.Figure
    """
    # Aggregate and pivot to create the impact × consequence count matrix
    heatmap_data = (
        df.groupby([imp_c, con_c])
        .size()
        .reset_index(name="Count")
    )
    pivot_data = (
        heatmap_data
        .pivot(index=imp_c, columns=con_c, values="Count")
        .fillna(0)
    )

    # Reorder rows from most to least severe impact
    impact_order = ["HIGH", "MODERATE", "MODIFIER", "LOW"]
    pivot_data = pivot_data.reindex(
        [i for i in impact_order if i in pivot_data.index]
    )

    fig = go.Figure(data=go.Heatmap(
        z             = pivot_data.values,
        x             = pivot_data.columns,
        y             = pivot_data.index,
        colorscale    = "YlOrRd",          # Yellow (low) → red (high)
        text          = pivot_data.values,
        texttemplate  = "%{text}",         # Display count value in each cell
        textfont      = {"size": 10},
        hoverongaps   = False,
        colorbar      = dict(title="Variant Count"),
    ))

    fig.update_layout(
        title       = "Impact vs Consequence Heatmap",
        xaxis_title = "Functional Consequence Type (VEP)",
        yaxis_title = "Variant Impact Level (VEP)",
        height      = 500,
        font        = dict(size=11),
    )

    # Rotate X-axis labels so long consequence names do not overlap
    fig.update_xaxes(tickangle=-45)

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 5 — GENE MUTATION FREQUENCY BY COHORT
# ══════════════════════════════════════════════════════════════════════════════

def create_gene_coverage_timeline(df, sym_c, sam_c, coh_c):
    """
    Grouped bar chart showing how many unique samples carry at least one
    variant in each of the top 10 most frequently mutated genes, split by
    cancer cohort.

    This plot answers the question: "Which genes are mutated in the most
    patients, and does the pattern differ between breast and endometrial cancer?"

    Parameters
    ----------
    df    : pd.DataFrame — full annotated variant table
    sym_c : str          — gene symbol column name
    sam_c : str          — sample identifier column name
    coh_c : str          — cohort column name

    Returns
    -------
    plotly.graph_objects.Figure
    """
    # Identify the top 10 genes by total variant row count across all samples
    top_genes = df[sym_c].value_counts().head(10).index

    # For each top gene × cohort combination, count how many unique samples
    # have at least one variant in that gene
    timeline_data = []
    for gene in top_genes:
        gene_df = df[df[sym_c] == gene]
        for cohort in df[coh_c].unique():
            n_samples = gene_df[gene_df[coh_c] == cohort][sam_c].nunique()
            timeline_data.append({
                "Gene":         gene,
                "Cohort":       cohort,
                "Sample_Count": n_samples,
            })

    timeline_df = pd.DataFrame(timeline_data)

    fig = px.bar(
        timeline_df,
        x                       = "Gene",
        y                       = "Sample_Count",
        color                   = "Cohort",
        title                   = "Top 10 Genes: Number of Affected Samples by Cohort",
        labels                  = {
            "Sample_Count": "Number of Samples with Variant",
            "Gene":         "Gene",
        },
        barmode                 = "group",  # Side-by-side bars per cohort
        height                  = 500,
        color_discrete_sequence = px.colors.qualitative.Bold,
    )

    fig.update_layout(font=dict(size=12))

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 6 — 3D VARIANT EXPLORER
# ══════════════════════════════════════════════════════════════════════════════

def create_3d_scatter(df, pos_c, sym_c, imp_c):
    """
    Three-dimensional interactive scatter plot for multi-dimensional variant
    exploration.

    Axes
    ----
    X — Genomic position on chromosome 17 (bp)
    Y — gnomAD NFE population allele frequency
        Shows whether each variant is rare or common in the general population.
    Z — Sequencing read depth (DP) at the variant site
        DP (Depth of coverage) is the number of sequencing reads that covered
        the genomic position where the variant was called. Higher depth
        indicates a more reliable, better-supported variant call.
        A variant with DP = 5 is much less certain than one with DP = 200.

    Colour — VEP impact level (HIGH = 0, MODERATE = 1, MODIFIER = 2, LOW = 3)
             Mapped to a red (HIGH) → green (LOW) colour scale.

    Interaction
    -----------
    Click and drag to rotate the plot in 3D. Scroll to zoom. Hover over any
    point to see the gene name, position, population frequency, and read depth.

    Performance and data handling
    ------------------------------
    Rendering thousands of 3D points can be slow in a web browser. A random
    sample of up to 500 points is used to keep the dashboard responsive.
    DP values above the 99th percentile are capped to prevent outliers from
    compressing the visible Z-axis range for the majority of points.

    Fallback behaviour
    ------------------
    If no DP column is found in the input data, the function returns an empty
    placeholder figure with an explanatory title rather than crashing.

    Parameters
    ----------
    df    : pd.DataFrame — full annotated variant table
    pos_c : str          — genomic position column name
    sym_c : str          — gene symbol column name
    imp_c : str          — VEP impact column name

    Returns
    -------
    plotly.graph_objects.Figure
    """
    # Look for a DP column (case-insensitive); return placeholder if absent
    dp_col = next((c for c in df.columns if c.upper() == "DP"), None)
    if dp_col is None:
        print("  WARNING: No DP (read depth) column found — 3D scatter will be a placeholder.")
        fig = go.Figure()
        fig.update_layout(
            title  = "3D Variant Explorer (unavailable — DP column not found in input data)",
            height = 700,
        )
        return fig

    # Prepare a clean copy with only the columns needed for this plot
    df_3d = df[[pos_c, "gnomADe_NFE_AF", dp_col, imp_c, sym_c]].copy()
    df_3d["Read_Depth"]    = pd.to_numeric(df_3d[dp_col],          errors="coerce")
    df_3d["gnomADe_NFE_AF"] = pd.to_numeric(df_3d["gnomADe_NFE_AF"], errors="coerce")

    # Drop rows where any axis value is missing
    df_3d = df_3d.dropna(subset=["Read_Depth", "gnomADe_NFE_AF"])

    # Cap DP at the 99th percentile to prevent extreme outliers from
    # collapsing the visible Z-axis range for the majority of variants
    depth_cap = df_3d["Read_Depth"].quantile(0.99)
    df_3d["Read_Depth"] = df_3d["Read_Depth"].clip(upper=depth_cap)

    # Sample up to 500 points for browser rendering performance
    df_sample = df_3d.sample(min(500, len(df_3d)), random_state=42)

    fig = go.Figure(data=[go.Scatter3d(
        x    = df_sample[pos_c],
        y    = df_sample["gnomADe_NFE_AF"],
        z    = df_sample["Read_Depth"],
        mode = "markers",
        marker=dict(
            size  = 5,
            # Map impact categories to integers for the colour scale:
            # 0 = HIGH (red) → 3 = LOW (green), reflecting severity order
            color     = df_sample[imp_c].map(
                {"HIGH": 0, "MODERATE": 1, "MODIFIER": 2, "LOW": 3}
            ),
            colorscale= "RdYlGn_r",   # Red (HIGH) → green (LOW)
            showscale = True,
            colorbar  = dict(
                title    = "Impact",
                tickvals = [0, 1, 2, 3],
                ticktext = ["HIGH", "MODERATE", "MODIFIER", "LOW"],
            ),
            line = dict(color="black", width=0.5),
        ),
        text         = df_sample[sym_c],
        hovertemplate= (
            "<b>%{text}</b><br>"
            "Position: %{x:,}<br>"
            "gnomAD NFE AF: %{y:.4f}<br>"
            "Read Depth: %{z:.0f}×"
            "<extra></extra>"
        ),
    )])

    fig.update_layout(
        title = "3D Variant Explorer: Position × Population Frequency × Read Depth",
        scene = dict(
            xaxis_title = "Genomic Position (Chr17)",
            yaxis_title = "gnomAD NFE Allele Frequency",
            zaxis_title = "Read Depth — DP (×)",
        ),
        height = 700,
        font   = dict(size=11),
    )

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 7 — INTERACTIVE DATA TABLE
# ══════════════════════════════════════════════════════════════════════════════

def create_interactive_table(df, sym_c, pos_c, imp_c, con_c):
    """
    Scrollable, formatted table showing key columns for the first 100 variants.

    The table is limited to 100 rows to keep the dashboard file size and
    browser rendering time manageable. It provides a direct way to read
    individual variant details — gene, position, protein change, impact,
    consequence, gnomAD frequency, sample, cohort, and tissue — without
    requiring the user to open the Excel source file.

    Parameters
    ----------
    df    : pd.DataFrame — full annotated variant table
    sym_c : str          — gene symbol column name
    pos_c : str          — genomic position column name
    imp_c : str          — VEP impact column name
    con_c : str          — VEP consequence column name

    Returns
    -------
    plotly.graph_objects.Figure
    """
    # Select the most informative columns for display
    table_df = df[[
        sym_c, pos_c, "HGVSp", imp_c, con_c,
        "gnomADe_NFE_AF", "Sample", "Cohort", "Tissue",
    ]].copy()

    # Round the gnomAD frequency to 4 decimal places for readable display
    table_df["gnomADe_NFE_AF"] = table_df["gnomADe_NFE_AF"].round(4)

    # Limit to first 100 rows for performance
    table_df = table_df.head(100)

    fig = go.Figure(data=[go.Table(
        header=dict(
            values     = list(table_df.columns),
            fill_color = "#3498db",
            font       = dict(color="white", size=12),
            align      = "left",
        ),
        cells=dict(
            values     = [table_df[col] for col in table_df.columns],
            # Alternate row shading (light grey / white) for readability
            fill_color = [["#ecf0f1", "white"] * (len(table_df) // 2 + 1)],
            align      = "left",
            font       = dict(size=11),
        ),
    )])

    fig.update_layout(
        title  = "Variant Data Table (first 100 rows — scroll to explore)",
        height = 600,
    )

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# SNP DATA LOADER
# ══════════════════════════════════════════════════════════════════════════════

def load_snp_data():
    """
    Load the SNP summary table produced by script 11.

    The SNP panels (plots 8–10) depend on this file. If it does not exist —
    for example because script 11 has not yet been run — this function returns
    None gracefully so the dashboard can still be generated without the SNP
    section rather than crashing.

    Returns
    -------
    pd.DataFrame or None
        The SNP summary table, or None if the file is not found.
    """
    snp_file = os.path.join(BASE_PATH, "11_Master_Unique_SNP_Summary.xlsx")

    if not os.path.exists(snp_file):
        print(f"  WARNING: SNP summary not found at {snp_file}")
        print("  Run script 11 first to enable SNP panels (plots 8–10).")
        return None

    df_snp = pd.read_excel(snp_file)
    df_snp.columns = df_snp.columns.str.strip()  # Remove any accidental whitespace
    print(f"  SNP data loaded: {len(df_snp)} unique SNPs")
    return df_snp


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 8 — SNP FREQUENCY BENCHMARKING
# ══════════════════════════════════════════════════════════════════════════════

def create_snp_benchmarking(df_snp):
    """
    Interactive scatter plot benchmarking observed SNP carrier frequencies in
    this study against gnomAD NFE population reference frequencies.

    This is the interactive version of the static identity plots from script 11.

    The data from script 11 is stored in wide format (one row per SNP, one
    column per cohort × tissue frequency). This function first reshapes
    (melts) it to long format so each row represents one SNP in one group,
    enabling colour and symbol encoding by tissue and cohort.

    Reading the plot
    ----------------
    X-axis : gnomAD NFE allele frequency (population reference)
    Y-axis : Observed carrier frequency in this study (%)
    Diagonal reference line : y = x × 100 — perfect agreement with population

    Points above the diagonal → SNP more common in our cohort than in the
      general European population (possible enrichment)
    Points below the diagonal → SNP less common (possible depletion)
    Points on the diagonal   → frequency matches population expectation

    Parameters
    ----------
    df_snp : pd.DataFrame — SNP summary table from script 11

    Returns
    -------
    plotly.graph_objects.Figure
    """
    # Detect cohort × tissue frequency columns dynamically
    # (e.g. "Breast_Tumour_Frequency_%", "Endometrium_Healthy_Frequency_%")
    freq_cols = [c for c in df_snp.columns if c.endswith("Frequency_%")]

    if not freq_cols or "gnomAD_NFE_AF" not in df_snp.columns:
        fig = go.Figure()
        fig.update_layout(
            title  = "SNP Benchmarking (required columns not found in SNP summary file)",
            height = 600,
        )
        return fig

    # Reshape from wide to long format:
    # Each row becomes one SNP × one cohort-tissue group
    id_vars = [c for c in ["Variant_ID", "SYMBOL", "Consequence", "IMPACT", "gnomAD_NFE_AF"]
               if c in df_snp.columns]
    df_long = df_snp.melt(
        id_vars   = id_vars,
        value_vars = freq_cols,
        var_name  = "Group",
        value_name= "Study_Frequency_%",
    )

    # Extract cohort and tissue from the column name
    # e.g. "Breast_Tumour_Frequency_%" → Cohort="Breast", Tissue="Tumour"
    df_long["Cohort"] = df_long["Group"].str.split("_").str[0]
    df_long["Tissue"] = df_long["Group"].str.split("_").str[1]

    # Exclude rows where the study frequency is zero (SNP absent from that group)
    df_long = df_long[df_long["Study_Frequency_%"] > 0].copy()

    fig = px.scatter(
        df_long,
        x          = "gnomAD_NFE_AF",
        y          = "Study_Frequency_%",
        color      = "Tissue",
        symbol     = "Cohort",
        hover_data = {
            "Variant_ID":        True,
            "SYMBOL":            True,
            "Consequence":       True,
            "IMPACT":            True,
            "gnomAD_NFE_AF":     ":.4f",
            "Study_Frequency_%": ":.1f",
            "Group":             False,   # Hide the raw column name from the tooltip
        },
        color_discrete_map = {"Tumour": "#e74c3c", "Healthy": "#3498db"},
        title      = "SNP Frequency Benchmarking: Study Cohort vs. gnomAD NFE Population",
        labels     = {
            "gnomAD_NFE_AF":     "gnomAD NFE Allele Frequency (population reference)",
            "Study_Frequency_%": "Carrier Frequency in Study (%)",
        },
        height     = 600,
    )

    # Diagonal reference line: y = x × 100 (population frequency scaled to %)
    x_max = df_long["gnomAD_NFE_AF"].max() * 1.05
    fig.add_trace(go.Scatter(
        x         = [0, x_max],
        y         = [0, x_max * 100],
        mode      = "lines",
        line      = dict(color="grey", dash="dash", width=1),
        name      = "Population baseline (y = x)",
        hoverinfo = "skip",   # Do not show tooltip on the reference line
    ))

    fig.update_layout(
        hovermode    = "closest",
        legend_title = "Tissue / Cohort",
        font         = dict(size=11),
    )

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 9 — SNP CARRIER FREQUENCY HEATMAP
# ══════════════════════════════════════════════════════════════════════════════

def create_snp_heatmap(df_snp):
    """
    Heatmap showing carrier frequencies for the top 40 SNPs across all
    cohort × tissue groups.

    Reading the heatmap
    -------------------
    Rows    : The 40 SNPs with the highest mean carrier frequency across groups
              (labelled as "rsID  (GENE)" or "GENE:HGVSp  (GENE)")
    Columns : The four cohort × tissue groups
              (e.g. "Breast Tumour", "Breast Healthy", "Endometrium Tumour", etc.)
    Cells   : Carrier frequency (%) — colour and numeric label

    Pan-cohort SNPs (common across all four columns) appear uniformly coloured.
    Group-specific SNPs appear as a single bright cell against pale neighbours.
    Hovering over a cell shows the full SNP details including gnomAD frequency
    and consequence type.

    Parameters
    ----------
    df_snp : pd.DataFrame — SNP summary table from script 11

    Returns
    -------
    plotly.graph_objects.Figure
    """
    freq_cols = [c for c in df_snp.columns if c.endswith("Frequency_%")]

    if not freq_cols:
        fig = go.Figure()
        fig.update_layout(
            title  = "SNP Heatmap (frequency columns not found)",
            height = 600,
        )
        return fig

    # Select the 40 SNPs with the highest mean carrier frequency across all groups
    df_snp["Mean_Freq"] = df_snp[freq_cols].mean(axis=1)
    top_snps = df_snp.nlargest(40, "Mean_Freq").copy()

    # Build display labels: "rs12345  (GSDMB)"
    label_col = "Variant_ID" if "Variant_ID" in top_snps.columns else top_snps.index
    top_snps["Label"] = (
        top_snps[label_col].astype(str) + "  (" + top_snps["SYMBOL"].astype(str) + ")"
        if "SYMBOL" in top_snps.columns
        else top_snps[label_col].astype(str)
    )

    z_data  = top_snps[freq_cols].values
    # Tidy column labels: "Breast_Tumour_Frequency_%" → "Breast Tumour"
    x_labels = [c.replace("_Frequency_%", "").replace("_", " ") for c in freq_cols]
    y_labels = top_snps["Label"].tolist()

    # Build rich hover text matrix (one string per cell)
    hover_text = []
    for _, row in top_snps.iterrows():
        row_text = []
        for col, x_lbl in zip(freq_cols, x_labels):
            row_text.append(
                f"SNP: {row.get('Variant_ID', '')}<br>"
                f"Gene: {row.get('SYMBOL', '')}<br>"
                f"Group: {x_lbl}<br>"
                f"Carrier frequency: {row[col]:.1f}%<br>"
                f"gnomAD NFE AF: {row.get('gnomAD_NFE_AF', 'N/A')}<br>"
                f"Consequence: {row.get('Consequence', 'N/A')}"
            )
        hover_text.append(row_text)

    fig = go.Figure(data=go.Heatmap(
        z             = z_data,
        x             = x_labels,
        y             = y_labels,
        colorscale    = "YlOrRd",   # Yellow (low frequency) → red (high frequency)
        text          = [[f"{v:.1f}%" for v in row] for row in z_data],
        texttemplate  = "%{text}",
        textfont      = dict(size=9),
        hovertext     = hover_text,
        hovertemplate = "%{hovertext}<extra></extra>",
        colorbar      = dict(title="Carrier<br>Frequency (%)"),
        zmin          = 0,
    ))

    fig.update_layout(
        title       = "Top 40 SNPs: Carrier Frequency Across Cohort × Tissue Groups",
        xaxis_title = "Cohort — Tissue Group",
        yaxis_title = "SNP (rsID / Gene)",
        height      = 900,
        font        = dict(size=11),
        yaxis       = dict(tickfont=dict(size=9)),
    )

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 10 — SNP CONSEQUENCE BREAKDOWN
# ══════════════════════════════════════════════════════════════════════════════

def create_snp_consequence_breakdown(df_snp):
    """
    Stacked horizontal bar chart showing unique SNP counts by consequence type,
    with bars split by gene.

    Reading the chart
    -----------------
    Each bar represents one VEP consequence type (e.g. intron_variant,
    missense_variant, synonymous_variant). Bar length = total number of unique
    SNPs with that consequence. Bar segments are coloured by gene, showing
    which genes contribute SNPs to each consequence category.

    Clicking a gene name in the legend toggles that gene's contribution
    on/off across all bars — useful for isolating GSDMB-specific SNPs vs.
    contributions from neighbouring genes.

    Parameters
    ----------
    df_snp : pd.DataFrame — SNP summary table from script 11

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if "Consequence" not in df_snp.columns or "SYMBOL" not in df_snp.columns:
        fig = go.Figure()
        fig.update_layout(
            title  = "SNP Consequence Breakdown (required columns not found)",
            height = 500,
        )
        return fig

    # Count unique SNPs per consequence × gene combination
    con_counts = (
        df_snp.groupby(["Consequence", "SYMBOL"])
        .size()
        .reset_index(name="Count")
    )

    # Order consequence types by total count (ascending) so the longest bars
    # appear at the top of the horizontal chart
    con_order = (
        con_counts.groupby("Consequence")["Count"]
        .sum()
        .sort_values(ascending=True)
        .index.tolist()
    )

    fig = px.bar(
        con_counts,
        x              = "Count",
        y              = "Consequence",
        color          = "SYMBOL",
        orientation    = "h",   # Horizontal bars — better for long consequence names
        title          = "SNP Distribution by Functional Consequence and Gene",
        labels         = {
            "Count":       "Number of Unique SNPs",
            "Consequence": "Functional Consequence Type (VEP)",
            "SYMBOL":      "Gene",
        },
        category_orders         = {"Consequence": con_order},
        height                  = 600,
        color_discrete_sequence = px.colors.qualitative.Bold,
    )

    fig.update_layout(
        barmode      = "stack",      # Stacked bars: gene contributions are additive
        legend_title = "Gene",
        font         = dict(size=11),
        hovermode    = "closest",
    )

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# MAIN ORCHESTRATION FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def create_dashboard():
    """
    Orchestrate the full dashboard creation pipeline.

    Execution order
    ---------------
    1.  Load the annotated variant table.
    2.  Standardise tissue and cohort labels.
    3.  Attempt to load the SNP summary from script 11 (optional).
    4.  Generate all 7 core interactive figures (plots 1–7).
    5.  If SNP data is available, generate plots 8–10.
    6.  Assemble all figures plus CSS-styled HTML scaffolding into a single
        HTML string.
    7.  Write the HTML string to the output file.
    """
    print("=" * 70)
    print("SCRIPT 14 — Interactive Visualisation Dashboard")
    print("=" * 70)
    print(f"\n  Input : {INPUT_FILE}")
    print(f"  Output: {OUTPUT_HTML}")
    print(f"  Technology: Plotly (self-contained interactive HTML)\n")
    print("=" * 70 + "\n")

    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: Input file not found: {INPUT_FILE}")
        print("       Please update 'BASE_PATH' in the CONFIGURATION block.")
        return

    # ── Load and standardise main variant table ───────────────────────────────
    df = pd.read_excel(INPUT_FILE, sheet_name="Biological_Annotations")
    print(f"Loaded {len(df)} rows from sheet 'Biological_Annotations'")

    sym_c = find_col(df, "SYMBOL")
    pos_c = find_col(df, "Pos")
    imp_c = find_col(df, "IMPACT")
    con_c = find_col(df, "CONSEQUENCE")
    coh_c = find_col(df, "COHORT")
    tis_c = find_col(df, "TISSUE")

    df[tis_c] = (
        df[tis_c].astype(str).str.strip()
        .replace({"Tumor": "Tumour", "Normal": "Healthy", "Control": "Healthy"})
    )
    df[coh_c] = df[coh_c].astype(str).str.strip()

    print(f"  Variants: {len(df):,} | Samples: {df['Sample'].nunique()} | "
          f"Genes: {df[sym_c].nunique()}\n")

    # ── Load SNP data (optional — script 11 output) ───────────────────────────
    print("Loading SNP summary data (from script 11)...")
    df_snp = load_snp_data()

    # ── Generate core figures (plots 1–7) ─────────────────────────────────────
    print("\nGenerating interactive plots:")
    print("  [1/7] Genomic landscape scatter...")
    fig1 = create_interactive_scatter(df, sym_c, pos_c, imp_c, coh_c, tis_c)

    print("  [2/7] Variant hierarchy sunburst...")
    fig2 = create_variant_sunburst(df, sym_c, imp_c, con_c)

    print("  [3/7] Variant burden box plots...")
    fig3 = create_cohort_comparison_box(df, coh_c, tis_c)

    print("  [4/7] Impact vs Consequence heatmap...")
    fig4 = create_impact_consequence_heatmap(df, imp_c, con_c)

    print("  [5/7] Gene mutation frequency bar chart...")
    fig5 = create_gene_coverage_timeline(df, sym_c, "Sample", coh_c)

    print("  [6/7] 3D variant explorer...")
    fig6 = create_3d_scatter(df, pos_c, sym_c, imp_c)

    print("  [7/7] Interactive data table...")
    fig7 = create_interactive_table(df, sym_c, pos_c, imp_c, con_c)

    # ── Generate SNP figures (plots 8–10, conditional) ────────────────────────
    snp_section_html   = ""
    snp_stat_card_html = ""

    if df_snp is not None:
        n_snps = len(df_snp)

        # SNP count stat card for the summary header
        snp_stat_card_html = f"""
                <div class="stat-card">
                    <div class="stat-number">{n_snps:,}</div>
                    <div class="stat-label">Established SNPs (gnomAD NFE &gt;1%)</div>
                </div>"""

        print("  [8/10] SNP frequency benchmarking...")
        fig8 = create_snp_benchmarking(df_snp)

        print("  [9/10] SNP carrier frequency heatmap...")
        fig9 = create_snp_heatmap(df_snp)

        print("  [10/10] SNP consequence breakdown...")
        fig10 = create_snp_consequence_breakdown(df_snp)

        snp_section_html = f"""
            <div class="section-divider">
                <h2 class="section-heading">&#128202; SNP Analysis (Script 11)</h2>
                <p class="section-desc">
                    Common variants defined as gnomAD NFE allele frequency &gt;1%.
                    Carrier frequencies are calculated per cohort–tissue group.
                </p>
            </div>
            <div class="plot-section">
                <h2 class="plot-title">8. SNP Frequency Benchmarking vs. gnomAD Population</h2>
                {fig8.to_html(include_plotlyjs=False, div_id='plot8')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">9. SNP Carrier Frequency Heatmap (Top 40)</h2>
                {fig9.to_html(include_plotlyjs=False, div_id='plot9')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">10. SNP Distribution by Functional Consequence</h2>
                {fig10.to_html(include_plotlyjs=False, div_id='plot10')}
            </div>"""

    # ── Assemble HTML document ────────────────────────────────────────────────
    # The first figure uses include_plotlyjs='cdn' to load the Plotly JavaScript
    # library from a Content Delivery Network. All subsequent figures use
    # include_plotlyjs=False to reuse the already-loaded library, keeping the
    # overall file size manageable.
    print("\nAssembling dashboard HTML...")

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>GSDMB Genomics Analysis Dashboard</title>
        <meta charset="utf-8">
        <style>
            body {{
                font-family: 'Arial', sans-serif;
                margin: 0;
                padding: 20px;
                background-color: #f5f5f5;
            }}
            .header {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 30px;
                text-align: center;
                border-radius: 10px;
                margin-bottom: 30px;
                box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            }}
            .header h1 {{ margin: 0; font-size: 2.5em; }}
            .header p  {{ margin: 10px 0 0 0; font-size: 1.2em; opacity: 0.9; }}
            .container {{
                max-width: 1400px;
                margin: 0 auto;
            }}
            .plot-section {{
                background: white;
                padding: 20px;
                margin-bottom: 30px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            .plot-title {{
                color: #333;
                border-bottom: 3px solid #667eea;
                padding-bottom: 10px;
                margin-bottom: 20px;
                font-size: 1.5em;
            }}
            .stats-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 20px;
                margin-bottom: 30px;
            }}
            .stat-card {{
                background: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                text-align: center;
            }}
            .stat-number {{
                font-size: 2.5em;
                font-weight: bold;
                color: #667eea;
            }}
            .stat-label {{
                color: #666;
                margin-top: 10px;
            }}
            .section-divider {{
                background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
                color: white;
                padding: 20px 30px;
                border-radius: 8px;
                margin-bottom: 20px;
                margin-top: 10px;
            }}
            .section-heading {{
                margin: 0 0 8px 0;
                font-size: 1.6em;
            }}
            .section-desc {{
                margin: 0;
                font-size: 0.95em;
                opacity: 0.92;
            }}
            .footer {{
                text-align: center;
                padding: 20px;
                color: #666;
                border-top: 1px solid #ddd;
                margin-top: 40px;
            }}
        </style>
    </head>
    <body>
        <div class="container">

            <!-- ── Header ── -->
            <div class="header">
                <h1>&#129516; GSDMB Genomics Analysis Dashboard</h1>
                <p>Interactive Variant Exploration — All plots support zoom, pan, and hover</p>
                <p style="font-size: 0.9em; margin-top: 15px;">
                    Dataset: {len(df):,} variants &nbsp;|&nbsp;
                    {df['Sample'].nunique()} samples &nbsp;|&nbsp;
                    {df[sym_c].nunique()} genes
                </p>
            </div>

            <!-- ── Summary stat cards ── -->
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-number">{len(df):,}</div>
                    <div class="stat-label">Total Variant Rows</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{df['Sample'].nunique()}</div>
                    <div class="stat-label">Unique Samples</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{df[sym_c].nunique()}</div>
                    <div class="stat-label">Genes Analysed</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{(df[imp_c] == 'HIGH').sum()}</div>
                    <div class="stat-label">HIGH Impact Variants</div>
                </div>
                {snp_stat_card_html}
            </div>

            <!-- ── Core plots (1–7) ── -->
            <div class="plot-section">
                <h2 class="plot-title">1. Genomic Landscape Explorer</h2>
                {fig1.to_html(include_plotlyjs='cdn', div_id='plot1')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">2. Variant Hierarchy (Sunburst)</h2>
                {fig2.to_html(include_plotlyjs=False, div_id='plot2')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">3. Variant Burden by Cohort and Tissue</h2>
                {fig3.to_html(include_plotlyjs=False, div_id='plot3')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">4. Impact vs Consequence Matrix</h2>
                {fig4.to_html(include_plotlyjs=False, div_id='plot4')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">5. Gene Mutation Frequency by Cohort</h2>
                {fig5.to_html(include_plotlyjs=False, div_id='plot5')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">6. 3D Variant Explorer (Position × Population Frequency × Read Depth)</h2>
                {fig6.to_html(include_plotlyjs=False, div_id='plot6')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">7. Variant Data Table (first 100 rows)</h2>
                {fig7.to_html(include_plotlyjs=False, div_id='plot7')}
            </div>

            <!-- ── SNP panels (8–10, conditional) ── -->
            {snp_section_html}

            <!-- ── Footer ── -->
            <div class="footer">
                <p><strong>GSDMB Genomics Analysis Pipeline</strong></p>
                <p>Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p style="margin-top: 10px; font-size: 0.9em; color: #999;">
                    Interactive visualisations powered by Plotly &nbsp;|&nbsp;
                    All plots support zoom, pan, hover, and legend filtering
                </p>
            </div>

        </div>
    </body>
    </html>
    """

    # ── Write HTML to file ────────────────────────────────────────────────────
    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html_content)

    file_size_mb = os.path.getsize(OUTPUT_HTML) / 1024 / 1024
    n_plots      = 10 if df_snp is not None else 7

    print(f"\n{'=' * 70}")
    print("Dashboard created successfully.")
    print(f"  Output file   : {OUTPUT_HTML}")
    print(f"  File size     : {file_size_mb:.2f} MB")
    print(f"  Total plots   : {n_plots}")
    print(f"\n  Open {OUTPUT_HTML} in any web browser to explore.")
    print(f"{'=' * 70}\n")


# ── Script entry point ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    create_dashboard()
