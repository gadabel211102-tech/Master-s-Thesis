#!/usr/bin/env python3
"""
Generates a self-contained interactive HTML dashboard for exploring the
annotated variant dataset produced by scripts 07–13. All plots are rendered
with Plotly and embedded in a single HTML file that can be opened in any
web browser without additional software.

The dashboard includes up to ten interactive plots:
  1–7: Core variant analysis (genomic landscape, hierarchy, burden, heatmaps,
       gene frequency, 3D explorer, data table)
  8–10: SNP-specific panels loaded from script 11 output, if available

Plots 2–10 load Plotly from the CDN reference injected by plot 1, keeping the
total file size manageable.
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import os
import warnings
warnings.filterwarnings('ignore')

# --- Configuration -----------------------------------------------------------
BASE_PATH   = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
INPUT_FILE  = os.path.join(BASE_PATH, "GSDMB_Annotated_Report_Fixed.xlsx")
OUTPUT_HTML = os.path.join(BASE_PATH, "14_Interactive_Dashboard.html")


def find_col(df: pd.DataFrame, target: str):
    """Case-insensitive column lookup — VEP field names vary in capitalisation between cache versions."""
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None


def create_interactive_scatter(df, sym_c, pos_c, imp_c, coh_c, tis_c):
    """
    Plot 1: Genomic position vs gnomAD NFE allele frequency, faceted by
    cohort and tissue type. Each point is one variant; colour indicates gene
    and marker shape indicates VEP impact class. Log scale on the y-axis
    separates the dense cluster of common variants from rare outliers.
    """
    fig = px.scatter(
        df,
        x=pos_c, y='gnomAD_NFE_AF',
        color=sym_c, symbol=imp_c,
        facet_row=coh_c, facet_col=tis_c,
        hover_data=['HGVSp', 'Consequence', 'Sample', 'gnomAD_NFE_Source'],
        title='Interactive Genomic Landscape',
        labels={
            pos_c: 'Genomic Position (Chr17)',
            'gnomAD_NFE_AF': 'gnomAD NFE AF',
            sym_c: 'Gene', imp_c: 'Impact'
        },
        height=800,
        color_discrete_sequence=px.colors.qualitative.Set3
    )
    fig.update_layout(hovermode='closest', font=dict(size=11), title_font_size=16)
    fig.update_xaxes(title_text="Position on Chr17")
    fig.update_yaxes(type="log", title_text="gnomAD NFE AF (log)")
    return fig


def create_variant_sunburst(df, sym_c, imp_c, con_c):
    """
    Plot 2: Sunburst chart showing the hierarchical breakdown of variants
    by gene -> VEP impact -> consequence. Segment size is proportional to
    variant count; colour intensity also encodes count. This provides an
    at-a-glance overview of where the functional burden lies across genes.
    """
    sunburst_data = df.groupby([sym_c, imp_c, con_c]).size().reset_index(name='Count')
    fig = px.sunburst(
        sunburst_data,
        path=[sym_c, imp_c, con_c], values='Count',
        title='Variant Hierarchy: Gene -> Impact -> Consequence',
        color='Count', color_continuous_scale='RdYlBu_r',
        height=700
    )
    fig.update_traces(textinfo='label+percent parent')
    fig.update_layout(font=dict(size=12))
    return fig


def create_cohort_comparison_box(df, coh_c, tis_c):
    """
    Plot 3: Box plot of per-sample variant counts grouped by cohort and
    tissue type. Individual sample points are overlaid (points='all') so
    that outlier samples are visible alongside the distributional summary.
    """
    sample_counts = df.groupby(['Sample', coh_c, tis_c]).size().reset_index(name='Variant_Count')
    fig = px.box(
        sample_counts,
        x=coh_c, y='Variant_Count', color=tis_c,
        title='Variant Burden Distribution by Cohort and Tissue',
        labels={coh_c: 'Cancer Type', 'Variant_Count': 'Number of Variants', tis_c: 'Tissue Type'},
        points='all', hover_data=['Sample'],
        height=500,
        color_discrete_map={'Tumour': '#e74c3c', 'Healthy': '#3498db'}
    )
    fig.update_layout(boxmode='group', font=dict(size=12))
    return fig


def create_impact_consequence_heatmap(df, imp_c, con_c):
    """
    Plot 4: Heatmap of VEP impact category (rows) vs functional consequence
    (columns), with cell values showing variant counts. Impact rows are
    ordered from most to least severe (HIGH -> MODIFIER) to match the
    ordering convention used throughout the analysis.
    """
    heatmap_data = df.groupby([imp_c, con_c]).size().reset_index(name='Count')
    pivot_data = heatmap_data.pivot(index=imp_c, columns=con_c, values='Count').fillna(0)

    impact_order = ['HIGH', 'MODERATE', 'MODIFIER', 'LOW']
    pivot_data = pivot_data.reindex([i for i in impact_order if i in pivot_data.index])

    fig = go.Figure(data=go.Heatmap(
        z=pivot_data.values,
        x=pivot_data.columns, y=pivot_data.index,
        colorscale='YlOrRd',
        text=pivot_data.values, texttemplate='%{text}',
        textfont={"size": 10}, hoverongaps=False,
        colorbar=dict(title="Count")
    ))
    fig.update_layout(
        title='Impact vs Consequence Heatmap',
        xaxis_title='Functional Consequence', yaxis_title='Variant Impact',
        height=500, font=dict(size=11)
    )
    fig.update_xaxes(tickangle=-45)
    return fig


def create_gene_coverage_timeline(df, sym_c, sam_c, coh_c):
    """
    Plot 5: Grouped bar chart showing the number of distinct samples carrying
    at least one variant in each of the top 10 most frequently mutated genes,
    split by cohort. Highlights genes with cohort-specific enrichment patterns.
    """
    top_genes = df[sym_c].value_counts().head(10).index
    timeline_data = []
    for gene in top_genes:
        gene_data = df[df[sym_c] == gene]
        for cohort in df[coh_c].unique():
            n_samples = gene_data[gene_data[coh_c] == cohort][sam_c].nunique()
            timeline_data.append({'Gene': gene, 'Cohort': cohort, 'Sample_Count': n_samples})

    fig = px.bar(
        pd.DataFrame(timeline_data),
        x='Gene', y='Sample_Count', color='Cohort',
        title='Top 10 Genes: Sample Coverage by Cohort',
        labels={'Sample_Count': 'Number of Samples', 'Gene': 'Gene'},
        barmode='group', height=500,
        color_discrete_sequence=px.colors.qualitative.Bold
    )
    fig.update_layout(font=dict(size=12))
    return fig


def create_3d_scatter(df, pos_c, sym_c, imp_c):
    """
    Plot 6: 3D scatter exploring the relationship between genomic position
    (x), gnomAD NFE population frequency (y), and sequencing depth (z).
    A maximum of 500 points are sampled for browser performance. Depth is
    capped at the 99th percentile to prevent extreme outliers from compressing
    the visible range. Falls back to an empty placeholder if DP is absent.
    """
    dp_col = next((c for c in df.columns if c.upper() == 'DP'), None)
    if dp_col is None:
        print("  Warning: no DP column found — 3D plot will be skipped.")
        fig = go.Figure()
        fig.update_layout(title='3D Variant Explorer (unavailable — DP column not found)', height=700)
        return fig

    df_3d = df[[pos_c, 'gnomAD_NFE_AF', dp_col, imp_c, sym_c, 'gnomAD_NFE_Source']].copy()
    df_3d['Read_Depth']    = pd.to_numeric(df_3d[dp_col], errors='coerce')
    df_3d['gnomAD_NFE_AF'] = pd.to_numeric(df_3d['gnomAD_NFE_AF'], errors='coerce')
    df_3d = df_3d.dropna(subset=['Read_Depth', 'gnomAD_NFE_AF'])

    if df_3d.empty:
        fig = go.Figure()
        fig.update_layout(title='3D Variant Explorer (unavailable — no valid AF/DP values)', height=700)
        return fig

    # Cap depth at the 99th percentile; extreme outliers (e.g. amplicon spikes)
    # otherwise collapse the visible z-range and obscure the main distribution
    depth_cap = df_3d['Read_Depth'].quantile(0.99)
    df_3d['Read_Depth'] = df_3d['Read_Depth'].clip(upper=depth_cap)

    # Subsample to keep the HTML file size manageable and browser rendering fast
    df_sample = df_3d.sample(min(500, len(df_3d)), random_state=42)

    fig = go.Figure(data=[go.Scatter3d(
        x=df_sample[pos_c],
        y=df_sample['gnomAD_NFE_AF'],
        z=df_sample['Read_Depth'],
        mode='markers',
        marker=dict(
            size=5,
            color=df_sample[imp_c].map({'HIGH': 0, 'MODERATE': 1, 'MODIFIER': 2, 'LOW': 3}),
            colorscale='RdYlGn_r', showscale=True,
            colorbar=dict(title="Impact", tickvals=[0, 1, 2, 3],
                          ticktext=['HIGH', 'MODERATE', 'MODIFIER', 'LOW']),
            line=dict(color='black', width=0.5)
        ),
        text=df_sample[sym_c],
        hovertemplate=(
            '<b>%{text}</b><br>'
            'Position: %{x:,}<br>'
            'gnomAD NFE AF: %{y:.4f}<br>'
            'Read Depth: %{z:.0f}x'
            '<extra></extra>'
        )
    )])
    fig.update_layout(
        title='3D Variant Explorer: Position x Population Frequency x Read Depth',
        scene=dict(xaxis_title='Genomic Position (Chr17)',
                   yaxis_title='gnomAD NFE AF',
                   zaxis_title='Read Depth (DP)'),
        height=700, font=dict(size=11)
    )
    return fig


def create_interactive_table(df, sym_c, pos_c, imp_c, con_c):
    """
    Plot 7: Scrollable data table showing the top 100 variants with the most
    relevant columns. Plotly's go.Table renders in the browser with alternating
    row shading; the 100-row limit keeps the embedded HTML size reasonable.
    """
    table_cols = [sym_c, pos_c, 'HGVSp', imp_c, con_c,
                  'gnomAD_NFE_AF', 'gnomAD_NFE_Source', 'Sample', 'Cohort', 'Tissue']
    table_df = df[table_cols].copy()
    table_df['gnomAD_NFE_AF'] = pd.to_numeric(table_df['gnomAD_NFE_AF'], errors='coerce').round(4)
    table_df = table_df.head(100)

    fig = go.Figure(data=[go.Table(
        header=dict(values=list(table_df.columns),
                    fill_color='#3498db', font=dict(color='white', size=12), align='left'),
        cells=dict(values=[table_df[col] for col in table_df.columns],
                   fill_color=[['#ecf0f1', 'white'] * (len(table_df)//2 + 1)],
                   align='left', font=dict(size=11))
    )])
    fig.update_layout(title='Top 100 Variants (Scrollable Table)', height=600)
    return fig


def load_snp_data():
    """
    Load the SNP summary produced by script 11. Returns None if the file is
    absent so that the dashboard still renders without the SNP panels.
    """
    snp_file = os.path.join(BASE_PATH, "11_Master_Unique_SNP_Summary.xlsx")
    if not os.path.exists(snp_file):
        print(f"  Warning: SNP summary not found at {snp_file}")
        print("  Run script 11 first to enable SNP panels.")
        return None
    df_snp = pd.read_excel(snp_file)
    df_snp.columns = df_snp.columns.str.strip()
    print(f"  SNP data loaded: {len(df_snp)} unique SNPs")
    return df_snp


def create_snp_benchmarking(df_snp: pd.DataFrame):
    """
    Plot 8: Scatter plot comparing each SNP's carrier frequency in the study
    cohort against its gnomAD NFE reference frequency. The dashed diagonal
    (y = x * 100) marks the expected position if the study population perfectly
    matches the European reference panel. Points above the line indicate
    enrichment relative to the reference; points below indicate depletion.
    """
    freq_cols = [c for c in df_snp.columns if c.endswith('Frequency_%')]
    if not freq_cols or 'gnomAD_NFE_AF' not in df_snp.columns:
        fig = go.Figure()
        fig.update_layout(title='SNP Benchmarking (required columns not found)', height=600)
        return fig

    id_vars = [c for c in ['Variant_ID', 'SYMBOL', 'Consequence', 'IMPACT',
                            'gnomAD_NFE_AF', 'gnomAD_NFE_Source'] if c in df_snp.columns]
    df_long = df_snp.melt(id_vars=id_vars, value_vars=freq_cols,
                          var_name='Group', value_name='Study_Frequency_%')
    df_long['Cohort'] = df_long['Group'].str.split('_').str[0]
    df_long['Tissue'] = df_long['Group'].str.split('_').str[1]
    df_long = df_long[df_long['Study_Frequency_%'] > 0].copy()

    fig = px.scatter(
        df_long,
        x='gnomAD_NFE_AF', y='Study_Frequency_%',
        color='Tissue', symbol='Cohort',
        hover_data={'Variant_ID': True, 'SYMBOL': True, 'Consequence': True,
                    'IMPACT': True, 'gnomAD_NFE_AF': ':.4f',
                    'gnomAD_NFE_Source': True, 'Study_Frequency_%': ':.1f', 'Group': False},
        color_discrete_map={'Tumour': '#e74c3c', 'Healthy': '#3498db'},
        title='SNP Frequency Benchmarking: Study Cohort vs gnomAD Population',
        labels={'gnomAD_NFE_AF': 'gnomAD NFE AF', 'Study_Frequency_%': 'Carrier Frequency in Study (%)'},
        height=600
    )

    # Diagonal reference line: points on this line match the gnomAD population exactly
    x_range = [0, df_long['gnomAD_NFE_AF'].max() * 1.05]
    fig.add_trace(go.Scatter(
        x=x_range, y=[v * 100 for v in x_range],
        mode='lines', line=dict(color='grey', dash='dash', width=1),
        name='Population baseline (y = x)', hoverinfo='skip'
    ))
    fig.update_layout(hovermode='closest', legend_title='Tissue / Cohort', font=dict(size=11))
    return fig


def create_snp_heatmap(df_snp: pd.DataFrame):
    """
    Plot 9: Carrier frequency heatmap for the 40 SNPs with the highest mean
    carrier frequency across all groups. Rows are SNPs labelled with rsID and
    gene symbol; columns are cohort-tissue groups. Rich hover text includes
    gnomAD AF and consequence for each cell.
    """
    freq_cols = [c for c in df_snp.columns if c.endswith('Frequency_%')]
    if not freq_cols:
        fig = go.Figure()
        fig.update_layout(title='SNP Heatmap (frequency columns not found)', height=600)
        return fig

    df_snp = df_snp.copy()
    df_snp['Mean_Freq'] = df_snp[freq_cols].mean(axis=1)
    top_snps = df_snp.nlargest(40, 'Mean_Freq').copy()

    label_col = 'Variant_ID' if 'Variant_ID' in top_snps.columns else top_snps.index
    top_snps['Label'] = (
        top_snps[label_col].astype(str) + '  (' + top_snps['SYMBOL'].astype(str) + ')'
        if 'SYMBOL' in top_snps.columns else top_snps[label_col].astype(str)
    )

    z_data   = top_snps[freq_cols].values
    x_labels = [c.replace('_Frequency_%', '').replace('_', ' ') for c in freq_cols]
    y_labels = top_snps['Label'].tolist()

    hover_text = []
    for _, row in top_snps.iterrows():
        row_text = []
        for col, x_lbl in zip(freq_cols, x_labels):
            row_text.append(
                f"SNP: {row.get('Variant_ID', '')}<br>"
                f"Gene: {row.get('SYMBOL', '')}<br>"
                f"Group: {x_lbl}<br>"
                f"Carrier freq: {row[col]:.1f}%<br>"
                f"gnomAD NFE AF: {row.get('gnomAD_NFE_AF', 'N/A')}<br>"
                f"NFE source: {row.get('gnomAD_NFE_Source', 'N/A')}<br>"
                f"Consequence: {row.get('Consequence', 'N/A')}"
            )
        hover_text.append(row_text)

    fig = go.Figure(data=go.Heatmap(
        z=z_data, x=x_labels, y=y_labels,
        colorscale='YlOrRd',
        text=[[f'{v:.1f}%' for v in row] for row in z_data],
        texttemplate='%{text}', textfont=dict(size=9),
        hovertext=hover_text, hovertemplate='%{hovertext}<extra></extra>',
        colorbar=dict(title='Carrier<br>Frequency (%)'), zmin=0
    ))
    fig.update_layout(
        title='Top 40 SNPs: Carrier Frequency Across Cohort Groups',
        xaxis_title='Cohort / Tissue Group', yaxis_title='SNP (rsID / Gene)',
        height=900, font=dict(size=11), yaxis=dict(tickfont=dict(size=9))
    )
    return fig


def create_snp_consequence_breakdown(df_snp: pd.DataFrame):
    """
    Plot 10: Horizontal stacked bar chart of SNP counts by functional
    consequence, coloured by gene. Consequences are sorted in ascending order
    of total count so the most common appear at the bottom of the y-axis.
    """
    if 'Consequence' not in df_snp.columns or 'SYMBOL' not in df_snp.columns:
        fig = go.Figure()
        fig.update_layout(title='SNP Consequence Breakdown (columns not found)', height=500)
        return fig

    con_counts  = (df_snp.groupby(['Consequence', 'SYMBOL']).size()
                   .reset_index(name='Count'))
    con_order   = (con_counts.groupby('Consequence')['Count']
                   .sum().sort_values(ascending=True).index.tolist())

    fig = px.bar(
        con_counts,
        x='Count', y='Consequence', color='SYMBOL', orientation='h',
        title='SNP Distribution by Functional Consequence and Gene',
        labels={'Count': 'Number of Unique SNPs',
                'Consequence': 'Functional Consequence', 'SYMBOL': 'Gene'},
        category_orders={'Consequence': con_order},
        height=600, color_discrete_sequence=px.colors.qualitative.Bold
    )
    fig.update_layout(barmode='stack', legend_title='Gene',
                      font=dict(size=11), hovermode='closest')
    return fig


def create_dashboard():
    """
    Orchestrate all plot functions and assemble the outputs into a single
    self-contained HTML file. Plot 1 injects the Plotly CDN script tag
    (include_plotlyjs='cdn'); all subsequent plots omit it (include_plotlyjs=False)
    to avoid embedding the same ~3 MB library multiple times.
    """
    print("=" * 70)
    print("SCRIPT 14: INTERACTIVE VISUALISATION DASHBOARD")
    print("=" * 70)
    print(f"\n  Input:  {INPUT_FILE}")
    print(f"  Output: {OUTPUT_HTML}")
    print("\n" + "=" * 70 + "\n")

    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: input file not found: {INPUT_FILE}"); return

    df = pd.read_excel(INPUT_FILE, sheet_name='Biological_Annotations')

    sym_c        = find_col(df, 'SYMBOL')
    pos_c        = find_col(df, 'Pos')
    imp_c        = find_col(df, 'IMPACT')
    con_c        = find_col(df, 'CONSEQUENCE')
    coh_c        = find_col(df, 'COHORT')
    tis_c        = find_col(df, 'TISSUE')
    exome_nfe_c  = find_col(df, 'gnomADe_NFE_AF')
    genome_nfe_c = find_col(df, 'gnomADg_NFE_AF')

    required = [sym_c, pos_c, imp_c, con_c, coh_c, tis_c, exome_nfe_c, genome_nfe_c]
    if any(c is None for c in required):
        print("ERROR: one or more required columns are missing from the input file.")
        print(f"  SYMBOL={sym_c}, Pos={pos_c}, IMPACT={imp_c}, CONSEQUENCE={con_c}, "
              f"COHORT={coh_c}, TISSUE={tis_c}, gnomADe_NFE_AF={exome_nfe_c}, "
              f"gnomADg_NFE_AF={genome_nfe_c}")
        return

    # Build combined gnomAD NFE AF column (exome preferred; genome as fallback)
    df['gnomAD_NFE_AF'] = df[exome_nfe_c].combine_first(df[genome_nfe_c])
    df['gnomAD_NFE_Source'] = np.where(
        df[exome_nfe_c].notna(), 'Exome_NFE',
        np.where(df[genome_nfe_c].notna(), 'Genome_NFE', 'Missing')
    )
    df[tis_c] = df[tis_c].astype(str).str.strip().replace(
        {'Tumor': 'Tumour', 'Normal': 'Healthy', 'Control': 'Healthy'}
    )
    df[coh_c] = df[coh_c].astype(str).str.strip()

    print(f"Total variants: {len(df)}\n")
    print("Loading SNP summary data...")
    df_snp = load_snp_data()

    print("Building interactive plots:")
    print("  [1/7] Genomic landscape scatter...")
    fig1 = create_interactive_scatter(df, sym_c, pos_c, imp_c, coh_c, tis_c)
    print("  [2/7] Variant hierarchy sunburst...")
    fig2 = create_variant_sunburst(df, sym_c, imp_c, con_c)
    print("  [3/7] Cohort comparison box plots...")
    fig3 = create_cohort_comparison_box(df, coh_c, tis_c)
    print("  [4/7] Impact-Consequence heatmap...")
    fig4 = create_impact_consequence_heatmap(df, imp_c, con_c)
    print("  [5/7] Gene mutation frequency bars...")
    fig5 = create_gene_coverage_timeline(df, sym_c, 'Sample', coh_c)
    print("  [6/7] 3D variant explorer...")
    fig6 = create_3d_scatter(df, pos_c, sym_c, imp_c)
    print("  [7/7] Scrollable data table...")
    fig7 = create_interactive_table(df, sym_c, pos_c, imp_c, con_c)

    snp_section_html    = ''
    snp_stat_card_html  = ''
    if df_snp is not None:
        n_snps = len(df_snp)
        snp_stat_card_html = f"""
                <div class="stat-card">
                    <div class="stat-number">{n_snps:,}</div>
                    <div class="stat-label">Established SNPs (gnomAD NFE AF &gt;1%)</div>
                </div>"""
        print("  [8/10] SNP frequency benchmarking...")
        fig8 = create_snp_benchmarking(df_snp)
        print("  [9/10] SNP carrier frequency heatmap...")
        fig9 = create_snp_heatmap(df_snp)
        print("  [10/10] SNP consequence breakdown...")
        fig10 = create_snp_consequence_breakdown(df_snp)

        snp_section_html = f"""
            <div class="section-divider">
                <h2 class="section-heading">SNP Analysis (Script 11)</h2>
                <p class="section-desc">
                    Established SNPs defined as variants with gnomAD NFE allele frequency above 1%.
                    Carrier frequencies are calculated per cohort and tissue group across all samples.
                </p>
            </div>
            <div class="plot-section">
                <h2 class="plot-title">8. SNP Frequency Benchmarking vs gnomAD Population</h2>
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

    # --- Assemble HTML -------------------------------------------------------
    # Plot 1 is rendered with include_plotlyjs='cdn', which injects a <script> tag
    # pointing to the Plotly CDN. All other plots use include_plotlyjs=False to
    # avoid embedding the library multiple times and bloating the output file.
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Genomics Analysis Dashboard</title>
        <meta charset="utf-8">
        <style>
            body {{
                font-family: 'Arial', sans-serif;
                margin: 0; padding: 20px;
                background-color: #f5f5f5;
            }}
            .header {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white; padding: 30px; text-align: center;
                border-radius: 10px; margin-bottom: 30px;
                box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            }}
            .header h1 {{ margin: 0; font-size: 2.5em; }}
            .header p  {{ margin: 10px 0 0 0; font-size: 1.2em; opacity: 0.9; }}
            .container {{ max-width: 1400px; margin: 0 auto; }}
            .plot-section {{
                background: white; padding: 20px; margin-bottom: 30px;
                border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            .plot-title {{
                color: #333; border-bottom: 3px solid #667eea;
                padding-bottom: 10px; margin-bottom: 20px; font-size: 1.5em;
            }}
            .stats-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 20px; margin-bottom: 30px;
            }}
            .stat-card {{
                background: white; padding: 20px; border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1); text-align: center;
            }}
            .stat-number {{ font-size: 2.5em; font-weight: bold; color: #667eea; }}
            .stat-label  {{ color: #666; margin-top: 10px; }}
            .section-divider {{
                background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
                color: white; padding: 20px 30px; border-radius: 8px;
                margin-bottom: 20px; margin-top: 10px;
            }}
            .section-heading {{ margin: 0 0 8px 0; font-size: 1.6em; }}
            .section-desc    {{ margin: 0; font-size: 0.95em; opacity: 0.92; }}
            .footer {{
                text-align: center; padding: 20px; color: #666;
                border-top: 1px solid #ddd; margin-top: 40px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>Genomics Analysis Interactive Dashboard</h1>
                <p>Comprehensive Variant Analysis and Visualisation</p>
                <p style="font-size: 0.9em; margin-top: 15px;">
                    Dataset: {len(df):,} variants | {df['Sample'].nunique()} samples |
                    {df[sym_c].nunique()} genes
                </p>
            </div>

            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-number">{len(df):,}</div>
                    <div class="stat-label">Total Variants</div>
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
                    <div class="stat-label">High Impact Variants</div>
                </div>
                {snp_stat_card_html}
            </div>

            <div class="plot-section">
                <h2 class="plot-title">1. Genomic Landscape Explorer</h2>
                {fig1.to_html(include_plotlyjs='cdn', div_id='plot1')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">2. Variant Hierarchy (Sunburst)</h2>
                {fig2.to_html(include_plotlyjs=False, div_id='plot2')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">3. Variant Burden by Cohort</h2>
                {fig3.to_html(include_plotlyjs=False, div_id='plot3')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">4. Impact vs Consequence Matrix</h2>
                {fig4.to_html(include_plotlyjs=False, div_id='plot4')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">5. Gene Mutation Frequency</h2>
                {fig5.to_html(include_plotlyjs=False, div_id='plot5')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">6. 3D Variant Explorer</h2>
                {fig6.to_html(include_plotlyjs=False, div_id='plot6')}
            </div>
            <div class="plot-section">
                <h2 class="plot-title">7. Variant Data Table (Top 100)</h2>
                {fig7.to_html(include_plotlyjs=False, div_id='plot7')}
            </div>

            {snp_section_html}

            <div class="footer">
                <p><strong>GSDMB Variant Analysis Pipeline</strong></p>
                <p>Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p style="margin-top: 10px; font-size: 0.9em; color: #999;">
                    Interactive visualisations powered by Plotly |
                    All plots support zoom, pan, and hover for details
                </p>
            </div>
        </div>
    </body>
    </html>
    """

    with open(OUTPUT_HTML, 'w', encoding='utf-8') as f:
        f.write(html_content)

    n_plots = 10 if df_snp is not None else 7
    print(f"\n{'=' * 70}")
    print(f"Dashboard created: {OUTPUT_HTML}")
    print(f"File size: {os.path.getsize(OUTPUT_HTML) / 1024 / 1024:.2f} MB")
    print(f"Total interactive plots: {n_plots}")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    create_dashboard()
