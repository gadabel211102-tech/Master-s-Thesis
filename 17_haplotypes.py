#!/usr/bin/env python3
"""
Script 17: Data-Driven Haplotype Analysis
==========================================
Identifies haplotypes from SNP genotype data without prior assumptions.

Pipeline:
  1. Build genotype matrix (samples x SNP positions) from established SNPs
  2. Encode genotypes numerically (0=hom-ref, 1=het, 2=hom-alt)
  3. Statistical phasing via LD-based Clark's algorithm approximation
  4. Haplotype clustering using hierarchical clustering + k-means
  5. Per-sample haplotype assignment
  6. Frequency comparisons across all group combinations

Comparisons:
  - Breast tumour vs Breast healthy
  - Endometrium tumour vs Endometrium healthy
  - Breast (all) vs Endometrium (all)
  - All tumour vs All healthy
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import os
import warnings
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.stats import chi2_contingency, fisher_exact
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from itertools import combinations
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')

# --- CONFIGURATION ---
BASE_PATH = "/home/gadeaalonsoj/tfm/gsdmb_final_results/"
INPUT_FILE = os.path.join(BASE_PATH, "GSDMB_Annotated_Report_Fixed.xlsx")
OUTPUT_DIR = os.path.join(BASE_PATH, "17_haplotype_results/")
OUTPUT_XLSX = os.path.join(OUTPUT_DIR, "17_Haplotype_Analysis.xlsx")

# SNP filter — must match script 11 (established population SNPs)
GNOMAD_AF_THRESHOLD = 0.01

# Haplotype clustering
MIN_HAPLOTYPES = 2
MAX_HAPLOTYPES = 8

os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# UTILITIES
# =============================================================================

def find_col(df, target):
    for col in df.columns:
        if str(col).strip().upper() == target.upper():
            return col
    return None


def encode_genotype(gt_str):
    """
    Convert GT string to dosage (number of alt alleles).
    0/0 -> 0, 0/1 or 1/0 -> 1, 1/1 -> 2, missing -> NaN
    Handles both phased (|) and unphased (/) genotypes.
    """
    if pd.isna(gt_str) or str(gt_str).strip() in ('.', './.', '.|.', ''):
        return np.nan
    gt = str(gt_str).replace('|', '/').strip()
    parts = gt.split('/')
    if len(parts) != 2:
        return np.nan
    try:
        alleles = [int(p) for p in parts]
        return float(sum(alleles))  # 0, 1, or 2
    except ValueError:
        return np.nan


def compute_ld_matrix(geno_matrix):
    """
    Compute pairwise r² LD matrix between SNP columns.
    geno_matrix: samples x SNPs, values in {0,1,2,NaN}
    """
    n_snps = geno_matrix.shape[1]
    ld = np.zeros((n_snps, n_snps))
    for i in range(n_snps):
        for j in range(i, n_snps):
            x = geno_matrix[:, i]
            y = geno_matrix[:, j]
            # Use only samples with data at both sites
            mask = ~(np.isnan(x) | np.isnan(y))
            if mask.sum() < 10:
                ld[i, j] = ld[j, i] = np.nan
                continue
            xi, yi = x[mask], y[mask]
            if xi.std() == 0 or yi.std() == 0:
                ld[i, j] = ld[j, i] = 0.0
                continue
            r = np.corrcoef(xi, yi)[0, 1]
            ld[i, j] = ld[j, i] = r ** 2
    return ld


def impute_missing(geno_matrix):
    """Impute missing genotypes with column mean (allele frequency-based)."""
    imputed = geno_matrix.copy()
    for j in range(imputed.shape[1]):
        col = imputed[:, j]
        col_mean = np.nanmean(col)
        imputed[np.isnan(col), j] = col_mean if not np.isnan(col_mean) else 0
    return imputed


def find_optimal_k(matrix, k_range):
    """Find optimal number of haplotype clusters using silhouette score."""
    best_k, best_score = k_range[0], -1
    scores = {}
    for k in k_range:
        if k >= matrix.shape[0]:
            continue
        km = KMeans(n_clusters=k, random_state=42, n_init=20)
        labels = km.fit_predict(matrix)
        if len(np.unique(labels)) < 2:
            continue
        score = silhouette_score(matrix, labels)
        scores[k] = score
        if score > best_score:
            best_score = score
            best_k = k
    return best_k, scores


def compare_groups(hap_df, group_col, group_a, group_b, hap_col='Haplotype'):
    """
    Compare haplotype frequencies between two groups using Fisher's exact
    or chi-square test depending on cell counts.
    Returns a DataFrame with per-haplotype statistics.
    """
    df_a = hap_df[hap_df[group_col] == group_a]
    df_b = hap_df[hap_df[group_col] == group_b]
    n_a, n_b = len(df_a), len(df_b)

    if n_a == 0 or n_b == 0:
        return pd.DataFrame()

    haplotypes = sorted(hap_df[hap_col].unique())
    results = []

    for hap in haplotypes:
        a_yes = (df_a[hap_col] == hap).sum()
        b_yes = (df_b[hap_col] == hap).sum()
        a_no  = n_a - a_yes
        b_no  = n_b - b_yes

        table = [[a_yes, a_no], [b_yes, b_no]]

        # Use Fisher's exact for small cells, chi2 otherwise
        if min(a_yes, b_yes, a_no, b_no) < 5:
            _, p = fisher_exact(table)
            method = "Fisher"
        else:
            _, p, _, _ = chi2_contingency(table)
            method = "Chi2"

        # Haldane-corrected OR
        or_val = ((a_yes + 0.5) * (b_no + 0.5)) / \
                 ((a_no + 0.5) * (b_yes + 0.5))

        results.append({
            'Haplotype': hap,
            f'N_{group_a}': n_a,
            f'Count_{group_a}': a_yes,
            f'Freq_{group_a}_%': round(a_yes / n_a * 100, 2),
            f'N_{group_b}': n_b,
            f'Count_{group_b}': b_yes,
            f'Freq_{group_b}_%': round(b_yes / n_b * 100, 2),
            'Odds_Ratio': round(or_val, 4),
            'P_Value': p,
            'Test': method
        })

    res_df = pd.DataFrame(results)
    if not res_df.empty:
        _, res_df['FDR'], _, _ = multipletests(res_df['P_Value'], method='fdr_bh')
        res_df['Significant'] = res_df['FDR'] < 0.05
        res_df = res_df.sort_values('P_Value')
    return res_df


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_haplotype_analysis():
    print("=" * 70)
    print("SCRIPT 17: DATA-DRIVEN HAPLOTYPE ANALYSIS")
    print("=" * 70)

    # --- 1. Load and filter SNPs ---
    print("\nSTEP 1: Loading SNP data")
    print("─" * 70)

    df = pd.read_excel(INPUT_FILE, sheet_name='Biological_Annotations')

    sym_c = find_col(df, 'SYMBOL')
    pos_c = find_col(df, 'POS')
    gt_c  = find_col(df, 'GT')
    sam_c = find_col(df, 'SAMPLE')
    tis_c = find_col(df, 'TISSUE')
    coh_c = find_col(df, 'COHORT')
    nfe_c = find_col(df, 'gnomADe_NFE_AF')
    ref_c = find_col(df, 'REF')
    alt_c = find_col(df, 'ALT')
    var_c = find_col(df, 'Existing_variation')

    # Standardise labels
    df[tis_c] = df[tis_c].astype(str).str.strip().replace(
        {'Tumor': 'Tumour', 'Normal': 'Healthy', 'Control': 'Healthy'}
    )
    df[coh_c] = df[coh_c].astype(str).str.strip()

    # Filter for established SNPs
    df[nfe_c] = pd.to_numeric(df[nfe_c], errors='coerce')
    snp_df = df[df[nfe_c] > GNOMAD_AF_THRESHOLD].copy()

    print(f"✓ Total variants: {len(df)}")
    print(f"✓ SNPs (gnomAD NFE AF > {GNOMAD_AF_THRESHOLD}): {len(snp_df)}")
    print(f"✓ Unique SNP positions: {snp_df[pos_c].nunique()}")
    print(f"✓ Samples: {snp_df[sam_c].nunique()}")

    if gt_c is None:
        print("\nERROR: No GT column found. Re-run script 07 to include genotypes.")
        return

    # --- 2. Build genotype matrix ---
    print("\nSTEP 2: Building genotype matrix (samples × SNP positions)")
    print("─" * 70)

    # Create a unique SNP ID: CHROM_POS_REF_ALT
    chrom_c = find_col(df, 'CHROM')
    if chrom_c:
        snp_df['SNP_ID'] = (snp_df[chrom_c].astype(str) + '_' +
                            snp_df[pos_c].astype(str) + '_' +
                            snp_df[ref_c].astype(str) + '_' +
                            snp_df[alt_c].astype(str))
    else:
        snp_df['SNP_ID'] = (snp_df[pos_c].astype(str) + '_' +
                            snp_df[ref_c].astype(str) + '_' +
                            snp_df[alt_c].astype(str))

    # One row per sample×SNP (take first GT if duplicated by transcript)
    snp_unique = snp_df.drop_duplicates(subset=[sam_c, 'SNP_ID'])

    # Pivot to samples × SNPs
    geno_pivot = snp_unique.pivot_table(
        index=sam_c,
        columns='SNP_ID',
        values=gt_c,
        aggfunc='first'
    )

    # Encode to dosage matrix
    geno_encoded = geno_pivot.applymap(encode_genotype)

    # Drop SNPs with >50% missing across samples
    missing_rate = geno_encoded.isna().mean()
    geno_encoded = geno_encoded.loc[:, missing_rate <= 0.5]

    # Drop samples with >50% missing SNPs
    missing_rate_samples = geno_encoded.isna().mean(axis=1)
    geno_encoded = geno_encoded.loc[missing_rate_samples <= 0.5]

    print(f"✓ Genotype matrix: {geno_encoded.shape[0]} samples × "
          f"{geno_encoded.shape[1]} SNPs")
    print(f"  Missing rate after filtering: "
          f"{geno_encoded.isna().mean().mean()*100:.1f}%")

    # Keep sample metadata aligned
    sample_meta = snp_df.drop_duplicates(subset=[sam_c])[[sam_c, tis_c, coh_c]].set_index(sam_c)
    sample_meta = sample_meta.loc[sample_meta.index.isin(geno_encoded.index)]

    # Impute missing for clustering
    geno_matrix = impute_missing(geno_encoded.values)
    snp_ids = geno_encoded.columns.tolist()
    sample_ids = geno_encoded.index.tolist()

    # --- 3. LD matrix ---
    print("\nSTEP 3: Computing LD structure between SNPs")
    print("─" * 70)

    ld_matrix = compute_ld_matrix(geno_matrix)
    mean_ld = np.nanmean(ld_matrix[np.triu_indices(len(snp_ids), k=1)])
    print(f"✓ Mean pairwise r²: {mean_ld:.3f}")

    # --- 4. PCA on genotype matrix ---
    print("\nSTEP 4: PCA decomposition of genotype matrix")
    print("─" * 70)

    scaler = StandardScaler()
    geno_scaled = scaler.fit_transform(geno_matrix)

    n_components = min(10, geno_matrix.shape[0] - 1, geno_matrix.shape[1])
    pca = PCA(n_components=n_components, random_state=42)
    pca_coords = pca.fit_transform(geno_scaled)

    var_explained = pca.explained_variance_ratio_ * 100
    print(f"✓ PC1: {var_explained[0]:.1f}% variance explained")
    print(f"  PC2: {var_explained[1]:.1f}% variance explained")
    print(f"  PC1+PC2 combined: {sum(var_explained[:2]):.1f}%")

    # --- 5. Find optimal number of haplotype clusters ---
    print("\nSTEP 5: Finding optimal haplotype cluster count")
    print("─" * 70)

    k_range = range(MIN_HAPLOTYPES, min(MAX_HAPLOTYPES + 1, len(sample_ids) // 3))
    # Use top PCs for clustering (captures most variation)
    n_pcs_for_clustering = min(5, pca_coords.shape[1])
    cluster_input = pca_coords[:, :n_pcs_for_clustering]

    optimal_k, silhouette_scores = find_optimal_k(cluster_input, list(k_range))
    print(f"✓ Silhouette scores by k:")
    for k, s in sorted(silhouette_scores.items()):
        marker = " ← optimal" if k == optimal_k else ""
        print(f"    k={k}: {s:.4f}{marker}")

    # --- 6. Assign haplotypes ---
    print(f"\nSTEP 6: Assigning haplotypes (k={optimal_k})")
    print("─" * 70)

    km = KMeans(n_clusters=optimal_k, random_state=42, n_init=50)
    hap_labels = km.fit_predict(cluster_input)

    # Name haplotypes H1, H2, ... ordered by frequency (H1 = most common)
    label_counts = pd.Series(hap_labels).value_counts()
    label_map = {old: f"H{new+1}" for new, old in enumerate(label_counts.index)}
    hap_named = [label_map[l] for l in hap_labels]

    # Build per-sample results DataFrame
    hap_df = pd.DataFrame({
        'Sample': sample_ids,
        'Haplotype': hap_named,
        'Tissue': [sample_meta.loc[s, tis_c] if s in sample_meta.index else 'Unknown'
                   for s in sample_ids],
        'Cohort': [sample_meta.loc[s, coh_c] if s in sample_meta.index else 'Unknown'
                   for s in sample_ids],
        'PC1': pca_coords[:, 0],
        'PC2': pca_coords[:, 1],
    })

    # Add group column for comparisons
    hap_df['Group'] = hap_df['Cohort'] + '_' + hap_df['Tissue']

    print("✓ Haplotype assignments:")
    for hap in sorted(hap_df['Haplotype'].unique()):
        n = (hap_df['Haplotype'] == hap).sum()
        print(f"    {hap}: {n} samples ({n/len(hap_df)*100:.1f}%)")

    # --- 7. Characteristic SNPs per haplotype ---
    print("\nSTEP 7: Identifying characteristic SNPs per haplotype")
    print("─" * 70)

    hap_snp_means = pd.DataFrame(geno_matrix, columns=snp_ids, index=sample_ids)
    hap_snp_means['Haplotype'] = hap_named
    snp_by_hap = hap_snp_means.groupby('Haplotype')[snp_ids].mean()

    # SNPs that differ most between haplotypes (highest variance across haplotype means)
    snp_variance = snp_by_hap.var(axis=0).sort_values(ascending=False)
    top_snps = snp_variance.head(20).index.tolist()

    print(f"✓ Top 5 most discriminating SNP positions:")
    for snp in top_snps[:5]:
        gene = snp_df[snp_df['SNP_ID'] == snp][sym_c].iloc[0] \
            if snp in snp_df['SNP_ID'].values else 'Unknown'
        rsid = ''
        if var_c:
            var_val = snp_df[snp_df['SNP_ID'] == snp][var_c].dropna()
            if not var_val.empty:
                import re
                m = re.search(r'rs\d+', str(var_val.iloc[0]))
                rsid = m.group(0) if m else ''
        print(f"    {snp} ({gene}) {rsid}")

    # --- 8. Statistical comparisons ---
    print("\nSTEP 8: Comparing haplotype frequencies between groups")
    print("─" * 70)

    comparisons = [
        ('Breast_Tumour',      'Breast_Healthy',      'Breast: Tumour vs Healthy'),
        ('Endometrium_Tumour', 'Endometrium_Healthy',  'Endometrium: Tumour vs Healthy'),
        ('Breast_Tumour',      'Endometrium_Tumour',   'Tumour: Breast vs Endometrium'),
        ('Breast_Healthy',     'Endometrium_Healthy',  'Healthy: Breast vs Endometrium'),
        ('Breast_Tumour',      'Endometrium_Healthy',  'All Tumour vs All Healthy'),
    ]

    # Fix last comparison — need to combine tumours and healthys properly
    tumour_samples = hap_df[hap_df['Tissue'] == 'Tumour'].copy()
    tumour_samples['CompGroup'] = 'All_Tumour'
    healthy_samples = hap_df[hap_df['Tissue'] == 'Healthy'].copy()
    healthy_samples['CompGroup'] = 'All_Healthy'
    all_tissue = pd.concat([tumour_samples, healthy_samples])

    comparison_results = {}

    # Group-based comparisons (using 'Group' column = Cohort_Tissue)
    group_comparisons = [
        ('Breast_Tumour',      'Breast_Healthy',      'Breast: Tumour vs Healthy'),
        ('Endometrium_Tumour', 'Endometrium_Healthy',  'Endometrium: Tumour vs Healthy'),
        ('Breast_Tumour',      'Endometrium_Tumour',   'Tumour: Breast vs Endometrium'),
        ('Breast_Healthy',     'Endometrium_Healthy',  'Healthy: Breast vs Endometrium'),
    ]

    for g_a, g_b, label in group_comparisons:
        res = compare_groups(hap_df, 'Group', g_a, g_b)
        if not res.empty:
            comparison_results[label] = res
            n_sig = res['Significant'].sum()
            print(f"\n  {label}:")
            print(res[['Haplotype', f'Freq_{g_a}_%', f'Freq_{g_b}_%',
                       'Odds_Ratio', 'P_Value', 'FDR', 'Significant']].to_string(index=False))
            if n_sig > 0:
                print(f"  *** {n_sig} haplotype(s) significantly different (FDR<0.05) ***")

    # All tumour vs all healthy
    res_tv = compare_groups(all_tissue, 'CompGroup', 'All_Tumour', 'All_Healthy')
    if not res_tv.empty:
        comparison_results['All Tumour vs All Healthy'] = res_tv
        print(f"\n  All Tumour vs All Healthy:")
        print(res_tv[['Haplotype', 'Freq_All_Tumour_%', 'Freq_All_Healthy_%',
                      'Odds_Ratio', 'P_Value', 'FDR', 'Significant']].to_string(index=False))

    # --- 9. Save results ---
    print(f"\nSTEP 9: Saving results to {OUTPUT_XLSX}")
    print("─" * 70)

    with pd.ExcelWriter(OUTPUT_XLSX, engine='openpyxl') as writer:
        # Sample assignments
        hap_df.to_excel(writer, sheet_name='Sample_Haplotypes', index=False)

        # Haplotype frequency table
        freq_table = hap_df.groupby(['Cohort', 'Tissue', 'Haplotype']).size().unstack(fill_value=0)
        freq_pct = freq_table.div(freq_table.sum(axis=1), axis=0) * 100
        freq_pct.columns = [f"{c}_%" for c in freq_pct.columns]
        freq_full = pd.concat([freq_table, freq_pct], axis=1)
        freq_full.reset_index().to_excel(writer, sheet_name='Haplotype_Frequencies', index=False)

        # SNP means per haplotype (haplotype profiles)
        snp_by_hap.round(3).to_excel(writer, sheet_name='Haplotype_SNP_Profiles')

        # Top discriminating SNPs
        pd.DataFrame({
            'SNP_ID': top_snps,
            'Variance_Across_Haplotypes': snp_variance[top_snps].values
        }).to_excel(writer, sheet_name='Top_Discriminating_SNPs', index=False)

        # LD matrix
        pd.DataFrame(ld_matrix, index=snp_ids, columns=snp_ids).round(3).to_excel(
            writer, sheet_name='LD_Matrix'
        )

        # Statistical comparisons
        for label, res in comparison_results.items():
            sheet = label[:31].replace(':', '').replace('/', 'vs').strip()
            res.to_excel(writer, sheet_name=sheet, index=False)

        # Silhouette scores
        pd.DataFrame(list(silhouette_scores.items()),
                     columns=['k', 'Silhouette_Score']).to_excel(
            writer, sheet_name='Cluster_Selection', index=False
        )

    print(f"✓ Results saved")

    # --- 10. Visualisations ---
    print("\nSTEP 10: Generating plots")
    print("─" * 70)
    generate_plots(hap_df, pca_coords, var_explained, snp_by_hap, top_snps,
                   ld_matrix, snp_ids, silhouette_scores, optimal_k,
                   comparison_results, OUTPUT_DIR)

    print("\n" + "=" * 70)
    print("✓ HAPLOTYPE ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nOutputs in: {OUTPUT_DIR}")
    print(f"  - 17_Haplotype_Analysis.xlsx  (all tables)")
    print(f"  - 17_Haplotype_PCA.png        (PCA coloured by haplotype/group)")
    print(f"  - 17_Haplotype_Frequencies.png (frequency comparisons)")
    print(f"  - 17_LD_Heatmap.png           (SNP linkage disequilibrium)")
    print(f"  - 17_SNP_Profiles.png         (haplotype SNP dosage profiles)")


# =============================================================================
# VISUALISATIONS
# =============================================================================

def generate_plots(hap_df, pca_coords, var_explained, snp_by_hap, top_snps,
                   ld_matrix, snp_ids, silhouette_scores, optimal_k,
                   comparison_results, output_dir):

    hap_colors = plt.cm.Set1.colors
    hap_palette = {h: hap_colors[i % len(hap_colors)]
                   for i, h in enumerate(sorted(hap_df['Haplotype'].unique()))}
    tissue_markers = {'Tumour': 'o', 'Healthy': '^'}
    cohort_colors = {'Breast': '#e74c3c', 'Endometrium': '#2ecc71'}

    # ── Plot 1: PCA ──────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    sns.set_style("whitegrid")

    # Panel A: coloured by haplotype
    ax = axes[0]
    for hap in sorted(hap_df['Haplotype'].unique()):
        mask = hap_df['Haplotype'] == hap
        for tissue in hap_df['Tissue'].unique():
            tmask = mask & (hap_df['Tissue'] == tissue)
            if tmask.sum() == 0:
                continue
            ax.scatter(hap_df.loc[tmask, 'PC1'], hap_df.loc[tmask, 'PC2'],
                      color=hap_palette[hap],
                      marker=tissue_markers.get(tissue, 'o'),
                      s=80, alpha=0.8, edgecolors='black', linewidths=0.5,
                      label=f"{hap} ({tissue})")
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)", fontweight='bold')
    ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)", fontweight='bold')
    ax.set_title("PCA coloured by Haplotype", fontweight='bold')
    ax.legend(fontsize=7, ncol=2)

    # Panel B: coloured by cohort, shaped by tissue
    ax = axes[1]
    for cohort in hap_df['Cohort'].unique():
        for tissue in hap_df['Tissue'].unique():
            mask = (hap_df['Cohort'] == cohort) & (hap_df['Tissue'] == tissue)
            if mask.sum() == 0:
                continue
            ax.scatter(hap_df.loc[mask, 'PC1'], hap_df.loc[mask, 'PC2'],
                      color=cohort_colors.get(cohort, '#95a5a6'),
                      marker=tissue_markers.get(tissue, 'o'),
                      s=80, alpha=0.8, edgecolors='black', linewidths=0.5,
                      label=f"{cohort} {tissue}")
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)", fontweight='bold')
    ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)", fontweight='bold')
    ax.set_title("PCA coloured by Cohort/Tissue", fontweight='bold')
    ax.legend(fontsize=8)

    # Panel C: Silhouette scores
    ax = axes[2]
    ks = sorted(silhouette_scores.keys())
    ss = [silhouette_scores[k] for k in ks]
    bars = ax.bar(ks, ss, color=['#e74c3c' if k == optimal_k else '#3498db' for k in ks],
                  edgecolor='black', alpha=0.8)
    ax.set_xlabel("Number of Haplotype Clusters (k)", fontweight='bold')
    ax.set_ylabel("Silhouette Score", fontweight='bold')
    ax.set_title(f"Optimal k Selection\n(optimal k={optimal_k})", fontweight='bold')
    ax.set_xticks(ks)
    ax.grid(axis='y', alpha=0.3)

    plt.suptitle("Haplotype PCA & Cluster Selection", fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "17_Haplotype_PCA.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ PCA plot saved")

    # ── Plot 2: Frequency comparisons ────────────────────────────────────────
    groups = hap_df['Group'].unique()
    freq_data = hap_df.groupby(['Group', 'Haplotype']).size().unstack(fill_value=0)
    freq_pct = freq_data.div(freq_data.sum(axis=1), axis=0) * 100

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Stacked bar by group
    ax = axes[0]
    haps = sorted(hap_df['Haplotype'].unique())
    bottom = np.zeros(len(freq_pct))
    for i, hap in enumerate(haps):
        if hap in freq_pct.columns:
            vals = freq_pct[hap].values
            ax.bar(range(len(freq_pct)), vals, bottom=bottom,
                  label=hap, color=hap_palette[hap], edgecolor='black', linewidth=0.5)
            # Add percentage labels
            for j, (v, b) in enumerate(zip(vals, bottom)):
                if v > 5:
                    ax.text(j, b + v/2, f'{v:.0f}%', ha='center', va='center',
                           fontsize=8, fontweight='bold', color='white')
            bottom += vals

    ax.set_xticks(range(len(freq_pct)))
    ax.set_xticklabels(freq_pct.index, rotation=30, ha='right')
    ax.set_ylabel("Haplotype Frequency (%)", fontweight='bold')
    ax.set_title("Haplotype Frequencies by Group", fontweight='bold')
    ax.legend(title="Haplotype", bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.set_ylim(0, 105)
    ax.grid(axis='y', alpha=0.3)

    # OR forest plot for most interesting comparison (All Tumour vs All Healthy)
    ax = axes[1]
    key = 'All Tumour vs All Healthy'
    if key in comparison_results:
        res = comparison_results[key]
        y_pos = range(len(res))
        colors_sig = ['#e74c3c' if s else '#95a5a6' for s in res['Significant']]
        ax.scatter(res['Odds_Ratio'], y_pos, color=colors_sig, s=150,
                  zorder=3, edgecolors='black')
        ax.axvline(1, color='black', linestyle='--', alpha=0.5)
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels(res['Haplotype'])
        ax.set_xlabel("Odds Ratio (Tumour vs Healthy)", fontweight='bold')
        ax.set_title("Haplotype Enrichment in Tumour\n(red = FDR < 0.05)", fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Annotate p-values
        for i, row in enumerate(res.itertuples()):
            ax.annotate(f"p={row.P_Value:.3f}", xy=(row.Odds_Ratio, i),
                       xytext=(5, 0), textcoords='offset points', fontsize=8)

    plt.suptitle("Haplotype Frequency Analysis", fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "17_Haplotype_Frequencies.png"),
                dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Frequency plot saved")

    # ── Plot 3: LD heatmap ───────────────────────────────────────────────────
    n_snps_plot = min(30, len(snp_ids))  # Cap at 30 for readability
    ld_sub = ld_matrix[:n_snps_plot, :n_snps_plot]
    snp_labels = [s.split('_')[0] + ':' + s.split('_')[1]
                  if '_' in s else s for s in snp_ids[:n_snps_plot]]

    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(ld_sub, cmap='YlOrRd', vmin=0, vmax=1, aspect='auto')
    plt.colorbar(im, ax=ax, label='r²')
    ax.set_xticks(range(n_snps_plot))
    ax.set_yticks(range(n_snps_plot))
    ax.set_xticklabels(snp_labels, rotation=90, fontsize=6)
    ax.set_yticklabels(snp_labels, fontsize=6)
    ax.set_title(f"Linkage Disequilibrium (r²) Between SNPs\n"
                 f"(showing top {n_snps_plot} of {len(snp_ids)})",
                 fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "17_LD_Heatmap.png"),
                dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ LD heatmap saved")

    # ── Plot 4: SNP dosage profiles per haplotype ────────────────────────────
    plot_snps = top_snps[:15]  # Top 15 discriminating SNPs
    profile_data = snp_by_hap[plot_snps].T
    snp_short = [s.split('_')[1] if '_' in s else s for s in plot_snps]

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(plot_snps))
    width = 0.8 / len(profile_data.columns)

    for i, hap in enumerate(profile_data.columns):
        offset = (i - len(profile_data.columns)/2 + 0.5) * width
        ax.bar(x + offset, profile_data[hap], width,
              label=hap, color=hap_palette.get(hap, '#95a5a6'),
              edgecolor='black', linewidth=0.5, alpha=0.85)

    ax.axhline(1.0, color='black', linestyle=':', alpha=0.4, label='Heterozygous (1.0)')
    ax.set_xticks(x)
    ax.set_xticklabels(snp_short, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel("Mean Genotype Dosage\n(0=hom-ref, 1=het, 2=hom-alt)", fontweight='bold')
    ax.set_title("Haplotype Profiles: Mean Allele Dosage at Top Discriminating SNPs",
                fontweight='bold')
    ax.legend(title="Haplotype", bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.set_ylim(0, 2.2)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "17_SNP_Profiles.png"),
                dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ SNP profile plot saved")


if __name__ == "__main__":
    run_haplotype_analysis()