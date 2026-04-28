"""
Script 09: GSDMB-Focused Landscape Mapping
=========================================

Purpose
-------
Generate split genomic landscapes restricted to GSDMB variants so the study
gene can be reviewed without mixing established common SNPs and lower-frequency
variants in the same main figure.

Methodological note
-------------------
Annotation labels come from the stage-07 annotation workbook, while explicit
WT / Het ALT / Hom ALT states come from stage 05b forced genotypes. The
resulting figures therefore use callable denominators instead of deriving
WT from annotation-derived subtraction.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from annotation_genotype_bridge import build_allele_genotype_long, summarise_variant_panel_genotypes
from landscape_mapping_utils import (
    BASE_IMPACT_ORDER,
    DEFAULT_GENOTYPE_COLORS,
    apply_variant_order,
    build_family_manifest,
    plot_genotype_landscape,
    plot_variant_landscape,
)
from pipeline_utils import (
    attach_amplicon_warning_columns,
    build_amplicon_warning_lookup,
    common_nfe_variant_mask,
    ensure_directory,
    get_paths,
    get_thresholds,
)
from pipeline_validation import validate_file_exists

PATHS = get_paths()
THRESHOLDS = get_thresholds()
ANNOTATION_FILE = PATHS['annotated_report']
FORCED_GENOTYPE_FILE = PATHS.get('forced_genotypes_dir', PATHS['results_dir'] / '05b_forced_genotypes') / 'GSDMB_Forced_Genotypes_Union_Sites.xlsx'
output_dir = ensure_directory(PATHS.get('gsdmb_landscape_dir', PATHS['results_dir'] / '09_gsdmb_landscape'))
output_common_image = output_dir / 'GSDMB_Target_Gene_Common_SNP_Landscape.png'
output_rare_image = output_dir / 'GSDMB_Target_Gene_Lower_Frequency_Variant_Landscape.png'
output_common_genotype_image = output_dir / 'GSDMB_Target_Gene_Common_SNP_Genotype_Landscape.png'
output_rare_genotype_image = output_dir / 'GSDMB_Target_Gene_Lower_Frequency_Variant_Genotype_Landscape.png'
output_snp_only_image = output_dir / 'GSDMB_Target_Gene_SNP_Only_Landscape.png'
output_snp_only_genotype_image = output_dir / 'GSDMB_Target_Gene_SNP_Only_Genotype_Landscape.png'
output_excel = output_dir / 'GSDMB_Target_Gene_Variant_Landscape_Split.xlsx'
comparison_significance_file = PATHS.get('landscape_comparison_dir', PATHS['results_dir'] / '09b_landscape_comparison') / 'Landscape_Frequency_Significance_Summary.csv'
LABEL_MODE = 'high_only'
FORCED_GENOTYPE_LABELS = {'rs11078928', 'rs869402'}
IMPACT_ORDER = BASE_IMPACT_ORDER
IMPACT_COLORS = {'HIGH': '#D62728', 'MODERATE': '#FF7F0E', 'LOW': '#2CA02C', 'MODIFIER': '#1F77B4'}


def label_impacts() -> list[str]:
    return ['HIGH'] if LABEL_MODE == 'high_only' else ['HIGH', 'MODERATE']


def build_variant_label(row: pd.Series) -> str:
    rsid = str(row.get('rsID', '') or '').strip()
    if rsid and rsid.lower() != 'nan':
        return rsid
    alt = str(row.get('ALT', '') or '').strip()
    return f"GSDMB:{int(row['POS'])}:{alt}" if alt else str(row.get('Variant_Key', ''))


def select_high_effect_labels(summary: pd.DataFrame) -> list[str]:
    ranked = (
        summary.loc[summary['Impact_Severity'].astype(str).str.upper() == 'HIGH', ['Variant_Label', 'Carrier_Percentage_Callables']]
        .assign(Variant_Label=lambda df: df['Variant_Label'].astype(str).str.strip())
        .query("Variant_Label != '' and Variant_Label != 'nan'")
        .groupby('Variant_Label', as_index=False)['Carrier_Percentage_Callables']
        .max()
        .sort_values(['Carrier_Percentage_Callables', 'Variant_Label'], ascending=[False, True], kind='stable')
    )
    return ranked['Variant_Label'].tolist()


def select_standout_label(summary: pd.DataFrame) -> list[str]:
    ranked = (
        summary[['Variant_Label', 'Carrier_Percentage_Callables']]
        .assign(Variant_Label=lambda df: df['Variant_Label'].astype(str).str.strip())
        .query("Variant_Label != '' and Variant_Label != 'nan'")
        .groupby('Variant_Label', as_index=False)['Carrier_Percentage_Callables']
        .max()
        .sort_values(['Carrier_Percentage_Callables', 'Variant_Label'], ascending=[False, True], kind='stable')
        .reset_index(drop=True)
    )
    if ranked.empty:
        return []
    return [str(ranked.loc[0, 'Variant_Label'])]


def load_comparison_significant_labels(map_name: str) -> set[str]:
    if not comparison_significance_file.exists():
        return set()
    significance_df = pd.read_csv(comparison_significance_file)
    map_df = significance_df[(significance_df['Map'] == map_name) & (significance_df['FDR_BH'].astype(float) < 0.05)]
    labels = map_df['rsID'].astype(str).str.strip()
    return {label for label in labels if label and label.lower() != 'nan'}


def prepare_variant_summary(long_df: pd.DataFrame) -> tuple[pd.DataFrame, dict[tuple[str, str], int]]:
    long_df = long_df[long_df['Gene_Symbol'].astype(str).str.contains('GSDMB', case=False, na=False)].copy()
    sample_size_lookup = long_df.groupby(['Cohort', 'Tissue'])['Sample'].nunique().to_dict()
    summary = summarise_variant_panel_genotypes(long_df, ['Cohort', 'Tissue', 'Variant_Key'])
    summary['Impact_Severity'] = summary['Impact_Severity'].astype(str).str.upper()
    summary['Is_Common_SNP_By_NFE'] = common_nfe_variant_mask(summary, THRESHOLDS['min_nfe_af'])
    summary['Variant_Class'] = summary['Is_Common_SNP_By_NFE'].map({True: 'Established_SNP', False: 'Lower_Frequency_Mutation'})
    summary['Variant_Label'] = summary.apply(build_variant_label, axis=1)
    summary['Labelled_On_Figure'] = summary['Impact_Severity'].isin(label_impacts()) & summary['Variant_Label'].astype(str).ne('')
    warning_lookup = build_amplicon_warning_lookup(summary, variant_col='Variant_Key', chrom_col='CHROM', pos_col='POS')
    summary = attach_amplicon_warning_columns(summary, warning_lookup, variant_col='Variant_Key')
    summary = summary.sort_values(['Cohort', 'Tissue', 'CHROM', 'POS', 'REF', 'ALT'], kind='stable').reset_index(drop=True)
    return summary, sample_size_lookup


def build_genotype_points(summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in summary.iterrows():
        for genotype_label, count_col, pct_col in [
            ('WT', 'N_WT', 'WT_Percentage_Callables'),
            ('Het', 'N_Het_ALT', 'Het_Percentage_Callables'),
            ('Hom', 'N_Hom_ALT', 'Hom_Percentage_Callables'),
        ]:
            rows.append({
                **row.to_dict(),
                'Genotype_Category': genotype_label,
                'Genotype_Count': int(row[count_col]),
                'Genotype_Percentage': row[pct_col],
            })
    out = pd.DataFrame(rows)
    out['Genotype_Category'] = pd.Categorical(out['Genotype_Category'], categories=['WT', 'Het', 'Hom'], ordered=True)
    return out


def select_family(summary: pd.DataFrame, genotype_points: pd.DataFrame, family_name: str, family_note: str, mask: pd.Series) -> dict[str, pd.DataFrame]:
    family_summary = summary[mask].copy()

    # For the GSDMB-only figures we still define the variant family once, using
    # the same SNP-only and lower-frequency rules as script 08, then reuse that
    # exact ordered manifest across the paired carrier and genotype figures.
    manifest, plot_variant_keys = build_family_manifest(family_summary, family_name, family_note)
    plotted_summary = apply_variant_order(family_summary[family_summary['Variant_Key'].isin(plot_variant_keys)].copy(), manifest)
    plotted_genotypes = apply_variant_order(genotype_points[genotype_points['Variant_Key'].isin(plot_variant_keys)].copy(), manifest)
    return {
        'summary': plotted_summary,
        'genotypes': plotted_genotypes,
        'manifest': manifest,
    }


def export_workbook(common_family: dict[str, pd.DataFrame], rare_family: dict[str, pd.DataFrame], sample_size_lookup: dict[tuple[str, str], int]) -> None:
    sample_totals = pd.DataFrame([
        {'Cohort': cohort, 'Tissue': tissue, 'Panel_Sample_Count': count}
        for (cohort, tissue), count in sorted(sample_size_lookup.items())
    ])
    family_summary = pd.concat([
        common_family['manifest'].assign(Figure_Family='GSDMB common SNP / SNP-only'),
        rare_family['manifest'].assign(Figure_Family='GSDMB lower-frequency'),
    ], ignore_index=True)
    summary_df = pd.DataFrame([
        {'Output': 'Annotation source', 'Value': str(ANNOTATION_FILE), 'Notes': 'Variant labels and functional interpretation fields come from the stage-07 annotation workbook.'},
        {'Output': 'Genotype source', 'Value': str(FORCED_GENOTYPE_FILE), 'Notes': 'Explicit WT / Het ALT / Hom ALT states come from stage 05b forced genotypes.'},
        {'Output': 'Gene restriction', 'Value': 'Gene_Symbol contains GSDMB', 'Notes': 'This script is restricted to GSDMB-only variants before family splitting.'},
        {'Output': 'SNP definition', 'Value': f"gnomAD NFE AF > {THRESHOLDS['min_nfe_af']:.2f}", 'Notes': 'This is the explicit SNP-only family used in both GSDMB carrier and genotype figures.'},
        {'Output': 'Lower-frequency definition', 'Value': f"gnomAD NFE AF <= {THRESHOLDS['min_nfe_af']:.2f} or missing", 'Notes': 'Missing gnomAD NFE values remain in the lower-frequency family.'},
        {'Output': 'Variant-family ordering', 'Value': 'CHROM, POS, REF, ALT, Variant_Key', 'Notes': 'The same ordered manifest is reused across comparable GSDMB-only figures.'},
        {'Output': 'Callable exclusion rule', 'Value': 'Exclude only variants with zero callable samples across all panels in that family', 'Notes': 'These variants remain documented in the manifest with Included_In_Plots = False.'},
        {'Output': 'SNP-only filenames', 'Value': f'{output_snp_only_image.name}; {output_snp_only_genotype_image.name}', 'Notes': 'Alias files added so the SNP-only outputs are explicit for reporting.'},
    ])
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        summary_df.to_excel(writer, sheet_name='summary', index=False)
        sample_totals.to_excel(writer, sheet_name='sample_totals', index=False)
        common_family['summary'].to_excel(writer, sheet_name='common_snp_points', index=False)
        rare_family['summary'].to_excel(writer, sheet_name='lower_freq_points', index=False)
        common_family['genotypes'].to_excel(writer, sheet_name='common_snp_genotypes', index=False)
        rare_family['genotypes'].to_excel(writer, sheet_name='lower_freq_genotypes', index=False)
        common_family['manifest'].to_excel(writer, sheet_name='common_snp_manifest', index=False)
        rare_family['manifest'].to_excel(writer, sheet_name='lower_freq_manifest', index=False)
        family_summary.to_excel(writer, sheet_name='family_manifest_all', index=False)


def generate_gsdmb_maps() -> None:
    print('--- Generating split GSDMB-only landscape maps from explicit forced genotypes ---')
    validate_file_exists(ANNOTATION_FILE, 'Script 09 annotation workbook')
    validate_file_exists(FORCED_GENOTYPE_FILE, 'Script 09 forced genotype workbook')

    long_df = build_allele_genotype_long(ANNOTATION_FILE, FORCED_GENOTYPE_FILE)
    variant_summary, sample_size_lookup = prepare_variant_summary(long_df)
    genotype_points = build_genotype_points(variant_summary)

    common_family = select_family(
        variant_summary,
        genotype_points,
        family_name='Established_SNP',
        family_note=f'gnomAD NFE AF > {THRESHOLDS["min_nfe_af"]:.2f}',
        mask=variant_summary['Is_Common_SNP_By_NFE'],
    )
    rare_family = select_family(
        variant_summary,
        genotype_points,
        family_name='Lower_Frequency_Mutation',
        family_note=f'gnomAD NFE AF <= {THRESHOLDS["min_nfe_af"]:.2f} or missing',
        mask=~variant_summary['Is_Common_SNP_By_NFE'],
    )

    comparison_labels = load_comparison_significant_labels('GSDMB (script 9)')
    common_annotation_labels = set(select_high_effect_labels(common_family['summary'])) | comparison_labels
    rare_annotation_labels = set(select_standout_label(rare_family['summary'])) | comparison_labels
    common_variant_annotation_mask = common_family['summary']['Variant_Label'].isin(common_annotation_labels)
    rare_variant_annotation_mask = rare_family['summary']['Variant_Label'].isin(rare_annotation_labels)
    common_genotype_annotation_mask = (
        common_family['genotypes']['Variant_Label'].isin(common_annotation_labels)
        & common_family['genotypes']['Genotype_Category'].isin(['Het', 'Hom'])
    )
    rare_genotype_annotation_mask = (
        rare_family['genotypes']['Variant_Label'].isin(rare_annotation_labels)
        & rare_family['genotypes']['Genotype_Category'].isin(['Het', 'Hom'])
    )

    plot_variant_landscape(
        common_family['summary'],
        sample_size_lookup,
        'GSDMB landscape: SNP-only established variants',
        [output_common_image, output_snp_only_image],
        hue_col='Impact_Severity',
        palette=IMPACT_COLORS,
        legend_title='Impact',
        annotation_mask=common_variant_annotation_mask,
        impact_order=IMPACT_ORDER,
        right_pad=0.88,
        highlight_top_variant=False,
    )
    plot_variant_landscape(
        rare_family['summary'],
        sample_size_lookup,
        'GSDMB landscape: lower-frequency non-SNP variants',
        output_rare_image,
        hue_col='Impact_Severity',
        palette=IMPACT_COLORS,
        legend_title='Impact',
        annotation_mask=rare_variant_annotation_mask,
        impact_order=IMPACT_ORDER,
        right_pad=0.90,
        highlight_top_variant=False,
    )
    plot_genotype_landscape(
        common_family['genotypes'],
        sample_size_lookup,
        'GSDMB genotype landscape: SNP-only established variants (WT / Het / Hom)',
        [output_common_genotype_image, output_snp_only_genotype_image],
        genotype_colors=DEFAULT_GENOTYPE_COLORS,
        annotation_mask=common_genotype_annotation_mask,
    )
    plot_genotype_landscape(
        rare_family['genotypes'],
        sample_size_lookup,
        'GSDMB genotype landscape: lower-frequency non-SNP variants (WT / Het / Hom)',
        output_rare_genotype_image,
        genotype_colors=DEFAULT_GENOTYPE_COLORS,
        annotation_mask=rare_genotype_annotation_mask,
    )

    export_workbook(common_family, rare_family, sample_size_lookup)
    print(f'SUCCESS: split GSDMB workbook saved to: {output_excel}')


if __name__ == '__main__':
    generate_gsdmb_maps()
