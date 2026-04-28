"""Script 09b: Normal-versus-Tumour Landscape Comparison
=========================================================

Purpose
-------
Formalise the follow-up review of the landscape plots from scripts 08 and 09.
The script compares normal and tumour carrier frequencies, identifies
normal-only and tumour-only variants, and checks whether the visual patterns
remain after restricting the comparison to higher-quality samples.

Methodological note
-------------------
Frequencies now use stage-05b explicit allele-level genotype states. Annotation
fields come from stage 07 only for biological labels and interpretation. The
comparison therefore uses callable denominators rather than inferring wild-type
status from annotation absence.
"""

from __future__ import annotations

from math import comb
import re
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from annotation_genotype_bridge import build_allele_genotype_long
from pipeline_utils import ensure_directory, get_paths
from pipeline_validation import validate_file_exists

PATHS = get_paths()
ANNOTATION_FILE = PATHS['annotated_report']
FORCED_GENOTYPE_FILE = PATHS.get('forced_genotypes_dir', PATHS['results_dir'] / '05b_forced_genotypes') / 'GSDMB_Forced_Genotypes_Union_Sites.xlsx'
QC_DIR = PATHS['qc_visualisation_dir']
OUTPUT_DIR = ensure_directory(PATHS.get('landscape_comparison_dir', PATHS['results_dir'] / '09b_landscape_comparison'))
WORKBOOK_PATH = OUTPUT_DIR / 'Landscape_Normal_vs_Tumour_with_QC_Significance_Genes_UPDATED.xlsx'
SUMMARY_CSV = OUTPUT_DIR / 'Landscape_Normal_vs_Tumour_Summary.csv'
GSDMB_TESTS_CSV = OUTPUT_DIR / 'GSDMB_Frequency_Tests.csv'
GLOBAL_TESTS_CSV = OUTPUT_DIR / 'Global_Frequency_Tests.csv'
SIGNIFICANCE_SUMMARY_CSV = OUTPUT_DIR / 'Landscape_Frequency_Significance_Summary.csv'

HEADER_FILL = PatternFill('solid', fgColor='1F4E78')
HEADER_FONT = Font(color='FFFFFF', bold=True)
SECTION_FILL = PatternFill('solid', fgColor='D9EAF7')


def normalize_sample_name(value: object) -> str:
    text = str(value).strip().lower()
    text = re.sub(r'[^a-z0-9]+', '_', text)
    return re.sub(r'_+', '_', text).strip('_')


def fisher_two_sided(normal_carriers: int, normal_noncarriers: int, tumour_carriers: int, tumour_noncarriers: int) -> float | None:
    if min(normal_noncarriers, tumour_noncarriers, normal_carriers, tumour_carriers) < 0:
        return None
    row1 = normal_carriers + normal_noncarriers
    row2 = tumour_carriers + tumour_noncarriers
    col1 = normal_carriers + tumour_carriers
    total = row1 + row2
    if row1 == 0 or row2 == 0 or total == 0:
        return None

    def hypergeom_prob(x_val: int) -> float:
        return comb(row1, x_val) * comb(row2, col1 - x_val) / comb(total, col1)

    min_x = max(0, col1 - row2)
    max_x = min(row1, col1)
    observed = hypergeom_prob(normal_carriers)
    p_value = 0.0
    for x_val in range(min_x, max_x + 1):
        prob = hypergeom_prob(x_val)
        if prob <= observed + 1e-12:
            p_value += prob
    return min(p_value, 1.0)


def benjamini_hochberg(p_values: list[float | None]) -> list[float | None]:
    valid = [(idx, p) for idx, p in enumerate(p_values) if p is not None and not pd.isna(p)]
    adjusted: list[float | None] = [None] * len(p_values)
    if not valid:
        return adjusted
    order = sorted(valid, key=lambda item: item[1])
    previous = 1.0
    n_tests = len(order)
    for rank in range(n_tests, 0, -1):
        idx, p_value = order[rank - 1]
        current = min(previous, p_value * n_tests / rank)
        previous = current
        adjusted[idx] = current
    return adjusted


def load_genotype_long() -> pd.DataFrame:
    validate_file_exists(ANNOTATION_FILE, 'Script 09b annotation workbook')
    validate_file_exists(FORCED_GENOTYPE_FILE, 'Script 09b forced genotype workbook')
    long_df = build_allele_genotype_long(ANNOTATION_FILE, FORCED_GENOTYPE_FILE)
    long_df = long_df[long_df['Tissue'].isin(['Healthy', 'Tumour'])].copy()
    long_df['sample_key'] = long_df['Sample'].map(normalize_sample_name)
    return long_df


def load_qc_audits() -> pd.DataFrame:
    audit_map = [
        ('Breast', 'Healthy', QC_DIR / 'CoverageAudit_Breast_Control.csv'),
        ('Breast', 'Tumour', QC_DIR / 'CoverageAudit_Breast_Tumour.csv'),
        ('Endometrium', 'Healthy', QC_DIR / 'CoverageAudit_Endometrium_Control.csv'),
        ('Endometrium', 'Tumour', QC_DIR / 'CoverageAudit_Endometrium_Tumour.csv'),
    ]
    frames: list[pd.DataFrame] = []
    for cohort, tissue, path in audit_map:
        validate_file_exists(path, f'Script 09b QC audit {cohort}/{tissue}')
        frame = pd.read_csv(path)
        sample_col = next(col for col in frame.columns if col.strip().lower() == 'sample')
        zero_col = next(col for col in frame.columns if col.strip().lower() == 'zeros')
        sub = frame[[sample_col, zero_col]].copy()
        sub.columns = ['sample_raw', 'zero_amplicons']
        sub['sample_key'] = sub['sample_raw'].map(normalize_sample_name)
        sub['Cohort'] = cohort
        sub['Tissue'] = tissue
        sub['high_qc'] = sub['zero_amplicons'].fillna(999).eq(0)
        frames.append(sub)
    return pd.concat(frames, ignore_index=True)


def prepare_scope_table(long_df: pd.DataFrame, scope: str) -> pd.DataFrame:
    if scope == 'GSDMB':
        scope_df = long_df[long_df['Gene_Symbol'].astype(str).str.contains('GSDMB', case=False, na=False)].copy()
    else:
        scope_df = long_df.copy()
    return scope_df.drop_duplicates(['Sample', 'Variant_Key'])


def _first_nonempty(series: pd.Series):
    for value in series:
        if pd.notna(value) and str(value).strip() not in {'', 'nan'}:
            return value
    return None


def build_frequency_tests(scope_df: pd.DataFrame) -> pd.DataFrame:
    metadata_cols = [
        'Variant_Key', 'Locus_Key', 'rsID', 'Gene_Symbol', 'Gene_ID', 'CHROM', 'POS', 'REF', 'ALT',
        'Impact_Severity', 'Detailed_Consequence', 'Functional_Class', 'Location_Class',
        'Feature_type', 'Representative_Transcript', 'Representative_Annotation_Basis',
        'gnomAD_NFE_AF_combined', 'gnomAD_NFE_Source', 'Multi_Allelic_Union_Site', 'N_Alt_Alleles',
    ]
    rows: list[dict[str, object]] = []
    grouped = scope_df.groupby(['Cohort', 'Variant_Key'], sort=False)
    pooled_normal = scope_df[scope_df['Tissue'] == 'Healthy'].copy()
    for (cohort, variant_key), group in grouped:
        normal = pooled_normal[pooled_normal['Variant_Key'] == variant_key]
        tumour = group[group['Tissue'] == 'Tumour']
        normal_callable = normal[normal['Callable_Allele_Genotype']]
        tumour_callable = tumour[tumour['Callable_Allele_Genotype']]
        normal_carriers = int(normal_callable['Variant_Is_Carrier'].sum())
        tumour_carriers = int(tumour_callable['Variant_Is_Carrier'].sum())
        normal_n = int(normal_callable['Sample'].nunique())
        tumour_n = int(tumour_callable['Sample'].nunique())
        row = {
            'Cohort': cohort,
            'Comparison_Label': f'{cohort} tumour vs pooled normal',
            'Comparison_Logic': f'{cohort} tumour compared with pooled normal controls (breast normal + endometrium normal)',
            'Normal_Label': 'Pooled normal controls',
            'Tumour_Label': f'{cohort} tumour',
            'Normal_Source': 'Breast normal + endometrium normal',
            'Tumour_Source': f'{cohort} tumour',
            'Variant_Key': variant_key,
            'Normal_Carriers': normal_carriers,
            'Tumour_Carriers': tumour_carriers,
            'Normal_N': normal_n,
            'Tumour_N': tumour_n,
            'Normal_NonCallable': int((~normal['Callable_Allele_Genotype']).sum()),
            'Tumour_NonCallable': int((~tumour['Callable_Allele_Genotype']).sum()),
            'Normal_Frequency_pct': (normal_carriers / normal_n * 100) if normal_n else None,
            'Tumour_Frequency_pct': (tumour_carriers / tumour_n * 100) if tumour_n else None,
        }
        row['Tumour_minus_Normal_pct'] = None if row['Normal_Frequency_pct'] is None or row['Tumour_Frequency_pct'] is None else row['Tumour_Frequency_pct'] - row['Normal_Frequency_pct']
        row['Variant_Status'] = 'shared' if normal_carriers > 0 and tumour_carriers > 0 else ('tumour_only' if tumour_carriers > 0 else 'normal_only')
        row['Direction'] = 'Not estimable' if row['Tumour_minus_Normal_pct'] is None else ('Lower in tumour' if row['Tumour_minus_Normal_pct'] < -1e-9 else ('Higher in tumour' if row['Tumour_minus_Normal_pct'] > 1e-9 else 'Same frequency'))
        row['p_value'] = fisher_two_sided(normal_carriers, normal_n - normal_carriers, tumour_carriers, tumour_n - tumour_carriers)
        for col in metadata_cols:
            row[col] = _first_nonempty(group[col]) if col in group.columns else None
        rows.append(row)

    tests = pd.DataFrame(rows)
    tests['FDR_BH'] = None
    for _, indices in tests.groupby('Cohort').groups.items():
        idxs = list(indices)
        tests.loc[idxs, 'FDR_BH'] = benjamini_hochberg(tests.loc[idxs, 'p_value'].tolist())
    tests['Raw_p_lt_0.05'] = tests['p_value'].apply(lambda value: False if value is None or pd.isna(value) else value < 0.05)
    tests['FDR_lt_0.05'] = tests['FDR_BH'].apply(lambda value: False if value is None or pd.isna(value) else value < 0.05)
    preferred_order = [
        'Cohort', 'Comparison_Label', 'Comparison_Logic', 'Normal_Label', 'Tumour_Label', 'Normal_Source', 'Tumour_Source',
        'Gene_Symbol', 'Gene_ID', 'Variant_Key', 'Locus_Key', 'CHROM', 'POS', 'REF', 'ALT', 'rsID',
        'Impact_Severity', 'Detailed_Consequence', 'Functional_Class', 'Location_Class', 'Feature_type',
        'Normal_Carriers', 'Tumour_Carriers', 'Normal_N', 'Tumour_N', 'Normal_NonCallable', 'Tumour_NonCallable',
        'Normal_Frequency_pct', 'Tumour_Frequency_pct', 'Tumour_minus_Normal_pct', 'Variant_Status', 'Direction',
        'p_value', 'FDR_BH', 'Raw_p_lt_0.05', 'FDR_lt_0.05', 'Representative_Transcript',
        'Representative_Annotation_Basis', 'gnomAD_NFE_AF_combined', 'gnomAD_NFE_Source',
        'Multi_Allelic_Union_Site', 'N_Alt_Alleles',
    ]
    return tests[preferred_order]


def build_variant_status_summary(tests_df: pd.DataFrame) -> pd.DataFrame:
    return (
        tests_df.groupby(['Cohort', 'Impact_Severity', 'Variant_Status'])
        .size()
        .reset_index(name='n_variants')
        .sort_values(['Cohort', 'Impact_Severity', 'Variant_Status'])
    )


def build_direction_summary(tests_df: pd.DataFrame, direction_col: str, out_col: str) -> pd.DataFrame:
    shared = tests_df[tests_df['Variant_Status'] == 'shared'].copy()
    return (
        shared.groupby(['Cohort', 'Impact_Severity', direction_col])
        .size()
        .reset_index(name=out_col)
        .sort_values(['Cohort', 'Impact_Severity', direction_col])
    )


def build_qc_sensitivity(scope_df: pd.DataFrame, qc_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    aligned = scope_df.merge(
        qc_df[['sample_key', 'Cohort', 'Tissue', 'zero_amplicons', 'high_qc']],
        on=['sample_key', 'Cohort', 'Tissue'],
        how='inner',
    )
    high_qc = aligned[aligned['high_qc']].copy()
    qc_counts = (
        aligned.groupby(['Cohort', 'Tissue'])
        .apply(
            lambda group: pd.Series({
                'aligned_samples': group['sample_key'].nunique(),
                'high_qc_samples': group.loc[group['high_qc'], 'sample_key'].nunique(),
                'mean_zero_amplicons': group[['sample_key', 'zero_amplicons']].drop_duplicates()['zero_amplicons'].mean(),
            }),
            include_groups=False,
        )
        .reset_index()
    )
    if high_qc.empty:
        return pd.DataFrame(columns=['Cohort', 'Variant_Key']), qc_counts

    pooled_normal = high_qc[high_qc['Tissue'] == 'Healthy'].copy()
    rows: list[dict[str, object]] = []
    for (cohort, variant_key), group in high_qc.groupby(['Cohort', 'Variant_Key'], sort=False):
        normal = pooled_normal[pooled_normal['Variant_Key'] == variant_key]
        tumour = group[group['Tissue'] == 'Tumour']
        normal_callable = normal[normal['Callable_Allele_Genotype']]
        tumour_callable = tumour[tumour['Callable_Allele_Genotype']]
        normal_n = int(normal_callable['Sample'].nunique())
        tumour_n = int(tumour_callable['Sample'].nunique())
        normal_carriers = int(normal_callable['Variant_Is_Carrier'].sum())
        tumour_carriers = int(tumour_callable['Variant_Is_Carrier'].sum())
        normal_pct = (normal_carriers / normal_n * 100) if normal_n else None
        tumour_pct = (tumour_carriers / tumour_n * 100) if tumour_n else None
        delta = None if normal_pct is None or tumour_pct is None else tumour_pct - normal_pct
        rows.append({
            'Cohort': cohort,
            'Variant_Key': variant_key,
            'Normal_N_HighQC': normal_n,
            'Tumour_N_HighQC': tumour_n,
            'Normal_Carriers_HighQC': normal_carriers,
            'Tumour_Carriers_HighQC': tumour_carriers,
            'Normal_Frequency_pct_HighQC': normal_pct,
            'Tumour_Frequency_pct_HighQC': tumour_pct,
            'Tumour_minus_Normal_pct_HighQC': delta,
            'Direction_HighQC': 'Not estimable' if delta is None else ('Lower in tumour' if delta < -1e-9 else ('Higher in tumour' if delta > 1e-9 else 'Same frequency')),
        })
    return pd.DataFrame(rows), qc_counts


def build_summary(scope_label: str, status_df: pd.DataFrame, tests_df: pd.DataFrame, qc_counts_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort in sorted(tests_df['Cohort'].dropna().unique()):
        cohort_status = status_df[status_df['Cohort'] == cohort].copy()
        cohort_tests = tests_df[tests_df['Cohort'] == cohort].copy()
        pivot = cohort_status.pivot_table(index='Impact_Severity', columns='Variant_Status', values='n_variants', fill_value=0)
        shared = cohort_tests[cohort_tests['Variant_Status'] == 'shared']
        lower_count = int((shared['Direction'] == 'Lower in tumour').sum())
        higher_count = int((shared['Direction'] == 'Higher in tumour').sum())
        raw_hits = cohort_tests[cohort_tests['Raw_p_lt_0.05']].sort_values('p_value')
        hit_text = ', '.join(
            f"{(row.rsID or row.Variant_Key)} (p={row.p_value:.4g})"
            for row in raw_hits.head(5).itertuples(index=False)
            if row.p_value is not None and not pd.isna(row.p_value)
        ) or 'None'
        qc_rows = qc_counts_df[qc_counts_df['Cohort'] == cohort]
        qc_text = '; '.join(
            f"{row.Tissue}: aligned n={int(row.aligned_samples)}, high-QC n={int(row.high_qc_samples)}"
            for row in qc_rows.itertuples(index=False)
        )
        best_row = cohort_tests.sort_values('p_value', na_position='last').iloc[0]
        rows.append({
            'Map': scope_label,
            'Cohort': cohort,
            'Tumour_only_variants': int((cohort_tests['Variant_Status'] == 'tumour_only').sum()),
            'Normal_only_variants': int((cohort_tests['Variant_Status'] == 'normal_only').sum()),
            'Shared_variants': int((cohort_tests['Variant_Status'] == 'shared').sum()),
            'Shared_lower_in_tumour': lower_count,
            'Shared_higher_in_tumour': higher_count,
            'Raw_p_lt_0.05': int(cohort_tests['Raw_p_lt_0.05'].sum()),
            'FDR_lt_0.05': int(cohort_tests['FDR_lt_0.05'].sum()),
            'Best_rsID': best_row['rsID'] if pd.notna(best_row['rsID']) else best_row['Variant_Key'],
            'Best_variant_status': best_row['Variant_Status'],
            'Best_direction': best_row['Direction'],
            'Best_p_value': best_row['p_value'],
            'Best_FDR_BH': best_row['FDR_BH'],
            'Example_nominal_hits': hit_text,
            'Tumour_only_breakdown': '; '.join(
                f"{impact}: {int(pivot.loc[impact, 'tumour_only'])}"
                for impact in pivot.index
                if 'tumour_only' in pivot.columns and int(pivot.loc[impact, 'tumour_only']) > 0
            ) or 'None',
            'Normal_only_breakdown': '; '.join(
                f"{impact}: {int(pivot.loc[impact, 'normal_only'])}"
                for impact in pivot.index
                if 'normal_only' in pivot.columns and int(pivot.loc[impact, 'normal_only']) > 0
            ) or 'None',
            'QC_alignment': qc_text or 'Not available',
        })
    return pd.DataFrame(rows)


def style_workbook(path: Path) -> None:
    workbook = load_workbook(path)
    for sheet in workbook.worksheets:
        if sheet.title == 'Guide':
            sheet.column_dimensions['A'].width = 22
            sheet.column_dimensions['B'].width = 110
            for row in sheet.iter_rows(min_row=1, max_row=sheet.max_row, min_col=1, max_col=2):
                row[0].font = Font(bold=True)
                row[0].fill = SECTION_FILL
                row[0].alignment = Alignment(vertical='top', wrap_text=True)
                row[1].alignment = Alignment(vertical='top', wrap_text=True)
            continue
        sheet.freeze_panes = 'A2'
        for cell in sheet[1]:
            cell.fill = HEADER_FILL
            cell.font = HEADER_FONT
            cell.alignment = Alignment(wrap_text=True, horizontal='center', vertical='center')
        for col_cells in sheet.columns:
            width = max(len(str(cell.value)) if cell.value is not None else 0 for cell in col_cells) + 2
            sheet.column_dimensions[get_column_letter(col_cells[0].column)].width = min(max(width, 12), 38)
        for row in sheet.iter_rows(min_row=2):
            for cell in row:
                cell.alignment = Alignment(wrap_text=True, vertical='top')
    workbook.save(path)


def main() -> None:
    print('======================================================================')
    print('SCRIPT 09B: NORMAL-VERSUS-TUMOUR LANDSCAPE COMPARISON')
    print('======================================================================')
    print(f'Annotation workbook: {ANNOTATION_FILE}')
    print(f'Forced genotype workbook: {FORCED_GENOTYPE_FILE}')
    print(f'QC audit directory: {QC_DIR}')
    print(f'Output directory: {OUTPUT_DIR}')
    print('======================================================================')

    long_df = load_genotype_long()
    qc_audits = load_qc_audits()

    gsdmb_scope = prepare_scope_table(long_df, 'GSDMB')
    global_scope = prepare_scope_table(long_df, 'GLOBAL')

    gsdmb_tests = build_frequency_tests(gsdmb_scope)
    global_tests = build_frequency_tests(global_scope)

    gsdmb_status = build_variant_status_summary(gsdmb_tests)
    global_status = build_variant_status_summary(global_tests)
    gsdmb_direction = build_direction_summary(gsdmb_tests, 'Direction', 'n_shared_variants')
    global_direction = build_direction_summary(global_tests, 'Direction', 'n_shared_variants')

    gsdmb_high_qc, gsdmb_qc_counts = build_qc_sensitivity(gsdmb_scope, qc_audits)
    global_high_qc, global_qc_counts = build_qc_sensitivity(global_scope, qc_audits)
    gsdmb_tests = gsdmb_tests.merge(gsdmb_high_qc, on=['Cohort', 'Variant_Key'], how='left')
    global_tests = global_tests.merge(global_high_qc, on=['Cohort', 'Variant_Key'], how='left')

    gsdmb_qc_source = gsdmb_tests[gsdmb_tests['Direction_HighQC'].notna()].copy()
    gsdmb_qc_source['Direction_QC'] = gsdmb_qc_source['Direction_HighQC']
    global_qc_source = global_tests[global_tests['Direction_HighQC'].notna()].copy()
    global_qc_source['Direction_QC'] = global_qc_source['Direction_HighQC']
    gsdmb_qc_direction = build_direction_summary(gsdmb_qc_source, 'Direction_QC', 'n_shared_variants_high_qc')
    global_qc_direction = build_direction_summary(global_qc_source, 'Direction_QC', 'n_shared_variants_high_qc')

    summary = pd.concat([
        build_summary('GSDMB (script 9)', gsdmb_status, gsdmb_tests, gsdmb_qc_counts),
        build_summary('Global (script 8)', global_status, global_tests, global_qc_counts),
    ], ignore_index=True)

    def label_significant_row(row: pd.Series) -> str:
        rs_part = row['rsID'] if pd.notna(row.get('rsID')) else row['Variant_Key']
        gene_part = f"{row['Gene_Symbol']} " if pd.notna(row.get('Gene_Symbol')) and str(row.get('Gene_Symbol')).strip() else ''
        comparison_part = f"{row['Comparison_Label']}: " if pd.notna(row.get('Comparison_Label')) and str(row.get('Comparison_Label')).strip() else ''
        p_val = row['p_value']
        fdr_val = row['FDR_BH']
        return f"{comparison_part}{gene_part}{rs_part} (p={p_val:.4g}, FDR={fdr_val:.4g})" if p_val is not None and fdr_val is not None and not pd.isna(p_val) and not pd.isna(fdr_val) else f"{comparison_part}{gene_part}{rs_part}"

    significance_counts_rows: list[dict[str, object]] = []
    significance_detail_rows: list[dict[str, object]] = []
    detail_columns = [
        'Map', 'Cohort', 'Comparison_Label', 'Comparison_Logic', 'Normal_Label', 'Tumour_Label', 'Normal_Source', 'Tumour_Source',
        'Variants_tested', 'Raw_p_lt_0.05', 'FDR_lt_0.05',
        'Significance_Level', 'Also_FDR_significant', 'Gene_Symbol', 'Gene_ID', 'rsID', 'Variant_Key',
        'POS', 'Impact_Severity', 'Detailed_Consequence', 'Functional_Class', 'Location_Class',
        'Variant_Status', 'Direction', 'Normal_Carriers', 'Tumour_Carriers', 'Normal_N', 'Tumour_N',
        'Normal_NonCallable', 'Tumour_NonCallable', 'Normal_Frequency_pct', 'Tumour_Frequency_pct',
        'Tumour_minus_Normal_pct', 'p_value', 'FDR_BH',
    ]

    for map_label, frame in [('GSDMB (script 9)', gsdmb_tests), ('Global (script 8)', global_tests)]:
        for cohort, sub in frame.groupby('Cohort'):
            nominal_hits = sub[sub['Raw_p_lt_0.05']].sort_values(['p_value', 'FDR_BH', 'POS'], na_position='last')
            fdr_hits = sub[sub['FDR_lt_0.05']].sort_values(['FDR_BH', 'p_value', 'POS'], na_position='last')
            counts = {
                'Map': map_label,
                'Cohort': cohort,
                'Comparison_Label': f'{cohort} tumour vs pooled normal',
                'Comparison_Logic': f'{cohort} tumour compared with pooled normal controls (breast normal + endometrium normal)',
                'Normal_Label': 'Pooled normal controls',
                'Tumour_Label': f'{cohort} tumour',
                'Normal_Source': 'Breast normal + endometrium normal',
                'Tumour_Source': f'{cohort} tumour',
                'Variants_tested': len(sub),
                'Raw_p_lt_0.05': int(sub['Raw_p_lt_0.05'].sum()),
                'FDR_lt_0.05': int(sub['FDR_lt_0.05'].sum()),
                'Shared_variants': int((sub['Variant_Status'] == 'shared').sum()),
                'Tumour_only_variants': int((sub['Variant_Status'] == 'tumour_only').sum()),
                'Normal_only_variants': int((sub['Variant_Status'] == 'normal_only').sum()),
                'Nominal_significant_variants': '; '.join(label_significant_row(row) for _, row in nominal_hits.iterrows()) or 'None',
                'FDR_significant_variants': '; '.join(label_significant_row(row) for _, row in fdr_hits.iterrows()) or 'None',
            }
            significance_counts_rows.append(counts)
            if nominal_hits.empty:
                significance_detail_rows.append({
                    **{key: counts[key] for key in ['Map', 'Cohort', 'Comparison_Label', 'Variants_tested', 'Raw_p_lt_0.05', 'FDR_lt_0.05']},
                    'Comparison_Logic': counts['Comparison_Logic'],
                    'Normal_Label': counts['Normal_Label'],
                    'Tumour_Label': counts['Tumour_Label'],
                    'Normal_Source': counts['Normal_Source'],
                    'Tumour_Source': counts['Tumour_Source'],
                    'Significance_Level': 'None',
                    'Also_FDR_significant': False,
                    'Gene_Symbol': None,
                    'Gene_ID': None,
                    'rsID': None,
                    'Variant_Key': None,
                    'POS': None,
                    'Impact_Severity': None,
                    'Detailed_Consequence': None,
                    'Functional_Class': None,
                    'Location_Class': None,
                    'Variant_Status': None,
                    'Direction': None,
                    'Normal_Carriers': None,
                    'Tumour_Carriers': None,
                    'Normal_N': None,
                    'Tumour_N': None,
                    'Normal_NonCallable': None,
                    'Tumour_NonCallable': None,
                    'Normal_Frequency_pct': None,
                    'Tumour_Frequency_pct': None,
                    'Tumour_minus_Normal_pct': None,
                    'p_value': None,
                    'FDR_BH': None,
                })
                continue
            for _, row in nominal_hits.iterrows():
                detail = {
                    **{key: counts[key] for key in ['Map', 'Cohort', 'Comparison_Label', 'Variants_tested', 'Raw_p_lt_0.05', 'FDR_lt_0.05']},
                    'Comparison_Logic': counts['Comparison_Logic'],
                    'Normal_Label': counts['Normal_Label'],
                    'Tumour_Label': counts['Tumour_Label'],
                    'Normal_Source': counts['Normal_Source'],
                    'Tumour_Source': counts['Tumour_Source'],
                    'Significance_Level': 'Nominal (p < 0.05)',
                    'Also_FDR_significant': bool(row['FDR_lt_0.05']),
                }
                for col in detail_columns:
                    if col in detail:
                        continue
                    detail[col] = row[col] if col in row.index else None
                significance_detail_rows.append(detail)

    significance_counts = pd.DataFrame(significance_counts_rows)
    significance_summary = pd.DataFrame(significance_detail_rows)[detail_columns]
    qc_counts = pd.concat([
        gsdmb_qc_counts.assign(Map='GSDMB (script 9)'),
        global_qc_counts.assign(Map='Global (script 8)'),
    ], ignore_index=True)

    guide = pd.DataFrame([
        ['Purpose', 'Summary of the visual patterns seen in scripts 8 and 9, with supporting variant counts, frequency tests, QC sensitivity checks, and a detailed significant-variant overview.'],
        ['How to read frequencies', 'Breast tumour and endometrium tumour are each compared against pooled normal controls (breast normal + endometrium normal). Frequencies are carrier frequencies among callable genotype states only.'],
        ['Annotation source', 'Biological labels such as gene, consequence, functional class, and impact come from the stage-07 annotation workbook.'],
        ['Genotype source', 'Carrier and non-carrier states come from the stage-05b forced-genotype matrix, not from annotation-derived subtraction.'],
        ['Variant status', 'shared = present in both tissues; tumour_only = seen only in tumour; normal_only = seen only in normal.'],
        ['Direction', 'Tumour_minus_Normal_pct < 0 means the variant is less frequent in tumour than in normal, among callable samples.'],
        ['QC sensitivity', 'High-QC comparisons use only aligned samples with zero_amplicons = 0 in the coverage audit generated by script 03.'],
        ['Frequency significance', 'p_value is a two-sided Fisher exact test on carrier vs callable non-carrier counts with pooled normal controls as the baseline. FDR_BH is Benjamini-Hochberg correction within each cohort and map.'],
    ])

    gsdmb_key = gsdmb_tests.sort_values(['Cohort', 'p_value', 'Tumour_minus_Normal_pct'], na_position='last').groupby('Cohort', group_keys=False).head(12)
    global_key = global_tests.sort_values(['Cohort', 'p_value', 'Tumour_minus_Normal_pct'], na_position='last').groupby('Cohort', group_keys=False).head(12)

    with pd.ExcelWriter(WORKBOOK_PATH, engine='openpyxl') as writer:
        guide.to_excel(writer, sheet_name='Guide', index=False, header=False)
        summary.to_excel(writer, sheet_name='Summary', index=False)
        significance_summary.to_excel(writer, sheet_name='Significance_Overview', index=False)
        significance_counts.to_excel(writer, sheet_name='Significance_Counts', index=False)
        qc_counts.to_excel(writer, sheet_name='QC_Sample_Counts', index=False)
        gsdmb_status.to_excel(writer, sheet_name='GSDMB_Variant_Status', index=False)
        gsdmb_direction.to_excel(writer, sheet_name='GSDMB_Shared_Direction', index=False)
        gsdmb_key.to_excel(writer, sheet_name='GSDMB_Key_with_p_FDR', index=False)
        gsdmb_tests.to_excel(writer, sheet_name='GSDMB_All_with_p_FDR', index=False)
        gsdmb_qc_direction.to_excel(writer, sheet_name='GSDMB_QC_Sensitivity', index=False)
        global_status.to_excel(writer, sheet_name='Global_Variant_Status', index=False)
        global_direction.to_excel(writer, sheet_name='Global_Shared_Direction', index=False)
        global_key.to_excel(writer, sheet_name='Global_Key_with_p_FDR', index=False)
        global_tests.to_excel(writer, sheet_name='Global_All_with_p_FDR', index=False)
        global_qc_direction.to_excel(writer, sheet_name='Global_QC_Sensitivity', index=False)

    style_workbook(WORKBOOK_PATH)
    summary.to_csv(SUMMARY_CSV, index=False)
    gsdmb_tests.to_csv(GSDMB_TESTS_CSV, index=False)
    global_tests.to_csv(GLOBAL_TESTS_CSV, index=False)
    significance_summary.to_csv(SIGNIFICANCE_SUMMARY_CSV, index=False)

    print('Saved workbook:', WORKBOOK_PATH)
    print('Saved summary CSV:', SUMMARY_CSV)
    print('Saved GSDMB tests CSV:', GSDMB_TESTS_CSV)
    print('Saved global tests CSV:', GLOBAL_TESTS_CSV)
    print('Saved significance summary CSV:', SIGNIFICANCE_SUMMARY_CSV)


if __name__ == '__main__':
    main()
