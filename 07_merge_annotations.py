"""
Script 07: Consolidation of VEP-Annotated VCFs
==============================================

Purpose
-------
Read every sample-level VEP-annotated VCF, expand transcript-level consequence
records, and merge the resulting annotations into a single Excel workbook for
all downstream descriptive and inferential analyses.

Methodological note
-------------------
The resulting workbook is intentionally annotation-focused. It preserves the
transcript-expanded VEP rows plus prioritised annotation summaries, but it is
not a complete sample-by-locus genotype matrix. Explicit WT / Het / Hom /
callability logic lives in script 05b and must be sourced from the forced-
genotype workbook rather than inferred here.

Output
------
- ``Biological_Annotations``: transcript-expanded raw annotation rows with
  stable keys and derived annotation helper columns.
- ``Sample_Variant_Annotations``: one prioritised annotation row per
  sample/variant allele.
- ``Variant_Annotations``: one prioritised annotation row per unique variant
  allele across the cohort, suitable for downstream joins.
- ``README``: scope and sheet guide for the workbook.
"""

from __future__ import annotations

import glob
import gzip
import os
import subprocess
from pathlib import Path

import pandas as pd

from pipeline_utils import find_col, get_paths
from pipeline_validation import print_validation_summary, validate_required_columns
from variant_annotation_utils import (
    choose_representative_annotation_row,
    derive_annotation_columns,
    describe_annotation_basis,
)

PATHS = get_paths()
results_dir = str(PATHS["annotated_variants_dir"])
os.makedirs(results_dir, exist_ok=True)
output_file = os.path.join(results_dir, "GSDMB_Annotated_Variants.xlsx")
tfm_root = str(PATHS["project_root"])


def get_vep_headers(vcf_path: str) -> list[str] | None:
    """Extract the CSQ field names directly from one VCF header."""
    try:
        opener = gzip.open if vcf_path.endswith('.gz') else open
        with opener(vcf_path, 'rt') as handle:
            for line in handle:
                if line.startswith('##INFO=<ID=CSQ'):
                    header_part = line.split('Format: ')[1].split('"')[0]
                    return header_part.strip().split('|')
    except Exception as exc:  # pragma: no cover - defensive logging for malformed VCFs
        print(f"Debug: header extraction failed for {vcf_path}: {exc}")
    return None


def bcftools_query_vcf(vcf_path: str) -> list[dict[str, object]]:
    """Extract explicit VCF fields plus transcript-expanded CSQ records."""
    try:
        hdr = subprocess.check_output(
            ['bcftools', 'view', '-h', vcf_path], stderr=subprocess.DEVNULL
        ).decode()
    except subprocess.CalledProcessError:
        print(f"  WARNING: could not read header for {vcf_path}")
        return []

    has_gt = '##FORMAT=<ID=GT,' in hdr
    has_dp = '##FORMAT=<ID=DP,' in hdr
    has_af = '##FORMAT=<ID=AF,' in hdr
    has_ao = '##FORMAT=<ID=AO,' in hdr
    has_ro = '##FORMAT=<ID=RO,' in hdr
    has_iaf = '##INFO=<ID=AF,' in hdr
    has_idp = '##INFO=<ID=DP,' in hdr

    fmt_parts = ['%CHROM', '%POS', '%REF', '%ALT', '%QUAL', '%FILTER']
    if has_dp:
        fmt_parts.append('[%DP]')
    elif has_idp:
        fmt_parts.append('%INFO/DP')
    else:
        fmt_parts.append('.')

    ao_ro_mode = has_ao and has_ro and not has_af
    if has_af:
        fmt_parts.append('[%AF]')
    elif ao_ro_mode:
        fmt_parts.append('[%AO]')
        fmt_parts.append('[%RO]')
    elif has_iaf:
        fmt_parts.append('%INFO/AF')
    else:
        fmt_parts.append('.')

    fmt_parts.append('[%GT]' if has_gt else '.')
    fmt_parts.append('%INFO/CSQ')
    fmt_string = '\t'.join(fmt_parts) + '\n'

    try:
        raw = subprocess.check_output(
            ['bcftools', 'query', '-f', fmt_string, vcf_path],
            stderr=subprocess.DEVNULL,
        ).decode()
    except subprocess.CalledProcessError as exc:
        print(f"  WARNING: bcftools query failed for {vcf_path}: {exc}")
        return []

    base_cols = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER']
    meta_cols = base_cols + (['DP', 'AO', 'RO', 'GT', 'CSQ_RAW'] if ao_ro_mode else ['DP', 'AF', 'GT', 'CSQ_RAW'])
    vep_fields = get_vep_headers(vcf_path)

    rows: list[dict[str, object]] = []
    for line in raw.splitlines():
        if not line.strip():
            continue
        parts = line.split('	', len(meta_cols) - 1)
        if len(parts) < len(meta_cols):
            parts += ['.'] * (len(meta_cols) - len(parts))
        record = dict(zip(meta_cols, parts))

        if ao_ro_mode:
            try:
                ao = float(record.get('AO', 0))
                ro = float(record.get('RO', 0))
                record['AF'] = f"{ao / (ao + ro):.4f}" if (ao + ro) > 0 else '.'
            except (ValueError, TypeError):
                record['AF'] = '.'
            del record['AO'], record['RO']

        csq_raw = record.pop('CSQ_RAW', '.')
        if csq_raw in ('.', ''):
            if vep_fields:
                for field in vep_fields:
                    record[field] = '.'
            rows.append(record)
            continue

        for entry in csq_raw.split(','):
            row = record.copy()
            values = entry.split('|')
            if vep_fields:
                for i, field in enumerate(vep_fields):
                    row[field] = values[i] if i < len(values) else '.'
            else:
                for i, value in enumerate(values):
                    row[f'CSQ_{i}'] = value
            rows.append(row)
    return rows


def build_sample_variant_annotations(raw_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse transcript rows to one representative annotation per sample/variant allele."""
    rows: list[dict[str, object]] = []
    group_cols = ['Sample', 'Cohort', 'Tissue', 'Variant_Key']
    for _, group in raw_df.groupby(group_cols, sort=False):
        rep = choose_representative_annotation_row(group)
        rows.append({
            'Sample': rep.get('Sample', ''),
            'Cohort': rep.get('Cohort', ''),
            'Tissue': rep.get('Tissue', ''),
            'Variant_Key': rep.get('Variant_Key', ''),
            'Locus_Key': rep.get('Locus_Key', ''),
            'rsID': rep.get('rsID', ''),
            'CHROM': rep.get('CHROM', ''),
            'POS': rep.get('POS', None),
            'REF': rep.get('REF', ''),
            'ALT': rep.get('ALT', ''),
            'Existing_variation': rep.get('Existing_variation', ''),
            'Gene_Symbol': rep.get('SYMBOL', ''),
            'Gene_ID': rep.get('Gene', ''),
            'Detailed_Consequence': rep.get('Detailed_Consequence', ''),
            'Consequence_Raw': rep.get('Consequence_Raw', ''),
            'Location_Class': rep.get('Location_Class', ''),
            'Functional_Class': rep.get('Functional_Class', ''),
            'Impact_Severity': rep.get('Impact_Severity', ''),
            'gnomAD_NFE_AF_combined': rep.get('gnomAD_NFE_AF_combined', None),
            'gnomAD_NFE_Source': rep.get('gnomAD_NFE_Source', ''),
            'SIFT_Raw': rep.get('SIFT_Raw', ''),
            'SIFT_Prediction': rep.get('SIFT_Prediction', ''),
            'PolyPhen_Raw': rep.get('PolyPhen_Raw', ''),
            'PolyPhen_Prediction': rep.get('PolyPhen_Prediction', ''),
            'HGVSc': rep.get('HGVSc', ''),
            'HGVSp': rep.get('HGVSp', ''),
            'Feature': rep.get('Feature', ''),
            'Feature_type': rep.get('Feature_type', ''),
            'CANONICAL': rep.get('CANONICAL', ''),
            'MANE_SELECT': rep.get('MANE_SELECT', ''),
            'PICK': rep.get('PICK', ''),
            'Source_GT': rep.get('Source_GT', ''),
            'Source_GT_Class': rep.get('Source_GT_Class', ''),
            'Transcript_Rows_For_Sample_Variant': int(len(group)),
            'Representative_Annotation_Basis': describe_annotation_basis(rep),
            'Workbook_Scope': 'annotation_summary_per_sample_variant',
            'Genotype_Interpretation_Note': 'Source_GT / Source_GT_Class come from the originating VCF row only and do not form a complete genotype matrix.',
        })
    out = pd.DataFrame(rows)
    return out.sort_values(['Cohort', 'Tissue', 'Sample', 'CHROM', 'POS', 'ALT'], kind='stable').reset_index(drop=True)


def build_variant_annotations(raw_df: pd.DataFrame, sample_variant_df: pd.DataFrame) -> pd.DataFrame:
    """Build one representative annotation row per unique variant allele."""
    rows: list[dict[str, object]] = []
    for _, group in raw_df.groupby(['Variant_Key'], sort=False):
        rep = choose_representative_annotation_row(group)
        sample_variants = sample_variant_df[sample_variant_df['Variant_Key'] == rep.get('Variant_Key', '')]
        rows.append({
            'Variant_Key': rep.get('Variant_Key', ''),
            'Locus_Key': rep.get('Locus_Key', ''),
            'rsID': rep.get('rsID', ''),
            'CHROM': rep.get('CHROM', ''),
            'POS': rep.get('POS', None),
            'REF': rep.get('REF', ''),
            'ALT': rep.get('ALT', ''),
            'Existing_variation': rep.get('Existing_variation', ''),
            'Gene_Symbol': rep.get('SYMBOL', ''),
            'Gene_ID': rep.get('Gene', ''),
            'Detailed_Consequence': rep.get('Detailed_Consequence', ''),
            'Consequence_Raw': rep.get('Consequence_Raw', ''),
            'Location_Class': rep.get('Location_Class', ''),
            'Functional_Class': rep.get('Functional_Class', ''),
            'Impact_Severity': rep.get('Impact_Severity', ''),
            'gnomAD_NFE_AF_combined': rep.get('gnomAD_NFE_AF_combined', None),
            'gnomAD_NFE_Source': rep.get('gnomAD_NFE_Source', ''),
            'SIFT_Raw': rep.get('SIFT_Raw', ''),
            'SIFT_Prediction': rep.get('SIFT_Prediction', ''),
            'PolyPhen_Raw': rep.get('PolyPhen_Raw', ''),
            'PolyPhen_Prediction': rep.get('PolyPhen_Prediction', ''),
            'HGVSc': rep.get('HGVSc', ''),
            'HGVSp': rep.get('HGVSp', ''),
            'Feature': rep.get('Feature', ''),
            'Feature_type': rep.get('Feature_type', ''),
            'Representative_Transcript': rep.get('Feature', ''),
            'Representative_Annotation_Basis': describe_annotation_basis(rep),
            'Transcript_Rows_For_Variant': int(len(group)),
            'Samples_Annotated': int(sample_variants['Sample'].nunique()),
            'Cohorts_Observed': '; '.join(sorted(sample_variants['Cohort'].dropna().astype(str).unique())),
            'Tissues_Observed': '; '.join(sorted(sample_variants['Tissue'].dropna().astype(str).unique())),
            'Workbook_Scope': 'annotation_summary_per_variant_allele',
            'Join_Use': 'Join to stage 05b genotype summaries by Variant_Key; use Locus_Key only when working at locus level.',
        })
    out = pd.DataFrame(rows)
    return out.sort_values(['CHROM', 'POS', 'REF', 'ALT'], kind='stable').reset_index(drop=True)


def build_readme() -> pd.DataFrame:
    return pd.DataFrame([
        ['Workbook purpose', 'Annotation-focused workbook built from VEP-annotated VCF rows. It is useful for gene, transcript, consequence, impact, and predictor labels.'],
        ['Critical scope note', 'This workbook is not a complete sample-by-locus genotype matrix. WT / Het / Hom / callable denominator logic must come from stage 05b forced genotypes.'],
        ['Sheet: Biological_Annotations', 'Transcript-expanded VEP rows with stable join keys plus derived annotation helper columns.'],
        ['Sheet: Sample_Variant_Annotations', 'One representative annotation row per sample and variant allele. Keeps source GT as a raw VCF field only.'],
        ['Sheet: Variant_Annotations', 'One representative annotation row per unique variant allele, suitable for downstream joining and figure labels.'],
        ['Stable join keys', 'Variant_Key = CHROM:POS:REF:ALT and Locus_Key = CHROM:POS:REF. Variant_Key is the preferred bridge to allele-specific genotype summaries.'],
        ['Transcript prioritisation', 'Representative annotations prefer PICK, then CANONICAL, then MANE_SELECT, then worst IMPACT / worst consequence.'],
        ['Functional summaries', 'Detailed_Consequence, Functional_Class, Location_Class, Impact_Severity, SIFT_Prediction, and PolyPhen_Prediction are derived helper fields layered onto the VEP output.'],
    ], columns=['Field', 'Description'])


def merge_vep_to_excel() -> None:
    print('--- Starting VEP consolidation (annotation-focused workbook) ---')
    vcf_files = glob.glob(os.path.join(tfm_root, '**/*.vep.vcf.gz'), recursive=True)
    if not vcf_files:
        print('ERROR: No *.vep.vcf.gz files found under', tfm_root)
        return

    print(f'Found {len(vcf_files)} VEP-annotated VCF(s)')
    all_variants: list[pd.DataFrame] = []

    for vcf_path in sorted(vcf_files):
        parts = vcf_path.split(os.sep)
        sample_name = os.path.basename(vcf_path).replace('.vep.vcf.gz', '')
        cohort_type = 'Unknown'
        tissue_type = 'Unknown'
        try:
            sample_dir_idx = next(i for i, part in enumerate(parts) if part == sample_name)
            if parts[sample_dir_idx - 1].lower() == 'dna_calls':
                tissue_type = parts[sample_dir_idx - 2].capitalize()
                cohort_type = parts[sample_dir_idx - 3].capitalize()
            else:
                tissue_type = parts[sample_dir_idx - 1].capitalize()
                cohort_type = parts[sample_dir_idx - 2].capitalize()
        except (StopIteration, IndexError):
            pass

        print(f'  Processing: {cohort_type}/{tissue_type}/{sample_name}')
        rows = bcftools_query_vcf(vcf_path)
        if not rows:
            print(f'    WARNING: no variants extracted from {vcf_path}')
            continue

        df = pd.DataFrame(rows)
        df['Sample'] = sample_name
        df['Cohort'] = cohort_type
        df['Tissue'] = tissue_type
        all_variants.append(df)

    if not all_variants:
        print('No variants found in any VCF.')
        return

    final_df = pd.concat(all_variants, ignore_index=True)
    final_df = derive_annotation_columns(final_df)

    impact_col = find_col(final_df, 'IMPACT')
    if impact_col:
        impact_order = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
        final_df['_impact_rank'] = final_df[impact_col].astype(str).str.upper().map(impact_order).fillna(4)
        final_df = final_df.sort_values(['_impact_rank', 'CHROM', 'POS', 'REF', 'ALT', 'Sample'], kind='stable').drop(columns=['_impact_rank'])

    priority = [
        'Cohort', 'Tissue', 'Sample', 'Variant_Key', 'Locus_Key', 'rsID',
        'CHROM', 'POS', 'REF', 'ALT', 'GT', 'Source_GT_Class',
        'SYMBOL', 'Gene', 'Feature', 'Feature_type',
        'Consequence', 'Consequence_Raw', 'Detailed_Consequence', 'Location_Class', 'Functional_Class',
        'IMPACT', 'Impact_Severity', 'gnomAD_NFE_AF_combined', 'gnomAD_NFE_Source', 'SIFT_Raw', 'SIFT_Prediction', 'PolyPhen_Raw', 'PolyPhen_Prediction',
        'HGVSc', 'HGVSp', 'Existing_variation', 'Annotation_Row_Scope',
    ]
    priority = [col for col in priority if col in final_df.columns]
    final_df = final_df[priority + [col for col in final_df.columns if col not in priority]]

    validate_required_columns(final_df, ['Cohort', 'Tissue', 'Sample', 'CHROM', 'POS', 'REF', 'ALT', 'Variant_Key', 'Locus_Key'], 'Script 07 merged annotation output')
    print_validation_summary(final_df, 'Sample', 'Script 07 merged annotations', ['Cohort', 'Tissue'])

    sample_variant_df = build_sample_variant_annotations(final_df)
    variant_annotations_df = build_variant_annotations(final_df, sample_variant_df)
    readme_df = build_readme()

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        readme_df.to_excel(writer, sheet_name='README', index=False)
        final_df.to_excel(writer, sheet_name='Biological_Annotations', index=False)
        sample_variant_df.to_excel(writer, sheet_name='Sample_Variant_Annotations', index=False)
        variant_annotations_df.to_excel(writer, sheet_name='Variant_Annotations', index=False)

    print(f"SUCCESS: Consolidated {len(final_df)} transcript-expanded rows into {output_file}")
    print(f"         Sample-level prioritised rows: {len(sample_variant_df)}")
    print(f"         Unique variant-level annotation rows: {len(variant_annotations_df)}")


if __name__ == '__main__':
    merge_vep_to_excel()
