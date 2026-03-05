"""
Reads all VEP-annotated VCFs produced by script 06 and consolidates them into
a single Excel workbook sheet ('Biological_Annotations'), with one row per
variant-transcript combination.

VCF parsing is done via bcftools query rather than by reading the VCF in Python.
This avoids the pipe-delimiter collision in the CSQ field: VEP encodes all
annotation sub-fields as a single pipe-delimited string inside the INFO column,
which cannot be split with a simple sed or string-split on the raw VCF text
without corrupting multi-allelic or multi-transcript records.

Output columns are sorted so that clinically relevant fields (gene symbol,
consequence, impact, HGVSp, allele frequency) appear first.
"""

import pandas as pd
import glob
import os
import gzip
import subprocess


# --- Configuration -----------------------------------------------------------
results_dir = "/home/gadeaalonsoj/tfm/gsdmb_final_results"
output_file = os.path.join(results_dir, "GSDMB_Annotated_Report_Fixed.xlsx")
tfm_root    = "/home/gadeaalonsoj/tfm"


def get_vep_headers(vcf_path: str):
    """
    Parse the VCF header to extract the ordered list of sub-field names embedded
    in the CSQ INFO tag. VEP writes these as a pipe-delimited 'Format:' string
    in the ##INFO=<ID=CSQ,...> header line. Knowing the field order is necessary
    to correctly unpack each CSQ entry during variant parsing.
    Returns a list of field name strings, or None if the header line is absent.
    """
    try:
        opener = gzip.open if vcf_path.endswith('.gz') else open
        with opener(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('##INFO=<ID=CSQ'):
                    header_part = line.split('Format: ')[1].split('"')[0]
                    return header_part.strip().split('|')
    except Exception as e:
        print(f"  Warning: could not extract CSQ header from {vcf_path}: {e}")
    return None


def find_col(df: pd.DataFrame, target_name: str):
    """
    Case-insensitive column lookup. VEP field names vary in capitalisation
    between cache versions, so direct string matching is unreliable.
    """
    for col in df.columns:
        if col.upper() == target_name.upper():
            return col
    return None


def bcftools_query_vcf(vcf_path: str) -> list:
    """
    Extract variant data from a VEP-annotated VCF using bcftools query and
    return a list of dicts, one per variant-transcript row.

    The FORMAT fields (GT, DP, AF) are extracted as explicit, separate columns
    rather than relying on the CSQ string. This matters because CSQ is a single
    pipe-delimited blob; pulling FORMAT values out of it via string splitting is
    fragile and error-prone. bcftools query handles multi-sample FORMAT columns
    and missing tags gracefully with the [...] per-sample syntax.

    Allele frequency (AF) extraction logic, in priority order:
      1. FORMAT/AF  — direct per-sample AF if the caller emits it
      2. AO / RO    — alt and ref observation counts; AF derived as AO/(AO+RO)
                      (TVC emits AO/RO rather than AF in some configurations)
      3. INFO/AF    — site-level AF as a fallback
      4. '.'        — placeholder if none of the above are present
    """
    try:
        hdr = subprocess.check_output(
            ['bcftools', 'view', '-h', vcf_path], stderr=subprocess.DEVNULL
        ).decode()
    except subprocess.CalledProcessError:
        print(f"  Warning: could not read VCF header for {vcf_path}")
        return []

    # Inspect the header to determine which FORMAT/INFO tags are available
    has_gt  = '##FORMAT=<ID=GT,' in hdr
    has_dp  = '##FORMAT=<ID=DP,' in hdr
    has_af  = '##FORMAT=<ID=AF,' in hdr
    has_ao  = '##FORMAT=<ID=AO,' in hdr
    has_ro  = '##FORMAT=<ID=RO,' in hdr
    has_iaf = '##INFO=<ID=AF,'   in hdr
    has_idp = '##INFO=<ID=DP,'   in hdr

    # Build the bcftools query format string dynamically based on available tags.
    # CSQ is kept as a single raw field and split later — if we asked bcftools to
    # split it, the pipe characters would be ambiguous with column delimiters.
    fmt_parts = ['%CHROM', '%POS', '%REF', '%ALT', '%QUAL', '%FILTER']

    # Depth: FORMAT/DP is per-sample; INFO/DP is the pooled total across samples
    if has_dp:
        fmt_parts.append('[%DP]')
    elif has_idp:
        fmt_parts.append('%INFO/DP')
    else:
        fmt_parts.append('.')

    # Allele fraction
    ao_ro_mode = has_ao and has_ro and not has_af   # TVC-style: no direct AF field
    if has_af:
        fmt_parts.append('[%AF]')
    elif ao_ro_mode:
        fmt_parts.append('[%AO]')
        fmt_parts.append('[%RO]')
    elif has_iaf:
        fmt_parts.append('%INFO/AF')
    else:
        fmt_parts.append('.')

    if has_gt:
        fmt_parts.append('[%GT]')
    else:
        fmt_parts.append('.')

    fmt_parts.append('%INFO/CSQ')   # must be last; contains internal pipe delimiters

    fmt_string = '\t'.join(fmt_parts) + '\n'

    try:
        raw = subprocess.check_output(
            ['bcftools', 'query', '-f', fmt_string, vcf_path],
            stderr=subprocess.DEVNULL
        ).decode()
    except subprocess.CalledProcessError as e:
        print(f"  Warning: bcftools query failed for {vcf_path}: {e}")
        return []

    # Column names matching the format string built above
    base_cols = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER']
    if ao_ro_mode:
        meta_cols = base_cols + ['DP', 'AO', 'RO', 'GT', 'CSQ_RAW']
    else:
        meta_cols = base_cols + ['DP', 'AF', 'GT', 'CSQ_RAW']

    vep_fields = get_vep_headers(vcf_path)

    rows = []
    for line in raw.splitlines():
        if not line.strip():
            continue

        # Split on tab only up to the expected number of columns; this preserves
        # any tabs that might appear inside the CSQ field itself
        parts = line.split('\t', len(meta_cols) - 1)
        if len(parts) < len(meta_cols):
            parts += ['.'] * (len(meta_cols) - len(parts))   # pad missing trailing fields
        record = dict(zip(meta_cols, parts))

        # Compute AF from AO and RO if in TVC-style mode
        if ao_ro_mode:
            try:
                ao = float(record.get('AO', 0))
                ro = float(record.get('RO', 0))
                record['AF'] = f"{ao / (ao + ro):.4f}" if (ao + ro) > 0 else '.'
            except (ValueError, TypeError):
                record['AF'] = '.'
            del record['AO'], record['RO']

        # Expand the CSQ field: VEP writes one comma-separated entry per transcript,
        # and each transcript entry is a pipe-delimited list of annotation values.
        # We emit one output row per transcript so that every consequence is visible.
        csq_raw = record.pop('CSQ_RAW', '.')
        if csq_raw in ('.', ''):
            # No annotation (e.g. intergenic variant with no CSQ) — keep one row
            if vep_fields:
                for f in vep_fields:
                    record[f] = '.'
            rows.append(record)
            continue

        for entry in csq_raw.split(','):
            row    = record.copy()
            values = entry.split('|')
            if vep_fields:
                for i, field in enumerate(vep_fields):
                    row[field] = values[i] if i < len(values) else '.'
            else:
                # Fallback when header parsing failed: use positional column names
                for i, v in enumerate(values):
                    row[f'CSQ_{i}'] = v
            rows.append(row)

    return rows


def merge_vep_to_excel():
    """
    Discover all VEP-annotated VCFs under tfm_root, extract and expand their
    variants, add cohort/tissue/sample metadata, and write the consolidated
    result to an Excel sheet sorted by variant impact severity.
    """
    print("Starting VEP annotation consolidation...")

    vcf_files = glob.glob(os.path.join(tfm_root, "**/*.vep.vcf.gz"), recursive=True)
    if not vcf_files:
        print(f"ERROR: no *.vep.vcf.gz files found under {tfm_root}")
        return

    print(f"Found {len(vcf_files)} VEP-annotated VCF(s)")
    all_variants = []

    for vcf_path in sorted(vcf_files):
        parts       = vcf_path.split(os.sep)
        sample_name = os.path.basename(vcf_path).replace('.vep.vcf.gz', '')
        cohort_type = 'Unknown'
        tissue_type = 'Unknown'

        # Infer cohort and tissue labels from the directory structure.
        # Two layouts are supported:
        #   .../tfm/<cohort>/<tissue>/dna_calls/<sample>/<sample>.vep.vcf.gz
        #   .../tfm/<cohort>/<tissue>/<sample>/<sample>.vep.vcf.gz
        try:
            sample_dir_idx = next(i for i, p in enumerate(parts) if p == sample_name)
            if parts[sample_dir_idx - 1].lower() == 'dna_calls':
                tissue_type = parts[sample_dir_idx - 2].capitalize()
                cohort_type = parts[sample_dir_idx - 3].capitalize()
            else:
                tissue_type = parts[sample_dir_idx - 1].capitalize()
                cohort_type = parts[sample_dir_idx - 2].capitalize()
        except (StopIteration, IndexError):
            pass   # labels remain 'Unknown' if the path layout is unexpected

        print(f"  Processing: {cohort_type}/{tissue_type}/{sample_name}")
        rows = bcftools_query_vcf(vcf_path)

        if not rows:
            print(f"    Warning: no variants extracted from {vcf_path}")
            continue

        df = pd.DataFrame(rows)
        df['Sample'] = sample_name
        df['Cohort'] = cohort_type
        df['Tissue'] = tissue_type
        all_variants.append(df)

    if not all_variants:
        print("No variants found across any VCF.")
        return

    final_df = pd.concat(all_variants, ignore_index=True)

    # --- Sort by impact severity ---------------------------------------------
    # VEP impact categories rank variants from most to least likely to affect
    # protein function: HIGH (frameshift, stop-gained) > MODERATE (missense) >
    # LOW (synonymous) > MODIFIER (intergenic, intronic)
    impact_col      = find_col(final_df, 'IMPACT')
    symbol_col      = find_col(final_df, 'SYMBOL')
    consequence_col = find_col(final_df, 'CONSEQUENCE')
    hgvsp_col       = find_col(final_df, 'HGVSp')
    variant_id_col  = find_col(final_df, 'Existing_variation')

    if impact_col:
        impact_order = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
        final_df['_impact_rank'] = final_df[impact_col].map(impact_order).fillna(4)
        final_df = final_df.sort_values('_impact_rank').drop(columns='_impact_rank')

    # --- Column ordering -----------------------------------------------------
    # Place the most clinically relevant fields first for readability in Excel
    priority = ['Cohort', 'Tissue', 'Sample', 'CHROM', 'POS', 'REF', 'ALT',
                'QUAL', 'FILTER', 'DP', 'AF', 'GT']
    for col in [symbol_col, hgvsp_col, consequence_col, impact_col, variant_id_col]:
        if col and col not in priority:
            priority.append(col)
    other_cols = [c for c in final_df.columns if c not in priority]
    final_df = final_df[priority + other_cols]

    # --- Write Excel ---------------------------------------------------------
    # If the file already exists, replace the sheet rather than appending a
    # duplicate; this allows re-running the script without manually deleting output
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    if not os.path.exists(output_file):
        final_df.to_excel(output_file, sheet_name='Biological_Annotations', index=False)
        print(f"Created: {output_file}")
    else:
        with pd.ExcelWriter(output_file, engine='openpyxl', mode='a',
                            if_sheet_exists='replace') as writer:
            final_df.to_excel(writer, sheet_name='Biological_Annotations', index=False)
        print(f"Updated existing file: {output_file}")

    print(f"Consolidated {len(final_df)} variant-transcript rows into 'Biological_Annotations'.")
    print(f"Key columns present: "
          f"REF={('REF' in final_df.columns)}, ALT={('ALT' in final_df.columns)}, "
          f"GT={('GT' in final_df.columns)}, DP={('DP' in final_df.columns)}, "
          f"AF={('AF' in final_df.columns)}")


if __name__ == "__main__":
    merge_vep_to_excel()
