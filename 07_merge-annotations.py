import pandas as pd
import glob
import os
import gzip
import subprocess

# --- Configuration ---
results_dir = "/home/gadeaalonsoj/tfm/gsdmb_final_results"
output_file = os.path.join(results_dir, "GSDMB_Annotated_Report_Fixed.xlsx")
tfm_root = "/home/gadeaalonsoj/tfm"


def get_vep_headers(vcf_path):
    """Extracts the exact CSQ field names directly from the VCF header."""
    try:
        opener = gzip.open if vcf_path.endswith('.gz') else open
        with opener(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('##INFO=<ID=CSQ'):
                    header_part = line.split('Format: ')[1].split('"')[0]
                    return header_part.strip().split('|')
    except Exception as e:
        print(f"Debug: Header extraction failed for {vcf_path}: {e}")
    return None


def find_col(df, target_name):
    """Finds a column name in a DataFrame case-insensitively."""
    for col in df.columns:
        if col.upper() == target_name.upper():
            return col
    return None


def bcftools_query_vcf(vcf_path):
    """
    Extract per-variant data from a VEP-annotated VCF using bcftools query.
    Returns a list of dicts, one per variant x transcript (CSQ split properly).
    GT, DP, AF, REF, ALT extracted as dedicated fields - NOT via CSQ pipe-splitting.
    """
    try:
        hdr = subprocess.check_output(
            ['bcftools', 'view', '-h', vcf_path], stderr=subprocess.DEVNULL
        ).decode()
    except subprocess.CalledProcessError:
        print(f"  WARNING: could not read header for {vcf_path}")
        return []

    has_gt  = '##FORMAT=<ID=GT,'  in hdr
    has_dp  = '##FORMAT=<ID=DP,'  in hdr
    has_af  = '##FORMAT=<ID=AF,'  in hdr
    has_ao  = '##FORMAT=<ID=AO,'  in hdr
    has_ro  = '##FORMAT=<ID=RO,'  in hdr
    has_iaf = '##INFO=<ID=AF,'    in hdr
    has_idp = '##INFO=<ID=DP,'    in hdr

    # Build format string: extract CSQ as single raw field, everything else explicit
    fmt_parts = ['%CHROM', '%POS', '%REF', '%ALT', '%QUAL', '%FILTER']

    # Depth: prefer FORMAT/DP, fall back to INFO/DP
    if has_dp:
        fmt_parts.append('[%DP]')
    elif has_idp:
        fmt_parts.append('%INFO/DP')
    else:
        fmt_parts.append('.')

    # Allele fraction: prefer FORMAT/AF, derive from AO/RO if needed, fall back to INFO/AF
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

    # Genotype
    if has_gt:
        fmt_parts.append('[%GT]')
    else:
        fmt_parts.append('.')

    # CSQ last - keep intact so we split on | ourselves
    fmt_parts.append('%INFO/CSQ')

    fmt_string = '\t'.join(fmt_parts) + '\n'

    try:
        raw = subprocess.check_output(
            ['bcftools', 'query', '-f', fmt_string, vcf_path],
            stderr=subprocess.DEVNULL
        ).decode()
    except subprocess.CalledProcessError as e:
        print(f"  WARNING: bcftools query failed for {vcf_path}: {e}")
        return []

    # Determine column layout
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
        parts = line.split('\t', len(meta_cols) - 1)
        if len(parts) < len(meta_cols):
            parts += ['.'] * (len(meta_cols) - len(parts))
        record = dict(zip(meta_cols, parts))

        # Derive AF from AO/RO if needed
        if ao_ro_mode:
            try:
                ao = float(record.get('AO', 0))
                ro = float(record.get('RO', 0))
                record['AF'] = f"{ao / (ao + ro):.4f}" if (ao + ro) > 0 else '.'
            except (ValueError, TypeError):
                record['AF'] = '.'
            del record['AO'], record['RO']

        # Expand CSQ: each transcript is a comma-delimited entry
        csq_raw = record.pop('CSQ_RAW', '.')
        if csq_raw in ('.', ''):
            if vep_fields:
                for f in vep_fields:
                    record[f] = '.'
            rows.append(record)
            continue

        csq_entries = csq_raw.split(',')
        for entry in csq_entries:
            row = record.copy()
            values = entry.split('|')
            if vep_fields:
                for i, field in enumerate(vep_fields):
                    row[field] = values[i] if i < len(values) else '.'
            else:
                for i, v in enumerate(values):
                    row[f'CSQ_{i}'] = v
            rows.append(row)

    return rows


def merge_vep_to_excel():
    print("--- Starting VEP Consolidation (VCF-direct version) ---")

    # Find all VEP-annotated VCFs written by script 06
    vcf_files = glob.glob(os.path.join(tfm_root, "**/*.vep.vcf.gz"), recursive=True)

    # BUG FIX: old version returned immediately after the fallback search even
    # if it found files. Now we only bail if truly nothing is found.
    if not vcf_files:
        print("ERROR: No *.vep.vcf.gz files found under", tfm_root)
        return

    print(f"Found {len(vcf_files)} VEP-annotated VCF(s)")

    all_variants = []

    for vcf_path in sorted(vcf_files):
        parts = vcf_path.split(os.sep)

        # Expected layout: .../tfm/<cohort>/<tissue>/dna_calls/<sample>/<sample>.vep.vcf.gz
        # or: .../tfm/<cohort>/<tissue>/<sample>/<sample>.vep.vcf.gz
        sample_name = os.path.basename(vcf_path).replace('.vep.vcf.gz', '')
        cohort_type = 'Unknown'
        tissue_type = 'Unknown'

        try:
            sample_dir_idx = next(
                i for i, p in enumerate(parts) if p == sample_name
            )
            # Expected: .../cohort/tissue/dna_calls/sample/sample.vep.vcf.gz
            # so parts[sample_dir_idx - 1] == 'dna_calls'
            #    parts[sample_dir_idx - 2] == tissue  (e.g. 'normal', 'tumour')
            #    parts[sample_dir_idx - 3] == cohort  (e.g. 'breast', 'endometrium')
            if parts[sample_dir_idx - 1].lower() == 'dna_calls':
                tissue_type = parts[sample_dir_idx - 2].capitalize()
                cohort_type = parts[sample_dir_idx - 3].capitalize()
            else:
                # Fallback for flatter layouts without dna_calls folder
                tissue_type = parts[sample_dir_idx - 1].capitalize()
                cohort_type = parts[sample_dir_idx - 2].capitalize()
        except (StopIteration, IndexError):
            pass

        print(f"  Processing: {cohort_type}/{tissue_type}/{sample_name}")
        rows = bcftools_query_vcf(vcf_path)

        if not rows:
            print(f"    WARNING: no variants extracted from {vcf_path}")
            continue

        df = pd.DataFrame(rows)
        df['Sample'] = sample_name
        df['Cohort'] = cohort_type
        df['Tissue'] = tissue_type
        all_variants.append(df)

    if not all_variants:
        print("No variants found in any VCF.")
        return

    final_df = pd.concat(all_variants, ignore_index=True)

    # --- Column ordering ---
    impact_col      = find_col(final_df, 'IMPACT')
    symbol_col      = find_col(final_df, 'SYMBOL')
    consequence_col = find_col(final_df, 'CONSEQUENCE')
    hgvsp_col       = find_col(final_df, 'HGVSp')
    variant_id_col  = find_col(final_df, 'Existing_variation')

    # Sort by impact
    if impact_col:
        impact_order = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
        final_df['_impact_rank'] = final_df[impact_col].map(impact_order).fillna(4)
        final_df = final_df.sort_values('_impact_rank').drop('_impact_rank', axis=1)

    # Most useful columns first, then everything else
    priority = ['Cohort', 'Tissue', 'Sample', 'CHROM', 'POS', 'REF', 'ALT',
                'QUAL', 'FILTER', 'DP', 'AF', 'GT']
    for col in [symbol_col, hgvsp_col, consequence_col, impact_col, variant_id_col]:
        if col and col not in priority:
            priority.append(col)

    other_cols = [c for c in final_df.columns if c not in priority]
    final_df = final_df[priority + other_cols]

    # --- Write Excel ---
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    if not os.path.exists(output_file):
        final_df.to_excel(output_file, sheet_name='Biological_Annotations', index=False)
        print(f"--- Created new report: {output_file} ---")
    else:
        with pd.ExcelWriter(output_file, engine='openpyxl', mode='a',
                            if_sheet_exists='replace') as writer:
            final_df.to_excel(writer, sheet_name='Biological_Annotations', index=False)
        print(f"--- Updated existing report: {output_file} ---")

    print(f"SUCCESS: Consolidated {len(final_df)} variant-transcript rows into "
          f"'Biological_Annotations'.")
    print(f"Columns available for QC: "
          f"REF={'REF' in final_df.columns}, ALT={'ALT' in final_df.columns}, "
          f"GT={'GT' in final_df.columns}, DP={'DP' in final_df.columns}, "
          f"AF={'AF' in final_df.columns}")


if __name__ == "__main__":
    merge_vep_to_excel()