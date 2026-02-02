import pandas as pd
import glob
import os
import gzip

# --- Configuration ---
results_dir = "/home/gadeaalonsoj/tfm/gsdmb_final_results"
output_file = os.path.join(results_dir, "GSDMB_Annotated_Report_Fixed.xlsx")
tfm_root = "/home/gadeaalonsoj/tfm"

def get_vep_headers(vcf_path):
    """Extracts the exact CSQ field names directly from the VCF header."""
    try:
        with gzip.open(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('##INFO=<ID=CSQ'):
                    header_part = line.split('Format: ')[1].split('"')[0]
                    return header_part.split('|')
    except Exception as e:
        print(f"Debug: Header extraction failed: {e}")
        return None
    return None

def find_col(df, target_name):
    """Finds a column name in a DataFrame case-insensitively."""
    for col in df.columns:
        if col.upper() == target_name.upper():
            return col
    return None

def merge_vep_to_excel():
    print("--- Starting VEP Consolidation (Robust Case-Insensitive Version) ---")
    
    vep_files = glob.glob(os.path.join(tfm_root, "**/*_vep_summary.tsv"), recursive=True)
    if not vep_files:
        print("ERROR: No _vep_summary.tsv files found.")
        return

    vcf_search = glob.glob(os.path.join(tfm_root, "**/*.vep.vcf.gz"), recursive=True)
    if not vcf_search:
        print("ERROR: No VCF files found to extract headers.")
        return
    
    vep_fields = get_vep_headers(vcf_search[0])
    all_variants = []
    
    for f in vep_files:
        parts = f.split(os.sep)
        try:
            idx = parts.index('dna_calls')
            tissue_type = parts[idx-1]
            cohort_type = parts[idx-2]
            sample_name = parts[idx+1]
            
            df = pd.read_csv(f, sep='\t', header=None)
            
            base_cols = ['Chr', 'Pos', 'Ref', 'Alt', 'Qual', 'Filter']
            if vep_fields:
                df.columns = base_cols + vep_fields[:df.shape[1]-len(base_cols)]
            else:
                df.columns = base_cols + [f'Col_{i}' for i in range(df.shape[1]-len(base_cols))]
            
            # Add our custom metadata
            df['Sample'] = sample_name
            df['Cohort'] = cohort_type.capitalize()
            df['Tissue'] = tissue_type.capitalize()
            all_variants.append(df)
                
        except (ValueError, IndexError):
            continue

    if all_variants:
        final_df = pd.concat(all_variants, ignore_index=True)
        
        # Identify the actual names used by Ensembl for key columns
        impact_col = find_col(final_df, 'IMPACT')
        symbol_col = find_col(final_df, 'SYMBOL')
        consequence_col = find_col(final_df, 'CONSEQUENCE')
        hgvsp_col = find_col(final_df, 'HGVSp')
        variant_id_col = find_col(final_df, 'Existing_variation')

        # Sort by Impact if column exists
        if impact_col:
            impact_order = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
            final_df['impact_rank'] = final_df[impact_col].map(impact_order).fillna(4)
            final_df = final_df.sort_values('impact_rank').drop('impact_rank', axis=1)

        # Build priority list using ONLY columns that actually exist
        priority = ['Cohort', 'Tissue', 'Sample']
        for col in [symbol_col, hgvsp_col, consequence_col, impact_col, variant_id_col]:
            if col: priority.append(col)
        
        # Final column arrangement
        other_cols = [c for c in final_df.columns if c not in priority]
        final_df = final_df[priority + other_cols]

        # Write to Excel
        with pd.ExcelWriter(output_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            final_df.to_excel(writer, sheet_name='Biological_Annotations', index=False)
            
        print(f"SUCCESS: Consolidated {len(final_df)} variants into 'Biological_Annotations'.")
    else:
        print("No variants found matching the folder criteria.")

if __name__ == "__main__":
    merge_vep_to_excel()