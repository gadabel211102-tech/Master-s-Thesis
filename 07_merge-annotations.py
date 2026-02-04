import pandas as pd
import glob
import os
import gzip

# --- Configuration ---
# Set output path and root directory for the project
results_dir = "/home/gadeaalonsoj/tfm/gsdmb_final_results"
output_file = os.path.join(results_dir, "GSDMB_Annotated_Report_Fixed.xlsx")
tfm_root = "/home/gadeaalonsoj/tfm"

def get_vep_headers(vcf_path):
    """Extracts the exact CSQ field names directly from the VCF header."""
    try:
        # Open the gzipped VCF as plain text
        with gzip.open(vcf_path, 'rt') as f:
            for line in f:
                # Look for the CSQ INFO line which defines the field order
                if line.startswith('##INFO=<ID=CSQ'):
                    # Parse the string after "Format: " to get field names
                    header_part = line.split('Format: ')[1].split('"')[0]
                    return header_part.split('|')
    except Exception as e:
        print(f"Debug: Header extraction failed: {e}")
        return None
    return None

def find_col(df, target_name):
    """Finds a column name in a DataFrame case-insensitively."""
    # Useful because VEP/BCFtools output capitalisation can vary
    for col in df.columns:
        if col.upper() == target_name.upper():
            return col
    return None

def merge_vep_to_excel():
    print("--- Starting VEP Consolidation (Robust Case-Insensitive Version) ---")
    
    # Recursively find all TSV summary files generated in previous steps
    vep_files = glob.glob(os.path.join(tfm_root, "**/*_vep_summary.tsv"), recursive=True)
    if not vep_files:
        print("ERROR: No _vep_summary.tsv files found.")
        return

    # Find VCFs to extract the header format for correct column naming
    vcf_search = glob.glob(os.path.join(tfm_root, "**/*.vep.vcf.gz"), recursive=True)
    if not vcf_search:
        print("ERROR: No VCF files found to extract headers.")
        return
    
    # Get the field names from the first VCF found
    vep_fields = get_vep_headers(vcf_search[0])
    all_variants = []
    
    for f in vep_files:
        parts = f.split(os.sep)
        try:
            # Parse the directory path to assign cohort, tissue, and sample metadata
            idx = parts.index('dna_calls')
            tissue_type = parts[idx-1]
            cohort_type = parts[idx-2]
            sample_name = parts[idx+1]
            
            # Read the current TSV file
            df = pd.read_csv(f, sep='\t', header=None)
            
            # Apply standard VCF columns + VEP fields as headers
            base_cols = ['Chr', 'Pos', 'Ref', 'Alt', 'Qual', 'Filter']
            if vep_fields:
                df.columns = base_cols + vep_fields[:df.shape[1]-len(base_cols)]
            else:
                # Fallback naming if header extraction fails
                df.columns = base_cols + [f'Col_{i}' for i in range(df.shape[1]-len(base_cols))]
            
            # Add the parsed metadata to the dataframe
            df['Sample'] = sample_name
            df['Cohort'] = cohort_type.capitalize()
            df['Tissue'] = tissue_type.capitalize()
            all_variants.append(df)
                
        except (ValueError, IndexError):
            continue

    if all_variants:
        # Combine all sample dataframes into one large table
        final_df = pd.concat(all_variants, ignore_index=True)
        
        # Locate important columns regardless of specific capitalisation
        impact_col = find_col(final_df, 'IMPACT')
        symbol_col = find_col(final_df, 'SYMBOL')
        consequence_col = find_col(final_df, 'CONSEQUENCE')
        hgvsp_col = find_col(final_df, 'HGVSp')
        variant_id_col = find_col(final_df, 'Existing_variation')

        # Custom sorting logic: prioritise HIGH impact variants at the top
        if impact_col:
            impact_order = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
            final_df['impact_rank'] = final_df[impact_col].map(impact_order).fillna(4)
            final_df = final_df.sort_values('impact_rank').drop('impact_rank', axis=1)

        # Define the column order for the final report
        priority = ['Cohort', 'Tissue', 'Sample']
        for col in [symbol_col, hgvsp_col, consequence_col, impact_col, variant_id_col]:
            if col: priority.append(col)
        
        # Move priority columns to the front, followed by everything else
        other_cols = [c for c in final_df.columns if c not in priority]
        final_df = final_df[priority + other_cols]

        # Append/Write the final table to the Excel workbook
        with pd.ExcelWriter(output_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            final_df.to_excel(writer, sheet_name='Biological_Annotations', index=False)
            
        print(f"SUCCESS: Consolidated {len(final_df)} variants into 'Biological_Annotations'.")
        print("No variants found matching the folder criteria.")

if __name__ == "__main__":
    merge_vep_to_excel()
