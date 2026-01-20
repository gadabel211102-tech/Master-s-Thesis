import pandas as pd
import glob
import os
import re

# --- Configuration ---
results_dir = "/home/gadeaalonsoj/tfm/gsdmb_final_results"
bed_file_path = "IAD255368_167_Submitted.bed"
output_file = os.path.join(results_dir, "GSDMB_Annotated_Report_Fixed.xlsx")
panel_files = glob.glob(os.path.join(results_dir, "GSDMB_*_Full_Panel.csv"))

def clean_name(name):
    """Removes invalid Excel characters and limits length."""
    return re.sub(r'[\\/*?:\[\]]', '', name)[:31]

def load_bed_annotations(path):
    """Loads BED file and combines ID and Gene Name columns."""
    # Loading col 0 (chr), 1 (start), 2 (end), 4 (RS ID), 5 (Gene Name)
    bed_df = pd.read_csv(path, sep='\t', comment='t', header=None, 
                         usecols=[0, 1, 2, 4, 5],
                         names=['chr', 'start', 'end', 'rs_id', 'gene'])
    
    # Create a single 'annotation' column: use rs_id if it's not '.', otherwise use gene
    bed_df['annotation'] = bed_df.apply(
        lambda row: row['rs_id'] if row['rs_id'] != '.' else row['gene'], axis=1
    )
    return bed_df[['chr', 'start', 'end', 'annotation']]

if not panel_files:
    print("Error: No 'Full_Panel.csv' files found.")
else:
    annotations_df = load_bed_annotations(bed_file_path)
    
    # Using 'xlsxwriter' engine often solves 'found a problem' errors better than openpyxl
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        for file in panel_files:
            cohort_name = os.path.basename(file).replace("GSDMB_", "").replace("_Full_Panel.csv", "")
            safe_sheet_name = clean_name(cohort_name)
            
            df = pd.read_csv(file)

            # Ensure coordinates are same type for merging
            for col in ['start', 'end']:
                df[col] = df[col].astype(int)
                annotations_df[col] = annotations_df[col].astype(int)

            # Left merge ensures we keep all depth data even if BED match is missing
            df_merged = pd.merge(df, annotations_df, on=['chr', 'start', 'end'], how='left')

            # Identify numeric sample columns
            metadata = ['chr', 'start', 'end', 'id', 'annotation']
            sample_cols = [c for c in df_merged.columns if c not in metadata]

            # Calculate Median and get Top 10 worst
            df_merged['median_depth'] = df_merged[sample_cols].median(axis=1)
            worst_10 = df_merged.sort_values('median_depth').head(10).copy()

            # Final clean export to Excel
            final_cols = ['chr', 'start', 'end', 'annotation', 'median_depth']
            # Only use columns that exist to prevent errors
            valid_cols = [c for c in final_cols if c in worst_10.columns]
            
            worst_10[valid_cols].to_excel(writer, sheet_name=safe_sheet_name, index=False)
            print(f"Sheet '{safe_sheet_name}' added successfully.")

    print(f"\nReport completed: {output_file}")
