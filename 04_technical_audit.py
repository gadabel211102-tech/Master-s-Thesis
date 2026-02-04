import pandas as pd
import glob
import os
import re

# --- Configuration ---
results_dir = "/home/gadeaalonsoj/tfm/gsdmb_final_results"
bed_file_path = "IAD255368_167_Submitted.bed"
output_file = os.path.join(results_dir, "GSDMB_Annotated_Report_Fixed.xlsx")
# Look for all "Full_Panel" CSVs produced by the previous steps
panel_files = glob.glob(os.path.join(results_dir, "GSDMB_*_Full_Panel.csv"))

def clean_name(name):
    """Removes invalid Excel characters and limits length."""
    # Excel sheet names can't have certain characters (like / or :) and must be < 31 chars
    return re.sub(r'[\\/*?:\[\]]', '', name)[:31]

def load_bed_annotations(path):
    """Loads BED file and combines ID and Gene Name columns."""
    # We pull specific columns: 0(Chr), 1(Start), 2(End), 4(rsID/SNP), 5(Gene)
    bed_df = pd.read_csv(path, sep='\t', comment='t', header=None, 
                         usecols=[0, 1, 2, 4, 5],
                         names=['chr', 'start', 'end', 'rs_id', 'gene'])
    
    # SMART ANNOTATION: If there is an rsID (like rs1234), use it. 
    # If the rsID is just a dot ('.'), use the Gene Name instead.
    bed_df['annotation'] = bed_df.apply(
        lambda row: row['rs_id'] if row['rs_id'] != '.' else row['gene'], axis=1
    )
    return bed_df[['chr', 'start', 'end', 'annotation']]

if not panel_files:
    print("Error: No 'Full_Panel.csv' files found.")
else:
    # Load the "key" (the BED file) that explains what each coordinate actually is
    annotations_df = load_bed_annotations(bed_file_path)
    
    # Create the Excel file using 'xlsxwriter' (more stable for large genetic datasets)
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        for file in panel_files:
            # Extract cohort name from filename (e.g., GSDMB_Control_Full_Panel.csv -> Control)
            cohort_name = os.path.basename(file).replace("GSDMB_", "").replace("_Full_Panel.csv", "")
            safe_sheet_name = clean_name(cohort_name)
            
            df = pd.read_csv(file)

            # DATA CLEANING: Ensure coordinates are integers so the merge doesn't fail
            for col in ['start', 'end']:
                df[col] = df[col].astype(int)
                annotations_df[col] = annotations_df[col].astype(int)

            # THE MERGE: This combines the Depth data with the Annotation data 
            # by matching the Chromosome, Start, and End positions.
            df_merged = pd.merge(df, annotations_df, on=['chr', 'start', 'end'], how='left')

            # Identify which columns are actually sample sequencing data
            metadata = ['chr', 'start', 'end', 'id', 'annotation']
            sample_cols = [c for c in df_merged.columns if c not in metadata]

            # CALCULATE STATISTICS:
            # Find the median depth for every single amplicon across all samples
            df_merged['median_depth'] = df_merged[sample_cols].median(axis=1)
            
            # THE TROUBLESHOOTER: Sort by lowest depth and pick the 10 amplicons that performed worst
            worst_10 = df_merged.sort_values('median_depth').head(10).copy()

            # Prepare columns for the final Excel sheet
            final_cols = ['chr', 'start', 'end', 'annotation', 'median_depth']
            valid_cols = [c for c in final_cols if c in worst_10.columns]
            
            # Save this cohort's Top 10 worst amplicons to a dedicated tab in Excel
            worst_10[valid_cols].to_excel(writer, sheet_name=safe_sheet_name, index=False)
            print(f"Sheet '{safe_sheet_name}' added successfully.")

    print(f"\nReport completed: {output_file}")
