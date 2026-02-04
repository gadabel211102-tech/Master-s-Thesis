#!/usr/bin/env bash  # Shebang: Ensures the script runs in the Bash shell environment
set -euo pipefail    # Strict mode: Exit on Error, Unset vars, or Pipe failures

# Usage: Instructions for the user on how to provide the input directory and output filename
# Example: Shows a real-world path example for the script to run correctly

QC_DIR="$1"          # Assigns the 1st command-line argument to QC_DIR (the input folder)
OUT="$2"             # Assigns the 2nd command-line argument to OUT (the result file)

QC_SUM="${QC_DIR}/qc_summary.tsv"  # Constructs path to the summary table (Tab Separated Values)
PASS_DIR="${QC_DIR}/pass_bams"      # Constructs path to the folder containing the actual data files

[[ -f "$QC_SUM" ]]      || { echo "ERROR: qc_summary.tsv not found: $QC_SUM"; exit 1; }  # Check if file exists; if not, exit
[[ -d "$PASS_DIR" ]]    || { echo "ERROR: pass_bams/ not found: $PASS_DIR"; exit 1; }   # Check if directory exists; if not, exit

echo "[INFO] Reading: $QC_SUM"      # Prints status message to the console
echo "[INFO] Looking for PASS samples..."

# Extract PASS sample names using 'awk' 
SAMPLES=$(awk -F'\t' '  # -F"\t" sets the field separator to a "Tab" character
  NR==1 {               # "NR==1" means: Process the first line (the header)
    for(i=1;i<=NF;i++) h[$i]=i  # Loop through every column and store its name/index in the "h" array
    if(!("sample" in h) || !("status" in h)) { # If "sample" or "status" columns are missing...
      print "ERROR: Missing sample/status columns in qc_summary.tsv" >"/dev/stderr" # ...print error to stderr
      exit 2            # Exit with code 2 (misconfiguration)
    }
    next                # Finished with header, skip to the next line of the file
  }
  $h["status"]=="PASS" { print $h["sample"] } # If the "status" column equals "PASS", print the sample name
' "$QC_SUM")            # Feed the summary file into the awk command above

# Build manifest using actual files in pass_bams/
echo "[INFO] Writing manifest → $OUT" # Status update
: > "$OUT"              # The ":" is a "null" command; combined with ">", it empties/creates the output file

for sample in $SAMPLES; do  # Start a loop for every sample name stored in the SAMPLES variable
  bam="${PASS_DIR}/${sample}.bam" # Construct the expected path to the .bam file
  if [[ -f "$bam" ]]; then  # If (-f) the file actually exists on the hard drive...
    echo "$bam" >> "$OUT"   # ...append (>>) the path to our manifest file
  else                      # If the sample is "PASS" in the table but the file is missing...
    echo "[WARN] Missing BAM for PASS sample: $sample" >&2 # ...warn the user in stderr
  fi                        # End 'if' block
done                        # End 'for' loop

echo "[INFO] Manifest done." # Completion message
echo "[INFO] Number of BAMs: $(wc -l < "$OUT")" # Use 'wc -l' (word count - lines) to count entries in the manifest
