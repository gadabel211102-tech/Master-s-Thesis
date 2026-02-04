#!/usr/bin/env bash  # Shebang: Runs the script using the Bash interpreter
set -euo pipefail    # Strict mode: Exit on Error, Unset variables, or Pipe failures

# ---------------------------------------------------------------
# Zero‑coverage analysis per amplicon, per sample, per cohort (FAST)
# ---------------------------------------------------------------
# Usage:
#    02.2_zero_base_coverage_fast.sh <ROOT_DIR> <OUT_DIR>
# ---------------------------------------------------------------

ROOT="$1"            # First argument: The main folder containing the data
OUT="$2"             # Second argument: Where to save the final results

# Define paths for the BED file (genomic regions) and BAM folders (sequencing data)
BED="${ROOT}/dna_qc/targets.clean.bed"
PASS="${ROOT}/dna_qc/pass_bams"
FAIL="${ROOT}/dna_qc/fail_bams"

mkdir -p "$OUT"      # Creates the output directory; -p prevents error if it already exists
TMP="$OUT/tmp"       # Path for a temporary folder to store intermediate files
mkdir -p "$TMP"      # Creates the temporary folder

# Collect samples: Finds all .bam files in both PASS and FAIL folders, sorts them, and saves to array BAMS
mapfile -t BAMS < <(find "$PASS" "$FAIL" -type f -name "*.bam" | sort)

if [[ ${#BAMS[@]} -eq 0 ]]; then  # If the count of BAMS (${#BAMS[@]}) is equal to 0...
    echo "ERROR: no BAMs found in $ROOT"
    exit 1           # Stop the script with an error
fi

echo "[INFO] Found ${#BAMS[@]} samples"
echo "[INFO] Using BED: $BED"

# Read BED: Loads the first 4 columns of the BED file into a list (array) called BED_ROWS
mapfile -t BED_ROWS < <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4}' "$BED")

# Build output header: Initialises the top row of the spreadsheet
HEADER="chr\tstart\tend\tannotation"
for BAM in "${BAMS[@]}"; do      # Loop through each BAM file to get sample names
    S=$(basename "$BAM" .bam)    # Strip the path and extension to get just the sample name (e.g., 'SampleA')
    HEADER="${HEADER}\t${S}"     # Add the sample name as a new column in the header
done

# Create a cohort name based on directory names and define the final output file path
COHORT_NAME=$(basename "$(dirname "$ROOT")")_$(basename "$ROOT")
OUTFILE="${OUT}/zero_cov_${COHORT_NAME}.tsv"
echo -e "$HEADER" > "$OUTFILE"   # Write the header to the file (-e enables tab interpretation)

echo "[INFO] Output → $OUTFILE"
echo "[INFO] Temporary mosdepth directory → $TMP"

# ---------------------------------------------------------------
# STEP 1 — Compute per-base mosdepth for each sample ONCE
# ---------------------------------------------------------------

declare -A PERBASE   # Creates an "Associative Array" (like a dictionary) to map samples to files

echo "[INFO] Running per-base mosdepth (one per sample)..."

for BAM in "${BAMS[@]}"; do      # For every BAM file...
    S=$(basename "$BAM" .bam)    # Get the sample name
    PREFIX="${TMP}/${S}"         # Define the naming prefix for the output

    if [[ ! -f "${PREFIX}.per-base.bed.gz" ]]; then  # If the mosdepth output doesn't exist yet...
        mosdepth --threads 4 \   # Run mosdepth using 4 CPU threads
            --fast-mode \        # Optimising for speed (skips some slower stats)
            "${PREFIX}" \        # Output filename prefix
            "$BAM" > /dev/null 2>&1 # Run silently (hide normal output and errors)
    fi

    PERBASE["$S"]="${PREFIX}.per-base.bed.gz" # Save the file path into our dictionary
done

echo "[INFO] All per-base files ready."

# ---------------------------------------------------------------
# STEP 2 — Iterate over amplicons and count zero bases quickly
# ---------------------------------------------------------------

echo "[INFO] Computing zero-base counts..."

for ROW in "${BED_ROWS[@]}"; do  # Loop through every genomic region in the BED file
    # Split the tab-separated ROW into 4 specific variables
    IFS=$'\t' read -r chr start end annot <<< "$ROW"
    line="${chr}\t${start}\t${end}\t${annot}" # Start building the output row with region info

    for BAM in "${BAMS[@]}"; do  # Now, for this specific region, check every sample
        S=$(basename "$BAM" .bam)
        PB="${PERBASE[$S]}"      # Get the path to the mosdepth file for this sample

        # Count depth==0 bases for this amplicon:
        # zcat: opens compressed file; awk: checks if region matches AND depth ($4) is 0
        z=$(zcat "$PB" | awk -v c="$chr" -v s="$start" -v e="$end" '
            $1==c && $2>=s && $3<=e && $4==0 {count++} 
            END {print count+0}') # Print the count (adding 0 ensures "0" is printed if empty)
        
        line="${line}\t${z}"     # Append the zero-count for this sample to the row
    done

    echo -e "$line" >> "$OUTFILE" # Write the completed row for this region to the output file
done

echo "[SUCCESS] FAST zero-base analysis completed."
echo "[NOTE] Per-base files stored in $TMP."
