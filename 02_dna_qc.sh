#!/usr/bin/env bash
# Ensures the script stops immediately if any command fails (-e), 
# if a variable is missing (-u), or if a pipe command fails (-o pipefail).
set -euo pipefail

# -------------------- Final Thresholds --------------------
# These are the "Rules of the Lab." If a sample doesn't hit these numbers, it fails.
MIN_TOTAL_READS=120000      # Must have at least 120k reads
MIN_MAPPED_PCT=90           # 90% of reads must align to the human genome
MIN_ON_TARGET_PCT=75        # 75% of those reads must hit our specific panel
MIN_MEAN_COVERAGE=200       # Average depth must be at least 200x
UNIFORMITY_THRESHOLD=100    # We want to see how many regions hit 100x
MIN_UNIFORMITY_PCT=95       # 95% of regions should be above that threshold
LOW_DEPTH_ALERT=50          # Anything below 50x is considered a "gap"

# 

# -------------------- Default Values ----------------------
# Setting up empty placeholders for the inputs you provide via the terminal.
INPUT_DIR=""
OUT_DIR=""
REF_FA=""
BED=""
THREADS="${THREADS:-4}" # Use 4 computer cores by default unless told otherwise

# Standard "How to use this script" message
usage() {
    cat <<EOF
Usage: $0 -i INPUT_DIR -o OUT_DIR -r REF_FASTA -b TARGETS_BED [-t THREADS]
EOF
}

# Logic to read the flags (like -i or -o) you type when running the script
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT_DIR="$2"; shift 2 ;;
        -o) OUT_DIR="$2"; shift 2 ;;
        -r) REF_FA="$2"; shift 2 ;;
        -b) BED="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        *) usage; exit 2 ;;
    esac
done

# If any required input is missing, show the usage message and quit.
[[ -z "$INPUT_DIR" || -z "$OUT_DIR" || -z "$REF_FA" || -z "$BED" ]] && { usage; exit 1; }

# Create the folder structure for the results
mkdir -p "${OUT_DIR}/pass_bams" "${OUT_DIR}/fail_bams" "${OUT_DIR}/stats"

# ==============================================================================
# STEP 1 — CLEAN BED PREPARATION
# ==============================================================================
# "mosdepth" (the tool we use for coverage) is very picky. It only wants 4 columns.
# This part "cleans" your target file to ensure it is exactly: chr, start, end, id.
BED_CLEAN="${OUT_DIR}/targets.clean.bed"

echo "[*] Preparing clean 4-column BED (chr/start/end/id)..."
awk '
BEGIN {FS=OFS="\t"}
NF >= 4 {
    chr=$1
    # Standardise chromosome names: if it says "17", change it to "chr17"
    if (chr=="17") chr="chr17"
    print chr, $2, $3, $4
}
' "$BED" > "$BED_CLEAN"

BED_CALC="${BED_CLEAN}"

# ==============================================================================
# STEP 2 — SAMPLE PROCESSING
# ==============================================================================

# Find all BAM files in the input folder and put them in a list (array)
mapfile -t BAMS < <(find "$INPUT_DIR" -maxdepth 2 -type f -name "*.bam" | sort)
NUM_BAMS=${#BAMS[@]}

if [[ "$NUM_BAMS" -eq 0 ]]; then
    echo "ERROR: No BAM files found in $INPUT_DIR"
    exit 1
fi

echo "[*] Found $NUM_BAMS samples."

# Create the header for our summary spreadsheet
SUMMARY="${OUT_DIR}/qc_summary.tsv"
echo -e "sample\ttotal_reads\tmapped_pct\ton_target_pct\tmean_cov\tuniformity_100x\tworst_amplicon\tstatus\tfail_reasons" \
    > "$SUMMARY"

# Start a loop: Process every BAM file one by one
for (( i=0; i<NUM_BAMS; i++ )); do
    BAM="${BAMS[$i]}"
    SAMPLE=$(basename "$BAM" .bam)
    SDIR="${OUT_DIR}/stats/${SAMPLE}"
    echo "[ $((i+1)) / $NUM_BAMS ] Processing: $SAMPLE"

    mkdir -p "$SDIR"

    # -------------------- BAM METRICS --------------------
    # Use 'samtools' to count how many reads are in the file and how many mapped
    FLAGSTAT=$(samtools flagstat -@ "$THREADS" "$BAM")
    TOTAL_READS=$(echo "$FLAGSTAT" | awk '/in total/ {print $1}')
    MAPPED_READS=$(echo "$FLAGSTAT" | awk '/ mapped / && !/primary/ {print $1}')

    # Calculate the % of mapped reads
    MAPPED_PCT=$(awk -v m="$MAPPED_READS" -v t="$TOTAL_READS" \
        'BEGIN{if(t>0){printf "%.2f", (m/t)*100}else{print "0"}}')

    # -------------------- ON-TARGET % --------------------
    # Count reads that physically overlap with our BED targets (on-target)
    ON_TARGET_READS=$(samtools view -@ "$THREADS" -F 4 -c -L "$BED_CALC" "$BAM")
    ON_TARGET_PCT=$(awk -v o="$ON_TARGET_READS" -v m="$MAPPED_READS" \
        'BEGIN{if(m>0){printf "%.2f", (o/m)*100}else{print "0"}}')

    # -------------------- MOSDEPTH -----------------------
    # The industry standard tool for calculating depth of coverage.
    # It generates a compressed file (.regions.bed.gz) with depth for every region.
    mosdepth --threads "$THREADS" \
             --no-per-base \
             --by "$BED_CALC" \
             --fast-mode \
             "${SDIR}/${SAMPLE}" \
             "$BAM"

    REG="${SDIR}/${SAMPLE}.regions.bed.gz"

    # Mean coverage calculation: (Total bases covered) / (Total length of targets)
    MEAN_COV=$(zcat "$REG" | \
        awk 'BEGIN{sum=0;bp=0} {len=$3-$2; sum+=len*$4; bp+=len} END{
             if (bp>0) printf "%.2f", sum/bp; else print "0"}')

    # Uniformity calculation: % of amplicons that have at least 100x depth
    UNIFORMITY=$(zcat "$REG" | \
        awk -v thr="$UNIFORMITY_THRESHOLD" \
            'BEGIN{c=0;tot=0} {tot++; if($4>=thr)c++} END{
             if(tot>0) printf "%.2f", (c/tot)*100; else print "0"}')

    # Find the "Worst Amplicon": the single lowest coverage value in the whole sample
    WORST_COV=$(zcat "$REG" | \
        awk 'NR==1{min=$4} {if($4<min)min=$4} END{printf "%.2f", min}')
   # -------------------- FAIL LOGIC --------------------
    # Start with the assumption that the sample is "PASS"
    STATUS="PASS"
    FAILS=() # An empty list (array) to store specific reasons if it fails

    # Check 1: Total Raw Reads
    # If total reads are less than 120,000, mark as FAIL and record the reason.
    [[ "$TOTAL_READS" -lt "$MIN_TOTAL_READS" ]] && STATUS="FAIL" && FAILS+=("Low reads")

    # Check 2: Mapping % (using 'bc' because Bash cannot handle decimals/floats naturally)
    (( $(echo "$MAPPED_PCT < $MIN_MAPPED_PCT" | bc -l) )) && STATUS="FAIL" && FAILS+=("Low mapping")

    # Check 3: On-Target %
    (( $(echo "$ON_TARGET_PCT < $MIN_ON_TARGET_PCT" | bc -l) )) && STATUS="FAIL" && FAILS+=("Low on target")

    # Check 4: Mean Coverage (e.g., must be > 200x)
    (( $(echo "$MEAN_COV < $MIN_MEAN_COVERAGE" | bc -l) )) && STATUS="FAIL" && FAILS+=("Low mean coverage")

    # Check 5: Uniformity (e.g., 95% of targets must have at least 100x depth)
    (( $(echo "$UNIFORMITY < $MIN_UNIFORMITY_PCT" | bc -l) )) && STATUS="FAIL" && FAILS+=("Low uniformity")

    # -------------------- WRITE QC SUMMARY --------------------
    # Use 'printf' to format the output into clean columns. 
    # This appends one row per sample to the 'qc_summary.tsv' file.
    printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\n" \
        "$SAMPLE" \
        "$TOTAL_READS" \
        "$MAPPED_PCT" \
        "$ON_TARGET_PCT" \
        "$MEAN_COV" \
        "$UNIFORMITY" \
        "$WORST_COV" \
        "$STATUS" \
        "$(IFS=','; echo "${FAILS[*]-}")" \ # Joins all failure reasons with a comma
        >> "$SUMMARY"

    # -------------------- MOVE BAMS ---------------------
    # Decide the destination folder based on the PASS/FAIL status.
    # This is excellent for downstream automation (e.g., only running variant calling on Pass BAMs).
    DEST=$([[ "$STATUS" == "PASS" ]] && echo "${OUT_DIR}/pass_bams" || echo "${OUT_DIR}/fail_bams")
    
    # Copy the BAM file and its index (.bai) to the Pass or Fail folder.
    # '-p' preserves the original timestamps.
    cp -p "$BAM" "$DEST/"
    [[ -f "${BAM}.bai" ]] && cp -p "${BAM}.bai" "$DEST/"
done

echo "[SUCCESS] DNA QC completed successfully."
