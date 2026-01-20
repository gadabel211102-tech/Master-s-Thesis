
#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Script: 02_dna_qc.sh (CLEAN, FIXED)
# Description:
#   - Runs full DNA-panel QC
#   - Uses *clean* BED (4 columns, chr/start/end/id)
#   - Produces valid mosdepth `.regions.bed.gz` files
#   - Avoids BED corruption, paste mismatches, and malformed outputs
#   - Compatible with 03_qc_analysis.py
# ==============================================================================

# -------------------- Final Thresholds --------------------
MIN_TOTAL_READS=120000
MIN_MAPPED_PCT=90
MIN_ON_TARGET_PCT=75
MIN_MEAN_COVERAGE=200
UNIFORMITY_THRESHOLD=100
MIN_UNIFORMITY_PCT=95
LOW_DEPTH_ALERT=50

# -------------------- Default Values ----------------------
INPUT_DIR=""
OUT_DIR=""
REF_FA=""
BED=""
THREADS="${THREADS:-4}"

usage() {
    cat <<EOF
Usage: $0 -i INPUT_DIR -o OUT_DIR -r REF_FASTA -b TARGETS_BED [-t THREADS]
EOF
}

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

[[ -z "$INPUT_DIR" || -z "$OUT_DIR" || -z "$REF_FA" || -z "$BED" ]] && { usage; exit 1; }

mkdir -p "${OUT_DIR}/pass_bams" "${OUT_DIR}/fail_bams" "${OUT_DIR}/stats"

# ==============================================================================
# STEP 1 — CLEAN BED PREPARATION (NO ANNOTATION PASTING)
# ==============================================================================

# We expect the BED to look like:
#   chr17   start   end   rsID
# Using ONLY columns 1–4 (mosdepth ignores columns >4)
BED_CLEAN="${OUT_DIR}/targets.clean.bed"

echo "[*] Preparing clean 4-column BED (chr/start/end/id)..."
awk '
BEGIN {FS=OFS="\t"}
NF >= 4 {
    # Ensure chr-prefix normalization ONLY
    chr=$1
    if (chr=="17") chr="chr17"
    if (chr=="chr17") chr="chr17"

    print chr, $2, $3, $4
}
' "$BED" > "$BED_CLEAN"

# BED for mosdepth calculation
BED_CALC="${BED_CLEAN}"

# ==============================================================================
# STEP 2 — SAMPLE PROCESSING
# ==============================================================================

mapfile -t BAMS < <(find "$INPUT_DIR" -maxdepth 2 -type f -name "*.bam" | sort)
NUM_BAMS=${#BAMS[@]}

if [[ "$NUM_BAMS" -eq 0 ]]; then
    echo "ERROR: No BAM files found in $INPUT_DIR"
    exit 1
fi

echo "[*] Found $NUM_BAMS samples."

SUMMARY="${OUT_DIR}/qc_summary.tsv"
echo -e "sample\ttotal_reads\tmapped_pct\ton_target_pct\tmean_cov\tuniformity_100x\tworst_amplicon\tstatus\tfail_reasons" \
    > "$SUMMARY"

for (( i=0; i<NUM_BAMS; i++ )); do
    BAM="${BAMS[$i]}"
    SAMPLE=$(basename "$BAM" .bam)
    SDIR="${OUT_DIR}/stats/${SAMPLE}"
    echo "[ $((i+1)) / $NUM_BAMS ] Processing: $SAMPLE"

    rm -rf "$SDIR"
    mkdir -p "$SDIR"

    # -------------------- BAM METRICS --------------------
    FLAGSTAT=$(samtools flagstat -@ "$THREADS" "$BAM")
    TOTAL_READS=$(echo "$FLAGSTAT" | awk '/in total/ {print $1}')
    MAPPED_READS=$(echo "$FLAGSTAT" | awk '/ mapped / && !/primary/ {print $1}')

    MAPPED_PCT=$(awk -v m="$MAPPED_READS" -v t="$TOTAL_READS" \
        'BEGIN{if(t>0){printf "%.2f", (m/t)*100}else{print "0"}}')

    # -------------------- ON-TARGET % --------------------
    ON_TARGET_READS=$(samtools view -@ "$THREADS" -F 4 -c -L "$BED_CALC" "$BAM")
    ON_TARGET_PCT=$(awk -v o="$ON_TARGET_READS" -v m="$MAPPED_READS" \
        'BEGIN{if(m>0){printf "%.2f", (o/m)*100}else{print "0"}}')

    # -------------------- MOSDEPTH -----------------------
    mosdepth --threads "$THREADS" \
             --no-per-base \
             --by "$BED_CALC" \
             --fast-mode \
             "${SDIR}/${SAMPLE}" \
             "$BAM"

    REG="${SDIR}/${SAMPLE}.regions.bed.gz"

    # Mean coverage
    MEAN_COV=$(zcat "$REG" | \
        awk 'BEGIN{sum=0;bp=0} {len=$3-$2; sum+=len*$4; bp+=len} END{
             if (bp>0) printf "%.2f", sum/bp; else print "0"}')

    # Fraction of amplicons ≥100x
    UNIFORMITY=$(zcat "$REG" | \
        awk -v thr="$UNIFORMITY_THRESHOLD" \
            'BEGIN{c=0;tot=0} {tot++; if($4>=thr)c++} END{
             if(tot>0) printf "%.2f", (c/tot)*100; else print "0"}')

    # Minimum depth across all amplicons
    WORST_COV=$(zcat "$REG" | \
        awk 'NR==1{min=$4} {if($4<min)min=$4} END{printf "%.2f", min}')

    # -------------------- FAIL LOGIC --------------------
    STATUS="PASS"
    FAILS=()

    [[ "$TOTAL_READS" -lt "$MIN_TOTAL_READS" ]] && STATUS="FAIL" && FAILS+=("Low reads")
    (( $(echo "$MAPPED_PCT < $MIN_MAPPED_PCT" | bc -l) )) && STATUS="FAIL" && FAILS+=("Low mapping")
    (( $(echo "$ON_TARGET_PCT < $MIN_ON_TARGET_PCT" | bc -l) )) && STATUS="FAIL" && FAILS+=("Low on target")
    (( $(echo "$MEAN_COV < $MIN_MEAN_COVERAGE" | bc -l) )) && STATUS="FAIL" && FAILS+=("Low mean coverage")
    (( $(echo "$UNIFORMITY < $MIN_UNIFORMITY_PCT" | bc -l) )) && STATUS="FAIL" && FAILS+=("Low uniformity")

    # -------------------- WRITE QC SUMMARY --------------------
    printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\n" \
        "$SAMPLE" \
        "$TOTAL_READS" \
        "$MAPPED_PCT" \
        "$ON_TARGET_PCT" \
        "$MEAN_COV" \
        "$UNIFORMITY" \
        "$WORST_COV" \
        "$STATUS" \
        "$(IFS=','; echo "${FAILS[*]-}")" \
        >> "$SUMMARY"

    # -------------------- MOVE BAMS ---------------------
    DEST=$([[ "$STATUS" == "PASS" ]] && echo "${OUT_DIR}/pass_bams" || echo "${OUT_DIR}/fail_bams")
    cp -p "$BAM" "$DEST/"
    [[ -f "${BAM}.bai" ]] && cp -p "${BAM}.bai" "$DEST/"
done

echo "[SUCCESS] DNA QC completed successfully."

