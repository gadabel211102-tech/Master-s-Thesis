#!/usr/bin/env bash
# Post-QC analysis run after 02_dna_qc.sh. Performs two tasks:
#   1. Generates a manifest of BAM paths for all samples that passed QC
#   2. Computes per-amplicon zero- and low-coverage metrics across all samples
# Usage: ./02b_more-qc.sh -q QC_DIR -m MANIFEST_OUT -z ZERO_COV_OUT [-t THREADS]
set -euo pipefail

# Print the line number of any failure to help with debugging
trap 'echo "[ERROR] Script failed at line $LINENO" >&2' ERR

# --- Input parameters --------------------------------------------------------
QC_DIR=""
MANIFEST_OUT=""
ZERO_COV_OUT=""
THREADS="${THREADS:-4}"

usage() {
    cat <<USAGE
Usage: $0 -q QC_DIR -m MANIFEST_OUT -z ZERO_COV_OUT [-t THREADS]
  -q  Output directory from 02_dna_qc.sh (must contain qc_summary.tsv,
      pass_bams/, fail_bams/, targets.sorted.bed)
  -m  Output path for the PASS manifest (one BAM path per line)
  -z  Output directory for zero-coverage results
  -t  CPU threads for mosdepth (default: 4)
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -q) QC_DIR="$2";        shift 2 ;;
        -m) MANIFEST_OUT="$2";  shift 2 ;;
        -z) ZERO_COV_OUT="$2";  shift 2 ;;
        -t) THREADS="$2";       shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1"; usage; exit 2 ;;
    esac
done

[[ -z "$QC_DIR" || -z "$MANIFEST_OUT" || -z "$ZERO_COV_OUT" ]] && { usage; exit 1; }

# Locate expected outputs from the previous QC step
QC_SUMMARY="${QC_DIR}/qc_summary.tsv"
PASS_DIR="${QC_DIR}/pass_bams"
FAIL_DIR="${QC_DIR}/fail_bams"
BED_FILE="${QC_DIR}/targets.sorted.bed"

[[ -f "$QC_SUMMARY" ]] || { echo "ERROR: qc_summary.tsv not found: $QC_SUMMARY" >&2; exit 1; }
[[ -d "$PASS_DIR"   ]] || { echo "ERROR: pass_bams/ not found: $PASS_DIR" >&2; exit 1; }
[[ -f "$BED_FILE"   ]] || { echo "ERROR: targets.sorted.bed not found: $BED_FILE" >&2; exit 1; }
[[ -d "$FAIL_DIR"   ]] || echo "WARNING: fail_bams/ not found — PASS samples only will be analysed"

mkdir -p "$(dirname "$MANIFEST_OUT")" "$ZERO_COV_OUT"

# =============================================================================
# Part 1: Generate manifest of PASS sample BAM paths
# =============================================================================
echo ""
echo "========================================"
echo "Part 1: Generating PASS manifest"
echo "========================================"

# Parse qc_summary.tsv dynamically so column order is not assumed;
# the header row is used to build a lookup table mapping column names to indices
PASS_SAMPLES=$(awk -F'\t' '
    NR==1 {
        for(i=1; i<=NF; i++) header[$i] = i
        if(!("sample" in header) || !("status" in header)) {
            print "ERROR: qc_summary.tsv is missing required columns" > "/dev/stderr"
            exit 2
        }
        next
    }
    $header["status"] == "PASS" { print $header["sample"] }
' "$QC_SUMMARY")

PASS_COUNT=$(echo "$PASS_SAMPLES" | grep -c . || echo 0)
echo "[*] Found $PASS_COUNT PASS sample(s)"

# Initialise an empty manifest; it will be populated or remain empty if no PASS samples exist
: > "$MANIFEST_OUT"

FOUND=0; MISSING=0
mapfile -t PASS_SAMPLES_ARRAY <<< "$PASS_SAMPLES"

for sample in "${PASS_SAMPLES_ARRAY[@]}"; do
    [[ -z "$sample" ]] && continue    # skip blank lines from empty input
    bam="${PASS_DIR}/${sample}.bam"
    if [[ -f "$bam" ]]; then
        echo "$bam" >> "$MANIFEST_OUT"
        (( FOUND++ )) || true
    else
        echo "WARNING: BAM not found for PASS sample: $sample (expected: $bam)" >&2
        (( MISSING++ )) || true
    fi
done

echo "[*] Manifest written: $MANIFEST_OUT ($FOUND found, $MISSING missing)"

# =============================================================================
# Part 2: Zero- and low-coverage analysis per amplicon per sample
# =============================================================================
echo ""
echo "========================================"
echo "Part 2: Zero-coverage analysis"
echo "========================================"

# Collect BAMs from PASS and (if present) FAIL directories so that coverage
# gaps can be assessed across the full cohort, not just passing samples
SEARCH_DIRS=("$PASS_DIR")
[[ -d "$FAIL_DIR" ]] && SEARCH_DIRS+=("$FAIL_DIR")

mapfile -t ALL_BAMS < <(find "${SEARCH_DIRS[@]}" -type f -name "*.bam" 2>/dev/null | sort | grep -v '^$')
[[ ${#ALL_BAMS[@]} -gt 0 ]] || { echo "ERROR: no BAM files found in search directories" >&2; exit 1; }
echo "[*] Found ${#ALL_BAMS[@]} BAM(s) for analysis"

# Load target regions from the cleaned BED file (cols: chr, start, end, annotation)
mapfile -t BED_ROWS < <(awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4}' "$BED_FILE") || true
[[ ${#BED_ROWS[@]} -gt 0 ]] || { echo "ERROR: BED file is empty: $BED_FILE" >&2; exit 1; }
NUM_REGIONS=${#BED_ROWS[@]}
echo "[*] Target regions: $NUM_REGIONS"

# Build TSV header: fixed coordinate columns followed by four metrics per sample
HEADER="chr\tstart\tend\tannotation\tamplicon_length_bp"
for BAM in "${ALL_BAMS[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)
    # For each sample: absolute base counts and percentages for both thresholds
    HEADER="${HEADER}\t${SAMPLE}_zero_bases\t${SAMPLE}_pct_zero\t${SAMPLE}_low_bases\t${SAMPLE}_pct_low"
done

# Derive a cohort label from the directory structure (e.g. breast_tumour/dna_qc -> breast_tumour)
COHORT_NAME=$(basename "$(dirname "$QC_DIR")")
ZERO_COV_FILE="${ZERO_COV_OUT}/zero_cov_${COHORT_NAME}.tsv"
echo -e "$HEADER" > "$ZERO_COV_FILE"
echo "[*] Output file: $ZERO_COV_FILE"

# Temporary directory for per-base mosdepth outputs; retained after completion
# so that reruns can skip samples already processed
TMP_DIR="${ZERO_COV_OUT}/tmp_mosdepth"
mkdir -p "$TMP_DIR"

LOW_COV_THRESHOLD=200   # bases below this depth are flagged as low-coverage

# --- Step 2.1: Run per-base mosdepth for each sample -------------------------
echo ""
echo "[*] Step 1/2: Running per-base mosdepth..."

declare -A PERBASE_FILES
TOTAL_SAMPLES=${#ALL_BAMS[@]}; CURRENT=0

for BAM in "${ALL_BAMS[@]}"; do
    CURRENT=$((CURRENT + 1))
    SAMPLE=$(basename "$BAM" .bam)
    PREFIX="${TMP_DIR}/${SAMPLE}"
    PERBASE_OUTPUT="${PREFIX}.per-base.bed.gz"

    if [[ -f "$PERBASE_OUTPUT" ]]; then
        echo "  [$CURRENT/$TOTAL_SAMPLES] $SAMPLE — already computed, skipping"
        PERBASE_FILES["$SAMPLE"]="$PERBASE_OUTPUT"
        continue
    fi

    echo "  [$CURRENT/$TOTAL_SAMPLES] $SAMPLE — computing per-base depth..."
    # Unlike script 02 which uses --no-per-base for speed, here we need base-level
    # resolution to count exactly how many bases fall below each coverage threshold
    mosdepth --threads "$THREADS" \
             --fast-mode \
             "${PREFIX}" \
             "$BAM" > /dev/null 2>&1
    # Output: ${PREFIX}.per-base.bed.gz (run-length encoded: chr, start, end, depth)

    PERBASE_FILES["$SAMPLE"]="$PERBASE_OUTPUT"
done

echo "[*] Per-base coverage files ready"

# --- Step 2.2: Count zero- and low-coverage bases per amplicon ---------------
echo ""
echo "[*] Step 2/2: Counting zero- and sub-${LOW_COV_THRESHOLD}x bases per amplicon..."

REGION_NUM=0

for ROW in "${BED_ROWS[@]}"; do
    REGION_NUM=$((REGION_NUM + 1))
    IFS=$'\t' read -r chr start end annot <<< "$ROW"
    amplicon_len=$(( end - start ))

    (( REGION_NUM % 10 == 0 )) && \
        echo "  Region $REGION_NUM/$NUM_REGIONS (${chr}:${start}-${end}, ${amplicon_len} bp)..."

    line="${chr}\t${start}\t${end}\t${annot}\t${amplicon_len}"

    for BAM in "${ALL_BAMS[@]}"; do
        SAMPLE=$(basename "$BAM" .bam)
        PERBASE_FILE="${PERBASE_FILES[$SAMPLE]}"

        # mosdepth per-base output is run-length encoded: each record represents a
        # consecutive run of bases all at the same depth. To count bases accurately
        # we must compute the overlap of each run with the amplicon window and
        # multiply by depth — simply counting records would overcount short runs
        results=$(zcat "$PERBASE_FILE" | awk \
            -v chr_target="$chr" \
            -v amp_start="$start" \
            -v amp_end="$end" \
            -v amp_len="$amplicon_len" \
            -v low_thr="$LOW_COV_THRESHOLD" '
            BEGIN { zero_bases=0; low_bases=0 }

            # Restrict to records on this chromosome that overlap the amplicon window
            $1 == chr_target && $3 > amp_start && $2 < amp_end {
                # Clip run boundaries to the amplicon edges to avoid counting
                # bases outside the target region
                clip_start = ($2 > amp_start) ? $2 : amp_start
                clip_end   = ($3 < amp_end)   ? $3 : amp_end
                bases      = clip_end - clip_start

                if ($4 == 0)       zero_bases += bases   # col 4 = depth for this run
                if ($4 < low_thr)  low_bases  += bases
            }

            END {
                zero_pct = (amp_len > 0) ? (zero_bases / amp_len) * 100 : 0
                low_pct  = (amp_len > 0) ? (low_bases  / amp_len) * 100 : 0
                printf "%d\t%.2f\t%d\t%.2f\n", zero_bases, zero_pct, low_bases, low_pct
            }
        ')

        line="${line}\t${results}"
    done

    echo -e "$line" >> "$ZERO_COV_FILE"
done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "========================================"
echo "Analysis complete."
echo "  Manifest:           $MANIFEST_OUT ($FOUND PASS samples)"
echo "  Zero-coverage file: $ZERO_COV_FILE ($NUM_REGIONS regions, ${#ALL_BAMS[@]} samples)"
echo "  Temporary files:    $TMP_DIR (safe to delete once results are verified)"
echo "========================================"
