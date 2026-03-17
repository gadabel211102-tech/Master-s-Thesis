#!/usr/bin/env bash
# ==============================================================================
# Script: 02b_more-qc.sh 
# ==============================================================================
# Purpose:
#   Post-QC analysis combining two tasks:
#   1. Generate manifest files listing all samples that passed QC
#   2. Perform zero-coverage analysis to identify poorly covered bases per amplicon
#   
# This script should be run AFTER 02_dna_qc.sh completes
#
# ==============================================================================

# Exit immediately if any command fails, treat unset variables as errors,
# and propagate errors through pipes
set -euo pipefail

# Add error trap to show where failures occur
trap 'echo "[ERROR] Script failed at line $LINENO with exit code $?" >&2' ERR

# ==============================================================================
# CONFIGURATION & USAGE
# ==============================================================================

usage() {
    cat <<EOF
Usage: $0 -q QC_DIR -m MANIFEST_OUT -z ZERO_COV_OUT [-t THREADS]

Required Arguments:
  -q QC_DIR           Path to dna_qc directory (output from 02_dna_qc.sh)
                      Must contain: qc_summary.tsv, pass_bams/, fail_bams/, targets.sorted.bed
  
  -m MANIFEST_OUT     Output path for manifest file (list of PASS BAMs)
                      Example: /path/to/manifests/pass_breast_tumour.txt
  
  -z ZERO_COV_OUT     Output directory for zero-coverage analysis results
                      Example: /path/to/zero_cov_analysis/

Optional Arguments:
  -t THREADS          Number of CPU threads for mosdepth (default: 4)

Example:
  $0 -q /data/breast/tumour/dna_qc \\
     -m /data/manifests/pass_breast_tumour.txt \\
     -z /data/zero_cov_results \\
     -t 8

Output Files:
  1. Manifest file listing all PASS BAM paths (one per line)
  2. Zero-coverage TSV with per-amplicon, per-sample metrics:
       - Number and % of bases with 0x coverage
       - Number and % of bases with <200x coverage
EOF
}

# ==============================================================================
# PARSE COMMAND LINE ARGUMENTS
# ==============================================================================

QC_DIR=""
MANIFEST_OUT=""
ZERO_COV_OUT=""
THREADS="${THREADS:-4}"  # Default to 4 threads if not specified

while [[ $# -gt 0 ]]; do
    case "$1" in
        -q) QC_DIR="$2"; shift 2 ;;
        -m) MANIFEST_OUT="$2"; shift 2 ;;
        -z) ZERO_COV_OUT="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: Unknown option: $1"; usage; exit 2 ;;
    esac
done

# ==============================================================================
# VALIDATE REQUIRED INPUTS
# ==============================================================================

if [[ -z "$QC_DIR" || -z "$MANIFEST_OUT" || -z "$ZERO_COV_OUT" ]]; then
    echo "ERROR: Missing required arguments"
    usage
    exit 1
fi

# Define expected input files/directories from QC step
QC_SUMMARY="${QC_DIR}/qc_summary.tsv"
PASS_DIR="${QC_DIR}/pass_bams"
FAIL_DIR="${QC_DIR}/fail_bams"
BED_FILE="${QC_DIR}/targets.sorted.bed"

# Verify all required inputs exist
if [[ ! -f "$QC_SUMMARY" ]]; then
    echo "ERROR: qc_summary.tsv not found: $QC_SUMMARY"
    exit 1
fi

if [[ ! -d "$PASS_DIR" ]]; then
    echo "ERROR: pass_bams/ directory not found: $PASS_DIR"
    exit 1
fi

# FAIL_DIR is optional - warn if missing but don't exit
if [[ ! -d "$FAIL_DIR" ]]; then
    echo "WARNING: fail_bams/ directory not found: $FAIL_DIR"
    echo "         Will analyze PASS samples only"
fi

if [[ ! -f "$BED_FILE" ]]; then
    echo "ERROR: targets.sorted.bed not found: $BED_FILE"
    exit 1
fi

# Create output directories
mkdir -p "$(dirname "$MANIFEST_OUT")"
mkdir -p "$ZERO_COV_OUT"

# ==============================================================================
# PART 1: GENERATE MANIFEST FILE FOR PASS SAMPLES
# ==============================================================================
echo ""
echo "========================================"
echo "PART 1: GENERATING PASS MANIFEST"
echo "========================================"
echo ""

echo "[INFO] Reading QC summary: $QC_SUMMARY"
echo "[INFO] Looking for PASS samples..."

# --------------------------------------------------
# Extract sample names that passed QC
# --------------------------------------------------
# Strategy:
# 1. Parse header to find column indices for "sample" and "status"
# 2. For data rows, extract sample name where status == "PASS"

PASS_SAMPLES=$(awk -F'\t' '
  # Process header row
  NR==1 {
    # Build hash table mapping column names to indices
    for(i=1; i<=NF; i++) {
      header[$i] = i
    }
    
    # Verify required columns exist
    if(!("sample" in header) || !("status" in header)) {
      print "ERROR: Missing sample/status columns in qc_summary.tsv" >"/dev/stderr"
      exit 2
    }
    next  # Skip to next line (start processing data)
  }
  
  # Process data rows
  # Extract sample name if status is PASS
  $header["status"] == "PASS" {
    print $header["sample"]
  }
' "$QC_SUMMARY")

# Count how many PASS samples were found
PASS_COUNT=$(echo "$PASS_SAMPLES" | grep -c . || echo 0)

if [[ $PASS_COUNT -eq 0 ]]; then
    echo "[WARN] No PASS samples found in QC summary"
    echo "[INFO] Creating empty manifest"
    : > "$MANIFEST_OUT"  # Create empty file
    echo "[INFO] Manifest written: $MANIFEST_OUT (0 samples)"
else
    echo "[INFO] Found $PASS_COUNT PASS samples"
fi

# --------------------------------------------------
# Build manifest file with full BAM paths
# --------------------------------------------------
echo "[INFO] Building manifest file..."
echo "[INFO] Output path: $MANIFEST_OUT"

# Initialize empty manifest file
: > "$MANIFEST_OUT"

# Counter for tracking
FOUND=0
MISSING=0

# Handle case where PASS_SAMPLES is empty
if [[ -z "$PASS_SAMPLES" ]]; then
    echo "[DEBUG] PASS_SAMPLES variable is empty"
else
    echo "[DEBUG] PASS_SAMPLES contains:"
    echo "$PASS_SAMPLES" | head -5
fi

# Convert PASS_SAMPLES string to array for proper iteration
mapfile -t PASS_SAMPLES_ARRAY <<< "$PASS_SAMPLES"

echo "[DEBUG] Array has ${#PASS_SAMPLES_ARRAY[@]} elements"

# For each PASS sample, verify BAM exists and add to manifest
for sample in "${PASS_SAMPLES_ARRAY[@]}"; do
    # Skip empty lines
    if [[ -z "$sample" ]]; then
        echo "[DEBUG] Skipping empty sample entry"
        continue
    fi
    
    bam="${PASS_DIR}/${sample}.bam"
    
    echo "[DEBUG] Checking BAM: $bam"
    
    if [[ -f "$bam" ]]; then
        # BAM file exists - add to manifest
        echo "$bam" >> "$MANIFEST_OUT"
        FOUND=$((FOUND + 1))
        echo "[DEBUG] Found BAM $FOUND: $sample"
    else
        # BAM file missing - warn but continue
        echo "[WARN] Missing BAM for PASS sample: $sample (expected: $bam)" >&2
        MISSING=$((MISSING + 1))
    fi
done

echo "[DEBUG] Loop completed. FOUND=$FOUND, MISSING=$MISSING"

echo ""
echo "[SUCCESS] Manifest generation completed"
echo "  Total PASS samples:    $PASS_COUNT"
echo "  BAMs found:            $FOUND"
echo "  BAMs missing:          $MISSING"
echo "  Manifest file:         $MANIFEST_OUT"
echo ""
echo "[DEBUG] *** PART 1 COMPLETE - About to start Part 2 ***"
echo "[DEBUG] set -e is: $(set +o | grep errexit)"
echo "[DEBUG] set -u is: $(set +o | grep nounset)"
echo "[DEBUG] set -o pipefail is: $(set +o | grep pipefail)"
echo "[DEBUG] Part 1 completed successfully, proceeding to Part 2..."

# ==============================================================================
# PART 2: ZERO-COVERAGE ANALYSIS
# ==============================================================================
echo ""
echo "========================================"
echo "PART 2: ZERO-COVERAGE ANALYSIS"
echo "========================================"
echo ""
echo "[DEBUG] *** ENTERING PART 2 ***"
echo "[DEBUG] Current directory: $(pwd)"
echo "[DEBUG] Script still running with PID: $$"
echo "[DEBUG] Starting Part 2 - zero coverage analysis"

# 1. Verify PASS directory exists
if [[ ! -d "$PASS_DIR" ]]; then
    echo "[ERROR] Pass directory missing: $PASS_DIR"
    exit 1
fi

# 2. Collect BAMs from directories that exist
echo "[INFO] Collecting BAM files from PASS directory..."
echo "[DEBUG] PASS_DIR: $PASS_DIR"

# First check what's in PASS directory
echo "[DEBUG] Contents of PASS_DIR:"
ls -la "$PASS_DIR" | head -20 || echo "[WARN] Could not list PASS_DIR"

# Build list of directories to search
SEARCH_DIRS=("$PASS_DIR")

# Add FAIL_DIR if it exists
if [[ -d "$FAIL_DIR" ]]; then
    echo "[INFO] FAIL directory found, including it in analysis"
    echo "[DEBUG] FAIL_DIR: $FAIL_DIR"
    echo "[DEBUG] Contents of FAIL_DIR:"
    ls -la "$FAIL_DIR" | head -20 || echo "[WARN] Could not list FAIL_DIR"
    SEARCH_DIRS+=("$FAIL_DIR")
else
    echo "[WARN] FAIL directory not found (will analyze PASS samples only): $FAIL_DIR"
fi

# Collect BAMs from all search directories
BAM_LIST=$(find "${SEARCH_DIRS[@]}" -type f -name "*.bam" 2>/dev/null | sort) || true

echo "[DEBUG] Number of BAMs found by find: $(echo "$BAM_LIST" | grep -c . || echo 0)"

if [[ -z "$BAM_LIST" ]]; then
    echo "[ERROR] No BAM files found in search directories:"
    for dir in "${SEARCH_DIRS[@]}"; do
        echo "  - $dir"
    done
    exit 1
fi

echo "[DEBUG] BAM_LIST contents (first 5 lines):"
echo "$BAM_LIST" | head -5

# Filter out empty lines and create array
mapfile -t ALL_BAMS < <(echo "$BAM_LIST" | grep -v '^$')

echo "[DEBUG] Array size after mapfile: ${#ALL_BAMS[@]}"
echo "[DEBUG] First BAM in array: ${ALL_BAMS[0]}"

# 3. Create temporary directory
TMP_DIR="${ZERO_COV_OUT}/tmp_mosdepth"
mkdir -p "$TMP_DIR"

echo "[INFO] Found ${#ALL_BAMS[@]} total samples for analysis"
echo "[INFO] Using BED file: $BED_FILE"

# 4. Read BED file into array
echo "[INFO] Reading target regions from BED file..."
mapfile -t BED_ROWS < <(awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4}' "$BED_FILE") || true

if [[ ${#BED_ROWS[@]} -eq 0 ]]; then
    echo "[ERROR] BED file is empty or formatted incorrectly: $BED_FILE"
    exit 1
fi

NUM_REGIONS=${#BED_ROWS[@]}
echo "[INFO] Found $NUM_REGIONS target regions"

# --------------------------------------------------
# Build output file header
# --------------------------------------------------

# Start with genomic coordinates, annotation, and amplicon length
HEADER="chr\tstart\tend\tannotation\tamplicon_length_bp"

# Add four columns per sample:
#   <sample>_zero_bases   - number of bases with 0x coverage
#   <sample>_pct_zero     - percentage of amplicon with 0x coverage
#   <sample>_low_bases    - number of bases with <200x coverage (includes zero)
#   <sample>_pct_low      - percentage of amplicon with <200x coverage
for BAM in "${ALL_BAMS[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)
    HEADER="${HEADER}\t${SAMPLE}_zero_bases\t${SAMPLE}_pct_zero\t${SAMPLE}_low_bases\t${SAMPLE}_pct_low"
done

# Determine cohort name from directory structure
# Format: parent_directory_current_directory
# Example: /data/breast/tumour/dna_qc â†’ breast_tumour
echo "[DEBUG] QC_DIR: $QC_DIR"
echo "[DEBUG] dirname QC_DIR: $(dirname "$QC_DIR")"
echo "[DEBUG] basename dirname QC_DIR: $(basename "$(dirname "$QC_DIR")")"

PARENT_DIR=$(basename "$(dirname "$QC_DIR")")
COHORT_NAME="${PARENT_DIR}"

echo "[DEBUG] COHORT_NAME: $COHORT_NAME"

# Create output file
ZERO_COV_FILE="${ZERO_COV_OUT}/zero_cov_${COHORT_NAME}.tsv"
echo -e "$HEADER" > "$ZERO_COV_FILE"

echo "[INFO] Output file: $ZERO_COV_FILE"
echo "[INFO] Temporary mosdepth directory: $TMP_DIR"
echo ""

# ==============================================================================
# STEP 2.1: Run per-base mosdepth for each sample (one-time computation)
# ==============================================================================
echo "[INFO] STEP 1/2: Running per-base coverage analysis with mosdepth..."
echo ""

# Associative array to store mosdepth output paths
declare -A PERBASE_FILES

TOTAL_SAMPLES=${#ALL_BAMS[@]}
CURRENT=0

for BAM in "${ALL_BAMS[@]}"; do
    CURRENT=$((CURRENT + 1))
    SAMPLE=$(basename "$BAM" .bam)
    PREFIX="${TMP_DIR}/${SAMPLE}"
    PERBASE_OUTPUT="${PREFIX}.per-base.bed.gz"
    
    # Check if this sample's per-base file already exists (resume capability)
    if [[ -f "$PERBASE_OUTPUT" ]]; then
        echo "  [$CURRENT/$TOTAL_SAMPLES] $SAMPLE - Already computed (skipping)"
        PERBASE_FILES["$SAMPLE"]="$PERBASE_OUTPUT"
        continue
    fi
    
    echo "  [$CURRENT/$TOTAL_SAMPLES] $SAMPLE - Computing per-base coverage..."
    
    # Run mosdepth in per-base mode
    # --threads: Use multiple CPUs
    # --fast-mode: Skip some internal checks for speed
    # Output: ${PREFIX}.per-base.bed.gz with format: chr, start, end, coverage
    mosdepth --threads "$THREADS" \
             --fast-mode \
             "${PREFIX}" \
             "$BAM" > /dev/null 2>&1
    
    # Store the output file path for later use
    PERBASE_FILES["$SAMPLE"]="$PERBASE_OUTPUT"
done

echo ""
echo "[SUCCESS] All per-base coverage files ready"
echo ""

# ==============================================================================
# STEP 2.2: Count zero-coverage and low-coverage bases per amplicon per sample
# ==============================================================================
echo "[INFO] STEP 2/2: Counting zero-coverage and sub-200x bases per amplicon..."
echo ""

# Coverage threshold for low-coverage metric
LOW_COV_THRESHOLD=200

# Progress tracking
REGION_NUM=0

# Iterate over each target region (amplicon)
for ROW in "${BED_ROWS[@]}"; do
    REGION_NUM=$((REGION_NUM + 1))
    
    # Parse BED row into variables
    IFS=$'\t' read -r chr start end annot <<< "$ROW"
    
    # Calculate amplicon length in bases
    amplicon_len=$(( end - start ))
    
    # Progress indicator every 10 regions
    if (( REGION_NUM % 10 == 0 )); then
        echo "  Processing region $REGION_NUM/$NUM_REGIONS (${chr}:${start}-${end}, ${amplicon_len}bp)..."
    fi
    
    # Start building the output line with genomic coordinates and amplicon length
    line="${chr}\t${start}\t${end}\t${annot}\t${amplicon_len}"
    
    # For each sample, count bases with 0x coverage and bases below 200x in this region
    for BAM in "${ALL_BAMS[@]}"; do
        SAMPLE=$(basename "$BAM" .bam)
        PERBASE_FILE="${PERBASE_FILES[$SAMPLE]}"
        
        # Count zero-coverage and low-coverage bases using AWK
        #
        # mosdepth per-base format (run-length encoded):
        #   col1: chr
        #   col2: start (0-based)
        #   col3: end   (0-based, exclusive)
        #   col4: coverage depth
        #
        # Each record represents a RUN of consecutive bases all at the same depth.
        # To get true base counts we must sum (overlap_end - overlap_start) for
        # each run, clipped to the amplicon boundaries, NOT simply count records.
        #
        # We clip each run to the amplicon region before counting to handle
        # runs that partially overlap the amplicon boundary.
        results=$(zcat "$PERBASE_FILE" | awk \
            -v chr_target="$chr" \
            -v amp_start="$start" \
            -v amp_end="$end" \
            -v amp_len="$amplicon_len" \
            -v low_thr="$LOW_COV_THRESHOLD" '
            BEGIN {
                zero_bases = 0
                low_bases  = 0
            }
            # Only process records on the correct chromosome that overlap this amplicon
            $1 == chr_target && $3 > amp_start && $2 < amp_end {

                # Clip the run to the amplicon boundaries
                # (handles partial overlaps at amplicon edges)
                clip_start = ($2 > amp_start) ? $2 : amp_start
                clip_end   = ($3 < amp_end)   ? $3 : amp_end
                bases      = clip_end - clip_start  # number of bases in this clipped run

                # Accumulate zero-coverage bases
                if ($4 == 0) {
                    zero_bases += bases
                }

                # Accumulate sub-200x bases (includes zero-coverage bases)
                if ($4 < low_thr) {
                    low_bases += bases
                }
            }
            END {
                # Calculate percentages (guard against zero-length amplicons)
                if (amp_len > 0) {
                    zero_pct = (zero_bases / amp_len) * 100
                    low_pct  = (low_bases  / amp_len) * 100
                } else {
                    zero_pct = 0
                    low_pct  = 0
                }
                # Output: zero_bases, pct_zero, low_bases, pct_low
                printf "%d\t%.2f\t%d\t%.2f\n", zero_bases, zero_pct, low_bases, low_pct
            }
        ')
        
        # Append metrics to the line
        line="${line}\t${results}"
    done
    
    # Write completed line to output file
    echo -e "$line" >> "$ZERO_COV_FILE"
done

echo ""
echo "[SUCCESS] Zero-coverage analysis completed"
echo ""

# ==============================================================================
# CLEANUP & FINAL SUMMARY
# ==============================================================================

echo ""
echo "========================================"
echo "ANALYSIS COMPLETE - SUMMARY"
echo "========================================"
echo ""
echo "PART 1 - Manifest:"
echo "  File:         $MANIFEST_OUT"
echo "  PASS samples: $FOUND"
echo ""
echo "PART 2 - Zero Coverage:"
echo "  File:         $ZERO_COV_FILE"
echo "  Regions:      $NUM_REGIONS"
echo "  Samples:      ${#ALL_BAMS[@]}"
echo "  Metrics/sample: zero_bases, pct_zero, low_bases (<${LOW_COV_THRESHOLD}x), pct_low"
echo ""
echo "Temporary per-base files stored in: $TMP_DIR"
echo "  (You can safely delete this directory to save space)"
echo ""
echo "Next steps:"
echo "  1. Review zero-coverage report for problematic amplicons"
echo "  2. Use manifest file for downstream variant calling"
echo "  3. Consider re-designing amplicons with high zero-coverage rates"
echo ""