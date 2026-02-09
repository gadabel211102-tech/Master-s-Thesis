#!/usr/bin/env bash
# ==============================================================================
# Script: 02_post_qc_analysis.sh (UNIFIED & FULLY ANNOTATED)
# ==============================================================================
# Purpose:
#   Post-QC analysis combining two tasks:
#   1. Generate manifest files listing all samples that passed QC
#   2. Perform zero-coverage analysis to identify poorly covered bases per amplicon
#
# This script should be run AFTER 02_dna_qc.sh completes
#
# Author: [Your name]
# Date: February 2026
# ==============================================================================

# Exit immediately if any command fails, treat unset variables as errors,
# and propagate errors through pipes
set -euo pipefail

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
  2. Zero-coverage TSV showing bases with 0x coverage per amplicon per sample
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

if [[ ! -d "$FAIL_DIR" ]]; then
    echo "ERROR: fail_bams/ directory not found: $FAIL_DIR"
    exit 1
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

# For each PASS sample, verify BAM exists and add to manifest
for sample in $PASS_SAMPLES; do
    bam="${PASS_DIR}/${sample}.bam"
    
    if [[ -f "$bam" ]]; then
        # BAM file exists - add to manifest
        echo "$bam" >> "$MANIFEST_OUT"
        ((FOUND++))
    else
        # BAM file missing - warn but continue
        echo "[WARN] Missing BAM for PASS sample: $sample (expected: $bam)" >&2
        ((MISSING++))
    fi
done

echo ""
echo "[SUCCESS] Manifest generation completed"
echo "  Total PASS samples:    $PASS_COUNT"
echo "  BAMs found:            $FOUND"
echo "  BAMs missing:          $MISSING"
echo "  Manifest file:         $MANIFEST_OUT"
echo ""

# ==============================================================================
# PART 2: ZERO-COVERAGE ANALYSIS
# ==============================================================================
echo ""
echo "========================================"
echo "PART 2: ZERO-COVERAGE ANALYSIS"
echo "========================================"
echo ""

# --------------------------------------------------
# Setup for zero-coverage analysis
# --------------------------------------------------

# Create temporary directory for per-base mosdepth outputs
TMP_DIR="${ZERO_COV_OUT}/tmp_mosdepth"
mkdir -p "$TMP_DIR"

echo "[INFO] Collecting BAM files from PASS and FAIL directories..."

# Collect all BAM files (both PASS and FAIL) for comprehensive analysis
mapfile -t ALL_BAMS < <(find "$PASS_DIR" "$FAIL_DIR" -type f -name "*.bam" | sort)

# Verify we found BAMs
if [[ ${#ALL_BAMS[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found in PASS or FAIL directories"
    exit 1
fi

echo "[INFO] Found ${#ALL_BAMS[@]} total samples (PASS + FAIL)"
echo "[INFO] Using BED file: $BED_FILE"

# --------------------------------------------------
# Read BED file into array
# --------------------------------------------------
# Extract only chr, start, end, annotation (first 4 columns)
echo "[INFO] Reading target regions from BED file..."

mapfile -t BED_ROWS < <(awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4}' "$BED_FILE")

NUM_REGIONS=${#BED_ROWS[@]}
echo "[INFO] Found $NUM_REGIONS target regions"

# --------------------------------------------------
# Build output file header
# --------------------------------------------------

# Start with genomic coordinates and annotation
HEADER="chr\tstart\tend\tannotation"

# Add each sample name as a column
for BAM in "${ALL_BAMS[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)
    HEADER="${HEADER}\t${SAMPLE}"
done

# Determine cohort name from directory structure
# Format: parent_directory_current_directory
# Example: /data/breast/tumour/dna_qc → breast_tumour
COHORT_NAME=$(basename "$(dirname "$QC_DIR")")_$(basename "$(dirname "$QC_DIR/dummy")")

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
    ((CURRENT++))
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
# STEP 2.2: Count zero-coverage bases per amplicon per sample
# ==============================================================================
echo "[INFO] STEP 2/2: Counting zero-coverage bases per amplicon..."
echo ""

# Progress tracking
REGION_NUM=0

# Iterate over each target region (amplicon)
for ROW in "${BED_ROWS[@]}"; do
    ((REGION_NUM++))
    
    # Parse BED row into variables
    IFS=$'\t' read -r chr start end annot <<< "$ROW"
    
    # Progress indicator every 10 regions
    if (( REGION_NUM % 10 == 0 )); then
        echo "  Processing region $REGION_NUM/$NUM_REGIONS (${chr}:${start}-${end})..."
    fi
    
    # Start building the output line with genomic coordinates
    line="${chr}\t${start}\t${end}\t${annot}"
    
    # For each sample, count bases with 0x coverage in this region
    for BAM in "${ALL_BAMS[@]}"; do
        SAMPLE=$(basename "$BAM" .bam)
        PERBASE_FILE="${PERBASE_FILES[$SAMPLE]}"
        
        # Count zero-coverage bases using AWK
        # Logic:
        # - Extract only records matching this amplicon (chr and overlapping coordinates)
        # - Count records where coverage (column 4) == 0
        # - Return count (default to 0 if no matches)
        zero_count=$(zcat "$PERBASE_FILE" | awk \
            -v chr_target="$chr" \
            -v start_target="$start" \
            -v end_target="$end" '
            BEGIN {
                count = 0
            }
            # Check if this base falls within our target region
            $1 == chr_target && $2 >= start_target && $3 <= end_target && $4 == 0 {
                count++
            }
            END {
                print count + 0  # Ensure numeric output (0 if no matches)
            }
        ')
        
        # Append zero-count to the line
        line="${line}\t${zero_count}"
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
echo ""
echo "Temporary per-base files stored in: $TMP_DIR"
echo "  (You can safely delete this directory to save space)"
echo ""
echo "Next steps:"
echo "  1. Review zero-coverage report for problematic amplicons"
echo "  2. Use manifest file for downstream variant calling"
echo "  3. Consider re-designing amplicons with high zero-coverage rates"
echo ""
