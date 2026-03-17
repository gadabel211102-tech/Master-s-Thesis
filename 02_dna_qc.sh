#!/usr/bin/env bash
# ==============================================================================
# Script: 02_dna_qc.sh (FULLY ANNOTATED VERSION)
# ==============================================================================
# Purpose:
#   Quality control pipeline for targeted DNA sequencing (Ion Torrent)
#   - Validates and indexes BAM files
#   - Calculates coverage metrics using mosdepth
#   - Identifies poorly performing amplicons
#   - Filters samples into PASS/FAIL categories
# ==============================================================================

# Exit immediately if any command fails, treat unset variables as errors,
# and propagate errors through pipes
set -euo pipefail

# ==============================================================================
# CONFIGURATION - QC THRESHOLDS
# ==============================================================================
# These thresholds determine whether a sample passes or fails QC

MIN_TOTAL_READS=120000      # Minimum total reads required per sample
MIN_MAPPED_PCT=90           # Minimum % of reads that map to reference
MIN_ON_TARGET_PCT=75        # Minimum % of mapped reads that hit target regions
MIN_MEAN_COVERAGE=200       # Minimum mean depth across all target regions
UNIFORMITY_THRESHOLD=100    # Depth threshold for uniformity calculation
MIN_UNIFORMITY_PCT=95       # Minimum % of amplicons above uniformity threshold

# ==============================================================================
# INPUT PARAMETERS - Command Line Arguments
# ==============================================================================
INPUT_DIR=""    # Directory containing input BAM files
OUT_DIR=""      # Output directory for QC results
REF_FA=""       # Reference genome FASTA file
BED=""          # BED file defining target regions (amplicons)
THREADS="${THREADS:-4}"  # Number of CPU threads (default: 4, can be set via environment)

# ==============================================================================
# FUNCTION: Display usage information
# ==============================================================================
usage() {
    cat <<EOF
Usage: $0 -i INPUT_DIR -o OUT_DIR -r REF_FASTA -b TARGETS_BED [-t THREADS]

Required Arguments:
  -i INPUT_DIR    Directory containing input BAM files
  -o OUT_DIR      Output directory for QC results
  -r REF_FASTA    Reference genome FASTA file (e.g., hg38.fa)
  -b TARGETS_BED  BED file with target regions (amplicon coordinates)

Optional Arguments:
  -t THREADS      Number of CPU threads to use (default: 4)

Example:
  $0 -i /data/bams -o /data/qc_results -r /ref/hg38.fa -b /ref/targets.bed -t 8
EOF
}

# ==============================================================================
# PARSE COMMAND LINE ARGUMENTS
# ==============================================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT_DIR="$2"; shift 2 ;;   # Capture input directory
        -o) OUT_DIR="$2"; shift 2 ;;     # Capture output directory
        -r) REF_FA="$2"; shift 2 ;;      # Capture reference FASTA
        -b) BED="$2"; shift 2 ;;         # Capture BED file
        -t) THREADS="$2"; shift 2 ;;     # Capture thread count
        *) usage; exit 2 ;;              # Unknown option - show help and exit
    esac
done

# ==============================================================================
# VALIDATE REQUIRED INPUTS
# ==============================================================================
# Check that all required parameters were provided
if [[ -z "$INPUT_DIR" || -z "$OUT_DIR" || -z "$REF_FA" || -z "$BED" ]]; then
    echo "ERROR: Missing required arguments"
    usage
    exit 1
fi

# ==============================================================================
# CREATE OUTPUT DIRECTORY STRUCTURE
# ==============================================================================
# Create subdirectories for organizing results:
#   - pass_bams/  : BAMs that passed QC
#   - fail_bams/  : BAMs that failed QC
#   - stats/      : Per-sample coverage statistics
mkdir -p "${OUT_DIR}/pass_bams" "${OUT_DIR}/fail_bams" "${OUT_DIR}/stats"

# ==============================================================================
# STEP 1: BED FILE CLEANING AND NORMALIZATION
# ==============================================================================
echo "[*] Step 1: Cleaning and normalizing BED file..."

# Output path for the cleaned BED file
BED_CALC="${OUT_DIR}/targets.sorted.bed"

# BED file cleaning logic:
# 1. Remove carriage returns (\r) that cause parsing issues
# 2. Skip 'track' and 'browser' header lines
# 3. Keep only lines starting with valid chromosome names (chr1-22, X, Y, M or 1-22)
# 4. Require at least 3 columns (chr, start, end)
# 5. Normalize chromosome names to include 'chr' prefix
# 6. Sort by chromosome and start position for efficient processing
# **CRITICAL FIX**: Preserve ALL columns (including rsID in column 5 and gene in column 6)
awk 'BEGIN{FS="\t"; OFS="\t"} {
    # Remove Windows-style line endings
    gsub(/\r/,""); 
    
    # Skip track and browser lines
    if ($1 ~ /^(track|browser)/) next;
    
    # Check if line starts with valid chromosome identifier
    if ($1 ~ /^(chr|[0-9XYM])/) {
        # Ensure minimum 3 columns present
        if (NF >= 3) {
            # Add "chr" prefix if missing (e.g., "1" -> "chr1")
            if ($1 !~ /^chr/) $1="chr"$1;
            # Output ALL columns to preserve rsID and gene annotations
            print $0
        }
    }
}' "$BED" | sort -k1,1 -k2,2n > "$BED_CALC"

# Verify the cleaned BED file is not empty
if [[ ! -s "$BED_CALC" ]]; then
    echo "ERROR: Cleaned BED file is empty. Check input BED format."
    exit 1
fi

echo "[*] BED file cleaned and sorted: $BED_CALC"

# ==============================================================================
# STEP 2: DISCOVER INPUT BAM FILES
# ==============================================================================
echo "[*] Step 2: Discovering BAM files..."

# Find all BAM files in input directory (max 2 levels deep) and sort them
mapfile -t BAMS < <(find "$INPUT_DIR" -maxdepth 2 -type f -name "*.bam" | sort)
NUM_BAMS=${#BAMS[@]}

# Verify we found at least one BAM file
if [[ $NUM_BAMS -eq 0 ]]; then
    echo "ERROR: No BAM files found in $INPUT_DIR"
    exit 1
fi

echo "[*] Found $NUM_BAMS BAM files to process"

# ==============================================================================
# STEP 3: INITIALIZE SUMMARY TABLE
# ==============================================================================
# Create TSV file with header for QC summary report
SUMMARY="${OUT_DIR}/qc_summary.tsv"
echo -e "sample\ttotal_reads\tmapped_pct\ton_target_pct\tmean_cov\tuniformity_100x\tworst_amplicon\tstatus\tfail_reasons" > "$SUMMARY"

# Column descriptions:
#   sample          - Sample name (from BAM filename)
#   total_reads     - Total number of reads in BAM
#   mapped_pct      - Percentage of reads mapped to reference
#   on_target_pct   - Percentage of mapped reads hitting target regions
#   mean_cov        - Weighted mean coverage across all amplicons
#   uniformity_100x - Percentage of amplicons with >=100x coverage
#   worst_amplicon  - Amplicon with lowest coverage (ID:depth)
#   status          - PASS or FAIL
#   fail_reasons    - Specific QC failures (if any)

# ==============================================================================
# STEP 4: PROCESS EACH BAM FILE
# ==============================================================================
echo "[*] Step 3: Processing samples..."

for (( i=0; i<NUM_BAMS; i++ )); do
    # --------------------------------------------------
    # 4.1: Setup for current sample
    # --------------------------------------------------
    BAM="${BAMS[$i]}"                           # Full path to BAM file
    SAMPLE=$(basename "$BAM" .bam)              # Extract sample name (remove .bam extension)
    SDIR="${OUT_DIR}/stats/${SAMPLE}"           # Sample-specific stats directory
    
    mkdir -p "$SDIR"  # Create output directory for this sample

    echo ""
    echo "================================================"
    echo "Processing sample $((i+1))/$NUM_BAMS: $SAMPLE"
    echo "================================================"

    # --------------------------------------------------
    # 4.2: Calculate basic BAM statistics
    # --------------------------------------------------
    echo "  [1/7] Calculating BAM statistics..."
    
    # Run samtools flagstat to get read counts
    FLAGSTAT=$(samtools flagstat -@ "$THREADS" "$BAM")
    
    # Extract total reads (first line of flagstat output)
    TOTAL_READS=$(echo "$FLAGSTAT" | awk '/in total/ {print $1}')
    
    # Extract mapped reads (exclude "primary mapped" line to avoid double-counting)
    MAPPED_READS=$(echo "$FLAGSTAT" | awk '/ mapped / && !/primary/ {print $1}')
    
    # Calculate mapping percentage
    # Uses awk for floating-point arithmetic: (mapped/total) * 100
    MAPPED_PCT=$(awk -v m="$MAPPED_READS" -v t="$TOTAL_READS" \
        'BEGIN{if(t>0) printf "%.2f", (m/t)*100; else print "0"}')

    # --------------------------------------------------
    # 4.3: Calculate on-target metrics
    # --------------------------------------------------
    echo "  [2/7] Calculating on-target percentage..."
    
    # Count reads overlapping target regions
    # -F 4: exclude unmapped reads
    # -c: count only (don't output reads)
    # -L: restrict to regions in BED file
    ON_TARGET_READS=$(samtools view -@ "$THREADS" -F 4 -c -L "$BED_CALC" "$BAM")
    
    # Calculate on-target percentage: (on-target / mapped) * 100
    ON_TARGET_PCT=$(awk -v o="$ON_TARGET_READS" -v m="$MAPPED_READS" \
        'BEGIN{if(m>0) printf "%.2f", (o/m)*100; else print "0"}')

    # --------------------------------------------------
    # 4.4: Run mosdepth for per-region coverage
    # --------------------------------------------------
    echo "  [3/7] Running mosdepth coverage analysis..."
    
    # mosdepth flags:
    #   --threads       : Use multiple CPUs
    #   --no-per-base   : Don't calculate per-base depth (saves time/space)
    #   --by            : Calculate coverage for regions in BED file
    #   --fast-mode     : Skip per-base calculations (faster)
    mosdepth --threads "$THREADS" \
             --no-per-base \
             --by "$BED_CALC" \
             --fast-mode \
             "${SDIR}/${SAMPLE}" \
             "$BAM"
    
    # This creates several output files:
    #   ${SAMPLE}.regions.bed.gz       - Per-region coverage (chr, start, end, region_id, coverage)
    #   ${SAMPLE}.mosdepth.summary.txt - Overall statistics
    #   ${SAMPLE}.mosdepth.region.dist.txt - Coverage distribution

    # --------------------------------------------------
    # 4.5: Calculate mean coverage (weighted by region length)
    # --------------------------------------------------
    echo "  [4/7] Calculating mean coverage..."
    
    REG="${SDIR}/${SAMPLE}.regions.bed.gz"
    
    # Mosdepth regions.bed.gz format:
    # Column 1: chromosome
    # Column 2: start position
    # Column 3: end position
    # Column 4: region ID (often numeric)
    # Column 5: mean coverage for this region *** THIS IS WHAT WE NEED ***
    
    # Calculate weighted mean:
    # For each region: multiply (length * coverage), sum all, divide by total bases
    # This gives more weight to larger regions, which is biologically appropriate
    MEAN_COV=$(zcat "$REG" | awk '
        BEGIN {
            sum=0      # Running sum of (length * coverage)
            bp=0       # Total bases covered
        }
        {
            len = $3 - $2              # Calculate region length
            sum += len * $5            # Add (length * coverage) to sum
            bp += len                  # Add length to total bases
        }
        END {
            if(bp > 0) 
                printf "%.2f", sum/bp  # Weighted mean = total_coverage / total_bases
            else 
                print "0"
        }
    ')

    # --------------------------------------------------
    # 4.6: Calculate uniformity (% of amplicons above threshold)
    # --------------------------------------------------
    echo "  [5/7] Calculating coverage uniformity..."
    
    # Count how many amplicons have coverage >= UNIFORMITY_THRESHOLD (100x)
    # This indicates how evenly covered the target regions are
    UNIFORMITY=$(zcat "$REG" | awk -v thr="$UNIFORMITY_THRESHOLD" '
        BEGIN {
            pass_count=0   # Amplicons passing threshold
            total=0        # Total amplicons
        }
        {
            total++
            if($5 >= thr)  # Check if coverage (col 5) >= threshold
                pass_count++
        }
        END {
            if(total > 0) 
                printf "%.2f", (pass_count/total)*100
            else 
                print "0"
        }
    ')

    # --------------------------------------------------
    # 4.7: Identify worst-performing amplicon
    # --------------------------------------------------
    echo "  [6/7] Identifying worst amplicon..."
    
    # Extract region_id and coverage, sort by coverage (ascending), take first entry
    # Format: region_id:coverage (e.g., "557279:123.45")
    WORST_AMP=$(zcat "$REG" | awk '{
            # Print: region_id:coverage (formatted to 2 decimal places)
            print $4 ":" sprintf("%.2f", $5)
        }' | \
        sort -t: -k2,2n |  # Sort by coverage (numeric, after the colon)
        head -1)            # Take the lowest value
    
    # Handle edge case where no regions were found
    [[ -z "$WORST_AMP" ]] && WORST_AMP="N/A"

    # --------------------------------------------------
    # 4.8: Determine PASS/FAIL status and reasons
    # --------------------------------------------------
    echo "  [7/7] Evaluating QC status..."
    
    STATUS="PASS"        # Start optimistic
    FAIL_REASONS=""      # Empty reasons list

    # Check each QC threshold and build a list of failures
    # bc -l: use bc calculator for floating-point comparison
    
    # Threshold 1: Total reads
    if (( $(echo "$TOTAL_READS < $MIN_TOTAL_READS" | bc -l) )); then
        STATUS="FAIL"
        FAIL_REASONS="${FAIL_REASONS}Low total reads; "
    fi

    # Threshold 2: Mapping rate
    if (( $(echo "$MAPPED_PCT < $MIN_MAPPED_PCT" | bc -l) )); then
        STATUS="FAIL"
        FAIL_REASONS="${FAIL_REASONS}Low mapping rate; "
    fi

    # Threshold 3: On-target percentage
    if (( $(echo "$ON_TARGET_PCT < $MIN_ON_TARGET_PCT" | bc -l) )); then
        STATUS="FAIL"
        FAIL_REASONS="${FAIL_REASONS}Low on target; "
    fi

    # Threshold 4: Mean coverage
    if (( $(echo "$MEAN_COV < $MIN_MEAN_COVERAGE" | bc -l) )); then
        STATUS="FAIL"
        FAIL_REASONS="${FAIL_REASONS}Low coverage; "
    fi

    # Threshold 5: Uniformity
    if (( $(echo "$UNIFORMITY < $MIN_UNIFORMITY_PCT" | bc -l) )); then
        STATUS="FAIL"
        FAIL_REASONS="${FAIL_REASONS}Poor uniformity; "
    fi

    # Clean up the fail reasons string (remove trailing "; ")
    FAIL_REASONS=$(echo "$FAIL_REASONS" | sed 's/; $//')
    
    # If sample passed, ensure fail_reasons is empty (not just whitespace)
    [[ "$STATUS" == "PASS" ]] && FAIL_REASONS=""

    # --------------------------------------------------
    # 4.9: Print progress summary
    # --------------------------------------------------
    echo ""
    echo "  Summary for $SAMPLE:"
    echo "    Total Reads:    $TOTAL_READS"
    echo "    Mapped:         ${MAPPED_PCT}%"
    echo "    On Target:      ${ON_TARGET_PCT}%"
    echo "    Mean Coverage:  ${MEAN_COV}x"
    echo "    Uniformity:     ${UNIFORMITY}%"
    echo "    Worst Amplicon: ${WORST_AMP}"
    echo "    STATUS:         $STATUS"
    [[ -n "$FAIL_REASONS" ]] && echo "    Fail Reasons:   $FAIL_REASONS"

    # --------------------------------------------------
    # 4.10: Write results to summary table
    # --------------------------------------------------
    printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\n" \
        "$SAMPLE" \
        "$TOTAL_READS" \
        "$MAPPED_PCT" \
        "$ON_TARGET_PCT" \
        "$MEAN_COV" \
        "$UNIFORMITY" \
        "$WORST_AMP" \
        "$STATUS" \
        "$FAIL_REASONS" >> "$SUMMARY"
    
    # --------------------------------------------------
    # 4.11: Copy BAM to appropriate folder (PASS or FAIL)
    # --------------------------------------------------
    # Determine destination based on QC status
    if [[ "$STATUS" == "PASS" ]]; then
        DEST="${OUT_DIR}/pass_bams"
    else
        DEST="${OUT_DIR}/fail_bams"
    fi
    
    # Copy BAM file (preserving permissions and timestamps)
    cp -p "$BAM" "$DEST/"
    
    # Also copy the BAM index if it exists
    if [[ -f "${BAM}.bai" ]]; then
        cp -p "${BAM}.bai" "$DEST/"
    fi
    
    echo "  → Copied to: $DEST"

done  # End of per-sample loop

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
echo ""
echo "========================================"
echo "DNA QC PIPELINE COMPLETED SUCCESSFULLY"
echo "========================================"
echo ""
echo "Summary file:   $SUMMARY"
echo "PASS samples:   ${OUT_DIR}/pass_bams/"
echo "FAIL samples:   ${OUT_DIR}/fail_bams/"
echo "Detailed stats: ${OUT_DIR}/stats/"
echo ""

# Count PASS and FAIL samples
PASS_COUNT=$(grep -c "PASS" "$SUMMARY" || echo 0)
FAIL_COUNT=$(grep -c "FAIL" "$SUMMARY" || echo 0)

echo "Results: $PASS_COUNT PASS, $FAIL_COUNT FAIL (out of $NUM_BAMS total)"
echo ""