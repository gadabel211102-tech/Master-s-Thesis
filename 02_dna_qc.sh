#!/usr/bin/env bash
# Quality control pipeline for targeted DNA sequencing data (Ion Torrent).
# Computes per-sample coverage metrics and classifies samples as PASS or FAIL.
# Usage: ./02_dna_qc.sh -i INPUT_DIR -o OUT_DIR -r REF_FASTA -b TARGETS_BED [-t THREADS]
set -euo pipefail

# --- QC thresholds -----------------------------------------------------------
# These values define the minimum acceptable performance for a sample to pass.
# Thresholds were selected based on standard practice for targeted amplicon panels.
MIN_TOTAL_READS=120000      # samples with fewer reads lack statistical power for variant calling
MIN_MAPPED_PCT=90           # high mapping rate expected given a targeted, amplicon-based design
MIN_ON_TARGET_PCT=75        # reflects expected enrichment efficiency for the panel used
MIN_MEAN_COVERAGE=200       # minimum depth required for reliable heterozygous variant detection
UNIFORMITY_THRESHOLD=100    # depth at which an amplicon is considered adequately covered
MIN_UNIFORMITY_PCT=95       # percentage of amplicons that must meet the uniformity threshold

# --- Input parameters --------------------------------------------------------
INPUT_DIR=""
OUT_DIR=""
REF_FA=""
BED=""
THREADS="${THREADS:-4}"

usage() {
    cat <<USAGE
Usage: $0 -i INPUT_DIR -o OUT_DIR -r REF_FASTA -b TARGETS_BED [-t THREADS]
  -i  Directory containing input BAM files
  -o  Output directory for QC results
  -r  Reference genome FASTA (e.g. hg38.fa)
  -b  BED file defining target amplicon coordinates
  -t  CPU threads (default: 4)
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT_DIR="$2"; shift 2 ;;
        -o) OUT_DIR="$2";   shift 2 ;;
        -r) REF_FA="$2";    shift 2 ;;
        -b) BED="$2";       shift 2 ;;
        -t) THREADS="$2";   shift 2 ;;
        *) usage; exit 2 ;;
    esac
done

[[ -z "$INPUT_DIR" || -z "$OUT_DIR" || -z "$REF_FA" || -z "$BED" ]] && { usage; exit 1; }

# Create output subdirectories for passing samples, failing samples, and per-sample stats
mkdir -p "${OUT_DIR}/pass_bams" "${OUT_DIR}/fail_bams" "${OUT_DIR}/stats"

# --- Step 1: Clean and normalise BED file ------------------------------------
echo "[*] Step 1: Cleaning and normalising BED file..."

BED_CALC="${OUT_DIR}/targets.sorted.bed"

# Strip Windows line endings, remove UCSC track/browser header lines, enforce
# chr-prefixed chromosome names, require at least 3 columns, and sort by
# coordinate — all columns (including rsID and gene in cols 5-6) are preserved
awk 'BEGIN{FS="\t"; OFS="\t"} {
    gsub(/\r/,"");
    if ($1 ~ /^(track|browser)/) next;
    if ($1 ~ /^(chr|[0-9XYM])/) {
        if (NF >= 3) {
            if ($1 !~ /^chr/) $1="chr"$1;
            print $0
        }
    }
}' "$BED" | sort -k1,1 -k2,2n > "$BED_CALC"

[[ -s "$BED_CALC" ]] || { echo "ERROR: cleaned BED file is empty — check input format"; exit 1; }
echo "[*] BED file cleaned: $BED_CALC"

# --- Step 2: Discover BAM files ----------------------------------------------
echo "[*] Step 2: Discovering BAM files..."

mapfile -t BAMS < <(find "$INPUT_DIR" -maxdepth 2 -type f -name "*.bam" | sort)
NUM_BAMS=${#BAMS[@]}
[[ $NUM_BAMS -gt 0 ]] || { echo "ERROR: no BAM files found in $INPUT_DIR"; exit 1; }
echo "[*] Found $NUM_BAMS BAM file(s)"

# --- Step 3: Initialise summary table ----------------------------------------
# One row per sample; columns record each QC metric and the final PASS/FAIL verdict
SUMMARY="${OUT_DIR}/qc_summary.tsv"
echo -e "sample\ttotal_reads\tmapped_pct\ton_target_pct\tmean_cov\tuniformity_100x\tworst_amplicon\tstatus\tfail_reasons" > "$SUMMARY"

# --- Step 4: Process each sample ---------------------------------------------
echo "[*] Step 4: Processing samples..."

for (( i=0; i<NUM_BAMS; i++ )); do

    BAM="${BAMS[$i]}"
    SAMPLE=$(basename "$BAM" .bam)     # derive sample name from filename
    SDIR="${OUT_DIR}/stats/${SAMPLE}"
    mkdir -p "$SDIR"

    echo ""
    echo "================================================"
    echo "Processing sample $((i+1))/$NUM_BAMS: $SAMPLE"
    echo "================================================"

    # [1/7] Basic read counts via samtools flagstat
    echo "  [1/7] Calculating BAM statistics..."
    FLAGSTAT=$(samtools flagstat -@ "$THREADS" "$BAM")
    TOTAL_READS=$(echo "$FLAGSTAT" | awk '/in total/ {print $1}')
    # Exclude the "primary mapped" line to avoid double-counting supplementary alignments
    MAPPED_READS=$(echo "$FLAGSTAT" | awk '/ mapped / && !/primary/ {print $1}')
    MAPPED_PCT=$(awk -v m="$MAPPED_READS" -v t="$TOTAL_READS" \
        'BEGIN{if(t>0) printf "%.2f", (m/t)*100; else print "0"}')

    # [2/7] On-target rate: proportion of mapped reads overlapping the target BED
    echo "  [2/7] Calculating on-target percentage..."
    # -F 4 excludes unmapped reads; -L restricts counting to target regions; -c returns a count only
    ON_TARGET_READS=$(samtools view -@ "$THREADS" -F 4 -c -L "$BED_CALC" "$BAM")
    ON_TARGET_PCT=$(awk -v o="$ON_TARGET_READS" -v m="$MAPPED_READS" \
        'BEGIN{if(m>0) printf "%.2f", (o/m)*100; else print "0"}')

    # [3/7] Per-region depth via mosdepth
    echo "  [3/7] Running mosdepth coverage analysis..."
    # --no-per-base and --fast-mode skip per-base output, reducing runtime and disk usage;
    # --by computes mean depth for each region defined in the target BED file
    mosdepth --threads "$THREADS" \
             --no-per-base \
             --by "$BED_CALC" \
             --fast-mode \
             "${SDIR}/${SAMPLE}" \
             "$BAM"
    # Produces: ${SAMPLE}.regions.bed.gz (cols: chr, start, end, region_id, mean_depth)

    # [4/7] Weighted mean coverage across all amplicons
    echo "  [4/7] Calculating mean coverage..."
    REG="${SDIR}/${SAMPLE}.regions.bed.gz"
    # Weighted by region length so that larger amplicons contribute proportionally to the mean
    MEAN_COV=$(zcat "$REG" | awk '
        BEGIN { sum=0; bp=0 }
        {
            len  = $3 - $2
            sum += len * $5    # length x depth for this region
            bp  += len
        }
        END { if(bp>0) printf "%.2f", sum/bp; else print "0" }
    ')

    # [5/7] Uniformity: fraction of amplicons meeting the minimum depth threshold
    echo "  [5/7] Calculating coverage uniformity..."
    UNIFORMITY=$(zcat "$REG" | awk -v thr="$UNIFORMITY_THRESHOLD" '
        BEGIN { pass=0; total=0 }
        {
            total++
            if($5 >= thr) pass++    # col 5 = mean depth for this amplicon
        }
        END { if(total>0) printf "%.2f", (pass/total)*100; else print "0" }
    ')

    # [6/7] Identify worst-performing amplicon (lowest mean depth)
    echo "  [6/7] Identifying worst amplicon..."
    # Format output as region_id:depth, sort ascending by depth, return the minimum
    WORST_AMP=$(zcat "$REG" | awk '{print $4 ":" sprintf("%.2f",$5)}' | \
        sort -t: -k2,2n | head -1)
    [[ -z "$WORST_AMP" ]] && WORST_AMP="N/A"

    # [7/7] Evaluate all thresholds; accumulate failure reasons if any are breached
    echo "  [7/7] Evaluating QC status..."
    STATUS="PASS"
    FAIL_REASONS=""

    # bc -l is used for floating-point comparisons, which bash cannot perform natively
    (( $(echo "$TOTAL_READS < $MIN_TOTAL_READS"     | bc -l) )) && STATUS="FAIL" && FAIL_REASONS+="Low total reads; "
    (( $(echo "$MAPPED_PCT < $MIN_MAPPED_PCT"       | bc -l) )) && STATUS="FAIL" && FAIL_REASONS+="Low mapping rate; "
    (( $(echo "$ON_TARGET_PCT < $MIN_ON_TARGET_PCT" | bc -l) )) && STATUS="FAIL" && FAIL_REASONS+="Low on-target rate; "
    (( $(echo "$MEAN_COV < $MIN_MEAN_COVERAGE"      | bc -l) )) && STATUS="FAIL" && FAIL_REASONS+="Low coverage; "
    (( $(echo "$UNIFORMITY < $MIN_UNIFORMITY_PCT"   | bc -l) )) && STATUS="FAIL" && FAIL_REASONS+="Poor uniformity; "

    FAIL_REASONS=$(echo "$FAIL_REASONS" | sed 's/; $//')
    [[ "$STATUS" == "PASS" ]] && FAIL_REASONS=""

    # Print per-sample summary to stdout
    echo ""
    echo "  Summary for $SAMPLE:"
    echo "    Total reads:    $TOTAL_READS"
    echo "    Mapped:         ${MAPPED_PCT}%"
    echo "    On-target:      ${ON_TARGET_PCT}%"
    echo "    Mean coverage:  ${MEAN_COV}x"
    echo "    Uniformity:     ${UNIFORMITY}%"
    echo "    Worst amplicon: ${WORST_AMP}"
    echo "    Status:         $STATUS"
    [[ -n "$FAIL_REASONS" ]] && echo "    Fail reasons:   $FAIL_REASONS"

    # Append this sample's metrics to the summary TSV
    printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\n" \
        "$SAMPLE" "$TOTAL_READS" "$MAPPED_PCT" "$ON_TARGET_PCT" \
        "$MEAN_COV" "$UNIFORMITY" "$WORST_AMP" "$STATUS" "$FAIL_REASONS" >> "$SUMMARY"

    # Copy BAM (and index if present) to the appropriate pass/fail directory
    if [[ "$STATUS" == "PASS" ]]; then
        DEST="${OUT_DIR}/pass_bams"
    else
        DEST="${OUT_DIR}/fail_bams"
    fi
    cp -p "$BAM" "$DEST/"
    [[ -f "${BAM}.bai" ]] && cp -p "${BAM}.bai" "$DEST/"
    echo "  -> Copied to: $DEST"

done

# --- Final summary -----------------------------------------------------------
PASS_COUNT=$(grep -c "PASS" "$SUMMARY" || echo 0)
FAIL_COUNT=$(grep -c "FAIL" "$SUMMARY" || echo 0)

echo ""
echo "========================================"
echo "QC pipeline complete."
echo "  Summary:          $SUMMARY"
echo "  PASS samples:     ${OUT_DIR}/pass_bams/"
echo "  FAIL samples:     ${OUT_DIR}/fail_bams/"
echo "  Per-sample stats: ${OUT_DIR}/stats/"
echo "  Results: $PASS_COUNT PASS, $FAIL_COUNT FAIL (of $NUM_BAMS total)"
echo "========================================"
