#!/usr/bin/env bash
# ==============================================================================
# Script: 02_dna_qc.sh (UPDATED: HARMONISED COORDINATES)
# Description: DNA Panel QC with RS-ID/Gene Annotation.
# Fix: Forces "chr17" prefix for cross-cohort compatibility in Script 3.
# ==============================================================================

set -euo pipefail

# --- Final Thresholds ---
MIN_TOTAL_READS=120000        
MIN_MAPPED_PCT=90
MIN_ON_TARGET_PCT=75
MIN_MEAN_COVERAGE=200         
UNIFORMITY_THRESHOLD=100      
MIN_UNIFORMITY_PCT=95         
LOW_DEPTH_ALERT=50            

# --- Default Values ---
INPUT_DIR="" ; OUT_DIR="" ; REF_FA="" ; BED="" ; THREADS="${THREADS:-4}"

usage(){
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

[[ -z "${INPUT_DIR}" || -z "${OUT_DIR}" || -z "${REF_FA}" || -z "${BED}" ]] && { usage; exit 1; }

mkdir -p "${OUT_DIR}/pass_bams" "${OUT_DIR}/fail_bams" "${OUT_DIR}/stats"

# ---------------------------------------------------------------------------
# Step 1: BED Preparation 
# ---------------------------------------------------------------------------
# We ensure the chromosome column is ALWAYS "chr17" to avoid 17 vs chr17 mismatches
BED_ANNOTATED="${OUT_DIR}/targets.annotated.bed"
grep -v "^#" "${BED}" | awk 'BEGIN{FS=OFS="\t"} {
    # Fix chromosome naming at the source
    chrom = ($1 == "17" || $1 == "chr17") ? "chr17" : $1;
    label = ($5 == "." || $5 == "") ? $4 : $4"|"$5;
    if(NF>=6 && $6 != "." && $6 != ""){ label = label"|"$6 }
    print chrom,$2,$3,label
}' > "${BED_ANNOTATED}"

# Clean 3-column BED for calculation
BED_CALC="${OUT_DIR}/targets.calc.bed"
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' "${BED_ANNOTATED}" > "${BED_CALC}"

COHORT_TEMP="${OUT_DIR}/stats/all_low_amps.tmp"
: > "${COHORT_TEMP}"

# ---------------------------------------------------------------------------
# Step 2: Sample Processing Loop
# ---------------------------------------------------------------------------
mapfile -t BAMS < <(find "${INPUT_DIR}" -maxdepth 2 -type f -name "*.bam" | sort)
NUM_BAMS=${#BAMS[@]}

if [ "$NUM_BAMS" -eq 0 ]; then
    echo "ERROR: No BAM files found in ${INPUT_DIR}"
    exit 1
fi

echo "[*] Found $NUM_BAMS samples. Starting harmonised cohort analysis..."

SUMMARY="${OUT_DIR}/qc_summary.tsv"
echo -e "sample\ttotal_reads\tmapped_pct\ton_target_pct\tmean_cov\tuniformity_100x\tworst_amplicon\tstatus\tfail_reasons" > "${SUMMARY}"

for ((i=0; i<NUM_BAMS; i++)); do
  BAM="${BAMS[$i]}"
  SAMPLE=$(basename "${BAM}" .bam)
  SDIR="${OUT_DIR}/stats/${SAMPLE}"
  
  echo "[$(($i+1))/$NUM_BAMS] Analysing: $SAMPLE"
  rm -rf "${SDIR}" && mkdir -p "${SDIR}"

  FLAGSTAT=$(samtools flagstat -@ "${THREADS}" "${BAM}")
  TOTAL_READS=$(echo "$FLAGSTAT" | awk '/in total/ {print $1}')
  MAPPED_READS=$(echo "$FLAGSTAT" | awk '/ mapped / && !/primary/ {print $1}')
  MAPPED_PCT=$(awk -v m="${MAPPED_READS}" -v t="${TOTAL_READS}" 'BEGIN{if(t>0){printf "%.2f",(m/t)*100}else{print "0"}}')
  
  ON_TARGET_READS=$(samtools view -@ "${THREADS}" -F 4 -c -L "${BED_CALC}" "${BAM}")
  ON_TARGET_PCT=$(awk -v o="${ON_TARGET_READS}" -v m="${MAPPED_READS}" 'BEGIN{if(m>0){printf "%.2f",(o/m)*100}else{print "0"}}')

  # Coverage calculation (mosdepth will now use the harmonised BED_CALC)
  mosdepth --threads "${THREADS}" --no-per-base --by "${BED_CALC}" --fast-mode "${SDIR}/${SAMPLE}" "${BAM}"
  REG_BEDGZ="${SDIR}/${SAMPLE}.regions.bed.gz"

  MEAN_COV=$(zcat "${REG_BEDGZ}" | awk 'BEGIN{sum=0;bp=0}{len=$3-$2;sum+=len*$4;bp+=len}END{if(bp>0){printf "%.2f",sum/bp}else{print "0"}}')
  UNIFORMITY=$(zcat "${REG_BEDGZ}" | awk -v t="$UNIFORMITY_THRESHOLD" 'BEGIN{c=0;tot=0}{tot++; if($4>=t)c++}END{if(tot>0)printf "%.2f",(c/tot)*100; else print "0"}')
  WORST_COV=$(zcat "${REG_BEDGZ}" | awk 'NR==1{min=$4}{if($4<min)min=$4}END{printf "%.2f",min}')

  paste <(zcat "${REG_BEDGZ}") <(awk '{print $4}' "${BED_ANNOTATED}") | \
  awk -v t="$LOW_DEPTH_ALERT" '$4 < t {print $5}' >> "${COHORT_TEMP}"

  STATUS="PASS"; FAILS=()
  [[ "${TOTAL_READS}" -lt "${MIN_TOTAL_READS}" ]] && STATUS="FAIL" && FAILS+=("low_reads")
  (( $(awk 'BEGIN{print ('"${MAPPED_PCT}"'<'"${MIN_MAPPED_PCT}"')}') )) && STATUS="FAIL" && FAILS+=("low_mapping")
  (( $(awk 'BEGIN{print ('"${ON_TARGET_PCT}"'<'"${MIN_ON_TARGET_PCT}"')}') )) && STATUS="FAIL" && FAILS+=("low_on_target")
  (( $(awk 'BEGIN{print ('"${MEAN_COV}"'<'"${MIN_MEAN_COVERAGE}"')}') )) && STATUS="FAIL" && FAILS+=("low_mean_cov")
  (( $(awk 'BEGIN{print ('"${UNIFORMITY}"'<'"${MIN_UNIFORMITY_PCT}"')}') )) && STATUS="FAIL" && FAILS+=("low_uniformity")

  printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\n" \
    "${SAMPLE}" "${TOTAL_READS}" "${MAPPED_PCT}" "${ON_TARGET_PCT}" "${MEAN_COV}" "${UNIFORMITY}" "${WORST_COV}" \
    "${STATUS}" "$(IFS=','; echo "${FAILS[*]-}")" >> "${SUMMARY}"

  DEST=$([[ "${STATUS}" == "PASS" ]] && echo "${OUT_DIR}/pass_bams" || echo "${OUT_DIR}/fail_bams")
  cp -p "${BAM}" "${DEST}/" && [[ -f "${BAM}.bai" ]] && cp -p "${BAM}.bai" "${DEST}/"
done


# Step 3: Global Failure Report (Old behaviour + coordinates + percentage)
echo "[*] Generating Global Failure Report using ${LOW_DEPTH_ALERT}x threshold..."

GLOBAL_REPORT="${OUT_DIR}/cohort_amplicon_failure_report.txt"
TEMP_IDS=$(mktemp)

# Total number of samples in cohort (already defined earlier)
TOTAL_SAMPLES="${NUM_BAMS}"

# 1. Extract amplicon IDs with depth < LOW_DEPTH_ALERT
find "${OUT_DIR}/stats/" -name "*.regions.bed.gz" -print0 \
| xargs -0 zcat \
| awk -v thr="${LOW_DEPTH_ALERT}" 'NR==FNR {id[$1":"$2":"$3]=$4; next}
    {
      key=$1":"$2":"$3;
      if($4 < thr && key in id) print id[key];
    }' "${BED_ANNOTATED}" - \
| sort \
| uniq -c \
> "${TEMP_IDS}"

# 2. If nothing failed, skip
if [[ ! -s "${TEMP_IDS}" ]]; then
    echo "No regions found below ${LOW_DEPTH_ALERT}x. Report not created."
    rm -f "${TEMP_IDS}"
    exit 0
fi

# 3. Header
echo -e "Chr\tStart\tEnd\tID\tFail_Count\tFail_Percent" > "${GLOBAL_REPORT}"

# 4. Match IDs back to BED_ANNOTATED, calculate percentage
while read -r COUNT AMP_ID; do
    awk -v id="${AMP_ID}" -v cnt="${COUNT}" -v total="${TOTAL_SAMPLES}" '
        $4 == id {
            perc = (cnt / total) * 100;
            printf "%s\t%s\t%s\t%s\t%s\t%.2f\n", $1, $2, $3, id, cnt, perc;
        }
    ' "${BED_ANNOTATED}" >> "${GLOBAL_REPORT}"
done < "${TEMP_IDS}"

echo "[SUCCESS] Global failure report saved to: ${GLOBAL_REPORT}"
rm -f "${TEMP_IDS}"
