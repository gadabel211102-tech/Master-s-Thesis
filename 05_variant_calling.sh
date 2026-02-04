#!/usr/bin/env bash
# Exit on error (-e), unset variables (-u), or pipe failures (-o pipefail)
set -euo pipefail

# --- Defaults ---
# Initialise file paths and default numerical values for the pipeline
MANIFEST=""; REF_FA=""; BED=""; OUT_DIR=""; THREADS="${THREADS:-4}"
TVC_JSON=""; HOTSPOT_VCF=""; ERROR_MOTIF=""
MIN_DP="${MIN_DP:-50}"; MIN_QUAL="${MIN_QUAL:-20}"

# --- Helpers ---
# Gets current date/time for logging
_now(){ date '+%Y-%m-%d %H:%M:%S'; }
# Prints a message with a timestamp
log(){ printf '[%s] %s\n' "$(_now)" "$*"; }
# Converts seconds into a readable HH:MM:SS format
fmt_dur(){ local s="$1"; printf '%02d:%02d:%02d' $((s/3600)) $(((s%3600)/60)) $((s%60)); }
# Calculates a future timestamp (Estimated Time of Arrival)
eta_ts(){ local s="$1"; date -d "@$(( $(date +%s) + s ))" '+%Y-%m-%d %H:%M:%S'; }
# Prints a horizontal line across the terminal for visual clarity
hr(){ printf '%*s' "${COLUMNS:-100}" '' | tr ' ' '-'; echo; }
# Draws a visual progress bar based on a percentage input
progress_bar(){ local pct="$1"; local width=30; local filled=$(( pct * width / 100 )); local empty=$(( width - filled ));
  printf '\r[%s%s] %3d%%' "$(printf '%0.s#' $(seq 1 $filled))" "$(printf '%0.s-' $(seq 1 $empty))" "${pct}"; }
# Computes a running average to predict remaining processing time
running_avg(){ local avg="$1" n="$2" new="$3"; avg=${avg//[^0-9]/}; n=${n//[^0-9]/}; new=${new//[^0-9]/}; echo $(( (avg * n + new) / (n + 1) )); }

# --- Step Tracking ---
STEP_NAME=""; STEP_T0=0
# Starts a timer and logs the beginning of a specific pipeline stage
step_start(){ STEP_NAME="$1"; STEP_T0=$(date +%s); log "▶ ${STEP_NAME}…"; }
# Ends the timer, logs duration, and returns the elapsed time
step_done(){ local t1=$(date +%s); local dt=$(( t1 - STEP_T0 )); log "✔ ${STEP_NAME} done ($(fmt_dur ${dt}))"; echo "${dt}"; }

CURRENT_SAMPLE=""
# Triggers an error message if any command fails during execution
trap 'echo; log "✖ Mistake while processing sample: ${CURRENT_SAMPLE:-<none>}"; exit 1' ERR

# Prints the help manual for the script
usage(){ cat <<'USAGE'
Usage: 05_variant_calling.sh -m MANIFEST -r REF.fa -b TARGETS.bed -o OUT_DIR [options]
...
USAGE
}

# --- Parse args ---
# Iterate through command line flags and assign them to variables
while (( $# > 0 )); do
  case "$1" in
    -m) MANIFEST="$2"; shift 2 ;;
    -r) REF_FA="$2";  shift 2 ;;
    -b) BED="$2";      shift 2 ;;
    -o) OUT_DIR="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    --min-dp)   MIN_DP="$2";   shift 2 ;;
    --min-qual) MIN_QUAL="$2"; shift 2 ;;
    --tvc-json) TVC_JSON="$2"; shift 2 ;;
    --hotspot-vcf) HOTSPOT_VCF="$2"; shift 2 ;;
    --error-motif) ERROR_MOTIF="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 2 ;;
  esac
done

# --- Basic checks ---
# Ensure required arguments are not empty strings
if [[ -z "${MANIFEST}" || -z "${REF_FA}" || -z "${BED}" || -z "${OUT_DIR}" ]]; then
  usage; exit 1
fi
# Verify that all required bioinformatics tools are installed on the system
for tool in bcftools samtools tabix bgzip gzip awk sed grep; do
  command -v "${tool}" >/dev/null || { echo "ERROR: ${tool} not found"; exit 1; }
done
# Specifically check for the Torrent Variant Caller (TVC) binary
if ! command -v tvc >/dev/null && ! command -v variant_caller_pipeline.py >/dev/null; then
  echo "ERROR: TVC or wrapper not found"; exit 1
fi
# Ensure the provided files actually exist on the disk
[[ -f "${REF_FA}" ]] || { echo "ERROR: reference FASTA not found"; exit 1; }
[[ -f "${BED}"    ]] || { echo "ERROR: BED not found"; exit 1; }
[[ -f "${MANIFEST}" ]] || { echo "ERROR: manifest not found"; exit 1; }

# --- Log configuration ---
# Create output and log folders
mkdir -p "${OUT_DIR}" "${OUT_DIR}/logs"
# Log the current configuration to the terminal
log "Output dir: ${OUT_DIR}"
log "Manifest:   ${MANIFEST}"
log "Reference:  ${REF_FA}"
log "Targets BED:${BED}"
log "Threads:    ${THREADS}"
log "Filters:    MIN_DP=${MIN_DP}, MIN_QUAL=${MIN_QUAL}"
hr

# --- Reference index ---
# Create a .fai index if it doesn't exist (required for samtools)
if [[ ! -f "${REF_FA}.fai" ]]; then
  step_start "Index reference"
  samtools faidx "${REF_FA}"
  step_done >/dev/null
fi

# --------------------------- BED normalisation ------------------
# Set path for the cleaned BED file
BED_NORM="${OUT_DIR}/targets.cleaned.bed"
step_start "Normalize BED -> ${BED_NORM}"
# Clean BED: remove metadata/comments, remove Windows line endings, keep first 6 columns
awk 'BEGIN{OFS="\t"} !/^track/ && !/^browser/ && !/^#/ { if(NF>=3){print $1,$2,$3,$4,$5,$6} }' "${BED}" | \
  sed 's/\r$//' > "${BED_NORM}"
step_done >/dev/null

# Load manifest into array; ignore comment lines and empty whitespace
mapfile -t BAMS < <(grep -v -E '^[[:space:]]*#' "${MANIFEST}" | sed '/^[[:space:]]*$/d')
# Stop if no BAM files were parsed from the manifest
[[ ${#BAMS[@]} -gt 0 ]] || { echo "ERROR: empty manifest"; exit 1; }
TOTAL=${#BAMS[@]}
log "Samples to process: ${TOTAL}"
hr

# Initialise arrays and counters for performance tracking and ETA
SAMPLE_VCFS=(); IDX=0; S_DONE=0; AVG_SAMPLE=0; AVG_TVC=0; N_TVC=0; AVG_NORM=0; N_NORM=0; AVG_FILTER=0; N_FILTER=0

# Start processing each BAM file found in the manifest
for BAM in "${BAMS[@]}"; do
  IDX=$((IDX+1))
  # Warn and skip if the specific BAM file cannot be found on disk
  if [[ ! -f "${BAM}" ]]; then
    log "WARN: BAM missing, skipping: ${BAM}"; continue
  fi
  
  # Determine sample name and create a filename-safe version
  SAMPLE=$(basename "$BAM" .bam)
  SAFE_SAMPLE="$(printf '%s' "${SAMPLE}" | tr -c 'A-Za-z0-9._-' '_')"
  SDIR="${OUT_DIR}/${SAFE_SAMPLE}"

  mkdir -p "${SDIR}/logs"

  # Create links for BAMs with "clean" names to avoid TVC tool errors with special characters
  bam_dir=$(dirname -- "${BAM}")
  bam_base=$(basename -- "${BAM}")
  safe_base="$(printf '%s' "${bam_base}" | tr -c 'A-Za-z0-9._-' '_')"
  SAFE_BAM="${bam_dir}/${safe_base}"
  if [[ "${SAFE_BAM}" != "${BAM}" && ! -e "${SAFE_BAM}" ]]; then
    ln -s "${BAM}" "${SAFE_BAM}"
    # Link corresponding index (.bai) files if they exist
    if   [[ -f "${BAM}.bai" ]]; then ln -s "${BAM}.bai" "${SAFE_BAM}.bai" 2>/dev/null || true
    elif [[ -f "${BAM%.bam}.bai" ]]; then ln -s "${BAM%.bam}.bai" "${SAFE_BAM}.bai" 2>/dev/null || true
    fi
  fi
  TVC_BAM="${SAFE_BAM}"

  # Build a BAM index if it is missing (required for variant calling)
  if [[ ! -f "${TVC_BAM}.bai" ]]; then
    step_start "Index BAM (${SAMPLE})"
    samtools index -@ "${THREADS}" "${TVC_BAM}"
    step_done >/dev/null
  fi

  # Log progress and calculate ETA based on the running average of previous samples
  echo; log "===== Sample ${IDX}/${TOTAL}: ${SAMPLE} ====="
  if [[ ${S_DONE} -gt 0 ]]; then
    REMAINING=$((TOTAL-(IDX-1)))
    est_secs=$((AVG_SAMPLE*REMAINING))
    log "ETA before start: ~$(fmt_dur ${est_secs}) (ETA $(eta_ts ${est_secs}))"
  fi
  log "Collecting timings on first sample…"
  progress_bar 0; printf '\n'; SAMPLE_T0=$(date +%s)

  # --- TVC call ---
  step_start "TVC (variant calling)"
  RAW_VCF=""
  # Choose between the Python wrapper script or the direct TVC binary
  if command -v variant_caller_pipeline.py >/dev/null; then
    # Construct command for Python wrapper
    TVC_CMD=(variant_caller_pipeline.py --input-bam "${TVC_BAM}" --reference-fasta "${REF_FA}" --region-bed "${BED_NORM}" --output-dir "${SDIR}" --num-threads "${THREADS}")
    [[ -n "${TVC_JSON}"    ]] && TVC_CMD+=(--parameters-file "${TVC_JSON}")
    [[ -n "${HOTSPOT_VCF}" ]] && TVC_CMD+=(--input-vcf       "${HOTSPOT_VCF}")
    [[ -n "${ERROR_MOTIF}" ]] && TVC_CMD+=(--error-motif     "${ERROR_MOTIF}")
    # Execute and log output
    "${TVC_CMD[@]}" > "${SDIR}/logs/tvc.stdout.log" 2> "${SDIR}/logs/tvc.stderr.log"
    RAW_VCF="${SDIR}/TSVC_variants.vcf"
    [[ ! -f "${RAW_VCF}" ]] && RAW_VCF="${SDIR}/variants.vcf"
  else
    # Construct command for direct binary call
    TVC_CMD=(tvc --input-bam "${TVC_BAM}" --reference "${REF_FA}" --output-dir "${SDIR}" --num-threads "${THREADS}" --target-file "${BED_NORM}" --output-vcf "small_variants.vcf")
    "${TVC_CMD[@]}" > "${SDIR}/logs/tvc.stdout.log" 2> "${SDIR}/logs/tvc.stderr.log"
    RAW_VCF="${SDIR}/TSVC_variants.vcf"
    [[ ! -f "${RAW_VCF}" ]] && RAW_VCF="${SDIR}/small_variants.vcf"
  fi
  # Verify that the variant calling produced an output file
  [[ ! -f "${RAW_VCF}" ]] && { log "ERROR: VCF not found after TVC"; exit 1; }
  tvc_dt=$(step_done); progress_bar 33; printf '\n'
  # Update the average time taken for the TVC step
  AVG_TVC=$(running_avg "${AVG_TVC}" "${N_TVC}" "${tvc_dt}"); N_TVC=$((N_TVC+1))

  # --- VCF Header Correction ---
  step_start "Validate/repair VCF header (reheader + optional rename)"
  # Check for ##contig lines; reheader using FASTA index if they are missing
  if ! bcftools view -h "${RAW_VCF}" | grep -q '^##contig'; then
    log "VCF header lacks ##contig. Reinserting contigs using ${REF_FA}.fai..."
    bcftools reheader --fai "${REF_FA}.fai" "${RAW_VCF}" -o "${SDIR}/variants.hdr.vcf" 
    bgzip -f "${SDIR}/variants.hdr.vcf"
    tabix -p vcf "${SDIR}/variants.hdr.vcf.gz"
    RAW_VCF="${SDIR}/variants.hdr.vcf.gz"
  fi
  # Determine if Reference and VCF use the same "chr" prefix convention
  HAS_CHR_FA=$(awk -F'\t' 'NR==1{print ($1 ~ /^chr/)?1:0}' "${REF_FA}.fai")
  FIRST_VAR=$(bcftools view -H "${RAW_VCF}" | head -n1 || true)
  HAS_CHR_VCF=0; [[ -n "${FIRST_VAR}" && "${FIRST_VAR%%\t*}" =~ ^chr ]] && HAS_CHR_VCF=1
  
  # Harmonise contig names if a mismatch exists
  if [[ -n "${FIRST_VAR}" ]]; then
    if [[ "${HAS_CHR_FA}" -eq 1 && "${HAS_CHR_VCF}" -eq 0 ]]; then
      # Add "chr" prefix to VCF records
      awk -F'\t' '{g=$1; sub(/^chr/,"",g); print g "\t" $1}' "${REF_FA}.fai" > "${SDIR}/chr_map.txt"
      bcftools annotate --rename-chrs "${SDIR}/chr_map.txt" "${RAW_VCF}" -Oz -o "${SDIR}/variants.renamed.vcf.gz"
      tabix -p vcf "${SDIR}/variants.renamed.vcf.gz"
      RAW_VCF="${SDIR}/variants.renamed.vcf.gz"
    elif [[ "${HAS_CHR_FA}" -eq 0 && "${HAS_CHR_VCF}" -eq 1 ]]; then
      # Remove "chr" prefix from VCF records
      awk -F'\t' '{g=$1; print $1 "\t" gensub(/^chr/,"","",g)}' "${REF_FA}.fai" > "${SDIR}/chr_map.txt"
      bcftools annotate --rename-chrs "${SDIR}/chr_map.txt" "${RAW_VCF}" -Oz -o "${SDIR}/variants.renamed.vcf.gz"
      tabix -p vcf "${SDIR}/variants.renamed.vcf.gz"
      RAW_VCF="${SDIR}/variants.renamed.vcf.gz"
    else
      log "Nomenclatura de contigs coherente entre FASTA y VCF."
    fi
  else
    log "VCF sin variantes; saltamos rename-chrs. (Header ya reinyectado si faltaba)"
  fi
  # Ensure the VCF is not corrupt and contains the required contig headers
  if ! bcftools view -h "${RAW_VCF}" | grep -q '^##contig'; then
    log "ERROR: header aún sin ##contig tras reheader. Abortando para evitar corrupción."; exit 1
  fi
  step_done

 # --- Normalise + compress ---
  step_start "Normalise + compress"
  # Split multi-allelic sites and check variants against the reference FASTA
  bcftools norm -f "${REF_FA}" -m -both "${RAW_VCF}" -Oz -o "${SDIR}/variants.norm.vcf.gz"
  # Index the newly normalised VCF
  tabix -p vcf "${SDIR}/variants.norm.vcf.gz"
  # Log duration and update the running average for normalisation time
  norm_dt=$(step_done); progress_bar 66; printf '\n'
  AVG_NORM=$(running_avg "${AVG_NORM}" "${N_NORM}" "${norm_dt}"); N_NORM=$((N_NORM+1))

  # --- Unified filtering ---
  step_start "Filter (PASS; DP>=${MIN_DP}; QUAL>=${MIN_QUAL})"
  # Filter for high-quality variants (QUAL), depth (DP), and the 'PASS' status
  bcftools view -i "MAX(FMT/DP)>=${MIN_DP} && QUAL>=${MIN_QUAL} && FILTER='PASS'" \
    "${SDIR}/variants.norm.vcf.gz" -Oz -o "${SDIR}/variants.filtered.vcf.gz"
  tabix -p vcf "${SDIR}/variants.filtered.vcf.gz"
  # Complete the per-sample progress bar and update timing
  filter_dt=$(step_done); progress_bar 100; printf '\n'
  AVG_FILTER=$(running_avg "${AVG_FILTER}" "${N_FILTER}" "${filter_dt}"); N_FILTER=$((N_FILTER+1))

  # --- Export manual review candidates ---
  # Define helper functions to check if specific INFO or FORMAT tags exist in the header
  has_info_tag(){ bcftools view -h "$1" | grep -q "^##INFO=<ID=$2,"; }
  has_fmt_tag(){  bcftools view -h "$1" | grep -q "^##FORMAT=<ID=$2,"; }
  VCF_HDR="${SDIR}/variants.norm.vcf.gz"
  INFO_TAGS=(HP Indel STB STBP HS SB) # Sequencing bias tags for manual inspection

  # 1) Per-sample export: Extract metrics for review
  step_start "Export manual review (per-sample metrics)"
  # Determine if AO (Allele Observation) and RO (Reference Observation) tags are available
  AO_RO_FMT=0; has_fmt_tag "${VCF_HDR}" AO && has_fmt_tag "${VCF_HDR}" RO && AO_RO_FMT=1 || true
  AF_FMT=0;    has_fmt_tag "${VCF_HDR}" AF && AF_FMT=1 || true
  
  if [[ "${AO_RO_FMT}" -eq 1 ]]; then
    # Extract raw counts to a temporary file
    bcftools query -f "[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%SAMPLE\t%DP\t%AO\t%RO\n]" \
      "${VCF_HDR}" > "${SDIR}/manual_review_samples.raw.tsv" || true
    # Calculate Allele Frequency (AF) manually using awk (AO / (AO + RO))
    awk -F'\t' 'BEGIN{OFS="\t"} {
      chrom=$1; pos=$2; ref=$3; alt=$4; qual=$5; filt=$6; sample=$7; dp=$8; ao=$9; ro=$10;
      if (ao ~ /^[0-9]+$/ && ro ~ /^[0-9]+$/ && (ao+ro)>0) { af = ao/(ao+ro) } else { af = "." }
      print chrom, pos, ref, alt, qual, filt, sample, dp, af
    }' "${SDIR}/manual_review_samples.raw.tsv" > "${SDIR}/manual_review_samples.tsv" || true
    rm -f "${SDIR}/manual_review_samples.raw.tsv"
  else
    # Fallback: export existing AF tag or just the DP if AF is missing
    if [[ "${AF_FMT}" -eq 1 ]]; then
      bcftools query -f "[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%SAMPLE\t%DP\t%AF\n]" \
        "${VCF_HDR}" > "${SDIR}/manual_review_samples.tsv" || true
    else
      bcftools query -f "[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%SAMPLE\t%DP\n]" \
        "${VCF_HDR}" > "${SDIR}/manual_review_samples.tmp.tsv" || true
      awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,"."}' \
        "${SDIR}/manual_review_samples.tmp.tsv" > "${SDIR}/manual_review_samples.tsv" || true
      rm -f "${SDIR}/manual_review_samples.tmp.tsv"
    fi
  fi
  step_done

  # 2) Site-level export: Dynamically build format string based on available INFO tags
  step_start "Export manual review (site-level INFO tags)"
  FMT_STR="%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER"
  if bcftools view -h "${VCF_HDR}" | grep -q '^##INFO=<ID=DP,'; then FMT_STR+="\t%INFO/DP"; else FMT_STR+="\t."; fi
  if bcftools view -h "${VCF_HDR}" | grep -q '^##INFO=<ID=AF,'; then FMT_STR+="\t%INFO/AF"; else FMT_STR+="\t."; fi
  # Loop through potential bias tags and add them to the TSV columns if they exist
  for tag in "${INFO_TAGS[@]}"; do
    if bcftools view -h "${VCF_HDR}" | grep -q "^##INFO=<ID=${tag},"; then FMT_STR+="\t%INFO/${tag}"; else FMT_STR+="\t."; fi
  done
  FMT_STR+="\n"
  bcftools query -f "${FMT_STR}" "${VCF_HDR}" > "${SDIR}/manual_review_site.tsv" || true
  step_done

  # --- Summary & bookkeeping ---
  # Count the number of variants that passed all filters
  N_PASS=$(bcftools view -H "${SDIR}/variants.filtered.vcf.gz" | wc -l || echo 0)
  log "Summary ${SAMPLE}: ${N_PASS} passing variants"; hr

  # Calculate final sample timing and update the cohort-wide list of VCFs
  SAMPLE_T1=$(date +%s); SAMPLE_DT=$((SAMPLE_T1-SAMPLE_T0))
  AVG_SAMPLE=$(running_avg "${AVG_SAMPLE}" "${S_DONE}" "${SAMPLE_DT}"); S_DONE=$((S_DONE+1))
  SAMPLE_VCFS+=("${SDIR}/variants.filtered.vcf.gz")

done

# --- Cohort merge & exports ---
# Only perform cohort analysis if more than one sample was processed
if [[ ${#SAMPLE_VCFS[@]} -gt 1 ]]; then
  COHORT_DIR="${OUT_DIR}/cohort"; mkdir -p "${COHORT_DIR}"
  step_start "Merge cohort"
  # Combine all individual VCFs into a single cohort-wide file
  bcftools merge "${SAMPLE_VCFS[@]}" -Oz -o "${COHORT_DIR}/cohort.filtered.vcf.gz" \
    > "${COHORT_DIR}/merge.log" 2>&1 || log "WARN: merge failed"
  # If merge succeeded, index and export genotype/metric tables
  if [[ -f "${COHORT_DIR}/cohort.filtered.vcf.gz" ]]; then
    tabix -p vcf "${COHORT_DIR}/cohort.filtered.vcf.gz"
    step_done
    # Export a matrix of genotypes across all samples
    bcftools query -f'%CHROM\t%POS\t%ID[\t%GT]\n' "${COHORT_DIR}/cohort.filtered.vcf.gz" > "${COHORT_DIR}/cohort_genotypes.tsv" || true
    # Export site-level quality and depth metrics for the whole cohort
    bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL[\t%DP]\n' "${COHORT_DIR}/cohort.filtered.vcf.gz" > "${COHORT_DIR}/site_metrics.tsv" || true
  else
    step_done
  fi
fi

log "Pipeline completed. Outputs in: ${OUT_DIR}"
