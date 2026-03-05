#!/usr/bin/env bash
# Ion Torrent targeted variant calling pipeline: TVC -> normalise -> filter -> export.
# Processes each BAM in a manifest independently, then merges into a cohort-level VCF.
# Usage: ./05_variant_calling.sh -m MANIFEST -r REF.fa -b TARGETS.bed -o OUT_DIR [options]
set -euo pipefail

# --- Defaults ----------------------------------------------------------------
MANIFEST=""; REF_FA=""; BED=""; OUT_DIR=""; THREADS="${THREADS:-4}"
TVC_JSON=""; HOTSPOT_VCF=""; ERROR_MOTIF=""
MIN_DP="${MIN_DP:-50}"     # minimum per-sample read depth to retain a variant call
MIN_QUAL="${MIN_QUAL:-20}" # minimum Phred-scaled variant quality score

# --- Helper functions --------------------------------------------------------

_now(){ date '+%Y-%m-%d %H:%M:%S'; }
log(){  printf '[%s] %s\n' "$(_now)" "$*"; }

# Convert a duration in seconds to HH:MM:SS
fmt_dur(){ local s="$1"; printf '%02d:%02d:%02d' $((s/3600)) $(((s%3600)/60)) $((s%60)); }

# Compute an absolute ETA timestamp given remaining seconds
eta_ts(){ local s="$1"; date -d "@$(( $(date +%s) + s ))" '+%Y-%m-%d %H:%M:%S'; }

# Print a full-width horizontal rule
hr(){ printf '%*s' "${COLUMNS:-100}" '' | tr ' ' '-'; echo; }

# ASCII progress bar: takes a percentage (0-100) and prints in place
progress_bar(){
  local pct="$1" width=30
  local filled=$(( pct * width / 100 ))
  local empty=$(( width - filled ))
  printf '\r[%s%s] %3d%%' \
    "$(printf '%0.s#' $(seq 1 $filled))" \
    "$(printf '%0.s-' $(seq 1 $empty))" \
    "${pct}"
}

# Incremental running average: avoids storing all timings; returns updated mean
# Args: current_avg  n_samples_so_far  new_value
running_avg(){
  local avg="$1" n="$2" new="$3"
  avg=${avg//[^0-9]/}; n=${n//[^0-9]/}; new=${new//[^0-9]/}
  echo $(( (avg * n + new) / (n + 1) ))
}

# Step timing: call step_start before a block and step_done after; step_done
# prints elapsed time and returns the duration in seconds for ETA tracking
STEP_NAME=""; STEP_T0=0
step_start(){ STEP_NAME="$1"; STEP_T0=$(date +%s); log ">> ${STEP_NAME}..."; }
step_done(){
  local t1=$(date +%s)
  local dt=$(( t1 - STEP_T0 ))
  log "-- ${STEP_NAME} done ($(fmt_dur ${dt}))"
  echo "${dt}"
}

# On any unhandled error, report which sample was being processed before exiting
CURRENT_SAMPLE=""
trap 'echo; log "ERROR: pipeline failed while processing sample: ${CURRENT_SAMPLE:-<none>}"; exit 1' ERR

# --- Usage -------------------------------------------------------------------
usage(){ cat <<'USAGE'
Usage: 05_variant_calling.sh -m MANIFEST -r REF.fa -b TARGETS.bed -o OUT_DIR [options]

Required:
  -m FILE          Manifest: one BAM path per line (lines beginning with # are ignored)
  -r FILE          Reference FASTA (indexed .fai will be created if absent)
  -b FILE          BED file defining target regions
  -o DIR           Output directory

Options:
  -t INT           Threads (default: env THREADS or 4)
  --min-dp INT     Minimum per-sample depth to retain a call (default: env MIN_DP or 50)
  --min-qual INT   Minimum QUAL score (default: env MIN_QUAL or 20)
  --tvc-json FILE  TVC parameters JSON (optional; uses TVC defaults if omitted)
  --hotspot-vcf FILE  Hotspot/whitelist VCF passed to TVC (optional)
  --error-motif FILE  Error motif file for TVC (optional)
  -h, --help       Show this message
USAGE
}

# --- Argument parsing --------------------------------------------------------
while (( $# > 0 )); do
  case "$1" in
    -m) MANIFEST="$2"; shift 2 ;;
    -r) REF_FA="$2";   shift 2 ;;
    -b) BED="$2";      shift 2 ;;
    -o) OUT_DIR="$2";  shift 2 ;;
    -t) THREADS="$2";  shift 2 ;;
    --min-dp)      MIN_DP="$2";      shift 2 ;;
    --min-qual)    MIN_QUAL="$2";    shift 2 ;;
    --tvc-json)    TVC_JSON="$2";    shift 2 ;;
    --hotspot-vcf) HOTSPOT_VCF="$2"; shift 2 ;;
    --error-motif) ERROR_MOTIF="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

# --- Input validation --------------------------------------------------------
if [[ -z "${MANIFEST}" || -z "${REF_FA}" || -z "${BED}" || -z "${OUT_DIR}" ]]; then
  usage; exit 1
fi

# Verify all required command-line tools are available before starting
for tool in bcftools samtools tabix bgzip gzip awk sed grep; do
  command -v "${tool}" >/dev/null || { echo "ERROR: ${tool} not found in PATH" >&2; exit 1; }
done

# TVC can be invoked either as the Python wrapper or the compiled binary
if ! command -v tvc >/dev/null && ! command -v variant_caller_pipeline.py >/dev/null; then
  echo "ERROR: neither tvc nor variant_caller_pipeline.py found in PATH" >&2; exit 1
fi

[[ -f "${REF_FA}"   ]] || { echo "ERROR: reference FASTA not found: ${REF_FA}" >&2;   exit 1; }
[[ -f "${BED}"      ]] || { echo "ERROR: BED file not found: ${BED}" >&2;             exit 1; }
[[ -f "${MANIFEST}" ]] || { echo "ERROR: manifest not found: ${MANIFEST}" >&2;        exit 1; }

# --- Initialise output directories and log run parameters -------------------
mkdir -p "${OUT_DIR}" "${OUT_DIR}/logs"
log "Output dir:  ${OUT_DIR}"
log "Manifest:    ${MANIFEST}"
log "Reference:   ${REF_FA}"
log "Targets BED: ${BED}"
log "Threads:     ${THREADS}"
log "Filters:     MIN_DP=${MIN_DP}, MIN_QUAL=${MIN_QUAL}"
hr

# --- Index reference if needed ----------------------------------------------
# A .fai index is required by bcftools norm and TVC for random-access to the FASTA
if [[ ! -f "${REF_FA}.fai" ]]; then
  step_start "Indexing reference FASTA"
  samtools faidx "${REF_FA}"
  step_done >/dev/null
fi

# --- Normalise BED file -----------------------------------------------------
BED_NORM="${OUT_DIR}/targets.cleaned.bed"
step_start "Normalising BED -> ${BED_NORM}"
# Strip UCSC track/browser headers and Windows line endings; retain up to 6 columns
awk 'BEGIN{OFS="\t"} !/^track/ && !/^browser/ && !/^#/ { if(NF>=3){print $1,$2,$3,$4,$5,$6} }' "${BED}" | \
  sed 's/\r$//' > "${BED_NORM}"
step_done >/dev/null

# --- Load manifest ----------------------------------------------------------
# Strip comment lines and blank lines; result is an array of BAM paths
mapfile -t BAMS < <(grep -v -E '^[[:space:]]*#' "${MANIFEST}" | sed '/^[[:space:]]*$/d')
[[ ${#BAMS[@]} -gt 0 ]] || { echo "ERROR: manifest contains no valid BAM paths" >&2; exit 1; }
TOTAL=${#BAMS[@]}
log "Samples to process: ${TOTAL}"
hr

# --- Per-sample ETA tracking state ------------------------------------------
# Running averages are maintained per pipeline stage to give accurate ETAs
# as the run progresses across samples
SAMPLE_VCFS=()
IDX=0; S_DONE=0; AVG_SAMPLE=0
AVG_TVC=0; N_TVC=0
AVG_NORM=0; N_NORM=0
AVG_FILTER=0; N_FILTER=0

# =============================================================================
# Per-sample processing loop
# =============================================================================
for BAM in "${BAMS[@]}"; do
  IDX=$((IDX+1))
  if [[ ! -f "${BAM}" ]]; then
    log "WARNING: BAM not found, skipping: ${BAM}"; continue
  fi

  SAMPLE=$(basename "$BAM" .bam)
  # Sanitise sample name for use as a directory name (replace special characters)
  SAFE_SAMPLE="$(printf '%s' "${SAMPLE}" | tr -c 'A-Za-z0-9._-' '_')"
  SDIR="${OUT_DIR}/${SAFE_SAMPLE}"
  mkdir -p "${SDIR}/logs"

  # TVC requires the BAM filename to contain only safe characters; if the original
  # filename contains spaces or special characters, create a symlink with a safe name
  bam_dir=$(dirname -- "${BAM}")
  bam_base=$(basename -- "${BAM}")
  safe_base="$(printf '%s' "${bam_base}" | tr -c 'A-Za-z0-9._-' '_')"
  SAFE_BAM="${bam_dir}/${safe_base}"
  if [[ "${SAFE_BAM}" != "${BAM}" && ! -e "${SAFE_BAM}" ]]; then
    ln -s "${BAM}" "${SAFE_BAM}"
    # Also symlink the index file so TVC can find it alongside the safe-named BAM
    if   [[ -f "${BAM}.bai" ]];       then ln -s "${BAM}.bai"       "${SAFE_BAM}.bai" 2>/dev/null || true
    elif [[ -f "${BAM%.bam}.bai" ]];  then ln -s "${BAM%.bam}.bai"  "${SAFE_BAM}.bai" 2>/dev/null || true
    fi
  fi
  TVC_BAM="${SAFE_BAM}"

  # Index the BAM if no .bai exists (required for TVC random access)
  if [[ ! -f "${TVC_BAM}.bai" ]]; then
    step_start "Indexing BAM (${SAMPLE})"
    samtools index -@ "${THREADS}" "${TVC_BAM}"
    step_done >/dev/null
  fi

  echo; log "===== Sample ${IDX}/${TOTAL}: ${SAMPLE} ====="
  if [[ ${S_DONE} -gt 0 ]]; then
    REMAINING=$(( TOTAL - (IDX - 1) ))
    est_secs=$(( AVG_SAMPLE * REMAINING ))
    log "Estimated time remaining: ~$(fmt_dur ${est_secs}) (finish ~$(eta_ts ${est_secs}))"
  fi
  progress_bar 0; printf '\n'
  SAMPLE_T0=$(date +%s)

  # --- Step 1: Variant calling with TVC ------------------------------------
  # TVC (Torrent Variant Caller) is Thermo Fisher's Ion Torrent-specific caller;
  # it models platform error profiles (homopolymer indels) that generic callers
  # such as GATK do not account for. The pipeline supports both the Python wrapper
  # (variant_caller_pipeline.py) and the compiled binary (tvc), preferring the
  # wrapper when available as it handles more pre/post-processing automatically.
  step_start "Variant calling (TVC)"
  RAW_VCF=""
  if command -v variant_caller_pipeline.py >/dev/null; then
    TVC_CMD=(variant_caller_pipeline.py
      --input-bam        "${TVC_BAM}"
      --reference-fasta  "${REF_FA}"
      --region-bed       "${BED_NORM}"
      --output-dir       "${SDIR}"
      --num-threads      "${THREADS}")
    [[ -n "${TVC_JSON}"    ]] && TVC_CMD+=(--parameters-file "${TVC_JSON}")
    [[ -n "${HOTSPOT_VCF}" ]] && TVC_CMD+=(--input-vcf       "${HOTSPOT_VCF}")
    [[ -n "${ERROR_MOTIF}" ]] && TVC_CMD+=(--error-motif     "${ERROR_MOTIF}")
    # Log the exact command used for reproducibility
    echo "CMD: ${TVC_CMD[*]}" > "${SDIR}/logs/tvc.cmd.log"
    "${TVC_CMD[@]}" > "${SDIR}/logs/tvc.stdout.log" 2> "${SDIR}/logs/tvc.stderr.log"
    RAW_VCF="${SDIR}/TSVC_variants.vcf"
    [[ ! -f "${RAW_VCF}" ]] && RAW_VCF="${SDIR}/variants.vcf"   # fallback output name
  else
    TVC_CMD=(tvc
      --input-bam   "${TVC_BAM}"
      --reference   "${REF_FA}"
      --output-dir  "${SDIR}"
      --num-threads "${THREADS}"
      --target-file "${BED_NORM}"
      --output-vcf  "small_variants.vcf")
    echo "CMD: ${TVC_CMD[*]}" > "${SDIR}/logs/tvc.cmd.log"
    "${TVC_CMD[@]}" > "${SDIR}/logs/tvc.stdout.log" 2> "${SDIR}/logs/tvc.stderr.log"
    RAW_VCF="${SDIR}/TSVC_variants.vcf"
    [[ ! -f "${RAW_VCF}" ]] && RAW_VCF="${SDIR}/small_variants.vcf"
  fi
  [[ ! -f "${RAW_VCF}" ]] && { log "ERROR: no VCF produced by TVC"; exit 1; }
  tvc_dt=$(step_done); progress_bar 33; printf '\n'
  AVG_TVC=$(running_avg "${AVG_TVC}" "${N_TVC}" "${tvc_dt}"); N_TVC=$((N_TVC+1))

  # --- Step 2: Validate and repair VCF header ------------------------------
  # Some TVC versions omit ##contig lines from the VCF header; bcftools norm
  # (used in step 3) requires these lines to resolve chromosome lengths.
  # Additionally, TVC may write chromosome names without the 'chr' prefix even
  # when the reference uses 'chr'-prefixed names, causing a mismatch that
  # prevents downstream tools from finding variants — this is corrected here.
  step_start "Validating/repairing VCF header"

  # Re-inject ##contig lines from the .fai index if they are absent
  if ! bcftools view -h "${RAW_VCF}" | grep -q '^##contig'; then
    log "VCF header missing ##contig lines — reinserting from reference index..."
    bcftools reheader --fai "${REF_FA}.fai" "${RAW_VCF}" -o "${SDIR}/variants.hdr.vcf"
    bgzip -f "${SDIR}/variants.hdr.vcf"
    tabix -p vcf "${SDIR}/variants.hdr.vcf.gz"
    RAW_VCF="${SDIR}/variants.hdr.vcf.gz"
  fi

  # Check whether the reference and VCF use the same chromosome naming convention
  HAS_CHR_FA=$(awk -F'\t' 'NR==1{print ($1 ~ /^chr/)?1:0}' "${REF_FA}.fai")
  FIRST_VAR=$(bcftools view -H "${RAW_VCF}" | head -n1 || true)
  HAS_CHR_VCF=0
  [[ -n "${FIRST_VAR}" && "${FIRST_VAR%%	*}" =~ ^chr ]] && HAS_CHR_VCF=1

  if [[ -n "${FIRST_VAR}" ]]; then
    if [[ "${HAS_CHR_FA}" -eq 1 && "${HAS_CHR_VCF}" -eq 0 ]]; then
      # Reference uses 'chr1' style but VCF uses '1' style — add prefix to VCF
      log "Chromosome naming mismatch: FASTA uses 'chr*', VCF does not. Renaming VCF contigs..."
      awk -F'\t' '{g=$1; sub(/^chr/,"",g); print g "\t" $1}' "${REF_FA}.fai" > "${SDIR}/chr_map.txt"
      bcftools annotate --rename-chrs "${SDIR}/chr_map.txt" "${RAW_VCF}" \
        -Oz -o "${SDIR}/variants.renamed.vcf.gz"
      tabix -p vcf "${SDIR}/variants.renamed.vcf.gz"
      RAW_VCF="${SDIR}/variants.renamed.vcf.gz"
    elif [[ "${HAS_CHR_FA}" -eq 0 && "${HAS_CHR_VCF}" -eq 1 ]]; then
      # Reference uses '1' style but VCF uses 'chr1' style — strip prefix from VCF
      log "Chromosome naming mismatch: VCF uses 'chr*', FASTA does not. Stripping 'chr' prefix..."
      awk -F'\t' '{g=$1; print $1 "\t" gensub(/^chr/,"","",g)}' "${REF_FA}.fai" > "${SDIR}/chr_map.txt"
      bcftools annotate --rename-chrs "${SDIR}/chr_map.txt" "${RAW_VCF}" \
        -Oz -o "${SDIR}/variants.renamed.vcf.gz"
      tabix -p vcf "${SDIR}/variants.renamed.vcf.gz"
      RAW_VCF="${SDIR}/variants.renamed.vcf.gz"
    else
      log "Chromosome naming consistent between FASTA and VCF — no renaming needed."
    fi
  else
    log "VCF contains no variants — skipping chromosome rename."
  fi

  # Final check: abort rather than pass a malformed header to downstream tools
  if ! bcftools view -h "${RAW_VCF}" | grep -q '^##contig'; then
    log "ERROR: VCF header still missing ##contig lines after repair — aborting to prevent corruption."
    exit 1
  fi
  step_done

  # --- Step 3: Left-align and split multi-allelic sites --------------------
  # bcftools norm left-aligns indels against the reference (canonical representation
  # required for consistent variant matching across samples) and splits
  # multi-allelic records into separate lines (-m -both) so that each ALT allele
  # can be filtered and annotated independently in downstream steps.
  step_start "Left-normalise and compress"
  bcftools norm -f "${REF_FA}" -m -both "${RAW_VCF}" \
    -Oz -o "${SDIR}/variants.norm.vcf.gz"
  tabix -p vcf "${SDIR}/variants.norm.vcf.gz"
  norm_dt=$(step_done); progress_bar 66; printf '\n'
  AVG_NORM=$(running_avg "${AVG_NORM}" "${N_NORM}" "${norm_dt}"); N_NORM=$((N_NORM+1))

  # --- Step 4: Hard filter -------------------------------------------------
  # Retains only calls that:
  #   - carry the PASS filter flag (TVC internal quality thresholds met)
  #   - have depth >= MIN_DP in at least one sample (MAX across FORMAT/DP)
  #   - have a Phred-scaled quality score >= MIN_QUAL
  # Using MAX(FMT/DP) rather than INFO/DP ensures the threshold applies per-sample,
  # which is important in multi-sample VCFs where INFO/DP reflects pooled depth.
  step_start "Hard filtering (PASS; DP>=${MIN_DP}; QUAL>=${MIN_QUAL})"
  bcftools view \
    -i "MAX(FMT/DP)>=${MIN_DP} && QUAL>=${MIN_QUAL} && FILTER='PASS'" \
    "${SDIR}/variants.norm.vcf.gz" \
    -Oz -o "${SDIR}/variants.filtered.vcf.gz"
  tabix -p vcf "${SDIR}/variants.filtered.vcf.gz"
  filter_dt=$(step_done); progress_bar 100; printf '\n'
  AVG_FILTER=$(running_avg "${AVG_FILTER}" "${N_FILTER}" "${filter_dt}"); N_FILTER=$((N_FILTER+1))

  # --- Step 5: Export tables for manual review ----------------------------
  # The normalised (pre-filter) VCF is used here so that borderline calls just
  # below the filter thresholds remain visible for manual inspection.
  # Helper functions check the VCF header before querying a tag — not all TVC
  # versions emit the same FORMAT/INFO fields.
  has_info_tag(){ bcftools view -h "$1" | grep -q "^##INFO=<ID=$2,"; }
  has_fmt_tag(){  bcftools view -h "$1" | grep -q "^##FORMAT=<ID=$2,"; }
  VCF_HDR="${SDIR}/variants.norm.vcf.gz"
  # Site-level bias indicators that indicate potential false positives when elevated
  INFO_TAGS=(HP Indel STB STBP HS SB)

  # Per-sample export: allele frequency is derived from AO/RO (alt/ref observations)
  # when available; falls back to FORMAT/AF if present, or '.' if neither is in the header
  step_start "Exporting per-sample review table"
  AO_RO_FMT=0; has_fmt_tag "${VCF_HDR}" AO && has_fmt_tag "${VCF_HDR}" RO && AO_RO_FMT=1 || true
  AF_FMT=0;    has_fmt_tag "${VCF_HDR}" AF && AF_FMT=1 || true

  if [[ "${AO_RO_FMT}" -eq 1 ]]; then
    # AF = AO / (AO + RO): proportion of reads supporting the alternate allele
    bcftools query -f "[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%SAMPLE\t%DP\t%AO\t%RO\n]" \
      "${VCF_HDR}" > "${SDIR}/manual_review_samples.raw.tsv" || true
    awk -F'\t' 'BEGIN{OFS="\t"} {
      ao=$9; ro=$10
      if (ao ~ /^[0-9]+$/ && ro ~ /^[0-9]+$/ && (ao+ro)>0) { af = ao/(ao+ro) } else { af = "." }
      print $1,$2,$3,$4,$5,$6,$7,$8,af
    }' "${SDIR}/manual_review_samples.raw.tsv" > "${SDIR}/manual_review_samples.tsv" || true
    rm -f "${SDIR}/manual_review_samples.raw.tsv"
  elif [[ "${AF_FMT}" -eq 1 ]]; then
    bcftools query -f "[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%SAMPLE\t%DP\t%AF\n]" \
      "${VCF_HDR}" > "${SDIR}/manual_review_samples.tsv" || true
  else
    # No AF-related fields available — output a placeholder column
    bcftools query -f "[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%SAMPLE\t%DP\n]" \
      "${VCF_HDR}" > "${SDIR}/manual_review_samples.tmp.tsv" || true
    awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,"."}' \
      "${SDIR}/manual_review_samples.tmp.tsv" > "${SDIR}/manual_review_samples.tsv" || true
    rm -f "${SDIR}/manual_review_samples.tmp.tsv"
  fi
  step_done

  # Site-level export: include bias tags only when present in the header to
  # avoid bcftools errors on VCFs that lack them; absent tags are written as '.'
  step_start "Exporting site-level review table"
  FMT_STR="%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER"
  if bcftools view -h "${VCF_HDR}" | grep -q '^##INFO=<ID=DP,'; then FMT_STR+="\t%INFO/DP"; else FMT_STR+="\t."; fi
  if bcftools view -h "${VCF_HDR}" | grep -q '^##INFO=<ID=AF,'; then FMT_STR+="\t%INFO/AF"; else FMT_STR+="\t."; fi
  for tag in "${INFO_TAGS[@]}"; do
    if bcftools view -h "${VCF_HDR}" | grep -q "^##INFO=<ID=${tag},"; then
      FMT_STR+="\t%INFO/${tag}"
    else
      FMT_STR+="\t."
    fi
  done
  FMT_STR+="\n"
  bcftools query -f "${FMT_STR}" "${VCF_HDR}" > "${SDIR}/manual_review_site.tsv" || true
  step_done

  # --- Per-sample summary --------------------------------------------------
  N_PASS=$(bcftools view -H "${SDIR}/variants.filtered.vcf.gz" | wc -l || echo 0)
  log "Summary — ${SAMPLE}: ${N_PASS} variants passing all filters"
  hr

  SAMPLE_T1=$(date +%s); SAMPLE_DT=$((SAMPLE_T1-SAMPLE_T0))
  AVG_SAMPLE=$(running_avg "${AVG_SAMPLE}" "${S_DONE}" "${SAMPLE_DT}"); S_DONE=$((S_DONE+1))
  SAMPLE_VCFS+=("${SDIR}/variants.filtered.vcf.gz")

done

# =============================================================================
# Cohort-level merge
# =============================================================================
# When more than one sample has been processed, merge all filtered VCFs into a
# single multi-sample VCF. This enables cohort-wide genotype comparisons and
# is required as input for downstream annotation (script 06) and mapping (script 08).
if [[ ${#SAMPLE_VCFS[@]} -gt 1 ]]; then
  COHORT_DIR="${OUT_DIR}/cohort"
  mkdir -p "${COHORT_DIR}"
  step_start "Merging cohort VCFs"
  bcftools merge "${SAMPLE_VCFS[@]}" \
    -Oz -o "${COHORT_DIR}/cohort.filtered.vcf.gz" \
    > "${COHORT_DIR}/merge.log" 2>&1 || log "WARNING: cohort merge failed — check ${COHORT_DIR}/merge.log"

  if [[ -f "${COHORT_DIR}/cohort.filtered.vcf.gz" ]]; then
    tabix -p vcf "${COHORT_DIR}/cohort.filtered.vcf.gz"
    step_done
    # Export genotype matrix: one row per site, one column per sample
    bcftools query -f '%CHROM\t%POS\t%ID[\t%GT]\n' \
      "${COHORT_DIR}/cohort.filtered.vcf.gz" > "${COHORT_DIR}/cohort_genotypes.tsv" || true
    # Export site-level depth summary across samples
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL[\t%DP]\n' \
      "${COHORT_DIR}/cohort.filtered.vcf.gz" > "${COHORT_DIR}/site_metrics.tsv" || true
  else
    step_done
  fi
fi

log "Pipeline complete. Outputs in: ${OUT_DIR}"
