
#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------------
# Ion Torrent targeted variant calling pipeline (TVC → normalise → filter → export)
# Corrected & hardened version
# ------------------------------------------------------------------------------------

# --------------------------- Defaults ---------------------------
MANIFEST=""; REF_FA=""; BED=""; OUT_DIR=""; THREADS="${THREADS:-4}"
TVC_JSON=""; HOTSPOT_VCF=""; ERROR_MOTIF=""
MIN_DP="${MIN_DP:-50}"; MIN_QUAL="${MIN_QUAL:-20}"

# --------------------------- Helpers ----------------------------
_now(){ date '+%Y-%m-%d %H:%M:%S'; }
log(){ printf '[%s] %s\n' "$(_now)" "$*"; }
fmt_dur(){ local s="$1"; printf '%02d:%02d:%02d' $((s/3600)) $(((s%3600)/60)) $((s%60)); }
eta_ts(){ local s="$1"; date -d "@$(( $(date +%s) + s ))" '+%Y-%m-%d %H:%M:%S'; }
hr(){ printf '%*s' "${COLUMNS:-100}" '' | tr ' ' '-'; echo; }
progress_bar(){ local pct="$1"; local width=30; local filled=$(( pct * width / 100 )); local empty=$(( width - filled ));
  printf '\r[%s%s] %3d%%' "$(printf '%0.s#' $(seq 1 $filled))" "$(printf '%0.s-' $(seq 1 $empty))" "${pct}"; }
# Running average
running_avg(){ local avg="$1" n="$2" new="$3"; avg=${avg//[^0-9]/}; n=${n//[^0-9]/}; new=${new//[^0-9]/}; echo $(( (avg * n + new) / (n + 1) )); }

STEP_NAME=""; STEP_T0=0
step_start(){ STEP_NAME="$1"; STEP_T0=$(date +%s); log "▶ ${STEP_NAME}…"; }
step_done(){ local t1=$(date +%s); local dt=$(( t1 - STEP_T0 )); log "✔ ${STEP_NAME} done ($(fmt_dur ${dt}))"; echo "${dt}"; }

CURRENT_SAMPLE=""
trap 'echo; log "✖ Mistake while processing sample: ${CURRENT_SAMPLE:-<none>}"; exit 1' ERR

usage(){ cat <<'USAGE'
Usage: 05_variant_calling.sh -m MANIFEST -r REF.fa -b TARGETS.bed -o OUT_DIR [options]

Required arguments:
  -m FILE              Manifest with one BAM path per line (comments with # allowed)
  -r FILE              Reference FASTA (indexed with .fai; will be created if missing)
  -b FILE              BED targets file
  -o DIR               Output directory

Options:
  -t INT               Threads (default: env THREADS or 4)
  --min-dp INT         Minimum per-sample depth (default: env MIN_DP or 50)
  --min-qual INT       Minimum QUAL (default: env MIN_QUAL or 20)
  --tvc-json FILE      TVC parameters JSON (optional)
  --hotspot-vcf FILE   Hotspot/whitelist VCF to feed TVC (optional)
  --error-motif FILE   Error motif file for TVC (optional)
  -h, --help           Show this message
USAGE
}

# --------------------------- Parse args -------------------------
while (( $# > 0 )); do
  case "$1" in
    -m) MANIFEST="$2"; shift 2 ;;
    -r) REF_FA="$2";  shift 2 ;;
    -b) BED="$2";     shift 2 ;;
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

# --------------------------- Basic checks -----------------------
if [[ -z "${MANIFEST}" || -z "${REF_FA}" || -z "${BED}" || -z "${OUT_DIR}" ]]; then
  usage; exit 1
fi
for tool in bcftools samtools tabix bgzip gzip awk sed grep; do
  command -v "${tool}" >/dev/null || { echo "ERROR: ${tool} not found"; exit 1; }
done
if ! command -v tvc >/dev/null && ! command -v variant_caller_pipeline.py >/dev/null; then
  echo "ERROR: TVC or wrapper not found"; exit 1
fi
[[ -f "${REF_FA}" ]] || { echo "ERROR: reference FASTA not found"; exit 1; }
[[ -f "${BED}"    ]] || { echo "ERROR: BED not found"; exit 1; }
[[ -f "${MANIFEST}" ]] || { echo "ERROR: manifest not found"; exit 1; }

# --------------------------- Log configuration ------------------
mkdir -p "${OUT_DIR}" "${OUT_DIR}/logs"
log "Output dir: ${OUT_DIR}"
log "Manifest:   ${MANIFEST}"
log "Reference:  ${REF_FA}"
log "Targets BED:${BED}"
log "Threads:    ${THREADS}"
log "Filters:    MIN_DP=${MIN_DP}, MIN_QUAL=${MIN_QUAL}"
hr

# --------------------------- Reference index --------------------
if [[ ! -f "${REF_FA}.fai" ]]; then
  step_start "Index reference"
  samtools faidx "${REF_FA}"
  step_done >/dev/null
fi

# --------------------------- BED normalisation ------------------
BED_NORM="${OUT_DIR}/targets.cleaned.bed"
step_start "Normalize BED -> ${BED_NORM}"
# Remove track/browser/comment lines, strip CRs, keep up to 6 cols
awk 'BEGIN{OFS="\t"} !/^track/ && !/^browser/ && !/^#/ { if(NF>=3){print $1,$2,$3,$4,$5,$6} }' "${BED}" | \
  sed 's/\r$//' > "${BED_NORM}"
step_done >/dev/null

# --------------------------- Load manifest ----------------------
mapfile -t BAMS < <(grep -v -E '^[[:space:]]*#' "${MANIFEST}" | sed '/^[[:space:]]*$/d')
[[ ${#BAMS[@]} -gt 0 ]] || { echo "ERROR: empty manifest"; exit 1; }
TOTAL=${#BAMS[@]}
log "Samples to process: ${TOTAL}"
hr

# --------------------------- ETA state --------------------------
SAMPLE_VCFS=(); IDX=0; S_DONE=0; AVG_SAMPLE=0; AVG_TVC=0; N_TVC=0; AVG_NORM=0; N_NORM=0; AVG_FILTER=0; N_FILTER=0

# --------------------------- Per-sample loop --------------------
for BAM in "${BAMS[@]}"; do
  IDX=$((IDX+1))
  if [[ ! -f "${BAM}" ]]; then
    log "WARN: BAM missing, skipping: ${BAM}"; continue
  fi
  
  SAMPLE=$(basename "$BAM" .bam)
  SAFE_SAMPLE="$(printf '%s' "${SAMPLE}" | tr -c 'A-Za-z0-9._-' '_')"
  SDIR="${OUT_DIR}/${SAFE_SAMPLE}"

  mkdir -p "${SDIR}/logs"

  # ---- sanitise BAM filename (basename only) ----
  bam_dir=$(dirname -- "${BAM}")
  bam_base=$(basename -- "${BAM}")
  safe_base="$(printf '%s' "${bam_base}" | tr -c 'A-Za-z0-9._-' '_')"
  SAFE_BAM="${bam_dir}/${safe_base}"
  if [[ "${SAFE_BAM}" != "${BAM}" && ! -e "${SAFE_BAM}" ]]; then
    ln -s "${BAM}" "${SAFE_BAM}"
    if   [[ -f "${BAM}.bai" ]]; then ln -s "${BAM}.bai" "${SAFE_BAM}.bai" 2>/dev/null || true
    elif [[ -f "${BAM%.bam}.bai" ]]; then ln -s "${BAM%.bam}.bai" "${SAFE_BAM}.bai" 2>/dev/null || true
    fi
  fi
  TVC_BAM="${SAFE_BAM}"

  # Ensure the BAM we will use has its index
  if [[ ! -f "${TVC_BAM}.bai" ]]; then
    step_start "Index BAM (${SAMPLE})"
    samtools index -@ "${THREADS}" "${TVC_BAM}"
    step_done >/dev/null
  fi

  echo; log "===== Sample ${IDX}/${TOTAL}: ${SAMPLE} ====="
  if [[ ${S_DONE} -gt 0 ]]; then
    REMAINING=$((TOTAL-(IDX-1)))
    est_secs=$((AVG_SAMPLE*REMAINING))
    log "ETA before start: ~$(fmt_dur ${est_secs}) (ETA $(eta_ts ${est_secs}))"
  fi
  log "Collecting timings on first sample…"
  progress_bar 0; printf '\n'; SAMPLE_T0=$(date +%s)

  # ----------------- TVC call -----------------
  step_start "TVC (variant calling)"
  RAW_VCF=""
  if command -v variant_caller_pipeline.py >/dev/null; then
    TVC_CMD=(variant_caller_pipeline.py \
      --input-bam "${TVC_BAM}" \
      --reference-fasta "${REF_FA}" \
      --region-bed "${BED_NORM}" \
      --output-dir "${SDIR}" \
      --num-threads "${THREADS}")
    [[ -n "${TVC_JSON}"    ]] && TVC_CMD+=(--parameters-file "${TVC_JSON}")
    [[ -n "${HOTSPOT_VCF}" ]] && TVC_CMD+=(--input-vcf       "${HOTSPOT_VCF}")
    [[ -n "${ERROR_MOTIF}" ]] && TVC_CMD+=(--error-motif     "${ERROR_MOTIF}")
    echo "TVC WRAPPER CMD: ${TVC_CMD[*]}" > "${SDIR}/logs/tvc.cmd.log"
    "${TVC_CMD[@]}" > "${SDIR}/logs/tvc.stdout.log" 2> "${SDIR}/logs/tvc.stderr.log"
    RAW_VCF="${SDIR}/TSVC_variants.vcf"
    [[ ! -f "${RAW_VCF}" ]] && RAW_VCF="${SDIR}/variants.vcf"
  else
    TVC_CMD=(tvc \
      --input-bam "${TVC_BAM}" \
      --reference "${REF_FA}" \
      --output-dir "${SDIR}" \
      --num-threads "${THREADS}" \
      --target-file "${BED_NORM}" \
      --output-vcf "small_variants.vcf")
    echo "TVC CMD: ${TVC_CMD[*]}" > "${SDIR}/logs/tvc.cmd.log"
    "${TVC_CMD[@]}" > "${SDIR}/logs/tvc.stdout.log" 2> "${SDIR}/logs/tvc.stderr.log"
    RAW_VCF="${SDIR}/TSVC_variants.vcf"
    [[ ! -f "${RAW_VCF}" ]] && RAW_VCF="${SDIR}/small_variants.vcf"
  fi
  [[ ! -f "${RAW_VCF}" ]] && { log "ERROR: VCF not found after TVC"; exit 1; }
  tvc_dt=$(step_done); progress_bar 33; printf '\n'
  AVG_TVC=$(running_avg "${AVG_TVC}" "${N_TVC}" "${tvc_dt}"); N_TVC=$((N_TVC+1))

  # ------------- Validate/repair VCF header ----------------
  step_start "Validate/repair VCF header (reheader + optional rename)"
  # 1) Re-inject contigs if missing
  if ! bcftools view -h "${RAW_VCF}" | grep -q '^##contig'; then
    log "VCF header lacks ##contig. Reinserting contigs using ${REF_FA}.fai..."
    bcftools reheader --fai "${REF_FA}.fai" "${RAW_VCF}" -o "${SDIR}/variants.hdr.vcf" 
    bgzip -f "${SDIR}/variants.hdr.vcf"
    tabix -p vcf "${SDIR}/variants.hdr.vcf.gz"
    RAW_VCF="${SDIR}/variants.hdr.vcf.gz"
  fi
  # 2) Detect contig naming scheme (chr vs no chr); use .fai and first variant line
  HAS_CHR_FA=$(awk -F'\t' 'NR==1{print ($1 ~ /^chr/)?1:0}' "${REF_FA}.fai")
  FIRST_VAR=$(bcftools view -H "${RAW_VCF}" | head -n1 || true)
  HAS_CHR_VCF=0; [[ -n "${FIRST_VAR}" && "${FIRST_VAR%%\t*}" =~ ^chr ]] && HAS_CHR_VCF=1
  # 3) Rename if mismatch
  if [[ -n "${FIRST_VAR}" ]]; then
    if [[ "${HAS_CHR_FA}" -eq 1 && "${HAS_CHR_VCF}" -eq 0 ]]; then
      log "FASTA usa 'chr*' y VCF no. Renombrando contigs en VCF -> 'chr*'..."
      awk -F'\t' '{g=$1; sub(/^chr/,"",g); print g "\t" $1}' "${REF_FA}.fai" > "${SDIR}/chr_map.txt"  # 1 -> chr1
      bcftools annotate --rename-chrs "${SDIR}/chr_map.txt" "${RAW_VCF}" -Oz -o "${SDIR}/variants.renamed.vcf.gz"
      tabix -p vcf "${SDIR}/variants.renamed.vcf.gz"
      RAW_VCF="${SDIR}/variants.renamed.vcf.gz"
    elif [[ "${HAS_CHR_FA}" -eq 0 && "${HAS_CHR_VCF}" -eq 1 ]]; then
      log "FASTA no usa 'chr*' y VCF sí. Renombrando contigs en VCF -> sin 'chr'..."
      awk -F'\t' '{g=$1; print $1 "\t" gensub(/^chr/,"","",g)}' "${REF_FA}.fai" > "${SDIR}/chr_map.txt"  # chr1 -> 1
      bcftools annotate --rename-chrs "${SDIR}/chr_map.txt" "${RAW_VCF}" -Oz -o "${SDIR}/variants.renamed.vcf.gz"
      tabix -p vcf "${SDIR}/variants.renamed.vcf.gz"
      RAW_VCF="${SDIR}/variants.renamed.vcf.gz"
    else
      log "Nomenclatura de contigs coherente entre FASTA y VCF."
    fi
  else
    log "VCF sin variantes; saltamos rename-chrs. (Header ya reinyectado si faltaba)"
  fi
  # 4) Final header validation
  if ! bcftools view -h "${RAW_VCF}" | grep -q '^##contig'; then
    log "ERROR: header aún sin ##contig tras reheader. Abortando para evitar corrupción."; exit 1
  fi
  step_done

  # ----------------- Normalize + compress -----------------
  step_start "Normalize + compress"
  bcftools norm -f "${REF_FA}" -m -both "${RAW_VCF}" -Oz -o "${SDIR}/variants.norm.vcf.gz"
  tabix -p vcf "${SDIR}/variants.norm.vcf.gz"
  norm_dt=$(step_done); progress_bar 66; printf '\n'
  AVG_NORM=$(running_avg "${AVG_NORM}" "${N_NORM}" "${norm_dt}"); N_NORM=$((N_NORM+1))

  # ----------------- Unified filtering --------------------
  step_start "Filter (PASS; DP>=${MIN_DP}; QUAL>=${MIN_QUAL})"
  bcftools view -i "MAX(FMT/DP)>=${MIN_DP} && QUAL>=${MIN_QUAL} && FILTER='PASS'" \
    "${SDIR}/variants.norm.vcf.gz" -Oz -o "${SDIR}/variants.filtered.vcf.gz"
  tabix -p vcf "${SDIR}/variants.filtered.vcf.gz"
  filter_dt=$(step_done); progress_bar 100; printf '\n'
  AVG_FILTER=$(running_avg "${AVG_FILTER}" "${N_FILTER}" "${filter_dt}"); N_FILTER=$((N_FILTER+1))

  # ----------------- Export manual review candidates --------
  # Helper funcs to detect header tags (take file arg)
  has_info_tag(){ bcftools view -h "$1" | grep -q "^##INFO=<ID=$2,"; }
  has_fmt_tag(){  bcftools view -h "$1" | grep -q "^##FORMAT=<ID=$2,"; }
  VCF_HDR="${SDIR}/variants.norm.vcf.gz"
  INFO_TAGS=(HP Indel STB STBP HS SB) # potential site-level bias tags

  # 1) Per-sample export (derive AF from AO/RO if present)
  step_start "Export manual review (per-sample metrics)"
  AO_RO_FMT=0; has_fmt_tag "${VCF_HDR}" AO && has_fmt_tag "${VCF_HDR}" RO && AO_RO_FMT=1 || true
  AF_FMT=0;    has_fmt_tag "${VCF_HDR}" AF && AF_FMT=1 || true
  if [[ "${AO_RO_FMT}" -eq 1 ]]; then
    bcftools query -f "[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%SAMPLE\t%DP\t%AO\t%RO\n]" \
      "${VCF_HDR}" > "${SDIR}/manual_review_samples.raw.tsv" || true
    awk -F'\t' 'BEGIN{OFS="\t"} {
      chrom=$1; pos=$2; ref=$3; alt=$4; qual=$5; filt=$6; sample=$7; dp=$8; ao=$9; ro=$10;
      if (ao ~ /^[0-9]+$/ && ro ~ /^[0-9]+$/ && (ao+ro)>0) { af = ao/(ao+ro) } else { af = "." }
      print chrom, pos, ref, alt, qual, filt, sample, dp, af
    }' "${SDIR}/manual_review_samples.raw.tsv" > "${SDIR}/manual_review_samples.tsv" || true
    rm -f "${SDIR}/manual_review_samples.raw.tsv"
  else
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

  # 2) Site-level export (INFO tags only if present)
  step_start "Export manual review (site-level INFO tags)"
  FMT_STR="%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER"
  if bcftools view -h "${VCF_HDR}" | grep -q '^##INFO=<ID=DP,'; then FMT_STR+="\t%INFO/DP"; else FMT_STR+="\t."; fi
  if bcftools view -h "${VCF_HDR}" | grep -q '^##INFO=<ID=AF,'; then FMT_STR+="\t%INFO/AF"; else FMT_STR+="\t."; fi
  for tag in "${INFO_TAGS[@]}"; do
    if bcftools view -h "${VCF_HDR}" | grep -q "^##INFO=<ID=${tag},"; then FMT_STR+="\t%INFO/${tag}"; else FMT_STR+="\t."; fi
  done
  FMT_STR+="\n"
  bcftools query -f "${FMT_STR}" "${VCF_HDR}" > "${SDIR}/manual_review_site.tsv" || true
  step_done

  # ----------------- Summary & bookkeeping -----------------
  N_PASS=$(bcftools view -H "${SDIR}/variants.filtered.vcf.gz" | wc -l || echo 0)
  log "Summary ${SAMPLE}: ${N_PASS} passing variants"
  hr

  SAMPLE_T1=$(date +%s); SAMPLE_DT=$((SAMPLE_T1-SAMPLE_T0))
  AVG_SAMPLE=$(running_avg "${AVG_SAMPLE}" "${S_DONE}" "${SAMPLE_DT}"); S_DONE=$((S_DONE+1))
  SAMPLE_VCFS+=("${SDIR}/variants.filtered.vcf.gz")

done

# --------------------------- Cohort merge & exports -------------
if [[ ${#SAMPLE_VCFS[@]} -gt 1 ]]; then
  COHORT_DIR="${OUT_DIR}/cohort"; mkdir -p "${COHORT_DIR}"
  step_start "Merge cohort"
  # Merge VCFs across samples
  bcftools merge "${SAMPLE_VCFS[@]}" -Oz -o "${COHORT_DIR}/cohort.filtered.vcf.gz" \
    > "${COHORT_DIR}/merge.log" 2>&1 || log "WARN: merge failed"
  if [[ -f "${COHORT_DIR}/cohort.filtered.vcf.gz" ]]; then
    tabix -p vcf "${COHORT_DIR}/cohort.filtered.vcf.gz"
    step_done
    bcftools query -f'%CHROM\t%POS\t%ID[\t%GT]\n' "${COHORT_DIR}/cohort.filtered.vcf.gz" > "${COHORT_DIR}/cohort_genotypes.tsv" || true
    bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL[\t%DP]\n' "${COHORT_DIR}/cohort.filtered.vcf.gz" > "${COHORT_DIR}/site_metrics.tsv" || true
  else
    step_done
  fi
fi

log "Pipeline completed. Outputs in: ${OUT_DIR}"
