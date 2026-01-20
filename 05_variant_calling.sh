
#!/usr/bin/env bash
set -euo pipefail
   
# -----------------------------------------------------------------------------
# Ion Torrent targeted variant calling pipeline (TVC → normalize → filter → export)
# -----------------------------------------------------------------------------

MANIFEST="" ; REF_FA="" ; BED="" ; OUT_DIR="" ; THREADS="${THREADS:-4}"
TVC_JSON="" ; HOTSPOT_VCF="" ; ERROR_MOTIF=""
MIN_DP="${MIN_DP:-50}" ; MIN_QUAL="${MIN_QUAL:-20}"

usage(){
  cat <<'EOF'
Usage: 03_variant_calling.sh -m PASS_MANIFEST -r REF_FASTA -b TARGETS_BED -o OUT_DIR
       [-t THREADS] [--min-dp N] [--min-qual Q] [--tvc-json FILE]
       [--hotspot-vcf FILE] [--error-motif FILE]
EOF
}

# ---------- Helpers ----------
_now() { date '+%Y-%m-%d %H:%M:%S'; }
log(){ printf '[%s] %s\n' "$(_now)" "$*"; }
fmt_dur(){ local s="$1"; printf '%02d:%02d:%02d' $((s/3600)) $(((s%3600)/60)) $((s%60)); }
eta_ts(){ local s="$1"; date -d "@$(( $(date +%s) + s ))" '+%Y-%m-%d %H:%M:%S'; }
hr(){ printf '%*s\n' "${COLUMNS:-100}" '' | tr ' ' '-'; }
progress_bar(){ local pct="$1"; local width=30; local filled=$(( pct * width / 100 )); local empty=$(( width - filled )); printf '\r[%s%s] %3d%%' "$(printf '%0.s#' $(seq 1 $filled))" "$(printf '%0.s-' $(seq 1 $empty))" "${pct}"; }

# Media corrida robusta (solo números)
running_avg(){
  local avg="$1" n="$2" new="$3"
  avg=${avg//[^0-9]/}; n=${n//[^0-9]/}; new=${new//[^0-9]/}
  echo $(( (avg * n + new) / (n + 1) ))
}

STEP_NAME=""; STEP_T0=0
step_start(){ STEP_NAME="$1"; STEP_T0=$(date +%s); log "▶ ${STEP_NAME}…"; }
step_done(){
  local t1=$(date +%s); local dt=$(( t1 - STEP_T0 ))
  log "✔ ${STEP_NAME} done ($(fmt_dur ${dt}))"
  echo "${dt}"
}

CURRENT_SAMPLE=""
trap 'echo; log "✖ Error mientras procesaba la muestra: ${CURRENT_SAMPLE:-<none>}"; exit 1' ERR

# ---------- Parse args ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -m) MANIFEST="$2"; shift 2 ;;
    -r) REF_FA="$2"; shift 2 ;;
    -b) BED="$2"; shift 2 ;;
    -o) OUT_DIR="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    --min-dp) MIN_DP="$2"; shift 2 ;;
    --min-qual) MIN_QUAL="$2"; shift 2 ;;
    --tvc-json) TVC_JSON="$2"; shift 2 ;;
    --hotspot-vcf) HOTSPOT_VCF="$2"; shift 2 ;;
    --error-motif) ERROR_MOTIF="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 2 ;;
  esac
done

# ---------- Basic checks ----------
[[ -z "${MANIFEST}" || -z "${REF_FA}" || -z "${BED}" || -z "${OUT_DIR}" ]] && { usage; exit 1; }

for tool in bcftools samtools tabix gzip bgzip; do
  command -v "${tool}" >/dev/null || { echo "ERROR: ${tool} not found"; exit 1; }
done
command -v tvc >/dev/null || command -v variant_caller_pipeline.py >/dev/null || { echo "ERROR: TVC or wrapper not found"; exit 1; }

[[ -f "${REF_FA}" ]] || { echo "ERROR: reference FASTA not found"; exit 1; }
[[ -f "${BED}" ]] || { echo "ERROR: BED not found"; exit 1; }
[[ -f "${MANIFEST}" ]] || { echo "ERROR: manifest not found"; exit 1; }

# ---------- Log configuration ----------
mkdir -p "${OUT_DIR}" "${OUT_DIR}/logs"
log "Output dir: ${OUT_DIR}"
log "Manifest:   ${MANIFEST}"
log "Reference:  ${REF_FA}"
log "Targets BED:${BED}"
log "Threads:    ${THREADS}"
log "Filters:    MIN_DP=${MIN_DP}, MIN_QUAL=${MIN_QUAL}"
hr

# ---------- Reference index ----------
[[ -f "${REF_FA}.fai" ]] || { step_start "Index reference"; samtools faidx "${REF_FA}"; step_done >/dev/null; }

# ---------- BED normalization ----------
BED_NORM="${OUT_DIR}/targets.cleaned.bed"
step_start "Normalize BED -> ${BED_NORM}"
awk 'BEGIN{OFS="\t"} /^track|^browser|^#/ {next} NF>=3 {print $1,$2,$3,$4,$5,$6}' "${BED}" | sed 's/\r$//' > "${BED_NORM}"
step_done >/dev/null

# ---------- Load manifest ----------
mapfile -t BAMS < "${MANIFEST}"
[[ ${#BAMS[@]} -gt 0 ]] || { echo "ERROR: empty manifest"; exit 1; }
TOTAL=${#BAMS[@]}
log "Samples to process: ${TOTAL}"
hr

# ---------- ETA state ----------
SAMPLE_VCFS=(); IDX=0; S_DONE=0; AVG_SAMPLE=0; AVG_TVC=0; N_TVC=0; AVG_NORM=0; N_NORM=0; AVG_FILTER=0; N_FILTER=0

# ---------- Per-sample loop ----------
for BAM in "${BAMS[@]}"; do
  IDX=$((IDX+1))
  [[ -f "${BAM}" ]] || { log "WARN: BAM missing, skipping: ${BAM}"; continue; }

  SAMPLE=$(basename "${BAM}" .bam)
  CURRENT_SAMPLE="${SAMPLE}"
  SDIR="${OUT_DIR}/${SAMPLE}"
  mkdir -p "${SDIR}/logs"

  [[ -f "${BAM}.bai" ]] || { step_start "Index BAM (${SAMPLE})"; samtools index "${BAM}"; step_done >/dev/null; }

  echo; log "===== Sample ${IDX}/${TOTAL}: ${SAMPLE} ====="
  [[ ${S_DONE} -gt 0 ]] && { REMAINING=$((TOTAL-(IDX-1))); est_secs=$((AVG_SAMPLE*REMAINING)); log "ETA before start: ~$(fmt_dur ${est_secs}) (ETA $(eta_ts ${est_secs}))"; } || log "Collecting timings on first sample…"

  progress_bar 0; printf '\n'; SAMPLE_T0=$(date +%s)

  # --- TVC call ---
  step_start "TVC (variant calling)"
  RAW_VCF=""
  if command -v variant_caller_pipeline.py >/dev/null; then
    TVC_CMD=(variant_caller_pipeline.py
      --input-bam       "${BAM}"
      --reference-fasta "${REF_FA}"
      --region-bed      "${BED_NORM}"
      --output-dir      "${SDIR}"
      --num-threads     "${THREADS}"
    )
    [[ -n "${TVC_JSON}"    ]] && TVC_CMD+=(--parameters-file "${TVC_JSON}")
    [[ -n "${HOTSPOT_VCF}" ]] && TVC_CMD+=(--input-vcf "${HOTSPOT_VCF}")
    [[ -n "${ERROR_MOTIF}" ]] && TVC_CMD+=(--error-motif "${ERROR_MOTIF}")
    echo "TVC WRAPPER CMD: ${TVC_CMD[*]}" > "${SDIR}/logs/tvc.cmd.log"
    "${TVC_CMD[@]}" > "${SDIR}/logs/tvc.stdout.log" 2> "${SDIR}/logs/tvc.stderr.log"
    RAW_VCF="${SDIR}/TSVC_variants.vcf"
    [[ -f "${RAW_VCF}" ]] || RAW_VCF="${SDIR}/variants.vcf"
  else
    TVC_CMD=(tvc
      --input-bam    "${BAM}"
      --reference    "${REF_FA}"
      --output-dir   "${SDIR}"
      --num-threads  "${THREADS}"
      --target-file  "${BED_NORM}"
      --output-vcf   "small_variants.vcf"
    )
    echo "TVC CMD: ${TVC_CMD[*]}" > "${SDIR}/logs/tvc.cmd.log"
    "${TVC_CMD[@]}" > "${SDIR}/logs/tvc.stdout.log" 2> "${SDIR}/logs/tvc.stderr.log"
    RAW_VCF="${SDIR}/TSVC_variants.vcf"
    [[ -f "${RAW_VCF}" ]] || RAW_VCF="${SDIR}/small_variants.vcf"
  fi

  [[ -f "${RAW_VCF}" ]] || { log "ERROR: VCF not found after TVC"; exit 1; }
  tvc_dt=$(step_done); progress_bar 33; printf '\n'
  AVG_TVC=$(running_avg "${AVG_TVC}" "${N_TVC}" "${tvc_dt}"); N_TVC=$((N_TVC+1))

  # --- VALIDATE/REPAIR VCF HEADER (reheader + optional rename) ---
  step_start "Validate/repair VCF header (reheader + optional rename)"
  # 1) Reinyectar contigs si faltan
  if ! bcftools view -h "${RAW_VCF}" | grep -q "^##contig"; then
    log "VCF header lacks ##contig. Reinserting contigs using ${REF_FA}.fai..."
    bcftools reheader --fai "${REF_FA}.fai" "${RAW_VCF}" -o "${SDIR}/variants.hdr.vcf"
    bgzip -f "${SDIR}/variants.hdr.vcf"
    tabix -p vcf "${SDIR}/variants.hdr.vcf.gz"
    RAW_VCF="${SDIR}/variants.hdr.vcf.gz"
  fi

  # 2) Detectar nomenclatura (chr vs sin chr)
  HAS_CHR_FA=$(grep -q "^>chr" "${REF_FA}" && echo 1 || echo 0)
  FIRST_VAR=$(bcftools view -H "${RAW_VCF}" | head -n 1 || true)
  HAS_CHR_VCF=$(echo "${FIRST_VAR}" | grep -q "^chr" && echo 1 || echo 0)

  # 3) Renombrado si hay desajuste
  if [[ -n "${FIRST_VAR}" ]]; then
    if [[ "${HAS_CHR_FA}" -eq 1 && "${HAS_CHR_VCF}" -eq 0 ]]; then
      log "FASTA usa 'chr*' y VCF no. Renombrando contigs en VCF -> 'chr*'..."
      awk '{g=$1; sub(/^chr/,"",g); print g "\t" $1}' "${REF_FA}.fai" > "${SDIR}/chr_map.txt"   # 1 -> chr1
      bcftools annotate --rename-chrs "${SDIR}/chr_map.txt" \
        "${RAW_VCF}" -Oz -o "${SDIR}/variants.renamed.vcf.gz"
      tabix -p vcf "${SDIR}/variants.renamed.vcf.gz"
      RAW_VCF="${SDIR}/variants.renamed.vcf.gz"
    elif [[ "${HAS_CHR_FA}" -eq 0 && "${HAS_CHR_VCF}" -eq 1 ]]; then
      log "FASTA no usa 'chr*' y VCF sí. Renombrando contigs en VCF -> sin 'chr'..."
      awk '{g=$1; print $1 "\t" gensub(/^chr/,"","",g)}' "${REF_FA}.fai" > "${SDIR}/chr_map.txt" # chr1 -> 1
      bcftools annotate --rename-chrs "${SDIR}/chr_map.txt" \
        "${RAW_VCF}" -Oz -o "${SDIR}/variants.renamed.vcf.gz"
      tabix -p vcf "${SDIR}/variants.renamed.vcf.gz"
      RAW_VCF="${SDIR}/variants.renamed.vcf.gz"
    else
      log "Nomenclatura de contigs coherente entre FASTA y VCF."
    fi
  else
    log "VCF sin variantes; saltamos rename-chrs. (Header ya reinyectado si faltaba)"
  fi

  # 4) Validación final del header
  if ! bcftools view -h "${RAW_VCF}" | grep -q "^##contig"; then
    log "ERROR: header aún sin ##contig tras reheader. Abortando para evitar corrupción."
    exit 1
  fi
  step_done

  # --- Normalize + compress ---
  step_start "Normalize + compress"
  bcftools norm -f "${REF_FA}" -m -both "${RAW_VCF}" -Oz -o "${SDIR}/variants.norm.vcf.gz"
  tabix -p vcf "${SDIR}/variants.norm.vcf.gz"
  norm_dt=$(step_done); progress_bar 66; printf '\n'
  AVG_NORM=$(running_avg "${AVG_NORM}" "${N_NORM}" "${norm_dt}"); N_NORM=$((N_NORM+1))

  # --- Unified filtering ---
  step_start "Filter (PASS; DP>=${MIN_DP}; QUAL>=${MIN_QUAL})"
  # Desambiguar DP: usar profundidad por muestra (FORMAT) y hacer el filtro robusto para 1..N muestras
  bcftools view -i "MAX(FMT/DP)>=${MIN_DP} && QUAL>=${MIN_QUAL} && FILTER='PASS'" \
    "${SDIR}/variants.norm.vcf.gz" -Oz -o "${SDIR}/variants.filtered.vcf.gz"
  tabix -p vcf "${SDIR}/variants.filtered.vcf.gz"
  filter_dt=$(step_done); progress_bar 100; printf '\n'
  AVG_FILTER=$(running_avg "${AVG_FILTER}" "${N_FILTER}" "${filter_dt}"); N_FILTER=$((N_FILTER+1))

  # --- Export manual review candidates (robusto y con AF derivada si falta) ---
  # Definimos helpers para detectar tags en el header
  has_info_tag() { bcftools view -h "${VCF_HDR}" | grep -q "^##INFO=<ID=$1,"; }
  has_fmt_tag()  { bcftools view -h "${VCF_HDR}" | grep -q "^##FORMAT=<ID=$1,"; }

  VCF_HDR="${SDIR}/variants.norm.vcf.gz"
  INFO_TAGS=(HP Indel STB STBP HS SB)  # posibles etiquetas de sesgo/homopolímero según versiones

  # 1) Export por-muestra (una fila por muestra y variante)
  step_start "Export manual review (per-sample metrics)"
  AO_RO_FMT=0; has_fmt_tag AO && has_fmt_tag RO && AO_RO_FMT=1
  AF_FMT=0;    has_fmt_tag AF && AF_FMT=1

  if [[ "${AO_RO_FMT}" -eq 1 ]]; then
    # Exportar CHROM POS REF ALT QUAL FILTER SAMPLE DP AO RO por muestra
    bcftools query -f "[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%SAMPLE\t%DP\t%AO\t%RO\n]" \
      "${VCF_HDR}" > "${SDIR}/manual_review_samples.raw.tsv" || true
    # Derivar AF = AO/(AO+RO); si AO/RO son '.', dejamos AF='.'
    awk -F'\t' 'BEGIN{OFS="\t"}
      {
        chrom=$1; pos=$2; ref=$3; alt=$4; qual=$5; filt=$6; sample=$7; dp=$8; ao=$9; ro=$10;
        if (ao ~ /^[0-9]+$/ && ro ~ /^[0-9]+$/ && (ao+ro)>0) { af = ao / (ao + ro) } else { af = "." }
        print chrom, pos, ref, alt, qual, filt, sample, dp, af
      }' "${SDIR}/manual_review_samples.raw.tsv" > "${SDIR}/manual_review_samples.tsv" || true
    rm -f "${SDIR}/manual_review_samples.raw.tsv"
  else
    # Si no hay AO/RO, intentar con AF en FORMAT; si tampoco, dejar AF='.'
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

  # 2) Export por-sitio (INFO tags solo si existen)
  step_start "Export manual review (site-level INFO tags)"
  FMT_STR="%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER"
  # Añadir INFO/DP si existe (métrica a nivel sitio)
  if bcftools view -h "${VCF_HDR}" | grep -q "^##INFO=<ID=DP,"; then
    FMT_STR+="\t%INFO/DP"
  else
    FMT_STR+="\t."
  fi
  # Añadir INFO/AF si existe (suele venir de ciertos pipelines; TVC a veces no lo pone)
  if bcftools view -h "${VCF_HDR}" | grep -q "^##INFO=<ID=AF,"; then
    FMT_STR+="\t%INFO/AF"
  else
    FMT_STR+="\t."
  fi
  # Añadir etiquetas de revisión si existen
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

  # --- Summary & separador visual ---
  N_PASS=$(bcftools view -H "${SDIR}/variants.filtered.vcf.gz" | wc -l || echo 0)
  log "Summary ${SAMPLE}: ${N_PASS} passing variants"
  hr

  # Timing per sample
  SAMPLE_T1=$(date +%s); SAMPLE_DT=$((SAMPLE_T1-SAMPLE_T0))
  AVG_SAMPLE=$(running_avg "${AVG_SAMPLE}" "${S_DONE}" "${SAMPLE_DT}"); S_DONE=$((S_DONE+1))

  SAMPLE_VCFS+=("${SDIR}/variants.filtered.vcf.gz")
done

# ---------- Cohort merge & exports ----------
if [[ ${#SAMPLE_VCFS[@]} -gt 1 ]]; then
  COHORT_DIR="${OUT_DIR}/cohort"; mkdir -p "${COHORT_DIR}"

  step_start "Merge cohort"
  bcftools merge "${SAMPLE_VCFS[@]}" -Oz -o "${COHORT_DIR}/cohort.filtered.vcf.gz" \
    > "${COHORT_DIR}/merge.log" 2>&1 || log "WARN: merge failed"
  [[ -f "${COHORT_DIR}/cohort.filtered.vcf.gz" ]] && tabix -p vcf "${COHORT_DIR}/cohort.filtered.vcf.gz"
  step_done

  bcftools query -f'%CHROM\t%POS\t%ID[\t%GT]\n' "${COHORT_DIR}/cohort.filtered.vcf.gz" \
    > "${COHORT_DIR}/cohort_genotypes.tsv" || true

  bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL[\t%DP]\n' "${COHORT_DIR}/cohort.filtered.vcf.gz" \
    > "${COHORT_DIR}/site_metrics.tsv" || true
fi

log "Pipeline completed. Outputs in: ${OUT_DIR}"

