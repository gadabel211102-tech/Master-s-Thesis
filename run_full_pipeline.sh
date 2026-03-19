#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

THREADS="${THREADS:-8}"
BAM_ENV_NAME="${BAM_ENV_NAME:-bam-steps}"
VEP_ENV_NAME="${VEP_ENV_NAME:-vep_env}"
REPO_VEP_ENV_BIN="${ROOT_DIR}/miniconda3/envs/${VEP_ENV_NAME}/bin"
TFM_ENV_PYTHON_CANDIDATE_1="${HOME}/tfm_env/bin/python"
TFM_ENV_PYTHON_CANDIDATE_2="${ROOT_DIR}/tfm_env/bin/python"
if [[ -x "${TFM_ENV_PYTHON_CANDIDATE_1}" ]]; then
  TFM_ENV_PYTHON="${TFM_ENV_PYTHON_CANDIDATE_1}"
else
  TFM_ENV_PYTHON="${TFM_ENV_PYTHON_CANDIDATE_2}"
fi
PYTHON_SOURCE="unset"
if [[ -n "${PYTHON_BIN:-}" ]]; then
  PYTHON_BIN="${PYTHON_BIN}"
  PYTHON_SOURCE="explicit"
elif [[ -n "${VIRTUAL_ENV:-}" ]] && [[ -x "${VIRTUAL_ENV}/bin/python" ]]; then
  PYTHON_BIN="${VIRTUAL_ENV}/bin/python"
  PYTHON_SOURCE="active_venv"
elif [[ -x "${TFM_ENV_PYTHON}" ]]; then
  PYTHON_BIN="${TFM_ENV_PYTHON}"
  PYTHON_SOURCE="tfm_env"
elif command -v python >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python)"
  PYTHON_SOURCE="path_python"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python3)"
  PYTHON_SOURCE="path_python3"
else
  PYTHON_BIN="${TFM_ENV_PYTHON}"
  PYTHON_SOURCE="missing"
fi
R_BIN="${R_BIN:-Rscript}"
REF_FA="${REF_FA:-${ROOT_DIR}/ref_alt/hg38_canonical.fa}"
BED_FILE="${BED_FILE:-${ROOT_DIR}/dna_bed/IAD255368_167_Submitted.bed}"
VCF_1000G="${VCF_1000G:-${ROOT_DIR}/ref/ref_panel_chr17_gsdmb_GRCh38.vcf.gz}"
PANEL_1000G="${PANEL_1000G:-/home/gadeaalonsoj/1000g_grch38/integrated_call_samples_v3.20130502.ALL.panel}"
ANALYSIS_ROOT="${ANALYSIS_ROOT:-${ROOT_DIR}/analysis_results}"
MIN_DP="${MIN_DP:-100}"
MIN_QUAL="${MIN_QUAL:-20}"

DRY_RUN=0
INCLUDE_1000G=1
FROM_STEP="01"
TO_STEP="19"
TO_STEP_EXPLICIT=0

STEP_ORDER=(01 02 02b 02c 03 04 05 06 07 07b 08 09 09b 10 11 12 13 14 15 15R 16 17 18 19)
PYTHON_STEPS=(02c 03 04 07 07b 08 09 09b 10 11 12 13 14 16 17 18 19)
BAM_STEPS=(01 02 02b 05 15)

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Run the GSDMB pipeline in the documented order, with the expected environment
switches handled automatically.

Options:
  --from STEP         First step to run (default: ${FROM_STEP})
  --to STEP           Last step to run (default: ${TO_STEP})
  --dry-run           Print commands without executing them
  -h, --help          Show this help

Supported step labels:
  ${STEP_ORDER[*]}

Environment overrides:
  THREADS             Threads for shell stages (default: ${THREADS})
  MIN_DP              Variant-calling minimum depth (default: ${MIN_DP})
  MIN_QUAL            Variant-calling minimum QUAL (default: ${MIN_QUAL})
  BAM_ENV_NAME        micromamba env for BAM/QC/calling/phasing (default: ${BAM_ENV_NAME})
  VEP_ENV_NAME        conda env for VEP annotation (default: ${VEP_ENV_NAME})
  PYTHON_BIN          Python interpreter for analysis scripts (default: ${PYTHON_BIN})
  R_BIN               Rscript executable (default: ${R_BIN})
  REF_FA              Reference FASTA (default: ${REF_FA})
  BED_FILE            Target BED (default: ${BED_FILE})
  VCF_1000G           1000 Genomes phased VCF/BCF (default: ${VCF_1000G})
  PANEL_1000G         1000 Genomes population panel (default: ${PANEL_1000G})
  ANALYSIS_ROOT       Consolidated analysis output root (default: ${ANALYSIS_ROOT})

Examples:
  ./run_full_pipeline.sh
  PYTHON_BIN=/path/to/venv/bin/python ./run_full_pipeline.sh
  PYTHON_BIN=/path/to/venv/bin/python ./run_full_pipeline.sh --from 07
EOF
}

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

announce_step() {
  local label="$1"
  shift
  log "STEP ${label} - $*"
}

quote_args() {
  local out=""
  local arg
  for arg in "$@"; do
    out+=" $(printf '%q' "$arg")"
  done
  printf '%s\n' "${out# }"
}

run_in_root() {
  log "(cd ${ROOT_DIR} && $(quote_args "$@"))"
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    (
      cd "${ROOT_DIR}"
      "$@"
    )
  fi
}

run_bam_script() {
  local script_name="$1"
  shift
  run_in_root micromamba run -n "${BAM_ENV_NAME}" bash "${ROOT_DIR}/${script_name}" "$@"
}

run_python_script() {
  local script_name="$1"
  shift
  run_in_root "${PYTHON_BIN}" "${ROOT_DIR}/${script_name}" "$@"
}

run_r_script() {
  local script_name="$1"
  shift
  run_in_root "${R_BIN}" "${ROOT_DIR}/${script_name}" "$@"
}

run_vep_script() {
  if [[ -x "${REPO_VEP_ENV_BIN}/vep" ]]; then
    run_in_root env PATH="${REPO_VEP_ENV_BIN}:${PATH}" bash "${ROOT_DIR}/06_annotation.sh"
  elif command -v conda >/dev/null 2>&1; then
    run_in_root conda run -n "${VEP_ENV_NAME}" bash "${ROOT_DIR}/06_annotation.sh"
  elif command -v micromamba >/dev/null 2>&1 && micromamba env list | awk 'NR > 1 {print $1}' | grep -Fx "${VEP_ENV_NAME}" >/dev/null; then
    run_in_root micromamba run -n "${VEP_ENV_NAME}" bash "${ROOT_DIR}/06_annotation.sh"
  elif [[ -x "${ROOT_DIR}/ensembl-vep/vep" ]]; then
    run_in_root env PATH="${ROOT_DIR}/ensembl-vep:${PATH}" bash "${ROOT_DIR}/06_annotation.sh"
  else
    die "No usable VEP runtime found. Expected repo-local ${REPO_VEP_ENV_BIN}/vep, conda env '${VEP_ENV_NAME}', micromamba env '${VEP_ENV_NAME}', or local ${ROOT_DIR}/ensembl-vep/vep"
  fi
}

step_index() {
  local target="$1"
  local i
  for i in "${!STEP_ORDER[@]}"; do
    if [[ "${STEP_ORDER[$i]}" == "${target}" ]]; then
      echo "${i}"
      return 0
    fi
  done
  return 1
}

should_run_step() {
  local target="$1"
  local idx
  idx="$(step_index "${target}")" || return 1
  [[ "${idx}" -ge "${FROM_INDEX}" && "${idx}" -le "${TO_INDEX}" ]]
}

should_run_any() {
  local step
  for step in "$@"; do
    if should_run_step "${step}"; then
      return 0
    fi
  done
  return 1
}

ensure_prereqs() {
  if should_run_any "${BAM_STEPS[@]}"; then
    command -v micromamba >/dev/null 2>&1 || die "micromamba is required for BAM-stage execution."
  fi

  if should_run_step 06; then
    if [[ -x "${REPO_VEP_ENV_BIN}/vep" ]]; then
      :
    elif command -v conda >/dev/null 2>&1; then
      :
    elif command -v micromamba >/dev/null 2>&1 && micromamba env list | awk 'NR > 1 {print $1}' | grep -Fx "${VEP_ENV_NAME}" >/dev/null; then
      :
    elif [[ -x "${ROOT_DIR}/ensembl-vep/vep" ]]; then
      :
    else
      die "No usable VEP runtime found. Expected repo-local ${REPO_VEP_ENV_BIN}/vep, conda env '${VEP_ENV_NAME}', micromamba env '${VEP_ENV_NAME}', or local ${ROOT_DIR}/ensembl-vep/vep"
    fi
  fi

  if should_run_any "${PYTHON_STEPS[@]}"; then
    [[ -x "${PYTHON_BIN}" ]] || die "Python interpreter not executable: ${PYTHON_BIN}"
    if [[ "${PYTHON_SOURCE}" == "path_python" || "${PYTHON_SOURCE}" == "path_python3" ]]; then
      [[ -x "${TFM_ENV_PYTHON}" ]] || die "tfm_env is not usable (${TFM_ENV_PYTHON} is missing or broken). Activate the real analysis venv first, or run with PYTHON_BIN=/path/to/your/python."
    fi
  fi

  if should_run_step 15R; then
    command -v "${R_BIN}" >/dev/null 2>&1 || die "R executable not found on PATH: ${R_BIN}"
  fi

  if should_run_any 02 05; then
    [[ -f "${REF_FA}" ]] || die "Reference FASTA not found: ${REF_FA}"
    [[ -f "${BED_FILE}" ]] || die "Target BED not found: ${BED_FILE}"
  fi

  if should_run_any 02b 05; then
    mkdir -p "${ROOT_DIR}/manifests"
  fi

  if should_run_step 19; then
    [[ -f "${VCF_1000G}" ]] || die "1000 Genomes VCF/BCF not found: ${VCF_1000G}"
    [[ -f "${PANEL_1000G}" ]] || die "1000 Genomes panel not found: ${PANEL_1000G}"
  fi
}

cohort_input_dir() {
  local cohort="$1"
  local tissue="$2"
  echo "${ROOT_DIR}/${cohort}/${tissue}/dna"
}

cohort_qc_dir() {
  local cohort="$1"
  local tissue="$2"
  echo "${ROOT_DIR}/${cohort}/${tissue}/dna_qc"
}

cohort_zero_cov_dir() {
  local cohort="$1"
  local tissue="$2"
  echo "${ROOT_DIR}/${cohort}/${tissue}/zero_coverage"
}

cohort_manifest() {
  local cohort="$1"
  local tissue="$2"
  echo "${ROOT_DIR}/manifests/${cohort}-${tissue}-pass_manifest.txt"
}

cohort_calls_dir() {
  local cohort="$1"
  local tissue="$2"
  echo "${ROOT_DIR}/${cohort}/${tissue}/dna_calls"
}

run_qc_loops() {
  local cohort
  local tissue
  local input_dir
  local qc_dir
  local zero_cov_dir
  local manifest

  for cohort in breast endometrium; do
    for tissue in normal tumour; do
      input_dir="$(cohort_input_dir "${cohort}" "${tissue}")"
      qc_dir="$(cohort_qc_dir "${cohort}" "${tissue}")"
      zero_cov_dir="$(cohort_zero_cov_dir "${cohort}" "${tissue}")"
      manifest="$(cohort_manifest "${cohort}" "${tissue}")"

      if should_run_any 01 02; then
        [[ -d "${input_dir}" ]] || die "Input DNA directory not found: ${input_dir}"
      fi

      if should_run_step 01; then
        announce_step "01" "Checking BAM integrity and indexing for ${cohort}/${tissue}"
        run_bam_script 01_check_and_index.sh -i "${input_dir}" -t "${THREADS}"
      fi

      if should_run_step 02; then
        announce_step "02" "Running DNA QC for ${cohort}/${tissue}"
        run_bam_script 02_dna_qc.sh -i "${input_dir}" -o "${qc_dir}" -r "${REF_FA}" -b "${BED_FILE}" -t "${THREADS}"
      fi

      if should_run_step 02b; then
        [[ -d "${qc_dir}" ]] || die "QC directory not found for ${cohort}/${tissue}: ${qc_dir}"
        announce_step "02b" "Building PASS manifest and zero-coverage summaries for ${cohort}/${tissue}"
        run_bam_script 02b_more-qc.sh -q "${qc_dir}" -m "${manifest}" -z "${zero_cov_dir}" -t "${THREADS}"
      fi
    done
  done
}

run_variant_calling_loop() {
  local cohort
  local tissue
  local manifest
  local calls_dir

  for cohort in breast endometrium; do
    for tissue in normal tumour; do
      manifest="$(cohort_manifest "${cohort}" "${tissue}")"
      calls_dir="$(cohort_calls_dir "${cohort}" "${tissue}")"
      [[ -f "${manifest}" ]] || die "Manifest not found for ${cohort}/${tissue}: ${manifest}"
      announce_step "05" "Calling variants for ${cohort}/${tissue}"
      run_bam_script 05_variant_calling.sh -m "${manifest}" -r "${REF_FA}" -b "${BED_FILE}" -o "${calls_dir}" -t "${THREADS}" --min-dp "${MIN_DP}" --min-qual "${MIN_QUAL}"
    done
  done
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --from) FROM_STEP="$2"; shift 2 ;;
      --to) TO_STEP="$2"; TO_STEP_EXPLICIT=1; shift 2 ;;
      --dry-run) DRY_RUN=1; shift ;;
      -h|--help) usage; exit 0 ;;
      *) die "Unknown option: $1" ;;
    esac
  done
}

parse_args "$@"

FROM_INDEX="$(step_index "${FROM_STEP}")" || die "Unknown --from step: ${FROM_STEP}"
TO_INDEX="$(step_index "${TO_STEP}")" || die "Unknown --to step: ${TO_STEP}"
[[ "${FROM_INDEX}" -le "${TO_INDEX}" ]] || die "--from must be earlier than or equal to --to"

ensure_prereqs

log "Pipeline root: ${ROOT_DIR}"
log "Running steps: ${FROM_STEP} -> ${TO_STEP}"
log "Threads: ${THREADS}"
log "BAM env: ${BAM_ENV_NAME}"
log "VEP env: ${VEP_ENV_NAME}"
if should_run_any "${PYTHON_STEPS[@]}"; then log "Python: ${PYTHON_BIN}"; fi
if should_run_step 15R; then log "R: ${R_BIN}"; fi
if should_run_step 19; then
  log "1000G VCF: ${VCF_1000G}"
  log "1000G panel: ${PANEL_1000G}"
fi

if should_run_any 01 02 02b; then run_qc_loops; fi
if should_run_step 02c; then announce_step "02c" "Building zero-coverage Excel report"; run_python_script 02c_zerocovinfo.py --zero "${ROOT_DIR}/breast/normal/zero_coverage/zero_cov_normal.tsv" "${ROOT_DIR}/breast/tumour/zero_coverage/zero_cov_tumour.tsv" "${ROOT_DIR}/endometrium/normal/zero_coverage/zero_cov_normal.tsv" "${ROOT_DIR}/endometrium/tumour/zero_coverage/zero_cov_tumour.tsv" --qc "${ROOT_DIR}/breast/normal/dna_qc/qc_summary.tsv" "${ROOT_DIR}/breast/tumour/dna_qc/qc_summary.tsv" "${ROOT_DIR}/endometrium/normal/dna_qc/qc_summary.tsv" "${ROOT_DIR}/endometrium/tumour/dna_qc/qc_summary.tsv" --output "${ANALYSIS_ROOT}/02c_zero_coverage_report/GSDMB_Zero_Coverage_Report.xlsx"; fi
if should_run_step 03; then announce_step "03" "Generating QC visualisations"; run_python_script 03_qc_visualisation.py; fi
if should_run_step 04; then announce_step "04" "Running technical audit"; run_python_script 04_technical_audit.py both --bed "${ROOT_DIR}/dna_bed/IAD255368_167_Submitted.bed" --fasta "${ROOT_DIR}/ref_alt/hg38_alt.fa" --output "${ANALYSIS_ROOT}/04_technical_audit"; fi
if should_run_step 05; then run_variant_calling_loop; fi
if should_run_step 06; then announce_step "06" "Annotating variants with VEP"; run_vep_script; fi
if should_run_step 07; then announce_step "07" "Merging annotated VCFs into the canonical workbook"; run_python_script 07_merge_annotations.py; fi
if should_run_step 07b; then announce_step "07b" "Running variant-level QC checks"; run_python_script 07b_variant_qc.py; fi
if should_run_step 08; then announce_step "08" "Generating the global variant landscape"; run_python_script 08_mapping.py; fi
if should_run_step 09; then announce_step "09" "Generating the GSDMB-only landscape"; run_python_script 09_gsdmb_only.py; fi
if should_run_step 09b; then announce_step "09b" "Comparing normal-versus-tumour landscape patterns with QC and significance checks"; run_python_script 09b_landscape_comparison.py; fi
if should_run_step 10; then announce_step "10" "Summarising descriptive variant statistics"; run_python_script 10_variant_stats.py; fi
if should_run_step 11; then announce_step "11" "Identifying common SNPs and benchmarking frequencies"; run_python_script 11_SNPs.py; fi
if should_run_step 12; then announce_step "12" "Running SNP enrichment testing"; run_python_script 12_stats_enrichment.py; fi
if should_run_step 13; then announce_step "13" "Running permutation-based SNP testing"; run_python_script 13_permutation_testing.py; fi
if should_run_step 14; then announce_step "14" "Building the interactive dashboard"; run_python_script 14_interactive_dashboard.py; fi
if should_run_step 15; then announce_step "15" "Phasing haplotypes with the cohort VCF set"; run_bam_script 15_haplotypes.sh; fi
if should_run_step 15R; then announce_step "15R" "Running haplotype statistics in R"; run_r_script 15_haplotype_stats.R; fi
if should_run_step 16; then announce_step "16" "Harmonising clinical workbooks"; run_python_script 16_excel_harmonisation.py; fi
if should_run_step 17; then announce_step "17" "Running SNP-clinical association analysis"; run_python_script 17_snp_association.py; fi
if should_run_step 18; then announce_step "18" "Running haplotype-clinical association analysis"; run_python_script 18_haplotype_association.py; fi
if should_run_step 19; then
  announce_step "19" "Comparing study haplotypes against 1000 Genomes"
  run_python_script 19_1000g_haplotype_comparison.py --vcf "${VCF_1000G}" --panel "${PANEL_1000G}"
fi

log "Pipeline run completed."
