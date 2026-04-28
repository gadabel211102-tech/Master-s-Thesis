#!/usr/bin/env bash
#
# Thesis launcher for the numbered GSDMB workflow stages.
# The file stays intentionally explicit so supervisors and external readers can
# see the execution order, runtime split, and optional branches without needing
# to reverse-engineer helper wrappers.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_PIPELINE_OVERRIDES="${ROOT_DIR}/pipeline_config.local.sh"
if [[ -f "${LOCAL_PIPELINE_OVERRIDES}" ]]; then
  # shellcheck source=/dev/null
  source "${LOCAL_PIPELINE_OVERRIDES}"
fi

# Runtime discovery prefers an explicitly chosen interpreter, then an active
# environment, then the repo-local thesis environment.
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
PANEL_1000G="${PANEL_1000G:-}"
ANALYSIS_ROOT="${ANALYSIS_ROOT:-${ROOT_DIR}/analysis_results}"
MIN_DP="${MIN_DP:-100}"
MIN_QUAL="${MIN_QUAL:-20}"
FEMALE_ONLY_SAMPLE_FILE="${FEMALE_ONLY_SAMPLE_FILE:-${ROOT_DIR}/analysis_results/15_haplotype_phasing/female_by_design_study_samples.tsv}"
COLLAB_CORE_XLSX="${COLLAB_CORE_XLSX:-}"
CROSSVAL_1000G_POPULATION="${CROSSVAL_1000G_POPULATION:-EUR}"


DRY_RUN=0
INCLUDE_1000G=1
PROFILE="${PROFILE:-core}"
FROM_STEP="01"
TO_STEP="26"

# Step groups mirror the workflow classification in README.md and
# script_classification_table.tsv.
STEP_ORDER=(01 02 02b 02c 03 04 04bV 04bP 05 05b 06 07 07b 08 09 09b 10 11 12 13 14 15 15R 16 17 18 19 19b 20 20b 21 22 23 24 25 26)
CORE_STEPS=(01 02 02b 05 05b 06 07 11 12 15 15R 16 17 18 19b 20 20b 23)
EXPLORATORY_VALIDATION_STEPS=(09b 13 19 22 24 25 26)
UTILITY_REPORTING_STEPS=(02c 03 04 04bV 04bP 07b 08 09 10 14 21)
PYTHON_STEPS=(02c 03 04 04bV 04bP 05b 07 07b 08 09 09b 10 11 12 13 14 16 17 18 19 19b 20 20b 21 22 23 24 25 26)
BAM_STEPS=(01 02 02b 05 15)

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Run the GSDMB thesis workflow in the documented order, with environment
switches handled automatically. The default profile now runs the reproducible
core pipeline only.

Options:
  --profile NAME      Workflow profile: core (default), exploratory, utility, full
  --from STEP         First step to run (default: ${FROM_STEP})
  --to STEP           Last step to run (default: ${TO_STEP})
  --dry-run           Print commands without executing them
  -h, --help          Show this help

Supported step labels:
  ${STEP_ORDER[*]}

Profile contents:
  core                01 02 02b 05 05b 06 07 11 12 15 15R 16 17 18 19b 20 20b 23
  exploratory         09b 13 19 22 24 25 26
  utility             02c 03 04 04bV 04bP 07b 08 09 10 14 21
  full                all classified steps in step order

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
  FEMALE_ONLY_SAMPLE_FILE  Optional explicit female-by-design study-sample list for the stage-15 female rerun (default: ${FEMALE_ONLY_SAMPLE_FILE})
  COLLAB_CORE_XLSX    Optional collaborator workbook for stage 21 core haplotype comparison (default: ${COLLAB_CORE_XLSX})
  CROSSVAL_1000G_POPULATION  Primary 1000 Genomes female reference population for stage 24 (default: ${CROSSVAL_1000G_POPULATION})


Examples:
  ./run_full_pipeline.sh
  ./run_full_pipeline.sh --profile exploratory
  ./run_full_pipeline.sh --profile utility --from 04 --to 14
  ./run_full_pipeline.sh --profile full
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

step_in_set() {
  local target="$1"
  shift
  local step
  for step in "$@"; do
    if [[ "${step}" == "${target}" ]]; then
      return 0
    fi
  done
  return 1
}

normalise_profile() {
  case "$1" in
    core) echo "core" ;;
    exploratory|exploratory-validation|validation) echo "exploratory" ;;
    utility|utility-reporting|reporting) echo "utility" ;;
    full|all) echo "full" ;;
    *) return 1 ;;
  esac
}

step_selected_by_profile() {
  local target="$1"
  case "${PROFILE}" in
    core) step_in_set "${target}" "${CORE_STEPS[@]}" ;;
    exploratory) step_in_set "${target}" "${EXPLORATORY_VALIDATION_STEPS[@]}" ;;
    utility) step_in_set "${target}" "${UTILITY_REPORTING_STEPS[@]}" ;;
    full) step_in_set "${target}" "${STEP_ORDER[@]}" ;;
    *) return 1 ;;
  esac
}

should_run_step() {
  local target="$1"
  local idx
  idx="$(step_index "${target}")" || return 1
  [[ "${idx}" -ge "${FROM_INDEX}" && "${idx}" -le "${TO_INDEX}" ]] || return 1
  step_selected_by_profile "${target}"
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
  # Only validate the toolchains required by the selected profile / step range.
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

  if should_run_any 19 24; then
    [[ -n "${VCF_1000G}" ]] || die "VCF_1000G is not set. Provide it via env or pipeline_config.local.sh"
    [[ -n "${PANEL_1000G}" ]] || die "PANEL_1000G is not set. Provide it via env or pipeline_config.local.sh"
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

  # The DNA QC stages iterate over the fixed cohort/tissue matrix that underpins
  # the whole thesis dataset.
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

  # The DNA QC stages iterate over the fixed cohort/tissue matrix that underpins
  # the whole thesis dataset.
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
      --profile) PROFILE="$2"; shift 2 ;;
      --from) FROM_STEP="$2"; shift 2 ;;
      --to) TO_STEP="$2"; shift 2 ;;
      --dry-run) DRY_RUN=1; shift ;;
      -h|--help) usage; exit 0 ;;
      *) die "Unknown option: $1" ;;
    esac
  done
}

parse_args "$@"

PROFILE="$(normalise_profile "${PROFILE}")" || die "Unknown --profile value: ${PROFILE}"
FROM_INDEX="$(step_index "${FROM_STEP}")" || die "Unknown --from step: ${FROM_STEP}"
TO_INDEX="$(step_index "${TO_STEP}")" || die "Unknown --to step: ${TO_STEP}"
[[ "${FROM_INDEX}" -le "${TO_INDEX}" ]] || die "--from must be earlier than or equal to --to"

ensure_prereqs

log "Pipeline root: ${ROOT_DIR}"
log "Profile: ${PROFILE}"
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
if [[ -n "${FEMALE_ONLY_SAMPLE_FILE}" ]]; then
  log "Female-only study sample file: ${FEMALE_ONLY_SAMPLE_FILE}"
fi
if [[ -n "${COLLAB_CORE_XLSX}" ]]; then
  log "Stage-21 collaborator workbook: ${COLLAB_CORE_XLSX}"
fi
log "Cross-validation 1000G population: ${CROSSVAL_1000G_POPULATION}"

# The execution block below is intentionally linear so the numbered thesis steps
# remain easy to cross-reference in the manuscript and repository docs.
if should_run_any 01 02 02b; then run_qc_loops; fi
if should_run_step 02c; then announce_step "02c" "Building zero-coverage Excel report"; run_python_script 02c_zerocovinfo.py --zero "${ROOT_DIR}/breast/normal/zero_coverage/zero_cov_normal.tsv" "${ROOT_DIR}/breast/tumour/zero_coverage/zero_cov_tumour.tsv" "${ROOT_DIR}/endometrium/normal/zero_coverage/zero_cov_normal.tsv" "${ROOT_DIR}/endometrium/tumour/zero_coverage/zero_cov_tumour.tsv" --qc "${ROOT_DIR}/breast/normal/dna_qc/qc_summary.tsv" "${ROOT_DIR}/breast/tumour/dna_qc/qc_summary.tsv" "${ROOT_DIR}/endometrium/normal/dna_qc/qc_summary.tsv" "${ROOT_DIR}/endometrium/tumour/dna_qc/qc_summary.tsv" --output "${ANALYSIS_ROOT}/02c_zero_coverage_report/GSDMB_Zero_Coverage_Report.xlsx"; fi
if should_run_step 03; then announce_step "03" "Generating QC visualisations"; run_python_script 03_qc_visualisation.py; fi
if should_run_step 04; then announce_step "04" "Running technical audit"; run_python_script 04_technical_audit.py both --bed "${ROOT_DIR}/dna_bed/IAD255368_167_Submitted.bed" --fasta "${ROOT_DIR}/ref_alt/hg38_alt.fa" --output "${ANALYSIS_ROOT}/04_technical_audit"; fi
if should_run_step 04bV; then announce_step "04bV" "Building flagged amplicon failure visualisations"; run_python_script 04b_amplicon_failure_visualisation.py; fi
if should_run_step 04bP; then announce_step "04bP" "Building the flagged amplicon overview panel"; run_python_script 04b_amplicon_failure_panel.py; fi
if should_run_step 05; then run_variant_calling_loop; fi
if should_run_step 05b; then announce_step "05b" "Force genotyping all cohort-observed loci across every sample"; run_python_script 05b_force_genotype_union_sites.py --project-root "${ROOT_DIR}" --reference "${REF_FA}" --output-dir "${ANALYSIS_ROOT}/05b_forced_genotypes" --threads "${THREADS}" --callable-dp "${MIN_DP}"; fi
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
if should_run_step 15R; then
  announce_step "15R" "Running haplotype statistics in R"
  run_r_script 15_haplotype_stats.R
  if [[ -n "${FEMALE_ONLY_SAMPLE_FILE}" && -f "${FEMALE_ONLY_SAMPLE_FILE}" ]]; then
    announce_step "15R" "Running female-by-design study-sample haplotype rerun"
    run_in_root env TFM_ANALYSIS_LABEL="FemaleOnly_StudySamples" TFM_SAMPLE_INCLUDE_FILE="${FEMALE_ONLY_SAMPLE_FILE}" TFM_ANALYSIS_CONTEXT="Study-sample rerun restricted to the female-by-design cohort set; the 1000 Genomes reference comparison remains mixed-sex." "${R_BIN}" "${ROOT_DIR}/15_haplotype_stats.R"
  else
    log "Skipping female-by-design stage-15 rerun: FEMALE_ONLY_SAMPLE_FILE is not available."
  fi
fi
if should_run_step 16; then announce_step "16" "Harmonising clinical workbooks"; run_python_script 16_excel_harmonisation.py; fi
if should_run_step 17; then announce_step "17" "Running SNP-clinical association analysis"; run_python_script 17_snp_association.py; fi
if should_run_step 18; then announce_step "18" "Running haplotype-clinical association analysis"; run_python_script 18_haplotype_association.py; fi
if should_run_step 19; then
  announce_step "19" "Comparing study haplotypes against 1000 Genomes"
  run_python_script 19_1000g_haplotype_comparison.py --vcf "${VCF_1000G}" --panel "${PANEL_1000G}"
fi
if should_run_step 19b; then announce_step "19b" "Running standalone objective-2 RNA QC"; run_python_script 19b_rna_qc.py; fi
if should_run_step 20; then announce_step "20" "Running workbook-first objective-2 isoform-expression association analysis"; run_python_script 20_objective2_isoform_expression_association.py; fi
if should_run_step 20b; then announce_step "20b" "Quantifying panel-gene RNA expression from strict and exploratory RNA BAMs"; run_python_script 20b_panel_gene_expression.py; fi
if should_run_step 21; then
  if [[ -n "${COLLAB_CORE_XLSX}" && -f "${COLLAB_CORE_XLSX}" ]]; then
    announce_step "21" "Building core haplotype comparison summaries"
    run_python_script 21_core_haplotype_comparison.py --collab-xlsx "${COLLAB_CORE_XLSX}"
  else
    log "Skipping optional stage 21: collaborator workbook not found at ${COLLAB_CORE_XLSX}"
  fi
fi
if should_run_step 22; then announce_step "22" "Building haplotype-first interpretation outputs"; run_python_script 22_haplotype_first_interpretation.py; fi
if should_run_step 23; then announce_step "23" "Integrating Excel RNA and BAM-derived RNA with significant SNP and haplotype signals"; run_python_script 23_rna_integration.py; fi

if should_run_step 24; then
  announce_step "24" "Running external SNP and haplotype cross-validation for breast and endometrial tumours against female-only 1000 Genomes"
  stage24_args=(--vcf "${VCF_1000G}" --panel "${PANEL_1000G}" --reference-population "${CROSSVAL_1000G_POPULATION}")
  run_python_script 24_external_cross_validation.py "${stage24_args[@]}"
fi
if should_run_step 25; then announce_step "25" "Running the focused rs11078928 / rs869402 haplotype analysis"; run_python_script 25_rs11078928_rs869402_haplotype_focus.py; fi
if should_run_step 26; then announce_step "26" "Building the final SNP functional interpretation and biological-context layer"; run_python_script 26_snp_functional_interpretation.py; fi

if should_run_any 02 02b 19b; then
  announce_step "QC-SUMMARY" "Refreshing the combined DNA/RNA sample-count summary"
  run_python_script pipeline_sample_summary.py --project-root "${ROOT_DIR}" --analysis-root "${ANALYSIS_ROOT}"
fi

log "Pipeline run completed."












