#!/usr/bin/env bash
# Validates BAM file integrity and generates indices prior to variant calling.
# Usage: ./01_check_and_index.sh -i INPUT_DIR [-t THREADS]
set -euo pipefail  # exit on error, treat unset variables as errors, fail on pipe errors

INPUT_DIR=""
THREADS="${THREADS:-4}"  # default to 4 threads if not set in the environment

usage() { echo "Usage: $0 -i INPUT_DIR [-t THREADS]" >&2; }

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }  # timestamped logging

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_DIR="$2"; shift 2 ;;   # path to directory containing BAM files
    -t) THREADS="$2";   shift 2 ;;   # optional: number of parallel threads for indexing
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

# Validate required inputs before proceeding
[[ -z "${INPUT_DIR}" ]]  && { usage; exit 1; }
[[ -d "${INPUT_DIR}" ]]  || { echo "ERROR: directory not found: ${INPUT_DIR}" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }

# Collect all BAM files in INPUT_DIR; non-recursive and sorted for reproducibility
mapfile -t BAMS < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*.bam" | sort)
[[ ${#BAMS[@]} -gt 0 ]] || { echo "ERROR: no BAM files found in ${INPUT_DIR}" >&2; exit 1; }

log "Found ${#BAMS[@]} BAM file(s) — threads: ${THREADS}"

PASS=0; SKIP=0  # counters for final summary

for BAM in "${BAMS[@]}"; do
  log "Checking: $(basename "${BAM}")"

  # samtools quickcheck validates the EOF marker and basic record structure;
  # corrupt or truncated BAMs are skipped here to prevent downstream failures
  if ! samtools quickcheck -v "${BAM}" 2>&1; then
    log "WARNING: $(basename "${BAM}") failed integrity check — skipping"
    (( SKIP++ )) || true; continue
  fi

  # A BAI index file is required by downstream tools (e.g. GATK, bcftools) to
  # enable random access to genomic regions without reading the entire file;
  # only generated if absent, as re-indexing an unchanged BAM is redundant
  if [[ -f "${BAM}.bai" ]]; then
    log "Index exists: $(basename "${BAM}").bai — skipping"
  else
    samtools index -@ "${THREADS}" "${BAM}"  # -@ sets the number of threads
    log "Indexed: $(basename "${BAM}").bai"
  fi

  (( PASS++ )) || true
done

log "Done — passed: ${PASS} | skipped: ${SKIP}"
