
#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Script 01: BAM Integrity Check and Index Generation
# =============================================================================
# Purpose:
#   Verify that each BAM file is readable and create a BAM index (.bai) when it
#   is missing. This is the first gate in the pipeline, ensuring that downstream
#   QC and variant calling operate on structurally valid alignment files.
#
# Inputs:
#   - A directory containing BAM files.
#   - The number of threads to use for samtools indexing.
#
# Outputs:
#   - Existing BAM files are left in place.
#   - Missing .bai index files are created alongside their BAM files.
#
# Analytical note:
#   This script does not assess biological quality. It only confirms technical
#   readability and indexing status so the subsequent QC stage starts from a
#   consistent set of alignment files.
# =============================================================================

# Usage: ./01_check_and_index.sh -i /path/to/bams -t 4
INPUT_DIR=""
THREADS="${THREADS:-4}"

usage() {
  echo "Usage: $0 -i INPUT_DIR [-t THREADS]" >&2
}

# Parse command-line arguments before any file-system work is attempted.
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_DIR="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 2 ;;
  esac
done

[[ -z "${INPUT_DIR}" ]] && { usage; exit 1; }
command -v samtools >/dev/null || { echo "samtools not found"; exit 1; }

# Enumerate the BAM inputs once so the script can fail early if the folder
# is empty instead of doing partial work.
mapfile -t BAMS < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*.bam" | sort)
[[ ${#BAMS[@]} -gt 0 ]] || { echo "No BAMs in ${INPUT_DIR}"; exit 1; }

# Each BAM is validated with samtools quickcheck before indexing. Files
# that fail integrity checks are reported and skipped rather than stopping the
# entire batch.
for BAM in "${BAMS[@]}"; do
  echo ">> Checking ${BAM}"
  if ! samtools quickcheck -v "${BAM}"; then
    echo "!! ${BAM} failed samtools quickcheck — skip"; continue
  fi

  if [[ -f "${BAM}.bai" ]]; then
    echo "== Index exists: ${BAM}.bai"
  else
    echo ">> Indexing ${BAM}"
    samtools index -@ "${THREADS}" "${BAM}"
  fi
done

echo "Done."
