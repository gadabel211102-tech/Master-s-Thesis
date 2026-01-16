
#!/usr/bin/env bash
set -euo pipefail

# Usage: ./01_check_and_index.sh -i /path/to/bams -t 4
INPUT_DIR=""
THREADS="${THREADS:-4}"

usage() {
  echo "Usage: $0 -i INPUT_DIR [-t THREADS]" >&2
}

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

mapfile -t BAMS < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*.bam" | sort)
[[ ${#BAMS[@]} -gt 0 ]] || { echo "No BAMs in ${INPUT_DIR}"; exit 1; }

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
