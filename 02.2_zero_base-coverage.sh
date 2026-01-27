#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------
# Zero‑coverage analysis per amplicon, per sample, per cohort (FAST)
# ---------------------------------------------------------------
# Usage:
#   02.2_zero_base_coverage_fast.sh <ROOT_DIR> <OUT_DIR>
#
# ROOT_DIR must contain:
#   dna_qc/
#      ├── pass_bams/
#      ├── fail_bams/
#      └── targets.clean.bed
# ---------------------------------------------------------------

ROOT="$1"
OUT="$2"

BED="${ROOT}/dna_qc/targets.clean.bed"
PASS="${ROOT}/dna_qc/pass_bams"
FAIL="${ROOT}/dna_qc/fail_bams"

mkdir -p "$OUT"
TMP="$OUT/tmp"
mkdir -p "$TMP"

# Collect samples
mapfile -t BAMS < <(find "$PASS" "$FAIL" -type f -name "*.bam" | sort)

if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "ERROR: no BAMs found in $ROOT"
    exit 1
fi

echo "[INFO] Found ${#BAMS[@]} samples"
echo "[INFO] Using BED: $BED"

# Read BED
mapfile -t BED_ROWS < <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4}' "$BED")

# Build output header
HEADER="chr\tstart\tend\tannotation"
for BAM in "${BAMS[@]}"; do
    S=$(basename "$BAM" .bam)
    HEADER="${HEADER}\t${S}"
done

COHORT_NAME=$(basename "$(dirname "$ROOT")")_$(basename "$ROOT")
OUTFILE="${OUT}/zero_cov_${COHORT_NAME}.tsv"
echo -e "$HEADER" > "$OUTFILE"

echo "[INFO] Output → $OUTFILE"
echo "[INFO] Temporary mosdepth directory → $TMP"

# ---------------------------------------------------------------
# STEP 1 — Compute per-base mosdepth for each sample ONCE
# ---------------------------------------------------------------

declare -A PERBASE

echo "[INFO] Running per-base mosdepth (one per sample)..."

for BAM in "${BAMS[@]}"; do
    S=$(basename "$BAM" .bam)
    PREFIX="${TMP}/${S}"

    if [[ ! -f "${PREFIX}.per-base.bed.gz" ]]; then
        mosdepth --threads 4 \
            --fast-mode \
            "${PREFIX}" \
            "$BAM" > /dev/null 2>&1
    fi

    PERBASE["$S"]="${PREFIX}.per-base.bed.gz"
done

echo "[INFO] All per-base files ready."

# ---------------------------------------------------------------
# STEP 2 — Iterate over amplicons and count zero bases quickly
# ---------------------------------------------------------------

echo "[INFO] Computing zero-base counts..."

for ROW in "${BED_ROWS[@]}"; do
    IFS=$'\t' read -r chr start end annot <<< "$ROW"
    line="${chr}\t${start}\t${end}\t${annot}"

    for BAM in "${BAMS[@]}"; do
        S=$(basename "$BAM" .bam)
        PB="${PERBASE[$S]}"

        # Count depth==0 bases for this amplicon
        z=$(zcat "$PB" | awk -v c="$chr" -v s="$start" -v e="$end" '
            $1==c && $2>=s && $3<=e && $4==0 {count++} 
            END {print count+0}')
        
        line="${line}\t${z}"
    done

    echo -e "$line" >> "$OUTFILE"
done

echo "[SUCCESS] FAST zero-base analysis completed."
echo "[NOTE] Per-base files stored in $TMP."