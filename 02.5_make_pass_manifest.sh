
#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./make_manifest.sh /path/to/cohort/dna_qc /path/to/output_manifest.txt
#
# Example:
#   ./make_manifest.sh /home/gadeaalonsoj/tfm/breast/tumour/dna_qc \
#       /home/gadeaalonsoj/tfm/manifests/pass_breast_tumour.txt

QC_DIR="$1"
OUT="$2"

QC_SUM="${QC_DIR}/qc_summary.tsv"
PASS_DIR="${QC_DIR}/pass_bams"

[[ -f "$QC_SUM" ]]      || { echo "ERROR: qc_summary.tsv not found: $QC_SUM"; exit 1; }
[[ -d "$PASS_DIR" ]]    || { echo "ERROR: pass_bams/ not found: $PASS_DIR"; exit 1; }

echo "[INFO] Reading: $QC_SUM"
echo "[INFO] Looking for PASS samples..."

# Extract PASS sample names
SAMPLES=$(awk -F'\t' '
  NR==1 {
    for(i=1;i<=NF;i++) h[$i]=i
    if(!("sample" in h) || !("status" in h)) {
      print "ERROR: Missing sample/status columns in qc_summary.tsv" >"/dev/stderr"
      exit 2
    }
    next
  }
  $h["status"]=="PASS" { print $h["sample"] }
' "$QC_SUM")

# Build manifest using actual files in pass_bams/
echo "[INFO] Writing manifest → $OUT"
: > "$OUT"

for sample in $SAMPLES; do
  bam="${PASS_DIR}/${sample}.bam"
  if [[ -f "$bam" ]]; then
    echo "$bam" >> "$OUT"
  else
    echo "[WARN] Missing BAM for PASS sample: $sample" >&2
  fi
done

echo "[INFO] Manifest done."
echo "[INFO] Number of BAMs: $(wc -l < "$OUT")"
