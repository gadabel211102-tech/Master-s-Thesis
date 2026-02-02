#!/usr/bin/env bash
set -euo pipefail

# --- WSL PERL FIXES (The "Secret Sauce") ---
unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
export PERL5LIB=""

# --- Configuration ---
ROOT_DIR="/home/gadeaalonsoj/tfm" 
REF_FA="/home/gadeaalonsoj/.vep/homo_sapiens/112_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
CACHE_DIR="/home/gadeaalonsoj/.vep"
THREADS=4

# 1. Ensure we are in the VEP environment
# If you are running this from VS Code terminal, make sure you ran 'conda activate vep_env' first
log(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }

log "Starting global annotation for all cohorts..."

# 2. Find every filtered VCF in any subfolder (breast, endometrium, etc.)
find "${ROOT_DIR}" -name "variants.filtered.vcf.gz" | while read -r VCF; do
    
    # Get folder names for logging
    SAMPLE_DIR=$(dirname "${VCF}")
    SAMPLE_NAME=$(basename "${SAMPLE_DIR}")
    
    # Extract the cohort name (e.g., breast or endometrium)
    # This assumes path like: tfm/breast/tumour/SampleA/...
    COHORT=$(echo "$VCF" | awk -F'/' '{print $(NF-3)}') 
    
    OUT_VCF="${SAMPLE_DIR}/${SAMPLE_NAME}.vep.vcf.gz"
    
    log "-------------------------------------------------------"
    log "Cohort: ${COHORT} | Sample: ${SAMPLE_NAME}"
    log "-------------------------------------------------------"

    # 3. Run VEP
    vep -i "${VCF}" \
        -o "${OUT_VCF}" \
        --format vcf --vcf --compress_output bgzip --force_overwrite \
        --dir_cache "${CACHE_DIR}" \
        --fasta "${REF_FA}" \
        --species homo_sapiens --assembly GRCh38 --offline \
        --everything --pick --fork "${THREADS}"

    tabix -f -p vcf "${OUT_VCF}"

    # 4. Generate the summary TSV for your final report
    # This extracts: Chrom, Pos, Ref, Alt, Gene, Consequence, Protein Change, and Impact
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%CSQ\n' "${OUT_VCF}" \
        | sed 's/|/\t/g' \
        > "${SAMPLE_DIR}/${SAMPLE_NAME}_vep_summary.tsv"
done

log "COMPLETED: All cohorts annotated."