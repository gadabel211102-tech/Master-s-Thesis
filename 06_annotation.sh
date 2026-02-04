#!/usr/bin/env bash
# Exit script if a command fails (-e), a variable is unset (-u), or a pipe fails (-o pipefail)
set -euo pipefail

# --- WSL PERL FIXES ---
# Clear Perl environment variables to prevent library conflicts within WSL/Conda
unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
export PERL5LIB=""

# --- Configuration ---
# Set the base project directory and paths to genomic reference data
ROOT_DIR="/home/gadeaalonsoj/tfm" 
REF_FA="/home/gadeaalonsoj/.vep/homo_sapiens/112_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
CACHE_DIR="/home/gadeaalonsoj/.vep"
THREADS=4

# Define a logging function that prints messages with a timestamp
log(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }

log "Starting global annotation for all cohorts..."

# 2. Find every filtered VCF in any subfolder and loop through them
find "${ROOT_DIR}" -name "variants.filtered.vcf.gz" | while read -r VCF; do
    
    # Identify the directory and sample name based on the file path
    SAMPLE_DIR=$(dirname "${VCF}")
    SAMPLE_NAME=$(basename "${SAMPLE_DIR}")
    
    # Parse the directory structure to identify the cohort name (e.g., breast, endometrium)
    COHORT=$(echo "$VCF" | awk -F'/' '{print $(NF-3)}') 
    
    # Define the output file path for the annotated VCF
    OUT_VCF="${SAMPLE_DIR}/${SAMPLE_NAME}.vep.vcf.gz"
    
    log "-------------------------------------------------------"
    log "Cohort: ${COHORT} | Sample: ${SAMPLE_NAME}"
    log "-------------------------------------------------------"

    # 3. Run VEP (Variant Effect Predictor)
    vep -i "${VCF}" \
        -o "${OUT_VCF}" \
        --format vcf --vcf --compress_output bgzip --force_overwrite \
        --dir_cache "${CACHE_DIR}" \
        --fasta "${REF_FA}" \
        --species homo_sapiens --assembly GRCh38 --offline \
        --everything --pick --fork "${THREADS}"

    # Create a tabix index for the annotated VCF (required for fast retrieval)
    tabix -f -p vcf "${OUT_VCF}"

    # 4. Generate the summary TSV for your final report
    # Extract variant data and the CSQ (Consequence) field
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%CSQ\n' "${OUT_VCF}" \
        | sed 's/|/\t/g' \
        > "${SAMPLE_DIR}/${SAMPLE_NAME}_vep_summary.tsv"
done

log "COMPLETED: All cohorts annotated."
