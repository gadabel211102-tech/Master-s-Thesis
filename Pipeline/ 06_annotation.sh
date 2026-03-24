#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Script 06: Functional Annotation with Ensembl VEP
# =============================================================================
# Purpose:
#   Annotate each filtered VCF with transcript consequence, impact prediction
#   and external reference resources such as gnomAD, CADD, dbNSFP and COSMIC.
#
# Inputs:
#   - variants.filtered.vcf.gz files produced by the variant-calling stage.
#   - A local VEP cache, reference FASTA and plugin resources.
#
# Outputs:
#   - One *.vep.vcf.gz file per sample plus a matching vep.log file.
#
# Methodological note:
#   Annotation is performed sample by sample but with a shared resource bundle,
#   ensuring that all cohorts are interpreted under the same transcript and
#   database context before the tables are merged in script 07.
# =============================================================================

# --- WSL PERL FIXES (The "Secret Sauce") ---
unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
export PERL5LIB=""

# --- Global resources and reference data ---
ROOT_DIR="/home/gadeaalonsoj/tfm" 
REF_FA="/home/gadeaalonsoj/tfm/.vep/homo_sapiens/112_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
CACHE_DIR="/home/gadeaalonsoj/tfm/.vep"
THREADS=4

# 1. Ensure we are in the VEP environment
# If you are running this from VS Code terminal, make sure you ran 'conda activate vep_env' first
log(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }

log "Starting global annotation for all cohorts..."
log ""
log "Enhanced Annotations Included:"
log "  - CADD: Combined Annotation Dependent Depletion (deleteriousness score)"
log "  - SIFT: Sorting Intolerant From Tolerant (functional prediction)"
log "  - PolyPhen2: Polymorphism Phenotyping v2 (pathogenicity prediction)"
log "  - REVEL: Rare Exome Variant Ensemble Learner"
log "  - DANN: Deep Neural Network pathogenicity score"
log "  - gnomAD: Population allele frequencies"
log "  - COSMIC: Catalogue of Somatic Mutations in Cancer"
log ""

# Iterate through every cohort/tissue/sample VCF produced upstream so that
# annotation is applied uniformly across the entire study.
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

    # Run VEP with a broad annotation profile so the downstream statistical
    # scripts can access consequence class, protein change, population frequency
    # context and pathogenicity predictors from one harmonised source.
    vep -i "${VCF}" \
        -o "${OUT_VCF}" \
        --format vcf --vcf --compress_output bgzip --force_overwrite \
        --dir_cache "${CACHE_DIR}" \
        --fasta "${REF_FA}" \
        --species homo_sapiens --assembly GRCh38 --offline \
        --everything --pick --fork "${THREADS}" \
        --plugin CADD,"${CACHE_DIR}/CADD_GRCh38_whole_genome_SNVs.tsv.gz" \
        --plugin dbNSFP,"${CACHE_DIR}/dbNSFP5.3.1a_grch38.gz",SIFT_pred,Polyphen2_HDIV_pred,CADD_phred,REVEL_score,DANN_score \
        --pubmed --domains --regulatory --variant_class \
        2>&1 | tee "${SAMPLE_DIR}/vep.log"
    
    # NOTE: TSV summary removed — script 07 reads VCF files directly.
    # The previous sed 's/|/\t/g' approach was corrupting VEP CSQ fields.
done

log "COMPLETED: All cohorts annotated."
