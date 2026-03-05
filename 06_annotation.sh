#!/usr/bin/env bash
# Runs VEP (Variant Effect Predictor) on every filtered VCF produced by script 05,
# annotating each variant with functional consequence, pathogenicity scores, and
# population frequency data. Output is a bgzip-compressed VCF per sample.
set -euo pipefail

# --- WSL Perl environment fix ------------------------------------------------
# Under Windows Subsystem for Linux, Conda can leave PERL5LIB pointing at Windows
# paths that are invisible inside the Linux filesystem, causing VEP to fail at
# import. Unsetting these variables forces Perl to use only the paths defined
# within the active Conda environment (vep_env).
unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
export PERL5LIB=""

# --- Configuration -----------------------------------------------------------
# NOTE: activate the vep_env Conda environment before running this script:
#   conda activate vep_env
ROOT_DIR="/home/gadeaalonsoj/tfm"
# GRCh38 reference FASTA used by VEP for HGVS notation and sequence-based plugins
REF_FA="/home/gadeaalonsoj/tfm/.vep/homo_sapiens/112_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
# Offline VEP cache directory (pre-downloaded; avoids network calls during the run)
CACHE_DIR="/home/gadeaalonsoj/tfm/.vep"
THREADS=4

log(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }

log "Starting variant annotation for all cohorts..."
log ""
log "Annotation sources:"
log "  CADD       - Combined Annotation Dependent Depletion (deleteriousness)"
log "  SIFT       - Sorting Intolerant From Tolerant (missense impact)"
log "  PolyPhen2  - Polymorphism Phenotyping v2 (pathogenicity)"
log "  REVEL      - Rare Exome Variant Ensemble Learner (rare variant score)"
log "  DANN       - Deep neural network deleteriousness score"
log "  gnomAD     - Population allele frequencies (via --everything)"
log "  COSMIC     - Catalogue of Somatic Mutations in Cancer (via --everything)"
log ""

# --- Discover all filtered VCFs and annotate each ----------------------------
# find locates variants.filtered.vcf.gz in any subdirectory under ROOT_DIR,
# covering all cohorts (breast, endometrium, etc.) without hardcoding paths
find "${ROOT_DIR}" -name "variants.filtered.vcf.gz" | while read -r VCF; do

    SAMPLE_DIR=$(dirname "${VCF}")
    SAMPLE_NAME=$(basename "${SAMPLE_DIR}")

    # Extract the cohort label from the directory path; assumes the structure:
    # <ROOT_DIR>/<cohort>/<tissue_type>/<sample>/variants.filtered.vcf.gz
    COHORT=$(echo "${VCF}" | awk -F'/' '{print $(NF-3)}')

    OUT_VCF="${SAMPLE_DIR}/${SAMPLE_NAME}.vep.vcf.gz"

    log "-------------------------------------------------------"
    log "Cohort: ${COHORT} | Sample: ${SAMPLE_NAME}"
    log "-------------------------------------------------------"

    # Run VEP with functional annotations and pathogenicity plugins.
    #
    # Key flags:
    #   --everything        enables all standard annotation fields including
    #                       gnomAD frequencies, COSMIC, regulatory regions, etc.
    #   --pick              selects one consequence per variant (most severe
    #                       transcript) to avoid duplicated rows in downstream tables
    #   --offline           uses the local cache only; no Ensembl REST calls
    #   --plugin CADD       appends pre-computed CADD phred scores from a tabix-
    #                       indexed file; scores >= 20 indicate the top 1% most
    #                       deleterious variants in the genome
    #   --plugin dbNSFP     retrieves SIFT, PolyPhen2, REVEL, and DANN predictions
    #                       from a single pre-indexed database file; these ensemble
    #                       scores are important for classifying missense variants
    #   --pubmed            links to relevant PubMed publications where available
    #   --domains           annotates protein domain overlap (Pfam, PROSITE, etc.)
    #   --regulatory        reports overlap with regulatory features (promoters,
    #                       enhancers) from Ensembl Regulatory Build
    #   --variant_class     adds SO sequence ontology class (SNV, insertion, etc.)
    #
    # Output is written as bgzip-compressed VCF so tabix can index it for script 07.
    # stderr and stdout are both captured to vep.log for debugging.
    #
    # NOTE: TSV export via sed 's/|/\t/g' was removed — pipe-delimited CSQ fields
    # cannot be safely split with sed; script 07 parses the VCF CSQ field directly.
    vep \
        -i "${VCF}" \
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

done

log "Annotation complete for all cohorts."
