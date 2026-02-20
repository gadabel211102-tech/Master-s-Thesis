#!/bin/bash
# =============================================================================
# Script 15: Merge all sample VCFs → Phase with BEAGLE → Extract phased GTs
# =============================================================================
#
# PURPOSE:
#   Constructs a cohort-wide haplotype phasing pipeline for the GSDMB locus
#   on chromosome 17. Single-sample VCFs from all four cohort groups (breast
#   normal, breast tumour, endometrium normal, endometrium tumour) are merged
#   into a single multi-sample VCF, filtered to biallelic chr17 SNPs, and
#   statistically phased using BEAGLE 5.4.
#
# STUDY DESIGN NOTE:
#   Samples are unpaired (tumour and healthy tissue from different patients).
#   Phasing is performed at the population level across all groups jointly,
#   which maximises the number of samples available to BEAGLE for LD-based
#   phasing. Only QC-passed samples are present in dna_calls/ directories
#   (script 02 copies only PASS BAMs forward), so no additional QC filtering
#   is required here.
#
# REFERENCE PANEL PREPARATION (run once before this script):
#   # Download 1000 Genomes Phase 3 chr17 (GRCh38 / hg38)
#   wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz
#   wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
#
#   # Subset to GSDMB locus only (speeds up BEAGLE significantly)
#   mkdir -p ~/tfm/ref
#   bcftools view \
#       --regions chr17:39800000-39990000 \
#       --output-type z \
#       --output ~/tfm/ref/ref_panel_chr17_gsdmb_GRCh38.vcf.gz \
#       CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz
#   bcftools index -t ~/tfm/ref/ref_panel_chr17_gsdmb_GRCh38.vcf.gz
#
#   # Download BEAGLE genetic map (GRCh38)
#   wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
#   unzip plink.GRCh38.map.zip -d ~/tfm/ref/
#
# OUTPUTS:
#   gsdmb_final_results/19_haplotype_phased/
#     cohort_merged.vcf.gz     — multi-sample merged VCF (all variants)
#     cohort_snps.vcf.gz       — biallelic chr17 SNPs only (pre-phasing)
#     cohort_phased.vcf.gz     — BEAGLE-phased haplotypes
#     phased_genotypes.tsv     — phased GT matrix (SNPs x samples) for R
#     sample_metadata.tsv      — sample group assignments for R
#     sample_order.txt         — ordered sample list matching TSV columns
#
# USAGE:
#   cd ~/tfm
#   bash 15_haplotypes.sh
#
# DEPENDENCIES:
#   bcftools >= 1.15, Java >= 11, BEAGLE 5.4
#
# Date:   2026-02-12
# =============================================================================

set -euo pipefail
# set -e        : exit immediately if any command returns non-zero
# set -u        : treat unset variables as errors (prevents silent bugs)
# set -o pipefail: pipeline fails if any stage within it fails

# =============================================================================
# CONFIGURATION
# Edit these paths if your directory structure differs
# =============================================================================

BASE="/home/gadeaalonsoj/tfm"
RESULTS="${BASE}/gsdmb_final_results"
HAPLO_DIR="${RESULTS}/19_haplotype_phased"
BEAGLE="${BASE}/beagle.jar"
THREADS=4

# 1000G Phase 3 reference panel (chr17:39800000-39950000 subset)
# See header comments for how to prepare this file
REF_PANEL="${BASE}/ref/ref_panel_chr17_gsdmb_GRCh38.vcf.gz"

# BEAGLE genetic map for GRCh37 chromosome 17
GENETIC_MAP="${BASE}/ref/plink.chr17.GRCh38.map"
# BEAGLE phasing parameters
NE=10000       # Effective population size; 10,000 is standard for European ancestry
ITERATIONS=20  # More iterations than BEAGLE default (12) for improved accuracy

mkdir -p "${HAPLO_DIR}"

echo "======================================================================"
echo "SCRIPT 15: BEAGLE PHASING PIPELINE"
echo "======================================================================"
echo ""
echo "Configuration:"
echo "  Base directory:  ${BASE}"
echo "  Output:          ${HAPLO_DIR}"
echo "  BEAGLE threads:  ${THREADS}"
echo "  ne:              ${NE}"
echo "  Iterations:      ${ITERATIONS}"
echo ""

# =============================================================================
# PRE-FLIGHT CHECKS
# Verify all required files exist before starting so the pipeline fails
# cleanly at the top rather than partway through a long run.
# =============================================================================
echo "PRE-FLIGHT: Checking required files"
echo "----------------------------------------------------------------------"

# BEAGLE jar is mandatory
if [ ! -f "${BEAGLE}" ]; then
    echo "ERROR: BEAGLE jar not found at ${BEAGLE}"
    echo "       Download from: https://faculty.washington.edu/browning/beagle/beagle.html"
    exit 1
fi
echo "  ✓ BEAGLE jar found"

# Reference panel is optional but strongly recommended
# Script will continue without it but phasing accuracy will be reduced
if [ ! -f "${REF_PANEL}" ]; then
    echo "  ⚠ WARNING: Reference panel not found at ${REF_PANEL}"
    echo "    Phasing will proceed without a reference panel (cohort LD only)."
    echo "    See script header for download instructions."
    USE_REF=false
else
    echo "  ✓ Reference panel found: ${REF_PANEL}"
    USE_REF=true
fi

# Genetic map is optional but recommended
if [ ! -f "${GENETIC_MAP}" ]; then
    echo "  ⚠ WARNING: Genetic map not found at ${GENETIC_MAP}"
    echo "    Phasing will assume uniform recombination rate across chr17."
    USE_MAP=false
else
    echo "  ✓ Genetic map found: ${GENETIC_MAP}"
    USE_MAP=true
fi

echo ""

# =============================================================================
# STEP 1: Collect all filtered VCFs and rename sample columns
# =============================================================================
# Each single-sample VCF from script 05 has its sample column named
# generically (often "SAMPLE" or the IonCode barcode). When bcftools merge
# combines them, every sample column must be uniquely named or the merge
# will fail. We rename each sample to its folder name (which is guaranteed
# unique) using bcftools reheader before merging.
#
# dna_calls/ directories contain only QC-passed samples because script 02
# copies BAMs to pass_bams/ only, and script 05 processes pass_bams/.
# No additional QC filtering is therefore required at this step.
# =============================================================================
echo "STEP 1: Renaming samples and collecting VCFs"
echo "----------------------------------------------------------------------"

SAMPLE_LIST="${HAPLO_DIR}/sample_list.txt"
RENAMED_DIR="${HAPLO_DIR}/renamed_vcfs"
mkdir -p "${RENAMED_DIR}"
> "${SAMPLE_LIST}"   # Truncate or create the sample list file

for cohort in breast endometrium; do
    for tissue in normal tumour; do
        CALLS_DIR="${BASE}/${cohort}/${tissue}/dna_calls"

        # Skip cohort/tissue combinations whose directory does not exist
        if [ ! -d "${CALLS_DIR}" ]; then
            echo "  ⚠ Not found: ${CALLS_DIR} — skipping"
            continue
        fi

        for sample_dir in "${CALLS_DIR}"/*/; do
            sample=$(basename "${sample_dir}")
            vcf="${sample_dir}/variants.filtered.vcf.gz"

            # Skip samples that have no filtered VCF (e.g. variant calling failed)
            if [ ! -f "${vcf}" ]; then
                echo "  ⚠ No VCF for ${sample} — skipping"
                continue
            fi

            out_vcf="${RENAMED_DIR}/${sample}.vcf.gz"

            # Write new sample name list (one name per line — bcftools reheader -s format)
            echo "${sample}" > "${HAPLO_DIR}/rename_${sample}.txt"

            # Rename sample column and recompress; -t creates a .tbi tabix index
            bcftools reheader \
                -s "${HAPLO_DIR}/rename_${sample}.txt" \
                -o "${out_vcf}" \
                "${vcf}"
            bcftools index -t "${out_vcf}"

            echo "${out_vcf}" >> "${SAMPLE_LIST}"
            echo "  ✓ ${cohort}/${tissue}/${sample}"
        done
    done
done

N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo ""
echo "✓ ${N_SAMPLES} samples prepared"

# Abort early if no samples were found — nothing to phase
if [ "${N_SAMPLES}" -eq 0 ]; then
    echo "ERROR: No sample VCFs found. Check that dna_calls/ directories"
    echo "       exist and contain variants.filtered.vcf.gz files."
    exit 1
fi

# =============================================================================
# STEP 2: Merge all single-sample VCFs into one cohort-wide VCF
# =============================================================================
# bcftools merge combines multiple single-sample VCFs into one multi-sample
# VCF with one column per sample.
#
# --missing-to-ref fills positions where a sample has no call with the
# reference genotype (0/0) rather than missing (./.). This is appropriate
# for targeted amplicon data where good coverage at a position with no
# variant call reliably indicates the reference allele. It is not appropriate
# for whole-genome data with variable coverage.
# =============================================================================
echo ""
echo "STEP 2: Merging into cohort VCF"
echo "----------------------------------------------------------------------"

MERGED_VCF="${HAPLO_DIR}/cohort_merged.vcf.gz"

bcftools merge \
    --file-list "${SAMPLE_LIST}" \
    --missing-to-ref \
    --output-type z \
    --output "${MERGED_VCF}" \
    --threads "${THREADS}"

bcftools index -t "${MERGED_VCF}"

N_VARS=$(bcftools view -H "${MERGED_VCF}" | wc -l)
echo "✓ Merged VCF: ${N_VARS} variants across ${N_SAMPLES} samples"

# =============================================================================
# STEP 3: Filter to biallelic SNPs on chr17
# =============================================================================
# We restrict to:
#   chr17 only       — our panel targets the GSDMB locus on chr17
#   SNPs only        — indels and MNPs excluded; haplo.stats requires SNPs
#   Biallelic sites  — multiallelic sites excluded for simplicity; rare in
#                      a targeted amplicon panel
#
# IMPORTANT: gnomAD NFE AF > 1% filtering is intentionally NOT applied here.
# We phase all chr17 biallelic SNPs to give BEAGLE maximum LD information,
# which improves phasing accuracy. The gnomAD frequency filter is applied
# downstream in 15_haplo_stats.r using positions from the annotated Excel
# report (matching the exact SNP set from script 11).
# =============================================================================
echo ""
echo "STEP 3: Filtering to biallelic SNPs on chr17"
echo "----------------------------------------------------------------------"

SNP_VCF="${HAPLO_DIR}/cohort_snps.vcf.gz"

bcftools view \
    --regions chr17 \
    --type snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --output-type z \
    --output "${SNP_VCF}" \
    "${MERGED_VCF}"

bcftools index -t "${SNP_VCF}"

N_SNPS=$(bcftools view -H "${SNP_VCF}" | wc -l)
echo "✓ Filtered SNP VCF: ${N_SNPS} biallelic SNPs on chr17"

# Common failure: contig named "17" instead of "chr17"
if [ "${N_SNPS}" -eq 0 ]; then
    echo "ERROR: No SNPs remain after filtering."
    echo "       Check your contig naming convention:"
    echo "       Run: bcftools view -H ${MERGED_VCF} | head -3 | cut -f1"
    echo "       If contigs are named '17' not 'chr17', change --regions to '17'"
    exit 1
fi

# =============================================================================
# STEP 4: Statistical phasing with BEAGLE 5.4
# =============================================================================
# BEAGLE uses a hidden Markov model to phase genotypes probabilistically,
# exploiting linkage disequilibrium patterns in the data.
#
# With 1000G reference panel (recommended):
#   - Borrows haplotype structure from ~2,500 phased European samples
#   - Substantially improves accuracy over cohort-only phasing, especially
#     for rarer variants and smaller cohorts such as ours (~138 samples)
#
# With genetic map (recommended):
#   - Uses chromosome-specific recombination rates from HapMap
#   - Improves phasing at recombination hotspots
#
# Key parameters:
#   ne=10000     : effective population size, standard for European ancestry
#   iterations=20: more EM iterations than default (12) for better accuracy
# =============================================================================
echo ""
echo "STEP 4: Statistical phasing with BEAGLE 5.4"
echo "----------------------------------------------------------------------"

PHASED_PREFIX="${HAPLO_DIR}/cohort_phased"

# Build BEAGLE command dynamically, adding ref and map only if files exist
BEAGLE_CMD="java -Xmx4g -jar ${BEAGLE} \
    gt=${SNP_VCF} \
    out=${PHASED_PREFIX} \
    nthreads=${THREADS} \
    iterations=${ITERATIONS} \
    ne=${NE}"

if [ "${USE_REF}" = true ]; then
    BEAGLE_CMD="${BEAGLE_CMD} ref=${REF_PANEL}"
    echo "  Using 1000G reference panel for improved phasing accuracy"
else
    echo "  ⚠ Running without reference panel — cohort LD only"
fi

if [ "${USE_MAP}" = true ]; then
    BEAGLE_CMD="${BEAGLE_CMD} map=${GENETIC_MAP}"
    echo "  Using genetic map for recombination-aware phasing"
else
    echo "  ⚠ Running without genetic map — uniform recombination assumed"
fi

echo ""
eval "${BEAGLE_CMD}"

bcftools index -t "${PHASED_PREFIX}.vcf.gz"

N_PHASED=$(bcftools view -H "${PHASED_PREFIX}.vcf.gz" | wc -l)
echo ""
echo "✓ Phasing complete: ${N_PHASED} SNPs phased across ${N_SAMPLES} samples"
echo "  Output: ${PHASED_PREFIX}.vcf.gz"

# =============================================================================
# STEP 5: Extract phased genotypes to TSV for R
# =============================================================================
# The R script reads a plain TSV rather than a VCF directly, which avoids
# a VCF parsing dependency in R and makes the data straightforward to inspect.
#
# Output format (tab-separated):
#   CHROM   POS       REF  ALT  sample1  sample2  ...
#   chr17   39905000  A    G    0|0      0|1      1|1  ...
#
# The '|' separator in genotypes indicates phased haplotypes:
#   0|1 = ref allele on haplotype 1, alt allele on haplotype 2
#   1|0 = alt on haplotype 1, ref on haplotype 2
# These are distinct haplotype assignments and treated differently by haplo.stats.
# =============================================================================
echo ""
echo "STEP 5: Extracting phased genotypes to TSV"
echo "----------------------------------------------------------------------"

PHASED_TSV="${HAPLO_DIR}/phased_genotypes.tsv"

# Save the ordered list of sample names exactly as they appear in the phased VCF
# The order here must match the column order in the TSV header
bcftools query -l "${PHASED_PREFIX}.vcf.gz" > "${HAPLO_DIR}/sample_order.txt"

# Extract one row per variant with all sample genotypes
# %GT retrieves the phased genotype field (e.g. 0|1) for each sample
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
    "${PHASED_PREFIX}.vcf.gz" \
    > "${PHASED_TSV}"

# Prepend a header row using sample names from sample_order.txt
HEADER="CHROM\tPOS\tREF\tALT\t$(paste -s -d'\t' "${HAPLO_DIR}/sample_order.txt")"
{ echo -e "${HEADER}"; cat "${PHASED_TSV}"; } > "${PHASED_TSV}.tmp"
mv "${PHASED_TSV}.tmp" "${PHASED_TSV}"

N_ROWS=$(( $(wc -l < "${PHASED_TSV}") - 1 ))  # Subtract header line
echo "✓ Phased genotype table: ${N_ROWS} SNPs × ${N_SAMPLES} samples"
echo "  Saved to: ${PHASED_TSV}"

# =============================================================================
# STEP 6: Create sample metadata file for R
# =============================================================================
# Maps each sample to its cohort, tissue, and group label. Used by
# 15_haplo_stats.r to define case/control groups for association testing.
#
# Columns:
#   Sample  — matches column names in phased_genotypes.tsv
#   Cohort  — Breast or Endometrium
#   Tissue  — Tumour or Healthy
#   Group   — combined label e.g. Breast_Tumour (used for grouped plots)
# =============================================================================
echo ""
echo "STEP 6: Creating sample metadata file"
echo "----------------------------------------------------------------------"

META_TSV="${HAPLO_DIR}/sample_metadata.tsv"
echo -e "Sample\tCohort\tTissue\tGroup" > "${META_TSV}"

for cohort in breast endometrium; do
    for tissue in normal tumour; do
        # Map directory names to display labels used throughout the pipeline
        tissue_label=$([ "${tissue}" = "normal" ] && echo "Healthy" || echo "Tumour")
        cohort_label="${cohort^}"  # Parameter expansion: capitalise first letter

        CALLS_DIR="${BASE}/${cohort}/${tissue}/dna_calls"
        if [ ! -d "${CALLS_DIR}" ]; then continue; fi

        for sample_dir in "${CALLS_DIR}"/*/; do
            sample=$(basename "${sample_dir}")
            echo -e "${sample}\t${cohort_label}\t${tissue_label}\t${cohort_label}_${tissue_label}"
        done
    done
done >> "${META_TSV}"

N_META=$(( $(wc -l < "${META_TSV}") - 1 ))
echo "✓ Metadata written for ${N_META} samples"

# Sanity check: metadata count should match the phased VCF sample count
# A mismatch indicates samples were skipped in Step 1 or Step 5
if [ "${N_META}" -ne "${N_SAMPLES}" ]; then
    echo "  ⚠ WARNING: Metadata has ${N_META} samples but VCF has ${N_SAMPLES}."
    echo "    Check sample_list.txt and the dna_calls/ directories for discrepancies."
fi

# =============================================================================
# COMPLETION SUMMARY
# =============================================================================
echo ""
echo "======================================================================"
echo "✓ PHASING PIPELINE COMPLETE"
echo "======================================================================"
echo ""
echo "Output files in: ${HAPLO_DIR}/"
echo "  cohort_merged.vcf.gz     — multi-sample merged VCF (all variants)"
echo "  cohort_snps.vcf.gz       — biallelic chr17 SNPs (pre-phasing input)"
echo "  cohort_phased.vcf.gz     — BEAGLE-phased haplotypes (all SNPs)"
echo "  phased_genotypes.tsv     — phased GT matrix (SNPs × samples)"
echo "  sample_metadata.tsv      — sample group assignments for R"
echo "  sample_order.txt         — sample column order matching TSV"
echo ""
echo "Phasing configuration:"
echo "  Reference panel used: ${USE_REF}"
echo "  Genetic map used:     ${USE_MAP}"
echo ""
echo "NOTE: gnomAD NFE AF > 1%% filtering is applied in 15_haplo_stats.r."
echo "      The phased VCF contains all biallelic chr17 SNPs."
echo ""
echo "Next step: Rscript 15_haplo_stats.r"
echo ""
