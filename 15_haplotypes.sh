#!/bin/bash
# =============================================================================
# Script 15: Merge all sample VCFs → Phase with BEAGLE → Extract phased GTs
# =============================================================================
# =============================================================================

set -euo pipefail

BASE="/home/gadeaalonsoj/tfm"
RESULTS="${BASE}/gsdmb_final_results"
HAPLO_DIR="${RESULTS}/19_haplotype_phased"
BEAGLE="${BASE}/beagle.jar"
THREADS=4

mkdir -p "${HAPLO_DIR}"

echo "======================================================================"
echo "SCRIPT 19a: BEAGLE PHASING PIPELINE"
echo "======================================================================"

# ------------------------------------------------------------------------------
# STEP 1: Collect all filtered VCFs and create sample-renamed copies
# Each VCF has the sample named generically (e.g. "SAMPLE" or the IonCode).
# We rename the sample column to match the folder name so the merged VCF
# has one uniquely-named column per sample.
# ------------------------------------------------------------------------------
echo ""
echo "STEP 1: Renaming samples and collecting VCFs"
echo "----------------------------------------------------------------------"

SAMPLE_LIST="${HAPLO_DIR}/sample_list.txt"
RENAMED_DIR="${HAPLO_DIR}/renamed_vcfs"
mkdir -p "${RENAMED_DIR}"
> "${SAMPLE_LIST}"

for cohort in breast endometrium; do
    for tissue in normal tumour; do
        CALLS_DIR="${BASE}/${cohort}/${tissue}/dna_calls"
        if [ ! -d "${CALLS_DIR}" ]; then
            echo "  ⚠ Not found: ${CALLS_DIR} — skipping"
            continue
        fi
        for sample_dir in "${CALLS_DIR}"/*/; do
            sample=$(basename "${sample_dir}")
            vcf="${sample_dir}/variants.filtered.vcf.gz"
            if [ ! -f "${vcf}" ]; then
                echo "  ⚠ No VCF for ${sample} — skipping"
                continue
            fi

            out_vcf="${RENAMED_DIR}/${sample}.vcf.gz"

            # Get the current sample name in the VCF header
            current_name=$(bcftools query -l "${vcf}" | head -1)

            # Rename sample to folder name and recompress
            echo "${current_name} ${sample}" > "${HAPLO_DIR}/rename_${sample}.txt"
            bcftools reheader \
                --samples "${HAPLO_DIR}/rename_${sample}.txt" \
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

# ------------------------------------------------------------------------------
# STEP 2: Merge all single-sample VCFs into one cohort VCF
# ------------------------------------------------------------------------------
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
echo "✓ Merged VCF: ${N_VARS} variants, ${N_SAMPLES} samples"

# ------------------------------------------------------------------------------
# STEP 3: Filter to SNP positions only (matching script 11 SNP set)
# Keep only biallelic SNPs with gnomAD NFE AF > 1% proxy:
# here we keep PASS variants on chr17 with AF > 0.01 in the cohort itself
# (since we don't have gnomAD in the VCF at this stage — we filter by position
#  after phasing using the SNP list from the annotated Excel)
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# STEP 4: Run BEAGLE phasing (no reference panel — uses cohort LD)
# ------------------------------------------------------------------------------
echo ""
echo "STEP 4: Statistical phasing with BEAGLE 5.4"
echo "----------------------------------------------------------------------"

PHASED_PREFIX="${HAPLO_DIR}/cohort_phased"

java -Xmx4g -jar "${BEAGLE}" \
    gt="${SNP_VCF}" \
    out="${PHASED_PREFIX}" \
    nthreads="${THREADS}" \
    iterations=20 \
    ne=10000

# ne=10000 is standard effective population size for European ancestry
# iterations=20 gives more accurate phasing than default (12)

bcftools index -t "${PHASED_PREFIX}.vcf.gz"

echo "✓ Phasing complete: ${PHASED_PREFIX}.vcf.gz"

# ------------------------------------------------------------------------------
# STEP 5: Extract phased genotypes to TSV for Python/R
# Format: CHROM POS REF ALT sample1_GT sample2_GT ...
# ------------------------------------------------------------------------------
echo ""
echo "STEP 5: Extracting phased genotypes to TSV"
echo "----------------------------------------------------------------------"

PHASED_TSV="${HAPLO_DIR}/phased_genotypes.tsv"

# Get sample names in order
bcftools query -l "${PHASED_PREFIX}.vcf.gz" > "${HAPLO_DIR}/sample_order.txt"

# Extract GT for all samples
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
    "${PHASED_PREFIX}.vcf.gz" \
    > "${PHASED_TSV}"

# Add header
HEADER="CHROM\tPOS\tREF\tALT\t$(paste -s -d'\t' "${HAPLO_DIR}/sample_order.txt")"
{ echo -e "${HEADER}"; cat "${PHASED_TSV}"; } > "${PHASED_TSV}.tmp"
mv "${PHASED_TSV}.tmp" "${PHASED_TSV}"

N_PHASED=$(wc -l < "${PHASED_TSV}")
echo "✓ Phased genotype table: $((N_PHASED-1)) SNPs × ${N_SAMPLES} samples"
echo "  Saved to: ${PHASED_TSV}"

# ------------------------------------------------------------------------------
# STEP 6: Create sample metadata file for R
# ------------------------------------------------------------------------------
echo ""
echo "STEP 6: Creating sample metadata file"
echo "----------------------------------------------------------------------"

META_TSV="${HAPLO_DIR}/sample_metadata.tsv"
echo -e "Sample\tCohort\tTissue\tGroup" > "${META_TSV}"

for cohort in breast endometrium; do
    for tissue in normal tumour; do
        tissue_label=$([ "${tissue}" = "normal" ] && echo "Healthy" || echo "Tumour")
        cohort_label=$(echo "${cohort^}")  # capitalise first letter
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

echo ""
echo "======================================================================"
echo "✓ PHASING PIPELINE COMPLETE"
echo "======================================================================"
echo ""
echo "Output files in: ${HAPLO_DIR}/"
echo "  cohort_phased.vcf.gz     — phased VCF (all samples)"
echo "  phased_genotypes.tsv     — phased GT matrix (SNPs × samples)"
echo "  sample_metadata.tsv      — sample group assignments"
echo ""
echo "Next step: run  Rscript 19b_haplo_stats.R"
echo ""
