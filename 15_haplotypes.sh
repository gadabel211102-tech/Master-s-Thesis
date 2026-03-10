#!/bin/bash
# Constructs a cohort-wide haplotype phasing pipeline for the GSDMB locus
# on chromosome 17. Single-sample VCFs from all four cohort groups (breast
# normal, breast tumour, endometrium normal, endometrium tumour) are merged
# into a single multi-sample VCF, filtered to biallelic chr17 SNPs, and
# statistically phased using BEAGLE 5.4.
#
# Samples are unpaired (tumour and healthy tissue from different patients).
# Phasing is performed jointly across all groups to maximise the number of
# samples available to BEAGLE for LD-based inference. Only QC-passed samples
# reach the dna_calls/ directories (script 02 copies only PASS BAMs forward),
# so no additional QC filtering is required here.
#
# REFERENCE PANEL PREPARATION (run once before this script):
#   # Download 1000 Genomes Phase 3 chr17 (GRCh38)
#   wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz
#   wget ...vcf.gz.tbi
#
#   # Subset to GSDMB locus (chr17:39800000-39990000) to speed up BEAGLE
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
# OUTPUTS (gsdmb_final_results/19_haplotype_phased/):
#   cohort_merged.vcf.gz   — multi-sample merged VCF (all variants)
#   cohort_snps.vcf.gz     — biallelic chr17 SNPs only (pre-phasing input)
#   cohort_phased.vcf.gz   — BEAGLE-phased haplotypes
#   phased_genotypes.tsv   — phased GT matrix (SNPs x samples) for R
#   sample_metadata.tsv    — sample group assignments for R
#   sample_order.txt       — ordered sample list matching TSV columns
#
# USAGE:
#   cd ~/tfm && bash 15_haplotypes.sh
#
# DEPENDENCIES:
#   bcftools >= 1.15, Java >= 11, BEAGLE 5.4

set -euo pipefail
# -e: exit on any error; -u: treat unset variables as errors;
# -o pipefail: fail the pipeline if any stage within it fails

# --- Configuration -----------------------------------------------------------

BASE="/home/gadeaalonsoj/tfm"
RESULTS="${BASE}/gsdmb_final_results"
HAPLO_DIR="${RESULTS}/19_haplotype_phased"
BEAGLE="${BASE}/beagle.jar"
THREADS=4

# 1000G Phase 3 reference panel subset to the GSDMB locus — see header for preparation
REF_PANEL="${BASE}/ref/ref_panel_chr17_gsdmb_GRCh38.vcf.gz"

# BEAGLE genetic map for GRCh38 chromosome 17
GENETIC_MAP="${BASE}/ref/plink.chr17.GRCh38.map"

NE=10000       # Effective population size; 10,000 is standard for European ancestry
ITERATIONS=20  # 20 EM iterations (above BEAGLE default of 12) for improved accuracy

mkdir -p "${HAPLO_DIR}"

echo "======================================================================"
echo "SCRIPT 15: BEAGLE PHASING PIPELINE"
echo "======================================================================"
echo ""
echo "  Base:       ${BASE}"
echo "  Output:     ${HAPLO_DIR}"
echo "  Threads:    ${THREADS}"
echo "  ne:         ${NE}"
echo "  Iterations: ${ITERATIONS}"
echo ""

# --- Pre-flight checks -------------------------------------------------------
# Verify required files before starting so the pipeline fails cleanly at the
# top rather than partway through a long run.
echo "Pre-flight checks"
echo "----------------------------------------------------------------------"

if [ ! -f "${BEAGLE}" ]; then
    echo "ERROR: BEAGLE jar not found at ${BEAGLE}"
    echo "       Download from: https://faculty.washington.edu/browning/beagle/beagle.html"
    exit 1
fi
echo "  BEAGLE jar: OK"

# Reference panel is optional but strongly recommended; without it, BEAGLE
# phases using only cohort-internal LD, which is less accurate at small n
if [ ! -f "${REF_PANEL}" ]; then
    echo "  Warning: reference panel not found at ${REF_PANEL}"
    echo "    Phasing will use cohort LD only. See script header for download instructions."
    USE_REF=false
else
    echo "  Reference panel: OK"
    USE_REF=true
fi

# Genetic map is optional; without it, BEAGLE assumes uniform recombination rate
if [ ! -f "${GENETIC_MAP}" ]; then
    echo "  Warning: genetic map not found at ${GENETIC_MAP}"
    echo "    Phasing will assume uniform recombination rate across chr17."
    USE_MAP=false
else
    echo "  Genetic map: OK"
    USE_MAP=true
fi
echo ""

# =============================================================================
# STEP 1: Collect sample VCFs and rename sample columns before merging
# =============================================================================
# Each single-sample VCF from script 05 has its sample column named
# generically (often "SAMPLE" or the IonCode barcode). When bcftools merge
# combines them, every sample column must be uniquely named or the merge will
# fail. We rename each sample to its folder name (guaranteed unique) using
# bcftools reheader before merging.
echo "STEP 1: Renaming samples and collecting VCFs"
echo "----------------------------------------------------------------------"

SAMPLE_LIST="${HAPLO_DIR}/sample_list.txt"
RENAMED_DIR="${HAPLO_DIR}/renamed_vcfs"
mkdir -p "${RENAMED_DIR}"
> "${SAMPLE_LIST}"   # Truncate or create the file

for cohort in breast endometrium; do
    for tissue in normal tumour; do
        CALLS_DIR="${BASE}/${cohort}/${tissue}/dna_calls"

        if [ ! -d "${CALLS_DIR}" ]; then
            echo "  Not found: ${CALLS_DIR} — skipping"
            continue
        fi

        for sample_dir in "${CALLS_DIR}"/*/; do
            sample=$(basename "${sample_dir}")
            vcf="${sample_dir}/variants.filtered.vcf.gz"

            if [ ! -f "${vcf}" ]; then
                echo "  No VCF for ${sample} — skipping"
                continue
            fi

            out_vcf="${RENAMED_DIR}/${sample}.vcf.gz"

            # Write a one-line rename file; bcftools reheader -s expects one
            # new sample name per line, matching the sample column order in the VCF
            echo "${sample}" > "${HAPLO_DIR}/rename_${sample}.txt"

            bcftools reheader \
                -s "${HAPLO_DIR}/rename_${sample}.txt" \
                -o "${out_vcf}" \
                "${vcf}"
            bcftools index -t "${out_vcf}"

            echo "${out_vcf}" >> "${SAMPLE_LIST}"
            echo "  ${cohort}/${tissue}/${sample}"
        done
    done
done

N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo ""
echo "  ${N_SAMPLES} samples prepared"

if [ "${N_SAMPLES}" -eq 0 ]; then
    echo "ERROR: No sample VCFs found. Check that dna_calls/ directories"
    echo "       exist and contain variants.filtered.vcf.gz files."
    exit 1
fi

# =============================================================================
# STEP 2: Merge all single-sample VCFs into one cohort-wide multi-sample VCF
# =============================================================================
# --missing-to-ref fills positions where a sample has no call with the
# reference genotype (0/0) rather than missing (./.). This is appropriate
# for targeted amplicon data: good coverage at a position with no variant call
# reliably indicates the reference allele. It is NOT appropriate for
# whole-genome data with variable coverage, where absence of a call may
# reflect genuine missing data rather than a confident reference call.
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
echo "  Merged VCF: ${N_VARS} variants across ${N_SAMPLES} samples"

# =============================================================================
# STEP 3: Filter to biallelic SNPs on chr17
# =============================================================================
# We restrict to:
#   chr17 only      — the panel targets the GSDMB locus on chr17
#   SNPs only       — indels and MNPs excluded; haplo.stats requires SNPs
#   Biallelic sites — multiallelic sites excluded; rare in a targeted panel
#
# IMPORTANT: gnomAD NFE AF > 1% filtering is intentionally NOT applied here.
# We phase all biallelic chr17 SNPs to give BEAGLE the maximum number of
# variants for LD-based inference — more variants improve phasing accuracy
# even if some will be excluded from the downstream association analysis.
# The frequency filter is applied in 15_haplo_stats.r using positions from
# the annotated Excel report, matching the SNP set from script 11.
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
echo "  Biallelic chr17 SNPs: ${N_SNPS}"

# A count of zero usually means contigs are named "17" rather than "chr17"
if [ "${N_SNPS}" -eq 0 ]; then
    echo "ERROR: No SNPs remain after filtering."
    echo "       Check contig naming: bcftools view -H ${MERGED_VCF} | head -3 | cut -f1"
    echo "       If contigs are named '17' not 'chr17', change --regions to '17'"
    exit 1
fi

# =============================================================================
# STEP 4: Statistical phasing with BEAGLE 5.4
# =============================================================================
# BEAGLE uses a hidden Markov model to phase genotypes probabilistically,
# exploiting linkage disequilibrium patterns across the input samples.
#
# With 1000G reference panel (recommended):
#   Borrows haplotype structure from ~2,500 phased European samples.
#   Substantially improves accuracy for smaller cohorts such as ours (~138
#   samples), particularly for rarer variants where cohort-internal LD alone
#   is insufficient to resolve phase.
#
# With genetic map (recommended):
#   Uses chromosome-specific recombination rates rather than assuming a
#   uniform rate, improving accuracy at recombination hotspots.
#
# BEAGLE command is built dynamically so it degrades gracefully if the
# reference panel or genetic map is absent.
echo ""
echo "STEP 4: Statistical phasing with BEAGLE 5.4"
echo "----------------------------------------------------------------------"

PHASED_PREFIX="${HAPLO_DIR}/cohort_phased"

BEAGLE_CMD="java -Xmx4g -jar ${BEAGLE} \
    gt=${SNP_VCF} \
    out=${PHASED_PREFIX} \
    nthreads=${THREADS} \
    iterations=${ITERATIONS} \
    ne=${NE}"

if [ "${USE_REF}" = true ]; then
    BEAGLE_CMD="${BEAGLE_CMD} ref=${REF_PANEL}"
    echo "  Using 1000G reference panel"
else
    echo "  Warning: running without reference panel — cohort LD only"
fi

if [ "${USE_MAP}" = true ]; then
    BEAGLE_CMD="${BEAGLE_CMD} map=${GENETIC_MAP}"
    echo "  Using genetic map for recombination-aware phasing"
else
    echo "  Warning: running without genetic map — uniform recombination assumed"
fi

echo ""
eval "${BEAGLE_CMD}"

bcftools index -t "${PHASED_PREFIX}.vcf.gz"

N_PHASED=$(bcftools view -H "${PHASED_PREFIX}.vcf.gz" | wc -l)
echo ""
echo "  Phasing complete: ${N_PHASED} SNPs phased across ${N_SAMPLES} samples"
echo "  Output: ${PHASED_PREFIX}.vcf.gz"

# =============================================================================
# STEP 5: Extract phased genotypes to TSV for R
# =============================================================================
# The R script reads a plain TSV rather than parsing the VCF directly,
# which avoids a VCF dependency in R and makes the data straightforward
# to inspect manually.
#
# Output format (tab-separated):
#   CHROM   POS       ID       REF  ALT  sample1  sample2  ...
#   chr17   39905000  rs12345  A    G    0|0      0|1      1|1  ...
#
# The '|' separator indicates phased haplotypes:
#   0|1 = ref allele on haplotype 1, alt on haplotype 2
#   1|0 = alt on haplotype 1, ref on haplotype 2
# These are distinct haplotype assignments treated differently by haplo.stats.
# %ID is included so the R script can cross-reference rsIDs from the annotated
# Excel report without requiring a positional join.
echo ""
echo "STEP 5: Extracting phased genotypes to TSV"
echo "----------------------------------------------------------------------"

PHASED_TSV="${HAPLO_DIR}/phased_genotypes.tsv"

# Save sample order from the phased VCF; must match the TSV column order
bcftools query -l "${PHASED_PREFIX}.vcf.gz" > "${HAPLO_DIR}/sample_order.txt"

# Extract one row per variant: positions, rsID, alleles, then all sample GTs
bcftools query \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
    "${PHASED_PREFIX}.vcf.gz" \
    > "${PHASED_TSV}"

# Prepend a header row; paste -s -d'\t' joins the sample names into one
# tab-separated line matching the column order set by bcftools query -l
HEADER="CHROM\tPOS\tID\tREF\tALT\t$(paste -s -d'\t' "${HAPLO_DIR}/sample_order.txt")"
{ echo -e "${HEADER}"; cat "${PHASED_TSV}"; } > "${PHASED_TSV}.tmp"
mv "${PHASED_TSV}.tmp" "${PHASED_TSV}"

N_ROWS=$(( $(wc -l < "${PHASED_TSV}") - 1 ))   # Subtract header line
echo "  Phased genotype table: ${N_ROWS} SNPs x ${N_SAMPLES} samples"
echo "  Saved to: ${PHASED_TSV}"

# =============================================================================
# STEP 6: Create sample metadata file for R
# =============================================================================
# Maps each sample to its cohort, tissue, and combined group label.
# Used by 15_haplo_stats.r to define case/control groups for association
# testing. Sample names must match the column headers in phased_genotypes.tsv.
echo ""
echo "STEP 6: Creating sample metadata file"
echo "----------------------------------------------------------------------"

META_TSV="${HAPLO_DIR}/sample_metadata.tsv"
echo -e "Sample\tCohort\tTissue\tGroup" > "${META_TSV}"

for cohort in breast endometrium; do
    for tissue in normal tumour; do
        tissue_label=$([ "${tissue}" = "normal" ] && echo "Healthy" || echo "Tumour")
        cohort_label="${cohort^}"   # Bash parameter expansion: capitalise first letter

        CALLS_DIR="${BASE}/${cohort}/${tissue}/dna_calls"
        [ ! -d "${CALLS_DIR}" ] && continue

        for sample_dir in "${CALLS_DIR}"/*/; do
            sample=$(basename "${sample_dir}")
            echo -e "${sample}\t${cohort_label}\t${tissue_label}\t${cohort_label}_${tissue_label}"
        done
    done
done >> "${META_TSV}"

N_META=$(( $(wc -l < "${META_TSV}") - 1 ))
echo "  Metadata written for ${N_META} samples"

# A mismatch here indicates samples were skipped in Step 1 or failed phasing
if [ "${N_META}" -ne "${N_SAMPLES}" ]; then
    echo "  Warning: metadata has ${N_META} samples but phased VCF has ${N_SAMPLES}."
    echo "    Check sample_list.txt and dna_calls/ directories for discrepancies."
fi

# --- Summary -----------------------------------------------------------------
echo ""
echo "======================================================================"
echo "PHASING PIPELINE COMPLETE"
echo "======================================================================"
echo ""
echo "Output files in: ${HAPLO_DIR}/"
echo "  cohort_merged.vcf.gz   — multi-sample merged VCF (all variants)"
echo "  cohort_snps.vcf.gz     — biallelic chr17 SNPs (pre-phasing input)"
echo "  cohort_phased.vcf.gz   — BEAGLE-phased haplotypes"
echo "  phased_genotypes.tsv   — phased GT matrix (SNPs x samples)"
echo "  sample_metadata.tsv    — sample group assignments for R"
echo "  sample_order.txt       — sample column order matching TSV"
echo ""
echo "  Reference panel used: ${USE_REF}"
echo "  Genetic map used:     ${USE_MAP}"
echo ""
echo "Note: gnomAD NFE AF > 1%% filtering is applied in 15_haplo_stats.r."
echo "      The phased VCF contains all biallelic chr17 SNPs."
echo ""
echo "Next step: Rscript 15_haplo_stats.r"
echo ""
