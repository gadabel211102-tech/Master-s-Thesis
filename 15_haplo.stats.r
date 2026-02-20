#!/usr/bin/env Rscript
# =============================================================================
# Script 15b: Haplotype Analysis with haplo.stats
# =============================================================================
#
# PURPOSE:
#   Uses BEAGLE-phased genotypes (from 15_haplotypes.sh) to:
#     1. Filter to established SNPs (gnomAD NFE AF > 1%) matching script 11
#     2. Estimate global haplotype frequencies using the EM algorithm
#     3. Test haplotype-phenotype associations across five group comparisons
#     4. Export results to Excel and generate publication-quality plots
#
# STUDY DESIGN NOTE:
#   Samples are unpaired (tumour and healthy tissue from different patients).
#   Association testing therefore compares haplotype frequencies between
#   groups at the population level — not within individuals. This is
#   appropriate and is the standard approach for population-level
#   haplotype association studies.
#
# STATISTICAL APPROACH:
#   - Haplotype frequencies estimated by EM algorithm (haplo.em)
#   - Per-haplotype association tested by Fisher's exact test (carrier status)
#   - Odds ratios with 95% CI estimated by logistic regression (dosage model)
#   - Multiple testing correction: Benjamini-Hochberg FDR within each comparison
#   - Minimum haplotype frequency threshold applied to exclude very rare
#     haplotypes that cannot be reliably estimated from sample sizes
#
# FIX (2026-02-20):
#   Previously, haplotypes in the per-comparison association tests were
#   renumbered independently (H1, H2, … in alphabetical/binary string order),
#   causing the labels in the forest plot and Fisher's test figures to be
#   completely inconsistent with the global H1–H9 labels in the frequency
#   plots. The fix maps every per-comparison haplotype string back to its
#   global ID (H1–H9 etc.) before labelling. Haplotypes that are present in
#   a comparison subset but were too rare to appear in the global EM are
#   labelled "Hrare_<binary_string>" so they remain interpretable.
#
# INPUTS (produced by 15_haplotypes.sh):
#   phased_genotypes.tsv   — phased GT matrix (SNPs x samples)
#   sample_metadata.tsv    — sample group assignments
#   GSDMB_Annotated_Report_Fixed.xlsx — for gnomAD NFE AF filtering
#
# OUTPUTS:
#   19_Haplotype_Results.xlsx  — full results (SNPs used, frequencies,
#                                 association tables, global test summary)
#   19_Global_Haplotype_Frequencies.png  — frequency bar chart
#   19_Haplotype_Freqs_By_Group.png      — stacked bar by group
#   19_Score_Test_Results.png            — p-value bar chart by comparison
#   19_OR_Forest_Plot.png                — forest plot (All Tumour vs Healthy)
#   19_Haplotype_Composition.png         — tile figure: REF/ALT base at each
#                                           SNP for every haplotype (rows =
#                                           haplotypes, columns = SNPs)
#
# USAGE:
#   Rscript 15_haplo_stats.r
#
# DEPENDENCIES:
#   haplo.stats, ggplot2, dplyr, tidyr, openxlsx
#   Install: install.packages(c("haplo.stats","ggplot2","dplyr","tidyr","openxlsx"))
#
# Author: Genomics Analysis Pipeline
# Date:   2026-02-20
# =============================================================================

suppressPackageStartupMessages({
  library(haplo.stats)   # EM haplotype frequency estimation and association tests
  library(ggplot2)       # Plotting
  library(dplyr)         # Data manipulation
  library(tidyr)         # Data reshaping
  library(openxlsx)      # Excel export
})

# =============================================================================
# CONFIGURATION
# =============================================================================

HAPLO_DIR  <- "/home/gadeaalonsoj/tfm/gsdmb_final_results/19_haplotype_phased"
GENO_FILE  <- file.path(HAPLO_DIR, "phased_genotypes.tsv")
META_FILE  <- file.path(HAPLO_DIR, "sample_metadata.tsv")
OUT_DIR    <- "/home/gadeaalonsoj/tfm/gsdmb_final_results/19_haplo_stats_results"
OUT_XLSX   <- file.path(OUT_DIR, "19_Haplotype_Results.xlsx")

# Annotated report used to extract established SNP positions (gnomAD NFE AF > 1%)
# This ensures the same SNP set as script 11 is used for haplotype analysis
ANNOTATED_REPORT <- "/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_Annotated_Report_Fixed.xlsx"

# Minimum haplotype frequency to report
# At ~138 samples (276 chromosomes), 2% corresponds to ~5-6 chromosomes
# Haplotypes below this threshold cannot be reliably estimated by EM
MIN_HAPLOTYPE_FREQ <- 0.02

# Minimum number of copies of a haplotype required in a comparison group
# to include it in association testing (prevents spurious results from
# haplotypes present in only 1-2 samples)
MIN_HAP_COPIES <- 3

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("======================================================================\n")
cat("SCRIPT 15b: HAPLOTYPE ANALYSIS WITH haplo.stats\n")
cat("======================================================================\n\n")
cat(sprintf("  Minimum haplotype frequency: %.0f%%\n", MIN_HAPLOTYPE_FREQ * 100))
cat(sprintf("  Minimum haplotype copies:    %d\n\n", MIN_HAP_COPIES))

# =============================================================================
# STEP 1: Load phased genotype data and sample metadata
# =============================================================================
cat("STEP 1: Loading phased genotype data\n")
cat("----------------------------------------------------------------------\n")

# check.names = FALSE preserves sample names that contain special characters
geno_raw <- read.table(GENO_FILE, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
meta     <- read.table(META_FILE, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)

cat(sprintf("✓ Genotype matrix loaded: %d SNPs × %d samples\n",
            nrow(geno_raw), ncol(geno_raw) - 4))
cat(sprintf("✓ Metadata loaded: %d samples\n", nrow(meta)))

# Verify input files are not empty
if (nrow(geno_raw) == 0) stop("Genotype file is empty. Check 15_haplotypes.sh completed successfully.")
if (nrow(meta) == 0)     stop("Metadata file is empty. Check 15_haplotypes.sh Step 6.")

# =============================================================================
# STEP 2: Filter to established SNPs (gnomAD NFE AF > 1%)
# =============================================================================
# The phased VCF contains all chr17 biallelic SNPs to maximise LD information
# available to BEAGLE. Here we restrict to the established SNP set (gnomAD
# NFE AF > 1%) that was used in script 11, so that haplotype analysis is
# based on the same biologically meaningful variants.
#
# Positions are matched between the phased TSV and the annotated Excel report.
# Any established SNPs not found in the phased file are reported as warnings
# (they may have been filtered by BEAGLE or bcftools due to missing data).
# =============================================================================
cat("\nSTEP 2: Filtering to established SNPs (gnomAD NFE AF > 1%)\n")
cat("----------------------------------------------------------------------\n")

if (!file.exists(ANNOTATED_REPORT)) {
  stop(sprintf("Annotated report not found:\n  %s\n  Run script 07 first.",
               ANNOTATED_REPORT))
}

cat(sprintf("  Reading: %s\n", basename(ANNOTATED_REPORT)))
annot <- read.xlsx(ANNOTATED_REPORT, sheet = "Biological_Annotations")

# Case-insensitive column finder — robust to capitalisation differences
find_col_ci <- function(df, pattern) {
  grep(pattern, colnames(df), ignore.case = TRUE, value = TRUE)[1]
}

pos_col <- find_col_ci(annot, "^pos$|^position$")
nfe_col <- find_col_ci(annot, "gnomad.*nfe.*af|nfe.*af")

if (is.na(pos_col)) {
  stop("Cannot find POS column in annotated report.\n",
       "Available columns: ", paste(colnames(annot), collapse = ", "))
}
if (is.na(nfe_col)) {
  stop("Cannot find gnomAD NFE AF column in annotated report.\n",
       "Available columns: ", paste(colnames(annot), collapse = ", "))
}

cat(sprintf("  POS column:           '%s'\n", pos_col))
cat(sprintf("  gnomAD NFE AF column: '%s'\n", nfe_col))

# Apply gnomAD NFE > 1% filter (same as script 11)
# suppressWarnings handles non-numeric values such as '.' for missing entries
annot[[nfe_col]]  <- suppressWarnings(as.numeric(annot[[nfe_col]]))
snps_established  <- annot[!is.na(annot[[nfe_col]]) & annot[[nfe_col]] > 0.01, ]
established_positions <- unique(as.integer(snps_established[[pos_col]]))
established_positions <- established_positions[!is.na(established_positions)]

cat(sprintf("  Established SNPs (gnomAD NFE > 1%%): %d positions\n",
            length(established_positions)))

# Filter phased matrix to established SNP positions
keep     <- geno_raw$POS %in% established_positions
n_before <- nrow(geno_raw)
geno_raw <- geno_raw[keep, ]
n_after  <- nrow(geno_raw)

cat(sprintf("  SNPs in phased file before filter: %d\n", n_before))
cat(sprintf("  SNPs remaining after filter:       %d\n", n_after))

if (n_after == 0) {
  stop(
    "No SNPs remain after filtering.\n",
    "Check that POS values in phased_genotypes.tsv match those in\n",
    "the annotated report (gnomAD NFE AF > 0.01 variants).\n",
    "First positions in phased file:    ",
    paste(head(as.integer(geno_raw$POS)), collapse = ", "), "\n",
    "First established positions found: ",
    paste(head(established_positions), collapse = ", ")
  )
}

# Report established SNPs not found in the phased file
# (may have been filtered by BEAGLE or bcftools due to missing data)
missing_from_phased <- established_positions[!established_positions %in% geno_raw$POS]
if (length(missing_from_phased) > 0) {
  cat(sprintf("  ⚠ %d established SNP(s) absent from phased file:\n",
              length(missing_from_phased)))
  cat(sprintf("    %s\n", paste(missing_from_phased, collapse = ", ")))
  cat("    These were likely filtered by BEAGLE (low call rate) or bcftools.\n")
}

cat(sprintf("✓ Proceeding with %d SNPs\n", n_after))

# =============================================================================
# STEP 3: Build SNP info and parse phased genotypes into allele matrices
# =============================================================================
cat("\nSTEP 3: Parsing phased genotypes into allele matrices\n")
cat("----------------------------------------------------------------------\n")

# SNP labels in format "POS_REF>ALT" for readable output
snp_info   <- geno_raw[, c("CHROM", "POS", "REF", "ALT")]
snp_labels <- paste0(geno_raw$POS, "_", geno_raw$REF, ">", geno_raw$ALT)

# Sample columns: everything after the first 4 metadata columns
sample_cols <- colnames(geno_raw)[-(1:4)]
n_samples   <- length(sample_cols)
n_snps      <- nrow(geno_raw)

# Parse phased genotype strings (e.g. "0|1") into separate allele vectors
# BEAGLE always outputs phased '|' genotypes; unphased '/' may appear
# if a site was impossible to phase (rare in practice after BEAGLE)
parse_phased_gt <- function(gt_vec) {
  allele1 <- sapply(gt_vec, function(g) {
    if (is.na(g) || g %in% c(".", ".|.", "./.", ".")) return(NA_integer_)
    parts <- strsplit(g, "[|/]")[[1]]
    suppressWarnings(as.integer(parts[1]))
  })
  allele2 <- sapply(gt_vec, function(g) {
    if (is.na(g) || g %in% c(".", ".|.", "./.", ".")) return(NA_integer_)
    parts <- strsplit(g, "[|/]")[[1]]
    suppressWarnings(as.integer(parts[2]))
  })
  list(a1 = allele1, a2 = allele2)
}

# Initialise allele matrices: rows = samples, columns = SNPs
allele1_mat <- matrix(NA_integer_, nrow = n_samples, ncol = n_snps)
allele2_mat <- matrix(NA_integer_, nrow = n_samples, ncol = n_snps)
rownames(allele1_mat) <- sample_cols
rownames(allele2_mat) <- sample_cols
colnames(allele1_mat) <- snp_labels
colnames(allele2_mat) <- snp_labels

# Fill matrices: for each SNP, parse all sample genotypes
for (j in seq_len(n_snps)) {
  gt_vec <- as.character(geno_raw[j, sample_cols])
  parsed <- parse_phased_gt(gt_vec)
  allele1_mat[, j] <- parsed$a1
  allele2_mat[, j] <- parsed$a2
}

# haplo.stats uses 1-based allele coding: 1 = reference, 2 = alternate
# Our matrices use 0-based coding (0 = ref, 1 = alt), so add 1
# NA values are preserved through arithmetic
allele1_mat_hap <- allele1_mat + 1L
allele2_mat_hap <- allele2_mat + 1L

cat(sprintf("✓ Allele matrices built: %d samples × %d SNPs\n", n_samples, n_snps))

# Report missing data rate
n_missing <- sum(is.na(allele1_mat))
total     <- n_samples * n_snps
cat(sprintf("  Missing genotypes: %d / %d (%.1f%%)\n",
            n_missing, total, 100 * n_missing / total))

# Align metadata to the sample order in the genotype matrix
# Samples in the metadata but not in the genotype matrix receive NA rows
meta_aligned <- meta[match(sample_cols, meta$Sample), ]
missing_meta <- sum(is.na(meta_aligned$Sample))
if (missing_meta > 0) {
  cat(sprintf("  ⚠ %d samples in genotype matrix have no metadata — excluded\n",
              missing_meta))
}

# Keep only samples with valid metadata
valid_samples <- !is.na(meta_aligned$Sample)
meta_aligned  <- meta_aligned[valid_samples, ]
a1 <- allele1_mat_hap[valid_samples, , drop = FALSE]
a2 <- allele2_mat_hap[valid_samples, , drop = FALSE]

cat(sprintf("✓ %d samples retained after metadata alignment\n", nrow(meta_aligned)))
cat(sprintf("  Cohorts: %s\n", paste(sort(unique(meta_aligned$Cohort)), collapse = ", ")))
cat(sprintf("  Tissues: %s\n", paste(sort(unique(meta_aligned$Tissue)), collapse = ", ")))

# =============================================================================
# STEP 4: Build genotype matrix for haplo.stats
# =============================================================================
# haplo.em() requires a data frame with alternating allele1/allele2 columns
# per SNP: (snp1_a1, snp1_a2, snp2_a1, snp2_a2, ...). This is constructed
# from the two allele matrices.
# =============================================================================

# Helper: interleave allele1 and allele2 matrices into haplo.stats format
build_geno_matrix <- function(a1_mat, a2_mat) {
  n    <- ncol(a1_mat)
  cols <- vector("list", n * 2)
  for (j in seq_len(n)) {
    cols[[2 * j - 1]] <- a1_mat[, j]   # allele 1 for SNP j
    cols[[2 * j]]     <- a2_mat[, j]   # allele 2 for SNP j
  }
  as.data.frame(cols)
}

# =============================================================================
# STEP 5: Estimate global haplotype frequencies (full cohort, EM algorithm)
# =============================================================================
cat("\nSTEP 4: Estimating global haplotype frequencies (EM algorithm)\n")
cat("----------------------------------------------------------------------\n")
cat("  Fitting haplo.em across all samples...\n")

geno_mat_full <- build_geno_matrix(a1, a2)

# haplo.em fits haplotype frequencies by the Expectation-Maximisation
# algorithm across the full cohort. miss.val = c(0, NA) treats both
# 0-coded alleles (reference in 0-based coding) and NA as missing —
# note we have already converted to 1-based coding so 0 should not appear,
# but it is retained as a safety net.
# min.posterior = 0.001 suppresses haplotypes with very low posterior probability
hap_em_global <- haplo.em(
  geno        = geno_mat_full,
  locus.label = snp_labels,
  miss.val    = c(0, NA),
  control     = haplo.em.control(min.posterior = 0.001)
)

# Build global frequency table
global_freqs            <- as.data.frame(hap_em_global$haplotype)
global_freqs$Frequency  <- hap_em_global$hap.prob
# Estimated count = frequency × total chromosomes (2 per diploid sample)
global_freqs$Count      <- round(hap_em_global$hap.prob * nrow(a1) * 2)
global_freqs            <- global_freqs[order(-global_freqs$Frequency), ]

# Apply minimum frequency filter
global_freqs <- global_freqs[global_freqs$Frequency >= MIN_HAPLOTYPE_FREQ, ]
global_freqs$Haplotype_ID <- paste0("H", seq_len(nrow(global_freqs)))

cat(sprintf("✓ %d haplotypes estimated (freq ≥ %.0f%%)\n",
            nrow(global_freqs), MIN_HAPLOTYPE_FREQ * 100))
cat("\nGlobal haplotype frequencies:\n")
print(global_freqs[, c("Haplotype_ID", "Frequency", "Count")], row.names = FALSE)

# =============================================================================
# Build global haplotype string lookup (used by association testing below)
# =============================================================================
# Convert the 1-based haplo.stats allele coding back to 0/1 binary strings
# so we can match per-comparison haplotype strings to global IDs.
#
# haplo.stats codes: 1 = REF, 2 = ALT  →  subtract 1 → 0 = REF, 1 = ALT
# paste(collapse="") produces e.g. "001010..." for a 15-SNP haplotype
#
# NOTE: global_freqs SNP columns may be factors after as.data.frame() +
# reordering, so we coerce explicitly to a plain numeric matrix first.
# =============================================================================
global_snp_matrix <- matrix(
  as.numeric(as.matrix(global_freqs[, snp_labels, drop = FALSE])),
  nrow = nrow(global_freqs),
  ncol = length(snp_labels)
)
global_hap_strings <- apply(
  global_snp_matrix - 1L,
  1,
  paste, collapse = ""
)
# Named vector: binary string → global haplotype ID (e.g. "001..." → "H3")
string_to_global_id <- setNames(global_freqs$Haplotype_ID, global_hap_strings)

cat(sprintf("\n✓ Global haplotype string lookup built (%d entries)\n",
            length(string_to_global_id)))

# =============================================================================
# STEP 6: Association testing — five group comparisons
# =============================================================================
# For each comparison we:
#   1. Subset samples to the two groups being compared
#   2. Convert haplotypes to binary carrier status per sample
#   3. Run Fisher's exact test (carrier vs non-carrier)
#   4. Estimate OR with 95% CI using logistic regression (dosage coding)
#   5. Apply BH FDR correction within each comparison
#
# FIX: Haplotypes are now mapped back to their global IDs (H1, H2, …) using
# string_to_global_id before being labelled. This ensures consistent labelling
# across all figures. Haplotypes not in the global set are labelled "Hrare".
#
# Comparisons:
#   a. Breast Tumour vs Breast Healthy        — within-cohort, breast
#   b. Endometrium Tumour vs Endometrium Healthy — within-cohort, endometrium
#   c. Breast vs Endometrium (Tumour only)    — cross-cancer, tumour
#   d. Breast vs Endometrium (Healthy only)   — cross-cancer, healthy tissue
#   e. All Tumour vs All Healthy              — combined cross-cohort
# =============================================================================
cat("\n\nSTEP 5: Haplotype association testing\n")
cat("----------------------------------------------------------------------\n")

comparisons <- list(
  list(name        = "Breast: Tumour vs Healthy",
       filter      = meta_aligned$Cohort == "Breast",
       outcome_col = "Tissue",
       case        = "Tumour",
       control     = "Healthy"),

  list(name        = "Endometrium: Tumour vs Healthy",
       filter      = meta_aligned$Cohort == "Endometrium",
       outcome_col = "Tissue",
       case        = "Tumour",
       control     = "Healthy"),

  list(name        = "Tumour: Breast vs Endometrium",
       filter      = meta_aligned$Tissue == "Tumour",
       outcome_col = "Cohort",
       case        = "Breast",
       control     = "Endometrium"),

  list(name        = "Healthy: Breast vs Endometrium",
       filter      = meta_aligned$Tissue == "Healthy",
       outcome_col = "Cohort",
       case        = "Breast",
       control     = "Endometrium"),

  list(name        = "All: Tumour vs Healthy",
       filter      = rep(TRUE, nrow(meta_aligned)),
       outcome_col = "Tissue",
       case        = "Tumour",
       control     = "Healthy")
)

cat(sprintf("  Cohort values in metadata: %s\n",
            paste(sort(unique(meta_aligned$Cohort)), collapse = ", ")))
cat(sprintf("  Tissue values in metadata: %s\n",
            paste(sort(unique(meta_aligned$Tissue)), collapse = ", ")))

# Helper: convert allele matrix (1-based) to haplotype strings (0/1 per SNP)
# e.g. a sample with ref at all SNPs = "00000"; alt at first SNP = "10000"
make_hap_strings <- function(a_mat) {
  apply(a_mat - 1L, 1, paste, collapse = "")
}

# Main association function
# Returns a data frame of per-haplotype results, or NULL if too few samples.
#
# FIXED: now accepts string_to_global_id to map haplotype strings back to the
# globally consistent H1, H2, … labels rather than re-numbering locally.
run_haplo_assoc <- function(a1_sub, a2_sub, y, snp_labels, comp_name,
                            case_label, control_label, n_case, n_control,
                            string_to_global_id) {

  # Build haplotype strings for each chromosome per sample
  hap1_strings <- make_hap_strings(a1_sub)
  hap2_strings <- make_hap_strings(a2_sub)

  # Count observed copies of each haplotype across both chromosomes
  all_haps    <- table(c(hap1_strings, hap2_strings))

  # Only test haplotypes with sufficient copies to be reliable
  common_haps <- names(all_haps[all_haps >= MIN_HAP_COPIES])
  n_haps      <- length(common_haps)
  cat(sprintf("    Haplotypes with >= %d copies: %d\n", MIN_HAP_COPIES, n_haps))

  if (n_haps < 2) {
    cat("    WARNING: Too few haplotypes — skipping comparison\n")
    return(NULL)
  }

  results <- list()

  for (h in seq_along(common_haps)) {
    hap <- common_haps[h]

    # --- Map haplotype string to its global ID ---
    # If this haplotype string matches one of the globally estimated haplotypes,
    # use that label (e.g. "H1", "H3"). Otherwise assign the next sequential
    # number continuing from where the global haplotypes left off, stored in
    # rare_hap_registry so the same rare haplotype gets the same label if it
    # appears in multiple comparisons.
    if (hap %in% names(string_to_global_id)) {
      hap_label <- string_to_global_id[hap]
    } else if (hap %in% names(rare_hap_registry)) {
      hap_label <- rare_hap_registry[[hap]]
    } else {
      next_id             <- paste0("H", length(string_to_global_id) +
                                         length(rare_hap_registry) + 1)
      rare_hap_registry[[hap]] <<- next_id   # <<- updates parent env
      hap_label           <- next_id
      cat(sprintf("    ⚠ Rare haplotype assigned label '%s'\n", hap_label))
    }

    # Dosage = number of chromosomes carrying this haplotype (0, 1, or 2)
    dosage <- as.integer(hap1_strings == hap) +
              as.integer(hap2_strings == hap)

    # Skip if dosage has no variance (all samples have same value — uninformative)
    if (var(dosage) == 0) next

    # Carrier status: 1 if the sample carries at least one copy, 0 otherwise
    carrier <- as.integer(dosage >= 1)

    # Haplotype frequency in cases and controls (per chromosome)
    freq_case    <- mean(dosage[y == 1]) / 2
    freq_control <- mean(dosage[y == 0]) / 2

    # Fisher's exact test on 2x2 carrier/non-carrier × case/control table
    tbl  <- table(carrier, y)
    if (nrow(tbl) < 2) next
    fish <- fisher.test(tbl)

    # Logistic regression for OR with 95% CI (dosage model: 0, 1, or 2 copies)
    # More powerful than carrier-only model as it uses full dosage information
    fit <- tryCatch(
      suppressWarnings(glm(y ~ dosage, family = binomial())),
      error = function(e) NULL
    )

    if (!is.null(fit) && nrow(summary(fit)$coefficients) >= 2) {
      beta  <- coef(fit)["dosage"]
      se    <- summary(fit)$coefficients["dosage", "Std. Error"]
      or    <- exp(beta)
      or_lo <- exp(beta - 1.96 * se)
      or_hi <- exp(beta + 1.96 * se)
      p_glm <- summary(fit)$coefficients["dosage", "Pr(>|z|)"]
    } else {
      # Logistic regression failed to converge (e.g. perfect separation)
      # Fall back to Fisher OR; CI reported as NA — do NOT impute values
      or    <- as.numeric(fish$estimate)
      or_lo <- NA_real_
      or_hi <- NA_real_
      p_glm <- fish$p.value
      cat(sprintf("    ⚠ Logistic regression failed for %s — OR from Fisher, CI = NA\n",
                  hap_label))
    }

    # Build readable haplotype definition showing alt-allele positions
    snp_pos  <- sapply(strsplit(snp_labels, "_"), `[`, 1)
    hap_bits <- strsplit(hap, "")[[1]]
    alt_pos  <- snp_pos[hap_bits == "1"]
    hap_def  <- if (length(alt_pos) == 0) {
      paste0(hap_label, "(ref)")      # All-reference haplotype
    } else {
      paste0(hap_label, "(", paste(alt_pos, collapse = "/"), ")")
    }

    results[[h]] <- data.frame(
      Haplotype      = hap_label,          # FIX: globally consistent ID
      Hap_Definition = hap_def,
      Global_Freq    = round(sum(all_haps[hap]) / sum(all_haps), 4),
      Freq_Case      = round(freq_case, 4),
      Freq_Control   = round(freq_control, 4),
      OR             = round(or, 4),
      OR_95CI_Lo     = round(or_lo, 4),   # NA if logistic regression failed
      OR_95CI_Hi     = round(or_hi, 4),   # NA if logistic regression failed
      P_Fisher       = fish$p.value,
      P_Logistic     = p_glm,
      Comparison     = comp_name,
      Case           = case_label,
      Control        = control_label,
      N_Case         = n_case,
      N_Control      = n_control,
      stringsAsFactors = FALSE
    )
  }

  if (length(results) == 0) return(NULL)

  res_df <- do.call(rbind, results)

  # BH FDR correction applied within this comparison only
  res_df$FDR                <- p.adjust(res_df$P_Fisher, method = "BH")
  res_df$Significant_Fisher <- res_df$P_Fisher < 0.05
  res_df$Significant_FDR    <- res_df$FDR < 0.05
  res_df                    <- res_df[order(res_df$P_Fisher), ]

  # Best individual haplotype p-value (note: this is NOT an omnibus test —
  # it is the minimum Fisher p across haplotypes within this comparison)
  best_p <- min(res_df$P_Fisher)
  cat(sprintf("    Best individual haplotype p (Fisher): %.4f | FDR: %.4f\n",
              best_p, min(res_df$FDR)))

  list(table = res_df, best_haplotype_p = best_p)
}

# Registry to assign consistent sequential IDs to rare haplotypes across
# all comparisons (e.g. H10, H11, ...) continuing from the global H1-H9 set
rare_hap_registry <- list()

all_results   <- list()
score_results <- list()

for (comp in comparisons) {
  cat(sprintf("\n  Running: %s\n", comp$name))

  # Subset to samples relevant to this comparison
  idx <- which(comp$filter &
               meta_aligned[[comp$outcome_col]] %in% c(comp$case, comp$control))

  meta_sub  <- meta_aligned[idx, ]
  a1_sub    <- a1[idx, , drop = FALSE]
  a2_sub    <- a2[idx, , drop = FALSE]
  n_case    <- sum(meta_sub[[comp$outcome_col]] == comp$case)
  n_control <- sum(meta_sub[[comp$outcome_col]] == comp$control)

  cat(sprintf("    %s (case): %d | %s (control): %d\n",
              comp$case, n_case, comp$control, n_control))

  # Require at least 3 samples per group
  if (n_case < 3 || n_control < 3) {
    cat("    WARNING: Too few samples — skipping\n")
    next
  }

  # Binary outcome vector: 1 = case, 0 = control
  y <- as.integer(meta_sub[[comp$outcome_col]] == comp$case)

  # Pass string_to_global_id so haplotypes get globally consistent labels
  res <- run_haplo_assoc(a1_sub, a2_sub, y, snp_labels,
                         comp$name, comp$case, comp$control,
                         n_case, n_control,
                         string_to_global_id)

  if (!is.null(res)) {
    score_results[[comp$name]] <- list(
      table            = res$table,
      best_haplotype_p = res$best_haplotype_p
    )
    all_results[[comp$name]] <- res$table
    cat(sprintf("    Saved %d haplotype results\n", nrow(res$table)))
  }
}

# =============================================================================
# STEP 7: Save results to Excel
# =============================================================================
cat("\n\nSTEP 6: Saving results to Excel\n")
cat("----------------------------------------------------------------------\n")

wb <- createWorkbook()

# Sheet 1: Which SNPs were used — provides full analytical provenance
snp_used_df <- data.frame(
  POS       = geno_raw$POS,
  REF       = geno_raw$REF,
  ALT       = geno_raw$ALT,
  SNP_Label = snp_labels
)
addWorksheet(wb, "SNPs_Used")
writeData(wb, "SNPs_Used", snp_used_df)

# Sheet 2: Global haplotype frequencies (full cohort EM estimate)
addWorksheet(wb, "Global_Haplotype_Freqs")
writeData(wb, "Global_Haplotype_Freqs", global_freqs)

# Sheets 3+: Per-comparison association results
for (comp_name in names(all_results)) {
  # Excel sheet names cannot contain colons, capped at 31 characters
  sheet_name <- substr(gsub(":", " -", comp_name), 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, all_results[[comp_name]])
}

# Final sheet: Summary of best p-values per comparison
# NOTE: best_haplotype_p is the minimum Fisher p across haplotypes within
# each comparison. It is NOT an omnibus global test p-value. An omnibus
# test (e.g. haplo.score) could be added in future work.
summary_df <- data.frame(
  Comparison           = names(score_results),
  Best_Haplotype_P     = sapply(score_results, function(x) x$best_haplotype_p),
  Significant_Fisher   = sapply(score_results, function(x) any(x$table$Significant_Fisher)),
  Significant_FDR      = sapply(score_results, function(x) any(x$table$Significant_FDR))
)
addWorksheet(wb, "Summary_Best_P_Values")
writeData(wb, "Summary_Best_P_Values", summary_df)

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("✓ Results saved to: %s\n", OUT_XLSX))

# =============================================================================
# STEP 8: Plots
# =============================================================================
cat("\nSTEP 7: Generating plots\n")
cat("----------------------------------------------------------------------\n")

# --- Plot 1: Global haplotype frequency bar chart ---
p1 <- ggplot(global_freqs,
             aes(x = reorder(Haplotype_ID, -Frequency),
                 y = Frequency * 100)) +
  geom_col(fill = "#3498db", colour = "black", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Frequency * 100)),
            vjust = -0.4, size = 3.5) +
  labs(
    title    = "Global Haplotype Frequencies",
    subtitle = sprintf("EM estimation across %d samples ",
                       nrow(a1), n_snps),
    x = "Haplotype",
    y = "Frequency (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "19_Global_Haplotype_Frequencies.png"),
       p1, width = 8, height = 5, dpi = 300)
cat("  ✓ Global frequency plot saved\n")

# --- Plot 2: Haplotype frequencies by group (stacked bar) ---
# Runs EM separately within each group to get group-specific frequencies.
# FIX: Group-level haplotypes are also mapped back to global IDs so that the
# same haplotype has the same colour and label across all four group bars.
group_freq_list <- list()
for (grp in unique(meta_aligned$Group)) {
  idx <- which(meta_aligned$Group == grp)
  if (length(idx) < 3) next    # Skip groups with fewer than 3 samples
  a1_g   <- a1[idx, , drop = FALSE]
  a2_g   <- a2[idx, , drop = FALSE]
  geno_g <- build_geno_matrix(a1_g, a2_g)
  tryCatch({
    em_g <- haplo.em(geno = geno_g, locus.label = snp_labels,
                     miss.val = c(0, NA),
                     control  = haplo.em.control(min.posterior = 0.001))
    hap_df           <- as.data.frame(em_g$haplotype)
    hap_df$Frequency <- em_g$hap.prob
    hap_df$Group     <- grp

    # Map each group-level haplotype string to its global ID
    # Coerce to numeric matrix first — columns may be factors after as.data.frame()
    grp_snp_matrix <- matrix(
      as.numeric(as.matrix(hap_df[, snp_labels, drop = FALSE])),
      nrow = nrow(hap_df),
      ncol = length(snp_labels)
    )
    grp_hap_strings <- apply(grp_snp_matrix - 1L, 1, paste, collapse = "")
    hap_df$Haplotype_ID <- sapply(grp_hap_strings, function(s) {
      if (s %in% names(string_to_global_id)) {
        string_to_global_id[s]
      } else if (s %in% names(rare_hap_registry)) {
        rare_hap_registry[[s]]
      } else {
        next_id <- paste0("H", length(string_to_global_id) +
                                length(rare_hap_registry) + 1)
        rare_hap_registry[[s]] <<- next_id
        next_id
      }
    })

    group_freq_list[[grp]] <- hap_df
  }, error = function(e) {
    cat(sprintf("  ⚠ EM failed for group %s: %s\n", grp, conditionMessage(e)))
    NULL
  })
}

if (length(group_freq_list) > 0) {
  group_freq_df <- bind_rows(lapply(names(group_freq_list), function(grp) {
    df <- group_freq_list[[grp]]
    df <- df[df$Frequency >= MIN_HAPLOTYPE_FREQ, c("Group", "Haplotype_ID", "Frequency")]

    # Collapse any haplotype not in the global H1-H9 set into "Other"
    # These are subgroup-specific EM noise haplotypes not reliably estimated
    global_ids <- global_freqs$Haplotype_ID
    df$Haplotype_ID <- ifelse(df$Haplotype_ID %in% global_ids,
                              df$Haplotype_ID, "Other")

    # Sum "Other" frequencies together per group
    df <- df %>%
      group_by(Group, Haplotype_ID) %>%
      summarise(Frequency = sum(Frequency), .groups = "drop")
    df
  }))

  # Order levels: H1-H9 by global frequency, then Other last
  hap_level_order <- c(global_freqs$Haplotype_ID, "Other")
  group_freq_df$Haplotype_ID <- factor(group_freq_df$Haplotype_ID,
                                        levels = hap_level_order)

  # Colour palette: Set1 for H1-H9, grey for Other
  n_global  <- nrow(global_freqs)
  pal_cols  <- RColorBrewer::brewer.pal(max(n_global, 3), "Set1")[seq_len(n_global)]
  fill_cols <- c(setNames(pal_cols, global_freqs$Haplotype_ID), "Other" = "#cccccc")

  p2 <- ggplot(group_freq_df,
               aes(x = Group, y = Frequency * 100, fill = Haplotype_ID)) +
    geom_col(colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    labs(
      title = "Haplotype Frequencies by Sample Group",
      subtitle = "Haplotypes below global 2% frequency threshold grouped as 'Other'",
      x     = NULL,
      y     = "Frequency (%)",
      fill  = "Haplotype"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(face = "bold")) +
    ylim(0, 105)

  ggsave(file.path(OUT_DIR, "19_Haplotype_Freqs_By_Group.png"),
         p2, width = 10, height = 6, dpi = 300)
  cat("  ✓ Group frequency plot saved\n")
}

# --- Plot 3: Association results — p-value bar chart ---
if (length(score_results) > 0) {
  score_plot_data <- bind_rows(lapply(names(score_results), function(cn) {
    tbl    <- score_results[[cn]]$table
    best_p <- score_results[[cn]]$best_haplotype_p
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    data.frame(
      Comparison   = cn,
      Best_P       = best_p,
      Haplotype    = tbl$Haplotype,
      OR           = tbl$OR,
      OR_Lo        = tbl$OR_95CI_Lo,
      OR_Hi        = tbl$OR_95CI_Hi,
      P_Value      = tbl$P_Fisher,
      FDR          = tbl$FDR,
      Freq_Case    = tbl$Freq_Case,
      Freq_Control = tbl$Freq_Control,
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(score_plot_data) && nrow(score_plot_data) > 0) {

    # Order haplotype factor levels globally so x-axis is consistent across facets
    rare_ids         <- unlist(rare_hap_registry)
    hap_levels_assoc <- c(global_freqs$Haplotype_ID, rare_ids)
    score_plot_data$Haplotype <- factor(score_plot_data$Haplotype,
                                         levels = hap_levels_assoc)
    score_plot_data$Comparison <- factor(score_plot_data$Comparison,
                                          levels = names(score_results))

    # Three-level significance category
    score_plot_data$Sig_Category <- ifelse(
      score_plot_data$FDR < 0.05, "FDR < 0.05",
      ifelse(score_plot_data$P_Value < 0.05, "Fisher p < 0.05", "Not significant")
    )
    score_plot_data$Sig_Category <- factor(
      score_plot_data$Sig_Category,
      levels = c("Not significant", "Fisher p < 0.05", "FDR < 0.05")
    )

    sig_colours <- c(
      "Not significant"  = "#95a5a6",
      "Fisher p < 0.05"  = "#f39c12",
      "FDR < 0.05"       = "#e74c3c"
    )

    # Bar chart: only haplotypes with sufficient copies shown per facet
    p3a <- ggplot(score_plot_data,
                  aes(x = Haplotype,
                      y = -log10(P_Value),
                      fill = Sig_Category)) +
      geom_col(colour = "black", alpha = 0.85) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                 colour = "red", alpha = 0.7) +
      scale_fill_manual(
        values = sig_colours,
        name   = "Significance",
        drop   = FALSE
      ) +
      facet_wrap(~ Comparison, scales = "free_x", ncol = 2) +
      labs(
        title    = "Haplotype Association: Fisher's Exact Test",
        subtitle = "Red dashed line = nominal p < 0.05  |  Orange = Fisher p < 0.05  |  Red = FDR < 0.05",
        x        = "Haplotype",
        y        = "-log10(p-value)"
      ) +
      theme_bw(base_size = 11) +
      theme(strip.text  = element_text(face = "bold", size = 7),
            plot.title  = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(OUT_DIR, "19_Score_Test_Results.png"),
           p3a, width = 14, height = 10, dpi = 300)
    cat("  ✓ Association p-value plot saved\n")
  }
}

# --- Plot 4: Haplotype composition tile figure ---
# Shows the actual allele (REF/ALT base) at each SNP for every global haplotype.
# Only haplotypes in global_freqs (freq >= MIN_HAPLOTYPE_FREQ) are displayed.
# Rows = haplotypes (ordered by descending frequency), columns = SNPs.
# Red tiles = alternative allele; blue tiles = reference allele.
cat("  Generating haplotype composition tile figure...\n")

tryCatch({

  # Build a long-format data frame: one row per haplotype × SNP combination
  hap_snp_rows <- list()
  for (i in seq_len(nrow(global_freqs))) {
    hap_id   <- global_freqs$Haplotype_ID[i]
    hap_freq <- global_freqs$Frequency[i]
    hap_count <- global_freqs$Count[i]

    for (j in seq_len(n_snps)) {
      # haplo.stats codes: 1 = REF, 2 = ALT
      allele_code <- global_freqs[i, snp_labels[j]]
      is_alt      <- !is.na(allele_code) && allele_code == 2
      base        <- if (is_alt) snp_info$ALT[j] else snp_info$REF[j]

      hap_snp_rows[[length(hap_snp_rows) + 1]] <- data.frame(
        Haplotype  = hap_id,
        Frequency  = hap_freq,
        Count      = hap_count,
        SNP_Label  = snp_labels[j],
        POS        = snp_info$POS[j],
        REF        = snp_info$REF[j],
        ALT        = snp_info$ALT[j],
        Allele     = as.integer(allele_code),
        Base       = base,
        Is_Alt     = is_alt,
        stringsAsFactors = FALSE
      )
    }
  }

  tile_df <- do.call(rbind, hap_snp_rows)

  # Order haplotypes by descending frequency (H1 at top)
  hap_order       <- global_freqs$Haplotype_ID   # already sorted by -Frequency
  tile_df$Haplotype <- factor(tile_df$Haplotype, levels = rev(hap_order))

  # SNP column order: left to right by genomic position
  snp_order <- snp_labels[order(snp_info$POS)]
  tile_df$SNP_Label <- factor(tile_df$SNP_Label, levels = snp_order)

  # Y-axis label: "H1  24.8% (n=68)"
  hap_label_map <- setNames(
    sprintf("%s  %.1f%% (n=%d)",
            global_freqs$Haplotype_ID,
            global_freqs$Frequency * 100,
            global_freqs$Count),
    global_freqs$Haplotype_ID
  )
  tile_df$Hap_Label <- factor(
    hap_label_map[as.character(tile_df$Haplotype)],
    levels = rev(hap_label_map[hap_order])
  )

  # ------------------------------------------------------------------
  # Build POS -> Gene lookup from the annotated report (SYMBOL column)
  # ------------------------------------------------------------------
  sym_col       <- find_col_ci(annot, "^symbol$|^gene$|^gene_symbol$")
  pos_col_annot <- find_col_ci(annot, "^pos$|^position$")

  snp_gene_map <- setNames(rep("Unknown", n_snps), as.character(snp_info$POS))

  if (!is.na(sym_col) && !is.na(pos_col_annot)) {
    annot_genes <- annot[, c(pos_col_annot, sym_col), drop = FALSE]
    annot_genes[[pos_col_annot]] <- as.integer(annot_genes[[pos_col_annot]])
    annot_genes <- annot_genes[!is.na(annot_genes[[pos_col_annot]]), ]
    for (pos_val in unique(snp_info$POS)) {
      hits <- annot_genes[annot_genes[[pos_col_annot]] == pos_val, sym_col]
      hits <- hits[!is.na(hits) & nchar(trimws(hits)) > 0]
      if (length(hits) > 0)
        snp_gene_map[as.character(pos_val)] <- trimws(hits[1])
    }
  } else {
    cat("  ⚠ SYMBOL column not found in annotated report — gene labels will show 'Unknown'\n")
  }

  # X-axis label: three-line "POS\nREF>ALT\nGENE"
  snp_xlab <- setNames(
    paste0(snp_info$POS, "\n",
           snp_info$REF, ">", snp_info$ALT, "\n",
           snp_gene_map[as.character(snp_info$POS)]),
    snp_labels
  )
  tile_df$SNP_XLabel <- factor(
    snp_xlab[as.character(tile_df$SNP_Label)],
    levels = snp_xlab[snp_order]
  )

  # Colour-code the x-axis tick labels by gene for quick visual grouping
  gene_colours <- c(
    "GSDMB"  = "#1A5276",
    "GSDMA"  = "#6C3483",
    "ORMDL3" = "#1E8449",
    "Unknown"= "#7F8C8D"
  )
  axis_tick_cols <- sapply(levels(tile_df$SNP_XLabel), function(lbl) {
    gene <- strsplit(lbl, "\n")[[1]][3]
    col  <- gene_colours[gene]
    if (is.na(col)) "#7F8C8D" else col
  })

  p4 <- ggplot(tile_df,
               aes(x     = SNP_XLabel,
                   y     = Hap_Label,
                   fill  = Is_Alt,
                   label = Base)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(size = 3.8, fontface = "bold",
              colour = ifelse(tile_df$Is_Alt[order(tile_df$SNP_XLabel, tile_df$Hap_Label)],
                              "white", "#1a1a2e")) +
    scale_fill_manual(
      values = c("FALSE" = "#B8D4E8", "TRUE" = "#E8453C"),
      labels = c("Reference allele", "Alternative allele"),
      name   = NULL,
      drop   = FALSE
    ) +
    labs(
      title    = "GSDMB Locus Haplotype Composition",
      subtitle = sprintf(
        "%d established SNPs (gnomAD NFE AF > 1%%) | EM estimation across %d samples | Haplotypes shown: %.1f%% of all chromosomes",
        n_snps,
        nrow(a1),
        sum(global_freqs$Frequency) * 100
      ),
      x = "SNP Position  |  Change  |  Gene",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 8.5, colour = "grey40"),
      axis.text.x     = element_text(size = 7, family = "mono",
                                     angle = 0, hjust = 0.5, lineheight = 0.95,
                                     colour = axis_tick_cols),
      axis.title.x    = element_text(size = 8.5, colour = "grey40", face = "italic"),
      axis.text.y     = element_text(size = 10, face = "bold", family = "mono"),
      legend.position = "bottom",
      legend.key.size = unit(0.5, "cm"),
      panel.grid      = element_blank()
    )

  # Adjust figure height dynamically based on number of haplotypes
  fig_height <- max(3.5, nrow(global_freqs) * 0.6 + 2.5)
  fig_width  <- max(10,  n_snps * 1.0 + 2.0)

  ggsave(file.path(OUT_DIR, "19_Haplotype_Composition.png"),
         p4, width = fig_width, height = fig_height, dpi = 300)
  cat("  ✓ Haplotype composition tile figure saved\n")

}, error = function(e) {
  cat(sprintf("  ⚠ Haplotype composition plot failed: %s\n", conditionMessage(e)))
})

# =============================================================================
# COMPLETION SUMMARY
# =============================================================================
cat("\n======================================================================\n")
cat("✓ HAPLOTYPE ANALYSIS COMPLETE\n")
cat("======================================================================\n\n")
cat(sprintf("Outputs in: %s\n\n", OUT_DIR))
cat("  19_Haplotype_Results.xlsx            — full results workbook\n")
cat("  19_Global_Haplotype_Frequencies.png  — global EM frequency bar chart\n")
cat("  19_Haplotype_Freqs_By_Group.png      — frequency by sample group\n")
cat("  19_Score_Test_Results.png            — association p-value bar chart\n")
cat("  19_Haplotype_Composition.png         — tile figure: allele at each SNP per haplotype\n\n")

cat("SUMMARY OF BEST P-VALUES PER COMPARISON:\n")
print(summary_df, row.names = FALSE)
cat("\n")
cat("NOTE: 'Best_Haplotype_P' is the minimum Fisher p across individual\n")
cat("haplotypes within each comparison. It is not an omnibus global test.\n")
cat("Interpret with caution given multiple testing across haplotypes.\n\n")
