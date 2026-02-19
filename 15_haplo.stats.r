#!/usr/bin/env Rscript
# =============================================================================
# Script 19b: Haplotype Analysis with haplo.stats
# =============================================================================
# Uses BEAGLE-phased genotypes to:
#   1. Filter to established SNPs (gnomAD NFE AF > 1%) from script 11
#   2. Estimate haplotype frequencies (EM algorithm)
#   3. Test haplotype-phenotype associations for all group comparisons
#   4. Export results + plots
#
# Run: Rscript 17_haplo_stats.r
# =============================================================================

suppressPackageStartupMessages({
  library(haplo.stats)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

# --- CONFIGURATION ---
HAPLO_DIR  <- "/home/gadeaalonsoj/tfm/gsdmb_final_results/19_haplotype_phased"
GENO_FILE  <- file.path(HAPLO_DIR, "phased_genotypes.tsv")
META_FILE  <- file.path(HAPLO_DIR, "sample_metadata.tsv")
OUT_DIR    <- "/home/gadeaalonsoj/tfm/gsdmb_final_results/19_haplo_stats_results"
OUT_XLSX   <- file.path(OUT_DIR, "19_Haplotype_Results.xlsx")

# Main annotated report — used to extract established SNP positions (gnomAD NFE AF > 1%)
# This matches the exact filter used in script 11
ANNOTATED_REPORT <- "/home/gadeaalonsoj/tfm/gsdmb_final_results/GSDMB_Annotated_Report_Fixed.xlsx"

# Minimum haplotype frequency to report (filter very rare haplotypes)
MIN_HAPLOTYPE_FREQ <- 0.02

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("======================================================================\n")
cat("SCRIPT 19b: HAPLOTYPE ANALYSIS WITH haplo.stats\n")
cat("======================================================================\n\n")

# =============================================================================
# 1. Load data
# =============================================================================
cat("STEP 1: Loading phased genotype data\n")
cat("----------------------------------------------------------------------\n")

geno_raw  <- read.table(GENO_FILE, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
meta      <- read.table(META_FILE, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)

cat(sprintf("✓ Genotype matrix loaded: %d SNPs × %d samples\n",
            nrow(geno_raw), ncol(geno_raw) - 4))
cat(sprintf("✓ Metadata: %d samples\n", nrow(meta)))

# =============================================================================
# 2. Filter to established SNPs from script 11 (gnomAD NFE AF > 1%)
# This removes the ~124 off-target / noise SNPs and keeps only the
# biologically meaningful ones identified in the annotation pipeline.
# =============================================================================
cat("\nSTEP 2: Filtering to established SNPs from script 11\n")
cat("----------------------------------------------------------------------\n")

if (!file.exists(ANNOTATED_REPORT)) {
  stop(sprintf("Annotated report not found: %s\n  Run script 07 first.",
               ANNOTATED_REPORT))
}

cat(sprintf("  Reading from: %s\n", basename(ANNOTATED_REPORT)))
annot <- read.xlsx(ANNOTATED_REPORT, sheet = "Biological_Annotations")

# Case-insensitive column finder helper
find_col_ci <- function(df, pattern) {
  grep(pattern, colnames(df), ignore.case = TRUE, value = TRUE)[1]
}

pos_col <- find_col_ci(annot, "^pos$|^position$")
nfe_col <- find_col_ci(annot, "gnomad.*nfe.*af|nfe.*af")

if (is.na(pos_col)) {
  stop("Cannot find POS column in GSDMB_Annotated_Report_Fixed.xlsx.\n",
       "Available columns: ", paste(colnames(annot), collapse = ", "))
}
if (is.na(nfe_col)) {
  stop("Cannot find gnomAD NFE AF column in GSDMB_Annotated_Report_Fixed.xlsx.\n",
       "Available columns: ", paste(colnames(annot), collapse = ", "))
}

cat(sprintf("  POS column:          '%s'\n", pos_col))
cat(sprintf("  gnomAD NFE AF column:'%s'\n", nfe_col))

# Apply the same filter as script 11: gnomAD NFE AF > 1%
annot[[nfe_col]] <- suppressWarnings(as.numeric(annot[[nfe_col]]))
snps_established <- annot[!is.na(annot[[nfe_col]]) & annot[[nfe_col]] > 0.01, ]

# Extract unique positions
established_positions <- unique(as.integer(snps_established[[pos_col]]))
established_positions <- established_positions[!is.na(established_positions)]

cat(sprintf("  Established SNPs from script 11: %d positions\n",
            length(established_positions)))

# Filter phased genotype matrix to only these positions
keep <- geno_raw$POS %in% established_positions
n_before <- nrow(geno_raw)
geno_raw  <- geno_raw[keep, ]
n_after   <- nrow(geno_raw)

cat(sprintf("  SNPs before filter: %d\n", n_before))
cat(sprintf("  SNPs after filter:  %d\n", n_after))

if (n_after == 0) {
  stop("No SNPs remain after filtering.\n",
       "Check that POS values in the phased genotype TSV match those in\n",
       "GSDMB_Annotated_Report_Fixed.xlsx (gnomAD NFE AF > 0.01 variants).\n",
       "First few positions in phased file:    ",
       paste(head(as.integer(geno_raw$POS)), collapse = ", "), "\n",
       "First few established positions found: ",
       paste(head(established_positions), collapse = ", "))
}

# Report any established SNPs missing from the phased file
missing_from_phased <- established_positions[
  !established_positions %in% geno_raw$POS
]
if (length(missing_from_phased) > 0) {
  cat(sprintf("  ⚠ %d established SNP(s) not found in phased file (may have been\n",
              length(missing_from_phased)))
  cat("    filtered out by BEAGLE or bcftools):\n")
  cat(sprintf("    %s\n", paste(missing_from_phased, collapse = ", ")))
}

cat(sprintf("✓ Proceeding with %d SNPs\n", n_after))

# =============================================================================
# 3. Build SNP info and sample columns
# =============================================================================
snp_info   <- geno_raw[, c("CHROM", "POS", "REF", "ALT")]
snp_labels <- paste0(geno_raw$POS, "_", geno_raw$REF, ">", geno_raw$ALT)

# Sample columns (everything after CHROM POS REF ALT)
sample_cols <- colnames(geno_raw)[-(1:4)]

# =============================================================================
# 4. Build allele matrix for haplo.stats
# haplo.stats needs two columns per SNP: allele1 and allele2 (phased)
# Input GT format from BEAGLE: "0|1", "1|0", "0|0", "1|1"
# =============================================================================
cat("\nSTEP 3: Parsing phased genotypes into allele matrix\n")
cat("----------------------------------------------------------------------\n")

parse_phased_gt <- function(gt_vec) {
  allele1 <- sapply(gt_vec, function(g) {
    if (is.na(g) || g %in% c(".", ".|.")) return(NA)
    as.integer(strsplit(g, "\\|")[[1]][1])
  })
  allele2 <- sapply(gt_vec, function(g) {
    if (is.na(g) || g %in% c(".", ".|.")) return(NA)
    as.integer(strsplit(g, "\\|")[[1]][2])
  })
  list(a1 = as.integer(allele1), a2 = as.integer(allele2))
}

n_samples <- length(sample_cols)
n_snps    <- nrow(geno_raw)

allele1_mat <- matrix(NA_integer_, nrow = n_samples, ncol = n_snps)
allele2_mat <- matrix(NA_integer_, nrow = n_samples, ncol = n_snps)
rownames(allele1_mat) <- sample_cols
rownames(allele2_mat) <- sample_cols
colnames(allele1_mat) <- snp_labels
colnames(allele2_mat) <- snp_labels

for (j in seq_len(n_snps)) {
  gt_vec <- as.character(geno_raw[j, sample_cols])
  parsed <- parse_phased_gt(gt_vec)
  allele1_mat[, j] <- parsed$a1
  allele2_mat[, j] <- parsed$a2
}

# haplo.stats uses 1-based allele coding: 1 = ref, 2 = alt
allele1_mat_hap <- allele1_mat + 1L   # 0 -> 1 (ref), 1 -> 2 (alt)
allele2_mat_hap <- allele2_mat + 1L

cat(sprintf("✓ Allele matrices built: %d samples × %d SNPs\n",
            n_samples, n_snps))

# Align metadata to sample order
meta_aligned  <- meta[match(sample_cols, meta$Sample), ]
missing_meta  <- sum(is.na(meta_aligned$Sample))
if (missing_meta > 0) {
  cat(sprintf("  ⚠ %d samples have no metadata — will be excluded\n", missing_meta))
}
valid_samples <- !is.na(meta_aligned$Sample)
meta_aligned  <- meta_aligned[valid_samples, ]
a1 <- allele1_mat_hap[valid_samples, , drop = FALSE]
a2 <- allele2_mat_hap[valid_samples, , drop = FALSE]

# =============================================================================
# 5. Estimate global haplotype frequencies (full cohort)
# =============================================================================
cat("\nSTEP 4: Estimating global haplotype frequencies (EM algorithm)\n")
cat("----------------------------------------------------------------------\n")

build_geno_matrix <- function(a1_mat, a2_mat) {
  n    <- ncol(a1_mat)
  cols <- vector("list", n * 2)
  for (j in seq_len(n)) {
    cols[[2*j - 1]] <- a1_mat[, j]
    cols[[2*j]]     <- a2_mat[, j]
  }
  as.data.frame(cols)
}

geno_mat_full <- build_geno_matrix(a1, a2)

hap_em_global <- haplo.em(
  geno        = geno_mat_full,
  locus.label = snp_labels,
  miss.val    = c(0, NA),
  control     = haplo.em.control(min.posterior = 0.001)
)

global_freqs <- hap_em_global$haplotype
global_freqs <- cbind(global_freqs,
                      Frequency = hap_em_global$hap.prob,
                      Count     = round(hap_em_global$hap.prob * nrow(a1) * 2))
global_freqs <- as.data.frame(global_freqs)
global_freqs <- global_freqs[order(-global_freqs$Frequency), ]
global_freqs <- global_freqs[global_freqs$Frequency >= MIN_HAPLOTYPE_FREQ, ]
global_freqs$Haplotype_ID <- paste0("H", seq_len(nrow(global_freqs)))

cat(sprintf("✓ %d haplotypes estimated (freq ≥ %.0f%%)\n",
            nrow(global_freqs), MIN_HAPLOTYPE_FREQ * 100))
cat("\nGlobal haplotype frequencies:\n")
print(global_freqs[, c("Haplotype_ID", "Frequency", "Count")], row.names = FALSE)

# =============================================================================
# 6. Association testing for all group comparisons
# =============================================================================
cat("\n\nSTEP 5: Haplotype association testing\n")
cat("----------------------------------------------------------------------\n")

comparisons <- list(
  list(name        = "Breast_Tumour_vs_Healthy",
       filter      = meta_aligned$Cohort == "Breast",
       outcome_col = "Tissue",
       case        = "Tumour", control = "Healthy"),

  list(name        = "Endometrium_Tumour_vs_Healthy",
       filter      = meta_aligned$Cohort == "Endometrium",
       outcome_col = "Tissue",
       case        = "Tumour", control = "Healthy"),

  list(name        = "Breast_vs_Endometrium_Tumour",
       filter      = meta_aligned$Tissue == "Tumour",
       outcome_col = "Cohort",
       case        = "Breast", control = "Endometrium"),

  list(name        = "Breast_vs_Endometrium_Healthy",
       filter      = meta_aligned$Tissue == "Healthy",
       outcome_col = "Cohort",
       case        = "Breast", control = "Endometrium"),

  list(name        = "All_Tumour_vs_Healthy",
       filter      = rep(TRUE, nrow(meta_aligned)),
       outcome_col = "Tissue",
       case        = "Tumour", control = "Healthy")
)
all_results   <- list()
score_results <- list()

cat(sprintf("  Unique Cohort values: %s\n",
            paste(sort(unique(meta_aligned$Cohort)), collapse = ", ")))
cat(sprintf("  Unique Tissue values: %s\n",
            paste(sort(unique(meta_aligned$Tissue)), collapse = ", ")))

# Association testing: haplo.em frequencies + Fisher's exact test per haplotype.
# We compute haplotype carrier status per sample directly from the phased
# allele matrices (no reliance on haplo.em internal object structure).

# Convert allele matrices back to haplotype strings for direct comparison
# a1/a2 are coded 1=ref, 2=alt — convert to 0/1 strings
make_hap_strings <- function(a_mat) {
  # a_mat: samples x SNPs, values 1 or 2 (1=ref, 2=alt)
  apply(a_mat - 1L, 1, paste, collapse = "")
}

run_haplo_assoc <- function(a1_sub, a2_sub, y, snp_labels, comp_name,
                             case_label, control_label, n_case, n_control) {

  # Get haplotype string for each chromosome (2 per sample)
  hap1_strings <- make_hap_strings(a1_sub)  # length = n_samples
  hap2_strings <- make_hap_strings(a2_sub)  # length = n_samples

  # All observed haplotypes
  all_haps <- table(c(hap1_strings, hap2_strings))
  # Keep only those with >= 3 copies
  common_haps <- names(all_haps[all_haps >= 3])
  n_haps <- length(common_haps)
  cat(sprintf("    Distinct haplotypes with >= 3 copies: %d\n", n_haps))

  if (n_haps < 2) {
    cat("    WARNING: Too few haplotypes -- skipping\n")
    return(NULL)
  }

  results <- list()

  for (h in seq_along(common_haps)) {
    hap <- common_haps[h]

    # Dosage = number of chromosomes carrying this haplotype (0, 1, or 2)
    dosage <- as.integer(hap1_strings == hap) +
              as.integer(hap2_strings == hap)

    # Skip if no variance
    if (var(dosage) == 0) next

    # Carrier status (dosage >= 1)
    carrier <- as.integer(dosage >= 1)

    # Frequencies in cases vs controls
    freq_case    <- mean(dosage[y == 1]) / 2
    freq_control <- mean(dosage[y == 0]) / 2

    # Fisher's exact on carrier vs non-carrier
    tbl <- table(carrier, y)
    if (nrow(tbl) < 2) next
    fish <- fisher.test(tbl)

    # Also logistic regression for OR with CI
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
      or <- fish$estimate; or_lo <- NA; or_hi <- NA; p_glm <- fish$p.value
    }

    # Haplotype label: position>allele format
    snp_pos  <- sapply(strsplit(snp_labels, "_"), `[`, 1)
    hap_bits <- strsplit(hap, "")[[1]]
    hap_name <- paste0("H", h, "(",
                       paste(snp_pos[hap_bits == "1"], collapse = "/"), ")")

    results[[h]] <- data.frame(
      Haplotype      = paste0("H", h),
      Hap_Definition = hap_name,
      Global_Freq    = round(sum(all_haps[hap]) / sum(all_haps), 4),
      Freq_Case      = round(freq_case, 4),
      Freq_Control   = round(freq_control, 4),
      OR             = round(or, 4),
      OR_95CI_Lo     = round(or_lo, 4),
      OR_95CI_Hi     = round(or_hi, 4),
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
  res_df$FDR         <- p.adjust(res_df$P_Fisher, method = "BH")
  res_df$Significant <- res_df$FDR < 0.05
  res_df             <- res_df[order(res_df$P_Fisher), ]

  global_p <- min(res_df$P_Fisher)
  cat(sprintf("    Best haplotype p (Fisher): %.4f | FDR: %.4f\n",
              global_p, min(res_df$FDR)))

  list(table = res_df, global_p = global_p)
}

for (comp in comparisons) {
  cat(sprintf("\n  Running: %s\n", comp$name))

  idx       <- which(comp$filter &
                     meta_aligned[[comp$outcome_col]] %in% c(comp$case, comp$control))
  meta_sub  <- meta_aligned[idx, ]
  a1_sub    <- a1[idx, , drop = FALSE]
  a2_sub    <- a2[idx, , drop = FALSE]
  n_case    <- sum(meta_sub[[comp$outcome_col]] == comp$case)
  n_control <- sum(meta_sub[[comp$outcome_col]] == comp$control)

  cat(sprintf("    %s (case): %d | %s (control): %d\n",
              comp$case, n_case, comp$control, n_control))

  if (n_case < 3 || n_control < 3) {
    cat("    WARNING: Too few samples -- skipping\n")
    next
  }

  y <- as.integer(meta_sub[[comp$outcome_col]] == comp$case)

  res <- run_haplo_assoc(a1_sub, a2_sub, y, snp_labels,
                         comp$name, comp$case, comp$control,
                         n_case, n_control)

  if (!is.null(res)) {
    score_results[[comp$name]] <- list(
      table    = res$table,
      global_p = res$global_p
    )
    all_results[[comp$name]] <- res$table
    cat(sprintf("    Saved %d haplotype results\n", nrow(res$table)))
  }
}

# 7. Save results
# =============================================================================
cat("\n\nSTEP 6: Saving results\n")
cat("----------------------------------------------------------------------\n")

wb <- createWorkbook()

# SNP filter provenance — which positions were used
snp_used_df <- data.frame(
  POS       = geno_raw$POS,
  REF       = geno_raw$REF,
  ALT       = geno_raw$ALT,
  SNP_Label = snp_labels
)
addWorksheet(wb, "SNPs_Used")
writeData(wb, "SNPs_Used", snp_used_df)

# Global frequencies
addWorksheet(wb, "Global_Haplotype_Freqs")
writeData(wb, "Global_Haplotype_Freqs", global_freqs)

# Per-comparison results
for (comp_name in names(all_results)) {
  sheet_name <- substr(gsub("_", " ", comp_name), 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, all_results[[comp_name]])
}

# Summary of global p-values
summary_df <- data.frame(
  Comparison     = names(score_results),
  Global_Score_P = sapply(score_results, function(x) x$global_p),
  Significant    = sapply(score_results, function(x) x$global_p < 0.05)
)
addWorksheet(wb, "Summary_Global_Tests")
writeData(wb, "Summary_Global_Tests", summary_df)

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("✓ Results saved to: %s\n", OUT_XLSX))

# =============================================================================
# 8. Plots
# =============================================================================
cat("\nSTEP 7: Generating plots\n")
cat("----------------------------------------------------------------------\n")

# Plot 1: Global haplotype frequency bar chart
p1 <- ggplot(global_freqs, aes(x = reorder(Haplotype_ID, -Frequency),
                                y = Frequency * 100)) +
  geom_col(fill = "#3498db", colour = "black", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Frequency * 100)),
            vjust = -0.4, size = 3.5) +
  labs(title    = "Global Haplotype Frequencies",
       subtitle = sprintf("EM estimation across %d samples | %d established SNPs",
                          nrow(a1), n_snps),
       x = "Haplotype", y = "Frequency (%)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "19_Global_Haplotype_Frequencies.png"),
       p1, width = 8, height = 5, dpi = 300)
cat("  ✓ Global frequency plot saved\n")

# Plot 2: Haplotype frequencies by group (stacked bar)
group_freq_list <- list()
for (grp in unique(meta_aligned$Group)) {
  idx <- which(meta_aligned$Group == grp)
  if (length(idx) < 3) next
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
    group_freq_list[[grp]] <- hap_df
  }, error = function(e) NULL)
}

if (length(group_freq_list) > 0) {
  group_freq_df <- bind_rows(lapply(names(group_freq_list), function(grp) {
    df             <- group_freq_list[[grp]]
    df             <- df[order(-df$Frequency), ]
    df$Rank        <- seq_len(nrow(df))
    df$Haplotype_ID <- paste0("H", df$Rank)
    df[df$Frequency >= MIN_HAPLOTYPE_FREQ,
       c("Group", "Haplotype_ID", "Frequency")]
  }))

  p2 <- ggplot(group_freq_df,
               aes(x = Group, y = Frequency * 100, fill = Haplotype_ID)) +
    geom_col(colour = "black", linewidth = 0.3) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Haplotype Frequencies by Sample Group",
         x = NULL, y = "Frequency (%)", fill = "Haplotype") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(face = "bold")) +
    ylim(0, 105)

  ggsave(file.path(OUT_DIR, "19_Haplotype_Freqs_By_Group.png"),
         p2, width = 10, height = 6, dpi = 300)
  cat("  ✓ Group frequency plot saved\n")
}

# Plot 3: Haplotype association results — OR forest plot + p-value bar chart
if (length(score_results) > 0) {
  score_plot_data <- bind_rows(lapply(names(score_results), function(cn) {
    tbl      <- score_results[[cn]]$table
    global_p <- score_results[[cn]]$global_p
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    data.frame(
      Comparison  = cn,
      Global_P    = global_p,
      Haplotype   = tbl$Haplotype,
      OR          = tbl$OR,
      OR_Lo       = tbl$OR_95CI_Lo,
      OR_Hi       = tbl$OR_95CI_Hi,
      P_Value     = tbl$P_Fisher,
      FDR         = tbl$FDR,
      Freq_Case   = tbl$Freq_Case,
      Freq_Control = tbl$Freq_Control,
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(score_plot_data) && nrow(score_plot_data) > 0) {
    score_plot_data$OR_Lo <- ifelse(is.na(score_plot_data$OR_Lo),
                                    score_plot_data$OR * 0.5,
                                    score_plot_data$OR_Lo)
    score_plot_data$OR_Hi <- ifelse(is.na(score_plot_data$OR_Hi),
                                    score_plot_data$OR * 2.0,
                                    score_plot_data$OR_Hi)

    # Panel A: -log10 p-value bar chart
    p3a <- ggplot(score_plot_data,
                  aes(x = Haplotype, y = -log10(P_Value),
                      fill = FDR < 0.05)) +
      geom_col(colour = "black", alpha = 0.85) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                 colour = "red", alpha = 0.7) +
      scale_fill_manual(values = c("FALSE" = "#95a5a6", "TRUE" = "#e74c3c"),
                        labels = c("FDR ≥ 0.05", "FDR < 0.05"),
                        name   = "Significance") +
      facet_wrap(~ Comparison, scales = "free_x", ncol = 2) +
      labs(title    = "Haplotype Association: Fisher Exact Test",
           subtitle = "Red dashed line = nominal p < 0.05 | Red fill = FDR < 0.05",
           x = "Haplotype", y = "-log10(p-value)") +
      theme_bw(base_size = 11) +
      theme(strip.text  = element_text(face = "bold", size = 7),
            plot.title  = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(OUT_DIR, "19_Score_Test_Results.png"),
           p3a, width = 14, height = 10, dpi = 300)
    cat("  ✓ Association p-value plot saved\n")

    # Panel B: OR forest plot for All_Tumour_vs_Healthy
    all_comp <- score_plot_data[
      score_plot_data$Comparison == "All_Tumour_vs_Healthy", ]

    if (nrow(all_comp) > 0) {
      all_comp <- all_comp[order(all_comp$OR), ]
      all_comp$y <- seq_len(nrow(all_comp))

      p3b <- ggplot(all_comp, aes(x = OR, y = y, colour = FDR < 0.05)) +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = OR_Lo, xmax = OR_Hi), height = 0.3) +
        geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
        scale_y_continuous(breaks = all_comp$y, labels = all_comp$Haplotype) +
        scale_colour_manual(values = c("FALSE" = "#95a5a6", "TRUE" = "#e74c3c"),
                            labels = c("FDR ≥ 0.05", "FDR < 0.05"),
                            name   = "Significance") +
        labs(title    = "Haplotype Odds Ratios: All Tumour vs Healthy",
             subtitle = "Error bars = 95% CI | Dashed line = OR 1.0",
             x = "Odds Ratio", y = "Haplotype") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold"))

      ggsave(file.path(OUT_DIR, "19_OR_Forest_Plot.png"),
             p3b, width = 8, height = 6, dpi = 300)
      cat("  ✓ OR forest plot saved\n")
    }
  }
}


cat("\n======================================================================\n")
cat("✓ HAPLO.STATS ANALYSIS COMPLETE\n")
cat("======================================================================\n\n")
cat(sprintf("Outputs in: %s\n", OUT_DIR))
cat("  19_Haplotype_Results.xlsx  (includes SNPs_Used sheet for traceability)\n")
cat("  19_Global_Haplotype_Frequencies.png\n")
cat("  19_Haplotype_Freqs_By_Group.png\n")
cat("  19_Score_Test_Results.png\n\n")

cat("SUMMARY OF GLOBAL SCORE TESTS:\n")
print(summary_df, row.names = FALSE)
cat("\n")
