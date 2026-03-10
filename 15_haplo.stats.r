#!/usr/bin/env Rscript
# Haplotype analysis for the GSDMB locus (chr17:39.8-39.99 Mb).
#
# Three complementary levels of analysis are run in sequence:
#   1. Full-region haplotypes — all whitelist SNPs treated as a single block.
#      Provides the broadest overview of haplotype structure across the locus.
#   2. LD-block haplotypes — the locus is partitioned into internally coherent
#      blocks using a Gabriel-style algorithm (all pairwise r² >= threshold).
#      Shorter blocks have fewer SNPs and therefore more phasing power.
#   3. Gene-specific haplotypes — SNPs are further grouped by gene (GSDMB here),
#      optionally sub-partitioned into LD blocks within the gene.
#
# For each region and level, the script estimates haplotype frequencies via
# EM (haplo.em), tests tumour vs healthy carrier differences by Fisher's exact
# test and Firth-penalised logistic regression, and applies BH-FDR correction
# within each region.
#
# Input files (produced by scripts 11 and 15_haplotypes.sh):
#   phased_genotypes.tsv       — phased GT matrix from BEAGLE
#   sample_metadata.tsv        — cohort/tissue labels
#   11_Master_Unique_SNP_Summary.xlsx — whitelist defining the target SNP set
#   GSDMB_Annotated_Report_Fixed.xlsx — gene/position annotation for each SNP
#
# Usage:  Rscript 15_haplo_stats.r

suppressPackageStartupMessages({
  library(haplo.stats)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(readxl)
  library(stringr)
  library(scales)
  library(forcats)
  for (pkg in c("patchwork", "logistf", "RColorBrewer", "gridExtra")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  # survminer depends on ggpubr -> ggrepel; install with dependencies=TRUE so the
  # chain resolves. If installation fails in a restricted environment (e.g. conda R)
  # the function make_km_plot exits silently and no KM plots are produced.
  for (pkg in c("ggrepel", "ggpubr", "survival", "survminer")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      tryCatch(
        install.packages(pkg, repos = "https://cloud.r-project.org",
                         dependencies = TRUE),
        warning = function(w) invisible(NULL),
        error   = function(e) invisible(NULL)
      )
    }
  }
  library(patchwork)
  library(logistf)
})

# Explicitly namespace dplyr verbs to avoid conflicts with base R and other packages
# (plyr, MASS, and stats all export functions with overlapping names)
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
rename    <- dplyr::rename
arrange   <- dplyr::arrange
distinct  <- dplyr::distinct
pull      <- dplyr::pull
group_by  <- dplyr::group_by
ungroup   <- dplyr::ungroup
summarise <- dplyr::summarise
left_join <- dplyr::left_join
inner_join <- dplyr::inner_join
bind_rows <- dplyr::bind_rows

# --- Global ggplot theme -----------------------------------------------------
# A clean publication theme applied to all figures for visual consistency
pub_theme <- theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    plot.title        = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle     = element_text(size = 9, colour = "grey40", hjust = 0),
    axis.title        = element_text(size = 10, face = "bold"),
    axis.text         = element_text(size = 9, colour = "black"),
    axis.line         = element_line(linewidth = 0.4),
    axis.ticks        = element_line(linewidth = 0.4),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey95", colour = "grey70"),
    strip.text        = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 8),
    legend.title      = element_text(size = 9, face = "bold"),
    legend.key.size   = unit(0.45, "cm"),
    plot.margin       = margin(8, 10, 6, 8)
  )

theme_set(pub_theme)

# Fixed colour map for the four cohort-tissue groups; used in cohort frequency plots
GROUP_COLS <- c(
  "Breast_Healthy"      = "#4393C3",
  "Breast_Tumour"       = "#D6604D",
  "Endometrium_Healthy" = "#74ADD1",
  "Endometrium_Tumour"  = "#A50026",
  "All_Healthy"         = "#92C5DE",
  "All_Tumour"          = "#F4A582"
)

# Wrapper around ggsave with consistent DPI and units
SAVE <- function(file, plot, width = 10, height = 6) {
  ggsave(file, plot, width = width, height = height, dpi = 300,
         units = "in", bg = "white")
  cat("  Saved:", basename(file), "\n")
}

# Excel sheet names must be <= 31 characters and cannot contain :\/? etc.
safe_sheet <- function(x) {
  x <- gsub("[:/\\?*\\[\\]]", "_", x)
  substr(x, 1, 31)
}

# --- Configuration -----------------------------------------------------------
BASE_DIR  <- "/home/gadeaalonsoj/tfm/gsdmb_final_results"
HAPLO_DIR <- file.path(BASE_DIR, "19_haplotype_phased")
OUT_DIR   <- file.path(BASE_DIR, "19_haplo_stats_results")

GENO_FILE <- file.path(HAPLO_DIR, "phased_genotypes.tsv")
META_FILE <- file.path(HAPLO_DIR, "sample_metadata.tsv")
WHITELIST <- file.path(BASE_DIR, "11_Master_Unique_SNP_Summary.xlsx")    # target SNP set
ANNOTATED <- file.path(BASE_DIR, "GSDMB_Annotated_Report_Fixed.xlsx")    # gene/position annotations

OUT_XLSX  <- file.path(OUT_DIR, "19_Haplotype_Results_v7_blocks_and_genes.xlsx")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Frequency thresholds:
#   MAIN_MIN_FREQ  — haplotypes below 1% excluded from the full-region overview;
#                    rare haplotypes are too infrequent for reliable carrier comparison
#   SUPP_MIN_FREQ  — looser threshold for LD-block and gene analyses where fewer SNPs
#                    per block means more distinct haplotypes at low frequency
MAIN_MIN_FREQ        <- 0.01
SUPP_MIN_FREQ        <- 0.005
MIN_HAP_COPIES       <- 3    # haplotype must be observed >= 3 times (across both alleles) to be tested
MIN_SAMPLES_PER_ARM  <- 3    # minimum tumour or healthy samples required to run any test
MIN_GROUP_SAMPLES_EM <- 3

# LD block parameters
RUN_BLOCK_ANALYSIS   <- TRUE
LD_BLOCK_THRESHOLD   <- 0.60   # minimum pairwise r² for all pairs within a candidate block
MIN_BLOCK_SNPS       <- 2
MAX_BLOCK_SNPS       <- 12     # cap prevents blocks so long that EM becomes unstable

# Gene-specific analysis parameters
RUN_GENE_SPECIFIC           <- TRUE
GENE_LIST_SUPP              <- c("GSDMB")
GENE_SPECIFIC_USE_LD_BLOCKS <- TRUE   # sub-partition gene SNPs into LD blocks before phasing

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

parse_phased_call <- function(x) {
  # Parse a single phased genotype string (e.g. "0|1" or "1/0") into a
  # length-2 integer vector of allele indices. Returns NA for missing calls.
  x <- trimws(as.character(x))
  if (is.na(x) || x == "" || x == "." || x == "./." || x == ".|.") {
    return(c(NA_integer_, NA_integer_))
  }
  bits <- strsplit(x, "[|/]")[[1]]
  if (length(bits) != 2) return(c(NA_integer_, NA_integer_))
  out <- suppressWarnings(as.integer(bits))
  if (any(is.na(out))) return(c(NA_integer_, NA_integer_))
  out
}

build_allele_matrices <- function(geno_df, sample_cols) {
  # Convert the raw genotype character matrix (one cell per sample per SNP)
  # into two integer allele matrices: a1 (haplotype 1) and a2 (haplotype 2).
  # Allele codes follow haplo.stats convention: 1 = REF, 2 = ALT.
  #
  # The matrix is transposed so rows are samples and columns are SNPs, matching
  # the orientation expected by haplo.em. Parsing is vectorised using sub()
  # regex capture groups rather than looping over each cell for speed.
  n_samples <- length(sample_cols)
  n_snps    <- nrow(geno_df)

  raw_mat <- t(as.matrix(geno_df[, sample_cols, drop = FALSE]))  # samples x SNPs

  parse_allele <- function(x, allele_idx) {
    x   <- trimws(as.character(x))
    bad <- is.na(x) | x == "" | x == "." | x == "./." | x == ".|."
    # Capture haplotype 1 or 2 from phased format "A|B"
    a <- if (allele_idx == 1L) {
      sub("([0-9]+)[|/]([0-9]+)", "\\1", x)
    } else {
      sub("([0-9]+)[|/]([0-9]+)", "\\2", x)
    }
    out <- suppressWarnings(as.integer(a)) + 1L  # shift 0/1 VCF coding to 1/2 haplo.stats coding
    out[bad | is.na(out)] <- NA_integer_
    out
  }

  a1 <- matrix(parse_allele(raw_mat, 1L), nrow = n_samples, ncol = n_snps,
               dimnames = list(sample_cols, geno_df$rsID_clean))
  a2 <- matrix(parse_allele(raw_mat, 2L), nrow = n_samples, ncol = n_snps,
               dimnames = list(sample_cols, geno_df$rsID_clean))

  list(a1 = a1, a2 = a2)
}

build_geno_matrix <- function(a1_mat, a2_mat) {
  # Interleave haplotype 1 and haplotype 2 columns into the paired-column
  # format required by haplo.em: SNP1_hap1, SNP1_hap2, SNP2_hap1, SNP2_hap2...
  stopifnot(nrow(a1_mat) == nrow(a2_mat), ncol(a1_mat) == ncol(a2_mat))
  out <- vector("list", ncol(a1_mat) * 2)
  for (j in seq_len(ncol(a1_mat))) {
    out[[2 * j - 1]] <- a1_mat[, j]
    out[[2 * j]]     <- a2_mat[, j]
  }
  as.data.frame(out, check.names = FALSE)
}

make_binary_hap_strings <- function(a_mat) {
  # Collapse each row of an allele matrix into a binary string (e.g. "0110")
  # where 0 = REF allele and 1 = ALT allele. Used to identify and count
  # haplotype combinations without relying on EM output.
  apply(a_mat - 1L, 1, paste, collapse = "")
}

make_em_freq_table <- function(a1_mat, a2_mat, snp_meta, min_freq) {
  # Estimate haplotype frequencies using the EM algorithm implemented in
  # haplo.em. EM is necessary here because phasing (even from BEAGLE) is
  # probabilistic — EM integrates over phase uncertainty to give population-
  # level frequency estimates rather than just counting observed haplotypes.
  #
  # miss.val = c(0, NA): both 0 and NA are treated as missing data codes.
  # min.posterior = 0.001: haplo-sample pairs with posterior probability below
  # this threshold are excluded, reducing noise from unlikely phase assignments.
  geno_mat <- build_geno_matrix(a1_mat, a2_mat)
  em_fit <- haplo.em(
    geno        = geno_mat,
    locus.label = snp_meta$rsID_clean,
    miss.val    = c(0, NA),
    control     = haplo.em.control(min.posterior = 0.001)
  )

  freqs <- as.data.frame(em_fit$haplotype)
  names(freqs) <- snp_meta$rsID_clean
  freqs$Frequency <- as.numeric(em_fit$hap.prob)
  freqs$Count     <- round(freqs$Frequency * nrow(a1_mat) * 2)  # expected haplotype copies in pool
  freqs <- freqs %>%
    arrange(desc(Frequency)) %>%
    mutate(Haplotype_ID = paste0("H", seq_len(n())))  # H1 = most frequent

  # Convert the allele-coded matrix (1/2) back to binary (0/1) for display
  snp_matrix <- matrix(
    as.numeric(as.matrix(freqs[, snp_meta$rsID_clean, drop = FALSE])),
    nrow = nrow(freqs), ncol = nrow(snp_meta)
  )
  freqs$Allele_String <- apply(snp_matrix, 1, function(x) paste0(x - 1L, collapse = ""))
  freqs$Alt_Alleles   <- rowSums(snp_matrix - 1L)  # number of ALT positions in each haplotype

  # Human-readable description: list the rsIDs that carry the ALT allele
  annot_strings <- apply(snp_matrix, 1, function(row) {
    alt_idx <- which(row == 2)
    if (length(alt_idx) == 0) return("Reference-like")
    paste0(snp_meta$rsID_clean[alt_idx], collapse = "; ")
  })
  freqs$Variant_Content <- annot_strings

  shown <- freqs %>% filter(Frequency >= min_freq)
  list(all = freqs, shown = shown, fit = em_fit)
}

compute_ld_matrix <- function(a1_mat, a2_mat, snp_meta) {
  # Compute pairwise r² (the square of the Pearson correlation between allele
  # dosages) for all SNP pairs. Dosage at each site is the count of ALT alleles
  # carried across both haplotypes (0, 1, or 2), which is a standard proxy for
  # LD estimation from phased data. Pairs with fewer than 3 complete observations
  # are assigned NA to avoid spuriously high r² from minimal data.
  dosage <- (a1_mat - 1L) + (a2_mat - 1L)
  n_snps <- ncol(dosage)
  ld <- matrix(NA_real_, nrow = n_snps, ncol = n_snps,
               dimnames = list(snp_meta$rsID_clean, snp_meta$rsID_clean))
  diag(ld) <- 1
  for (i in seq_len(n_snps)) {
    for (j in i:n_snps) {
      xi   <- dosage[, i]
      xj   <- dosage[, j]
      keep <- complete.cases(xi, xj)
      if (sum(keep) < 3) {
        r2 <- NA_real_
      } else {
        r  <- suppressWarnings(cor(xi[keep], xj[keep]))
        r2 <- ifelse(is.na(r), NA_real_, r^2)
      }
      ld[i, j] <- r2
      ld[j, i] <- r2  # matrix is symmetric
    }
  }
  ld
}

define_ld_blocks <- function(ld_matrix, snp_meta, threshold = 0.6,
                             min_snps = 2, max_snps = 12) {
  # Gabriel-style LD block definition: a window of SNPs i:j is accepted as a
  # block only if ALL pairwise r² values within it are >= threshold (i.e. the
  # entire block is internally coherent, not just adjacent pairs).
  #
  # The algorithm scans left to right. At each starting SNP i it tries the
  # largest candidate window first (greedily maximising block length), then
  # shrinks until it finds a valid block or moves on. SNPs not assigned to
  # any block are skipped — they may be isolated variants in weak LD with
  # their neighbours.
  snp_names <- rownames(ld_matrix)
  n         <- length(snp_names)
  if (is.null(snp_names) || n < min_snps) return(list())

  # Check whether all lower-triangle r² values in a candidate window exceed threshold
  block_coherent <- function(name_idx) {
    nms   <- snp_names[name_idx]
    sub   <- ld_matrix[nms, nms, drop = FALSE]
    pairs <- sub[lower.tri(sub)]
    length(pairs) > 0 && all(!is.na(pairs) & pairs >= threshold)
  }

  blocks   <- list()
  block_id <- 1
  i        <- 1

  while (i <= n) {
    if (i + min_snps - 1L > n) break   # not enough SNPs left to form a block

    best_end <- NA
    max_j    <- min(n, i + max_snps - 1L)

    # Try longest window first; take first (longest) coherent window found
    for (j in max_j:(i + min_snps - 1L)) {
      if (block_coherent(i:j)) {
        best_end <- j
        break
      }
    }

    if (!is.na(best_end)) {
      blocks[[paste0("Block", block_id)]] <- i:best_end
      block_id <- block_id + 1
      i <- best_end + 1  # next block starts immediately after current one
    } else {
      i <- i + 1  # this SNP cannot start a valid block; advance
    }
  }
  blocks
}

run_assoc_for_region <- function(region_name, a1_mat, a2_mat, region_meta,
                                 metadata_df, min_copies = 3,
                                 comparison_set = c("Breast", "Endometrium", "All")) {
  # Test each haplotype for carrier frequency differences between tumour and
  # healthy tissue using two complementary methods:
  #
  # (1) Fisher's exact test on a 2x2 table (carrier / non-carrier x tumour / healthy).
  #     Chosen because sample sizes are small and cell counts are sometimes near zero.
  #
  # (2) Firth-penalised logistic regression (logistf). Standard logistic regression
  #     fails or gives infinite estimates when a haplotype is absent in one group
  #     (complete separation). Firth's penalisation shrinks the likelihood to avoid
  #     this, making it the recommended method for rare binary outcomes in small samples.
  #     Falls back to standard GLM if logistf fails, and to Fisher OR if GLM also fails.
  #
  # BH-FDR correction is applied once across ALL haplotype-comparison combinations
  # within a region (not per comparison group), which is more conservative but
  # avoids inflating significance by treating the three comparison groups as independent.

  # Count haplotype occurrences across both alleles for all samples
  hap1 <- make_binary_hap_strings(a1_mat)
  hap2 <- make_binary_hap_strings(a2_mat)
  hap_counts  <- sort(table(c(hap1, hap2)), decreasing = TRUE)
  tested_haps <- names(hap_counts[hap_counts >= min_copies])

  if (length(tested_haps) == 0) return(NULL)

  results <- list()

  comp_defs <- list(
    list(label = "Breast",      keep = metadata_df$Cohort == "Breast"),
    list(label = "Endometrium", keep = metadata_df$Cohort == "Endometrium"),
    list(label = "All",         keep = rep(TRUE, nrow(metadata_df)))
  )
  comp_defs <- comp_defs[sapply(comp_defs, function(x) x$label %in% comparison_set)]

  for (comp in comp_defs) {
    idx <- which(comp$keep & metadata_df$Tissue %in% c("Healthy", "Tumour"))
    if (length(idx) == 0) next

    meta_sub  <- metadata_df[idx, , drop = FALSE]
    y         <- as.integer(meta_sub$Tissue == "Tumour")   # 1 = case, 0 = control
    n_case    <- sum(y == 1)
    n_control <- sum(y == 0)
    if (n_case < MIN_SAMPLES_PER_ARM || n_control < MIN_SAMPLES_PER_ARM) next

    hap1_sub <- hap1[idx]
    hap2_sub <- hap2[idx]

    comp_rows <- list()
    for (h in tested_haps) {
      # Dosage = number of copies of haplotype h carried by each sample (0, 1, or 2)
      dosage <- as.integer(hap1_sub == h) + as.integer(hap2_sub == h)
      if (var(dosage) == 0) next  # no variation — cannot test

      carrier <- as.integer(dosage >= 1)
      tbl     <- table(carrier, y)
      if (nrow(tbl) < 2 || ncol(tbl) < 2) next  # degenerate table — skip

      fish <- fisher.test(tbl)
      freq_case    <- mean(dosage[y == 1], na.rm = TRUE) / 2  # carrier freq = mean dosage / 2
      freq_control <- mean(dosage[y == 0], na.rm = TRUE) / 2

      # Primary model: Firth-penalised logistic regression on dosage (additive model)
      fit <- tryCatch(
        logistf::logistf(y ~ dosage),
        error = function(e) NULL
      )

      if (!is.null(fit) && "dosage" %in% names(coef(fit))) {
        beta  <- coef(fit)["dosage"]
        ci    <- as.numeric(suppressWarnings(confint(fit, parm = "dosage")))
        or    <- exp(beta)
        or_lo <- exp(ci[1])
        or_hi <- exp(ci[2])
        p_mod <- fit$prob["dosage"]
      } else {
        # Fallback 1: standard logistic regression
        fit_glm <- tryCatch(suppressWarnings(glm(y ~ dosage, family = binomial())),
                            error = function(e) NULL)
        if (!is.null(fit_glm) && nrow(summary(fit_glm)$coefficients) >= 2) {
          beta  <- coef(fit_glm)["dosage"]
          se    <- summary(fit_glm)$coefficients["dosage", "Std. Error"]
          or    <- exp(beta)
          or_lo <- exp(beta - 1.96 * se)
          or_hi <- exp(beta + 1.96 * se)
          p_mod <- summary(fit_glm)$coefficients["dosage", "Pr(>|z|)"]
        } else {
          # Fallback 2: use Fisher OR directly (no CI from regression available)
          or    <- as.numeric(fish$estimate)
          or_lo <- NA_real_
          or_hi <- NA_real_
          p_mod <- fish$p.value
        }
      }

      # Decode the binary haplotype string back to a list of ALT-carrying rsIDs
      hap_bits <- strsplit(h, "")[[1]]
      alt_idx  <- which(hap_bits == "1")
      hap_def  <- if (length(alt_idx) == 0) {
        "Reference-like"
      } else {
        paste(region_meta$rsID_clean[alt_idx], collapse = "; ")
      }

      comp_rows[[length(comp_rows) + 1]] <- data.frame(
        Region           = region_name,
        Comparison       = comp$label,
        Haplotype_String = h,
        Haplotype_Def    = hap_def,
        Haplotype_Copies = as.integer(hap_counts[h]),
        Global_Freq      = as.integer(hap_counts[h]) / sum(hap_counts),
        Freq_Tumour      = freq_case,
        Freq_Healthy     = freq_control,
        OR               = or,
        OR_95CI_Lo       = or_lo,
        OR_95CI_Hi       = or_hi,
        P_Fisher         = fish$p.value,
        P_Model          = p_mod,
        N_Tumour         = n_case,
        N_Healthy        = n_control,
        stringsAsFactors = FALSE
      )
    }

    if (length(comp_rows) > 0) {
      results[[comp$label]] <- bind_rows(comp_rows) %>% arrange(P_Fisher)
    }
  }

  if (length(results) == 0) return(NULL)
  out <- bind_rows(results)

  # BH-FDR across all haplotype-comparison tests in this region jointly
  out$FDR                <- p.adjust(out$P_Fisher, method = "BH")
  out$Significant_Fisher <- out$P_Fisher < 0.05
  out$Significant_FDR    <- out$FDR < 0.05

  # Assign Haplotype_IDs consistently by global frequency so H1 is always the
  # most common haplotype regardless of which comparison group it was first seen in
  freq_rank <- out %>%
    distinct(Haplotype_String, Global_Freq) %>%
    arrange(desc(Global_Freq)) %>%
    mutate(Haplotype_ID = paste0("H", seq_len(n())))
  out <- out %>% left_join(freq_rank %>% select(Haplotype_String, Haplotype_ID),
                           by = "Haplotype_String")
  out
}

build_region_summary <- function(region_name, region_type, snp_meta, freq_obj, assoc_df = NULL) {
  out <- data.frame(
    Region       = region_name,
    Region_Type  = region_type,
    N_SNPs       = nrow(snp_meta),
    SNPs         = paste(snp_meta$rsID_clean, collapse = ", "),
    Genes        = paste(unique(snp_meta$Gene), collapse = ", "),
    Haps_All     = nrow(freq_obj$all),
    Haps_Shown   = nrow(freq_obj$shown),
    Top_Hap      = ifelse(nrow(freq_obj$shown) > 0, freq_obj$shown$Haplotype_ID[1], NA),
    Top_Hap_Freq = ifelse(nrow(freq_obj$shown) > 0, freq_obj$shown$Frequency[1], NA),
    N_Tests      = if (!is.null(assoc_df) && nrow(assoc_df) > 0) nrow(assoc_df) else 0L,
    stringsAsFactors = FALSE
  )
  if (!is.null(assoc_df) && nrow(assoc_df) > 0) {
    best <- assoc_df %>% arrange(P_Fisher) %>% slice(1)
    out$Best_Comparison <- best$Comparison
    out$Best_Hap_Def    <- best$Haplotype_Def
    out$Best_P_Fisher   <- best$P_Fisher
    out$Best_FDR        <- best$FDR
  }
  out
}
# =============================================================================
# FIGURE FUNCTIONS
# =============================================================================

make_ld_plot <- function(ld_matrix, snp_meta, out_file, ld_blocks = NULL) {
  # LD heatmap: r² between all pairs of SNPs, with optional block boundary
  # rectangles overlaid. Blocks are drawn as red outlines on the grid.
  ld_df <- as.data.frame(as.table(ld_matrix), stringsAsFactors = FALSE)
  names(ld_df) <- c("SNP1", "SNP2", "R2")
  pos_map <- setNames(snp_meta$POS_int, snp_meta$rsID_clean)
  ld_df$POS1 <- pos_map[ld_df$SNP1]
  ld_df$POS2 <- pos_map[ld_df$SNP2]

  snp_levels <- snp_meta$rsID_clean

  # Convert block index vectors to rectangle coordinates on the grid
  block_rects <- NULL
  if (!is.null(ld_blocks) && length(ld_blocks) > 0) {
    block_rects <- bind_rows(lapply(ld_blocks, function(idx) {
      data.frame(
        xmin = idx[1] - 0.5,
        xmax = idx[length(idx)] + 0.5,
        ymin = length(snp_levels) - idx[length(idx)] + 0.5,
        ymax = length(snp_levels) - idx[1] + 1.5
      )
    }))
  }

  p <- ggplot(ld_df, aes(x = factor(SNP1, levels = snp_levels),
                         y = factor(SNP2, levels = rev(snp_levels)),
                         fill = R2)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    { if (!is.null(block_rects))
        geom_rect(data = block_rects,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE,
                  colour = "#D6604D", fill = NA, linewidth = 0.8)
      else list() } +
    scale_fill_gradientn(
      colours  = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
      na.value = "grey95",
      limits   = c(0, 1),
      labels   = percent_format(accuracy = 1),
      name     = expression(r^2)
    ) +
    labs(title    = "Pairwise linkage disequilibrium across selected SNPs",
         subtitle = expression(paste("Pairwise ", r^2, " estimated from phased dosages | red boxes = LD blocks")),
         x = NULL, y = NULL) +
    theme(axis.text.x    = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          axis.text.y    = element_text(size = 7),
          legend.position = "right")
  SAVE(out_file, p,
       width  = max(8, nrow(snp_meta) * 0.40 + 2),
       height = max(7, nrow(snp_meta) * 0.35 + 1))
}

make_freq_plot <- function(freq_df, title, subtitle, out_file) {
  # Simple bar chart of global haplotype frequencies, labelled with percentages
  if (nrow(freq_df) == 0) return(invisible(NULL))
  p <- ggplot(freq_df, aes(x = fct_reorder(Haplotype_ID, Frequency, .desc = TRUE),
                           y = Frequency * 100)) +
    geom_col(fill = "#4393C3", colour = "black", alpha = 0.85) +
    geom_text(aes(label = sprintf("%.1f%%", Frequency * 100)), vjust = -0.35, size = 3.2) +
    labs(title = title, subtitle = subtitle, x = "Haplotype", y = "Frequency (%)") +
    coord_cartesian(ylim = c(0, max(freq_df$Frequency * 100) * 1.15))
  SAVE(out_file, p, width = max(7, nrow(freq_df) * 0.7), height = 5)
}

# -----------------------------------------------------------------------------
# Forest plot: OR + 95% CI per haplotype, faceted by comparison group
# Y-axis uses short IDs (H1, H2...); OR and p-value shown as text labels.
# CI whiskers are capped at [0.05, 20] for readability; raw values are
# preserved in the data table.
# -----------------------------------------------------------------------------
make_forest_plot <- function(assoc_df, region_name, out_file) {
  if (is.null(assoc_df) || nrow(assoc_df) == 0) return(invisible(NULL))

  plot_df <- assoc_df %>%
    filter(!is.na(OR), !is.na(OR_95CI_Lo), !is.na(OR_95CI_Hi)) %>%
    mutate(
      label      = sprintf("%s (%.0f%%)", Haplotype_ID, Global_Freq * 100),
      sig_label  = case_when(
        FDR < 0.05       ~ "FDR < 0.05",
        P_Fisher < 0.05  ~ "p < 0.05",
        TRUE             ~ "ns"
      ),
      OR_lo_plot = pmax(OR_95CI_Lo, 0.05),   # cap for display only
      OR_hi_plot = pmin(OR_95CI_Hi, 20),
      or_text    = sprintf("OR %.2f [%.2f\u2013%.2f]  p=%s",
                           OR, OR_95CI_Lo, OR_95CI_Hi,
                           ifelse(P_Fisher < 0.001,
                                  sprintf("%.2e", P_Fisher),
                                  sprintf("%.3f", P_Fisher)))
    )

  if (nrow(plot_df) == 0) return(invisible(NULL))

  x_text_pos <- 10^(log10(20) * 1.05)

  p <- ggplot(plot_df, aes(x = OR, y = fct_reorder(label, OR),
                           colour = sig_label, shape = sig_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = OR_lo_plot, xmax = OR_hi_plot),
                   height = 0.25, linewidth = 0.6) +
    geom_point(size = 3) +
    geom_text(aes(x = OR_hi_plot, label = or_text),
              hjust = -0.07, size = 2.8, colour = "grey30") +
    scale_colour_manual(values = c("FDR < 0.05" = "#D6604D",
                                   "p < 0.05"   = "#F4A582",
                                   "ns"         = "grey55"),
                        name = "Significance") +
    scale_shape_manual(values = c("FDR < 0.05" = 18,
                                  "p < 0.05"   = 17,
                                  "ns"         = 16),
                       name = "Significance") +
    scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10),
                  labels = c("0.1", "0.25", "0.5", "1", "2", "4", "10")) +
    coord_cartesian(clip = "off") +
    facet_wrap(~Comparison, ncol = 1, scales = "free_y") +
    labs(title    = sprintf("Forest plot \u2014 %s", region_name),
         subtitle = "Firth-penalised logistic regression | OR (95% CI) | log10 x-axis | label = Haplotype ID (global freq%)",
         x = "Odds ratio (log scale)", y = NULL) +
    theme(axis.text.y     = element_text(size = 8.5, face = "bold"),
          strip.text      = element_text(face = "bold"),
          legend.position = "right",
          plot.margin     = margin(6, 160, 6, 8))

  n_haps <- length(unique(plot_df$label))
  n_comp <- length(unique(plot_df$Comparison))
  SAVE(out_file, p,
       width  = 11,
       height = max(5, n_haps * 0.5 * n_comp + 2.5))
}

# -----------------------------------------------------------------------------
# Allele composition tile grid: one row per haplotype, one column per SNP.
# Red tiles = ALT allele; white = REF. Frequency labels on the right margin.
# Rows are ordered by descending global frequency (H1 at top).
# -----------------------------------------------------------------------------
make_haplotype_composition_plot <- function(freq_df, snp_meta_sub, region_name, out_file) {
  if (nrow(freq_df) == 0) return(invisible(NULL))

  snp_cols <- intersect(snp_meta_sub$rsID_clean, names(freq_df))
  if (length(snp_cols) < 2) return(invisible(NULL))

  hap_mat <- freq_df %>%
    select(Haplotype_ID, Frequency, all_of(snp_cols)) %>%
    arrange(desc(Frequency))

  tile_df <- hap_mat %>%
    pivot_longer(cols = all_of(snp_cols),
                 names_to = "SNP", values_to = "Allele_Code") %>%
    mutate(
      Allele       = ifelse(Allele_Code == 2, "Alt", "Ref"),
      SNP          = factor(SNP, levels = snp_cols),
      Haplotype_ID = factor(Haplotype_ID, levels = rev(hap_mat$Haplotype_ID)),
      Freq_label   = sprintf("%.1f%%", Frequency * 100)
    )

  gene_map  <- setNames(snp_meta_sub$Gene, snp_meta_sub$rsID_clean)
  tile_df$Gene <- gene_map[as.character(tile_df$SNP)]

  p <- ggplot(tile_df, aes(x = SNP, y = Haplotype_ID, fill = Allele)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_manual(values = c("Alt" = "#D6604D", "Ref" = "#F7F7F7"),
                      name = "Allele") +
    geom_text(data = hap_mat %>%
                mutate(Haplotype_ID = factor(Haplotype_ID, levels = rev(hap_mat$Haplotype_ID))),
              aes(x = length(snp_cols) + 0.6,
                  y = Haplotype_ID,
                  label = sprintf("%.1f%%", Frequency * 100)),
              inherit.aes = FALSE, size = 3, hjust = 0) +
    coord_cartesian(clip = "off") +
    labs(title    = sprintf("Haplotype allele composition \u2014 %s", region_name),
         subtitle = "Red = alt allele carried | rows ordered by global frequency",
         x = NULL, y = "Haplotype") +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7.5),
      axis.text.y     = element_text(size = 8),
      plot.margin     = margin(8, 50, 6, 8),
      legend.position = "top"
    )

  n_haps <- nrow(hap_mat)
  n_snps <- length(snp_cols)
  SAVE(out_file, p,
       width  = max(7, n_snps * 0.55 + 2),
       height = max(4, n_haps * 0.38 + 2))
}

# -----------------------------------------------------------------------------
# Tumour vs Healthy frequency bar chart with a haplotype definition table below.
# Short Haplotype_IDs (H1, H2...) are used on the x-axis; the definition table
# maps each ID to its constituent ALT SNPs and Fisher p-value. Combined into a
# single figure using patchwork.
# -----------------------------------------------------------------------------
make_freq_comparison_plot <- function(assoc_df, region_name, out_file,
                                      comparison = "All") {
  if (is.null(assoc_df) || nrow(assoc_df) == 0) return(invisible(NULL))

  sub <- assoc_df %>%
    filter(Comparison == comparison) %>%
    arrange(desc(Global_Freq)) %>%
    mutate(
      sig_mark  = case_when(
        FDR < 0.05      ~ "**",
        P_Fisher < 0.05 ~ "*",
        TRUE            ~ ""
      ),
      Hap_Label = Haplotype_ID
    )

  if (nrow(sub) == 0) return(invisible(NULL))

  plot_df <- sub %>%
    select(Hap_Label, Haplotype_ID, Freq_Tumour, Freq_Healthy,
           P_Fisher, FDR, sig_mark, Global_Freq) %>%
    pivot_longer(cols = c(Freq_Tumour, Freq_Healthy),
                 names_to = "Group", values_to = "Freq") %>%
    mutate(
      Group     = recode(Group, "Freq_Tumour" = "Tumour", "Freq_Healthy" = "Healthy"),
      Group     = factor(Group, levels = c("Healthy", "Tumour")),
      Hap_Label = factor(Hap_Label, levels = unique(sub$Hap_Label))
    )

  # Significance markers placed above the taller bar in each pair
  sig_df  <- plot_df %>%
    group_by(Hap_Label, sig_mark) %>%
    summarise(y_max = max(Freq, na.rm = TRUE), .groups = "drop") %>%
    filter(sig_mark != "")

  y_ceil <- max(plot_df$Freq * 100, na.rm = TRUE) * 1.18

  bar_plot <- ggplot(plot_df,
                     aes(x = Hap_Label, y = Freq * 100, fill = Group)) +
    geom_col(position = position_dodge(0.72), width = 0.62,
             colour = "black", linewidth = 0.3, alpha = 0.88) +
    { if (nrow(sig_df) > 0)
        geom_text(data = sig_df,
                  aes(x = Hap_Label, y = y_max * 100 + y_ceil * 0.04,
                      label = sig_mark),
                  inherit.aes = FALSE, size = 5, colour = "#A50026")
      else list() } +
    scale_fill_manual(values = c("Healthy" = "#4393C3", "Tumour" = "#D6604D"),
                      name = NULL) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_ceil)) +
    labs(title    = sprintf("Haplotype frequency: Tumour vs Healthy \u2014 %s (%s)",
                            region_name, comparison),
         subtitle = "* p\u202f<\u202f0.05 (Fisher)   ** FDR\u202f<\u202f0.05   bars show carrier frequency",
         x = NULL, y = "Haplotype frequency (%)") +
    theme(axis.text.x    = element_text(size = 9, face = "bold"),
          legend.position = "top",
          plot.margin    = margin(6, 8, 2, 8))

  # Definition table rendered as a ggplot grob so it can be combined with patchwork
  tbl_df <- sub %>%
    transmute(
      ID         = Haplotype_ID,
      `Global%`  = sprintf("%.1f%%", Global_Freq * 100),
      `Alt SNPs` = ifelse(nchar(Haplotype_Def) > 60,
                          paste0(substr(Haplotype_Def, 1, 57), "\u2026"),
                          Haplotype_Def),
      `p (Fisher)` = ifelse(P_Fisher < 0.001,
                             sprintf("%.2e", P_Fisher),
                             sprintf("%.3f", P_Fisher)),
      FDR_col      = ifelse(FDR < 0.001,
                            sprintf("%.2e", FDR),
                            sprintf("%.3f", FDR))
    ) %>%
    rename(FDR = FDR_col)

  tbl_plot <- ggplot() +
    theme_void() +
    annotation_custom(
      gridExtra::tableGrob(
        tbl_df, rows = NULL,
        theme = gridExtra::ttheme_minimal(
          core    = list(fg_params = list(cex = 0.72, hjust = 0, x = 0.02)),
          colhead = list(fg_params = list(cex = 0.75, fontface = "bold",
                                          hjust = 0, x = 0.02)),
          padding = unit(c(2, 4), "mm")
        )
      )
    )

  combined <- patchwork::wrap_plots(bar_plot, tbl_plot,
                                    ncol    = 1,
                                    heights = c(3, max(1, nrow(tbl_df) * 0.28)))
  n_haps <- nrow(sub)
  SAVE(out_file, combined,
       width  = max(7, n_haps * 0.85 + 2),
       height = max(6, n_haps * 0.28 + 5))
}

# -----------------------------------------------------------------------------
# Association overview bubble plot: -log10(best Fisher p) per region.
# Bubble size encodes the number of haplotypes tested in that region; colour
# encodes region type (Full / LD_Block / Gene). Two reference lines are drawn:
# the nominal p = 0.05 threshold and the Bonferroni-corrected threshold across
# the total number of haplotype-comparison tests run.
# -----------------------------------------------------------------------------
make_association_overview_plot <- function(summary_df, out_file) {
  plot_df <- summary_df %>%
    filter(!is.na(Best_P_Fisher)) %>%
    mutate(
      neg_log_p   = -log10(Best_P_Fisher),
      neg_log_fdr = -log10(pmax(Best_FDR, 1e-6)),
      Region_Type = factor(Region_Type, levels = c("Full", "LD_Block", "Gene")),
      sig_label   = case_when(
        Best_FDR < 0.05      ~ "FDR < 0.05",
        Best_P_Fisher < 0.05 ~ "p < 0.05",
        TRUE                 ~ "ns"
      ),
      Region_short = gsub("Gene_GSDMB_", "GSDMB\n", gsub("LD_", "", Region))
    )

  if (nrow(plot_df) == 0) return(invisible(NULL))

  # Bonferroni correction: divide 0.05 by the total number of individual tests
  n_tests <- if ("N_Tests" %in% names(plot_df)) sum(plot_df$N_Tests, na.rm = TRUE) else sum(plot_df$Haps_Shown)
  n_tests <- max(n_tests, 1L)

  p <- ggplot(plot_df, aes(x = fct_reorder(Region_short, neg_log_p, .desc = TRUE),
                           y = neg_log_p,
                           size   = Haps_Shown,
                           colour = Region_Type,
                           shape  = sig_label)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "grey60", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05 / n_tests),
               linetype = "dotted", colour = "#D6604D", linewidth = 0.5) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(3, 9), name = "Haplotypes\ntested") +
    scale_colour_manual(
      values = c("Full" = "#666666", "LD_Block" = "#7B3294", "Gene" = "#1B9E77"),
      name   = "Region type") +
    scale_shape_manual(
      values = c("FDR < 0.05" = 18, "p < 0.05" = 17, "ns" = 16),
      name   = "Significance") +
    annotate("text", x = Inf, y = -log10(0.05) + 0.1,
             label = "p = 0.05", hjust = 1.1, size = 3, colour = "grey50") +
    annotate("text", x = Inf, y = -log10(0.05 / n_tests) + 0.1,
             label = sprintf("Bonferroni (n=%d tests)", n_tests),
             hjust = 1.1, size = 3, colour = "#D6604D") +
    labs(title    = "Association overview across all haplotype regions",
         subtitle = "-log10(best Fisher p) per region | bubble size = haplotypes tested",
         x = NULL, y = expression(-log[10](p))) +
    theme(axis.text.x    = element_text(angle = 40, hjust = 1, size = 7.5),
          legend.position = "right")

  SAVE(out_file, p,
       width  = max(10, nrow(plot_df) * 0.65 + 2),
       height = 6)
}

make_block_ld_plot <- function(ld_matrix_full, block_idx, snp_meta_full,
                               block_name, out_file) {
  # Extract the LD sub-matrix for a single block and pass to make_ld_plot
  sub_snps <- snp_meta_full$rsID_clean[block_idx]
  sub_ld   <- ld_matrix_full[sub_snps, sub_snps, drop = FALSE]
  sub_meta <- snp_meta_full[block_idx, , drop = FALSE]
  make_ld_plot(sub_ld, sub_meta, out_file)
}

# -----------------------------------------------------------------------------
# Stacked bar chart: haplotype distribution across the four cohort-tissue groups.
# Each bar segment represents one haplotype; "Other" aggregates haplotypes below
# the frequency threshold. A consistent colour palette is used across all cohort
# frequency plots so H1 is always the same colour.
# -----------------------------------------------------------------------------
make_cohort_freq_plot <- function(a1_mat, a2_mat, freq_df, metadata_df,
                                  region_name, out_file, min_freq = 0.01,
                                  hap_palette = NULL) {
  shown_haps <- freq_df %>% filter(Frequency >= min_freq) %>% pull(Haplotype_ID)
  if (length(shown_haps) == 0) return(invisible(NULL))

  hap1 <- make_binary_hap_strings(a1_mat)
  hap2 <- make_binary_hap_strings(a2_mat)

  # Map binary allele strings back to Haplotype_IDs
  str_to_id <- setNames(freq_df$Haplotype_ID, freq_df$Allele_String)

  n    <- nrow(metadata_df)
  hid1 <- str_to_id[hap1]
  hid2 <- str_to_id[hap2]
  hid1[is.na(hid1) | !(hid1 %in% shown_haps)] <- "Other"
  hid2[is.na(hid2) | !(hid2 %in% shown_haps)] <- "Other"
  group_vec <- paste(metadata_df$Cohort, metadata_df$Tissue, sep = "\n")

  # Build long-form table vectorised (rep() avoids a row-by-row loop)
  rows <- data.frame(
    Haplotype = c(hid1, hid2),
    Group     = rep(group_vec, 2L),
    stringsAsFactors = FALSE
  )

  plot_df <- bind_rows(rows) %>%
    group_by(Group, Haplotype) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Group) %>%
    mutate(Freq = n / sum(n)) %>%
    ungroup() %>%
    mutate(
      Haplotype = factor(Haplotype,
                         levels = c(shown_haps[shown_haps %in% unique(Haplotype)], "Other")),
      Group     = factor(Group, levels = sort(unique(Group)))
    )

  n_haps <- length(shown_haps)
  if (!is.null(hap_palette) && all(shown_haps %in% names(hap_palette))) {
    hap_pal <- hap_palette[c(shown_haps, "Other")]
  } else {
    # Paired palette (max 12 colours); hue_pal() fills in any additional haplotypes
    hap_pal <- setNames(
      c(RColorBrewer::brewer.pal(min(n_haps, 12), "Paired")[seq_len(min(n_haps, 12))],
        if (n_haps > 12) scales::hue_pal()(n_haps - 12) else character(0),
        "grey75"),
      c(shown_haps, "Other")
    )
  }

  p <- ggplot(plot_df, aes(x = Group, y = Freq * 100, fill = Haplotype)) +
    geom_col(colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = hap_pal, name = "Haplotype") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 102)) +
    labs(title    = sprintf("Haplotype distribution by cohort \u2014 %s", region_name),
         subtitle = "Stacked bars show haplotype frequency within each group",
         x = NULL, y = "Frequency (%)") +
    theme(axis.text.x    = element_text(size = 9),
          legend.position = "right")

  SAVE(out_file, p, width = max(7, length(unique(plot_df$Group)) * 1.2 + 3), height = 5.5)
}

# -----------------------------------------------------------------------------
# Manhattan-style plot: -log10(Fisher p) per haplotype plotted at the median
# genomic position of its ALT-carrying SNPs. Point size encodes global
# haplotype frequency; colour encodes gene. Reference-like haplotypes (no ALT
# SNPs) are plotted at the median position of all locus SNPs.
# -----------------------------------------------------------------------------
make_manhattan_hap_plot <- function(assoc_df, snp_meta, out_file,
                                    comparison = "All") {
  if (is.null(assoc_df) || nrow(assoc_df) == 0) return(invisible(NULL))

  sub <- assoc_df %>% filter(Comparison == comparison, !is.na(P_Fisher))
  if (nrow(sub) == 0) return(invisible(NULL))

  get_median_pos <- function(hap_def, snp_meta) {
    if (hap_def == "Reference-like") return(median(snp_meta$POS_int))
    snp_ids <- trimws(strsplit(hap_def, ";")[[1]])
    pos     <- snp_meta$POS_int[snp_meta$rsID_clean %in% snp_ids]
    if (length(pos) == 0) return(median(snp_meta$POS_int))
    median(pos)
  }

  sub <- sub %>%
    rowwise() %>%
    mutate(
      med_pos    = get_median_pos(Haplotype_Def, snp_meta),
      neg_log_p  = -log10(P_Fisher),
      gene_label = {
        snp_ids <- trimws(strsplit(Haplotype_Def, ";")[[1]])
        genes   <- snp_meta$Gene[snp_meta$rsID_clean %in% snp_ids]
        if (length(genes) == 0) "Unknown" else paste(unique(genes), collapse = "/")
      }
    ) %>%
    ungroup()

  gene_lvls <- unique(snp_meta$Gene)
  gene_cols <- setNames(scales::hue_pal()(length(gene_lvls)), gene_lvls)

  p <- ggplot(sub, aes(x = med_pos / 1e6, y = neg_log_p,
                       colour = gene_label, size = Global_Freq * 100)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_point(alpha = 0.8) +
    scale_colour_manual(values = gene_cols, name = "Gene") +
    scale_size_continuous(range = c(2, 7), name = "Global freq (%)") +
    labs(title    = sprintf("Haplotype association by genomic position \u2014 %s comparison", comparison),
         subtitle = "-log10(Fisher p) | point size = global haplotype frequency | coloured by gene",
         x = "Genomic position (Mb)", y = expression(-log[10](p))) +
    theme(legend.position = "right")

  SAVE(out_file, p, width = 10, height = 5.5)
}

# -----------------------------------------------------------------------------
# Kaplan-Meier survival plot: carrier vs non-carrier of a specified haplotype.
# Requires os_months/os_event or pfs_months/pfs_event columns in the metadata,
# which are only present if clinical data were merged upstream. The function
# exits silently if survival packages or the necessary columns are absent.
# Log-rank p-value is computed manually and added as a subtitle rather than
# using survminer's built-in pval argument, to control decimal formatting.
# -----------------------------------------------------------------------------
make_km_plot <- function(a1_mat, a2_mat, freq_df, metadata_df,
                         region_name, out_dir,
                         hap_rank  = 1,
                         time_col  = "os_months", event_col = "os_event") {

  if (!requireNamespace("survival",  quietly = TRUE)) return(invisible(NULL))
  if (!requireNamespace("survminer", quietly = TRUE)) return(invisible(NULL))
  if (!all(c(time_col, event_col) %in% names(metadata_df))) return(invisible(NULL))

  top_hap_str <- freq_df$Allele_String[hap_rank]
  top_hap_id  <- freq_df$Haplotype_ID[hap_rank]
  if (is.na(top_hap_str)) return(invisible(NULL))

  hap1 <- make_binary_hap_strings(a1_mat)
  hap2 <- make_binary_hap_strings(a2_mat)

  # Carrier = sample carries at least one copy of the target haplotype
  carrier      <- as.integer(hap1 == top_hap_str) + as.integer(hap2 == top_hap_str)
  carrier_flag <- factor(ifelse(carrier >= 1, "Carrier", "Non-carrier"),
                         levels = c("Non-carrier", "Carrier"))

  surv_df <- data.frame(
    time    = as.numeric(metadata_df[[time_col]]),
    event   = as.integer(metadata_df[[event_col]]),
    carrier = carrier_flag
  ) %>% filter(!is.na(time), !is.na(event), time >= 0)

  # Require at least 10 subjects and both carrier groups to be represented
  if (nrow(surv_df) < 10 || length(unique(surv_df$carrier)) < 2) return(invisible(NULL))

  fit       <- survival::survfit(survival::Surv(time, event) ~ carrier, data = surv_df)
  pval      <- survival::survdiff(survival::Surv(time, event) ~ carrier, data = surv_df)
  pval_text <- sprintf("Log-rank p = %.3f", 1 - pchisq(pval$chisq, df = 1))

  p <- survminer::ggsurvplot(
    fit,
    data         = surv_df,
    pval         = FALSE,
    conf.int     = TRUE,
    risk.table   = TRUE,
    palette      = c("#4393C3", "#D6604D"),
    legend.labs  = c("Non-carrier", sprintf("Carrier (%s)", top_hap_id)),
    title        = sprintf("Survival by %s carrier status \u2014 %s", top_hap_id, region_name),
    subtitle     = pval_text,
    xlab         = ifelse(grepl("pfs", time_col, ignore.case = TRUE),
                          "Time (months) \u2014 PFS", "Time (months) \u2014 OS"),
    ylab         = "Survival probability",
    ggtheme      = pub_theme,
    risk.table.height = 0.25
  )

  safe_nm   <- gsub("[^A-Za-z0-9_]", "_", region_name)
  safe_time <- gsub("[^A-Za-z0-9]", "_", time_col)
  out_file  <- file.path(out_dir,
                         sprintf("19_KM_%s_%s_%s.png", safe_nm, top_hap_id, safe_time))

  png(out_file, width = 2800, height = 2400, res = 300)
  print(p)
  dev.off()
  cat("  Saved:", basename(out_file), "\n")
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================
cat("\n======================================================================\n")
cat("SCRIPT 15: HAPLOTYPE ANALYSIS\n")
cat("======================================================================\n\n")
cat(sprintf("  Main frequency threshold:       %.1f%%\n", MAIN_MIN_FREQ * 100))
cat(sprintf("  Supplementary frequency cutoff: %.1f%%\n", SUPP_MIN_FREQ * 100))
cat(sprintf("  LD block threshold (r^2):       %.2f\n\n", LD_BLOCK_THRESHOLD))

if (!requireNamespace("survminer", quietly = TRUE)) {
  cat("  Note: survminer not available \u2014 KM plots will be skipped.\n")
  cat("  To enable: conda install -c conda-forge r-survminer\n\n")
}

# --- Step 1: Load phased genotypes and annotations ---------------------------
cat("Step 1: Loading phased genotypes and annotations\n")
cat("----------------------------------------------------------------------\n")

# The whitelist defines which SNPs to include: those with gnomAD NFE AF > 1%
# from script 11. Only these SNPs enter the haplotype analysis.
target_rsids <- read_excel(WHITELIST) %>%
  pull(Variant_ID) %>%
  as.character() %>%
  trimws() %>%
  unique()

# Build a per-rsID annotation map (position, gene, alleles) from the main
# annotated report. distinct() keeps one row per rsID because VEP outputs one
# row per transcript; we only need positional and gene-level information here.
annot <- read_excel(ANNOTATED, sheet = "Biological_Annotations")
annotation_map <- annot %>%
  mutate(
    rsID_clean = str_extract(Existing_variation, "rs[0-9]+"),
    POS_int    = as.integer(gsub("[^0-9]", "", as.character(POS))),
    Gene       = ifelse(is.na(SYMBOL) | trimws(SYMBOL) == "", "Unknown", trimws(SYMBOL))
  ) %>%
  filter(rsID_clean %in% target_rsids) %>%
  select(POS_int, rsID_clean, Gene, REF, ALT) %>%
  distinct(rsID_clean, .keep_all = TRUE) %>%
  arrange(POS_int)

geno_raw <- read.table(GENO_FILE, header = TRUE, sep = "\t", check.names = FALSE, quote = "")
meta_raw <- read.table(META_FILE, header = TRUE, sep = "\t", check.names = FALSE, quote = "")

stopifnot(all(c("Sample", "Cohort", "Tissue") %in% names(meta_raw)))

# The VCF ID field from BEAGLE is sometimes formatted as "CHROM:POS:REF:ALT"
# rather than an rsID. Extract the numeric POS component to join against the
# annotation map by genomic position.
geno_raw <- geno_raw %>%
  mutate(extracted_pos = as.integer(str_split_fixed(ID, ":", 4)[, 2]))

# inner_join on POS_int retains only SNPs present in both the phased VCF and
# the whitelist annotation map. Where the join produces duplicate REF/ALT
# columns (.x from VCF, .y from annotation_map), the curated annotation values
# (.y) are kept because they have been verified against the reference genome.
geno_filtered <- geno_raw %>%
  inner_join(annotation_map, by = c("extracted_pos" = "POS_int")) %>%
  mutate(POS_int = extracted_pos) %>%
  { if ("REF.y" %in% names(.)) rename(., REF = REF.y, ALT = ALT.y) else . } %>%
  select(-any_of(c("REF.x", "ALT.x"))) %>%
  arrange(POS_int) %>%
  distinct(rsID_clean, .keep_all = TRUE)   # guard against many-to-many collisions at shared positions

# Align metadata rows to the sample column order in the genotype table;
# the strict equality check below ensures no mismatch propagates silently
sample_cols <- intersect(names(geno_filtered), meta_raw$Sample)
meta <- meta_raw %>%
  filter(Sample %in% sample_cols) %>%
  mutate(Group = paste(Cohort, Tissue, sep = "_"))

sample_cols <- sample_cols[sample_cols %in% meta$Sample]
meta <- meta[match(sample_cols, meta$Sample), , drop = FALSE]
stopifnot(all(sample_cols == meta$Sample))

cat(sprintf("  Genotype table: %d SNPs x %d samples\n", nrow(geno_filtered), length(sample_cols)))
cat(sprintf("  Metadata:       %d samples retained\n", nrow(meta)))
cat(sprintf("  Genes covered:  %s\n", paste(unique(geno_filtered$Gene), collapse = ", ")))

# --- Step 2: Build phased allele matrices ------------------------------------
cat("\nStep 2: Building phased allele matrices\n")
cat("----------------------------------------------------------------------\n")

am <- build_allele_matrices(geno_filtered, sample_cols)
a1 <- am$a1
a2 <- am$a2

# snp_meta is the per-SNP metadata table used throughout the analysis;
# it is kept in genomic position order to ensure consistent matrix indexing
snp_meta <- geno_filtered %>%
  select(POS_int, rsID_clean, Gene, REF, ALT) %>%
  distinct(rsID_clean, .keep_all = TRUE) %>%
  arrange(POS_int)

n_samples     <- nrow(a1)
n_snps        <- ncol(a1)
n_missing     <- sum(is.na(a1) | is.na(a2))
total_alleles <- 2 * n_samples * n_snps

cat(sprintf("  Allele matrices: %d samples x %d SNPs\n", n_samples, n_snps))
cat(sprintf("  Missing alleles: %d / %d (%.2f%%)\n",
            n_missing, total_alleles, 100 * n_missing / total_alleles))

# --- Step 3: Full-region haplotype analysis ----------------------------------
cat("\nStep 3: Full-region haplotypes\n")
cat("----------------------------------------------------------------------\n")

full_freq  <- make_em_freq_table(a1, a2, snp_meta, min_freq = MAIN_MIN_FREQ)
full_assoc <- run_assoc_for_region("Full_Region", a1, a2, snp_meta, meta,
                                   min_copies = MIN_HAP_COPIES)

# Build a single consistent haplotype colour palette from the full-region
# frequency table so that H1, H2... are always the same colour across all
# cohort-frequency plots regardless of which region is being plotted
.all_hap_ids <- full_freq$shown$Haplotype_ID
.n_haps      <- length(.all_hap_ids)
.hap_palette <- if (.n_haps > 0) {
  setNames(
    c(RColorBrewer::brewer.pal(min(.n_haps, 12), "Paired")[seq_len(min(.n_haps, 12))],
      if (.n_haps > 12) scales::hue_pal()(.n_haps - 12) else character(0)),
    .all_hap_ids
  )
} else {
  character(0)
}
.hap_palette["Other"] <- "grey75"

cat(sprintf("  Full-region: %d haplotypes total | %d shown at >= %.1f%%\n",
            nrow(full_freq$all), nrow(full_freq$shown), MAIN_MIN_FREQ * 100))
if (nrow(full_freq$shown) > 0) {
  cat(sprintf("  Top haplotype: %s (%.1f%%)\n",
              full_freq$shown$Haplotype_ID[1], full_freq$shown$Frequency[1] * 100))
}

# --- Step 4: LD analysis and block definition --------------------------------
cat("\nStep 4: Linkage disequilibrium and LD blocks\n")
cat("----------------------------------------------------------------------\n")

ld_matrix <- compute_ld_matrix(a1, a2, snp_meta)
ld_blocks  <- if (RUN_BLOCK_ANALYSIS) {
  define_ld_blocks(ld_matrix, snp_meta,
                   threshold = LD_BLOCK_THRESHOLD,
                   min_snps  = MIN_BLOCK_SNPS,
                   max_snps  = MAX_BLOCK_SNPS)
} else {
  list()
}

if (length(ld_blocks) == 0) {
  cat("  No LD blocks met the configured criteria\n")
} else {
  cat(sprintf("  LD blocks defined: %d\n", length(ld_blocks)))
  for (nm in names(ld_blocks)) {
    idx <- ld_blocks[[nm]]
    cat(sprintf("  %s: %d SNPs | %s\n",
                nm, length(idx), paste(snp_meta$rsID_clean[idx], collapse = ", ")))
  }
}

block_freq_tables  <- list()
block_assoc_tables <- list()
block_summary_rows <- list()

if (length(ld_blocks) > 0) {
  for (nm in names(ld_blocks)) {
    idx         <- ld_blocks[[nm]]
    region_meta <- snp_meta[idx, , drop = FALSE]
    region_name <- sprintf("LD_%s", nm)

    freq_obj <- make_em_freq_table(a1[, idx, drop = FALSE], a2[, idx, drop = FALSE],
                                   region_meta, min_freq = SUPP_MIN_FREQ)
    assoc_df <- run_assoc_for_region(region_name,
                                     a1[, idx, drop = FALSE], a2[, idx, drop = FALSE],
                                     region_meta, meta, min_copies = MIN_HAP_COPIES)

    block_freq_tables[[region_name]]  <- freq_obj$shown
    block_assoc_tables[[region_name]] <- assoc_df
    block_summary_rows[[region_name]] <- build_region_summary(region_name, "LD_Block",
                                                              region_meta, freq_obj, assoc_df)
  }
}

# --- Step 5: Gene-specific haplotypes ----------------------------------------
cat("\nStep 5: Gene-specific haplotypes\n")
cat("----------------------------------------------------------------------\n")

gene_freq_tables  <- list()
gene_assoc_tables <- list()
gene_summary_rows <- list()

if (RUN_GENE_SPECIFIC) {
  genes_to_run <- intersect(GENE_LIST_SUPP, unique(snp_meta$Gene))
  if (length(genes_to_run) == 0) {
    cat("  None of the requested genes were found in the SNP panel\n")
  } else {
    for (gene in genes_to_run) {
      gene_idx <- which(snp_meta$Gene == gene)
      if (length(gene_idx) < 2) {
        cat(sprintf("  Skipping %s: fewer than 2 SNPs\n", gene))
        next
      }

      # Optionally sub-partition gene SNPs into LD blocks before phasing;
      # this reduces the number of SNPs per EM call and improves power
      # where clear block structure exists within the gene
      if (GENE_SPECIFIC_USE_LD_BLOCKS) {
        gene_ld     <- ld_matrix[gene_idx, gene_idx, drop = FALSE]
        gene_blocks <- define_ld_blocks(gene_ld, snp_meta[gene_idx, , drop = FALSE],
                                        threshold = LD_BLOCK_THRESHOLD,
                                        min_snps  = MIN_BLOCK_SNPS,
                                        max_snps  = MAX_BLOCK_SNPS)
        if (length(gene_blocks) == 0) gene_blocks <- list(All = seq_along(gene_idx))
      } else {
        gene_blocks <- list(All = seq_along(gene_idx))
      }

      for (sub_nm in names(gene_blocks)) {
        sub_local_idx <- gene_blocks[[sub_nm]]
        sub_idx       <- gene_idx[sub_local_idx]
        region_meta   <- snp_meta[sub_idx, , drop = FALSE]
        region_name   <- if (sub_nm == "All") sprintf("Gene_%s", gene) else
                                              sprintf("Gene_%s_%s", gene, sub_nm)

        freq_obj <- make_em_freq_table(a1[, sub_idx, drop = FALSE],
                                       a2[, sub_idx, drop = FALSE],
                                       region_meta, min_freq = SUPP_MIN_FREQ)
        assoc_df <- run_assoc_for_region(region_name,
                                         a1[, sub_idx, drop = FALSE],
                                         a2[, sub_idx, drop = FALSE],
                                         region_meta, meta, min_copies = MIN_HAP_COPIES)

        gene_freq_tables[[region_name]]  <- freq_obj$shown
        gene_assoc_tables[[region_name]] <- assoc_df
        gene_summary_rows[[region_name]] <- build_region_summary(region_name, "Gene",
                                                                 region_meta, freq_obj, assoc_df)
        cat(sprintf("  %s: %d SNPs | %d shown haplotypes\n",
                    region_name, nrow(region_meta), nrow(freq_obj$shown)))
      }
    }
  }
}

# --- Step 6: Save results workbook -------------------------------------------
cat("\nStep 6: Saving results workbook\n")
cat("----------------------------------------------------------------------\n")

wb <- createWorkbook()

addWorksheet(wb, "SNPs_Used");       writeData(wb, "SNPs_Used",       snp_meta)
addWorksheet(wb, "Full_Global_Freqs"); writeData(wb, "Full_Global_Freqs", full_freq$shown)

if (!is.null(full_assoc) && nrow(full_assoc) > 0) {
  addWorksheet(wb, "Full_Associations")
  writeData(wb, "Full_Associations", full_assoc)
}

# LD matrix in long format for easier downstream handling in R or Python
ld_long <- as.data.frame(as.table(ld_matrix), stringsAsFactors = FALSE)
names(ld_long) <- c("SNP1", "SNP2", "R2")
addWorksheet(wb, "LD_Matrix_Long")
writeData(wb, "LD_Matrix_Long", ld_long)

if (length(ld_blocks) > 0) {
  ld_block_df <- bind_rows(lapply(names(ld_blocks), function(nm) {
    idx <- ld_blocks[[nm]]
    data.frame(Block = nm, SNP_Order = seq_along(idx),
               rsID = snp_meta$rsID_clean[idx], POS = snp_meta$POS_int[idx],
               Gene = snp_meta$Gene[idx], stringsAsFactors = FALSE)
  }))
  addWorksheet(wb, "LD_Blocks")
  writeData(wb, "LD_Blocks", ld_block_df)
}

for (nm in names(block_freq_tables)) {
  addWorksheet(wb, safe_sheet(paste0(nm, "_Freqs")))
  writeData(wb, safe_sheet(paste0(nm, "_Freqs")), block_freq_tables[[nm]])
}
for (nm in names(block_assoc_tables)) {
  if (!is.null(block_assoc_tables[[nm]]) && nrow(block_assoc_tables[[nm]]) > 0) {
    addWorksheet(wb, safe_sheet(paste0(nm, "_Assoc")))
    writeData(wb, safe_sheet(paste0(nm, "_Assoc")), block_assoc_tables[[nm]])
  }
}
for (nm in names(gene_freq_tables)) {
  addWorksheet(wb, safe_sheet(paste0(nm, "_Freqs")))
  writeData(wb, safe_sheet(paste0(nm, "_Freqs")), gene_freq_tables[[nm]])
}
for (nm in names(gene_assoc_tables)) {
  if (!is.null(gene_assoc_tables[[nm]]) && nrow(gene_assoc_tables[[nm]]) > 0) {
    addWorksheet(wb, safe_sheet(paste0(nm, "_Assoc")))
    writeData(wb, safe_sheet(paste0(nm, "_Assoc")), gene_assoc_tables[[nm]])
  }
}

summary_parts <- list(build_region_summary("Full_Region", "Full", snp_meta, full_freq, full_assoc))
if (length(block_summary_rows) > 0) summary_parts[[length(summary_parts) + 1]] <- bind_rows(block_summary_rows)
if (length(gene_summary_rows)  > 0) summary_parts[[length(summary_parts) + 1]] <- bind_rows(gene_summary_rows)
summary_rows <- bind_rows(summary_parts)
addWorksheet(wb, "Region_Summary")
writeData(wb, "Region_Summary", summary_rows)

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("  Workbook saved: %s\n", OUT_XLSX))

# --- Step 7: Generate figures ------------------------------------------------
cat("\nStep 7: Generating figures\n")
cat("----------------------------------------------------------------------\n")

make_freq_plot(
  full_freq$shown,
  title    = "Full-region global haplotype frequencies",
  subtitle = sprintf("%d SNPs across %d samples | shown haplotypes >= %.1f%%",
                     n_snps, n_samples, MAIN_MIN_FREQ * 100),
  out_file = file.path(OUT_DIR, "19_FullRegion_Haplotype_Frequencies.png")
)

make_ld_plot(ld_matrix, snp_meta, ld_blocks = ld_blocks,
             out_file = file.path(OUT_DIR, "19_LD_Heatmap_FullRegion.png"))

if (length(block_summary_rows) > 0) {
  block_plot_df <- bind_rows(block_summary_rows) %>%
    mutate(Region = factor(Region, levels = Region[order(N_SNPs, decreasing = TRUE)]))

  p_blocks <- ggplot(block_plot_df, aes(x = Region, y = Haps_Shown)) +
    geom_col(fill = "#7B3294", colour = "black", alpha = 0.85) +
    geom_text(aes(label = paste0(N_SNPs, " SNPs")), vjust = -0.35, size = 3.2) +
    labs(title    = "LD-block haplotype complexity",
         subtitle = sprintf("Shown haplotypes at >= %.1f%% frequency", SUPP_MIN_FREQ * 100),
         x = NULL, y = "Shown haplotypes") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  SAVE(file.path(OUT_DIR, "19_LD_Block_Haplotype_Counts.png"), p_blocks, width = 8, height = 5)
}

if (length(gene_summary_rows) > 0) {
  gene_plot_df <- bind_rows(gene_summary_rows) %>%
    mutate(Region = factor(Region, levels = Region[order(Haps_Shown, decreasing = TRUE)]))

  p_genes <- ggplot(gene_plot_df, aes(x = Region, y = Haps_Shown)) +
    geom_col(fill = "#1B9E77", colour = "black", alpha = 0.85) +
    geom_text(aes(label = paste0(N_SNPs, " SNPs")), vjust = -0.35, size = 3.2) +
    labs(title    = "Gene-specific haplotype complexity",
         subtitle = sprintf("Shown haplotypes at >= %.1f%% frequency", SUPP_MIN_FREQ * 100),
         x = NULL, y = "Shown haplotypes") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  SAVE(file.path(OUT_DIR, "19_Gene_Haplotype_Counts.png"), p_genes, width = 8, height = 5)
}

make_association_overview_plot(summary_rows,
                               out_file = file.path(OUT_DIR, "19_Association_Overview_Bubble.png"))

if (!is.null(full_assoc) && nrow(full_assoc) > 0) {
  make_forest_plot(full_assoc, region_name = "Full_Region",
                   out_file = file.path(OUT_DIR, "19_Forest_FullRegion.png"))
}

make_haplotype_composition_plot(freq_df = full_freq$shown, snp_meta_sub = snp_meta,
                                region_name = "Full_Region",
                                out_file = file.path(OUT_DIR, "19_Composition_FullRegion.png"))

if (!is.null(full_assoc) && nrow(full_assoc) > 0) {
  make_freq_comparison_plot(full_assoc, region_name = "Full_Region",
                            out_file   = file.path(OUT_DIR, "19_FreqComparison_FullRegion.png"),
                            comparison = "All")
}

if (!is.null(full_assoc) && nrow(full_assoc) > 0) {
  make_manhattan_hap_plot(full_assoc, snp_meta,
                          out_file   = file.path(OUT_DIR, "19_Manhattan_Haplotypes_All.png"),
                          comparison = "All")
  for (.comp in c("Breast", "Endometrium")) {
    if (.comp %in% full_assoc$Comparison) {
      make_manhattan_hap_plot(full_assoc, snp_meta,
                              out_file   = file.path(OUT_DIR, sprintf("19_Manhattan_Haplotypes_%s.png", .comp)),
                              comparison = .comp)
    }
  }
}

make_cohort_freq_plot(a1_mat = a1, a2_mat = a2, freq_df = full_freq$shown,
                      metadata_df = meta, region_name = "Full_Region",
                      out_file    = file.path(OUT_DIR, "19_CohortFreq_FullRegion.png"),
                      min_freq    = MAIN_MIN_FREQ, hap_palette = .hap_palette)

# KM plots: run for the most frequent haplotype plus any nominally significant ones
if (nrow(full_freq$shown) > 0) {
  km_hap_ranks <- unique(c(
    1L,
    if (!is.null(full_assoc) && nrow(full_assoc) > 0)
      which(full_freq$shown$Allele_String %in%
              full_assoc$Haplotype_String[full_assoc$Significant_Fisher])
    else integer(0)
  ))
  for (.tc in c("os_months", "pfs_months")) {
    .ec <- sub("months", "event", .tc)
    for (.rank in km_hap_ranks) {
      make_km_plot(a1, a2, full_freq$shown, meta,
                   region_name = "Full_Region", out_dir = OUT_DIR,
                   hap_rank = .rank, time_col = .tc, event_col = .ec)
    }
  }
}

# Per-gene KM plots
if (length(gene_freq_tables) > 0) {
  for (.gnm in names(gene_freq_tables)) {
    .gfreq <- gene_freq_tables[[.gnm]]
    if (is.null(.gfreq) || nrow(.gfreq) == 0) next
    .gsnp_ids <- trimws(strsplit(gene_summary_rows[[.gnm]]$SNPs, ",")[[1]])
    .gidx     <- which(snp_meta$rsID_clean %in% .gsnp_ids)
    if (length(.gidx) < 2) next
    .g_assoc       <- gene_assoc_tables[[.gnm]]
    km_hap_ranks_g <- unique(c(
      1L,
      if (!is.null(.g_assoc) && nrow(.g_assoc) > 0)
        which(.gfreq$Allele_String %in%
                .g_assoc$Haplotype_String[.g_assoc$Significant_Fisher])
      else integer(0)
    ))
    for (.tc in c("os_months", "pfs_months")) {
      .ec <- sub("months", "event", .tc)
      for (.rank in km_hap_ranks_g) {
        make_km_plot(a1[, .gidx, drop = FALSE], a2[, .gidx, drop = FALSE],
                     .gfreq, meta, region_name = .gnm, out_dir = OUT_DIR,
                     hap_rank = .rank, time_col = .tc, event_col = .ec)
      }
    }
  }
}

# Per-gene-block figures: forest, composition, frequency comparison, LD heatmaps
if (length(gene_assoc_tables) > 0) {
  for (region_nm in names(gene_assoc_tables)) {
    assoc_df_g <- gene_assoc_tables[[region_nm]]
    freq_df_g  <- gene_freq_tables[[region_nm]]
    safe_nm    <- gsub("[^A-Za-z0-9_]", "_", region_nm)

    if (!is.null(assoc_df_g) && nrow(assoc_df_g) > 0) {
      make_forest_plot(assoc_df_g, region_name = region_nm,
                       out_file = file.path(OUT_DIR, sprintf("19_Forest_%s.png", safe_nm)))
    }

    if (!is.null(freq_df_g) && nrow(freq_df_g) > 0) {
      gene_block_meta <- tryCatch({
        snp_ids <- trimws(strsplit(gene_summary_rows[[region_nm]]$SNPs, ",")[[1]])
        snp_meta[snp_meta$rsID_clean %in% snp_ids, , drop = FALSE]
      }, error = function(e) NULL)
      if (!is.null(gene_block_meta) && nrow(gene_block_meta) > 0) {
        make_haplotype_composition_plot(freq_df = freq_df_g, snp_meta_sub = gene_block_meta,
                                        region_name = region_nm,
                                        out_file = file.path(OUT_DIR, sprintf("19_Composition_%s.png", safe_nm)))
      }
    }

    if (!is.null(assoc_df_g) && nrow(assoc_df_g) > 0) {
      for (.comp in unique(assoc_df_g$Comparison)) {
        .comp_safe <- gsub("[^A-Za-z0-9]", "_", .comp)
        make_freq_comparison_plot(assoc_df_g, region_name = region_nm,
                                  out_file   = file.path(OUT_DIR, sprintf("19_FreqComparison_%s_%s.png", safe_nm, .comp_safe)),
                                  comparison = .comp)
        make_forest_plot(assoc_df_g %>% filter(Comparison == .comp),
                         region_name = sprintf("%s (%s)", region_nm, .comp),
                         out_file    = file.path(OUT_DIR, sprintf("19_Forest_%s_%s.png", safe_nm, .comp_safe)))
      }
    }
  }
}

# Per-GSDMB-block LD heatmaps
if (RUN_GENE_SPECIFIC && exists("ld_matrix")) {
  gsdmb_idx <- which(snp_meta$Gene == "GSDMB")
  if (length(gsdmb_idx) >= 2) {
    gsdmb_ld <- ld_matrix[snp_meta$rsID_clean[gsdmb_idx],
                           snp_meta$rsID_clean[gsdmb_idx], drop = FALSE]
    gsdmb_meta         <- snp_meta[gsdmb_idx, , drop = FALSE]
    gsdmb_blocks_local <- define_ld_blocks(gsdmb_ld, gsdmb_meta,
                                           threshold = LD_BLOCK_THRESHOLD,
                                           min_snps  = MIN_BLOCK_SNPS,
                                           max_snps  = MAX_BLOCK_SNPS)
    for (blk_nm in names(gsdmb_blocks_local)) {
      blk_idx <- gsdmb_idx[gsdmb_blocks_local[[blk_nm]]]
      make_block_ld_plot(ld_matrix_full = ld_matrix, block_idx = blk_idx,
                         snp_meta_full  = snp_meta, block_name = blk_nm,
                         out_file = file.path(OUT_DIR, sprintf("19_LD_GSDMB_%s.png", blk_nm)))
    }
  }
}

# --- Completion summary ------------------------------------------------------
cat("\n======================================================================\n")
cat("HAPLOTYPE ANALYSIS COMPLETE\n")
cat("======================================================================\n\n")
cat(sprintf("Outputs in: %s\n\n", OUT_DIR))

out_files <- sort(list.files(OUT_DIR, full.names = FALSE))
cat(sprintf("Files created (%d total):\n", length(out_files)))
for (f in out_files) cat(sprintf("  %s\n", f))

cat("\nRegion summary:\n")
print(summary_rows, row.names = FALSE)
cat("\n")
