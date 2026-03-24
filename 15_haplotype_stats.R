#!/usr/bin/env Rscript
# =============================================================================
# Script 15: GSDMB Haplotype Analysis
# Full-region + LD-block + gene-specific haplotypes
# =============================================================================
# Purpose:
#   Perform downstream statistical analysis of phased haplotypes, including
#   frequency estimation, tumour-versus-control comparison and clinical
#   association modelling for the GSDMB locus.
#
# Inputs:
#   - Phased genotype matrices derived from the BEAGLE workflow.
#   - Harmonised clinical metadata and SNP backbone definitions.
#
# Outputs:
#   - Publication-style figures, Excel summaries and modelling tables that
#     describe haplotype prevalence and association signals.
# =============================================================================

suppressPackageStartupMessages({
  required_pkgs <- c(
    "haplo.stats", "ggplot2", "dplyr", "tidyr", "openxlsx", "readxl",
    "stringr", "scales", "forcats", "patchwork", "RColorBrewer",
    "gridExtra"
  )
  missing_required <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_required) > 0) {
    stop(
      paste0(
        "Missing required R packages: ",
        paste(missing_required, collapse = ", "),
        ". Run Rscript /home/gadeaalonsoj/tfm/install_haplotype_r_dependencies.R before rerunning stage 15R."
      )
    )
  }

  library(haplo.stats)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(readxl)
  library(stringr)
  library(scales)
  library(forcats)
  library(patchwork)

  HAS_LOGISTF <- requireNamespace("logistf", quietly = TRUE)
  HAS_SURVIVAL <- requireNamespace("survival", quietly = TRUE)
  HAS_SURVMINER <- HAS_SURVIVAL && requireNamespace("survminer", quietly = TRUE)
})

# Explicitly prefer dplyr verbs
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

# =============================================================================
# 0. GLOBAL THEME
# =============================================================================
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

GROUP_COLS <- c(
  "Breast_Healthy"      = "#F4CAE4",
  "Breast_Tumour"       = "#CC79A7",
  "Endometrium_Healthy" = "#B3E2CD",
  "Endometrium_Tumour"  = "#009E73",
  "All_Healthy"         = "#56B4E9",
  "All_Tumour"          = "#D55E00"
)

SAVE <- function(file, plot, width = 10, height = 6) {
  ggsave(file, plot, width = width, height = height, dpi = 300,
         units = "in", bg = "white")
  cat("  Saved:", basename(file), "\n")
}

safe_sheet <- function(x) {
  x <- gsub("[:/\\?*\\[\\]]", "_", x)
  substr(x, 1, 31)
}

format_threshold_label <- function(threshold) {
  sprintf("r2_%03d", round(threshold * 100))
}

# =============================================================================
# 1. CONFIGURATION
# =============================================================================
BASE_DIR  <- "/home/gadeaalonsoj/tfm/analysis_results"
HAPLO_DIR <- file.path(BASE_DIR, "15_haplotype_phasing")
OUT_DIR   <- file.path(BASE_DIR, "15_haplotype_statistics")

GENO_FILE <- file.path(HAPLO_DIR, "phased_genotypes.tsv")
META_FILE <- file.path(HAPLO_DIR, "sample_metadata.tsv")
WHITELIST <- file.path(BASE_DIR, "11_common_snps", "GSDMB_Common_SNP_Frequency_Summary.xlsx")
ANNOTATED <- file.path(BASE_DIR, "07_annotated_variants", "GSDMB_Annotated_Variants.xlsx")

OUT_XLSX  <- file.path(OUT_DIR, "GSDMB_Haplotype_Results.xlsx")
MAIN_FIG_DIR <- file.path(OUT_DIR, "main_figures")
SUPP_FIG_DIR <- file.path(OUT_DIR, "supplementary_figures")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(MAIN_FIG_DIR)) dir.create(MAIN_FIG_DIR, recursive = TRUE)
if (!dir.exists(SUPP_FIG_DIR)) dir.create(SUPP_FIG_DIR, recursive = TRUE)

# Main filters
MAIN_MIN_FREQ        <- 0.01   # full-region global haplotypes shown from 1%
SUPP_MIN_FREQ        <- 0.005  # LD-block and gene-specific analyses can be looser
MIN_HAP_COPIES       <- 3
MIN_SAMPLES_PER_ARM  <- 3
MIN_GROUP_SAMPLES_EM <- 3

# LD block settings
RUN_BLOCK_ANALYSIS   <- TRUE
LD_BLOCK_THRESHOLDS  <- c(0.60, 0.80)
MIN_BLOCK_SNPS       <- 2
MAX_BLOCK_SNPS       <- 12

# Gene-specific settings
RUN_GENE_SPECIFIC    <- TRUE
GENE_LIST_SUPP       <- c("GSDMB")
GENE_SPECIFIC_USE_LD_BLOCKS <- TRUE
GENERATE_SUPPLEMENTARY_FIGURES <- FALSE
MAIN_TOP_FULL_HAPS <- 10
MAIN_TOP_COMPARISON_HAPS <- 8
MAX_KEY_BLOCKS_MAIN <- 3

# =============================================================================
# 2. HELPERS
# =============================================================================


parse_phased_call <- function(x) {
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
  n_samples <- length(sample_cols)
  n_snps    <- nrow(geno_df)

  # Extract the raw genotype character matrix once (samples × SNPs after transpose)
  raw_mat <- t(as.matrix(geno_df[, sample_cols, drop = FALSE]))  # n_samples × n_snps

  # Vectorised parse: split on | or /, extract first and second allele
  parse_allele <- function(x, allele_idx) {
    x  <- trimws(as.character(x))
    bad <- is.na(x) | x == "" | x == "." | x == "./." | x == ".|."
    # fast split: use sub to extract each haplotype
    a <- if (allele_idx == 1L) {
      sub("([0-9]+)[|/]([0-9]+)", "\\1", x)
    } else {
      sub("([0-9]+)[|/]([0-9]+)", "\\2", x)
    }
    out <- suppressWarnings(as.integer(a)) + 1L
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
  stopifnot(nrow(a1_mat) == nrow(a2_mat), ncol(a1_mat) == ncol(a2_mat))
  out <- vector("list", ncol(a1_mat) * 2)
  for (j in seq_len(ncol(a1_mat))) {
    out[[2 * j - 1]] <- a1_mat[, j]
    out[[2 * j]]     <- a2_mat[, j]
  }
  as.data.frame(out, check.names = FALSE)
}

make_binary_hap_strings <- function(a_mat) {
  apply(a_mat - 1L, 1, paste, collapse = "")
}

make_em_freq_table <- function(a1_mat, a2_mat, snp_meta, min_freq) {
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
  freqs$Count     <- round(freqs$Frequency * nrow(a1_mat) * 2)
  freqs <- freqs %>%
    arrange(desc(Frequency)) %>%
    mutate(Haplotype_ID = paste0("H", seq_len(n())))

  snp_matrix <- matrix(
    as.numeric(as.matrix(freqs[, snp_meta$rsID_clean, drop = FALSE])),
    nrow = nrow(freqs), ncol = nrow(snp_meta)
  )
  freqs$Allele_String <- apply(snp_matrix, 1, function(x) paste0(x - 1L, collapse = ""))
  freqs$Alt_Alleles   <- rowSums(snp_matrix - 1L)

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
  dosage <- (a1_mat - 1L) + (a2_mat - 1L)
  n_snps <- ncol(dosage)
  ld <- matrix(NA_real_, nrow = n_snps, ncol = n_snps,
               dimnames = list(snp_meta$rsID_clean, snp_meta$rsID_clean))
  diag(ld) <- 1
  for (i in seq_len(n_snps)) {
    for (j in i:n_snps) {
      xi <- dosage[, i]
      xj <- dosage[, j]
      keep <- complete.cases(xi, xj)
      if (sum(keep) < 3) {
        r2 <- NA_real_
      } else {
        r <- suppressWarnings(cor(xi[keep], xj[keep]))
        r2 <- ifelse(is.na(r), NA_real_, r^2)
      }
      ld[i, j] <- r2
      ld[j, i] <- r2
    }
  }
  ld
}

define_ld_blocks <- function(ld_matrix, snp_meta, threshold = 0.6,
                             min_snps = 2, max_snps = 12) {
  # Gabriel-style block definition: a candidate window i:j forms a block only if
  # ALL pairwise r2 within it are >= threshold.
  # Uses rownames of ld_matrix for indexing so this is safe on any sub-matrix.
  snp_names <- rownames(ld_matrix)
  n         <- length(snp_names)
  if (is.null(snp_names) || n < min_snps) return(list())

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
    # Need at least min_snps remaining SNPs to form any block
    if (i + min_snps - 1L > n) break

    best_end <- NA
    max_j    <- min(n, i + max_snps - 1L)

    # Scan from largest window down to smallest; take first (longest) valid block
    for (j in max_j:(i + min_snps - 1L)) {
      if (block_coherent(i:j)) {
        best_end <- j
        break
      }
    }

    if (!is.na(best_end)) {
      blocks[[paste0("Block", block_id)]] <- i:best_end
      block_id <- block_id + 1
      i <- best_end + 1
    } else {
      i <- i + 1
    }
  }
  blocks
}

run_assoc_for_region <- function(region_name, a1_mat, a2_mat, region_meta,
                                 metadata_df, min_copies = 3,
                                 comparison_set = c("Breast", "Endometrium", "All")) {
  hap1 <- make_binary_hap_strings(a1_mat)
  hap2 <- make_binary_hap_strings(a2_mat)
  hap_counts <- sort(table(c(hap1, hap2)), decreasing = TRUE)
  tested_haps <- names(hap_counts[hap_counts >= min_copies])

  if (length(tested_haps) == 0) return(NULL)

  results <- list()

  comp_defs <- list(
    list(label = "Breast",      keep = (metadata_df$Cohort == "Breast" & metadata_df$Tissue == "Tumour") | metadata_df$Tissue == "Healthy"),
    list(label = "Endometrium", keep = (metadata_df$Cohort == "Endometrium" & metadata_df$Tissue == "Tumour") | metadata_df$Tissue == "Healthy"),
    list(label = "All",         keep = rep(TRUE, nrow(metadata_df)))
  )
  comp_defs <- comp_defs[sapply(comp_defs, function(x) x$label %in% comparison_set)]

  for (comp in comp_defs) {
    idx <- which(comp$keep & metadata_df$Tissue %in% c("Healthy", "Tumour"))
    if (length(idx) == 0) next

    meta_sub <- metadata_df[idx, , drop = FALSE]
    y <- as.integer(meta_sub$Tissue == "Tumour")
    n_case <- sum(y == 1)
    n_control <- sum(y == 0)
    if (n_case < MIN_SAMPLES_PER_ARM || n_control < MIN_SAMPLES_PER_ARM) next

    hap1_sub <- hap1[idx]
    hap2_sub <- hap2[idx]

    comp_rows <- list()
    for (h in tested_haps) {
      dosage <- as.integer(hap1_sub == h) + as.integer(hap2_sub == h)
      if (var(dosage) == 0) next

      carrier <- as.integer(dosage >= 1)
      carrier_total <- sum(carrier == 1, na.rm = TRUE)
      noncarrier_total <- sum(carrier == 0, na.rm = TRUE)
      if (carrier_total < min_copies || noncarrier_total < MIN_SAMPLES_PER_ARM) next
      tbl <- table(carrier, y)
      if (nrow(tbl) < 2 || ncol(tbl) < 2) next

      fish <- fisher.test(tbl)
      freq_case <- mean(dosage[y == 1], na.rm = TRUE) / 2
      freq_control <- mean(dosage[y == 0], na.rm = TRUE) / 2

      fit <- if (HAS_LOGISTF) {
        tryCatch(
          logistf::logistf(y ~ dosage),
          error = function(e) NULL
        )
      } else {
        NULL
      }

      if (!is.null(fit) && "dosage" %in% names(coef(fit))) {
        beta  <- coef(fit)["dosage"]
        ci    <- as.numeric(suppressWarnings(confint(fit, parm = "dosage")))
        or    <- exp(beta)
        or_lo <- exp(ci[1])
        or_hi <- exp(ci[2])
        p_mod <- fit$prob["dosage"]
      } else {
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
          or    <- as.numeric(fish$estimate)
          or_lo <- NA_real_
          or_hi <- NA_real_
          p_mod <- fish$p.value
        }
      }

      hap_bits <- strsplit(h, "")[[1]]
      alt_idx <- which(hap_bits == "1")
      hap_def <- if (length(alt_idx) == 0) {
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
      comp_df <- bind_rows(comp_rows) %>%
        arrange(P_Fisher)
      results[[comp$label]] <- comp_df
    }
  }

  if (length(results) == 0) return(NULL)
  out <- bind_rows(results)
  # Single BH correction across ALL haplotype-comparison tests in this region.
  out$FDR                    <- p.adjust(out$P_Fisher, method = "BH")
  out$Significant_Fisher     <- out$P_Fisher < 0.05
  out$Significant_FDR        <- out$FDR < 0.05
  # Add short Haplotype_ID by ranking on Global_Freq (H1 = most frequent).
  # Use all haplotypes across all comparisons so ranks are consistent.
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

make_ld_plot <- function(ld_matrix, snp_meta, out_file, ld_blocks = NULL) {
  ld_df <- as.data.frame(as.table(ld_matrix), stringsAsFactors = FALSE)
  names(ld_df) <- c("SNP1", "SNP2", "R2")
  pos_map <- setNames(snp_meta$POS_int, snp_meta$rsID_clean)
  ld_df$POS1 <- pos_map[ld_df$SNP1]
  ld_df$POS2 <- pos_map[ld_df$SNP2]

  snp_levels <- snp_meta$rsID_clean

  # Build block boundary rectangles if block info provided
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
                  colour = "#E69F00", fill = NA, linewidth = 0.8)
      else list() } +
    scale_fill_gradientn(
      colours  = c("#00204C", "#2E5EAA", "#4E9F9F", "#A5B85C", "#F0E442"),
      na.value = "grey95",
      limits   = c(0, 1),
      labels   = percent_format(accuracy = 1),
      name     = expression(r^2)
    ) +
    labs(title    = "LD Matrix Across the SNP Panel",
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
  if (nrow(freq_df) == 0) return(invisible(NULL))
  # Canonical order: H1 → H2 → … (already ranked by descending frequency in make_em_freq_table)
  hap_levels <- freq_df %>% arrange(desc(Frequency)) %>% pull(Haplotype_ID)
  freq_df <- freq_df %>%
    mutate(Haplotype_ID = factor(Haplotype_ID, levels = hap_levels))
  p <- ggplot(freq_df, aes(x = Haplotype_ID,
                           y = Frequency * 100)) +
    geom_col(fill = "#0072B2", colour = "black", alpha = 0.85) +
    geom_text(aes(label = sprintf("%.1f%%", Frequency * 100)), vjust = -0.35, size = 3.2) +
    labs(title = title, subtitle = subtitle, x = "Haplotype", y = "Frequency (%)") +
    coord_cartesian(ylim = c(0, max(freq_df$Frequency * 100) * 1.15))
  SAVE(out_file, p, width = max(7, nrow(freq_df) * 0.7), height = 5)
}

# -----------------------------------------------------------------------------
# FIGURE: Forest plot — OR + 95% CI per haplotype, faceted by comparison group
# Y-axis uses short IDs (H1, H2…); OR and p-value shown as text on each row.
# -----------------------------------------------------------------------------
make_forest_plot <- function(assoc_df, region_name, out_file) {
  if (is.null(assoc_df) || nrow(assoc_df) == 0) return(invisible(NULL))

  plot_df <- assoc_df %>%
    filter(!is.na(OR), !is.na(OR_95CI_Lo), !is.na(OR_95CI_Hi)) %>%
    mutate(
      # Short y-axis label: ID + global frequency
      label      = sprintf("%s (%.0f%%)", Haplotype_ID, Global_Freq * 100),
      sig_label  = case_when(
        FDR < 0.05       ~ "FDR < 0.05",
        P_Fisher < 0.05  ~ "p < 0.05",
        TRUE             ~ "ns"
      ),
      OR_lo_plot = pmax(OR_95CI_Lo, 0.05),
      OR_hi_plot = pmin(OR_95CI_Hi, 20),
      # Annotation text: OR [lo–hi]  p=xxx
      or_text    = sprintf("OR %.2f [%.2f–%.2f]  p=%s",
                           OR,
                           OR_95CI_Lo,
                           OR_95CI_Hi,
                           ifelse(P_Fisher < 0.001,
                                  sprintf("%.2e", P_Fisher),
                                  sprintf("%.3f", P_Fisher)))
    )

  if (nrow(plot_df) == 0) return(invisible(NULL))

  # Canonical y-axis order: H1 at top → HN at bottom (i.e. ascending so ggplot puts H1 on top)
  # Build label order from Global_Freq so it is consistent with H-numbering
  hap_order <- plot_df %>%
    distinct(Haplotype_ID, Global_Freq, label) %>%
    arrange(Global_Freq) %>%   # ascending → H1 ends up at top of y-axis
    pull(label)
  plot_df <- plot_df %>%
    mutate(label = factor(label, levels = hap_order))

  # x-axis range for placing text labels just past the right edge
  x_text_pos <- 10^(log10(20) * 1.05)

  p <- ggplot(plot_df, aes(x = OR, y = label,
                           colour = sig_label, shape = sig_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = OR_lo_plot, xmax = OR_hi_plot),
                   height = 0.25, linewidth = 0.6) +
    geom_point(size = 3) +
    geom_text(aes(x = OR_hi_plot, label = or_text),
              hjust = -0.07, size = 2.8, colour = "grey30") +
    scale_colour_manual(values = c("FDR < 0.05" = "#D55E00",
                                   "p < 0.05"   = "#E69F00",
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
    labs(title    = sprintf("Haplotype Association: %s", region_name),
         subtitle = "Firth-penalised logistic regression | OR (95% CI) | log10 x-axis | label = Haplotype ID (global freq%)",
         x = "Odds ratio (log scale)", y = NULL) +
    theme(axis.text.y  = element_text(size = 8.5, face = "bold"),
          strip.text   = element_text(face = "bold"),
          legend.position = "right",
          plot.margin  = margin(6, 160, 6, 8))

  n_haps <- length(unique(plot_df$label))
  n_comp <- length(unique(plot_df$Comparison))
  SAVE(out_file, p,
       width  = 11,
       height = max(5, n_haps * 0.5 * n_comp + 2.5))
}

# -----------------------------------------------------------------------------
# FIGURE: Haplotype allele composition tile grid
# Each row = one haplotype; each column = one SNP; filled = alt allele carried
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
      Allele     = ifelse(Allele_Code == 2, "Alt", "Ref"),
      SNP        = factor(SNP, levels = snp_cols),
      Haplotype_ID = factor(Haplotype_ID,
                            levels = rev(hap_mat$Haplotype_ID)),
      Freq_label = sprintf("%.1f%%", Frequency * 100)
    )

  # Gene annotation strip along x-axis
  gene_map <- setNames(snp_meta_sub$Gene, snp_meta_sub$rsID_clean)
  tile_df$Gene <- gene_map[as.character(tile_df$SNP)]

  p <- ggplot(tile_df, aes(x = SNP, y = Haplotype_ID, fill = Allele)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_manual(values = c("Alt" = "#D55E00", "Ref" = "#F2F2F2"),
                      name = "Allele") +
    # Frequency labels on the right
    geom_text(data = hap_mat %>%
                mutate(Haplotype_ID = factor(Haplotype_ID, levels = rev(hap_mat$Haplotype_ID))),
              aes(x = length(snp_cols) + 0.6,
                  y = Haplotype_ID,
                  label = sprintf("%.1f%%", Frequency * 100)),
              inherit.aes = FALSE, size = 3, hjust = 0) +
    coord_cartesian(clip = "off") +
    labs(title    = sprintf("Haplotype Allelic Composition: %s", region_name),
         subtitle = "Red = alt allele carried | rows ordered by global frequency",
         x = NULL, y = "Haplotype") +
    theme(
      axis.text.x    = element_text(angle = 45, hjust = 1, size = 7.5),
      axis.text.y    = element_text(size = 8),
      plot.margin    = margin(8, 50, 6, 8),
      legend.position = "top"
    )

  n_haps <- nrow(hap_mat)
  n_snps <- length(snp_cols)
  SAVE(out_file, p,
       width  = max(7, n_snps * 0.55 + 2),
       height = max(4, n_haps * 0.38 + 2))
}

# -----------------------------------------------------------------------------
# FIGURE: Tumour vs Healthy haplotype frequency grouped bar chart
# X-axis uses short IDs (H1, H2…); a definition table is printed below.
# -----------------------------------------------------------------------------
make_freq_comparison_plot <- function(assoc_df, region_name, out_file,
                                      comparison = "All") {
  if (is.null(assoc_df) || nrow(assoc_df) == 0) return(invisible(NULL))

  sub <- assoc_df %>%
    filter(Comparison == comparison) %>%
    arrange(desc(Global_Freq)) %>%
    mutate(
      sig_mark = case_when(
        FDR < 0.05      ~ "**",
        P_Fisher < 0.05 ~ "*",
        TRUE            ~ ""
      ),
      # Short label for x-axis
      Hap_Label = Haplotype_ID
    )

  if (nrow(sub) == 0) return(invisible(NULL))

  # Canonical x-axis order: H1 → H2 → … (descending frequency, already in sub after arrange)
  hap_x_levels <- sub %>% pull(Hap_Label)   # stable: sub is already sorted by desc(Global_Freq)

  plot_df <- sub %>%
    select(Hap_Label, Haplotype_ID, Freq_Tumour, Freq_Healthy,
           P_Fisher, FDR, sig_mark, Global_Freq) %>%
    pivot_longer(cols = c(Freq_Tumour, Freq_Healthy),
                 names_to = "Group", values_to = "Freq") %>%
    mutate(
      Group     = recode(Group, "Freq_Tumour" = "Tumour", "Freq_Healthy" = "Pooled control"),
      Group     = factor(Group, levels = c("Pooled control", "Tumour")),
      Hap_Label = factor(Hap_Label, levels = hap_x_levels)
    )

  # Significance marker positions (above the taller bar)
  sig_df <- plot_df %>%
    group_by(Hap_Label, sig_mark) %>%
    summarise(y_max = max(Freq, na.rm = TRUE), .groups = "drop") %>%
    filter(sig_mark != "")

  bar_plot <- ggplot(plot_df,
                     aes(x = Hap_Label, y = Freq * 100, fill = Group)) +
    geom_col(position = position_dodge(0.72), width = 0.62,
             colour = "black", linewidth = 0.3, alpha = 0.88) +
    { if (nrow(sig_df) > 0)
        geom_text(data = sig_df,
                  aes(x = Hap_Label, y = pmin(y_max * 100 + 3, 98),
                      label = sig_mark),
                  inherit.aes = FALSE, size = 5, colour = "#D55E00")
      else list() } +
    scale_fill_manual(values = c("Pooled control" = "#0072B2", "Tumour" = "#D55E00"),
                      name = NULL) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      labels = function(x) sprintf("%d%%", as.integer(x))
    ) +
    labs(title    = sprintf("Haplotype Carrier Frequency: %s (%s)",
                            region_name, comparison),
         subtitle = "* p\u202f<\u202f0.05 (Fisher)   ** FDR\u202f<\u202f0.05   bars show tumour versus pooled-control carrier frequency",
         x = NULL, y = "Haplotype frequency (%)") +
    theme(axis.text.x   = element_text(size = 9, face = "bold"),
          legend.position = "top",
          plot.margin   = margin(6, 8, 2, 8))

  # ── Definition table ────────────────────────────────────────────────────────
  tbl_df <- sub %>%
    transmute(
      ID       = Haplotype_ID,
      `Global%`  = sprintf("%.1f%%", Global_Freq * 100),
      `Alt SNPs` = ifelse(nchar(Haplotype_Def) > 60,
                          paste0(substr(Haplotype_Def, 1, 57), "…"),
                          Haplotype_Def),
      `p (Fisher)` = ifelse(P_Fisher < 0.001,
                             sprintf("%.2e", P_Fisher),
                             sprintf("%.3f", P_Fisher)),
      FDR_col  = ifelse(FDR < 0.001,
                        sprintf("%.2e", FDR),
                        sprintf("%.3f", FDR))
    ) %>%
    rename(FDR = FDR_col)

  # Render as a ggplot table using annotate
  tbl_plot <- ggplot() +
    theme_void() +
    annotation_custom(
      gridExtra::tableGrob(
        tbl_df,
        rows = NULL,
        theme = gridExtra::ttheme_minimal(
          core    = list(fg_params = list(cex = 0.72, hjust = 0, x = 0.02)),
          colhead = list(fg_params = list(cex = 0.75, fontface = "bold",
                                          hjust = 0, x = 0.02)),
          padding = unit(c(2, 4), "mm")
        )
      )
    )

  combined <- patchwork::wrap_plots(bar_plot, tbl_plot,
                                     ncol = 1,
                                     heights = c(3, max(1, nrow(tbl_df) * 0.28)))

  n_haps <- nrow(sub)
  SAVE(out_file, combined,
       width  = max(7, n_haps * 0.85 + 2),
       height = max(6, n_haps * 0.28 + 5))
}

# -----------------------------------------------------------------------------
# FIGURE: Haplotype dosage composition by comparison arm
# For each tested haplotype, show the percentage of samples carrying 0, 1, or 2
# copies within the comparison arms. This complements the allele tile grid by
# revealing whether a signal is driven by non-carriers, heterozygous carriers,
# or homozygous carriers.
# -----------------------------------------------------------------------------
make_haplotype_dosage_plot <- function(a1_mat, a2_mat, assoc_df, metadata_df,
                                       region_name, out_file,
                                       comparison = "All") {
  if (is.null(assoc_df) || nrow(assoc_df) == 0) return(invisible(NULL))

  sub <- assoc_df %>%
    filter(Comparison == comparison) %>%
    arrange(desc(Global_Freq))

  if (nrow(sub) == 0) return(invisible(NULL))

  keep_idx <- switch(
    comparison,
    "Breast"      = (metadata_df$Cohort == "Breast" & metadata_df$Tissue == "Tumour") | metadata_df$Tissue == "Healthy",
    "Endometrium" = (metadata_df$Cohort == "Endometrium" & metadata_df$Tissue == "Tumour") | metadata_df$Tissue == "Healthy",
    metadata_df$Tissue %in% c("Healthy", "Tumour")
  )

  idx <- which(keep_idx)
  if (length(idx) == 0) return(invisible(NULL))

  meta_sub <- metadata_df[idx, , drop = FALSE]
  hap1_sub <- make_binary_hap_strings(a1_mat[idx, , drop = FALSE])
  hap2_sub <- make_binary_hap_strings(a2_mat[idx, , drop = FALSE])

  group_counts <- meta_sub %>%
    mutate(Group_Label = ifelse(Tissue == "Healthy", "Pooled control", "Tumour")) %>%
    count(Group_Label, name = "N_Group")

  dosage_rows <- lapply(seq_len(nrow(sub)), function(i) {
    hap_string <- sub$Haplotype_String[i]
    dosage <- as.integer(hap1_sub == hap_string) + as.integer(hap2_sub == hap_string)
    data.frame(
      Sample        = meta_sub$Sample,
      Group_Label   = ifelse(meta_sub$Tissue == "Healthy", "Pooled control", "Tumour"),
      Haplotype_ID  = sub$Haplotype_ID[i],
      Global_Freq   = sub$Global_Freq[i],
      Dosage_Class  = factor(
        dosage,
        levels = c(0, 1, 2),
        labels = c("Non-carrier", "Heterozygous carrier", "Homozygous carrier")
      ),
      stringsAsFactors = FALSE
    )
  })

  plot_df <- bind_rows(dosage_rows) %>%
    count(Group_Label, Haplotype_ID, Global_Freq, Dosage_Class, name = "N") %>%
    group_by(Group_Label, Haplotype_ID) %>%
    mutate(Percent = 100 * N / sum(N)) %>%
    ungroup() %>%
    left_join(group_counts, by = "Group_Label") %>%
    mutate(
      Group_Label = sprintf("%s (n=%d)", Group_Label, N_Group),
      Haplotype_ID = factor(Haplotype_ID, levels = sub$Haplotype_ID),
      Dosage_Class = factor(Dosage_Class, levels = c("Non-carrier", "Heterozygous carrier", "Homozygous carrier"))
    )

  p <- ggplot(plot_df, aes(x = Haplotype_ID, y = Percent, fill = Dosage_Class)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.3) +
    geom_text(
      data = plot_df %>% filter(Percent >= 8),
      aes(label = sprintf("%.0f%%", Percent)),
      position = position_stack(vjust = 0.5),
      size = 3,
      colour = "white",
      fontface = "bold"
    ) +
    facet_wrap(~ Group_Label, ncol = 1) +
    scale_fill_manual(
      values = c("Non-carrier" = "#D9D9D9", "Heterozygous carrier" = "#56B4E9", "Homozygous carrier" = "#D55E00"),
      name = "Haplotype status"
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = c(0, 0),
      breaks = seq(0, 100, by = 20),
      labels = function(x) sprintf("%.0f%%", x)
    ) +
    labs(
      title = sprintf("Haplotype Dosage Composition: %s (%s)", region_name, comparison),
      subtitle = "Bars are proportional within each arm and compare cohort tumours against the pooled healthy control arm",
      x = NULL,
      y = "Samples within comparison arm (%)"
    ) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(size = 9, face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  n_haps <- nrow(sub)
  SAVE(out_file, p,
       width = max(7, n_haps * 0.75 + 2),
       height = 6.4)
}

# -----------------------------------------------------------------------------
# FIGURE: Association -log10(p) bubble plot across all regions
# Bubble size = N haplotypes tested; colour = region type
# -----------------------------------------------------------------------------
make_association_overview_plot <- function(summary_df, out_file) {
  plot_df <- summary_df %>%
    filter(!is.na(Best_P_Fisher)) %>%
    mutate(
      neg_log_p  = -log10(Best_P_Fisher),
      neg_log_fdr = -log10(pmax(Best_FDR, 1e-6)),
      Region_Type = factor(Region_Type, levels = c("Full", "LD_Block", "Gene")),
      sig_label  = case_when(
        Best_FDR < 0.05  ~ "FDR < 0.05",
        Best_P_Fisher < 0.05 ~ "p < 0.05",
        TRUE             ~ "ns"
      ),
      Region_short = gsub("Gene_GSDMB_", "GSDMB\n", gsub("LD_", "", Region))
    )

  if (nrow(plot_df) == 0) return(invisible(NULL))

  # Bonferroni threshold over total haplotype-comparison tests (sum of Haps_Shown
  # across all regions times the number of comparison groups used).
  # Fall back to number of shown haplotypes if N_Tests is absent.
  n_tests <- if ("N_Tests" %in% names(plot_df)) sum(plot_df$N_Tests, na.rm = TRUE) else sum(plot_df$Haps_Shown)
  n_tests <- max(n_tests, 1L)

  p <- ggplot(plot_df, aes(x = fct_reorder(Region_short, neg_log_p, .desc = TRUE),
                           y = neg_log_p,
                           size = Haps_Shown,
                           colour = Region_Type,
                           shape = sig_label)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "grey60", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05 / n_tests),
               linetype = "dotted", colour = "#D55E00", linewidth = 0.5) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(3, 9), name = "Haplotypes\ntested") +
    scale_colour_manual(
      values = c("Full" = "#666666", "LD_Block" = "#CC79A7", "Gene" = "#009E73"),
      name = "Region type") +
    scale_shape_manual(
      values = c("FDR < 0.05" = 18, "p < 0.05" = 17, "ns" = 16),
      name = "Significance") +
    annotate("text", x = Inf, y = -log10(0.05) + 0.1,
             label = "p = 0.05", hjust = 1.1, size = 3, colour = "grey50") +
    annotate("text", x = Inf,
             y = -log10(0.05 / n_tests) + 0.1,
             label = sprintf("Bonferroni (n=%d tests)", n_tests),
             hjust = 1.1, size = 3, colour = "#D55E00") +
    labs(title    = "Regional Haplotype Association Overview",
         subtitle = "-log10(best Fisher p) per region | bubble size = haplotypes tested",
         x = NULL, y = expression(-log[10](p))) +
    theme(axis.text.x   = element_text(angle = 40, hjust = 1, size = 7.5),
          legend.position = "right")

  SAVE(out_file, p,
       width  = max(10, nrow(plot_df) * 0.65 + 2),
       height = 6)
}


make_haplotype_atlas_plot <- function(freq_df, snp_meta_sub, out_file,
                                      top_n = MAIN_TOP_FULL_HAPS,
                                      max_snps = 12) {
  if (nrow(freq_df) == 0) return(invisible(NULL))

  visible <- freq_df %>% arrange(desc(Frequency)) %>% slice_head(n = top_n)
  if (nrow(visible) == 0) return(invisible(NULL))

  snp_cols <- intersect(snp_meta_sub$rsID_clean, names(visible))
  if (length(snp_cols) == 0) return(invisible(NULL))
  variability <- sapply(snp_cols, function(col) dplyr::n_distinct(visible[[col]]))
  selected <- names(sort(variability, decreasing = TRUE))
  selected <- selected[seq_len(min(length(selected), max_snps))]
  if (length(selected) == 0) selected <- snp_cols[seq_len(min(length(snp_cols), max_snps))]

  bar_df <- visible %>%
    transmute(Haplotype_ID = factor(Haplotype_ID, levels = rev(Haplotype_ID)),
              FrequencyPct = Frequency * 100,
              Count = Count)

  other_freq <- sum(freq_df$Frequency, na.rm = TRUE) - sum(visible$Frequency, na.rm = TRUE)
  other_count <- sum(freq_df$Count, na.rm = TRUE) - sum(visible$Count, na.rm = TRUE)
  if (other_freq > 0) {
    bar_df <- bind_rows(
      bar_df,
      data.frame(Haplotype_ID = factor("Other", levels = c(levels(bar_df$Haplotype_ID), "Other")),
                 FrequencyPct = other_freq * 100,
                 Count = other_count)
    )
  }

  bar_plot <- ggplot(bar_df, aes(x = FrequencyPct, y = Haplotype_ID)) +
    geom_col(fill = "#0B5C8C", alpha = 0.9) +
    geom_text(aes(label = sprintf("%.1f%%  (n=%s)", FrequencyPct, Count)), hjust = -0.05, size = 3.1) +
    coord_cartesian(xlim = c(0, max(bar_df$FrequencyPct) * 1.25), clip = "off") +
    labs(title = "Full-region haplotype atlas",
         subtitle = "Top full-region haplotypes by global frequency; rare remainder pooled as Other",
         x = "Global frequency (%)", y = NULL) +
    theme(plot.margin = margin(8, 20, 8, 8))

  tile_df <- visible %>%
    select(Haplotype_ID, Frequency, all_of(selected)) %>%
    arrange(desc(Frequency)) %>%
    pivot_longer(cols = all_of(selected), names_to = "SNP", values_to = "Allele_Code") %>%
    mutate(
      Allele = ifelse(Allele_Code == 2, "Alt", "Ref"),
      Haplotype_ID = factor(Haplotype_ID, levels = rev(visible$Haplotype_ID)),
      SNP = factor(SNP, levels = selected)
    )

  tile_plot <- ggplot(tile_df, aes(x = SNP, y = Haplotype_ID, fill = Allele)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_manual(values = c("Alt" = "#D55E00", "Ref" = "#F2F2F2"), name = "Allele") +
    labs(title = NULL, x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
          legend.position = "top")

  combined <- patchwork::wrap_plots(bar_plot, tile_plot, ncol = 2, widths = c(1.1, 1.8))
  SAVE(out_file, combined, width = max(11, length(selected) * 0.55 + 6), height = max(5.5, nrow(visible) * 0.45 + 2.5))
}

make_full_region_comparison_plot <- function(a1_mat, a2_mat, freq_df, metadata_df, out_file,
                                             top_n = MAIN_TOP_COMPARISON_HAPS) {
  if (nrow(freq_df) == 0) return(invisible(NULL))
  shown <- freq_df %>% arrange(desc(Frequency)) %>% slice_head(n = top_n)
  if (nrow(shown) == 0) return(invisible(NULL))
  hap1 <- make_binary_hap_strings(a1_mat)
  hap2 <- make_binary_hap_strings(a2_mat)
  group_idx <- list(
    "Pooled control" = which(metadata_df$Tissue == "Healthy"),
    "Breast tumour" = which(metadata_df$Cohort == "Breast" & metadata_df$Tissue == "Tumour"),
    "Endometrium tumour" = which(metadata_df$Cohort == "Endometrium" & metadata_df$Tissue == "Tumour")
  )
  rows <- bind_rows(lapply(seq_len(nrow(shown)), function(i) {
    hap <- shown$Allele_String[i]
    dosage <- as.integer(hap1 == hap) + as.integer(hap2 == hap)
    bind_rows(lapply(names(group_idx), function(group_name) {
      idx <- group_idx[[group_name]]
      data.frame(Haplotype_ID = shown$Haplotype_ID[i],
                 Group = group_name,
                 FrequencyPct = if (length(idx) > 0) mean(dosage[idx], na.rm = TRUE) * 50 else NA_real_)
    }))
  }))
  p <- ggplot(rows, aes(x = Haplotype_ID, y = FrequencyPct, colour = Group, group = Group)) +
    geom_line(position = position_dodge(width = 0.25), alpha = 0.7) +
    geom_point(position = position_dodge(width = 0.25), size = 3) +
    scale_colour_manual(values = c("Pooled control" = "#666666", "Breast tumour" = "#B55D6A", "Endometrium tumour" = "#3D7EA6")) +
    labs(title = "Full-region cohort comparison",
         subtitle = "Main full-region haplotypes across pooled controls and tumour cohorts",
         x = NULL, y = "Estimated haplotype frequency (%)") +
    theme(legend.position = "top")
  SAVE(out_file, p, width = max(8.5, nrow(shown) * 0.8 + 2), height = 5.5)
}

make_key_block_overview_plot <- function(summary_df, block_assoc_tables, gene_assoc_tables, out_file,
                                         max_regions = MAX_KEY_BLOCKS_MAIN) {
  rank_df <- summary_df %>%
    filter(Region_Type != "Full", !is.na(Best_P_Fisher), Haps_Shown > 0) %>%
    mutate(GenePriority = ifelse(Region_Type == "Gene", 1, 0)) %>%
    arrange(Best_FDR, Best_P_Fisher, desc(Haps_Shown), desc(GenePriority)) %>%
    slice_head(n = max_regions)
  if (nrow(rank_df) == 0) return(invisible(NULL))

  rows <- bind_rows(lapply(rank_df$Region, function(region_nm) {
    assoc_df <- if (region_nm %in% names(gene_assoc_tables)) gene_assoc_tables[[region_nm]] else block_assoc_tables[[region_nm]]
    if (is.null(assoc_df) || nrow(assoc_df) == 0) return(NULL)
    assoc_df %>%
      filter(Comparison %in% c("Breast", "Endometrium")) %>%
      group_by(Comparison) %>%
      arrange(desc(Global_Freq), P_Fisher) %>%
      slice_head(n = 5) %>%
      ungroup() %>%
      transmute(Region = region_nm,
                Haplotype_ID = Haplotype_ID,
                Group = ifelse(Comparison == "Breast", "Breast tumour", "Endometrium tumour"),
                FrequencyPct = Freq_Tumour * 100)
  }))
  if (is.null(rows) || nrow(rows) == 0) return(invisible(NULL))
  p <- ggplot(rows, aes(x = Haplotype_ID, y = FrequencyPct, colour = Group, group = Group)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 2.8) +
    facet_wrap(~Region, scales = "free_x") +
    scale_colour_manual(values = c("Breast tumour" = "#B55D6A", "Endometrium tumour" = "#3D7EA6")) +
    labs(title = "Key-block supplementary highlights",
         subtitle = "Top-ranked haplotypes from the strongest non-full-region blocks",
         x = NULL, y = "Tumour frequency (%)") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 35, hjust = 1))
  SAVE(out_file, p, width = 12, height = max(5.5, nrow(rank_df) * 2.1))
}

# -----------------------------------------------------------------------------
# FIGURE: Per-block LD heatmap (for GSDMB gene blocks only)
# -----------------------------------------------------------------------------
make_block_ld_plot <- function(ld_matrix_full, block_idx, snp_meta_full,
                               block_name, out_file) {
  sub_snps <- snp_meta_full$rsID_clean[block_idx]
  sub_ld   <- ld_matrix_full[sub_snps, sub_snps, drop = FALSE]
  sub_meta <- snp_meta_full[block_idx, , drop = FALSE]
  make_ld_plot(sub_ld, sub_meta, out_file)
}

# -----------------------------------------------------------------------------
# FIGURE: Haplotype frequency stacked bar — cohort breakdown
# Shows how each haplotype is distributed across Breast/Endometrium × Healthy/Tumour
# -----------------------------------------------------------------------------
make_cohort_freq_plot <- function(a1_mat, a2_mat, freq_df, metadata_df,
                                  region_name, out_file, min_freq = 0.01,
                                  hap_palette = NULL) {
  # Canonical order: H1 → HN by descending frequency, then "Other" last.
  # freq_df rows are already ranked by Frequency desc (from make_em_freq_table).
  shown_haps <- freq_df %>%
    filter(Frequency >= min_freq) %>%
    arrange(desc(Frequency)) %>%
    pull(Haplotype_ID)
  if (length(shown_haps) == 0) return(invisible(NULL))

  hap1 <- make_binary_hap_strings(a1_mat)
  hap2 <- make_binary_hap_strings(a2_mat)

  # Map binary strings -> Haplotype_ID via freq_df
  str_to_id <- setNames(freq_df$Haplotype_ID, freq_df$Allele_String)

  # Vectorised: build long-form table directly without row-by-row appending
  n <- nrow(metadata_df)
  hid1 <- str_to_id[hap1]
  hid2 <- str_to_id[hap2]
  hid1[is.na(hid1) | !(hid1 %in% shown_haps)] <- "Other"
  hid2[is.na(hid2) | !(hid2 %in% shown_haps)] <- "Other"
  group_vec <- paste(metadata_df$Cohort, metadata_df$Tissue, sep = "\n")

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
      Group = factor(Group, levels = sort(unique(Group)))
    )

  # Colour palette: use global palette if provided for cross-figure consistency
  n_haps  <- length(shown_haps)
  if (!is.null(hap_palette) && all(shown_haps %in% names(hap_palette))) {
    hap_pal <- hap_palette[c(shown_haps, "Other")]
  } else {
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
    labs(title    = sprintf("Haplotype Frequency by Cohort: %s", region_name),
         subtitle = "Stacked bars show haplotype frequency within each group",
         x = NULL, y = "Frequency (%)") +
    theme(axis.text.x    = element_text(size = 9),
          legend.position = "right")

  SAVE(out_file, p, width = max(7, length(unique(plot_df$Group)) * 1.2 + 3), height = 5.5)
}

# -----------------------------------------------------------------------------
# FIGURE: Manhattan-style plot — p-values per haplotype coloured by gene
# (full-region association results plotted by SNP position)
# -----------------------------------------------------------------------------
make_manhattan_hap_plot <- function(assoc_df, snp_meta, out_file,
                                    comparison = "All") {
  if (is.null(assoc_df) || nrow(assoc_df) == 0) return(invisible(NULL))

  sub <- assoc_df %>% filter(Comparison == comparison, !is.na(P_Fisher))
  if (nrow(sub) == 0) return(invisible(NULL))

  # For each haplotype, get the median position of its alt SNPs
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
      med_pos   = get_median_pos(Haplotype_Def, snp_meta),
      neg_log_p = -log10(P_Fisher),
      gene_label = {
        snp_ids <- trimws(strsplit(Haplotype_Def, ";")[[1]])
        genes   <- snp_meta$Gene[snp_meta$rsID_clean %in% snp_ids]
        if (length(genes) == 0) "Unknown" else paste(unique(genes), collapse = "/")
      }
    ) %>%
    ungroup()

  gene_lvls <- unique(snp_meta$Gene)
  gene_cols <- setNames(
    scales::hue_pal()(length(gene_lvls)),
    gene_lvls
  )

  p <- ggplot(sub, aes(x = med_pos / 1e6, y = neg_log_p,
                       colour = gene_label, size = Global_Freq * 100)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_point(alpha = 0.8) +
    scale_colour_manual(values = gene_cols, name = "Gene") +
    scale_size_continuous(range = c(2, 7), name = "Global freq (%)") +
    labs(title    = sprintf("Haplotype Association Position: %s", comparison),
         subtitle = "-log10(Fisher p) | point size = global haplotype frequency | coloured by gene",
         x = "Genomic position (Mb)", y = expression(-log[10](p))) +
    theme(legend.position = "right")

  SAVE(out_file, p, width = 10, height = 5.5)
}


# -----------------------------------------------------------------------------
# FIGURE: Kaplan-Meier survival — carrier vs non-carrier of top haplotype
# Requires: survival, survminer packages
# Runs for both OS and PFS where available; handles metadata with os_months /
# pfs_months columns if present (merged in from clinical data).
# If the metadata does not contain survival columns this function skips silently.
# -----------------------------------------------------------------------------
make_km_plot <- function(a1_mat, a2_mat, freq_df, metadata_df,
                         region_name, out_dir,
                         hap_rank = 1,
                         time_col = "os_months", event_col = "os_event") {

  if (!HAS_SURVIVAL || !HAS_SURVMINER) return(invisible(NULL))
  if (!all(c(time_col, event_col) %in% names(metadata_df))) return(invisible(NULL))

  top_hap_str <- freq_df$Allele_String[hap_rank]
  top_hap_id  <- freq_df$Haplotype_ID[hap_rank]
  if (is.na(top_hap_str)) return(invisible(NULL))

  hap1 <- make_binary_hap_strings(a1_mat)
  hap2 <- make_binary_hap_strings(a2_mat)

  carrier <- as.integer(hap1 == top_hap_str) + as.integer(hap2 == top_hap_str)
  carrier_flag <- factor(ifelse(carrier >= 1, "Carrier", "Non-carrier"),
                         levels = c("Non-carrier", "Carrier"))

  surv_df <- data.frame(
    time    = as.numeric(metadata_df[[time_col]]),
    event   = as.integer(metadata_df[[event_col]]),
    carrier = carrier_flag
  ) %>% filter(!is.na(time), !is.na(event), time >= 0)

  if (nrow(surv_df) < 10 || length(unique(surv_df$carrier)) < 2) return(invisible(NULL))

  fit  <- survival::survfit(survival::Surv(time, event) ~ carrier, data = surv_df)
  pval <- survival::survdiff(survival::Surv(time, event) ~ carrier, data = surv_df)
  pval_text <- sprintf("Log-rank p = %.3f", 1 - pchisq(pval$chisq, df = 1))

  p <- survminer::ggsurvplot(
    fit,
    data         = surv_df,
    pval         = FALSE,
    conf.int     = TRUE,
    risk.table   = TRUE,
    palette      = c("#0072B2", "#D55E00"),
    legend.labs  = c("Non-carrier", sprintf("Carrier (%s)", top_hap_id)),
    title        = sprintf("Survival by %s carrier status — %s", top_hap_id, region_name),
    subtitle     = pval_text,
    xlab         = ifelse(grepl("pfs", time_col, ignore.case = TRUE),
                          "Time (months) — PFS", "Time (months) — OS"),
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
# 3. DATA LOADING
# =============================================================================
cat("\n======================================================================\n")
cat("SCRIPT 15: HAPLOTYPE ANALYSIS WITH FULL-REGION + LD BLOCKS + GENES\n")
cat("======================================================================\n\n")
cat(sprintf("  Main frequency threshold:       %.1f%%\n", MAIN_MIN_FREQ * 100))
cat(sprintf("  Supplementary frequency cutoff: %.1f%%\n", SUPP_MIN_FREQ * 100))
cat(sprintf("  LD block thresholds (adjacent): %s\n\n", paste(sprintf("r^2 >= %.2f", LD_BLOCK_THRESHOLDS), collapse = " and ")))

if (!requireNamespace("survminer", quietly = TRUE)) {
  cat("  NOTE: survminer not available — KM plots will be skipped.\n")
  cat("  To enable, run: conda install -c conda-forge r-survminer\n\n")
}

cat("STEP 1: Loading phased genotypes and annotations\n")
cat("----------------------------------------------------------------------\n")

target_rsids <- read_excel(WHITELIST) %>%
  pull(Variant_ID) %>%
  as.character() %>%
  trimws() %>%
  unique()

annot <- read_excel(ANNOTATED, sheet = "Biological_Annotations")
annotation_map <- annot %>%
  mutate(
    rsID_clean = str_extract(Existing_variation, "rs[0-9]+"),
    POS_int    = as.integer(gsub("[^0-9]", "", as.character(POS))),
    REF        = toupper(trimws(as.character(REF))),
    ALT        = toupper(trimws(as.character(ALT))),
    Gene       = ifelse(is.na(SYMBOL) | trimws(SYMBOL) == "", "Unknown", trimws(SYMBOL))
  ) %>%
  filter(rsID_clean %in% target_rsids) %>%
  select(POS_int, rsID_clean, Gene, REF, ALT) %>%
  distinct(POS_int, REF, ALT, .keep_all = TRUE) %>%
  arrange(POS_int)

geno_raw <- read.table(GENO_FILE, header = TRUE, sep = "	", check.names = FALSE, quote = "") %>%
  mutate(
    POS_int = as.integer(POS),
    REF     = toupper(trimws(as.character(REF))),
    ALT     = toupper(trimws(as.character(ALT)))
  )
meta_raw <- read.table(META_FILE, header = TRUE, sep = "	", check.names = FALSE, quote = "")

stopifnot(all(c("Sample", "Cohort", "Tissue") %in% names(meta_raw)))

geno_filtered <- geno_raw %>%
  inner_join(annotation_map, by = c("POS_int", "REF", "ALT")) %>%
  arrange(POS_int) %>%
  distinct(rsID_clean, .keep_all = TRUE)   # guard against any residual many:many collisions

sample_cols <- intersect(names(geno_filtered), meta_raw$Sample)
meta <- meta_raw %>%
  filter(Sample %in% sample_cols) %>%
  mutate(Group = paste(Cohort, Tissue, sep = "_"))

sample_cols <- sample_cols[sample_cols %in% meta$Sample]
meta <- meta[match(sample_cols, meta$Sample), , drop = FALSE]
stopifnot(all(sample_cols == meta$Sample))

cat(sprintf("✓ Genotype table loaded: %d SNPs × %d samples\n", nrow(geno_filtered), length(sample_cols)))
cat(sprintf("✓ Metadata loaded:       %d samples retained\n", nrow(meta)))
cat(sprintf("✓ Unique genes covered:  %s\n", paste(unique(geno_filtered$Gene), collapse = ", ")))

# =============================================================================
# 4. BUILD ALLELE MATRICES
# =============================================================================
cat("\nSTEP 2: Building phased allele matrices\n")
cat("----------------------------------------------------------------------\n")

am <- build_allele_matrices(geno_filtered, sample_cols)
a1 <- am$a1
a2 <- am$a2
snp_meta <- geno_filtered %>%
  select(POS_int, rsID_clean, Gene, REF, ALT) %>%
  distinct(rsID_clean, .keep_all = TRUE) %>%   # belt-and-braces; geno_filtered already deduped
  arrange(POS_int)

n_samples <- nrow(a1)
n_snps    <- ncol(a1)
n_missing <- sum(is.na(a1) | is.na(a2))
total_alleles <- 2 * n_samples * n_snps

cat(sprintf("✓ Allele matrices built: %d samples × %d SNPs\n", n_samples, n_snps))
cat(sprintf("  Missing alleles: %d / %d (%.2f%%)\n",
            n_missing, total_alleles, 100 * n_missing / total_alleles))

# =============================================================================
# 5. FULL-REGION MAIN ANALYSIS
# =============================================================================
cat("\nSTEP 3: Full-region haplotypes\n")
cat("----------------------------------------------------------------------\n")

full_freq <- make_em_freq_table(a1, a2, snp_meta, min_freq = MAIN_MIN_FREQ)
full_assoc <- run_assoc_for_region("Full_Region", a1, a2, snp_meta, meta,
                                   min_copies = MIN_HAP_COPIES)

# Build a GLOBAL consistent haplotype colour palette (used by all cohort freq plots)
# so H1 is always the same colour across every figure.
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

cat(sprintf("✓ Full-region haplotypes estimated: %d total | %d shown at >= %.1f%%\n",
            nrow(full_freq$all), nrow(full_freq$shown), MAIN_MIN_FREQ * 100))
if (nrow(full_freq$shown) > 0) {
  cat(sprintf("  Top haplotype: %s (%.1f%%)\n",
              full_freq$shown$Haplotype_ID[1], full_freq$shown$Frequency[1] * 100))
}

# =============================================================================
# 6. LD ANALYSIS + BLOCKS
# =============================================================================
cat("\nSTEP 4: Linkage disequilibrium and LD blocks\n")
cat("----------------------------------------------------------------------\n")

ld_matrix <- compute_ld_matrix(a1, a2, snp_meta)
ld_blocks_by_threshold <- list()

if (RUN_BLOCK_ANALYSIS) {
  for (thr in LD_BLOCK_THRESHOLDS) {
    thr_label <- format_threshold_label(thr)
    thr_blocks <- define_ld_blocks(
      ld_matrix, snp_meta,
      threshold = thr,
      min_snps = MIN_BLOCK_SNPS,
      max_snps = MAX_BLOCK_SNPS
    )
    ld_blocks_by_threshold[[thr_label]] <- list(
      threshold = thr,
      blocks = thr_blocks
    )

    if (length(thr_blocks) == 0) {
      cat(sprintf("  No LD blocks met the configured criteria at %s\n",
                  sprintf("r^2 >= %.2f", thr)))
    } else {
      cat(sprintf("? LD blocks defined at %s: %d\n",
                  sprintf("r^2 >= %.2f", thr), length(thr_blocks)))
      for (nm in names(thr_blocks)) {
        idx <- thr_blocks[[nm]]
        cat(sprintf("  %s (%s): %d SNPs | %s\n",
                    nm, thr_label, length(idx),
                    paste(snp_meta$rsID_clean[idx], collapse = ", ")))
      }
    }
  }
}

block_freq_tables <- list()
block_assoc_tables <- list()
block_summary_rows <- list()

if (length(ld_blocks_by_threshold) > 0) {
  for (thr_label in names(ld_blocks_by_threshold)) {
    thr <- ld_blocks_by_threshold[[thr_label]]$threshold
    thr_blocks <- ld_blocks_by_threshold[[thr_label]]$blocks
    if (length(thr_blocks) == 0) {
      next
    }

    for (nm in names(thr_blocks)) {
      idx <- thr_blocks[[nm]]
      region_meta <- snp_meta[idx, , drop = FALSE]
      region_name <- sprintf("LD_Block_%s_%s", thr_label, nm)

      freq_obj <- make_em_freq_table(a1[, idx, drop = FALSE], a2[, idx, drop = FALSE],
                                     region_meta, min_freq = SUPP_MIN_FREQ)
      assoc_df <- run_assoc_for_region(region_name,
                                       a1[, idx, drop = FALSE],
                                       a2[, idx, drop = FALSE],
                                       region_meta, meta,
                                       min_copies = MIN_HAP_COPIES)

      block_freq_tables[[region_name]] <- freq_obj$shown
      block_assoc_tables[[region_name]] <- assoc_df
      summary_row <- build_region_summary(region_name, "LD_Block",
                                          region_meta, freq_obj, assoc_df)
      summary_row$LD_Threshold <- thr
      summary_row$Threshold_Label <- thr_label
      block_summary_rows[[region_name]] <- summary_row
    }
  }
}

# =============================================================================
# 7. GENE-SPECIFIC ANALYSIS
# =============================================================================
cat("\nSTEP 5: Gene-specific haplotypes\n")
cat("----------------------------------------------------------------------\n")

gene_freq_tables <- list()
gene_assoc_tables <- list()
gene_summary_rows <- list()

if (RUN_GENE_SPECIFIC) {
  genes_available <- unique(snp_meta$Gene)
  genes_to_run <- intersect(GENE_LIST_SUPP, genes_available)
  if (length(genes_to_run) == 0) {
    cat("  None of the requested genes were found in the SNP panel\n")
  } else {
    for (gene in genes_to_run) {
      gene_idx <- which(snp_meta$Gene == gene)
      if (length(gene_idx) < 2) {
        cat(sprintf("  Skipping %s: fewer than %d SNPs\n", gene, MIN_BLOCK_SNPS))
        next
      }

      gene_thresholds <- if (GENE_SPECIFIC_USE_LD_BLOCKS) LD_BLOCK_THRESHOLDS else NA_real_

      for (thr in gene_thresholds) {
        if (GENE_SPECIFIC_USE_LD_BLOCKS) {
          thr_label <- format_threshold_label(thr)
          gene_ld <- ld_matrix[gene_idx, gene_idx, drop = FALSE]
          gene_blocks <- define_ld_blocks(gene_ld, snp_meta[gene_idx, , drop = FALSE],
                                          threshold = thr,
                                          min_snps = MIN_BLOCK_SNPS,
                                          max_snps = MAX_BLOCK_SNPS)
          if (length(gene_blocks) == 0) {
            gene_blocks <- list(All = seq_along(gene_idx))
          }
        } else {
          thr_label <- "all_snps"
          gene_blocks <- list(All = seq_along(gene_idx))
        }

        for (sub_nm in names(gene_blocks)) {
          sub_local_idx <- gene_blocks[[sub_nm]]
          sub_idx <- gene_idx[sub_local_idx]
          region_meta <- snp_meta[sub_idx, , drop = FALSE]
          region_name <- if (sub_nm == "All") {
            sprintf("Gene_%s_%s_All", gene, thr_label)
          } else {
            sprintf("Gene_%s_%s_%s", gene, thr_label, sub_nm)
          }

          freq_obj <- make_em_freq_table(a1[, sub_idx, drop = FALSE],
                                         a2[, sub_idx, drop = FALSE],
                                         region_meta, min_freq = SUPP_MIN_FREQ)
          assoc_df <- run_assoc_for_region(region_name,
                                           a1[, sub_idx, drop = FALSE],
                                           a2[, sub_idx, drop = FALSE],
                                           region_meta, meta,
                                           min_copies = MIN_HAP_COPIES)

          gene_freq_tables[[region_name]] <- freq_obj$shown
          gene_assoc_tables[[region_name]] <- assoc_df
          summary_row <- build_region_summary(region_name, "Gene",
                                              region_meta, freq_obj, assoc_df)
          summary_row$LD_Threshold <- if (is.na(thr)) NA_real_ else thr
          summary_row$Threshold_Label <- thr_label
          gene_summary_rows[[region_name]] <- summary_row

          cat(sprintf("  %s: %d SNPs | %d shown haplotypes\n",
                      region_name, nrow(region_meta), nrow(freq_obj$shown)))
        }
      }
    }
  }
}

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================
cat("\nSTEP 6: Saving results workbook\n")
cat("----------------------------------------------------------------------\n")

wb <- createWorkbook()

addWorksheet(wb, "SNPs_Used")
writeData(wb, "SNPs_Used", snp_meta)

addWorksheet(wb, "Full_Global_Freqs")
writeData(wb, "Full_Global_Freqs", full_freq$shown)

if (!is.null(full_assoc) && nrow(full_assoc) > 0) {
  addWorksheet(wb, "Full_Associations")
  writeData(wb, "Full_Associations", full_assoc)
}

ld_long <- as.data.frame(as.table(ld_matrix), stringsAsFactors = FALSE)
names(ld_long) <- c("SNP1", "SNP2", "R2")
addWorksheet(wb, "LD_Matrix_Long")
writeData(wb, "LD_Matrix_Long", ld_long)

if (length(ld_blocks_by_threshold) > 0) {
  ld_block_df <- bind_rows(lapply(names(ld_blocks_by_threshold), function(thr_label) {
    thr <- ld_blocks_by_threshold[[thr_label]]$threshold
    thr_blocks <- ld_blocks_by_threshold[[thr_label]]$blocks
    if (length(thr_blocks) == 0) {
      return(NULL)
    }

    bind_rows(lapply(names(thr_blocks), function(nm) {
      idx <- thr_blocks[[nm]]
      data.frame(
        LD_Threshold = thr,
        Threshold_Label = thr_label,
        Block = nm,
        SNP_Order = seq_along(idx),
        rsID = snp_meta$rsID_clean[idx],
        POS  = snp_meta$POS_int[idx],
        Gene = snp_meta$Gene[idx],
        stringsAsFactors = FALSE
      )
    }))
  }))
  if (!is.null(ld_block_df) && nrow(ld_block_df) > 0) {
    addWorksheet(wb, "LD_Blocks")
    writeData(wb, "LD_Blocks", ld_block_df)
  }
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
if (length(gene_summary_rows) > 0) summary_parts[[length(summary_parts) + 1]] <- bind_rows(gene_summary_rows)
summary_rows <- bind_rows(summary_parts)
addWorksheet(wb, "Region_Summary")
writeData(wb, "Region_Summary", summary_rows)

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("✓ Workbook saved: %s\n", OUT_XLSX))

# =============================================================================
# 9. PLOTS
# =============================================================================
cat("\nSTEP 7: Generating figures\n")
cat("----------------------------------------------------------------------\n")
# ?? Main figure set ??????????????????????????????????????????????????????????
make_haplotype_atlas_plot(
  full_freq$all,
  snp_meta,
  out_file = file.path(MAIN_FIG_DIR, "Haplotype_Atlas_FullRegion.png")
)

make_full_region_comparison_plot(
  a1_mat = a1,
  a2_mat = a2,
  freq_df = full_freq$shown,
  metadata_df = meta,
  out_file = file.path(MAIN_FIG_DIR, "Haplotype_Cohort_Comparison_FullRegion.png")
)

make_association_overview_plot(
  summary_rows,
  out_file = file.path(MAIN_FIG_DIR, "Haplotype_Association_Overview.png")
)

make_key_block_overview_plot(
  summary_rows,
  block_assoc_tables,
  gene_assoc_tables,
  out_file = file.path(MAIN_FIG_DIR, "Haplotype_KeyBlocks_Overview.png")
)

make_haplotype_composition_plot(
  freq_df      = full_freq$shown,
  snp_meta_sub = snp_meta,
  region_name  = "Full_Region",
  out_file     = file.path(MAIN_FIG_DIR, "Haplotype_Composition_FullRegion.png")
)
if (!is.null(full_assoc) && nrow(full_assoc) > 0) {
  make_freq_comparison_plot(
    full_assoc,
    region_name = "Full_Region",
    out_file    = file.path(MAIN_FIG_DIR, "Haplotype_FreqComparison_FullRegion.png"),
    comparison  = "All"
  )
}

if (GENERATE_SUPPLEMENTARY_FIGURES) {
  if (length(ld_blocks_by_threshold) > 0) {
    for (thr_label in names(ld_blocks_by_threshold)) {
      make_ld_plot(
        ld_matrix,
        snp_meta,
        ld_blocks = ld_blocks_by_threshold[[thr_label]]$blocks,
        out_file  = file.path(SUPP_FIG_DIR, sprintf("19_LD_Heatmap_FullRegion_%s.png", thr_label))
      )
    }
  } else {
    make_ld_plot(ld_matrix, snp_meta, out_file = file.path(SUPP_FIG_DIR, "19_LD_Heatmap_FullRegion.png"))
  }

  if (length(block_summary_rows) > 0) {
    block_plot_df <- bind_rows(block_summary_rows) %>% mutate(Region = factor(Region, levels = Region[order(N_SNPs, decreasing = TRUE)]))
    p_blocks <- ggplot(block_plot_df, aes(x = Region, y = Haps_Shown)) +
      geom_col(fill = "#CC79A7", colour = "black", alpha = 0.85) +
      geom_text(aes(label = paste0(N_SNPs, " SNPs")), vjust = -0.35, size = 3.2) +
      labs(title = "Haplotype complexity across LD blocks", x = NULL, y = "Shown haplotypes") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    SAVE(file.path(SUPP_FIG_DIR, "19_LD_Block_Haplotype_Counts.png"), p_blocks, width = 8, height = 5)
  }

  if (length(gene_summary_rows) > 0) {
    gene_plot_df <- bind_rows(gene_summary_rows) %>% mutate(Region = factor(Region, levels = Region[order(Haps_Shown, decreasing = TRUE)]))
    p_genes <- ggplot(gene_plot_df, aes(x = Region, y = Haps_Shown)) +
      geom_col(fill = "#009E73", colour = "black", alpha = 0.85) +
      geom_text(aes(label = paste0(N_SNPs, " SNPs")), vjust = -0.35, size = 3.2) +
      labs(title = "Gene-specific haplotype complexity", x = NULL, y = "Shown haplotypes") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    SAVE(file.path(SUPP_FIG_DIR, "19_Gene_Haplotype_Counts.png"), p_genes, width = 8, height = 5)
  }

  if (!is.null(full_assoc) && nrow(full_assoc) > 0) {
    make_forest_plot(full_assoc, region_name = "Full_Region", out_file = file.path(SUPP_FIG_DIR, "19_Forest_FullRegion.png"))
    for (.comp in unique(full_assoc$Comparison)) {
      .comp_safe <- gsub("[^A-Za-z0-9]", "_", .comp)
      make_haplotype_dosage_plot(
        a1_mat      = a1,
        a2_mat      = a2,
        assoc_df    = full_assoc,
        metadata_df = meta,
        region_name = "Full_Region",
        out_file    = file.path(SUPP_FIG_DIR, sprintf("19_DosageComposition_FullRegion_%s.png", .comp_safe)),
        comparison  = .comp
      )
    }
    make_manhattan_hap_plot(full_assoc, snp_meta, out_file = file.path(SUPP_FIG_DIR, "19_Manhattan_Haplotypes_All.png"), comparison = "All")
  }
}

# =============================================================================
# 10. COMPLETION SUMMARY
# =============================================================================
cat("\n======================================================================\n")
cat("✓ HAPLOTYPE ANALYSIS COMPLETE\n")
cat("======================================================================\n\n")
cat(sprintf("Outputs in: %s\n\n", OUT_DIR))
cat("Key changes in this version:\n")
cat("  - Full-region haplotypes retained as the main overview\n")
cat("  - LD heatmap with block-boundary overlays\n")
cat("  - LD blocks defined with non-greedy maximal-window algorithm\n")
cat("  - Gene-specific haplotypes analysed without sliding windows\n")
cat("  - Single global BH-FDR applied per region (no per-comparison FDR)\n")
cat("  - KM plots generated for top + significantly-associated haplotypes\n\n")

# Dynamically list all output files actually created
out_files <- sort(list.files(OUT_DIR, full.names = FALSE))
cat(sprintf("Files created (%d total):\n", length(out_files)))
for (f in out_files) cat(sprintf("  %s\n", f))

cat("\nRegion summary:\n")
print(summary_rows, row.names = FALSE)
cat("\n")
