#!/usr/bin/env Rscript
# =============================================================================
# Script 15: GSDMB Haplotype Analysis 
# =============================================================================

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
  for (pkg in c("patchwork", "logistf", "genetics")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(patchwork)
  library(logistf)
  library(genetics)   # for pairwise r² (LD); no LDheatmap needed
})

# Explicitly re-attach dplyr verbs so they win over genetics/MASS masking
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
left_join  <- dplyr::left_join
inner_join <- dplyr::inner_join
bind_rows  <- dplyr::bind_rows

# =============================================================================
# 0. GLOBAL THEME  (publication-ready: Nature/Genetics style)
# =============================================================================
pub_theme <- theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    plot.title        = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle     = element_text(size = 9,  color = "grey40", hjust = 0),
    axis.title        = element_text(size = 10, face = "bold"),
    axis.text         = element_text(size = 9,  color = "black"),
    axis.line         = element_line(linewidth = 0.4),
    axis.ticks        = element_line(linewidth = 0.4),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey95", color = "grey70"),
    strip.text        = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 8),
    legend.title      = element_text(size = 9, face = "bold"),
    legend.key.size   = unit(0.45, "cm"),
    plot.margin       = margin(8, 10, 6, 8)
  )

theme_set(pub_theme)

# Cohort/group colour palette (accessible, print-safe)
GROUP_COLS <- c(
  "Breast_Healthy"      = "#4393C3",
  "Breast_Tumour"       = "#D6604D",
  "Endometrium_Healthy" = "#74ADD1",
  "Endometrium_Tumour"  = "#A50026",
  "All_Healthy"         = "#92C5DE",
  "All_Tumour"          = "#F4A582"
)
COHORT_COLS  <- c("Breast" = "#4393C3", "Endometrium" = "#D6604D", "All" = "#878787")
TISSUE_COLS  <- c("Healthy" = "#4393C3", "Tumour" = "#D6604D")

SAVE <- function(file, plot, width = 10, height = 6) {
  ggsave(file, plot, width = width, height = height, dpi = 300,
         units = "in", bg = "white")
  cat("  Saved:", basename(file), "\n")
}

# =============================================================================
# 1. CONFIGURATION & PATHS
# =============================================================================
BASE_DIR  <- "/home/gadeaalonsoj/tfm/gsdmb_final_results"
HAPLO_DIR <- file.path(BASE_DIR, "19_haplotype_phased")
OUT_DIR   <- file.path(BASE_DIR, "19_haplo_stats_results")

GENO_FILE  <- file.path(HAPLO_DIR, "phased_genotypes.tsv")
META_FILE  <- file.path(HAPLO_DIR, "sample_metadata.tsv")
WHITELIST  <- file.path(BASE_DIR, "11_Master_Unique_SNP_Summary.xlsx")
ANNOTATED  <- file.path(BASE_DIR, "GSDMB_Annotated_Report_Fixed.xlsx")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# =============================================================================
# 2. DATA LOADING & ANNOTATION
# =============================================================================
cat("\nSTEP 2: Loading and annotating SNPs...\n")

target_rsids <- read_excel(WHITELIST) %>% pull(Variant_ID) %>% trimws()

annotation_map <- read_excel(ANNOTATED, sheet = "Biological_Annotations") %>%
  mutate(
    rsID_clean = str_extract(Existing_variation, "rs[0-9]+"),
    POS_int    = as.integer(gsub("[^0-9]", "", as.character(POS)))
  ) %>%
  filter(rsID_clean %in% target_rsids) %>%
  select(POS_int, rsID_clean, Gene = SYMBOL, REF, ALT) %>%
  arrange(rsID_clean, desc(Gene)) %>%
  distinct(rsID_clean, .keep_all = TRUE) %>%
  arrange(POS_int)

geno_raw <- read.table(GENO_FILE, header = TRUE, sep = "\t", check.names = FALSE)
metadata <- read.table(META_FILE,  header = TRUE, sep = "\t", check.names = FALSE)
common_samples <- intersect(colnames(geno_raw), metadata$Sample)

stopifnot(all(c("Sample", "Cohort", "Tissue") %in% colnames(metadata)))

geno_filtered <- geno_raw %>%
  mutate(extracted_pos = as.integer(str_split_fixed(ID, ":", 4)[, 2])) %>%
  inner_join(annotation_map, by = c("extracted_pos" = "POS_int"))

n_snps     <- nrow(geno_filtered)
snp_labels <- geno_filtered$rsID_clean

# Carrier counts for annotation labels
annotation_map$Carrier_Count <- sapply(annotation_map$rsID_clean, function(rsid) {
  row_data <- geno_filtered %>% filter(rsID_clean == rsid)
  sum(grepl("1", unlist(row_data[, common_samples])))
})

annotation_map <- annotation_map %>%
  mutate(Detailed_Label = paste0(
    rsID_clean, "\n(", Gene, ", ", REF, ">", ALT, ")\n[n=", Carrier_Count, "]"
  ))

ordered_labels <- annotation_map$Detailed_Label

# =============================================================================
# 3. BUILD GENOTYPE MATRIX (all common samples)
# =============================================================================
cat("STEP 3: Building genotype matrix...\n")

mat_all <- matrix(NA, nrow = length(common_samples), ncol = n_snps * 2,
                  dimnames = list(common_samples, NULL))

for (i in seq_along(common_samples)) {
  alleles <- as.numeric(unlist(strsplit(
    as.character(geno_filtered[, common_samples[i]]), "[|/]"
  )))
  mat_all[i, ] <- alleles + 1   # haplo.stats requires 1-based coding
}

# =============================================================================
# 4. GLOBAL HAPLOTYPE EM  + 1% FREQUENCY FILTER
# =============================================================================
cat("STEP 4: Global EM and frequency filtering (≥1%)...\n")

geno_setup_all <- setupGeno(mat_all, locus.label = snp_labels)
global_em      <- haplo.em(geno_setup_all)

freq_df <- as.data.frame(global_em$haplotype)
colnames(freq_df) <- snp_labels
freq_df$em_index  <- seq_len(nrow(freq_df))
freq_df$Frequency <- as.numeric(global_em$hap.prob)

freq_df <- freq_df %>%
  filter(Frequency >= 0.01) %>%
  arrange(desc(Frequency)) %>%
  mutate(
    Haplotype = paste0("H", seq_len(n())),
    Haplotype = factor(Haplotype, levels = Haplotype),
    Allele_String = apply(
      across(all_of(snp_labels)), 1,
      function(x) paste(ifelse(x == 1, "Ref", "Alt"), collapse = "-")
    )
  )

em_to_haplabel    <- setNames(as.character(freq_df$Haplotype), freq_df$em_index)
common_haplotypes <- as.character(freq_df$Haplotype)
cat(sprintf("  Haplotypes passing ≥1%% filter: %d\n", length(common_haplotypes)))

# =============================================================================
# REFERENCE HAPLOTYPE SELECTION
# H1 (most frequent) is used as the statistical reference for all GLM
# comparisons. This maximises power by providing the largest reference group,
# giving the most precise odds ratio estimates.
#
# NOTE ON ANCESTRAL STATE: Direct comparison to the reference genome is not
# possible because no sample in the dataset carries the all-reference haplotype
# on both chromosomes with certainty. As a transparency measure, the script
# identifies and reports the haplotype with the fewest alt alleles (closest
# empirical proxy for the reference genome sequence) so it can be noted in the
# methods. This is NOT used as the statistical reference.
# =============================================================================
freq_df <- freq_df %>%
  mutate(
    across(all_of(snp_labels), ~ as.numeric(as.character(.))),
    n_alt_alleles = rowSums(across(all_of(snp_labels), ~ . - 1))
  )

# H1 = most frequent = statistical reference
ref_haplo_global <- as.character(freq_df$Haplotype[1])

# Identify all-reference proxy for reporting purposes only
allref_candidates <- freq_df %>% filter(n_alt_alleles == min(n_alt_alleles))
allref_row        <- allref_candidates %>% arrange(desc(Frequency)) %>% slice(1)
allref_haplo      <- as.character(allref_row$Haplotype)

cat(sprintf("  Statistical reference haplotype : %s  (most frequent, freq=%.1f%%)\n",
            ref_haplo_global, freq_df$Frequency[1] * 100))
cat(sprintf("  All-ref proxy haplotype (report): %s  (n_alt=%d, freq=%.1f%%) — noted in methods, not used as reference\n",
            allref_haplo, min(freq_df$n_alt_alleles), allref_row$Frequency * 100))

# =============================================================================
# 5. STATISTICAL MODELS  — Tumour vs Healthy within each comparison group
# =============================================================================
cat("STEP 5: Running association models...\n")

comparisons <- list(
  list(label = "Breast",      cohorts = "Breast",                  name = "Breast"),
  list(label = "Endometrium", cohorts = "Endometrium",             name = "Endometrium"),
  list(label = "All",         cohorts = c("Breast","Endometrium"), name = "All")
)

run_comparison <- function(comp) {
  sub_meta <- metadata %>%
    filter(Cohort %in% comp$cohorts, Sample %in% common_samples,
           Tissue %in% c("Healthy", "Tumour"))

  n_normal <- sum(sub_meta$Tissue == "Healthy")
  n_tumour <- sum(sub_meta$Tissue == "Tumour")

  cat(sprintf("  %s: n=%d  (Healthy=%d, Tumour=%d)\n",
              comp$label, nrow(sub_meta), n_normal, n_tumour))

  if (nrow(sub_meta) < 10) {
    warning(sprintf("  SKIPPING %s: too few total samples (%d).", comp$label, nrow(sub_meta)))
    return(NULL)
  }
  if (n_normal == 0 || n_tumour == 0) {
    warning(sprintf(
      "  SKIPPING %s: both Healthy and Tumour required (Healthy=%d, Tumour=%d).",
      comp$label, n_normal, n_tumour))
    return(NULL)
  }

  sub_samples <- sub_meta$Sample
  sub_mat     <- mat_all[sub_samples, , drop = FALSE]
  y           <- as.numeric(sub_meta$Tissue == "Tumour")

  geno_sub <- setupGeno(sub_mat, locus.label = snp_labels)

  sc <- tryCatch(
    haplo.score(y, geno_sub, trait.type = "binomial",
                min.count = ceiling(0.01 * length(y))),
    error = function(e) { message("  Score test failed: ", e$message); NULL }
  )

  n_haps_g     <- length(global_em$hap.prob)
  hap_labels_g <- em_to_haplabel[as.character(seq_len(n_haps_g))]
  global_to_local <- match(common_samples, sub_samples)

  dosage_mat <- matrix(0, nrow = nrow(sub_mat), ncol = length(common_haplotypes),
                       dimnames = list(sub_samples, common_haplotypes))

  em_row_g  <- global_em$indx.subj
  em_hap1_g <- global_em$hap1code
  em_hap2_g <- global_em$hap2code
  em_post_g <- global_em$post

  for (k in seq_along(em_row_g)) {
    global_i <- em_row_g[k]
    local_i  <- global_to_local[global_i]
    if (is.na(local_i)) next
    h1   <- em_hap1_g[k]; h2 <- em_hap2_g[k]; post <- em_post_g[k]
    lbl1 <- if (!is.na(h1) && h1 >= 1 && h1 <= n_haps_g) hap_labels_g[h1] else NA
    lbl2 <- if (!is.na(h2) && h2 >= 1 && h2 <= n_haps_g) hap_labels_g[h2] else NA
    if (!is.na(lbl1) && lbl1 %in% common_haplotypes)
      dosage_mat[local_i, lbl1] <- dosage_mat[local_i, lbl1] + post
    if (!is.na(lbl2) && lbl2 %in% common_haplotypes)
      dosage_mat[local_i, lbl2] <- dosage_mat[local_i, lbl2] + post
  }

  ref_haplo_label <- ref_haplo_global
  test_haplos     <- setdiff(common_haplotypes, ref_haplo_label)
  if (length(test_haplos) == 0) return(NULL)

  results <- lapply(test_haplos, function(hname) {
    df_fit <- data.frame(y = y, hap = dosage_mat[, hname])
    if (var(df_fit$hap) == 0) return(NULL)
    fit <- tryCatch(
      logistf(y ~ hap, data = df_fit, firth = TRUE, pl = FALSE,
              control = logistf.control(maxit = 500, maxstep = 10)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    beta <- coef(fit)["hap"]
    se   <- sqrt(diag(vcov(fit)))["hap"]
    p    <- fit$prob["hap"]
    data.frame(
      Haplotype  = hname, Beta = beta, SE = se, p_glm = p,
      OR = exp(beta), Lower = exp(beta - 1.96 * se),
      Upper = exp(beta + 1.96 * se), Comparison = comp$label,
      row.names = NULL
    )
  })

  res <- bind_rows(results)
  if (nrow(res) == 0) return(NULL)

  if (!is.null(sc)) {
    score_p <- data.frame(
      Haplotype = paste0("H", seq_along(sc$score.haplo.p)),
      Score_P   = sc$score.haplo.p
    )
    res <- left_join(res, score_p, by = "Haplotype")
  } else {
    res$Score_P <- NA_real_
  }
  res
}

all_stats_list <- lapply(comparisons, run_comparison)
all_stats_raw  <- bind_rows(all_stats_list)
cat(sprintf("  Rows returned by GLM: %d\n", nrow(all_stats_raw)))

# =============================================================================
# 5b. BUILD all_stats WITH BH-FDR SIGNIFICANCE LABELS
# BH correction applied within each comparison group separately.
# sig_label reflects FDR-adjusted p, not raw Score_P.
# =============================================================================
if (nrow(all_stats_raw) > 0) {
  all_stats <- all_stats_raw %>%
    mutate(
      Comparison = factor(Comparison, levels = c("Breast", "Endometrium", "All")),
      Haplotype  = factor(Haplotype,  levels = common_haplotypes)
    ) %>%
    group_by(Comparison) %>%
    mutate(
      fdr = p.adjust(Score_P, method = "BH"),
      sig_label = case_when(
        !is.na(fdr) & fdr < 0.001 ~ "***",
        !is.na(fdr) & fdr < 0.01  ~ "**",
        !is.na(fdr) & fdr < 0.05  ~ "*",
        TRUE                       ~ ""
      )
    ) %>%
    ungroup()
  cat(sprintf("  Comparisons with results: %s\n",
              paste(unique(as.character(all_stats$Comparison)), collapse = ", ")))
  print(head(all_stats %>% select(Haplotype, Comparison, OR, Score_P, fdr), 6))
} else {
  all_stats <- all_stats_raw
  cat("  WARNING: No GLM results.\n")
  print(warnings())
}

# =============================================================================
# STEP 5c. EXTRACT PHASED HAPLOTYPE ALLELE MATRICES FOR LD COMPUTATION
#
# mat_all rows = samples, cols interleaved: SNP1_allele1, SNP1_allele2,
# SNP2_allele1, SNP2_allele2, ... (1-based: 1=ref, 2=alt)
#
# For LD we treat each phased chromosome as an independent observation.
# This is the correct approach — using phased data from BEAGLE gives accurate
# r² because phase ambiguity (the main cause of LD deflation from unphased
# dosage) is already resolved.
# =============================================================================

# Extract the two haploid allele matrices from mat_all
# allele1_mat[i, j] = first phased allele of sample i at SNP j (1 or 2)
# allele2_mat[i, j] = second phased allele of sample i at SNP j (1 or 2)
allele1_mat <- mat_all[, seq(1, n_snps * 2, by = 2), drop = FALSE]  # odd cols
allele2_mat <- mat_all[, seq(2, n_snps * 2, by = 2), drop = FALSE]  # even cols
colnames(allele1_mat) <- snp_labels
colnames(allele2_mat) <- snp_labels

# Stack: 2*n_samples haploid rows x n_snps cols (0=ref, 1=alt, 0-based)
# Each row is one phased chromosome — fully resolved by BEAGLE
hap_mat <- rbind(allele1_mat - 1L, allele2_mat - 1L)  # convert to 0/1

# =============================================================================
# 6. FIGURE 0 — LD Heatmap (r² between SNPs in the haplotype analysis)
#
# WHY: Shows the haplotype block structure. High r² between SNPs explains
# why they co-segregate into consistent haplotypes. This is the biological
# justification for doing haplotype analysis at all rather than single-SNP.
# Computed from BEAGLE-phased haplotypes for accurate r² (no phase ambiguity).
# =============================================================================
cat("STEP 6: Figure 0 — LD Heatmap...\n")

tryCatch({

  # ── Step 1: compute r² directly from phased haplotype matrix ─────────────
  # r²(A,B) = [cov(A,B)]² / [var(A) * var(B)]
  # where A and B are 0/1 vectors across all 2*n_samples phased chromosomes.
  # This is exact and requires no external LD package.
  compute_r2_phased <- function(hap_mat) {
    n_snps_h <- ncol(hap_mat)
    r2_mat   <- matrix(NA_real_, nrow = n_snps_h, ncol = n_snps_h,
                       dimnames = list(colnames(hap_mat), colnames(hap_mat)))
    # allele frequencies (p = freq of alt allele = 1)
    p <- colMeans(hap_mat, na.rm = TRUE)
    for (i in seq_len(n_snps_h)) {
      for (j in seq(i, n_snps_h)) {
        if (i == j) { r2_mat[i, j] <- 1; next }
        # only use chromosomes with no missing data at either SNP
        ok  <- !is.na(hap_mat[, i]) & !is.na(hap_mat[, j])
        a   <- hap_mat[ok, i]
        b   <- hap_mat[ok, j]
        pa  <- mean(a); pb <- mean(b)
        if (pa == 0 | pa == 1 | pb == 0 | pb == 1) { r2_mat[i,j] <- r2_mat[j,i] <- NA; next }
        # D = observed AB haplotype frequency minus expected
        D   <- mean(a * b) - pa * pb
        r2  <- D^2 / (pa * (1 - pa) * pb * (1 - pb))
        r2_mat[i, j] <- r2_mat[j, i] <- round(r2, 4)
      }
    }
    r2_mat
  }

  r2_mat <- compute_r2_phased(hap_mat)
  cat(sprintf("  r² matrix computed: %d x %d SNPs, %d phased chromosomes\n",
              ncol(r2_mat), ncol(r2_mat), nrow(hap_mat)))
  cat(sprintf("  r² range (off-diagonal): %.3f – %.3f\n",
              min(r2_mat[row(r2_mat) != col(r2_mat)], na.rm = TRUE),
              max(r2_mat[row(r2_mat) != col(r2_mat)], na.rm = TRUE)))

  # ── Step 3: axis labels — rsID + gene + position ──────────────────────────
  ax_labels <- annotation_map %>%
    mutate(label = paste0(rsID_clean, "\n", Gene,
                          " (", round(POS_int / 1e3, 1), " kb)")) %>%
    pull(label, name = rsID_clean)

  # ── Step 4: tidy long-format, LOWER triangle + diagonal ───────────────────
  # Lower triangle keeps genomic order on both axes; easy to read
  r2_df <- as.data.frame(as.table(r2_mat), stringsAsFactors = FALSE) %>%
    rename(SNP1 = Var1, SNP2 = Var2, r2 = Freq) %>%
    mutate(
      SNP1  = factor(SNP1, levels = snp_labels),
      SNP2  = factor(SNP2, levels = snp_labels),
      r2    = as.numeric(r2),
      keep  = as.integer(SNP1) >= as.integer(SNP2),
      # show value if >= 0.10; blank for diagonal (self-comparisons)
      r2_label = case_when(
        SNP1 == SNP2 ~ "",
        r2   >= 0.10 ~ sprintf("%.2f", r2),
        TRUE         ~ ""
      ),
      # white text on dark cells, black on light cells
      txt_col = ifelse(r2 >= 0.55, "white", "black")
    ) %>%
    filter(keep)

  # ── Step 5: colour palette ─────────────────────────────────────────────────
  ld_pal <- colorRampPalette(
    c("#FFFFFF", "#FFF0E0", "#FDCC8A", "#FC8D59", "#E34A33", "#8B0000")
  )(101)

  # ── Step 6: square heatmap ────────────────────────────────────────────────
  fig0 <- ggplot(r2_df, aes(x = SNP2, y = SNP1, fill = r2)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = r2_label, colour = txt_col),
              size = 2.6, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradientn(
      colours = ld_pal,
      limits  = c(0, 1),
      breaks  = c(0, 0.25, 0.5, 0.75, 1),
      labels  = c("0", "0.25", "0.50", "0.75", "1.00"),
      name    = expression(r^2),
      guide   = guide_colorbar(
        barheight = 8, barwidth = 1.2,
        title.position = "top", title.hjust = 0.5,
        frame.colour = "grey40", ticks.colour = "grey40"
      )
    ) +
    scale_x_discrete(labels = ax_labels[snp_labels], position = "bottom") +
    scale_y_discrete(labels = ax_labels[snp_labels], limits = rev) +
    coord_fixed() +
    labs(
      title    = "Linkage Disequilibrium Between SNPs",
      subtitle = sprintf(
        "Pairwise r² from BEAGLE-phased haplotypes | %d samples (%d chromosomes) | Lower triangle | r² ≥ 0.10 labelled",
        length(common_samples), nrow(hap_mat)),
      x = NULL, y = NULL
    ) +
    theme(
      axis.text.x     = element_text(size = 7.5, angle = 45, hjust = 1,
                                     lineheight = 0.85),
      axis.text.y     = element_text(size = 7.5, lineheight = 0.85),
      axis.ticks      = element_blank(),
      axis.line       = element_blank(),
      panel.grid      = element_blank(),
      legend.position = "right",
      plot.margin     = margin(6, 10, 6, 8)
    )

  fig_sz <- max(7, n_snps * 0.75 + 2)
  SAVE(file.path(OUT_DIR, "19_Fig0_LD_Heatmap.png"), fig0,
       width = fig_sz, height = fig_sz - 1)

}, error = function(e) {
  cat(sprintf("  WARNING: LD heatmap failed (%s) — skipping.\n", e$message))
})

# =============================================================================
# 7. FIGURE 1 — Haplotype Composition
#    1a: SNP Alt-Allele Carrier Counts (bar) — QC/methods context
#    1b: Haplotype × SNP allele tile — the Rosetta Stone of the analysis
#        IMPROVED: added allele-string summary column on right margin
# =============================================================================
cat("STEP 7: Figure 1 — Haplotype Composition...\n")

tile_data <- freq_df %>%
  select(Haplotype, Frequency, Allele_String, all_of(snp_labels)) %>%
  pivot_longer(cols = all_of(snp_labels), names_to = "rsID_clean", values_to = "Val") %>%
  left_join(annotation_map, by = "rsID_clean") %>%
  mutate(
    Allele     = factor(as.numeric(Val) - 1, levels = c(0, 1), labels = c("Ref", "Alt")),
    Haplotype  = factor(Haplotype, levels = rev(common_haplotypes)),
    Freq_Label = percent(Frequency, accuracy = 0.1)
  )

# ── Fig 1a: SNP carrier count bar ─────────────────────────────────────────────
fig1a <- annotation_map %>%
  mutate(
    Short_Label = paste0(rsID_clean, "\n(", Gene, ")"),
    Short_Label = factor(Short_Label, levels = annotation_map %>%
                           mutate(Short_Label = paste0(rsID_clean, "\n(", Gene, ")")) %>%
                           pull(Short_Label))
  ) %>%
  ggplot(aes(x = Short_Label, y = Carrier_Count)) +
  geom_col(fill = "#2166AC", width = 0.65) +
  geom_text(aes(label = Carrier_Count), vjust = -0.4, size = 3,
            fontface = "bold", color = "#2166AC") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     breaks = scales::pretty_breaks(5)) +
  labs(
    title    = "Alt-Allele Carrier Counts per SNP",
    subtitle = "Number of samples carrying ≥1 alternate allele at each variant",
    x        = "Variant (rsID, Gene)", y = "Carriers (n)"
  ) +
  theme(
    axis.text.x        = element_text(angle = 40, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3)
  )

SAVE(file.path(OUT_DIR, "19_Fig1a_SNP_Carrier_Counts.png"), fig1a, width = 12, height = 5)

# ── Fig 1b: Haplotype allele composition tile ─────────────────────────────────
# IMPROVED: allele-string summary shown on right margin so each haplotype
# can be described in a single line (e.g. "Ref-Alt-Ref-Ref-Alt")
tile_data2 <- tile_data %>%
  mutate(
    Compact_Label = paste0(rsID_clean, "\n[", Gene, ", n=", Carrier_Count, "]"),
    Compact_Label = factor(Compact_Label,
      levels = annotation_map %>%
        mutate(l = paste0(rsID_clean, "\n[", Gene, ", n=", Carrier_Count, "]")) %>%
        pull(l))
  )

ordered_compact <- levels(tile_data2$Compact_Label)
n_cols <- length(ordered_compact)

fig1b <- tile_data2 %>%
  ggplot(aes(x = Compact_Label, y = Haplotype, fill = Allele)) +
  geom_tile(color = "white", linewidth = 0.45) +
  # Frequency label
  geom_text(
    data = tile_data2 %>% distinct(Haplotype, Frequency, Freq_Label),
    aes(x = n_cols + 0.9, y = Haplotype, label = Freq_Label, fill = NULL),
    size = 2.7, hjust = 0, color = "grey25", fontface = "bold"
  ) +
  # Allele string summary (e.g. Ref-Alt-Ref…) — compact haplotype fingerprint
  geom_text(
    data = tile_data2 %>% distinct(Haplotype, Allele_String),
    aes(x = n_cols + 2.6, y = Haplotype, label = Allele_String, fill = NULL),
    size = 2.2, hjust = 0, color = "grey40", family = "mono"
  ) +
  scale_fill_manual(values = c("Ref" = "#E8E8E8", "Alt" = "#C0392B"), name = "Allele") +
  scale_x_discrete(expand = expansion(add = c(0.5, 5.5))) +
  labs(
    title    = "GSDMB Haplotype Allele Composition",
    subtitle = sprintf(
      "Each row = one haplotype (H1 most frequent); red = alternate allele | %d haplotypes ≥1%% global freq.",
      length(common_haplotypes)),
    x = "Variant (rsID, Gene, n carriers)", y = "Haplotype"
  ) +
  theme(
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 7.5),
    panel.grid      = element_blank(),
    legend.position = "right",
    plot.margin     = margin(6, 14, 6, 8)
  )

SAVE(file.path(OUT_DIR, "19_Fig1b_Haplotype_Composition.png"), fig1b, width = 14, height = 8)

fig1_combined <- fig1a / fig1b + plot_layout(heights = c(1, 2.8))
SAVE(file.path(OUT_DIR, "19_Fig1_Haplotype_Composition.png"), fig1_combined, width = 14, height = 14)

# =============================================================================
# 8. FIGURE 2 — Global Haplotype Frequencies  (SIMPLIFIED)
#    CHANGED: removed per-cohort lollipop panels (redundant with Fig 3).
#    Now a single clean lollipop with frequency % and estimated count labels.
#    Points sized by frequency; reference haplotype H1 highlighted in red.
# =============================================================================
cat("STEP 8: Figure 2 — Global Haplotype Frequencies...\n")

hap_order <- rev(common_haplotypes)   # H1 at top when coord_flip applied

fig2 <- freq_df %>%
  mutate(
    Haplotype  = factor(as.character(Haplotype), levels = hap_order),
    n_chromos  = round(Frequency * length(common_samples) * 2),
    is_ref     = as.character(Haplotype) == ref_haplo_global,
    point_col  = ifelse(is_ref, "#e74c3c", "#2166AC"),
    freq_label = sprintf("%s  (n≈%d)", percent(Frequency, accuracy = 0.1), n_chromos)
  ) %>%
  ggplot(aes(x = Haplotype, y = Frequency)) +
  geom_segment(aes(xend = Haplotype, y = 0, yend = Frequency, colour = is_ref),
               linewidth = 0.9) +
  geom_point(aes(size = Frequency, colour = is_ref), alpha = 0.9) +
  geom_text(aes(label = freq_label), hjust = -0.12, size = 2.7, color = "grey25") +
  scale_colour_manual(values = c("FALSE" = "#2166AC", "TRUE" = "#e74c3c"),
                      labels = c("FALSE" = "Other", "TRUE" = "Reference (H1)"),
                      name = NULL) +
  scale_size_continuous(range = c(2, 9), guide = "none") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.5)),
                     breaks = scales::pretty_breaks(5)) +
  coord_flip() +
  labs(
    title    = "Global GSDMB Haplotype Frequencies",
    subtitle = sprintf(
      "EM-estimated frequencies across all %d samples (%d chromosomes)  |  %s = statistical reference (most frequent, shown in red)  |  %s = all-ref proxy (fewest alt alleles)",
      length(common_samples), length(common_samples) * 2, ref_haplo_global, allref_haplo),
    x = "Haplotype", y = "Global Frequency"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position    = "top",
    legend.justification = "left"
  )

SAVE(file.path(OUT_DIR, "19_Fig2_Global_Frequencies.png"), fig2, width = 10, height = 6.5)

# =============================================================================
# 9. FIGURE 3 — Haplotype Frequencies: Healthy vs Tumour
#    CHANGED from grouped bars to paired dot-line plot.
#    Each haplotype has two dots (Healthy, Tumour) connected by a line.
#    Line colour encodes direction: red = higher in tumour, blue = lower.
#    Faceted by Breast | Endometrium | All (pooled).
#    Much easier to read frequency shifts than grouped bars.
# =============================================================================
cat("STEP 9: Figure 3 — Stratified Frequencies (Healthy vs Tumour)...\n")

compute_stratum_freqs <- function(tissue_val, cohort_val, label) {
  sub_meta <- metadata %>%
    filter(
      Sample %in% common_samples,
      Tissue == tissue_val,
      if (!is.null(cohort_val)) Cohort %in% cohort_val else TRUE
    )
  if (nrow(sub_meta) < 5) return(NULL)
  sub_mat <- mat_all[sub_meta$Sample, , drop = FALSE]
  em <- tryCatch(
    haplo.em(setupGeno(sub_mat, locus.label = snp_labels)),
    error = function(e) NULL
  )
  if (is.null(em)) return(NULL)
  hap_df <- as.data.frame(em$haplotype)
  colnames(hap_df) <- snp_labels
  global_patterns  <- apply(freq_df[, snp_labels], 1, paste, collapse = "-")
  stratum_patterns <- apply(hap_df[, snp_labels],  1, paste, collapse = "-")
  hap_df$Haplotype <- as.character(freq_df$Haplotype)[match(stratum_patterns, global_patterns)]
  hap_df$Frequency <- as.numeric(em$hap.prob)
  hap_df$Group     <- label
  hap_df %>% filter(Haplotype %in% common_haplotypes) %>%
    select(Haplotype, Frequency, Group)
}

strat_freqs <- bind_rows(
  compute_stratum_freqs("Healthy", "Breast",                  "Breast_Healthy"),
  compute_stratum_freqs("Tumour",  "Breast",                  "Breast_Tumour"),
  compute_stratum_freqs("Healthy", "Endometrium",             "Endometrium_Healthy"),
  compute_stratum_freqs("Tumour",  "Endometrium",             "Endometrium_Tumour"),
  compute_stratum_freqs("Healthy", c("Breast","Endometrium"), "All_Healthy"),
  compute_stratum_freqs("Tumour",  c("Breast","Endometrium"), "All_Tumour")
) %>%
  mutate(
    Haplotype = factor(Haplotype, levels = common_haplotypes),
    Group     = factor(Group, levels = names(GROUP_COLS)),
    Tissue    = ifelse(grepl("Healthy", Group), "Healthy", "Tumour"),
    Cohort    = case_when(
      grepl("^Breast",      Group) ~ "Breast",
      grepl("^Endometrium", Group) ~ "Endometrium",
      TRUE                         ~ "All (Pooled)"
    ),
    Cohort    = factor(Cohort, levels = c("Breast", "Endometrium", "All (Pooled)"))
  )

# Wide format for drawing connecting lines between Healthy and Tumour dots
strat_wide <- strat_freqs %>%
  select(Haplotype, Tissue, Frequency, Cohort) %>%
  pivot_wider(names_from = Tissue, values_from = Frequency) %>%
  mutate(
    direction = case_when(
      Tumour > Healthy ~ "Higher in Tumour",
      Tumour < Healthy ~ "Lower in Tumour",
      TRUE             ~ "Equal"
    ),
    delta = Tumour - Healthy
  )

# Long format for dots
strat_long <- strat_freqs %>%
  mutate(Tissue = factor(Tissue, levels = c("Healthy", "Tumour")))

fig3 <- ggplot() +
  # Connecting lines coloured by direction
  geom_segment(
    data = strat_wide,
    aes(x = Haplotype, xend = Haplotype,
        y = Healthy,   yend = Tumour,
        colour = direction),
    linewidth = 0.7, alpha = 0.7
  ) +
  # Dots for each tissue
  geom_point(
    data = strat_long,
    aes(x = Haplotype, y = Frequency,
        shape = Tissue, fill = Tissue),
    size = 2.8, stroke = 0.4, colour = "white", alpha = 0.95
  ) +
  scale_colour_manual(
    values = c("Higher in Tumour" = "#D6604D",
               "Lower in Tumour"  = "#4393C3",
               "Equal"            = "grey60"),
    name = "Frequency shift"
  ) +
  scale_shape_manual(values = c("Healthy" = 21, "Tumour" = 24), name = "Tissue") +
  scale_fill_manual(values = TISSUE_COLS, name = "Tissue") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0.05, 0.1))) +
  facet_wrap(~ Cohort, ncol = 1, scales = "free_y") +
  labs(
    title    = "Haplotype Frequencies: Healthy vs Tumour",
    subtitle = "EM-estimated per stratum  |  Lines connect Healthy → Tumour for each haplotype  |  Colour = direction of shift",
    x        = "Haplotype", y = "EM-Estimated Frequency"
  ) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.box      = "horizontal"
  )

SAVE(file.path(OUT_DIR, "19_Fig3_Stratified_Frequencies.png"), fig3, width = 11, height = 12)

# =============================================================================
# 10. FIGURE 4 — Forest Plot (Odds Ratios, Tumour vs Healthy)
#     All v3 improvements:
#       • Left frequency strip (lollipop; H1 in red)
#       • H1 reference plotted as diamond at OR=1
#       • Stars on point (not CI end), BH-FDR corrected
#       • Arrow heads for clipped upper CIs
#       • Y-axis labels only on left strip
# =============================================================================
cat("STEP 10: Figure 4 — Forest Plot...\n")

ref_haplo <- ref_haplo_global
ref_freq  <- percent(freq_df$Frequency[freq_df$Haplotype == ref_haplo_global], accuracy = 0.1)

if (nrow(all_stats) == 0) {
  cat("  SKIPPING Fig4: no GLM results available.\n")
} else {

  lo_vals <- all_stats$Lower[is.finite(all_stats$Lower) & all_stats$Lower > 0]
  hi_vals <- all_stats$Upper[is.finite(all_stats$Upper) & all_stats$Upper > 0]
  x_min   <- max(0.05, min(lo_vals, na.rm = TRUE) * 0.7)
  x_max   <- min(50,   max(hi_vals, na.rm = TRUE) * 1.4)

  all_breaks <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32)
  use_breaks <- all_breaks[all_breaks >= x_min & all_breaks <= x_max]
  if (!1 %in% use_breaks) use_breaks <- sort(c(1, use_breaks))

  all_hap_levels <- unique(c(ref_haplo, rev(levels(all_stats$Haplotype))))

  ref_row <- expand.grid(
    Haplotype  = ref_haplo,
    Comparison = levels(all_stats$Comparison),
    stringsAsFactors = FALSE
  ) %>%
    mutate(OR = 1, Lower = NA_real_, Upper = NA_real_,
           sig_label = "", fdr = NA_real_, is_ref = TRUE)

  plot_data <- all_stats %>%
    mutate(is_ref = FALSE) %>%
    bind_rows(ref_row) %>%
    mutate(
      Haplotype  = factor(Haplotype,  levels = rev(all_hap_levels)),
      Comparison = factor(Comparison, levels = c("Breast", "Endometrium", "All")),
      ci_clipped = !is_ref & !is.na(Upper) & Upper > x_max
    )

  clipped_data <- plot_data %>% filter(ci_clipped)

  forest_main <- ggplot(plot_data,
    aes(x = OR, y = Haplotype, colour = Comparison, shape = Comparison)) +

    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey45", linewidth = 0.5) +

    geom_errorbarh(
      data = plot_data %>% filter(!is_ref, !ci_clipped, !is.na(Lower), !is.na(Upper)),
      aes(xmin = pmax(Lower, x_min * 0.9), xmax = pmin(Upper, x_max * 1.05)),
      height = 0.28, linewidth = 0.55, alpha = 0.85
    ) +

    geom_errorbarh(
      data = plot_data %>% filter(!is_ref, ci_clipped, !is.na(Lower), !is.na(Upper)),
      aes(xmin = pmax(Lower, x_min * 0.9), xmax = x_max * 0.98),
      height = 0.28, linewidth = 0.55, alpha = 0.85
    ) +

    geom_segment(
      data = clipped_data,
      aes(x = x_max * 0.96, xend = x_max * 1.04,
          y = Haplotype,    yend = Haplotype),
      arrow     = arrow(length = unit(0.12, "cm"), type = "open"),
      linewidth = 0.55, alpha = 0.85
    ) +

    geom_point(data = plot_data %>% filter(!is_ref), size = 2.4, stroke = 0.5) +

    geom_point(
      data = plot_data %>% filter(is_ref),
      aes(x = 1, y = Haplotype), shape = 18, size = 3.8,
      colour = "black", inherit.aes = FALSE
    ) +

    geom_text(
      data = plot_data %>% filter(is_ref, Comparison == "Breast"),
      aes(x = 1, y = Haplotype, label = "ref."),
      colour = "black", size = 2.4, vjust = -1.2, hjust = 0.5,
      inherit.aes = FALSE
    ) +

    geom_text(
      data = plot_data %>% filter(!is_ref, sig_label != ""),
      aes(label = sig_label), vjust = -0.75, hjust = 0.5, size = 3.0,
      show.legend = FALSE
    ) +

    scale_x_log10(
      limits = c(x_min, x_max), breaks = use_breaks,
      labels = ifelse(use_breaks == as.integer(use_breaks),
                      as.character(as.integer(use_breaks)),
                      as.character(use_breaks))
    ) +
    scale_colour_manual(values = COHORT_COLS, name = "Comparison") +
    scale_shape_manual(values = c("Breast" = 16, "Endometrium" = 17, "All" = 15),
                       name = "Comparison") +

    facet_wrap(~ Comparison, ncol = length(unique(all_stats$Comparison))) +
    labs(
      title    = "Haplotype Association with Tumour (vs Healthy Reference)",
      subtitle = sprintf(
        "Reference: %s (%s global freq.)  |  ◆ = reference  |  Stars = BH-FDR: * <0.05  ** <0.01  *** <0.001  |  ► CI truncated",
        ref_haplo, ref_freq),
      x = "Odds Ratio (log scale)", y = NULL
    ) +
    theme(
      legend.position    = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      axis.text.x        = element_text(size = 9),
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      strip.text         = element_text(size = 11, face = "bold")
    )

  freq_strip_data <- freq_df %>%
    mutate(
      Haplotype = factor(as.character(Haplotype), levels = rev(all_hap_levels)),
      is_ref    = as.character(Haplotype) == ref_haplo_global,
      freq_pct  = Frequency * 100
    )

  freq_strip <- ggplot(freq_strip_data,
    aes(x = freq_pct, y = Haplotype, colour = is_ref)) +
    geom_segment(aes(x = 0, xend = freq_pct, yend = Haplotype), linewidth = 0.6) +
    geom_point(size = 2.2) +
    geom_text(aes(label = sprintf("%.1f%%", freq_pct)),
              hjust = -0.3, size = 2.4, colour = "grey25") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.6)),
                       breaks = scales::pretty_breaks(3)) +
    scale_colour_manual(values = c("FALSE" = "grey55", "TRUE" = "#e74c3c"),
                        guide = "none") +
    labs(title = "Global\nfreq.", subtitle = " ",
         x = "Frequency (%)", y = "Haplotype") +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      axis.text.y        = element_text(size = 9, face = "bold"),
      axis.title.x       = element_text(size = 8.5, face = "bold"),
      axis.title.y       = element_text(size = 9,   face = "bold"),
      plot.title         = element_text(size = 9,   face = "bold", hjust = 0.5),
      plot.subtitle      = element_text(size = 7,   colour = "white")
    )

  fig4 <- freq_strip + forest_main + plot_layout(widths = c(0.45, 3))

  SAVE(file.path(OUT_DIR, "19_Fig4_Forest_Plot.png"), fig4, width = 17, height = 7)
}

# =============================================================================
# 11. FIGURE 5 — FDR Heatmap across comparisons
#     IMPROVED vs v2/v3:
#       • Uses BH-FDR (not raw Score_P) for both colour and text labels
#       • Cells with FDR < 0.05 show the FDR value AND a significance symbol
#       • Non-significant cells show "ns" in small grey text (not blank)
#       • Colour scale anchored to actual data range, not a fixed ceiling
#       • Added annotation row showing haplotype global frequency
# =============================================================================
cat("STEP 11: Figure 5 — FDR Heatmap...\n")

if (nrow(all_stats) == 0) {
  cat("  SKIPPING Fig5: no GLM results available.\n")
} else {
  fdr_data <- all_stats %>%
    select(Haplotype, Comparison, fdr, sig_label) %>%
    mutate(
      neglog10_fdr = -log10(pmax(fdr, 1e-10)),
      cell_label   = case_when(
        !is.na(fdr) & fdr < 0.001 ~ paste0(formatC(fdr, digits = 1, format = "e"), "\n***"),
        !is.na(fdr) & fdr < 0.01  ~ paste0(formatC(fdr, digits = 1, format = "e"), "\n**"),
        !is.na(fdr) & fdr < 0.05  ~ paste0(formatC(fdr, digits = 1, format = "e"), "\n*"),
        !is.na(fdr)                ~ "ns",
        TRUE                       ~ "NA"
      ),
      text_colour = ifelse(!is.na(fdr) & fdr < 0.05, "white", "grey60")
    )

  max_val <- max(fdr_data$neglog10_fdr[is.finite(fdr_data$neglog10_fdr)], na.rm = TRUE)
  if (is.infinite(max_val) || is.na(max_val)) max_val <- 3

  # Add global frequency annotation strip at the bottom
  freq_anno <- freq_df %>%
    mutate(
      Comparison   = "Global\nFreq.",
      neglog10_fdr = NA_real_,
      cell_label   = percent(Frequency, accuracy = 0.1),
      text_colour  = "grey20",
      sig_label    = ""
    ) %>%
    select(Haplotype, Comparison, neglog10_fdr, cell_label, text_colour, sig_label)

  plot_data5 <- bind_rows(fdr_data, freq_anno) %>%
    mutate(
      Comparison = factor(Comparison, levels = c("Breast", "Endometrium", "All", "Global\nFreq."))
    )

  fig5 <- ggplot(plot_data5, aes(x = Comparison, y = Haplotype, fill = neglog10_fdr)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = cell_label, colour = text_colour),
              size = 2.4, lineheight = 0.9) +
    scale_colour_identity() +
    scale_fill_gradientn(
      colors   = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
      values   = rescale(c(0, max_val * 0.2, max_val * 0.5,
                           max_val * 0.75, max_val)),
      name     = expression(-log[10](FDR)),
      na.value = "grey93",
      limits   = c(0, max_val),
      guide    = guide_colorbar(barheight = 8, barwidth = 1.2)
    ) +
    # Vertical separator before the frequency annotation column
    geom_vline(xintercept = 3.5, colour = "grey60", linewidth = 0.8, linetype = "dashed") +
    labs(
      title    = "Haplotype Association: BH-FDR Significance Heatmap",
      subtitle = "Colour = −log₁₀(FDR)  |  Text: FDR value + * symbol where FDR < 0.05; 'ns' = not significant  |  Right column = global haplotype frequency",
      x        = "Comparison (Tumour vs Healthy)", y = "Haplotype"
    ) +
    theme(
      panel.grid  = element_blank(),
      axis.line   = element_blank(),
      axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5)
    )

  SAVE(file.path(OUT_DIR, "19_Fig5_FDR_Heatmap.png"), fig5, width = 8, height = 7)
}

# =============================================================================
# 11b. FIGURE 6 — Raw P-value Heatmap
#      Companion to Fig 5. Uses uncorrected p-values (p_glm from Firth logistic
#      regression) so the reader can see nominal signals that are obscured by
#      FDR correction in a small-sample exploratory analysis.
#      Threshold for labelling: p < 0.05 (nominal, uncorrected).
#      Cells labelled with raw p-value + * where p < 0.05; 'ns' otherwise.
# =============================================================================
cat("STEP 11b: Figure 6 — Raw P-value Heatmap...\n")

if (nrow(all_stats) == 0) {
  cat("  SKIPPING Fig6: no GLM results available.\n")
} else {
  praw_data <- all_stats %>%
    select(Haplotype, Comparison, p_glm) %>%
    mutate(
      neglog10_p   = -log10(pmax(p_glm, 1e-10)),
      cell_label   = case_when(
        !is.na(p_glm) & p_glm < 0.001 ~ paste0(formatC(p_glm, digits = 1, format = "e"), "\n***"),
        !is.na(p_glm) & p_glm < 0.01  ~ paste0(formatC(p_glm, digits = 1, format = "e"), "\n**"),
        !is.na(p_glm) & p_glm < 0.05  ~ paste0(formatC(p_glm, digits = 1, format = "e"), "\n*"),
        !is.na(p_glm)                  ~ "ns",
        TRUE                           ~ "NA"
      ),
      text_colour = ifelse(!is.na(p_glm) & p_glm < 0.05, "white", "grey60")
    )

  max_val_p <- max(praw_data$neglog10_p[is.finite(praw_data$neglog10_p)], na.rm = TRUE)
  if (is.infinite(max_val_p) || is.na(max_val_p)) max_val_p <- 3

  # Add global frequency annotation strip
  freq_anno_p <- freq_df %>%
    mutate(
      Comparison   = "Global\nFreq.",
      neglog10_p   = NA_real_,
      cell_label   = percent(Frequency, accuracy = 0.1),
      text_colour  = "grey20"
    ) %>%
    select(Haplotype, Comparison, neglog10_p, cell_label, text_colour)

  plot_data6 <- bind_rows(praw_data, freq_anno_p) %>%
    mutate(
      Comparison = factor(Comparison,
                          levels = c("Breast", "Endometrium", "All", "Global\nFreq."))
    )

  fig6 <- ggplot(plot_data6,
                 aes(x = Comparison, y = Haplotype, fill = neglog10_p)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = cell_label, colour = text_colour),
              size = 2.4, lineheight = 0.9) +
    scale_colour_identity() +
    scale_fill_gradientn(
      colors   = c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"),
      values   = rescale(c(0, max_val_p * 0.2, max_val_p * 0.5,
                           max_val_p * 0.75, max_val_p)),
      name     = expression(-log[10](p)),
      na.value = "grey93",
      limits   = c(0, max_val_p),
      guide    = guide_colorbar(barheight = 8, barwidth = 1.2)
    ) +
    geom_vline(xintercept = 3.5, colour = "grey60", linewidth = 0.8,
               linetype = "dashed") +
    labs(
      title    = "Haplotype Association: Nominal P-value Heatmap (Uncorrected)",
      subtitle = paste0(
        "Colour = −log₁₀(raw p)  |  Firth logistic regression p-values, NO multiple-testing correction  |",
        "  * p<0.05  ** p<0.01  *** p<0.001  |  'ns' = not significant  |  Right column = global freq."
      ),
      x = "Comparison (Tumour vs Healthy)", y = "Haplotype"
    ) +
    theme(
      panel.grid  = element_blank(),
      axis.line   = element_blank(),
      axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5)
    )

  SAVE(file.path(OUT_DIR, "19_Fig6_RawP_Heatmap.png"), fig6, width = 8, height = 7)
}

# =============================================================================
# 12. EXPORT RESULTS TO EXCEL
# =============================================================================
cat("STEP 12: Exporting results...\n")

write.xlsx(
  list(
    Global_Frequencies = freq_df %>%
      select(Haplotype, Frequency, Allele_String, all_of(snp_labels)),
    Stratified_Freqs   = strat_freqs,
    Association_Stats  = if (nrow(all_stats) > 0)
      all_stats %>%
        select(any_of(c("Haplotype", "Comparison", "OR", "Lower", "Upper",
                        "Beta", "SE", "p_glm", "Score_P", "fdr", "sig_label")))
      else data.frame(Note = "No GLM results")
  ),
  file.path(OUT_DIR, "19_Haplotype_Results_v5.xlsx"),
  overwrite = TRUE
)

# =============================================================================
# 13. SUMMARY TO CONSOLE
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("✓ ANALYSIS COMPLETE — v5\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  SNPs in analysis     : %d\n", n_snps))
cat(sprintf("  Common samples       : %d\n", length(common_samples)))
cat(sprintf("  Haplotypes (≥1%%)    : %d\n", length(common_haplotypes)))
cat(sprintf("  Statistical reference  : %s\n", ref_haplo_global))
cat(sprintf("  All-ref proxy (report): %s\n", allref_haplo))
cat("\n  Global haplotype freqs (top 5):\n")
print(head(freq_df %>% select(Haplotype, Frequency, Allele_String), 5))
cat("\n  Figures saved to:", OUT_DIR, "\n")
cat("    19_Fig0_LD_Heatmap.png\n")
cat("    19_Fig1a_SNP_Carrier_Counts.png\n")
cat("    19_Fig1b_Haplotype_Composition.png\n")
cat("    19_Fig1_Haplotype_Composition.png   (combined)\n")
cat("    19_Fig2_Global_Frequencies.png\n")
cat("    19_Fig3_Stratified_Frequencies.png\n")
cat("    19_Fig4_Forest_Plot.png\n")
cat("    19_Fig5_FDR_Heatmap.png\n")
cat("    19_Fig6_RawP_Heatmap.png\n")
cat("    19_Haplotype_Results_v5.xlsx\n\n")
