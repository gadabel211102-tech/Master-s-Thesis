#!/usr/bin/env Rscript
# =============================================================================
# Script 15: GSDMB Haplotype Analysis — v2 (Publication-Ready)
# Comparisons:
#   (1) Breast Healthy vs Breast Tumour
#   (2) Endometrium Healthy vs Endometrium Tumour
#   (3) All Healthy vs All Tumour (pooled)
#   (4) Global haplotype frequencies
# Haplotypes filtered to ≥1% global frequency.
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
  if (!requireNamespace("patchwork", quietly = TRUE))
    install.packages("patchwork", repos = "https://cloud.r-project.org")
  library(patchwork)
  if (!requireNamespace("logistf", quietly = TRUE))
    install.packages("logistf", repos = "https://cloud.r-project.org")
  library(logistf)
})

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
COHORT_COLS <- c("Breast" = "#4393C3", "Endometrium" = "#D6604D", "All" = "#878787")

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

# Confirm expected columns exist in metadata
stopifnot(all(c("Sample", "Cohort", "Tissue") %in% colnames(metadata)))
# Expected Tissue values: "Healthy" and "Tumour"
# Expected Cohort values: "Breast"  and "Endometrium"

geno_filtered <- geno_raw %>%
  mutate(extracted_pos = as.integer(str_split_fixed(ID, ":", 4)[, 2])) %>%
  inner_join(annotation_map, by = c("extracted_pos" = "POS_int"))

n_snps      <- nrow(geno_filtered)
snp_labels  <- geno_filtered$rsID_clean

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
# Store original haplo.em internal index (1-based row number) for dosage mapping
freq_df$em_index   <- seq_len(nrow(freq_df))
freq_df$Frequency  <- as.numeric(global_em$hap.prob)

# Apply ≥1% filter and sort by descending frequency BEFORE assigning labels
# → H1 = most frequent, H2 = second most frequent, etc.
freq_df <- freq_df %>%
  filter(Frequency >= 0.01) %>%
  arrange(desc(Frequency)) %>%
  mutate(
    Haplotype = paste0("H", seq_len(n())),   # H1 = most frequent
    Haplotype = factor(Haplotype, levels = Haplotype),
    Allele_String = apply(
      across(all_of(snp_labels)), 1,
      function(x) paste(ifelse(x == 1, "Ref", "Alt"), collapse = "-")
    )
  )

# Build lookup: em_index -> new Haplotype label (for dosage matrix mapping)
em_to_haplabel <- setNames(as.character(freq_df$Haplotype), freq_df$em_index)
# em_to_haplabel[as.character(em_index)] gives the new label, NA if filtered out

common_haplotypes <- as.character(freq_df$Haplotype)
cat(sprintf("  Haplotypes passing ≥1%% filter: %d\n", length(common_haplotypes)))

# =============================================================================
# 5. STATISTICAL MODELS  — Tumour vs Healthy within each comparison group
# =============================================================================
cat("STEP 5: Running association models...\n")

# Define the three comparisons
comparisons <- list(
  list(label = "Breast",      cohorts = "Breast",                name = "Breast"),
  list(label = "Endometrium", cohorts = "Endometrium",           name = "Endometrium"),
  list(label = "All",         cohorts = c("Breast","Endometrium"), name = "All")
)

run_comparison <- function(comp) {
  sub_meta <- metadata %>%
    filter(Cohort %in% comp$cohorts, Sample %in% common_samples,
           Tissue %in% c("Healthy", "Tumour"))

  n_normal <- sum(sub_meta$Tissue == "Healthy")
  n_tumour <- sum(sub_meta$Tissue == "Tumour")

  cat(sprintf(
    "  %s: n=%d  (Healthy=%d, Tumour=%d)\n",
    comp$label, nrow(sub_meta), n_normal, n_tumour
  ))

  if (nrow(sub_meta) < 10) {
    warning(sprintf("  SKIPPING %s: too few total samples (%d).", comp$label, nrow(sub_meta)))
    return(NULL)
  }
  if (n_normal == 0 || n_tumour == 0) {
    warning(sprintf(
      "  SKIPPING %s: both Healthy and Tumour samples required (Healthy=%d, Tumour=%d).\n  This cohort may be tumour-only — GLM and score test cannot run.",
      comp$label, n_normal, n_tumour
    ))
    return(NULL)
  }

  sub_samples <- sub_meta$Sample
  sub_mat     <- mat_all[sub_samples, , drop = FALSE]
  y           <- as.numeric(sub_meta$Tissue == "Tumour")   # 0 = Healthy, 1 = Tumour

  geno_sub <- setupGeno(sub_mat, locus.label = snp_labels)

  # ── Score test (haplo.stats) ──────────────────────────────────────────────
  sc <- tryCatch(
    haplo.score(y, geno_sub, trait.type = "binomial",
                min.count = ceiling(0.01 * length(y))),
    error = function(e) { message("  Score test failed: ", e$message); NULL }
  )

  # ── Build dosage matrix using GLOBAL EM haplotype definitions ─────────────
  # This ensures haplotype labels are consistent with common_haplotypes.
  # We run haplo.em on this stratum's geno but map results to global haplotype
  # allele strings so indices are comparable.
  #
  # Strategy: use global_em haplotype allele table.
  # For each subject, their most likely haplotype pair (from global EM on full
  # dataset) is already captured in global_em$row/hap1/hap2/post.
  # Subset to subjects in this comparison stratum.

  n_sub      <- nrow(sub_mat)
  n_haps_g     <- length(global_em$hap.prob)
  # Map each haplo.em internal index to our frequency-ranked label (NA if <1%)
  hap_labels_g <- em_to_haplabel[as.character(seq_len(n_haps_g))]

  # Map: global sample index -> local row in dosage_mat (NA if not in stratum)
  # global_em$row is 1-based index into common_samples
  # We need: for each entry in global_em, what row is it in sub_samples?
  # sub_samples are a subset of common_samples, so:
  #   global_sample_name = common_samples[global_em$row[k]]
  #   local_i = which(sub_samples == global_sample_name)
  global_to_local <- match(common_samples, sub_samples)
  # global_to_local[i] = row in sub_samples for common_samples[i], NA if not present

  dosage_mat <- matrix(0, nrow = n_sub, ncol = length(common_haplotypes),
                       dimnames = list(sub_samples, common_haplotypes))

  em_row_g  <- global_em$indx.subj   # 1-based subject index into common_samples
  em_hap1_g <- global_em$hap1code    # first haplotype index in pair
  em_hap2_g <- global_em$hap2code    # second haplotype index in pair
  em_post_g <- global_em$post        # posterior probability for this pair

  for (k in seq_along(em_row_g)) {
    global_i <- em_row_g[k]               # 1-based index into common_samples
    local_i  <- global_to_local[global_i] # row in sub_samples (NA if not in stratum)
    if (is.na(local_i)) next
    h1   <- em_hap1_g[k]
    h2   <- em_hap2_g[k]
    post <- em_post_g[k]
    lbl1 <- if (!is.na(h1) && h1 >= 1 && h1 <= n_haps_g) hap_labels_g[h1] else NA
    lbl2 <- if (!is.na(h2) && h2 >= 1 && h2 <= n_haps_g) hap_labels_g[h2] else NA
    if (!is.na(lbl1) && lbl1 %in% common_haplotypes)
      dosage_mat[local_i, lbl1] <- dosage_mat[local_i, lbl1] + post
    if (!is.na(lbl2) && lbl2 %in% common_haplotypes)
      dosage_mat[local_i, lbl2] <- dosage_mat[local_i, lbl2] + post
  }

  # ── Firth penalised logistic regression (handles complete separation) ─────
  # Reference = globally most frequent haplotype (H1 after freq_df ordering)
  # Reference = H1 (most frequent globally); test all other common haplotypes
  ref_haplo_label <- as.character(freq_df$Haplotype[1])
  test_haplos     <- setdiff(common_haplotypes, ref_haplo_label)


  if (length(test_haplos) == 0) return(NULL)

  results <- lapply(test_haplos, function(hname) {
    df_fit <- data.frame(y = y, hap = dosage_mat[, hname])
    # Skip if haplotype dosage has zero variance in this stratum
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
      Haplotype  = hname,
      Beta       = beta,
      SE         = se,
      p_glm      = p,
      OR         = exp(beta),
      Lower      = exp(beta - 1.96 * se),
      Upper      = exp(beta + 1.96 * se),
      Comparison = comp$label,
      row.names  = NULL
    )
  })

  res <- bind_rows(results)
  if (nrow(res) == 0) return(NULL)

  # ── Attach score test p-values ────────────────────────────────────────────
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
if (nrow(all_stats_raw) > 0) {
  all_stats <- all_stats_raw %>%
    mutate(
      Comparison = factor(Comparison, levels = c("Breast","Endometrium","All")),
      Haplotype  = factor(Haplotype,  levels = common_haplotypes),
      sig_label  = case_when(
        !is.na(Score_P) & Score_P < 0.001 ~ "***",
        !is.na(Score_P) & Score_P < 0.01  ~ "**",
        !is.na(Score_P) & Score_P < 0.05  ~ "*",
        TRUE                               ~ ""
      )
    )
  cat(sprintf("  Comparisons with results: %s\n",
              paste(unique(as.character(all_stats$Comparison)), collapse = ", ")))
  print(head(all_stats %>% select(Haplotype, Comparison, OR, Score_P), 6))
} else {
  all_stats <- all_stats_raw
  cat("  WARNING: No GLM results — check the 50+ warnings above.\n")
  cat("  Running warnings() to diagnose:\n")
  print(warnings())
}

# =============================================================================
# 6. FIGURE 1 — Haplotype Composition Tile Plot  (IMPROVED)
#    • Fig 1a: SNP Alt-Allele Carrier Counts per variant (standalone bar chart)
#    • Fig 1b: Haplotype × SNP allele tile grid with global frequency annotation
# =============================================================================
cat("STEP 6: Figure 1 — Haplotype Composition...\n")

tile_data <- freq_df %>%
  select(Haplotype, Frequency, all_of(snp_labels)) %>%
  pivot_longer(cols = all_of(snp_labels), names_to = "rsID_clean", values_to = "Val") %>%
  left_join(annotation_map, by = "rsID_clean") %>%
  mutate(
    Allele     = factor(as.numeric(Val) - 1, levels = c(0, 1), labels = c("Ref", "Alt")),
    Haplotype  = factor(Haplotype, levels = rev(common_haplotypes)),
    Freq_Label = percent(Frequency, accuracy = 0.1)
  )

# ── Fig 1a: SNP carrier count bar (now a standalone figure with clear labelling) ──
# "Carrier" = any sample carrying ≥1 alt allele at that SNP
fig1a <- annotation_map %>%
  mutate(
    Short_Label = paste0(rsID_clean, "\n(", Gene, ")"),
    Short_Label = factor(Short_Label, levels = annotation_map %>%
                           mutate(Short_Label = paste0(rsID_clean, "\n(", Gene, ")")) %>%
                           pull(Short_Label))
  ) %>%
  ggplot(aes(x = Short_Label, y = Carrier_Count)) +
  geom_col(fill = "#2166AC", width = 0.65) +
  geom_text(aes(label = Carrier_Count), vjust = -0.4, size = 3, fontface = "bold",
            color = "#2166AC") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     breaks = scales::pretty_breaks(5)) +
  labs(
    title    = "Alt-Allele Carrier Counts per GSDMB/GSDMA SNP",
    subtitle = "Number of samples carrying ≥1 alternate allele at each variant",
    x        = "Variant (rsID, Gene)",
    y        = "Carriers (n)"
  ) +
  theme(
    axis.text.x        = element_text(angle = 40, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3)
  )

SAVE(file.path(OUT_DIR, "19_Fig1a_SNP_Carrier_Counts.png"), fig1a, width = 12, height = 5)

# ── Fig 1b: Haplotype allele composition tile (cleaner, no embedded bar) ──
# Add Alt-allele count column to give per-SNP context via x-axis subtitle
snp_alt_counts <- annotation_map %>%
  mutate(x_label = paste0(rsID_clean, "\n[", Gene, ", n=", Carrier_Count, "]")) %>%
  pull(x_label, name = Detailed_Label)  # maps old Detailed_Label → compact label

tile_data2 <- tile_data %>%
  mutate(
    Compact_Label = paste0(rsID_clean, "\n[", Gene, ", n=", Carrier_Count, "]"),
    Compact_Label = factor(Compact_Label,
      levels = annotation_map %>%
        mutate(l = paste0(rsID_clean, "\n[", Gene, ", n=", Carrier_Count, "]")) %>%
        pull(l))
  )

ordered_compact <- levels(tile_data2$Compact_Label)

fig1b <- tile_data2 %>%
  ggplot(aes(x = Compact_Label, y = Haplotype, fill = Allele)) +
  geom_tile(color = "white", linewidth = 0.45) +
  # Global frequency label on right
  geom_text(
    data = tile_data2 %>% distinct(Haplotype, Frequency, Freq_Label),
    aes(x = length(ordered_compact) + 0.75, y = Haplotype,
        label = Freq_Label, fill = NULL),
    size = 2.7, hjust = 0, color = "grey25", fontface = "bold"
  ) +
  # Haplotype row label (left, inside)
  scale_fill_manual(
    values = c("Ref" = "#E8E8E8", "Alt" = "#C0392B"),
    name   = "Allele"
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 2.0))) +
  labs(
    title    = "GSDMB Haplotype Allele Composition",
    subtitle = sprintf(
      "Each row = one haplotype (H1 most frequent); red = alternate allele | n=%d haplotypes ≥1%% global freq.",
      length(common_haplotypes)
    ),
    x = "Variant (rsID, Gene, n carriers)",
    y = "Haplotype"
  ) +
  theme(
    axis.text.x    = element_text(angle = 40, hjust = 1, size = 7.5),
    panel.grid     = element_blank(),
    legend.position = "right",
    plot.margin    = margin(6, 12, 6, 8)
  )

SAVE(file.path(OUT_DIR, "19_Fig1b_Haplotype_Composition.png"), fig1b, width = 14, height = 8)

# Keep a combined version for convenience (carrier bar on top, tile below)
fig1_combined <- fig1a / fig1b + plot_layout(heights = c(1, 2.8))
SAVE(file.path(OUT_DIR, "19_Fig1_Haplotype_Composition.png"), fig1_combined, width = 14, height = 14)

# =============================================================================
# 7. FIGURE 2 — Haplotype Frequencies: Global + Per-Cohort  (IMPROVED)
#    Left panel : global lollipop (H1–H15, descending frequency)
#    Right panels: per-cohort × tissue faceted dot plot (Breast / Endometrium)
#    Both panels share the same y-axis haplotype order.
# =============================================================================
cat("STEP 7: Figure 2 — Global + Per-Cohort Haplotype Frequencies...\n")

# ── Haplotype level order: H1 at top (most frequent) ─────────────────────────
hap_order <- rev(common_haplotypes)   # rev so H1 ends up at top in coord_flip

# ── Panel A: Global lollipop ─────────────────────────────────────────────────
global_plot_df <- freq_df %>%
  mutate(Haplotype = factor(Haplotype, levels = hap_order))

p_global <- ggplot(global_plot_df, aes(x = Haplotype, y = Frequency)) +
  geom_segment(aes(xend = Haplotype, y = 0, yend = Frequency),
               color = "grey55", linewidth = 0.8) +
  geom_point(aes(size = Frequency), color = "#2166AC", alpha = 0.88) +
  geom_text(aes(label = percent(Frequency, accuracy = 0.1)),
            hjust = -0.3, size = 2.6, color = "grey25") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.25)),
                     breaks = scales::pretty_breaks(4)) +
  scale_size_continuous(range = c(1.5, 8), guide = "none") +
  coord_flip() +
  labs(title = "Global", subtitle = "All samples\ncombined",
       x = NULL, y = "Frequency") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    plot.title    = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey45"),
    axis.text.y   = element_text(size = 8, face = "bold")
  )

# ── Panel B: Per-cohort × tissue faceted dot plot ────────────────────────────
# Cohort-level frequencies (Breast Healthy, Breast Tumour,
#                           Endometrium Healthy, Endometrium Tumour)
# strat_freqs is built in STEP 8 — so we need to run panel B AFTER strat_freqs.
# We defer the full figure assembly to after Step 8.
# (Placeholder — actual fig2 is assembled below step 8.)

# Store global panel for later assembly
.fig2_global_panel <- p_global

# =============================================================================
# 8. FIGURE 3 — Haplotype Frequencies Stratified by Group
#    (estimated from GLM intercept + betas, compared Healthy vs Tumour)
# =============================================================================
cat("STEP 8: Figure 3 — Stratified Frequencies by Group...\n")

# Compute haplotype frequencies per Tissue × Cohort using EM on each stratum
compute_stratum_freqs <- function(tissue_val, cohort_val, label) {
  sub_meta <- metadata %>%
    filter(
      Sample %in% common_samples,
      Tissue == tissue_val,
      if (!is.null(cohort_val)) Cohort %in% cohort_val else TRUE
    )
  if (nrow(sub_meta) < 5) return(NULL)
  sub_mat <- mat_all[sub_meta$Sample, , drop = FALSE]
  em      <- tryCatch(
    haplo.em(setupGeno(sub_mat, locus.label = snp_labels)),
    error = function(e) NULL
  )
  if (is.null(em)) return(NULL)
  hap_df <- as.data.frame(em$haplotype)
  colnames(hap_df) <- snp_labels
  # Map em_index to frequency-ranked global label
  hap_df$em_index  <- seq_len(nrow(hap_df))
  # Match allele patterns to global freq_df to get consistent labels
  global_patterns  <- apply(freq_df[, snp_labels], 1, paste, collapse="-")
  stratum_patterns <- apply(hap_df[, snp_labels], 1, paste, collapse="-")
  hap_df$Haplotype <- as.character(freq_df$Haplotype)[match(stratum_patterns, global_patterns)]
  hap_df$Frequency  <- as.numeric(em$hap.prob)
  hap_df$Group      <- label
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
    )
  )

fig3 <- ggplot(strat_freqs,
               aes(x = Haplotype, y = Frequency, fill = Group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = GROUP_COLS, name = "Group",
                    labels = function(x) gsub("_", " ", x)) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~Cohort, ncol = 1, scales = "free_y") +
  labs(
    title    = "Haplotype Frequencies: Healthy vs Tumour",
    subtitle = "Stratified by cohort; EM-estimated frequencies per stratum",
    x        = "Haplotype", y = "Frequency"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

SAVE(file.path(OUT_DIR, "19_Fig3_Stratified_Frequencies.png"), fig3, width = 11, height = 12)

# =============================================================================
# 8b. FIGURE 2 — ASSEMBLY (Global + Per-Cohort lollipop panels)
#     Now that strat_freqs is ready, build the cohort panels and combine.
# =============================================================================
cat("STEP 8b: Figure 2 — Assembling global + cohort frequency plot...\n")

cohort_plot_df <- strat_freqs %>%
  filter(Cohort != "All (Pooled)") %>%
  mutate(
    Haplotype = factor(Haplotype, levels = hap_order),
    Tissue    = factor(Tissue, levels = c("Healthy", "Tumour")),
    Cohort    = factor(Cohort, levels = c("Breast", "Endometrium"))
  )

TISSUE_COLS <- c("Healthy" = "#4393C3", "Tumour" = "#D6604D")

make_cohort_panel <- function(cohort_name) {
  df <- cohort_plot_df %>% filter(Cohort == cohort_name)
  df_h <- df %>% filter(Tissue == "Healthy")
  df_t <- df %>% filter(Tissue == "Tumour")
  ggplot(df, aes(x = Haplotype, y = Frequency)) +
    # Healthy stems + points (nudged slightly left)
    geom_segment(data = df_h,
      aes(xend = Haplotype, y = 0, yend = Frequency),
      color = TISSUE_COLS["Healthy"], linewidth = 0.6, alpha = 0.55,
      position = position_nudge(x = -0.18)
    ) +
    geom_point(data = df_h, aes(shape = Tissue),
      color = TISSUE_COLS["Healthy"],
      position = position_nudge(x = -0.18),
      size = 2.6, alpha = 0.88
    ) +
    # Tumour stems + points (nudged slightly right)
    geom_segment(data = df_t,
      aes(xend = Haplotype, y = 0, yend = Frequency),
      color = TISSUE_COLS["Tumour"], linewidth = 0.6, alpha = 0.55,
      position = position_nudge(x = 0.18)
    ) +
    geom_point(data = df_t, aes(shape = Tissue),
      color = TISSUE_COLS["Tumour"],
      position = position_nudge(x = 0.18),
      size = 2.6, alpha = 0.88
    ) +
    scale_shape_manual(values = c("Healthy" = 16, "Tumour" = 17), name = "Tissue") +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.18)),
                       breaks = scales::pretty_breaks(4)) +
    coord_flip() +
    labs(title = cohort_name, subtitle = "EM freq. per\nstratum",
         x = NULL, y = "Frequency") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      axis.text.y        = element_blank(),   # shared with left panel
      axis.ticks.y       = element_blank(),
      plot.title         = element_text(size = 10, face = "bold"),
      plot.subtitle      = element_text(size = 8, color = "grey45"),
      legend.position    = "right"
    )
}

p_breast  <- make_cohort_panel("Breast")
p_endo    <- make_cohort_panel("Endometrium")

fig2 <- (.fig2_global_panel | p_breast | p_endo) +
  plot_layout(widths = c(1.4, 1, 1), guides = "collect") +
  plot_annotation(
    title    = "Haplotype Frequencies — Global and Per-Cohort (GSDMB)",
    subtitle = "Haplotypes ≥1% global frequency shown; cohort panels show Healthy vs Tumour EM estimates",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  )

SAVE(file.path(OUT_DIR, "19_Fig2_Global_Frequencies.png"), fig2, width = 13, height = 6.5)

# =============================================================================
# 9. FIGURE 4 — Forest Plot (Odds Ratios, Tumour vs Healthy)
# =============================================================================
cat("STEP 9: Figure 4 — Forest Plot...\n")

# Reference haplotype label (most frequent)
ref_haplo <- as.character(freq_df$Haplotype[1])
ref_freq  <- percent(freq_df$Frequency[1], accuracy = 0.1)

if (nrow(all_stats) == 0) {
  cat("  SKIPPING Fig4: no GLM results available.\n")
} else {
  # Compute sensible x-axis limits from actual data (exclude infinite/extreme)
  or_vals <- all_stats$OR[is.finite(all_stats$OR) & all_stats$OR > 0]
  lo_vals <- all_stats$Lower[is.finite(all_stats$Lower) & all_stats$Lower > 0]
  hi_vals <- all_stats$Upper[is.finite(all_stats$Upper) & all_stats$Upper > 0]
  x_min   <- max(0.05, min(lo_vals, na.rm = TRUE) * 0.7)
  x_max   <- min(50,   max(hi_vals, na.rm = TRUE) * 1.4)

  # Nice log10 breaks within the actual data range
  all_breaks <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32)
  use_breaks <- all_breaks[all_breaks >= x_min & all_breaks <= x_max]
  if (1 %in% use_breaks == FALSE) use_breaks <- sort(c(1, use_breaks))

  fig4 <- ggplot(all_stats,
                 aes(x = OR, y = Haplotype, color = Comparison, shape = Comparison)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_errorbar(aes(xmin = pmax(Lower, x_min * 0.9),
                      xmax = pmin(Upper, x_max * 1.1)),
                  width = 0.3, linewidth = 0.6, alpha = 0.85,
                  orientation = "y") +
    geom_point(size = 2.5, stroke = 0.5) +
    geom_text(aes(label = sig_label, x = pmin(Upper, x_max) * 1.05),
              hjust = -0.2, size = 3, show.legend = FALSE) +
    scale_x_log10(
      limits = c(x_min, x_max),
      breaks = use_breaks,
      labels = ifelse(use_breaks == as.integer(use_breaks),
                      as.character(as.integer(use_breaks)),
                      as.character(use_breaks))
    ) +
    scale_color_manual(values = COHORT_COLS, name = "Comparison") +
    scale_shape_manual(values = c("Breast" = 16, "Endometrium" = 17, "All" = 15),
                       name = "Comparison") +
    facet_wrap(~Comparison, ncol = length(unique(all_stats$Comparison))) +
    labs(
      title    = "Haplotype Association with Tumour (vs Healthy Reference)",
      subtitle = sprintf(
        "Reference haplotype: %s (%s global freq.) | * p<0.05  ** p<0.01  *** p<0.001",
        ref_haplo, ref_freq),
      x = "Odds Ratio (log scale)", y = "Haplotype"
    ) +
    theme(
      legend.position  = "none",
      panel.grid.major.y = element_blank(),
      axis.text.x      = element_text(size = 9, angle = 0, hjust = 0.5),
      strip.text       = element_text(size = 11, face = "bold")
    )
  SAVE(file.path(OUT_DIR, "19_Fig4_Forest_Plot.png"), fig4, width = 14, height = 7)
}

# =============================================================================
# 10. FIGURE 5 — Score Test  –log10(p) Heatmap across comparisons
# =============================================================================
cat("STEP 10: Figure 5 — Score Test Heatmap...\n")

if (nrow(all_stats) == 0) {
  cat("  SKIPPING Fig5: no GLM results available.\n")
} else {
  score_data    <- all_stats %>%
    select(Haplotype, Comparison, Score_P) %>%
    mutate(neglog10p = -log10(pmax(Score_P, 1e-10)))
  sig_threshold <- -log10(0.05)
  max_p         <- max(score_data$neglog10p, na.rm = TRUE)
  if (is.infinite(max_p) || is.na(max_p)) max_p <- sig_threshold * 3

  fig5 <- ggplot(score_data,
                 aes(x = Comparison, y = Haplotype, fill = neglog10p)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(!is.na(Score_P) & Score_P < 0.05,
                                 formatC(Score_P, digits = 2, format = "e"), "")),
              size = 2.5, color = "white", fontface = "bold") +
    scale_fill_gradientn(
      colors   = c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404"),
      # Anchor colour stops to actual data range so differences are visible
      values   = rescale(c(0,
                           max_p * 0.25,
                           max_p * 0.55,
                           max_p * 0.80,
                           max_p)),
      name     = expression(-log[10](p)),
      na.value = "grey85",
      limits   = c(0, max_p),
      guide    = guide_colorbar(barheight = 8, barwidth = 1.2)
    ) +
    labs(
      title    = "Haplotype Score Test Results",
      subtitle = "Text shows p-values for associations with p < 0.05",
      x        = "Comparison (Tumour vs Healthy)", y = "Haplotype"
    ) +
    theme(panel.grid = element_blank(), axis.line = element_blank())

  SAVE(file.path(OUT_DIR, "19_Fig5_ScoreTest_Heatmap.png"), fig5, width = 7, height = 6)
}

# =============================================================================
# 11. EXPORT RESULTS TO EXCEL
# =============================================================================
cat("STEP 11: Exporting results...\n")

write.xlsx(
  list(
    Global_Frequencies   = freq_df %>% select(Haplotype, Frequency, all_of(snp_labels)),
    Stratified_Freqs     = strat_freqs,
    Association_Stats    = if (nrow(all_stats) > 0)
      all_stats %>% select(any_of(c("Haplotype","Comparison","OR","Lower","Upper","Beta","SE","p_glm","Score_P","sig_label")))
      else data.frame(Note = "No GLM results — complete separation in all comparisons")
  ),
  file.path(OUT_DIR, "19_Haplotype_Results_v2.xlsx"),
  overwrite = TRUE
)

# =============================================================================
# 12. SUMMARY TO CONSOLE
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("✓ ANALYSIS COMPLETE\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  SNPs in analysis     : %d\n", n_snps))
cat(sprintf("  Common samples       : %d\n", length(common_samples)))
cat(sprintf("  Haplotypes (≥1%%)    : %d\n", length(common_haplotypes)))
cat(sprintf("  Global haplotype freqs (top 5):\n"))
print(head(freq_df %>% select(Haplotype, Frequency, Allele_String), 5))
cat("\n  Figures saved to   :", OUT_DIR, "\n")
cat("    19_Fig1a_SNP_Carrier_Counts.png\n")
cat("    19_Fig1b_Haplotype_Composition.png\n")
cat("    19_Fig1_Haplotype_Composition.png  (combined)\n")
cat("    19_Fig2_Global_Frequencies.png\n")
cat("    19_Fig3_Stratified_Frequencies.png\n")
cat("    19_Fig4_Forest_Plot.png\n")
cat("    19_Fig5_ScoreTest_Heatmap.png\n")
cat("    19_Haplotype_Results_v2.xlsx\n\n")
