#!/usr/bin/env Rscript

# ==============================================================================
# Script 13B: AMRFinderPlus Gene Heatmap and Clade Enrichment Statistics
#
# Inputs:
#   1. results/13_amrfinderplus/master_AMRFinder_report.tsv
#   2. results/12d_clade_analysis/all_preharvest_704_clades.csv
#   3. results/12d_clade_analysis/clade_summary_counts.csv
#
# Outputs:
#   1. AMR gene presence/absence matrix
#   2. AMR gene prevalence by clade
#   3. Fisher exact test for each AMR gene by clade
#   4. BH-FDR corrected enrichment/depletion statistics
#   5. Heatmap clustered by clade
#   6. Enrichment dotplot
#
# Default analysis level:
#   gene-level using "Element symbol"
#
# To run allele-level instead:
#   AMRFINDER_FEATURE_LEVEL=allele Rscript scripts/13B_amrfinder_clade_heatmap_enrichment.R
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(tibble)
})

message("=== SCRIPT 13B: AMRFINDERPLUS CLADE HEATMAP + FDR ENRICHMENT ===")

# ==============================================================================
# 1. Paths and parameters
# ==============================================================================

base_dir <- "/home/samuelajulo/epi/Ecoli_epi"

amr_file <- file.path(
  base_dir,
  "results",
  "13_amrfinderplus",
  "master_AMRFinder_report.tsv"
)

clade_file <- file.path(
  base_dir,
  "results",
  "12d_clade_analysis",
  "all_preharvest_704_clades.csv"
)

clade_counts_file <- file.path(
  base_dir,
  "results",
  "12d_clade_analysis",
  "clade_summary_counts.csv"
)

metadata_file <- file.path(
  base_dir,
  "results",
  "12c_all_preharvest_combined",
  "target_metadata_704_preharvest_non_dairy.csv"
)

metadata_fallback <- "/home/samuelajulo/epi/data6a.csv"

out_dir <- file.path(
  base_dir,
  "results",
  "13_amrfinderplus",
  "13B_clade_heatmap_enrichment"
)

table_dir <- file.path(out_dir, "tables")
matrix_dir <- file.path(out_dir, "matrices")
figure_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

FEATURE_LEVEL <- Sys.getenv("AMRFINDER_FEATURE_LEVEL", unset = "gene")
FEATURE_LEVEL <- tolower(FEATURE_LEVEL)

if (!FEATURE_LEVEL %in% c("gene", "allele")) {
  stop("[ERROR] AMRFINDER_FEATURE_LEVEL must be 'gene' or 'allele'.", call. = FALSE)
}

MAX_HEATMAP_FEATURES <- 100
MIN_GLOBAL_PREVALENCE_FOR_HEATMAP <- 0.01
FDR_CUTOFF <- 0.05

if (!file.exists(amr_file)) {
  stop("[ERROR] AMRFinderPlus report not found: ", amr_file, call. = FALSE)
}

if (!file.exists(clade_file)) {
  stop("[ERROR] Clade file not found: ", clade_file, call. = FALSE)
}

if (!file.exists(clade_counts_file)) {
  stop("[ERROR] Clade counts file not found: ", clade_counts_file, call. = FALSE)
}

metadata_used <- if (file.exists(metadata_file)) {
  metadata_file
} else {
  metadata_fallback
}

message("[INFO] AMRFinderPlus file: ", amr_file)
message("[INFO] Clade file:         ", clade_file)
message("[INFO] Feature level:      ", FEATURE_LEVEL)
message("[INFO] Output directory:   ", out_dir)

# ==============================================================================
# 2. Helpers
# ==============================================================================

normalize_tip_id <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  x <- str_remove(x, "\\.gz$")
  x <- str_remove(x, "\\.(fasta|fna|fa|gff|gff3|faa|ffn|fa|aln|tsv|csv)$")
  accession_like <- str_extract(x, "GC[AF]_[0-9]+\\.[0-9]+")
  ifelse(is.na(accession_like), x, accession_like)
}

extract_year_num <- function(x) {
  x <- as.character(x)
  out <- str_extract(x, "(18|19|20)[0-9]{2}")
  suppressWarnings(as.integer(out))
}

pick_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}

safe_fisher <- function(tab, alternative = "greater") {
  tryCatch(
    fisher.test(tab, alternative = alternative)$p.value,
    error = function(e) NA_real_
  )
}

theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(linewidth = 0.25, color = "grey88"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "grey40"),
      strip.text = element_text(face = "bold")
    )
}

save_plot <- function(plot, filename, width = 10, height = 6) {
  ggsave(file.path(figure_dir, paste0(filename, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(figure_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 600)
}

# ==============================================================================
# 3. Read clades and metadata
# ==============================================================================

clades <- read_csv(clade_file, show_col_types = FALSE) %>%
  mutate(
    Tip_ID = normalize_tip_id(Tip_ID),
    Clade = as.character(Clade)
  ) %>%
  distinct(Tip_ID, .keep_all = TRUE)

clade_counts <- read_csv(clade_counts_file, show_col_types = FALSE) %>%
  mutate(Clade = as.character(Clade))

if (!"Isolate_Count" %in% names(clade_counts)) {
  clade_counts <- clades %>%
    count(Clade, name = "Isolate_Count")
}

clade_counts <- clade_counts %>%
  arrange(desc(Isolate_Count), Clade)

clade_levels <- clade_counts$Clade

metadata_clean <- tibble(Tip_ID = clades$Tip_ID)

if (file.exists(metadata_used)) {
  metadata <- read_csv(metadata_used, show_col_types = FALSE)

  assembly_col <- pick_first_existing_col(
    metadata,
    c("Assembly", "Tip_ID", "Isolate", "Strain", "Sample", "sample")
  )

  year_col <- pick_first_existing_col(
    metadata,
    c("Collection_year", "Collection year", "collection_year", "Year", "year", "Collection date", "collection_date")
  )

  if (!is.na(assembly_col)) {
    metadata_clean <- metadata %>%
      mutate(
        Tip_ID = normalize_tip_id(.data[[assembly_col]]),
        Collection_year_num = if (!is.na(year_col)) extract_year_num(.data[[year_col]]) else NA_integer_
      ) %>%
      distinct(Tip_ID, .keep_all = TRUE)
  }
}

sample_annotation <- clades %>%
  left_join(clade_counts, by = "Clade") %>%
  left_join(metadata_clean, by = "Tip_ID") %>%
  mutate(
    Clade = factor(Clade, levels = clade_levels)
  ) %>%
  arrange(Clade, Tip_ID)

write_csv(
  sample_annotation,
  file.path(table_dir, "13B_sample_clade_metadata_annotation.csv")
)

# ==============================================================================
# 4. Read AMRFinderPlus report and define features
# ==============================================================================

amr <- read_tsv(amr_file, show_col_types = FALSE)

required_cols <- c("Sample", "Element symbol")
missing_required <- setdiff(required_cols, names(amr))

if (length(missing_required) > 0) {
  stop("[ERROR] Missing required AMRFinder columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

amr_clean <- amr %>%
  mutate(
    Tip_ID = normalize_tip_id(.data[["Sample"]]),
    Element_symbol_clean = str_squish(as.character(.data[["Element symbol"]])),
    Element_name_clean = if ("Element name" %in% names(.)) str_squish(as.character(.data[["Element name"]])) else NA_character_,
    Class_clean = if ("Class" %in% names(.)) str_squish(as.character(.data[["Class"]])) else NA_character_,
    Subclass_clean = if ("Subclass" %in% names(.)) str_squish(as.character(.data[["Subclass"]])) else NA_character_,
    Type_clean = if ("Type" %in% names(.)) str_squish(as.character(.data[["Type"]])) else NA_character_,
    Subtype_clean = if ("Subtype" %in% names(.)) str_squish(as.character(.data[["Subtype"]])) else NA_character_,
    AMR_feature = case_when(
      FEATURE_LEVEL == "allele" & !is.na(Element_name_clean) & Element_name_clean != "" ~
        paste(Element_symbol_clean, Element_name_clean, sep = " | "),
      TRUE ~ Element_symbol_clean
    )
  ) %>%
  filter(
    !is.na(Tip_ID),
    Tip_ID %in% sample_annotation$Tip_ID,
    !is.na(AMR_feature),
    AMR_feature != "",
    AMR_feature != "-"
  )

amr_samples_found <- sort(unique(amr_clean$Tip_ID))
amr_samples_missing <- setdiff(sample_annotation$Tip_ID, amr_samples_found)
amr_samples_extra <- setdiff(unique(normalize_tip_id(amr$Sample)), sample_annotation$Tip_ID)

write_lines(amr_samples_missing, file.path(log_dir, "13B_samples_in_clades_without_amrfinder_hits.txt"))
write_lines(amr_samples_extra, file.path(log_dir, "13B_samples_in_amrfinder_not_in_12D_clades.txt"))

message("[INFO] AMRFinder rows in target clade set: ", nrow(amr_clean))
message("[INFO] Target genomes with >=1 AMRFinder hit: ", length(amr_samples_found))
message("[INFO] Target genomes with no AMRFinder hit: ", length(amr_samples_missing))

feature_info <- amr_clean %>%
  group_by(AMR_feature) %>%
  summarise(
    Element_symbols = paste(sort(unique(Element_symbol_clean[!is.na(Element_symbol_clean)])), collapse = ";"),
    Element_names = paste(sort(unique(Element_name_clean[!is.na(Element_name_clean)])), collapse = ";"),
    Classes = paste(sort(unique(Class_clean[!is.na(Class_clean)])), collapse = ";"),
    Subclasses = paste(sort(unique(Subclass_clean[!is.na(Subclass_clean)])), collapse = ";"),
    Types = paste(sort(unique(Type_clean[!is.na(Type_clean)])), collapse = ";"),
    Subtypes = paste(sort(unique(Subtype_clean[!is.na(Subtype_clean)])), collapse = ";"),
    n_raw_amrfinder_rows = n(),
    n_genomes_with_feature = n_distinct(Tip_ID),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genomes_with_feature), AMR_feature)

write_csv(feature_info, file.path(table_dir, "13B_amrfinder_feature_info.csv"))

presence_events <- amr_clean %>%
  distinct(Tip_ID, AMR_feature)

all_samples <- sample_annotation$Tip_ID
all_features <- sort(unique(presence_events$AMR_feature))

if (length(all_features) == 0) {
  stop("[ERROR] No AMRFinder features detected among the 12D clade genomes.", call. = FALSE)
}

presence_matrix <- matrix(
  0L,
  nrow = length(all_samples),
  ncol = length(all_features),
  dimnames = list(all_samples, all_features)
)

event_row <- match(presence_events$Tip_ID, all_samples)
event_col <- match(presence_events$AMR_feature, all_features)
presence_matrix[cbind(event_row, event_col)] <- 1L

write_csv(
  as.data.frame(presence_matrix) %>% rownames_to_column("Tip_ID"),
  file.path(matrix_dir, "13B_amrfinder_presence_absence_matrix.csv.gz")
)

# ==============================================================================
# 5. AMR burden per genome
# ==============================================================================

sample_burden <- tibble(
  Tip_ID = rownames(presence_matrix),
  AMR_feature_count = rowSums(presence_matrix)
) %>%
  left_join(sample_annotation, by = "Tip_ID") %>%
  arrange(Clade, Tip_ID)

write_csv(sample_burden, file.path(table_dir, "13B_sample_amr_feature_burden_by_clade.csv"))

# ==============================================================================
# 6. Clade prevalence and Fisher/FDR statistics
# ==============================================================================

message("[INFO] Running Fisher exact tests for each AMR feature by clade...")

stat_rows <- vector("list", length(clade_levels) * length(all_features))
k <- 1

for (cl in clade_levels) {
  in_samples <- sample_annotation %>%
    filter(as.character(Clade) == cl) %>%
    pull(Tip_ID)

  out_samples <- sample_annotation %>%
    filter(as.character(Clade) != cl) %>%
    pull(Tip_ID)

  in_idx <- match(in_samples, rownames(presence_matrix))
  out_idx <- match(out_samples, rownames(presence_matrix))

  n_in <- length(in_idx)
  n_out <- length(out_idx)

  for (feature in all_features) {
    feature_idx <- match(feature, colnames(presence_matrix))

    a <- sum(presence_matrix[in_idx, feature_idx], na.rm = TRUE)
    b <- n_in - a
    c <- sum(presence_matrix[out_idx, feature_idx], na.rm = TRUE)
    d <- n_out - c

    tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)

    p_enrich <- safe_fisher(tab, alternative = "greater")
    p_deplete <- safe_fisher(tab, alternative = "less")
    p_two <- safe_fisher(tab, alternative = "two.sided")

    odds_ratio_haldane <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

    stat_rows[[k]] <- tibble(
      Clade = cl,
      AMR_feature = feature,
      N_in_clade = n_in,
      N_out_clade = n_out,
      In_clade_present = a,
      In_clade_absent = b,
      Out_clade_present = c,
      Out_clade_absent = d,
      In_clade_percent = round(100 * a / n_in, 3),
      Out_clade_percent = round(100 * c / n_out, 3),
      Prevalence_difference_percent = round((100 * a / n_in) - (100 * c / n_out), 3),
      Odds_ratio_Haldane = round(odds_ratio_haldane, 5),
      Fisher_p_enrichment = p_enrich,
      Fisher_p_depletion = p_deplete,
      Fisher_p_two_sided = p_two
    )

    k <- k + 1
  }
}

stats <- bind_rows(stat_rows) %>%
  mutate(
    FDR_enrichment_all_tests = p.adjust(Fisher_p_enrichment, method = "BH"),
    FDR_depletion_all_tests = p.adjust(Fisher_p_depletion, method = "BH"),
    FDR_two_sided_all_tests = p.adjust(Fisher_p_two_sided, method = "BH")
  ) %>%
  group_by(Clade) %>%
  mutate(
    FDR_enrichment_within_clade = p.adjust(Fisher_p_enrichment, method = "BH"),
    FDR_depletion_within_clade = p.adjust(Fisher_p_depletion, method = "BH"),
    FDR_two_sided_within_clade = p.adjust(Fisher_p_two_sided, method = "BH")
  ) %>%
  ungroup() %>%
  left_join(feature_info, by = "AMR_feature") %>%
  mutate(
    Enriched_FDR05_all_tests =
      FDR_enrichment_all_tests <= FDR_CUTOFF &
      Prevalence_difference_percent > 0,

    Depleted_FDR05_all_tests =
      FDR_depletion_all_tests <= FDR_CUTOFF &
      Prevalence_difference_percent < 0,

    Enriched_FDR05_within_clade =
      FDR_enrichment_within_clade <= FDR_CUTOFF &
      Prevalence_difference_percent > 0,

    Depleted_FDR05_within_clade =
      FDR_depletion_within_clade <= FDR_CUTOFF &
      Prevalence_difference_percent < 0
  ) %>%
  arrange(Clade, FDR_enrichment_all_tests, desc(Prevalence_difference_percent), AMR_feature)

write_csv(
  stats,
  file.path(table_dir, "13B_amrfinder_gene_clade_fisher_fdr_manifest.csv")
)

write_csv(
  stats %>%
    filter(Enriched_FDR05_all_tests) %>%
    arrange(Clade, FDR_enrichment_all_tests, desc(Prevalence_difference_percent)),
  file.path(table_dir, "13B_amrfinder_gene_clade_enriched_FDR05_all_tests.csv")
)

write_csv(
  stats %>%
    filter(Depleted_FDR05_all_tests) %>%
    arrange(Clade, FDR_depletion_all_tests, Prevalence_difference_percent),
  file.path(table_dir, "13B_amrfinder_gene_clade_depleted_FDR05_all_tests.csv")
)

prevalence_by_clade <- stats %>%
  select(
    Clade,
    AMR_feature,
    N_in_clade,
    In_clade_present,
    In_clade_percent,
    Out_clade_present,
    Out_clade_percent,
    Prevalence_difference_percent,
    Classes,
    Subclasses,
    Types,
    Subtypes
  ) %>%
  arrange(Clade, desc(In_clade_percent), AMR_feature)

write_csv(
  prevalence_by_clade,
  file.path(table_dir, "13B_amrfinder_gene_prevalence_by_clade.csv")
)

# ==============================================================================
# 7. Feature and sample ordering for heatmap
# ==============================================================================

global_feature_summary <- tibble(
  AMR_feature = colnames(presence_matrix),
  Global_present = colSums(presence_matrix),
  Global_percent = round(100 * colSums(presence_matrix) / nrow(presence_matrix), 3)
) %>%
  left_join(
    stats %>%
      group_by(AMR_feature) %>%
      summarise(
        best_FDR_enrichment_all_tests = min(FDR_enrichment_all_tests, na.rm = TRUE),
        max_prevalence_difference_percent = max(Prevalence_difference_percent, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "AMR_feature"
  ) %>%
  arrange(desc(Global_present), best_FDR_enrichment_all_tests, AMR_feature)

write_csv(
  global_feature_summary,
  file.path(table_dir, "13B_amrfinder_global_feature_summary.csv")
)

selected_features <- global_feature_summary %>%
  filter(
    Global_percent >= MIN_GLOBAL_PREVALENCE_FOR_HEATMAP * 100 |
      best_FDR_enrichment_all_tests <= FDR_CUTOFF
  ) %>%
  arrange(desc(Global_present), best_FDR_enrichment_all_tests, AMR_feature) %>%
  slice_head(n = MAX_HEATMAP_FEATURES) %>%
  pull(AMR_feature)

if (length(selected_features) == 0) {
  selected_features <- global_feature_summary %>%
    slice_head(n = min(MAX_HEATMAP_FEATURES, n())) %>%
    pull(AMR_feature)
}

heatmap_matrix <- presence_matrix[, selected_features, drop = FALSE]

if (ncol(heatmap_matrix) > 2) {
  feature_order <- colnames(heatmap_matrix)[hclust(dist(t(heatmap_matrix)))$order]
} else {
  feature_order <- colnames(heatmap_matrix)
}

sample_order <- c()

for (cl in clade_levels) {
  ids <- sample_annotation %>%
    filter(as.character(Clade) == cl) %>%
    pull(Tip_ID)

  ids <- ids[ids %in% rownames(heatmap_matrix)]

  if (length(ids) > 2 && ncol(heatmap_matrix) > 1) {
    submat <- heatmap_matrix[ids, , drop = FALSE]
    hc <- tryCatch(hclust(dist(submat)), error = function(e) NULL)
    if (!is.null(hc)) {
      ids <- ids[hc$order]
    }
  }

  sample_order <- c(sample_order, ids)
}

write_csv(
  tibble(Order = seq_along(sample_order), Tip_ID = sample_order) %>%
    left_join(sample_annotation, by = "Tip_ID"),
  file.path(table_dir, "13B_amrfinder_heatmap_sample_order_clustered_within_clade.csv")
)

write_csv(
  tibble(Order = seq_along(feature_order), AMR_feature = feature_order) %>%
    left_join(feature_info, by = "AMR_feature"),
  file.path(table_dir, "13B_amrfinder_heatmap_feature_order_clustered.csv")
)

# ==============================================================================
# 8. Figures
# ==============================================================================

message("[INFO] Generating figures...")

heatmap_long <- tidyr::expand_grid(
  Tip_ID = sample_order,
  AMR_feature = feature_order
) %>%
  mutate(
    Presence = as.integer(
      heatmap_matrix[
        cbind(match(Tip_ID, rownames(heatmap_matrix)), match(AMR_feature, colnames(heatmap_matrix)))
      ]
    )
  ) %>%
  left_join(sample_annotation %>% select(Tip_ID, Clade), by = "Tip_ID") %>%
  mutate(
    Tip_ID = factor(Tip_ID, levels = rev(sample_order)),
    AMR_feature = factor(AMR_feature, levels = feature_order),
    Clade = factor(Clade, levels = clade_levels)
  )

p_heatmap <- ggplot(heatmap_long, aes(x = AMR_feature, y = Tip_ID, fill = Presence)) +
  geom_tile(color = NA) +
  facet_grid(Clade ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "white", high = "black", breaks = c(0, 1), labels = c("Absent", "Present")) +
  labs(
    title = "AMRFinderPlus gene presence/absence clustered by clade",
    subtitle = paste0("Rows grouped by 12D clade and clustered within clade; feature level = ", FEATURE_LEVEL),
    x = "AMR feature",
    y = "Genomes",
    fill = "Presence"
  ) +
  theme_publication(base_size = 9) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
    panel.spacing.y = unit(0.08, "lines")
  )

save_plot(
  p_heatmap,
  "13B_amrfinder_presence_absence_heatmap_clustered_by_clade",
  width = 14,
  height = 10
)

p_prev <- stats %>%
  filter(AMR_feature %in% feature_order) %>%
  mutate(
    AMR_feature = factor(AMR_feature, levels = feature_order),
    Clade = factor(Clade, levels = clade_levels)
  ) %>%
  ggplot(aes(x = AMR_feature, y = Clade, fill = In_clade_percent)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "white", high = "black", labels = function(x) paste0(x, "%")) +
  labs(
    title = "AMRFinderPlus gene prevalence by clade",
    subtitle = "Cell values represent percent of genomes in each clade carrying the AMR feature",
    x = "AMR feature",
    y = "Clade",
    fill = "Prevalence"
  ) +
  theme_publication(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7)
  )

save_plot(
  p_prev,
  "13B_amrfinder_clade_prevalence_heatmap",
  width = 14,
  height = 5.5
)

enrichment_plot_df <- stats %>%
  filter(Enriched_FDR05_all_tests) %>%
  group_by(Clade) %>%
  arrange(FDR_enrichment_all_tests, desc(Prevalence_difference_percent), .by_group = TRUE) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  mutate(
    AMR_feature = fct_reorder(AMR_feature, Prevalence_difference_percent),
    Clade = factor(Clade, levels = clade_levels),
    neg_log10_FDR = -log10(pmax(FDR_enrichment_all_tests, 1e-300))
  )

if (nrow(enrichment_plot_df) > 0) {
  p_enrich <- ggplot(
    enrichment_plot_df,
    aes(
      x = Clade,
      y = AMR_feature,
      size = Prevalence_difference_percent,
      fill = neg_log10_FDR
    )
  ) +
    geom_point(shape = 21, color = "black", alpha = 0.9) +
    scale_size_continuous(range = c(2, 8)) +
    scale_fill_gradient(low = "white", high = "black") +
    labs(
      title = "AMRFinderPlus genes enriched by clade",
      subtitle = "Fisher exact test, BH-FDR corrected across all gene-by-clade tests",
      x = "Clade",
      y = "AMR feature",
      size = "Prevalence difference (%)",
      fill = "-log10(FDR)"
    ) +
    theme_publication(base_size = 10)

  save_plot(
    p_enrich,
    "13B_amrfinder_clade_enrichment_dotplot_FDR05",
    width = 9,
    height = max(5.5, 0.22 * nrow(enrichment_plot_df))
  )
}

p_burden <- sample_burden %>%
  mutate(Clade = factor(Clade, levels = clade_levels)) %>%
  ggplot(aes(x = Clade, y = AMR_feature_count)) +
  geom_boxplot(outlier.alpha = 0.3, linewidth = 0.35) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "AMRFinderPlus AMR feature burden by clade",
    subtitle = "Each point is one genome; burden is the number of AMRFinderPlus features detected",
    x = "Clade",
    y = "AMR feature count per genome"
  ) +
  theme_publication(base_size = 11)

save_plot(
  p_burden,
  "13B_amrfinder_amr_feature_burden_by_clade",
  width = 8,
  height = 5.5
)

# ==============================================================================
# 9. Manifest
# ==============================================================================

manifest <- tibble(
  output_type = c(
    "input_amrfinder_report",
    "input_clade_file",
    "input_clade_counts",
    "feature_level",
    "main_presence_matrix",
    "main_fdr_manifest",
    "enriched_fdr05_table",
    "depleted_fdr05_table",
    "prevalence_by_clade",
    "sample_burden",
    "figures"
  ),
  path_or_value = c(
    amr_file,
    clade_file,
    clade_counts_file,
    FEATURE_LEVEL,
    file.path(matrix_dir, "13B_amrfinder_presence_absence_matrix.csv.gz"),
    file.path(table_dir, "13B_amrfinder_gene_clade_fisher_fdr_manifest.csv"),
    file.path(table_dir, "13B_amrfinder_gene_clade_enriched_FDR05_all_tests.csv"),
    file.path(table_dir, "13B_amrfinder_gene_clade_depleted_FDR05_all_tests.csv"),
    file.path(table_dir, "13B_amrfinder_gene_prevalence_by_clade.csv"),
    file.path(table_dir, "13B_sample_amr_feature_burden_by_clade.csv"),
    figure_dir
  )
)

write_csv(manifest, file.path(out_dir, "13B_output_manifest.csv"))

# ==============================================================================
# 10. Terminal summary
# ==============================================================================

message("")
message("==========================================================")
message(" SCRIPT 13B COMPLETED SUCCESSFULLY")
message("==========================================================")
message("[OUTPUT] Tables:   ", table_dir)
message("[OUTPUT] Matrices: ", matrix_dir)
message("[OUTPUT] Figures:  ", figure_dir)
message("")
message("=== Input summary ===")
print(tibble(
  target_genomes = nrow(sample_annotation),
  amrfinder_rows_in_target_set = nrow(amr_clean),
  genomes_with_amrfinder_hits = length(amr_samples_found),
  genomes_without_amrfinder_hits = length(amr_samples_missing),
  amr_features_detected = length(all_features),
  feature_level = FEATURE_LEVEL
))
message("")
message("=== Enriched AMR genes by clade, FDR <= 0.05 ===")
print(
  stats %>%
    filter(Enriched_FDR05_all_tests) %>%
    count(Clade, name = "n_enriched_features") %>%
    arrange(Clade)
)
message("")
message("[DONE] AMRFinderPlus clade heatmap and FDR enrichment analysis finished.")
