#!/usr/bin/env Rscript

# ==============================================================================
# Script 12E: Pairwise SNP Distance, Temporal Clade Dynamics, and Clade-Enriched SNPs
#
# Purpose:
#   1. Load Panaroo core genome alignment from Script 12C
#   2. Extract A/C/G/T variable core SNP sites
#   3. Compute pairwise SNP distances among all 704 non-dairy pre-harvest genomes
#   4. Merge SNP distances with Script 12D clade assignments
#   5. Summarize:
#        - global pairwise SNP distances
#        - within-clade SNP distances
#        - between-clade SNP distances
#        - top 3 major clades
#        - top 5 overview clades
#   6. Screen clade temporal behavior:
#        - persistent-like
#        - emerging-like
#        - re-emergence / unsampled-persistence-like
#        - sporadic/minority
#   7. Identify clade-enriched SNP alleles:
#        - strict rule:
#            allele count >=40 within clade AND <1% outside clade
#        - prevalence rule:
#            allele prevalence >=80% within clade AND <1% outside clade
#
# Notes:
#   - SNP distances are computed directly from Panaroo core alignment.
#   - Only A/C/G/T sites are used for SNP distances.
#   - Gaps, N, and ambiguous bases are ignored pairwise.
#   - Clade labels come from Script 12D.
# ==============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(forcats)
})

message("=== SCRIPT 12E: PAIRWISE SNP DISTANCE + CLADE DYNAMICS + ENRICHED SNPS ===")

# ==============================================================================
# 1. Paths and parameters
# ==============================================================================

base_dir <- "/home/samuelajulo/epi/Ecoli_epi"

alignment_filtered <- file.path(
  base_dir,
  "results",
  "12c_all_preharvest_combined",
  "panaroo",
  "core_gene_alignment_filtered.aln"
)

alignment_unfiltered <- file.path(
  base_dir,
  "results",
  "12c_all_preharvest_combined",
  "panaroo",
  "core_gene_alignment.aln"
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

out_dir <- file.path(base_dir, "results", "12e_pairwise_snp_distance")
table_dir <- file.path(out_dir, "tables")
matrix_dir <- file.path(out_dir, "matrices")
figure_dir <- file.path(out_dir, "figures")
enriched_dir <- file.path(out_dir, "clade_enriched_snps")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enriched_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

EXPECTED_N <- 704
TOP_N_MAJOR_CLADES <- 3
TOP_N_OVERVIEW_CLADES <- 5
BLOCK_SIZE <- 50000

STRICT_MIN_IN_CLADE_COUNT <- 40
STRICT_MAX_OUT_CLADE_PREV <- 0.01

PREV_MIN_IN_CLADE_PREV <- 0.80
PREV_MAX_OUT_CLADE_PREV <- 0.01

bases <- c("A", "C", "G", "T")

alignment_file <- if (file.exists(alignment_filtered)) {
  alignment_filtered
} else {
  alignment_unfiltered
}

if (!file.exists(alignment_file)) {
  stop(
    "[ERROR] Core alignment file not found. Expected one of:\n",
    alignment_filtered, "\n",
    alignment_unfiltered,
    call. = FALSE
  )
}

if (!file.exists(clade_file)) {
  stop("[ERROR] 12D clade file not found: ", clade_file, call. = FALSE)
}

if (!file.exists(clade_counts_file)) {
  stop("[ERROR] 12D clade counts file not found: ", clade_counts_file, call. = FALSE)
}

metadata_used <- if (file.exists(metadata_file)) {
  metadata_file
} else {
  metadata_fallback
}

if (!file.exists(metadata_used)) {
  stop("[ERROR] Metadata file not found: ", metadata_used, call. = FALSE)
}

message("[INFO] Alignment file:      ", alignment_file)
message("[INFO] Clade file:          ", clade_file)
message("[INFO] Clade counts file:   ", clade_counts_file)
message("[INFO] Metadata file:       ", metadata_used)
message("[INFO] Output directory:    ", out_dir)

# ==============================================================================
# 2. Helper functions
# ==============================================================================

normalize_tip_id <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  x <- str_remove(x, "\\.gz$")
  x <- str_remove(x, "\\.(fasta|fna|fa|gff|gff3|aln)$")
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

summarise_snp <- function(df) {
  df <- df %>%
    filter(!is.na(SNP_distance), is.finite(SNP_distance))

  if (nrow(df) == 0) {
    return(tibble(
      n_pairs = 0,
      min_snp = NA_real_,
      q1_snp = NA_real_,
      median_snp = NA_real_,
      mean_snp = NA_real_,
      q3_snp = NA_real_,
      max_snp = NA_real_,
      iqr_snp = NA_real_,
      p95_snp = NA_real_,
      p99_snp = NA_real_
    ))
  }

  tibble(
    n_pairs = nrow(df),
    min_snp = min(df$SNP_distance, na.rm = TRUE),
    q1_snp = as.numeric(quantile(df$SNP_distance, 0.25, na.rm = TRUE)),
    median_snp = median(df$SNP_distance, na.rm = TRUE),
    mean_snp = mean(df$SNP_distance, na.rm = TRUE),
    q3_snp = as.numeric(quantile(df$SNP_distance, 0.75, na.rm = TRUE)),
    max_snp = max(df$SNP_distance, na.rm = TRUE),
    iqr_snp = IQR(df$SNP_distance, na.rm = TRUE),
    p95_snp = as.numeric(quantile(df$SNP_distance, 0.95, na.rm = TRUE)),
    p99_snp = as.numeric(quantile(df$SNP_distance, 0.99, na.rm = TRUE))
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

safe_percent <- function(num, denom) {
  ifelse(is.na(denom) | denom == 0, NA_real_, 100 * num / denom)
}

# ==============================================================================
# 3. Read clades and metadata
# ==============================================================================

message("[INFO] Reading 12D clade assignments...")

clades <- read_csv(clade_file, show_col_types = FALSE) %>%
  mutate(
    Tip_ID = normalize_tip_id(Tip_ID),
    Clade = as.character(Clade)
  ) %>%
  distinct(Tip_ID, .keep_all = TRUE)

clade_counts <- read_csv(clade_counts_file, show_col_types = FALSE) %>%
  mutate(Clade = as.character(Clade)) %>%
  arrange(desc(Isolate_Count), Clade)

major_clades <- clade_counts %>%
  slice_head(n = TOP_N_MAJOR_CLADES) %>%
  pull(Clade)

overview_clades <- clade_counts %>%
  slice_head(n = TOP_N_OVERVIEW_CLADES) %>%
  pull(Clade)

message("[INFO] Number of clades from 12D: ", nrow(clade_counts))
message("[INFO] Top 3 major clades: ", paste(major_clades, collapse = ", "))
message("[INFO] Top 5 overview clades: ", paste(overview_clades, collapse = ", "))

write_csv(
  clade_counts,
  file.path(table_dir, "12E_input_clade_counts_from_12D.csv")
)

write_csv(
  tibble(Rank = seq_along(major_clades), Major_Clade = major_clades),
  file.path(table_dir, "12E_top3_major_clades.csv")
)

metadata <- read_csv(metadata_used, show_col_types = FALSE)

assembly_col <- pick_first_existing_col(
  metadata,
  c("Assembly", "Tip_ID", "Isolate", "Strain", "Sample", "sample")
)

if (is.na(assembly_col)) {
  stop("[ERROR] Could not identify Assembly/Tip_ID column in metadata.", call. = FALSE)
}

year_col <- pick_first_existing_col(
  metadata,
  c("Collection_year", "Collection year", "collection_year", "Year", "year", "Collection date", "collection_date")
)

metadata_clean <- metadata %>%
  mutate(
    Tip_ID = normalize_tip_id(.data[[assembly_col]]),
    Collection_year_num = if (!is.na(year_col)) extract_year_num(.data[[year_col]]) else NA_integer_
  ) %>%
  distinct(Tip_ID, .keep_all = TRUE)

sample_annotation <- clades %>%
  left_join(metadata_clean, by = "Tip_ID")

write_csv(
  sample_annotation,
  file.path(table_dir, "12E_sample_clade_metadata_annotation.csv")
)

# ==============================================================================
# 4. Read core alignment
# ==============================================================================

message("[INFO] Reading core alignment with ape::read.dna...")

dna <- read.dna(alignment_file, format = "fasta", as.character = TRUE)

if (!is.matrix(dna)) {
  dna <- as.matrix(dna)
}

raw_ids <- rownames(dna)
clean_ids <- normalize_tip_id(raw_ids)

if (anyDuplicated(clean_ids) > 0) {
  dup_ids <- clean_ids[duplicated(clean_ids)]
  stop(
    "[ERROR] Duplicate normalized alignment IDs detected. Example: ",
    paste(head(unique(dup_ids), 10), collapse = ", "),
    call. = FALSE
  )
}

rownames(dna) <- clean_ids
dna <- dna[order(rownames(dna)), , drop = FALSE]

n_genomes <- nrow(dna)
alignment_length <- ncol(dna)

message("[INFO] Genomes in alignment: ", n_genomes)
message("[INFO] Alignment length: ", alignment_length)

if (n_genomes != EXPECTED_N) {
  warning("[WARN] Expected ", EXPECTED_N, " genomes, but alignment has ", n_genomes, ".")
}

missing_clade <- setdiff(rownames(dna), clades$Tip_ID)
missing_alignment <- setdiff(clades$Tip_ID, rownames(dna))

if (length(missing_clade) > 0) {
  write_lines(missing_clade, file.path(log_dir, "alignment_tips_missing_from_12D_clades.txt"))
  stop("[ERROR] Some alignment tips are missing from 12D clades. See logs.", call. = FALSE)
}

if (length(missing_alignment) > 0) {
  write_lines(missing_alignment, file.path(log_dir, "12D_clade_tips_missing_from_alignment.txt"))
  stop("[ERROR] Some 12D clade tips are missing from alignment. See logs.", call. = FALSE)
}

dna <- toupper(dna)

# ==============================================================================
# 5. Identify core SNP sites
# ==============================================================================

message("[INFO] Identifying variable A/C/G/T SNP sites...")

var_sites <- logical(alignment_length)
block_starts <- seq(1, alignment_length, by = BLOCK_SIZE)

for (block_start in block_starts) {
  block_end <- min(block_start + BLOCK_SIZE - 1, alignment_length)
  idx <- block_start:block_end

  block <- dna[, idx, drop = FALSE]
  encoded_block <- matrix(
    match(as.vector(block), bases),
    nrow = n_genomes,
    ncol = length(idx)
  )

  allele_presence_count <-
    (colSums(encoded_block == 1, na.rm = TRUE) > 0) +
    (colSums(encoded_block == 2, na.rm = TRUE) > 0) +
    (colSums(encoded_block == 3, na.rm = TRUE) > 0) +
    (colSums(encoded_block == 4, na.rm = TRUE) > 0)

  var_sites[idx] <- allele_presence_count >= 2

  message("[INFO] SNP scan block ", block_start, "-", block_end, " complete.")
}

snp_site_idx <- which(var_sites)
n_snp_sites <- length(snp_site_idx)

message("[INFO] Core SNP sites detected: ", n_snp_sites)

snp_site_summary <- tibble(
  alignment_file = alignment_file,
  n_genomes = n_genomes,
  alignment_length_bp = alignment_length,
  core_snp_sites_acgt = n_snp_sites,
  percent_sites_variable_acgt = round((n_snp_sites / alignment_length) * 100, 4),
  ambiguous_or_gap_rule = "Pairwise SNP distances count only A/C/G/T mismatches; gaps, N, and ambiguous bases are ignored pairwise."
)

write_csv(
  snp_site_summary,
  file.path(table_dir, "12E_snp_site_summary.csv")
)

if (n_snp_sites == 0) {
  stop("[ERROR] No SNP sites detected in the alignment.", call. = FALSE)
}

# ==============================================================================
# 6. Build SNP-only encoded matrix
# ==============================================================================

message("[INFO] Building SNP-only encoded A/C/G/T matrix...")

snp_dna <- dna[, snp_site_idx, drop = FALSE]

encoded_snp <- matrix(
  match(as.vector(snp_dna), bases),
  nrow = n_genomes,
  ncol = n_snp_sites
)

rownames(encoded_snp) <- rownames(dna)
colnames(encoded_snp) <- as.character(snp_site_idx)

valid_snp <- !is.na(encoded_snp)

# ==============================================================================
# 7. Compute pairwise SNP distance matrix
# ==============================================================================

message("[INFO] Computing pairwise SNP distance matrix in blocks...")

same_counts <- matrix(0, nrow = n_genomes, ncol = n_genomes)
valid_counts <- matrix(0, nrow = n_genomes, ncol = n_genomes)

snp_block_starts <- seq(1, n_snp_sites, by = BLOCK_SIZE)

for (block_start in snp_block_starts) {
  block_end <- min(block_start + BLOCK_SIZE - 1, n_snp_sites)
  idx <- block_start:block_end

  encoded_block <- encoded_snp[, idx, drop = FALSE]
  valid_block <- !is.na(encoded_block)

  valid_counts <- valid_counts + tcrossprod(valid_block * 1)

  for (a in seq_along(bases)) {
    allele_mat <- (!is.na(encoded_block)) & (encoded_block == a)
    same_counts <- same_counts + tcrossprod(allele_mat * 1)
  }

  message("[INFO] Pairwise block ", block_start, "-", block_end, " of ", n_snp_sites, " SNP sites complete.")
}

snp_matrix <- valid_counts - same_counts
diag(snp_matrix) <- 0
rownames(snp_matrix) <- rownames(dna)
colnames(snp_matrix) <- rownames(dna)

valid_site_matrix <- valid_counts
diag(valid_site_matrix) <- n_snp_sites
rownames(valid_site_matrix) <- rownames(dna)
colnames(valid_site_matrix) <- rownames(dna)

message("[INFO] Pairwise SNP matrix complete.")

write_csv(
  as.data.frame(snp_matrix) %>% tibble::rownames_to_column("Tip_ID"),
  file.path(matrix_dir, "12E_pairwise_snp_distance_matrix_all_704.csv.gz")
)

write_csv(
  as.data.frame(valid_site_matrix) %>% tibble::rownames_to_column("Tip_ID"),
  file.path(matrix_dir, "12E_pairwise_valid_site_count_matrix_all_704.csv.gz")
)

# ==============================================================================
# 8. Long-form pairwise table
# ==============================================================================

message("[INFO] Building long-form pairwise SNP table...")

pair_idx <- which(lower.tri(snp_matrix), arr.ind = TRUE)

pairwise_long <- tibble(
  Tip_ID_1 = rownames(snp_matrix)[pair_idx[, "row"]],
  Tip_ID_2 = colnames(snp_matrix)[pair_idx[, "col"]],
  SNP_distance = as.integer(round(snp_matrix[pair_idx])),
  Valid_sites_compared = as.integer(round(valid_site_matrix[pair_idx]))
) %>%
  left_join(
    clades %>% select(Tip_ID, Clade) %>% rename(Clade_1 = Clade),
    by = c("Tip_ID_1" = "Tip_ID")
  ) %>%
  left_join(
    clades %>% select(Tip_ID, Clade) %>% rename(Clade_2 = Clade),
    by = c("Tip_ID_2" = "Tip_ID")
  ) %>%
  left_join(
    metadata_clean %>% select(Tip_ID, Collection_year_num) %>% rename(Collection_year_1 = Collection_year_num),
    by = c("Tip_ID_1" = "Tip_ID")
  ) %>%
  left_join(
    metadata_clean %>% select(Tip_ID, Collection_year_num) %>% rename(Collection_year_2 = Collection_year_num),
    by = c("Tip_ID_2" = "Tip_ID")
  ) %>%
  mutate(
    Pair_type = if_else(Clade_1 == Clade_2, "Within-clade", "Between-clade"),
    Clade_pair_A = pmin(Clade_1, Clade_2),
    Clade_pair_B = pmax(Clade_1, Clade_2),
    Clade_pair = paste(Clade_pair_A, Clade_pair_B, sep = " vs "),
    Year_gap = abs(Collection_year_1 - Collection_year_2),
    Year_gap_bin = case_when(
      is.na(Year_gap) ~ "Missing year",
      Year_gap == 0 ~ "Same year",
      Year_gap <= 3 ~ "1-3 years",
      Year_gap <= 9 ~ "4-9 years",
      Year_gap >= 10 ~ ">=10 years",
      TRUE ~ "Missing year"
    ),
    Year_gap_bin = factor(
      Year_gap_bin,
      levels = c("Same year", "1-3 years", "4-9 years", ">=10 years", "Missing year")
    )
  )

pairwise_quality_summary <- tibble(
  total_pairwise_comparisons_before_filter = nrow(pairwise_long),
  nonfinite_snp_distances = sum(is.na(pairwise_long$SNP_distance) | !is.finite(pairwise_long$SNP_distance)),
  zero_valid_site_comparisons = sum(pairwise_long$Valid_sites_compared <= 0, na.rm = TRUE)
)

write_csv(
  pairwise_quality_summary,
  file.path(table_dir, "12E_pairwise_distance_quality_summary.csv")
)

pairwise_long <- pairwise_long %>%
  filter(
    !is.na(SNP_distance),
    is.finite(SNP_distance),
    Valid_sites_compared > 0
  )

write_csv(
  pairwise_long,
  file.path(table_dir, "12E_pairwise_snp_distances_all_pairs.csv.gz")
)

message("[INFO] Pairwise comparisons after filtering: ", nrow(pairwise_long))

# ==============================================================================
# 9. SNP distance summaries
# ==============================================================================

message("[INFO] Summarizing SNP distances...")

global_summary <- pairwise_long %>%
  summarise_snp() %>%
  mutate(
    comparison_scope = "All genome pairs",
    n_genomes = n_genomes,
    n_clades = n_distinct(clades$Clade),
    .before = 1
  )

within_summary <- pairwise_long %>%
  filter(Pair_type == "Within-clade") %>%
  group_by(Clade_1) %>%
  group_modify(~summarise_snp(.x)) %>%
  ungroup() %>%
  rename(Clade = Clade_1) %>%
  left_join(clade_counts, by = "Clade") %>%
  relocate(Clade, Isolate_Count, n_pairs) %>%
  arrange(desc(Isolate_Count), Clade)

between_summary <- pairwise_long %>%
  filter(Pair_type == "Between-clade") %>%
  group_by(Clade_pair_A, Clade_pair_B, Clade_pair) %>%
  group_modify(~summarise_snp(.x)) %>%
  ungroup() %>%
  arrange(median_snp, mean_snp)

near_summary <- pairwise_long %>%
  filter(Pair_type == "Within-clade") %>%
  group_by(Clade_1) %>%
  summarise(
    n_pairs = n(),
    pairs_le_5_snp = sum(SNP_distance <= 5, na.rm = TRUE),
    pairs_le_10_snp = sum(SNP_distance <= 10, na.rm = TRUE),
    pairs_le_25_snp = sum(SNP_distance <= 25, na.rm = TRUE),
    pairs_le_50_snp = sum(SNP_distance <= 50, na.rm = TRUE),
    percent_pairs_le_5_snp = round(100 * pairs_le_5_snp / n_pairs, 2),
    percent_pairs_le_10_snp = round(100 * pairs_le_10_snp / n_pairs, 2),
    percent_pairs_le_25_snp = round(100 * pairs_le_25_snp / n_pairs, 2),
    percent_pairs_le_50_snp = round(100 * pairs_le_50_snp / n_pairs, 2),
    .groups = "drop"
  ) %>%
  rename(Clade = Clade_1) %>%
  left_join(clade_counts, by = "Clade") %>%
  relocate(Clade, Isolate_Count, n_pairs) %>%
  arrange(desc(Isolate_Count), Clade)

write_csv(global_summary, file.path(table_dir, "12E_global_pairwise_snp_summary.csv"))
write_csv(within_summary, file.path(table_dir, "12E_within_clade_pairwise_snp_summary_all_clades.csv"))
write_csv(between_summary, file.path(table_dir, "12E_between_clade_pairwise_snp_summary_all_clades.csv"))
write_csv(near_summary, file.path(table_dir, "12E_near_identical_within_clade_pair_counts.csv"))

# ==============================================================================
# 10. Top 3 major clades and top 5 overview clades
# ==============================================================================

major_pairwise_within <- pairwise_long %>%
  filter(
    Pair_type == "Within-clade",
    Clade_1 %in% major_clades
  )

major_within_summary <- major_pairwise_within %>%
  group_by(Clade_1) %>%
  group_modify(~summarise_snp(.x)) %>%
  ungroup() %>%
  rename(Clade = Clade_1) %>%
  left_join(clade_counts, by = "Clade") %>%
  relocate(Clade, Isolate_Count, n_pairs) %>%
  arrange(match(Clade, major_clades))

major_between_summary <- pairwise_long %>%
  filter(
    Clade_pair_A %in% major_clades,
    Clade_pair_B %in% major_clades,
    Clade_pair_A != Clade_pair_B
  ) %>%
  group_by(Clade_pair_A, Clade_pair_B, Clade_pair) %>%
  group_modify(~summarise_snp(.x)) %>%
  ungroup() %>%
  arrange(median_snp, mean_snp)

major_sample_metadata <- sample_annotation %>%
  filter(Clade %in% major_clades)

production_summary <- if ("Production_system" %in% names(major_sample_metadata)) {
  major_sample_metadata %>%
    count(Clade, Production_system, name = "n") %>%
    group_by(Clade) %>%
    mutate(
      clade_total = sum(n),
      percent = round(100 * n / clade_total, 2)
    ) %>%
    ungroup()
} else {
  tibble()
}

temporal_summary <- major_sample_metadata %>%
  filter(!is.na(Collection_year_num)) %>%
  count(Clade, Collection_year_num, name = "n") %>%
  group_by(Clade) %>%
  mutate(
    clade_total_with_year = sum(n),
    percent = round(100 * n / clade_total_with_year, 2)
  ) %>%
  ungroup()

overview_pair_summary <- pairwise_long %>%
  filter(Clade_pair_A %in% overview_clades, Clade_pair_B %in% overview_clades) %>%
  group_by(Clade_pair_A, Clade_pair_B) %>%
  summarise(
    n_pairs = n(),
    median_snp = median(SNP_distance, na.rm = TRUE),
    mean_snp = mean(SNP_distance, na.rm = TRUE),
    min_snp = min(SNP_distance, na.rm = TRUE),
    max_snp = max(SNP_distance, na.rm = TRUE),
    .groups = "drop"
  )

overview_matrix <- overview_pair_summary %>%
  select(Clade_pair_A, Clade_pair_B, median_snp) %>%
  pivot_wider(
    names_from = Clade_pair_B,
    values_from = median_snp
  )

write_csv(major_pairwise_within, file.path(table_dir, "12E_top3_major_clades_within_pairwise_snp_distances.csv.gz"))
write_csv(major_within_summary, file.path(table_dir, "12E_top3_major_clades_within_snp_summary.csv"))
write_csv(major_between_summary, file.path(table_dir, "12E_top3_major_clades_between_snp_summary.csv"))
write_csv(production_summary, file.path(table_dir, "12E_top3_major_clades_production_system_composition.csv"))
write_csv(temporal_summary, file.path(table_dir, "12E_top3_major_clades_temporal_counts.csv"))
write_csv(overview_pair_summary, file.path(table_dir, "12E_top5_clades_pairwise_median_snp_summary.csv"))
write_csv(overview_matrix, file.path(matrix_dir, "12E_top5_clades_median_snp_matrix.csv"))

# ==============================================================================
# 11. Temporal persistence / emergence / re-emergence screen
# ==============================================================================

message("[INFO] Building clade temporal persistence/emergence screen...")

year_presence_by_clade <- sample_annotation %>%
  filter(!is.na(Collection_year_num)) %>%
  count(Clade, Collection_year_num, name = "n") %>%
  arrange(Clade, Collection_year_num)

write_csv(
  year_presence_by_clade,
  file.path(table_dir, "12E_all_clades_yearly_genome_counts.csv")
)

clade_temporal_summary <- year_presence_by_clade %>%
  group_by(Clade) %>%
  summarise(
    first_year = min(Collection_year_num, na.rm = TRUE),
    last_year = max(Collection_year_num, na.rm = TRUE),
    year_span = last_year - first_year,
    n_years_detected = n_distinct(Collection_year_num),
    total_with_year = sum(n),
    detected_years = paste(sort(unique(Collection_year_num)), collapse = ";"),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    longest_detection_gap = {
      yrs <- as.integer(unlist(strsplit(detected_years, ";")))
      if (length(yrs) <= 1) 0 else max(diff(sort(yrs))) - 1
    }
  ) %>%
  ungroup() %>%
  left_join(clade_counts, by = "Clade") %>%
  left_join(
    within_summary %>%
      select(
        Clade,
        within_n_pairs = n_pairs,
        within_median_snp = median_snp,
        within_mean_snp = mean_snp,
        within_min_snp = min_snp,
        within_max_snp = max_snp,
        within_iqr_snp = iqr_snp,
        within_p95_snp = p95_snp
      ),
    by = "Clade"
  ) %>%
  left_join(
    near_summary %>%
      select(
        Clade,
        pairs_le_5_snp,
        pairs_le_10_snp,
        pairs_le_25_snp,
        pairs_le_50_snp,
        percent_pairs_le_5_snp,
        percent_pairs_le_10_snp,
        percent_pairs_le_25_snp,
        percent_pairs_le_50_snp
      ),
    by = "Clade"
  ) %>%
  mutate(
    clade_size_class = if_else(Clade %in% major_clades, "Major clade", "Minority clade"),
    temporal_pattern_screen = case_when(
      Isolate_Count < 10 | n_years_detected <= 2 ~ "Sporadic/minority detection",
      first_year >= 2018 & n_years_detected >= 3 ~ "Emerging-like detection",
      longest_detection_gap >= 5 & n_years_detected >= 3 ~ "Re-emergence or unsampled-persistence-like",
      year_span >= 5 & n_years_detected >= 4 ~ "Persistent-like detection",
      TRUE ~ "Limited temporal evidence"
    ),
    interpretation_caution = "Pattern is based on public genome sampling; absence in a year means not observed, not true biological absence."
  ) %>%
  arrange(desc(Isolate_Count), Clade)

top3_temporal_screen <- clade_temporal_summary %>%
  filter(Clade %in% major_clades) %>%
  mutate(Clade = factor(Clade, levels = major_clades)) %>%
  arrange(Clade)

within_pairwise_temporal_gap_summary <- pairwise_long %>%
  filter(Pair_type == "Within-clade") %>%
  group_by(Clade_1, Year_gap_bin) %>%
  group_modify(~summarise_snp(.x)) %>%
  ungroup() %>%
  rename(Clade = Clade_1) %>%
  left_join(clade_counts, by = "Clade") %>%
  relocate(Clade, Isolate_Count, Year_gap_bin, n_pairs) %>%
  arrange(desc(Isolate_Count), Clade, Year_gap_bin)

top3_within_pairwise_temporal_gap_summary <- within_pairwise_temporal_gap_summary %>%
  filter(Clade %in% major_clades)

write_csv(clade_temporal_summary, file.path(table_dir, "12E_clade_temporal_persistence_emergence_screen_all_clades.csv"))
write_csv(top3_temporal_screen, file.path(table_dir, "12E_top3_major_clades_temporal_persistence_emergence_screen.csv"))
write_csv(within_pairwise_temporal_gap_summary, file.path(table_dir, "12E_within_clade_snp_distance_by_temporal_gap_summary.csv"))
write_csv(top3_within_pairwise_temporal_gap_summary, file.path(table_dir, "12E_top3_major_clades_snp_distance_by_temporal_gap_summary.csv"))

# ==============================================================================
# 12. Clade-enriched SNP alleles
# ==============================================================================

message("[INFO] Identifying clade-enriched SNP alleles...")

alignment_clade_map <- tibble(
  Tip_ID = rownames(encoded_snp)
) %>%
  left_join(clades %>% select(Tip_ID, Clade), by = "Tip_ID")

if (any(is.na(alignment_clade_map$Clade))) {
  stop("[ERROR] Missing clade labels for some alignment tips.", call. = FALSE)
}

all_clades <- sort(unique(alignment_clade_map$Clade))
enriched_rows <- list()
row_counter <- 1

for (cl in all_clades) {
  message("[INFO] Enriched SNP scan for ", cl)

  in_idx <- which(alignment_clade_map$Clade == cl)
  out_idx <- which(alignment_clade_map$Clade != cl)

  in_valid <- colSums(valid_snp[in_idx, , drop = FALSE])
  out_valid <- colSums(valid_snp[out_idx, , drop = FALSE])
  global_valid <- in_valid + out_valid

  for (allele_code in seq_along(bases)) {
    allele <- bases[allele_code]

    allele_mat <- (!is.na(encoded_snp)) & (encoded_snp == allele_code)

    in_count <- colSums(allele_mat[in_idx, , drop = FALSE])
    out_count <- colSums(allele_mat[out_idx, , drop = FALSE])
    global_count <- in_count + out_count

    in_prev <- ifelse(in_valid > 0, in_count / in_valid, NA_real_)
    out_prev <- ifelse(out_valid > 0, out_count / out_valid, NA_real_)
    global_prev <- ifelse(global_valid > 0, global_count / global_valid, NA_real_)

    strict_hit <- in_count >= STRICT_MIN_IN_CLADE_COUNT &
      !is.na(out_prev) &
      out_prev < STRICT_MAX_OUT_CLADE_PREV

    prevalence_hit <- !is.na(in_prev) &
      in_prev >= PREV_MIN_IN_CLADE_PREV &
      !is.na(out_prev) &
      out_prev < PREV_MAX_OUT_CLADE_PREV

    hit <- strict_hit | prevalence_hit

    if (any(hit)) {
      hit_idx <- which(hit)

      enriched_rows[[row_counter]] <- tibble(
        Clade = cl,
        Alignment_position = snp_site_idx[hit_idx],
        SNP_site_rank = hit_idx,
        Enriched_allele = allele,
        In_clade_allele_count = as.integer(in_count[hit_idx]),
        In_clade_total_valid = as.integer(in_valid[hit_idx]),
        In_clade_prevalence = round(in_prev[hit_idx], 5),
        Out_clade_allele_count = as.integer(out_count[hit_idx]),
        Out_clade_total_valid = as.integer(out_valid[hit_idx]),
        Out_clade_prevalence = round(out_prev[hit_idx], 5),
        Global_allele_count = as.integer(global_count[hit_idx]),
        Global_total_valid = as.integer(global_valid[hit_idx]),
        Global_allele_prevalence = round(global_prev[hit_idx], 5),
        Strict_count40_outlt1pct = strict_hit[hit_idx],
        Prevalence80_outlt1pct = prevalence_hit[hit_idx],
        Rule = case_when(
          strict_hit[hit_idx] & prevalence_hit[hit_idx] ~ "strict_and_prevalence",
          strict_hit[hit_idx] ~ "strict_count40_outlt1pct",
          prevalence_hit[hit_idx] ~ "prevalence80_outlt1pct",
          TRUE ~ "none"
        )
      )

      row_counter <- row_counter + 1
    }
  }
}

if (length(enriched_rows) > 0) {
  enriched_snps_all <- bind_rows(enriched_rows) %>%
    arrange(Clade, Alignment_position, Enriched_allele)
} else {
  enriched_snps_all <- tibble(
    Clade = character(),
    Alignment_position = integer(),
    SNP_site_rank = integer(),
    Enriched_allele = character(),
    In_clade_allele_count = integer(),
    In_clade_total_valid = integer(),
    In_clade_prevalence = numeric(),
    Out_clade_allele_count = integer(),
    Out_clade_total_valid = integer(),
    Out_clade_prevalence = numeric(),
    Global_allele_count = integer(),
    Global_total_valid = integer(),
    Global_allele_prevalence = numeric(),
    Strict_count40_outlt1pct = logical(),
    Prevalence80_outlt1pct = logical(),
    Rule = character()
  )
}

strict_enriched <- enriched_snps_all %>%
  filter(Strict_count40_outlt1pct)

prevalence_enriched <- enriched_snps_all %>%
  filter(Prevalence80_outlt1pct)

enriched_summary <- clade_counts %>%
  left_join(
    strict_enriched %>%
      count(Clade, name = "strict_enriched_snp_count"),
    by = "Clade"
  ) %>%
  left_join(
    prevalence_enriched %>%
      count(Clade, name = "prevalence80_enriched_snp_count"),
    by = "Clade"
  ) %>%
  mutate(
    strict_enriched_snp_count = replace_na(strict_enriched_snp_count, 0L),
    prevalence80_enriched_snp_count = replace_na(prevalence80_enriched_snp_count, 0L),
    strict_rule = paste0("in-clade allele count >= ", STRICT_MIN_IN_CLADE_COUNT, "; out-of-clade prevalence < ", STRICT_MAX_OUT_CLADE_PREV * 100, "%"),
    prevalence_rule = paste0("in-clade allele prevalence >= ", PREV_MIN_IN_CLADE_PREV * 100, "%; out-of-clade prevalence < ", PREV_MAX_OUT_CLADE_PREV * 100, "%")
  ) %>%
  arrange(desc(Isolate_Count), Clade)

top3_enriched_summary <- enriched_summary %>%
  filter(Clade %in% major_clades) %>%
  mutate(Clade = factor(Clade, levels = major_clades)) %>%
  arrange(Clade)

write_csv(enriched_snps_all, file.path(enriched_dir, "12E_clade_enriched_snps_all_rules_all_clades.csv.gz"))
write_csv(strict_enriched, file.path(enriched_dir, "12E_clade_enriched_snps_strict_count40_outlt1pct_all_clades.csv.gz"))
write_csv(prevalence_enriched, file.path(enriched_dir, "12E_clade_enriched_snps_prevalence80_outlt1pct_all_clades.csv.gz"))
write_csv(enriched_summary, file.path(enriched_dir, "12E_clade_enriched_snp_summary_by_clade.csv"))
write_csv(top3_enriched_summary, file.path(enriched_dir, "12E_top3_major_clades_enriched_snp_summary.csv"))

write_csv(enriched_summary, file.path(table_dir, "12E_clade_enriched_snp_summary_by_clade.csv"))
write_csv(top3_enriched_summary, file.path(table_dir, "12E_top3_major_clades_enriched_snp_summary.csv"))

# ==============================================================================
# 13. Figures
# ==============================================================================

message("[INFO] Generating figures...")

plot_within_all <- pairwise_long %>%
  filter(Pair_type == "Within-clade") %>%
  left_join(clade_counts, by = c("Clade_1" = "Clade")) %>%
  mutate(
    Clade_plot = fct_reorder(Clade_1, Isolate_Count, .desc = TRUE)
  ) %>%
  ggplot(aes(x = Clade_plot, y = SNP_distance)) +
  geom_boxplot(outlier.alpha = 0.25, linewidth = 0.35) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Within-clade pairwise SNP distances",
    subtitle = "All clades from the 704 non-dairy cattle pre-harvest core-genome alignment",
    x = "Clade",
    y = "Pairwise SNP distance"
  ) +
  theme_publication(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(plot_within_all, "12E_within_clade_snp_distance_boxplot_all_clades", width = 12, height = 6.5)

plot_major <- major_pairwise_within %>%
  mutate(Clade_1 = factor(Clade_1, levels = major_clades)) %>%
  ggplot(aes(x = Clade_1, y = SNP_distance)) +
  geom_violin(trim = FALSE, alpha = 0.45, linewidth = 0.25) +
  geom_boxplot(width = 0.18, outlier.alpha = 0.25, linewidth = 0.35) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Pairwise SNP distances within the top 3 major clades",
    subtitle = "Distances computed directly from Panaroo core SNP sites",
    x = "Major clade",
    y = "Pairwise SNP distance"
  ) +
  theme_publication(base_size = 12)

save_plot(plot_major, "12E_top3_major_clades_within_snp_distance_violin_boxplot", width = 9, height = 6)

major_comparison_data <- pairwise_long %>%
  filter(
    Clade_pair_A %in% major_clades,
    Clade_pair_B %in% major_clades
  ) %>%
  mutate(
    Comparison = if_else(
      Clade_pair_A == Clade_pair_B,
      paste0(Clade_pair_A, " within"),
      paste0(Clade_pair_A, " vs ", Clade_pair_B)
    )
  )

if (nrow(major_comparison_data) > 0) {
  comparison_levels <- major_comparison_data %>%
    group_by(Comparison) %>%
    summarise(median_snp = median(SNP_distance, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_snp) %>%
    pull(Comparison)

  plot_between_major <- major_comparison_data %>%
    mutate(Comparison = factor(Comparison, levels = comparison_levels)) %>%
    ggplot(aes(x = Comparison, y = SNP_distance)) +
    geom_boxplot(outlier.alpha = 0.25, linewidth = 0.35) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Within- and between-clade SNP distances for the top 3 major clades",
      subtitle = "Within-clade and between-clade genetic distances from the core alignment",
      x = "Comparison",
      y = "Pairwise SNP distance"
    ) +
    theme_publication(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(plot_between_major, "12E_top3_major_clades_within_between_snp_distance_boxplot", width = 10, height = 6)
}

plot_top5_heatmap <- overview_pair_summary %>%
  mutate(
    Clade_pair_A = factor(Clade_pair_A, levels = overview_clades),
    Clade_pair_B = factor(Clade_pair_B, levels = overview_clades)
  ) %>%
  ggplot(aes(x = Clade_pair_A, y = Clade_pair_B, fill = median_snp)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = round(median_snp, 0)), size = 3) +
  scale_fill_gradient(low = "white", high = "black", labels = comma) +
  labs(
    title = "Median pairwise SNP distance among the top 5 clades",
    subtitle = "Top 3 major clades plus 2 minority/overview clades by isolate count",
    x = "Clade",
    y = "Clade",
    fill = "Median SNPs"
  ) +
  theme_publication(base_size = 11)

save_plot(plot_top5_heatmap, "12E_top5_clades_median_snp_distance_heatmap", width = 7.5, height = 6.5)

if (nrow(production_summary) > 0) {
  plot_prod <- production_summary %>%
    mutate(Clade = factor(Clade, levels = major_clades)) %>%
    ggplot(aes(x = Clade, y = n, fill = Production_system)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2) +
    geom_text(
      aes(label = paste0(n, " (", percent, "%)")),
      position = position_stack(vjust = 0.5),
      size = 3
    ) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Production-system composition of the top 3 major clades",
      subtitle = "Beef vs cattle-unspecified among non-dairy pre-harvest genomes",
      x = "Major clade",
      y = "Genome count",
      fill = "Production system"
    ) +
    theme_publication(base_size = 11)

  save_plot(plot_prod, "12E_top3_major_clades_production_system_composition", width = 8.5, height = 5.5)
}

if (nrow(year_presence_by_clade) > 0) {
  plot_top3_year_counts <- year_presence_by_clade %>%
    filter(Clade %in% major_clades) %>%
    mutate(Clade = factor(Clade, levels = major_clades)) %>%
    ggplot(aes(x = Collection_year_num, y = n)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = n), vjust = -0.25, size = 2.6) +
    facet_wrap(~ Clade, ncol = 1, scales = "free_y") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "Observed temporal distribution of the top 3 major clades",
      subtitle = "Used to screen for persistent-like, emerging-like, and re-emergence-like detection patterns",
      x = "Collection year",
      y = "Genome count"
    ) +
    theme_publication(base_size = 11)

  save_plot(plot_top3_year_counts, "12E_top3_major_clades_observed_yearly_counts", width = 10, height = 8)

  plot_top5_year_counts <- year_presence_by_clade %>%
    filter(Clade %in% overview_clades) %>%
    mutate(Clade = factor(Clade, levels = overview_clades)) %>%
    ggplot(aes(x = Collection_year_num, y = n)) +
    geom_col(width = 0.8) +
    facet_wrap(~ Clade, ncol = 1, scales = "free_y") +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Observed temporal distribution of the top 5 clades",
      subtitle = "Top 3 major clades plus 2 minority/overview clades",
      x = "Collection year",
      y = "Genome count"
    ) +
    theme_publication(base_size = 11)

  save_plot(plot_top5_year_counts, "12E_top5_clades_observed_yearly_counts", width = 10, height = 9)
}

plot_top3_gap <- pairwise_long %>%
  filter(
    Pair_type == "Within-clade",
    Clade_1 %in% major_clades,
    !is.na(Year_gap_bin),
    Year_gap_bin != "Missing year"
  ) %>%
  mutate(
    Clade_1 = factor(Clade_1, levels = major_clades),
    Year_gap_bin = factor(
      as.character(Year_gap_bin),
      levels = c("Same year", "1-3 years", "4-9 years", ">=10 years")
    )
  ) %>%
  ggplot(aes(x = Year_gap_bin, y = SNP_distance)) +
  geom_boxplot(outlier.alpha = 0.25, linewidth = 0.35) +
  facet_wrap(~ Clade_1, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Within-clade SNP distance by temporal separation",
    subtitle = "Low SNP distances across long year gaps support persistence-like or unsampled-continuity patterns",
    x = "Temporal separation between genome pair",
    y = "Pairwise SNP distance"
  ) +
  theme_publication(base_size = 11)

save_plot(plot_top3_gap, "12E_top3_major_clades_snp_distance_by_temporal_gap", width = 10, height = 8)

plot_top3_snp_yeargap_scatter <- pairwise_long %>%
  filter(
    Pair_type == "Within-clade",
    Clade_1 %in% major_clades,
    !is.na(Year_gap)
  ) %>%
  mutate(Clade_1 = factor(Clade_1, levels = major_clades)) %>%
  ggplot(aes(x = Year_gap, y = SNP_distance)) +
  geom_point(alpha = 0.18, size = 0.7) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.7) +
  facet_wrap(~ Clade_1, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Relationship between temporal gap and SNP distance in major clades",
    subtitle = "Each point is one within-clade genome pair",
    x = "Absolute year gap between genome pair",
    y = "Pairwise SNP distance"
  ) +
  theme_publication(base_size = 11)

save_plot(plot_top3_snp_yeargap_scatter, "12E_top3_major_clades_snp_distance_vs_year_gap", width = 10, height = 8)

plot_enriched_summary <- enriched_summary %>%
  filter(strict_enriched_snp_count > 0 | prevalence80_enriched_snp_count > 0) %>%
  mutate(Clade = fct_reorder(Clade, strict_enriched_snp_count, .desc = TRUE)) %>%
  select(Clade, strict_enriched_snp_count, prevalence80_enriched_snp_count) %>%
  pivot_longer(
    cols = c(strict_enriched_snp_count, prevalence80_enriched_snp_count),
    names_to = "Rule",
    values_to = "Enriched_SNP_count"
  ) %>%
  mutate(
    Rule = recode(
      Rule,
      strict_enriched_snp_count = "Strict: count >=40, outside <1%",
      prevalence80_enriched_snp_count = "Prevalence: inside >=80%, outside <1%"
    )
  ) %>%
  ggplot(aes(x = Clade, y = Enriched_SNP_count, fill = Rule)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Clade-enriched SNP alleles",
    subtitle = "Candidate clade-defining SNP alleles from the Panaroo core SNP alignment",
    x = "Clade",
    y = "Enriched SNP allele count",
    fill = "Rule"
  ) +
  theme_publication(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(plot_enriched_summary, "12E_clade_enriched_snp_summary_by_clade", width = 11, height = 6)

plot_top3_enriched <- top3_enriched_summary %>%
  mutate(Clade = factor(Clade, levels = major_clades)) %>%
  select(Clade, strict_enriched_snp_count, prevalence80_enriched_snp_count) %>%
  pivot_longer(
    cols = c(strict_enriched_snp_count, prevalence80_enriched_snp_count),
    names_to = "Rule",
    values_to = "Enriched_SNP_count"
  ) %>%
  mutate(
    Rule = recode(
      Rule,
      strict_enriched_snp_count = "Strict: count >=40, outside <1%",
      prevalence80_enriched_snp_count = "Prevalence: inside >=80%, outside <1%"
    )
  ) %>%
  ggplot(aes(x = Clade, y = Enriched_SNP_count, fill = Rule)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(
    aes(label = Enriched_SNP_count),
    position = position_dodge(width = 0.8),
    vjust = -0.25,
    size = 3
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Clade-enriched SNP alleles in the top 3 major clades",
    subtitle = "Strict and prevalence-based candidate clade-defining SNP allele rules",
    x = "Major clade",
    y = "Enriched SNP allele count",
    fill = "Rule"
  ) +
  theme_publication(base_size = 11)

save_plot(plot_top3_enriched, "12E_top3_major_clades_enriched_snp_summary", width = 9, height = 5.5)

# ==============================================================================
# 14. Final manifest and terminal report
# ==============================================================================

manifest <- tibble(
  output_type = c(
    "input_alignment",
    "input_clade_table",
    "input_clade_counts",
    "input_metadata",
    "table_dir",
    "matrix_dir",
    "figure_dir",
    "enriched_snp_dir",
    "log_dir"
  ),
  path = c(
    alignment_file,
    clade_file,
    clade_counts_file,
    metadata_used,
    table_dir,
    matrix_dir,
    figure_dir,
    enriched_dir,
    log_dir
  )
)

write_csv(manifest, file.path(out_dir, "12E_output_manifest.csv"))

message("")
message("==========================================================")
message(" SCRIPT 12E COMPLETED SUCCESSFULLY")
message("==========================================================")
message("[OUTPUT] Tables:        ", table_dir)
message("[OUTPUT] Matrices:      ", matrix_dir)
message("[OUTPUT] Figures:       ", figure_dir)
message("[OUTPUT] Enriched SNPs: ", enriched_dir)
message("")
message("=== SNP site summary ===")
print(snp_site_summary)
message("")
message("=== Pairwise quality summary ===")
print(pairwise_quality_summary)
message("")
message("=== Global pairwise SNP summary ===")
print(global_summary)
message("")
message("=== Top 3 major clade within-SNP summary ===")
print(major_within_summary)
message("")
message("=== Top 3 temporal persistence/emergence screen ===")
print(top3_temporal_screen)
message("")
message("=== Top 3 enriched SNP summary ===")
print(top3_enriched_summary)
message("")
message("[DONE] 12E pairwise SNP, temporal dynamics, and enriched SNP analysis finished.")
