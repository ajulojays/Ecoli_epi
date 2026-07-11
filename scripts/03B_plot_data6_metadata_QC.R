#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 03B: Plot metadata QC distributions for data6
# ============================================================
#
# Purpose:
#   Plot N50 and contig-number distributions for genomes in data6.
#
# Default input priority:
#   1. DATA6_CSV environment variable if supplied
#   2. /home/samuelajulo/epi/data6.csv
#   3. ~/epi/marker_screen/results/02_O157H7_cattle_focus/tables/data6_cattle_O157H7_only.csv
#   4. ~/epi/marker_screen/results/02_O157H7_cattle_focus/tables/data6_cattle_with_ectyper_metadata.csv
#
# Output:
#   ~/epi/marker_screen/results/03_O157H7_cattle_QC/tables/data6_metadata_QC_summary.csv
#   ~/epi/marker_screen/results/03_O157H7_cattle_QC/figures/metadata_QC/
#
# Figures:
#   data6_N50_distribution.pdf/png
#   data6_contig_count_distribution.pdf/png
#   data6_N50_vs_contigs_scatter.pdf/png
#   data6_metadata_QC_multipanel.pdf/png
#
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(tidyr)
})

wd <- path.expand(Sys.getenv("ECOLI_EPI_WORKDIR", unset = "~/epi/marker_screen"))

candidate_inputs <- c(
  Sys.getenv("DATA6_CSV", unset = ""),
  "/home/samuelajulo/epi/data6.csv",
  file.path(wd, "results", "02_O157H7_cattle_focus", "tables", "data6_cattle_O157H7_only.csv"),
  file.path(wd, "results", "02_O157H7_cattle_focus", "tables", "data6_cattle_with_ectyper_metadata.csv")
)

candidate_inputs <- candidate_inputs[candidate_inputs != ""]

data6_file <- candidate_inputs[file.exists(candidate_inputs)][1]

if (is.na(data6_file) || length(data6_file) == 0) {
  stop(
    "No data6 input file found. Tried:\n",
    paste(candidate_inputs, collapse = "\n")
  )
}

base_dir <- file.path(wd, "results", "03_O157H7_cattle_QC")
table_dir <- file.path(base_dir, "tables")
figure_dir <- file.path(base_dir, "figures", "metadata_QC")

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

message("=== Script 03B started: ", Sys.time(), " ===")
message("[INFO] data6 input: ", data6_file)
message("[INFO] table output: ", table_dir)
message("[INFO] figure output: ", figure_dir)

data6_raw <- read_csv(data6_file, show_col_types = FALSE)

names_original <- names(data6_raw)
names_clean <- names_original %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("[^A-Za-z0-9_]", "_") %>%
  str_replace_all("_+", "_") %>%
  str_replace_all("^_|_$", "")

names(data6_raw) <- names_clean

find_col <- function(possible_names, available_names) {
  possible_clean <- possible_names %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[^A-Za-z0-9_]", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")

  hit <- possible_clean[possible_clean %in% available_names]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}

assembly_col <- find_col(
  c("Assembly", "assembly", "Assembly_Accession", "accession"),
  names(data6_raw)
)

n50_col <- find_col(
  c("N50", "n50", "Assembly_N50", "assembly_n50"),
  names(data6_raw)
)

contigs_col <- find_col(
  c("Contigs", "contigs", "Number_of_contigs", "number_of_contigs", "contig_count", "Contig_count"),
  names(data6_raw)
)

if (is.na(n50_col)) {
  stop(
    "Could not find N50 column. Available columns:\n",
    paste(names(data6_raw), collapse = ", ")
  )
}

if (is.na(contigs_col)) {
  stop(
    "Could not find Contigs column. Available columns:\n",
    paste(names(data6_raw), collapse = ", ")
  )
}

message("[INFO] Assembly column: ", assembly_col)
message("[INFO] N50 column: ", n50_col)
message("[INFO] Contigs column: ", contigs_col)

data6 <- data6_raw %>%
  mutate(
    Assembly_clean = if (!is.na(assembly_col)) as.character(.data[[assembly_col]]) else as.character(row_number()),
    N50_numeric = suppressWarnings(as.numeric(str_replace_all(as.character(.data[[n50_col]]), ",", ""))),
    Contigs_numeric = suppressWarnings(as.numeric(str_replace_all(as.character(.data[[contigs_col]]), ",", ""))),
    N50_Mb = N50_numeric / 1000000
  )

qc_summary <- tibble(
  metric = c(
    "n_genomes",
    "n_missing_N50",
    "n_missing_contigs",
    "N50_min",
    "N50_Q1",
    "N50_median",
    "N50_mean",
    "N50_Q3",
    "N50_max",
    "contigs_min",
    "contigs_Q1",
    "contigs_median",
    "contigs_mean",
    "contigs_Q3",
    "contigs_max",
    "n_N50_below_100kb",
    "n_contigs_above_300"
  ),
  value = c(
    nrow(data6),
    sum(is.na(data6$N50_numeric)),
    sum(is.na(data6$Contigs_numeric)),
    min(data6$N50_numeric, na.rm = TRUE),
    quantile(data6$N50_numeric, 0.25, na.rm = TRUE),
    median(data6$N50_numeric, na.rm = TRUE),
    mean(data6$N50_numeric, na.rm = TRUE),
    quantile(data6$N50_numeric, 0.75, na.rm = TRUE),
    max(data6$N50_numeric, na.rm = TRUE),
    min(data6$Contigs_numeric, na.rm = TRUE),
    quantile(data6$Contigs_numeric, 0.25, na.rm = TRUE),
    median(data6$Contigs_numeric, na.rm = TRUE),
    mean(data6$Contigs_numeric, na.rm = TRUE),
    quantile(data6$Contigs_numeric, 0.75, na.rm = TRUE),
    max(data6$Contigs_numeric, na.rm = TRUE),
    sum(data6$N50_numeric < 100000, na.rm = TRUE),
    sum(data6$Contigs_numeric > 300, na.rm = TRUE)
  )
)

write_csv(qc_summary, file.path(table_dir, "data6_metadata_QC_summary.csv"))

data6_qc_table <- data6 %>%
  mutate(
    metadata_qc_flag = case_when(
      is.na(N50_numeric) | is.na(Contigs_numeric) ~ "missing_metadata",
      N50_numeric < 100000 ~ "N50_below_100kb",
      Contigs_numeric > 300 ~ "contigs_above_300",
      TRUE ~ "metadata_pass"
    )
  )

write_csv(data6_qc_table, file.path(table_dir, "data6_metadata_QC_per_genome.csv"))

theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size - 1),
      plot.caption = element_text(size = base_size - 2, color = "gray30")
    )
}

n_genomes <- nrow(data6_qc_table)

p_n50 <- ggplot(data6_qc_table, aes(x = N50_Mb)) +
  geom_histogram(bins = 30, color = "gray20", linewidth = 0.2) +
  geom_vline(xintercept = 0.1, linetype = "dashed", linewidth = 0.6) +
  labs(
    title = "N50 distribution of data6 genomes",
    subtitle = "Dashed line marks N50 = 100 kb",
    x = "N50 (Mb)",
    y = "Number of genomes",
    caption = paste0("n = ", n_genomes, " genomes")
  ) +
  theme_publication(base_size = 12)

ggsave(file.path(figure_dir, "data6_N50_distribution.pdf"), p_n50, width = 8, height = 6)
ggsave(file.path(figure_dir, "data6_N50_distribution.png"), p_n50, width = 8, height = 6, dpi = 600)

p_contigs <- ggplot(data6_qc_table, aes(x = Contigs_numeric)) +
  geom_histogram(bins = 30, color = "gray20", linewidth = 0.2) +
  geom_vline(xintercept = 300, linetype = "dashed", linewidth = 0.6) +
  labs(
    title = "Contig count distribution of data6 genomes",
    subtitle = "Dashed line marks 300 contigs",
    x = "Number of contigs",
    y = "Number of genomes",
    caption = paste0("n = ", n_genomes, " genomes")
  ) +
  theme_publication(base_size = 12)

ggsave(file.path(figure_dir, "data6_contig_count_distribution.pdf"), p_contigs, width = 8, height = 6)
ggsave(file.path(figure_dir, "data6_contig_count_distribution.png"), p_contigs, width = 8, height = 6, dpi = 600)

p_scatter <- ggplot(data6_qc_table, aes(x = Contigs_numeric, y = N50_Mb)) +
  geom_point(alpha = 0.75, size = 2.1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = 300, linetype = "dashed", linewidth = 0.6) +
  labs(
    title = "N50 versus contig count in data6 genomes",
    subtitle = "Dashed lines mark N50 = 100 kb and contigs = 300",
    x = "Number of contigs",
    y = "N50 (Mb)",
    caption = paste0("n = ", n_genomes, " genomes")
  ) +
  theme_publication(base_size = 12)

ggsave(file.path(figure_dir, "data6_N50_vs_contigs_scatter.pdf"), p_scatter, width = 8, height = 6)
ggsave(file.path(figure_dir, "data6_N50_vs_contigs_scatter.png"), p_scatter, width = 8, height = 6, dpi = 600)

pdf(file.path(figure_dir, "data6_metadata_QC_multipanel.pdf"), width = 8, height = 6)
print(p_n50)
print(p_contigs)
print(p_scatter)
dev.off()

png(file.path(figure_dir, "data6_metadata_QC_multipanel.png"), width = 8, height = 6, units = "in", res = 600)
print(p_n50)
print(p_contigs)
print(p_scatter)
dev.off()

message("")
message("=== data6 metadata QC summary ===")
print(qc_summary, n = Inf)

message("")
message("[DONE] Tables written to:")
message("  ", table_dir)

message("")
message("[DONE] Figures written to:")
message("  ", figure_dir)

message("=== Script 03B completed: ", Sys.time(), " ===")
