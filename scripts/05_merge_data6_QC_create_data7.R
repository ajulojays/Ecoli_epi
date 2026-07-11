#!/usr/bin/env Rscript

# ============================================================
# Script 05
# Purpose:
#   1. Read cattle-focused data6
#   2. Read FASTA assembly QC metrics from Script 04
#   3. Read CheckM2 completeness and contamination results
#   4. Merge QC metrics with metadata
#   5. Create data7 using strict assembly-level QC filters
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

message("=== SCRIPT 05: MERGE QC + CREATE DATA7 STARTED ===")

# =========================
# Paths
# =========================
workdir <- Sys.getenv("ECOLI_EPI_WORKDIR", unset = "~/epi/marker_screen")
workdir <- path.expand(workdir)

data6_file <- Sys.getenv(
  "DATA6_CSV",
  unset = file.path(workdir, "results/02_O157H7_cattle_focus/tables/data6_cattle_O157H7_only.csv")
)

qc04_dir <- Sys.getenv(
  "QC04_DIR",
  unset = file.path(workdir, "results/04_assembly_QC_checkm2")
)

qc05_dir <- Sys.getenv(
  "QC05_DIR",
  unset = file.path(workdir, "results/05_merge_QC_create_data7")
)

out_dir <- file.path(qc05_dir, "tables")
fig_dir <- file.path(qc05_dir, "figures")
log_dir <- file.path(qc05_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

fasta_qc_file <- file.path(qc04_dir, "tables/data6_fasta_assembly_stats.tsv")
checkm2_file  <- file.path(qc04_dir, "checkm2/quality_report.tsv")

# =========================
# QC thresholds
# =========================
max_checkm2_contamination <- as.numeric(Sys.getenv("MAX_CHECKM2_CONTAMINATION", unset = "1.0"))
min_checkm2_completeness  <- as.numeric(Sys.getenv("MIN_CHECKM2_COMPLETENESS", unset = "95.0"))
min_total_bp              <- as.numeric(Sys.getenv("MIN_TOTAL_BP", unset = "4500000"))
max_total_bp              <- as.numeric(Sys.getenv("MAX_TOTAL_BP", unset = "6000000"))

# GC is tracked but not used as hard exclusion by default
min_gc_percent <- Sys.getenv("MIN_GC_PERCENT", unset = "")
max_gc_percent <- Sys.getenv("MAX_GC_PERCENT", unset = "")

min_gc_percent <- ifelse(min_gc_percent == "", NA_real_, as.numeric(min_gc_percent))
max_gc_percent <- ifelse(max_gc_percent == "", NA_real_, as.numeric(max_gc_percent))

# =========================
# Output files
# =========================
data7_all_csv      <- file.path(out_dir, "data7_O157H7_cattle_QC_all.csv")
data7_pass_csv     <- file.path(out_dir, "data7_O157H7_cattle_QC_pass.csv")
data7_fail_csv     <- file.path(out_dir, "data7_O157H7_cattle_QC_fail.csv")
data7_pass_ids_txt <- file.path(out_dir, "data7_O157H7_cattle_QC_pass_assemblies.txt")
summary_txt        <- file.path(out_dir, "data7_assembly_QC_summary.txt")
plot_pdf           <- file.path(fig_dir, "data7_assembly_QC_plots.pdf")

data7_all_rds  <- file.path(out_dir, "data7_O157H7_cattle_QC_all.rds")
data7_pass_rds <- file.path(out_dir, "data7_O157H7_cattle_QC_pass.rds")
data7_fail_rds <- file.path(out_dir, "data7_O157H7_cattle_QC_fail.rds")

# =========================
# Helper functions
# =========================
standardize_names <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_to_lower()
}

normalize_assembly_key <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  x <- str_remove(x, "\\.gz$")
  x <- str_remove(x, "\\.(fasta|fna|fa)$")
  accession_like <- str_extract(x, "GC[AF]_[0-9]+\\.[0-9]+")
  ifelse(is.na(accession_like), x, accession_like)
}

pick_col <- function(df, candidates, label) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    stop(
      paste0(
        "Could not find ", label,
        ". Available columns: ",
        paste(names(df), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  hit[1]
}

range_fail <- function(x, lower = NA_real_, upper = NA_real_) {
  fail <- rep(FALSE, length(x))

  if (!is.na(lower)) {
    fail <- fail | is.na(x) | x < lower
  }

  if (!is.na(upper)) {
    fail <- fail | is.na(x) | x > upper
  }

  fail
}

safe_metric_summary <- function(x, label) {
  tibble(
    metric = c(
      paste0(label, "_min"),
      paste0(label, "_Q1"),
      paste0(label, "_median"),
      paste0(label, "_mean"),
      paste0(label, "_Q3"),
      paste0(label, "_max")
    ),
    value = c(
      suppressWarnings(min(x, na.rm = TRUE)),
      suppressWarnings(as.numeric(quantile(x, 0.25, na.rm = TRUE))),
      suppressWarnings(median(x, na.rm = TRUE)),
      suppressWarnings(mean(x, na.rm = TRUE)),
      suppressWarnings(as.numeric(quantile(x, 0.75, na.rm = TRUE))),
      suppressWarnings(max(x, na.rm = TRUE))
    )
  ) %>%
    mutate(value = ifelse(is.infinite(value), NA_real_, value))
}

# =========================
# Input checks
# =========================
if (!file.exists(data6_file)) {
  stop("[ERROR] data6 file not found: ", data6_file, call. = FALSE)
}

if (!file.exists(fasta_qc_file)) {
  stop("[ERROR] FASTA QC table not found. Run Script 04 first: ", fasta_qc_file, call. = FALSE)
}

if (!file.exists(checkm2_file)) {
  stop("[ERROR] CheckM2 report not found. Run Script 04 first: ", checkm2_file, call. = FALSE)
}

# =========================
# Load data
# =========================
message("[INFO] Reading data6: ", data6_file)
data6 <- read_csv(data6_file, show_col_types = FALSE)

message("[INFO] Reading FASTA QC: ", fasta_qc_file)
fasta_qc <- read_tsv(fasta_qc_file, show_col_types = FALSE)

message("[INFO] Reading CheckM2 report: ", checkm2_file)
checkm2_raw <- read_tsv(checkm2_file, show_col_types = FALSE)

assembly_col <- pick_col(
  data6,
  c("Assembly", "assembly", "assembly_accession", "accession"),
  "assembly accession column in data6"
)

data6 <- data6 %>%
  mutate(assembly_key = normalize_assembly_key(.data[[assembly_col]]))

names(fasta_qc) <- standardize_names(names(fasta_qc))

fasta_key_col <- pick_col(
  fasta_qc,
  c("assembly_key", "assembly", "name", "genome"),
  "assembly key column in FASTA QC table"
)

fasta_qc <- fasta_qc %>%
  mutate(assembly_key = normalize_assembly_key(.data[[fasta_key_col]])) %>%
  mutate(
    total_bp = as.numeric(total_bp),
    gc_percent = as.numeric(gc_percent),
    n_contigs = as.numeric(n_contigs),
    contig_n50 = as.numeric(contig_n50),
    longest_contig = as.numeric(longest_contig),
    n_bases = as.numeric(n_bases)
  ) %>%
  select(
    assembly_key,
    original_fasta,
    total_bp,
    gc_percent,
    n_contigs,
    contig_n50,
    longest_contig,
    n_bases
  ) %>%
  distinct(assembly_key, .keep_all = TRUE)

names(checkm2_raw) <- standardize_names(names(checkm2_raw))

checkm2_id_col <- pick_col(
  checkm2_raw,
  c("name", "assembly_key", "genome", "bin_id", "sample"),
  "CheckM2 ID column"
)

checkm2_completeness_col <- pick_col(
  checkm2_raw,
  c("completeness", "checkm2_completeness"),
  "CheckM2 completeness column"
)

checkm2_contamination_col <- pick_col(
  checkm2_raw,
  c("contamination", "checkm2_contamination"),
  "CheckM2 contamination column"
)

checkm2_qc <- checkm2_raw %>%
  mutate(
    assembly_key = normalize_assembly_key(.data[[checkm2_id_col]]),
    checkm2_completeness = as.numeric(.data[[checkm2_completeness_col]]),
    checkm2_contamination = as.numeric(.data[[checkm2_contamination_col]])
  ) %>%
  select(
    assembly_key,
    checkm2_completeness,
    checkm2_contamination
  ) %>%
  distinct(assembly_key, .keep_all = TRUE)

# =========================
# Merge QC metrics
# =========================
data7_all <- data6 %>%
  left_join(fasta_qc, by = "assembly_key") %>%
  left_join(checkm2_qc, by = "assembly_key")

# =========================
# Apply QC gates
# =========================
data7_all <- data7_all %>%
  mutate(
    qc_missing_fasta = is.na(total_bp),
    qc_missing_checkm2 = is.na(checkm2_contamination) | is.na(checkm2_completeness),

    qc_checkm2_contamination_fail = is.na(checkm2_contamination) |
      checkm2_contamination >= max_checkm2_contamination,

    qc_checkm2_completeness_fail = is.na(checkm2_completeness) |
      checkm2_completeness < min_checkm2_completeness,

    qc_size_fail = range_fail(total_bp, lower = min_total_bp, upper = max_total_bp),
    qc_gc_fail = range_fail(gc_percent, lower = min_gc_percent, upper = max_gc_percent),

    qc_exclusion_reason = case_when(
      qc_missing_fasta ~ "missing_fasta_qc",
      qc_missing_checkm2 ~ "missing_checkm2",
      qc_checkm2_contamination_fail ~ "checkm2_contamination_ge_1",
      qc_checkm2_completeness_fail ~ "checkm2_completeness_lt_95",
      qc_size_fail ~ "genome_size_outside_4.5_6.0Mb",
      qc_gc_fail ~ "gc_outlier",
      TRUE ~ "pass"
    ),

    qc_pass = qc_exclusion_reason == "pass"
  )

data7_pass <- data7_all %>% filter(qc_pass)
data7_fail <- data7_all %>% filter(!qc_pass)

# =========================
# Export tables
# =========================
write_csv(data7_all, data7_all_csv)
write_csv(data7_pass, data7_pass_csv)
write_csv(data7_fail, data7_fail_csv)

write_lines(data7_pass$assembly_key, data7_pass_ids_txt)

saveRDS(data7_all, data7_all_rds)
saveRDS(data7_pass, data7_pass_rds)
saveRDS(data7_fail, data7_fail_rds)

# =========================
# Summary report
# =========================
exclusion_counts <- data7_all %>%
  count(qc_exclusion_reason, sort = TRUE)

summary_metrics <- bind_rows(
  tibble(metric = "n_data6_input", value = nrow(data6)),
  tibble(metric = "n_merged", value = nrow(data7_all)),
  tibble(metric = "n_data7_pass", value = nrow(data7_pass)),
  tibble(metric = "n_data7_fail", value = nrow(data7_fail)),
  tibble(metric = "max_checkm2_contamination_allowed", value = max_checkm2_contamination),
  tibble(metric = "min_checkm2_completeness_required", value = min_checkm2_completeness),
  tibble(metric = "min_total_bp_allowed", value = min_total_bp),
  tibble(metric = "max_total_bp_allowed", value = max_total_bp),
  safe_metric_summary(data7_all$total_bp, "total_bp"),
  safe_metric_summary(data7_all$gc_percent, "gc_percent"),
  safe_metric_summary(data7_all$n_contigs, "fasta_n_contigs"),
  safe_metric_summary(data7_all$contig_n50, "fasta_contig_n50"),
  safe_metric_summary(data7_all$checkm2_completeness, "checkm2_completeness"),
  safe_metric_summary(data7_all$checkm2_contamination, "checkm2_contamination")
)

write_csv(summary_metrics, file.path(out_dir, "data7_assembly_QC_summary_metrics.csv"))

sink(summary_txt)
cat("DATA7 ASSEMBLY QC SUMMARY\n")
cat("=========================\n\n")
cat("Input data6 rows:       ", nrow(data6), "\n")
cat("Merged rows:            ", nrow(data7_all), "\n")
cat("QC pass rows:           ", nrow(data7_pass), "\n")
cat("QC fail rows:           ", nrow(data7_fail), "\n\n")

cat("QC thresholds\n")
cat("-------------\n")
cat("CheckM2 contamination:  <", max_checkm2_contamination, "%\n")
cat("CheckM2 completeness:   >=", min_checkm2_completeness, "%\n")
cat("Genome size:            ", min_total_bp, "to", max_total_bp, "bp\n")
cat("GC hard filter:         ", ifelse(is.na(min_gc_percent), "none", min_gc_percent), "to", ifelse(is.na(max_gc_percent), "none", max_gc_percent), "\n\n")

cat("Exclusion counts\n")
cat("----------------\n")
print(exclusion_counts)

cat("\nSummary metrics\n")
cat("---------------\n")
print(summary_metrics)
sink()

# =========================
# Plots
# =========================
pdf(plot_pdf, width = 8, height = 6)

print(
  ggplot(data7_all, aes(x = checkm2_contamination)) +
    geom_histogram(bins = 50, color = "white") +
    geom_vline(xintercept = max_checkm2_contamination, linetype = "dashed", linewidth = 1) +
    theme_bw() +
    labs(
      title = "CheckM2 contamination distribution",
      x = "Contamination (%)",
      y = "Number of assemblies"
    )
)

print(
  ggplot(data7_all, aes(x = checkm2_completeness)) +
    geom_histogram(bins = 50, color = "white") +
    geom_vline(xintercept = min_checkm2_completeness, linetype = "dashed", linewidth = 1) +
    theme_bw() +
    labs(
      title = "CheckM2 completeness distribution",
      x = "Completeness (%)",
      y = "Number of assemblies"
    )
)

print(
  ggplot(data7_all, aes(x = total_bp / 1e6)) +
    geom_histogram(bins = 50, color = "white") +
    geom_vline(xintercept = c(min_total_bp / 1e6, max_total_bp / 1e6), linetype = "dashed", linewidth = 1) +
    theme_bw() +
    labs(
      title = "FASTA-derived genome size distribution",
      x = "Total assembly size (Mb)",
      y = "Number of assemblies"
    )
)

print(
  ggplot(data7_all, aes(x = gc_percent)) +
    geom_histogram(bins = 50, color = "white") +
    theme_bw() +
    labs(
      title = "FASTA-derived GC content distribution",
      x = "GC (%)",
      y = "Number of assemblies"
    )
)

print(
  ggplot(data7_all, aes(x = qc_exclusion_reason)) +
    geom_bar() +
    coord_flip() +
    theme_bw() +
    labs(
      title = "Assembly QC outcome",
      x = "QC outcome",
      y = "Number of assemblies"
    )
)

print(
  ggplot(data7_all, aes(x = checkm2_contamination, y = checkm2_completeness)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = max_checkm2_contamination, linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = min_checkm2_completeness, linetype = "dashed", linewidth = 1) +
    theme_bw() +
    labs(
      title = "CheckM2 completeness versus contamination",
      x = "Contamination (%)",
      y = "Completeness (%)"
    )
)

dev.off()

message("[DONE] Script 05 complete.")
message("[OUTPUT] data7 all:  ", data7_all_csv)
message("[OUTPUT] data7 pass: ", data7_pass_csv)
message("[OUTPUT] data7 fail: ", data7_fail_csv)
message("[OUTPUT] summary:    ", summary_txt)
message("[OUTPUT] plots:      ", plot_pdf)
