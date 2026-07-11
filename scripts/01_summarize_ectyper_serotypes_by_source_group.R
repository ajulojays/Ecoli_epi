#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 01: Summarize ECtyper serotypes by source group
# ============================================================
#
# Purpose:
#   This script processes the full ECtyper output table, joins it to the
#   curated metadata table, summarizes serotypes by source group, and
#   generates publication-ready figures.
#
# Main inputs:
#
#   1. Metadata:
#      ~/epi/marker_screen/data5.csv
#
#   2. ECtyper combined output:
#      ~/epi/marker_screen/ectyper_results/ectyper_output_all.tsv
#
# Required metadata columns:
#   Assembly
#   source_group
#
# Required ECtyper columns:
#   Name
#   O-type
#   H-type
#   Serotype
#
# Join key:
#   ECtyper Name = metadata Assembly
#
# Main output tables:
#
#   current_ectyper_joined_metadata.csv
#   current_serotype_summary_by_source_group.csv
#   current_serotype_counts_by_source_group.csv
#   current_Otype_counts_by_source_group.csv
#   current_Htype_counts_by_source_group.csv
#   current_top10_serotypes_by_source_group.csv
#   current_serotype_plot_groups_by_source_group.csv
#   current_temporal_serotype_counts_by_source_group.csv
#   current_temporal_serotype_percent_by_source_group.csv
#
# Main figure outputs:
#
#   figures/pie_top10_serotypes_<source_group>.pdf
#   figures/pie_top10_serotypes_<source_group>.png
#   figures/temporal_stacked_counts_top10_serotypes_<source_group>.pdf
#   figures/temporal_stacked_counts_top10_serotypes_<source_group>.png
#   figures/temporal_stacked_percent_top10_serotypes_<source_group>.pdf
#   figures/temporal_stacked_percent_top10_serotypes_<source_group>.png
#
# Grouping logic for figures:
#
#   For each source group independently:
#     - Identify the top 10 called serotypes.
#     - Keep those top 10 serotypes as individual categories.
#     - Collapse all other called serotypes into "Other serotypes".
#     - Keep genomes with no serotype as "No serotype call".
#
# This avoids dropping incomplete calls and prevents long legends with
# hundreds of rare serotypes.
#
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(forcats)
})

# ============================================================
# 1. Paths
# ============================================================
#
# ECOLI_EPI_WORKDIR allows another user to run this script outside Samuel's
# exact directory structure.
#
# Default:
#   ~/epi/marker_screen
#
# Example:
#   ECOLI_EPI_WORKDIR=/path/to/workdir Rscript scripts/01_...
#
# ============================================================

wd <- path.expand(Sys.getenv("ECOLI_EPI_WORKDIR", unset = "~/epi/marker_screen"))

metadata_file <- file.path(wd, "data5.csv")
ectyper_file  <- file.path(wd, "ectyper_results", "ectyper_output_all.tsv")

out_dir <- file.path(wd, "ectyper_results")
fig_dir <- file.path(out_dir, "figures")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

out_joined <- file.path(out_dir, "current_ectyper_joined_metadata.csv")
out_summary <- file.path(out_dir, "current_serotype_summary_by_source_group.csv")
out_serotype_counts <- file.path(out_dir, "current_serotype_counts_by_source_group.csv")
out_otype_counts <- file.path(out_dir, "current_Otype_counts_by_source_group.csv")
out_htype_counts <- file.path(out_dir, "current_Htype_counts_by_source_group.csv")
out_top10 <- file.path(out_dir, "current_top10_serotypes_by_source_group.csv")
out_plot_groups <- file.path(out_dir, "current_serotype_plot_groups_by_source_group.csv")
out_temporal_counts <- file.path(out_dir, "current_temporal_serotype_counts_by_source_group.csv")
out_temporal_percent <- file.path(out_dir, "current_temporal_serotype_percent_by_source_group.csv")

message("=== ECtyper serotype summary started: ", Sys.time(), " ===")
message("[INFO] Working directory: ", wd)
message("[INFO] Metadata file: ", metadata_file)
message("[INFO] ECtyper file: ", ectyper_file)
message("[INFO] Figure directory: ", fig_dir)

# ============================================================
# 2. Safety checks
# ============================================================

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

if (!file.exists(ectyper_file)) {
  stop("ECtyper output file not found: ", ectyper_file)
}

# ============================================================
# 3. Load metadata
# ============================================================

metadata <- read_csv(metadata_file, show_col_types = FALSE)

required_metadata_cols <- c("Assembly", "source_group")
missing_metadata_cols <- setdiff(required_metadata_cols, names(metadata))

if (length(missing_metadata_cols) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_metadata_cols, collapse = ", ")
  )
}

metadata <- metadata %>%
  mutate(
    Assembly = str_squish(as.character(Assembly)),
    source_group = str_squish(as.character(source_group))
  ) %>%
  filter(!is.na(Assembly), Assembly != "") %>%
  distinct(Assembly, .keep_all = TRUE)

message("[INFO] Metadata assemblies loaded: ", nrow(metadata))

# ============================================================
# 4. Load ECtyper output
# ============================================================

ectyper <- read_tsv(ectyper_file, show_col_types = FALSE)

required_ectyper_cols <- c("Name", "O-type", "H-type", "Serotype")
missing_ectyper_cols <- setdiff(required_ectyper_cols, names(ectyper))

if (length(missing_ectyper_cols) > 0) {
  stop(
    "Missing required ECtyper columns: ",
    paste(missing_ectyper_cols, collapse = ", ")
  )
}

ectyper <- ectyper %>%
  mutate(Name = str_squish(as.character(Name))) %>%
  filter(!is.na(Name), Name != "")

message("[INFO] ECtyper rows loaded: ", nrow(ectyper))

# Check duplicate ECtyper rows by genome accession.
# Duplicates should not happen in the completed workflow, but this report is
# useful if retry scripts accidentally append the same genome twice.

duplicate_names <- ectyper %>%
  count(Name, name = "n") %>%
  filter(n > 1)

if (nrow(duplicate_names) > 0) {
  warning("[WARN] Duplicate ECtyper Name values detected: ", nrow(duplicate_names))
  write_csv(duplicate_names, file.path(out_dir, "current_duplicate_ectyper_names.csv"))
} else {
  message("[INFO] No duplicate ECtyper Name values detected.")
}

# ============================================================
# 5. Helper functions
# ============================================================

clean_call <- function(x) {
  x <- str_squish(as.character(x))

  case_when(
    is.na(x) ~ "No_call",
    x == "" ~ "No_call",
    x %in% c("-", "NA", "N/A", "nan", "NaN", "None", "NULL") ~ "No_call",
    TRUE ~ x
  )
}

make_clean_serotype <- function(serotype, otype, htype) {
  serotype_clean <- clean_call(serotype)
  otype_clean <- clean_call(otype)
  htype_clean <- clean_call(htype)

  case_when(
    serotype_clean != "No_call" ~ serotype_clean,
    otype_clean != "No_call" & htype_clean != "No_call" ~ paste0(otype_clean, ":", htype_clean),
    otype_clean != "No_call" & htype_clean == "No_call" ~ paste0(otype_clean, ":No_H_call"),
    otype_clean == "No_call" & htype_clean != "No_call" ~ paste0("No_O_call:", htype_clean),
    TRUE ~ "No_serotype_call"
  )
}

safe_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

clean_year <- function(x) {
  x <- as.character(x)
  y <- str_extract(x, "\\d{4}")
  suppressWarnings(as.integer(y))
}

# ============================================================
# 6. Join ECtyper output to metadata
# ============================================================

joined <- ectyper %>%
  left_join(metadata, by = c("Name" = "Assembly")) %>%
  mutate(
    source_group_clean = case_when(
      is.na(source_group) ~ "Unknown_source_group",
      str_squish(as.character(source_group)) == "" ~ "Unknown_source_group",
      TRUE ~ str_squish(as.character(source_group))
    ),
    O_type_clean = clean_call(`O-type`),
    H_type_clean = clean_call(`H-type`),
    Serotype_clean = make_clean_serotype(
      serotype = Serotype,
      otype = `O-type`,
      htype = `H-type`
    ),
    serotype_call_status = case_when(
      Serotype_clean == "No_serotype_call" ~ "no_O_or_H_call",
      str_detect(Serotype_clean, "No_H_call") ~ "partial_O_only",
      str_detect(Serotype_clean, "No_O_call") ~ "partial_H_only",
      TRUE ~ "complete_or_reconstructed_call"
    )
  )

# Add a cleaned year column if Collection_year exists.
# If Collection_year is absent, temporal plots will be skipped.

if ("Collection_year" %in% names(joined)) {
  joined <- joined %>%
    mutate(Collection_year_clean = clean_year(Collection_year))
} else {
  joined <- joined %>%
    mutate(Collection_year_clean = NA_integer_)
}

message("[INFO] Joined rows: ", nrow(joined))

unmatched_metadata <- joined %>%
  filter(source_group_clean == "Unknown_source_group") %>%
  nrow()

message("[INFO] Rows without matched metadata/source group: ", unmatched_metadata)

write_csv(joined, out_joined)

# ============================================================
# 7. Summary tables
# ============================================================

summary_by_source <- joined %>%
  group_by(source_group_clean) %>%
  summarise(
    total_done = n(),
    serotype_called = sum(Serotype_clean != "No_serotype_call"),
    no_serotype_call = sum(Serotype_clean == "No_serotype_call"),
    partial_O_only = sum(serotype_call_status == "partial_O_only"),
    partial_H_only = sum(serotype_call_status == "partial_H_only"),
    percent_no_serotype_call = round(100 * no_serotype_call / total_done, 2),
    unique_serotypes_called = n_distinct(Serotype_clean[Serotype_clean != "No_serotype_call"]),
    .groups = "drop"
  ) %>%
  arrange(source_group_clean)

write_csv(summary_by_source, out_summary)

serotype_counts_by_source <- joined %>%
  count(source_group_clean, Serotype_clean, name = "n") %>%
  group_by(source_group_clean) %>%
  mutate(
    source_group_total = sum(n),
    percent_within_source_group = round(100 * n / source_group_total, 2)
  ) %>%
  ungroup() %>%
  arrange(source_group_clean, desc(n), Serotype_clean)

write_csv(serotype_counts_by_source, out_serotype_counts)

otype_counts_by_source <- joined %>%
  count(source_group_clean, O_type_clean, name = "n") %>%
  group_by(source_group_clean) %>%
  mutate(
    source_group_total = sum(n),
    percent_within_source_group = round(100 * n / source_group_total, 2)
  ) %>%
  ungroup() %>%
  arrange(source_group_clean, desc(n), O_type_clean)

write_csv(otype_counts_by_source, out_otype_counts)

htype_counts_by_source <- joined %>%
  count(source_group_clean, H_type_clean, name = "n") %>%
  group_by(source_group_clean) %>%
  mutate(
    source_group_total = sum(n),
    percent_within_source_group = round(100 * n / source_group_total, 2)
  ) %>%
  ungroup() %>%
  arrange(source_group_clean, desc(n), H_type_clean)

write_csv(htype_counts_by_source, out_htype_counts)

# ============================================================
# 8. Top 10 serotypes per source group for figures
# ============================================================
#
# Top 10 is calculated independently within each source group.
# "No_serotype_call" is not allowed to consume a top-10 slot.
#
# Plot groups:
#   - top 10 called serotypes remain individual labels
#   - all other called serotypes become "Other serotypes"
#   - no serotype calls become "No serotype call"
#
# ============================================================

top10_serotypes <- joined %>%
  filter(Serotype_clean != "No_serotype_call") %>%
  count(source_group_clean, Serotype_clean, name = "n") %>%
  group_by(source_group_clean) %>%
  arrange(desc(n), Serotype_clean, .by_group = TRUE) %>%
  mutate(rank_within_source_group = row_number()) %>%
  filter(rank_within_source_group <= 10) %>%
  ungroup()

write_csv(top10_serotypes, out_top10)

joined_plot <- joined %>%
  left_join(
    top10_serotypes %>%
      transmute(
        source_group_clean,
        Serotype_clean,
        is_top10_serotype = TRUE
      ),
    by = c("source_group_clean", "Serotype_clean")
  ) %>%
  mutate(
    is_top10_serotype = if_else(is.na(is_top10_serotype), FALSE, is_top10_serotype),
    Serotype_plot_group = case_when(
      Serotype_clean == "No_serotype_call" ~ "No serotype call",
      is_top10_serotype ~ Serotype_clean,
      TRUE ~ "Other serotypes"
    )
  )

plot_group_counts <- joined_plot %>%
  count(source_group_clean, Serotype_plot_group, name = "n") %>%
  group_by(source_group_clean) %>%
  mutate(
    source_group_total = sum(n),
    percent_within_source_group = round(100 * n / source_group_total, 2)
  ) %>%
  ungroup() %>%
  arrange(source_group_clean, desc(n), Serotype_plot_group)

write_csv(plot_group_counts, out_plot_groups)

# ============================================================
# 9. Publication-ready pie charts
# ============================================================

theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size - 1),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(size = base_size - 2, color = "gray30")
    )
}

make_pie_plot <- function(df_source, source_name) {
  df_plot <- df_source %>%
    arrange(desc(n), Serotype_plot_group) %>%
    mutate(
      Serotype_plot_group = fct_reorder(Serotype_plot_group, n),
      pct = 100 * n / sum(n),
      slice_label = if_else(
        pct >= 5,
        paste0(round(pct, 1), "%"),
        ""
      )
    )

  ggplot(df_plot, aes(x = "", y = n, fill = Serotype_plot_group)) +
    geom_col(width = 1, color = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = slice_label),
      position = position_stack(vjust = 0.5),
      size = 3.5,
      color = "black"
    ) +
    labs(
      title = paste0("Top 10 ECtyper serotypes in ", source_name),
      subtitle = "Other called serotypes collapsed as 'Other serotypes'",
      fill = "Serotype group",
      caption = paste0("n = ", sum(df_plot$n), " genomes")
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 9),
      plot.caption = element_text(hjust = 0.5, size = 9)
    )
}

source_groups <- sort(unique(plot_group_counts$source_group_clean))

message("[INFO] Creating pie charts for source groups: ", paste(source_groups, collapse = ", "))

for (src in source_groups) {
  df_src <- plot_group_counts %>%
    filter(source_group_clean == src)

  p <- make_pie_plot(df_src, src)

  base <- file.path(fig_dir, paste0("pie_top10_serotypes_", safe_filename(src)))

  ggsave(paste0(base, ".pdf"), p, width = 8, height = 6, device = cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = 8, height = 6, dpi = 600)
}

# Multi-page PDF containing all source-group pie charts.

pdf(file.path(fig_dir, "pie_top10_serotypes_all_source_groups.pdf"), width = 8, height = 6)
for (src in source_groups) {
  df_src <- plot_group_counts %>%
    filter(source_group_clean == src)

  print(make_pie_plot(df_src, src))
}
dev.off()

# ============================================================
# 10. Temporal stacked bar plots
# ============================================================
#
# Uses Collection_year_clean.
# If year is missing for some genomes, those genomes are excluded from
# temporal plots but retained in all non-temporal summaries.
#
# Two temporal versions are created:
#
#   1. Count stacked bars.
#   2. Percent stacked bars.
#
# ============================================================

temporal_data <- joined_plot %>%
  filter(!is.na(Collection_year_clean)) %>%
  filter(Collection_year_clean >= 1900, Collection_year_clean <= 2100)

if (nrow(temporal_data) == 0) {
  warning("[WARN] No usable Collection_year values found. Temporal plots skipped.")
} else {

  temporal_counts <- temporal_data %>%
    count(source_group_clean, Collection_year_clean, Serotype_plot_group, name = "n") %>%
    group_by(source_group_clean, Collection_year_clean) %>%
    mutate(
      year_source_total = sum(n),
      percent_within_year_source = round(100 * n / year_source_total, 2)
    ) %>%
    ungroup() %>%
    arrange(source_group_clean, Collection_year_clean, desc(n), Serotype_plot_group)

  temporal_percent <- temporal_counts %>%
    mutate(prop_within_year_source = n / year_source_total)

  write_csv(temporal_counts, out_temporal_counts)
  write_csv(temporal_percent, out_temporal_percent)

  make_temporal_count_plot <- function(df_source, source_name) {
    ggplot(
      df_source,
      aes(
        x = factor(Collection_year_clean),
        y = n,
        fill = Serotype_plot_group
      )
    ) +
      geom_col(width = 0.85, color = "gray20", linewidth = 0.1) +
      labs(
        title = paste0("Temporal distribution of top ECtyper serotypes in ", source_name),
        subtitle = "Top 10 serotypes shown individually; remaining called serotypes collapsed",
        x = "Collection year",
        y = "Number of genomes",
        fill = "Serotype group"
      ) +
      theme_publication(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right"
      )
  }

  make_temporal_percent_plot <- function(df_source, source_name) {
    ggplot(
      df_source,
      aes(
        x = factor(Collection_year_clean),
        y = prop_within_year_source,
        fill = Serotype_plot_group
      )
    ) +
      geom_col(width = 0.85, color = "gray20", linewidth = 0.1) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      labs(
        title = paste0("Temporal serotype composition in ", source_name),
        subtitle = "Percent stacked by collection year",
        x = "Collection year",
        y = "Percent of genomes",
        fill = "Serotype group"
      ) +
      theme_publication(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right"
      )
  }

  temporal_sources <- sort(unique(temporal_counts$source_group_clean))

  message("[INFO] Creating temporal stacked bar plots for source groups: ", paste(temporal_sources, collapse = ", "))

  for (src in temporal_sources) {
    df_src <- temporal_percent %>%
      filter(source_group_clean == src)

    p_count <- make_temporal_count_plot(df_src, src)
    p_percent <- make_temporal_percent_plot(df_src, src)

    base_count <- file.path(fig_dir, paste0("temporal_stacked_counts_top10_serotypes_", safe_filename(src)))
    base_percent <- file.path(fig_dir, paste0("temporal_stacked_percent_top10_serotypes_", safe_filename(src)))

    ggsave(paste0(base_count, ".pdf"), p_count, width = 12, height = 7, device = cairo_pdf)
    ggsave(paste0(base_count, ".png"), p_count, width = 12, height = 7, dpi = 600)

    ggsave(paste0(base_percent, ".pdf"), p_percent, width = 12, height = 7, device = cairo_pdf)
    ggsave(paste0(base_percent, ".png"), p_percent, width = 12, height = 7, dpi = 600)
  }

  # Multi-page PDFs.

  pdf(file.path(fig_dir, "temporal_stacked_counts_top10_serotypes_all_source_groups.pdf"), width = 12, height = 7)
  for (src in temporal_sources) {
    df_src <- temporal_percent %>%
      filter(source_group_clean == src)

    print(make_temporal_count_plot(df_src, src))
  }
  dev.off()

  pdf(file.path(fig_dir, "temporal_stacked_percent_top10_serotypes_all_source_groups.pdf"), width = 12, height = 7)
  for (src in temporal_sources) {
    df_src <- temporal_percent %>%
      filter(source_group_clean == src)

    print(make_temporal_percent_plot(df_src, src))
  }
  dev.off()
}

# ============================================================
# 11. Terminal report
# ============================================================

message("")
message("=== Summary by source group ===")
print(summary_by_source, n = Inf)

message("")
message("=== Top 10 serotypes per source group ===")
print(top10_serotypes, n = Inf)

message("")
message("=== Figure grouping counts ===")
print(plot_group_counts, n = Inf)

message("")
message("[DONE] Wrote tables:")
message("  ", out_joined)
message("  ", out_summary)
message("  ", out_serotype_counts)
message("  ", out_otype_counts)
message("  ", out_htype_counts)
message("  ", out_top10)
message("  ", out_plot_groups)
message("  ", out_temporal_counts)
message("  ", out_temporal_percent)

message("")
message("[DONE] Wrote figures to:")
message("  ", fig_dir)

message("=== ECtyper serotype summary completed: ", Sys.time(), " ===")
