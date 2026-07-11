#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 01: ECtyper serotype summary and figures
# ============================================================
#
# This script belongs to:
#   results/01_serotype_summary
#
# It reads Script 00 output from:
#   results/00_ectyper_batch_run/tables/ectyper_output_all.tsv
#
# It writes Script 01 tables to:
#   results/01_serotype_summary/tables
#
# It writes Script 01 figures to:
#   results/01_serotype_summary/figures
#
# Figure behavior:
#   1. Top 10 serotypes are calculated separately per source group.
#   2. Other called serotypes are collapsed as Other serotypes.
#   3. No serotype calls are retained as No serotype call.
#   4. Temporal plots group years <=2000 together.
#   5. Temporal plots display vertical n labels on each year bar.
#   6. Separate O157:H7 versus non-O157:H7 pie charts are produced.
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

wd <- path.expand(Sys.getenv("ECOLI_EPI_WORKDIR", unset = "~/epi/marker_screen"))

metadata_file <- file.path(wd, "data5.csv")

ectyper_file_new <- file.path(
  wd,
  "results",
  "00_ectyper_batch_run",
  "tables",
  "ectyper_output_all.tsv"
)

ectyper_file_old <- file.path(
  wd,
  "ectyper_results",
  "ectyper_output_all.tsv"
)

ectyper_file <- if (file.exists(ectyper_file_new)) {
  ectyper_file_new
} else {
  ectyper_file_old
}

script01_dir <- file.path(wd, "results", "01_serotype_summary")
table_dir <- file.path(script01_dir, "tables")
figure_dir <- file.path(script01_dir, "figures")

fig_top10_dir <- file.path(figure_dir, "top10_serotype_pies")
fig_o157_dir <- file.path(figure_dir, "O157H7_focus_pies")
fig_temporal_count_dir <- file.path(figure_dir, "temporal_stacked_counts")
fig_temporal_percent_dir <- file.path(figure_dir, "temporal_stacked_percent")
fig_multipage_dir <- file.path(figure_dir, "multipage_pdfs")

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_top10_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_o157_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_temporal_count_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_temporal_percent_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_multipage_dir, showWarnings = FALSE, recursive = TRUE)

out_joined <- file.path(table_dir, "current_ectyper_joined_metadata.csv")
out_summary <- file.path(table_dir, "current_serotype_summary_by_source_group.csv")
out_serotype_counts <- file.path(table_dir, "current_serotype_counts_by_source_group.csv")
out_otype_counts <- file.path(table_dir, "current_Otype_counts_by_source_group.csv")
out_htype_counts <- file.path(table_dir, "current_Htype_counts_by_source_group.csv")
out_top10 <- file.path(table_dir, "current_top10_serotypes_by_source_group.csv")
out_plot_groups <- file.path(table_dir, "current_serotype_plot_groups_by_source_group.csv")
out_temporal_counts <- file.path(table_dir, "current_temporal_serotype_counts_by_source_group.csv")
out_temporal_percent <- file.path(table_dir, "current_temporal_serotype_percent_by_source_group.csv")
out_o157_focus <- file.path(table_dir, "current_O157H7_focus_by_source_group.csv")

message("=== Script 01 started: ", Sys.time(), " ===")
message("[INFO] Working directory: ", wd)
message("[INFO] Metadata file: ", metadata_file)
message("[INFO] ECtyper file: ", ectyper_file)
message("[INFO] Script 01 table directory: ", table_dir)
message("[INFO] Script 01 figure directory: ", figure_dir)

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

if (!file.exists(ectyper_file)) {
  stop("ECtyper output file not found: ", ectyper_file)
}

# ============================================================
# 2. Load inputs
# ============================================================

metadata <- read_csv(metadata_file, show_col_types = FALSE)

required_metadata_cols <- c("Assembly", "source_group")
missing_metadata_cols <- setdiff(required_metadata_cols, names(metadata))

if (length(missing_metadata_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_metadata_cols, collapse = ", "))
}

metadata <- metadata %>%
  mutate(
    Assembly = str_squish(as.character(Assembly)),
    source_group = str_squish(as.character(source_group))
  ) %>%
  filter(!is.na(Assembly), Assembly != "") %>%
  distinct(Assembly, .keep_all = TRUE)

ectyper <- read_tsv(ectyper_file, show_col_types = FALSE)

required_ectyper_cols <- c("Name", "O-type", "H-type", "Serotype")
missing_ectyper_cols <- setdiff(required_ectyper_cols, names(ectyper))

if (length(missing_ectyper_cols) > 0) {
  stop("Missing required ECtyper columns: ", paste(missing_ectyper_cols, collapse = ", "))
}

ectyper <- ectyper %>%
  mutate(Name = str_squish(as.character(Name))) %>%
  filter(!is.na(Name), Name != "")

message("[INFO] Metadata assemblies: ", nrow(metadata))
message("[INFO] ECtyper rows: ", nrow(ectyper))

duplicate_names <- ectyper %>%
  count(Name, name = "n") %>%
  filter(n > 1)

if (nrow(duplicate_names) > 0) {
  warning("[WARN] Duplicate ECtyper Name values detected: ", nrow(duplicate_names))
  write_csv(duplicate_names, file.path(table_dir, "current_duplicate_ectyper_names.csv"))
} else {
  message("[INFO] No duplicate ECtyper Name values detected.")
}

# ============================================================
# 3. Helper functions
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

make_temporal_year_group <- function(y) {
  case_when(
    is.na(y) ~ NA_character_,
    y <= 2000 ~ "<=2000",
    TRUE ~ as.character(y)
  )
}

make_temporal_year_factor <- function(x) {
  vals <- unique(as.character(x[!is.na(x)]))
  numeric_years <- suppressWarnings(as.integer(vals[vals != "<=2000"]))
  numeric_years <- sort(numeric_years[!is.na(numeric_years)])
  level_order <- c("<=2000", as.character(numeric_years))
  factor(as.character(x), levels = level_order[level_order %in% vals])
}

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

# ============================================================
# 4. Join and clean
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

if ("Collection_year" %in% names(joined)) {
  joined <- joined %>%
    mutate(
      Collection_year_clean = clean_year(Collection_year),
      Collection_year_temporal_group = make_temporal_year_group(Collection_year_clean)
    )
} else {
  joined <- joined %>%
    mutate(
      Collection_year_clean = NA_integer_,
      Collection_year_temporal_group = NA_character_
    )
}

write_csv(joined, out_joined)

message("[INFO] Joined rows: ", nrow(joined))
message("[INFO] Rows without matched metadata/source group: ", sum(joined$source_group_clean == "Unknown_source_group"))

# ============================================================
# 5. Summary tables
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
# 6. Top 10 serotype grouping
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
# 7. Top 10 serotype pie charts
# ============================================================

make_top10_pie_plot <- function(df_source, source_name) {
  df_plot <- df_source %>%
    arrange(desc(n), Serotype_plot_group) %>%
    mutate(
      Serotype_plot_group = fct_reorder(Serotype_plot_group, n),
      pct = 100 * n / sum(n),
      slice_label = if_else(pct >= 5, paste0(round(pct, 1), "%"), "")
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
      subtitle = "Other called serotypes collapsed as Other serotypes",
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

for (src in source_groups) {
  df_src <- plot_group_counts %>%
    filter(source_group_clean == src)

  p <- make_top10_pie_plot(df_src, src)

  base <- file.path(fig_top10_dir, paste0("pie_top10_serotypes_", safe_filename(src)))

  ggsave(paste0(base, ".pdf"), p, width = 8, height = 6)
  ggsave(paste0(base, ".png"), p, width = 8, height = 6, dpi = 600)
}

pdf(file.path(fig_multipage_dir, "pie_top10_serotypes_all_source_groups.pdf"), width = 8, height = 6)
for (src in source_groups) {
  df_src <- plot_group_counts %>%
    filter(source_group_clean == src)

  print(make_top10_pie_plot(df_src, src))
}
dev.off()

# ============================================================
# 8. O157:H7 focus pie charts
# ============================================================

o157_focus_counts <- joined %>%
  mutate(
    O157H7_focus_group = case_when(
      Serotype_clean == "O157:H7" ~ "O157:H7",
      Serotype_clean == "No_serotype_call" ~ "No serotype call",
      TRUE ~ "Non-O157:H7"
    )
  ) %>%
  count(source_group_clean, O157H7_focus_group, name = "n") %>%
  group_by(source_group_clean) %>%
  mutate(
    source_group_total = sum(n),
    percent_within_source_group = round(100 * n / source_group_total, 2)
  ) %>%
  ungroup() %>%
  arrange(source_group_clean, desc(n), O157H7_focus_group)

write_csv(o157_focus_counts, out_o157_focus)

make_o157_pie_plot <- function(df_source, source_name) {
  df_plot <- df_source %>%
    arrange(desc(n), O157H7_focus_group) %>%
    mutate(
      O157H7_focus_group = fct_reorder(O157H7_focus_group, n),
      pct = 100 * n / sum(n),
      slice_label = paste0(round(pct, 1), "%\n", "n=", n)
    )

  ggplot(df_plot, aes(x = "", y = n, fill = O157H7_focus_group)) +
    geom_col(width = 1, color = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = slice_label),
      position = position_stack(vjust = 0.5),
      size = 3.4,
      color = "black"
    ) +
    labs(
      title = paste0("O157:H7 focus in ", source_name),
      subtitle = "O157:H7, non-O157:H7, and no serotype calls",
      fill = "Serotype category",
      caption = paste0("n = ", sum(df_plot$n), " genomes")
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      plot.caption = element_text(hjust = 0.5, size = 9)
    )
}

for (src in source_groups) {
  df_src <- o157_focus_counts %>%
    filter(source_group_clean == src)

  p <- make_o157_pie_plot(df_src, src)

  base <- file.path(fig_o157_dir, paste0("pie_O157H7_focus_", safe_filename(src)))

  ggsave(paste0(base, ".pdf"), p, width = 7, height = 6)
  ggsave(paste0(base, ".png"), p, width = 7, height = 6, dpi = 600)
}

pdf(file.path(fig_multipage_dir, "pie_O157H7_focus_all_source_groups.pdf"), width = 7, height = 6)
for (src in source_groups) {
  df_src <- o157_focus_counts %>%
    filter(source_group_clean == src)

  print(make_o157_pie_plot(df_src, src))
}
dev.off()

# ============================================================
# 9. Temporal stacked bar plots
# ============================================================

temporal_data <- joined_plot %>%
  filter(!is.na(Collection_year_clean)) %>%
  filter(Collection_year_clean >= 1900, Collection_year_clean <= 2100) %>%
  mutate(
    Collection_year_temporal_group = make_temporal_year_group(Collection_year_clean)
  )

if (nrow(temporal_data) == 0) {
  warning("[WARN] No usable Collection_year values found. Temporal plots skipped.")
} else {

  temporal_counts <- temporal_data %>%
    count(source_group_clean, Collection_year_temporal_group, Serotype_plot_group, name = "n") %>%
    group_by(source_group_clean, Collection_year_temporal_group) %>%
    mutate(
      year_source_total = sum(n),
      percent_within_year_source = round(100 * n / year_source_total, 2),
      prop_within_year_source = n / year_source_total
    ) %>%
    ungroup() %>%
    mutate(
      Collection_year_temporal_group = make_temporal_year_factor(Collection_year_temporal_group)
    ) %>%
    arrange(source_group_clean, Collection_year_temporal_group, desc(n), Serotype_plot_group)

  temporal_percent <- temporal_counts

  write_csv(temporal_counts, out_temporal_counts)
  write_csv(temporal_percent, out_temporal_percent)

  make_temporal_count_plot <- function(df_source, source_name) {
    total_df <- df_source %>%
      distinct(Collection_year_temporal_group, year_source_total)

    ggplot(
      df_source,
      aes(
        x = Collection_year_temporal_group,
        y = n,
        fill = Serotype_plot_group
      )
    ) +
      geom_col(width = 0.85, color = "gray20", linewidth = 0.1) +
      geom_text(
        data = total_df,
        aes(
          x = Collection_year_temporal_group,
          y = year_source_total,
          label = paste0("n=", year_source_total)
        ),
        inherit.aes = FALSE,
        angle = 90,
        vjust = -0.35,
        hjust = 0.5,
        size = 3.0
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      coord_cartesian(clip = "off") +
      labs(
        title = paste0("Temporal distribution of top ECtyper serotypes in ", source_name),
        subtitle = "Years <=2000 grouped; top 10 serotypes shown individually",
        x = "Collection year",
        y = "Number of genomes",
        fill = "Serotype group"
      ) +
      theme_publication(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        plot.margin = margin(t = 20, r = 20, b = 10, l = 10)
      )
  }

  make_temporal_percent_plot <- function(df_source, source_name) {
    total_df <- df_source %>%
      distinct(Collection_year_temporal_group, year_source_total) %>%
      mutate(label_y = 1.02)

    ggplot(
      df_source,
      aes(
        x = Collection_year_temporal_group,
        y = prop_within_year_source,
        fill = Serotype_plot_group
      )
    ) +
      geom_col(width = 0.85, color = "gray20", linewidth = 0.1) +
      geom_text(
        data = total_df,
        aes(
          x = Collection_year_temporal_group,
          y = label_y,
          label = paste0("n=", year_source_total)
        ),
        inherit.aes = FALSE,
        angle = 90,
        vjust = -0.20,
        hjust = 0.5,
        size = 3.0
      ) +
      scale_y_continuous(
        labels = percent_format(accuracy = 1),
        limits = c(0, 1.12),
        expand = expansion(mult = c(0, 0))
      ) +
      coord_cartesian(clip = "off") +
      labs(
        title = paste0("Temporal serotype composition in ", source_name),
        subtitle = "Years <=2000 grouped; percent stacked by collection year",
        x = "Collection year",
        y = "Percent of genomes",
        fill = "Serotype group"
      ) +
      theme_publication(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        plot.margin = margin(t = 20, r = 20, b = 10, l = 10)
      )
  }

  temporal_sources <- sort(unique(temporal_counts$source_group_clean))

  for (src in temporal_sources) {
    df_src <- temporal_percent %>%
      filter(source_group_clean == src)

    p_count <- make_temporal_count_plot(df_src, src)
    p_percent <- make_temporal_percent_plot(df_src, src)

    base_count <- file.path(fig_temporal_count_dir, paste0("temporal_stacked_counts_top10_serotypes_", safe_filename(src)))
    base_percent <- file.path(fig_temporal_percent_dir, paste0("temporal_stacked_percent_top10_serotypes_", safe_filename(src)))

    ggsave(paste0(base_count, ".pdf"), p_count, width = 13, height = 7)
    ggsave(paste0(base_count, ".png"), p_count, width = 13, height = 7, dpi = 600)

    ggsave(paste0(base_percent, ".pdf"), p_percent, width = 13, height = 7)
    ggsave(paste0(base_percent, ".png"), p_percent, width = 13, height = 7, dpi = 600)
  }

  pdf(file.path(fig_multipage_dir, "temporal_stacked_counts_top10_serotypes_all_source_groups.pdf"), width = 13, height = 7)
  for (src in temporal_sources) {
    df_src <- temporal_percent %>%
      filter(source_group_clean == src)

    print(make_temporal_count_plot(df_src, src))
  }
  dev.off()

  pdf(file.path(fig_multipage_dir, "temporal_stacked_percent_top10_serotypes_all_source_groups.pdf"), width = 13, height = 7)
  for (src in temporal_sources) {
    df_src <- temporal_percent %>%
      filter(source_group_clean == src)

    print(make_temporal_percent_plot(df_src, src))
  }
  dev.off()
}

# ============================================================
# 10. Terminal report
# ============================================================

message("")
message("=== Summary by source group ===")
print(summary_by_source, n = Inf)

message("")
message("=== Top 10 serotypes per source group ===")
print(top10_serotypes, n = Inf)

message("")
message("=== O157:H7 focus by source group ===")
print(o157_focus_counts, n = Inf)

message("")
message("[DONE] Tables written to:")
message("  ", table_dir)

message("")
message("[DONE] Figures written to:")
message("  ", figure_dir)

message("=== Script 01 completed: ", Sys.time(), " ===")
