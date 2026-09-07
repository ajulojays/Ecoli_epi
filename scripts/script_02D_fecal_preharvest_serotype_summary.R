#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 02D: Fecal Pre-Harvest ECtyper Serotype Summary & Figures
# ============================================================
#
# Target Cohorts:
#   1. Beef - Pre-harvest
#   2. Unspecified / Cattle-unspecified - Pre-harvest 
#      (strictly feces/manure matrices, excluding dairy or cow annotations)
#
# Output Directory:
#   results/02D_fecal_preharvest_serotype_summary
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
# 1. Paths and Setup
# ============================================================

BASE_DIR <- "/home/samuelajulo/epi/Ecoli_epi"
METADATA_CSV <- "/home/samuelajulo/epi/data6a.csv"

ectyper_file_new <- file.path(
  BASE_DIR,
  "results",
  "00_ectyper_batch_run",
  "tables",
  "ectyper_output_all.tsv"
)

ectyper_file_old <- file.path(
  "/home/samuelajulo/epi/marker_screen",
  "ectyper_results",
  "ectyper_output_all.tsv"
)

ectyper_file <- if (file.exists(ectyper_file_new)) {
  ectyper_file_new
} else {
  ectyper_file_old
}

script02d_dir <- file.path(BASE_DIR, "results", "02D_fecal_preharvest_serotype_summary")
table_dir <- file.path(script02d_dir, "tables")
figure_dir <- file.path(script02d_dir, "figures")

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

message("=== Script 02D started: ", Sys.time(), " ===")
message("[INFO] Base directory: ", BASE_DIR)
message("[INFO] Metadata file: ", METADATA_CSV)
message("[INFO] ECtyper file: ", ectyper_file)
message("[INFO] Output directory: ", script02d_dir)

if (!file.exists(METADATA_CSV)) stop("Metadata file not found: ", METADATA_CSV)
if (!file.exists(ectyper_file)) stop("ECtyper output file not found: ", ectyper_file)

# ============================================================
# 2. Load Inputs and Apply Targeted Filtering
# ============================================================

metadata <- read_csv(METADATA_CSV, show_col_types = FALSE)
ectyper <- read_tsv(ectyper_file, show_col_types = FALSE)

metadata <- metadata %>%
  mutate(
    Assembly = str_squish(as.character(Assembly)),
    Production_system = str_squish(as.character(Production_system)),
    Supply_chain_stage = str_squish(as.character(Supply_chain_stage))
  ) %>%
  filter(!is.na(Assembly), Assembly != "") %>%
  distinct(Assembly, .keep_all = TRUE)

ectyper <- ectyper %>%
  mutate(Name = str_squish(as.character(Name))) %>%
  filter(!is.na(Name), Name != "")

# Identify available text columns for isolation source / matrix checking
text_cols <- names(metadata)[sapply(metadata, is.character)]

# Combine text content per row for robust keyword searching
metadata <- metadata %>%
  rowwise() %>%
  mutate(combined_text = tolower(paste(c_across(all_of(text_cols)), collapse = " "))) %>%
  ungroup()

# Apply targeted filtering criteria
filtered_metadata <- metadata %>%
  filter(
    (str_detect(tolower(Production_system), "beef") & str_detect(tolower(Supply_chain_stage), "pre")) |
    ((str_detect(tolower(Production_system), "unspecified") | str_detect(tolower(Production_system), "cattle-unspecified")) & 
     str_detect(tolower(Supply_chain_stage), "pre") &
     (str_detect(combined_text, "fec") | str_detect(combined_text, "manure")) &
     !str_detect(combined_text, "dairy") & 
     !str_detect(combined_text, "cow"))
  ) %>%
  mutate(
    cohort_group = case_when(
      str_detect(tolower(Production_system), "beef") ~ "Beef_Preharvest",
      TRUE ~ "Unspecified_Fecal_Preharvest"
    )
  )

message("[INFO] Total assemblies in metadata: ", nrow(metadata))
message("[INFO] Assemblies meeting targeted pre-harvest filter: ", nrow(filtered_metadata))
print(table(filtered_metadata$cohort_group, useNA = "ifany"))

# ============================================================
# 3. Helper Functions
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
# 4. Join and Clean
# ============================================================

joined <- ectyper %>%
  inner_join(filtered_metadata, by = c("Name" = "Assembly")) %>%
  mutate(
    source_group_clean = cohort_group,
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

year_col <- names(joined)[str_detect(tolower(names(joined)), "year|date")][1]
if (!is.na(year_col)) {
  joined <- joined %>%
    mutate(
      Collection_year_clean = clean_year(.data[[year_col]]),
      Collection_year_temporal_group = make_temporal_year_group(Collection_year_clean)
    )
} else {
  joined <- joined %>%
    mutate(
      Collection_year_clean = NA_integer_,
      Collection_year_temporal_group = NA_character_
    )
}

write_csv(joined, file.path(table_dir, "fecal_preharvest_joined_metadata.csv"))
message("[INFO] Joined filtered rows: ", nrow(joined))

# ============================================================
# 5. Summary Tables & Top 10 Serotypes
# ============================================================

summary_by_source <- joined %>%
  group_by(source_group_clean) %>%
  summarise(
    total_done = n(),
    serotype_called = sum(Serotype_clean != "No_serotype_call"),
    no_serotype_call = sum(Serotype_clean == "No_serotype_call"),
    percent_no_serotype_call = round(100 * no_serotype_call / total_done, 2),
    unique_serotypes_called = n_distinct(Serotype_clean[Serotype_clean != "No_serotype_call"]),
    .groups = "drop"
  )

write_csv(summary_by_source, file.path(table_dir, "fecal_preharvest_summary_by_cohort.csv"))

top10_serotypes <- joined %>%
  filter(Serotype_clean != "No_serotype_call") %>%
  count(source_group_clean, Serotype_clean, name = "n") %>%
  group_by(source_group_clean) %>%
  arrange(desc(n), Serotype_clean, .by_group = TRUE) %>%
  mutate(rank_within_source_group = row_number()) %>%
  filter(rank_within_source_group <= 10) %>%
  ungroup()

write_csv(top10_serotypes, file.path(table_dir, "fecal_preharvest_top10_serotypes.csv"))

joined_plot <- joined %>%
  left_join(
    top10_serotypes %>% transmute(source_group_clean, Serotype_clean, is_top10 = TRUE),
    by = c("source_group_clean", "Serotype_clean")
  ) %>%
  mutate(
    is_top10 = if_else(is.na(is_top10), FALSE, is_top10),
    Serotype_plot_group = case_when(
      Serotype_clean == "No_serotype_call" ~ "No serotype call",
      is_top10 ~ Serotype_clean,
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
  ungroup()

write_csv(plot_group_counts, file.path(table_dir, "fecal_preharvest_plot_groups.csv"))

# ============================================================
# 6. Top 10 Serotype Pie Charts
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
    geom_text(aes(label = slice_label), position = position_stack(vjust = 0.5), size = 3.5, color = "black") +
    labs(
      title = paste0("Top 10 Serotypes: ", source_name),
      subtitle = "Filtered Fecal Pre-Harvest Cohort",
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
  df_src <- plot_group_counts %>% filter(source_group_clean == src)
  p <- make_top10_pie_plot(df_src, src)
  base <- file.path(fig_top10_dir, paste0("pie_top10_", safe_filename(src)))
  ggsave(paste0(base, ".pdf"), p, width = 8, height = 6)
  ggsave(paste0(base, ".png"), p, width = 8, height = 6, dpi = 600)
}

# ============================================================
# 7. O157:H7 Focus Pie Charts
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
  ungroup()

write_csv(o157_focus_counts, file.path(table_dir, "fecal_preharvest_O157H7_focus.csv"))

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
    geom_text(aes(label = slice_label), position = position_stack(vjust = 0.5), size = 3.4, color = "black") +
    labs(
      title = paste0("O157:H7 Focus: ", source_name),
      subtitle = "O157:H7 vs Non-O157:H7 vs No Call",
      fill = "Category",
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
  df_src <- o157_focus_counts %>% filter(source_group_clean == src)
  p <- make_o157_pie_plot(df_src, src)
  base <- file.path(fig_o157_dir, paste0("pie_O157H7_focus_", safe_filename(src)))
  ggsave(paste0(base, ".pdf"), p, width = 7, height = 6)
  ggsave(paste0(base, ".png"), p, width = 7, height = 6, dpi = 600)
}

# ============================================================
# 8. Temporal Stacked Bar Plots (<=2000 Grouping)
# ============================================================

temporal_data <- joined_plot %>%
  filter(!is.na(Collection_year_clean), Collection_year_clean >= 1900, Collection_year_clean <= 2100) %>%
  mutate(Collection_year_temporal_group = make_temporal_year_group(Collection_year_clean))

if (nrow(temporal_data) > 0) {
  temporal_counts <- temporal_data %>%
    count(source_group_clean, Collection_year_temporal_group, Serotype_plot_group, name = "n") %>%
    group_by(source_group_clean, Collection_year_temporal_group) %>%
    mutate(
      year_source_total = sum(n),
      prop_within_year_source = n / year_source_total
    ) %>%
    ungroup() %>%
    mutate(Collection_year_temporal_group = make_temporal_year_factor(Collection_year_temporal_group))

  write_csv(temporal_counts, file.path(table_dir, "fecal_preharvest_temporal_counts.csv"))

  make_temporal_count_plot <- function(df_source, source_name) {
    total_df <- df_source %>% distinct(Collection_year_temporal_group, year_source_total)

    ggplot(df_source, aes(x = Collection_year_temporal_group, y = n, fill = Serotype_plot_group)) +
      geom_col(width = 0.85, color = "gray20", linewidth = 0.1) +
      geom_text(
        data = total_df,
        aes(x = Collection_year_temporal_group, y = year_source_total, label = paste0("n=", year_source_total)),
        inherit.aes = FALSE, angle = 90, vjust = -0.35, hjust = 0.5, size = 3.0
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      coord_cartesian(clip = "off") +
      labs(
        title = paste0("Temporal Distribution: ", source_name),
        subtitle = "Years <=2000 grouped; Top 10 serotypes shown",
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

  for (src in unique(temporal_counts$source_group_clean)) {
    df_src <- temporal_counts %>% filter(source_group_clean == src)
    p_count <- make_temporal_count_plot(df_src, src)
    base_count <- file.path(fig_temporal_count_dir, paste0("temporal_counts_", safe_filename(src)))
    ggsave(paste0(base_count, ".pdf"), p_count, width = 13, height = 7)
    ggsave(paste0(base_count, ".png"), p_count, width = 13, height = 7, dpi = 600)
  }
}

message("")
message("=== Summary by Cohort Group ===")
print(summary_by_source)

message("")
message("[DONE] Script 02D completed successfully.")
message("Tables saved to: ", table_dir)
message("Figures saved to: ", figure_dir)
message("=== Script 02D completed: ", Sys.time(), " ===")
