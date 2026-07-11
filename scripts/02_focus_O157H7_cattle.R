#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 02: Cattle O157:H7-focused dataset and temporal analysis
# ============================================================
#
# Purpose:
#   Create data6, a cattle-focused dataset derived from data5 joined with
#   all available ECtyper metadata, then summarize and plot the temporal
#   proportion of O157:H7 compared with all other cattle serotypes.
#
# Input:
#   Metadata:
#     ~/epi/marker_screen/data5.csv
#
#   ECtyper output from Script 00:
#     ~/epi/marker_screen/results/00_ectyper_batch_run/tables/ectyper_output_all.tsv
#
# Output folder:
#   ~/epi/marker_screen/results/02_O157H7_cattle_focus/
#
# Main outputs:
#   tables/data6_cattle_with_ectyper_metadata.csv
#   tables/data6_cattle_O157H7_only.csv
#   tables/data6_cattle_serotype_counts.csv
#   tables/data6_cattle_temporal_O157H7_summary.csv
#   tables/data6_cattle_temporal_O157H7_long_counts.csv
#
# Figures:
#   figures/cattle_O157H7_temporal_proportion.pdf
#   figures/cattle_O157H7_temporal_proportion.png
#   figures/cattle_O157H7_temporal_proportion_line.pdf
#   figures/cattle_O157H7_temporal_proportion_line.png
#   figures/cattle_O157H7_temporal_counts.pdf
#   figures/cattle_O157H7_temporal_counts.png
#   figures/cattle_O157H7_temporal_percent_stacked.pdf
#   figures/cattle_O157H7_temporal_percent_stacked.png
#
# Temporal rule:
#   Years <=2000 are grouped together as "<=2000".
#   Years after 2000 are shown individually.
#
# O157:H7 logic:
#   O157:H7 = clean ECtyper serotype exactly equal to "O157:H7".
#   Non-O157:H7 = any other called serotype.
#   No serotype call = retained separately.
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

script02_dir <- file.path(wd, "results", "02_O157H7_cattle_focus")
table_dir <- file.path(script02_dir, "tables")
figure_dir <- file.path(script02_dir, "figures")

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

out_data6 <- file.path(table_dir, "data6_cattle_with_ectyper_metadata.csv")
out_data6_o157 <- file.path(table_dir, "data6_cattle_O157H7_only.csv")
out_summary <- file.path(table_dir, "data6_cattle_temporal_O157H7_summary.csv")
out_serotype_counts <- file.path(table_dir, "data6_cattle_serotype_counts.csv")
out_temporal_long <- file.path(table_dir, "data6_cattle_temporal_O157H7_long_counts.csv")

message("=== Script 02 started: ", Sys.time(), " ===")
message("[INFO] Working directory: ", wd)
message("[INFO] Metadata file: ", metadata_file)
message("[INFO] ECtyper file: ", ectyper_file)
message("[INFO] Script 02 output directory: ", script02_dir)

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

if (!file.exists(ectyper_file)) {
  stop("ECtyper output file not found: ", ectyper_file)
}

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

message("[INFO] data5 unique assemblies: ", nrow(metadata))

ectyper <- read_tsv(ectyper_file, show_col_types = FALSE)

required_ectyper_cols <- c("Name", "O-type", "H-type", "Serotype")
missing_ectyper_cols <- setdiff(required_ectyper_cols, names(ectyper))

if (length(missing_ectyper_cols) > 0) {
  stop("Missing required ECtyper columns: ", paste(missing_ectyper_cols, collapse = ", "))
}

ectyper <- ectyper %>%
  mutate(Name = str_squish(as.character(Name))) %>%
  filter(!is.na(Name), Name != "") %>%
  distinct(Name, .keep_all = TRUE)

message("[INFO] ECtyper unique genome rows: ", nrow(ectyper))

joined <- metadata %>%
  left_join(ectyper, by = c("Assembly" = "Name")) %>%
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
    O157H7_status = case_when(
      Serotype_clean == "O157:H7" ~ "O157:H7",
      Serotype_clean == "No_serotype_call" ~ "No serotype call",
      TRUE ~ "Non-O157:H7"
    ),
    is_O157H7 = Serotype_clean == "O157:H7"
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

data6 <- joined %>%
  filter(source_group_clean == "Cattle")

data6_o157 <- data6 %>%
  filter(is_O157H7)

write_csv(data6, out_data6)
write_csv(data6_o157, out_data6_o157)

message("[INFO] data6 cattle rows: ", nrow(data6))
message("[INFO] data6 O157:H7 cattle rows: ", nrow(data6_o157))

cattle_serotype_counts <- data6 %>%
  count(Serotype_clean, name = "n") %>%
  mutate(
    total_cattle = sum(n),
    percent_of_cattle = round(100 * n / total_cattle, 3)
  ) %>%
  arrange(desc(n), Serotype_clean)

write_csv(cattle_serotype_counts, out_serotype_counts)

temporal_data <- data6 %>%
  filter(!is.na(Collection_year_clean)) %>%
  filter(Collection_year_clean >= 1900, Collection_year_clean <= 2100) %>%
  mutate(
    Collection_year_temporal_group = make_temporal_year_group(Collection_year_clean)
  )

if (nrow(temporal_data) == 0) {
  warning("[WARN] No usable Collection_year values found. Temporal outputs skipped.")
} else {

  temporal_summary <- temporal_data %>%
    group_by(Collection_year_temporal_group) %>%
    summarise(
      total_cattle_genomes = n(),
      n_O157H7 = sum(O157H7_status == "O157:H7"),
      n_non_O157H7 = sum(O157H7_status == "Non-O157:H7"),
      n_no_serotype_call = sum(O157H7_status == "No serotype call"),
      n_serotype_called = n_O157H7 + n_non_O157H7,
      pct_O157H7_among_all_cattle = round(100 * n_O157H7 / total_cattle_genomes, 3),
      pct_O157H7_among_serotype_called = if_else(
        n_serotype_called > 0,
        round(100 * n_O157H7 / n_serotype_called, 3),
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      Collection_year_temporal_group = make_temporal_year_factor(Collection_year_temporal_group)
    ) %>%
    arrange(Collection_year_temporal_group) %>%
    rowwise() %>%
    mutate(
      prop_O157H7_among_all_cattle = n_O157H7 / total_cattle_genomes,
      ci = list(stats::binom.test(n_O157H7, total_cattle_genomes)$conf.int),
      prop_O157H7_ci_low = ci[[1]][1],
      prop_O157H7_ci_high = ci[[1]][2],
      pct_O157H7_ci_low = round(100 * prop_O157H7_ci_low, 3),
      pct_O157H7_ci_high = round(100 * prop_O157H7_ci_high, 3)
    ) %>%
    select(-ci) %>%
    ungroup()

  write_csv(temporal_summary, out_summary)

  temporal_long_counts <- temporal_data %>%
    count(Collection_year_temporal_group, O157H7_status, name = "n") %>%
    group_by(Collection_year_temporal_group) %>%
    mutate(
      year_total = sum(n),
      prop = n / year_total,
      percent = 100 * prop
    ) %>%
    ungroup() %>%
    mutate(
      Collection_year_temporal_group = make_temporal_year_factor(Collection_year_temporal_group),
      O157H7_status = factor(
        O157H7_status,
        levels = c("O157:H7", "Non-O157:H7", "No serotype call")
      )
    )

  write_csv(temporal_long_counts, out_temporal_long)

  p_prop <- ggplot(
    temporal_summary,
    aes(
      x = Collection_year_temporal_group,
      y = pct_O157H7_among_all_cattle,
      group = 1
    )
  ) +
    geom_col(width = 0.75, color = "gray20", linewidth = 0.2) +
    geom_text(
      aes(
        label = paste0(
          round(pct_O157H7_among_all_cattle, 1),
          "%\n",
          "n=",
          n_O157H7,
          "/",
          total_cattle_genomes
        )
      ),
      angle = 90,
      vjust = -0.25,
      hjust = 0.5,
      size = 3.0
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.22))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Temporal proportion of O157:H7 among cattle E. coli genomes",
      subtitle = "Years <=2000 grouped; denominator is all cattle genomes in data6",
      x = "Collection year",
      y = "O157:H7 proportion among cattle genomes",
      caption = "Labels show O157:H7 count / total cattle genomes per temporal group."
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = margin(t = 25, r = 20, b = 10, l = 10)
    )

  ggsave(
    file.path(figure_dir, "cattle_O157H7_temporal_proportion.pdf"),
    p_prop,
    width = 13,
    height = 7
  )

  ggsave(
    file.path(figure_dir, "cattle_O157H7_temporal_proportion.png"),
    p_prop,
    width = 13,
    height = 7,
    dpi = 600
  )

  p_line <- ggplot(
    temporal_summary,
    aes(
      x = as.integer(Collection_year_temporal_group),
      y = prop_O157H7_among_all_cattle,
      group = 1
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.8) +
    geom_text(
      aes(
        label = paste0(
          n_O157H7,
          "/",
          total_cattle_genomes
        )
      ),
      angle = 90,
      vjust = -0.55,
      hjust = 0.5,
      size = 2.8
    ) +
    scale_x_continuous(
      breaks = as.integer(temporal_summary$Collection_year_temporal_group),
      labels = as.character(temporal_summary$Collection_year_temporal_group)
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1.05),
      expand = expansion(mult = c(0.02, 0.12))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Temporal trend of O157:H7 among cattle E. coli genomes",
      subtitle = "Points show O157:H7 proportion by collection year",
      x = "Collection year",
      y = "O157:H7 proportion among cattle genomes",
      caption = "Years <=2000 grouped. Labels show O157:H7 count / total cattle genomes per temporal group."
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = margin(t = 25, r = 20, b = 10, l = 10)
    )

  ggsave(
    file.path(figure_dir, "cattle_O157H7_temporal_proportion_line.pdf"),
    p_line,
    width = 13,
    height = 7
  )

  ggsave(
    file.path(figure_dir, "cattle_O157H7_temporal_proportion_line.png"),
    p_line,
    width = 13,
    height = 7,
    dpi = 600
  )

  year_totals <- temporal_long_counts %>%
    distinct(Collection_year_temporal_group, year_total)

  p_counts <- ggplot(
    temporal_long_counts,
    aes(
      x = Collection_year_temporal_group,
      y = n,
      fill = O157H7_status
    )
  ) +
    geom_col(width = 0.85, color = "gray20", linewidth = 0.15) +
    geom_text(
      data = year_totals,
      aes(
        x = Collection_year_temporal_group,
        y = year_total,
        label = paste0("n=", year_total)
      ),
      inherit.aes = FALSE,
      angle = 90,
      vjust = -0.95,
      hjust = 0.5,
      size = 3.0
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.26))) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Temporal counts of O157:H7 and non-O157:H7 cattle E. coli",
      subtitle = "Years <=2000 grouped",
      x = "Collection year",
      y = "Number of cattle genomes",
      fill = "Serotype category"
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      plot.margin = margin(t = 25, r = 20, b = 10, l = 10)
    )

  ggsave(
    file.path(figure_dir, "cattle_O157H7_temporal_counts.pdf"),
    p_counts,
    width = 13,
    height = 7
  )

  ggsave(
    file.path(figure_dir, "cattle_O157H7_temporal_counts.png"),
    p_counts,
    width = 13,
    height = 7,
    dpi = 600
  )

  p_percent <- ggplot(
    temporal_long_counts,
    aes(
      x = Collection_year_temporal_group,
      y = prop,
      fill = O157H7_status
    )
  ) +
    geom_col(width = 0.85, color = "gray20", linewidth = 0.15) +
    geom_text(
      data = year_totals %>% mutate(label_y = 1.08),
      aes(
        x = Collection_year_temporal_group,
        y = label_y,
        label = paste0("n=", year_total)
      ),
      inherit.aes = FALSE,
      angle = 90,
      vjust = -0.65,
      hjust = 0.5,
      size = 3.0
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1.20),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Temporal composition of O157:H7 and non-O157:H7 cattle E. coli",
      subtitle = "Years <=2000 grouped; percent stacked by temporal group",
      x = "Collection year",
      y = "Percent of cattle genomes",
      fill = "Serotype category"
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      plot.margin = margin(t = 25, r = 20, b = 10, l = 10)
    )

  ggsave(
    file.path(figure_dir, "cattle_O157H7_temporal_percent_stacked.pdf"),
    p_percent,
    width = 13,
    height = 7
  )

  ggsave(
    file.path(figure_dir, "cattle_O157H7_temporal_percent_stacked.png"),
    p_percent,
    width = 13,
    height = 7,
    dpi = 600
  )
}

message("")
message("=== Script 02 summary ===")
message("[INFO] data6 cattle rows: ", nrow(data6))
message("[INFO] O157:H7 cattle rows: ", nrow(data6_o157))

message("")
message("=== Top cattle serotypes ===")
print(head(cattle_serotype_counts, 20), n = 20)

if (exists("temporal_summary")) {
  message("")
  message("=== Temporal O157:H7 summary ===")
  print(temporal_summary, n = Inf)
}

message("")
message("[DONE] Tables written to:")
message("  ", table_dir)

message("")
message("[DONE] Figures written to:")
message("  ", figure_dir)

message("=== Script 02 completed: ", Sys.time(), " ===")
