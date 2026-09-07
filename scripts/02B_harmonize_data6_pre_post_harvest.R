#!/usr/bin/env Rscript

# ============================================================
# Script 02B
# Purpose:
#   Harmonize "Isolation source" in data6.csv into:
#     1. pre-harvest and harvest
#     2. post-harvest
#
#   Create data6a.csv with an added harmonized source-stage column.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

message("=== SCRIPT 02B: HARMONIZE DATA6 INTO PRE/HARVEST VS POST-HARVEST ===")

# =========================
# Paths
# =========================
data6_file <- "/home/samuelajulo/epi/data6.csv"

out_dir <- "/home/samuelajulo/epi/Ecoli_epi/results/02B_data6_isolation_source_summary"
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

data6a_root_out <- "/home/samuelajulo/epi/data6a.csv"
data6a_repo_out <- file.path(table_dir, "data6a_with_pre_post_harvest_source.csv")

mapping_out <- file.path(table_dir, "data6a_isolation_source_pre_post_harvest_map.csv")
counts_out <- file.path(table_dir, "data6a_pre_post_harvest_counts.csv")
original_counts_out <- file.path(table_dir, "data6_original_isolation_source_counts.csv")

plot_pdf <- file.path(fig_dir, "data6a_pre_post_harvest_barplot.pdf")
plot_png <- file.path(fig_dir, "data6a_pre_post_harvest_barplot.png")

source_plot_pdf <- file.path(fig_dir, "data6_original_isolation_source_barplot_all.pdf")
source_plot_png <- file.path(fig_dir, "data6_original_isolation_source_barplot_all.png")

# =========================
# Input checks
# =========================
if (!file.exists(data6_file)) {
  stop("[ERROR] data6 file not found: ", data6_file, call. = FALSE)
}

data6 <- read_csv(data6_file, show_col_types = FALSE)

source_col <- "Isolation source"

if (!source_col %in% names(data6)) {
  stop(
    "[ERROR] Column not found: ", source_col,
    "\nAvailable columns:\n",
    paste(names(data6), collapse = "\n"),
    call. = FALSE
  )
}

# =========================
# Helper functions
# =========================
clean_source <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("\\([^)]*\\)", " ") %>%
    str_replace_all("[_\\-/]+", " ") %>%
    str_replace_all("[^a-z0-9 ]+", " ") %>%
    str_squish()
}

assign_pre_post_harvest <- function(x) {
  case_when(
    is.na(x) | x == "" ~ "pre-harvest and harvest",

    # Post-harvest: food, retail, dairy, processed or consumer-facing products
    str_detect(x, "ground beef|hamburger|burger|beef patty|lean ground beef") ~ "post-harvest",
    str_detect(x, "beef trim|trimming|trim\\b|boneless beef|raw beef") ~ "post-harvest",
    str_detect(x, "\\bbeef\\b") & !str_detect(x, "\\bcalf\\b|\\bveal\\b|\\bcattle\\b|\\bbovine\\b|\\bcow\\b|\\bsteer\\b|\\bheifer\\b") ~ "post-harvest",
    str_detect(x, "\\bmilk\\b|\\bcheese\\b|retail") ~ "post-harvest",

    # Pre-harvest and harvest: animal, farm, feedlot, feces, manure, hide, carcass, rectal
    str_detect(x, "\\bmanure\\b|\\bdung\\b|\\bsoil\\b") ~ "pre-harvest and harvest",
    str_detect(x, "\\bfeces\\b|\\bfecal\\b|\\brectal\\b") ~ "pre-harvest and harvest",
    str_detect(x, "\\bhide\\b|\\bcarcass\\b") ~ "pre-harvest and harvest",
    str_detect(x, "\\bbovine\\b|\\bcattle\\b|\\bcalf\\b|\\bveal\\b|\\bcow\\b|\\bsteer\\b|\\bheifer\\b|\\banimal\\b|\\bfeedlot\\b") ~ "pre-harvest and harvest",

    # Conservative fallback for this cattle-focused data6
    TRUE ~ "pre-harvest and harvest"
  )
}

# =========================
# Create data6a
# =========================
data6a <- data6 %>%
  mutate(
    Isolation_source_original = .data[[source_col]],
    Isolation_source_clean = clean_source(.data[[source_col]]),
    Isolation_source_stage = assign_pre_post_harvest(Isolation_source_clean)
  )

write_csv(data6a, data6a_root_out)
write_csv(data6a, data6a_repo_out)

# =========================
# Mapping/audit table
# =========================
mapping_table <- data6a %>%
  count(
    Isolation_source_original,
    Isolation_source_clean,
    Isolation_source_stage,
    name = "n",
    sort = TRUE
  ) %>%
  mutate(percent = 100 * n / sum(n))

write_csv(mapping_table, mapping_out)

# =========================
# Counts
# =========================
stage_counts <- data6a %>%
  count(Isolation_source_stage, name = "n", sort = TRUE) %>%
  mutate(
    percent = 100 * n / sum(n),
    Isolation_source_stage = factor(
      Isolation_source_stage,
      levels = rev(Isolation_source_stage)
    )
  )

write_csv(stage_counts, counts_out)

original_counts <- data6a %>%
  mutate(
    Isolation_source_original = as.character(Isolation_source_original),
    Isolation_source_original = str_squish(Isolation_source_original),
    Isolation_source_original = if_else(
      is.na(Isolation_source_original) | Isolation_source_original == "",
      "Missing/blank",
      Isolation_source_original
    )
  ) %>%
  count(Isolation_source_original, name = "n", sort = TRUE) %>%
  mutate(
    percent = 100 * n / sum(n),
    Isolation_source_wrapped = str_wrap(Isolation_source_original, width = 45),
    Isolation_source_wrapped = factor(
      Isolation_source_wrapped,
      levels = rev(unique(Isolation_source_wrapped))
    )
  )

write_csv(
  original_counts %>%
    select(Isolation_source = Isolation_source_original, n, percent),
  original_counts_out
)

# =========================
# Plot stage counts
# =========================
p_stage <- ggplot(stage_counts, aes(x = Isolation_source_stage, y = n)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw(base_size = 13) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Pre-harvest/harvest versus post-harvest distribution in data6a",
    subtitle = paste0("n = ", comma(sum(stage_counts$n)), " assemblies"),
    x = "Harmonized isolation source stage",
    y = "Number of assemblies"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(plot_pdf, p_stage, width = 8, height = 4.5)
ggsave(plot_png, p_stage, width = 8, height = 4.5, dpi = 300)

# =========================
# Plot original sources too
# =========================
n_categories <- nrow(original_counts)
source_plot_height <- max(6, min(40, 0.28 * n_categories + 3))

p_source <- ggplot(original_counts, aes(x = Isolation_source_wrapped, y = n)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw(base_size = 12) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Original isolation source distribution in data6",
    subtitle = paste0("All unique values before harmonization; n = ", comma(sum(original_counts$n)), " assemblies"),
    x = "Original isolation source",
    y = "Number of assemblies"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank()
  )

ggsave(source_plot_pdf, p_source, width = 11, height = source_plot_height, limitsize = FALSE)
ggsave(source_plot_png, p_source, width = 11, height = source_plot_height, dpi = 300, limitsize = FALSE)

message("[DONE] data6a created with Isolation_source_stage column.")
message("[OUTPUT] Root data6a: ", data6a_root_out)
message("[OUTPUT] Repo data6a: ", data6a_repo_out)
message("[OUTPUT] Mapping table: ", mapping_out)
message("[OUTPUT] Stage counts: ", counts_out)
message("[OUTPUT] Stage plot PDF: ", plot_pdf)

message("")
message("Stage counts:")
print(stage_counts)

message("")
message("Mapping preview:")
print(head(mapping_table, 30))
