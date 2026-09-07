#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

message("=== SCRIPT 06: DATA6 ISOLATION SOURCE BAR PLOT ===")

data6_file <- "/home/samuelajulo/epi/data6.csv"

out_dir <- path.expand("~/epi/marker_screen/results/06_data6_isolation_source_summary")
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

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

source_counts <- data6 %>%
  mutate(
    Isolation_source_clean = as.character(.data[[source_col]]),
    Isolation_source_clean = str_squish(Isolation_source_clean),
    Isolation_source_clean = if_else(
      is.na(Isolation_source_clean) | Isolation_source_clean == "",
      "Missing/blank",
      Isolation_source_clean
    )
  ) %>%
  count(Isolation_source_clean, name = "n", sort = TRUE) %>%
  mutate(
    percent = 100 * n / sum(n),
    Isolation_source_wrapped = str_wrap(Isolation_source_clean, width = 45),
    Isolation_source_wrapped = factor(
      Isolation_source_wrapped,
      levels = rev(unique(Isolation_source_wrapped))
    )
  )

write_csv(
  source_counts %>%
    select(Isolation_source = Isolation_source_clean, n, percent),
  file.path(table_dir, "data6_isolation_source_counts.csv")
)

n_categories <- nrow(source_counts)

plot_height <- max(6, min(40, 0.28 * n_categories + 3))

p <- ggplot(source_counts, aes(x = Isolation_source_wrapped, y = n)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw(base_size = 12) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Isolation source distribution in data6",
    subtitle = paste0("All unique values in the Isolation source column; n = ", comma(sum(source_counts$n)), " assemblies"),
    x = "Isolation source",
    y = "Number of assemblies"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "data6_isolation_source_barplot_all.pdf"),
  plot = p,
  width = 11,
  height = plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(fig_dir, "data6_isolation_source_barplot_all.png"),
  plot = p,
  width = 11,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE
)

message("[DONE] Isolation source bar plot complete.")
message("[OUTPUT] Counts table: ", file.path(table_dir, "data6_isolation_source_counts.csv"))
message("[OUTPUT] PDF: ", file.path(fig_dir, "data6_isolation_source_barplot_all.pdf"))
message("[OUTPUT] PNG: ", file.path(fig_dir, "data6_isolation_source_barplot_all.png"))

message("")
message("Top isolation sources:")
print(head(source_counts %>% select(Isolation_source_clean, n, percent), 20))
