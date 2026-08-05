#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 02E
# Temporal production-system and supply-chain-stage plots
# for data6a.csv
# ============================================================
#
# Goals:
#   1. Yearly count plots for Beef, Dairy, Cattle-unspecified
#      with all years <=2000 grouped as "<=2000".
#
#   2. Four-year binned plots:
#      <=2000, 2001-2004, 2005-2008, 2009-2012, etc.
#      The final bin stops at the final year present in the dataset.
#      Example: if final year is 2026, final bin is 2025-2026,
#      not 2025-2028.
#
#   3. Repeat across supply-chain subcategories:
#      Pre-harvest, Peri-harvest, Post-harvest.
#
#   4. Focus only on Pre-harvest.
#
#   5. Focus only on Beef and Cattle-unspecified Pre-harvest.
#
#   6. Add total pie charts for:
#      - Production system
#      - Supply-chain stage
#      - Production-stage group
#      - Pre-harvest production system
#      - Pre-harvest Beef vs Cattle-unspecified
#      - Source matrix class, when available
#
# Inputs:
#   /home/samuelajulo/epi/data6a.csv
#
# Outputs:
#   /home/samuelajulo/epi/Ecoli_epi/results/02E_data6a_temporal_production_stage
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

message("=== SCRIPT 02E: TEMPORAL PRODUCTION/STAGE PLOTS STARTED ===")
message("[TIME] ", Sys.time())

# ============================================================
# 1. Paths
# ============================================================

data6a_file <- "/home/samuelajulo/epi/data6a.csv"

base_dir <- "/home/samuelajulo/epi/Ecoli_epi"
out_dir <- file.path(base_dir, "results", "02E_data6a_temporal_production_stage")
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")

fig_yearly_dir <- file.path(fig_dir, "yearly_counts")
fig_4yr_dir <- file.path(fig_dir, "four_year_bins")
fig_focus_dir <- file.path(fig_dir, "focused_panels")
fig_pie_dir <- file.path(fig_dir, "total_piecharts")
chart_csv_dir <- file.path(out_dir, "chart_csv_files")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_yearly_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_4yr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_focus_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_pie_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chart_csv_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data6a_file)) {
  stop("[ERROR] data6a not found: ", data6a_file, call. = FALSE)
}

# ============================================================
# 2. Helper functions
# ============================================================

extract_year_num <- function(x) {
  x_chr <- as.character(x)

  year_chr <- str_extract(x_chr, "(18|19|20)\\d{2}")
  year_num <- suppressWarnings(as.integer(year_chr))

  le2000 <- str_detect(
    str_to_lower(x_chr),
    "<=\\s*2000|≤\\s*2000|before\\s*2000|pre\\s*2000"
  )

  year_num[le2000] <- 2000L
  year_num
}

pick_year_column <- function(df) {
  candidates <- c(
    "Collection_year",
    "Collection year",
    "collection_year",
    "Year",
    "year",
    "Collection date",
    "Collection_date",
    "collection_date",
    "collection date",
    "Date",
    "date",
    "Isolation date",
    "isolation_date",
    "Isolate collection date",
    "Collection_year_temporal_group"
  )

  present <- candidates[candidates %in% names(df)]

  if (length(present) == 0) {
    pattern_cols <- names(df)[str_detect(str_to_lower(names(df)), "year|date")]
    present <- pattern_cols
  }

  if (length(present) == 0) {
    stop(
      "[ERROR] No year/date column found.\nAvailable columns:\n",
      paste(names(df), collapse = "\n"),
      call. = FALSE
    )
  }

  scores <- sapply(present, function(col) {
    y <- extract_year_num(df[[col]])
    sum(!is.na(y) & y >= 1900 & y <= 2100)
  })

  best <- present[which.max(scores)]

  if (max(scores) == 0) {
    stop(
      "[ERROR] Year/date columns were found, but no valid years could be parsed.\nCandidate columns:\n",
      paste(present, collapse = ", "),
      call. = FALSE
    )
  }

  best
}

standardize_production <- function(x) {
  x <- str_squish(as.character(x))
  xl <- str_to_lower(x)

  case_when(
    str_detect(xl, "^beef$|beef") ~ "Beef",
    str_detect(xl, "^dairy$|dairy") ~ "Dairy",
    str_detect(xl, "cattle.*unspecified|unspecified") ~ "Cattle-unspecified",
    is.na(x) | x == "" ~ "Cattle-unspecified",
    TRUE ~ x
  )
}

standardize_stage <- function(x) {
  x <- str_squish(as.character(x))
  xl <- str_to_lower(x)

  case_when(
    str_detect(xl, "post") ~ "Post-harvest",
    str_detect(xl, "peri|per\\s*harvest|harvest-interface|interface") ~ "Peri-harvest",
    str_detect(xl, "pre") ~ "Pre-harvest",
    is.na(x) | x == "" ~ "Pre-harvest",
    TRUE ~ x
  )
}

standardize_matrix_class <- function(x) {
  x <- str_squish(as.character(x))
  xl <- str_to_lower(x)

  case_when(
    is.na(x) | x == "" ~ "Other/unknown",
    str_detect(xl, "feces|fecal|rectal") ~ "Feces/rectal",
    str_detect(xl, "manure|environment|soil|dung") ~ "Manure/environment",
    str_detect(xl, "hide") ~ "Hide",
    str_detect(xl, "carcass") ~ "Carcass",
    str_detect(xl, "meat|product|beef|veal|trim|ground|burger") ~ "Meat/product",
    str_detect(xl, "milk|dairy|cheese") ~ "Milk/dairy product",
    str_detect(xl, "calf") ~ "Calf-associated non-fecal",
    str_detect(xl, "cow") ~ "Cow-associated non-fecal",
    str_detect(xl, "bovine|cattle") ~ "Bovine/cattle non-fecal",
    TRUE ~ x
  )
}

make_year_plot_label <- function(y) {
  case_when(
    is.na(y) ~ NA_character_,
    y <= 2000 ~ "<=2000",
    TRUE ~ as.character(y)
  )
}

make_4yr_bin <- function(y, final_year) {
  out <- rep(NA_character_, length(y))

  valid <- !is.na(y)
  out[valid & y <= 2000] <- "<=2000"

  valid_after <- valid & y > 2000
  start <- 2001 + 4 * floor((y[valid_after] - 2001) / 4)
  end <- pmin(start + 3, final_year)

  # If the final bin contains only one year, label it as ">=YEAR"
  # Example: 2025-2025 becomes >=2025
  out[valid_after] <- ifelse(
    start == end,
    paste0(">=", start),
    paste0(start, "-", end)
  )

  out
}

make_year_levels <- function(y) {
  y <- y[!is.na(y)]
  y <- y[y >= 1900 & y <= 2100]

  if (length(y) == 0) {
    return(character(0))
  }

  max_year <- max(y)

  after_2000 <- if (max_year > 2000) {
    as.character(seq(2001, max_year, by = 1))
  } else {
    character(0)
  }

  if (any(y <= 2000)) {
    c("<=2000", after_2000)
  } else {
    after_2000
  }
}

make_4yr_levels <- function(y) {
  y <- y[!is.na(y)]
  y <- y[y >= 1900 & y <= 2100]

  if (length(y) == 0) {
    return(character(0))
  }

  max_year <- max(y)

  after_2000 <- character(0)

  if (max_year > 2000) {
    starts <- seq(
      2001,
      2001 + 4 * floor((max_year - 2001) / 4),
      by = 4
    )

    ends <- pmin(starts + 3, max_year)

    # If the final bin contains only one year, label it as ">=YEAR"
    # Example: 2025-2025 becomes >=2025
    after_2000 <- ifelse(
      starts == ends,
      paste0(">=", starts),
      paste0(starts, "-", ends)
    )
  }

  if (any(y <= 2000)) {
    c("<=2000", after_2000)
  } else {
    after_2000
  }
}

safe_filename <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

theme_nature <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      plot.caption = element_text(size = base_size - 2, color = "gray30", hjust = 0),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      legend.position = "top",
      legend.justification = "left",
      legend.box.just = "left",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      panel.spacing = grid::unit(1.1, "lines")
    )
}

theme_pie_publication <- function(base_size = 11) {
  theme_void(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      plot.caption = element_text(size = base_size - 2, color = "gray30", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.position = "right"
    )
}

save_plot_all <- function(plot, path_no_ext, width, height) {
  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

  ggsave(
    filename = paste0(path_no_ext, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = pdf_device
  )

  ggsave(
    filename = paste0(path_no_ext, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )

  ggsave(
    filename = paste0(path_no_ext, ".tiff"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}

make_stacked_bar <- function(
  df,
  x_col,
  title,
  subtitle,
  x_label,
  y_label = "Counts",
  show_total_labels = FALSE,
  facet_col = NULL,
  rotate_x = TRUE,
  keep_all_production_levels = TRUE
) {
  if (keep_all_production_levels) {
    fill_values <- c(
      "Beef" = "#D55E00",
      "Dairy" = "#0072B2",
      "Cattle-unspecified" = "#999999"
    )

    legend_breaks <- c("Beef", "Dairy", "Cattle-unspecified")

    df <- df %>%
      mutate(
        Production_system = factor(
          as.character(Production_system),
          levels = legend_breaks
        )
      )
  } else {
    fill_values <- c(
      "Beef" = "#D55E00",
      "Cattle-unspecified" = "#999999"
    )

    legend_breaks <- c("Beef", "Cattle-unspecified")

    df <- df %>%
      filter(as.character(Production_system) %in% legend_breaks) %>%
      mutate(
        Production_system = factor(
          as.character(Production_system),
          levels = legend_breaks
        )
      ) %>%
      droplevels()
  }

  p <- ggplot(
    df,
    aes(x = .data[[x_col]], y = n, fill = Production_system)
  ) +
    geom_col(width = 0.82, color = "white", linewidth = 0.15) +
    scale_fill_manual(
      values = fill_values,
      breaks = legend_breaks,
      limits = legend_breaks,
      drop = TRUE
    ) +
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0, 0.14))
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label,
      fill = "Production system"
    ) +
    theme_nature(base_size = 11)

  if (rotate_x) {
    p <- p +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      )
  }

  if (!is.null(facet_col)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)), ncol = 1, scales = "free_y")
  }

  if (show_total_labels) {
    group_vars <- c(x_col)

    if (!is.null(facet_col)) {
      group_vars <- c(group_vars, facet_col)
    }

    total_df <- df %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(total_n = sum(n), .groups = "drop") %>%
      filter(total_n > 0)

    p <- p +
      geom_text(
        data = total_df,
        aes(
          x = .data[[x_col]],
          y = total_n,
          label = paste0("n=", total_n)
        ),
        inherit.aes = FALSE,
        vjust = -0.35,
        size = 3.0
      )
  }

  p
}

make_pie_plot <- function(
  df,
  category_col,
  title,
  subtitle,
  legend_title,
  min_label_percent = 4
) {
  df_plot <- df %>%
    filter(!is.na(.data[[category_col]]), .data[[category_col]] != "") %>%
    filter(n > 0) %>%
    mutate(
      category_clean = as.character(.data[[category_col]])
    ) %>%
    group_by(category_clean) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    arrange(desc(n), category_clean) %>%
    mutate(
      total = sum(n),
      percent = 100 * n / total,
      percent_label = paste0(round(percent, 1), "%"),
      slice_label = if_else(
        percent >= min_label_percent,
        paste0("n=", n, "\n", percent_label),
        ""
      ),
      legend_label = paste0(category_clean, " (n=", n, ", ", percent_label, ")"),
      legend_label = factor(legend_label, levels = legend_label)
    )

  ggplot(df_plot, aes(x = "", y = n, fill = legend_label)) +
    geom_col(width = 1, color = "white", linewidth = 0.35) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = slice_label),
      position = position_stack(vjust = 0.5),
      size = 3.2,
      lineheight = 0.9,
      color = "black"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      fill = legend_title,
      caption = paste0("Total n = ", comma(sum(df_plot$n)))
    ) +
    theme_pie_publication(base_size = 11) +
    guides(fill = guide_legend(ncol = 1))
}

# ============================================================
# 3. Read and clean data
# ============================================================

data6a <- read_csv(data6a_file, show_col_types = FALSE)

required_cols <- c("Production_system", "Supply_chain_stage")
missing_cols <- setdiff(required_cols, names(data6a))

if (length(missing_cols) > 0) {
  stop(
    "[ERROR] Missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    "\nRun updated 02B first.\nAvailable columns:\n",
    paste(names(data6a), collapse = "\n"),
    call. = FALSE
  )
}

year_col <- pick_year_column(data6a)

message("[INFO] Using year/date column: ", year_col)

data6a_clean <- data6a %>%
  mutate(
    Collection_year_source_column = year_col,
    Collection_year_raw = as.character(.data[[year_col]]),
    Collection_year_num = extract_year_num(Collection_year_raw),
    Production_system = standardize_production(Production_system),
    Supply_chain_stage = standardize_stage(Supply_chain_stage),
    Production_stage_group = paste(Production_system, Supply_chain_stage, sep = " | ")
  )

if ("Source_matrix_class" %in% names(data6a_clean)) {
  data6a_clean <- data6a_clean %>%
    mutate(Source_matrix_class_clean = standardize_matrix_class(Source_matrix_class))
} else if ("Source_matrix_class_02D" %in% names(data6a_clean)) {
  data6a_clean <- data6a_clean %>%
    mutate(Source_matrix_class_clean = standardize_matrix_class(Source_matrix_class_02D))
} else {
  data6a_clean <- data6a_clean %>%
    mutate(Source_matrix_class_clean = NA_character_)
}

invalid_years <- data6a_clean %>%
  filter(
    is.na(Collection_year_num) |
      Collection_year_num < 1900 |
      Collection_year_num > 2100
  )

valid_data <- data6a_clean %>%
  filter(
    !is.na(Collection_year_num),
    Collection_year_num >= 1900,
    Collection_year_num <= 2100,
    Production_system %in% c("Beef", "Dairy", "Cattle-unspecified"),
    Supply_chain_stage %in% c("Pre-harvest", "Peri-harvest", "Post-harvest")
  )

final_year <- max(valid_data$Collection_year_num, na.rm = TRUE)

valid_data <- valid_data %>%
  mutate(
    Year_plot = make_year_plot_label(Collection_year_num),
    Year_4yr_bin = make_4yr_bin(Collection_year_num, final_year = final_year)
  )

year_levels <- make_year_levels(valid_data$Collection_year_num)
bin4_levels <- make_4yr_levels(valid_data$Collection_year_num)

production_levels <- c("Beef", "Dairy", "Cattle-unspecified")
stage_levels <- c("Pre-harvest", "Peri-harvest", "Post-harvest")

valid_data <- valid_data %>%
  mutate(
    Year_plot = factor(Year_plot, levels = year_levels),
    Year_4yr_bin = factor(Year_4yr_bin, levels = bin4_levels),
    Production_system = factor(Production_system, levels = production_levels),
    Supply_chain_stage = factor(Supply_chain_stage, levels = stage_levels)
  )

write_csv(
  data6a_clean,
  file.path(table_dir, "data6a_temporal_input_with_clean_years.csv")
)

write_csv(
  invalid_years,
  file.path(table_dir, "data6a_missing_or_invalid_years.csv")
)

message("[INFO] Total rows in data6a: ", nrow(data6a))
message("[INFO] Valid rows for temporal plotting: ", nrow(valid_data))
message("[INFO] Missing/invalid year rows: ", nrow(invalid_years))
message("[INFO] Final year detected: ", final_year)
message("[INFO] Four-year bins: ", paste(bin4_levels, collapse = ", "))

# ============================================================
# 4. Count tables
# ============================================================

yearly_prod <- valid_data %>%
  count(Year_plot, Production_system, name = "n") %>%
  complete(
    Year_plot = factor(year_levels, levels = year_levels),
    Production_system = factor(production_levels, levels = production_levels),
    fill = list(n = 0)
  ) %>%
  group_by(Year_plot) %>%
  mutate(
    year_total = sum(n),
    percent_within_year = if_else(year_total > 0, round(100 * n / year_total, 2), 0)
  ) %>%
  ungroup()

bin4_prod <- valid_data %>%
  count(Year_4yr_bin, Production_system, name = "n") %>%
  complete(
    Year_4yr_bin = factor(bin4_levels, levels = bin4_levels),
    Production_system = factor(production_levels, levels = production_levels),
    fill = list(n = 0)
  ) %>%
  group_by(Year_4yr_bin) %>%
  mutate(
    bin_total = sum(n),
    percent_within_bin = if_else(bin_total > 0, round(100 * n / bin_total, 2), 0)
  ) %>%
  ungroup()

yearly_stage <- valid_data %>%
  count(Year_plot, Supply_chain_stage, Production_system, name = "n") %>%
  complete(
    Year_plot = factor(year_levels, levels = year_levels),
    Supply_chain_stage = factor(stage_levels, levels = stage_levels),
    Production_system = factor(production_levels, levels = production_levels),
    fill = list(n = 0)
  ) %>%
  group_by(Year_plot, Supply_chain_stage) %>%
  mutate(
    year_stage_total = sum(n),
    percent_within_year_stage = if_else(
      year_stage_total > 0,
      round(100 * n / year_stage_total, 2),
      0
    )
  ) %>%
  ungroup()

bin4_stage <- valid_data %>%
  count(Year_4yr_bin, Supply_chain_stage, Production_system, name = "n") %>%
  complete(
    Year_4yr_bin = factor(bin4_levels, levels = bin4_levels),
    Supply_chain_stage = factor(stage_levels, levels = stage_levels),
    Production_system = factor(production_levels, levels = production_levels),
    fill = list(n = 0)
  ) %>%
  group_by(Year_4yr_bin, Supply_chain_stage) %>%
  mutate(
    bin_stage_total = sum(n),
    percent_within_bin_stage = if_else(
      bin_stage_total > 0,
      round(100 * n / bin_stage_total, 2),
      0
    )
  ) %>%
  ungroup()

preharvest_yearly <- valid_data %>%
  filter(Supply_chain_stage == "Pre-harvest") %>%
  count(Year_plot, Production_system, name = "n") %>%
  complete(
    Year_plot = factor(year_levels, levels = year_levels),
    Production_system = factor(production_levels, levels = production_levels),
    fill = list(n = 0)
  ) %>%
  group_by(Year_plot) %>%
  mutate(
    year_total = sum(n),
    percent_within_year = if_else(year_total > 0, round(100 * n / year_total, 2), 0)
  ) %>%
  ungroup()

preharvest_bin4 <- valid_data %>%
  filter(Supply_chain_stage == "Pre-harvest") %>%
  count(Year_4yr_bin, Production_system, name = "n") %>%
  complete(
    Year_4yr_bin = factor(bin4_levels, levels = bin4_levels),
    Production_system = factor(production_levels, levels = production_levels),
    fill = list(n = 0)
  ) %>%
  group_by(Year_4yr_bin) %>%
  mutate(
    bin_total = sum(n),
    percent_within_bin = if_else(bin_total > 0, round(100 * n / bin_total, 2), 0)
  ) %>%
  ungroup()

preharvest_focus_products <- c("Beef", "Cattle-unspecified")

preharvest_focus_yearly <- valid_data %>%
  filter(
    Supply_chain_stage == "Pre-harvest",
    Production_system %in% preharvest_focus_products
  ) %>%
  count(Year_plot, Production_system, name = "n") %>%
  complete(
    Year_plot = factor(year_levels, levels = year_levels),
    Production_system = factor(preharvest_focus_products, levels = preharvest_focus_products),
    fill = list(n = 0)
  ) %>%
  group_by(Year_plot) %>%
  mutate(
    year_total = sum(n),
    percent_within_year = if_else(year_total > 0, round(100 * n / year_total, 2), 0)
  ) %>%
  ungroup()

preharvest_focus_bin4 <- valid_data %>%
  filter(
    Supply_chain_stage == "Pre-harvest",
    Production_system %in% preharvest_focus_products
  ) %>%
  count(Year_4yr_bin, Production_system, name = "n") %>%
  complete(
    Year_4yr_bin = factor(bin4_levels, levels = bin4_levels),
    Production_system = factor(preharvest_focus_products, levels = preharvest_focus_products),
    fill = list(n = 0)
  ) %>%
  group_by(Year_4yr_bin) %>%
  mutate(
    bin_total = sum(n),
    percent_within_bin = if_else(bin_total > 0, round(100 * n / bin_total, 2), 0)
  ) %>%
  ungroup()

# ============================================================
# 5. Total category count tables for pie charts
# ============================================================

pie_total_production <- valid_data %>%
  count(Production_system, name = "n") %>%
  mutate(
    total = sum(n),
    percent = round(100 * n / total, 2)
  ) %>%
  arrange(desc(n))

pie_total_stage <- valid_data %>%
  count(Supply_chain_stage, name = "n") %>%
  mutate(
    total = sum(n),
    percent = round(100 * n / total, 2)
  ) %>%
  arrange(desc(n))

pie_total_production_stage <- valid_data %>%
  count(Production_stage_group, name = "n") %>%
  mutate(
    total = sum(n),
    percent = round(100 * n / total, 2)
  ) %>%
  arrange(desc(n))

pie_preharvest_production <- valid_data %>%
  filter(Supply_chain_stage == "Pre-harvest") %>%
  count(Production_system, name = "n") %>%
  mutate(
    total = sum(n),
    percent = round(100 * n / total, 2)
  ) %>%
  arrange(desc(n))

pie_preharvest_focus <- valid_data %>%
  filter(
    Supply_chain_stage == "Pre-harvest",
    Production_system %in% preharvest_focus_products
  ) %>%
  count(Production_system, name = "n") %>%
  mutate(
    total = sum(n),
    percent = round(100 * n / total, 2)
  ) %>%
  arrange(desc(n))

pie_matrix_class <- valid_data %>%
  filter(!is.na(Source_matrix_class_clean), Source_matrix_class_clean != "") %>%
  count(Source_matrix_class_clean, name = "n") %>%
  mutate(
    total = sum(n),
    percent = round(100 * n / total, 2)
  ) %>%
  arrange(desc(n))

write_csv(yearly_prod, file.path(table_dir, "data6a_yearly_counts_by_production_system.csv"))
write_csv(bin4_prod, file.path(table_dir, "data6a_4yr_counts_by_production_system.csv"))
write_csv(yearly_stage, file.path(table_dir, "data6a_yearly_counts_by_production_stage.csv"))
write_csv(bin4_stage, file.path(table_dir, "data6a_4yr_counts_by_production_stage.csv"))
write_csv(preharvest_yearly, file.path(table_dir, "data6a_preharvest_yearly_counts.csv"))
write_csv(preharvest_bin4, file.path(table_dir, "data6a_preharvest_4yr_counts.csv"))
write_csv(preharvest_focus_yearly, file.path(table_dir, "data6a_preharvest_beef_cattle_unspecified_yearly_counts.csv"))
write_csv(preharvest_focus_bin4, file.path(table_dir, "data6a_preharvest_beef_cattle_unspecified_4yr_counts.csv"))

write_csv(pie_total_production, file.path(table_dir, "data6a_pie_total_by_production_system.csv"))
write_csv(pie_total_stage, file.path(table_dir, "data6a_pie_total_by_supply_chain_stage.csv"))
write_csv(pie_total_production_stage, file.path(table_dir, "data6a_pie_total_by_production_stage_group.csv"))
write_csv(pie_preharvest_production, file.path(table_dir, "data6a_pie_preharvest_by_production_system.csv"))
write_csv(pie_preharvest_focus, file.path(table_dir, "data6a_pie_preharvest_beef_cattle_unspecified.csv"))
write_csv(pie_matrix_class, file.path(table_dir, "data6a_pie_total_by_source_matrix_class.csv"))

# ============================================================
# 5B. Chart-specific CSV files
# ============================================================

write_chart_csv <- function(df, filename) {
  write_csv(df, file.path(chart_csv_dir, filename))
}

write_chart_csv(
  yearly_prod,
  "01_yearly_counts_by_production_system.csv"
)

write_chart_csv(
  bin4_prod,
  "02_four_year_counts_by_production_system.csv"
)

write_chart_csv(
  yearly_stage,
  "03_yearly_counts_by_production_system_and_stage.csv"
)

write_chart_csv(
  bin4_stage,
  "04_four_year_counts_by_production_system_and_stage.csv"
)

write_chart_csv(
  preharvest_yearly,
  "05_preharvest_yearly_counts_by_production_system.csv"
)

write_chart_csv(
  preharvest_bin4,
  "06_preharvest_four_year_counts_by_production_system.csv"
)

write_chart_csv(
  preharvest_focus_yearly,
  "07_preharvest_beef_cattle_unspecified_yearly_counts.csv"
)

write_chart_csv(
  preharvest_focus_bin4,
  "08_preharvest_beef_cattle_unspecified_four_year_counts.csv"
)

write_chart_csv(
  pie_total_production,
  "09_pie_total_by_production_system.csv"
)

write_chart_csv(
  pie_total_stage,
  "10_pie_total_by_supply_chain_stage.csv"
)

write_chart_csv(
  pie_total_production_stage,
  "11_pie_total_by_production_stage_group.csv"
)

write_chart_csv(
  pie_preharvest_production,
  "12_pie_preharvest_by_production_system.csv"
)

write_chart_csv(
  pie_preharvest_focus,
  "13_pie_preharvest_beef_cattle_unspecified.csv"
)

write_chart_csv(
  pie_matrix_class,
  "14_pie_total_by_source_matrix_class.csv"
)

chart_csv_manifest <- tibble::tribble(
  ~chart_id, ~figure_prefix, ~csv_file, ~description, ~filter_used,
  "01", "01_yearly_counts_by_production_system", "01_yearly_counts_by_production_system.csv", "Annual counts by Beef, Dairy, and Cattle-unspecified", "All valid temporal records",
  "02", "02_four_year_counts_by_production_system", "02_four_year_counts_by_production_system.csv", "Four-year binned counts by Beef, Dairy, and Cattle-unspecified", "All valid temporal records",
  "03", "03_yearly_counts_by_production_system_and_stage", "03_yearly_counts_by_production_system_and_stage.csv", "Annual counts by production system faceted by supply-chain stage", "All valid temporal records",
  "04", "04_four_year_counts_by_production_system_and_stage", "04_four_year_counts_by_production_system_and_stage.csv", "Four-year binned counts by production system faceted by supply-chain stage", "All valid temporal records",
  "05", "05_preharvest_yearly_counts_by_production_system", "05_preharvest_yearly_counts_by_production_system.csv", "Annual pre-harvest counts by production system", "Supply_chain_stage == Pre-harvest",
  "06", "06_preharvest_four_year_counts_by_production_system", "06_preharvest_four_year_counts_by_production_system.csv", "Four-year pre-harvest counts by production system", "Supply_chain_stage == Pre-harvest",
  "07", "07_preharvest_beef_cattle_unspecified_yearly_counts", "07_preharvest_beef_cattle_unspecified_yearly_counts.csv", "Annual pre-harvest counts for Beef and Cattle-unspecified only", "Pre-harvest; Beef and Cattle-unspecified only",
  "08", "08_preharvest_beef_cattle_unspecified_four_year_counts", "08_preharvest_beef_cattle_unspecified_four_year_counts.csv", "Four-year pre-harvest counts for Beef and Cattle-unspecified only", "Pre-harvest; Beef and Cattle-unspecified only",
  "09", "09_pie_total_by_production_system", "09_pie_total_by_production_system.csv", "Total distribution by production system", "All valid temporal records",
  "10", "10_pie_total_by_supply_chain_stage", "10_pie_total_by_supply_chain_stage.csv", "Total distribution by supply-chain stage", "All valid temporal records",
  "11", "11_pie_total_by_production_stage_group", "11_pie_total_by_production_stage_group.csv", "Total distribution by production-stage group", "All valid temporal records",
  "12", "12_pie_preharvest_by_production_system", "12_pie_preharvest_by_production_system.csv", "Pre-harvest distribution by production system", "Supply_chain_stage == Pre-harvest",
  "13", "13_pie_preharvest_beef_cattle_unspecified", "13_pie_preharvest_beef_cattle_unspecified.csv", "Pre-harvest distribution for Beef and Cattle-unspecified only", "Pre-harvest; Beef and Cattle-unspecified only",
  "14", "14_pie_total_by_source_matrix_class", "14_pie_total_by_source_matrix_class.csv", "Total distribution by source matrix class", "All valid temporal records with matrix class"
)

write_csv(
  chart_csv_manifest,
  file.path(chart_csv_dir, "00_chart_csv_manifest.csv")
)

write_csv(
  chart_csv_manifest,
  file.path(table_dir, "data6a_chart_csv_manifest.csv")
)


# ============================================================
# 6. Bar figures
# ============================================================

p_yearly_prod <- make_stacked_bar(
  df = yearly_prod,
  x_col = "Year_plot",
  title = "Annual O157:H7 cattle-associated counts by production system",
  subtitle = paste0("Years <=2000 collapsed; n = ", comma(nrow(valid_data)), " assemblies"),
  x_label = "Collection year",
  y_label = "Counts",
  show_total_labels = FALSE,
  rotate_x = TRUE
)

save_plot_all(
  p_yearly_prod,
  file.path(fig_yearly_dir, "01_yearly_counts_by_production_system"),
  width = 11,
  height = 5.5
)

p_4yr_prod <- make_stacked_bar(
  df = bin4_prod,
  x_col = "Year_4yr_bin",
  title = "Four-year O157:H7 cattle-associated counts by production system",
  subtitle = paste0("Years <=2000 collapsed; final bin ends at ", final_year, "; n = ", comma(nrow(valid_data)), " assemblies"),
  x_label = "Collection-year interval",
  y_label = "Counts",
  show_total_labels = TRUE,
  rotate_x = TRUE
)

save_plot_all(
  p_4yr_prod,
  file.path(fig_4yr_dir, "02_four_year_counts_by_production_system"),
  width = 9,
  height = 5.5
)

p_yearly_stage <- make_stacked_bar(
  df = yearly_stage,
  x_col = "Year_plot",
  title = "Annual counts by production system and supply-chain stage",
  subtitle = "Faceted by stage; years <=2000 collapsed",
  x_label = "Collection year",
  y_label = "Counts",
  show_total_labels = FALSE,
  facet_col = "Supply_chain_stage",
  rotate_x = TRUE
)

save_plot_all(
  p_yearly_stage,
  file.path(fig_yearly_dir, "03_yearly_counts_by_production_system_and_stage"),
  width = 11,
  height = 8.5
)

p_4yr_stage <- make_stacked_bar(
  df = bin4_stage,
  x_col = "Year_4yr_bin",
  title = "Four-year counts by production system and supply-chain stage",
  subtitle = paste0("Faceted by stage; final bin ends at ", final_year),
  x_label = "Collection-year interval",
  y_label = "Counts",
  show_total_labels = TRUE,
  facet_col = "Supply_chain_stage",
  rotate_x = TRUE
)

save_plot_all(
  p_4yr_stage,
  file.path(fig_4yr_dir, "04_four_year_counts_by_production_system_and_stage"),
  width = 10,
  height = 8.5
)

p_pre_yearly <- make_stacked_bar(
  df = preharvest_yearly,
  x_col = "Year_plot",
  title = "Annual pre-harvest counts by production system",
  subtitle = "Filtered to Supply_chain_stage = Pre-harvest; years <=2000 collapsed",
  x_label = "Collection year",
  y_label = "Counts",
  show_total_labels = FALSE,
  rotate_x = TRUE
)

save_plot_all(
  p_pre_yearly,
  file.path(fig_focus_dir, "05_preharvest_yearly_counts_by_production_system"),
  width = 11,
  height = 5.5
)

p_pre_4yr <- make_stacked_bar(
  df = preharvest_bin4,
  x_col = "Year_4yr_bin",
  title = "Four-year pre-harvest counts by production system",
  subtitle = paste0("Filtered to Supply_chain_stage = Pre-harvest; final bin ends at ", final_year),
  x_label = "Collection-year interval",
  y_label = "Counts",
  show_total_labels = TRUE,
  rotate_x = TRUE
)

save_plot_all(
  p_pre_4yr,
  file.path(fig_focus_dir, "06_preharvest_four_year_counts_by_production_system"),
  width = 9,
  height = 5.5
)

p_pre_focus_yearly <- make_stacked_bar(
  df = preharvest_focus_yearly,
  x_col = "Year_plot",
  title = "Annual pre-harvest counts: Beef and cattle-unspecified",
  subtitle = "Filtered to Supply_chain_stage = Pre-harvest; Dairy excluded from this panel",
  x_label = "Collection year",
  y_label = "Counts",
  show_total_labels = FALSE,
  rotate_x = TRUE,
  keep_all_production_levels = FALSE
)

save_plot_all(
  p_pre_focus_yearly,
  file.path(fig_focus_dir, "07_preharvest_beef_cattle_unspecified_yearly_counts"),
  width = 11,
  height = 5.5
)

p_pre_focus_4yr <- make_stacked_bar(
  df = preharvest_focus_bin4,
  x_col = "Year_4yr_bin",
  title = "Four-year pre-harvest counts: Beef and cattle-unspecified",
  subtitle = paste0("Filtered to Supply_chain_stage = Pre-harvest; final bin ends at ", final_year),
  x_label = "Collection-year interval",
  y_label = "Counts",
  show_total_labels = TRUE,
  rotate_x = TRUE,
  keep_all_production_levels = FALSE
)

save_plot_all(
  p_pre_focus_4yr,
  file.path(fig_focus_dir, "08_preharvest_beef_cattle_unspecified_four_year_counts"),
  width = 9,
  height = 5.5
)

# ============================================================
# 7. Pie figures
# ============================================================

p_pie_production <- make_pie_plot(
  df = pie_total_production,
  category_col = "Production_system",
  title = "Total distribution by production system",
  subtitle = "All valid temporal assemblies",
  legend_title = "Production system"
)

save_plot_all(
  p_pie_production,
  file.path(fig_pie_dir, "09_pie_total_by_production_system"),
  width = 7.5,
  height = 5.5
)

p_pie_stage <- make_pie_plot(
  df = pie_total_stage,
  category_col = "Supply_chain_stage",
  title = "Total distribution by supply-chain stage",
  subtitle = "All valid temporal assemblies",
  legend_title = "Supply-chain stage"
)

save_plot_all(
  p_pie_stage,
  file.path(fig_pie_dir, "10_pie_total_by_supply_chain_stage"),
  width = 7.5,
  height = 5.5
)

p_pie_production_stage <- make_pie_plot(
  df = pie_total_production_stage,
  category_col = "Production_stage_group",
  title = "Total distribution by production-stage group",
  subtitle = "Production system combined with supply-chain stage",
  legend_title = "Production-stage group",
  min_label_percent = 3
)

save_plot_all(
  p_pie_production_stage,
  file.path(fig_pie_dir, "11_pie_total_by_production_stage_group"),
  width = 9,
  height = 6.5
)

p_pie_preharvest_production <- make_pie_plot(
  df = pie_preharvest_production,
  category_col = "Production_system",
  title = "Pre-harvest distribution by production system",
  subtitle = "Filtered to Supply_chain_stage = Pre-harvest",
  legend_title = "Production system"
)

save_plot_all(
  p_pie_preharvest_production,
  file.path(fig_pie_dir, "12_pie_preharvest_by_production_system"),
  width = 7.5,
  height = 5.5
)

p_pie_preharvest_focus <- make_pie_plot(
  df = pie_preharvest_focus,
  category_col = "Production_system",
  title = "Pre-harvest distribution: Beef and cattle-unspecified",
  subtitle = "Filtered to Pre-harvest; Dairy excluded from this panel",
  legend_title = "Production system"
)

save_plot_all(
  p_pie_preharvest_focus,
  file.path(fig_pie_dir, "13_pie_preharvest_beef_cattle_unspecified"),
  width = 7.5,
  height = 5.5
)

if (nrow(pie_matrix_class) > 0) {
  p_pie_matrix <- make_pie_plot(
    df = pie_matrix_class,
    category_col = "Source_matrix_class_clean",
    title = "Total distribution by source matrix class",
    subtitle = "Source-matrix classification from data6a",
    legend_title = "Source matrix class",
    min_label_percent = 3
  )

  save_plot_all(
    p_pie_matrix,
    file.path(fig_pie_dir, "14_pie_total_by_source_matrix_class"),
    width = 9,
    height = 6.5
  )
}

# ============================================================
# 8. Console summary
# ============================================================

message("")
message("=== Year column used ===")
message(year_col)

message("")
message("=== Final year detected ===")
message(final_year)

message("")
message("=== Four-year bins used ===")
message(paste(bin4_levels, collapse = ", "))

message("")
message("=== Valid temporal rows by production system ===")
print(valid_data %>% count(Production_system, sort = TRUE), n = Inf)

message("")
message("=== Valid temporal rows by supply-chain stage ===")
print(valid_data %>% count(Supply_chain_stage, sort = TRUE), n = Inf)

message("")
message("=== Valid temporal rows by production-stage group ===")
print(valid_data %>% count(Production_stage_group, sort = TRUE), n = Inf)

message("")
message("=== Pie totals: production system ===")
print(pie_total_production, n = Inf)

message("")
message("=== Pie totals: supply-chain stage ===")
print(pie_total_stage, n = Inf)

message("")
message("=== Pie totals: production-stage group ===")
print(pie_total_production_stage, n = Inf)

message("")
message("=== Pie totals: pre-harvest production system ===")
print(pie_preharvest_production, n = Inf)

message("")
message("=== Pie totals: pre-harvest Beef + Cattle-unspecified ===")
print(pie_preharvest_focus, n = Inf)

message("")
message("=== Missing/invalid year rows ===")
print(nrow(invalid_years))

message("")
message("[DONE] Script 02E completed successfully.")
message("[OUTPUT] Tables: ", table_dir)
message("[OUTPUT] Figures: ", fig_dir)
message("[OUTPUT] Pie charts: ", fig_pie_dir)
message("[OUTPUT] Chart CSV files: ", chart_csv_dir)
message("[TIME] ", Sys.time())
