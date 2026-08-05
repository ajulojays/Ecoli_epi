#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 02D
# Fecal/Manure Pre-Harvest O157:H7 Source Summary & Figures
# ============================================================
#
# Target cohorts:
#   1. Beef_Preharvest
#   2. Unspecified_Feces_Rectal_Preharvest
#   3. Unspecified_Manure_Environment_Preharvest
#
# Key logic:
#   - data6a is already O157:H7-focused.
#   - This script therefore summarizes source/matrix/temporal patterns,
#     not serotype diversity.
#   - Manure/feces/rectal are retained as risk-relevant matrices.
#   - Generic bovine/calf/cow without fecal/manure evidence are excluded
#     from the unspecified fecal/manure comparison.
#   - Explicit dairy signal only: dairy / milk / cheese / dairy cow.
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

message("=== SCRIPT 02D STARTED ===")
message("[TIME] ", Sys.time())

# ============================================================
# 1. Paths
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

fig_o157_dir <- file.path(figure_dir, "O157H7_focus_pies")
fig_temporal_count_dir <- file.path(figure_dir, "temporal_stacked_counts")
fig_temporal_percent_dir <- file.path(figure_dir, "temporal_stacked_percent")
fig_original_source_dir <- file.path(figure_dir, "original_source_barplots")
fig_multipage_dir <- file.path(figure_dir, "multipage_pdfs")

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_o157_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_temporal_count_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_temporal_percent_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_original_source_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_multipage_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(METADATA_CSV)) {
  stop("[ERROR] Metadata file not found: ", METADATA_CSV, call. = FALSE)
}

if (!file.exists(ectyper_file)) {
  stop("[ERROR] ECtyper output file not found: ", ectyper_file, call. = FALSE)
}

message("[INFO] Metadata file: ", METADATA_CSV)
message("[INFO] ECtyper file: ", ectyper_file)
message("[INFO] Output directory: ", script02d_dir)

# ============================================================
# 2. Helper functions
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

normalize_tip_id <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  x <- str_remove(x, "\\.gz$")
  x <- str_remove(x, "\\.(fasta|fna|fa|gff|gff3|txt|tsv|csv)$")
  accession_like <- str_extract(x, "GC[AF]_[0-9]+\\.[0-9]+")
  ifelse(is.na(accession_like), x, accession_like)
}

collapse_available_cols <- function(df, cols) {
  cols <- cols[cols %in% names(df)]

  if (length(cols) == 0) {
    return(rep("", nrow(df)))
  }

  apply(
    df[, cols, drop = FALSE],
    1,
    function(row) {
      row <- as.character(row)
      row[is.na(row)] <- ""
      paste(row, collapse = " ")
    }
  )
}

clean_text_for_filter <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[_\\-/]+", " ") %>%
    str_replace_all("[^a-z0-9 ]+", " ") %>%
    str_squish()
}

find_col <- function(df, candidates, label) {
  exact <- candidates[candidates %in% names(df)]
  if (length(exact) > 0) {
    return(exact[1])
  }

  clean_names <- names(df) %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "")

  clean_candidates <- candidates %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "")

  idx <- match(clean_candidates, clean_names)
  idx <- idx[!is.na(idx)]

  if (length(idx) > 0) {
    return(names(df)[idx[1]])
  }

  stop(
    "[ERROR] Could not find ECtyper column for: ", label,
    "\nTried: ", paste(candidates, collapse = ", "),
    "\nAvailable ECtyper columns:\n",
    paste(names(df), collapse = "\n"),
    call. = FALSE
  )
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
# 3. Load inputs
# ============================================================

metadata <- read_csv(METADATA_CSV, show_col_types = FALSE)
ectyper_raw <- read_tsv(ectyper_file, show_col_types = FALSE)

message("")
message("=== ECtyper columns detected ===")
print(names(ectyper_raw))

name_col <- find_col(
  ectyper_raw,
  c("Name", "name", "Assembly", "assembly", "Filename", "filename", "Sample", "sample"),
  "Name / assembly ID"
)

o_col <- find_col(
  ectyper_raw,
  c("O-type", "O_type", "O type", "Otype", "O", "O antigen", "O-antigen", "O_antigen"),
  "O-type"
)

h_col <- find_col(
  ectyper_raw,
  c("H-type", "H_type", "H type", "Htype", "H", "H antigen", "H-antigen", "H_antigen"),
  "H-type"
)

serotype_col <- find_col(
  ectyper_raw,
  c("Serotype", "serotype", "Serotype_clean", "Predicted serotype", "Predicted_serotype"),
  "Serotype"
)

message("[INFO] ECtyper Name column: ", name_col)
message("[INFO] ECtyper O column: ", o_col)
message("[INFO] ECtyper H column: ", h_col)
message("[INFO] ECtyper Serotype column: ", serotype_col)

ectyper <- ectyper_raw %>%
  transmute(
    ectyper_Name = str_squish(as.character(.data[[name_col]])),
    Assembly_key = normalize_tip_id(ectyper_Name),
    ectyper_O_type = .data[[o_col]],
    ectyper_H_type = .data[[h_col]],
    ectyper_Serotype = .data[[serotype_col]]
  ) %>%
  filter(!is.na(Assembly_key), Assembly_key != "") %>%
  distinct(Assembly_key, .keep_all = TRUE)

required_metadata_cols <- c("Assembly", "Production_system", "Supply_chain_stage")
missing_metadata_cols <- setdiff(required_metadata_cols, names(metadata))

if (length(missing_metadata_cols) > 0) {
  stop(
    "[ERROR] Missing required metadata column(s): ",
    paste(missing_metadata_cols, collapse = ", "),
    "\nAvailable metadata columns:\n",
    paste(names(metadata), collapse = "\n"),
    call. = FALSE
  )
}

metadata <- metadata %>%
  mutate(
    Assembly = str_squish(as.character(Assembly)),
    Assembly_key = normalize_tip_id(Assembly),
    Production_system = str_squish(as.character(Production_system)),
    Supply_chain_stage = str_squish(as.character(Supply_chain_stage))
  ) %>%
  filter(!is.na(Assembly_key), Assembly_key != "") %>%
  distinct(Assembly_key, .keep_all = TRUE)

message("[INFO] Metadata rows after Assembly_key deduplication: ", nrow(metadata))
message("[INFO] ECtyper rows after Assembly_key deduplication: ", nrow(ectyper))

# ============================================================
# 4. Matrix-aware targeted filtering
# ============================================================

source_matrix_cols <- c(
  "Isolation source",
  "Isolation_source_original",
  "Isolation_source_clean",
  "Food origin",
  "food_origin",
  "IFSAC category",
  "IFSAC_category",
  "Source",
  "source"
)

context_cols <- c(
  "Isolation source",
  "Isolation_source_original",
  "Isolation_source_clean",
  "Host",
  "host",
  "Food origin",
  "food_origin",
  "IFSAC category",
  "IFSAC_category",
  "Source",
  "source",
  "Location",
  "location",
  "geo_loc_name",
  "Geo_loc_name",
  "Collected by",
  "collected_by",
  "SRA Center",
  "SRA_Center",
  "BioProject",
  "Submitter",
  "Organization",
  "Center",
  "Lab"
)

metadata <- metadata %>%
  mutate(
    Isolation_source_for_audit = case_when(
      "Isolation source" %in% names(.) ~ as.character(.data[["Isolation source"]]),
      "Isolation_source_original" %in% names(.) ~ as.character(.data[["Isolation_source_original"]]),
      TRUE ~ NA_character_
    ),

    source_matrix_text = clean_text_for_filter(
      collapse_available_cols(., source_matrix_cols)
    ),

    context_text = clean_text_for_filter(
      collapse_available_cols(., context_cols)
    ),

    Source_matrix_class_02D = case_when(
      str_detect(source_matrix_text, "\\bfeces\\b|\\bfecal\\b|\\brectal\\b") ~ "Feces/rectal",
      str_detect(source_matrix_text, "\\bmanure\\b|\\bdung\\b|\\bsoil\\b") ~ "Manure/environment",
      str_detect(source_matrix_text, "\\bhide\\b") ~ "Hide",
      str_detect(source_matrix_text, "\\bcarcass\\b") ~ "Carcass",
      str_detect(source_matrix_text, "ground beef|raw ground beef|lean ground beef|beef patty|beef trim|trimming|boneless beef|raw beef|hamburger|burger") ~ "Meat/product",
      str_detect(source_matrix_text, "\\bmilk\\b|\\bcheese\\b") ~ "Milk/dairy product",
      str_detect(source_matrix_text, "\\bcalf\\b") ~ "Calf-associated non-fecal",
      str_detect(source_matrix_text, "\\bcow\\b") ~ "Cow-associated non-fecal",
      str_detect(source_matrix_text, "\\bbovine\\b|\\bcattle\\b") ~ "Bovine/cattle non-fecal",
      TRUE ~ "Other/unknown"
    ),

    is_fecal_manure_matrix = Source_matrix_class_02D %in% c(
      "Feces/rectal",
      "Manure/environment"
    ),

    has_explicit_dairy_signal = str_detect(
      context_text,
      "\\bdairy\\b|\\bmilk\\b|\\bcheese\\b|dairy cow"
    )
  )

unspecified_preharvest_audit <- metadata %>%
  filter(
    Production_system == "Cattle-unspecified",
    Supply_chain_stage == "Pre-harvest"
  ) %>%
  count(
    Isolation_source_for_audit,
    Source_matrix_class_02D,
    is_fecal_manure_matrix,
    has_explicit_dairy_signal,
    name = "n",
    sort = TRUE
  )

write_csv(
  unspecified_preharvest_audit,
  file.path(table_dir, "02D_unspecified_preharvest_matrix_audit_all.csv")
)

filtered_metadata <- metadata %>%
  filter(
    (
      Production_system == "Beef" &
        Supply_chain_stage == "Pre-harvest"
    ) |
    (
      Production_system == "Cattle-unspecified" &
        Supply_chain_stage == "Pre-harvest" &
        is_fecal_manure_matrix &
        !has_explicit_dairy_signal
    )
  ) %>%
  mutate(
    cohort_group = case_when(
      Production_system == "Beef" ~ "Beef_Preharvest",
      Production_system == "Cattle-unspecified" &
        Source_matrix_class_02D == "Feces/rectal" ~ "Unspecified_Feces_Rectal_Preharvest",
      Production_system == "Cattle-unspecified" &
        Source_matrix_class_02D == "Manure/environment" ~ "Unspecified_Manure_Environment_Preharvest",
      TRUE ~ "Other"
    ),
    Production_system_basis_02D = if ("Production_system_basis" %in% names(.)) {
      as.character(Production_system_basis)
    } else {
      NA_character_
    },
    Supply_chain_stage_basis_02D = if ("Supply_chain_stage_basis" %in% names(.)) {
      as.character(Supply_chain_stage_basis)
    } else {
      NA_character_
    }
  )

excluded_unspecified_preharvest <- metadata %>%
  filter(
    Production_system == "Cattle-unspecified",
    Supply_chain_stage == "Pre-harvest",
    !(
      is_fecal_manure_matrix &
        !has_explicit_dairy_signal
    )
  ) %>%
  count(
    Isolation_source_for_audit,
    Source_matrix_class_02D,
    is_fecal_manure_matrix,
    has_explicit_dairy_signal,
    name = "n",
    sort = TRUE
  )

included_unspecified_preharvest <- filtered_metadata %>%
  filter(Production_system == "Cattle-unspecified") %>%
  count(
    cohort_group,
    Isolation_source_for_audit,
    Source_matrix_class_02D,
    name = "n",
    sort = TRUE
  )

write_csv(
  excluded_unspecified_preharvest,
  file.path(table_dir, "02D_excluded_unspecified_preharvest_non_fecal_non_manure.csv")
)

write_csv(
  included_unspecified_preharvest,
  file.path(table_dir, "02D_included_unspecified_fecal_manure_preharvest.csv")
)

write_csv(
  filtered_metadata,
  file.path(table_dir, "fecal_preharvest_filtered_metadata_before_ectyper_join.csv")
)

filter_audit <- filtered_metadata %>%
  count(
    cohort_group,
    Production_system,
    Supply_chain_stage,
    Production_system_basis_02D,
    Supply_chain_stage_basis_02D,
    Source_matrix_class_02D,
    Isolation_source_for_audit,
    is_fecal_manure_matrix,
    has_explicit_dairy_signal,
    name = "n",
    sort = TRUE
  )

filter_summary <- filtered_metadata %>%
  count(
    cohort_group,
    Production_system,
    Supply_chain_stage,
    Source_matrix_class_02D,
    name = "n",
    sort = TRUE
  )

write_csv(
  filter_audit,
  file.path(table_dir, "fecal_preharvest_filter_audit_table.csv")
)

write_csv(
  filter_summary,
  file.path(table_dir, "fecal_preharvest_filter_summary.csv")
)

# ============================================================
# 5. Join ECtyper and clean calls
# ============================================================

joined <- ectyper %>%
  inner_join(filtered_metadata, by = "Assembly_key") %>%
  mutate(
    source_group_clean = cohort_group,

    O_type_clean = clean_call(ectyper_O_type),
    H_type_clean = clean_call(ectyper_H_type),

    Serotype_clean = make_clean_serotype(
      serotype = ectyper_Serotype,
      otype = ectyper_O_type,
      htype = ectyper_H_type
    ),

    serotype_call_status = case_when(
      Serotype_clean == "No_serotype_call" ~ "no_O_or_H_call",
      str_detect(Serotype_clean, "No_H_call") ~ "partial_O_only",
      str_detect(Serotype_clean, "No_O_call") ~ "partial_H_only",
      TRUE ~ "complete_or_reconstructed_call"
    )
  )

message("[INFO] Joined filtered rows with ECtyper: ", nrow(joined))

write_csv(
  joined,
  file.path(table_dir, "fecal_preharvest_joined_metadata_ectyper.csv")
)

join_check <- tibble(
  metric = c(
    "n_metadata_total",
    "n_filtered_metadata",
    "n_ectyper_total",
    "n_joined_filtered_ectyper",
    "n_filtered_without_ectyper"
  ),
  value = c(
    nrow(metadata),
    nrow(filtered_metadata),
    nrow(ectyper),
    nrow(joined),
    nrow(filtered_metadata) - nrow(joined)
  )
)

write_csv(
  join_check,
  file.path(table_dir, "fecal_preharvest_join_check.csv")
)

# ============================================================
# 6. Temporal year handling
# ============================================================

if ("Collection_year_temporal_group" %in% names(joined)) {
  joined <- joined %>%
    mutate(
      Collection_year_temporal_group_existing = as.character(Collection_year_temporal_group)
    )
}

year_candidates <- c(
  "Collection_year",
  "Collection year",
  "collection_year",
  "Year",
  "year",
  "Collection date",
  "collection_date",
  "Date",
  "date"
)

year_col <- year_candidates[year_candidates %in% names(joined)][1]

if (!is.na(year_col)) {
  joined <- joined %>%
    mutate(
      Collection_year_clean = clean_year(.data[[year_col]]),
      Collection_year_temporal_group = make_temporal_year_group(Collection_year_clean)
    )
} else if ("Collection_year_temporal_group_existing" %in% names(joined)) {
  joined <- joined %>%
    mutate(
      Collection_year_clean = clean_year(Collection_year_temporal_group_existing),
      Collection_year_temporal_group = Collection_year_temporal_group_existing
    )
} else {
  joined <- joined %>%
    mutate(
      Collection_year_clean = NA_integer_,
      Collection_year_temporal_group = NA_character_
    )
}

write_csv(
  joined,
  file.path(table_dir, "fecal_preharvest_joined_metadata_ectyper.csv")
)

# ============================================================
# 7. Summary tables
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

write_csv(
  summary_by_source,
  file.path(table_dir, "fecal_preharvest_summary_by_cohort.csv")
)

serotype_counts <- joined %>%
  count(source_group_clean, Source_matrix_class_02D, Serotype_clean, name = "n", sort = TRUE) %>%
  group_by(source_group_clean) %>%
  mutate(
    source_group_total = sum(n),
    percent_within_source_group = round(100 * n / source_group_total, 2)
  ) %>%
  ungroup()

write_csv(
  serotype_counts,
  file.path(table_dir, "fecal_preharvest_all_serotype_counts.csv")
)

plot_group_counts <- joined %>%
  mutate(
    Serotype_plot_group = case_when(
      Serotype_clean == "O157:H7" ~ "O157:H7",
      Serotype_clean == "No_serotype_call" ~ "No serotype call",
      TRUE ~ "Other serotype"
    )
  ) %>%
  count(source_group_clean, Serotype_plot_group, name = "n") %>%
  group_by(source_group_clean) %>%
  mutate(
    source_group_total = sum(n),
    percent_within_source_group = round(100 * n / source_group_total, 2)
  ) %>%
  ungroup()

write_csv(
  plot_group_counts,
  file.path(table_dir, "fecal_preharvest_plot_groups.csv")
)

# ============================================================
# 8. Original source barplot tables and figures
# ============================================================

source_counts <- joined %>%
  mutate(
    Original_source = str_squish(as.character(Isolation_source_for_audit)),
    Original_source = if_else(
      is.na(Original_source) | Original_source == "",
      "Missing/blank",
      Original_source
    )
  ) %>%
  count(
    source_group_clean,
    Source_matrix_class_02D,
    Serotype_clean,
    Original_source,
    name = "n",
    sort = TRUE
  ) %>%
  group_by(source_group_clean, Source_matrix_class_02D) %>%
  mutate(
    group_total = sum(n),
    percent_within_group = round(100 * n / group_total, 2)
  ) %>%
  ungroup() %>%
  arrange(source_group_clean, Source_matrix_class_02D, desc(n), Original_source)

overall_source_counts <- joined %>%
  mutate(
    Original_source = str_squish(as.character(Isolation_source_for_audit)),
    Original_source = if_else(
      is.na(Original_source) | Original_source == "",
      "Missing/blank",
      Original_source
    )
  ) %>%
  count(Original_source, name = "n", sort = TRUE) %>%
  mutate(
    total = sum(n),
    percent_overall = round(100 * n / total, 2)
  )

cohort_matrix_counts <- joined %>%
  count(
    source_group_clean,
    Source_matrix_class_02D,
    name = "n",
    sort = TRUE
  ) %>%
  mutate(
    total = sum(n),
    percent_overall = round(100 * n / total, 2)
  )

write_csv(
  source_counts,
  file.path(table_dir, "02D_included_original_source_counts_by_cohort_matrix.csv")
)

write_csv(
  overall_source_counts,
  file.path(table_dir, "02D_included_original_source_counts_overall.csv")
)

write_csv(
  cohort_matrix_counts,
  file.path(table_dir, "02D_included_cohort_matrix_counts.csv")
)

p_overall_sources <- overall_source_counts %>%
  mutate(
    Original_source = fct_reorder(Original_source, n)
  ) %>%
  ggplot(aes(x = n, y = Original_source)) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = paste0("n=", n, " (", percent_overall, "%)")),
    hjust = -0.05,
    size = 3.3
  ) +
  scale_x_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Original source distribution among included 02D records",
    subtitle = paste0("n = ", comma(sum(overall_source_counts$n)), " records"),
    x = "Number of records",
    y = "Original source"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(
  file.path(fig_original_source_dir, "02D_original_source_barplot_overall.pdf"),
  p_overall_sources,
  width = 10,
  height = 7
)

ggsave(
  file.path(fig_original_source_dir, "02D_original_source_barplot_overall.png"),
  p_overall_sources,
  width = 10,
  height = 7,
  dpi = 600
)

source_counts_plot <- source_counts %>%
  mutate(
    cohort_matrix_group = paste(
      source_group_clean,
      Source_matrix_class_02D,
      sep = " | "
    ),
    Original_source_plot = paste0(Original_source, "   [", cohort_matrix_group, "]"),
    Original_source_plot = fct_reorder(Original_source_plot, n)
  )

p_by_group_sources <- ggplot(
  source_counts_plot,
  aes(x = n, y = Original_source_plot)
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = paste0("n=", n, " (", percent_within_group, "%)")),
    hjust = -0.05,
    size = 3.0
  ) +
  scale_x_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "Original sources within each included 02D cohort/matrix group",
    subtitle = paste0("n = ", comma(sum(source_counts$n)), " records"),
    x = "Number of records",
    y = "Original source [cohort | matrix class]"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(
  file.path(fig_original_source_dir, "02D_original_source_barplot_by_cohort_matrix.pdf"),
  p_by_group_sources,
  width = 14,
  height = 9
)

ggsave(
  file.path(fig_original_source_dir, "02D_original_source_barplot_by_cohort_matrix.png"),
  p_by_group_sources,
  width = 14,
  height = 9,
  dpi = 600
)

p_facet_sources <- source_counts %>%
  mutate(
    cohort_matrix_group = paste(
      source_group_clean,
      Source_matrix_class_02D,
      sep = " | "
    ),
    Original_source = fct_reorder(Original_source, n)
  ) %>%
  ggplot(aes(x = n, y = Original_source)) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = paste0("n=", n)),
    hjust = -0.05,
    size = 3.0
  ) +
  facet_wrap(~ cohort_matrix_group, scales = "free_y", ncol = 1) +
  scale_x_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Original source distribution by included 02D cohort/matrix group",
    subtitle = paste0("n = ", comma(sum(source_counts$n)), " records"),
    x = "Number of records",
    y = "Original source"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(
  file.path(fig_original_source_dir, "02D_original_source_barplot_faceted_by_cohort_matrix.pdf"),
  p_facet_sources,
  width = 12,
  height = 12
)

ggsave(
  file.path(fig_original_source_dir, "02D_original_source_barplot_faceted_by_cohort_matrix.png"),
  p_facet_sources,
  width = 12,
  height = 12,
  dpi = 600
)

# ============================================================
# 9. O157:H7 focus pie charts
# ============================================================

o157_focus_counts <- joined %>%
  mutate(
    O157_focus_group = case_when(
      Serotype_clean == "O157:H7" ~ "O157:H7",
      Serotype_clean == "No_serotype_call" ~ "No serotype call",
      TRUE ~ "Non-O157:H7"
    )
  ) %>%
  count(source_group_clean, O157_focus_group, name = "n") %>%
  group_by(source_group_clean) %>%
  mutate(
    source_group_total = sum(n),
    percent_within_source_group = round(100 * n / source_group_total, 2)
  ) %>%
  ungroup()

write_csv(
  o157_focus_counts,
  file.path(table_dir, "fecal_preharvest_O157H7_focus.csv")
)

make_o157_pie_plot <- function(df_source, source_name) {
  df_plot <- df_source %>%
    arrange(desc(n), O157_focus_group) %>%
    mutate(
      O157_focus_group = fct_reorder(O157_focus_group, n),
      pct = 100 * n / sum(n),
      slice_label = paste0(round(pct, 1), "%\n", "n=", n)
    )

  ggplot(df_plot, aes(x = "", y = n, fill = O157_focus_group)) +
    geom_col(width = 1, color = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = slice_label),
      position = position_stack(vjust = 0.5),
      size = 3.4,
      color = "black"
    ) +
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

source_groups <- sort(unique(plot_group_counts$source_group_clean))

for (src in source_groups) {
  df_src <- o157_focus_counts %>% filter(source_group_clean == src)

  if (nrow(df_src) > 0) {
    p <- make_o157_pie_plot(df_src, src)
    base <- file.path(fig_o157_dir, paste0("pie_O157H7_focus_", safe_filename(src)))

    ggsave(paste0(base, ".pdf"), p, width = 7, height = 6)
    ggsave(paste0(base, ".png"), p, width = 7, height = 6, dpi = 600)
  }
}

# ============================================================
# 10. Temporal stacked bar plots
# ============================================================

temporal_data <- joined %>%
  mutate(
    Serotype_plot_group = case_when(
      Serotype_clean == "O157:H7" ~ "O157:H7",
      Serotype_clean == "No_serotype_call" ~ "No serotype call",
      TRUE ~ "Other serotype"
    )
  ) %>%
  filter(
    !is.na(Collection_year_clean),
    Collection_year_clean >= 1900,
    Collection_year_clean <= 2100
  ) %>%
  mutate(
    Collection_year_temporal_group = make_temporal_year_group(Collection_year_clean)
  )

if (nrow(temporal_data) > 0) {
  temporal_counts <- temporal_data %>%
    count(
      source_group_clean,
      Collection_year_temporal_group,
      Serotype_plot_group,
      name = "n"
    ) %>%
    group_by(source_group_clean, Collection_year_temporal_group) %>%
    mutate(
      year_source_total = sum(n),
      prop_within_year_source = n / year_source_total,
      percent_within_year_source = 100 * prop_within_year_source
    ) %>%
    ungroup() %>%
    mutate(
      Collection_year_temporal_group = make_temporal_year_factor(Collection_year_temporal_group)
    )

  write_csv(
    temporal_counts,
    file.path(table_dir, "fecal_preharvest_temporal_counts.csv")
  )

  make_temporal_count_plot <- function(df_source, source_name) {
    total_df <- df_source %>%
      distinct(Collection_year_temporal_group, year_source_total)

    ggplot(df_source, aes(x = Collection_year_temporal_group, y = n, fill = Serotype_plot_group)) +
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
        title = paste0("Temporal Distribution: ", source_name),
        subtitle = "Years <=2000 grouped",
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
      distinct(Collection_year_temporal_group, year_source_total)

    ggplot(df_source, aes(x = Collection_year_temporal_group, y = percent_within_year_source, fill = Serotype_plot_group)) +
      geom_col(width = 0.85, color = "gray20", linewidth = 0.1) +
      geom_text(
        data = total_df,
        aes(
          x = Collection_year_temporal_group,
          y = 100,
          label = paste0("n=", year_source_total)
        ),
        inherit.aes = FALSE,
        angle = 90,
        vjust = -0.35,
        hjust = 0.5,
        size = 3.0
      ) +
      scale_y_continuous(
        labels = function(x) paste0(x, "%"),
        limits = c(0, 112),
        breaks = seq(0, 100, 20)
      ) +
      coord_cartesian(clip = "off") +
      labs(
        title = paste0("Temporal Percentage Distribution: ", source_name),
        subtitle = "Years <=2000 grouped",
        x = "Collection year",
        y = "Percent within year",
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

    p_percent <- make_temporal_percent_plot(df_src, src)
    base_percent <- file.path(fig_temporal_percent_dir, paste0("temporal_percent_", safe_filename(src)))

    ggsave(paste0(base_percent, ".pdf"), p_percent, width = 13, height = 7)
    ggsave(paste0(base_percent, ".png"), p_percent, width = 13, height = 7, dpi = 600)
  }
} else {
  message("[WARN] No valid temporal data available for temporal plots.")
}

# ============================================================
# 11. Multipage PDF summary
# ============================================================

multipage_pdf <- file.path(fig_multipage_dir, "fecal_preharvest_source_summary_all_plots.pdf")

pdf(multipage_pdf, width = 11, height = 8.5)

print(p_overall_sources)
print(p_by_group_sources)
print(p_facet_sources)

if (length(source_groups) > 0) {
  for (src in source_groups) {
    df_src_o157 <- o157_focus_counts %>% filter(source_group_clean == src)

    if (nrow(df_src_o157) > 0) {
      print(make_o157_pie_plot(df_src_o157, src))
    }
  }
} else {
  plot.new()
  text(0.5, 0.5, "No filtered data available for plotting")
}

dev.off()

# ============================================================
# 12. Final summary
# ============================================================

message("")
message("=== Summary by Cohort Group ===")
print(summary_by_source)

message("")
message("=== Filter Summary ===")
print(filter_summary, n = Inf)

message("")
message("=== Join Check ===")
print(join_check)

message("")
message("=== Included unspecified rows ===")
print(included_unspecified_preharvest, n = Inf)

message("")
message("=== Excluded unspecified rows ===")
print(excluded_unspecified_preharvest, n = Inf)

message("")
message("=== Original source counts by cohort/matrix ===")
print(source_counts, n = Inf)

message("")
message("[DONE] Script 02D completed successfully.")
message("[OUTPUT] Tables saved to: ", table_dir)
message("[OUTPUT] Figures saved to: ", figure_dir)
message("[OUTPUT] Original source plots saved to: ", fig_original_source_dir)
message("[OUTPUT] Multipage PDF: ", multipage_pdf)
message("[TIME] ", Sys.time())
