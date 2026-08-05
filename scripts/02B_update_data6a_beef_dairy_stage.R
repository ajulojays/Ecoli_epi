#!/usr/bin/env Rscript

# ============================================================
# Script 02B
# Purpose:
#   Update data6.csv into data6a.csv with harmonized:
#     1. Production_system
#     2. Supply_chain_stage
#     3. Production_stage_group
#     4. Source_matrix_class
#     5. Production_system_basis
#     6. Supply_chain_stage_basis
#
# Key classification decisions:
#   - Animal-Calf-Bob Veal        -> Beef | Peri-harvest
#   - Formula-fed Veal            -> Beef | Peri-harvest
#   - Feedlot / finishing cattle  -> Beef | Pre-harvest
#   - Feces/manure/rectal         -> Pre-harvest
#   - Hide/carcass                -> Peri-harvest
#   - Ground beef/trim/raw beef   -> Post-harvest
#   - Milk/cheese/dairy           -> Dairy | Post-harvest unless dairy animal/manure/feces
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

message("=== SCRIPT 02B: UPDATE DATA6A WITH BEEF/DAIRY/STAGE/MATRIX CLASS ===")

# =========================
# Paths
# =========================
data6_file <- "/home/samuelajulo/epi/data6.csv"
data6a_root_out <- "/home/samuelajulo/epi/data6a.csv"

out_dir <- "/home/samuelajulo/epi/Ecoli_epi/results/02B_data6_isolation_source_summary"
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

data6a_repo_out <- file.path(table_dir, "data6a_with_beef_dairy_supply_chain_stage.csv")
mapping_out <- file.path(table_dir, "data6a_source_harmonization_mapping.csv")

production_counts_out <- file.path(table_dir, "data6a_production_system_counts.csv")
stage_counts_out <- file.path(table_dir, "data6a_supply_chain_stage_counts.csv")
matrix_counts_out <- file.path(table_dir, "data6a_source_matrix_class_counts.csv")
combined_counts_out <- file.path(table_dir, "data6a_production_stage_group_counts.csv")
basis_counts_out <- file.path(table_dir, "data6a_classification_basis_counts.csv")
original_counts_out <- file.path(table_dir, "data6_original_isolation_source_counts.csv")

production_plot_pdf <- file.path(fig_dir, "data6a_production_system_barplot.pdf")
production_plot_png <- file.path(fig_dir, "data6a_production_system_barplot.png")

stage_plot_pdf <- file.path(fig_dir, "data6a_supply_chain_stage_barplot.pdf")
stage_plot_png <- file.path(fig_dir, "data6a_supply_chain_stage_barplot.png")

matrix_plot_pdf <- file.path(fig_dir, "data6a_source_matrix_class_barplot.pdf")
matrix_plot_png <- file.path(fig_dir, "data6a_source_matrix_class_barplot.png")

combined_plot_pdf <- file.path(fig_dir, "data6a_production_stage_group_barplot.pdf")
combined_plot_png <- file.path(fig_dir, "data6a_production_stage_group_barplot.png")

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

collapse_metadata_text <- function(df) {
  metadata_cols <- c(
    "BioProject",
    "Collected by",
    "SRA Center",
    "Organization",
    "Submitter",
    "Center",
    "Lab",
    "Location",
    "Food origin",
    "Host",
    "IFSAC category",
    "Source",
    "source",
    "geo_loc_name",
    "Geo_loc_name"
  )

  metadata_cols <- metadata_cols[metadata_cols %in% names(df)]

  if (length(metadata_cols) == 0) {
    return(rep("", nrow(df)))
  }

  apply(
    df[, metadata_cols, drop = FALSE],
    1,
    function(row) {
      row <- as.character(row)
      row[is.na(row)] <- ""
      paste(row, collapse = " ")
    }
  )
}

assign_source_matrix_class <- function(x) {
  case_when(
    is.na(x) | x == "" ~ "Other/unknown",

    str_detect(x, "\\bfeces\\b|\\bfecal\\b|\\brectal\\b") ~ "Feces/rectal",
    str_detect(x, "\\bmanure\\b|\\bdung\\b|\\bsoil\\b") ~ "Manure/environment",

    str_detect(x, "\\bhide\\b") ~ "Hide",
    str_detect(x, "\\bcarcass\\b") ~ "Carcass",

    str_detect(x, "ground beef|raw ground beef|lean ground beef|beef patty") ~ "Meat/product",
    str_detect(x, "beef trim|trimming|trim\\b|boneless beef|raw beef|hamburger|burger") ~ "Meat/product",
    str_detect(x, "boneless veal|raw veal|veal meat|veal trim|retail veal") ~ "Meat/product",

    str_detect(x, "\\bmilk\\b|\\bcheese\\b") ~ "Milk/dairy product",

    str_detect(x, "bob veal|formula fed veal|formula-fed veal") ~ "Veal animal category",
    str_detect(x, "\\bcalf\\b") ~ "Calf-associated non-fecal",
    str_detect(x, "\\bcow\\b") ~ "Cow-associated non-fecal",
    str_detect(x, "\\bbovine\\b|\\bcattle\\b") ~ "Bovine/cattle non-fecal",

    TRUE ~ "Other/unknown"
  )
}

assign_production_system <- function(source, meta) {
  case_when(
    is.na(source) | source == "" ~ "Cattle-unspecified",

    # Strong metadata-assisted beef rules
    str_detect(meta, "\\bt\\s*m\\s*arthur\\b") & str_detect(source, "\\bhide\\b") ~ "Beef",
    str_detect(meta, "\\busmarc\\b|\\bu\\s*s\\s*marc\\b|meat animal research center|roman l hruska|hruska") ~ "Beef",
    str_detect(meta, "\\busda\\s*fsis\\b|\\bfsis\\b") & str_detect(source, "\\bheifer\\b") ~ "Beef",

    # Feedlot and finishing cattle are beef
    str_detect(source, "\\bfeedlot\\b|finishing cattle|cattle on feed|feed yard|feedyard") ~ "Beef",

    # Veal-associated labels are beef/veal-associated
    str_detect(source, "bob veal|formula fed veal|formula-fed veal|\\bveal\\b") ~ "Beef",

    # Dairy-specific labels
    str_detect(source, "\\bdairy\\b|\\bmilk\\b|\\bcheese\\b") ~ "Dairy",

    # Beef food/product labels
    str_detect(source, "ground beef|raw ground beef|lean ground beef|beef patty") ~ "Beef",
    str_detect(source, "\\bbeef\\b|hamburger|burger|beef trim|trimming|boneless beef|raw beef") ~ "Beef",

    # Beef-associated animal classes
    str_detect(source, "\\bsteer\\b|beef calf") ~ "Beef",

    # Harvest/processing labels usually beef chain
    str_detect(source, "\\bcarcass\\b") ~ "Beef",

    # Ambiguous cattle labels remain unspecified
    str_detect(source, "\\bbovine\\b|\\bcattle\\b|\\bcow\\b|\\bcalf\\b|\\bheifer\\b|\\banimal\\b") ~ "Cattle-unspecified",

    TRUE ~ "Cattle-unspecified"
  )
}

assign_supply_chain_stage <- function(source, meta) {
  case_when(
    is.na(source) | source == "" ~ "Pre-harvest",

    # Veal animal-class labels are harvest-interface / peri-harvest
    str_detect(source, "animal calf bob veal|bob veal|formula fed veal|formula-fed veal") ~ "Peri-harvest",

    # Feedlot and finishing cattle are pre-harvest
    str_detect(source, "\\bfeedlot\\b|finishing cattle|cattle on feed|feed yard|feedyard") ~ "Pre-harvest",

    # Strong metadata-assisted rules
    str_detect(meta, "\\bt\\s*m\\s*arthur\\b") & str_detect(source, "\\bhide\\b") ~ "Peri-harvest",
    str_detect(meta, "\\busmarc\\b|\\bu\\s*s\\s*marc\\b|meat animal research center|roman l hruska|hruska") ~ "Pre-harvest",
    str_detect(meta, "\\busda\\s*fsis\\b|\\bfsis\\b") & str_detect(source, "\\bheifer\\b") ~ "Pre-harvest",

    # Post-harvest food/product labels
    str_detect(source, "ground beef|raw ground beef|lean ground beef|beef patty|hamburger|burger") ~ "Post-harvest",
    str_detect(source, "beef trim|trimming|trim\\b|boneless beef|raw beef") ~ "Post-harvest",
    str_detect(source, "boneless veal|raw veal|veal meat|veal trim|retail veal") ~ "Post-harvest",
    str_detect(source, "\\bmilk\\b|\\bcheese\\b|retail") ~ "Post-harvest",

    # Peri-harvest slaughter/processing interface
    str_detect(source, "\\bhide\\b|\\bcarcass\\b") ~ "Peri-harvest",

    # Pre-harvest matrices / live animal environments
    str_detect(source, "\\bmanure\\b|\\bdung\\b|\\bsoil\\b") ~ "Pre-harvest",
    str_detect(source, "\\bfeces\\b|\\bfecal\\b|\\brectal\\b") ~ "Pre-harvest",
    str_detect(source, "\\bfarm\\b") ~ "Pre-harvest",

    # Generic live-animal labels without harvest/product evidence
    str_detect(source, "\\bbovine\\b|\\bcattle\\b|\\bcalf\\b|\\bcow\\b|\\bsteer\\b|\\bheifer\\b|\\banimal\\b") ~ "Pre-harvest",

    TRUE ~ "Pre-harvest"
  )
}

assign_production_basis <- function(source, meta, prod) {
  case_when(
    str_detect(meta, "\\bt\\s*m\\s*arthur\\b") & str_detect(source, "\\bhide\\b") ~ "org_and_hide_beef_peri_harvest",
    str_detect(meta, "\\busmarc\\b|\\bu\\s*s\\s*marc\\b|meat animal research center|roman l hruska|hruska") ~ "org_prior_USMARC_meat_animal_research",
    str_detect(meta, "\\busda\\s*fsis\\b|\\bfsis\\b") & str_detect(source, "\\bheifer\\b") ~ "org_prior_USDA_FSIS_heifer",
    str_detect(source, "\\bfeedlot\\b|finishing cattle|cattle on feed|feed yard|feedyard") ~ "feedlot_rule",
    str_detect(source, "bob veal|formula fed veal|formula-fed veal|\\bveal\\b") ~ "veal_animal_category_rule",
    prod == "Dairy" & str_detect(source, "\\bdairy\\b|\\bmilk\\b|\\bcheese\\b") ~ "direct_source_text_dairy",
    prod == "Beef" & str_detect(source, "ground beef|raw ground beef|lean ground beef|beef patty|\\bbeef\\b|hamburger|burger|beef trim|trimming|boneless beef|raw beef|\\bsteer\\b|beef calf|\\bcarcass\\b") ~ "direct_source_text_beef",
    prod == "Cattle-unspecified" & str_detect(source, "\\bfeces\\b|\\bfecal\\b|\\brectal\\b|\\bmanure\\b") ~ "unresolved_cattle_fecal_manure",
    prod == "Cattle-unspecified" & str_detect(source, "\\bcalf\\b") ~ "unresolved_calf",
    prod == "Cattle-unspecified" & str_detect(source, "\\bcow\\b") ~ "unresolved_cow",
    prod == "Cattle-unspecified" ~ "unresolved_cattle_unspecified",
    TRUE ~ "unresolved"
  )
}

assign_stage_basis <- function(source, meta, stage) {
  case_when(
    str_detect(source, "animal calf bob veal|bob veal|formula fed veal|formula-fed veal") ~ "veal_harvest_interface_peri_harvest",
    str_detect(meta, "\\bt\\s*m\\s*arthur\\b") & str_detect(source, "\\bhide\\b") ~ "org_and_hide_peri_harvest",
    str_detect(meta, "\\busmarc\\b|\\bu\\s*s\\s*marc\\b|meat animal research center|roman l hruska|hruska") ~ "org_prior_pre_harvest",
    str_detect(meta, "\\busda\\s*fsis\\b|\\bfsis\\b") & str_detect(source, "\\bheifer\\b") ~ "org_prior_pre_harvest",
    str_detect(source, "\\bfeedlot\\b|finishing cattle|cattle on feed|feed yard|feedyard") ~ "feedlot_pre_harvest",
    stage == "Post-harvest" ~ "direct_source_text_post_harvest",
    stage == "Peri-harvest" ~ "direct_source_text_peri_harvest",
    stage == "Pre-harvest" ~ "direct_source_text_pre_harvest",
    TRUE ~ "unresolved"
  )
}

# =========================
# Create data6a
# =========================
metadata_text <- collapse_metadata_text(data6)

data6a <- data6 %>%
  mutate(
    Isolation_source_original = .data[[source_col]],
    Isolation_source_clean = clean_source(.data[[source_col]]),
    Metadata_text_clean = clean_source(metadata_text),
    Source_matrix_class = assign_source_matrix_class(Isolation_source_clean),
    Production_system = assign_production_system(Isolation_source_clean, Metadata_text_clean),
    Supply_chain_stage = assign_supply_chain_stage(Isolation_source_clean, Metadata_text_clean),
    Production_stage_group = paste(Production_system, Supply_chain_stage, sep = " | "),
    Production_system_basis = assign_production_basis(
      Isolation_source_clean,
      Metadata_text_clean,
      Production_system
    ),
    Supply_chain_stage_basis = assign_stage_basis(
      Isolation_source_clean,
      Metadata_text_clean,
      Supply_chain_stage
    ),

    # Backward-compatible older column name
    Isolation_source_stage = case_when(
      Supply_chain_stage == "Post-harvest" ~ "post-harvest",
      Supply_chain_stage %in% c("Pre-harvest", "Peri-harvest") ~ "pre-harvest and harvest",
      TRUE ~ "pre-harvest and harvest"
    )
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
    Source_matrix_class,
    Production_system,
    Supply_chain_stage,
    Production_stage_group,
    Production_system_basis,
    Supply_chain_stage_basis,
    Isolation_source_stage,
    name = "n",
    sort = TRUE
  ) %>%
  mutate(percent = 100 * n / sum(n))

write_csv(mapping_table, mapping_out)

# =========================
# Counts
# =========================
production_counts <- data6a %>%
  count(Production_system, name = "n", sort = TRUE) %>%
  mutate(percent = 100 * n / sum(n))

stage_counts <- data6a %>%
  count(Supply_chain_stage, name = "n", sort = TRUE) %>%
  mutate(percent = 100 * n / sum(n))

matrix_counts <- data6a %>%
  count(Source_matrix_class, name = "n", sort = TRUE) %>%
  mutate(percent = 100 * n / sum(n))

combined_counts <- data6a %>%
  count(Production_stage_group, name = "n", sort = TRUE) %>%
  mutate(percent = 100 * n / sum(n))

basis_counts <- data6a %>%
  count(Production_system_basis, Supply_chain_stage_basis, name = "n", sort = TRUE) %>%
  mutate(percent = 100 * n / sum(n))

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

write_csv(production_counts, production_counts_out)
write_csv(stage_counts, stage_counts_out)
write_csv(matrix_counts, matrix_counts_out)
write_csv(combined_counts, combined_counts_out)
write_csv(basis_counts, basis_counts_out)

write_csv(
  original_counts %>%
    select(Isolation_source = Isolation_source_original, n, percent),
  original_counts_out
)

# =========================
# Plots
# =========================
p_production <- production_counts %>%
  mutate(Production_system = factor(Production_system, levels = rev(Production_system))) %>%
  ggplot(aes(x = Production_system, y = n)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw(base_size = 13) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Production system distribution in data6a",
    subtitle = paste0("n = ", comma(sum(production_counts$n)), " assemblies"),
    x = "Production system",
    y = "Number of assemblies"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(production_plot_pdf, p_production, width = 8, height = 4.5)
ggsave(production_plot_png, p_production, width = 8, height = 4.5, dpi = 300)

p_stage <- stage_counts %>%
  mutate(Supply_chain_stage = factor(Supply_chain_stage, levels = rev(Supply_chain_stage))) %>%
  ggplot(aes(x = Supply_chain_stage, y = n)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw(base_size = 13) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Supply-chain stage distribution in data6a",
    subtitle = paste0("n = ", comma(sum(stage_counts$n)), " assemblies"),
    x = "Supply-chain stage",
    y = "Number of assemblies"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(stage_plot_pdf, p_stage, width = 8, height = 4.5)
ggsave(stage_plot_png, p_stage, width = 8, height = 4.5, dpi = 300)

p_matrix <- matrix_counts %>%
  mutate(Source_matrix_class = factor(Source_matrix_class, levels = rev(Source_matrix_class))) %>%
  ggplot(aes(x = Source_matrix_class, y = n)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw(base_size = 12) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Source matrix class distribution in data6a",
    subtitle = paste0("n = ", comma(sum(matrix_counts$n)), " assemblies"),
    x = "Source matrix class",
    y = "Number of assemblies"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(matrix_plot_pdf, p_matrix, width = 9, height = 5.5)
ggsave(matrix_plot_png, p_matrix, width = 9, height = 5.5, dpi = 300)

p_combined <- combined_counts %>%
  mutate(Production_stage_group = factor(Production_stage_group, levels = rev(Production_stage_group))) %>%
  ggplot(aes(x = Production_stage_group, y = n)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw(base_size = 12) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Production system by supply-chain stage in data6a",
    subtitle = paste0("n = ", comma(sum(combined_counts$n)), " assemblies"),
    x = "Production-stage group",
    y = "Number of assemblies"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(combined_plot_pdf, p_combined, width = 10, height = 6)
ggsave(combined_plot_png, p_combined, width = 10, height = 6, dpi = 300)

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

# =========================
# Console summary
# =========================
message("[DONE] data6a created/updated.")
message("[OUTPUT] Root data6a: ", data6a_root_out)
message("[OUTPUT] Repo data6a: ", data6a_repo_out)
message("[OUTPUT] Mapping table: ", mapping_out)

message("")
message("Production system counts:")
print(production_counts)

message("")
message("Supply-chain stage counts:")
print(stage_counts)

message("")
message("Source matrix class counts:")
print(matrix_counts)

message("")
message("Production-stage group counts:")
print(combined_counts)

message("")
message("Classification basis counts:")
print(basis_counts)

message("")
message("Mapping preview:")
print(head(mapping_table, 50))
