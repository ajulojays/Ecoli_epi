#!/usr/bin/env Rscript

# ============================================================
# SCRIPT 01: Summarize ECtyper serotypes by source group
# ============================================================
#
# Purpose:
#   This script joins ECtyper output to curated metadata and summarizes
#   predicted serotypes by source group.
#
# Main inputs:
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
#   ECtyper Name  =  metadata Assembly
#
# Main outputs:
#   current_ectyper_joined_metadata.csv
#   current_serotype_summary_by_source_group.csv
#   current_serotype_counts_by_source_group.csv
#   current_Otype_counts_by_source_group.csv
#   current_Htype_counts_by_source_group.csv
#
# Missing serotype logic:
#   Some genomes may not have complete ECtyper calls.
#   These are retained and labeled explicitly:
#      No_serotype_call
#      O8:No_H_call
#      No_O_call:H7
#
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

# Allow users to override the large working directory.
# Default used in the original project:
#   ~/epi/marker_screen
wd <- path.expand(Sys.getenv("ECOLI_EPI_WORKDIR", unset = "~/epi/marker_screen"))

metadata_file <- file.path(wd, "data5.csv")
ectyper_file  <- file.path(wd, "ectyper_results", "ectyper_output_all.tsv")

out_dir <- file.path(wd, "ectyper_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_joined <- file.path(out_dir, "current_ectyper_joined_metadata.csv")
out_summary <- file.path(out_dir, "current_serotype_summary_by_source_group.csv")
out_serotype_counts <- file.path(out_dir, "current_serotype_counts_by_source_group.csv")
out_otype_counts <- file.path(out_dir, "current_Otype_counts_by_source_group.csv")
out_htype_counts <- file.path(out_dir, "current_Htype_counts_by_source_group.csv")

message("=== ECtyper serotype summary started: ", Sys.time(), " ===")
message("[INFO] Metadata file: ", metadata_file)
message("[INFO] ECtyper file: ", ectyper_file)

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

if (!file.exists(ectyper_file)) {
  stop("ECtyper output file not found: ", ectyper_file)
}

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

message("[INFO] Joined rows: ", nrow(joined))

unmatched_metadata <- joined %>%
  filter(source_group_clean == "Unknown_source_group") %>%
  nrow()

message("[INFO] Rows without matched metadata/source group: ", unmatched_metadata)

write_csv(joined, out_joined)

summary_by_source <- joined %>%
  group_by(source_group_clean) %>%
  summarise(
    total_done = n(),
    serotype_called = sum(Serotype_clean != "No_serotype_call"),
    no_serotype_call = sum(Serotype_clean == "No_serotype_call"),
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

message("")
message("=== Summary by source group ===")
print(summary_by_source, n = Inf)

message("")
message("=== Top 15 serotypes per source group ===")

top_serotypes <- serotype_counts_by_source %>%
  group_by(source_group_clean) %>%
  slice_max(order_by = n, n = 15, with_ties = FALSE) %>%
  ungroup()

print(top_serotypes, n = Inf)

message("")
message("[DONE] Wrote:")
message("  ", out_joined)
message("  ", out_summary)
message("  ", out_serotype_counts)
message("  ", out_otype_counts)
message("  ", out_htype_counts)

message("=== ECtyper serotype summary completed: ", Sys.time(), " ===")
