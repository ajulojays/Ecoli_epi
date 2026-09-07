#!/usr/bin/env Rscript

# ============================================================
# Script 02C
# Purpose:
#   1. Derive and normalize Location using 'Location', 'Food origin', and 'Collected by'
#   2. Standardize all US states to 2-letter postal codes (USA:XX)
#   3. Update data6a.csv with final normalized Location
#   4. Create iTOL annotation files:
#      - Collection_year_temporal_group (COLORSTRIP)
#      - Isolation_source_stage (COLORSTRIP)
#      - Location (COLORSTRIP)
#      - Original Isolation source (TEXT dataset)
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

message("=== SCRIPT 02C: CREATE iTOL ANNOTATIONS & UPDATE DATA6A LOCATION ===")

data6a_file <- "/home/samuelajulo/epi/data6a.csv"

out_dir <- "/home/samuelajulo/epi/Ecoli_epi/results/02C_itol_annotations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data6a_file)) {
  stop("[ERROR] data6a file not found: ", data6a_file, call. = FALSE)
}

data6a <- read_csv(data6a_file, show_col_types = FALSE)

required_cols <- c(
  "Collection_year_temporal_group",
  "Isolation_source_stage",
  "Isolation source",
  "Location",
  "Food origin",
  "Collected by"
)

missing_cols <- setdiff(required_cols, names(data6a))

if (length(missing_cols) > 0) {
  stop(
    "[ERROR] Missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    "\nAvailable columns:\n",
    paste(names(data6a), collapse = "\n"),
    call. = FALSE
  )
}

id_candidates <- c(
  "Assembly",
  "assembly",
  "assembly_key",
  "Assembly_accession",
  "assembly_accession",
  "Accession",
  "accession",
  "Name",
  "name"
)

id_col <- id_candidates[id_candidates %in% names(data6a)][1]

if (is.na(id_col)) {
  stop(
    "[ERROR] Could not find a tree ID column. Expected one of: ",
    paste(id_candidates, collapse = ", "),
    call. = FALSE
  )
}

message("[INFO] Using tree tip ID column: ", id_col)

normalize_tip_id <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  x <- str_remove(x, "\\.gz$")
  x <- str_remove(x, "\\.(fasta|fna|fa|gff|gff3)$")
  accession_like <- str_extract(x, "GC[AF]_[0-9]+\\.[0-9]+")
  ifelse(is.na(accession_like), x, accession_like)
}

clean_label <- function(x, missing_label) {
  x <- as.character(x)
  x <- str_squish(x)
  x <- if_else(is.na(x) | x == "", missing_label, x)
  x
}

sanitize_itol_text <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\t", " ")
  x <- str_replace_all(x, "\r|\n", " ")
  x <- str_squish(x)
  x
}

state_map <- c(
  "Alabama" = "AL", "Alaska" = "AK", "Arizona" = "AZ", "Arkansas" = "AR",
  "California" = "CA", "Colorado" = "CO", "Connecticut" = "CT", "Delaware" = "DE",
  "Florida" = "FL", "Georgia" = "GA", "Hawaii" = "HI", "Idaho" = "ID",
  "Illinois" = "IL", "Indiana" = "IN", "Iowa" = "IA", "Kansas" = "KS",
  "Kentucky" = "KY", "Louisiana" = "LA", "Maine" = "ME", "Maryland" = "MD",
  "Massachusetts" = "MA", "Michigan" = "MI", "Minnesota" = "MN", "Mississippi" = "MS",
  "Missouri" = "MO", "Montana" = "MT", "Nebraska" = "NE", "Nevada" = "NV",
  "New Hampshire" = "NH", "New Jersey" = "NJ", "New Mexico" = "NM", "New York" = "NY",
  "North Carolina" = "NC", "North Dakota" = "ND", "Ohio" = "OH", "Oklahoma" = "OK",
  "Oregon" = "OR", "Pennsylvania" = "PA", "Rhode Island" = "RI", "South Carolina" = "SC",
  "South Dakota" = "SD", "Tennessee" = "TN", "Texas" = "TX", "Utah" = "UT",
  "Vermont" = "VT", "Virginia" = "VA", "Washington" = "WA", "West Virginia" = "WV",
  "Wisconsin" = "WI", "Wyoming" = "WY", "Puerto Rico" = "PR"
)

normalize_state_code <- function(loc) {
  if (is.na(loc) || loc == "" || loc == "NA") return("NA")
  
  loc <- str_replace(loc, "^USA\\s*:\\s*", "USA:")
  
  if (str_detect(loc, "^USA:")) {
    sub_loc <- str_remove(loc, "^USA:")
    sub_loc_clean <- str_squish(sub_loc)
    
    for (st_name in names(state_map)) {
      if (str_detect(sub_loc_clean, regex(paste0("\\b", st_name, "\\b"), ignore_case = TRUE))) {
        return(paste0("USA:", state_map[[st_name]]))
      }
    }
    
    if (str_detect(sub_loc_clean, "^[A-Za-z]{2}$")) {
      return(paste0("USA:", str_to_upper(sub_loc_clean)))
    }
    
    return(paste0("USA:", sub_loc_clean))
  }
  
  return(loc)
}

derive_location <- function(loc, food_orig, collected_by) {
  loc_clean  <- str_squish(as.character(loc))
  food_clean <- str_squish(as.character(food_orig))
  coll_clean <- str_squish(as.character(collected_by))

  is_empty <- function(x) is.na(x) | x == "" | x == "NA"
  
  is_just_usa <- function(x) {
    if (is_empty(x)) return(FALSE)
    x_upper <- str_to_upper(x)
    x_upper %in% c("USA", "USA:", "USA :", "UNITED STATES", "US")
  }

  raw_loc <- "NA"

  if (!is_empty(loc_clean) && !is_just_usa(loc_clean)) {
    raw_loc <- loc_clean
  }
  else if (!is_empty(food_clean) && !is_just_usa(food_clean)) {
    raw_loc <- food_clean
  }
  else if (!is_empty(coll_clean)) {
    if (str_detect(coll_clean, regex("Kansas State University", ignore_case = TRUE))) {
      raw_loc <- "USA:KS"
    } else if (str_detect(coll_clean, regex("Nebraska", ignore_case = TRUE))) {
      raw_loc <- "USA:NE"
    } else if (is_just_usa(loc_clean) || is_just_usa(food_clean) || is_just_usa(coll_clean)) {
      raw_loc <- "USA"
    }
  }
  else if (is_just_usa(loc_clean) || is_just_usa(food_clean)) {
    raw_loc <- "USA"
  }

  final_loc <- normalize_state_code(raw_loc)
  return(final_loc)
}

# Apply derived & normalized Location directly to data6a
data6a$Location <- mapply(
  derive_location, 
  data6a[["Location"]], 
  data6a[["Food origin"]], 
  data6a[["Collected by"]], 
  USE.NAMES = FALSE
)

# Overwrite data6a.csv
write_csv(data6a, data6a_file)
message("[SUCCESS] Updated ", data6a_file, " with final normalized Location values.")

data6a_itol <- data6a %>%
  mutate(
    Tip_ID = normalize_tip_id(.data[[id_col]]),

    Collection_year_temporal_group = clean_label(
      Collection_year_temporal_group,
      "Unknown year"
    ),

    Isolation_source_stage = clean_label(
      Isolation_source_stage,
      "Unknown source stage"
    ),

    Isolation_source_original = clean_label(
      .data[["Isolation source"]],
      "Unknown isolation source"
    ),
    Isolation_source_original = sanitize_itol_text(Isolation_source_original),

    Location = sanitize_itol_text(Location)
  ) %>%
  distinct(Tip_ID, .keep_all = TRUE)

if (anyDuplicated(data6a_itol$Tip_ID) > 0) {
  duplicated_ids <- data6a_itol$Tip_ID[duplicated(data6a_itol$Tip_ID)]
  stop(
    "[ERROR] Duplicate Tip_ID values detected. First duplicates: ",
    paste(head(duplicated_ids, 20), collapse = ", "),
    call. = FALSE
  )
}

make_color_map <- function(values, fixed_colors = NULL, palette = "Dynamic") {
  values <- sort(unique(values))
  n_vals <- length(values)

  if (palette == "Distinct" || n_vals > 12) {
    # Generate evenly spaced, distinct HCL colors around the color wheel
    hues <- seq(15, 375, length.out = n_vals + 1)[1:n_vals]
    cols <- hcl(h = hues, c = 85, l = 55)
  } else {
    cols <- grDevices::hcl.colors(n_vals, palette = palette)
  }

  color_map <- tibble(
    value = values,
    color = cols
  )

  if (!is.null(fixed_colors)) {
    fixed_tbl <- tibble(
      value = names(fixed_colors),
      fixed_color = unname(fixed_colors)
    )

    color_map <- color_map %>%
      left_join(fixed_tbl, by = "value") %>%
      mutate(color = if_else(!is.na(fixed_color), fixed_color, color)) %>%
      select(value, color)
  }

  color_map
}

# 1. Isolation Stage Color Map
stage_fixed_colors <- c(
  "pre-harvest and harvest" = "#1B9E77",
  "post-harvest" = "#D95F02",
  "Unknown source stage" = "#7570B3"
)

stage_color_map <- make_color_map(
  data6a_itol$Isolation_source_stage,
  fixed_colors = stage_fixed_colors,
  palette = "Dark 3"
) %>%
  mutate(
    stage_order = match(value, c("pre-harvest and harvest", "post-harvest", "Unknown source stage"))
  ) %>%
  arrange(stage_order, value) %>%
  select(-stage_order)

# 2. Collection Year Color Map
year_values <- unique(data6a_itol$Collection_year_temporal_group)

year_order <- tibble(value = year_values) %>%
  mutate(
    num_part = suppressWarnings(as.numeric(str_extract(value, "[0-9]{4}"))),
    sort_key = case_when(
      value == "<=2000" ~ 0,
      value == "Unknown year" ~ 999999,
      !is.na(num_part) ~ num_part,
      TRUE ~ 999998
    )
  ) %>%
  arrange(sort_key, value) %>%
  pull(value)

year_color_map <- tibble(
  value = year_order,
  color = grDevices::hcl.colors(length(year_order), palette = "Viridis")
) %>%
  mutate(
    color = case_when(
      value == "<=2000" ~ "#666666",
      value == "Unknown year" ~ "#BDBDBD",
      TRUE ~ color
    )
  )

# 3. Location Color Map (Distinct 29-color spectrum generator)
location_color_map <- make_color_map(
  data6a_itol$Location,
  palette = "Distinct"
)

write_itol_colorstrip <- function(
  df,
  color_map,
  value_col,
  dataset_label,
  output_file,
  show_legend = TRUE
) {
  plot_df <- df %>%
    select(Tip_ID, value = all_of(value_col)) %>%
    left_join(color_map, by = "value") %>%
    transmute(
      Tip_ID = Tip_ID,
      color = color,
      label = value
    )

  legend_map <- color_map %>%
    filter(value %in% unique(plot_df$label)) %>%
    distinct(value, color)

  con <- file(output_file, open = "wt")

  writeLines("DATASET_COLORSTRIP", con)
  writeLines("SEPARATOR TAB", con)
  writeLines(paste("DATASET_LABEL", dataset_label, sep = "\t"), con)
  writeLines(paste("COLOR", "#000000", sep = "\t"), con)
  writeLines(paste("STRIP_WIDTH", "40", sep = "\t"), con)
  writeLines(paste("MARGIN", "5", sep = "\t"), con)
  writeLines(paste("BORDER_WIDTH", "1", sep = "\t"), con)
  writeLines(paste("BORDER_COLOR", "#000000", sep = "\t"), con)
  writeLines(paste("SHOW_INTERNAL", "0", sep = "\t"), con)

  if (show_legend) {
    writeLines(paste("LEGEND_TITLE", dataset_label, sep = "\t"), con)

    writeLines(
      paste(
        c("LEGEND_SHAPES", rep("1", nrow(legend_map))),
        collapse = "\t"
      ),
      con
    )

    writeLines(
      paste(
        c("LEGEND_COLORS", legend_map$color),
        collapse = "\t"
      ),
      con
    )

    writeLines(
      paste(
        c("LEGEND_LABELS", legend_map$value),
        collapse = "\t"
      ),
      con
    )
  }

  writeLines("DATA", con)

  apply(plot_df, 1, function(row) {
    writeLines(
      paste(
        row[["Tip_ID"]],
        row[["color"]],
        row[["label"]],
        sep = "\t"
      ),
      con
    )
  })

  close(con)

  invisible(plot_df)
}

write_itol_text_dataset <- function(
  df,
  value_col,
  dataset_label,
  output_file
) {
  text_df <- df %>%
    select(Tip_ID, text = all_of(value_col)) %>%
    mutate(
      text = sanitize_itol_text(text),
      position = "-1",
      color = "#000000",
      style = "normal",
      size_factor = "1",
      rotation = "0"
    )

  con <- file(output_file, open = "wt")

  writeLines("DATASET_TEXT", con)
  writeLines("SEPARATOR TAB", con)
  writeLines(paste("DATASET_LABEL", dataset_label, sep = "\t"), con)
  writeLines(paste("COLOR", "#000000", sep = "\t"), con)
  writeLines(paste("MARGIN", "5", sep = "\t"), con)
  writeLines(paste("SHOW_INTERNAL", "0", sep = "\t"), con)
  writeLines(paste("ALIGN_TO_TREE", "0", sep = "\t"), con)
  writeLines(paste("SIZE_FACTOR", "1", sep = "\t"), con)
  writeLines(paste("ROTATION", "0", sep = "\t"), con)
  writeLines("DATA", con)

  apply(text_df, 1, function(row) {
    writeLines(
      paste(
        row[["Tip_ID"]],
        row[["text"]],
        row[["position"]],
        row[["color"]],
        row[["style"]],
        row[["size_factor"]],
        row[["rotation"]],
        sep = "\t"
      ),
      con
    )
  })

  close(con)

  invisible(text_df)
}

# Output file paths
year_itol_file <- file.path(out_dir, "itol_collection_year_temporal_group_colorstrip.txt")
stage_itol_file <- file.path(out_dir, "itol_isolation_source_stage_colorstrip.txt")
location_itol_file <- file.path(out_dir, "itol_location_colorstrip.txt")
source_text_file <- file.path(out_dir, "itol_isolation_source_original_text.txt")

# Write Colorstrips
write_itol_colorstrip(
  df = data6a_itol,
  color_map = year_color_map,
  value_col = "Collection_year_temporal_group",
  dataset_label = "Collection year temporal group",
  output_file = year_itol_file,
  show_legend = TRUE
)

write_itol_colorstrip(
  df = data6a_itol,
  color_map = stage_color_map,
  value_col = "Isolation_source_stage",
  dataset_label = "Isolation source stage",
  output_file = stage_itol_file,
  show_legend = TRUE
)

write_itol_colorstrip(
  df = data6a_itol,
  color_map = location_color_map,
  value_col = "Location",
  dataset_label = "Location",
  output_file = location_itol_file,
  show_legend = TRUE
)

# Write Text Dataset
write_itol_text_dataset(
  df = data6a_itol,
  value_col = "Isolation_source_original",
  dataset_label = "Isolation source",
  output_file = source_text_file
)

# Save Master Data and Color Maps
write_csv(
  data6a_itol %>%
    select(
      Tip_ID,
      Collection_year_temporal_group,
      Isolation_source_stage,
      Isolation_source_original,
      Location
    ),
  file.path(out_dir, "data6a_itol_annotation_input.csv")
)

write_csv(
  year_color_map,
  file.path(out_dir, "itol_collection_year_temporal_group_color_map.csv")
)

write_csv(
  stage_color_map,
  file.path(out_dir, "itol_isolation_source_stage_color_map.csv")
)

write_csv(
  location_color_map,
  file.path(out_dir, "itol_location_color_map.csv")
)

# Save Summary Counts Tables
write_csv(
  data6a_itol %>% count(Collection_year_temporal_group, name = "n", sort = TRUE),
  file.path(out_dir, "data6a_collection_year_temporal_group_counts.csv")
)

write_csv(
  data6a_itol %>% count(Isolation_source_stage, name = "n", sort = TRUE),
  file.path(out_dir, "data6a_isolation_source_stage_counts.csv")
)

write_csv(
  data6a_itol %>% count(Isolation_source_original, name = "n", sort = TRUE),
  file.path(out_dir, "data6a_isolation_source_original_counts.csv")
)

write_csv(
  data6a_itol %>% count(Location, name = "n", sort = TRUE),
  file.path(out_dir, "data6a_location_counts.csv")
)

message("[DONE] iTOL annotation files created and data6a.csv updated.")
message("[OUTPUT] Colorstrip: ", location_itol_file)

message("")
message("Normalized Location counts:")
print(data6a_itol %>% count(Location, name = "n", sort = TRUE), n = Inf)
