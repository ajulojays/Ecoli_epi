#!/usr/bin/env Rscript

# ==============================================================================
# 12E Add-on: Official Clade-Defining SNP Manifest
#
# Criterion:
#   A SNP allele is clade-defining when:
#     1. allele prevalence inside the clade is >=40%
#     2. allele prevalence outside the clade is <1%
#
# Outputs:
#   - one manifest of all clade-defining SNP alleles
#   - summary by clade
#   - optional gene/impact annotation if a core-alignment position-to-gene map exists
#
# Important:
#   Gene/impact annotation requires a map from Panaroo core-alignment positions
#   to genes. If no map is present, gene and impact fields are retained but marked
#   as unmapped rather than guessed.
# ==============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(forcats)
})

message("=== 12E OFFICIAL CLADE-DEFINING SNP MANIFEST: 40% IN-CLADE, <1% OUT-CLADE ===")

# ==============================================================================
# 1. Paths
# ==============================================================================

base_dir <- "/home/samuelajulo/epi/Ecoli_epi"

alignment_filtered <- file.path(
  base_dir,
  "results",
  "12c_all_preharvest_combined",
  "panaroo",
  "core_gene_alignment_filtered.aln"
)

alignment_unfiltered <- file.path(
  base_dir,
  "results",
  "12c_all_preharvest_combined",
  "panaroo",
  "core_gene_alignment.aln"
)

clade_file <- file.path(
  base_dir,
  "results",
  "12d_clade_analysis",
  "all_preharvest_704_clades.csv"
)

clade_counts_file <- file.path(
  base_dir,
  "results",
  "12d_clade_analysis",
  "clade_summary_counts.csv"
)

metadata_file <- file.path(
  base_dir,
  "results",
  "12c_all_preharvest_combined",
  "target_metadata_704_preharvest_non_dairy.csv"
)

metadata_fallback <- "/home/samuelajulo/epi/data6a.csv"

out_dir <- file.path(
  base_dir,
  "results",
  "12e_pairwise_snp_distance",
  "clade_defining_snps_40pct"
)

fig_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

alignment_file <- if (file.exists(alignment_filtered)) {
  alignment_filtered
} else {
  alignment_unfiltered
}

if (!file.exists(alignment_file)) {
  stop(
    "[ERROR] Core alignment not found. Expected one of:\n",
    alignment_filtered, "\n",
    alignment_unfiltered,
    call. = FALSE
  )
}

if (!file.exists(clade_file)) {
  stop("[ERROR] Clade file not found: ", clade_file, call. = FALSE)
}

if (!file.exists(clade_counts_file)) {
  stop("[ERROR] Clade counts file not found: ", clade_counts_file, call. = FALSE)
}

metadata_used <- if (file.exists(metadata_file)) {
  metadata_file
} else {
  metadata_fallback
}

if (!file.exists(metadata_used)) {
  stop("[ERROR] Metadata file not found: ", metadata_used, call. = FALSE)
}

# Optional position-to-gene map candidates
gene_map_candidates <- c(
  file.path(out_dir, "core_alignment_position_gene_map.csv"),
  file.path(base_dir, "results", "12e_pairwise_snp_distance", "core_alignment_position_gene_map.csv"),
  file.path(base_dir, "results", "12c_all_preharvest_combined", "panaroo", "core_alignment_position_gene_map.csv"),
  file.path(base_dir, "results", "12c_all_preharvest_combined", "panaroo", "core_gene_alignment_position_gene_map.csv"),
  file.path(base_dir, "results", "12c_all_preharvest_combined", "panaroo", "core_alignment_gene_map.csv")
)

gene_map_file <- gene_map_candidates[file.exists(gene_map_candidates)][1]

# ==============================================================================
# 2. Criteria
# ==============================================================================

IN_CLADE_MIN_PREVALENCE <- 0.40
OUT_CLADE_MAX_PREVALENCE <- 0.01
EXPECTED_N <- 704
BLOCK_SIZE <- 50000
bases <- c("A", "C", "G", "T")

criterion_text <- paste0(
  "in-clade allele prevalence >= ", IN_CLADE_MIN_PREVALENCE * 100,
  "%; out-of-clade allele prevalence < ", OUT_CLADE_MAX_PREVALENCE * 100, "%"
)

message("[INFO] Alignment file: ", alignment_file)
message("[INFO] Clade file: ", clade_file)
message("[INFO] Metadata file: ", metadata_used)
message("[INFO] Output directory: ", out_dir)
message("[INFO] Official criterion: ", criterion_text)

if (!is.na(gene_map_file)) {
  message("[INFO] Gene map detected: ", gene_map_file)
} else {
  message("[WARN] No core-alignment position-to-gene map detected.")
  message("[WARN] Gene/impact columns will be retained but marked as unmapped.")
}

# ==============================================================================
# 3. Helpers
# ==============================================================================

normalize_tip_id <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  x <- str_remove(x, "\\.gz$")
  x <- str_remove(x, "\\.(fasta|fna|fa|gff|gff3|aln)$")
  accession_like <- str_extract(x, "GC[AF]_[0-9]+\\.[0-9]+")
  ifelse(is.na(accession_like), x, accession_like)
}

extract_year_num <- function(x) {
  x <- as.character(x)
  out <- str_extract(x, "(18|19|20)[0-9]{2}")
  suppressWarnings(as.integer(out))
}

pick_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}

safe_pct <- function(x) {
  round(100 * x, 3)
}

standardize_gene_map <- function(gmap) {
  names(gmap) <- make.names(names(gmap))

  pos_col <- pick_first_existing_col(
    gmap,
    c(
      "Alignment_position",
      "alignment_position",
      "Aln_position",
      "aln_position",
      "Position",
      "position",
      "Core_alignment_position",
      "core_alignment_position"
    )
  )

  if (is.na(pos_col)) {
    stop("[ERROR] Gene map exists but no alignment-position column was detected.", call. = FALSE)
  }

  out <- gmap %>%
    mutate(Alignment_position = as.integer(.data[[pos_col]]))

  wanted <- list(
    Gene = c("Gene", "gene", "Gene_ID", "gene_id", "Locus_tag", "locus_tag", "Name", "name"),
    Product = c("Product", "product", "Annotation", "annotation", "Description", "description"),
    Feature_type = c("Feature_type", "feature_type", "Type", "type"),
    Contig = c("Contig", "contig", "Sequence", "sequence", "Seqid", "seqid"),
    Gene_start = c("Gene_start", "gene_start", "Start", "start"),
    Gene_end = c("Gene_end", "gene_end", "End", "end"),
    Strand = c("Strand", "strand"),
    Codon_position = c("Codon_position", "codon_position"),
    Reference_codon = c("Reference_codon", "reference_codon", "Ref_codon", "ref_codon"),
    Alternate_codon = c("Alternate_codon", "alternate_codon", "Alt_codon", "alt_codon"),
    Reference_AA = c("Reference_AA", "reference_aa", "Ref_AA", "ref_aa"),
    Alternate_AA = c("Alternate_AA", "alternate_aa", "Alt_AA", "alt_aa"),
    Likely_impact = c("Likely_impact", "likely_impact", "Impact", "impact", "Effect", "effect")
  )

  for (new_col in names(wanted)) {
    old_col <- pick_first_existing_col(out, wanted[[new_col]])
    if (is.na(old_col)) {
      out[[new_col]] <- NA_character_
    } else {
      out[[new_col]] <- as.character(out[[old_col]])
    }
  }

  out %>%
    select(
      Alignment_position,
      Gene,
      Product,
      Feature_type,
      Contig,
      Gene_start,
      Gene_end,
      Strand,
      Codon_position,
      Reference_codon,
      Alternate_codon,
      Reference_AA,
      Alternate_AA,
      Likely_impact
    ) %>%
    distinct(Alignment_position, .keep_all = TRUE)
}

# ==============================================================================
# 4. Read clades, clade counts, metadata
# ==============================================================================

clades <- read_csv(clade_file, show_col_types = FALSE) %>%
  mutate(
    Tip_ID = normalize_tip_id(Tip_ID),
    Clade = as.character(Clade)
  ) %>%
  distinct(Tip_ID, .keep_all = TRUE)

clade_counts <- read_csv(clade_counts_file, show_col_types = FALSE) %>%
  mutate(Clade = as.character(Clade)) %>%
  arrange(desc(Isolate_Count), Clade)

if (!"Isolate_Count" %in% names(clade_counts)) {
  clade_counts <- clades %>%
    count(Clade, name = "Isolate_Count") %>%
    arrange(desc(Isolate_Count), Clade)
}

metadata <- read_csv(metadata_used, show_col_types = FALSE)

assembly_col <- pick_first_existing_col(
  metadata,
  c("Assembly", "Tip_ID", "Isolate", "Strain", "Sample", "sample")
)

if (is.na(assembly_col)) {
  stop("[ERROR] Could not identify Assembly/Tip_ID column in metadata.", call. = FALSE)
}

year_col <- pick_first_existing_col(
  metadata,
  c("Collection_year", "Collection year", "collection_year", "Year", "year", "Collection date", "collection_date")
)

metadata_clean <- metadata %>%
  mutate(
    Tip_ID = normalize_tip_id(.data[[assembly_col]]),
    Collection_year_num = if (!is.na(year_col)) extract_year_num(.data[[year_col]]) else NA_integer_
  ) %>%
  distinct(Tip_ID, .keep_all = TRUE)

sample_annotation <- clades %>%
  left_join(metadata_clean, by = "Tip_ID") %>%
  left_join(clade_counts, by = "Clade")

write_csv(
  sample_annotation,
  file.path(out_dir, "12E_40pct_sample_clade_metadata_annotation.csv")
)

# ==============================================================================
# 5. Read core alignment and extract SNP sites
# ==============================================================================

message("[INFO] Reading core alignment...")

dna <- read.dna(alignment_file, format = "fasta", as.character = TRUE)

if (!is.matrix(dna)) {
  dna <- as.matrix(dna)
}

clean_ids <- normalize_tip_id(rownames(dna))

if (anyDuplicated(clean_ids) > 0) {
  stop(
    "[ERROR] Duplicate normalized alignment IDs detected: ",
    paste(head(unique(clean_ids[duplicated(clean_ids)]), 10), collapse = ", "),
    call. = FALSE
  )
}

rownames(dna) <- clean_ids
dna <- dna[order(rownames(dna)), , drop = FALSE]
dna <- toupper(dna)

n_genomes <- nrow(dna)
alignment_length <- ncol(dna)

message("[INFO] Genomes in alignment: ", n_genomes)
message("[INFO] Alignment length: ", alignment_length)

if (n_genomes != EXPECTED_N) {
  warning("[WARN] Expected ", EXPECTED_N, " genomes but found ", n_genomes, ".")
}

missing_clade <- setdiff(rownames(dna), clades$Tip_ID)
missing_alignment <- setdiff(clades$Tip_ID, rownames(dna))

if (length(missing_clade) > 0) {
  write_lines(missing_clade, file.path(log_dir, "alignment_tips_missing_from_clades.txt"))
  stop("[ERROR] Alignment tips missing from clade file. See logs.", call. = FALSE)
}

if (length(missing_alignment) > 0) {
  write_lines(missing_alignment, file.path(log_dir, "clade_tips_missing_from_alignment.txt"))
  stop("[ERROR] Clade tips missing from alignment. See logs.", call. = FALSE)
}

message("[INFO] Identifying variable A/C/G/T SNP sites...")

var_sites <- logical(alignment_length)

for (block_start in seq(1, alignment_length, by = BLOCK_SIZE)) {
  block_end <- min(block_start + BLOCK_SIZE - 1, alignment_length)
  idx <- block_start:block_end

  block <- dna[, idx, drop = FALSE]
  encoded_block <- matrix(
    match(as.vector(block), bases),
    nrow = n_genomes,
    ncol = length(idx)
  )

  allele_presence_count <-
    (colSums(encoded_block == 1, na.rm = TRUE) > 0) +
    (colSums(encoded_block == 2, na.rm = TRUE) > 0) +
    (colSums(encoded_block == 3, na.rm = TRUE) > 0) +
    (colSums(encoded_block == 4, na.rm = TRUE) > 0)

  var_sites[idx] <- allele_presence_count >= 2

  message("[INFO] SNP scan block ", block_start, "-", block_end, " complete.")
}

snp_site_idx <- which(var_sites)
n_snp_sites <- length(snp_site_idx)

if (n_snp_sites == 0) {
  stop("[ERROR] No A/C/G/T SNP sites detected.", call. = FALSE)
}

message("[INFO] Core SNP sites detected: ", n_snp_sites)

snp_dna <- dna[, snp_site_idx, drop = FALSE]

encoded_snp <- matrix(
  match(as.vector(snp_dna), bases),
  nrow = n_genomes,
  ncol = n_snp_sites
)

rownames(encoded_snp) <- rownames(dna)
colnames(encoded_snp) <- as.character(snp_site_idx)

valid_snp <- !is.na(encoded_snp)

snp_site_summary <- tibble(
  alignment_file = alignment_file,
  n_genomes = n_genomes,
  alignment_length_bp = alignment_length,
  core_snp_sites_acgt = n_snp_sites,
  criterion = criterion_text
)

write_csv(
  snp_site_summary,
  file.path(out_dir, "12E_40pct_clade_defining_snp_site_summary.csv")
)

# ==============================================================================
# 6. Optional gene map
# ==============================================================================

if (!is.na(gene_map_file)) {
  gene_map <- read_csv(gene_map_file, show_col_types = FALSE) %>%
    standardize_gene_map()

  gene_map_status <- tibble(
    gene_map_status = "detected_and_loaded",
    gene_map_file = gene_map_file,
    n_positions_mapped = nrow(gene_map)
  )
} else {
  gene_map <- tibble(
    Alignment_position = integer(),
    Gene = character(),
    Product = character(),
    Feature_type = character(),
    Contig = character(),
    Gene_start = character(),
    Gene_end = character(),
    Strand = character(),
    Codon_position = character(),
    Reference_codon = character(),
    Alternate_codon = character(),
    Reference_AA = character(),
    Alternate_AA = character(),
    Likely_impact = character()
  )

  gene_map_status <- tibble(
    gene_map_status = "not_detected",
    gene_map_file = NA_character_,
    n_positions_mapped = 0L
  )
}

write_csv(
  gene_map_status,
  file.path(out_dir, "12E_40pct_gene_impact_mapping_status.csv")
)

template_map <- tibble(
  Alignment_position = snp_site_idx,
  Gene = NA_character_,
  Product = NA_character_,
  Feature_type = NA_character_,
  Contig = NA_character_,
  Gene_start = NA_character_,
  Gene_end = NA_character_,
  Strand = NA_character_,
  Codon_position = NA_character_,
  Reference_codon = NA_character_,
  Alternate_codon = NA_character_,
  Reference_AA = NA_character_,
  Alternate_AA = NA_character_,
  Likely_impact = NA_character_
)

write_csv(
  template_map,
  file.path(out_dir, "core_alignment_position_gene_map_TEMPLATE.csv")
)

# ==============================================================================
# 7. Identify clade-defining SNP alleles
# ==============================================================================

message("[INFO] Scanning clade-defining SNP alleles using official 40%/<1% criterion...")

alignment_clade_map <- tibble(
  Tip_ID = rownames(encoded_snp)
) %>%
  left_join(clades %>% select(Tip_ID, Clade), by = "Tip_ID") %>%
  left_join(clade_counts, by = "Clade")

all_clades <- sort(unique(alignment_clade_map$Clade))

manifest_rows <- list()
row_counter <- 1

for (cl in all_clades) {
  message("[INFO] Processing ", cl)

  in_idx <- which(alignment_clade_map$Clade == cl)
  out_idx <- which(alignment_clade_map$Clade != cl)

  in_valid <- colSums(valid_snp[in_idx, , drop = FALSE])
  out_valid <- colSums(valid_snp[out_idx, , drop = FALSE])
  global_valid <- in_valid + out_valid

  for (allele_code in seq_along(bases)) {
    allele <- bases[allele_code]

    allele_mat <- (!is.na(encoded_snp)) & (encoded_snp == allele_code)

    in_count <- colSums(allele_mat[in_idx, , drop = FALSE])
    out_count <- colSums(allele_mat[out_idx, , drop = FALSE])
    global_count <- in_count + out_count

    in_prev <- ifelse(in_valid > 0, in_count / in_valid, NA_real_)
    out_prev <- ifelse(out_valid > 0, out_count / out_valid, NA_real_)
    global_prev <- ifelse(global_valid > 0, global_count / global_valid, NA_real_)

    hit <- !is.na(in_prev) &
      !is.na(out_prev) &
      in_prev >= IN_CLADE_MIN_PREVALENCE &
      out_prev < OUT_CLADE_MAX_PREVALENCE

    if (any(hit)) {
      hit_idx <- which(hit)

      manifest_rows[[row_counter]] <- tibble(
        Clade = cl,
        Alignment_position = snp_site_idx[hit_idx],
        SNP_site_rank = hit_idx,
        Clade_defining_allele = allele,
        In_clade_allele_count = as.integer(in_count[hit_idx]),
        In_clade_total_valid = as.integer(in_valid[hit_idx]),
        In_clade_prevalence = round(in_prev[hit_idx], 5),
        In_clade_percent = round(100 * in_prev[hit_idx], 3),
        Out_clade_allele_count = as.integer(out_count[hit_idx]),
        Out_clade_total_valid = as.integer(out_valid[hit_idx]),
        Out_clade_prevalence = round(out_prev[hit_idx], 5),
        Out_clade_percent = round(100 * out_prev[hit_idx], 3),
        Global_allele_count = as.integer(global_count[hit_idx]),
        Global_total_valid = as.integer(global_valid[hit_idx]),
        Global_allele_prevalence = round(global_prev[hit_idx], 5),
        Global_allele_percent = round(100 * global_prev[hit_idx], 3),
        Criterion = criterion_text
      )

      row_counter <- row_counter + 1
    }
  }
}

if (length(manifest_rows) == 0) {
  clade_snp_manifest <- tibble(
    Clade = character(),
    Alignment_position = integer(),
    SNP_site_rank = integer(),
    Clade_defining_allele = character(),
    In_clade_allele_count = integer(),
    In_clade_total_valid = integer(),
    In_clade_prevalence = numeric(),
    In_clade_percent = numeric(),
    Out_clade_allele_count = integer(),
    Out_clade_total_valid = integer(),
    Out_clade_prevalence = numeric(),
    Out_clade_percent = numeric(),
    Global_allele_count = integer(),
    Global_total_valid = integer(),
    Global_allele_prevalence = numeric(),
    Global_allele_percent = numeric(),
    Criterion = character()
  )
} else {
  clade_snp_manifest <- bind_rows(manifest_rows)
}

# Add clade counts and optional gene/impact annotation
clade_snp_manifest <- clade_snp_manifest %>%
  left_join(clade_counts, by = "Clade") %>%
  left_join(gene_map, by = "Alignment_position") %>%
  mutate(
    Gene_mapping_status = case_when(
      !is.na(Gene) & Gene != "" ~ "mapped_to_gene",
      TRUE ~ "unmapped_no_core_alignment_position_to_gene_map"
    ),
    Likely_impact = case_when(
      !is.na(Likely_impact) & Likely_impact != "" ~ Likely_impact,
      Gene_mapping_status == "mapped_to_gene" ~ "genic SNP; coding effect not inferred in this run",
      TRUE ~ "unknown; no core-alignment position-to-gene map available"
    ),
    Impact_confidence = case_when(
      Gene_mapping_status == "mapped_to_gene" & !is.na(Likely_impact) & Likely_impact != "" ~ "from_position_gene_impact_map",
      Gene_mapping_status == "mapped_to_gene" ~ "gene_mapped_effect_not_inferred",
      TRUE ~ "not_assessed"
    ),
    SNP_manifest_id = paste(Clade, Alignment_position, Clade_defining_allele, sep = "_")
  ) %>%
  select(
    SNP_manifest_id,
    Clade,
    Isolate_Count,
    Alignment_position,
    SNP_site_rank,
    Clade_defining_allele,
    In_clade_allele_count,
    In_clade_total_valid,
    In_clade_prevalence,
    In_clade_percent,
    Out_clade_allele_count,
    Out_clade_total_valid,
    Out_clade_prevalence,
    Out_clade_percent,
    Global_allele_count,
    Global_total_valid,
    Global_allele_prevalence,
    Global_allele_percent,
    Gene_mapping_status,
    Gene,
    Product,
    Feature_type,
    Contig,
    Gene_start,
    Gene_end,
    Strand,
    Codon_position,
    Reference_codon,
    Alternate_codon,
    Reference_AA,
    Alternate_AA,
    Likely_impact,
    Impact_confidence,
    Criterion
  ) %>%
  arrange(Clade, Alignment_position, Clade_defining_allele)

# ==============================================================================
# 8. Summary outputs
# ==============================================================================

summary_by_clade <- clade_counts %>%
  left_join(
    clade_snp_manifest %>%
      count(Clade, name = "clade_defining_snp_allele_count"),
    by = "Clade"
  ) %>%
  left_join(
    clade_snp_manifest %>%
      filter(Gene_mapping_status == "mapped_to_gene") %>%
      count(Clade, name = "mapped_to_gene_count"),
    by = "Clade"
  ) %>%
  mutate(
    clade_defining_snp_allele_count = replace_na(clade_defining_snp_allele_count, 0L),
    mapped_to_gene_count = replace_na(mapped_to_gene_count, 0L),
    unmapped_count = clade_defining_snp_allele_count - mapped_to_gene_count,
    criterion = criterion_text
  ) %>%
  arrange(desc(Isolate_Count), Clade)

summary_by_gene <- clade_snp_manifest %>%
  filter(Gene_mapping_status == "mapped_to_gene") %>%
  count(Clade, Gene, Product, Likely_impact, name = "clade_defining_snp_count") %>%
  arrange(Clade, desc(clade_defining_snp_count), Gene)

write_csv(
  clade_snp_manifest,
  file.path(out_dir, "12E_clade_defining_snp_manifest_40pct_inclade_outlt1pct.csv.gz")
)

write_csv(
  summary_by_clade,
  file.path(out_dir, "12E_clade_defining_snp_summary_by_clade_40pct.csv")
)

write_csv(
  summary_by_gene,
  file.path(out_dir, "12E_clade_defining_snp_summary_by_gene_40pct.csv")
)

readme <- c(
  "12E official clade-defining SNP criterion",
  "",
  "Criterion:",
  paste0("A SNP allele is called clade-defining when in-clade allele prevalence is >= ", IN_CLADE_MIN_PREVALENCE * 100, "% and out-of-clade allele prevalence is < ", OUT_CLADE_MAX_PREVALENCE * 100, "%."),
  "",
  "Interpretation:",
  "These SNP alleles are candidate clade-defining genomic markers because they are common inside a clade and nearly absent outside that clade.",
  "",
  "Important caveat:",
  "This is a marker-enrichment definition, not a functional proof.",
  "",
  "Gene and impact annotation:",
  "Gene and likely-impact fields are populated only when a core-alignment position-to-gene map is available.",
  "Without that map, the manifest keeps the columns but marks positions as unmapped rather than guessing gene or impact.",
  "",
  "Main manifest:",
  "12E_clade_defining_snp_manifest_40pct_inclade_outlt1pct.csv.gz",
  "",
  "Required columns for an optional gene map:",
  "Alignment_position, Gene, Product, Feature_type, Contig, Gene_start, Gene_end, Strand, Codon_position, Reference_codon, Alternate_codon, Reference_AA, Alternate_AA, Likely_impact"
)

write_lines(
  readme,
  file.path(out_dir, "README_12E_clade_defining_snp_40pct_criterion.txt")
)

# ==============================================================================
# 9. Figure
# ==============================================================================

plot_df <- summary_by_clade %>%
  mutate(Clade = fct_reorder(Clade, Isolate_Count, .desc = TRUE))

p <- ggplot(plot_df, aes(x = Clade, y = clade_defining_snp_allele_count)) +
  geom_col(width = 0.78) +
  geom_text(
    aes(label = clade_defining_snp_allele_count),
    vjust = -0.25,
    size = 3.2
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Clade-defining SNP alleles",
    subtitle = "Official criterion: present in >=40% inside clade and <1% outside clade",
    x = "Clade",
    y = "Clade-defining SNP allele count"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(fig_dir, "12E_clade_defining_snp_counts_40pct.pdf"),
  p,
  width = 9,
  height = 5.5
)

ggsave(
  file.path(fig_dir, "12E_clade_defining_snp_counts_40pct.png"),
  p,
  width = 9,
  height = 5.5,
  dpi = 600
)

# ==============================================================================
# 10. Final report
# ==============================================================================

message("")
message("==========================================================")
message(" 12E CLADE-DEFINING SNP MANIFEST COMPLETE")
message("==========================================================")
message("[OUTPUT] ", out_dir)
message("")
message("=== SNP site summary ===")
print(snp_site_summary)
message("")
message("=== Gene/impact mapping status ===")
print(gene_map_status)
message("")
message("=== Clade-defining SNP summary ===")
print(summary_by_clade)
message("")
message("[DONE] Official 40% in-clade / <1% out-clade SNP manifest generated.")
