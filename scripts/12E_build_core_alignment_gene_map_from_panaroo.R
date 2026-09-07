#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

base_dir <- "/home/samuelajulo/epi/Ecoli_epi"

core_alignment <- file.path(
  base_dir,
  "results",
  "12c_all_preharvest_combined",
  "panaroo",
  "core_gene_alignment_filtered.aln"
)

if (!file.exists(core_alignment)) {
  core_alignment <- file.path(
    base_dir,
    "results",
    "12c_all_preharvest_combined",
    "panaroo",
    "core_gene_alignment.aln"
  )
}

aligned_gene_dir <- file.path(
  base_dir,
  "results",
  "12c_all_preharvest_combined",
  "panaroo",
  "aligned_gene_sequences"
)

out_dir <- file.path(
  base_dir,
  "results",
  "12e_pairwise_snp_distance"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_map <- file.path(out_dir, "core_alignment_position_gene_map.csv")
out_summary <- file.path(out_dir, "core_alignment_position_gene_map_summary.csv")

if (!file.exists(core_alignment)) {
  stop("[ERROR] Core alignment not found: ", core_alignment, call. = FALSE)
}

if (!dir.exists(aligned_gene_dir)) {
  stop("[ERROR] Panaroo aligned_gene_sequences directory not found: ", aligned_gene_dir, call. = FALSE)
}

message("[INFO] Core alignment: ", core_alignment)
message("[INFO] Aligned gene directory: ", aligned_gene_dir)

core_dna <- read.dna(core_alignment, format = "fasta", as.character = TRUE)
core_len <- ncol(as.matrix(core_dna))

message("[INFO] Core alignment length: ", core_len)

gene_files <- list.files(
  aligned_gene_dir,
  pattern = "\\.aln\\.fas$",
  full.names = TRUE
)

if (length(gene_files) == 0) {
  stop("[ERROR] No *.aln.fas files found in: ", aligned_gene_dir, call. = FALSE)
}

gene_files <- sort(gene_files)

message("[INFO] Gene alignment files found: ", length(gene_files))

rows <- list()
summary_rows <- list()

cursor <- 1L

for (i in seq_along(gene_files)) {
  f <- gene_files[i]

  gene_name <- basename(f) %>%
    str_remove("\\.aln\\.fas$")

  aln <- read.dna(f, format = "fasta", as.character = TRUE)
  aln <- as.matrix(aln)

  gene_len <- ncol(aln)
  start_pos <- cursor
  end_pos <- cursor + gene_len - 1L

  rows[[i]] <- tibble(
    Alignment_position = start_pos:end_pos,
    Gene = gene_name,
    Product = NA_character_,
    Feature_type = "Panaroo_core_gene_alignment_block",
    Contig = NA_character_,
    Gene_start = NA_character_,
    Gene_end = NA_character_,
    Strand = NA_character_,
    Codon_position = NA_character_,
    Reference_codon = NA_character_,
    Alternate_codon = NA_character_,
    Reference_AA = NA_character_,
    Alternate_AA = NA_character_,
    Likely_impact = "genic SNP; coding effect not inferred from block map"
  )

  summary_rows[[i]] <- tibble(
    Gene = gene_name,
    Gene_alignment_file = f,
    Block_start_alignment_position = start_pos,
    Block_end_alignment_position = end_pos,
    Gene_alignment_length = gene_len
  )

  cursor <- end_pos + 1L
}

gene_map <- bind_rows(rows)
gene_summary <- bind_rows(summary_rows)

mapped_len <- nrow(gene_map)

status <- tibble(
  core_alignment_file = core_alignment,
  aligned_gene_sequences_dir = aligned_gene_dir,
  n_gene_alignment_files = length(gene_files),
  core_alignment_length = core_len,
  mapped_gene_block_length = mapped_len,
  lengths_match = mapped_len == core_len,
  map_type = "gene-block map inferred from Panaroo aligned_gene_sequences sorted by filename",
  impact_note = "Maps SNPs to Panaroo core gene blocks. Does not infer synonymous/nonsynonymous impact."
)

write_csv(gene_map, out_map)
write_csv(gene_summary, out_summary)
write_csv(status, file.path(out_dir, "core_alignment_position_gene_map_build_status.csv"))

message("")
message("=== Build status ===")
print(status)

if (mapped_len != core_len) {
  message("")
  message("[WARNING] Mapped gene-block length does not equal core alignment length.")
  message("[WARNING] Do not use this map for gene annotation until ordering/length mismatch is resolved.")
} else {
  message("")
  message("[SUCCESS] Gene-block map created and length matches core alignment.")
  message("[OUTPUT] ", out_map)
}
