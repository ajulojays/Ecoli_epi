# Ecoli_epi

A reproducible, code-only workflow for large-scale *Escherichia coli* genome serotyping using ECtyper, followed by source-group serotype summarization.

This repository is designed for a livestock-associated *E. coli* comparative genomics workflow. The active analysis directory remains:

```text
~/epi/marker_screen/
```

This GitHub repository stores scripts and documentation only. It should not store downloaded genomes, large metadata files, ECtyper result tables, Prokka outputs, Panaroo outputs, CheckM2 outputs, or phylogenetic outputs.

---

## Main scripts

| Script | Purpose |
|---|---|
| `scripts/00_run_ectyper_data5.sh` | Original ECtyper batch workflow: download assemblies from NCBI, run ECtyper, save outputs, and clean temporary genomes after each batch. |
| `scripts/01_summarize_ectyper_serotypes_by_source_group.R` | R summary workflow: join ECtyper output to metadata and summarize serotypes, O-types, and H-types by source group. |

---

## Biological goal

The immediate goal is to assign ECtyper serotypes to a curated set of livestock-associated *E. coli* genomes.

The downstream goal is to select approximately 2,000 genomes for:

1. Genome QC.
2. Contamination screening.
3. Filtering genomes with contamination `<1%`.
4. Pangenome analysis.
5. Core genome alignment.
6. Core genome SNP calling.
7. Phylogenetic analysis.

---

## Input file expected by Script 00

```text
~/epi/marker_screen/data5.csv
```

Required column:

```text
Assembly
```

Useful metadata columns for Script 01:

```text
source_group
Collection_year
Collection_year_group
Isolation source
Host
N50
Contigs
```

---

## Output files from Script 00

```text
~/epi/marker_screen/ectyper_results/ectyper_output_all.tsv
~/epi/marker_screen/ectyper_results/ectyper_download_manifest_all.tsv
```

If 18,560 genomes are completed, expected line count is:

```text
18,561 lines
```

That means:

```text
18,560 genomes + 1 header
```

---

## Why ECtyper?

A custom two-marker BLAST screen was initially considered, but the marker FASTA did not pass positive-control validation. ECtyper is preferred because it provides a standardized *E. coli* serotyping workflow with O-type, H-type, serotype, species, MASH similarity, QC, and evidence fields.

---

## Why batch processing?

Batching prevents the workflow from depending on one massive download. It also keeps disk use manageable.

Each batch:

1. Downloads a subset of genomes.
2. Runs ECtyper.
3. Saves batch and combined outputs.
4. Deletes temporary genome files.
5. Moves to the next batch.

If a batch download fails, the batch name is saved to:

```text
~/epi/marker_screen/ectyper_logs/failed_batches.txt
```

---

## Why missing serotypes are retained

Some genomes may not receive a full serotype call. Script 01 does not drop them. Instead, it labels incomplete calls explicitly:

```text
No_serotype_call
O8:No_H_call
No_O_call:H7
```

This prevents hidden bias from removing genomes with incomplete typing.

