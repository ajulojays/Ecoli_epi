# Workflow notes

## Script 00

```text
scripts/00_run_ectyper_data5.sh
```

Script 00 is the original ECtyper batch workflow.

It does the following:

1. Activates the conda environment at `~/epi/envs/epi`.
2. Reads `~/epi/marker_screen/data5.csv`.
3. Extracts unique NCBI assembly accessions from the `Assembly` column.
4. Splits accessions into batch files.
5. Downloads each batch using NCBI `datasets`.
6. Unzips downloaded genome files.
7. Locates each assembly FASTA file.
8. Creates a manifest for each batch.
9. Creates a flat ECtyper FASTA input folder using symbolic links.
10. Runs ECtyper without `--verify`.
11. Saves the ECtyper batch output.
12. Appends batch output to the final combined table.
13. Compresses the full ECtyper output folder for each batch.
14. Deletes temporary downloaded genomes and workspace files.
15. Records failed batches.

## Script 01

```text
scripts/01_summarize_ectyper_serotypes_by_source_group.R
```

Script 01 summarizes ECtyper output by source group.

It does the following:

1. Reads `data5.csv`.
2. Reads `ectyper_output_all.tsv`.
3. Joins metadata and ECtyper results by `Assembly` / `Name`.
4. Cleans missing serotype calls.
5. Keeps genomes with incomplete serotype calls.
6. Writes a joined metadata table.
7. Writes serotype counts by source group.
8. Writes O-type counts by source group.
9. Writes H-type counts by source group.
10. Prints a terminal summary.

## Missing serotype labels

The R script uses explicit labels:

| Situation | Label |
|---|---|
| No O-type and no H-type | `No_serotype_call` |
| O-type present, H-type missing | `Otype:No_H_call` |
| O-type missing, H-type present | `No_O_call:Htype` |
| Complete call | `Otype:Htype` |


## Script 01 figure outputs

Script 01 now generates publication-ready figure outputs in:

```text
~/epi/marker_screen/ectyper_results/figures/
```

The figure logic is source-group specific.

For each source group:

1. The top 10 called serotypes are retained individually.
2. All other called serotypes are collapsed into `Other serotypes`.
3. Genomes without a serotype call are retained as `No serotype call`.

This produces:

```text
pie_top10_serotypes_<source_group>.pdf
pie_top10_serotypes_<source_group>.png
temporal_stacked_counts_top10_serotypes_<source_group>.pdf
temporal_stacked_counts_top10_serotypes_<source_group>.png
temporal_stacked_percent_top10_serotypes_<source_group>.pdf
temporal_stacked_percent_top10_serotypes_<source_group>.png
```

Multi-page combined PDFs are also created:

```text
pie_top10_serotypes_all_source_groups.pdf
temporal_stacked_counts_top10_serotypes_all_source_groups.pdf
temporal_stacked_percent_top10_serotypes_all_source_groups.pdf
```


## Script 01 figure update

The temporal stacked bar plots now group all collection years <=2000 as one temporal category.

Each temporal bar includes a vertical total genome count label:

```text
n=<number_of_genomes>
```

Script 01 also generates O157:H7-focused pie charts by source group. These compare:

```text
O157:H7
Non-O157:H7
No serotype call
```

The no-serotype group is retained separately to avoid incorrectly forcing untyped genomes into the non-O157:H7 category.


## Organized result folders

The working analysis directory now uses script-numbered result folders:

  results/00_ectyper_batch_run
  results/01_serotype_summary
  results/02_downstream_selection

Script 00 output belongs in results/00_ectyper_batch_run.

Script 01 output belongs in results/01_serotype_summary.

This prevents ECtyper batch files, summary tables, and figures from being mixed in one folder.
