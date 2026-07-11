# Run commands

## Activate existing environment

```bash
conda activate ~/epi/envs/epi
```

## Run Script 00

```bash
cd ~/epi/marker_screen
bash ~/epi/Ecoli_epi/scripts/00_run_ectyper_data5.sh
```

## Check outputs

```bash
wc -l ~/epi/marker_screen/ectyper_results/ectyper_output_all.tsv
wc -l ~/epi/marker_screen/ectyper_results/ectyper_download_manifest_all.tsv
cat ~/epi/marker_screen/ectyper_logs/failed_batches.txt
```

## Run Script 01

```bash
Rscript ~/epi/Ecoli_epi/scripts/01_summarize_ectyper_serotypes_by_source_group.R
```

## View summary

```bash
column -s, -t ~/epi/marker_screen/ectyper_results/current_serotype_summary_by_source_group.csv | less -S
```

## View serotype counts

```bash
column -s, -t ~/epi/marker_screen/ectyper_results/current_serotype_counts_by_source_group.csv | less -S
```

