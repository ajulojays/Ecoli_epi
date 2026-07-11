# Environment and setup guide

This guide explains how another user can clone, configure, and run the `Ecoli_epi` workflow from scratch.

The repository is code-only. It does not store downloaded genomes, ECtyper outputs, logs, pangenome files, CheckM2 outputs, Prokka annotations, or phylogenetic outputs.

Large files should remain in a separate working directory.

---

## 1. Clone the repository

```bash
git clone https://github.com/ajulojays/Ecoli_epi.git
cd Ecoli_epi
```

---

## 2. Recommended project layout

The recommended layout is:

```text
~/epi/
├── Ecoli_epi/              # GitHub repository: scripts and documentation only
├── envs/
│   └── epi/                # Conda environment
└── marker_screen/          # Large working analysis directory
    ├── data5.csv
    ├── ectyper_batches/
    ├── ectyper_work/
    ├── ectyper_results/
    └── ectyper_logs/
```

The GitHub repository should remain small. The working directory stores large analysis files.

---

## 3. Create the conda environment

Create the environment from the repository environment file:

```bash
conda env create -p ~/epi/envs/epi -f env/environment.yml
```

Activate it:

```bash
conda activate ~/epi/envs/epi
```

Confirm the environment works:

```bash
which ectyper
ectyper --version

which datasets
datasets --version

which Rscript
Rscript --version

which git
git --version
```

If using GitHub CLI:

```bash
which gh
gh --version
gh auth status
```

If GitHub CLI is not authenticated:

```bash
gh auth login
```

---

## 4. Required tools

The environment includes:

| Tool | Purpose |
|---|---|
| `ectyper` | Predicts *E. coli* O-type, H-type, serotype, species fields, and evidence fields. |
| `ncbi-datasets-cli` / `datasets` | Downloads genome assemblies from NCBI. |
| `unzip` | Extracts downloaded NCBI genome zip files. |
| `R` | Runs the serotype summary script. |
| `readr`, `dplyr`, `stringr`, `tidyr` | R packages used for joining, cleaning, and summarizing ECtyper output. |
| `git` | Version control. |
| `gh` | GitHub CLI for pushing or managing the repository. |

---

## 5. Input metadata

The main input file is:

```text
data5.csv
```

Place it in the large working directory:

```bash
mkdir -p ~/epi/marker_screen
cp /path/to/data5.csv ~/epi/marker_screen/data5.csv
```

Required columns:

| Column | Description |
|---|---|
| `Assembly` | NCBI assembly accession, for example `GCA_012836955.1`. |
| `source_group` | Curated source group, for example `Cattle`, `Poultry`, or `Swine`. |

Useful optional columns:

| Column | Description |
|---|---|
| `Collection_year` | Year of isolation or collection. |
| `Collection_year_group` | Binned year group. |
| `Isolation source` | Original NCBI isolation source text. |
| `Host` | Host metadata. |
| `N50` | Assembly N50. |
| `Contigs` | Number of contigs. |

Check the header:

```bash
head -n 1 ~/epi/marker_screen/data5.csv
```

---

## 6. Script 00: ECtyper batch workflow

Script:

```text
scripts/00_run_ectyper_data5.sh
```

Purpose:

1. Reads `data5.csv`.
2. Extracts unique NCBI assembly accessions from the `Assembly` column.
3. Splits accessions into batches.
4. Downloads genome FASTAs using NCBI `datasets`.
5. Unzips downloaded genome files.
6. Finds each `.fna` assembly FASTA.
7. Builds a per-batch download manifest.
8. Creates a flat ECtyper input folder using symbolic links.
9. Runs ECtyper without `--verify`.
10. Saves per-batch ECtyper output.
11. Appends each successful batch to one final table.
12. Compresses full ECtyper batch output.
13. Deletes temporary downloaded genomes after each batch.
14. Records failed batches.

---

## 7. Run a small test

Always run a small test before a full run.

```bash
cd ~/epi/Ecoli_epi

ECOLI_EPI_WORKDIR=~/epi/marker_screen \
CONDA_ENV_PATH=~/epi/envs/epi \
BATCH_SIZE=100 \
THREADS=8 \
MAX_GENOMES=10 \
bash scripts/00_run_ectyper_data5.sh
```

Expected outputs:

```text
~/epi/marker_screen/ectyper_results/ectyper_output_all.tsv
~/epi/marker_screen/ectyper_results/ectyper_download_manifest_all.tsv
```

For 10 genomes, expected line count is 11:

```text
10 genome rows + 1 header
```

Check:

```bash
wc -l ~/epi/marker_screen/ectyper_results/ectyper_output_all.tsv
wc -l ~/epi/marker_screen/ectyper_results/ectyper_download_manifest_all.tsv
```

---

## 8. Run the full ECtyper workflow

For a normal workstation:

```bash
cd ~/epi/Ecoli_epi

ECOLI_EPI_WORKDIR=~/epi/marker_screen \
CONDA_ENV_PATH=~/epi/envs/epi \
BATCH_SIZE=1000 \
THREADS=32 \
MAX_GENOMES=0 \
bash scripts/00_run_ectyper_data5.sh
```

For a large workstation:

```bash
cd ~/epi/Ecoli_epi

ECOLI_EPI_WORKDIR=~/epi/marker_screen \
CONDA_ENV_PATH=~/epi/envs/epi \
BATCH_SIZE=4000 \
THREADS=128 \
MAX_GENOMES=0 \
bash scripts/00_run_ectyper_data5.sh
```

`MAX_GENOMES=0` means use all genomes in `data5.csv`.

---

## 9. Runtime notes

Internet is required during the NCBI `datasets download` step.

Once a batch has downloaded, ECtyper itself does not need internet.

If a download fails, the failed batch is recorded in:

```text
~/epi/marker_screen/ectyper_logs/failed_batches.txt
```

Do not rerun Script 00 just to recover failed batches, because Script 00 cleans previous outputs at startup.

---

## 10. Output files from Script 00

Main combined outputs:

```text
~/epi/marker_screen/ectyper_results/ectyper_output_all.tsv
~/epi/marker_screen/ectyper_results/ectyper_download_manifest_all.tsv
```

Per-batch outputs:

```text
~/epi/marker_screen/ectyper_results/batch_aa_ectyper_output.tsv
~/epi/marker_screen/ectyper_results/batch_aa_download_manifest.tsv
~/epi/marker_screen/ectyper_results/batch_aa_ectyper_full_output.tar.gz
```

Logs:

```text
~/epi/marker_screen/ectyper_logs/batch_aa_datasets.log
~/epi/marker_screen/ectyper_logs/batch_aa_ectyper.log
~/epi/marker_screen/ectyper_logs/failed_batches.txt
~/epi/marker_screen/ectyper_logs/missing_fastas_all_batches.tsv
```

---

## 11. Script 01: Serotype summary by source group

Script:

```text
scripts/01_summarize_ectyper_serotypes_by_source_group.R
```

Purpose:

1. Reads `data5.csv`.
2. Reads `ectyper_output_all.tsv`.
3. Joins ECtyper output to metadata using `Name = Assembly`.
4. Keeps all typed genomes.
5. Cleans missing O-type, H-type, and serotype values.
6. Labels incomplete calls explicitly.
7. Summarizes serotypes by source group.
8. Summarizes O-types by source group.
9. Summarizes H-types by source group.
10. Writes joined per-genome metadata.

Run:

```bash
cd ~/epi/Ecoli_epi

ECOLI_EPI_WORKDIR=~/epi/marker_screen \
Rscript scripts/01_summarize_ectyper_serotypes_by_source_group.R
```

Outputs:

```text
~/epi/marker_screen/ectyper_results/current_ectyper_joined_metadata.csv
~/epi/marker_screen/ectyper_results/current_serotype_summary_by_source_group.csv
~/epi/marker_screen/ectyper_results/current_serotype_counts_by_source_group.csv
~/epi/marker_screen/ectyper_results/current_Otype_counts_by_source_group.csv
~/epi/marker_screen/ectyper_results/current_Htype_counts_by_source_group.csv
```

View summary:

```bash
column -s, -t ~/epi/marker_screen/ectyper_results/current_serotype_summary_by_source_group.csv | less -S
```

View serotype counts:

```bash
column -s, -t ~/epi/marker_screen/ectyper_results/current_serotype_counts_by_source_group.csv | less -S
```

---

## 12. Missing serotype handling

Some genomes may not receive a complete serotype call.

Script 01 does not drop these genomes.

It labels them as:

| Situation | Label |
|---|---|
| No O-type and no H-type | `No_serotype_call` |
| O-type present, H-type missing | `Otype:No_H_call` |
| O-type missing, H-type present | `No_O_call:Htype` |
| Complete call | `Otype:Htype` |

This prevents hidden bias from dropping genomes with incomplete typing.

---

## 13. Batch size guidance

| Computer type | Suggested `BATCH_SIZE` |
|---|---|
| Laptop or low disk | 250–500 |
| Normal workstation | 1000 |
| Large workstation | 2000–4000 |

Larger batches reduce overhead but increase temporary disk use and make download failures more costly.

---

## 14. Reproducibility checklist

Before running:

```bash
conda activate ~/epi/envs/epi
ectyper --version
datasets --version
Rscript --version
```

Check input:

```bash
head -n 1 ~/epi/marker_screen/data5.csv
```

Required columns:

```text
Assembly
source_group
```

Run test:

```bash
ECOLI_EPI_WORKDIR=~/epi/marker_screen \
CONDA_ENV_PATH=~/epi/envs/epi \
BATCH_SIZE=100 \
THREADS=8 \
MAX_GENOMES=10 \
bash scripts/00_run_ectyper_data5.sh
```

Only after the test passes should the full run be started.

