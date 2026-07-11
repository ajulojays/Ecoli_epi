#!/bin/bash
set -euo pipefail

# ============================================================
# SCRIPT 00: Batch ECtyper workflow for curated E. coli genomes
# ============================================================
#
# Purpose:
#   This script reads a curated metadata file, extracts NCBI assembly
#   accessions, downloads the corresponding genome FASTA files in batches,
#   runs ECtyper on each batch, saves the ECtyper output, and deletes
#   temporary genome files before moving to the next batch.
#
# Main input:
#   ~/epi/marker_screen/data5.csv
#
# Required column:
#   Assembly
#
# Main outputs:
#   ~/epi/marker_screen/ectyper_results/ectyper_output_all.tsv
#   ~/epi/marker_screen/ectyper_results/ectyper_download_manifest_all.tsv
#
# Key design:
#   The script processes genomes in batches so that disk use remains
#   controlled and failed NCBI downloads do not destroy already completed
#   results.
#
# Important:
#   ECtyper is run WITHOUT --verify because the upstream metadata curation
#   already restricted the dataset to confirmed E. coli genomes. Removing
#   --verify makes the workflow faster.
#
# ============================================================

echo "=== ECTYPER BATCH RUN STARTED: $(date) ==="

# ============================================================
# 1. Activate conda environment
# ============================================================

source ~/.bashrc
eval "$(conda shell.bash hook)"

# Allow users to override the conda environment path.
# Default used in the original project:
#   ~/epi/envs/epi
CONDA_ENV_PATH="${CONDA_ENV_PATH:-$HOME/epi/envs/epi}"

if [[ ! -d "$CONDA_ENV_PATH" ]]; then
  echo "[ERROR] Conda environment not found: $CONDA_ENV_PATH"
  echo "[HINT] Create it with:"
  echo "       conda env create -p $CONDA_ENV_PATH -f env/environment.yml"
  exit 1
fi

conda activate "$CONDA_ENV_PATH"

# ============================================================
# 2. Configuration
# ============================================================

# Allow users to override the large working directory.
# Default used in the original project:
#   ~/epi/marker_screen
WD="${ECOLI_EPI_WORKDIR:-$HOME/epi/marker_screen}"

METADATA_FILE="$WD/data5.csv"
ASSEMBLY_COLUMN="Assembly"

# Current full-run setting used in this project.
# 4000 genomes per batch gives fewer batches but larger temporary disk use.
BATCH_SIZE="${BATCH_SIZE:-4000}"

# Number of CPU threads for ECtyper.
THREADS="${THREADS:-128}"

# MAX_GENOMES controls test vs full run:
#   10   = quick test
#   100  = larger test
#   4000 = one full batch
#   0    = all genomes
MAX_GENOMES="${MAX_GENOMES:-0}"

BATCH_DIR="$WD/ectyper_batches"
WORK_DIR="$WD/ectyper_work"
RESULTS_DIR="$WD/ectyper_results"
LOGS="$WD/ectyper_logs"

mkdir -p "$BATCH_DIR" "$WORK_DIR" "$RESULTS_DIR" "$LOGS"

FINAL_OUTPUT="$RESULTS_DIR/ectyper_output_all.tsv"
FINAL_MANIFEST="$RESULTS_DIR/ectyper_download_manifest_all.tsv"
FAILED_BATCHES="$LOGS/failed_batches.txt"
MISSING_FASTA_LOG="$LOGS/missing_fastas_all_batches.tsv"

echo "[INFO] Working directory: $WD"
echo "[INFO] Metadata file: $METADATA_FILE"
echo "[INFO] Assembly column: $ASSEMBLY_COLUMN"
echo "[INFO] Batch size: $BATCH_SIZE"
echo "[INFO] Threads: $THREADS"
echo "[INFO] MAX_GENOMES: $MAX_GENOMES"

# ============================================================
# 3. Safety checks
# ============================================================

if [[ ! -f "$METADATA_FILE" ]]; then
  echo "[ERROR] Missing metadata file: $METADATA_FILE"
  exit 1
fi

if ! command -v datasets >/dev/null 2>&1; then
  echo "[ERROR] NCBI datasets command not found."
  echo "[HINT] Install with: conda install -c conda-forge -c bioconda ncbi-datasets-cli"
  exit 1
fi

if ! command -v ectyper >/dev/null 2>&1; then
  echo "[ERROR] ECtyper command not found."
  echo "[HINT] Install with: conda install -c conda-forge -c bioconda ectyper"
  exit 1
fi

# ============================================================
# 4. Clean old outputs
# ============================================================
#
# WARNING:
#   This deletes previous Script 00 outputs.
#   Do not run this main script to retry failed batches after a partial run.
#   Failed batches should be retried with a separate retry script.
#
# ============================================================

echo "[STEP 0] Cleaning previous ECtyper batch outputs..."

rm -f "$RESULTS_DIR"/batch_*_ectyper_output.tsv
rm -f "$RESULTS_DIR"/batch_*_download_manifest.tsv
rm -f "$RESULTS_DIR"/batch_*_ectyper_full_output.tar.gz
rm -f "$FINAL_OUTPUT"
rm -f "$FINAL_MANIFEST"
rm -f "$FAILED_BATCHES"
rm -f "$MISSING_FASTA_LOG"
rm -f "$BATCH_DIR"/batch_*
rm -f "$BATCH_DIR"/all_accessions.txt
rm -f "$BATCH_DIR"/all_accessions.full.txt
rm -rf "$WORK_DIR"/*

: > "$FAILED_BATCHES"
echo -e "batch\tassembly" > "$MISSING_FASTA_LOG"
echo -e "batch\tassembly\tstatus\tfna_file\tfile_size_bytes\tn_contigs\ttotal_bp" > "$FINAL_MANIFEST"

# ============================================================
# 5. Extract unique assembly accessions
# ============================================================

echo "[STEP 1] Extracting unique Assembly accessions from data5.csv..."

python3 - "$METADATA_FILE" "$ASSEMBLY_COLUMN" > "$BATCH_DIR/all_accessions.full.txt" <<'PY_ECTYPER_ACCESSIONS'
import csv
import sys

metadata_file = sys.argv[1]
assembly_col = sys.argv[2]

with open(metadata_file, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)

    if assembly_col not in reader.fieldnames:
        raise SystemExit(
            f"ERROR: Column '{assembly_col}' not found. Available columns: {reader.fieldnames}"
        )

    seen = set()

    for row in reader:
        acc = (row.get(assembly_col) or "").strip()

        if acc and acc not in seen:
            print(acc)
            seen.add(acc)
PY_ECTYPER_ACCESSIONS

TOTAL_GENOMES=$(wc -l < "$BATCH_DIR/all_accessions.full.txt")
echo "[INFO] Total unique assemblies in data5.csv: $TOTAL_GENOMES"

if [[ "$MAX_GENOMES" -gt 0 ]]; then
  echo "[INFO] Subset mode: using first $MAX_GENOMES genomes"
  head -n "$MAX_GENOMES" "$BATCH_DIR/all_accessions.full.txt" > "$BATCH_DIR/all_accessions.txt"
else
  echo "[INFO] Full mode: using all genomes"
  cp "$BATCH_DIR/all_accessions.full.txt" "$BATCH_DIR/all_accessions.txt"
fi

RUN_GENOMES=$(wc -l < "$BATCH_DIR/all_accessions.txt")
echo "[INFO] Genomes selected for this run: $RUN_GENOMES"

echo "[INFO] First few selected accessions:"
head "$BATCH_DIR/all_accessions.txt"

split -l "$BATCH_SIZE" "$BATCH_DIR/all_accessions.txt" "$BATCH_DIR/batch_"

echo "[INFO] Number of batches:"
ls "$BATCH_DIR"/batch_* | wc -l

# ============================================================
# 6. Sequential batch loop
# ============================================================

first_output=1

for BATCH_FILE in "$BATCH_DIR"/batch_*; do

  BATCH_NAME=$(basename "$BATCH_FILE")
  BATCH_WORK="$WORK_DIR/$BATCH_NAME"
  BATCH_FASTA_DIR="$BATCH_WORK/fasta_inputs"
  BATCH_ECTYPER_OUT="$BATCH_WORK/ectyper_out"
  BATCH_MANIFEST="$RESULTS_DIR/${BATCH_NAME}_download_manifest.tsv"
  BATCH_OUTPUT_COPY="$RESULTS_DIR/${BATCH_NAME}_ectyper_output.tsv"

  batch_start=$SECONDS

  echo "========================================"
  echo "[BATCH] Starting $BATCH_NAME at $(date)"
  echo "[BATCH] Assemblies:"
  wc -l "$BATCH_FILE"
  echo "========================================"

  rm -rf "$BATCH_WORK"
  mkdir -p "$BATCH_WORK" "$BATCH_FASTA_DIR" "$BATCH_ECTYPER_OUT"

  echo "[BATCH $BATCH_NAME] Downloading genomes..."

  if ! datasets download genome accession \
      --inputfile "$BATCH_FILE" \
      --include genome \
      --filename "$BATCH_WORK/temp_batch.zip" \
      >> "$LOGS/${BATCH_NAME}_datasets.log" 2>&1; then

      echo "[FAIL] datasets download failed: $BATCH_NAME"
      echo "$BATCH_NAME" >> "$FAILED_BATCHES"
      rm -rf "$BATCH_WORK"
      continue
  fi

  echo "[BATCH $BATCH_NAME] Unzipping..."
  unzip -q "$BATCH_WORK/temp_batch.zip" -d "$BATCH_WORK/temp_fastas"

  echo "[BATCH $BATCH_NAME] Creating manifest and ECtyper FASTA folder..."

  echo -e "batch\tassembly\tstatus\tfna_file\tfile_size_bytes\tn_contigs\ttotal_bp" > "$BATCH_MANIFEST"

  while read -r ACC; do

    fna_file=$(find "$BATCH_WORK/temp_fastas/ncbi_dataset/data/$ACC" -type f -name "*.fna" 2>/dev/null | head -n 1 || true)

    if [[ -z "${fna_file}" ]]; then
      fna_file=$(find "$BATCH_WORK/temp_fastas/ncbi_dataset/data" -type f -path "*$ACC*" -name "*.fna" 2>/dev/null | head -n 1 || true)
    fi

    if [[ -z "${fna_file}" ]]; then
      echo -e "$BATCH_NAME\t$ACC\tmissing_fasta\tNA\tNA\tNA\tNA" >> "$BATCH_MANIFEST"
      echo -e "$BATCH_NAME\t$ACC" >> "$MISSING_FASTA_LOG"
      continue
    fi

    file_size=$(stat -c%s "$fna_file")
    n_contigs=$(grep -c "^>" "$fna_file" || true)

    total_bp=$(awk '
      BEGIN {sum=0}
      /^>/ {next}
      {
        gsub(/[[:space:]]/, "", $0)
        sum += length($0)
      }
      END {print sum}
    ' "$fna_file")

    echo -e "$BATCH_NAME\t$ACC\tok\t$(basename "$fna_file")\t$file_size\t$n_contigs\t$total_bp" >> "$BATCH_MANIFEST"

    ln -s "$fna_file" "$BATCH_FASTA_DIR/${ACC}.fna"

  done < "$BATCH_FILE"

  tail -n +2 "$BATCH_MANIFEST" >> "$FINAL_MANIFEST"

  n_fasta=$(find "$BATCH_FASTA_DIR" -type l -name "*.fna" | wc -l)
  echo "[BATCH $BATCH_NAME] FASTA files prepared for ECtyper: $n_fasta"

  if [[ "$n_fasta" -eq 0 ]]; then
    echo "[FAIL] No FASTA files prepared for $BATCH_NAME"
    echo "$BATCH_NAME" >> "$FAILED_BATCHES"
    rm -rf "$BATCH_WORK"
    continue
  fi

  echo "[BATCH $BATCH_NAME] Running ECtyper without --verify..."

  ectyper_start=$SECONDS

  if ! ectyper \
      -i "$BATCH_FASTA_DIR" \
      -o "$BATCH_ECTYPER_OUT" \
      -c "$THREADS" \
      >> "$LOGS/${BATCH_NAME}_ectyper.log" 2>&1; then

      echo "[FAIL] ECtyper failed for $BATCH_NAME"
      echo "$BATCH_NAME" >> "$FAILED_BATCHES"
      rm -rf "$BATCH_WORK"
      continue
  fi

  ectyper_runtime=$((SECONDS - ectyper_start))
  echo "[BATCH $BATCH_NAME] ECtyper runtime: $((ectyper_runtime / 60)) min $((ectyper_runtime % 60)) sec"

  if [[ -f "$BATCH_ECTYPER_OUT/output.tsv" ]]; then
    cp "$BATCH_ECTYPER_OUT/output.tsv" "$BATCH_OUTPUT_COPY"
  else
    echo "[FAIL] ECtyper output.tsv missing for $BATCH_NAME"
    echo "$BATCH_NAME" >> "$FAILED_BATCHES"
    rm -rf "$BATCH_WORK"
    continue
  fi

  if [[ "$first_output" -eq 1 ]]; then
    awk -v batch="$BATCH_NAME" 'BEGIN{OFS="\t"} NR==1{print "batch",$0; next} {print batch,$0}' "$BATCH_OUTPUT_COPY" > "$FINAL_OUTPUT"
    first_output=0
  else
    awk -v batch="$BATCH_NAME" 'BEGIN{OFS="\t"} NR>1{print batch,$0}' "$BATCH_OUTPUT_COPY" >> "$FINAL_OUTPUT"
  fi

  tar -czf "$RESULTS_DIR/${BATCH_NAME}_ectyper_full_output.tar.gz" -C "$BATCH_WORK" "ectyper_out"

  echo "[BATCH $BATCH_NAME] Cleaning batch workspace..."
  rm -rf "$BATCH_WORK"

  batch_runtime=$((SECONDS - batch_start))
  echo "[BATCH] Finished $BATCH_NAME in $((batch_runtime / 60)) min $((batch_runtime % 60)) sec"
  echo "----------------------------------------"

done

echo "========================================"
echo "[FINAL] ECtyper run completed at $(date)"
echo "========================================"

echo "[INFO] Final ECtyper output:"
echo "$FINAL_OUTPUT"

if [[ -f "$FINAL_OUTPUT" ]]; then
  echo "[INFO] Rows in final ECtyper output:"
  wc -l "$FINAL_OUTPUT"
else
  echo "[WARN] No final ECtyper output file was created"
fi

echo "[INFO] Combined download manifest:"
echo "$FINAL_MANIFEST"
wc -l "$FINAL_MANIFEST"

echo "[INFO] Failed batches:"
cat "$FAILED_BATCHES" || true

echo "[INFO] Missing FASTAs:"
wc -l "$MISSING_FASTA_LOG"
head "$MISSING_FASTA_LOG"

echo "=== ECTYPER BATCH RUN COMPLETE: $(date) ==="
