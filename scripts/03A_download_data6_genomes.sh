#!/bin/bash
set -euo pipefail

# ============================================================
# SCRIPT 03A: Download genomes from data6 accession list
# ============================================================
#
# Purpose:
#   Download genome FASTA files for assemblies listed in data6.csv.
#
# Default input:
#   /home/samuelajulo/epi/data6.csv
#
# Required column:
#   Assembly
#
# Default output:
#   /home/samuelajulo/epi/marker_screen/results/03_O157H7_cattle_QC/fastas/raw
#
# Main outputs:
#   accession list
#   downloaded genome FASTAs
#   download manifest
#   failed accession list
#
# ============================================================

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV_PATH:-$HOME/epi/envs/epi}"

DATA6_CSV="${DATA6_CSV:-$HOME/epi/data6.csv}"

BASE_DIR="${BASE_DIR:-$HOME/epi/marker_screen/results/03_O157H7_cattle_QC}"
RAW_FASTA_DIR="$BASE_DIR/fastas/raw"
WORK_DIR="$BASE_DIR/download_work"
TABLE_DIR="$BASE_DIR/tables"
LOG_DIR="$BASE_DIR/logs"

ACCESSION_LIST="$TABLE_DIR/data6_assembly_accessions.txt"
MANIFEST="$TABLE_DIR/data6_download_manifest.tsv"
FAILED_LIST="$TABLE_DIR/data6_failed_downloads.txt"

BATCH_SIZE="${BATCH_SIZE:-200}"
THREADS="${THREADS:-16}"

mkdir -p "$RAW_FASTA_DIR" "$WORK_DIR" "$TABLE_DIR" "$LOG_DIR"

echo "=== SCRIPT 03A STARTED: $(date) ==="
echo "[INFO] DATA6_CSV: $DATA6_CSV"
echo "[INFO] BASE_DIR: $BASE_DIR"
echo "[INFO] RAW_FASTA_DIR: $RAW_FASTA_DIR"
echo "[INFO] BATCH_SIZE: $BATCH_SIZE"

if [[ ! -f "$DATA6_CSV" ]]; then
  echo "[ERROR] Input data6 CSV not found: $DATA6_CSV"
  exit 1
fi

if ! command -v datasets >/dev/null 2>&1; then
  echo "[ERROR] NCBI datasets command not found in current environment."
  echo "[HINT] Install with:"
  echo "       conda install -y -c conda-forge -c bioconda ncbi-datasets-cli"
  exit 1
fi

if ! command -v unzip >/dev/null 2>&1; then
  echo "[ERROR] unzip command not found."
  exit 1
fi

echo "[STEP 1] Extracting Assembly accessions from data6.csv..."

python3 - "$DATA6_CSV" "$ACCESSION_LIST" <<'PY'
import csv
import sys
from pathlib import Path

input_csv = Path(sys.argv[1])
output_txt = Path(sys.argv[2])

with input_csv.open(newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)

    if reader.fieldnames is None:
        raise SystemExit("ERROR: CSV has no header.")

    if "Assembly" not in reader.fieldnames:
        raise SystemExit(
            "ERROR: Required column 'Assembly' not found. "
            f"Available columns: {reader.fieldnames}"
        )

    seen = set()
    accessions = []

    for row in reader:
        acc = (row.get("Assembly") or "").strip()
        if acc and acc not in seen:
            accessions.append(acc)
            seen.add(acc)

output_txt.parent.mkdir(parents=True, exist_ok=True)

with output_txt.open("w") as out:
    for acc in accessions:
        out.write(acc + "\n")

print(f"[INFO] Unique accessions written: {len(accessions)}")
PY

TOTAL=$(wc -l < "$ACCESSION_LIST")
echo "[INFO] Total unique assemblies to download: $TOTAL"

if [[ "$TOTAL" -eq 0 ]]; then
  echo "[ERROR] No accessions found in $DATA6_CSV"
  exit 1
fi

echo -e "assembly\tstatus\tfasta_path\tfile_size_bytes\tn_contigs\ttotal_bp" > "$MANIFEST"
: > "$FAILED_LIST"

echo "[STEP 2] Splitting accession list into batches..."

rm -f "$WORK_DIR"/batch_*
split -l "$BATCH_SIZE" "$ACCESSION_LIST" "$WORK_DIR/batch_"

echo "[INFO] Number of batches:"
ls "$WORK_DIR"/batch_* | wc -l

echo "[STEP 3] Downloading genomes batch by batch..."

for BATCH_FILE in "$WORK_DIR"/batch_*; do
  BATCH_NAME=$(basename "$BATCH_FILE")
  BATCH_WORK="$WORK_DIR/${BATCH_NAME}_work"
  ZIP_FILE="$BATCH_WORK/${BATCH_NAME}.zip"
  EXTRACT_DIR="$BATCH_WORK/extracted"

  echo "========================================"
  echo "[BATCH] $BATCH_NAME started at $(date)"
  echo "[BATCH] Assemblies:"
  wc -l "$BATCH_FILE"
  echo "========================================"

  rm -rf "$BATCH_WORK"
  mkdir -p "$BATCH_WORK" "$EXTRACT_DIR"

  if ! datasets download genome accession \
      --inputfile "$BATCH_FILE" \
      --include genome \
      --filename "$ZIP_FILE" \
      > "$LOG_DIR/${BATCH_NAME}_datasets.log" 2>&1; then

    echo "[WARN] Batch download failed: $BATCH_NAME"
    cat "$BATCH_FILE" >> "$FAILED_LIST"
    rm -rf "$BATCH_WORK"
    continue
  fi

  echo "[BATCH $BATCH_NAME] Unzipping..."
  unzip -q "$ZIP_FILE" -d "$EXTRACT_DIR"

  echo "[BATCH $BATCH_NAME] Recovering FASTA files..."

  while read -r ACC; do
    [[ -z "$ACC" ]] && continue

    OUT_FASTA="$RAW_FASTA_DIR/${ACC}.fna"

    if [[ -s "$OUT_FASTA" ]]; then
      file_size=$(stat -c%s "$OUT_FASTA")
      n_contigs=$(grep -c "^>" "$OUT_FASTA" || true)
      total_bp=$(awk 'BEGIN{sum=0} /^>/{next} {gsub(/[[:space:]]/,"",$0); sum+=length($0)} END{print sum}' "$OUT_FASTA")
      echo -e "$ACC\talready_exists\t$OUT_FASTA\t$file_size\t$n_contigs\t$total_bp" >> "$MANIFEST"
      continue
    fi

    fna_file=$(find "$EXTRACT_DIR/ncbi_dataset/data/$ACC" -type f -name "*.fna" 2>/dev/null | head -n 1 || true)

    if [[ -z "$fna_file" ]]; then
      fna_file=$(find "$EXTRACT_DIR/ncbi_dataset/data" -type f -path "*$ACC*" -name "*.fna" 2>/dev/null | head -n 1 || true)
    fi

    if [[ -z "$fna_file" ]]; then
      echo "[WARN] Missing FASTA for $ACC"
      echo "$ACC" >> "$FAILED_LIST"
      echo -e "$ACC\tmissing_fasta\tNA\tNA\tNA\tNA" >> "$MANIFEST"
      continue
    fi

    cp "$fna_file" "$OUT_FASTA"

    file_size=$(stat -c%s "$OUT_FASTA")
    n_contigs=$(grep -c "^>" "$OUT_FASTA" || true)
    total_bp=$(awk 'BEGIN{sum=0} /^>/{next} {gsub(/[[:space:]]/,"",$0); sum+=length($0)} END{print sum}' "$OUT_FASTA")

    echo -e "$ACC\tok\t$OUT_FASTA\t$file_size\t$n_contigs\t$total_bp" >> "$MANIFEST"

  done < "$BATCH_FILE"

  echo "[BATCH $BATCH_NAME] Cleaning temporary files..."
  rm -rf "$BATCH_WORK"

  echo "[BATCH] $BATCH_NAME finished at $(date)"
done

echo "========================================"
echo "[FINAL] Download complete at $(date)"
echo "========================================"

echo "[INFO] FASTA files downloaded:"
find "$RAW_FASTA_DIR" -type f -name "*.fna" | wc -l

echo "[INFO] Manifest:"
echo "$MANIFEST"
wc -l "$MANIFEST"

echo "[INFO] Failed accession list:"
echo "$FAILED_LIST"
wc -l "$FAILED_LIST"

if [[ -s "$FAILED_LIST" ]]; then
  echo "[WARN] Some downloads failed. First few:"
  head "$FAILED_LIST"
else
  echo "[INFO] No failed downloads."
fi

echo "=== SCRIPT 03A COMPLETE: $(date) ==="
