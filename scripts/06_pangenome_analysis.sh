#!/usr/bin/env bash
# ==============================================================================
# SCRIPT 06: Pangenome Analysis (Bakta + Panaroo)
# ==============================================================================
set -euo pipefail

# Directories
PROJECT_DIR="$(pwd)"
ASSEMBLY_DIR="/home/samuelajulo/epi/marker_screen/results/03_O157H7_cattle_QC/fastas/raw"
BASE_OUT="${PROJECT_DIR}/results"
BAKTA_OUT="${BASE_OUT}/bakta"
PANAROO_OUT="${BASE_OUT}/panaroo"
PASS_LIST="/home/samuelajulo/epi/marker_screen/results/05_merge_QC_create_data7/tables/data7_O157H7_cattle_QC_pass_assemblies.txt"
GENOME_LIST="${BASE_OUT}/genome_list.txt"

mkdir -p "$BAKTA_OUT" "$PANAROO_OUT"

# Configuration
BAKTA_CPUS=8
PANAROO_CPUS=16
JOBS=8
export BAKTA_DB="$HOME/bakta_db/db"

echo "[INFO] Using Bakta DB: $BAKTA_DB"

# 1. BUILD TARGET LIST
awk -v dir="$ASSEMBLY_DIR" '{print dir "/" $1 ".fna"}' "$PASS_LIST" > "$GENOME_LIST"

# 2. RUN BAKTA (Parallelized)
export BAKTA_CPUS
export BAKTA_OUT
parallel -j "$JOBS" --bar --env BAKTA_DB --env BAKTA_CPUS --env BAKTA_OUT '
    genome="{}"
    prefix=$(basename "{}" .fna)
    outdir="${BAKTA_OUT}/${prefix}"
    mkdir -p "${outdir}"
    bakta --db "${BAKTA_DB}" --output "${outdir}" --prefix "${prefix}" \
          --threads "${BAKTA_CPUS}" --force "${genome}"
' :::: "$GENOME_LIST"

# 3. RUN PANAROO
echo "[INFO] Starting Panaroo..."
panaroo -i "${BAKTA_OUT}"/*/*.gff3 \
        -o "$PANAROO_OUT" \
        --clean-mode sensitive \
        --remove-invalid-genes \
        -t "$PANAROO_CPUS"

echo "[DONE] Pangenome analysis complete."
