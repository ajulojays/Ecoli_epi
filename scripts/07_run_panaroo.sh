#!/usr/bin/env bash
# ==============================================================================
# SCRIPT 07: Panaroo Pangenome Analysis and Core-Genome Alignment
# ==============================================================================

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BAKTA_OUT="/home/samuelajulo/epi/Ecoli_epi/results/bakta"
BASE_OUT="/home/samuelajulo/epi/Ecoli_epi/results"
PANAROO_OUT="${BASE_OUT}/panaroo"
GFF_LIST="${BASE_OUT}/panaroo_gff_list.txt"
LOG_DIR="${BASE_OUT}/logs"
TMP_DIR="/dev/shm/panaroo_tmp"

THREADS=64
CORE_THRESHOLD=0.95

mkdir -p "$BASE_OUT" "$LOG_DIR" "$TMP_DIR"

# ==============================================================================
# CHECK SOFTWARE
# ==============================================================================

command -v panaroo >/dev/null 2>&1 || {
    echo "[ERROR] Panaroo is not available in the active environment."
    exit 1
}

command -v mafft >/dev/null 2>&1 || {
    echo "[ERROR] MAFFT is not available in the active environment."
    exit 1
}

echo "[INFO] Panaroo version:"
panaroo --version

# ==============================================================================
# BUILD GFF3 INPUT LIST
# ==============================================================================

find "$BAKTA_OUT" \
    -mindepth 2 \
    -maxdepth 2 \
    -type f \
    -name "*.gff3" \
    -print0 |
    sort -z |
    tr '\0' '\n' > "$GFF_LIST"

GFF_COUNT=$(wc -l < "$GFF_LIST")

echo "[INFO] Bakta GFF3 files detected: $GFF_COUNT"

if (( GFF_COUNT < 2 )); then
    echo "[ERROR] At least two GFF3 files are required."
    exit 1
fi

# ==============================================================================
# CHECK COMPLETION
# ==============================================================================

PANGENOME_FILE="${PANAROO_OUT}/gene_presence_absence.csv"
CORE_ALIGNMENT="${PANAROO_OUT}/core_gene_alignment.aln"
FILTERED_ALIGNMENT="${PANAROO_OUT}/core_gene_alignment_filtered.aln"

if [[ -s "$PANGENOME_FILE" && \
      ( -s "$CORE_ALIGNMENT" || -s "$FILTERED_ALIGNMENT" ) ]]; then

    echo "[INFO] Panaroo pangenome and core alignment already exist."
    echo "[INFO] Pangenome table: $PANGENOME_FILE"

    if [[ -s "$FILTERED_ALIGNMENT" ]]; then
        echo "[INFO] Recommended alignment: $FILTERED_ALIGNMENT"
    else
        echo "[INFO] Core alignment: $CORE_ALIGNMENT"
    fi

    echo "[DONE] Nothing to run."
    exit 0
fi

# ==============================================================================
# PROTECT PARTIAL OR INCOMPATIBLE OUTPUT
# ==============================================================================

if [[ -d "$PANAROO_OUT" ]]; then
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    BACKUP_DIR="${PANAROO_OUT}_partial_${TIMESTAMP}"

    echo "[WARNING] Incomplete Panaroo output detected."
    echo "[INFO] Moving existing output to:"
    echo "       $BACKUP_DIR"

    mv "$PANAROO_OUT" "$BACKUP_DIR"
fi

mkdir -p "$PANAROO_OUT"

# ==============================================================================
# RUN PANAROO ONCE
# ==============================================================================

echo "[INFO] Starting Panaroo pangenome analysis."
echo "[INFO] Threads: $THREADS"
echo "[INFO] Core threshold: $CORE_THRESHOLD"
echo "[INFO] Temporary directory: $TMP_DIR"

panaroo \
    -i "$GFF_LIST" \
    -o "$PANAROO_OUT" \
    --clean-mode sensitive \
    --remove-invalid-genes \
    -a core \
    --aligner mafft \
    --core_threshold "$CORE_THRESHOLD" \
    --tmpdir "$TMP_DIR" \
    -t "$THREADS" \
    --verbose \
    2>&1 | tee "${LOG_DIR}/script07_panaroo.log"

# ==============================================================================
# VERIFY OUTPUT
# ==============================================================================

if [[ ! -s "$PANGENOME_FILE" ]]; then
    echo "[ERROR] Panaroo did not produce gene_presence_absence.csv."
    exit 1
fi

if [[ -s "$FILTERED_ALIGNMENT" ]]; then
    FINAL_ALIGNMENT="$FILTERED_ALIGNMENT"
elif [[ -s "$CORE_ALIGNMENT" ]]; then
    FINAL_ALIGNMENT="$CORE_ALIGNMENT"
else
    echo "[ERROR] Panaroo did not produce a core-genome alignment."
    exit 1
fi

echo
echo "[DONE] Panaroo analysis completed successfully."
echo "[INFO] Pangenome table: $PANGENOME_FILE"
echo "[INFO] Core alignment: $FINAL_ALIGNMENT"
echo "[INFO] Log file: ${LOG_DIR}/script07_panaroo.log"
