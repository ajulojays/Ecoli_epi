#!/usr/bin/env bash
# ==============================================================================
# SCRIPT 08: IQ-TREE 2 Phylogenetic Inference
# ==============================================================================
set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# Input from your previous Panaroo run
INPUT_FILE="/home/samuelajulo/epi/Ecoli_epi/results/panaroo/core_gene_alignment.aln"
OUTDIR="/home/samuelajulo/epi/Ecoli_epi/results/iqtree"
LOG_DIR="/home/samuelajulo/epi/Ecoli_epi/results/logs"

mkdir -p "$OUTDIR" "$LOG_DIR"

CORES=64
MEMORY="180G"

# ==============================================================================
# RUN INFERENCE
# ==============================================================================
echo "[INFO] Starting IQ-TREE 2: Model Selection + ML Tree + 1000 UFBoot..."
echo "[INFO] Input file: $INPUT_FILE"
echo "[INFO] Utilizing $CORES cores and $MEMORY memory limit."

# NOTE: +ASC removed because Panaroo alignments contain constant/invariant sites.
# If the run takes too long, change 'MFP' to 'GTR+F+G4' to skip model testing.
iqtree2 -s "$INPUT_FILE" \
        -pre "${OUTDIR}/Ecoli_core_panaroo" \
        -st DNA \
        -m MFP \
        -bb 1000 \
        -bnni \
        -nt "$CORES" \
        -mem "$MEMORY" \
        2>&1 | tee "${LOG_DIR}/script08_iqtree.log"

echo "[DONE] Analysis complete. Check tree file at: ${OUTDIR}/Ecoli_core_panaroo.treefile"
