#!/usr/bin/env bash
set -euo pipefail

# Collect the manuscript-critical inputs from a full/restored Ecoli_epi project.
# Usage:
#   bash collect_inputs.sh /path/to/Ecoli_epi
# Default:
#   $HOME/epi/Ecoli_epi

PROJECT_DIR="${1:-$HOME/epi/Ecoli_epi}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="$HERE/inputs"
mkdir -p "$INPUT_DIR"

copy_required () {
    local src="$1"
    local dst="$2"
    if [[ ! -f "$src" ]]; then
        echo "[MISSING] $src" >&2
        return 1
    fi
    cp -f "$src" "$dst"
    echo "[COPIED] $src -> $dst"
}

echo "============================================================"
echo " Collecting cgMLST persistence/clade inputs"
echo " Project: $PROJECT_DIR"
echo "============================================================"

missing=0

copy_required \
  "$PROJECT_DIR/results/12c_all_preharvest_combined/target_metadata_704_preharvest_non_dairy.csv" \
  "$INPUT_DIR/metadata_704.csv" || missing=1

copy_required \
  "$PROJECT_DIR/results/chewbbaca/cgmlst/core_95/cgMLST95.tsv" \
  "$INPUT_DIR/cgMLST95.tsv" || missing=1

# Prefer the explicit clade CSV if present. The bundled iTOL strip is a valid fallback.
if [[ -f "$PROJECT_DIR/results/12d_clade_analysis/all_preharvest_704_clades.csv" ]]; then
    cp -f "$PROJECT_DIR/results/12d_clade_analysis/all_preharvest_704_clades.csv" \
          "$INPUT_DIR/all_preharvest_704_clades.csv"
    echo "[COPIED] 704 clade CSV"
else
    echo "[INFO] all_preharvest_704_clades.csv not found; bundled iTOL clade strip will be used."
fi

# Tree is optional for computation, but useful for immediate iTOL visualization.
if [[ -f "$PROJECT_DIR/results/12c_all_preharvest_combined/iqtree/all_preharvest_704_core_phylo.treefile" ]]; then
    cp -f "$PROJECT_DIR/results/12c_all_preharvest_combined/iqtree/all_preharvest_704_core_phylo.treefile" \
          "$INPUT_DIR/all_preharvest_704_core_phylo.treefile"
    echo "[COPIED] 704 IQ-TREE treefile"
else
    echo "[INFO] 704 IQ-TREE treefile not found locally. Analysis can still run; tree visualization will require the tree later."
fi

echo
if [[ "$missing" -ne 0 ]]; then
    echo "[STOP] One or more REQUIRED files are missing." >&2
    echo "They were not committed in the working GitHub repo and must come from the restored/rescue project." >&2
    exit 2
fi

echo "[SUCCESS] Required inputs collected in: $INPUT_DIR"
