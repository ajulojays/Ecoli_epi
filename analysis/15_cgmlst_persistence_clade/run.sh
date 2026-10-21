#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

PYTHON="${PYTHON:-python3}"

if [[ ! -d .venv ]]; then
    echo "[INFO] Creating local Python environment..."
    "$PYTHON" -m venv .venv
fi

source .venv/bin/activate
python -m pip install --upgrade pip >/dev/null
python -m pip install -r requirements.txt

echo "[INFO] Running HC10/HC20/HC50 persistence-by-clade analysis..."
python 15_cgmlst_persistence_clade.py \
  --input-dir inputs \
  --output-dir outputs

echo
echo "[SUCCESS] Main outputs:"
echo "  outputs/cluster_temporal_summary_HC10_HC20_HC50.csv"
echo "  outputs/clade_temporal_genome_burden.csv"
echo "  outputs/clade_persistence_enrichment_cluster_level.csv"
echo "  outputs/itol/"
