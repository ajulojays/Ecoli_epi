#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# Fail early if the two manuscript-critical inputs are absent.
for f in inputs/metadata_704.csv inputs/cgMLST95.tsv; do
    if [[ ! -f "$f" ]]; then
        echo "[ERROR] Missing required input: $f" >&2
        echo >&2
        echo "First run ONE of:" >&2
        echo "  bash collect_inputs.sh /path/to/full/Ecoli_epi" >&2
        echo "or, if the full project is unavailable:" >&2
        echo "  bash recover_inputs_from_rescue.sh" >&2
        exit 2
    fi
done

# Prefer an already-installed Conda/Mamba because it avoids distro ensurepip issues.
if command -v micromamba >/dev/null 2>&1; then
    ENV_NAME="ecoli_persistence"
    if ! micromamba env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        echo "[INFO] Creating micromamba environment: $ENV_NAME"
        micromamba create -y -n "$ENV_NAME" -c conda-forge python=3.11 numpy pandas scipy
    fi
    RUNNER=(micromamba run -n "$ENV_NAME" python)

elif command -v conda >/dev/null 2>&1; then
    ENV_NAME="ecoli_persistence"
    if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        echo "[INFO] Creating conda environment: $ENV_NAME"
        conda create -y -n "$ENV_NAME" -c conda-forge python=3.11 numpy pandas scipy
    fi
    RUNNER=(conda run -n "$ENV_NAME" python)

else
    PYTHON="${PYTHON:-python3}"
    if [[ ! -d .venv ]]; then
        echo "[INFO] No Conda/Mamba detected; creating local Python venv..."
        if ! "$PYTHON" -m venv .venv; then
            echo >&2
            echo "[ERROR] Python venv support is missing." >&2
            echo "On Ubuntu/Debian run:" >&2
            echo "  sudo apt update && sudo apt install -y python3-venv" >&2
            echo "Then rerun: bash run.sh" >&2
            exit 3
        fi
    fi
    source .venv/bin/activate
    python -m pip install --upgrade pip >/dev/null
    python -m pip install -r requirements.txt
    RUNNER=(python)
fi

echo "[INFO] Running HC10/HC20/HC50 persistence-by-clade analysis..."
"${RUNNER[@]}" 15_cgmlst_persistence_clade.py \
  --input-dir inputs \
  --output-dir outputs

echo
echo "[SUCCESS] Main outputs:"
echo "  outputs/cluster_temporal_summary_HC10_HC20_HC50.csv"
echo "  outputs/clade_temporal_genome_burden.csv"
echo "  outputs/clade_persistence_enrichment_cluster_level.csv"
echo "  outputs/itol/"
