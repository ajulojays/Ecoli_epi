#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script 04
# Purpose:
#   1. Read downloaded public assembly FASTA files from data6
#   2. Calculate FASTA-based assembly size, GC%, contig count, N50
#   3. Prepare normalized FASTA names for CheckM2
#   4. Run CheckM2 completeness and contamination estimation
# ============================================================

echo "=== SCRIPT 04: ASSEMBLY QC + CHECKM2 STARTED: $(date) ==="

# =========================
# Environment
# =========================
source ~/.bashrc || true
eval "$(conda shell.bash hook)"

CONDA_ENV_PATH="${CONDA_ENV_PATH:-$HOME/epi/envs/epi}"
conda activate "$CONDA_ENV_PATH"

# =========================
# Project paths
# =========================
ECOLI_EPI_WORKDIR="${ECOLI_EPI_WORKDIR:-$HOME/epi/marker_screen}"

DATA6_CSV="${DATA6_CSV:-$ECOLI_EPI_WORKDIR/results/02_O157H7_cattle_focus/tables/data6_cattle_O157H7_only.csv}"

RAW_FASTA_DIR="${RAW_FASTA_DIR:-$ECOLI_EPI_WORKDIR/results/03_O157H7_cattle_QC/fastas/raw}"

QC04_DIR="${QC04_DIR:-$ECOLI_EPI_WORKDIR/results/04_assembly_QC_checkm2}"

TABLE_DIR="$QC04_DIR/tables"
LOG_DIR="$QC04_DIR/logs"
CHECKM2_INPUT_DIR="$QC04_DIR/fastas/checkm2_input"
CHECKM2_OUT_DIR="$QC04_DIR/checkm2"

THREADS="${THREADS:-64}"
RUN_CHECKM2="${RUN_CHECKM2:-true}"

mkdir -p "$TABLE_DIR" "$LOG_DIR" "$CHECKM2_INPUT_DIR" "$CHECKM2_OUT_DIR"

# =========================
# Sanity checks
# =========================
if [[ ! -f "$DATA6_CSV" ]]; then
  echo "[ERROR] data6 CSV not found: $DATA6_CSV"
  exit 1
fi

if [[ ! -d "$RAW_FASTA_DIR" ]]; then
  echo "[ERROR] Raw FASTA directory not found: $RAW_FASTA_DIR"
  exit 1
fi

if [[ -z "$(find "$RAW_FASTA_DIR" -type f \( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' -o -name '*.fna.gz' -o -name '*.fa.gz' -o -name '*.fasta.gz' \) -print -quit)" ]]; then
  echo "[ERROR] No FASTA files found in: $RAW_FASTA_DIR"
  exit 1
fi

# =========================
# Python FASTA parser
# =========================
PY_SCRIPT="$QC04_DIR/calc_assembly_stats_and_prepare_checkm2.py"

cat > "$PY_SCRIPT" <<'PY'
#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import shutil
from pathlib import Path

FASTA_EXTS = (
    ".fna", ".fa", ".fasta",
    ".fna.gz", ".fa.gz", ".fasta.gz"
)

def open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def get_assembly_key(path):
    name = Path(path).name

    match = re.search(r"(GC[AF]_[0-9]+\.[0-9]+)", name)
    if match:
        return match.group(1)

    if name.endswith(".gz"):
        name = name[:-3]

    for ext in [".fasta", ".fna", ".fa"]:
        if name.endswith(ext):
            name = name[:-len(ext)]
            break

    return name

def n50(lengths):
    if not lengths:
        return 0

    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    running = 0

    for length in lengths:
        running += length
        if running >= half:
            return length

    return 0

def fasta_stats(path):
    lengths = []
    current_len = 0
    gc = 0
    atgc = 0
    n_bases = 0

    with open_maybe_gzip(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                seq = line.upper()
                current_len += len(seq)
                gc += seq.count("G") + seq.count("C")
                atgc += seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")
                n_bases += seq.count("N")

    if current_len > 0:
        lengths.append(current_len)

    total_bp = sum(lengths)
    gc_percent = (gc / atgc * 100) if atgc > 0 else 0

    return {
        "assembly_key": get_assembly_key(path),
        "original_fasta": str(path),
        "total_bp": total_bp,
        "gc_percent": round(gc_percent, 4),
        "n_contigs": len(lengths),
        "contig_n50": n50(lengths),
        "longest_contig": max(lengths) if lengths else 0,
        "n_bases": n_bases
    }

def write_normalized_fasta(src, dst):
    src = Path(src)
    dst = Path(dst)

    if dst.exists() or dst.is_symlink():
        dst.unlink()

    if str(src).endswith(".gz"):
        with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        os.symlink(src.resolve(), dst)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--stats-out", required=True)
    parser.add_argument("--checkm2-input-dir", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    checkm2_input_dir = Path(args.checkm2_input_dir)
    checkm2_input_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = []
    for path in input_dir.rglob("*"):
        if path.is_file() and str(path).endswith(FASTA_EXTS):
            fasta_files.append(path)

    fasta_files = sorted(fasta_files)

    if not fasta_files:
        raise SystemExit(f"No FASTA files found in: {input_dir}")

    rows = []
    seen = set()

    for fasta in fasta_files:
        row = fasta_stats(fasta)
        key = row["assembly_key"]

        if key in seen:
            raise SystemExit(f"Duplicate assembly key detected: {key}")

        seen.add(key)
        rows.append(row)

        normalized_fasta = checkm2_input_dir / f"{key}.fna"
        write_normalized_fasta(fasta, normalized_fasta)

    columns = [
        "assembly_key",
        "original_fasta",
        "total_bp",
        "gc_percent",
        "n_contigs",
        "contig_n50",
        "longest_contig",
        "n_bases"
    ]

    with open(args.stats_out, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row[col]) for col in columns) + "\n")

    print(f"[INFO] FASTA files processed: {len(rows)}")
    print(f"[INFO] Stats written to: {args.stats_out}")
    print(f"[INFO] CheckM2 input folder: {checkm2_input_dir}")

if __name__ == "__main__":
    main()
PY

chmod +x "$PY_SCRIPT"

# =========================
# FASTA QC
# =========================
FASTA_STATS="$TABLE_DIR/data6_fasta_assembly_stats.tsv"

python3 "$PY_SCRIPT" \
  --input-dir "$RAW_FASTA_DIR" \
  --stats-out "$FASTA_STATS" \
  --checkm2-input-dir "$CHECKM2_INPUT_DIR"

# =========================
# Run CheckM2
# =========================
if [[ "$RUN_CHECKM2" == "true" ]]; then

  if ! command -v checkm2 >/dev/null 2>&1; then
    echo "[ERROR] checkm2 not found in PATH."
    echo "[INFO] Install it in your conda environment, then rerun Script 04."
    exit 1
  fi

  checkm2 predict \
    --input "$CHECKM2_INPUT_DIR" \
    --extension fna \
    --output-directory "$CHECKM2_OUT_DIR" \
    --threads "$THREADS" \
    --force \
    > "$LOG_DIR/checkm2_run.log" 2>&1

else
  echo "[INFO] RUN_CHECKM2=false, skipping CheckM2."
fi

rm -f "$PY_SCRIPT"

echo "============================================================"
echo "[DONE] Script 04 complete."
echo "[INFO] FASTA QC table: $FASTA_STATS"
echo "[INFO] CheckM2 report: $CHECKM2_OUT_DIR/quality_report.tsv"
echo "============================================================"
