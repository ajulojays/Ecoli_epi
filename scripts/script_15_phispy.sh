#!/usr/bin/env bash

# ==============================================================================
# script_15_phispy.sh
#
# High-speed Prokka annotation & ultra-parallel PhiSpy prophage identification
# pipeline for 704 non-dairy cattle pre-harvest E. coli assemblies.
#
# Environment: phage_env
# ==============================================================================

set -euo pipefail

# ==========================================
# 1. GLOBAL CONFIG & PATHS
# ==========================================
BASE_DIR="/home/samuelajulo/epi/Ecoli_epi"
METADATA_CSV="/home/samuelajulo/epi/data6a.csv"
RESULTS_DIR="${BASE_DIR}/results/15_phispy"

PROKKA_DIR="${RESULTS_DIR}/prokka_annotations"
PHISPY_DIR="${RESULTS_DIR}/phage_results"

PHAGE_ENV="phage_env"

# --- THREADRIPPER TUNING (128 CORES) ---
JOBS_PROKKA=64
THREADS_PROKKA=2

JOBS_PHISPY=128
THREADS_PHISPY=1

EXPECTED_METADATA_COUNT=704
MIN_FASTA_BYTES=1000

LOG_DIR="${RESULTS_DIR}/logs"
TS="$(date +%F_%H%M%S)"
ARCHIVE_DIR="${RESULTS_DIR}/_archive_${TS}"

MANIFEST="${RESULTS_DIR}/assemblies_manifest.tsv"
JOBLOG_PROKKA="${RESULTS_DIR}/prokka_joblog_${TS}.tsv"
JOBLOG_PHISPY="${RESULTS_DIR}/phispy_joblog_${TS}.tsv"
MASTER_SUMMARY="${RESULTS_DIR}/master_phage_summary.tsv"
TARGET_METADATA_CSV="${RESULTS_DIR}/target_metadata_beef_cattle_unspecified_preharvest.csv"

# ==========================================
# 2. SAFETY GUARDS FOR OUTPUT DIRECTORY
# ==========================================
if [[ -z "${RESULTS_DIR:-}" || "$RESULTS_DIR" == "/" ]]; then
    echo "[ERROR] RESULTS_DIR is unsafe: '$RESULTS_DIR' (refusing to proceed)"
    exit 1
fi

if [[ "$RESULTS_DIR" != "$BASE_DIR/"* ]]; then
    echo "[ERROR] RESULTS_DIR is not under BASE_DIR: '$RESULTS_DIR'"
    exit 1
fi

# ==========================================
# 3. CONDA INITIALIZATION
# ==========================================
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
fi

echo "=========================================================="
echo " 15: Prokka + PhiSpy Phage Pipeline for 704 E. coli Genomes"
echo "=========================================================="
echo "[INFO] Base directory:     ${BASE_DIR}"
echo "[INFO] Metadata:           ${METADATA_CSV}"
echo "[INFO] Results directory:  ${RESULTS_DIR}"
echo "[INFO] Conda Environment:  ${PHAGE_ENV}"
echo "[INFO] Prokka Parallel:    ${JOBS_PROKKA} jobs (${THREADS_PROKKA} threads/job)"
echo "[INFO] PhiSpy Parallel:    ${JOBS_PHISPY} jobs (${THREADS_PHISPY} thread/job)"
echo ""

# ==========================================
# 4. PRE-FLIGHT CHECKS
# ==========================================
echo "----------------------------------------------------------"
echo " Pre-flight checks"
echo "----------------------------------------------------------"

ERRORS=0

if [ ! -f "$METADATA_CSV" ]; then
    echo "[ERROR] Metadata file not found: $METADATA_CSV"
    ERRORS=$((ERRORS + 1))
else
    echo "[OK] Metadata file found."
fi

if ! conda env list | awk '{print $1}' | grep -qx "$PHAGE_ENV"; then
    echo "[ERROR] Conda environment not found: $PHAGE_ENV"
    ERRORS=$((ERRORS + 1))
else
    echo "[OK] Conda environment exists: $PHAGE_ENV"

    for tool in prokka phispy parallel; do
        if ! conda run -n "$PHAGE_ENV" which "$tool" >/dev/null 2>&1 && ! command -v "$tool" >/dev/null 2>&1; then
            echo "[ERROR] Binary '$tool' not found in $PHAGE_ENV or PATH"
            ERRORS=$((ERRORS + 1))
        else
            echo "[OK] '$tool' executable verified."
        fi
    done
fi

if [ "$ERRORS" -gt 0 ]; then
    echo "[FATAL] Pre-flight checks failed."
    exit 1
fi

echo "[SUCCESS] Pre-flight checks passed."
echo ""

# Activate Conda Env
conda activate "$PHAGE_ENV"

# ==========================================
# 5. ARCHIVE PREVIOUS OUTPUTS
# ==========================================
mkdir -p "$RESULTS_DIR" "$LOG_DIR" "$PROKKA_DIR" "$PHISPY_DIR"

if [ -d "$ARCHIVE_DIR" ]; then
    mkdir -p "$ARCHIVE_DIR" "$ARCHIVE_DIR/logs"
    shopt -s nullglob
    for f in "$RESULTS_DIR"/*.tsv "$RESULTS_DIR"/*.csv "$RESULTS_DIR"/*.log; do
        [[ -e "$f" ]] && mv -f "$f" "$ARCHIVE_DIR"/ || true
    done
    shopt -u nullglob
fi

# ==========================================
# 6. BUILD MANIFEST FOR EXACT 704 GENOMES
# ==========================================
echo "----------------------------------------------------------"
echo " Building Assembly Manifest (704 Target Genomes)"
echo "----------------------------------------------------------"

python3 - <<PY
import os
import re
import sys
import pandas as pd
from pathlib import Path

metadata_csv = "$METADATA_CSV"
base_dir = Path("$BASE_DIR")
manifest_file = "$MANIFEST"
target_metadata_csv = "$TARGET_METADATA_CSV"
expected_count = int("$EXPECTED_METADATA_COUNT")
min_bytes = int("$MIN_FASTA_BYTES")

df = pd.read_csv(metadata_csv)

id_col = None
for candidate in ["Assembly", "Isolate", "Strain", "Sample", "sample"]:
    if candidate in df.columns:
        id_col = candidate
        break
if id_col is None:
    id_col = df.columns[0]

def norm_text(x):
    if pd.isna(x): return ""
    return str(x).strip().lower()

def classify_stage(x):
    s = norm_text(x)
    if "post" in s: return "Post-harvest"
    if "peri" in s or "interface" in s: return "Peri-harvest"
    return "Pre-harvest"

def classify_prod(x):
    s = norm_text(x)
    if "beef" in s: return "Beef"
    if "dairy" in s: return "Dairy"
    return "Cattle-unspecified"

def normalize_id(x):
    x = str(x).strip()
    x = os.path.basename(x)
    x = re.sub(r"\.gz$", "", x)
    x = re.sub(r"\.(fasta|fna|fa|gff|gff3|txt|tsv|csv)$", "", x)
    m = re.search(r"GC[AF]_[0-9]+\.[0-9]+", x)
    return m.group(0) if m else x

df["stage"] = df["Supply_chain_stage"].apply(classify_stage)
df["prod"]  = df["Production_system"].apply(classify_prod)
df["sample_id"] = df[id_col].apply(lambda x: str(x).strip())
df["sample_key"] = df["sample_id"].apply(normalize_id)

subset = df[
    (df["stage"] == "Pre-harvest") &
    (df["prod"].isin(["Beef", "Cattle-unspecified"]))
].copy().drop_duplicates(subset=["sample_key"])

subset.to_csv(target_metadata_csv, index=False)
target_count = len(subset)
print(f"[INFO] Targeted Pre-Harvest Entries: {target_count}")

if target_count != expected_count:
    print(f"[ERROR] Expected {expected_count} genomes, found {target_count}.")
    sys.exit(1)

search_dirs = [
    base_dir / "assemblies",
    base_dir / "fasta",
    base_dir / "fna",
    base_dir / "shovill_assemblies",
    base_dir / "all_assemblies",
    base_dir / "results" / "bakta"
]

fasta_index = {}
for sdir in search_dirs:
    if not sdir.exists():
        continue
    for p in sdir.rglob("*"):
        if not p.is_file(): continue
        if p.suffix.lower() not in [".fna", ".fasta", ".fa"]: continue
        if p.stat().st_size < min_bytes: continue

        keys = {normalize_id(p.name), normalize_id(p.stem), normalize_id(p.parent.name)}
        m = re.search(r"GC[AF]_[0-9]+\.[0-9]+", str(p))
        if m: keys.add(m.group(0))

        for k in keys:
            if k and k not in fasta_index:
                fasta_index[k] = str(p)

manifest_entries = []
missing_keys = []

for _, row in subset.iterrows():
    s_key = row["sample_key"]
    fa_path = fasta_index.get(s_key)

    if fa_path:
        manifest_entries.append(f"{s_key}\t{fa_path}")
    else:
        missing_keys.append(s_key)

if missing_keys:
    print(f"[ERROR] Could not locate valid FASTA assemblies for {len(missing_keys)} samples.")
    print(f"[ERROR] First 5 missing keys: {missing_keys[:5]}")
    sys.exit(1)

with open(manifest_file, "w") as h:
    for line in manifest_entries:
        h.write(line + "\n")

print(f"[SUCCESS] Successfully built manifest with {len(manifest_entries)} FASTA files.")
PY

TOTAL_ASM=$(wc -l < "$MANIFEST" || echo 0)
echo "[INFO] FASTA assemblies mapped in manifest: $TOTAL_ASM"

if [ "$TOTAL_ASM" -ne "$EXPECTED_METADATA_COUNT" ]; then
    echo "[FATAL] Manifest count ($TOTAL_ASM) does not match expected ($EXPECTED_METADATA_COUNT)."
    exit 1
fi

# ==========================================
# 7. PHASE 1: HIGH-SPEED PROKKA ANNOTATION
# ==========================================
echo "----------------------------------------------------------"
echo " Phase 1: High-Speed Prokka Annotation (64 Jobs x 2 Threads)"
echo "----------------------------------------------------------"

export PROKKA_DIR
export LOG_DIR
export THREADS_PROKKA

cut -f1,2 "$MANIFEST" | tr '\n' '\0' | \
parallel -0 -j "$JOBS_PROKKA" --colsep '\t' --joblog "$JOBLOG_PROKKA" --eta '
    sample="{1}"
    fa="{2}"
    out_dir="$PROKKA_DIR/${sample}"
    log="$LOG_DIR/prokka_${sample}.log"

    if [ ! -f "$out_dir/${sample}.gbk" ]; then
        prokka --outdir "$out_dir" \
               --prefix "$sample" \
               --cpus "$THREADS_PROKKA" \
               --fast \
               --metagenome \
               --force \
               --quiet "$fa" >"$log" 2>&1
    fi
'

sync

# Verify Prokka Outputs
PROKKA_GENBANK_COUNT=$(find "$PROKKA_DIR" -type f -name "*.gbk" | wc -l || echo 0)
echo "[INFO] Prokka GenBank (.gbk) files generated: $PROKKA_GENBANK_COUNT / $EXPECTED_METADATA_COUNT"

if [ "$PROKKA_GENBANK_COUNT" -ne "$EXPECTED_METADATA_COUNT" ]; then
    echo "[WARN] Some Prokka annotations failed or are missing. Proceeding with available .gbk files."
fi

# ==========================================
# 8. PHASE 2: ULTRA-PARALLEL PHISPY
# ==========================================
echo "----------------------------------------------------------"
echo " Phase 2: Ultra-Parallel PhiSpy Prophage Detection (128 Jobs)"
echo "----------------------------------------------------------"

export PHISPY_DIR
export THREADS_PHISPY

find "$PROKKA_DIR" -type f -name "*.gbk" | tr '\n' '\0' | \
parallel -0 -j "$JOBS_PHISPY" --joblog "$JOBLOG_PHISPY" --eta '
    gbk="{}"
    sample=$(basename "$gbk" .gbk)
    out_dir="$PHISPY_DIR/${sample}"
    log="$LOG_DIR/phispy_${sample}.log"

    if [ ! -f "$out_dir/prophage.tbl" ]; then
        phispy -o "$out_dir" --threads "$THREADS_PHISPY" "$gbk" >"$log" 2>&1
    fi
'

sync

# ==========================================
# 9. PHASE 3: MASTER SUMMARY REPORT
# ==========================================
echo "----------------------------------------------------------"
echo " Phase 3: Building Master Phage Summary"
echo "----------------------------------------------------------"

echo -e "Sample\tphage_count" > "$MASTER_SUMMARY"

for tbl in "$PHISPY_DIR"/*/prophage.tbl; do
    if [ -f "$tbl" ]; then
        sample=$(basename "$(dirname "$tbl")")
        count=$(grep -v "^#" "$tbl" | wc -l || echo 0)
        echo -e "${sample}\t${count}" >> "$MASTER_SUMMARY"
    fi
done

TOTAL_PROCESSED=$(tail -n +2 "$MASTER_SUMMARY" | wc -l || echo 0)

echo ""
echo "=========================================================="
echo " ✅ Prokka + PhiSpy Pipeline Completed Successfully!"
echo "=========================================================="
echo " Results directory:  $RESULTS_DIR"
echo " Prokka Annotations: $PROKKA_DIR"
echo " PhiSpy Results:     $PHISPY_DIR"
echo " Master Summary:     $MASTER_SUMMARY ($TOTAL_PROCESSED samples)"
echo "=========================================================="
