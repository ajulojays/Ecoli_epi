#!/usr/bin/env bash

# ==============================================================================
# script_14_abricate.sh
#
# Abricate screening pipeline across multi-databases (resfinder, vfdb, plasmidfinder,
# ncbi, card) for all 704 non-dairy cattle pre-harvest E. coli assemblies.
#
# Environment: abricate_env
# ==============================================================================

set -euo pipefail

# ==========================================
# 1. GLOBAL CONFIG & PATHS
# ==========================================
BASE_DIR="/home/samuelajulo/epi/Ecoli_epi"
METADATA_CSV="/home/samuelajulo/epi/data6a.csv"
ABRICATE_BASE_OUT="${BASE_DIR}/results/14_abricate"

ABRICATE_ENV="abricate_env"
DATABASES=("resfinder" "vfdb" "plasmidfinder" "ncbi" "card")

# Parallel Jobs (Abricate is fast and lightweight per genome)
JOBS=32

EXPECTED_METADATA_COUNT=704
MIN_FASTA_BYTES=1000

LOG_DIR="${ABRICATE_BASE_OUT}/logs"
TS="$(date +%F_%H%M%S)"
ARCHIVE_DIR="${ABRICATE_BASE_OUT}/_archive_${TS}"

MANIFEST="${ABRICATE_BASE_OUT}/assemblies_manifest.tsv"
JOBLOG_BASE="${ABRICATE_BASE_OUT}/abricate_joblog"
TARGET_METADATA_CSV="${ABRICATE_BASE_OUT}/target_metadata_beef_cattle_unspecified_preharvest.csv"

# ==========================================
# 2. SAFETY GUARDS FOR OUTPUT DIRECTORY
# ==========================================
if [[ -z "${ABRICATE_BASE_OUT:-}" || "$ABRICATE_BASE_OUT" == "/" ]]; then
    echo "[ERROR] ABRICATE_BASE_OUT is unsafe: '$ABRICATE_BASE_OUT' (refusing to proceed)"
    exit 1
fi

if [[ "$ABRICATE_BASE_OUT" != "$BASE_DIR/"* ]]; then
    echo "[ERROR] ABRICATE_BASE_OUT is not under BASE_DIR: '$ABRICATE_BASE_OUT'"
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
echo " 14: Abricate Screening Pipeline for 704 E. coli Genomes"
echo "=========================================================="
echo "[INFO] Base directory:     ${BASE_DIR}"
echo "[INFO] Metadata:           ${METADATA_CSV}"
echo "[INFO] Results directory:  ${ABRICATE_BASE_OUT}"
echo "[INFO] Conda Environment:  ${ABRICATE_ENV}"
echo "[INFO] Databases:          ${DATABASES[*]}"
echo "[INFO] Parallel Jobs:      ${JOBS}"
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

if ! conda env list | awk '{print $1}' | grep -qx "$ABRICATE_ENV"; then
    echo "[ERROR] Conda environment not found: $ABRICATE_ENV"
    ERRORS=$((ERRORS + 1))
else
    echo "[OK] Conda environment exists: $ABRICATE_ENV"

    if ! conda run -n "$ABRICATE_ENV" which abricate >/dev/null 2>&1; then
        echo "[ERROR] 'abricate' binary not found in $ABRICATE_ENV"
        ERRORS=$((ERRORS + 1))
    else
        echo "[OK] 'abricate' binary verified."
    fi

    if ! command -v parallel >/dev/null 2>&1 && ! conda run -n "$ABRICATE_ENV" which parallel >/dev/null 2>&1; then
        echo "[ERROR] 'parallel' binary not found in PATH or $ABRICATE_ENV"
        ERRORS=$((ERRORS + 1))
    else
        echo "[OK] 'parallel' executable verified."
    fi
fi

if [ "$ERRORS" -gt 0 ]; then
    echo "[FATAL] Pre-flight checks failed."
    exit 1
fi

echo "[SUCCESS] Pre-flight checks passed."
echo ""

# Activate Conda Env
conda activate "$ABRICATE_ENV"

# ==========================================
# 5. ARCHIVE PREVIOUS OUTPUTS
# ==========================================
mkdir -p "$ABRICATE_BASE_OUT" "$LOG_DIR"

if [ -d "$ARCHIVE_DIR" ]; then
    mkdir -p "$ARCHIVE_DIR" "$ARCHIVE_DIR/logs"
    shopt -s nullglob
    for f in "$ABRICATE_BASE_OUT"/*.tsv "$ABRICATE_BASE_OUT"/*.txt "$ABRICATE_BASE_OUT"/*.log; do
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
# 7. RUN ABRICATE FOR EACH DATABASE
# ==========================================
for DB in "${DATABASES[@]}"; do
    echo "----------------------------------------------------"
    echo " Processing Database: $DB"
    echo "----------------------------------------------------"
    
    DB_OUT_DIR="$ABRICATE_BASE_OUT/$DB"
    DB_LOG_DIR="$LOG_DIR/$DB"
    mkdir -p "$DB_OUT_DIR" "$DB_LOG_DIR"
    
    JOB_LOG="${JOBLOG_BASE}_${DB}_${TS}.tsv"

    # Run Abricate in parallel reading directly from the manifest
    export DB_OUT_DIR
    export DB_LOG_DIR
    export DB

    cut -f1,2 "$MANIFEST" | tr '\n' '\0' | \
    parallel -0 -j "$JOBS" --colsep '\t' --joblog "$JOB_LOG" --eta '
        sample="{1}"
        fa="{2}"
        output="$DB_OUT_DIR/${sample}_${DB}.tsv"
        log="$DB_LOG_DIR/${sample}.log"

        if [ ! -f "$output" ]; then
            abricate --db "$DB" "$fa" > "$output" 2>"$log"
        fi
    '

    # Ensure background writes are flushed
    sync

    # ==========================================
    # 8. MERGE RESULTS INTO MASTER REPORT
    # ==========================================
    MASTER_REPORT="$DB_OUT_DIR/master_${DB}_report.tsv"
    echo "[INFO] Merging $DB results into $MASTER_REPORT..."

    FIRST_FILE=$(ls "$DB_OUT_DIR"/*.tsv 2>/dev/null | grep -v "master_" | head -n 1 || true)

    if [ -z "$FIRST_FILE" ]; then
        echo "[WARN] No individual result files found in $DB_OUT_DIR. Skipping merge."
        continue
    fi

    # Create Header with 'Sample' column prepended
    echo -e "Sample\t$(head -n 1 "$FIRST_FILE")" > "$MASTER_REPORT"

    # Append data lines
    for f in "$DB_OUT_DIR"/*.tsv; do
        filename=$(basename "$f")
        if [[ "$filename" != master_* ]]; then
            sample=$(echo "$filename" | sed "s/_${DB}.tsv//")
            tail -n +2 "$f" | sed "s/^/${sample}\t/" >> "$MASTER_REPORT"
        fi
    done

    echo "[SUCCESS] Master report created: $MASTER_REPORT"
done

echo ""
echo "=========================================================="
echo " ✅ Abricate Pipeline Completed Successfully!"
echo "=========================================================="
echo " Results directory: $ABRICATE_BASE_OUT"
echo " Manifest:          $MANIFEST"
echo "=========================================================="
