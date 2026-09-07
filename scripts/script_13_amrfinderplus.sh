#!/usr/bin/env bash

# ==============================================================================
# script_13_amrfinderplus.sh
#
# AMRFinderPlus screening pipeline for all 704 non-dairy cattle pre-harvest
# Escherichia coli genomes (Beef + Cattle-unspecified, Pre-harvest stage).
#
# Organism: Escherichia (includes point mutation detection)
# Environment: amrfinder_env
# ==============================================================================

set -euo pipefail

# ==========================================
# 1. GLOBAL CONFIG & PATHS
# ==========================================
BASE_DIR="/home/samuelajulo/epi/Ecoli_epi"
METADATA_CSV="/home/samuelajulo/epi/data6a.csv"
RESULTS_DIR="${BASE_DIR}/results/13_amrfinderplus"

AMRFINDER_ENV="amrfinder_env"
ORG="Escherichia"

# Threadripper tuning (JOBS * THREADS = 128)
JOBS=64
THREADS=2

EXPECTED_METADATA_COUNT=704
MIN_FASTA_BYTES=1000

LOG_DIR="${RESULTS_DIR}/logs"
TS="$(date +%F_%H%M%S)"
ARCHIVE_DIR="${RESULTS_DIR}/_archive_${TS}"

MASTER_REPORT="${RESULTS_DIR}/master_AMRFinder_report.tsv"
MANIFEST="${RESULTS_DIR}/assemblies_manifest.tsv"
JOBLOG="${RESULTS_DIR}/amrfinder_joblog_${TS}.tsv"
FAIL_LIST="${RESULTS_DIR}/failed_amrfinder_${TS}.txt"
HITCOUNT_TSV="${RESULTS_DIR}/amrfinder_hitcounts_${TS}.tsv"
MERGE_REPORT="${RESULTS_DIR}/master_build_report_${TS}.tsv"
TARGET_METADATA_CSV="${RESULTS_DIR}/target_metadata_beef_cattle_unspecified_preharvest.csv"

# ==========================================
# 2. CONDA INITIALIZATION
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
echo " 13: AMRFinderPlus Pipeline for 704 E. coli Genomes"
echo "=========================================================="
echo "[INFO] Base directory:     ${BASE_DIR}"
echo "[INFO] Metadata:           ${METADATA_CSV}"
echo "[INFO] Results directory:  ${RESULTS_DIR}"
echo "[INFO] Conda Environment:  ${AMRFINDER_ENV}"
echo "[INFO] Organism:           ${ORG}"
echo "[INFO] Parallel Jobs:      ${JOBS} (Threads/job: ${THREADS})"
echo ""

# ==========================================
# 3. PRE-FLIGHT CHECKS
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

if ! conda env list | awk '{print $1}' | grep -qx "$AMRFINDER_ENV"; then
    echo "[ERROR] Conda environment not found: $AMRFINDER_ENV"
    ERRORS=$((ERRORS + 1))
else
    echo "[OK] Conda environment exists: $AMRFINDER_ENV"

    if ! conda run -n "$AMRFINDER_ENV" which amrfinder >/dev/null 2>&1; then
        echo "[ERROR] 'amrfinder' binary not executable inside $AMRFINDER_ENV"
        ERRORS=$((ERRORS + 1))
    else
        echo "[OK] 'amrfinder' executable verified."
    fi

    if ! conda run -n "$AMRFINDER_ENV" which parallel >/dev/null 2>&1 && ! command -v parallel >/dev/null 2>&1; then
        echo "[ERROR] 'parallel' binary not found in PATH or $AMRFINDER_ENV"
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

# ==========================================
# 4. ARCHIVE PREVIOUS OUTPUTS
# ==========================================
mkdir -p "$RESULTS_DIR" "$LOG_DIR"

echo "[INFO] Archiving previous outputs if present..."
mkdir -p "$ARCHIVE_DIR" "$ARCHIVE_DIR/logs"

shopt -s nullglob
for f in "$RESULTS_DIR"/*.tsv "$RESULTS_DIR"/*.txt "$RESULTS_DIR"/*.log; do
    [[ -e "$f" ]] && mv -f "$f" "$ARCHIVE_DIR"/ || true
done
for f in "$LOG_DIR"/*; do
    [[ -e "$f" ]] && mv -f "$f" "$ARCHIVE_DIR/logs"/ || true
done
shopt -u nullglob

rm -rf "$LOG_DIR"
mkdir -p "$LOG_DIR"
: > "$FAIL_LIST"

# ==========================================
# 5. UPDATE AMRFINDERPLUS DATABASE
# ==========================================
echo "[INFO] Activating $AMRFINDER_ENV and updating AMRFinderPlus DB..."
conda activate "$AMRFINDER_ENV"
amrfinder -u

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

# Determine sample ID column
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

# Scan candidate FASTA directories
search_dirs = [
    base_dir / "assemblies",
    base_dir / "fasta",
    base_dir / "fna",
    base_dir / "results" / "bakta",
    base_dir / "gff_files"
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
# 7. RUN AMRFINDERPLUS IN PARALLEL
# ==========================================
echo "----------------------------------------------------------"
echo " Running AMRFinderPlus on 704 Assemblies"
echo "----------------------------------------------------------"

export OUT_DIR="$RESULTS_DIR"
export LOG_DIR
export THREADS
export ORG
export FAIL_LIST

cut -f1,2 "$MANIFEST" | tr '\n' '\0' | \
parallel -0 -j "$JOBS" --colsep '\t' --joblog "$JOBLOG" --eta '
  sample="{1}"
  fa="{2}"
  out="$OUT_DIR/${sample}.tsv"
  log="$LOG_DIR/${sample}.log"

  amrfinder -n "$fa" --organism "'"$ORG"'" --plus --threads "$THREADS" -o "$out" >"$log" 2>&1

  if [[ ! -s "$out" ]]; then
    echo "[FAIL] $sample (empty or missing output)" >> "'"$FAIL_LIST"'"
    exit 1
  fi
'

# ==========================================
# 8. BUILD MASTER REPORT & SUMMARY
# ==========================================
echo "----------------------------------------------------------"
echo " Building Master AMRFinderPlus Report"
echo "----------------------------------------------------------"

mapfile -t OUTS < <(cut -f1 "$MANIFEST" | awk '{print "'"$RESULTS_DIR"'/"$1".tsv"}')

# Determine dominant column count (header NF)
EXPECTED_NF="$(
  for f in "${OUTS[@]}"; do
    [[ -s "$f" ]] || continue
    awk -F'\t' 'NR==1{print NF; exit}' "$f"
  done | sort | uniq -c | sort -nr | awk 'NR==1{print $2}'
)"

if [[ -z "$EXPECTED_NF" ]]; then
    echo "[FATAL] Could not determine expected header NF from AMRFinderPlus outputs."
    exit 1
fi
echo "[INFO] Master report column count (excluding Sample): NF=$EXPECTED_NF"

# Header Template
HEADER_TSV=""
for f in "${OUTS[@]}"; do
  [[ -s "$f" ]] || continue
  nf="$(awk -F'\t' 'NR==1{print NF; exit}' "$f")"
  if [[ "$nf" == "$EXPECTED_NF" ]]; then
    HEADER_TSV="$f"
    break
  fi
done

TMP_MASTER="${MASTER_REPORT}.tmp"
echo -e "Sample\t$(head -n 1 "$HEADER_TSV")" > "$TMP_MASTER"

kept=0
skipped=0

for f in "${OUTS[@]}"; do
  [[ -s "$f" ]] || { ((skipped++)); continue; }
  nf="$(awk -F'\t' 'NR==1{print NF; exit}' "$f")"
  if [[ "$nf" != "$EXPECTED_NF" ]]; then
    ((skipped++))
    continue
  fi
  sample="$(basename "$f" .tsv)"
  tail -n +2 "$f" | awk -v s="$sample" 'BEGIN{OFS="\t"} {print s, $0}' >> "$TMP_MASTER"
  ((kept++))
done

mv -f "$TMP_MASTER" "$MASTER_REPORT"

# Hit counts summary
echo -e "Sample\tN_hits" > "$HITCOUNT_TSV"
for f in "${OUTS[@]}"; do
  [[ -s "$f" ]] || continue
  sample="$(basename "$f" .tsv)"
  n=$(( $(wc -l < "$f") - 1 )); (( n < 0 )) && n=0
  echo -e "${sample}\t${n}" >> "$HITCOUNT_TSV"
done

# Build Master Metrics Report
{
  echo -e "metric\tvalue"
  echo -e "organism\t$ORG"
  echo -e "expected_header_NF\t$EXPECTED_NF"
  echo -e "master_header_NF_including_Sample\t$((EXPECTED_NF+1))"
  echo -e "assemblies_total\t$TOTAL_ASM"
  echo -e "outputs_kept\t$kept"
  echo -e "outputs_skipped\t$skipped"
  echo -e "master_report_lines\t$(wc -l < "$MASTER_REPORT")"
} > "$MERGE_REPORT"

echo ""
echo "=========================================================="
echo " ✅ AMRFinderPlus Pipeline Completed Successfully!"
echo "=========================================================="
echo " Master Report:  $MASTER_REPORT"
echo " Manifest:       $MANIFEST"
echo " Hitcounts:      $HITCOUNT_TSV"
echo " Summary Build:  $MERGE_REPORT"
echo " Job Log:        $JOBLOG"
echo " Failed List:    $FAIL_LIST"
echo "=========================================================="
