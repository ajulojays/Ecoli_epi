#!/usr/bin/env bash
# ==============================================================================
# SCRIPT 08: High-Throughput Multi-Locus Sequence Typing (MLST)
# ==============================================================================

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================
BAKTA_OUT="/home/samuelajulo/epi/Ecoli_epi/results/bakta"
BASE_OUT="/home/samuelajulo/epi/Ecoli_epi/results"
MLST_OUT_DIR="${BASE_OUT}/mlst"
LOG_DIR="${BASE_OUT}/logs"
THREADS=64

# Output Files
RAW_MLST_OUT="${MLST_OUT_DIR}/mlst_raw_results.tsv"
CLEAN_MLST_OUT="${MLST_OUT_DIR}/mlst_summary_report.tsv"
MLST_LOG="${LOG_DIR}/script08_mlst.log"

# Create directories
mkdir -p "$MLST_OUT_DIR" "$LOG_DIR"

echo "=========================================================================="
echo " Starting MLST Profile Profiling for E. coli Assemblies"
echo "=========================================================================="

# ==============================================================================
# SANITY CHECKS
# ==============================================================================
command -v mlst >/dev/null 2>&1 || { echo "[ERROR] 'mlst' tool is not installed."; exit 1; }

# Find all FNA files inside the Bakta output folders
echo "[INFO] Scanning for assembly .fna files inside Bakta directories..."
mapfile -t FNA_FILES < <(find "$BAKTA_OUT" -mindepth 2 -maxdepth 2 -type f -name "*.fna" | sort)
FNA_COUNT=${#FNA_FILES[@]}

echo "[INFO] Found $FNA_COUNT assembly files."
if (( FNA_COUNT == 0 )); then
    echo "[ERROR] No assembly .fna files found in $BAKTA_OUT."
    exit 1
fi

# ==============================================================================
# RUN MLST PIPELINE
# ==============================================================================
echo "[INFO] Launching MLST across $THREADS threads..."

# Pass the array of files directly as arguments to mlst to bypass stdin issues
mlst --threads "$THREADS" --legacy --scheme ecoli "${FNA_FILES[@]}" > "$RAW_MLST_OUT" 2> "$MLST_LOG"

# ==============================================================================
# POST-PROCESSING: CREATE CLEAN RESULTS
# ==============================================================================
# The output has columns: [Filepath] [Scheme] [ST] [Allele1] [Allele2] ...
# For E. coli (Achtman), alleles are: adk fumC gyrB icd mdh purA recA
echo "[INFO] Post-processing output to format sample names..."

echo -e "Sample\tScheme\tST\tadk\tfumC\tgyrB\ticd\tmdh\tpurA\trecA" > "$CLEAN_MLST_OUT"

while IFS=$'\t' read -r filepath scheme st alleles; do
    # Extract only the file basename (e.g., sample_name.fna)
    filename=$(basename "$filepath")
    # Strip the .fna extension to get the clean Sample ID
    sample_id="${filename%.fna}"
    
    echo -e "${sample_id}\t${scheme}\t${st}\t${alleles}" >> "$CLEAN_MLST_OUT"
done < "$RAW_MLST_OUT"

# ==============================================================================
# SUMMARY METRICS
# ==============================================================================
TOTAL_TYPED=$(awk '$3 != "-" {print}' "$CLEAN_MLST_OUT" | wc -l)
(( TOTAL_TYPED-- )) # Subtract header row
NEW_STs=$(awk '$3 == "STnew" || $3 == "?" {print}' "$CLEAN_MLST_OUT" | wc -l)

echo "--------------------------------------------------------------------------"
echo " [SUCCESS] MLST execution completed!"
echo " Raw results: $RAW_MLST_OUT"
echo " Clean summary report: $CLEAN_MLST_OUT"
echo " Total genomes successfully assigned an ST: $TOTAL_TYPED / $FNA_COUNT"
echo " Genomes with novel/unassigned profiles (STnew/?): $NEW_STs"
echo "--------------------------------------------------------------------------"
