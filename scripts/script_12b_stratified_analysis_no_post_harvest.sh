#!/bin/bash
# script_12b_stratified_analysis_no_post_harvest.sh - Full Production Stratified Pangenome & Phylogeny Pipeline (Excluding Post-Harvest)
set -e

# ==========================================
# 1. GLOBAL PATHS & ENVIRONMENT DEFINITIONS
# ==========================================
BASE_DIR="/home/samuelajulo/epi/Ecoli_epi"
METADATA_CSV="/home/samuelajulo/epi/data6a.csv"
BAKTA_DIR="${BASE_DIR}/results/bakta"
RESULTS_DIR="${BASE_DIR}/results/stratified_analysis_without_post_harvest" # Dedicated outputs directory

PANAROO_ENV="panaroo_env"
IQTREE_ENV="iqtree_env"

THREADS=64
CORE_THRESHOLD=0.95

# Ensure Conda is initialized for the subshells
source ~/miniforge3/etc/profile.d/conda.sh

echo "=========================================================="
echo " Running Deep Pre-Flight Path & Dependency Checks..."
echo "=========================================================="
ERRORS=0

# Check Metadata File & Required Columns
if [ ! -f "$METADATA_CSV" ]; then
    echo "  [ERROR] Metadata file not found at: $METADATA_CSV"
    ERRORS=$((ERRORS + 1))
else
    echo "  [OK] Metadata file found."
    if ! python3 -c "import pandas as pd; df=pd.read_csv('$METADATA_CSV'); exit(0 if 'Production_system' in df.columns and 'Supply_chain_stage' in df.columns else 1)" 2>/dev/null; then
        echo "  [ERROR] Required columns 'Production_system' or 'Supply_chain_stage' are missing from $METADATA_CSV"
        ERRORS=$((ERRORS + 1))
    else
        echo "  [OK] Required metadata columns verified."
    fi
fi

# Check Bakta Output Directory
if [ ! -d "$BAKTA_DIR" ]; then
    echo "  [ERROR] Source Bakta directory does not exist: $BAKTA_DIR"
    ERRORS=$((ERRORS + 1))
else
    GFF_COUNT=$(find "$BAKTA_DIR" -mindepth 2 -maxdepth 2 -type f -name "*.gff3" | wc -l || echo 0)
    if [ "$GFF_COUNT" -eq 0 ]; then
        echo "  [ERROR] Bakta directory contains zero '.gff3' files."
        ERRORS=$((ERRORS + 1))
    else
        echo "  [OK] Bakta directory verified (Found $GFF_COUNT total .gff3 files)."
    fi
fi

# Verify panaroo_env & Tool Dependencies
if ! conda env list | grep -q "$PANAROO_ENV"; then
    echo "  [ERROR] Conda environment '$PANAROO_ENV' could not be found."
    ERRORS=$((ERRORS + 1))
else
    echo "  [OK] Conda environment '$PANAROO_ENV' exists."
    if ! conda run -n "$PANAROO_ENV" which panaroo >/dev/null 2>&1; then
        echo "  [ERROR] 'panaroo' binary is NOT executable inside '$PANAROO_ENV'."
        ERRORS=$((ERRORS + 1))
    else
        echo "  [OK] 'panaroo' tool installation verified."
    fi
    if ! conda run -n "$PANAROO_ENV" which mafft >/dev/null 2>&1; then
        echo "  [ERROR] 'mafft' alignment tool is NOT executable inside '$PANAROO_ENV'."
        ERRORS=$((ERRORS + 1))
    else
        echo "  [OK] 'mafft' alignment tool installation verified."
    fi
fi

# Verify iqtree_env & Tool Dependencies
if ! conda env list | grep -q "$IQTREE_ENV"; then
    echo "  [ERROR] Conda environment '$IQTREE_ENV' could not be found."
    ERRORS=$((ERRORS + 1))
else
    echo "  [OK] Conda environment '$IQTREE_ENV' exists."
    if ! conda run -n "$IQTREE_ENV" which iqtree >/dev/null 2>&1; then
        echo "  [ERROR] 'iqtree' binary is NOT executable inside '$IQTREE_ENV'."
        ERRORS=$((ERRORS + 1))
    else
        echo "  [OK] 'iqtree' tool installation verified."
    fi
fi

if [ "$ERRORS" -gt 0 ]; then
    echo -e "\n[FATAL] Pre-flight checks failed. Review errors above before running."
    exit 1
else
    echo -e "\n[SUCCESS] All checks passed! Launching production pipeline (excluding post-harvest)...\n"
fi

# ==========================================
# 2. CORE WORKFLOW LOOP
# ==========================================
COHORT_LIST=("Beef" "Dairy" "Cattle-unspecified")

for GROUP in "${COHORT_LIST[@]}"; do
    GRP_DIR=$(echo "$GROUP" | tr '[:upper:]' '[:lower:]' | tr ' ' '_' | tr '-' '_')
    
    GRP_OUT="${RESULTS_DIR}/${GRP_DIR}"
    PANAROO_OUT="${GRP_OUT}/panaroo"
    IQTREE_OUT="${GRP_OUT}/iqtree"
    GFF_LIST_FILE="${GRP_OUT}/panaroo_gff_list.txt"
    
    echo "----------------------------------------------------------"
    echo " PROCESSING COHORT: $GROUP (Excluding Post-Harvest Mode)"
    echo "----------------------------------------------------------"
    
    mkdir -p "$PANAROO_OUT"
    mkdir -p "$IQTREE_OUT"
    
    # Subset ALL sample IDs matching this specific production group, filtering out post-harvest stages
    echo "[INFO] Extracting assembly accession list for $GROUP (excluding post-harvest)..."
    SAMPLE_LIST=$(python3 - <<PYEOF
import pandas as pd
df = pd.read_csv("$METADATA_CSV")

id_col = 'Assembly'

subset = df[df['Production_system'].str.lower() == "$GROUP".lower()].copy()
if 'Supply_chain_stage' in subset.columns:
    subset = subset[subset['Supply_chain_stage'].str.lower().str.strip() != 'post-harvest']

for idx in subset[id_col].dropna():
    print(str(idx).strip())
PYEOF
)
    
    TOTAL_METADATA_COUNT=$(echo "$SAMPLE_LIST" | grep -v '^$' | wc -l || echo 0)
    if [ "$TOTAL_METADATA_COUNT" -eq 0 ]; then
        echo "[WARN] No samples matched the group $GROUP (excluding post-harvest) in metadata. Skipping."
        continue
    fi
    echo "[INFO] Isolated $TOTAL_METADATA_COUNT targeted entries. Generating file inventory..."

    # Generate the file path tracking list for Panaroo
    > "$GFF_LIST_FILE"
    while read -r SAMPLE_ID; do
        [ -z "$SAMPLE_ID" ] && continue
        SAMPLE_ID=$(echo "$SAMPLE_ID" | tr -d '\r')
        SRC_GFF="${BAKTA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}.gff3"
        
        if [ -f "$SRC_GFF" ]; then
            echo "$SRC_GFF" >> "$GFF_LIST_FILE"
        fi
    done <<< "$SAMPLE_LIST"

    MAPPED_GFF_COUNT=$(wc -l < "$GFF_LIST_FILE" || echo 0)
    echo "[INFO] Successfully mapped $MAPPED_GFF_COUNT / $TOTAL_METADATA_COUNT annotations on disk."

    if [ "$MAPPED_GFF_COUNT" -le 1 ]; then
        echo "[ERROR] Found $MAPPED_GFF_COUNT valid paths. Panaroo requires >= 2 files. Skipping."
        continue
    fi

    # Run Panaroo Pangenome & Core Gene Alignment Generation
    echo "[INFO] Activating ${PANAROO_ENV}..."
    conda activate "$PANAROO_ENV"
    
    echo "[INFO] Running Panaroo Core Genome Alignment..."
    panaroo -i "$GFF_LIST_FILE" \
            -o "$PANAROO_OUT" \
            --clean-mode sensitive \
            --remove-invalid-genes \
            -a core \
            --aligner mafft \
            --core_threshold "$CORE_THRESHOLD" \
            -t "$THREADS"
            
    # Pinpoint the final target alignment matrix
    if [ -s "${PANAROO_OUT}/core_gene_alignment_filtered.aln" ]; then
        CORE_ALIGNMENT="${PANAROO_OUT}/core_gene_alignment_filtered.aln"
    else
        CORE_ALIGNMENT="${PANAROO_OUT}/core_gene_alignment.aln"
    fi

    if [ ! -s "$CORE_ALIGNMENT" ]; then
        echo "[ERROR] Panaroo failed to output a core alignment file for $GROUP. Skipping tree step."
        continue
    fi
    echo "[SUCCESS] Found alignment file: $CORE_ALIGNMENT"

    # Run IQ-TREE Profiling (ModelFinder Plus enabled)
    echo "[INFO] Activating ${IQTREE_ENV}..."
    conda activate "$IQTREE_ENV"
    
    echo "[INFO] Running IQ-TREE Inference with ModelFinder..."
    cd "$IQTREE_OUT"
    
    iqtree -s "$CORE_ALIGNMENT" \
           -m MFP \
           -bb 1000 \
           -nt AUTO \
           -ntmax "$THREADS" \
           -pre "${GRP_DIR}_core_phylo"
           
    echo "[SUCCESS] Finished production workflow for cohort: $GROUP (excluding post-harvest)"
done

echo -e "\n=========================================================="
echo " Stratified Analysis Pipeline (Without Post-Harvest) Completed!"
echo " All cohorts fully resolved through Panaroo & IQ-TREE."
echo "=========================================================="
