#!/bin/bash
# script_11.sh - cgMLST Extraction & Metrics Summary Only (Skips QC)
set -e

BASE_DIR="/home/samuelajulo/epi/Ecoli_epi"
ALLELE_DIR="${BASE_DIR}/results/chewbbaca/allele_call/results"
ALLELE_TSV="${ALLELE_DIR}/results_alleles.tsv"
CGMLST_OUT="${BASE_DIR}/results/chewbbaca/cgmlst/core_95"

echo "=========================================================="
echo " Running cgMLST Core Extraction & Metrics Summary"
echo "=========================================================="

# ----------------------------------------------------------------
# Step 3: Extract the Core Genome (cgMLST at 95% using v3+ syntax)
# ----------------------------------------------------------------
echo -e "\n[INFO] Running Step 3: ExtractCgMLST (95% Threshold)..."
rm -rf "$CGMLST_OUT"
mkdir -p "$CGMLST_OUT"

chewBBACA.py ExtractCgMLST -i "$ALLELE_TSV" -o "$CGMLST_OUT" --t 0.95
echo "[SUCCESS] Core genome profile extracted successfully."

# ----------------------------------------------------------------
# Step 4: Run Run-Level Performance Summary
# ----------------------------------------------------------------
echo -e "\n[INFO] Running Step 4: Generating Metrics Summary..."
python3 - << 'PYEOF'
import pandas as pd
import sys

stats_path = "/home/samuelajulo/epi/Ecoli_epi/results/chewbbaca/allele_call/results/results_statistics.tsv"
alleles_path = "/home/samuelajulo/epi/Ecoli_epi/results/chewbbaca/allele_call/results/results_alleles.tsv"

try:
    stats = pd.read_csv(stats_path, sep="\t", index_col=0)
    alleles = pd.read_csv(alleles_path, sep="\t", index_col=0)
except Exception as e:
    print(f"[ERROR] Failed to read files: {e}")
    sys.exit(1)

total_genomes = len(stats)
total_loci = alleles.shape[1]

avg_alm = stats['ALM'].mean()
avg_asm = stats['ASM'].mean()
avg_missing = stats['LNF'].mean() + stats['INF'].mean()

print("\n" + "="*58)
print("            chewBBACA wgMLST Run Summary                ")
print("="*58)
print(f"Total Genomes Analyzed: {total_genomes}")
print(f"Total Loci in Schema:   {total_loci}")
print("-"*58)
print(f"Average exact allele matches (ALM) per genome: {avg_alm:.1f} ({ (avg_alm/total_loci)*100 :.1f}%)")
print(f"Average novel alleles found (ASM) per genome:  {avg_asm:.1f} ({ (avg_asm/total_loci)*100 :.1f}%)")
print(f"Average missing/incomplete loci per genome:    {avg_missing:.1f} ({ (avg_missing/total_loci)*100 :.1f}%)")
print("="*58 + "\n")
PYEOF

echo "=========================================================="
echo " Pipeline Finished Successfully!"
echo "=========================================================="
