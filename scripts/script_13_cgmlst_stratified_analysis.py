#!/usr/bin/env python3
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform

# Set styling for publication-quality charts
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['savefig.dpi'] = 300

# Global Path Definitions
BASE_DIR = "/home/samuelajulo/epi/Ecoli_epi"
METADATA_CSV = "/home/samuelajulo/epi/data6a.csv"
CGMLST_TSV = f"{BASE_DIR}/results/chewbbaca/cgmlst/core_95/cgMLST95.tsv"
OUTPUT_DIR = f"{BASE_DIR}/results/cgmlst_stratified"

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("==========================================================")
print(" Starting Stratified cgMLST Profiling & Visualization")
print("==========================================================")

# 1. Load Datasets
print("[INFO] Loading datasets...")
if not os.path.exists(METADATA_CSV) or not os.path.exists(CGMLST_TSV):
    print("[ERROR] Missing required input files. Verify your paths.")
    sys.exit(1)

meta_df = pd.read_csv(METADATA_CSV)
cgmlst_df = pd.read_csv(CGMLST_TSV, sep="\t", index_col=0)

# Clean index names (remove spaces or line breaks if present)
cgmlst_df.index = cgmlst_df.index.str.strip()
meta_df['Assembly'] = meta_df['Assembly'].astype(str).str.strip()

# Align metadata to include only genomes present in the cgMLST matrix
meta_df = meta_df[meta_df['Assembly'].isin(cgmlst_df.index)]
print(f"[SUCCESS] Aligned matrices. Ready to analyze {len(meta_df)} assemblies.")

# Define Cohorts
stratifications = {
    "Production_System": {
        "column": "Production_system",
        "groups": ["Beef", "Dairy", "Cattle-unspecified"]
    },
    "Supply_Chain_Stage": {
        "column": "Supply_chain_stage",
        "groups": ["Pre-harvest", "Peri-harvest", "Post-harvest"]
    }
}

# 2. Distance Matrix Calculation Function
def calculate_pairwise_distances(matrix_df):
    """Calculates pairwise hamming distances, accounting for chewBBACA missing values (-1)"""
    arr = matrix_df.values
    
    # Custom distance metric function to handle standard calls vs missing values (-1)
    def custom_hamming(u, v):
        # Only compare loci where BOTH genomes have a valid non-negative allele designation
        valid_mask = (u >= 0) & (v >= 0)
        if not np.any(valid_mask):
            return 1.0
        # Calculate fraction of mismatches among valid loci
        mismatches = np.sum(u[valid_mask] != v[valid_mask])
        return mismatches / np.sum(valid_mask)
    
    distances = pdist(arr, metric=custom_hamming)
    # Convert fractional distance back to total absolute locus differences
    total_loci = arr.shape[1]
    absolute_distances = distances * total_loci
    return squareform(absolute_distances)

# 3. Process Cohorts Loop
all_distance_data = []

for strat_name, config in stratifications.items():
    print(f"\n--- Processing Stratification Type: {strat_name} ---")
    col = config["column"]
    
    # Verify column existence
    if col not in meta_df.columns:
        print(f"[WARN] Column '{col}' not found in metadata. Skipping.")
        continue
        
    for group in config["groups"]:
        # Standardize matching strings
        group_meta = meta_df[meta_df[col].astype(str).str.lower() == group.lower()]
        group_ids = group_meta['Assembly'].tolist()
        
        subset_matrix = cgmlst_df.loc[cgmlst_df.index.isin(group_ids)]
        genome_count = len(subset_matrix)
        
        print(f"  Cohort: {group} -> Found {genome_count} genomes matching on disk.")
        
        if genome_count < 2:
            print(f"  [WARN] Not enough genomes for {group} to calculate distances. Skipping.")
            continue
            
        # Compute and save matrix
        dist_matrix = calculate_pairwise_distances(subset_matrix)
        dist_df = pd.DataFrame(dist_matrix, index=subset_matrix.index, columns=subset_matrix.index)
        
        group_clean_name = group.lower().replace(" ", "_").replace("-", "_")
        matrix_out_path = f"{OUTPUT_DIR}/{group_clean_name}_pairwise_locus_distances.tsv"
        dist_df.to_csv(matrix_out_path, sep="\t")
        
        # Save flattened vector for population metrics comparison
        flat_distances = dist_matrix[np.triu_indices(dist_matrix.shape[0], k=1)]
        for d in flat_distances:
            all_distance_data.append({
                "Stratification": strat_name,
                "Cohort": group,
                "Locus Difference": d
            })

# 4. Generate & Save Charts
print("\n[INFO] Generating comparative population plots...")
plot_df = pd.DataFrame(all_distance_data)

if not plot_df.empty:
    for strat_name in stratifications.keys():
        strat_plot_data = plot_df[plot_df['Stratification'] == strat_name]
        if strat_plot_data.empty:
            continue
            
        plt.figure(figsize=(12, 7))
        # Generate Boxplot distributions of distances
        ax = sns.boxplot(
            data=strat_plot_data, 
            x="Cohort", 
            y="Locus Difference", 
            palette="Set2",
            hue="Cohort",
            legend=False
        )
        
        plt.title(f"Pairwise Locus Variations by {strat_name.replace('_', ' ')}")
        plt.xlabel("Stratified Cohort Group")
        plt.ylabel("Number of Allele Differences (3,229 Total Loci)")
        
        chart_clean_name = strat_name.lower()
        plt.tight_layout()
        plt.savefig(f"{OUTPUT_DIR}/{chart_clean_name}_distance_distribution.png")
        plt.close()
        print(f"  [SUCCESS] Chart exported: {OUTPUT_DIR}/{chart_clean_name}_distance_distribution.png")

print("\n==========================================================")
print(" Stratified Master Profiling Completed Successfully!")
print(f" All files and plots saved to: {OUTPUT_DIR}")
print("==========================================================")
