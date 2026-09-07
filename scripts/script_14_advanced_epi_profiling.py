#!/usr/bin/env python3
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
import scipy.cluster.hierarchy as sch
import networkx as nx

sns.set_theme(style="white", context="talk")
plt.rcParams['savefig.dpi'] = 300

BASE_DIR = "/home/samuelajulo/epi/Ecoli_epi"
METADATA_CSV = "/home/samuelajulo/epi/data6a.csv"
CGMLST_TSV = f"{BASE_DIR}/results/chewbbaca/cgmlst/core_95/cgMLST95.tsv"
OUTPUT_DIR = f"{BASE_DIR}/results/cgmlst_advanced_epi"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("==========================================================")
print(" Running Unified Advanced Epi Profiling Pipeline (Script 14)")
print("==========================================================")

meta_df = pd.read_csv(METADATA_CSV)
cgmlst_df = pd.read_csv(CGMLST_TSV, sep="\t", index_col=0)

cgmlst_df.index = cgmlst_df.index.str.strip()
meta_df['Assembly'] = meta_df['Assembly'].astype(str).str.strip()
meta_filtered = meta_df[meta_df['Assembly'].isin(cgmlst_df.index)].copy()

# Dynamic Year Detection
year_col = None
for col in meta_filtered.columns:
    if 'year' in col.lower() or 'date' in col.lower():
        year_col = col
        break

print(f"[INFO] Detected year column: {year_col}")

if year_col:
    def extract_year(val):
        try:
            s = str(val).strip()
            if len(s) >= 4 and s[:4].isdigit():
                return int(s[:4])
            return np.nan
        except:
            return np.nan
    meta_filtered['Numeric_Year'] = meta_filtered[year_col].apply(extract_year)
else:
    meta_filtered['Numeric_Year'] = np.nan

def custom_hamming_distances(matrix_df):
    arr = matrix_df.values
    def custom_hamming(u, v):
        valid = (u >= 0) & (v >= 0)
        if not np.any(valid): return 1.0
        return np.sum(u[valid] != v[valid]) / np.sum(valid)
    distances = pdist(arr, metric=custom_hamming)
    return squareform(distances * arr.shape[1])

# Compute Global Distance Matrix and HierCC Single-Linkage Clustering
print("[INFO] Computing global cgMLST distance matrix for HierCC...")
global_dist_matrix = custom_hamming_distances(cgmlst_df)
linkage_matrix = sch.linkage(squareform(global_dist_matrix), method='single')

hiercc_levels = {'HC10': 10, 'HC50': 50}
clusters_dict = {'Assembly': cgmlst_df.index}
for level_name, threshold in hiercc_levels.items():
    labels = sch.fcluster(linkage_matrix, t=threshold, criterion='distance')
    clusters_dict[level_name] = [f"{level_name}_{c}" for c in labels]

clusters_df = pd.DataFrame(clusters_dict)
meta_filtered = pd.merge(meta_filtered, clusters_df, on='Assembly', how='left')

# Save Enriched Metadata
enriched_meta_path = f"{OUTPUT_DIR}/enriched_metadata_with_hiercc.csv"
meta_filtered.to_csv(enriched_meta_path, index=False)
print(f"[SUCCESS] Saved enriched metadata to: {enriched_meta_path}")

prod_systems = meta_filtered['Production_system'].dropna().astype(str).unique().tolist()

def sort_chronological(indices):
    special_first = [str(x) for x in indices if str(x) == '<=2000']
    middle = sorted([str(x) for x in indices if str(x) != '<=2000' and str(x) != '>=2011'])
    special_last = [str(x) for x in indices if str(x) == '>=2011']
    return special_first + middle + special_last

for prod_sys in prod_systems:
    prod_clean = str(prod_sys).lower().replace(" ", "_").replace("-", "_")
    base_prod_dir = f"{OUTPUT_DIR}/{prod_clean}"
    
    prod_meta_all = meta_filtered[meta_filtered['Production_system'].astype(str) == str(prod_sys)].copy()
    if prod_meta_all.empty: continue

    conditions = {
        "with_post_harvest": prod_meta_all,
        "without_post_harvest": prod_meta_all[prod_meta_all['Supply_chain_stage'].astype(str).str.lower() != 'post-harvest']
    }

    for cond_name, sub_meta in conditions.items():
        if sub_meta.empty: continue
        
        cond_dir = f"{base_prod_dir}/{cond_name}"
        os.makedirs(cond_dir, exist_ok=True)
        print(f"[INFO] Processing [{prod_sys}] -> Condition: {cond_name} -> {cond_dir}")

        # 1. Supply Chain Count Bar Plot
        if 'Supply_chain_stage' in sub_meta.columns:
            sc_counts = sub_meta['Supply_chain_stage'].value_counts()
            if not sc_counts.empty:
                sc_counts.index = [f"{idx} (n={sc_counts[idx]})" for idx in sc_counts.index]
                plt.figure(figsize=(10, 6))
                sc_counts.plot(kind='bar', color='steelblue', width=0.5)
                plt.title(f"[{prod_sys} - {cond_name}] Isolate Count per Supply Chain Stage", fontsize=12)
                plt.ylabel("Isolate Count", fontsize=12)
                plt.xlabel("Supply Chain Stage", fontsize=12)
                plt.xticks(rotation=25, ha='right', fontsize=10)
                plt.yticks(fontsize=10)
                plt.tight_layout()
                plt.savefig(f"{cond_dir}/count_per_supply_chain.png", bbox_inches='tight')
                plt.close()

        # 2. Stacked Proportion Plots for HC10 and HC50 (Supply Chain Stage) - Top 14 and Top 19
        for thresh in ['HC10', 'HC50']:
            st_counts = sub_meta.groupby(['Supply_chain_stage', thresh]).size().unstack(fill_value=0)
            if st_counts.empty: continue
            
            row_sums = st_counts.sum(axis=1).replace(0, 1)
            st_proportions = st_counts.div(row_sums, axis=0).fillna(0)
            
            stage_counts = sub_meta['Supply_chain_stage'].value_counts()
            st_proportions.index = [f"{idx} (n={stage_counts.get(idx, 0)})" for idx in st_proportions.index]

            for top_n in [14, 19]:
                if st_proportions.shape[1] > top_n:
                    top_cols = st_proportions.sum(axis=0).nlargest(top_n).index
                    other_cols = [c for c in st_proportions.columns if c not in top_cols]
                    df_plot = st_proportions[top_cols].copy()
                    df_plot['Others'] = st_proportions[other_cols].sum(axis=1)
                else:
                    df_plot = st_proportions.copy()

                plt.figure(figsize=(14, 8))
                df_plot.plot(kind='bar', stacked=True, cmap='tab20', width=0.6)
                plt.title(f"[{prod_sys} - {cond_name}] Lineage Dynamics ({thresh}) - Top {top_n}", fontsize=12)
                plt.ylabel("Proportional Scale", fontsize=12)
                plt.xlabel("Supply Chain Stage", fontsize=12)
                plt.legend(title=f"HierCC ({thresh})", bbox_to_anchor=(1.02, 1), loc='upper left', ncol=2, fontsize=10, title_fontsize=11)
                plt.xticks(rotation=25, ha='right', fontsize=10)
                plt.yticks(fontsize=10)
                plt.tight_layout()
                plt.savefig(f"{cond_dir}/supply_chain_proportions_{thresh.lower()}_top{top_n}.png", bbox_inches='tight')
                plt.close()

        # Helper function for year bins
        def get_year_bins(y, mode):
            if y <= 2000:
                return "<=2000"
            if mode == 'individual':
                return str(int(y))
            elif mode == '5year':
                start = 2001 + ((int(y) - 2001) // 5) * 5
                end = start + 4
                return f"{start}-{end}"
            elif mode == '10year':
                if y >= 2011:
                    return ">=2011"
                start = 2001 + ((int(y) - 2001) // 10) * 10
                end = start + 9
                return f"{start}-{end}"
            return str(int(y))

        # 3. Temporal Count Bar Plots (<=2000 first, >=2011 last)
        if year_col and not sub_meta['Numeric_Year'].dropna().empty:
            sub_meta_temp = sub_meta.copy()
            
            # Individual Years
            sub_meta_temp['Year_Individual'] = sub_meta_temp['Numeric_Year'].apply(lambda x: get_year_bins(x, 'individual') if pd.notnull(x) else 'Unknown')
            ind_counts = sub_meta_temp[sub_meta_temp['Year_Individual'] != 'Unknown']['Year_Individual'].value_counts()
            ind_counts = ind_counts.reindex(sort_chronological(ind_counts.index)).fillna(0)
            ind_counts.index = [f"{idx} (n={int(ind_counts[idx])})" for idx in ind_counts.index]

            plt.figure(figsize=(12, 6))
            ind_counts.plot(kind='bar', color='forestgreen', width=0.6)
            plt.title(f"[{prod_sys} - {cond_name}] Count per Year", fontsize=12)
            plt.ylabel("Isolate Count", fontsize=12)
            plt.xlabel("Year", fontsize=12)
            plt.xticks(rotation=45, ha='right', fontsize=10)
            plt.yticks(fontsize=10)
            plt.tight_layout()
            plt.savefig(f"{cond_dir}/count_per_year_individual.png", bbox_inches='tight')
            plt.close()

            # 5-Year Intervals
            sub_meta_temp['Year_5Yr'] = sub_meta_temp['Numeric_Year'].apply(lambda x: get_year_bins(x, '5year') if pd.notnull(x) else 'Unknown')
            y5_counts = sub_meta_temp[sub_meta_temp['Year_5Yr'] != 'Unknown']['Year_5Yr'].value_counts()
            y5_counts = y5_counts.reindex(sort_chronological(y5_counts.index)).fillna(0)
            y5_counts.index = [f"{idx} (n={int(y5_counts[idx])})" for idx in y5_counts.index]

            plt.figure(figsize=(10, 6))
            y5_counts.plot(kind='bar', color='darkorange', width=0.6)
            plt.title(f"[{prod_sys} - {cond_name}] Count per 5-Year Interval", fontsize=12)
            plt.ylabel("Isolate Count", fontsize=12)
            plt.xlabel("5-Year Interval", fontsize=12)
            plt.xticks(rotation=25, ha='right', fontsize=10)
            plt.yticks(fontsize=10)
            plt.tight_layout()
            plt.savefig(f"{cond_dir}/count_per_year_5interval.png", bbox_inches='tight')
            plt.close()

            # 10-Year Intervals (Merged >=2011)
            sub_meta_temp['Year_10Yr'] = sub_meta_temp['Numeric_Year'].apply(lambda x: get_year_bins(x, '10year') if pd.notnull(x) else 'Unknown')
            y10_counts = sub_meta_temp[sub_meta_temp['Year_10Yr'] != 'Unknown']['Year_10Yr'].value_counts()
            y10_counts = y10_counts.reindex(sort_chronological(y10_counts.index)).fillna(0)
            y10_counts.index = [f"{idx} (n={int(y10_counts[idx])})" for idx in y10_counts.index]

            plt.figure(figsize=(10, 6))
            y10_counts.plot(kind='bar', color='purple', width=0.6)
            plt.title(f"[{prod_sys} - {cond_name}] Count per 10-Year Interval (>=2011 merged)", fontsize=12)
            plt.ylabel("Isolate Count", fontsize=12)
            plt.xlabel("10-Year Interval", fontsize=12)
            plt.xticks(rotation=25, ha='right', fontsize=10)
            plt.yticks(fontsize=10)
            plt.tight_layout()
            plt.savefig(f"{cond_dir}/count_per_year_10interval.png", bbox_inches='tight')
            plt.close()

            # 3b. Stacked Proportion Plots for HC10 and HC50 across 10-Year Intervals (Top 14 and Top 19, >=2011 merged)
            sub_meta_temporal_hc = sub_meta.copy()
            sub_meta_temporal_hc['Year_10Yr'] = sub_meta_temporal_hc['Numeric_Year'].apply(lambda x: get_year_bins(x, '10year') if pd.notnull(x) else 'Unknown')
            sub_meta_temporal_hc = sub_meta_temporal_hc[sub_meta_temporal_hc['Year_10Yr'] != 'Unknown']
            
            for thresh in ['HC10', 'HC50']:
                yr_counts = sub_meta_temporal_hc.groupby(['Year_10Yr', thresh]).size().unstack(fill_value=0)
                if yr_counts.empty: continue
                
                row_sums = yr_counts.sum(axis=1).replace(0, 1)
                yr_proportions = yr_counts.div(row_sums, axis=0).fillna(0)
                yr_proportions = yr_proportions.reindex(sort_chronological(yr_proportions.index)).fillna(0)
                
                bin_counts = sub_meta_temporal_hc['Year_10Yr'].value_counts()
                yr_proportions.index = [f"{idx} (n={bin_counts.get(idx, 0)})" for idx in yr_proportions.index]

                for top_n in [14, 19]:
                    if yr_proportions.shape[1] > top_n:
                        top_cols = yr_proportions.sum(axis=0).nlargest(top_n).index
                        other_cols = [c for c in yr_proportions.columns if c not in top_cols]
                        df_plot_yr = yr_proportions[top_cols].copy()
                        df_plot_yr['Others'] = yr_proportions[other_cols].sum(axis=1)
                    else:
                        df_plot_yr = yr_proportions.copy()

                    plt.figure(figsize=(14, 8))
                    df_plot_yr.plot(kind='bar', stacked=True, cmap='tab20', width=0.6)
                    plt.title(f"[{prod_sys} - {cond_name}] Lineage Dynamics across 10-Year Intervals ({thresh}) - Top {top_n}", fontsize=12)
                    plt.ylabel("Proportional Scale", fontsize=12)
                    plt.xlabel("10-Year Interval", fontsize=12)
                    plt.legend(title=f"HierCC ({thresh})", bbox_to_anchor=(1.02, 1), loc='upper left', ncol=2, fontsize=10, title_fontsize=11)
                    plt.xticks(rotation=25, ha='right', fontsize=10)
                    plt.yticks(fontsize=10)
                    plt.tight_layout()
                    plt.savefig(f"{cond_dir}/year_10yr_proportions_{thresh.lower()}_top{top_n}.png", bbox_inches='tight')
                    plt.close()

        # 4. Heatmaps and GrapeTree MSTs for the subset
        g_ids = sub_meta['Assembly'].tolist()
        sub_matrix = cgmlst_df.loc[cgmlst_df.index.isin(g_ids)]
        
        if len(sub_matrix) >= 2:
            sub_dist_matrix = custom_hamming_distances(sub_matrix)
            
            if len(sub_matrix) <= 200:
                plt.figure(figsize=(12, 10))
                g = sns.clustermap(sub_dist_matrix, cmap="rocket_r", xticklabels=False, yticklabels=False, cbar_kws={'label': 'Allele Differences'})
                g.fig.suptitle(f"Transmission Cluster Map: {prod_sys} ({cond_name})", y=1.02, fontsize=12)
                g.savefig(f"{cond_dir}/{prod_clean}_{cond_name}_epi_heatmap.png", bbox_inches='tight')
                plt.close()

            G = nx.Graph()
            nodes = sub_matrix.index.tolist()
            for i in range(len(nodes)):
                for j in range(i+1, len(nodes)):
                    G.add_edge(nodes[i], nodes[j], weight=sub_dist_matrix[i, j])
            mst = nx.minimum_spanning_tree(G, weight='weight')
            
            def to_newick(g, node, parent=None):
                subtrees = []
                for neighbor in g.neighbors(node):
                    if neighbor != parent:
                        edge_w = g[node][neighbor]['weight']
                        subtrees.append(to_newick(g, neighbor, node) + f":{edge_w}")
                if subtrees:
                    return f"({','.join(subtrees)}){node}"
                return node

            with open(f"{cond_dir}/{prod_clean}_{cond_name}_grapetree_mst.nwk", "w") as f:
                f.write(to_newick(mst, nodes[0]) + ";")

print("\n==========================================================")
print(" Success! Script 14 updated with reduced title font size.")
print(f" Output root directory: {OUTPUT_DIR}")
print("==========================================================")
