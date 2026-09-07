#!/usr/bin/env python3

import os
import re
import math
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
sns.set_theme(style="white", context="talk")
plt.rcParams["savefig.dpi"] = 300

# ==============================================================================
# Script 14C: Binned cgMLST Temporal Dynamics within Core-SNP Clades
#
# Purpose:
#   Use broad, better-supported time periods instead of sparse exact years.
#
# Time bins:
#   <2010
#   2010-2015
#   >2016
#
# Main biological question:
#   Which cgMLST clusters within each core-SNP clade show evidence of:
#     - persistent-like repeated detection across broad time periods
#     - emerging-like late-period detection
#     - re-emergence-like / unsampled-persistence pattern
#     - temporally localized expansion
#     - sporadic/minority detection
# ==============================================================================

BASE_DIR = "/home/samuelajulo/epi/Ecoli_epi"

MAIN_CGMLST_LEVEL = os.environ.get("MAIN_CGMLST_LEVEL", "HC20").upper()
TOP_N_CGMLST_PLOT = os.environ.get("TOP_N_CGMLST_PLOT", "ALL").upper()

INPUT_ROOT = f"{BASE_DIR}/results/14B_cgmlst_clade_temporal_dynamics_704_{MAIN_CGMLST_LEVEL}"

# fallback for older HC10 folder
if not os.path.exists(INPUT_ROOT) and MAIN_CGMLST_LEVEL == "HC10":
    INPUT_ROOT = f"{BASE_DIR}/results/14B_cgmlst_clade_temporal_dynamics_704"

GENOME_ASSIGN = f"{INPUT_ROOT}/tables/14B_genome_level_temporal_category_assignments.csv"
CLUSTER_SUMMARY = f"{INPUT_ROOT}/tables/14B_cgmlst_cluster_temporal_category_summary.csv"

OUTPUT_DIR = f"{BASE_DIR}/results/14C_cgmlst_binned_temporal_dynamics_704_{MAIN_CGMLST_LEVEL}"
TABLE_DIR = f"{OUTPUT_DIR}/tables"
FIGURE_DIR = f"{OUTPUT_DIR}/figures"
DEBUG_DIR = f"{OUTPUT_DIR}/debug"

for d in [OUTPUT_DIR, TABLE_DIR, FIGURE_DIR, DEBUG_DIR]:
    os.makedirs(d, exist_ok=True)

print("==========================================================")
print(" Script 14C: Binned cgMLST Temporal Dynamics")
print("==========================================================")
print(f"[INFO] cgMLST level: {MAIN_CGMLST_LEVEL}")
print(f"[INFO] Input root: {INPUT_ROOT}")
print(f"[INFO] Output root: {OUTPUT_DIR}")

if not os.path.exists(GENOME_ASSIGN):
    raise SystemExit(f"[ERROR] Missing 14B genome assignment file: {GENOME_ASSIGN}")

if not os.path.exists(CLUSTER_SUMMARY):
    raise SystemExit(f"[ERROR] Missing 14B cluster summary file: {CLUSTER_SUMMARY}")

# ==============================================================================
# Helpers
# ==============================================================================

def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def accession_root(x):
    x = str(x)
    m = re.search(r"(GC[AF]_[0-9]+)(?:\.[0-9]+)?", x)
    return m.group(1) if m else x.strip()

def extract_year(x):
    if pd.isna(x):
        return np.nan
    m = re.search(r"(18|19|20)[0-9]{2}", str(x))
    return int(m.group(0)) if m else np.nan

def assign_time_bin(y):
    if pd.isna(y):
        return "Unknown"
    y = int(y)
    if y < 2010:
        return "<2010"
    if 2010 <= y <= 2015:
        return "2010-2015"
    if y >= 2016:
        return ">2016"
    return "Unknown"

BIN_LEVELS = ["<2010", "2010-2015", ">2016"]

def years_string(years):
    years = sorted(set([int(y) for y in years if not pd.isna(y)]))
    return ";".join(map(str, years))

def bin_string(bins):
    bins = [b for b in BIN_LEVELS if b in set(bins)]
    return ";".join(bins)

def amr_summary_string(vals):
    vals = [int(v) for v in vals if not pd.isna(v)]
    if len(vals) == 0:
        return "median=0; mean=0; range=0-0"
    return f"median={np.median(vals):.1f}; mean={np.mean(vals):.2f}; range={np.min(vals)}-{np.max(vals)}"

def classify_binned(row):
    n = int(row["N_genomes"])
    n_bins = int(row["N_time_bins_detected"])
    has_early = row["Has_<2010"]
    has_mid = row["Has_2010_2015"]
    has_late = row["Has_>2016"]

    if row["N_genomes_with_year"] == 0:
        return "Unknown year"

    if n <= 2:
        return "Sporadic/minority"

    if has_early and has_mid and has_late:
        return "Persistent-like across broad time periods"

    if has_early and (not has_mid) and has_late:
        return "Re-emergence-like / unsampled persistence"

    if (not has_early) and (not has_mid) and has_late:
        if n >= 5:
            return "Emerging-like late-period cluster"
        return "Late-period sporadic/minority"

    if (not has_early) and has_mid and has_late:
        return "Late-period repeated detection"

    if has_early and has_mid and (not has_late):
        return "Early-to-mid repeated detection; not observed late"

    if n_bins == 1 and n >= 10:
        return "Temporally localized expansion"

    if n_bins == 1:
        return "Sporadic/minority"

    if n_bins >= 2:
        return "Repeated detection, limited temporal evidence"

    return "Limited temporal evidence"

def evidence_level(row):
    n = int(row["N_genomes"])
    n_bins = int(row["N_time_bins_detected"])
    snp_available = not pd.isna(row.get("Median_within_cluster_SNP", np.nan))

    if n >= 10 and n_bins >= 2 and snp_available:
        return "Strong"
    if n >= 5 and n_bins >= 2 and snp_available:
        return "Moderate"
    if n >= 3:
        return "Limited"
    return "Exploratory"

def make_interpretation(row):
    clade = row["Clade"]
    cluster = row["cgMLST_cluster"]
    n = int(row["N_genomes"])
    cat = row["Temporal_category_binned"]
    bins = row["Time_bins_detected"]
    years = row["Exact_years_detected"]

    med_snp = row.get("Median_within_cluster_SNP", np.nan)
    snp_text = (
        f"median within-cluster SNP distance = {med_snp:.1f}"
        if not pd.isna(med_snp)
        else "within-cluster SNP distance not available"
    )

    if cat == "Persistent-like across broad time periods":
        return (
            f"{cluster} within {clade} was observed across all broad time periods "
            f"({bins}; exact years: {years}) among sampled genomes, with {snp_text}. "
            f"This supports a persistent-like or repeatedly circulating sublineage pattern."
        )

    if cat == "Re-emergence-like / unsampled persistence":
        return (
            f"{cluster} within {clade} was observed in an early period and again in the late period "
            f"({bins}; exact years: {years}) but not in the middle bin. This is consistent with "
            f"apparent re-emergence or unsampled persistence; {snp_text}."
        )

    if cat == "Emerging-like late-period cluster":
        return (
            f"{cluster} within {clade} was observed only in the late period ({bins}; exact years: {years}) "
            f"with {n} genomes, supporting an apparent late-period emergence pattern among sampled genomes; "
            f"{snp_text}."
        )

    if cat == "Temporally localized expansion":
        return (
            f"{cluster} within {clade} was concentrated in one broad time period ({bins}; exact years: {years}) "
            f"with {n} genomes, suggesting a temporally localized expansion among sampled genomes; {snp_text}."
        )

    if cat == "Sporadic/minority":
        return (
            f"{cluster} within {clade} had limited temporal representation ({bins}; exact years: {years}) "
            f"and is best interpreted as sporadic/minority detection among sampled genomes."
        )

    return (
        f"{cluster} within {clade} was observed in {bins} (exact years: {years}) with {n} genomes. "
        f"The pattern is best interpreted as {cat}; {snp_text}."
    )

# ==============================================================================
# Read inputs
# ==============================================================================

df = pd.read_csv(GENOME_ASSIGN)
cluster_14b = pd.read_csv(CLUSTER_SUMMARY)

id_col = pick_col(df, ["Accession_root", "Assembly", "Tip_ID", "Isolate", "Sample", "sample"])
if id_col is None:
    raise SystemExit("[ERROR] Could not detect genome ID column in 14B assignment file.")

if "Accession_root" not in df.columns:
    df["Accession_root"] = df[id_col].apply(accession_root)

year_col = pick_col(df, ["Numeric_Year", "Collection_year_num", "Collection_year", "Collection year", "Year", "year"])

if year_col is None:
    for c in df.columns:
        if "year" in c.lower() or "date" in c.lower():
            year_col = c
            break

if year_col is None:
    raise SystemExit("[ERROR] Could not detect year column in 14B assignment file.")

required_cols = ["Clade", "cgMLST_cluster"]
for c in required_cols:
    if c not in df.columns:
        raise SystemExit(f"[ERROR] Missing required column in 14B assignment file: {c}")

df["Year_clean"] = df[year_col].apply(extract_year)
df["Time_bin"] = df["Year_clean"].apply(assign_time_bin)
df["Time_bin"] = pd.Categorical(df["Time_bin"], categories=BIN_LEVELS + ["Unknown"], ordered=True)

df.to_csv(f"{TABLE_DIR}/14C_genome_level_assignments_with_time_bins.csv", index=False)

# ==============================================================================
# Overall bin counts
# ==============================================================================

overall_bin_counts = (
    df[df["Time_bin"] != "Unknown"]
    .groupby("Time_bin", observed=False)
    .agg(N_genomes=("Accession_root", "nunique"))
    .reset_index()
)

overall_bin_counts = overall_bin_counts[overall_bin_counts["Time_bin"].isin(BIN_LEVELS)].copy()
overall_bin_counts["Percent_of_704"] = round(100 * overall_bin_counts["N_genomes"] / df["Accession_root"].nunique(), 2)

overall_bin_counts.to_csv(f"{TABLE_DIR}/14C_overall_time_bin_counts_704.csv", index=False)

# ==============================================================================
# Cluster-bin count table
# ==============================================================================

cluster_bin_counts = (
    df[df["Time_bin"] != "Unknown"]
    .groupby(["Clade", "cgMLST_cluster", "Time_bin"], observed=False)
    .agg(N_genomes=("Accession_root", "nunique"))
    .reset_index()
)

cluster_bin_counts = cluster_bin_counts[cluster_bin_counts["Time_bin"].isin(BIN_LEVELS)].copy()
cluster_bin_counts.to_csv(f"{TABLE_DIR}/14C_cgmlst_cluster_time_bin_counts.csv", index=False)

exact_year_counts = (
    df.dropna(subset=["Year_clean"])
    .assign(Year_clean=lambda d: d["Year_clean"].astype(int))
    .groupby(["Clade", "cgMLST_cluster", "Year_clean"])
    .agg(N_genomes=("Accession_root", "nunique"))
    .reset_index()
    .sort_values(["Clade", "cgMLST_cluster", "Year_clean"])
)

exact_year_counts.to_csv(f"{TABLE_DIR}/14C_cgmlst_cluster_exact_year_counts.csv", index=False)

# ==============================================================================
# Binned cluster summary
# ==============================================================================

def summarize_cluster(g):
    years = sorted(set([int(y) for y in g["Year_clean"].dropna()]))
    bins = [str(x) for x in g["Time_bin"].dropna().astype(str) if x != "Unknown"]
    bin_set = set(bins)

    n = g["Accession_root"].nunique()
    n_with_year = g.dropna(subset=["Year_clean"])["Accession_root"].nunique()

    bin_counts = {
        b: g[g["Time_bin"].astype(str) == b]["Accession_root"].nunique()
        for b in BIN_LEVELS
    }

    amr_col = "AMRFinder_feature_count" if "AMRFinder_feature_count" in g.columns else None
    amr_vals = g[amr_col].fillna(0).astype(int).tolist() if amr_col else [0] * len(g)

    return pd.Series({
        "N_genomes": n,
        "N_genomes_with_year": n_with_year,
        "First_year": min(years) if years else np.nan,
        "Last_year": max(years) if years else np.nan,
        "Year_span": (max(years) - min(years)) if years else np.nan,
        "Exact_years_detected": years_string(years),
        "N_exact_years_detected": len(years),
        "Time_bins_detected": bin_string(bins),
        "N_time_bins_detected": len([b for b in BIN_LEVELS if b in bin_set]),
        "Has_<2010": "<2010" in bin_set,
        "Has_2010_2015": "2010-2015" in bin_set,
        "Has_>2016": ">2016" in bin_set,
        "Count_<2010": bin_counts["<2010"],
        "Count_2010_2015": bin_counts["2010-2015"],
        "Count_>2016": bin_counts[">2016"],
        "AMRFinder_feature_count_summary": amr_summary_string(amr_vals),
        "Median_AMRFinder_feature_count": float(np.median(amr_vals)) if amr_vals else 0,
        "Mean_AMRFinder_feature_count": float(np.mean(amr_vals)) if amr_vals else 0,
        "Max_AMRFinder_feature_count": int(np.max(amr_vals)) if amr_vals else 0
    })

cluster_summary = (
    df.groupby(["Clade", "cgMLST_cluster"], dropna=False)
    .apply(summarize_cluster)
    .reset_index()
)

snp_cols = [
    "Clade",
    "cgMLST_cluster",
    "Median_within_cluster_SNP",
    "Mean_within_cluster_SNP",
    "Max_within_cluster_SNP",
    "IQR_within_cluster_SNP"
]

available_snp_cols = [c for c in snp_cols if c in cluster_14b.columns]

if "Median_within_cluster_SNP" in cluster_14b.columns:
    cluster_summary = cluster_summary.merge(
        cluster_14b[available_snp_cols].drop_duplicates(["Clade", "cgMLST_cluster"]),
        on=["Clade", "cgMLST_cluster"],
        how="left"
    )
else:
    cluster_summary["Median_within_cluster_SNP"] = np.nan
    cluster_summary["Mean_within_cluster_SNP"] = np.nan
    cluster_summary["Max_within_cluster_SNP"] = np.nan
    cluster_summary["IQR_within_cluster_SNP"] = np.nan

cluster_summary["Temporal_category_binned"] = cluster_summary.apply(classify_binned, axis=1)
cluster_summary["Evidence_level_binned"] = cluster_summary.apply(evidence_level, axis=1)
cluster_summary["Interpretation_binned"] = cluster_summary.apply(make_interpretation, axis=1)

cluster_summary = cluster_summary.sort_values(
    ["Clade", "N_genomes", "cgMLST_cluster"],
    ascending=[True, False, True]
)

main_cols = [
    "Clade",
    "cgMLST_cluster",
    "N_genomes",
    "N_genomes_with_year",
    "First_year",
    "Last_year",
    "Year_span",
    "Exact_years_detected",
    "N_exact_years_detected",
    "Time_bins_detected",
    "N_time_bins_detected",
    "Count_<2010",
    "Count_2010_2015",
    "Count_>2016",
    "Median_within_cluster_SNP",
    "Mean_within_cluster_SNP",
    "Max_within_cluster_SNP",
    "IQR_within_cluster_SNP",
    "AMRFinder_feature_count_summary",
    "Median_AMRFinder_feature_count",
    "Mean_AMRFinder_feature_count",
    "Max_AMRFinder_feature_count",
    "Temporal_category_binned",
    "Evidence_level_binned",
    "Interpretation_binned"
]

cluster_summary[main_cols].to_csv(
    f"{TABLE_DIR}/14C_cgmlst_cluster_binned_temporal_category_summary.csv",
    index=False
)

# ==============================================================================
# Genome-level binned assignments
# ==============================================================================

genome_binned = df.merge(
    cluster_summary[[
        "Clade",
        "cgMLST_cluster",
        "N_genomes",
        "Time_bins_detected",
        "N_time_bins_detected",
        "Temporal_category_binned",
        "Evidence_level_binned",
        "Interpretation_binned"
    ]],
    on=["Clade", "cgMLST_cluster"],
    how="left",
    suffixes=("", "_binned_cluster")
)

genome_binned.to_csv(
    f"{TABLE_DIR}/14C_genome_level_binned_temporal_category_assignments.csv",
    index=False
)

# ==============================================================================
# Category summaries
# ==============================================================================

category_summary = (
    cluster_summary
    .groupby(["Clade", "Temporal_category_binned", "Evidence_level_binned"])
    .agg(
        N_cgMLST_clusters=("cgMLST_cluster", "nunique"),
        N_genomes=("N_genomes", "sum")
    )
    .reset_index()
    .sort_values(["Clade", "Temporal_category_binned", "Evidence_level_binned"])
)

category_summary.to_csv(
    f"{TABLE_DIR}/14C_binned_temporal_category_summary_by_clade.csv",
    index=False
)

# ==============================================================================
# Figures
# ==============================================================================

print("[INFO] Generating binned figures...")

# Overall bin count
plt.figure(figsize=(8, 5))
ax = sns.barplot(
    data=overall_bin_counts,
    x="Time_bin",
    y="N_genomes",
    color="grey"
)
for container in ax.containers:
    ax.bar_label(container, fmt="%d", fontsize=10)
plt.title("Genome counts by broad time period")
plt.xlabel("Time period")
plt.ylabel("Genome count")
plt.tight_layout()
plt.savefig(f"{FIGURE_DIR}/14C_overall_genome_counts_by_time_bin.png", bbox_inches="tight")
plt.savefig(f"{FIGURE_DIR}/14C_overall_genome_counts_by_time_bin.pdf", bbox_inches="tight")
plt.close()

# Select clusters for plot
plot_summary = cluster_summary.sort_values("N_genomes", ascending=False).copy()

if TOP_N_CGMLST_PLOT == "ALL":
    plot_clusters = plot_summary[["Clade", "cgMLST_cluster", "N_genomes", "Temporal_category_binned"]].copy()
else:
    top_n = int(TOP_N_CGMLST_PLOT)
    plot_clusters = plot_summary.head(top_n)[["Clade", "cgMLST_cluster", "N_genomes", "Temporal_category_binned"]].copy()

plot_clusters.to_csv(
    f"{TABLE_DIR}/14C_cgmlst_clusters_in_binned_heatmap_plot.csv",
    index=False
)

plot_df = cluster_bin_counts.merge(
    plot_clusters,
    on=["Clade", "cgMLST_cluster"],
    how="inner"
)

plot_df["Cluster_label"] = (
    plot_df["Clade"].astype(str)
    + " | "
    + plot_df["cgMLST_cluster"].astype(str)
    + " (n="
    + plot_df["N_genomes_y"].astype(int).astype(str)
    + "; "
    + plot_df["Temporal_category_binned"].astype(str)
    + ")"
)

row_order = (
    plot_df.groupby("Cluster_label")["N_genomes_x"]
    .sum()
    .sort_values(ascending=False)
    .index
    .tolist()
)

plot_df["Time_bin"] = pd.Categorical(plot_df["Time_bin"], categories=BIN_LEVELS, ordered=True)

heat_mat = (
    plot_df
    .pivot_table(
        index="Cluster_label",
        columns="Time_bin",
        values="N_genomes_x",
        fill_value=0,
        aggfunc="sum",
        observed=False
    )
)

heat_mat = heat_mat.reindex(row_order)
heat_mat = heat_mat[BIN_LEVELS]

heat_labels = heat_mat.astype(int).astype(str)
heat_labels = heat_labels.mask(heat_mat == 0, "")

fig_height = max(8, 0.45 * len(heat_mat))

plt.figure(figsize=(13, fig_height))
ax = sns.heatmap(
    heat_mat,
    cmap="Greys",
    linewidths=0.35,
    linecolor="white",
    cbar_kws={"label": "Genome count"},
    annot=heat_labels,
    fmt="",
    annot_kws={"fontsize": 8},
    yticklabels=True,
    xticklabels=True
)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=10)

if TOP_N_CGMLST_PLOT == "ALL":
    plot_title = f"All {MAIN_CGMLST_LEVEL} cgMLST clusters across broad time periods"
else:
    plot_title = f"Top {TOP_N_CGMLST_PLOT} {MAIN_CGMLST_LEVEL} cgMLST clusters across broad time periods"

plt.title(plot_title)
plt.xlabel("Broad time period")
plt.ylabel("Clade | cgMLST cluster (n; binned temporal category)")
plt.tight_layout()
plt.savefig(f"{FIGURE_DIR}/14C_cgmlst_cluster_binned_heatmap.png", bbox_inches="tight")
plt.savefig(f"{FIGURE_DIR}/14C_cgmlst_cluster_binned_heatmap.pdf", bbox_inches="tight")
plt.close()

# Bubble plot observed only
bubble_df = plot_df[plot_df["N_genomes_x"] > 0].copy()
bubble_df["Cluster_label"] = pd.Categorical(
    bubble_df["Cluster_label"],
    categories=row_order[::-1],
    ordered=True
)

plt.figure(figsize=(13, fig_height))
ax = sns.scatterplot(
    data=bubble_df,
    x="Time_bin",
    y="Cluster_label",
    size="N_genomes_x",
    hue="Temporal_category_binned",
    sizes=(80, 700),
    alpha=0.85,
    edgecolor="black",
    linewidth=0.3
)
plt.title(plot_title.replace("across broad time periods", "observed across broad time periods"))
plt.xlabel("Broad time period")
plt.ylabel("Clade | cgMLST cluster (n; binned temporal category)")
plt.legend(
    title="Temporal category / count",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0
)
plt.tight_layout()
plt.savefig(f"{FIGURE_DIR}/14C_cgmlst_cluster_binned_bubbleplot_observed_only.png", bbox_inches="tight")
plt.savefig(f"{FIGURE_DIR}/14C_cgmlst_cluster_binned_bubbleplot_observed_only.pdf", bbox_inches="tight")
plt.close()

plot_df.to_csv(
    f"{TABLE_DIR}/14C_binned_heatmap_plotting_data.csv",
    index=False
)

# ==============================================================================
# README and manifest
# ==============================================================================

readme = [
    "Script 14C: Binned cgMLST temporal dynamics",
    "",
    "Time bins used:",
    "<2010",
    "2010-2015",
    ">2016",
    "",
    "Why bins were used:",
    "Exact collection years were highly imbalanced. Broad bins improve interpretability and reduce overclaiming from sparsely sampled years.",
    "",
    "Main output:",
    "tables/14C_cgmlst_cluster_binned_temporal_category_summary.csv",
    "",
    "Genome-level output:",
    "tables/14C_genome_level_binned_temporal_category_assignments.csv",
    "",
    "Interpretation caution:",
    "These categories describe observed temporal patterns among public genomes. Absence in a bin means not observed, not confirmed biological absence."
]

with open(f"{OUTPUT_DIR}/README_14C_binned_temporal_dynamics.txt", "w") as f:
    f.write("\n".join(readme) + "\n")

manifest = pd.DataFrame({
    "output_type": [
        "overall_time_bin_counts",
        "cluster_time_bin_counts",
        "cluster_binned_summary",
        "genome_level_binned_assignment",
        "category_summary_by_clade",
        "binned_heatmap",
        "binned_bubbleplot",
        "readme"
    ],
    "path": [
        f"{TABLE_DIR}/14C_overall_time_bin_counts_704.csv",
        f"{TABLE_DIR}/14C_cgmlst_cluster_time_bin_counts.csv",
        f"{TABLE_DIR}/14C_cgmlst_cluster_binned_temporal_category_summary.csv",
        f"{TABLE_DIR}/14C_genome_level_binned_temporal_category_assignments.csv",
        f"{TABLE_DIR}/14C_binned_temporal_category_summary_by_clade.csv",
        f"{FIGURE_DIR}/14C_cgmlst_cluster_binned_heatmap.pdf",
        f"{FIGURE_DIR}/14C_cgmlst_cluster_binned_bubbleplot_observed_only.pdf",
        f"{OUTPUT_DIR}/README_14C_binned_temporal_dynamics.txt"
    ]
})
manifest.to_csv(f"{OUTPUT_DIR}/14C_output_manifest.csv", index=False)

# ==============================================================================
# Terminal report
# ==============================================================================

print("")
print("==========================================================")
print(" Script 14C completed successfully")
print("==========================================================")
print(f"[OUTPUT] {OUTPUT_DIR}")

print("\n=== Overall bin counts ===")
print(overall_bin_counts.to_string(index=False))

print("\n=== Binned temporal category summary by clade ===")
print(category_summary.to_string(index=False))

print("\n=== Top 25 cgMLST clusters by N_genomes ===")
print(
    cluster_summary.sort_values("N_genomes", ascending=False)[[
        "Clade",
        "cgMLST_cluster",
        "N_genomes",
        "Time_bins_detected",
        "Count_<2010",
        "Count_2010_2015",
        "Count_>2016",
        "Median_within_cluster_SNP",
        "Max_within_cluster_SNP",
        "Temporal_category_binned",
        "Evidence_level_binned"
    ]].head(25).to_string(index=False)
)

print("\n[DONE] Main table:")
print(f"{TABLE_DIR}/14C_cgmlst_cluster_binned_temporal_category_summary.csv")
