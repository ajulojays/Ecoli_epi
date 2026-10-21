# cgMLST persistence across livestock core-SNP clades

This module tests whether temporally persistent cgMLST lineages are concentrated in particular phylogenetic clades in the **704 non-dairy cattle, pre-harvest E. coli O157:H7 genomes** used by the 12C/12D/14B workflow.

## What it does

1. Loads the exact 704-genome metadata and 95% cgMLST allele matrix.
2. Reconstructs study-defined single-linkage **HC10, HC20, and HC50** clusters.
3. Applies the same temporal-category logic used in Script 14B:
   - Persistent-like
   - Emerging-like
   - Re-emergence-like / unsampled persistence
   - Repeated detection, limited temporal evidence
   - Sporadic/minority
4. Maps each genome and cluster to its 12D core-SNP clade.
5. Tests **cluster-level enrichment** of persistent-like lineages in each clade using Fisher's exact test with Benjamini-Hochberg FDR correction.
6. Reports descriptive genome burden by clade.
7. Generates iTOL color-strip files for:
   - core-SNP clade
   - HC10 / HC20 / HC50
   - temporal category at each HC level

## Important scope

The preserved 12D clade mapping contains 704 tips:
- Clade_001: 322
- Clade_002: 194
- Clade_003: 181
- Clade_004: 4
- Clade_005: 3

This is the exact 704-genome pre-harvest, non-dairy analysis set. It is **not** the full 1,029-genome livestock tree. If the scientific target is the full 1,029-genome tree, HC assignments and temporal categories should be recalculated for that population rather than treating the 704-genome clade labels as if they were defined on all 1,029 genomes.

## Inputs

The clade colorstrip is bundled in `inputs/`.

Two required source files were not committed to the working GitHub branch and therefore must be collected from a full/restored project copy:

- `results/12c_all_preharvest_combined/target_metadata_704_preharvest_non_dairy.csv`
- `results/chewbbaca/cgmlst/core_95/cgMLST95.tsv`

Optional:
- `results/12d_clade_analysis/all_preharvest_704_clades.csv`
- `results/12c_all_preharvest_combined/iqtree/all_preharvest_704_core_phylo.treefile`

## Fast run

From this module:

```bash
bash collect_inputs.sh /path/to/Ecoli_epi
bash run.sh
```

The collector renames/copies the required inputs into `inputs/`. The runner creates a local virtual environment, installs NumPy/Pandas/SciPy, and executes the analysis.

## Clone only this folder

Git itself does not clone a single directory directly, but sparse checkout does:

```bash
git clone --filter=blob:none --no-checkout https://github.com/ajulojays/Ecoli_epi.git
cd Ecoli_epi
git sparse-checkout init --cone
git sparse-checkout set analysis/15_cgmlst_persistence_clade
git checkout main
cd analysis/15_cgmlst_persistence_clade
```

Then collect the two required data inputs from the restored/full project and run:

```bash
bash collect_inputs.sh /path/to/restored/Ecoli_epi
bash run.sh
```

## Main outputs

- `outputs/genome_level_hc_assignments.csv`
- `outputs/cluster_temporal_summary_HC10_HC20_HC50.csv`
- `outputs/genome_temporal_assignments_HC10_HC20_HC50.csv`
- `outputs/clade_temporal_genome_burden.csv`
- `outputs/clade_persistence_enrichment_cluster_level.csv`
- `outputs/itol/itol_HC10_colorstrip.txt`
- `outputs/itol/itol_HC20_colorstrip.txt`
- `outputs/itol/itol_HC50_colorstrip.txt`
- `outputs/itol/itol_temporal_category_HC10.txt`
- `outputs/itol/itol_temporal_category_HC20.txt`
- `outputs/itol/itol_temporal_category_HC50.txt`
- `outputs/itol/itol_coreSNP_clades_colorstrip.txt`

## Statistical interpretation

The enrichment analysis treats the **cgMLST cluster** as the inferential unit. This avoids counting many genomes from one persistent lineage as independent observations. Genome counts and percentages are reported separately as descriptive epidemiologic burden.

A low-HC cluster or persistent-like assignment is evidence of genomic/temporal relatedness among sampled genomes; it is not by itself proof of continuous biological persistence or direct transmission. Years without detections represent unsampled/non-observed periods, not demonstrated absence.
