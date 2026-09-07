Script 14B: cgMLST cluster temporal dynamics within 12D core-SNP clades

What this analysis is designed to achieve:
1. Restrict analysis to the exact non-dairy cattle pre-harvest genome set used in 12C/12D.
2. Use 12D core-SNP clades as broad genomic lineages.
3. Use cgMLST HierCC-style clusters as finer sublineages within each clade.
4. Summarize temporal behavior of each cgMLST cluster using first year, last year, year span, years detected, detection gaps, and consecutive-year detection.
5. Merge within-cluster SNP distances from 12E to check genomic coherence.
6. Merge AMRFinderPlus feature counts to describe AMR burden per cluster.
7. Classify each cgMLST cluster as persistent-like, emerging-like, re-emergence-like / unsampled persistence, sporadic/minority, or limited.
8. Assign every genome to the temporal category of its cgMLST cluster.

Interpretation caution:
These are sampling-aware categories. Absence in a year means not observed in public genomes, not confirmed biological absence.

Main output table:
tables/14B_cgmlst_cluster_temporal_category_summary.csv

Genome-level assignment table:
tables/14B_genome_level_temporal_category_assignments.csv
