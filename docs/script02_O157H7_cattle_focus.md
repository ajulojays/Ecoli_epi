# Script 02: O157:H7 cattle focus

Script 02 creates data6, a cattle-only dataset containing the original data5 metadata joined to ECtyper metadata.

Input files:

- data5.csv
- results/00_ectyper_batch_run/tables/ectyper_output_all.tsv

Output folder:

- results/02_O157H7_cattle_focus/

Main output tables:

- tables/data6_cattle_with_ectyper_metadata.csv
- tables/data6_cattle_O157H7_only.csv
- tables/data6_cattle_serotype_counts.csv
- tables/data6_cattle_temporal_O157H7_summary.csv
- tables/data6_cattle_temporal_O157H7_long_counts.csv

Main figures:

- figures/cattle_O157H7_temporal_proportion.pdf
- figures/cattle_O157H7_temporal_proportion.png
- figures/cattle_O157H7_temporal_counts.pdf
- figures/cattle_O157H7_temporal_counts.png
- figures/cattle_O157H7_temporal_percent_stacked.pdf
- figures/cattle_O157H7_temporal_percent_stacked.png

Temporal grouping:

- <=2000
- 2001
- 2002
- 2003
- each year afterward

O157:H7 categories:

- O157:H7
- Non-O157:H7
- No serotype call

No serotype call is retained separately so untyped genomes are not incorrectly counted as non-O157:H7.

Additional figure update:

- The vertical n labels on temporal stacked bar plots are positioned higher above the bars.
- Script 02 also generates a temporal O157:H7 proportion line plot with 95% exact binomial confidence intervals.

Additional figure files:

- figures/cattle_O157H7_temporal_proportion_line_CI.pdf
- figures/cattle_O157H7_temporal_proportion_line_CI.png

The line plot uses:

- numerator: number of O157:H7 cattle genomes per temporal group
- denominator: all cattle genomes in data6 per temporal group
- CI method: exact binomial confidence interval from stats::binom.test

