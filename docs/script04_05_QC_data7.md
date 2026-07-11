# Script 04 and Script 05: Assembly QC and data7 creation

Script 04:
scripts/04_run_assembly_QC_checkm2.sh

Purpose:
- Reads downloaded public assembly FASTA files from data6
- Calculates FASTA-derived assembly size, GC%, contig count, N50, longest contig, and N bases
- Stages normalized FASTA files for CheckM2
- Runs CheckM2 completeness and contamination estimation

Default local output:
~/epi/marker_screen/results/04_assembly_QC_checkm2/

Script 05:
scripts/05_merge_data6_QC_create_data7.R

Purpose:
- Reads cattle O157:H7 data6
- Reads FASTA assembly QC metrics from Script 04
- Reads CheckM2 completeness and contamination results
- Merges QC metrics with metadata
- Creates data7 using strict assembly-level QC filters

Default QC thresholds:
- CheckM2 contamination < 1%
- CheckM2 completeness >= 95%
- Genome size between 4.5 Mb and 6.0 Mb

Default local output:
~/epi/marker_screen/results/05_merge_QC_create_data7/

GitHub tracking:
- Scripts are committed.
- Small summary result tables can be committed.
- Raw FASTA files and full CheckM2 working folders should remain local.
