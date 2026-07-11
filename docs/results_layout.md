# Results layout

Large analysis outputs are organized by script number under the working directory.

Working directory:

  ~/epi/marker_screen

Result root:

  ~/epi/marker_screen/results

Recommended layout:

  results/
    00_ectyper_batch_run/
      tables/
      batch_outputs/
      manifests/
      archives/
      logs/

    01_serotype_summary/
      tables/
      figures/
        top10_serotype_pies/
        O157H7_focus_pies/
        temporal_stacked_counts/
        temporal_stacked_percent/
        multipage_pdfs/

    02_downstream_selection/

Script 00 writes ECtyper run outputs under:

  results/00_ectyper_batch_run

Script 01 reads the final ECtyper table from:

  results/00_ectyper_batch_run/tables/ectyper_output_all.tsv

Script 01 writes tables and figures under:

  results/01_serotype_summary

This keeps raw serotyping outputs separate from downstream summary and plotting outputs.

## Script 02

Script 02 writes to:

results/02_O157H7_cattle_focus/

It creates data6, a cattle-only dataset containing the original data5 metadata joined to ECtyper metadata.

