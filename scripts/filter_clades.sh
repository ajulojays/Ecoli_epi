#!/bin/bash

# 1. Correct paths based on your terminal output
# UPDATE THIS LINE if data6e.csv is saved somewhere else!
METADATA="../data/data6e.csv" 

CLADES="../results/12d_clade_analysis/all_preharvest_704_clades.csv"
OUTPUT="../results/12d_clade_analysis/data6e_704_genomes_with_clades.csv"

echo "Filtering metadata down to 704 genomes and adding clades..."

# 2. Use awk with comma separators (-F',')
awk -F',' -v OFS=',' '
    # First pass: read the clade info file
    NR==FNR {
        if (FNR==1) next # Skip the header of the clades file
        clade_map[$1] = $2
        next
    }
    # Second pass: read data6e
    FNR==1 {
        # Print the original header and append "Clade" to the end
        print $0, "Clade"
        next
    }
    # For all other rows, check if Genome_ID (column 1) is in our 704 list
    {
        if ($1 in clade_map) {
            # If it is, print the whole row plus the matched Clade
            print $0, clade_map[$1]
        }
    }
' "$CLADES" "$METADATA" > "$OUTPUT"

echo "Done! Saved to $OUTPUT"
