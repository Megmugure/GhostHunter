#!/bin/bash

# Concatenate median TMRCA values from all .tmrca.txt files
# Author: mwanjiku
# Date: April 2025
# Project: Ghost Populations - Coalescent Time Modality Tests

# Change to working directory with .tmrca.txt files
cd /usr/scratch/userdata/mwanjiku/ghost-pop-gen/Models/Model4_fasta_files || {
    echo "Directory not found. Exiting."
    exit 1
}

# Output file name
OUTPUT_FILE="all_model4_median_tmrca_values.txt"

# Remove existing output to prevent appending duplicates
rm -f "$OUTPUT_FILE"

# Loop through files and extract 5th column (median)
for f in *.tmrca.txt; do
    echo "Extracting median TMRCA from $f"
    awk '{print $5}' "$f" >> "$OUTPUT_FILE"
done

# Confirm line count
echo "-----------------------------------------"
echo "Total number of median TMRCA values:"
wc -l "$OUTPUT_FILE"
echo "-----------------------------------------"
echo "All values saved to $OUTPUT_FILE"

