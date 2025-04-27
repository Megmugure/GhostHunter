#!/bin/bash

# Description:
# This script runs 3-population likelihood ratio tests using IMa3
# on all input files in the ./original.u_Files_3pop/ directory with the .u extension.

# Requirements:
# - IMa3 must be compiled and located at ../../IMa3/IMa3 relative to the script
# - Genealogy files must be present in ./geneology.ti_Files_3pop/ and named like: <input>.u.out.ti
# - Input files must have a .u extension
# - Output will be saved in ./LRT_outfiles_3pop/

# Create output directory if it doesn't exist
mkdir -p ./LRT_outfiles_3pop/

# Loop over all .u files in the specified input directory
for FILE in ./original.u_Files_3pop/*.u; do
    echo "Processing file: $FILE"

    # Get the base name without directory and extension
    BASENAME=$(basename "$FILE" .u)
    echo "Running 3-pop LRT on: $BASENAME"

    # Define expected .ti file
    TI_FILE="./geneology.ti_Files_3pop/${BASENAME}.u.out.ti"

    # Check if the genealogy .ti file exists
    if [[ ! -f "$TI_FILE" ]]; then
        echo "Skipping $BASENAME — missing genealogy file: $TI_FILE"
        continue
    fi

    # Run the IMa3 program with specified parameters using mpirun
    mpirun -np 1 ../../IMa3/IMa3 \
        -i "$FILE" \
        -r0 \
        -v "$TI_FILE" \
        -L 5000 \
        -m 5.5 \
        -q 100 \
        -t 5.5 \
        -o "./LRT_outfiles_3pop/${BASENAME}.3pop.LRT.out"

    echo "Finished 3-pop LRT: $BASENAME"
done

echo "3-pop LRT analysis completed for all files."
