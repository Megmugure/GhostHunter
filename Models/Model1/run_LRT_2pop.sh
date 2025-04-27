#!/bin/bash

# Description:
# This script runs 2-population likelihood ratio tests using IMa3
# on all input files in the ./original.u_Files_2pop/ directory with the .u extension.

# Requirements:
# - IMa3 must be compiled and located at ../../IMa3/IMa3 relative to the script
# - Input files must have a .u extension
# - Output will be saved in ./LRT_outfiles_2pop/

# Loop over all .u files in the specified input directory
for FILE in ./original.u_Files_2pop/*.u; do
    echo "Processing file: $FILE"

    # Strip the directory and .u extension to get a base name for the output
    BASENAME=$(basename "$FILE" .u)
    echo "Running 2-pop LRT on: $BASENAME"

    # Run the IMa3 program with specified parameters using mpirun
    mpirun -np 1 ../../IMa3/IMa3 \
        -i "$FILE" \
        -r0 \
        -v "./geneology.ti_Files_2pop/" \
        -L 5000 \
        -m 5.5 \
        -q 100 \
        -t 5.5 \
        -o "./LRT_outfiles_2pop/${BASENAME}.2pop.LRT.out"
done

echo "2-pop LRT analysis completed for all files."
