#!/bin/bash
# Usage: bash run_argweaver.sh input.fasta output.tmrca.txt

set -euo pipefail

FASTA=$1
OUTFILE=$2

# Remove existing stats file if any
rm -f "${FASTA}.arg.stats"

# Run Argweaver sampling
arg-sample --fasta "$FASTA" --output "${FASTA}.arg" --sample-step 100 --verbose 0 --overwrite

# Extract TMRCA
arg-extract-tmrca "${FASTA}.arg.%d.smc.gz" > "$OUTFILE"

