#!/bin/bash
# Script: run_argweaver.sh
# Author: Margaret Wanjiku
# Purpose:
#   Runs ARGweaver's arg-sample and arg-extract-tmrca for a given FASTA.
# Usage:
#   bash run_argweaver.sh input.fasta output.tmrca.txt

set -euo pipefail

FASTA=$1
OUTFILE=$2

# Remove any existing stats file to avoid reuse
rm -f "${FASTA}.arg.stats"

# Run ARGweaver to sample ARGs
arg-sample --fasta "$FASTA" --output "${FASTA}.arg" --sample-step 100 --verbose 0 --overwrite

# Extract TMRCA per site
arg-extract-tmrca "${FASTA}.arg.%d.smc.gz" > "$OUTFILE"
