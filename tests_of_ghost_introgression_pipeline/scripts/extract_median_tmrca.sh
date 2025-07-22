#!/bin/bash
# Script: extract_median_tmrca.sh
# Author: Margaret Wanjiku
# Purpose:
#   Extracts the median TMRCA per site (5th column) from ARGweaver output.
# Usage:
#   bash extract_median_tmrca.sh input.tmrca.txt output.txt

set -euo pipefail

INFILE=$1
OUTFILE=$2

# Extract column 5 (median TMRCA) using awk
awk '{print $5}' "$INFILE" > "$OUTFILE"
