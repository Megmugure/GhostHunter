#!/bin/bash
# Usage: bash extract_median_tmrca.sh input.tmrca.txt output.txt

set -euo pipefail

INFILE=$1
OUTFILE=$2

awk '{print $5}' "$INFILE" > "$OUTFILE"

