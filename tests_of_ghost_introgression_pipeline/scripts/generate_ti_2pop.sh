#!/bin/bash
# Script: generate_ti_2pop.sh
# Author: Margaret Wanjiku
# Purpose:
#   Generates .ti files from IMa3 2-pop input files.
# Usage:
#   bash generate_ti_2pop.sh input_dir output_dir ima3_path ncores burnin length interval chains hn ha hb

set -euo pipefail

INPUT_DIR=$1
OUTPUT_DIR=$2
IMA3=$3
NCORES=${4:-2}
BURNIN=${5:-10000}
LENGTH=${6:-1000}
INTERVAL=${7:-200}
CHAIN_SWAP=${8:-1}
HN=${9:-20}
HA=${10:-0.95}
HB=${11:-0.85}

mkdir -p "$OUTPUT_DIR"

for FILE in "$INPUT_DIR"/*.u; do
    BASENAME=$(basename "$FILE" .u_2pop.u)
    echo "→ Sampling genealogies for $BASENAME using $NCORES cores"

    mpirun -np "$NCORES" "$IMA3" \
        -i "$FILE" \
        -o "${OUTPUT_DIR}/${BASENAME}.u_2pop_null.out" \
        -q 100 -m 5.5 -t 5.5 \
        -b "$BURNIN" -L "$LENGTH" -d "$INTERVAL" -p "$CHAIN_SWAP" -r 245 \
        -hn "$HN" -ha "$HA" -hb "$HB"
done
