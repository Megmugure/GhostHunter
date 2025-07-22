#!/bin/bash
# Script: run_lrt_3pop.sh
# Author: Margaret Wanjiku
# Purpose:
#   Runs LRT for IMa3 3-population models.
# Usage:
#   bash run_lrt_3pop.sh input_dir ti_dir output_dir ima3_path ncores

set -euo pipefail

INPUT_DIR=$1
TI_DIR=$2
OUT_DIR=$3
IMA3=$4
NCORES=${5:-1}

mkdir -p "$OUT_DIR"

for FILE in "$INPUT_DIR"/*.u; do
    BASENAME=$(basename "$FILE" .u)
    TI_FILE="${TI_DIR}/${BASENAME}.u.out.ti"
    OUT_FILE="${OUT_DIR}/${BASENAME}.3pop.LRT.out"

    echo "→ Running 3pop LRT for $BASENAME"
    mpirun -np "$NCORES" "$IMA3" \
        -i "$FILE" -r 0 -v "$TI_FILE" \
        -L 5000 -m 5.5 -q 100 -t 5.5 \
        -o "$OUT_FILE"
done
