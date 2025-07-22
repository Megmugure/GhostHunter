#!/bin/bash
# Script: run_lrt_2pop.sh
# Author: Margaret Wanjiku
# Purpose:
#   Runs LRT (likelihood ratio test) for IMa3 2-pop models.
# Usage:
#   bash run_lrt_2pop.sh input_dir ti_dir output_dir model_file ima3_path ncores

set -euo pipefail

INPUT_DIR=$1
TI_DIR=$2
OUT_DIR=$3
MODEL_FILE=$4
IMA3=$5
NCORES=${6:-1}

mkdir -p "$OUT_DIR"

for FILE in "$INPUT_DIR"/*.u; do
    BASENAME=$(basename "$FILE" .u_2pop.u)
    TI_FILE="${TI_DIR}/${BASENAME}.u_2pop_null.out.ti"
    OUT_FILE="${OUT_DIR}/${BASENAME}.2pop.LRT.out"

    echo "→ Running LRT for $BASENAME"
    mpirun -np "$NCORES" "$IMA3" \
        -i "$FILE" -r 0 -v "$TI_FILE" -w "$MODEL_FILE" \
        -c 2 -L 5000 -m 5.5 -q 100 -t 5.5 \
        -o "$OUT_FILE"
done
