#!/bin/bash
# Robust parser for 2pop LRT results (from Snakemake directory layout)

set -euo pipefail

OUTPUT_FILE="results/ima3/All_LRT_results_2pop.csv"
echo -e "Model,Replicate,Filename,2LLR,DF,p-value" > "$OUTPUT_FILE"

for file in results/ima3/LRT_outfiles_2pop/*.LRT.out; do
    [[ -e "$file" ]] || continue

    model=$(echo "$file" | grep -oP 'model\d+')
    replicate=$(echo "$file" | grep -oP 'replicate\d+')
    filename=$(basename "$file")

    model_line=$(grep -P "^\s*2\s+-?\d+\.\d+\s+\d+\s+\d+\*?\s+\d+\.\d+" "$file" || true)

    if [[ -n "$model_line" ]]; then
        df=$(echo "$model_line" | awk '{gsub("\\*", "", $4); print $4}')
        llr=$(echo "$model_line" | awk '{print $5}')
        pval=$(python3 -c "import scipy.stats as s; print(round(s.chi2.sf($llr, $df), 6))")
        echo -e "$model,$replicate,$filename,$llr,$df,$pval" >> "$OUTPUT_FILE"
    else
        echo "⚠ No valid LRT summary found in $file" >&2
    fi
done

echo "✔ Parsed 2pop LRT results saved to $OUTPUT_FILE"
