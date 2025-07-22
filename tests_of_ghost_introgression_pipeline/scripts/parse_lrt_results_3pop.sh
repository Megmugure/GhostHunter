#!/bin/bash
# Script: parse_lrt_results_3pop.sh
# Author: Margaret Wanjiku
# Purpose:
#   Parses IMa3 3-pop LRT outputs and extracts LRT stats.
# Output:
#   results/ima3/All_LRT_results_3pop.csv

set -euo pipefail

OUTPUT_FILE="results/ima3/All_LRT_results_3pop.csv"
echo -e "Model,Replicate,Filename,2LLR,DF,p-value" > "$OUTPUT_FILE"

for file in results/ima3/LRT_outfiles_3pop/*.LRT.out; do
    [[ -e "$file" ]] || continue

    model=$(echo "$file" | grep -oP 'model\d+')
    replicate=$(echo "$file" | grep -oP 'replicate\d+')
    filename=$(basename "$file")

    while IFS= read -r line; do
        if [[ $line =~ [[:space:]]+([0-9]+\.[0-9]+)[[:space:]]+([0-9]+)[[:space:]]+([0-9]+)[[:space:]]+[0-9]+\.[0-9]+ ]]; then
            llr="${BASH_REMATCH[1]}"
            df="${BASH_REMATCH[2]}"
            pval=$(python3 -c "import scipy.stats as s; print(round(s.chi2.sf($llr, $df), 6))")
            echo -e "$model,$replicate,$filename,$llr,$df,$pval" >> "$OUTPUT_FILE"
            break
        fi
    done < "$file"
done

echo "Parsed 3pop LRT results saved to $OUTPUT_FILE"
