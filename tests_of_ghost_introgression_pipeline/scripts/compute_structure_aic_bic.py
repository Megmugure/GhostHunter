#!/usr/bin/env python3
"""
Script: compute_structure_aic_bic.py
Author: Margaret Wanjiku
Purpose:
    Calculates AIC and BIC values for STRUCTURE runs based on
    log-likelihoods, number of individuals, number of loci, and K.

Arguments:
    --indir: Directory containing STRUCTURE output files
    --strdir: Directory containing corresponding .str input files
    --outfile: Output CSV file (default: structure_aic_bic_summary.csv)

Output:
    CSV file containing model, replicate, K, lnL, I, A, p, AIC, BIC
"""

import os
import re
import pandas as pd
import numpy as np
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--indir", required=True, help="Path to STRUCTURE output directory")
parser.add_argument("--strdir", required=True, help="Path to .str input files directory")
parser.add_argument("--outfile", default="structure_aic_bic_summary.csv", help="Output CSV file")
args = parser.parse_args()

results = []

# Loop through model and replicate directories
for model in sorted(os.listdir(args.indir)):
    model_path = os.path.join(args.indir, model)
    if not os.path.isdir(model_path):
        continue

    for replicate in sorted(os.listdir(model_path)):
        rep_path = os.path.join(model_path, replicate)
        if not os.path.isdir(rep_path):
            continue

        replicate_num_match = re.search(r'\d+', replicate)
        model_num_match = re.search(r'\d+', model)
        if not replicate_num_match or not model_num_match:
            continue

        replicate_num = replicate_num_match.group()
        model_num = model_num_match.group()

        # Locate corresponding .str input file
        str_file = os.path.join(args.strdir, f"model{model_num}_replicate{replicate_num}_cleaned.str")
        if not os.path.exists(str_file):
            print(f"Missing .str file: {str_file}")
            continue

        try:
            with open(str_file) as f:
                lines = [line.strip() for line in f if line.strip()]
            I = len(lines)
            num_loci = len(lines[0].split()) - 1
            A = num_loci * (2 - 1)  # Assume haploid data
        except Exception as e:
            print(f"Error parsing {str_file}: {e}")
            continue

        # Loop over STRUCTURE output files
        for file in os.listdir(rep_path):
            match = re.match(r"structure_run_K(\d+)_f", file)
            if not match:
                continue

            K = int(match.group(1))
            output_file = os.path.join(rep_path, file)

            try:
                with open(output_file, 'r') as f:
                    lines = f.readlines()
                lnL_line = [l for l in lines if "Estimated Ln Prob of Data" in l][-1]
                lnL = float(re.search(r"-?\d+\.\d+", lnL_line).group())
            except Exception as e:
                print(f"Error reading {output_file}: {e}")
                continue

            # Calculate number of parameters
            p = I * (K - 1) + K * A
            aic = -2 * lnL + 2 * p
            bic = -2 * lnL + p * np.log(I)

            results.append({
                "model": f"model{model_num}",
                "replicate": int(replicate_num),
                "K": K,
                "lnL": lnL,
                "I": I,
                "A": A,
                "p": p,
                "AIC": aic,
                "BIC": bic
            })

# Save results
df = pd.DataFrame(results)

if not df.empty:
    df.sort_values(by=["model", "replicate", "K"], inplace=True)
    df.to_csv(args.outfile, index=False)
    print(f"Done! Output saved to {args.outfile}")
else:
    print("No data found. Check if STRUCTURE outputs and .str files were correctly matched.")
