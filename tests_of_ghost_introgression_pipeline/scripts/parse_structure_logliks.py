#!/usr/bin/env python3
"""
Script: parse_structure_logliks.py
Author: Margaret Wanjiku
Purpose:
    Walks through STRUCTURE output directories and extracts log-likelihoods
    from each STRUCTURE run. Saves a summary CSV for further analysis.

Arguments:
    --indir: Path to base STRUCTURE output directory
    --outfile: Output CSV file (default: structure_loglik_summary.csv)

Output:
    CSV file containing:
        model, replicate, K, log_likelihood, file
"""

import os
import re
import pandas as pd
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--indir", required=True, help="Path to STRUCTURE output base directory")
parser.add_argument("--outfile", default="structure_loglik_summary.csv", help="Output CSV file")
args = parser.parse_args()

output_rows = []

print("Scanning for STRUCTURE output files...")

# Walk through all files in the output directory
for root, dirs, files in os.walk(args.indir):
    for fname in files:
        # Match STRUCTURE output files like structure_run_K{K}_f
        if fname.startswith("structure_run_K") and fname.endswith("_f"):
            full_path = os.path.join(root, fname)

            # Extract model, replicate, K from file path
            match = re.search(r'model(\d+)/replicate(\d+)/structure_run_K(\d+)_f', full_path)
            if not match:
                print(f"Skipped (pattern mismatch): {full_path}")
                continue

            model = f"model{match.group(1)}"
            replicate = f"replicate{match.group(2)}"
            K = int(match.group(3))

            # Extract log-likelihood from file
            with open(full_path, 'r') as f:
                lines = f.readlines()

            logL = None
            for line in lines:
                if "Estimated Ln Prob of Data" in line:
                    try:
                        logL = float(line.strip().split("=")[-1])
                    except ValueError:
                        print(f"Could not parse likelihood in: {full_path}")
                    break

            if logL is not None:
                output_rows.append({
                    "model": model,
                    "replicate": replicate,
                    "K": K,
                    "log_likelihood": logL,
                    "file": full_path
                })

# Save to CSV
if not output_rows:
    print("No valid STRUCTURE output files parsed.")
else:
    df = pd.DataFrame(output_rows)
    df = df.sort_values(by=["model", "replicate", "K"])
    df.to_csv(args.outfile, index=False)
    print(f"Done. Parsed log-likelihoods written to {args.outfile} with {len(df)} rows.")
