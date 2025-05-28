#!/usr/bin/env python3
import os
import re
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--indir", required=True)
parser.add_argument("--outfile", default="structure_loglik_summary.csv")
args = parser.parse_args()

output_rows = []

print("🔍 Scanning for STRUCTURE output files...")

for root, dirs, files in os.walk(args.indir):
    for fname in files:
        if fname.startswith("structure_run_K") and fname.endswith("_f"):
            full_path = os.path.join(root, fname)

            match = re.search(r'model(\d+)/replicate(\d+)/structure_run_K(\d+)_f', full_path)
            if not match:
                print(f"Skipped (pattern mismatch): {full_path}")
                continue

            model = f"model{match.group(1)}"
            replicate = f"replicate{match.group(2)}"
            K = int(match.group(3))

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

if not output_rows:
    print("❌ No valid STRUCTURE output files parsed.")
else:
    df = pd.DataFrame(output_rows)
    df = df.sort_values(by=["model", "replicate", "K"])
    df.to_csv(args.outfile, index=False)
    print(f"✅ Done. Parsed log-likelihoods written to {args.outfile} with {len(df)} rows.")

