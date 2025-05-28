#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input CSV file from STRUCTURE log-likelihoods")
    parser.add_argument("--output", default="bootstrap_lrt_results.csv")
    parser.add_argument("--bootstraps", type=int, default=100, help="Number of bootstrap replicates")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df = df.dropna(subset=["model", "replicate", "K", "log_likelihood"])

    lrt_results = []

    for model in df["model"].unique():
        model_df = df[df["model"] == model]
        for replicate in model_df["replicate"].unique():
            rep_df = model_df[model_df["replicate"] == replicate]
            k_values = sorted(rep_df["K"].unique())

            for k in k_values:
                k_next = k + 1
                if k_next not in rep_df["K"].values:
                    continue

                row_k = rep_df[rep_df["K"] == k]
                row_k1 = rep_df[rep_df["K"] == k_next]

                if row_k.empty or row_k1.empty:
                    continue

                loglik_k = row_k["log_likelihood"].values[0]
                loglik_k1 = row_k1["log_likelihood"].values[0]
                T_obs = -2 * (loglik_k - loglik_k1)

                boot_stats = [
                    -2 * (np.random.choice(row_k["log_likelihood"].values) - 
                          np.random.choice(row_k["log_likelihood"].values))
                    for _ in range(args.bootstraps)
                ]

                p_val = np.mean(np.array(boot_stats) > T_obs)

                lrt_results.append({
                    "model": model,
                    "replicate": replicate,
                    "K0": k,
                    "K1": k_next,
                    "loglik_K0": loglik_k,
                    "loglik_K1": loglik_k1,
                    "T_obs": T_obs,
                    "p_value": p_val
                })

    pd.DataFrame(lrt_results).to_csv(args.output, index=False)
    print(f"Bootstrap LRT results written to: {args.output}")

if __name__ == "__main__":
    main()
