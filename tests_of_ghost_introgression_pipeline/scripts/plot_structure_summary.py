#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--loglik", required=True, help="Log-likelihood summary CSV")
    parser.add_argument("--aicbic", required=True, help="AIC/BIC summary CSV")
    parser.add_argument("--lrt", required=True, help="Bootstrap LRT results CSV")
    parser.add_argument("--outdir", default="results", help="Directory to save plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    loglik_df = pd.read_csv(args.loglik)
    aic_bic_df = pd.read_csv(args.aicbic)
    lrt_df = pd.read_csv(args.lrt)

    # Plot 1: Log-Likelihood vs K per Model
    plt.figure(figsize=(10, 6))
    for model in loglik_df["model"].unique():
        df = loglik_df[loglik_df["model"] == model]
        mean_ll = df.groupby("K")["log_likelihood"].mean()
        plt.plot(mean_ll.index, mean_ll.values, marker='o', label=model)

    plt.title("STRUCTURE Log-Likelihood vs K")
    plt.xlabel("Number of Clusters (K)")
    plt.ylabel("Mean Log-Likelihood")
    plt.legend(title="Model", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "plot_loglik_vs_k.png"))
    plt.close()

    # Plot 2: AIC and BIC across K per Model
    plt.figure(figsize=(12, 6))
    for model in aic_bic_df["model"].unique():
        df = aic_bic_df[aic_bic_df["model"] == model]
        mean_aic = df.groupby("K")["AIC"].mean()
        mean_bic = df.groupby("K")["BIC"].mean()
        plt.plot(mean_aic.index, mean_aic.values, linestyle='--', label=f"{model} AIC")
        plt.plot(mean_bic.index, mean_bic.values, linestyle='-', label=f"{model} BIC")

    plt.title("Model Selection using AIC and BIC")
    plt.xlabel("Number of Clusters (K)")
    plt.ylabel("Score")
    plt.legend(title="Criterion", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "plot_aic_bic_vs_k.png"))
    plt.close()

    # Plot 3: Bootstrap LRT p-values for K Comparisons
    plt.figure(figsize=(10, 6))
    mean_pvals = lrt_df.groupby(["K0", "K1"])["p_value"].mean().reset_index()
    plt.plot(mean_pvals["K1"], mean_pvals["p_value"], marker='o')
    plt.axhline(0.05, color='red', linestyle='--', label="p = 0.05 threshold")

    plt.title("Bootstrap LRT Mean p-values Across K")
    plt.xlabel("Tested K (K1)")
    plt.ylabel("Mean p-value")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "plot_bootstrap_lrt_pvalues.png"))
    plt.close()

    print("Summary plots saved to:", args.outdir)

if __name__ == "__main__":
    main()

