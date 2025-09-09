#!/usr/bin/env python3
"""
Purpose
-------
Generate quick QC plots and a small numerical summary for PPP windowed Weir–FST
(all windows) vs the uniformly sampled windows used downstream.

Usage
-----
python scripts/qc_fst_plots.py \
  --filtered data/fst/CEU_CHS.allchr.windowed.weir.filtered.fst \
  --sampled  data/bins/CEU_CHS.uniform3.N300.windowed.weir.fst \
  --outdir   figs \
  --chr      1

Inputs
------
- Filtered FST table (all chromosomes), columns:
  CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST, MEAN_FST
- Sampled windows table with the same leading columns.

Outputs
-------
- figs/fst_hist_all_vs_sampled.png
- figs/fst_ecdf_all_vs_sampled.png
- figs/fst_track_chr{chr}.png
- figs/sampled_windows_per_chrom.png
- figs/fst_vs_nvariants_hexbin.png
- figs/summary.txt   # counts and key quantiles

Notes
-----
- Tries to be resilient to slight header differences by renaming first 6 columns.
- Uses default matplotlib; no external style dependencies.
"""

import argparse, os, math
import pandas as pd
import matplotlib.pyplot as plt

def ecdf(x):
    """Return ECDF x,y for a 1D series (NaNs dropped)."""
    x = pd.Series(x).dropna().sort_values().values
    y = (pd.Series(range(1, len(x)+1)) / len(x)).values
    return x, y

def midpoints(df):
    """Midpoint (bp) of each window: (BIN_START + BIN_END)/2."""
    return (df["BIN_START"].astype(float) + df["BIN_END"].astype(float)) / 2.0

def main():
    p = argparse.ArgumentParser(description="QC plots for PPP windowed Fst + sampled windows")
    p.add_argument("--filtered", default="data/fst/CEU_CHS.allchr.windowed.weir.filtered.fst",
                   help="Filtered windowed Fst table (all chromosomes)")
    p.add_argument("--sampled", default="data/bins/CEU_CHS.uniform3.N300.windowed.weir.fst",
                   help="Sampled windows table")
    p.add_argument("--outdir", default="figs", help="Output figures directory")
    p.add_argument("--chr", default="1", help="Chromosome to draw a track for (e.g. 1)")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Read tables
    cols = ["CHROM","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST"]
    all_df = pd.read_csv(args.filtered, sep="\t")
    if list(all_df.columns[:6]) != cols:
        # handle if header slightly differs; try rename safely
        rename_map = dict(zip(all_df.columns[:6], cols))
        all_df = all_df.rename(columns=rename_map)
    samp_df = pd.read_csv(args.sampled, sep="\t")
    if list(samp_df.columns[:6]) != cols:
        rename_map = dict(zip(samp_df.columns[:6], cols))
        samp_df = samp_df.rename(columns=rename_map)

    # 1) Histogram of WEIGHTED_FST (all vs sampled)
    plt.figure()
    all_df["WEIGHTED_FST"].plot(kind="hist", bins=60, density=True, label="all", alpha=0.5)
    samp_df["WEIGHTED_FST"].plot(kind="hist", bins=60, density=True, label="sampled", alpha=0.5)
    plt.xlabel("Weighted FST per 50kb window")
    plt.ylabel("Density")
    plt.legend()
    plt.title("FST distribution: all windows vs sampled")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "fst_hist_all_vs_sampled.png"), dpi=200)
    plt.close()

    # 2) ECDF (all vs sampled)
    x_all, y_all = ecdf(all_df["WEIGHTED_FST"])
    x_s, y_s = ecdf(samp_df["WEIGHTED_FST"])
    plt.figure()
    plt.plot(x_all, y_all, label="all")
    plt.plot(x_s, y_s, label="sampled")
    plt.xlabel("Weighted FST")
    plt.ylabel("ECDF")
    plt.legend()
    plt.title("ECDF of FST: all windows vs sampled")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "fst_ecdf_all_vs_sampled.png"), dpi=200)
    plt.close()

    # 3) Manhattan-like track for one chromosome (all windows), highlight sampled ones
    chr_sel = str(args.chr)
    all_chr = all_df[all_df["CHROM"].astype(str) == chr_sel].copy()
    samp_chr = samp_df[samp_df["CHROM"].astype(str) == chr_sel].copy()
    if len(all_chr) > 0:
        plt.figure()
        plt.scatter(midpoints(all_chr)/1e6, all_chr["WEIGHTED_FST"], s=6, label="all")
        if len(samp_chr) > 0:
            plt.scatter(midpoints(samp_chr)/1e6, samp_chr["WEIGHTED_FST"], s=14, marker="x", label="sampled")
        plt.xlabel(f"Chromosome {chr_sel} position (Mb)")
        plt.ylabel("Weighted FST")
        plt.title(f"Windowed FST track (chr{chr_sel})")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, f"fst_track_chr{chr_sel}.png"), dpi=200)
        plt.close()

    # 4) Count sampled windows per chromosome (bar)
    samp_counts = samp_df.groupby("CHROM").size().reset_index(name="n")
    samp_counts = samp_counts.sort_values("CHROM", key=lambda s: s.astype(str).str.replace("chr","", regex=False).astype(int))
    plt.figure()
    plt.bar(samp_counts["CHROM"].astype(str), samp_counts["n"])
    plt.xlabel("Chromosome")
    plt.ylabel("# sampled windows")
    plt.title("Sampled windows per chromosome")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "sampled_windows_per_chrom.png"), dpi=200)
    plt.close()

    # 5) Relationship: N_VARIANTS vs WEIGHTED_FST (all windows)
    plt.figure()
    plt.hexbin(all_df["N_VARIANTS"], all_df["WEIGHTED_FST"], gridsize=50)
    plt.xlabel("N_VARIANTS (per 50kb window)")
    plt.ylabel("Weighted FST")
    plt.title("N_VARIANTS vs FST (all windows)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "fst_vs_nvariants_hexbin.png"), dpi=200)
    plt.close()

    # 6) Print a small text summary (helps Methods + sanity checks)
    q = all_df["WEIGHTED_FST"].quantile([0.33, 0.66]).values
    summary = {
        "n_all_windows": int(len(all_df)),
        "n_sampled": int(len(samp_df)),
        "mean_fst_all": float(all_df["WEIGHTED_FST"].mean()),
        "mean_fst_sampled": float(samp_df["WEIGHTED_FST"].mean()),
        "fst_33pct_all": float(q[0]),
        "fst_66pct_all": float(q[1]),
    }
    with open(os.path.join(args.outdir, "summary.txt"), "w") as fh:
        for k,v in summary.items():
            fh.write(f"{k}\t{v}\n")

if __name__ == "__main__":
    main()
