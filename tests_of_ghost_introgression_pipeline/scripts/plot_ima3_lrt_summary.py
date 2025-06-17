#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# CONFIGURATION
RESULTS_DIR = "results/ima3"
OUTDIR = RESULTS_DIR

# LOAD DATA
lrt_2pop = pd.read_csv(os.path.join(RESULTS_DIR, "All_LRT_results_2pop.csv"))
lrt_3pop = pd.read_csv(os.path.join(RESULTS_DIR, "All_LRT_results_3pop.csv"))

lrt_2pop["Test"] = "Ghost Population"
lrt_3pop["Test"] = "Ghost Gene Flow"

lrt_all = pd.concat([lrt_2pop, lrt_3pop], ignore_index=True)
lrt_all["Model"] = pd.Categorical(
    lrt_all["Model"], categories=sorted(lrt_all["Model"].unique()), ordered=True
)
lrt_all["-log10(p-value)"] = -np.log10(lrt_all["p-value"].replace(0, 1e-10))

# PLOT STYLE
sns.set_theme(style="whitegrid", context="talk", font_scale=1.1)

# Boxplot of LRT values
plt.figure(figsize=(14, 6))
sns.boxplot(data=lrt_all, x="Model", y="2LLR", hue="Test")
plt.axhline(y=9.21, color='gray', linestyle='--', linewidth=1, label="p=0.01 threshold (df=2)")
plt.title("LRT Values by Model and Test Type")
plt.ylabel("LRT Statistic (2LLR)")
plt.legend(title="Test Type")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/plot1_boxplot_2llr_by_model.png", dpi=300)
print("plot1_boxplot_2llr_by_model.png saved.")

# Barplot of -log10(p-value)
plt.figure(figsize=(14, 6))
sns.barplot(data=lrt_all, x="Model", y="-log10(p-value)", hue="Test", errorbar=None)
plt.axhline(y=-np.log10(0.05), color="red", linestyle="--", label="p = 0.05")
plt.title("Significance of LRTs by Model")
plt.ylabel("-log10(p-value)")
plt.legend(title="Test Type")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/plot2_barplot_logpval_by_model.png", dpi=300)
print("plot2_barplot_logpval_by_model.png saved.")

# FacetGrid: 2LLR by model and test
g = sns.catplot(
    data=lrt_all, kind="box",
    x="Test", y="2LLR", hue="Test", col="Model", col_wrap=3,
    height=4, aspect=1.1, palette="muted", legend=False
)
for ax in g.axes.ravel():
    ax.axhline(9.21, linestyle="--", color="gray", linewidth=1)
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle("LRT 2LLR by Test Type Across Models")
plt.savefig(f"{OUTDIR}/plot3_facetgrid_boxplot_by_model.png", dpi=300)
print("plot3_facetgrid_boxplot_by_model.png saved.")

# Stripplot: replicate-level p-values
plt.figure(figsize=(14, 6))
sns.stripplot(data=lrt_all, x="Model", y="-log10(p-value)", hue="Test", dodge=True, jitter=True, alpha=0.6)
plt.axhline(-np.log10(0.05), color="red", linestyle="--", label="p = 0.05")
plt.title("Per-Replicate LRT Significance")
plt.ylabel("-log10(p-value)")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUTDIR}/plot4_stripplot_logpval_by_replicate.png", dpi=300)
print("plot4_stripplot_logpval_by_replicate.png saved.")

# Heatmap: average -log10(p) per model & test
heatmap_df = lrt_all.groupby(["Model", "Test"], observed=True)["-log10(p-value)"].mean().unstack()
plt.figure(figsize=(8, 6))
sns.heatmap(heatmap_df, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': '-log10(p-value)'})
plt.title("Mean -log10(p-value) per Model and Test Type")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/plot5_heatmap_avg_logpval.png", dpi=300)
print("plot5_heatmap_avg_logpval.png saved.")

# Violin plot: distribution of LRT values
plt.figure(figsize=(14, 6))
sns.violinplot(data=lrt_all, x="Model", y="2LLR", hue="Test", split=True, inner="box", palette="Set2")
plt.axhline(9.21, color="gray", linestyle="--", linewidth=1)
plt.title("LRT Statistic (2LLR) Distribution with Threshold")
plt.ylabel("2LLR")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/plot6_violin_2llr_by_model.png", dpi=300)
print("plot6_violin_2llr_by_model.png saved.")

# Done
print("All plots generated and saved successfully.")

