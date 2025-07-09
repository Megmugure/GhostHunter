#!/usr/bin/env python3

"""
visualize_bootstrap_lrt.py

This script generatess visualizations from the STRUCTURE LRT bootstrap analysis.

Input:
- `bootstrap_lrt_results_all.csv` containing:
    - model: model name
    - replicate: replicate ID
    - K: number of clusters tested
    - K+1: next cluster
    - T_obs: observed likelihood ratio statistic
    - p_value: bootstrap p-value for K→K+1

Outputs:
1. Boxplot of T_obs by K and model
2. Heatmap of p-values by replicate and K transition
3. Stacked bar plot of p-value categories by model and K transition

Each plot is saved as a high-resolution PNG (300 dpi).
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Data loading and preprocessing

# Load CSV
df = pd.read_csv("bootstrap_lrt_results_all.csv")

# Ensure numeric types
df['K'] = df['K'].astype(int)
df['K+1'] = df['K+1'].astype(int)

# Create K transition label, e.g., "2→3"
df['K_pair'] = df['K'].astype(str) + "→" + df['K+1'].astype(str)

# Boxplot of T_obs by K 

plt.figure(figsize=(12, 6))

# Show boxplot of T_obs across K, grouped by model
sns.boxplot(data=df, x='K', y='T_obs', hue='model', showfliers=False)

# Add jittered individual points for visual detail
sns.stripplot(data=df, x='K', y='T_obs', hue='model', dodge=True, alpha=0.4, linewidth=0.5)

# Reference line at 0 for visual aid
plt.axhline(0, color='black', linestyle='--', linewidth=1)

# Labels and formatting
plt.title('Distribution of LRT Statistics (T_obs) by K Across Models', fontsize=14)
plt.ylabel('LRT Statistic (T_obs)', fontsize=12)
plt.xlabel('K', fontsize=12)
plt.legend(title='Model', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(True, linestyle=':', linewidth=0.5)

# Save the figure
plt.tight_layout()
plt.savefig("plot_tobs_by_K.png", dpi=300)
plt.show()

# Heatmap of p-values per replicate

# Pivot data to heatmap shape: replicates x K transitions
heatmap_df = df.pivot_table(index='replicate', columns='K_pair', values='p_value')

plt.figure(figsize=(16, 8))

# Show heatmap with p-values annotated
sns.heatmap(heatmap_df, annot=True, fmt=".2f", cmap='coolwarm', linewidths=0.5, cbar_kws={'label': 'p-value'})

# Formatting
plt.title('p-value Heatmap by Replicate and K Transition', fontsize=14)
plt.ylabel('Replicate', fontsize=12)
plt.xlabel('K → K+1', fontsize=12)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("heatmap_pvalues_by_replicate.png", dpi=300)
plt.show()

# Stacked Bar Plot of Significance Categories

# Categorize p-values based on significance
def classify_p(p):
    if p < 0.05:
        return 'Significant (p < 0.05)'
    elif p > 0.95:
        return 'Highly non-significant (p > 0.95)'
    else:
        return 'Ambiguous'

# Apply classification
df['p_class'] = df['p_value'].apply(classify_p)

# Count proportions per (model, K_pair)
summary = df.groupby(['model', 'K_pair'])['p_class'].value_counts(normalize=True).unstack().fillna(0)

# Reorder columns for consistent color interpretation
class_order = ['Significant (p < 0.05)', 'Ambiguous', 'Highly non-significant (p > 0.95)']
summary = summary[class_order]

# Plot as stacked bar
summary.plot(kind='bar', stacked=True, figsize=(16, 6), color=sns.color_palette("Set2"))

# Labels and formatting
plt.title("Proportion of p-value Categories by Model and K Transition", fontsize=14)
plt.ylabel("Proportion", fontsize=12)
plt.xlabel("Model / K→K+1", fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.legend(title="p-value Category", bbox_to_anchor=(1.02, 1), loc='upper left')

# Save
plt.tight_layout()
plt.savefig("pvalue_category_proportions.png", dpi=300)
plt.show()
