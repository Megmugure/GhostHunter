#!/usr/bin/env Rscript
# Script: run_modality_tests.R
# Author: Margaret Wanjiku
# Purpose:
#   Performs statistical modality testing on TMRCA values from ARGweaver output.
#   - Kolmogorov-Smirnov test vs exponential distribution
#   - Hartigan's Dip Test for multimodality
#   Also generates histogram plots and summary statistics.
#
# Usage:
#   Rscript run_modality_tests.R input_medians.txt summary.csv hist.pdf hist.png stats.txt

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]            # TMRCA median values
outfile_summary <- args[2]   # Summary CSV of statistical test results
outfile_pdf <- args[3]       # Histogram PDF output
outfile_png <- args[4]       # Histogram PNG output
outfile_stats <- args[5]     # Raw TMRCA stats summary

library(diptest)     # For Hartigan's dip test
library(ggplot2)     # For plotting

# Load median TMRCA values
tmrca <- scan(infile)

# Kolmogorov–Smirnov test against exponential distribution
ks_result <- ks.test(tmrca, "pexp", rate = 1/mean(tmrca))

# Hartigan's Dip Test for multimodality
dip_result <- dip.test(tmrca)

# Write raw summary stats to text
summary_stats <- summary(tmrca)
writeLines(capture.output(summary_stats), outfile_stats)

# Compile summary table for export
result_df <- data.frame(
  KS_D = ks_result$statistic,
  KS_p = ks_result$p.value,
  Dip_D = dip_result$statistic,
  Dip_p = dip_result$p.value,
  N = length(tmrca)
)
write.csv(result_df, file = outfile_summary, row.names = FALSE)

# Plot histogram to PDF
pdf(outfile_pdf, width = 8, height = 6)
hist(tmrca, breaks = 50, col = "skyblue", main = "TMRCA Distribution", xlab = "TMRCA")
dev.off()

# Plot histogram to PNG
png(outfile_png, width = 800, height = 600)
hist(tmrca, breaks = 50, col = "skyblue", main = "TMRCA Distribution", xlab = "TMRCA")
dev.off()
