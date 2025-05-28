#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile_summary <- args[2]
outfile_pdf <- args[3]
outfile_png <- args[4]
outfile_stats <- args[5]

library(diptest)
library(ggplot2)

tmrca <- scan(infile)

# KS Test vs exponential
ks_result <- ks.test(tmrca, "pexp", rate = 1/mean(tmrca))

# Dip Test
dip_result <- dip.test(tmrca)

# Summary stats
summary_stats <- summary(tmrca)
writeLines(capture.output(summary_stats), outfile_stats)

# Write test results
result_df <- data.frame(
  KS_D = ks_result$statistic,
  KS_p = ks_result$p.value,
  Dip_D = dip_result$statistic,
  Dip_p = dip_result$p.value,
  N = length(tmrca)
)
write.csv(result_df, file = outfile_summary, row.names = FALSE)

# Plot PDFs
pdf(outfile_pdf, width=8, height=6)
hist(tmrca, breaks=50, col="skyblue", main="TMRCA Distribution", xlab="TMRCA")
dev.off()

png(outfile_png, width=800, height=600)
hist(tmrca, breaks=50, col="skyblue", main="TMRCA Distribution", xlab="TMRCA")
dev.off()

