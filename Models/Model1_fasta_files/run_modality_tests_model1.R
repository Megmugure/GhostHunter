# Load required libraries
# diptest: for modality testing using Hartigan's Dip Test
# ggplot2: for clean, publication-quality plots
library(diptest)
library(ggplot2)

# Create a directory to store all results and output files
results_dir <- "results_model1"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Read in the list of median TMRCA values from a plain text file
# Each value should be on a new line
tmrca <- scan("all_model1_median_tmrca_values.txt")

# Generate basic descriptive statistics (min, max, mean, quartiles)
# This helps understand the distribution shape before any tests
summary_stats <- summary(tmrca)
print(summary_stats)

# Save the summary statistics to a text file for documentation
writeLines(capture.output(summary_stats), file.path(results_dir, "tmrca_summary_statistics.txt"))

# Perform a Kolmogorov–Smirnov (KS) test
# Null hypothesis: the TMRCA values follow an exponential distribution (as expected under a standard coalescent model)
# Note: ties in data may affect the accuracy of this test
ks_result <- ks.test(tmrca, "pexp", rate = 1 / mean(tmrca))

# Perform Hartigan's Dip Test for unimodality vs multimodality
# Null hypothesis: the distribution is unimodal (single peak)
dip_result <- dip.test(tmrca)

# Combine test results and save them to a text file
# Useful for reproducibility and manuscript documentation
results <- c(
  "Kolmogorov–Smirnov Test",
  paste("D =", ks_result$statistic),
  paste("p-value =", ks_result$p.value),
  "",
  "Hartigan's Dip Test",
  paste("D =", dip_result$statistic),
  paste("p-value =", dip_result$p.value)
)
writeLines(results, file.path(results_dir, "modality_test_results_model1.txt"))

# Create and save a basic histogram using base R
# This gives a quick visual of the distribution shape
pdf(file.path(results_dir, "tmrca_distribution_model1.pdf"), width = 8, height = 6)
hist(tmrca, breaks = 50, col = "skyblue", main = "TMRCA Distribution", xlab = "Median TMRCA")
dev.off()

# Save the same histogram as a PNG for use in presentations or web
png(file.path(results_dir, "tmrca_distribution_model1.png"), width = 800, height = 600)
hist(tmrca, breaks = 50, col = "skyblue", main = "TMRCA Distribution", xlab = "Median TMRCA")
dev.off()

# Create a ggplot2 histogram with cleaner formatting for publication
# This is useful for inclusion in a paper or figure supplement
tmrca_df <- data.frame(TMRCA = tmrca)

pdf(file.path(results_dir, "tmrca_distribution_ggplot_model1.pdf"), width = 8, height = 6)
ggplot(tmrca_df, aes(x = TMRCA)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(
    title = "TMRCA Distribution (ggplot)",
    x = "Median TMRCA",
    y = "Frequency"
  )
dev.off()
