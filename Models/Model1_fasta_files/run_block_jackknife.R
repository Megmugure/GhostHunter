# Load necessary libraries
library(diptest)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# The first argument is the directory with .tmrca.txt files
model_dir <- args[1]

# The second argument is the directory where results will be saved
results_dir <- args[2]

# Number of replicates for the jackknife
reps <- as.numeric(args[3])

set.seed(42)  # For reproducibility

# List all .tmrca.txt files in the model directory
all_files <- list.files(model_dir, pattern = "*.tmrca.txt", full.names = TRUE)
total <- length(all_files)
block_size <- ceiling(total * 0.10)  # Dropping 10% of the files

cat("Total files:", total, "\n")
cat("Block size (files dropped each replicate):", block_size, "\n")

# Initialize vectors to store test results
dip_stats <- numeric(reps)
dip_pvals <- numeric(reps)
ks_stats <- numeric(reps)
ks_pvals <- numeric(reps)

# Loop over the replicates
for (i in 1:reps) {
  cat("Running replicate", i, "\n")
  
  # Randomly drop block_size files (block jackknife)
  dropped <- sample(seq_along(all_files), block_size)
  kept_files <- all_files[-dropped]
  
  # Extract the 5th column (TMRCA values) from the kept files
  tmrca_values <- unlist(lapply(kept_files, function(f) {
    values <- tryCatch(scan(f, what = numeric(), quiet = TRUE), error = function(e) NA)
    if (length(values) >= 5) {
      return(values[5])  # Return the 5th value (TMRCA)
    } else {
      cat("Skipping file:", f, "due to insufficient data.\n")
      return(NA)  # Return NA if the file doesn't have enough values
    }
  }))
  
  # Remove any NA values (from files that had insufficient data)
  tmrca_values <- na.omit(tmrca_values)

  # Debugging: Print the number of retained files and their data
  cat("Retained files:", length(kept_files), "\n")
  cat("Number of TMRCA values:", length(tmrca_values), "\n")
  if (length(tmrca_values) > 0) {
    cat("First few TMRCA values for replicate", i, ":", head(tmrca_values), "\n")
  }

  # Check if the TMRCA values have variation (i.e., they are not all the same)
  if (length(tmrca_values) < 2) {
    cat("Skipping replicate", i, "due to insufficient data.\n")
    next
  }

  # Check if the values are all the same
  if (length(unique(tmrca_values)) == 1) {
    cat("Skipping replicate", i, "because TMRCA values are constant.\n")
    next
  }

  # KS test: Does the distribution match exponential?
  ks_test <- ks.test(tmrca_values, "pexp", rate = 1 / mean(tmrca_values))
  ks_stats[i] <- ks_test$statistic
  ks_pvals[i] <- ks_test$p.value

  # Dip test: Is it unimodal?
  dip_test <- dip.test(tmrca_values)
  dip_stats[i] <- dip_test$statistic
  dip_pvals[i] <- dip_test$p.value
}

# Save the results to a CSV file
results <- data.frame(
  Replicate = 1:reps,
  KS_Statistic = ks_stats,
  KS_PValue = ks_pvals,
  Dip_Statistic = dip_stats,
  Dip_PValue = dip_pvals
)

# Specify the output file path
output_file <- file.path(results_dir, "block_jackknife_results.csv")

# Write the results to a CSV file
write.csv(results, output_file, row.names = FALSE)

cat("Results saved to", results_dir, "\n")
