#!/bin/bash

# Script to extract 2LLR values, degrees of freedom (df), 
# compute p-values, and save results to a text file.

# Output file name
output_file="LRT_results.txt"

# Write the header to the output file
echo -e "Filename\t2LLR\t df\t p-value" > "$output_file"

# Function to process files in a directory
process_files_in_directory() {
    dir="$1"  # Directory to process

    # Loop through all files in the directory
    for file in "$dir"/*; do
        # Only process files with .LRT.out extension
        if [[ $file =~ \.LRT\.out$ ]]; then
            echo "Processing file: $file"

            # Read each line in the file
            while IFS= read -r line; do
                # Use regular expression to extract the 2LLR statistic and degrees of freedom (df)
                if [[ $line =~ [[:space:]]+([0-9]+\.[0-9]+)[[:space:]]+([0-9]+)[[:space:]]+([0-9]+)[[:space:]]+[0-9]+\.[0-9]+ ]]; then
                    llr="${BASH_REMATCH[1]}"  # Extracted 2LLR value
                    df="${BASH_REMATCH[2]}"   # Extracted degrees of freedom

                    # Compute the p-value using Python's SciPy package
                    p_value=$(python3 -c "import scipy.stats as stats; print(stats.chi2.sf($llr, $df))")

                    # Append results to the output file, including the filename
                    echo -e "$file\t$llr\t $df\t $p_value" >> "$output_file"
                fi
            done < "$file"
        fi
    done
}

# Process the 3pop and 2pop directories
process_files_in_directory "LRT_outfiles_3pop"
process_files_in_directory "LRT_outfiles_2pop"

# Print confirmation message
echo "Results saved to $output_file"
L