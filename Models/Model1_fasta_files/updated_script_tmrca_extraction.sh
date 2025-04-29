#!/bin/bash

#PBS -V
#PBS -l nodes=1:ppn=2:mpi
#PBS -N tmrca_extraction_1
#PBS -joe
#PBS -q batch
#PBS -l walltime=10:00:00

# Load the Miniconda3 environment (adjust the path based on your setup)
source /home3/mwanjiku/miniconda3/bin/activate argweaver_py2

# Make sure the necessary commands are in the PATH (in case they are not automatically found)
export PATH="/home3/mwanjiku/miniconda3/envs/argweaver_py2/bin:$PATH"

# Change directory to the working directory where the job was submitted
cd $PBS_O_WORKDIR

# Get the number of cores assigned to the job
NCORES=$(wc -w < $PBS_NODEFILE)

# Capture the current date and hostname for logging purposes
DATE=$(date)
STARTHOST=$(hostname)

# Output the start details for logging
echo "Running on host: $STARTHOST"
echo "Job submitted: $DATE"

# Define the directory containing the FASTA files (adjust the path accordingly)
DIR="/usr/scratch/userdata/mwanjiku/ghost-pop-gen/Models/Model1_fasta_files"
cd $DIR

# Function to process each FASTA file
process_file() {
    FILE=$1
    echo "Processing file: $FILE"
    
    # Delete any existing .arg.stats files to avoid the "file already exists" error
    if [ -f "${FILE}.arg.stats" ]; then
        echo "Deleting existing .arg.stats file for $FILE"
        rm -f ${FILE}.arg.stats
    fi
    
    # Generate the ARG for each FASTA file using arg-sample
    arg-sample --fasta $FILE --output ${FILE}.arg --verbose 0 --sample-step 100
    
    # Extract TMRCA for each ARG
    echo "Extracting TMRCA for $FILE"
    arg-extract-tmrca ${FILE}.arg.%d.smc.gz > ${FILE%.fasta}.tmrca.txt

    
    echo "TMRCA extraction completed for $FILE"
}

# Export the function to be used by parallel
export -f process_file

# Run the processing of all FASTA files in parallel, distributing the work across NCORES
/home3/mwanjiku/miniconda3/envs/argweaver_py2/bin/parallel -j $NCORES process_file ::: *.fasta

