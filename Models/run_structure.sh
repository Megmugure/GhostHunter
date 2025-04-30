#!/bin/bash
#PBS -V
#PBS -l nodes=1:ppn=8
#PBS -N STRUCTURE_Run
#PBS -joe
#PBS -q batch
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR
NCORES=$(wc -w < $PBS_NODEFILE)
echo "Running STRUCTURE on all models using $NCORES cores"

# Path to STRUCTURE binary (adjust if needed)
STRUCTURE_EXEC=./structure

# STRUCTURE parameters
BURNIN=100000
NUMREPS=500000
OUTPUT_BASE="structure_outputs"
mkdir -p "$OUTPUT_BASE"

# Range of K values to test
K_VALUES=(1 2 3 4 5 6)

# Loop over each model directory
for model_dir in Model1 Model2 Model3 Model4 Model5; do
    echo "Processing $model_dir..."

    # Loop over each replicate file in the model directory
    for replicate_file in ${model_dir}/*_cleaned.str; do
        filename=$(basename "$replicate_file")
        model_name=$(echo "$filename" | cut -d'_' -f1)
        replicate_num=$(echo "$filename" | grep -oP '(?<=replicate)\d+')

        # Get the number of loci and individuals from the replicate file
        num_fields=$(awk 'NR==1 {print NF}' "$replicate_file")
        num_loci=$((num_fields - 1))
        num_individuals=$(wc -l < "$replicate_file")

        # Create output directory for this replicate
        output_dir="${OUTPUT_BASE}/${model_name}/replicate${replicate_num}"
        mkdir -p "$output_dir"

        # Create the mainparams.txt file for the replicate
        mainparams_file="${output_dir}/mainparams.txt"
        cat > "$mainparams_file" <<EOL
#define MAXPOPS 10
#define NUMINDS ${num_individuals}
#define NUMLOCI ${num_loci}
#define PLOIDY 1
#define MISSING 0
#define ONEROWPERIND 1
#define LABEL 1
#define POPDATA 0
#define POPFLAG 0
#define LOCDATA 0
#define PHENOTYPE 0
#define EXTRACOLS 0
#define MARKERNAMES 0
#define RECESSIVEALLELES 0
#define MAPDISTANCES 0
#define PHASED 0
#define PHASEINFO 0
#define PHASEDINPUT 0
#define USEPHASEINFORECORDS 0
#define MARKOVPHASE 0
#define FASTPHASE 0
#define BURNIN ${BURNIN}
#define NUMREPS ${NUMREPS}
EOL

        # Create the extraparams.txt file for the replicate
        extraparams_file="${output_dir}/extraparams.txt"
        cat > "$extraparams_file" <<EOL
INFERALPHA 1
ALPHA 1.0
POPALPHAS 0
ALPHAMAX 10.0
FREQSCORR 1
ONEFST 0
USEPOPINFO 0
LOCPRIOR 0
PRINTLIKES 0
PRINTKLD 0
EOL

        # Loop over each K value
        for K in "${K_VALUES[@]}"; do
            output_prefix="${output_dir}/structure_run_K${K}"
            echo "Running STRUCTURE for ${filename} with K=${K}..."
            ${STRUCTURE_EXEC} -K ${K} -m "$mainparams_file" -e "$extraparams_file" -i "$replicate_file" -o "$output_prefix"
        done
    done
done

echo "STRUCTURE runs completed."
