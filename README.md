````markdown
# Ghost Population Detection Pipeline

This repository contains a Snakemake-based workflow to detect unsampled "ghost" populations in genomic datasets using demographic inference, statistical testing, and model selection.
The pipeline integrates STRUCTURE, IMa3, ARGweaver, and custom analysis scripts to evaluate signals of ghost introgression in population genomics data.


## Key Features

- **STRUCTURE-based inference** of admixture models across multiple K values
- **Likelihood Ratio Tests (LRTs)** using IMa3 for model comparison
- **Bootstrap testing** and AIC/BIC model selection
- **ARGweaver-based coalescent simulations** and TMRCA distribution analysis
- **Multimodality tests** (e.g., Hartigan’s Dip Test) to detect non-standard coalescent patterns
- Fully automated with **Snakemake**
- Reproducible environments with **Conda**


## Repository Structure

```text
ghost-pop-gen/
├── Snakefile                  # Main Snakemake pipeline
├── Snakefile.part3            # ARGweaver + modality analysis
├── environment.yml            # Main conda environment
├── config/                    # YAML config files and model specifications
│   ├── config.yaml
│   ├── model1.par
│   └── nested_models_2pop.txt
├── data/                      # Input files (FASTA, .u, .str)
│   ├── fasta/
│   ├── ima3_inputs_2pop/
│   ├── ima3_inputs_3pop/
│   └── structure_inputs/
├── envs/                      # Conda envs for specific tools
│   └── argweaver_py2.yaml
├── results/                   # Output files and visualizations
│   ├── structure_outputs/
│   ├── ima3/
│   ├── *.csv
│   ├── *.png
├── scripts/                   # R, Python, and Bash helper scripts
├── software/                  # Compiled tools (e.g., ARGweaver)
└── README.md                 
````


## Installation and Setup

1. Clone the repository:

```bash
git clone https://github.com/Megmugure/ghost-pop-gen.git
cd ghost-pop-gen
```

2. Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate ghost-pop-gen
```


## Running the Pipeline

To run the full workflow:

```bash
snakemake --cores 4
```

To perform a dry run:

```bash
snakemake -n
```

To run the ARGweaver + modality testing component separately:

```bash
snakemake -s Snakefile.part3 --cores 4
```

To generate a DAG (workflow graph):

```bash
snakemake --dag | dot -Tpng > dag.png
```

## Example Analysis Commands

```bash
# Run Kolmogorov-Smirnov test on TMRCA values
Rscript scripts/KS_tests.R data/input.tmrca

# Run LRT test
python scripts/LRT_test.py results/lrt_values.txt

# Run STRUCTURE bootstrap LRT
bash scripts/bootstrap_test.sh data/structure_data.str
```

## Citing This Work

If you use this pipeline in your research, please cite:

> (preprint link or DOI coming soon)

## License

This project is licensed under the MIT License.
See the [LICENSE](LICENSE) file for full details.


## Author

**Margaret Wanjiku**
[margaretwmugure@gmail.com](mailto:margaretwmugure@gmail.com)
[GitHub: Megmugure](https://github.com/Megmugure)
