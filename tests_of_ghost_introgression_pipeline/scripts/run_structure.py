#!/usr/bin/env python3
"""
Script: run_structure.py
Author: Margaret Wanjiku
Purpose:
    Wraps the STRUCTURE executable for a specific replicate and value of K.
    Generates parameter files and invokes STRUCTURE with appropriate flags.

Inputs:
    --input: .str file
    --outdir: output directory
    --structure_exec: path to STRUCTURE binary
    --K: number of clusters
    --burnin: number of burn-in iterations
    --numreps: number of MCMC reps after burn-in

Output:
    STRUCTURE output file (structure_run_K{K}_f)
"""

import os
import argparse
import subprocess

# Functions to build STRUCTURE parameter files
def build_mainparams(path, num_inds, num_loci, burnin, numreps):
    """Write mainparams file with STRUCTURE parameters"""
    with open(path, "w") as f:
        f.write(f"""#define MAXPOPS 10
#define NUMINDS {num_inds}
#define NUMLOCI {num_loci}
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
#define BURNIN {burnin}
#define NUMREPS {numreps}
""")


def build_extraparams(path):
    """Write extraparams file for STRUCTURE"""
    with open(path, "w") as f:
        f.write("""INFERALPHA 1
ALPHA 1.0
POPALPHAS 0
ALPHAMAX 10.0
FREQSCORR 1
ONEFST 0
USEPOPINFO 0
LOCPRIOR 0
PRINTLIKES 0
PRINTKLD 0
""")


# CLI interface
def main():
    parser = argparse.ArgumentParser(description="Run STRUCTURE for a given replicate and K.")
    parser.add_argument("--input", required=True, help="Input .str file")
    parser.add_argument("--outdir", required=True, help="Directory to save STRUCTURE output")
    parser.add_argument("--structure_exec", required=True, help="Path to STRUCTURE executable")
    parser.add_argument("--K", type=int, required=True, help="Number of clusters (K)")
    parser.add_argument("--burnin", type=int, default=100000, help="Number of burn-in iterations")
    parser.add_argument("--numreps", type=int, default=500000, help="Number of MCMC reps after burn-in")
    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.outdir, exist_ok=True)

    # Determine number of individuals and loci from input .str file
    with open(args.input) as f:
        lines = [l.strip() for l in f if l.strip()]
        num_inds = len(lines)
        num_loci = len(lines[0].split()) - 1  # First column is individual ID

    # Create parameter files
    mainparams = os.path.join(args.outdir, "mainparams.txt")
    extraparams = os.path.join(args.outdir, "extraparams.txt")
    build_mainparams(mainparams, num_inds, num_loci, args.burnin, args.numreps)
    build_extraparams(extraparams)

    # STRUCTURE output file prefix (STRUCTURE appends _f to it)
    output_prefix = os.path.join(args.outdir, f"structure_run_K{args.K}")

    # Build STRUCTURE command
    cmd = [
        args.structure_exec,
        "-K", str(args.K),
        "-m", mainparams,
        "-e", extraparams,
        "-i", args.input,
        "-o", output_prefix
    ]

    print("Running STRUCTURE:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"STRUCTURE run complete: {output_prefix}_f")

if __name__ == "__main__":
    main()
