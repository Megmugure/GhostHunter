import os
import argparse
import subprocess

def build_mainparams(path, num_inds, num_loci, burnin, numreps):
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

def main():
    parser = argparse.ArgumentParser(description="Run STRUCTURE for a given replicate and K.")
    parser.add_argument("--input", required=True, help="Input .str file")
    parser.add_argument("--outdir", required=True, help="Output directory for this STRUCTURE run")
    parser.add_argument("--structure_exec", required=True, help="Path to STRUCTURE executable")
    parser.add_argument("--K", type=int, required=True, help="Number of clusters (K)")
    parser.add_argument("--burnin", type=int, default=100000)
    parser.add_argument("--numreps", type=int, default=500000)

    args = parser.parse_args()

    input_file = args.input
    K = args.K
    outdir = args.outdir
    structure_exec = args.structure_exec
    burnin = args.burnin
    numreps = args.numreps

    os.makedirs(outdir, exist_ok=True)

    with open(input_file) as f:
        lines = [l.strip() for l in f if l.strip()]
        num_inds = len(lines)
        num_loci = len(lines[0].split()) - 1

    mainparams = os.path.join(outdir, "mainparams.txt")
    extraparams = os.path.join(outdir, "extraparams.txt")
    output_prefix = os.path.join(outdir, f"structure_run_K{K}")  # STRUCTURE appends _f

    build_mainparams(mainparams, num_inds, num_loci, burnin, numreps)
    build_extraparams(extraparams)

    cmd = [
        structure_exec,
        "-K", str(K),
        "-m", mainparams,
        "-e", extraparams,
        "-i", input_file,
        "-o", output_prefix
    ]

    print("🔧 Running STRUCTURE:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"STRUCTURE run complete: {output_prefix}_f")

if __name__ == "__main__":
    main()
