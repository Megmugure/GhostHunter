########################################
# CONFIG
########################################

configfile: "config/config.yaml"

import os
from glob import glob

########################################
# STRUCTURE SETTINGS
########################################

STRUCTURE_DIR = config["input_dirs"]["structure"]
KS = config.get("ks", [1,2,3,4,5,6])

STRUCTURE_FILES = glob(os.path.join(STRUCTURE_DIR, "*.str"))
STRUCTURE_SAMPLES = [
    os.path.splitext(os.path.basename(f))[0]
    for f in STRUCTURE_FILES
]

########################################
# IMA3 SETTINGS
########################################

IMA3_DIR = config["input_dirs"]["ima3"]

IMA3_FILES = glob(os.path.join(IMA3_DIR, "*.u"))
IMA3_SAMPLES = [
    os.path.splitext(os.path.basename(f))[0]
    for f in IMA3_FILES
]

########################################
# ARGWEAVER SETTINGS
########################################

ARG_CFG = config.get("argweaver", {})
RUN_ARGWEAVER = ARG_CFG.get("run", False)

FASTA_DIR = ARG_CFG.get("input_dir", "data/argweaver_inputs")

# If user explicitly provides samples, use them
if "samples" in ARG_CFG and ARG_CFG["samples"]:
    FASTA_SAMPLES = ARG_CFG["samples"]
else:
    FASTA_FILES = glob(os.path.join(FASTA_DIR, "*.fasta"))
    FASTA_SAMPLES = [
        os.path.splitext(os.path.basename(f))[0]
        for f in FASTA_FILES
    ] if RUN_ARGWEAVER else []

########################################
# FINAL TARGET
########################################

rule all:
    input:
        expand(
            "results/structure/{sample}/STRUCTURE_REPORT.txt",
            sample=STRUCTURE_SAMPLES
        ),
        expand(
            "results/ima3/{sample}/ima3_summary.txt",
            sample=IMA3_SAMPLES
        )
        +
        (
            expand(
                "results/argweaver/{sample}/{sample}.stats.txt",
                sample=FASTA_SAMPLES
            ) if RUN_ARGWEAVER else []
        )
        +
        expand(
    "results/argweaver/{sample}/{sample}.ghost_summary.txt",
    sample=FASTA_SAMPLES
        )

########################################
# STRUCTURE
########################################

rule run_structure:
    input:
        strfile=os.path.join(STRUCTURE_DIR, "{sample}.str")

    output:
        "results/structure/{sample}/K{K}/structure_run_K{K}_f"

    params:
        structure_exec=config["structure_exec"],
        outdir="results/structure/{sample}/K{K}"

    shell:
        """
        set -euo pipefail

        mkdir -p {params.outdir}

        python scripts/STRUCTURE_analysis.py \
            --input {input.strfile} \
            --outdir {params.outdir} \
            --structure_exec {params.structure_exec} \
            --K {wildcards.K}
        """


rule parse_lnprob:
    input:
        expand(
            "results/structure/{{sample}}/K{K}/structure_run_K{K}_f",
            K=KS
        )

    output:
        "results/structure/{sample}/lnprob.tsv"

    shell:
        """
        python scripts/parse_lnprob.py \
            --indir results/structure/{wildcards.sample} \
            --out {output}
        """


rule model_selection:
    input:
        lnprob="results/structure/{sample}/lnprob.tsv",
        strfile=os.path.join(STRUCTURE_DIR, "{sample}.str")

    output:
        evanno="results/structure/{sample}/summary/evanno.tsv",
        aicbic="results/structure/{sample}/summary/aic_bic.tsv"

    shell:
        """
        mkdir -p results/structure/{wildcards.sample}/summary

        I=$(wc -l < {input.strfile})
        L=$(awk '{{print NF-1; exit}}' {input.strfile})

        python scripts/compute_model_selection.py \
            --input {input.lnprob} \
            --outdir results/structure/{wildcards.sample}/summary \
            --I $I \
            --L $L
        """


rule infer_k:
    input:
        evanno="results/structure/{sample}/summary/evanno.tsv",
        aicbic="results/structure/{sample}/summary/aic_bic.tsv"

    output:
        "results/structure/{sample}/summary/bestK.txt"

    shell:
        """
        python scripts/infer_best_k.py \
            --evanno {input.evanno} \
            --aicbic {input.aicbic} \
            --out {output}
        """


rule report:
    input:
        evanno="results/structure/{sample}/summary/evanno.tsv",
        aicbic="results/structure/{sample}/summary/aic_bic.tsv",
        bestk="results/structure/{sample}/summary/bestK.txt"

    output:
        "results/structure/{sample}/STRUCTURE_REPORT.txt"

    shell:
        """
        python scripts/render_report.py \
            --evanno {input.evanno} \
            --aicbic {input.aicbic} \
            --bestk {input.bestk} \
            --out {output}
        """

########################################
# IMA3
########################################

rule run_ima3:
    input:
        ufile=os.path.join(IMA3_DIR, "{sample}.u")

    output:
        "results/ima3/{sample}/ima3.out"

    params:
        exec=config["ima3_exec"],
        burnin=config["ima3"]["burnin"],
        chain=config["ima3"]["chain_length"],
        migration=config["ima3"]["migration"],
        theta=config["ima3"]["theta"]

    shell:
        """
        set -euo pipefail

        mkdir -p results/ima3/{wildcards.sample}

        {params.exec} \
            -i {input.ufile} \
            -o {output} \
            -q100 \
            -m{params.migration} \
            -t{params.theta} \
            -b{params.burnin} \
            -L{params.chain} \
            -d200 \
            -p 2 \
            -r245 \
            -hn 20 \
            -ha 0.99 \
            -hb 0.9
        """


rule parse_ima3:
    input:
        "results/ima3/{sample}/ima3.out"

    output:
        "results/ima3/{sample}/ima3_summary.txt"

    shell:
        """
        python scripts/parse_ima3.py \
            --input {input} \
            --out {output}
        """

########################################
# ARGWEAVER
########################################

rule run_argweaver:
    input:
        fasta=lambda wc: os.path.join(
            FASTA_DIR,
            f"{wc.sample}.fasta"
        )

    output:
        smc="results/argweaver/{sample}/{sample}.arg.0.smc.gz"

    params:
        exe=ARG_CFG["executable"],
        popsize=ARG_CFG["popsize"],
        mutrate=ARG_CFG["mutrate"],
        recombrate=ARG_CFG["recombrate"],
        ntimes=ARG_CFG["ntimes"],
        maxtime=ARG_CFG["maxtime"],
        iters=ARG_CFG["iters"],
        step=ARG_CFG["sample_step"]

    shell:
        r"""
        mkdir -p results/argweaver/{wildcards.sample}

        {params.exe} \
            --fasta {input.fasta} \
            --output results/argweaver/{wildcards.sample}/{wildcards.sample}.arg \
            --popsize {params.popsize} \
            --mutrate {params.mutrate} \
            --recombrate {params.recombrate} \
            --ntimes {params.ntimes} \
            --maxtime {params.maxtime} \
            --iters {params.iters} \
            --sample-step {params.step} \
            --overwrite \
            --verbose 1
        """

rule extract_tmrca:
    input:
        "results/argweaver/{sample}/{sample}.arg.0.smc.gz"

    output:
        "results/argweaver/{sample}/{sample}.tmrca.tsv"

    params:
        extract=ARG_CFG["extract_exec"]

    shell:
        r"""
        {params.extract} \
            results/argweaver/{wildcards.sample}/{wildcards.sample}.arg.*.smc.gz \
            > {output}
        """

rule analyze_tmrca:
    input:
        tsv="results/argweaver/{sample}/{sample}.tmrca.tsv"

    output:
        csv="results/argweaver/{sample}/{sample}.summary.csv",
        png="results/argweaver/{sample}/{sample}.hist.png",
        txt="results/argweaver/{sample}/{sample}.ghost_summary.txt"

    shell:
        r"""
        set -euo pipefail

        /home3/mwanjiku/Tools/miniconda3/bin/conda run -n ghost-pop-gen-r \
            Rscript scripts/analyze_tmrca.R \
            {input.tsv} \
            {output.csv} \
            {output.png} \
            {output.txt}
        """