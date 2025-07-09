#!/usr/bin/env python3

"""
bootstrap_structure_lrt_all.py

This script performs a parametric bootstrap likelihood ratio test (LRT) on STRUCTURE output
to assess model support for different values of K (number of ancestral populations).

Steps:
1. Parses STRUCTURE output files to extract Q matrices and allele frequencies.
2. Simulates genotypes under the current K model.
3. Re-runs ADMIXTURE on simulated data at K and K+1.
4. Computes the test statistic T_obs and compares it to a bootstrapped null distribution.
5. Outputs per-model, per-replicate test results into CSV.

Author: [Your Name or Institution]
Date: [Optional]
"""

import os
import re
import numpy as np
import pandas as pd
from glob import glob
import subprocess

# STRUCTURE Parsing

def parse_structure_output(filepath):
    """
    Parses STRUCTURE output file to extract:
    - Q matrix: ancestry proportions (individuals x clusters)
    - Allele frequencies per locus per cluster
    """
    q_matrix = []
    freqs = {}
    reading_q = False
    reading_freqs = False
    current_locus = None
    locus_counter = 1
    cluster_id = 0  # STRUCTURE outputs one cluster block at a time

    with open(filepath, 'r') as file:
        lines = file.readlines()

    for line in lines:
        stripped = line.strip()

        if "Inferred ancestry of individuals" in stripped:
            reading_q = True
            continue

        if reading_q:
            if not stripped or "Label" in stripped:
                continue
            if "Estimated Allele Frequencies in each cluster" in stripped:
                reading_q = False
                reading_freqs = True
                continue
            try:
                parts = stripped.split(":")[-1].strip().split()
                q_matrix.append([float(x) for x in parts])
            except ValueError:
                continue

        if reading_freqs and "First column gives estimated ancestral frequencies" in stripped:
            continue

        if reading_freqs:
            if stripped.startswith("Locus"):
                current_locus = f"L{locus_counter}"
                locus_counter += 1
                continue
            if stripped.startswith("Values of parameters used"):
                break
            parts = stripped.split()
            if len(parts) >= 3 and parts[0].isdigit():
                allele = parts[0]
                try:
                    freq = float(parts[-1])
                    freqs[(cluster_id, current_locus, allele)] = freq
                except ValueError:
                    continue

    return np.array(q_matrix), freqs

# Genotype simulation

def convert_freqs_to_array(freqs, loci_ids, alleles, K):
    """
    Converts allele frequencies from dictionary to structured NumPy array of shape (K, L, A)
    """
    L = len(loci_ids)
    A = len(alleles)
    arr = np.zeros((K, L, A))
    for k in range(K):
        for l_idx, locus in enumerate(loci_ids):
            for a_idx, a in enumerate(alleles):
                key = (k, locus, a)
                if key in freqs:
                    arr[k, l_idx, a_idx] = freqs[key]
    return arr

def simulate_genotypes(I, L, allele_freqs, q_matrix):
    """
    Simulates diploid genotypes given:
    - I individuals
    - L loci
    - allele_freqs: shape (K, L, 2)
    - q_matrix: shape (I, K)
    
    Returns: array of shape (I, L, 2)
    """
    K = allele_freqs.shape[0]
    ploidy = 2
    genos = np.zeros((I, L, ploidy), dtype=int)
    for i in range(I):
        q_probs = q_matrix[i]
        if np.sum(q_probs) == 0:
            raise ValueError(f"Invalid q_matrix for individual {i}: sum=0")
        q_probs = np.clip(q_probs, 0, 1)
        q_probs /= np.sum(q_probs)
        for l in range(L):
            for p in range(ploidy):
                k = np.random.choice(K, p=q_probs)
                p1 = allele_freqs[k, l, 1]
                allele_probs = np.clip([1 - p1, p1], 0, 1)
                allele_probs /= np.sum(allele_probs)
                allele = np.random.choice([0, 1], p=allele_probs)
                genos[i, l, p] = allele
    return genos

# PLINK/ADMIXTURE I/O

def write_ped_map(genos, loci_ids, output_prefix):
    """
    Writes PLINK .ped and .map files for ADMIXTURE input
    """
    I, L, P = genos.shape
    with open(f"{output_prefix}.ped", 'w') as ped:
        for i in range(I):
            row = ['FAM1', f'ind{i+1}', '0', '0', '1', '1']
            for l in range(L):
                alleles = [str(genos[i, l, j] + 1) for j in range(P)]
                row.extend(alleles)
            ped.write(' '.join(row) + '\n')
    with open(f"{output_prefix}.map", 'w') as mapf:
        for idx, locus in enumerate(loci_ids):
            mapf.write(f"1\t{locus}\t0\t{idx+1}\n")

def convert_to_bed(ped_prefix):
    """
    Converts PLINK .ped/.map to binary .bed format using PLINK
    """
    cmd = ["plink", "--file", ped_prefix, "--make-bed", "--out", ped_prefix]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

def run_admixture(ped_prefix, K):
    """
    Runs ADMIXTURE on binary .bed file and extracts log-likelihood from the log output.
    Ensures all output is written inside the correct directory.
    """
    convert_to_bed(ped_prefix)

    dirname = os.path.dirname(ped_prefix)
    basename = os.path.basename(ped_prefix)
    log_file = os.path.join(dirname, f"{basename}.{K}.log")
    bed_file = f"{basename}.bed"

    cmd = ["admixture", "--cv", bed_file, str(K)]
    
    with open(log_file, "w") as lf:
        subprocess.run(cmd, cwd=dirname, stdout=lf, stderr=lf, check=True)

    with open(log_file, "r") as f:
        for line in f:
            if "Loglikelihood" in line:
                match = re.search(r"Loglikelihood:.*?([-\d.]+)", line)
                if match:
                    return float(match.group(1))

    raise RuntimeError(f"Loglikelihood not found in {log_file}")


# File helpers

def find_k_files(rep_dir):
    return sorted([f for f in os.listdir(rep_dir) if re.match(r'structure_run_K\d+_f$', f)])

def extract_k(filename):
    return int(re.findall(r'K(\d+)', filename)[0])

# Main routine

def main():
    """
    Main entry point for running bootstrap LRT analysis
    """
    root = "structure_outputs"
    bootstraps = 100
    loci_sample = 200
    alleles = ['1', '2']
    results = []

    # Create clean output directory for bootstrap simulations
    output_dir = "bootstrap_temp"
    os.makedirs(output_dir, exist_ok=True)

    print("Starting parametric bootstrap analysis...\n")

    for model_path in sorted(glob(f"{root}/model*/")):
        model = os.path.basename(os.path.normpath(model_path))

        for rep_path in sorted(glob(f"{model_path}/replicate*/")):
            replicate = os.path.basename(os.path.normpath(rep_path))
            k_files = find_k_files(rep_path)
            k_vals = sorted(set(extract_k(f) for f in k_files))

            if not k_vals:
                print(f"No K files found in {rep_path}")
                continue

            for k in k_vals:
                k1 = k + 1
                k_file = os.path.join(rep_path, f'structure_run_K{k}_f')
                k1_file = os.path.join(rep_path, f'structure_run_K{k1}_f')
                if not os.path.exists(k_file) or not os.path.exists(k1_file):
                    print(f"Skipping {rep_path} K={k}→{k1}: one or both files missing")
                    continue

                try:
                    q_k, f_k = parse_structure_output(k_file)
                    q_k1, f_k1 = parse_structure_output(k1_file)
                    print(f"Parsed {rep_path} K={k}: Q={q_k.shape}, Freqs={len(f_k)}")
                except Exception as e:
                    print(f"Parsing error in {rep_path}: {e}")
                    continue

                K_actual = q_k.shape[1]
                loci_ids = sorted(set(k[1] for k in f_k.keys()))
                if len(loci_ids) < loci_sample:
                    print(f"Skipping {rep_path}: only {len(loci_ids)} usable loci (need {loci_sample})")
                    continue

                loci_sampled = loci_ids[:loci_sample]
                I = q_k.shape[0]
                f_arr_k = convert_freqs_to_array(f_k, loci_sampled, alleles, K_actual)

                try:
                    sim_prefix = os.path.join(output_dir, "bootstrap_sim")
                    genos = simulate_genotypes(I, loci_sample, f_arr_k, q_k)
                    write_ped_map(genos, loci_sampled, sim_prefix)
                    ll_k = run_admixture(sim_prefix, k)
                    ll_k1 = run_admixture(sim_prefix, k1)
                    T_obs = -2 * (ll_k - ll_k1)
                except Exception as e:
                    print(f"ADMIXTURE failed on {model}/{replicate} K={k}: {e}")
                    continue

                bootstrap_Ts = []
                for b in range(bootstraps):
                    try:
                        g_b = simulate_genotypes(I, loci_sample, f_arr_k, q_k)
                        prefix = os.path.join(output_dir, f"bootstrap_b{b}")
                        write_ped_map(g_b, loci_sampled, prefix)
                        llb_k = run_admixture(prefix, k)
                        llb_k1 = run_admixture(prefix, k1)
                        T_b = -2 * (llb_k - llb_k1)
                        bootstrap_Ts.append(T_b)
                    except Exception as e:
                        print(f"Bootstrap {b} failed for {model}/{replicate}: {e}")
                        continue

                p_val = np.mean([t > T_obs for t in bootstrap_Ts])
                print(f"[{model}/{replicate}] K={k} → K+1={k1}: T_obs={T_obs:.2f}, p={p_val:.3f}")

                results.append({
                    'model': model,
                    'replicate': replicate,
                    'K': k,
                    'K+1': k1,
                    'T_obs': T_obs,
                    'p_value': p_val
                })

    df = pd.DataFrame(results)
    df.to_csv("bootstrap_lrt_results_all.csv", index=False)
    print("\nDone. Results written to 'bootstrap_lrt_results_all.csv'")

# Entry

if __name__ == "__main__":
    main()
