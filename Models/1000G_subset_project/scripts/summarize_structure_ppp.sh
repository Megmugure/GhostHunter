#!/usr/bin/env bash
# ------------------------------------------------------------------------------------
# Purpose:
#   Parse STRUCTURE outputs to summarize lnP(D) by K, compute Evanno ΔK and bestK,
#   compute AIC/BIC per K using p = I*(K-1) + K*L, pick best replicate per K,
#   and snapshot key counts. Produces tidy TSVs + small text files for downstream plots.
#
# Usage:
#   BASE=results/structure_ppp_sampled ./scripts/summarize_structure.sh
#
# Inputs:
#   - ${BASE}/K*/rep*/structure.stdout        # captured logs with "Estimated Ln Prob..."
#   - ${BASE}/ppp_sampled.thin.{bim,fam}      # for I, L counts
#   - data/panel/integrated_call_samples_v3.20130502.ALL.panel (optional CEU/CHS counts)
#   - data/subsets_fst_sampled/CEU_CHS.sampled.windows.vcf.gz (optional mapping to panel)
#
# Outputs (to ${BASE}/summary/):
#   - summary_lnprob.tsv
#   - evanno.tsv, bestK.txt
#   - aic_bic.tsv
#   - best_rep_per_K.tsv
#   - counts.txt
#   - params_snapshot.txt
#
# Dependencies:
#   - bash coreutils, awk, grep, find
#   - python (2.7+ or 3.x) for the Evanno/AIC/BIC block
#   - bcftools (optional) if you want CEU/CHS counts from VCF samples
# ------------------------------------------------------------------------------------

set -euo pipefail
export LC_ALL=C

# ---------- locations (match your PBS outputs) ----------
BASE="${BASE:-results/structure_ppp_sampled}"
OUT="${BASE}/summary"
PANEL="${PANEL:-data/panel/integrated_call_samples_v3.20130502.ALL.panel}"
mkdir -p "$OUT"

# Inputs produced by STRUCTURE_PPP.pbs (PLINK stage + STRUCTURE input)
BED_PREFIX="${BASE}/ppp_sampled.thin"       # .fam/.bim
STR_IN="${BASE}/ppp_sampled.thin_structure.recode.strct_in"
MAINP="${BASE}/mainparams.fix"
EXTRAP="${BASE}/extraparams.fix"

# ---------- sanity ----------
[[ -f "$BED_PREFIX.fam" && -f "$BED_PREFIX.bim" ]] || { echo "Missing PLINK outputs under $BASE"; exit 1; }
[[ -f "$STR_IN" ]] || { echo "Missing STRUCTURE .str input at $STR_IN"; exit 1; }

# ---------- counts ----------
I=$(wc -l < "$BED_PREFIX.fam")   # number of individuals
L=$(wc -l < "$BED_PREFIX.bim")   # number of loci (biallelic SNPs)

# Optional: CEU/CHS counts for Methods (using sampled VCF if present)
VCF_PATH="data/subsets_fst_sampled/CEU_CHS.sampled.windows.vcf.gz"
VCF_SAMPLES="$OUT/samples.in.vcf.txt"
CEU_CHS_MAP="$OUT/CEU_CHS.map.tsv"
if command -v bcftools >/dev/null 2>&1 && [[ -f "$VCF_PATH" ]]; then
  bcftools query -l "$VCF_PATH" > "$VCF_SAMPLES"
  if [[ -s "$PANEL" ]]; then
    awk '$2=="CEU" || $2=="CHS"{print $1"\t"$2}' "$PANEL" \
      | awk 'NR==FNR{a[$1]=$2; next} ($1 in a){print $1"\t"a[$1]}' - "$VCF_SAMPLES" \
      > "$CEU_CHS_MAP" || true
  fi
fi
N_CEU=$(awk '$2=="CEU"{c++} END{print c+0}' "$CEU_CHS_MAP" 2>/dev/null || echo 0)
N_CHS=$(awk '$2=="CHS"{c++} END{print c+0}' "$CEU_CHS_MAP" 2>/dev/null || echo 0)

{
  echo -e "metric\tvalue"
  echo -e "N_individuals\t$I"
  echo -e "N_loci\t$L"
  echo -e "N_CEU\t$N_CEU"
  echo -e "N_CHS\t$N_CHS"
} > "$OUT/counts.txt"

# ---------- collect lnProb across all runs ----------
SUMMARY="$OUT/summary_lnprob.tsv"
echo -e "K\trep\tlnprob\tstdout_path" > "$SUMMARY"

while IFS= read -r -d '' f; do
  K=$(echo "$f" | sed -n 's#.*/K\([0-9]\+\)/rep\([0-9]\+\)/structure.stdout$#\1#p')
  REP=$(echo "$f" | sed -n 's#.*/K\([0-9]\+\)/rep\([0-9]\+\)/structure.stdout$#\2#p')
  LP=$(grep -m1 -E "Estimated Ln Prob of Data" "$f" | awk '{print $NF}')
  [[ -z "$LP" ]] && LP="NA"
  echo -e "${K}\t${REP}\t${LP}\t${f}" >> "$SUMMARY"
done < <(find "$BASE" -type f -name structure.stdout -print0)

# ---------- choose a Python interpreter (prefer python3) ----------
PYTHON="${PYTHON:-$(command -v python3 || command -v python || true)}"
if [[ -z "${PYTHON}" ]]; then
  echo "Need python (2.7+ or 3.x). Set PYTHON=/path/to/python" >&2
  exit 1
fi

# ---------- Evanno ΔK, bestK, best rep per K, AIC/BIC ----------
"$PYTHON" - "$SUMMARY" "$OUT" "$I" "$L" <<'PY'
# This inline Python block:
#  - Reads summary_lnprob.tsv
#  - Computes mean/sd lnP(D) by K
#  - Computes L'(K), |L''(K)|, and Evanno ΔK = |L''(K)| / SD(lnP(D)_K)
#  - Writes evanno.tsv and bestK.txt
#  - Picks best replicate per K (max lnP(D))
#  - Computes AIC/BIC per K using p = I*(K-1) + K*L (I=indivs, L=loci)
from __future__ import print_function
import sys, csv, math, collections, os

def to_float(x):
    try: return float(x)
    except: return None

def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs)/float(len(xs)) if xs else float('nan')

def pstdev(xs):
    xs = [x for x in xs if x is not None]
    n = len(xs)
    if n <= 1:
        return float('nan')
    mu = mean(xs)
    var = sum((x - mu) ** 2 for x in xs) / float(n)
    return math.sqrt(var)

summary, outdir, I, L = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

byK = collections.defaultdict(list)
with open(summary) as fh:
    rdr = csv.DictReader(fh, delimiter='\t')
    for r in rdr:
        try:
            k = int(r['K']); rep = int(r['rep'])
            lp = to_float(r['lnprob'])
        except:
            continue
        if lp is not None:
            byK[k].append((rep, lp, r['stdout_path']))

Ks = sorted(byK.keys())
means = dict((K, mean([lp for (_, lp, _) in byK[K]])) for K in Ks)
sds   = dict((K, pstdev([lp for (_, lp, _) in byK[K]])) for K in Ks)

# L'(K) and |L''(K)|
L1 = {}
L2 = {}
for i, K in enumerate(Ks):
    if i > 0:
        prevK = Ks[i-1]
        L1[K] = means[K] - means[prevK]
for i, K in enumerate(Ks):
    if 0 < i < len(Ks)-1:
        L2[K] = abs(means[Ks[i+1]] - 2.0*means[K] + means[Ks[i-1]])

DeltaK = {}
for i, K in enumerate(Ks):
    if 0 < i < len(Ks)-1:
        sd = sds[K]
        if sd == sd and sd > 0:  # not NaN and >0
            DeltaK[K] = L2[K] / sd
        else:
            DeltaK[K] = float('nan')

# Write evanno table
with open(os.path.join(outdir, "evanno.tsv"), "w") as fw:
    fw.write("K\tmean_lnprob\tsd_lnprob\tL1\tabs_L2\tDeltaK\n")
    for K in Ks:
        fw.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
            K, means.get(K, 'NA'), sds.get(K, 'NA'),
            L1.get(K, 'NA'), L2.get(K, 'NA'), DeltaK.get(K, 'NA')
        ))

# Best K = argmax DeltaK (Evanno)
bestK = None
bestVal = None
for K, val in DeltaK.items():
    if val == val:  # not NaN
        if bestVal is None or val > bestVal:
            bestVal, bestK = val, K
with open(os.path.join(outdir, "bestK.txt"), "w") as fw:
    fw.write(str(bestK) if bestK is not None else "NA")

# Best replicate per K (highest lnprob)
with open(os.path.join(outdir, "best_rep_per_K.tsv"), "w") as fw:
    fw.write("K\tbest_rep\tbest_lnprob\tstdout_path\n")
    for K in Ks:
        rep, lp, path = max(byK[K], key=lambda t: t[1])
        fw.write("{}\t{}\t{}\t{}\n".format(K, rep, lp, path))

# AIC/BIC per K using manuscript p = I*(K-1) + K*L  (L≈#biallelic loci so A≈L)
def aic(ll, p): return -2.0*ll + 2.0*p
def bic(ll, p, n): return -2.0*ll + p*math.log(float(n))
with open(os.path.join(outdir, "aic_bic.tsv"), "w") as fw:
    fw.write("K\tmean_lnprob\tp_params\tAIC\tBIC\n")
    for K in Ks:
        p = I*(K-1) + K*L
        ll = means[K]
        fw.write("{}\t{}\t{}\t{}\t{}\n".format(K, ll, p, aic(ll,p), bic(ll,p,I)))
PY

# ---------- snapshot params (for Methods) ----------
{
  echo "[files]"
  [[ -f "$MAINP" ]] && echo "mainparams: $MAINP"
  [[ -f "$EXTRAP" ]] && echo "extraparams: $EXTRAP"
  echo "[counts]"
  echo "N_individuals (I): $I"
  echo "N_loci (L): $L"
} > "$OUT/params_snapshot.txt"

echo "[ok] Wrote:"
echo "  - $SUMMARY"
echo "  - $OUT/evanno.tsv  (bestK: $(cat "$OUT/bestK.txt" 2>/dev/null || echo NA))"
echo "  - $OUT/aic_bic.tsv"
echo "  - $OUT/best_rep_per_K.tsv"
echo "  - $OUT/counts.txt"
echo "  - $OUT/params_snapshot.txt"
