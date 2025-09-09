#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# ========= Config =========
PLINK_PREFIX="${PLINK_PREFIX:-results/structure_ppp_sampled/ppp_sampled.thin}"
OUTROOT="${OUTROOT:-results/admixture_bootstrap_grid}"
K_MIN="${K_MIN:-2}"
K_MAX="${K_MAX:-3}"          # inclusive upper bound; loop runs K0=K_MIN..K_MAX-1
B="${B:-100}"                # bootstrap replicates
CORES="${CORES:-8}"
SEED="${SEED:-12345}"
PYTHON="${PYTHON:-$(command -v python3 || command -v python || true)}"

# ========= Checks =========
need(){ command -v "$1" >/dev/null 2>&1 || { echo "Missing: $1" >&2; exit 1; }; }
need plink
need admixture
[[ -n "${PYTHON:-}" ]] || { echo "Need python3/python" >&2; exit 1; }
# numpy check (fail fast)
"$PYTHON" - <<'PY' >/dev/null 2>&1 || { echo "Python lacks NumPy (set PYTHON=...)" >&2; exit 1; }
import numpy as np
PY

# ========= Paths =========
mkdir -p "$OUTROOT"
OUTROOT="$(cd "$OUTROOT" && pwd)"
PLINK_DIR="$(cd "$(dirname "$PLINK_PREFIX")" && pwd)"
PLINK_BASE="$(basename "$PLINK_PREFIX")"
[[ -f "${PLINK_DIR}/${PLINK_BASE}.bed" && -f "${PLINK_DIR}/${PLINK_BASE}.bim" && -f "${PLINK_DIR}/${PLINK_BASE}.fam" ]] \
  || { echo "Missing PLINK files: ${PLINK_PREFIX}.{bed,bim,fam}" >&2; exit 1; }

LOG="$OUTROOT/run.log"

# ========= Lock (avoid duplicate runs) =========
LOCK="$OUTROOT/.bootstrap.lock"
exec 9>"$LOCK"
if ! flock -n 9; then
  echo "[$(date)] [lock] Another bootstrap run is already in progress: $LOCK" >> "$LOG"
  exit 1
fi
trap 'rm -f "$LOCK"' EXIT

echo "[$(date)] Bootstrap grid: K_MIN=$K_MIN K_MAX=$K_MAX B=$B CORES=$CORES" >> "$LOG"

# ========= Helper to parse log-likelihood =========
get_ll(){ awk '/Loglikelihood/ {ll=$2} END{print ll+0}' "$1"; }

# ========= Simulator (parametric under K) =========
# NOTE: now accepts an optional SEED as argv[6] so each replicate is unique.
SIMPY="$OUTROOT/simulate_from_QP.py"
cat > "$SIMPY" <<'PY'
from __future__ import print_function
import sys, numpy as np

if len(sys.argv) < 6:
    sys.stderr.write("usage: simulate_from_QP.py fam bim Q P out [seed]\n")
    sys.exit(2)

fam, bim, Qf, Pf, out = sys.argv[1:6]
seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0

Q = np.loadtxt(Qf)      # I x K
P = np.loadtxt(Pf)      # K x L  (ADMIXTURE .P)
# allow L x K
if P.shape[0] != Q.shape[1]:
    if P.shape[1] == Q.shape[1]:
        P = P.T
    else:
        sys.stderr.write("P/Q dimension mismatch\n"); sys.exit(2)
I, K = Q.shape
Kp, L = P.shape
assert Kp == K

# load fam/bim
ids = [ln.strip().split()[:2] for ln in open(fam)]
bimrows = [ln.strip().split() for ln in open(bim)]
A1 = [r[4] for r in bimrows]  # from .bim
A2 = [r[5] for r in bimrows]

# p_il = sum_k q_ik * p_lk
p = np.dot(Q, P)  # I x L
rng = np.random.RandomState(seed)
G = rng.binomial(2, np.clip(p, 1e-9, 1-1e-9))  # diploid genotypes

# write MAP/PED
with open(out + ".map","w") as m:
    for r in bimrows:
        m.write("{}\t{}\t{}\t{}\n".format(r[0], r[1], r[2], r[3]))
with open(out + ".ped","w") as ped:
    for i,(fid,iid) in enumerate(ids):
        fields = [fid, iid, "0","0","0","-9"]
        for l in range(L):
            g = G[i,l]
            if g==0: fields.extend([A2[l], A2[l]])
            elif g==1: fields.extend([A1[l], A2[l]])
            else: fields.extend([A1[l], A1[l]])
        ped.write(" ".join(fields)+"\n")
PY

# ========= Master summary =========
SUMMARY="$OUTROOT/summary.tsv"
echo -e "K0\tK1\tTobs\tLL_obs_K0\tLL_obs_K1\tB\tp_value\tdir" > "$SUMMARY"

# ========= Main loop over K0 =========
for K0 in $(seq "$K_MIN" $((K_MAX-1))); do
  K1=$((K0+1))
  DIR="$OUTROOT/K${K0}_vs_K${K1}"
  mkdir -p "$DIR"
  echo "[$(date)] Observed fits: K=${K0}, ${K1}" >> "$LOG"

  # Run ADMIXTURE in the PLINK directory, write logs to DIR
  pushd "$PLINK_DIR" >/dev/null
  admixture -j"$CORES" --seed "$SEED" "${PLINK_BASE}.bed" "$K0" > "$DIR/obs.K${K0}.log" 2>&1
  mv "${PLINK_BASE}.${K0}.P" "$DIR/obs.K${K0}.P"; mv "${PLINK_BASE}.${K0}.Q" "$DIR/obs.K${K0}.Q"
  admixture -j"$CORES" --seed "$SEED" "${PLINK_BASE}.bed" "$K1" > "$DIR/obs.K${K1}.log" 2>&1
  mv "${PLINK_BASE}.${K1}.P" "$DIR/obs.K${K1}.P"; mv "${PLINK_BASE}.${K1}.Q" "$DIR/obs.K${K1}.Q"
  popd >/dev/null

  LL0=$(get_ll "$DIR/obs.K${K0}.log"); LL1=$(get_ll "$DIR/obs.K${K1}.log")
  Tobs=$(python - <<PY
ll0=$LL0; ll1=$LL1
print(-2.0*(ll0-ll1))
PY
)
  echo "  [obs] K${K0}: $LL0  K${K1}: $LL1  => Tobs=$Tobs" >> "$LOG"
  echo "$Tobs" > "$DIR/Tobs.txt"

  TSV="$DIR/bootstrap.tsv"
  echo -e "rep\tLL_K0\tLL_K1\tT" > "$TSV"

  for rep in $(seq 1 "$B"); do
    RDIR="$DIR/rep$(printf "%04d" "$rep")"; mkdir -p "$RDIR"

    # simulate under K0 using observed P/Q with a UNIQUE seed each rep
    SIM_SEED=$((SEED + rep))
    "$PYTHON" "$SIMPY" "${PLINK_DIR}/${PLINK_BASE}.fam" "${PLINK_DIR}/${PLINK_BASE}.bim" \
      "$DIR/obs.K${K0}.Q" "$DIR/obs.K${K0}.P" "$RDIR/sim" "$SIM_SEED"

    # PLINK -> bed
    plink --noweb --file "$RDIR/sim" --make-bed --out "$RDIR/sim" >/dev/null 2>&1

    # fit admixture at K0 and K1 on the simulated dataset (also seed)
    admixture -j"$CORES" --seed "$((SEED + rep))" "$RDIR/sim.bed" "$K0" > "$RDIR/K${K0}.log" 2>&1
    admixture -j"$CORES" --seed "$((SEED + rep))" "$RDIR/sim.bed" "$K1" > "$RDIR/K${K1}.log" 2>&1

    ll0=$(get_ll "$RDIR/K${K0}.log"); ll1=$(get_ll "$RDIR/K${K1}.log")
    Tb=$(python - <<PY
ll0=$ll0; ll1=$ll1
print(-2.0*(ll0-ll1))
PY
)
    echo -e "${rep}\t${ll0}\t${ll1}\t${Tb}" >> "$TSV"
    echo "    [rep $rep/$B] T=$Tb" >> "$LOG"
  done

  pval=$(awk -v t="$Tobs" 'NR>1{n++; if($4>=t) ge++} END{if(n>0) printf("%.6f", ge/n); else print "NA"}' "$TSV")
  echo "$pval" > "$DIR/pvalue.txt"
  echo -e "${K0}\t${K1}\t${Tobs}\t${LL0}\t${LL1}\t${B}\t${pval}\t${DIR}" >> "$SUMMARY"
  echo "  [done] K${K0}->K${K1}: p=$pval" >> "$LOG"
done

echo "[$(date)] [summary] $SUMMARY" >> "$LOG"
