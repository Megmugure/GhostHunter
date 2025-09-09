#!/usr/bin/env bash
# ------------------------------------------------------------------------------------
# Purpose:
#   Compute windowed Weir–FST (PPP / py-popgen) for CEU vs CHS on your subsetted VCFs,
#   merge & filter windows, uniformly sample windows across FST bins (e.g., low/mid/high),
#   and extract those genomic regions back into per-chr and combined VCFs.
#
# Usage:
#   WIN=50000 STEP=50000 NBINS=3 SAMPLE_TOTAL=300 SEED=12345 PPP_ENV=ppp_env \
#     ./scripts/ppp_fst_binning_sample.sh
#
# Inputs:
#   - data/subsets/CEU_CHS.chr{1..22}.vcf.gz           (from step 1)
#   - data/panel/integrated_call_samples_v3.20130502.ALL.panel
#
# Outputs:
#   - data/fst/CEU_CHS.chrN.windowed.weir.fst          # per-chr tables
#   - data/fst/CEU_CHS.allchr.windowed.weir.fst        # merged
#   - data/fst/CEU_CHS.allchr.windowed.weir.filtered.fst
#   - data/bins/CEU_CHS.uniform${NBINS}.N${SAMPLE_TOTAL}.windowed.weir.fst
#   - data/bins/regions.chrN.tsv                       # sampled windows per chr
#   - data/subsets_fst_sampled/CEU_CHS.sampled.chrN.vcf.gz (+ .tbi)
#   - data/subsets_fst_sampled/CEU_CHS.sampled.windows.vcf.gz (+ .tbi)
#
# Dependencies:
#   - bcftools
#   - py-popgen (PPP) CLI available in env: vcf_calc.py, stat_sampler.py
#   - awk, sort, head/tail, GNU coreutils
#
# Key params (env):
#   - WIN: window size (bp); STEP: slide (bp). Default 50kb non-overlapping.
#   - NBINS: # of FST bins for uniform sampling (e.g., 3 => low/mid/high).
#   - SAMPLE_TOTAL: total # windows to sample; must be divisible by NBINS.
#   - SEED: optional random seed passed to PPP sampler for reproducibility.
#   - PPP_ENV / PPP_BIN: how to find PPP CLIs (defaults to conda env named ppp_env).
#
# Safety checks:
#   - Verifies PPP CLIs exist; verifies CEU/CHS lists not empty; ensures
#     SAMPLE_TOTAL % NBINS == 0; ensures your bcftools is NOT coming from PPP env.
# ------------------------------------------------------------------------------------

set -euo pipefail

############################################
# PPP: Windowed Fst binning + uniform sampling
# Inputs:
#   - data/subsets/CEU_CHS.chr{1..22}.vcf.gz (from your subsetting step)
#   - data/panel/integrated_call_samples_v3.20130502.ALL.panel
# Outputs:
#   - data/fst/CEU_CHS.chrN.windowed.weir.fst
#   - data/fst/CEU_CHS.allchr.windowed.weir.filtered.fst
#   - data/bins/CEU_CHS.uniform${NBINS}.N${SAMPLE_TOTAL}.windowed.weir.fst
#   - data/bins/regions.chrN.tsv
#   - data/subsets_fst_sampled/CEU_CHS.sampled.chrN.vcf.gz (+ .tbi)
#   - data/subsets_fst_sampled/CEU_CHS.sampled.windows.vcf.gz (+ .tbi)
############################################

# ---------- Config ----------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
PANEL_DIR="$ROOT/data/panel"
SUBSETS_DIR="$ROOT/data/subsets"
FST_DIR="$ROOT/data/fst"
BINS_DIR="$ROOT/data/bins"
EXTRACT_DIR="$ROOT/data/subsets_fst_sampled"
LOGS="$ROOT/logs"

# Windowing & sampling params
WIN="${WIN:-50000}"                    # window size (bp)
STEP="${STEP:-50000}"                  # step = WIN for non-overlapping
NBINS="${NBINS:-3}"                    # low/mid/high
SAMPLE_TOTAL="${SAMPLE_TOTAL:-300}"    # divisible by NBINS (e.g., 300 -> 100/bin)
SEED="${SEED:-}"                       # optional: SEED=12345 (PPP uses --random-seed)

# Use PPP from a dedicated env by absolute path (no conda run)
PPP_ENV="${PPP_ENV:-ppp_env}"
CONDA_BASE="$(conda info --base 2>/dev/null || echo "${HOME}/miniconda3")"
PPP_BIN="${PPP_BIN:-$CONDA_BASE/envs/$PPP_ENV/bin}"
VCF_CALC="${VCF_CALC:-$PPP_BIN/vcf_calc.py}"
STAT_SAMPLER="${STAT_SAMPLER:-$PPP_BIN/stat_sampler.py}"

mkdir -p "$FST_DIR" "$BINS_DIR" "$EXTRACT_DIR" "$LOGS"

# ---------- Helpers ----------
log(){ echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGS/ppp_fst.log" ; }
die(){ echo "$*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

# ---------- Dependency checks ----------
need bcftools
[[ -f "$VCF_CALC" ]] || die "PPP CLI not found: $VCF_CALC (install py-popgen in env '$PPP_ENV')"
[[ -f "$STAT_SAMPLER" ]] || die "PPP CLI not found: $STAT_SAMPLER (install py-popgen in env '$PPP_ENV')"

# (Optional) ensure we're not accidentally using bcftools from ppp_env
BCFTOOLS_PATH="$(command -v bcftools || true)"
if [[ "$BCFTOOLS_PATH" == "$CONDA_BASE/envs/$PPP_ENV/"* ]]; then
  die "bcftools is coming from $PPP_ENV ($BCFTOOLS_PATH). Activate your bcftools env (e.g., 'bcf_env') and rerun."
fi

# ---------- Sanity checks ----------
PANEL_FILE="$PANEL_DIR/integrated_call_samples_v3.20130502.ALL.panel"
[[ -f "$PANEL_FILE" ]] || die "Panel file not found: $PANEL_FILE"

have_any_vcf=false
for c in {1..22}; do
  [[ -f "$SUBSETS_DIR/CEU_CHS.chr${c}.vcf.gz" ]] && have_any_vcf=true && break
done
"$have_any_vcf" || die "No per-chromosome CEU_CHS subset VCFs found in $SUBSETS_DIR"

(( SAMPLE_TOTAL % NBINS == 0 )) || die "SAMPLE_TOTAL ($SAMPLE_TOTAL) must be divisible by NBINS ($NBINS)"

# ---------- Make CEU / CHS lists ----------
CEU="$PANEL_DIR/ceu_samples.txt"
CHS="$PANEL_DIR/chs_samples.txt"

awk '$2=="CEU"{print $1}' "$PANEL_FILE" | sort -u > "$CEU"
awk '$2=="CHS"{print $1}' "$PANEL_FILE" | sort -u > "$CHS"

CEU_N=$(wc -l < "$CEU" || echo 0)
CHS_N=$(wc -l < "$CHS" || echo 0)
(( CEU_N > 0 && CHS_N > 0 )) || die "Empty CEU/CHS list (CEU=$CEU_N, CHS=$CHS_N)"

log "CEU n=$CEU_N, CHS n=$CHS_N"
log "Params: WIN=$WIN STEP=$STEP NBINS=$NBINS SAMPLE_TOTAL=$SAMPLE_TOTAL SEED=${SEED:-<none>}"
log "Using PPP from: $PPP_BIN"
log "bcftools: $BCFTOOLS_PATH"

# ---------- 1) Compute windowed Fst per chromosome ----------
FST_FILES=()
for chr in {1..22}; do
  VCF="$SUBSETS_DIR/CEU_CHS.chr${chr}.vcf.gz"
  if [[ ! -f "$VCF" ]]; then
    log "Missing $VCF — skipping chr${chr}"
    continue
  fi

  out_prefix="$FST_DIR/CEU_CHS.chr${chr}"
  out_file="${out_prefix}.windowed.weir.fst"

  if [[ -s "$out_file" ]]; then
    log "Fst already exists: $(basename "$out_file") — skipping"
  else
    log "Computing windowed Fst for chr${chr}..."
    "$VCF_CALC" \
      --vcf "$VCF" \
      --calc-statistic windowed-weir-fst \
      --pop-file "$CEU" \
      --pop-file "$CHS" \
      --statistic-window-size "$WIN" \
      --statistic-window-step "$STEP" \
      --out-prefix "$out_prefix" \
      --overwrite
  fi

  [[ -s "$out_file" ]] && FST_FILES+=("$out_file")
done

((${#FST_FILES[@]} > 0)) || die "No per-chromosome Fst files were produced."

# ---------- Merge all chromosomes into one Fst table ----------
MERGED="$FST_DIR/CEU_CHS.allchr.windowed.weir.fst"
log "Merging $((${#FST_FILES[@]})) per-chr Fst tables -> $(basename "$MERGED")"
{
  head -n 1 "${FST_FILES[0]}"
  for f in "${FST_FILES[@]}"; do
    tail -n +2 "$f"
  done
} > "$MERGED"

# ---------- Filter out empty/invalid windows before sampling ----------
FILTERED="$FST_DIR/CEU_CHS.allchr.windowed.weir.filtered.fst"
# Expect columns: CHROM(1) BIN_START(2) BIN_END(3) N_VARIANTS(4) WEIGHTED_FST(5) MEAN_FST(6)
awk 'NR==1 || ($4+0 > 0 && $5 != "nan" && $5 != "NaN" && $5 != "NA" && $5 != ".")' \
  "$MERGED" > "$FILTERED"

# ---------- 2) Uniform sampling across Fst bins (low/mid/high) ----------
SAMPLED="$BINS_DIR/CEU_CHS.uniform${NBINS}.N${SAMPLE_TOTAL}.windowed.weir.fst"
log "Sampling windows uniformly across ${NBINS} Fst bins (N=${SAMPLE_TOTAL})..."
STAT_ARGS=( --statistic-file "$FILTERED"
            --calc-statistic windowed-weir-fst
            --sampling-scheme uniform
            --uniform-bins "$NBINS"
            --sample-size "$SAMPLE_TOTAL"
            --out "$SAMPLED"
            --overwrite )
# PPP 0.1.12 uses --random-seed
if [[ -n "${SEED}" ]]; then
  STAT_ARGS+=( --random-seed "$SEED" )
fi
"$STAT_SAMPLER" "${STAT_ARGS[@]}"

[[ -s "$SAMPLED" ]] || die "Sampling produced no output: $SAMPLED"

# ---------- 3) Convert sampled windows to per-chromosome region lists ----------
REGIONS_ALL="$BINS_DIR/CEU_CHS.sampled_regions.tsv"
awk 'NR>1 {print $1"\t"$2"\t"$3}' "$SAMPLED" > "$REGIONS_ALL"

# Split by chrom
log "Splitting sampled regions per chromosome..."
rm -f "$BINS_DIR"/regions.chr*.tsv || true
awk -v outdir="$BINS_DIR" '{print > (outdir "/regions.chr" $1 ".tsv")}' "$REGIONS_ALL"

# ---------- 4) Extract those windows from CEU+CHS per-chr VCFs ----------
TMP_LIST=()
for chr in {1..22}; do
  REG="$BINS_DIR/regions.chr${chr}.tsv"
  VCF="$SUBSETS_DIR/CEU_CHS.chr${chr}.vcf.gz"
  [[ -s "$REG" && -f "$VCF" ]] || continue

  OUT="$EXTRACT_DIR/CEU_CHS.sampled.chr${chr}.vcf.gz"
  if [[ -s "$OUT" ]]; then
    log "Already extracted sampled windows for chr${chr}"
  else
    log "Extracting sampled windows from chr${chr}..."
    bcftools view -R "$REG" -Oz -o "$OUT" "$VCF"
    bcftools index -f "$OUT"
  fi
  TMP_LIST+=("$OUT")
done

# Concatenate per-chr sampled VCFs into one file for downstream tools
if (( ${#TMP_LIST[@]} > 0 )); then
  COMBINED="$EXTRACT_DIR/CEU_CHS.sampled.windows.vcf.gz"
  log "Concatenating ${#TMP_LIST[@]} per-chr sampled VCFs -> $(basename "$COMBINED")"
  bcftools concat -a -Oz -o "$COMBINED" "${TMP_LIST[@]}"
  bcftools index -f "$COMBINED"
  log "Final sampled VCF: $COMBINED"
else
  log "No sampled regions matched existing chromosomes (nothing to extract)."
fi

# ---------- 5) Tiny post-run summary ----------
log "Summary of sampled windows:"
TOTAL_WINDOWS=$(( $(wc -l < "$SAMPLED") - 1 ))
log "  Total sampled windows: ${TOTAL_WINDOWS} (requested: ${SAMPLE_TOTAL})"
log "  By chromosome:"
awk 'NR>1{c[$1]++} END{for(k in c) printf("    chr%s: %d\n", k, c[k])}' "$SAMPLED" | sort -V | tee -a "$LOGS/ppp_fst.log" >/dev/null || true

log "✔️  PPP Fst binning + sampling complete."
