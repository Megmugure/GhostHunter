#!/usr/bin/env bash
# ------------------------------------------------------------------------------------
# Purpose:
#   Download 1000 Genomes Phase 3 (20130502) per-chromosome VCFs + panel, build a
#   CEU+CHS sample list, and subset each chr's VCF to those samples (chr1–22).
#
# Usage:
#   ./scripts/subset_1000G_CEU_CHS.sh                # process chr1..22
#   THREADS=8 MAKE_TBI=1 ./scripts/subset_1000G_CEU_CHS.sh 1 2 10   # only chr1,2,10
#
# Inputs (fetched if missing):
#   - 1000G panel: integrated_call_samples_v3.20130502.ALL.panel
#   - Per-chr VCFs: ALL.chrN.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz (+ .tbi)
#
# Outputs:
#   - data/panel/ceu_chs_samples.txt                # sample IDs for CEU+CHS
#   - data/subsets/CEU_CHS.chrN.vcf.gz (+ .tbi/.csi)# subsetted VCF per chromosome
#
# Dependencies:
#   - bcftools, wget, awk, sed, sort, comm
#
# Env vars (optional):
#   - THREADS: bcftools threads (default 0 => auto)
#   - MAKE_TBI: "1" => tabix index (.tbi), "0" => CSI (.csi). Default 1.
#   - WGET_OPTS: extra wget reliability flags (already set sane defaults)
#
# Notes:
#   - Prints overlap count between the CEU/CHS list and the VCF header to catch
#     sample-ID mismatches early.
#   - Skips any output file that already exists.
# ------------------------------------------------------------------------------------

set -euo pipefail
trap 'echo "Error on line $LINENO"; exit 1' ERR

# ================== Config ==================
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

# Project-relative paths (script is in project/scripts/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

RAW_DIR="$PROJECT_DIR/data/raw"         # downloaded VCFs + .tbi
PANEL_DIR="$PROJECT_DIR/data/panel"     # panel + sample lists
OUT_DIR="$PROJECT_DIR/data/subsets"     # subsetted outputs
LOG_DIR="$PROJECT_DIR/logs"

PANEL_FILE="$PANEL_DIR/integrated_call_samples_v3.20130502.ALL.panel"
SAMPLE_LIST="$PANEL_DIR/ceu_chs_samples.txt"

# Tuning
THREADS="${THREADS:-0}"      # set THREADS=8 ... to control
MAKE_TBI="${MAKE_TBI:-1}"    # 1 => tabix .tbi, 0 => default .csi
WGET_OPTS=(-q -c --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=20 --tries=5)

mkdir -p "$RAW_DIR" "$PANEL_DIR" "$OUT_DIR" "$LOG_DIR"

# ================== Helpers ==================
need() { command -v "$1" >/dev/null || { echo "Missing dependency: $1" >&2; exit 1; }; }
need bcftools
need wget

# ================== Inputs ===================
# Download panel if missing
if [[ ! -f "$PANEL_FILE" ]]; then
  echo "Downloading panel: $(basename "$PANEL_FILE")"
  wget "${WGET_OPTS[@]}" -O "$PANEL_FILE" \
    "$BASE_URL/$(basename "$PANEL_FILE")"
fi

# Build CEU+CHS sample list (trim CRLF, dedupe)
if [[ ! -f "$SAMPLE_LIST" ]]; then
  echo "Building CEU+CHS sample list -> $SAMPLE_LIST"
  # Column 1: sample ID ; Column 2: population code
  awk '$2=="CEU" || $2=="CHS" {print $1}' "$PANEL_FILE" \
    | sed 's/\r$//' \
    | awk 'length($0)>0' \
    | sort -u > "$SAMPLE_LIST"
fi

# Which chromosomes to process (args override)
if [[ $# -gt 0 ]]; then
  CHRS=("$@")
else
  mapfile -t CHRS < <(seq 1 22)
fi

# ================== Work =====================
for chr in "${CHRS[@]}"; do
  echo "Processing chromosome ${chr}..."

  VCF_BASENAME="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  VCF="$RAW_DIR/$VCF_BASENAME"
  TBI="$RAW_DIR/${VCF_BASENAME}.tbi"
  OUT_VCF="$OUT_DIR/CEU_CHS.chr${chr}.vcf.gz"

  # Fetch source VCF + index if missing
  if [[ ! -f "$VCF" ]]; then
    echo "Downloading $VCF_BASENAME"
    wget "${WGET_OPTS[@]}" -O "$VCF" "$BASE_URL/$VCF_BASENAME"
  fi
  if [[ ! -f "$TBI" ]]; then
    echo "Downloading ${VCF_BASENAME}.tbi"
    wget "${WGET_OPTS[@]}" -O "$TBI" "$BASE_URL/${VCF_BASENAME}.tbi"
  fi

  # Subset (skip if already present)
  if [[ ! -f "$OUT_VCF" ]]; then
    echo "  ✂️  Subsetting to CEU+CHS -> $(basename "$OUT_VCF")"
    # Show overlap count (helps catch mismatched sample IDs)
    OVERLAP=$(
      comm -12 <(sort -u "$SAMPLE_LIST") \
              <(bcftools query -l "$VCF" | sort -u) | wc -l
    )
    echo "Overlapping samples in VCF: $OVERLAP"

    # If overlap is zero, still run bcftools (will warn) but you likely want to fix the list.
    bcftools view \
      --samples-file "$SAMPLE_LIST" \
      --threads "$THREADS" \
      --output-type z \
      --output "$OUT_VCF" \
      "$VCF"

    # Index output: .tbi or .csi
    if [[ "$MAKE_TBI" == "1" ]]; then
      bcftools index -t --threads "$THREADS" "$OUT_VCF"
    else
      bcftools index    --threads "$THREADS" "$OUT_VCF"
    fi

    # Report final sample count
    FINAL_N=$(bcftools query -l "$OUT_VCF" | wc -l || echo 0)
    echo "Wrote $(basename "$OUT_VCF") with ${FINAL_N} samples"
  else
    echo "Already exists: $(basename "$OUT_VCF")"
  fi
done

echo -e "\nDone. Outputs in: $OUT_DIR"
