#!/bin/bash
#SBATCH --job-name=perm_setup
#SBATCH --cpus-per-task=64
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=permutation_setup_output.txt
#SBATCH --error=permutation_setup_output.txt

## Description
## One-time setup for the permutation test. Samples a random subset of the
## PRE-GENERATED 10 kbp window FASTAs used for the real miRanda scan (plus and
## minus strand), splits the sRNA FASTA into individual queries, and runs
## miRanda on the real (unshuffled) sequences in parallel.
##
## Sampling from the same windows as the real scan guarantees an identical
## search space (boundaries, coordinates, strand handling) — no genome
## re-windowing, bedtools or seqkit needed.
##
## Called by run_permutation.sh — not normally submitted directly.

## Requirements
## miranda (1.9), Python 3.6+, GNU parallel

## Inputs
## $1 : Scripts directory
## $2 : Plus-strand window FASTA (pre-generated; e.g. *_10000bp_windows.fa).
##      The minus-strand file is found by replacing '_windows.fa' with
##      '_windows_minus.fa'; both are pooled and sampled together.
## $3 : Small RNA FASTA
## $4 : Output directory
## $5 : Run mode ("tRF" or "miRNA")
## $6 : Score threshold (e.g. 80)
## $7 : Window fraction to sample (default: 0.20)
## $8 : Random seed (default: 42)

echo "Permutation setup started on $(date)"

if [ "$#" -lt 6 ]; then
  echo "Usage: $0 <scripts_dir> <windows_plus.fa> <sRNA.fa> <out_dir> <run_mode> <score_threshold> [fraction] [seed]" >&2
  exit 1
fi

SCRIPTS="$1"
WINDOWS_PLUS="$2"
SRNA_FA="$3"
OUTDIR="$4"
RUN_MODE="$5"
SCORE_THRESHOLD="$6"
FRACTION="${7:-0.20}"
SEED="${8:-42}"

WINDOWS_MINUS="${WINDOWS_PLUS/_windows.fa/_windows_minus.fa}"
if [ ! -f "$WINDOWS_MINUS" ]; then
  echo "ERROR: expected minus-strand windows at '${WINDOWS_MINUS}' (derived from" >&2
  echo "       '${WINDOWS_PLUS}'). Rename it to match, or edit this script." >&2
  exit 1
fi

N_CORES=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "${OUTDIR}"

# ── 1. Sample windows from the pooled plus + minus target files ───────────
# Draws uniformly from both strand files (canonical chromosomes, hardcoded),
# exactly as the real scan treats plus and minus as separate targets. No
# reverse-complementation: the sampled targets are a single pooled FASTA.

echo "[$(date)] Sampling windows (fraction=${FRACTION}, seed=${SEED})..."

python3 "${SCRIPTS}/miranda/permutation_test/sample_windows_from_fasta.py" \
  "${WINDOWS_PLUS}" \
  "${WINDOWS_MINUS}" \
  "${OUTDIR}/real_windows.fa" \
  --fraction "$FRACTION" \
  --seed "$SEED"

N_WINDOWS=$(grep -c "^>" "${OUTDIR}/real_windows.fa")
echo "[$(date)] Sampled ${N_WINDOWS} pooled window targets."

# ── 3. Split sRNA FASTA into individual files ─────────────────────────────

echo "[$(date)] Splitting sRNA FASTA into individual queries..."

SRNA_DIR="${OUTDIR}/sRNA_queries"
mkdir -p "$SRNA_DIR"

# The FASTA header (e.g. ">tDR-55:76-Ala-AGC-1|tRF_1 Sprinzl_position: 55..76")
# is written unchanged into each file, so miRanda still reports the full tDRnamer
# id in its output. Only the *filename* is sanitised, because tDRnamer ids contain
# '|' and ':' which are awkward on the shell.
awk '/^>/ {
  if (out) close(out)
  name = substr($1, 2)
  fname = name
  gsub(/[|:\/ ]/, "_", fname)
  out = "'"${SRNA_DIR}"'/" fname ".fa"
  print > out
  next
}
{ print >> out }' "$SRNA_FA"

N_SRNAS=$(ls "${SRNA_DIR}"/*.fa | wc -l)
echo "[$(date)] Split into ${N_SRNAS} individual sRNA files."

# ── 4. Run miRanda on real sequences (parallel across tRFs) ──────────────

echo "[$(date)] Running miRanda at threshold ${SCORE_THRESHOLD} (${N_CORES} cores)..."

if [ "$RUN_MODE" == "tRF" ]; then
  MIRANDA_FLAGS="-sc ${SCORE_THRESHOLD}.0 -en 0 -scale 1.0 -loose"
elif [ "$RUN_MODE" == "miRNA" ]; then
  MIRANDA_FLAGS="-sc ${SCORE_THRESHOLD}.0 -en 0"
else
  echo "ERROR: RUN_MODE must be 'tRF' or 'miRNA', got '${RUN_MODE}'" >&2
  exit 1
fi

REAL_DIR="${OUTDIR}/real_results"
mkdir -p "$REAL_DIR"

run_one_trf() {
  local SRNA_FILE="$1"
  local SRNA_NAME=$(basename "${SRNA_FILE}" .fa)

  # Single scan against the pooled target file (both strands are separate
  # records within it).
  miranda "$SRNA_FILE" "${OUTDIR}/real_windows.fa" \
    ${MIRANDA_FLAGS} \
    -out "${REAL_DIR}/result_${SRNA_NAME}" \
    -quiet
}
export -f run_one_trf
export OUTDIR MIRANDA_FLAGS REAL_DIR

parallel -j "$N_CORES" run_one_trf ::: "${SRNA_DIR}"/*.fa

echo "[$(date)] All miRanda jobs complete."

# ── 5. Count hits ─────────────────────────────────────────────────────────

HITS=$(grep -c "Scores for this hit:" "${REAL_DIR}"/result_* | awk -F: '{s+=$2} END {print s+0}')
echo "$HITS" > "${OUTDIR}/observed_hits.txt"

# Per-(tRF,score) histogram of the real 20% scan — a cross-check that the null
# (shuffled) score distribution matches the real one in the chance-dominated bulk.
python3 "${SCRIPTS}/miranda/permutation_test/parse_miranda_hits.py" \
  --out "${OUTDIR}/observed_score_histogram.csv" \
  "${REAL_DIR}"/result_*

echo "[$(date)] Observed hits at threshold ${SCORE_THRESHOLD}: ${HITS}"
echo "[$(date)] Permutation setup complete."
