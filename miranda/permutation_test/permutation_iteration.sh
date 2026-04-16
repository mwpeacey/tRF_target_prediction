#!/bin/bash
#SBATCH --job-name=perm_iter
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=48:00:00

## Description
## A single iteration of the permutation test. Shuffles the selected window
## sequences (dinucleotide-preserving), runs miRanda on both strands, and
## writes the total hit count to the results directory.
##
## Designed to be submitted as a SLURM array job by run_permutation.sh.
## The iteration number comes from SLURM_ARRAY_TASK_ID.

## Requirements
## miranda (1.9), seqkit (2.10.0)
## Python 3.6+ with biopython, ushuffle

## Inputs
## $1 : Scripts directory
## $2 : Small RNA FASTA
## $3 : Output directory (same as permutation_setup.sh)
## $4 : Run mode ("tRF" or "miRNA")

echo "Permutation iteration started on $(date)"

if [ "$#" -lt 4 ]; then
  echo "Usage: $0 <scripts_dir> <sRNA.fa> <out_dir> <run_mode>" >&2
  exit 1
fi

SCRIPTS="$1"
SRNA_FA="$2"
OUTDIR="$3"
RUN_MODE="$4"

ITER=${SLURM_ARRAY_TASK_ID}
if [ -z "$ITER" ]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID not set. Submit as an array job." >&2
  exit 1
fi

ITER_DIR="${OUTDIR}/iterations/iter_${ITER}"
mkdir -p "$ITER_DIR"

# ── 1. Shuffle window sequences ────────────────────────────────────────────

echo "[$(date)] Iteration ${ITER}: shuffling windows..."

python3 "${SCRIPTS}/miranda/permutation_test/shuffle_fasta.py" \
  "${OUTDIR}/real_windows.fa" \
  "${ITER_DIR}/shuffled_windows.fa" \
  --seed "$ITER" \
  --klet 2

# Generate minus strand from shuffled plus strand
seqkit seq -r -p "${ITER_DIR}/shuffled_windows.fa" \
  > "${ITER_DIR}/shuffled_windows_minus.fa"

# ── 2. Run miRanda ─────────────────────────────────────────────────────────

if [ "$RUN_MODE" == "tRF" ]; then
  SC=90.0
  MIRANDA_FLAGS="-sc ${SC} -en 0 -scale 1.0 -loose"
elif [ "$RUN_MODE" == "miRNA" ]; then
  SC=150.0
  MIRANDA_FLAGS="-sc ${SC} -en 0"
else
  echo "ERROR: RUN_MODE must be 'tRF' or 'miRNA'" >&2
  exit 1
fi

echo "[$(date)] Iteration ${ITER}: scanning plus strand..."
miranda "$SRNA_FA" "${ITER_DIR}/shuffled_windows.fa" \
  ${MIRANDA_FLAGS} \
  -out "${ITER_DIR}/result_plus" \
  -quiet

echo "[$(date)] Iteration ${ITER}: scanning minus strand..."
miranda "$SRNA_FA" "${ITER_DIR}/shuffled_windows_minus.fa" \
  ${MIRANDA_FLAGS} \
  -out "${ITER_DIR}/result_minus" \
  -quiet

# ── 3. Count hits ──────────────────────────────────────────────────────────

count_hits() {
  awk '
    /Scores for this hit:/ { scores_next = 1; next }
    scores_next { if (/^>/) count++; scores_next = 0 }
    END { print count+0 }
  ' "$1"
}

HITS_PLUS=$(count_hits "${ITER_DIR}/result_plus")
HITS_MINUS=$(count_hits "${ITER_DIR}/result_minus")
HITS_TOTAL=$(( HITS_PLUS + HITS_MINUS ))

echo "${ITER}	${HITS_TOTAL}" > "${OUTDIR}/iterations/hits_${ITER}.txt"

echo "[$(date)] Iteration ${ITER}: ${HITS_PLUS} (+) + ${HITS_MINUS} (-) = ${HITS_TOTAL} total hits"

# ── 4. Clean up large intermediate files ───────────────────────────────────

rm -f "${ITER_DIR}/result_plus" "${ITER_DIR}/result_minus"
rm -f "${ITER_DIR}/shuffled_windows.fa" "${ITER_DIR}/shuffled_windows_minus.fa"
rmdir "${ITER_DIR}" 2>/dev/null

echo "[$(date)] Iteration ${ITER} complete."
