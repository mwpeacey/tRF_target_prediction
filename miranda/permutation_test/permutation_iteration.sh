#!/bin/bash
#SBATCH --job-name=perm_iter
#SBATCH --cpus-per-task=64
#SBATCH --mem=16G
#SBATCH --time=48:00:00

## Description
## A single iteration of the permutation test. Shuffles the selected window
## sequences (dinucleotide-preserving), runs miRanda on both strands in
## parallel (one process per tRF), and counts hits.
##
## Designed to be submitted as a SLURM array job by run_permutation.sh.
## The iteration number comes from SLURM_ARRAY_TASK_ID.

## Requirements
## miranda (1.9), seqkit (2.10.0), GNU parallel
## Python 3.6+ with ushuffle

## Inputs
## $1 : Scripts directory
## $2 : Small RNA FASTA
## $3 : Output directory (same as permutation_setup.sh)
## $4 : Run mode ("tRF" or "miRNA")
## $5 : Score threshold (e.g. 80)

echo "Permutation iteration started on $(date)"

if [ "$#" -lt 5 ]; then
  echo "Usage: $0 <scripts_dir> <sRNA.fa> <out_dir> <run_mode> <score_threshold>" >&2
  exit 1
fi

SCRIPTS="$1"
SRNA_FA="$2"
OUTDIR="$3"
RUN_MODE="$4"
SCORE_THRESHOLD="$5"

N_CORES=${SLURM_CPUS_PER_TASK:-16}

ITER=${SLURM_ARRAY_TASK_ID}
if [ -z "$ITER" ]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID not set. Submit as an array job." >&2
  exit 1
fi

ITER_DIR="${OUTDIR}/iterations/iter_${ITER}"
mkdir -p "$ITER_DIR"

# ── 1. Shuffle window sequences ──────────────────────────────────────────

echo "[$(date)] Iteration ${ITER}: shuffling windows..."

python3 "${SCRIPTS}/miranda/permutation_test/shuffle_fasta.py" \
  "${OUTDIR}/real_windows.fa" \
  "${ITER_DIR}/shuffled_windows.fa" \
  --seed "$ITER" \
  --klet 2

seqkit seq -r -p "${ITER_DIR}/shuffled_windows.fa" \
  > "${ITER_DIR}/shuffled_windows_minus.fa"

# ── 2. Run miRanda (parallel across tRFs) ────────────────────────────────

if [ "$RUN_MODE" == "tRF" ]; then
  MIRANDA_FLAGS="-sc ${SCORE_THRESHOLD}.0 -en 0 -scale 1.0 -loose"
elif [ "$RUN_MODE" == "miRNA" ]; then
  MIRANDA_FLAGS="-sc ${SCORE_THRESHOLD}.0 -en 0"
else
  echo "ERROR: RUN_MODE must be 'tRF' or 'miRNA'" >&2
  exit 1
fi

SRNA_DIR="${OUTDIR}/sRNA_queries"
RESULT_DIR="${ITER_DIR}/results"
mkdir -p "$RESULT_DIR"

run_one_trf() {
  local SRNA_FILE="$1"
  local SRNA_NAME=$(basename "${SRNA_FILE}" .fa)

  miranda "$SRNA_FILE" "${ITER_DIR}/shuffled_windows.fa" \
    ${MIRANDA_FLAGS} \
    -out "${RESULT_DIR}/result_${SRNA_NAME}_plus" \
    -quiet

  miranda "$SRNA_FILE" "${ITER_DIR}/shuffled_windows_minus.fa" \
    ${MIRANDA_FLAGS} \
    -out "${RESULT_DIR}/result_${SRNA_NAME}_minus" \
    -quiet
}
export -f run_one_trf
export ITER_DIR MIRANDA_FLAGS RESULT_DIR

echo "[$(date)] Iteration ${ITER}: scanning with ${N_CORES} cores..."
parallel -j "$N_CORES" run_one_trf ::: "${SRNA_DIR}"/*.fa

# ── 3. Count hits ─────────────────────────────────────────────────────────

HITS=$(grep -c "Scores for this hit:" "${RESULT_DIR}"/result_* | awk -F: '{s+=$2} END {print s+0}')
echo -e "${ITER}\t${HITS}" > "${OUTDIR}/iterations/hits_${ITER}.txt"

echo "[$(date)] Iteration ${ITER} complete. Hits: ${HITS}"

# ── 4. Clean up ──────────────────────────────────────────────────────────

rm -rf "${RESULT_DIR}"
rm -f "${ITER_DIR}/shuffled_windows.fa" "${ITER_DIR}/shuffled_windows_minus.fa"
rmdir "${ITER_DIR}" 2>/dev/null

echo "[$(date)] Iteration ${ITER} done."
