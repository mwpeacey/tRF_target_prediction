#!/bin/bash
#SBATCH --job-name=perm_iter
#SBATCH --cpus-per-task=64
#SBATCH --mem=16G
#SBATCH --time=48:00:00

## Description
## A single iteration of the permutation test. Shuffles either the pooled window
## targets (dinucleotide-preserving, k=2) or the query sequences (mononucleotide-
## preserving, k=1), then runs miRanda (one process per tRF) against the single
## pooled target file, and counts hits. Both strands are already present as
## separate records in the pooled target file, so no reverse-complementation is
## performed.
##
## Designed to be submitted as a SLURM array job by run_permutation.sh.
## The iteration number comes from SLURM_ARRAY_TASK_ID.

## Requirements
## miranda (1.9), GNU parallel
## Python 3.6+ with ushuffle

## Inputs
## $1 : Scripts directory
## $2 : Small RNA FASTA
## $3 : Output directory (same as permutation_setup.sh)
## $4 : Run mode ("tRF" or "miRNA")
## $5 : Shuffle target ("genome" or "query")
## $6 : Score threshold (e.g. 80)

echo "Permutation iteration started on $(date)"

if [ "$#" -lt 6 ]; then
  echo "Usage: $0 <scripts_dir> <sRNA.fa> <out_dir> <run_mode> <shuffle_target> <score_threshold>" >&2
  exit 1
fi

SCRIPTS="$1"
SRNA_FA="$2"
OUTDIR="$3"
RUN_MODE="$4"
SHUFFLE_TARGET="$5"
SCORE_THRESHOLD="$6"

N_CORES=${SLURM_CPUS_PER_TASK:-16}

ITER=${SLURM_ARRAY_TASK_ID}
if [ -z "$ITER" ]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID not set. Submit as an array job." >&2
  exit 1
fi

ITER_DIR="${OUTDIR}/iterations/iter_${ITER}"
mkdir -p "$ITER_DIR"

# ── 1. Shuffle sequences ────────────────────────────────────────────────

if [ "$SHUFFLE_TARGET" == "genome" ]; then
  echo "[$(date)] Iteration ${ITER}: shuffling pooled window targets..."

  python3 "${SCRIPTS}/miranda/permutation_test/shuffle_fasta.py" \
    "${OUTDIR}/real_windows.fa" \
    "${ITER_DIR}/shuffled_windows.fa" \
    --seed "$ITER" \
    --klet 2

  TARGET="${ITER_DIR}/shuffled_windows.fa"
  QUERY_DIR="${OUTDIR}/sRNA_queries"

elif [ "$SHUFFLE_TARGET" == "query" ]; then
  echo "[$(date)] Iteration ${ITER}: shuffling query sequences..."

  QUERY_DIR="${ITER_DIR}/shuffled_queries"
  mkdir -p "$QUERY_DIR"

  for SRNA_FILE in "${OUTDIR}/sRNA_queries"/*.fa; do
    python3 "${SCRIPTS}/miranda/permutation_test/shuffle_fasta.py" \
      "$SRNA_FILE" \
      "${QUERY_DIR}/$(basename "$SRNA_FILE")" \
      --seed "$ITER" \
      --klet 1
  done

  TARGET="${OUTDIR}/real_windows.fa"
else
  echo "ERROR: SHUFFLE_TARGET must be 'genome' or 'query'" >&2
  exit 1
fi

# ── 2. Run miRanda (parallel across tRFs) ────────────────────────────────

if [ "$RUN_MODE" == "tRF" ]; then
  MIRANDA_FLAGS="-sc ${SCORE_THRESHOLD}.0 -en 0 -scale 1.0 -loose"
elif [ "$RUN_MODE" == "miRNA" ]; then
  MIRANDA_FLAGS="-sc ${SCORE_THRESHOLD}.0 -en 0"
else
  echo "ERROR: RUN_MODE must be 'tRF' or 'miRNA'" >&2
  exit 1
fi

RESULT_DIR="${ITER_DIR}/results"
mkdir -p "$RESULT_DIR"

run_one_trf() {
  local SRNA_FILE="$1"
  local SRNA_NAME=$(basename "${SRNA_FILE}" .fa)

  # Single scan against the pooled target file (both strands are separate
  # records within it).
  miranda "$SRNA_FILE" "$TARGET" \
    ${MIRANDA_FLAGS} \
    -out "${RESULT_DIR}/result_${SRNA_NAME}" \
    -quiet
}
export -f run_one_trf
export TARGET MIRANDA_FLAGS RESULT_DIR

echo "[$(date)] Iteration ${ITER}: scanning with ${N_CORES} cores..."
parallel -j "$N_CORES" run_one_trf ::: "${QUERY_DIR}"/*.fa

# ── 3. Count hits and record the per-(tRF,score) histogram ────────────────
# The count feeds the existing global z-test; the score histogram is the null
# score distribution consumed by karlin_altschul_significance.py --null.

HITS=$(grep -c "Scores for this hit:" "${RESULT_DIR}"/result_* | awk -F: '{s+=$2} END {print s+0}')
echo -e "${ITER}\t${HITS}" > "${OUTDIR}/iterations/hits_${ITER}.txt"

python3 "${SCRIPTS}/miranda/permutation_test/parse_miranda_hits.py" \
  --out "${OUTDIR}/iterations/scores_${ITER}.csv" \
  "${RESULT_DIR}"/result_*

echo "[$(date)] Iteration ${ITER} complete. Hits: ${HITS}"

# ── 4. Clean up ──────────────────────────────────────────────────────────

rm -rf "${RESULT_DIR}"
rm -f "${ITER_DIR}/shuffled_windows.fa"
rm -rf "${ITER_DIR}/shuffled_queries"
rmdir "${ITER_DIR}" 2>/dev/null

echo "[$(date)] Iteration ${ITER} done."
