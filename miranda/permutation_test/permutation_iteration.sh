#!/bin/bash
#SBATCH --job-name=perm_iter
#SBATCH --cpus-per-task=64
#SBATCH --mem=16G
#SBATCH --time=48:00:00

## Description
## A single iteration of the permutation test. Shuffles the selected window
## sequences (dinucleotide-preserving), runs miRanda on both strands in
## parallel (one process per tRF), and writes hit counts at multiple score
## thresholds.
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

echo "Permutation iteration started on $(date)"

if [ "$#" -lt 4 ]; then
  echo "Usage: $0 <scripts_dir> <sRNA.fa> <out_dir> <run_mode>" >&2
  exit 1
fi

SCRIPTS="$1"
SRNA_FA="$2"
OUTDIR="$3"
RUN_MODE="$4"

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

# Generate minus strand from shuffled plus strand
seqkit seq -r -p "${ITER_DIR}/shuffled_windows.fa" \
  > "${ITER_DIR}/shuffled_windows_minus.fa"

# ── 2. Run miRanda (parallel across tRFs) ────────────────────────────────

# Run at minimum threshold to capture all hits
if [ "$RUN_MODE" == "tRF" ]; then
  MIRANDA_FLAGS="-sc 70.0 -en 0 -scale 1.0 -loose"
elif [ "$RUN_MODE" == "miRNA" ]; then
  MIRANDA_FLAGS="-sc 150.0 -en 0"
else
  echo "ERROR: RUN_MODE must be 'tRF' or 'miRNA'" >&2
  exit 1
fi

# sRNA query files were split during setup
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

# ── 3. Count hits at multiple score thresholds ───────────────────────────

# Extract all scores
for result_file in "${RESULT_DIR}"/result_*; do
  awk '
    /Scores for this hit:/ { scores_next = 1; next }
    scores_next {
      if (/^>/) {
        split($0, fields, "\t")
        print fields[3]
      }
      scores_next = 0
    }
  ' "$result_file"
done > "${ITER_DIR}/scores.txt"

# Count at each threshold and write as a single tab-separated line:
# iteration \t hits_70 \t hits_75 \t hits_80 \t ...
LINE="${ITER}"
for CUTOFF in $(seq 70 5 165); do
  HITS=$(awk -v cutoff="$CUTOFF" '$1 >= cutoff { count++ } END { print count+0 }' "${ITER_DIR}/scores.txt")
  LINE="${LINE}\t${HITS}"
done
echo -e "$LINE" > "${OUTDIR}/iterations/hits_${ITER}.txt"

echo "[$(date)] Iteration ${ITER} complete. Hit counts by threshold:"
cat "${OUTDIR}/iterations/hits_${ITER}.txt"

# ── 4. Clean up large intermediate files ─────────────────────────────────

rm -rf "${RESULT_DIR}"
rm -f "${ITER_DIR}/shuffled_windows.fa" "${ITER_DIR}/shuffled_windows_minus.fa"
rm -f "${ITER_DIR}/scores.txt"
rmdir "${ITER_DIR}" 2>/dev/null

echo "[$(date)] Iteration ${ITER} complete."
