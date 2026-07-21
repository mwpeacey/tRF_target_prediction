#!/bin/bash
#SBATCH --job-name=run_permutation
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --output=run_permutation_output.txt
#SBATCH --error=run_permutation_output.txt

## Description
## Master script for the permutation test.
##
## Shuffle target controls what gets randomised each iteration:
##   "genome" — shuffle genome windows, run real tRFs (tests genome enrichment)
##   "query"  — shuffle tRF sequences, run against real genome (tests query specificity)
##
## 1. Submits permutation_setup.sh — selects windows, runs miRanda on real
##    (unshuffled) sequences at the given score threshold.
## 2. Submits permutation_iteration.sh as a SLURM array — one shuffle
##    iteration per task, each parallelised across tRFs.
## 3. Consolidates per-iteration hit counts into shuffled_hits.tsv.
##    Run permutation_zscore.R locally to compute Z-scores and histogram data.
##
## Usage:
##   sbatch run_permutation.sh <scripts_dir> <genome.fa> <sRNA.fa> \
##          <out_dir> <run_mode> <shuffle_target> <score_threshold> \
##          <n_iterations> [fraction] [seed] [qos]
##
## Example (shuffle genome, threshold 80, 1000 iterations, 20% of genome):
##   sbatch run_permutation.sh \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction \
##     /grid/schorn/home/mpeacey/genomes/mm10/mm10.fa \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/import/mm10_tRF3b.fasta \
##     /grid/schorn/home/mpeacey/permutation_test \
##     tRF genome 80 1000 0.20 42 slow_nice
##
## Example (shuffle queries, threshold 80, 1000 iterations, 20% of genome):
##   sbatch run_permutation.sh \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction \
##     /grid/schorn/home/mpeacey/genomes/mm10/mm10.fa \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/import/mm10_tRF3b.fasta \
##     /grid/schorn/home/mpeacey/permutation_test_query \
##     tRF query 80 1000 0.20 42 slow_nice

echo "Permutation test master script started on $(date)"

if [ "$#" -lt 8 ]; then
  echo "Usage: $0 <scripts_dir> <genome.fa> <sRNA.fa> <out_dir> <run_mode> <shuffle_target> <score_threshold> <n_iterations> [fraction] [seed] [qos]" >&2
  exit 1
fi

SCRIPTS="$1"
GENOME_FA="$2"
SRNA_FA="$3"
OUTDIR="$4"
RUN_MODE="$5"
SHUFFLE_TARGET="$6"
SCORE_THRESHOLD="$7"
N_ITER="$8"
FRACTION="${9:-0.20}"
SEED="${10:-42}"
QOS="${11:-slow_nice}"

if [ "$SHUFFLE_TARGET" != "genome" ] && [ "$SHUFFLE_TARGET" != "query" ]; then
  echo "ERROR: shuffle_target must be 'genome' or 'query', got '${SHUFFLE_TARGET}'" >&2
  exit 1
fi

mkdir -p "${OUTDIR}/iterations"

echo "Configuration:"
echo "  Scripts:         ${SCRIPTS}"
echo "  Genome:          ${GENOME_FA}"
echo "  sRNA FASTA:      ${SRNA_FA}"
echo "  Output:          ${OUTDIR}"
echo "  Run mode:        ${RUN_MODE}"
echo "  Shuffle target:  ${SHUFFLE_TARGET}"
echo "  Score threshold: ${SCORE_THRESHOLD}"
echo "  Iterations:      ${N_ITER}"
echo "  Fraction:        ${FRACTION}"
echo "  Seed:            ${SEED}"
echo "  QOS:             ${QOS}"

# ── Step 1: Submit setup ───────────────────────────────────────────────────

SETUP_JOB=$(sbatch --parsable \
  --qos="${QOS}" \
  --output="${OUTDIR}/permutation_setup_output.txt" \
  --error="${OUTDIR}/permutation_setup_output.txt" \
  "${SCRIPTS}/miranda/permutation_test/permutation_setup.sh" \
  "$SCRIPTS" "$GENOME_FA" "$SRNA_FA" "$OUTDIR" "$RUN_MODE" "$SCORE_THRESHOLD" "$FRACTION" "$SEED"
)

echo "Submitted setup job: ${SETUP_JOB}"

# ── Step 2: Submit iteration array (depends on setup) ─────────────────────

ITER_JOB=$(sbatch --parsable \
  --qos="${QOS}" \
  --dependency=afterok:${SETUP_JOB} \
  --array=1-${N_ITER} \
  --output="${OUTDIR}/iterations/perm_iter_%a_output.txt" \
  --error="${OUTDIR}/iterations/perm_iter_%a_output.txt" \
  "${SCRIPTS}/miranda/permutation_test/permutation_iteration.sh" \
  "$SCRIPTS" "$SRNA_FA" "$OUTDIR" "$RUN_MODE" "$SHUFFLE_TARGET" "$SCORE_THRESHOLD"
)

echo "Submitted iteration array job: ${ITER_JOB} (${N_ITER} tasks)"

# ── Step 3: Consolidate (depends on all iterations) ───────────────────────

CONSOLIDATE_JOB=$(sbatch --parsable \
  --qos="${QOS}" \
  --dependency=afterok:${ITER_JOB} \
  --cpus-per-task=1 \
  --mem=1G \
  --time=0:10:00 \
  --job-name=perm_consolidate \
  --output="${OUTDIR}/permutation_consolidate_output.txt" \
  --error="${OUTDIR}/permutation_consolidate_output.txt" \
  --wrap="echo -e 'iteration\thits' > ${OUTDIR}/shuffled_hits.tsv && cat ${OUTDIR}/iterations/hits_*.txt | sort -n >> ${OUTDIR}/shuffled_hits.tsv && echo 'Consolidation complete.'"
)

echo "Submitted consolidation job: ${CONSOLIDATE_JOB} (depends on ${ITER_JOB})"

echo ""
echo "Pipeline submitted. Monitor with:"
echo "  squeue -u \$USER"
echo ""
echo "Once complete, run the R analysis locally:"
echo "  Rscript permutation_zscore.R ${SCORE_THRESHOLD} ${N_ITER}"
echo ""
echo "Output files:"
echo "  ${OUTDIR}/observed_hits.txt    (single integer: real hit count)"
echo "  ${OUTDIR}/shuffled_hits.tsv    (iteration + hits per shuffle)"
