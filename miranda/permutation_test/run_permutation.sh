#!/bin/bash
#SBATCH --job-name=run_permutation
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --output=run_permutation_output.txt
#SBATCH --error=run_permutation_output.txt

## Description
## Master script for the genome-shuffling permutation test.
##
## 1. Submits permutation_setup.sh to select windows and run miRanda on real
##    (unshuffled) sequences (parallel across tRFs).
## 2. Once setup completes, submits permutation_iteration.sh as a SLURM array
##    job (one task per shuffle iteration, each parallelised across tRFs).
## 3. Once all iterations complete, consolidates hit counts into a single TSV.
##    Run permutation_zscore.R locally to compute Z-scores and p-values.
##
## Usage:
##   sbatch run_permutation.sh <scripts_dir> <genome.fa> <sRNA.fa> \
##          <out_dir> <run_mode> [n_iterations] [fraction] [score_cutoff] [seed]
##
## Example (development — 5% of genome, 10 iterations, score ≥ 90):
##   sbatch run_permutation.sh \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction \
##     /grid/schorn/home/mpeacey/genomes/mm10/mm10.fa \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/import/mm10_tRF3b.fasta \
##     /grid/schorn/home/mpeacey/permutation_test \
##     tRF 10 0.05 90 42
##
## Example (production — 20% of genome, 100 iterations, score ≥ 80):
##   sbatch run_permutation.sh \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction \
##     /grid/schorn/home/mpeacey/genomes/mm10/mm10.fa \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/import/mm10_tRF3b.fasta \
##     /grid/schorn/home/mpeacey/permutation_test \
##     tRF 100 0.20 80 42

echo "Permutation test master script started on $(date)"

if [ "$#" -lt 5 ]; then
  echo "Usage: $0 <scripts_dir> <genome.fa> <sRNA.fa> <out_dir> <run_mode> [n_iterations] [fraction] [score_cutoff] [seed]" >&2
  exit 1
fi

SCRIPTS="$1"
GENOME_FA="$2"
SRNA_FA="$3"
OUTDIR="$4"
RUN_MODE="$5"
N_ITER="${6:-100}"
FRACTION="${7:-0.05}"
SCORE_CUTOFF="${8:-90}"
SEED="${9:-42}"

mkdir -p "${OUTDIR}/iterations"

echo "Configuration:"
echo "  Scripts:    ${SCRIPTS}"
echo "  Genome:     ${GENOME_FA}"
echo "  sRNA FASTA: ${SRNA_FA}"
echo "  Output:     ${OUTDIR}"
echo "  Run mode:   ${RUN_MODE}"
echo "  Iterations: ${N_ITER}"
echo "  Fraction:   ${FRACTION}"
echo "  Score:      ${SCORE_CUTOFF}"
echo "  Seed:       ${SEED}"

# ── Step 1: Submit setup ───────────────────────────────────────────────────

SETUP_JOB=$(sbatch --parsable \
  --output="${OUTDIR}/permutation_setup_output.txt" \
  --error="${OUTDIR}/permutation_setup_output.txt" \
  "${SCRIPTS}/miranda/permutation_test/permutation_setup.sh" \
  "$SCRIPTS" "$GENOME_FA" "$SRNA_FA" "$OUTDIR" "$RUN_MODE" "$FRACTION" "$SCORE_CUTOFF" "$SEED"
)

echo "Submitted setup job: ${SETUP_JOB}"

# ── Step 2: Submit iteration array (depends on setup) ─────────────────────

ITER_JOB=$(sbatch --parsable \
  --dependency=afterok:${SETUP_JOB} \
  --array=1-${N_ITER} \
  --output="${OUTDIR}/iterations/perm_iter_%a_output.txt" \
  --error="${OUTDIR}/iterations/perm_iter_%a_output.txt" \
  "${SCRIPTS}/miranda/permutation_test/permutation_iteration.sh" \
  "$SCRIPTS" "$SRNA_FA" "$OUTDIR" "$RUN_MODE" "$SCORE_CUTOFF"
)

echo "Submitted iteration array job: ${ITER_JOB} (${N_ITER} tasks)"

# ── Step 3: Consolidate per-iteration hit counts (depends on all iterations) ─

CONSOLIDATE_JOB=$(sbatch --parsable \
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
echo "  Rscript ${SCRIPTS}/miranda/permutation_test/permutation_zscore.R ${OUTDIR}"
echo ""
echo "Output files:"
echo "  ${OUTDIR}/observed_hits.txt   (real hit count)"
echo "  ${OUTDIR}/shuffled_hits.tsv   (per-iteration shuffled counts)"
