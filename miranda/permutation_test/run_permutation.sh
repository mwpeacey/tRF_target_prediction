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
##    (unshuffled) sequences.
## 2. Once setup completes, submits permutation_iteration.sh as a SLURM array
##    job (one task per shuffle iteration).
## 3. Once all iterations complete, submits the R analysis script to compute
##    Z-scores and p-values.
##
## Usage:
##   sbatch run_permutation.sh <scripts_dir> <genome.fa> <sRNA.fa> \
##          <out_dir> <run_mode> [n_iterations] [seed]
##
## Example:
##   sbatch run_permutation.sh \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction \
##     /grid/schorn/home/mpeacey/genomes/mm10/mm10.fa \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/import/mm10_tRF3b.fasta \
##     /grid/schorn/home/mpeacey/permutation_test \
##     tRF 100 42

echo "Permutation test master script started on $(date)"

if [ "$#" -lt 5 ]; then
  echo "Usage: $0 <scripts_dir> <genome.fa> <sRNA.fa> <out_dir> <run_mode> [n_iterations] [seed]" >&2
  exit 1
fi

SCRIPTS="$1"
GENOME_FA="$2"
SRNA_FA="$3"
OUTDIR="$4"
RUN_MODE="$5"
N_ITER="${6:-100}"
SEED="${7:-42}"

mkdir -p "${OUTDIR}/iterations"

echo "Configuration:"
echo "  Scripts:    ${SCRIPTS}"
echo "  Genome:     ${GENOME_FA}"
echo "  sRNA FASTA: ${SRNA_FA}"
echo "  Output:     ${OUTDIR}"
echo "  Run mode:   ${RUN_MODE}"
echo "  Iterations: ${N_ITER}"
echo "  Seed:       ${SEED}"

# ── Step 1: Submit setup ───────────────────────────────────────────────────

SETUP_JOB=$(sbatch --parsable \
  --output="${OUTDIR}/permutation_setup_output.txt" \
  --error="${OUTDIR}/permutation_setup_output.txt" \
  "${SCRIPTS}/miranda/permutation_test/permutation_setup.sh" \
  "$SCRIPTS" "$GENOME_FA" "$SRNA_FA" "$OUTDIR" "$RUN_MODE" "$SEED"
)

echo "Submitted setup job: ${SETUP_JOB}"

# ── Step 2: Submit iteration array (depends on setup) ─────────────────────

ITER_JOB=$(sbatch --parsable \
  --dependency=afterok:${SETUP_JOB} \
  --array=1-${N_ITER} \
  --output="${OUTDIR}/iterations/perm_iter_%a_output.txt" \
  --error="${OUTDIR}/iterations/perm_iter_%a_output.txt" \
  "${SCRIPTS}/miranda/permutation_test/permutation_iteration.sh" \
  "$SCRIPTS" "$SRNA_FA" "$OUTDIR" "$RUN_MODE"
)

echo "Submitted iteration array job: ${ITER_JOB} (${N_ITER} tasks)"

# ── Step 3: Consolidate per-iteration hit counts, then run analysis ───────

ANALYSIS_JOB=$(sbatch --parsable \
  --dependency=afterok:${ITER_JOB} \
  --cpus-per-task=1 \
  --mem=4G \
  --time=1:00:00 \
  --job-name=perm_analysis \
  --output="${OUTDIR}/permutation_analysis_output.txt" \
  --error="${OUTDIR}/permutation_analysis_output.txt" \
  --wrap="echo -e 'iteration\thits' > ${OUTDIR}/shuffled_hits.tsv && cat ${OUTDIR}/iterations/hits_*.txt | sort -n >> ${OUTDIR}/shuffled_hits.tsv && Rscript ${SCRIPTS}/miranda/permutation_test/permutation_zscore.R ${OUTDIR}"
)

echo "Submitted analysis job: ${ANALYSIS_JOB} (depends on ${ITER_JOB})"

echo ""
echo "Pipeline submitted. Monitor with:"
echo "  squeue -u \$USER"
echo ""
echo "Results will be written to:"
echo "  ${OUTDIR}/observed_hits.txt       (real hit count)"
echo "  ${OUTDIR}/shuffled_hits.tsv       (per-iteration shuffled counts)"
echo "  ${OUTDIR}/permutation_results.csv (Z-score, p-values)"
