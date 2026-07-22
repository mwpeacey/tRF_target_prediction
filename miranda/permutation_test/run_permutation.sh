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
##   sbatch run_permutation.sh <scripts_dir> <windows_plus.fa> <sRNA.fa> \
##          <out_dir> <run_mode> <shuffle_target> <score_threshold> \
##          <n_iterations> [fraction] [seed] [qos] [chroms]
##
## <windows_plus.fa> is the PRE-GENERATED plus-strand window FASTA used for the
## real scan (e.g. *_10000bp_windows.fa); the minus-strand file is derived by
## replacing '_windows.fa' with '_windows_minus.fa'.
##
## Example (shuffle genome, threshold 70, 1000 iterations, 20% of windows,
## restricted to the canonical chromosomes):
##   sbatch run_permutation.sh \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction \
##     /grid/schorn/home/mpeacey/mouse_two_cell_transcriptome/target_prediction/GRCm38.primary_assembly.genome_10000bp_windows.fa \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/import/tDRnamer/mm10-tDR.fa \
##     /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/import/miranda/genome_shuffle_null \
##     tRF genome 70 1000 0.20 42 slow_nice \
##     chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX,chrY

echo "Permutation test master script started on $(date)"

if [ "$#" -lt 8 ]; then
  echo "Usage: $0 <scripts_dir> <windows_plus.fa> <sRNA.fa> <out_dir> <run_mode> <shuffle_target> <score_threshold> <n_iterations> [fraction] [seed] [qos] [chroms]" >&2
  exit 1
fi

SCRIPTS="$1"
WINDOWS_PLUS="$2"
SRNA_FA="$3"
OUTDIR="$4"
RUN_MODE="$5"
SHUFFLE_TARGET="$6"
SCORE_THRESHOLD="$7"
N_ITER="$8"
FRACTION="${9:-0.20}"
SEED="${10:-42}"
QOS="${11:-slow_nice}"
CHROMS="${12:-}"   # comma-separated chromosomes to restrict windows to (default: all)

if [ "$SHUFFLE_TARGET" != "genome" ] && [ "$SHUFFLE_TARGET" != "query" ]; then
  echo "ERROR: shuffle_target must be 'genome' or 'query', got '${SHUFFLE_TARGET}'" >&2
  exit 1
fi

mkdir -p "${OUTDIR}/iterations"

echo "Configuration:"
echo "  Scripts:         ${SCRIPTS}"
echo "  Windows (plus):  ${WINDOWS_PLUS}"
echo "  sRNA FASTA:      ${SRNA_FA}"
echo "  Output:          ${OUTDIR}"
echo "  Run mode:        ${RUN_MODE}"
echo "  Shuffle target:  ${SHUFFLE_TARGET}"
echo "  Score threshold: ${SCORE_THRESHOLD}"
echo "  Iterations:      ${N_ITER}"
echo "  Fraction:        ${FRACTION}"
echo "  Seed:            ${SEED}"
echo "  QOS:             ${QOS}"
echo "  Chromosomes:     ${CHROMS:-all}"

# ── Step 1: Submit setup ───────────────────────────────────────────────────

SETUP_JOB=$(sbatch --parsable \
  --qos="${QOS}" \
  --output="${OUTDIR}/permutation_setup_output.txt" \
  --error="${OUTDIR}/permutation_setup_output.txt" \
  "${SCRIPTS}/miranda/permutation_test/permutation_setup.sh" \
  "$SCRIPTS" "$WINDOWS_PLUS" "$SRNA_FA" "$OUTDIR" "$RUN_MODE" "$SCORE_THRESHOLD" "$FRACTION" "$SEED" "$CHROMS"
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
  --wrap="echo -e 'iteration\thits' > ${OUTDIR}/shuffled_hits.tsv && cat ${OUTDIR}/iterations/hits_*.txt | sort -n >> ${OUTDIR}/shuffled_hits.tsv && awk -F, 'FNR==1{next} {c[\$1 FS \$2]+=\$3} END{print \"tRF,alignment_score,count\"; for(k in c) print k FS c[k]}' ${OUTDIR}/iterations/scores_*.csv > ${OUTDIR}/shuffled_score_histogram.csv && echo 'Consolidation complete.'"
)

echo "Submitted consolidation job: ${CONSOLIDATE_JOB} (depends on ${ITER_JOB})"

echo ""
echo "Pipeline submitted. Monitor with:"
echo "  squeue -u \$USER"
echo ""
echo "Once complete:"
echo "  # Global enrichment z-test (aggregate hit count):"
echo "  Rscript permutation_zscore.R ${SCORE_THRESHOLD} ${N_ITER} ${OUTDIR}"
echo ""
echo "  # Per-prediction Karlin-Altschul E-value / FDR, using the shuffled null"
echo "  # to score your REAL whole-genome predictions:"
echo "  python3 miranda/permutation_test/karlin_altschul_significance.py \\"
echo "      import/miranda/miranda_output_${SCORE_THRESHOLD}.csv \\"
echo "      import/miranda/miranda_output_${SCORE_THRESHOLD}_ka.csv \\"
echo "      --null ${OUTDIR}/shuffled_score_histogram.csv \\"
echo "      --null-fraction ${FRACTION} --null-iterations ${N_ITER} \\"
echo "      --summary import/miranda/ka_summary.csv"
echo ""
echo "Output files:"
echo "  ${OUTDIR}/observed_hits.txt              (single integer: real hit count)"
echo "  ${OUTDIR}/shuffled_hits.tsv              (iteration + hits per shuffle)"
echo "  ${OUTDIR}/shuffled_score_histogram.csv   (tRF,score,count — the NULL distribution)"
