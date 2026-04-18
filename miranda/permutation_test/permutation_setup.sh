#!/bin/bash
#SBATCH --job-name=perm_setup
#SBATCH --cpus-per-task=64
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=permutation_setup_output.txt
#SBATCH --error=permutation_setup_output.txt

## Description
## One-time setup for the permutation test. Selects a random subset of
## non-overlapping 10 kbp windows from the genome, extracts their sequences
## (plus and minus strand), splits the sRNA FASTA into individual queries,
## and runs miRanda at the minimum score threshold on the real (unshuffled)
## sequences in parallel. Hit counts are recorded at multiple score
## thresholds (70, 75, 80, ... ) from a single miRanda run.
##
## This script is called by run_permutation.sh and should not normally be
## submitted directly.

## Requirements
## samtools (1.20), bedtools (2.31.1), seqkit (2.10.0), miranda (1.9)
## Python 3.6+, GNU parallel

## Inputs
## $1 : Scripts directory (e.g. /path/to/tRF_target_prediction)
## $2 : Genome FASTA (indexed with samtools faidx)
## $3 : Small RNA FASTA (all tRFs in one file)
## $4 : Output directory for permutation test
## $5 : Run mode ("tRF" or "miRNA")
## $6 : Genome fraction to sample (e.g. 0.05 for 5%, default: 0.05)
## $7 : Random seed (default: 42)

echo "Permutation setup started on $(date)"

if [ "$#" -lt 5 ]; then
  echo "Usage: $0 <scripts_dir> <genome.fa> <sRNA.fa> <out_dir> <run_mode> [fraction] [seed]" >&2
  exit 1
fi

SCRIPTS="$1"
GENOME_FA="$2"
SRNA_FA="$3"
OUTDIR="$4"
RUN_MODE="$5"
FRACTION="${6:-0.05}"
SEED="${7:-42}"

N_CORES=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "${OUTDIR}"

# ── 1. Select non-overlapping windows ─────────────────────────────────────

echo "[$(date)] Selecting windows (fraction=${FRACTION}, seed=${SEED})..."

python3 "${SCRIPTS}/miranda/permutation_test/select_windows.py" \
  "${GENOME_FA}" \
  "${OUTDIR}/selected_windows.bed" \
  --fraction "$FRACTION" \
  --window 10000 \
  --seed "$SEED"

N_WINDOWS=$(wc -l < "${OUTDIR}/selected_windows.bed")
echo "[$(date)] Selected ${N_WINDOWS} windows."

# ── 2. Extract real sequences ─────────────────────────────────────────────

echo "[$(date)] Extracting real window sequences..."

bedtools getfasta \
  -fi "$GENOME_FA" \
  -bed "${OUTDIR}/selected_windows.bed" \
  -name \
  -fo "${OUTDIR}/real_windows.fa"

seqkit seq -r -p "${OUTDIR}/real_windows.fa" \
  > "${OUTDIR}/real_windows_minus.fa"

# ── 3. Split sRNA FASTA into individual files ─────────────────────────────

echo "[$(date)] Splitting sRNA FASTA into individual queries..."

SRNA_DIR="${OUTDIR}/sRNA_queries"
mkdir -p "$SRNA_DIR"

awk '/^>/ {
  if (out) close(out)
  name = substr($1, 2)
  out = "'"${SRNA_DIR}"'/" name ".fa"
  print > out
  next
}
{ print >> out }' "$SRNA_FA"

N_SRNAS=$(ls "${SRNA_DIR}"/*.fa | wc -l)
echo "[$(date)] Split into ${N_SRNAS} individual sRNA files."

# ── 4. Run miRanda on real sequences (parallel across tRFs) ──────────────

# Run at minimum threshold (70 for tRF, 150 for miRNA) to capture all hits.
# Counts at higher thresholds are extracted from the same output.

echo "[$(date)] Running miRanda on real (unshuffled) windows (${N_CORES} cores)..."

if [ "$RUN_MODE" == "tRF" ]; then
  MIRANDA_FLAGS="-sc 70.0 -en 0 -scale 1.0 -loose"
elif [ "$RUN_MODE" == "miRNA" ]; then
  MIRANDA_FLAGS="-sc 150.0 -en 0"
else
  echo "ERROR: RUN_MODE must be 'tRF' or 'miRNA', got '${RUN_MODE}'" >&2
  exit 1
fi

REAL_DIR="${OUTDIR}/real_results"
mkdir -p "$REAL_DIR"

# Function run by GNU parallel: scan one tRF against both strands
run_one_trf() {
  local SRNA_FILE="$1"
  local SRNA_NAME=$(basename "${SRNA_FILE}" .fa)

  miranda "$SRNA_FILE" "${OUTDIR}/real_windows.fa" \
    ${MIRANDA_FLAGS} \
    -out "${REAL_DIR}/result_${SRNA_NAME}_plus" \
    -quiet

  miranda "$SRNA_FILE" "${OUTDIR}/real_windows_minus.fa" \
    ${MIRANDA_FLAGS} \
    -out "${REAL_DIR}/result_${SRNA_NAME}_minus" \
    -quiet
}
export -f run_one_trf
export OUTDIR MIRANDA_FLAGS REAL_DIR

parallel -j "$N_CORES" run_one_trf ::: "${SRNA_DIR}"/*.fa

echo "[$(date)] All miRanda jobs complete."

# ── 5. Count hits at multiple score thresholds ───────────────────────────

echo "[$(date)] Counting observed hits at multiple thresholds..."

# Extract all alignment scores from miRanda output into a single file.
# Score is field 6 on hit lines (lines starting with ">" after "Scores for").
for result_file in "${REAL_DIR}"/result_*; do
  awk '
    /Scores for this hit:/ { scores_next = 1; next }
    scores_next {
      if (/^>/) {
        # Split on tab; score is the 6th field in miRanda >line output
        split($0, fields, "\t")
        print fields[3]
      }
      scores_next = 0
    }
  ' "$result_file"
done > "${OUTDIR}/real_scores.txt"

# Count hits at each threshold
echo -e "cutoff\thits" > "${OUTDIR}/observed_hits.tsv"
for CUTOFF in $(seq 70 5 165); do
  HITS=$(awk -v cutoff="$CUTOFF" '$1 >= cutoff { count++ } END { print count+0 }' "${OUTDIR}/real_scores.txt")
  echo -e "${CUTOFF}\t${HITS}" >> "${OUTDIR}/observed_hits.tsv"
done

echo "[$(date)] Observed hits by threshold:"
cat "${OUTDIR}/observed_hits.tsv"
echo "[$(date)] Permutation setup complete."
