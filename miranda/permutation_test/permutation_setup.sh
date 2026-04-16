#!/bin/bash
#SBATCH --job-name=perm_setup
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=permutation_setup_output.txt
#SBATCH --error=permutation_setup_output.txt

## Description
## One-time setup for the permutation test. Selects a random 20% subset of
## non-overlapping 10 kbp windows from the genome, extracts their sequences
## (plus and minus strand), and runs miRanda on the real (unshuffled) sequences
## to obtain the observed hit count.
##
## This script is called by run_permutation.sh and should not normally be
## submitted directly.

## Requirements
## samtools (1.20), bedtools (2.31.1), seqkit (2.10.0), miranda (1.9)
## Python 3.6+ with biopython

## Inputs
## $1 : Scripts directory (e.g. /path/to/tRF_target_prediction)
## $2 : Genome FASTA (indexed with samtools faidx)
## $3 : Small RNA FASTA (all tRFs in one file)
## $4 : Output directory for permutation test
## $5 : Run mode ("tRF" or "miRNA")
## $6 : Random seed (default: 42)

echo "Permutation setup started on $(date)"

if [ "$#" -lt 5 ]; then
  echo "Usage: $0 <scripts_dir> <genome.fa> <sRNA.fa> <out_dir> <run_mode> [seed]" >&2
  exit 1
fi

SCRIPTS="$1"
GENOME_FA="$2"
SRNA_FA="$3"
OUTDIR="$4"
RUN_MODE="$5"
SEED="${6:-42}"

mkdir -p "${OUTDIR}"

# ── 1. Select 20% non-overlapping windows ──────────────────────────────────

echo "[$(date)] Selecting windows (seed=${SEED})..."

# Ensure genome is indexed
if [ ! -f "${GENOME_FA}.fai" ]; then
  samtools faidx "$GENOME_FA"
fi

python3 "${SCRIPTS}/miranda/permutation_test/select_windows.py" \
  "${GENOME_FA}.fai" \
  "${OUTDIR}/selected_windows.bed" \
  --fraction 0.2 \
  --window 10000 \
  --seed "$SEED"

N_WINDOWS=$(wc -l < "${OUTDIR}/selected_windows.bed")
echo "[$(date)] Selected ${N_WINDOWS} windows."

# ── 2. Extract real sequences ───────────────────────────────────────────────

echo "[$(date)] Extracting real window sequences..."

bedtools getfasta \
  -fi "$GENOME_FA" \
  -bed "${OUTDIR}/selected_windows.bed" \
  -name \
  -fo "${OUTDIR}/real_windows.fa"

seqkit seq -r -p "${OUTDIR}/real_windows.fa" \
  > "${OUTDIR}/real_windows_minus.fa"

# ── 3. Run miRanda on real sequences ───────────────────────────────────────

echo "[$(date)] Running miRanda on real (unshuffled) windows..."

if [ "$RUN_MODE" == "tRF" ]; then
  SC=90.0
  MIRANDA_FLAGS="-sc ${SC} -en 0 -scale 1.0 -loose"
elif [ "$RUN_MODE" == "miRNA" ]; then
  SC=150.0
  MIRANDA_FLAGS="-sc ${SC} -en 0"
else
  echo "ERROR: RUN_MODE must be 'tRF' or 'miRNA', got '${RUN_MODE}'" >&2
  exit 1
fi

# Plus strand
echo "[$(date)] Scanning plus strand..."
miranda "$SRNA_FA" "${OUTDIR}/real_windows.fa" \
  ${MIRANDA_FLAGS} \
  -out "${OUTDIR}/real_result_plus" \
  -quiet

# Minus strand
echo "[$(date)] Scanning minus strand..."
miranda "$SRNA_FA" "${OUTDIR}/real_windows_minus.fa" \
  ${MIRANDA_FLAGS} \
  -out "${OUTDIR}/real_result_minus" \
  -quiet

# ── 4. Count hits ──────────────────────────────────────────────────────────

echo "[$(date)] Counting observed hits..."

# Count hit lines in miRanda output (lines starting with ">" after "Scores for")
count_hits() {
  awk '
    /Scores for this hit:/ { scores_next = 1; next }
    scores_next { if (/^>/) count++; scores_next = 0 }
    END { print count+0 }
  ' "$1"
}

HITS_PLUS=$(count_hits "${OUTDIR}/real_result_plus")
HITS_MINUS=$(count_hits "${OUTDIR}/real_result_minus")
HITS_TOTAL=$(( HITS_PLUS + HITS_MINUS ))

echo "${HITS_TOTAL}" > "${OUTDIR}/observed_hits.txt"

echo "[$(date)] Observed hits: ${HITS_PLUS} (+) + ${HITS_MINUS} (-) = ${HITS_TOTAL} total"
echo "[$(date)] Permutation setup complete."
