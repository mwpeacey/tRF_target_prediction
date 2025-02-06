#!/bin/bash
#$ -cwd
#$ -pe threads 4 
#$ -l m_mem_free=4G

# Description: Runs miranda with optimized settings for miRNA or tRF mode.

# Requirements: miranda, GNU parallel (for parallelization)

# Inputs:
# $1: FASTA file containing query small RNA sequences
# $2: FASTA file containing reference sequences to be scanned
# $3: Run ID
# $4: Small RNA name (e.g., the miRNA or tRF name)
# $5: Run mode (miRNA or tRF)

set -e  # Exit on error
set -o pipefail  # Fail if any command in a pipeline fails

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <small_rna_fasta> <transcriptome_fasta> <run_id> <small_rna> <run_mode>"
  exit 1
fi

SMALL_RNA_FASTA=$1
TRANSCRIPTOME_FASTA=$2
RUN_ID=$3
sRNA=$4
RUN_MODE=$5

OUTPUT_DIR="/grid/schorn/home/mpeacey/tRF_targets/data/miranda_output/${RUN_ID}/${sRNA}"
mkdir -p "$OUTPUT_DIR"
SUMMARY_FILE="${OUTPUT_DIR}/${sRNA}_summary.txt"

echo "New run started on $(date)"
echo "Run ID: ${RUN_ID}"
echo "Small RNA: ${sRNA}"
echo "Reference transcriptome: ${TRANSCRIPTOME_FASTA}"
echo "Run mode: ${RUN_MODE}"

# Define miranda parameters
if [[ ${RUN_MODE} == 'tRF' ]]; then
    MIRANDA_PARAMS="-sc 75.0 -en -20.0 -scale 1.0 -loose"
elif [[ ${RUN_MODE} == 'miRNA' ]]; then
    MIRANDA_PARAMS="-sc 120.0"
else
    echo "Invalid run mode: ${RUN_MODE}"
    exit 1
fi

# Function to process each transcript
process_transcript() {
    transcript=$1
    temp_fasta=$(mktemp)
    
    # Extract transcript sequence safely
    awk -v id="$transcript" '
        /^>/ {if (found) exit; found=index($0, id) > 0} 
        found' "$TRANSCRIPTOME_FASTA" > "$temp_fasta"

    if [[ ! -s "$temp_fasta" ]]; then
        echo "Transcript ${transcript} not found in reference." >&2
        rm -f "$temp_fasta"
        return
    fi

    result_file="${OUTPUT_DIR}/result_${sRNA}_${transcript}"
    miranda "$SMALL_RNA_FASTA" "$temp_fasta" -out "$result_file" $MIRANDA_PARAMS

    # Collect summary information if there are hits
    if grep -q "Scores for this hit:" "$result_file"; then
        grep -A 1 "Scores for this hit:" "$result_file" | grep '>' >> "$SUMMARY_FILE"
    else
        rm -f "$result_file"
    fi

    rm -f "$temp_fasta"
}

export -f process_transcript
export SMALL_RNA_FASTA TRANSCRIPTOME_FASTA sRNA OUTPUT_DIR SUMMARY_FILE MIRANDA_PARAMS

# Run in parallel with 4 threads
cat "../transcript_list.txt" | parallel -j 4 process_transcript {}

echo "Finished on $(date)."

