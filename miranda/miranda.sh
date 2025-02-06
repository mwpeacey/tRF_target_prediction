#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=4G

# Description: Runs miranda with optimized settings for miRNA or tRF mode.

# Requirements: miranda

# Inputs:
# $1: FASTA file containing query small RNA sequences
# $2: FASTA file containing reference sequences to be scanned
# $3: Run ID
# $4: Small RNA name (e.g., the miRNA or tRF name)
# $5: Run mode (miRNA or tRF)

SMALL_RNA_FASTA=$1
TRANSCRIPTOME_FASTA=$2
RUN_ID=$3
sRNA=$4
RUN_MODE=$5

OUTPUT_DIR="/grid/schorn/home/mpeacey/tRF_targets/data/miranda_output/${RUN_ID}/${sRNA}"
mkdir -p "$OUTPUT_DIR"

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <small_rna_fasta> <transcriptome_fasta> <run_id> <small_rna> <run_mode>"
  exit 1
fi

echo "New run started on $(date)"
echo "Run ID: ${RUN_ID}"
echo "Small RNA: ${sRNA}"
echo "Reference transcriptome: ${TRANSCRIPTOME_FASTA}"
echo "Run mode: ${RUN_MODE}"

# Count total number of transcripts for progress tracking
transcript_number=$(wc -l < "../transcript_list.txt")
transcript_counter=1

# Iterate through each transcript in the list
while read -r transcript; do
  echo "Scanning transcript #${transcript_counter} of ${transcript_number}, ID ${transcript}"

  # Extract the specific transcript sequence and write to a temporary file
  sed -n "/${transcript}/,/>/p" "$TRANSCRIPTOME_FASTA" | head -n -1 > "temp_${transcript}.fasta"

  # Choose miranda parameters based on the run mode
  if [[ ${RUN_MODE} == 'tRF' ]]; then
    miranda "$SMALL_RNA_FASTA" "temp_${transcript}.fasta" -out "result_${sRNA}_${transcript}" -sc 70.0 -en -20.0 -scale 1.0 -loose
  elif [[ ${RUN_MODE} == 'miRNA' ]]; then
    miranda "$SMALL_RNA_FASTA" "temp_${transcript}.fasta" -out "result_${sRNA}_${transcript}" -sc 120.0
  else
    echo "Invalid run mode: ${RUN_MODE}"
    exit 1
  fi

  # Collect summary information if there are hits
  if grep -q "Scores for this hit:" "result_${sRNA}_${transcript}"; then
    grep -A 1 "Scores for this hit:" "result_${sRNA}_${transcript}" | grep '>' >> "${sRNA}_summary"
  fi

  # Remove result file if no hits are found
  grep -q "No Hits Found above Threshold" "result_${sRNA}_${transcript}" && rm -f "result_${sRNA}_${transcript}"

  # Clean up the temporary transcript file
  rm -f "temp_${transcript}.fasta"

  ((transcript_counter++))
done < "../transcript_list.txt"

echo "Finished on $(date)."
