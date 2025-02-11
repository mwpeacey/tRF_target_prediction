#!/bin/bash
#$ -cwd
#$ -pe threads 4  # Request 4 threads for parallel processing
#$ -l m_mem_free=32G

## Optimized miranda.sh for faster transcript processing

## Inputs
## $1 : fasta of transcripts to scan
## $2 : ID for the run
## $3 : sRNA name
## $4 : run mode ("miRNA" or "tRF")
## $5 : output directory

TRANSCRIPTOME_FASTA=$1
RUN_ID=$2
sRNA=$3
RUN_MODE=$4
OUTPUT_DIRECTORY=$5

cd ${OUTPUT_DIRECTORY}/${RUN_ID}/${sRNA}

echo "Processing started on $(date) with $(nproc) CPUs"

## Pre-split transcriptome for faster access
if [ ! -d "split_transcripts" ]; then
    mkdir split_transcripts
    awk '/^>/ {if (seq) print seq > filename; filename="split_transcripts/" substr($0,2) ".fa"; seq=""} {if (!/^>/) seq=seq $0} END {if (seq) print seq > filename}' ${TRANSCRIPTOME_FASTA}
fi

## Count total transcripts
total_transcripts=$(wc -l < ../transcript_list.txt)
echo "Total transcripts to scan: ${total_transcripts}"

## Function to process a transcript
process_transcript() {
    transcript=$1
    transcript_file="split_transcripts/${transcript}.fa"
    if [ ! -f "$transcript_file" ]; then
        echo "Warning: Transcript file $transcript_file not found. Skipping."
        return
    fi
    
    ## Track progress
    transcript_counter=$(grep -n "${transcript}" ../transcript_list.txt | cut -d: -f1)
    echo "Scanning transcript #${transcript_counter} of ${total_transcripts}, ID ${transcript}"
    
    ## Run miranda
    if [ ${RUN_MODE} == "tRF" ]; then
        miranda temp_${sRNA}.fasta "$transcript_file" -out "result_${sRNA}_${transcript}" -sc 75.0 -en -20.0 -scale 1.0 -loose
    else
        miranda temp_${sRNA}.fasta "$transcript_file" -out "result_${sRNA}_${transcript}" -sc 120.0
    fi
    
    ## Extract hits
    grep -A 1 "Scores for this hit:" "result_${sRNA}_${transcript}" | grep '>' >> ${sRNA}_summary
    
    ## Remove result file if no hits
    grep -l "No Hits Found above Threshold" "result_${sRNA}_${transcript}" | xargs -r rm -f
}

## Run miranda in parallel for multiple transcripts
export -f process_transcript
export sRNA RUN_MODE

parallel --jobs 4 process_transcript ::: $(cat ../transcript_list.txt)

rm temp_${sRNA}.fasta

echo "Processing finished on $(date)"
