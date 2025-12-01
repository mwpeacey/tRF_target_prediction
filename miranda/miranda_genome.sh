#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G

## Description
## Uses miRanda to predict tRF or miRNA targets in an input genome. Designed for use with miranda_wrapper.sh.

## Requirements
## miranda (1.9)

## Inputs
## $1 : directory containing fasta files, each of which represents a transcript to scan.
## $2 : ID for the run, used as the directory name in the output directory.
## $3 : sRNA name. Corresponds to one entry in a larger small RNA fasta file.
## $4 : determines the settings to use for miRanda. "miRNA" runs with standard seed-weighting,
##      whereas "tRF" runs with seed-weighting removed and the alignment score threshold adjusted
##      accordingly.
## $5 : output directory.

echo "New run started on $(date)"

GENOME_FASTA=$1
RUN_ID=$2
sRNA=$3
RUN_MODE=$4
OUTPUT_DIRECTORY=$5

MINUS_FA="${GENOME_FASTA%.*}_minus.${GENOME_FASTA##*.}"

echo "Output directory: ${OUTPUT_DIRECTORY}"
echo "Genome fasta: ${GENOME_FASTA}"
echo "Run ID: ${RUN_ID}"
echo "Small RNA: ${sRNA}"
echo "Run mode: ${RUN_MODE}"

cd ${OUTPUT_DIRECTORY}/${RUN_ID}/${sRNA}

echo "Scanning plus strand on $(date)."

if [ ${RUN_MODE} == 'tRF' ]; then

	miranda temp_${sRNA}.fasta \
	${GENOME_FASTA} \
        -out result_${sRNA}_plus \
        -sc 70.0 -en 0 -scale 1.0 -loose

elif [ ${RUN_MODE} == 'miRNA' ]; then

	miranda temp_${sRNA}.fasta \
        ${GENOME_FASTA} \
        -out result_${sRNA}_plus \
        -sc 150.0 -en 0

fi

grep -A1 "Scores for this hit:" result_${sRNA}_plus \
  | sort \
  | grep '>' \
  | awk -v strand="+" '{ print $0 "\t" strand }' \
  > summary_${sRNA}

echo "Scanning minus strand on $(date)."

if [ ${RUN_MODE} == 'tRF' ]; then

        miranda temp_${sRNA}.fasta \
        ${MINUS_FA} \
        -out result_${sRNA}_minus \
        -sc 70.0 -en 0 -scale 1.0 -loose

elif [ ${RUN_MODE} == 'miRNA' ]; then

        miranda temp_${sRNA}.fasta \
        ${MINUS_FA} \
        -out result_${sRNA}_minus \
        -sc 150.0 -en 0

fi

grep -A1 "Scores for this hit:" result_${sRNA}_minus \
  | sort \
  | grep '>' \
  | awk -v strand="-" '{ print $0 "\t" strand }' \
  >> summary_${sRNA}

echo "Finished on $(date)."

