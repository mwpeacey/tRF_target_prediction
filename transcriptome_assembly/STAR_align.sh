#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=32G

## Aligns adpater-trimmed reads with STAR. Settings allow for multi-mappers.

## Requirements
## STAR (v2.7.11a)

## Inputs
## $1 : full path to data directory.
## $2 : full path to directory containing adapter-trimmed reads.
## $3 : run identifier (used as output folder name in data directory).
## $4 : full path to STAR index
## $5 : Sample name
## $6 : Read 1 fastq file
## $7 : Read 2 fastq file

echo "###########################"
echo "New run started on $(date)"
echo "Analyzing data from {$1}"
echo "###########################"

DATA_DIRECTORY=$1
FASTQ_DIRECTORY=$2
RUN_NAME=$3
INDEX_DIRECTORY=$4
SAMPLE_NAME=$5
READ_1=$6
READ_2=$7

cd ${FASTQ_DIRECTORY}

STAR --runThreadN 4 --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${READ_1} ${READ_2} \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 200 \
     --chimMultimapNmax 100 \
     --outSAMstrandField intronMotif \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix ${DATA_DIRECTORY}/${RUN_NAME}/${SAMPLE_NAME}_ \
     --outFilterType BySJout

echo "Finished run on $(date)"

