#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=4G
#$ -N STAR_wrapper
#$ -o STAR_wrapper_output.txt
#$ -e STAR_wrapper_output.txt

## Wrapper for STAR_align.sh. Submits separate jobs for each sample.
## Aligns adpater-trimmed reads with STAR. Settings allow for multi-mappers.

## Requirements
## STAR (v2.7.11a)

## Inputs
## $1 : scripts root directory (e.g. /grid/schorn/home/mpeacey/scripts/tRF_target_prediction)
## $2 : full path to data directory.
## $3 : full path to directory containing adapter-trimmed reads.
## $4 : run identifier (used as output folder name in data directory).
## $5 : full path to STAR index

SCRIPTS=$1
DATA_DIRECTORY=$2
FASTQ_DIRECTORY=$3
RUN_NAME=$4
INDEX_DIRECTORY=$5

cd ${DATA_DIRECTORY}
mkdir ${RUN_NAME}

cd ${FASTQ_DIRECTORY}

for SAMPLE in *_1.fastq; do

	SAMPLE_NAME=`echo ${SAMPLE} | cut -d'_' -f 1`

	echo "${SAMPLE_NAME}"

	R1=${SAMPLE_NAME}_1.fastq
	R2=${SAMPLE_NAME}_2.fastq

	qsub -N ${SAMPLE_NAME}_STAR \
             -o ${DATA_DIRECTORY}/${RUN_NAME}/${SAMPLE_NAME}_STAR_output.txt \
             -e ${DATA_DIRECTORY}/${RUN_NAME}/${SAMPLE_NAME}_STAR_output.txt \
             ${SCRIPTS}/transcriptome_assembly/STAR_align.sh \
             ${DATA_DIRECTORY} ${FASTQ_DIRECTORY} ${RUN_NAME} ${INDEX_DIRECTORY} ${SAMPLE_NAME} ${R1} ${R2}

done
