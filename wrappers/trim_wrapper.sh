#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=4G
#$ -N trim_wrapper
#$ -o trim_wrapper_output.txt
#$ -e trim_wrapper_output.txt

## Description
## Removes adapters, poly(A) sequences, and low quality bases from paired RNA-seq reads.

## Requirements
## cutadapt v4.5

## Inputs
## $1 : scripts root directory (e.g. /grid/schorn/home/mpeacey/scripts/tRF_target_prediction)
## $2 : Full path to data directory
## $3 : include adapter trimming (T) or skip if adapters have already been trimmed (F).
##      If F then no adapter sequences need be specified. 
## $4 : Forward adapter sequence (e.g. Truseq: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA)
## $5 : Reverse adapter sequence (e.g. Truseq: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT)

SCRIPTS=$1
FASTQ_DIRECTORY=$2
ADAPTER=$3
FWD_ADAPTER=$4
REV_ADAPTER=$5

cd ${FASTQ_DIRECTORY}
mkdir cutadapt_processed

for SAMPLE in *_1.fastq; do

	SAMPLE_NAME=`echo ${SAMPLE} | cut -d'_' -f 1`

	if [ ${ADAPTER} = T ]
	then

		qsub -N ${SAMPLE_NAME}_trim \
		-o ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_trim_output.txt \
		-e ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_trim_output.txt \
		${SCRIPTS}/transcriptome_assembly/trim.sh \
		${FASTQ_DIRECTORY} ${SAMPLE_NAME} ${FWD_ADAPTER} ${REV_ADAPTER}

	else

		qsub -N ${SAMPLE_NAME}_trim \
                -o ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_trim_output.txt \
                -e ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_trim_output.txt \
                ${SCRIPTS}/transcriptome_assembly/trim_no_adapter.sh \
                ${FASTQ_DIRECTORY} ${SAMPLE_NAME}

	fi

done
