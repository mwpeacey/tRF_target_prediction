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
## $1 : Full path to data directory
## $2 : include adapter trimming (T) or skip if adapters have already been trimmed (F).
##      If F then no adapter sequences need be specified. 
## $3 : Forward adapter sequence (e.g. Truseq: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA)
## $4 : Reverse adapter sequence (e.g. Truseq: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT)

FASTQ_DIRECTORY=$1
ADAPTER=$2
FWD_ADAPTER=$3
REV_ADAPTER=$4

cd ${FASTQ_DIRECTORY}
mkdir cutadapt_processed

for SAMPLE in *_1.fastq; do

	SAMPLE_NAME=`echo ${SAMPLE} | cut -d'_' -f 1`

	if [ ${ADAPTER} = T ]
	then

		qsub -N ${SAMPLE_NAME}_trim \
		-o ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_trim_output.txt \
		-e ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_trim_output.txt \
		/grid/schorn/home/mpeacey/scripts/tRF_target_prediction/transcriptome_assembly/trim.sh \
		${FASTQ_DIRECTORY} ${SAMPLE_NAME} ${FWD_ADAPTER} ${REV_ADAPTER}

	else

		qsub -N ${SAMPLE_NAME}_trim \
                -o ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_trim_output.txt \
                -e ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_trim_output.txt \
                /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/transcriptome_assembly/trim_no_adapter.sh \
                ${FASTQ_DIRECTORY} ${SAMPLE_NAME}

	fi

done
