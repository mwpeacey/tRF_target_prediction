#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=4G
#$ -N QC_wrapper
#$ -o QC_wrapper_output.txt
#$ -e QC_wrapper_output.txt

## Requirements
## FastQC v0.12.1

## Inputs
## $1 : scripts root directory (e.g. /grid/schorn/home/mpeacey/scripts/tRF_target_prediction) 
## $2 : directory containing adapter-trimmed and quality-filtered fastq files 

SCRIPTS=$1
FASTQ_DIRECTORY=$2

cd ${FASTQ_DIRECTORY}
mkdir fastqc

for SAMPLE in *.fastq; do

	SAMPLE_NAME=`echo ${SAMPLE} | cut -d'.' -f 1`

	qsub -N ${SAMPLE_NAME}_QC \
	-o ${FASTQ_DIRECTORY}/fastqc/${SAMPLE_NAME}_QC_output.txt \
	-e ${FASTQ_DIRECTORY}/fastqc/${SAMPLE_NAME}_QC_output.txt \
	${SCRIPTS}/transcriptome_assembly/quality_control.sh \
	${SAMPLE} ${FASTQ_DIRECTORY}

done
