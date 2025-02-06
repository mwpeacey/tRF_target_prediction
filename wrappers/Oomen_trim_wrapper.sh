#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=4G
#$ -N Oomen_trim_wrapper
#$ -o Oomen_trim_wrapper_output.txt
#$ -e Oomen_trim_wrapper_output.txt

## Description
## For SMART-seq+5' data. Removes adapters and sorts reads by their origin (5 prime, internal, or three prime).
## Adapted from Oomen et al 2025: An atlas of transcription initiation reveals regulatory principles of gene and 
## transposable element expression in early mammalian development. Should be followed by trim_wrapper.sh with the
## "no adapter" setting, for consistent poly A and quality trimming. 

## Requirements
## conda environment "standard_RNA_seq"
## trimmomatic v0.39
## Biopython

## Inputs
## $1 : scripts root directory (e.g. /grid/schorn/home/mpeacey/scripts/tRF_target_prediction)
## $2 : Full path to data directory
## $3 : Full path to directory containing raw fastq

SCRIPTS=$1
DATA_DIRECTORY=$2
FASTQ_DIRECTORY=$3

cd ${FASTQ_DIRECTORY}

mkdir five
mkdir internal
mkdir three
mkdir all

for SAMPLE in *_1.fastq; do

	SAMPLE_NAME=`echo ${SAMPLE} | cut -d'_' -f 1`

	echo "Processing ${SAMPLE_NAME}"

	qsub -N ${SAMPLE_NAME}_Oomen_trim \
	-o ${FASTQ_DIRECTORY}/${SAMPLE_NAME}_Oomen_trim_output.txt \
	-e ${FASTQ_DIRECTORY}/${SAMPLE_NAME}_Oomen_trim_output.txt \
	${SCRIPTS}/transcriptome_assembly/Oomen_trim.sh \
	${SCRIPTS} ${DATA_DIRECTORY} ${FASTQ_DIRECTORY} ${SAMPLE_NAME}

done
