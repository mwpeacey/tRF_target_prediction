#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=4G
#$ -N stringtie_wrapper
#$ -o stringtie_wrapper_output.txt
#$ -e stringtie_wrapper_output.txt

## Description
## Assembles a GTF annotation from RNA-seq data using stringtie, either with default settings
## (STRICT = F) or with settings designed to reduce false positives (STRICT = T).

## Requirements
## stringtie 2.2.1 (conda environment stringtie)

## Inputs
## $1 : full path to data home directory.
## $2 : full path to directory containing .bam output from alignment.
## $3 : full path to GTF reference annotation file.
## $4 : read strandedness (0 for unstranded, 1 for second strand, 2 for first strand).
## $5 : Not in use.

DATA_DIRECTORY=$1
ALIGNMENT_DIRECTORY=$2
REFERENCE=$3
STRANDEDNESS=$4
STRICT=$5

cd ${DATA_DIRECTORY}
mkdir stringtie

cd ${ALIGNMENT_DIRECTORY}

for SAMPLE in *.bam; do

	SAMPLE_NAME=`echo ${SAMPLE} | cut -d'_' -f 1`

	if [ ${STRICT} = T ]
	then

	#qsub -N ${SAMPLE_NAME}_stringtie \
	#-o ${DATA_DIRECTORY}/stringtie/${SAMPLE_NAME}_stringtie_strict_output.txt \
	#-e ${DATA_DIRECTORY}/stringtie/${SAMPLE_NAME}_stringtie_strict_output.txt \
	#/grid/schorn/home/mpeacey/scripts/standard_RNA_seq/stringtie/stringtie_strict.sh \
	#${DATA_DIRECTORY} ${ALIGNMENT_DIRECTORY} ${REFERENCE} ${SAMPLE} ${SAMPLE_NAME} ${STRANDEDNESS}

	else

	qsub -N ${SAMPLE_NAME}_stringtie \
        -o ${DATA_DIRECTORY}/stringtie/${SAMPLE_NAME}_stringtie_output.txt \
        -e ${DATA_DIRECTORY}/stringtie/${SAMPLE_NAME}_stringtie_output.txt \
        /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/transcriptome_assembly/stringtie.sh \
        ${DATA_DIRECTORY} ${ALIGNMENT_DIRECTORY} ${REFERENCE} ${SAMPLE} ${SAMPLE_NAME} ${STRANDEDNESS}

	fi

done
