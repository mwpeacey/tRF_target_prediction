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
## stringtie 3.0.0 (conda environment stringtie)

## Inputs
## $1 : scripts root directory (e.g. /grid/schorn/home/mpeacey/scripts/tRF_target_prediction)
## $2 : full path to data home directory.
## $3 : full path to directory containing .bam output from alignment.
## $4 : full path to GTF reference annotation file.
## $5 : read strandedness (0 for unstranded, 1 for second strand, 2 for first strand).
## $6 : Not in use. Used to optimize parameters.

SCRIPTS=$1
DATA_DIRECTORY=$2
ALIGNMENT_DIRECTORY=$3
REFERENCE=$4
STRANDEDNESS=$5
STRICT=$6

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
	#${SCRIPTS}/transcriptome_assembly/stringtie_strict.sh \
	#${DATA_DIRECTORY} ${ALIGNMENT_DIRECTORY} ${REFERENCE} ${SAMPLE} ${SAMPLE_NAME} ${STRANDEDNESS}

	else

	qsub -N ${SAMPLE_NAME}_stringtie \
        -o ${DATA_DIRECTORY}/stringtie/${SAMPLE_NAME}_stringtie_output.txt \
        -e ${DATA_DIRECTORY}/stringtie/${SAMPLE_NAME}_stringtie_output.txt \
        ${SCRIPTS}/transcriptome_assembly/stringtie.sh \
        ${DATA_DIRECTORY} ${ALIGNMENT_DIRECTORY} ${REFERENCE} ${SAMPLE} ${SAMPLE_NAME} ${STRANDEDNESS}

	fi

done
