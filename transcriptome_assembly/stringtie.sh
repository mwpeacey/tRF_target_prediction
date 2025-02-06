#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=16G

## Description
## Assembles a GTF annotation from RNA-seq data using stringtie with conservative settings.
## Designed for use with stringtie_wrapper.sh

## Requirements
## stringtie v3.0.0 (conda environment stringtie)

echo "#########################"
echo "New run started on $(date)"
echo "#########################"

DATA_DIRECTORY=$1
ALIGNMENT_DIRECTORY=$2
REFERENCE=$3
SAMPLE=$4
SAMPLE_NAME=$5
STRANDEDNESS=$6
RUN_ID=$7

cd ${ALIGNMENT_DIRECTORY}

echo "Data directory : ${DATA_DIRECTORY}"
echo "Alignment directory : ${ALIGNMENT_DIRECTORY}"
echo "Sample name : ${SAMPLE_NAME}"
echo "Strandedness : ${STRANDEDNESS}"

if [ ${STRANDEDNESS} -eq 0 ]
then

	stringtie -o ${DATA_DIRECTORY}/stringtie/${RUN_ID}/${SAMPLE_NAME}_stringtie-output.gtf \
	 	  -G ${REFERENCE} \
	  	  -p 4 -c 2 -f 0.05 -j 2 -s 5 \
          	  ${SAMPLE}

elif [ ${STRANDEDNESS} -eq 1 ]
then

	stringtie -o ${DATA_DIRECTORY}/stringtie/${RUN_ID}/${SAMPLE_NAME}_stringtie-output.gtf \
                  -G ${REFERENCE} \
                  -p 4 -c 2 -f 0.05 -j 2 -s 5 \
	 	  --fr \
                  ${SAMPLE}

elif [ ${STRANDEDNESS} -eq 2 ]
then

	stringtie -o ${DATA_DIRECTORY}/stringtie/${RUN_ID}/${SAMPLE_NAME}_stringtie-output.gtf \
                  -G ${REFERENCE} \
                  -p 4 -c 2 -f 0.05 -j 2 -s 5 \
                  --rf \
                  ${SAMPLE}

fi

echo "Finished run on $(date)"
