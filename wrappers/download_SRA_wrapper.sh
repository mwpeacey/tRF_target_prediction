#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=4G
#$ -N download_SRA_wrapper
#$ -o download_SRA_wrapper_output.txt
#$ -e download_SRA_wrapper_output.txt

## Description
## Downloads all files specified in an SRA explorer metadata file (downloaded from https://sra-explorer.info).

## Requirements
## sra-tools (v3.0.9)

## Inputs
## $1 : .tsv format metadata or .txt format URLs (specify which in input 3)
## $2 : output directory
## $3 : the format of the input. Either "full_metadata" or "URLs"

INPUT=$1
OUTPUT=$2
FORMAT=$3

if [ ${FORMAT} == 'full_metadata' ]; then

	codes=$(cut -f1 ${INPUT} | tail -n +2)

elif [ ${FORMAT} == 'URLs' ]; then

	codes=$(awk -F '/' '{print $NF}' ${INPUT} | sed 's/\.sra$//')

fi

echo "${codes}"

for code in $codes; do

	qsub -N ${code}_download_SRA \
	-o ${OUTPUT}/${code}_SRA_download_output.txt \
	-e ${OUTPUT}/${code}_SRA_download_output.txt \
	/grid/schorn/home/mpeacey/scripts/tRF_target_prediction/transcriptome_assembly/download_SRA.sh \
	${code} ${OUTPUT}

done
