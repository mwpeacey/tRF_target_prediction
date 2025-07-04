#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G
#$ -o gather_summaries_output.txt
#$ -e gather_summaries_output.txt

## Description
## After miRanda has run, this script gathers summary files into one output folder.

## Requirements
## None

## Inputs
## $1 : miranda run directory

echo "New run started on $(date)"

RUN_DIRECTORY=$1
ID=$2

echo "Run directory: ${1}"

cd $RUN_DIRECTORY

mkdir summary_files_${ID}

for folder in */ ; do

	echo "$folder"

	cd "$folder"

	cp summary_* ${RUN_DIRECTORY}/summary_files_${ID}

	cd ${RUN_DIRECTORY}

done

cd summary_files_${ID}

find . -type f -empty -print -delete

echo "Finished on $(date)."

mv gather_summaries_output.txt ${RUN_DIRECTORY}
