#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=16G

## Description

## Requirements

## Inputs

echo "New run started on $(date)"

sRNA=$1
OUTPUT_DIRECTORY=$2

echo "Output directory: ${OUTPUT_DIRECTORY}"
echo "Small RNA: ${sRNA}"

cd ${OUTPUT_DIRECTORY}/${RUN_ID}

blastn -query ${sRNA} \
	-strand plus \
	-task blastn-short \
	-db /grid/schorn/home/mpeacey/index/BLAST/mm10/mm10_genome/mm10_genome \
	-num_alignments 100 \
	-out ${OUTPUT_DIRECTORY}/${sRNA##*/}.blastn.txt \
	-word_size 7 \
	-num_threads 4

echo "Finished on $(date)."

