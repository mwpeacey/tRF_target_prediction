#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=16G
#$ -o make_blastdb_output.txt
#$ -e make_blastdb_output.txt

## Description

## Requirements

## Inputs

echo "New run started on $(date)"

OUTPUT_DIRECTORY=$1
FASTA=$2
ID=$3

echo "Output directory: ${OUTPUT_DIRECTORY}"
echo "Fasta input: ${FASTA}"

cd ${OUTPUT_DIRECTORY}

makeblastdb -in ${FASTA} -dbtype nucl -out ${ID} -title ${ID}

echo "Finished on $(date)."

