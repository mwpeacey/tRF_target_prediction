#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G
#$ -N blast_ORFs
#$ -o blast_ORFs_output.txt
#$ -e blast_ORFs_output.txt

## Description
## Takes open reading frames found in assembled transcripts (by find_ORFs.sh) and blasts
## them against the mouse refseq database.

## Inputs
## $1 : input directory
## $2 : input prefix

## Requirements
## blast 2.16.0

echo "#########################"
echo "New run started on $(date)"
echo "Running on samples in {$1}"
echo "#########################"

INPUT_DIRECTORY=$1
INPUT_PREFIX=$2

cd ${INPUT_DIRECTORY}

blastp -query ${INPUT_PREFIX}_ORFs.faa \
       -db /grid/schorn/home/mpeacey/annotations/mm10/refseq/mouse_protein_combined.faa \
       -max_target_seqs 1 -outfmt 10 -out ${INPUT_PREFIX}_ORFs_blast.csv

echo "Finished run on $(date)"

mv /grid/schorn/home/mpeacey/scripts/standard_RNA_seq/stringtie/blast_ORFs_output.txt ${INPUT_DIRECTORY}

