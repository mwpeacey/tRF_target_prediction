#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=16G
#$ -N find_ORFs
#$ -o find_ORFs_output.txt
#$ -e find_ORFs_output.txt

## Description
## Finds open reading frames in transcripts in fasta format

## Inputs
## $1 : input file (fasta format)
## $2 : prefix for output files
## $3 : output directory

## Requirements
## orfipy 0.0.4

echo "#########################"
echo "New run started on $(date)"
echo "Running on samples in {$1}"
echo "#########################"

INPUT_DIRECTORY=$1
INPUT_PREFIX=$2

cd ${INPUT_DIRECTORY}

orfipy ${INPUT_PREFIX}.fa --bed ${INPUT_PREFIX}_ORFs.bed --pep ${INPUT_PREFIX}_ORFs.faa \
       --strand f --start ATG --min 300 --outdir ${INPUT_DIRECTORY} --procs 4

echo "Finished run on $(date)"

mv /grid/schorn/home/mpeacey/scripts/standard_RNA_seq/stringtie/find_ORFs_output.txt ${INPUT_DIRECTORY}

