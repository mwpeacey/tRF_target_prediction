#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=8G

## Description
## Runs quality control on trimmed RNA-seq reads.

## Inputs
## $1: directory containing cutadapt-processed reads in fastq format

## Requirements
## FastQC v0.12.1

echo "#########################"
echo "New run started on $(date)"
echo "Running on samples in {$1}"
echo "#########################"

SAMPLE=$1
FASTQ_DIRECTORY=$2

cd ${FASTQ_DIRECTORY}

fastqc ${SAMPLE} -o fastqc

echo "Finished run on $(date)"


