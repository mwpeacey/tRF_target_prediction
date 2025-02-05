#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G

## Description
## Removes poly(A) sequences and low quality bases (but not adapters) from paired RNA-seq reads.
## Filters for minimum read length of 25 nt.
## Designed for use with trim_wrapper.sh.

## Requirements
## cutadapt (v4.5)

## Inputs
## $1 : Full path to data directory
## $2 : Sample name

echo "#########################"
echo "New run started on $(date)"
echo "Running on samples in {$1}"
echo "#########################"

FASTQ_DIRECTORY=$1
SAMPLE_NAME=$2

cd ${FASTQ_DIRECTORY}

R1=${SAMPLE_NAME}_1.fastq
R2=${SAMPLE_NAME}_2.fastq

echo "Processing sample ${SAMPLE_NAME}"

cutadapt -o ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_1.fastq \
         -p ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_2.fastq \
         --poly-a --quality-cutoff=15,10 --minimum-length=25 \
         ${R1} ${R2}

echo "Finished run on $(date)"


