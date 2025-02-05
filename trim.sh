#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G

## Description
## Removes adapters, poly(A) sequences, and low quality bases from paired RNA-seq reads.
## Filters for minimum read length of 25 nt.
## Designed for use with trim_wrapper.sh.

## Requirements
## cutadapt (v4.5)

## Inputs
## $1 : Full path to data directory
## $2 : Sample name
## $3 : Forward adapter sequence (e.g. Truseq: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA)
## $4 : Reverse adapter sequence (e.g. Truseq: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT)

echo "#########################"
echo "New run started on $(date)"
echo "Running on samples in {$1}"
echo "#########################"

FASTQ_DIRECTORY=$1
SAMPLE_NAME=$2
FWD_ADAPTER=$3
REV_ADAPTER=$4

cd ${FASTQ_DIRECTORY}

R1=${SAMPLE_NAME}_1.fastq
R2=${SAMPLE_NAME}_2.fastq

echo "Processing sample ${SAMPLE_NAME}"

cutadapt -a ${FWD_ADAPTER} \
         -A ${REV_ADAPTER} \
         -o ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_1.fastq \
         -p ${FASTQ_DIRECTORY}/cutadapt_processed/${SAMPLE_NAME}_2.fastq \
         --poly-a --quality-cutoff=15,10 --minimum-length=25 \
         ${R1} ${R2}

echo "Finished run on $(date)"


