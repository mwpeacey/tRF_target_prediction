#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=8G

## Description
## Downloads the fastq file(s) associated with the input SRR code.
## Use with "download_SRA_wrapper.sh" to download all files specified in SRA explorer metadata files.

## Requirements
## sra-tools v3.0.9

## Inputs
## $1 : SRR code
## $2 : output directory

echo "###########################"
echo "New run started on $(date)"
echo "###########################"

code=$1
OUTPUT=$2

cd $OUTPUT

echo "Downloading SRR code: ${code}"

fastq-dump ${code} --split-3 --skip-technical

echo "###########################"
echo "Finished on $(date)."
echo "###########################"

