#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G
#$ -N gtf_to_fasta
#$ -o gtf_to_fasta_output.txt
#$ -e gtf_to_fasta_output.txt

## Description
## Extracts exonic transcript sequences from a GTF file using a genome fasta as reference.
## The resulting multi-fasta is split into individual fasta files for faster miRanda processing.

## Inputs
## $1 : scripts root directory (e.g. /grid/schorn/home/mpeacey/scripts/tRF_target_prediction)
## $2 : output directory
## $3 : GTF input *prefix* i.e. "A" if input is A.gtf
## $4 : reference genome fasta (e.g. mm10_STAR/GRCm38.primary_assembly.genome.fa)

## Requirements
## Conda environment "stringtie"
## gffread 0.12.7

echo "#########################"
echo "New run started on $(date)"
echo "Running on samples in {$1}"
echo "#########################"

SCRIPTS=$1
OUTPUT_DIRECTORY=$2
GTF=$3
GENOME_FASTA=$4

echo "Output directory      : ${OUTPUT_DIRECTORY}"
echo "GTF to convert        : ${GTF}.gtf"
echo "Genome fasta file     : ${GENOME_FASTA}"

cd ${OUTPUT_DIRECTORY}

# Extracts exonic DNA sequences from the GTF file.
echo "Extracting fasta..."
gffread -w ${GTF}.fa --w-nocds -g ${GENOME_FASTA} ${GTF}.gtf

# Split transcripts for faster miranda processing
echo "Splitting transcripts..."
mkdir -p split_transcripts

# Input FASTA file
input_fasta="${GTF}.fa"

awk '/^>/ {if (f) close(f); f="split_transcripts/" substr($0,2) ".fasta"; print > f; next} {print > f}' "$input_fasta"

echo "Finished run on $(date)"

mv ${SCRIPTS}/transcriptome_assembly/gtf_to_fasta_output.txt ${OUTPUT_DIRECTORY}

