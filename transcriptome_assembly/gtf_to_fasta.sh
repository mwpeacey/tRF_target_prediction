#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=32G
#$ -N gtf_to_fasta
#$ -o gtf_to_fasta_output.txt
#$ -e gtf_to_fasta_output.txt

## Description
## Assembles custom transcripts from STAR-aligned RNA-seq reads and extracts exonic DNA sequences from the resulting GTF file.

## Inputs
## $1 : alignment directory (e.g. /grid/schorn/home/mpeacey/tRF_targets/data/public_RNAseq/Abe_2015/STAR_alignment)

## Requirements
## stringtie
## bedtools

echo "#########################"
echo "New run started on $(date)"
echo "Running on samples in {$1}"
echo "#########################"

TRANSCRIPT_SET_ID=$1
STRINGTIE_DIRECTORY=$2
GENOME_FASTA=$3

echo "Genome fasta  file : ${GENOME_FASTA}"

cd ${STRINGTIE_DIRECTORY}

gffread -w ${TRANSCRIPT_SET_ID}_assembled_transcripts.fa -g ${GENOME_FASTA} stringtie_merged_filtered.gtf

echo "Finished run on $(date)"

mv /grid/schorn/home/mpeacey/scripts/standard_RNA_seq/stringtie/gtf_to_fasta_output.txt ${STRINGTIE_DIRECTORY}/gtf_to_fasta_output.txt
