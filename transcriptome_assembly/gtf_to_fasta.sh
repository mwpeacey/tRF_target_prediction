#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G
#$ -N gtf_to_fasta
#$ -o gtf_to_fasta_output.txt
#$ -e gtf_to_fasta_output.txt

## Description
## Extracts exonic transcript sequences from a GTF file using a genome fasta as reference.
## The resulting multi-fasta is split into individiual fasta files for fasta miRanda processing.

## Inputs
## $1 : scripts root directory (e.g. /grid/schorn/home/mpeacey/scripts/tRF_target_prediction)
## $2 : directory containing the stringtie output.
## $3 : reference genome fasta (e.g. mm10_STAR/GRCm38.primary_assembly.genome.fa)

## Requirements
## Conda environment "stringtie"
## gffread 0.12.7

echo "#########################"
echo "New run started on $(date)"
echo "Running on samples in {$1}"
echo "#########################"

SCRIPTS=$1
STRINGTIE_DIRECTORY=$2
GENOME_FASTA=$3

echo "Stringtie directory : ${STRINGTIE_DIRECTORY}"
echo "Genome fasta  file : ${GENOME_FASTA}"

cd ${STRINGTIE_DIRECTORY}

# Extracts exonic DNA sequences form the GTF file.
echo "Extracting fasta..."
gffread -w stringtie_merged_filtered.fa -g ${GENOME_FASTA} stringtie_merged_filtered.gtf

# Splits the transcriptome multi-fasta into individual fasta files for faster miRanda processing.
mkdir split_transcripts
awk '/^>/ {if (seq) print seq > filename; filename="split_transcripts/" substr($0,2) ".fa"; seq=""} {if (!/^>/) seq=seq $0} END {if (seq) print seq > filename}' stringtie_merged_filtered.fa

echo "Finished run on $(date)"

mv ${SCRIPTS}/transcriptome_assembly/gtf_to_fasta_output.txt ${STRINGTIE_DIRECTORY}

