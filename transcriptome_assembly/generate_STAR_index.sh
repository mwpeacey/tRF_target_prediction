#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=32G
#$ -N generate_STAR_index
#$ -o generate_STAR_index_output.txt
#$ -e generate_STAR_index_output.txt

## Requirements:
## STAR (v2.7.11a)

## Inputs:
## $1: Genome directory (e.g. /grid/schorn/home/mpeacey/index/STAR/mm10_STAR_index)
## $2: Genome FASTA (e.g. /grid/schorn/home/mpeacey/index/mm10/GRCm38.primary_assembly.genome.fa)
## $3: Gene annotation (e.g. /grid/schorn/home/mpeacey/index/mm10/gencode.vM10.primary_assembly.annotation.gtf)

echo "#####################################"
echo "New run started on $(date)"
echo "#####################################"

echo "Genome directory: $1"
echo "Genome FASTA: $2"
echo "Genome annotaiton: $3"

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir $1 \
--genomeFastaFiles $2 \
--sjdbGTFfile $3 \
--sjdbOverhang 100

mv generate_STAR_index_output.txt $1

echo "Finished run on $(date)"
