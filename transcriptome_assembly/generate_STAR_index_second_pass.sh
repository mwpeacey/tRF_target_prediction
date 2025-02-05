#!/bin/bash
#$ -cwd
#$ -pe threads 10
#$ -l m_mem_free=32G
#$ -N generate_STAR_index_second_pass
#$ -o generate_STAR_index_second_pass_output.txt
#$ -e generate_STAR_index_second_pass_output.txt

## Collects splice junctions from the first pass alignment, selects unique novel junctions supported by at least 3 reads
## (unique or multi-mapping) in at least 1 sample, then generates a STAR index using these new junctions as well as a
## reference annotation.

## Requirements
## STAR (v2.7.11a)

## Inputs
## $1 : full path to directory containing output of first-pass STAR alignment.
## $2 : full path to output directory.
## $3 : full path to genome fasta file.
## $4 : full path to GTF reference annotation file.

echo "###########################"
echo "New run started on $(date)"
echo "Analyzing data from {$1}"
echo "###########################"

SJ_DIRECTORY=$1
INDEX_DIRECTORY=$2
GENOME_FASTA=$3
REFERENCE=$4

cd ${SJ_DIRECTORY}

cat *SJ.out.tab | \
awk '{FS="\t";OFS="\t"} {if ($5 > 0 && ($7 > 2 || $8 > 2) && $6==0 && $4==1) print $1, $2, $3, "+"; else if ($5 > 0 && ($7 > 2 || $8 > 2) && $6==0 && $4==2) print $1, $2, $3, "-";}' | \
sort | uniq > SJ.filtered.tab

STAR --runThreadN 10 \
     --runMode genomeGenerate \
     --genomeDir ${INDEX_DIRECTORY} \
     --genomeFastaFiles ${GENOME_FASTA} \
     --sjdbGTFfile ${REFERENCE} \
     --sjdbFileChrStartEnd ${SJ_DIRECTORY}/SJ.filtered.tab \
     --sjdbOverhang 100

echo "Finished run on $(date)"

mv /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/generate_STAR_index_second_pass_output.txt ${INDEX_DIRECTORY}
