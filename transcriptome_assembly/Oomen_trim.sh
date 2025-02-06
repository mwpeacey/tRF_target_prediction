#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G

## Description
## For SMART-seq+5' data. Removes adapters and sorts reads by their origin (5 prime, internal, or three prime).
## Adapted from Oomen et al 2025: An atlas of transcription initiation reveals regulatory principles of gene and 
## transposable element expression in early mammalian development. Should be followed by trim_wrapper.sh with the  
## "no adapter" setting, for consistent poly A and quality trimming. Designed for use with "Oomen_trim_wrapper.sh".

## Requirements
## conda environment "standard_RNA_seq" 
## trimmomatic v0.39 
## Biopython 

## Inputs
## $1 : Full path to data directory 
## $2 : Full path to directory containing raw fastq 
## $3 : sample name

SCRIPTS=$1
DATA_DIRECTORY=$2
FASTQ_DIRECTORY=$3
SAMPLE_NAME=$4

module load EBModules
module load Biopython

cd ${FASTQ_DIRECTORY}

echo "Processing sample ${SAMPLE_NAME}"

trimmomatic PE ${SAMPLE_NAME}_1.fastq ${SAMPLE_NAME}_2.fastq -baseout ${SAMPLE_NAME}_clean.fastq -threads 1 \
               ILLUMINACLIP:${DATA_DIRECTORY}/adaptors.fa:2:30:10:1:true ILLUMINACLIP:${DATA_DIRECTORY}/other_adaptors.fa:2:30:10:1:true MINLEN:30

rm ${SAMPLE_NAME}_clean_1U.fastq
rm ${SAMPLE_NAME}_clean_2U.fastq

python ${SCRIPTS}/transcriptome_assembly/processTailedRNAseq.py ${SAMPLE_NAME}_clean_1P.fastq ${SAMPLE_NAME}_clean_2P.fastq \
       "AAGCAGTGGTATCAACGCAGAGTAC" "AAGCAGTGGTATCAACGCAGAGTACATGGG" "AAGCAGTGGTATCAACGCAGAGTACTTTTT" \
       ${SAMPLE_NAME}_five ${SAMPLE_NAME}_three ${SAMPLE_NAME}_internal

cat ${SAMPLE_NAME}_five_1.fastq >> ${SAMPLE_NAME}_all_1.fastq
cat ${SAMPLE_NAME}_five_2.fastq >> ${SAMPLE_NAME}_all_2.fastq
cat ${SAMPLE_NAME}_three_1.fastq >> ${SAMPLE_NAME}_all_1.fastq
cat ${SAMPLE_NAME}_three_2.fastq >> ${SAMPLE_NAME}_all_2.fastq
cat ${SAMPLE_NAME}_internal_1.fastq >> ${SAMPLE_NAME}_all_1.fastq
cat ${SAMPLE_NAME}_internal_2.fastq >> ${SAMPLE_NAME}_all_2.fastq

rm ${SAMPLE_NAME}_clean_1P.fastq
rm ${SAMPLE_NAME}_clean_2P.fastq

trimmomatic PE ${SAMPLE_NAME}_five_1.fastq ${SAMPLE_NAME}_five_2.fastq -baseout ${SAMPLE_NAME}_five_clean.fastq -threads 1 TRAILING:10 MINLEN:30
trimmomatic PE ${SAMPLE_NAME}_three_1.fastq ${SAMPLE_NAME}_three_2.fastq -baseout ${SAMPLE_NAME}_three_clean.fastq -threads 1 TRAILING:10 MINLEN:30
trimmomatic PE ${SAMPLE_NAME}_internal_1.fastq ${SAMPLE_NAME}_internal_2.fastq -baseout ${SAMPLE_NAME}_internal_clean.fastq -threads 1 TRAILING:10 MINLEN:30
trimmomatic PE ${SAMPLE_NAME}_all_1.fastq ${SAMPLE_NAME}_all_2.fastq -baseout ${SAMPLE_NAME}_all_clean.fastq -threads 1 TRAILING:10 MINLEN:30

rm ${SAMPLE_NAME}_five_1.fastq ${SAMPLE_NAME}_five_2.fastq
rm ${SAMPLE_NAME}_three_1.fastq ${SAMPLE_NAME}_three_2.fastq
rm ${SAMPLE_NAME}_internal_1.fastq ${SAMPLE_NAME}_internal_2.fastq
rm ${SAMPLE_NAME}_all_1.fastq ${SAMPLE_NAME}_all_2.fastq

mv ${SAMPLE_NAME}_five_clean_1P.fastq five/${SAMPLE_NAME}_1.fastq
mv ${SAMPLE_NAME}_five_clean_2P.fastq five/${SAMPLE_NAME}_2.fastq

mv ${SAMPLE_NAME}_three_clean_1P.fastq three/${SAMPLE_NAME}_1.fastq
mv ${SAMPLE_NAME}_three_clean_2P.fastq three/${SAMPLE_NAME}_2.fastq

mv ${SAMPLE_NAME}_internal_clean_1P.fastq internal/${SAMPLE_NAME}_1.fastq
mv ${SAMPLE_NAME}_internal_clean_2P.fastq internal/${SAMPLE_NAME}_2.fastq

mv ${SAMPLE_NAME}_all_clean_1P.fastq all/${SAMPLE_NAME}_1.fastq
mv ${SAMPLE_NAME}_all_clean_2P.fastq all/${SAMPLE_NAME}_2.fastq

rm ${SAMPLE_NAME}_*_clean_*U.fastq
