#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G
#$ -N align_seq
#$ -o align_seq_output.txt
#$ -e align_seq_output.txt

## Aligns a list of sequences (e.g. PBS, tRF) to the a genome with bowtie.

## Requirements
## bowtie, samtools, bedtools

## Inputs
## $1 : output directory
## $2 : FASTA file containing sequences to be aligned
## $3: Basename of bowtie index (e.g. mm10-M23)
## $4: index directory
## $5: run ID
## $6: max mismatches in alignment (up to 3)

OUTPUT_DIRECTORY=$1
SEQ=$2
BASENAME=$3
INDEX=$4
RUN_ID=$5
MAX_MISMATCH=$6

echo "###########################"
echo "New run with ID ${RUN_ID} started on $(date)."
echo "Output deposited in ${OUTPUT_DIRECTORY}."
echo "Aligning sequences in ${SEQ}"
echo "Allowing ${MAX_MISMATCH} mismatches"
echo "To genome with index basename ${BASENAME}"
echo "In folder ${INDEX}."
echo "###########################"

cd ${INDEX}
echo 'Aligning...'

bowtie -x ${BASENAME} \
-f ${SEQ} \
-S ${OUTPUT_DIRECTORY}/${RUN_ID}.sam \
-a \
-v ${MAX_MISMATCH}

cd ${OUTPUT_DIRECTORY}

echo 'Coverting to sorted bam...'

samtools view -bS ${RUN_ID}.sam > ${RUN_ID}.bam
samtools sort ${RUN_ID}.bam -o ${RUN_ID}.bam.sorted
samtools index ${RUN_ID}.bam.sorted

echo 'Converting to bed...'

bedtools bamtobed -i ${RUN_ID}.bam.sorted > ${RUN_ID}.bed

echo "Finished run on $(date)"
mv /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/bowtie/align_seq_output.txt ${OUTPUT_DIRECTORY}
