#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=16G
#$ -N stringtie_merge
#$ -o stringtie_merge_output.txt
#$ -e stringtie_merge_output.txt

## Description
## Merges a set of stringtie-generated GTFs together with a reference GTF.

## Inputs
## $1 : directory containing stringtie-generated GTFs to be merged.
## $2 : reference annotation GTF.

## Requirements
## stringtie 3.0.0
## gffcompare v0.12.6

echo "#########################"
echo "New run started on $(date)"
echo "Running on samples in {$1}"
echo "#########################"

STRINGTIE_DIRECTORY=$1
REFERENCE=$2

cd ${STRINGTIE_DIRECTORY}

counter=1
for file in *.gtf; do

	if [ $counter = 1 ]; then

		echo $file > file_list.txt

	else

		echo $file >> file_list.txt

	fi

	((counter=counter+1))

done

stringtie --merge -p 4 -G ${REFERENCE} -o stringtie_merged.gtf file_list.txt

# Filters the output to remove transcripts with no strand information.
awk -F "\t" '$7 != "." { print}' stringtie_merged.gtf > stringtie_merged_filtered.gtf

gffcompare -R -r ${REFERENCE} -o merged stringtie_merged_filtered.gtf

echo "Finished run on $(date)"

mv /grid/schorn/home/mpeacey/scripts/standard_RNA_seq/stringtie/stringtie_merge_output.txt ${STRINGTIE_DIRECTORY}

