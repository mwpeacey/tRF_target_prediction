#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=16G
#$ -N stringtie_merge
#$ -o stringtie_merge_output.txt
#$ -e stringtie_merge_output.txt

## Description
## Merges a set of stringtie-generated GTFs together with a reference GTF,
## and filters out transcripts without strand information.

## Inputs
## $1 : scripts root directory (e.g. /grid/schorn/home/mpeacey/scripts/tRF_target_prediction)
## $2 : directory containing stringtie-generated GTFs to be merged.
## $3 : reference annotation GTF (e.g. mm10_STAR/gencode.vM23.primary_assembly.annotation.gtf)

## Requirements
## Conda environment "stringtie"
## stringtie 3.0.0
## gffcompare 0.12.6

echo "#########################"
echo "New run started on $(date)"
echo "Running on samples in {$1}"
echo "#########################"

SCRIPTS=$1
STRINGTIE_DIRECTORY=$2
REFERENCE=$3

echo "Stringtie directory : ${STRINGTIE_DIRECTORY}"
echo "Reference annotation : ${REFERENCE}"

cd ${STRINGTIE_DIRECTORY}

# Gathers gtf annotations from individual samples
echo "Gathering GTFs..."
counter=1
for file in *.gtf; do

	if [ $counter = 1 ]; then

		echo $file > file_list.txt

	else

		echo $file >> file_list.txt

	fi

	((counter=counter+1))

done

# Runs stringtie in merge mode
echo "Merging transcriptome annotations..."
stringtie --merge -p 4 -G ${REFERENCE} -o stringtie_merged.gtf file_list.txt

# Filters the output to remove transcripts with no strand information.
echo "Filtering output..."
awk -F "\t" '$7 != "." { print}' stringtie_merged.gtf > stringtie_merged_filtered.gtf

# Generates statistics about novel transcript variants
echo "Generating statistics..."
gffcompare -R -r ${REFERENCE} -o merged stringtie_merged_filtered.gtf

echo "Finished run on $(date)"

mv ${SCRIPTS}/transcriptome_assembly/stringtie_merge_output.txt ${STRINGTIE_DIRECTORY}

