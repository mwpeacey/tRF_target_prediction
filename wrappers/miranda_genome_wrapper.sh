#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=4G
#$ -N miranda_wrapper
#$ -o miranda_wrapper_output.txt
#$ -e miranda_wrapper_output.txt

## Description
## Uses miRanda to predict tRF or miRNA targets in an input transcriptome.

## Requirements
## miRanda (1.9)
## seqkit ()

## Inputs
## $1 : scripts root directory (e.g. /grid/schorn/home/mpeacey/scripts/tRF_target_prediction)
## $2 : fasta of small RNAs (e.g. tRF3b sequences).
## $3 : directory containing fasta files, each of which represents a transcript to scan.
## $4 : output directory.
## $5 : ID for the run, used as the directory name in the output directory.
##	Specifiy the name of the previous run with $7.

SCRIPTS=$1
SMALL_RNA_FASTA=$2
GENOME_FASTA=$3
OUTPUT_DIRECTORY=$4
RUN_ID=$5
RUN_MODE=$6

cd ${OUTPUT_DIRECTORY}

mkdir -p ${RUN_ID}
cd ${RUN_ID}

if  [[ ! -f sRNA_list.txt ]]; then

	cat ${SMALL_RNA_FASTA} | grep '>' > temp_sRNA_list.txt

	i=1
	for line in $(cat temp_sRNA_list.txt); do

        	result=${line#*>}

        	if [ $i = 1 ]; then

                	echo $result > sRNA_list.txt

        	else

            		echo $result >> sRNA_list.txt

        	fi

		((i=i+1))

	done

fi

rm temp_sRNA_list.txt

for sRNA in $(cat sRNA_list.txt); do

	cd ${OUTPUT_DIRECTORY}/${RUN_ID}

	echo "Submitting job for sRNA ${sRNA}..."

	mkdir -p ${sRNA}
	cd ${sRNA}

	cat ${SMALL_RNA_FASTA} | sed -n "/${sRNA}/,/>/p" | head -2 > temp_${sRNA}.fasta

	qsub -N ${RUN_ID}_${sRNA}_miranda \
	-o ${OUTPUT_DIRECTORY}/${RUN_ID}/${sRNA}/${RUN_ID}_${sRNA}_miranda_output.txt \
	-e ${OUTPUT_DIRECTORY}/${RUN_ID}/${sRNA}/${RUN_ID}_${sRNA}_miranda_output.txt \
	${SCRIPTS}/miranda/miranda_genome.sh \
	${GENOME_FASTA} ${RUN_ID} \
	${sRNA} ${RUN_MODE} ${OUTPUT_DIRECTORY}

done
