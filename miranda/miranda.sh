#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=32G

## Description
## Uses miRanda to predict tRF or miRNA targets in an input transcriptome. Designed for use with miranda_wrapper.sh.

## Requirements
## miranda v1.9

## Inputs
## $1 : fasta of transcripts to scan.
## $2 : ID for the run, used as the directory name in the output directory.
## $3 : sRNA name. Corresponds to one entry in a larger small RNA fasta file.
## $4 : determines the settings to use for miRanda. "miRNA" runs with standard seed-weighting,
##      whereas "tRF" runs with seed-weighting removed and the alignment score threshold adjusted
##      accordingly.
## $5 : output directory.

echo "New run started on $(date)"

TRANSCRIPTOME_FASTA=$1
RUN_ID=$2
sRNA=$3
RUN_MODE=$4
OUTPUT_DIRECTORY=$5

echo "Output directory: ${OUTPUT_DIRECTORY}"
echo "Run ID: ${RUN_ID}"
echo "Small RNA: ${sRNA}"
echo "Reference transcriptome: ${TRANSCRIPTOME_FASTA}"
echo "Run mode: ${RUN_MODE}"

cd ${OUTPUT_DIRECTORY}/${RUN_ID}

echo "Iterating through transcripts..."

transcript_counter=1
transcript_number=$(echo "$(cat transcript_list.txt)" | wc -l)

cd ${OUTPUT_DIRECTORY}/${RUN_ID}/${sRNA}

for transcript in $(cat ../transcript_list.txt); do

	echo "Scanning transcript #${transcript_counter} of ${transcript_number}, ID ${transcript}"

	## Move the target transcript to a temporary fasta file

        ##cat ${TRANSCRIPTOME_FASTA} | sed -n "/${transcript}/,/>/p" | head -n -1 > temp_${transcript}.fasta

	## Run miranda using that temporary file

	if [ ${RUN_MODE} == 'tRF' ]; then

		miranda temp_${sRNA}.fasta \
			${OUTPUT_DIRECTORY}/${RUN_ID}/split_transcripts/${transcript}.fa \
			-out result_${sRNA}_${transcript} \
			-sc 75.0 -en -20.0 -scale 1.0 -loose

	fi

	if [ ${RUN_MODE} == 'miRNA' ]; then

		miranda temp_${sRNA}.fasta \
                        ${OUTPUT_DIRECTORY}/${RUN_ID}/split_transcripts/${transcript}.fa \
                        -out result_${sRNA}_${transcript} \
                        -sc 120

	fi

	# rm temp_${transcript}.fasta

	## If hits, move the summary line to a summary file

	grep -A 1 "Scores for this hit:" result_${sRNA}_${transcript} | sort | grep '>' > temp_summary_${transcript_counter}

	if [ $transcript_counter = 1 ]; then

	               	mv temp_summary_1 ${sRNA}_summary

	        else

	            	cat temp_summary_$transcript_counter >> ${sRNA}_summary

	        fi

	rm temp_summary_${transcript_counter}

	## Remove result file if no hits

	grep -l "No Hits Found above Threshold" result_${sRNA}_${transcript} | xargs -r rm -f

	((transcript_counter=transcript_counter+1))

done

rm temp_${sRNA}.fasta

echo "Finished on $(date)."

