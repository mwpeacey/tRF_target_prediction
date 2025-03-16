#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=8G

## Description
## Uses miRanda to predict tRF or miRNA targets in an input transcriptome. Designed for use with miranda_wrapper.sh.

## Requirements
## miranda v1.9

## Inputs
## $1 : directory containing fasta files, each of which represents a transcript to scan.
## $2 : ID for the run, used as the directory name in the output directory.
## $3 : sRNA name. Corresponds to one entry in a larger small RNA fasta file.
## $4 : determines the settings to use for miRanda. "miRNA" runs with standard seed-weighting,
##      whereas "tRF" runs with seed-weighting removed and the alignment score threshold adjusted
##      accordingly.
## $5 : output directory.

echo "New run started on $(date)"

TRANSCRIPTOME_DIRECTORY=$1
RUN_ID=$2
sRNA=$3
RUN_MODE=$4
OUTPUT_DIRECTORY=$5

echo "Output directory: ${OUTPUT_DIRECTORY}"
echo "Transcriptome directory: ${TRANSCRIPTOME_DIRECTORY}"
echo "Run ID: ${RUN_ID}"
echo "Small RNA: ${sRNA}"
echo "Run mode: ${RUN_MODE}"

cd ${OUTPUT_DIRECTORY}/${RUN_ID}

if [[ -f mouse_tRF3b_v2_${sRNA}_miranda_output.txt ]]; then

	echo "Output file found. Checking if finished..."

	if grep -q "Finished" mouse_tRF3b_v2_${sRNA}_miranda_output.txt; then

        	echo "Already finished. Exiting..."

	else

    		echo "Not finished. Iterating through transcripts..."

        	cd ${TRANSCRIPTOME_DIRECTORY}

        	transcript_counter=1
        	transcript_number=$(ls -1 | wc -l)

        	for transcript in *.fasta; do

                	TRANSCRIPT_NAME="${transcript%.fasta}"

                	if grep -q ${TRANSCRIPT_NAME} mouse_tRF3b_v2_${sRNA}_miranda_output.txt; then

                        	echo "${TRANSCRIPT_NAME} finished, skipping..."

                	else

				echo "Scanning transcript #${transcript_counter} of ${transcript_number}, ID ${TRANSCRIPT_NAME}"

                		## Run miranda using that temporary file

                		cd ${OUTPUT_DIRECTORY}/${RUN_ID}/${sRNA}

                		if [ ${RUN_MODE} == 'tRF' ]; then

                        		miranda temp_${sRNA}.fasta \
                                		${TRANSCRIPTOME_DIRECTORY}/${transcript} \
                                		-out result_${sRNA}_${TRANSCRIPT_NAME} \
                                		-sc 75.0 -en -20.0 -scale 1.0 -loose

                		fi

                		if [ ${RUN_MODE} == 'miRNA' ]; then

                        		miranda temp_${sRNA}.fasta \
                                		${TRANSCRIPTOME_DIRECTORY}/${transcript} \
                                		-out result_${sRNA}_${TRANSCRIPT_NAME} \
                                	-sc 120

                		fi

                		## If hits, move the summary line to a summary file

                		grep -A 1 "Scores for this hit:" result_${sRNA}_${TRANSCRIPT_NAME} | sort | grep '>' > temp_summary_${transcript_counter}

                		if [ $transcript_counter = 1 ]; then

                        		mv temp_summary_1 ${sRNA}_summary

                		else

                    			cat temp_summary_${transcript_counter} >> ${sRNA}_summary
                        		rm temp_summary_${transcript_counter}

                		fi

                		## Remove result file if no hits

                		grep -l "No Hits Found above Threshold" result_${sRNA}_${TRANSCRIPT_NAME} | xargs -r rm -f

			fi

                	((transcript_counter=transcript_counter+1))

        	done

	fi

else

	echo "No previous output found. Starting from scratch..."

	cd ${TRANSCRIPTOME_DIRECTORY}

	echo "Iterating through transcripts..."

	transcript_counter=1
	transcript_number=$(ls -1 | wc -l)

	for transcript in *.fasta; do
	
		TRANSCRIPT_NAME="${transcript%.fasta}"

		echo "Scanning transcript #${transcript_counter} of ${transcript_number}, ID ${TRANSCRIPT_NAME}"

		cd ${OUTPUT_DIRECTORY}/${RUN_ID}/${sRNA}

		if [ ${RUN_MODE} == 'tRF' ]; then

			miranda temp_${sRNA}.fasta \
				${TRANSCRIPTOME_DIRECTORY}/${transcript} \
				-out result_${sRNA}_${TRANSCRIPT_NAME} \
				-sc 75.0 -en -20.0 -scale 1.0 -loose

		fi

		if [ ${RUN_MODE} == 'miRNA' ]; then

                	miranda temp_${sRNA}.fasta \
                        	${TRANSCRIPTOME_DIRECTORY}/${transcript} \
                        	-out result_${sRNA}_${TRANSCRIPT_NAME} \
                        	-sc 120
		fi

		## If hits, move the summary line to a summary file

		grep -A 1 "Scores for this hit:" result_${sRNA}_${TRANSCRIPT_NAME} | sort | grep '>' > temp_summary_${transcript_counter}

		if [ $transcript_counter = 1 ]; then

	               	mv temp_summary_1 ${sRNA}_summary

	        else

	            	cat temp_summary_${transcript_counter} >> ${sRNA}_summary
			rm temp_summary_${transcript_counter}

	        fi

		## Remove result file if no hits

		grep -l "No Hits Found above Threshold" result_${sRNA}_${TRANSCRIPT_NAME} | xargs -r rm -f

		((transcript_counter=transcript_counter+1))

	done

fi

rm ${OUTPUT_DIRECTORY}/${RUN_ID}/${sRNA}/temp_${sRNA}.fasta

echo "Finished on $(date)."

