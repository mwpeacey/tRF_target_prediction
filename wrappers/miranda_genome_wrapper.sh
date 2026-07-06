#!/bin/bash
#SBATCH --job-name=miranda_wrapper
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=miranda_wrapper_output.txt
#SBATCH --error=miranda_wrapper_output.txt

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

	# tDRnamer headers look like:
	#   >tDR-55:76-Ala-AGC-1|tRF_1 Sprinzl_position: 55..76
	# Use the clean tRF_N id (after the pipe, before any space) as the internal
	# identifier: it is unique and safe for filenames / SLURM job names.
	grep '^>' "${SMALL_RNA_FASTA}" \
		| sed -E 's/^>.*\|(tRF_[0-9]+).*/\1/' \
		> sRNA_list.txt

fi

for sRNA in $(cat sRNA_list.txt); do

	cd ${OUTPUT_DIRECTORY}/${RUN_ID}

	echo "Submitting job for sRNA ${sRNA}..."

	mkdir -p ${sRNA}
	cd ${sRNA}

	# Pull this fragment's record and rebuild a clean, whitespace-free 2-line
	# fasta. The header is trimmed to its first token ("tDRname|tRF_N"), which
	# drops the "Sprinzl_position: ..." description so miRanda's Seq1 is exactly
	# "<tDRname>|<tRF_N>" and carries both ids through to the output. The pipe
	# match is anchored (|tRF_1 not |tRF_10) to avoid partial-id collisions.
	awk -v id="${sRNA}" '
		/^>/ { keep = ($0 ~ ("\\|" id "([[:space:]]|$)"))
		       if (keep) { split($0, a, /[[:space:]]/); print a[1]; next } }
		keep { print }
	' "${SMALL_RNA_FASTA}" | head -2 > temp_${sRNA}.fasta
	
	sbatch --job-name=${RUN_ID}_${sRNA}_miranda \
    	--output=${OUTPUT_DIRECTORY}/${RUN_ID}/${sRNA}/${RUN_ID}_${sRNA}_miranda_output.txt \
    	--error=${OUTPUT_DIRECTORY}/${RUN_ID}/${sRNA}/${RUN_ID}_${sRNA}_miranda_output.txt \
    	${SCRIPTS}/miranda/miranda_genome.sh \
    	${GENOME_FASTA} ${RUN_ID} \
    	${sRNA} ${RUN_MODE} ${OUTPUT_DIRECTORY}	

done
