#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=32G
#$ -o permutation_analysis_universe_output.txt
#$ -e permutation_analysis_universe_output.txt

## $1 : miRanda CSV (output of miranda_to_bed.sh)
## $2 : RepeatMasker CSV
## $3 : universe (.rds GRanges or .bed)
## $4 : output directory
## $5 : minimum alignment score (default 80)
## $6 : number of permutations (default 1000)

data_file=$1
rmsk_file=$2
universe_file=$3
output_directory=$4
cutoff=${5:-80}
ntimes=${6:-1000}

module load EBModules
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/R_scripts/permutation_analysis_universe.R \
	${data_file} ${rmsk_file} ${universe_file} ${output_directory} ${cutoff} ${ntimes}

mv permutation_analysis_universe_output.txt ${output_directory}
