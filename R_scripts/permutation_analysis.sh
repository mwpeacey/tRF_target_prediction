#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=32G
#$ -o permutation_analysis_output.txt
#$ -e permutation_analysis_output.txt

data_file=$1
rmsk_file=$2
#gtf_file=$3
min_cutoff=$3
output_directory=$4

module load EBModules
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/R_scripts/permutation_analysis.R ${data_file} ${rmsk_file} ${min_cutoff} ${output_directory}

mv permutation_analysis_output.txt ${output_directory}
