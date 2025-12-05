#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=32G
#$ -o permutation_analysis_output.txt
#$ -e permutation_analysis_output.txt

tRF_file=$1
miRNA_file=$2
five_UTR=$3
output_directory=$4

module load EBModules
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/R_scripts/permutation_analysis_miRNA.R ${tRF_file} ${miRNA_file} ${five_UTR} ${output_directory}

mv permutation_analysis_output.txt ${output_directory}
