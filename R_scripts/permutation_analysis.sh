#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=32G

data_file=$1
rmsk_file=$2
gtf_file=$3
min_cutoff=$4
output_directory=$5

module load EBModules
module load R/4.0.4-fosscuda-2020a

Rscript /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/R_scripts/permutation_analysis.R ${data_file} ${rmsk_file} ${gtf_file} ${min_cutoff} ${output_directory}
