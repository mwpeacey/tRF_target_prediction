#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=32G

data_file=$1
rmsk_file=$2
gtf_file=$3
min_cutoff=$4
output_directory=$5

#module load EBModules
#module load R/4.0.4-fosscuda-2020a
#module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

source /grid/schorn/home/mpeacey/.bashrc
conda activate r_genomic

Rscript /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/R_scripts/permutation_analysis.R ${data_file} ${rmsk_file} ${gtf_file} ${min_cutoff} ${output_directory}
