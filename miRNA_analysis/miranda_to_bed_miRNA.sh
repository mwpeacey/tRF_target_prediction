#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=32G
#$ -o miranda_to_bed_output.txt
#$ -e miranda_to_bed_output.txt

summary_directory=$1
output_directory=$2
score_cutoff=$3
utr5_annotation_rds=$4

module load EBModules
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/R_scripts/miranda_to_bed_miRNA.R ${summary_directory} ${output_directory} ${score_cutoff} ${utr5_annotation_rds}

mv miranda_to_bed_output.txt ${output_directory}
