#!/bin/bash
#SBATCH --job-name=miranda_to_bed
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=miranda_to_bed_output.txt
#SBATCH --error=miranda_to_bed_output.txt

summary_directory=$1
output_directory=$2
score_cutoff=$3

module load EBModules
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/R_scripts/miranda_to_bed.R ${summary_directory} ${output_directory} ${score_cutoff}

mv miranda_to_bed_output.txt ${output_directory}
