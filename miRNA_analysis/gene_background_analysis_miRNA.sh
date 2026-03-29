#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=32G
#$ -o gene_background_analysis_output.txt
#$ -e gene_background_analysis_output.txt

tRF_file=$1
miRNA_file=$2
five_UTR=$3
rmsk_file=$4
transcript_gtf=$5
target_gene=$6
output_directory=$7
universe_mode=${8:-utr_only}

module load EBModules
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /grid/schorn/home/mpeacey/scripts/tRF_target_prediction/miRNA_analysis/gene_background_analysis_miRNA.R \
  "${tRF_file}" "${miRNA_file}" "${five_UTR}" "${rmsk_file}" \
  "${transcript_gtf}" "${target_gene}" "${output_directory}" "${universe_mode}"

mv gene_background_analysis_output.txt "${output_directory}"
