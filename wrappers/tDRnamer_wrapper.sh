#!/bin/bash
#SBATCH --job-name=tDRnamer_wrapper
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=tDRnamer_wrapper_output.txt
#SBATCH --error=tDRnamer_wrapper_output.txt

## Description
## Generates tRF fasta/csv in R, names the tRFs with tDRnamer, then joins the
## tDRnamer names back into the csv in R. Bridges the R <-> bash handoff so the
## whole thing runs in one shot.

## Requirements
## R (with Biostrings, glue, stringr, dplyr, readr, seqinr, tidyverse)
## tDRnamer

## Inputs
## $1 : scripts root dir (e.g. /grid/schorn/home/mpeacey/scripts/tRF_target_prediction)
## $2 : genome (e.g. mm10)
## $3 : tRF type (tRF3a or tRF3b)
## $4 : tDRnamer reference db prefix
##      (e.g. /grid/schorn/home/mpeacey/software/tDRnamer/reference_database/mm10/mm10)

set -euo pipefail

SCRIPTS=$1
GENOME=$2
TRF_TYPE=$3
DB=$4

cd "${SCRIPTS}"

## 1. R: generate the tRF fasta + csv
## (generate_tRF_fasta.R currently hard-codes `genome` near the top; set it to
##  ${GENOME} there, or parameterise it to read commandArgs, before relying on
##  this step.)
Rscript miranda/generate_tRF_fasta.R

## 2. bash: name the tRFs with tDRnamer
mkdir -p import/tDRnamer
tDRnamer --mode seq \
  --seq "import/${GENOME}_${TRF_TYPE}.fasta" \
  --db  "${DB}" \
  -o    "import/tDRnamer/${GENOME}"

## 3. R: join the tDRnamer names back into the csv
Rscript miranda/add_tDRnamer_names.R "${GENOME}" "${TRF_TYPE}"

echo "Done. See import/${GENOME}_${TRF_TYPE}.csv (now with a tDRnamer column)."
