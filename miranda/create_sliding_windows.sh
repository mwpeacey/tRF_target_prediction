#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G
#$ -o create_sliding_windows_output.txt
#$ -e create_sliding_windows_output.txt

## Description
## Takes a genome fasta and splits it into sliding windows, for input into miRanda.

## Requirements
## samtools (1.20)
## bedtools (2.31.1)
## seqkit (2.10.0)

## Inputs
## 1 : Scripts directory
## 2 : Genome fasta
## 3 : Output directory

echo "New run started on $(date)"

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <scripts_dir> <genome.fa> <out_dir>" >&2
  exit 1
fi

SCRIPTS="$1"
GENOME_FA="$2"
OUTDIR="$3"
WINDOW=10000
OVERLAP=50
STEP=$(( WINDOW - OVERLAP ))

BASE=$(basename "${GENOME_FA%.*}")

OUT_FA="${OUTDIR}/${BASE}_${WINDOW}bp_windows.fa"

echo "[$(date)] Starting window split of $GENOME_FA → $OUT_FA"

samtools faidx "$GENOME_FA"

bedtools makewindows \
    -g "${GENOME_FA}.fai" \
    -w "$WINDOW" \
    -s "$STEP" \
  | awk -v OFS="\t" '{ printf("%s\t%s\t%s\t%s_%d\n",$1,$2,$3,$1,$2+1) }' \
  | bedtools getfasta \
      -fi "$GENOME_FA" \
      -bed - \
      -name \
      -fo "$OUT_FA"

seqkit seq -r -p "$OUT_FA" \
  > "${OUTDIR}/${BASE}_${WINDOW}bp_windows_minus.fa"

echo "Finished on $(date)."

mv ${SCRIPTS}/miranda/create_sliding_windows_output.txt ${OUTPUT_DIRECTORY}/create_sliding_windows_output.txt
