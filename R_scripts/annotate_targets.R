################################################################################
# annotate_targets.R
# Annotates genomic features of hit sites, including hit transcript biotype,
# the position of the hit relative to the CDS (for protein coding genes), and
# overlap with long terminal repeats.
################################################################################

library(tidyverse)
library(glue)
library(GenomicFeatures)

################################################################################
## Import data
################################################################################

# Load data

data = read.csv(file = 'import/miranda/miranda_output_70.csv', header = TRUE) %>%
  dplyr::rename(transcript_start = start, transcript_end = end, start = genome_start, end = genome_end) %>%
  dplyr::filter(alignment_score >= 80)

# Add gene information 

transcript_to_gene = as.data.frame(rtracklayer::import('import/transcriptome_assembly/stringtie_merged_filtered.gtf')) %>%
  dplyr::filter(type == 'transcript') %>%
  dplyr::rename(target_id = 'transcript_id') %>%
  dplyr::select(c('target_id', 'ref_gene_id', 'gene_name'))

data = merge(data, transcript_to_gene, by = 'target_id') %>%
  separate(ref_gene_id, into = c('ref_gene_id'))

# Add tRF information

tRF_infomation = read.csv('import/mm10_tRF3b.csv') %>%
  dplyr::rename(tRF = 'unique_name')

data$tRF = stringr::str_remove(data$tRF, '>')

data = merge(data, tRF_infomation, by = 'tRF')

################################################################################
## Annotate hit overlap with repeats
################################################################################

# Create GRanges object for LTRs

LTR = rtracklayer::readGFF(file = 'import/annotation_tables/mm10_rmsk_TE.gtf') %>%
  dplyr::filter(class_id == 'LTR')

subject = makeGRangesFromDataFrame(LTR, keep.extra.columns = TRUE)

# Create GRanges object for query data

query = GRanges(data)

# Find overlaps between query and LTRs

LTR_overlap = findOverlaps(query = query, subject = subject, type = 'any', ignore.strand = T)
LTR_overlap = as.data.frame(LTR_overlap)

# Permutation analysis
library(regioneR)

GTF = makeTxDbFromGFF(file = 'import/transcriptome_assembly/stringtie_merged_filtered.gtf')
transcriptome_background = reduce(exons(GTF))
seqlevelsStyle(transcriptome_background) = 'UCSC'

perm_results = overlapPermTest(
  A = subject,
  B = query,
  ntimes = 20,
  genome = 'mm10', 
  universe=transcriptome_background ,
  verbose = T)

