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
  dplyr::filter(alignment_score >= 70)

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

data = merge(data, tRF_infomation, by = 'tRF') %>%
  separate(col = source_tRNAs, sep = '-', into = c('A', 'B', 'C', 'D'), remove = F) %>%
  unite(tRNA_anticodon, c('B', 'C'), remove = F, sep = '-') %>%
  dplyr::select(-c('A', 'B', 'C', 'D'))

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

# Add LTR indicator to data

data = data %>%
  dplyr::mutate(LTR = dplyr::row_number() %in% LTR_overlap$queryHits)

# Initialize LTR family vector
LTR_family = vector(length = nrow(data))

# Assign LTR family based on overlaps
for (row in 1:nrow(data)) {
  if (data[row, 'LTR']) {
    subject_row = LTR_overlap$subjectHits[LTR_overlap$queryHits == row][1]
    LTR_family[row] = subject[subject_row,]$family_id
  } else {
    LTR_family[row] = NA
  }
}

data$LTR_family = LTR_family

################################################################################
## Heatmap
################################################################################

library(pheatmap)

input = dplyr::filter(data, LTR == T) %>%
  group_by(tRNA_anticodon, LTR_family) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = 'tRNA_anticodon', values_from = 'n') %>%
  as.data.frame()

input[is.na(input)] = 0

row.names(input) = input$LTR_family

input = dplyr::select(input, -c('LTR_family'))

input = log10(input+0.01)

pheatmap(input, cluster_col = T, cluster_row = T)

