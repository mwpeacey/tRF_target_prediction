################################################################################
# tRF_bowtie_alignment.R
################################################################################

library(tidyverse)
library(GenomicRanges)

# Import repeat annotation

repeats = rtracklayer::readGFF('import/annotation_tables/dm6_rmsk_TE.gtf') %>%
  dplyr::filter(gene_id == 'MDG1_I-int' | gene_id == 'MDG1_LTR')

repeats = rtracklayer::readGFF('import/annotation_tables/danRer11_rmsk_TE.gtf') %>%
  dplyr::filter(gene_id == 'GypsyDR1')

subject = makeGRangesFromDataFrame(repeats, 
                                   keep.extra.columns = T)

# Import bed input 

alignment = read.csv(file = 'import/bowtie/dm6_tRF3a.bed', header = F, sep = '\t') %>%
  dplyr::rename('seqnames' = 'V1', 
                'start' = 'V2', 
                'end' = 'V3',
                'info' = 'V4',
                'score' = 'V5',
                'strand' = 'V6') %>%
  mutate(strand = case_when(strand == '+' ~ '-', strand == '-' ~ '+'))

alignment = read.csv(file = 'import/bowtie/danRer11_tRF3a.bed', header = F, sep = '\t') %>%
  dplyr::rename('seqnames' = 'V1', 
                'start' = 'V2', 
                'end' = 'V3',
                'info' = 'V4',
                'score' = 'V5',
                'strand' = 'V6') %>%
  mutate(strand = case_when(strand == '+' ~ '-', strand == '-' ~ '+'))

# tRNA info

tRNAs_info = read.csv('import/mm10_tRF3a.csv') %>%
  dplyr::rename(info = unique_name) 

# Overlap

alignment = merge(alignment, tRNAs_info)

query = GRanges(alignment)

overlap = findOverlaps(query = query, 
                       subject = subject,
                       type = 'any')

overlap = as.data.frame(overlap)

tRF = vector()

for (i in 1:nrow(repeats)){
  
  tRNAs = dplyr::filter(overlap, subjectHits == i)$queryHits
  
  tRF[i] = paste(alignment[tRNAs, ]$info, collapse = ', ')
  
}

repeats$tRF = tRF

repeats = mutate(repeats, tRF_name = case_when(tRF == 'tRF_21' ~ 'Arg-TCG',
                                               T ~ 'None'))

repeats = mutate(repeats, tRF_name = case_when(tRF == 'tRF_2, tRF_5, tRF_6' | 
                                                 tRF == 'tRF_5, tRF_6' | 
                                                 tRF == 'tRF_113, tRF_2, tRF_5, tRF_6' |
                                                 tRF == 'tRF_1' ~ 'Ala-AGC',
                                               tRF == 'tRF_12, tRF_9' |
                                                 tRF == 'tRF_91' |
                                                 tRF == 'tRF_27' |
                                                 tRF == 'tRF_63, tRF_64' ~ 'other',
                                               tRF == 'tRF_125' ~ 'iMet-CAT',
                                               T ~ 'None'))

input = dplyr::filter(repeats, tRF_name != 'None') %>%
  group_by(tRF_name) %>%
  summarize(count = n())

ggplot(dplyr::filter(input), aes(x="", y=count, fill=tRF_name)) +
  geom_bar(stat="identity", width=1, color = 'black') +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position = 'top')





