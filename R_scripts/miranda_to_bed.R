################################################################################
# miranda_to_bed.R
# Combines miRanda summary files into a single output with genomic coordinates
# of hit sites. Writes to .csv and .bed format.
################################################################################

library(tidyverse)
library(glue)

args = commandArgs(TRUE)

summary_directory = args[1]
output_directory = args[2]
score_cutoff = as.numeric(args[3])

################################################################################
## Import data
################################################################################

filenames = list.files(summary_directory,
                       pattern="summary*", 
                       full.names=TRUE)

counter = 1
for (file in filenames){
  
  temp = read.csv(file = file, header = F, sep = '\t', colClasses = 'character')
  
  if (counter == 1){
    
    data = temp
    
  }
  
  else{
    
    data = bind_rows(data, temp)
    
  }
  
  counter = counter + 1
  
}

remove(temp)

################################################################################
## Reformat and filter
################################################################################

# Rename columns and enforce a post-run score cutoff.

data = data[, -which(names(data) %in% c("V9", "V10"))]

miranda_output = data %>%
  dplyr::rename('tRF' = V1, 
                'coordinates' = V2, 
                'alignment_score' = V3,
                'energy' = V4,
                'Z_score' = V5,
                'miRNA_position' = V6,
                'target_position' = V7,
                'alignment_length' = V8,
                'strand' = V11) %>% 
  dplyr::filter(alignment_score >= score_cutoff)

# Find start and end positions of the hit in the window

miranda_output = separate(miranda_output, 'target_position', c('start', 'end'))

miranda_output$start = as.numeric(miranda_output$start)
miranda_output$end = as.numeric(miranda_output$end)

# Find start and end coordinates of the window

miranda_output = separate(miranda_output, 'coordinates', c('A', 'B'), sep = '::') %>%
  separate('B', c('seqnames', 'C'), sep = ':') %>%
  separate('C', c('window_start', 'window_end'), sep = '-')

miranda_output$window_start = as.numeric(miranda_output$window_start)
miranda_output$window_end = as.numeric(miranda_output$window_end)

# Determine genomic coordinates of hit site

miranda_output = miranda_output %>%
  mutate(
    genomic_start = if_else(
      strand == "+",
      window_start + start,
      window_end   - end   + 1
    ),
    genomic_end   = if_else(
      strand == "+",
      window_start + end,
      window_end   - start + 1
    )
  )


# Remove duplicate hits that occur within window overlaps

miranda_output = distinct(miranda_output, tRF, seqnames, genomic_start, genomic_end, strand, .keep_all = TRUE)

# Cleanup 

miranda_output = dplyr::select(miranda_output, c('tRF', 
                                                 'seqnames', 
                                                 'alignment_score', 
                                                 'energy', 
                                                 'genomic_start', 
                                                 'genomic_end', 
                                                 'alignment_length', 
                                                 'strand'))

################################################################################
# Export
################################################################################

## Write csv 

write_csv(miranda_output, file = glue('{output_directory}/miranda_output_{score_cutoff}.csv'))

## Write bed

miranda_bed = miranda_output %>%
  dplyr::select(c('seqnames', 'genomic_start', 'genomic_end', 'tRF', 'strand', 'alignment_score')) %>%
  dplyr::rename(chrom = 'seqnames', chromStart = 'genomic_start', chromEnd = 'genomic_end', name = 'tRF', score = 'alignment_score')

rtracklayer::export.bed(con = glue("{output_directory}/miranda_output_{score_cutoff}.bed"), 
                        object = GenomicRanges::makeGRangesFromDataFrame(miranda_bed, keep.extra.columns = T),
                        ignore.strand = F)
