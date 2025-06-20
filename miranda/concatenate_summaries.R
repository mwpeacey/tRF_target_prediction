library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(glue)

import_directory = '~/tRF_target_prediction/import/human_tRF_rules'

################################################################################
## Import data
################################################################################

filenames = list.files(glue("{import_directory}"),
                       pattern="*summary", 
                       full.names=TRUE)

counter = 1
for (file in filenames){
  
  temp = read.csv(file = file, header = F, sep = '\t')
  
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

score_cutoff = 70

miranda_output = data %>%
  dplyr::rename('tRF' = 'V1', 
                'target_id' = 'V2', 
                'alignment_score' = 'V3',
                'energy' = 'V4',
                'Z-score' = 'V5',
                'miRNA_position' = 'V6',
                'target_position' = 'V7',
                'alignment_length' = 'V8') %>%
  dplyr::filter(alignment_score >= score_cutoff)

## Find start and end positions of the hit on the transcript

miranda_output = separate(miranda_output, 'target_position', c('start', 'end'))

## Temporary fix to remove weird target_id entries

#miranda_output = miranda_output[grepl('ENSMUS', miranda_output$target_id) | grepl('MSTRG', miranda_output$target_id), ]

################################################################################
# Adds genomic coordinates of hit sites
################################################################################

# Import GTF annotation as TxDb object

transcriptome_gtf = list.files(glue("{import_directory}"),
                       pattern="*.gtf", 
                       full.names=TRUE)

transcriptome_TxDb = makeTxDbFromGFF(file=transcriptome_gtf,
                                     organism='Homo sapiens')

exons_by_transcript = exonsBy(transcriptome_TxDb, "tx", use.names=TRUE)

# Filters miranda output to ensure all hit transcripts are in the GTF file.
# They may be absent if e.g. the GTF file has been filtered after the
# miranda run.

miranda_output = miranda_output[miranda_output$target_id %in% names(exons_by_transcript), ]

# Convert miranda hit coordinates to GRanges object

rng_tx = IRanges(start = as.numeric(miranda_output$start), 
                 end = as.numeric(miranda_output$end),
                 names = miranda_output$target_id)

miranda_GRanges = GRanges(seqnames=miranda_output$target_id, ranges=rng_tx, strand = "*")

# Map transcript coordinates to genome

genome_coordinates = mapFromTranscripts(miranda_GRanges, exons_by_transcript)
remove(miranda_GRanges)

# Combine with original miranda output

df = as.data.frame(genome_coordinates, 
                   row.names = genome_coordinates$xHits)

df$target_id = names(genome_coordinates)

miranda_output$genome_start = df$start
miranda_output$genome_end = df$end
miranda_output$seqnames = df$seqnames
miranda_output$strand = df$strand

remove(df)


# Bed files

bed_files = list.files(summary_directory,
                       pattern="*.bed", 
                       full.names=TRUE)

counter = 1
for (file in bed_files){
  
  temp = as.data.frame(import.bed(file))
  
  if (counter == 1){
    
    output = temp
    
  }
  
  else{
    
    output = bind_rows(output, temp)
    
  }
  
  counter = counter + 1
  
}

export.bed(con = glue("{summary_directory}/miranda_output.bed"), 
           object = GenomicRanges::makeGRangesFromDataFrame(output, keep.extra.columns = T),
           ignore.strand = F)

# csv files

csv_files = list.files(summary_directory,
                       pattern="*.csv", 
                       full.names=TRUE)

counter = 1
for (file in csv_files){
  
  temp = read_csv(file)
  
  if (counter == 1){
    
    output = temp
    
  }
  
  else{
    
    output = bind_rows(output, temp)
    
  }
  
  counter = counter + 1
  
}

write_csv(output, file = glue::glue("{summary_directory}/miranda_output.csv"))
