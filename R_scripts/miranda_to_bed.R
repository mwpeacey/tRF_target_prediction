################################################################################
# miranda_to_bed.R
# Combines miRanda summary files into a single output with genomic coordinates
# of hit sites. Writes to .bed format.
################################################################################

library(tidyverse)
library(GenomicFeatures)
library(glue)

################################################################################
## Import data
################################################################################

filenames = list.files(glue("import/miranda"),
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

# For each transcript, takes the hit with the highest alignment score (deprecated)

# miranda_output = miranda_output %>%
#  group_by(tRF, target_id) %>%
#  dplyr::slice(which.max(alignment_score))

# Find start and end positions of the hit on the transcript

miranda_output = separate(miranda_output, 'target_position', c('start', 'end'))

# Temporary fix to remove weird target_id entries

miranda_output = miranda_output[grepl('ENSMUS', miranda_output$target_id) | grepl('MSTRG', miranda_output$target_id), ]

################################################################################
# Adds genomic coordinates of hit sites
################################################################################

# Import GTF annotation as TxDb object

transcriptome_TxDb = makeTxDbFromGFF(file="import/transcriptome_assembly/stringtie_merged_filtered.gtf",
                                     dataSource='Stringtie assembled transcriptome combined with GENCODE vM23.',
                                     organism='Mus musculus')

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

# Filter to remove exon junction spanning hits

miranda_output = miranda_output %>%
  mutate(hit_length = as.numeric(genome_end) - as.numeric(genome_start) + 1) %>%
  filter(hit_length <= 40)

################################################################################
# Export
################################################################################

## Write bed

miranda_bed = miranda_output %>%
  dplyr::select(c('seqnames', 'genome_start', 'genome_end', 'tRF', 'strand')) %>%
  dplyr::rename(chrom = 'seqnames', chromStart = 'genome_start', chromEnd = 'genome_end', name = 'tRF') 

rtracklayer::export.bed(con = glue("import/miranda/miranda_output_{score_cutoff}.bed"), 
                        object = GenomicRanges::makeGRangesFromDataFrame(miranda_bed, keep.extra.columns = T),
                        ignore.strand = F)
