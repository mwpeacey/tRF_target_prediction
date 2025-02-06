## Generates a fasta file of unique 18nt (tRF3a) and 22nt (tRF3b) tRNA fragments.
## Input is a fasta of mature tRNA sequences (minus CCA) downloaded from GtRNAdb.

library(Biostrings)
library(glue)
library(stringr)
library(tidyverse)

genome = 'hg38'

tRNA_fasta = Biostrings::readRNAStringSet(glue("~/data/{genome}-mature-tRNAs.fa"),
                                     format = 'fasta')

tRNA_fasta = DNAStringSet(tRNA_fasta)

tRNA_df = data.frame(sequence = stringr::str_c(as.character(tRNA_fasta), 'CCA'), 
                  name = names(tRNA_fasta))

for (i in 1:nrow(tRNA_df)){
  
  tRNA_df[i, 'name'] = str_split(tRNA_df[i, 'name'], ' ')[[1]][1]
  
  tRNA_df[i, 'name'] = str_split(tRNA_df[i, 'name'], '_')[[1]][3]
  
  tRNA_df[i, 'tRF3a'] = str_sub(tRNA_df[i, 'sequence'], start = -18, end = -1)
  
  tRNA_df[i, 'tRF3b'] = str_sub(tRNA_df[i, 'sequence'], start = -22, end = -1)
  
}

###############################################################################
## Generate tRNA fasta
###############################################################################

seqinr::write.fasta(sequences = as.list(tRNA_df$sequence),
                    names = as.list(tRNA_df$name),
                    file.out = glue::glue('~/tRF_targets_new/{genome}_tRNA.fasta'),
                    as.string = T)

###############################################################################
## Generate tRF3a fasta
###############################################################################

tRF3_seq = vector()
tRF3a_name = vector()

for (i in 1:length(unique(tRNA_df$tRF3a))){
  
  tRF3_seq[i] = unique(tRNA_df$tRF3a)[i]
  
  tRF3a_name[i]= paste(c(filter(tRNA_df, tRF3a == tRF3_seq[i])$name), collapse="_")
  
}

df = data.frame(sequence = tRF3_seq,
                source_tRNAs = tRF3a_name)

unique_name = vector()

for (i in 1:nrow(df)){
  
  unique_name[i] = glue::glue("tRF_{i}")
  
}

df$unique_name = unique_name

# Export

seqinr::write.fasta(sequences = as.list(df$sequence),
                    names = as.list(df$unique_name),
                    file.out = glue::glue('~/tRF_targets_new/{genome}_tRF3a.fasta'),
                    as.string = T)

write_csv(df, file = glue('~/tRF_targets_new/{genome}_tRF3a.csv'))

###############################################################################
## Generate tRF3b fasta
###############################################################################

tRF3_seq = vector()
tRF3b_name = vector()

for (i in 1:length(unique(tRNA_df$tRF3b))){
  
  tRF3_seq[i] = unique(tRNA_df$tRF3b)[i]
  
  tRF3b_name[i]= paste(c(dplyr::filter(tRNA_df, tRF3b == tRF3_seq[i])$name), collapse="_")
  
}

df = data.frame(sequence = tRF3_seq,
                source_tRNAs = tRF3b_name)

unique_name = vector()

for (i in 1:nrow(df)){
  
  unique_name[i] = glue::glue("tRF_{i}")
  
}

df$unique_name = unique_name

# Export

seqinr::write.fasta(sequences = as.list(df$sequence),
                    names = as.list(df$unique_name),
                    file.out = glue::glue('~/tRF_targets_new/{genome}_tRF3b.fasta'),
                    as.string = T)

write_csv(df, file = glue('~/tRF_targets_new/{genome}_tRF3b.csv'))