################################################################################
# annotate_targets.R
# Annotates genomic features of hit sites, including hit transcript biotype,
# the position of the hit relative to the CDS (for protein coding genes), and
# overlap with long terminal repeats.
################################################################################

library(tidyverse)
library(glue)
library(GenomicFeatures)
library(pheatmap)
library(AnnotationHub)

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

# Filter for GENCODE transcripts only

#data = data[grepl('ENSMUS', data$target_id), ]

################################################################################
## Annotate hit overlap with transcripts
################################################################################

# Initialize annotation database
ah = AnnotationHub()
edb = ah[["AH75036"]]

query = GenomicRanges::GRanges(data)
GenomeInfoDb::seqlevelsStyle(query) = 'NCBI'

# Does the hit overlap a transcript?

subject = transcripts(edb)
GenomeInfoDb::seqlevelsStyle(subject) = 'NCBI'

transcript_overlap = findOverlaps(query = query, 
                                  subject = subject,
                                  type = 'any')

transcript_overlap = as.data.frame(transcript_overlap)

# Does it overlap an exon?

subject = exonsBy(edb, by = 'tx')

exon_overlap = findOverlaps(query = query, 
                            subject = subject,
                            type = 'any')

exon_overlap = as.data.frame(exon_overlap)

# Does it overlap a transcript with a CDS?

transcripts = transcripts(edb)

subject = exonsBy(edb, by = 'tx')

subject = subject[names(subject) %in% transcripts[transcripts$tx_biotype == 'protein_coding']$tx_id]

coding_overlap = findOverlaps(query = query, 
                              subject = subject,
                              type = 'any')

coding_overlap = as.data.frame(coding_overlap)

# Does it overlap the CDS itself?

subject = cdsBy(edb, by = 'tx')

CDS_overlap = findOverlaps(query = query, 
                           subject = subject,
                           type = 'any')

CDS_overlap = as.data.frame(CDS_overlap)

# Does it overlap the 5' UTR?

subject = fiveUTRsByTranscript(edb)

five_overlap = findOverlaps(query = query, 
                            subject = subject,
                            type = 'any')

five_overlap = as.data.frame(five_overlap)

# Does it overlap the 3' UTR?

subject = threeUTRsByTranscript(edb)

three_overlap = findOverlaps(query = query, 
                             subject = subject,
                             type = 'any')

three_overlap = as.data.frame(three_overlap)

## Annotate accordingly

data = data %>%
  mutate(intergenic = !(seq_len(n()) %in% transcript_overlap$queryHits)) %>%
  mutate(intronic = !intergenic & !(seq_len(n()) %in% exon_overlap$queryHits)) %>%
  mutate(non_coding = !intergenic & !intronic & !(seq_len(n()) %in% coding_overlap$queryHits)) %>%
  mutate(CDS = !intergenic & !intronic & !non_coding & (seq_len(n()) %in% CDS_overlap$queryHits)) %>%
  mutate(five_prime_UTR = !intergenic & !intronic & !non_coding & !CDS & (seq_len(n()) %in% five_overlap$queryHits)) %>%
  mutate(three_prime_UTR = !intergenic & !intronic & !non_coding & !CDS & !five_prime_UTR & (seq_len(n()) %in% three_overlap$queryHits)) %>%
  mutate(location = case_when(
    intergenic ~ "None",
    intronic ~ "Intron",
    non_coding ~ "Exon - Non-coding",
    CDS ~ "Exon - CDS",
    five_prime_UTR ~ "Exon - 5 'UTR",
    three_prime_UTR ~ "Exon - 3 'UTR",
    TRUE ~ "other"
  ))

# Does it overlap a transcript with a SCAN domain

mouse_scan_genes = read_csv('import/annotation_tables/mouse_scan_genes.csv')

subject = exonsBy(edb, by = 'tx')

subject = subject[names(subject) %in% transcripts[transcripts$gene_id %in% mouse_scan_genes$gene_id]$tx_id]

SCAN_overlap = findOverlaps(query = query, 
                            subject = subject,
                            type = 'any')

SCAN_overlap = as.data.frame(SCAN_overlap)

data = data %>%
  mutate(SCAN = seq_len(n()) %in% SCAN_overlap$queryHits)

################################################################################
## Annotate hit overlap with repeats
################################################################################

# Create GRanges object for LTRs

LTR = rtracklayer::readGFF(file = 'import/annotation_tables/mm10_rmsk_TE.gtf') %>%
  dplyr::filter(class_id == 'LTR')

subject = makeGRangesFromDataFrame(LTR, keep.extra.columns = TRUE)
end(subject[strand(subject) == '+']) = end(subject[strand(subject) == '+']) + 200
start(subject[strand(subject) == '-']) = start(subject[strand(subject) == '-']) - 200

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
## Collapse overlapping sites
################################################################################

unique_data = data %>%
  arrange(desc(alignment_score), energy) %>%
  group_by(seqnames, start, end, strand) %>%
  slice(1) %>%
  ungroup()

################################################################################
## Export
################################################################################

write_csv(data, file = 'import/miranda_output_annotated.csv')
write_csv(unique_data, file = 'import/miranda_output_unique_annotated.csv')

################################################################################
## Plots
################################################################################

library(patchwork)
library(tidyverse)

# Distribution of alignment scores

input = dplyr::filter(unique_data, alignment_score >= 70) %>%
  group_by(location) %>%
  summarize(n = n())

size = 1

permutation_analysis = read_csv(file = 'import/LTR_enrichment_by_cutoff.csv')
permutation_analysis$padj = p.adjust(permutation_analysis$p_value, method = "BH")

cutoff = 80

n = nrow(dplyr::filter(unique_data, alignment_score >= cutoff, LTR == T))

# Histogram
hist_plot = ggplot(unique_data, aes(x = alignment_score)) +
  geom_histogram(fill = '#7AAFD3', color = 'black', binwidth = 1, boundary = 0) +
  geom_vline(xintercept = 80, linetype = 'dashed') +
  coord_cartesian(xlim = c(74, 90)) +
  scale_y_continuous(limits = c(0, 25000),
                     expand = expansion(mult = c(0, .1))) +
  xlab('Alignment score') +
  ylab('Count') +
  annotate("text", x = 85, y = 20000, 
           label = paste0("n = ", n, "\n(score ≥ 80)"), 
           hjust = 1, size = 5) +
  theme_bw() + theme(
    axis.text.x = element_text(size = 14, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(size = 16, margin = margin(r = 7)),
    axis.line = element_line(size = size),
    axis.ticks = element_line(size = size, color = 'black'),
    panel.border = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.background = element_rect(color = "white", fill = "white", size = 1.5, linetype = "solid"), 
    plot.title = element_text(hjust = 0.5, size = 10)
  )

# Z-score line
z_plot = ggplot(permutation_analysis, aes(x = cutoff, y = z_score)) +
  geom_line(color = 'firebrick', size = 1) +
  geom_vline(xintercept = 80, linetype = 'dashed') +
  coord_cartesian(xlim = c(74, 90)) +
  xlab('Alignment score') +
  ylab('Z-score (LTR enrichment)') +
  theme_bw() + theme(
    axis.text.x = element_text(size = 14, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title.x = element_text(size = 16, margin = margin(t = 7)),
    axis.title.y.left = element_text(size = 16, margin = margin(r = 7)),
    axis.title.y.right = element_text(size = 16, margin = margin(l = 7)),
    axis.line = element_line(size = size),
    axis.ticks = element_line(size = size, color = 'black'),
    panel.border = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.background = element_rect(color = "white", fill = "white", size = 1.5, linetype = "solid"), 
    plot.title = element_text(hjust = 0.5, size = 10)
  )

# Combine them
plot = hist_plot / z_plot + plot_layout(heights = c(2, 1))

plot + theme_bw() + theme(
  axis.text.x = element_text(size = 14, color = 'black'),
  axis.text.y = element_text(size = 14, color = 'black'),
  axis.title.x = element_text(size = 16, margin = margin(t = 7)),
  axis.title.y.left = element_text(size = 16, margin = margin(r = 7)),
  axis.title.y.right = element_text(size = 16, margin = margin(l = 7)),
  axis.line = element_line(size = size),
  axis.ticks = element_line(size = size, color = 'black'),
  panel.border = element_blank(),
  legend.position = "none",
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text.x = element_text(size = 14, face = "bold"),
  strip.background = element_rect(color = "white", fill = "white", size = 1.5, linetype = "solid"), 
  plot.title = element_text(hjust = 0.5, size = 10)
)

# What position in a transcript is hit? (LTRs_)

input = dplyr::filter(unique_data, LTR == T, alignment_score >= 80) %>%
  group_by(location) %>%
  summarize(n = n())

# What position in a transcript is hit? (SCAN genes)

input = dplyr::filter(unique_data, alignment_score >= 80) %>%
  group_by(location) %>%
  summarize(n = n())
    
size = 0.3527778

plot = ggplot(input, aes(x = location, y= n)) +
  geom_bar(stat='identity', fill = 'grey', width = 0.75) +
  xlab('Target site location') +
  ylab('# of target sites') +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  geom_text(aes(label =n), vjust = -0.5, size = 2.75)

plot + theme_bw() + theme(
  axis.text.x = element_text(size = 8, color = 'black',  angle = 45, vjust = 0.5),
  axis.text.y = element_text(size = 8, color = 'black'),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 10, margin = margin(r = 5)),
  axis.line = element_line(size = size),
  axis.ticks = element_line(size = size, color = 'black'),
  panel.border = element_blank(),
  legend.position = "none",
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text.x = element_text(size = 14, face = "bold"),
  strip.background = element_rect(color = "white", fill = "white", size = 1.5, linetype = "solid"), 
  plot.title = element_text(hjust = 0.5, size = 10)
)

# Heatmap 

input = dplyr::filter(unique_data, LTR == T, LTR_family != 'LTR', alignment_score >= 75) %>%
  group_by(tRNA_anticodon, LTR_family) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = 'tRNA_anticodon', values_from = 'n') %>%
  as.data.frame()

my_palette = colorRampPalette(c("#f0f0f0", "#b30000"))(100)

input[is.na(input)] = 0

row.names(input) = input$LTR_family

input = dplyr::select(input, -c('LTR_family'))

pheatmap(input, cluster_col = T, cluster_row = T, scale = 'row', color = my_palette)

pheatmap(input, cluster_col = T, cluster_row = T)


