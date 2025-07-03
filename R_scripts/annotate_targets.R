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
  dplyr::rename(start = genomic_start, end = genomic_end) %>%
  dplyr::filter(alignment_score >= 70)

# Add tRF information

tRF_infomation = read.csv('import/mm10_tRF3b.csv') %>%
  dplyr::rename(tRF = 'unique_name')

data$tRF = stringr::str_remove(data$tRF, '>')

data = merge(data, tRF_infomation, by = 'tRF') %>%
  separate(col = source_tRNAs, sep = '-', into = c('A', 'B', 'C', 'D'), remove = F) %>%
  unite(tRNA_anticodon, c('B', 'C'), remove = F, sep = '-') %>%
  dplyr::select(-c('A', 'B', 'C', 'D'))

################################################################################
## Annotate hit overlap with transcripts
################################################################################

#Build GRanges of hits
gr_hits = makeGRangesFromDataFrame(
  data,
  seqnames.field   = "seqnames",
  start.field      = "start",
  end.field        = "end",
  strand.field     = "strand",
  keep.extra.columns = TRUE
)

# Load annotation 
ah  = AnnotationHub()
txdb = ah[["AH75036"]]

# Extract feature GRanges
gr_tx     = transcripts(txdb)
gr_exons  = exons(txdb)
gr_cds    = cdsBy(txdb, by="tx") %>% unlist()
gr_utr5   = fiveUTRsByTranscript(txdb) %>% unlist()
gr_utr3   = threeUTRsByTranscript(txdb) %>% unlist()

# Harmonize seqlevel styles
seqlevelsStyle(gr_hits)  = "UCSC"
seqlevelsStyle(gr_tx)    = "UCSC"
seqlevelsStyle(gr_exons) = "UCSC"
seqlevelsStyle(gr_cds)   = "UCSC"
seqlevelsStyle(gr_utr5)  = "UCSC"
seqlevelsStyle(gr_utr3)  = "UCSC"

gr_exonsByTx = exonsBy(txdb, by="tx")
seqlevelsStyle(gr_exonsByTx) = "UCSC"

gr_exon_tx <- unlist(gr_exonsByTx)
# capture the transcript ID in a metadata column
mcols(gr_exon_tx)$tx_id = rep(names(gr_exonsByTx),
                               elementNROWS(gr_exonsByTx))

# 2) map each hit to any overlapping exon → transcript
mcols(gr_hits)$hit_idx = seq_along(gr_hits)
ex_ol = findOverlaps(gr_hits, gr_exon_tx, type="any")

tx_map = tibble(
  hit_idx    = mcols(gr_hits)$hit_idx[queryHits(ex_ol)],
  transcript = mcols(gr_exon_tx)$tx_id[subjectHits(ex_ol)]
) %>%
  group_by(hit_idx) %>%
  summarize(transcript = paste(unique(transcript), collapse = ";"), .groups="drop")

# 3) Build your logical annotations as before
annot = tibble(
  hit_idx       = mcols(gr_hits)$hit_idx,
  overlaps_5utr  = overlapsAny(gr_hits, gr_utr5),
  overlaps_3utr  = overlapsAny(gr_hits, gr_utr3),
  overlaps_cds   = overlapsAny(gr_hits, gr_cds),
  overlaps_exon  = overlapsAny(gr_hits, gr_exons)  # unstranded exons
)

# 4) Stitch transcript names into annot, then join into your data
annot = annot %>%
  left_join(tx_map, by="hit_idx")

data <- data %>%
  bind_cols(annot %>% dplyr::select(-hit_idx)) %>%
  mutate(
    location = case_when(
      overlaps_5utr ~ "Exon - 5' UTR",
      overlaps_3utr ~ "Exon - 3' UTR",
      overlaps_cds  ~ "Exon - CDS",
      overlaps_exon ~ "Exon - non-coding",
      TRUE          ~ "Intergenic"
    )
  ) %>%
  dplyr::select(
    tRF, seqnames, alignment_score, energy,
    start, end, alignment_length,
    strand, location, transcript
  )

# Does it overlap a transcript with a SCAN domain?

mouse_scan_genes = read_csv('import/annotation_tables/mouse_scan_genes.csv')

subject = exonsBy(txdb, by = 'tx')

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
## Annotate overlap with chimeric transcripts 
################################################################################

# Read chimeras
chimeras = read_csv('~/Downloads/Oomen_chimeras.csv') %>%
  separate(TE_position, sep = ':', into = c('seqnames', 'position')) %>%
  separate(position, sep = '-', into = c('start', 'end')) %>%
  dplyr::select(c('TE_strand', 'seqnames', 'start', 'end')) %>%
  dplyr::rename(strand = 'TE_strand')

# Create GRanges object 
subject = makeGRangesFromDataFrame(chimeras, keep.extra.columns = TRUE)
seqlevelsStyle(subject) = 'NCBI'

# Create GRanges object for query data
query = GRanges(data)
seqlevelsStyle(query) = 'NCBI'

# Find overlaps between query and subject
chimera_overlap = findOverlaps(query = query, subject = subject, type = 'within', ignore.strand = T)
chimera_overlap = as.data.frame(chimera_overlap)

# Add chimera indicator to data
data = data %>%
  dplyr::mutate(chimera = seq_len(n()) %in% chimera_overlap$queryHits)

################################################################################
## Annotate Stringtie transcripts with ORFs
################################################################################

# Import open reading frames deleted in assembled transcripts

open_reading_frames = read.table('~/tRF_targets_new/mouse_embryo_tRF3b_v2/mouse_embryo_assembled_transcripts_ORFs.bed') %>%
  dplyr::rename(target_id = V1, ORF_start = V2, ORF_end = V3, ORF_info = V4) %>%
  tidyr::separate(ORF_info, sep = ';', into = c('ORF_ID', 'ORF_type', 'ORF_length', 'ORF_frame', 'other'))

open_reading_frames$ORF_ID = stringr::str_remove(open_reading_frames$ORF_ID, 'ID=')
open_reading_frames$ORF_type = stringr::str_remove(open_reading_frames$ORF_type, 'ORF_type=')
open_reading_frames$ORF_length = stringr::str_remove(open_reading_frames$ORF_length, 'ORF_len=')
open_reading_frames$ORF_frame = stringr::str_remove(open_reading_frames$ORF_frame, 'ORF_frame=')

# Add blast information if present

blast_output = read.csv(file = '~/tRF_targets_new/mouse_embryo_tRF3b_v2/mouse_embryo_assembled_transcripts_ORFs_blast.csv', header = F) %>%
  dplyr::rename(ORF_ID = V1, blast_match = V2, pident = V3, length = V4, mismatch = V5, gapopen = V6, qstart = V7, qend = V8,
                sstart = V9, send = V10, evalue = V11, bitscore = V12)

open_reading_frames_annotated = merge(open_reading_frames, blast_output, by = 'ORF_ID')

open_reading_frames_annotated_filtered = open_reading_frames_annotated %>%
  group_by(target_id) %>%
  slice_max(bitscore, with_ties = FALSE) %>%
  ungroup()

data = merge(data, open_reading_frames_annotated_filtered, by = 'target_id', all.x = T)

## Annotate tRF position

data = dplyr::mutate(data, position_relative_to_ORF = case_when(is.na(ORF_ID) ~ 'None',
                                                                transcript_end < ORF_start ~ '5_UTR',
                                                                transcript_start > ORF_end ~ '3_UTR',
                                                                T ~ 'CDS'))

################################################################################
## Collapse overlapping sites
################################################################################

unique_data = data %>%
  arrange(desc(alignment_score), energy) %>%
  group_by(seqnames, start, end, strand) %>%
  dplyr::slice(1) %>%
  ungroup()

################################################################################
## Export
################################################################################

write_csv(data, file = 'import/miranda_output_annotated.csv')
write_csv(unique_data, file = 'import/miranda_output_unique_annotated.csv')

data = read_csv('import/miranda_output_annotated.csv')
unique_data = read_csv('import/miranda_output_unique_annotated.csv')


## Overlap 

query = GRanges(unique_data)

CLIP = read.table('~/Downloads/GSE140838_CLIP_tagIP_vs_untagIP_peaksFC.bed') %>%
  dplyr::rename(seqnames = V1, start = V2, end = V3, p = V4, strand = V6) %>%
  dplyr::filter(p >= 1)

subject = GRanges(CLIP)

overlap = findOverlaps(query = query, subject = subject, type = 'any', ignore.strand = F)
overlap = as.data.frame(overlap)

unique_data  = unique_data %>%
  dplyr::mutate(AGO_peak = dplyr::row_number() %in% overlap$queryHits)

df = dplyr::filter(unique_data, AGO_peak == T, LTR == T)

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

# What position in a transcript is hit? (LTRs)

Leu_tRFs = c('tRF_76', 'tRF_77', 'tRF_78', 'tRF_79')

input = dplyr::filter(unique_data, LTR == T, alignment_score >= 75, tRF %in% Leu_tRFs) %>%
  group_by(location, LTR_family) %>%
  summarize(n = n())

# What position in a transcript is hit? (SCAN genes)

input = dplyr::filter(unique_data, SCAN == T, alignment_score >= 75) %>%
  group_by(location) %>%
  summarize(n = n())
    
size = 0.3527778

plot = ggplot(input, aes(x = location, y = n, fill = LTR_family)) +
  geom_bar(stat='identity', width = 0.75) +
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


