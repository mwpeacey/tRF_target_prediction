################################################################################
# annotate_targets.R
#
# Annotates genomic features of hit sites, including hit transcript biotype,
# the position of the hit relative to the CDS (for protein coding genes), and
# overlap with long terminal repeats.
################################################################################

library(tidyverse)
library(glue)
library(GenomicFeatures)
library(pheatmap)
library(AnnotationHub)
library(GenomicScores)

#data = read_csv('import/miranda/miranda_output_annotated.csv')
#unique_data = read_csv('import/miranda/miranda_output_unique_annotated.csv')

################################################################################
## Import data
################################################################################

# Load data

canonical_chromosomes = paste0("chr", c(1:19, "X", "Y"))

data = read.csv(file = 'import/miranda/miranda_output_70.csv', header = TRUE) %>%
  dplyr::rename(start = genomic_start, end = genomic_end) %>%
  dplyr::filter(alignment_score >= 70) %>%
  dplyr::filter(seqnames %in% canonical_chromosomes)

# Add tRF information

tRF_infomation = read.csv('import/mm10_tRF3b.csv') %>%
  dplyr::rename(tRF = 'unique_name')

data$tRF = stringr::str_remove(data$tRF, '>')

data = merge(data, tRF_infomation, by = 'tRF') %>%
  separate(col = source_tRNAs, sep = '-', into = c('A', 'B', 'C', 'D'), remove = F) %>%
  unite(tRNA_anticodon, c('B', 'C'), remove = F, sep = '-') %>%
  dplyr::select(-c('A', 'B', 'C', 'D'))

################################################################################
## Annotate hit overlap with GENCODE transcripts
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
ex_ol = findOverlaps(gr_hits, gr_exon_tx, type="within")

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
    gencode_location = case_when(
      overlaps_5utr ~ "Exon - 5' UTR",
      overlaps_3utr ~ "Exon - 3' UTR",
      overlaps_cds  ~ "Exon - CDS",
      overlaps_exon ~ "Exon - non-coding",
      TRUE          ~ "Intergenic"
    ) 
  ) %>%
  dplyr::rename(gencode_transcripts = transcript) %>%
  dplyr::select(
    tRF, seqnames, alignment_score, energy, tRNA_anticodon,
    start, end, alignment_length,
    strand, gencode_location, gencode_transcripts
  )

################################################################################
## Annotate overlap with StringTie transcripts
################################################################################

# Import stringtie info

txdb_stringtie = makeTxDbFromGFF('import/transcriptome_assembly/stringtie_merged_filtered.gtf')
transcripts_gr = transcripts(txdb_stringtie)

transcripts_df = data.frame(
  tx_id = transcripts_gr$tx_id,
  tx_name = transcripts_gr$tx_name,
  tx_chrom = as.character(seqnames(transcripts_gr)),
  tx_strand = as.character(strand(transcripts_gr)),
  tx_start = start(transcripts_gr),
  tx_end = end(transcripts_gr)
)

# Build "splicings" dataframe

exons_by_tx = exonsBy(txdb_stringtie, by = "tx", use.names = TRUE)

# Import open reading frame information 

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

# Map ORFs to genomic coordinates

orfs_tx_ranges = GRanges(
  seqnames = open_reading_frames_annotated_filtered$target_id,
  ranges = IRanges(
    start = open_reading_frames_annotated_filtered$ORF_start,
    end = open_reading_frames_annotated_filtered$ORF_end
  ),
  strand = open_reading_frames_annotated_filtered$V6,
  ORF_ID = open_reading_frames_annotated_filtered$ORF_ID
)

genomic_orfs = mapFromTranscripts(orfs_tx_ranges, exons_by_tx)

# Build splicings dataframe

exons_flat <- unlist(exons_by_tx, use.names = FALSE)
exons_df <- as.data.frame(exons_flat)
exons_df$tx_name <- rep(names(exons_by_tx), lengths(exons_by_tx))
exons_df <- left_join(exons_df, transcripts_df, by = "tx_name")

splicings_df <- exons_df %>%
  transmute(
    tx_id = as.integer(tx_id),
    exon_rank = as.integer(exon_rank),
    exon_start = start,
    exon_end = end,
    exon_chrom = tx_chrom,
    exon_strand = tx_strand
  )

# Genomic ORFs already in GRanges (genomic_orfs)
orf_df <- as.data.frame(genomic_orfs)
orf_df$ORF_ID <- mcols(orfs_tx_ranges)$ORF_ID[orf_df$transcriptsHits]  # backtrack

# Create GRanges from exon coordinates
exon_gr <- GRanges(
  seqnames = splicings_df$exon_chrom,
  ranges = IRanges(splicings_df$exon_start, splicings_df$exon_end),
  strand = splicings_df$exon_strand
)

# Find overlaps between exons and CDS
hits <- findOverlaps(exon_gr, genomic_orfs)

# For each overlap, get intersection region
cds_ranges <- pintersect(exon_gr[queryHits(hits)], genomic_orfs[subjectHits(hits)])

# Extract relevant data
cds_df <- as.data.frame(cds_ranges)
cds_df$tx_id <- splicings_df$tx_id[queryHits(hits)]
cds_df$exon_rank <- splicings_df$exon_rank[queryHits(hits)]

cds_df <- cds_df %>%
  dplyr::rename(cds_start = start, cds_end = end)

# Collapse overlapping CDS fragments per exon
cds_df_collapsed <- cds_df %>%
  group_by(tx_id, exon_rank) %>%
  summarise(
    cds_start = min(cds_start, na.rm = TRUE),
    cds_end = max(cds_end, na.rm = TRUE),
    .groups = "drop"
  )

# Join collapsed CDS back to splicings
splicings_with_cds <- left_join(
  splicings_df,
  cds_df_collapsed,
  by = c("tx_id", "exon_rank")
)

# Build TxDb

txdb_custom = txdb_custom <- makeTxDb(
  transcripts = transcripts_df,
  splicings = splicings_with_cds
)

#Build GRanges of hits

gr_hits = makeGRangesFromDataFrame(
  data,
  seqnames.field   = "seqnames",
  start.field      = "start",
  end.field        = "end",
  strand.field     = "strand",
  keep.extra.columns = TRUE
)

# Extract feature GRanges
gr_tx     = transcripts(txdb_custom)
gr_exons  = exons(txdb_custom)
gr_cds    = cdsBy(txdb_custom, by="tx") %>% unlist()
gr_utr5   = fiveUTRsByTranscript(txdb_custom) %>% unlist()
gr_utr3   = threeUTRsByTranscript(txdb_custom) %>% unlist()

gr_exonsByTx = exonsBy(txdb_custom, by="tx")
gr_exon_tx <- unlist(gr_exonsByTx)

# capture the transcript ID in a metadata column
mcols(gr_exon_tx)$tx_id = rep(names(gr_exonsByTx),
                              elementNROWS(gr_exonsByTx))

# 2) map each hit to any overlapping exon → transcript
mcols(gr_hits)$hit_idx = seq_along(gr_hits)
ex_ol = findOverlaps(gr_hits, gr_exon_tx, type="within")

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
    stringtie_location = case_when(
      overlaps_5utr ~ "Exon - 5' UTR",
      overlaps_3utr ~ "Exon - 3' UTR",
      overlaps_cds  ~ "Exon - CDS",
      overlaps_exon ~ "Exon - non-coding",
      TRUE          ~ "Intergenic"
    )
  ) %>%
  dplyr::rename(stringtie_transcripts = transcript) %>%
  dplyr::select(
    tRF, seqnames, alignment_score, energy, tRNA_anticodon,
    start, end, alignment_length,
    strand, gencode_location, gencode_transcripts,
    stringtie_location, stringtie_transcripts
  )

################################################################################
## Annotate hit overlap with repeats
################################################################################

# Load and filter LTR annotation
LTR = rtracklayer::readGFF('import/annotation_tables/mm10_rmsk_TE.gtf') %>%
  dplyr::filter(class_id == 'LTR')

LTR = LTR[!grepl('int', LTR$gene_id),]

# Extend LTR elements by 200bp downstream
subject = makeGRangesFromDataFrame(LTR, keep.extra.columns = TRUE)
end(subject[strand(subject) == '+']) = end(subject[strand(subject) == '+']) + 200
start(subject[strand(subject) == '-']) = start(subject[strand(subject) == '-']) - 200

# Query object from your data
query = GRanges(data)

# Find overlaps
hits = findOverlaps(query, subject, type = 'any', ignore.strand = TRUE)
hits_df = as.data.frame(hits)

# Extract the overlapping metadata
matched_LTR = subject[hits_df$subjectHits]
overlap_info = data.frame(
  queryHits = hits_df$queryHits,
  LTR_family = mcols(matched_LTR)$family_id,
  LTR_gene_id = mcols(matched_LTR)$gene_id
)

# Keep only the first match per query (or you can customize to collapse all hits)
overlap_info = overlap_info[!duplicated(overlap_info$queryHits), ]

# Add LTR columns to main data
data = data %>%
  dplyr::mutate(row_id = row_number()) %>%
  left_join(overlap_info, by = c("row_id" = "queryHits")) %>%
  dplyr::mutate(LTR = !is.na(LTR_family)) %>%
  dplyr::select(-row_id)  # clean up temporary column


################################################################################
## Annotate target site hits with phast con scores
################################################################################

phastCons = getGScores("phastCons60way.UCSC.mm10")

gr_hits = makeGRangesFromDataFrame(
  data,
  seqnames.field   = "seqnames",
  start.field      = "start",
  end.field        = "end",
  strand.field     = "strand",
  keep.extra.columns = TRUE
)

scores = gscores(phastCons, gr_hits)

mcols(scores)$phastCons_mean = mcols(scores)$default
mcols(scores)$default = NULL

data = as.data.frame(scores)

################################################################################
## Collapse overlapping sites
################################################################################

# Convert to GRanges
gr = GRanges(
  seqnames = data$seqnames,
  ranges = IRanges(start = data$start, end = data$end),
  strand = data$strand
)

# Reduce overlapping ranges
reduced_gr = reduce(gr)

# Find overlaps
hits = findOverlaps(reduced_gr, gr)
hits_df = as.data.frame(hits)

# Count how many original ranges map to each reduced range
tRF_counts = hits_df %>%
  group_by(queryHits) %>%
  summarise(collapsed_tRF_count = n(), .groups = "drop")

# Pick the best hit for each reduced range
best_hits = hits_df %>%
  group_by(queryHits) %>%
  slice_max(order_by = data$alignment_score[subjectHits], n = 1, with_ties = FALSE)

# Extract best original rows
unique_data = data[best_hits$subjectHits, ]

# Add the tRF count to unique_data
unique_data$collapsed_tRF_count = tRF_counts$collapsed_tRF_count[match(best_hits$queryHits, tRF_counts$queryHits)]

################################################################################
## Export
################################################################################

write_csv(data, file = 'import/miranda/miranda_output_annotated.csv')
write_csv(unique_data, file = 'import/miranda/miranda_output_unique_annotated.csv')

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