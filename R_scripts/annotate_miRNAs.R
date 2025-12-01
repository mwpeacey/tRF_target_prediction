library(tidyverse)
library(glue)
library(GenomicFeatures)
library(pheatmap)
library(AnnotationHub)

################################################################################
## Import data
################################################################################

# Load data

canonical_chromosomes = paste0("chr", c(1:19, "X", "Y"))

data = read.csv(file = '~/tRF_target_prediction/import/miranda/test/miranda_output_150.csv', header = TRUE) %>%
  dplyr::rename(start = genomic_start, end = genomic_end) %>%
  dplyr::filter(alignment_score >= 150) %>%
  dplyr::filter(seqnames %in% canonical_chromosomes)

################################################################################
## Annotate hit overlap with GENCODE transcripts
################################################################################

# Build GRanges of hits
gr_hits = GenomicRanges::makeGRangesFromDataFrame(
  data,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

# Load annotation 
ah = AnnotationHub::AnnotationHub()
txdb = ah[["AH75036"]]

# Extract feature GRanges
gr_tx    = GenomicFeatures::transcripts(txdb)
gr_exons = GenomicFeatures::exons(txdb)
gr_cds   = GenomicFeatures::cdsBy(txdb, by = "tx") %>% unlist()
gr_utr5  = GenomicFeatures::fiveUTRsByTranscript(txdb) %>% unlist()
gr_utr3  = GenomicFeatures::threeUTRsByTranscript(txdb) %>% unlist()

# Harmonize seqlevel styles
GenomeInfoDb::seqlevelsStyle(gr_hits)   = "UCSC"
GenomeInfoDb::seqlevelsStyle(gr_tx)     = "UCSC"
GenomeInfoDb::seqlevelsStyle(gr_exons)  = "UCSC"
GenomeInfoDb::seqlevelsStyle(gr_cds)    = "UCSC"
GenomeInfoDb::seqlevelsStyle(gr_utr5)   = "UCSC"
GenomeInfoDb::seqlevelsStyle(gr_utr3)   = "UCSC"

gr_exonsByTx = GenomicFeatures::exonsBy(txdb, by = "tx")
GenomeInfoDb::seqlevelsStyle(gr_exonsByTx) = "UCSC"

gr_exon_tx = unlist(gr_exonsByTx)
S4Vectors::mcols(gr_exon_tx)$tx_id = rep(
  names(gr_exonsByTx),
  S4Vectors::elementNROWS(gr_exonsByTx)
)

# Map each hit to overlapping exon → transcript
S4Vectors::mcols(gr_hits)$hit_idx = seq_along(gr_hits)
ex_ol = GenomicRanges::findOverlaps(gr_hits, gr_exon_tx, type = "within")

tx_map = dplyr::tibble(
  hit_idx = S4Vectors::mcols(gr_hits)$hit_idx[S4Vectors::queryHits(ex_ol)],
  transcript = S4Vectors::mcols(gr_exon_tx)$tx_id[S4Vectors::subjectHits(ex_ol)]
) %>%
  dplyr::group_by(hit_idx) %>%
  dplyr::summarize(transcript = paste(unique(transcript), collapse = ";"), .groups = "drop")

# Region type annotations
annot = dplyr::tibble(
  hit_idx        = S4Vectors::mcols(gr_hits)$hit_idx,
  overlaps_5utr  = IRanges::overlapsAny(gr_hits, gr_utr5),
  overlaps_3utr  = IRanges::overlapsAny(gr_hits, gr_utr3),
  overlaps_cds   = IRanges::overlapsAny(gr_hits, gr_cds),
  overlaps_exon  = IRanges::overlapsAny(gr_hits, gr_exons)
)

# Join annotations
annot = dplyr::left_join(annot, tx_map, by = "hit_idx")

data = data %>%
  dplyr::bind_cols(annot %>% dplyr::select(-hit_idx)) %>%
  dplyr::mutate(
    gencode_location = dplyr::case_when(
      overlaps_5utr ~ "Exon - 5' UTR",
      overlaps_3utr ~ "Exon - 3' UTR",
      overlaps_cds  ~ "Exon - CDS",
      overlaps_exon ~ "Exon - non-coding",
      TRUE          ~ "Intergenic"
    )
  ) %>%
  dplyr::rename(gencode_transcripts = transcript) %>%
  dplyr::select(-c('overlaps_5utr', 'overlaps_3utr', 'overlaps_cds',
                   'overlaps_exon'))

################################################################################
## Collapse overlapping sites
################################################################################

# Convert to GRanges
gr = GenomicRanges::GRanges(
  seqnames = data$seqnames,
  ranges = IRanges::IRanges(start = data$start, end = data$end),
  strand = data$strand
)

# Reduce overlapping regions
reduced_gr = GenomicRanges::reduce(gr)

# Find overlaps (mapping each reduced range to original hits)
hits = GenomicRanges::findOverlaps(reduced_gr, gr)
hits_df = tibble::as_tibble(hits)

# Add original alignment_score column for ranking
hits_df = hits_df %>%
  dplyr::mutate(
    alignment_score = data$alignment_score[subjectHits]
  )

# Identify best-scoring original hit for each reduced region
best_hits = hits_df %>%
  dplyr::group_by(queryHits) %>%
  dplyr::slice_max(order_by = alignment_score, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# Count total hits per reduced region
tRF_counts = hits_df %>%
  dplyr::count(queryHits, name = "collapsed_tRF_count")

# Merge best hit info with tRF counts
best_hits = dplyr::left_join(best_hits, tRF_counts, by = "queryHits")

# Extract rows from data
unique_data = data[best_hits$subjectHits, ]
unique_data$collapsed_tRF_count = best_hits$collapsed_tRF_count

################################################################################
## Export
################################################################################

write_csv(data, file = '~/tRF_target_prediction/import/miranda/test/miranda_miRNA_output_annotated.csv')
write_csv(unique_data, file = '~/tRF_target_prediction/import/miranda/test/miranda_output_miRNA_unique_annotated.csv')
