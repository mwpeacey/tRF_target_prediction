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
library(AnnotationHub)

#data = read_csv('import/miranda/miranda_output_annotated.csv')
#unique_data = read_csv('import/miranda/miranda_output_unique_annotated.csv')

################################################################################
## Import data
################################################################################

# Load data

canonical_chromosomes = paste0("chr", c(1:22, "X", "Y"))

data = read.csv(file = 'import/miranda/hg38/miranda_output_70.csv', header = TRUE) %>%
  dplyr::rename(start = genomic_start, end = genomic_end) %>%
  dplyr::filter(alignment_score >= 80) %>%
  dplyr::filter(seqnames %in% canonical_chromosomes)

# Add tRF information

tRF_infomation = read.csv('import/hg38_tRF3b.csv') %>%
  dplyr::rename(tRF = 'unique_name')

tRF_infomation = tRF_infomation %>%
  dplyr::mutate(
    source_tRNAs = source_tRNAs %>%
      str_split("_") %>% 
      map(~ .x %>%
            str_replace("-\\d+$", "") %>%  
            unique()                       
      ) %>%
      map_chr(~ paste(.x, collapse = "_"))
  )

data$tRF = stringr::str_remove(data$tRF, '>')

data = merge(data, tRF_infomation, by = 'tRF') %>%
  separate(col = source_tRNAs, sep = '-', into = c('A', 'B', 'C', 'D'), remove = F) %>%
  unite(tRNA_anticodon, c('B', 'C'), remove = F, sep = '-') %>%
  dplyr::select(-c('A', 'B', 'C', 'D'))

################################################################################
## Annotate hit overlap with repeats
################################################################################

# Load and filter LTR annotation
LTR = read.csv('import/annotation_tables/UCSC_hg38_LTR.csv', comment.char = '#', header = FALSE)
colnames(LTR) = c('bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns', 'genoName', 'genoStart', 'genoEnd', 'genoLeft', 'strand', 'repName', 'repClass', 'repFamily', 'repStart', 'repEnd', 'repLeft', 'id')
LTR = LTR[!grepl('int', LTR$repName),]

# Extend LTR elements by 200bp downstream
subject = makeGRangesFromDataFrame(LTR, 
                                   keep.extra.columns = TRUE,
                                   start.field = 'genoStart',
                                   end.field = 'genoEnd',
                                   seqnames.field = 'genoName')
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
  LTR_family = mcols(matched_LTR)$repFamily,
  LTR_gene_id = mcols(matched_LTR)$repName
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

query(ah, c("Homo sapiens", "EnsDb", "104"))
txdb = ah[["AH95744"]]

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

annot = dplyr::tibble(
  hit_idx        = S4Vectors::mcols(gr_hits)$hit_idx,
  overlaps_5utr  = IRanges::overlapsAny(gr_hits, gr_utr5, type = "within"),
  overlaps_3utr  = IRanges::overlapsAny(gr_hits, gr_utr3, type = "within"),
  overlaps_cds   = IRanges::overlapsAny(gr_hits, gr_cds,  type = "within"),
  overlaps_exon  = IRanges::overlapsAny(gr_hits, gr_exons,type = "within")
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
## Add gene names and gene IDs directly from TxDb
################################################################################

# Build mapping from transcript → gene_id
tx_to_gene = dplyr::tibble(
  transcript_id = gr_tx$tx_name,
  gene_id = gr_tx$gene_id
)

# Extract gene name from genes() metadata
genes_gr = GenomicFeatures::genes(txdb)

gene_metadata = S4Vectors::mcols(genes_gr) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name)

# Join gene_name to transcript info
tx_to_gene = dplyr::left_join(tx_to_gene, gene_metadata, by = "gene_id")

# Expand transcripts per row
data = data %>% dplyr::mutate(row_id = dplyr::row_number())

tx_long = data %>%
  dplyr::select(row_id, gencode_transcripts) %>%
  tidyr::separate_rows(gencode_transcripts, sep = ";")

# Join gene_id + gene_name
tx_long_annot = dplyr::left_join(
  tx_long,
  tx_to_gene,
  by = c("gencode_transcripts" = "transcript_id")
) %>%
  dplyr::filter(!is.na(gene_name))

# Collapse back to per-row
gene_info_per_row = tx_long_annot %>%
  dplyr::group_by(row_id) %>%
  dplyr::summarize(
    gencode_gene_id   = paste(unique(gene_id), collapse = ";"),
    gencode_gene_name = paste(unique(gene_name), collapse = ";"),
    .groups = "drop"
  )

# Join into data
data = data %>%
  dplyr::left_join(gene_info_per_row, by = "row_id") %>%
  dplyr::select(-row_id)

################################################################################
## Annotate Gag-like genes
################################################################################

#gag_genes = read_csv('import/annotation_tables/mouse_gag_genes_new.csv')

#data = dplyr::mutate(data, gag_gene = case_when(gencode_gene_name %in% gag_genes$gene_name ~ T,
#                                                T ~ F))

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
## Exon distance analysis
################################################################################

#Build GRanges of hits

gr_hits = makeGRangesFromDataFrame(
  unique_data,
  seqnames.field   = "seqnames",
  start.field      = "start",
  end.field        = "end",
  strand.field     = "strand",
  keep.extra.columns = TRUE
)

# Load annotation 

ah  = AnnotationHub()
txdb = ah[["AH75036"]]
gr_tx = transcripts(txdb)
gr_exons  = exonsBy(txdb, by = "tx")

# Harmonize seqlevel styles
seqlevelsStyle(gr_hits)  = "UCSC"
seqlevelsStyle(gr_tx)    = "UCSC"
seqlevelsStyle(gr_exons) = "UCSC"

# Filter to protein-coding transcripts

protein_coding_tx = gr_tx[gr_tx$tx_biotype == "protein_coding"]
protein_coding_ids = protein_coding_tx$tx_name

exons_by_tx_coding = gr_exons[names(gr_exons) %in% protein_coding_ids]

# Get TSS positions of protein-coding transcripts

tss_pc = promoters(protein_coding_tx, upstream = 1, downstream = 1)

# Find overlaps with LTRs

LTR = rtracklayer::readGFF('import/annotation_tables/mm10_rmsk_TE.gtf') %>%
  dplyr::filter(class_id == 'LTR')

ltr_gr = makeGRangesFromDataFrame(LTR, keep.extra.columns = TRUE)

hits = findOverlaps(tss_pc, ltr_gr, ignore.strand = FALSE)

ltr_initiated_tx = protein_coding_tx[queryHits(hits)]$tx_name

exons_by_tx_ltr_coding = gr_exons[names(gr_exons) %in% ltr_initiated_tx]

# Sort exons by their 5' end depending on strand

first_splice_sites_list = lapply(exons_by_tx_ltr_coding, function(exs) {
  strand_sign = as.character(strand(exs)[1])
  exs_sorted = if (strand_sign == "+") exs[order(start(exs))] else exs[order(end(exs), decreasing = TRUE)]
  
  first_exon = exs_sorted[1]
  pos = if (strand_sign == "+") end(first_exon) else start(first_exon)
  
  GRanges(
    seqnames = seqnames(first_exon),
    ranges = IRanges(start = pos, end = pos),
    strand = strand(first_exon),
    tx_name = names(exs)[1]
  )
})

first_splice_sites = Reduce(c, first_splice_sites_list)

# Compute nearest distance

hits = distanceToNearest(gr_hits, first_splice_sites)

# Attach distances back to your data

gr_hits$distance_to_5p_splice = NA_integer_
gr_hits$nearest_tx = NA_character_
gr_hits$distance_to_5p_splice[queryHits(hits)] = mcols(hits)$distance

unique_data = as.data.frame(gr_hits)

################################################################################
## Export
################################################################################

write_csv(data, file = 'import/miranda/hg38/hg38_80_output_annotated.csv')
write_csv(unique_data, file = 'import/miranda/hg38/hg38_80_output_unique_annotated.csv')

