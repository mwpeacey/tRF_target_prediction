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

canonical_chromosomes = paste0("chr", c(1:19, "X", "Y"))

data = read.csv(file = 'import/miranda/new/miranda_output_70.csv', header = TRUE) %>%
  dplyr::rename(start = genomic_start, end = genomic_end) %>%
  dplyr::filter(alignment_score >= 70) %>%
  dplyr::filter(seqnames %in% canonical_chromosomes)

# Add tRF information

tRF_infomation = read.csv('import/mm10_tRF3b.csv') %>%
  dplyr::rename(tRF = 'unique_name') %>%
  # tDRnamer name is taken from the miRanda query id below, so drop any
  # tDRnamer column here to avoid a duplicate/clashing column at the merge.
  dplyr::select(-dplyr::any_of('tDRnamer'))

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

# miRanda query names now arrive as ">tDR-55:76-Ala-AGC-1|tRF_1", carrying both
# the tDRnamer name and the original unique id joined by "|". Strip the ">",
# split on "|" into a tDRnamer column and the unique tRF id used as the key.
data = data %>%
  dplyr::mutate(
    tRF      = stringr::str_remove(tRF, '^>'),
    tDRnamer = stringr::str_remove(tRF, '\\|.*$'),   # e.g. tDR-55:76-Ala-AGC-1
    tRF      = stringr::str_remove(tRF, '^.*\\|')     # e.g. tRF_1
  )

data = merge(data, tRF_infomation, by = 'tRF') %>%
  separate(col = source_tRNAs, sep = '-', into = c('A', 'B', 'C', 'D'), remove = F) %>%
  unite(tRNA_anticodon, c('B', 'C'), remove = F, sep = '-') %>%
  dplyr::select(-c('A', 'B', 'C', 'D'))

################################################################################
## Annotate hit overlap with repeats
################################################################################

# Load and filter LTR annotation
LTR = read.csv('import/annotation_tables/UCSC_mm10_LTR.csv', header = T)
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

gag_genes = read_csv('import/annotation_tables/mouse_gag_genes_new.csv')

data = dplyr::mutate(data, gag_gene = case_when(gencode_gene_name %in% gag_genes$gene_name ~ T,
                                                T ~ F))

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

mcols(gr_exon_tx)$tx_id = rep(names(gr_exonsByTx),
                              elementNROWS(gr_exonsByTx))

mcols(gr_hits)$hit_idx = seq_along(gr_hits)
ex_ol = findOverlaps(gr_hits, gr_exon_tx, type="within")

tx_map = tibble(
  hit_idx = mcols(gr_hits)$hit_idx[queryHits(ex_ol)],
  tx_id   = as.integer(mcols(gr_exon_tx)$tx_id[subjectHits(ex_ol)])  
) %>%
  left_join(transcripts_df %>% dplyr::select(tx_id, tx_name), by = "tx_id") %>%
  group_by(hit_idx) %>%
  summarize(transcript = paste(unique(tx_name), collapse = ";"), .groups = "drop")

annot = tibble(
  hit_idx       = mcols(gr_hits)$hit_idx,
  overlaps_5utr  = overlapsAny(gr_hits, gr_utr5, type = "within"),
  overlaps_3utr  = overlapsAny(gr_hits, gr_utr3, type = "within"),
  overlaps_cds   = overlapsAny(gr_hits, gr_cds,  type = "within"),
  overlaps_exon  = overlapsAny(gr_hits, gr_exons,type = "within")
)

annot = annot %>%
  left_join(tx_map, by="hit_idx")

orf_match_table = open_reading_frames_annotated_filtered %>%
  dplyr::select(target_id, blast_match)

tx_to_protein = tx_map %>%
  separate_rows(transcript, sep = ";") %>%
  dplyr::mutate(tx_id = as.integer(transcript)) %>%
  left_join(transcripts_df %>% dplyr::select(tx_id, tx_name), by = "tx_id") %>%
  left_join(open_reading_frames_annotated_filtered %>% dplyr::select(target_id, blast_match),
            by = c("tx_name" = "target_id")) %>%
  dplyr::filter(!is.na(blast_match)) %>%
  group_by(hit_idx) %>%
  summarize(blast_match = paste(unique(blast_match), collapse = ";"), .groups = "drop")

data$hit_idx = seq_len(nrow(data))

data = data %>%
  bind_cols(annot %>% dplyr::select(-hit_idx)) %>%
  left_join(tx_to_protein, by = "hit_idx") %>%
  dplyr::mutate(
    stringtie_location = case_when(
      overlaps_5utr ~ "Exon - 5' UTR",
      overlaps_3utr ~ "Exon - 3' UTR",
      overlaps_cds  ~ "Exon - CDS",
      overlaps_exon ~ "Exon - non-coding",
      TRUE          ~ "Intergenic"
    )
  ) %>%
  dplyr::rename(stringtie_transcripts = transcript) %>%
  dplyr::select(-c('overlaps_5utr', 'overlaps_3utr', 'overlaps_cds',
                   'overlaps_exon', 'hit_idx'))

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

#imprinted_genes = c('Rtl1', 'Rasgrf1', 'Impact', 'Slc38a4', 'Kcnq1ot1', 
#                    'Mest', 'Snrpn', 'Cdh15', 'Peg3', 'Peg10')

imprinted_genes = read_csv('import/annotation_tables/geneimprint_mouse_imprinted_genes.csv')

imprinted_gene_names = dplyr::filter(imprinted_genes, Status == "Imprinted")$Gene

annotate_imprint_status <- function(gene_names, imprinted_gene_names) {
  purrr::map_chr(gene_names, function(x) {
    if (is.na(x) || !nzchar(x)) {
      return(NA_character_)
    }

    genes <- stringr::str_trim(stringr::str_split(x, ";")[[1]])
    genes <- genes[nzchar(genes)]

    if (length(genes) == 0L) {
      return(NA_character_)
    }

    if (any(genes %in% imprinted_gene_names)) {
      "Imprinted"
    } else {
      "Not imprinted"
    }
  })
}

data$imprint_status = annotate_imprint_status(data$gencode_gene_name, imprinted_gene_names)
unique_data$imprint_status = annotate_imprint_status(unique_data$gencode_gene_name, imprinted_gene_names)

data$is_imprinted = !is.na(data$imprint_status) & data$imprint_status == "Imprinted"
unique_data$is_imprinted = !is.na(unique_data$imprint_status) &
  unique_data$imprint_status == "Imprinted"

write_csv(data, file = 'import/miranda/miranda_output_annotated.csv')
write_csv(unique_data, file = 'import/miranda/miranda_output_unique_annotated.csv')

################################################################################
## Alignment visualization (HTML report)
################################################################################

#' Generate a browsable HTML report of miRanda alignments with colour-coded
#' match characters. Each target site is rendered as a card showing annotation
#' metadata and the reconstructed alignment block. A search bar allows filtering
#' by gene name, tRF, coordinates, or LTR family.
#'
#' @param df        A data frame containing at minimum: tRF, seqnames, start,
#'                  end, strand, alignment_score, energy, query_alignment,
#'                  match_string, ref_alignment. Optional columns (used when
#'                  present): gencode_gene_name, gencode_location, LTR,
#'                  LTR_family, LTR_gene_id, stringtie_location,
#'                  tRNA_anticodon, imprint_status, is_imprinted.
#' @param output_file  Path to write the HTML report.
#' @param title     Title displayed at the top of the report.
#' @param max_rows  Maximum number of rows to render. Set to Inf to include all.

generate_alignment_report <- function(df, output_file,
                                      title = "tRF Target Site Alignment Browser",
                                      max_rows = Inf) {

  # Check that alignment columns are present
  required_cols <- c("query_alignment", "match_string", "ref_alignment")
  if (!all(required_cols %in% names(df))) {
    warning("Alignment columns not found in data frame. Skipping report generation.")
    return(invisible(NULL))
  }

  # Truncate if necessary
  if (nrow(df) > max_rows) {
    message(sprintf("Truncating report to %d rows (of %d total). Increase max_rows to include more.",
                    max_rows, nrow(df)))
    df <- df[seq_len(max_rows), ]
  }

  # Helpers
  escape_html <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;",  x, fixed = TRUE)
    x <- gsub(">", "&gt;",  x, fixed = TRUE)
    x
  }

  safe_field <- function(df, row, col, fallback = "\u2014") {
    if (col %in% names(df) && !is.na(df[[col]][row])) as.character(df[[col]][row]) else fallback
  }

  # Normalize miRanda alignment rows for display. The query/reference lines
  # include labels plus flanking 3'/5' annotations, while the match line in our
  # CSV is just the base-pair markers. We therefore strip only the labels from
  # the query/reference rows, then pad the marker line to the first/last base.
  realign_display <- function(query_aln, match_str, ref_aln) {
    strip_label <- function(x, label) {
      sub(paste0("^\\s*", label, ":\\s*"), "", x)
    }

    extract_core <- function(x) {
      x %>%
        sub("^\\s*[35]'\\s*", "", .) %>%
        sub("\\s*[35]'\\s*$", "", .)
    }

    edge_lowercase_width <- function(x, side = c("left", "right")) {
      side <- match.arg(side)
      pattern <- if (side == "left") "^[a-z]+" else "[a-z]+$"
      match <- regmatches(x, regexpr(pattern, x))

      if (length(match) == 0L || identical(match, character(0)) || identical(match, "")) {
        return(0L)
      }

      nchar(match, type = "chars")
    }

    leading_space_width <- function(x) {
      match <- regmatches(x, regexpr("^\\s*", x))

      if (length(match) == 0L || identical(match, character(0)) || is.na(match)) {
        return(0L)
      }

      nchar(match, type = "chars")
    }

    trailing_space_width <- function(x) {
      match <- regmatches(x, regexpr("\\s*$", x))

      if (length(match) == 0L || identical(match, character(0)) || is.na(match)) {
        return(0L)
      }

      nchar(match, type = "chars")
    }

    normalize_match_string <- function(
      match_str,
      display_width,
      label_prefix_width,
      desired_left_pad,
      desired_right_pad
    ) {
      if (is.na(match_str)) {
        return("")
      }

      lead_spaces <- leading_space_width(match_str)
      if (label_prefix_width > 0L && lead_spaces >= label_prefix_width) {
        match_str <- substring(match_str, label_prefix_width + 1L)
      }

      match_width <- nchar(match_str, type = "chars")
      if (match_width == display_width) {
        return(match_str)
      }

      lead_spaces <- leading_space_width(match_str)
      trail_spaces <- trailing_space_width(match_str)

      if (match_width > display_width) {
        extra_width <- match_width - display_width
        trim_left <- min(extra_width, max(0L, lead_spaces - desired_left_pad))
        extra_width <- extra_width - trim_left

        trim_right <- min(extra_width, max(0L, trail_spaces - desired_right_pad))
        extra_width <- extra_width - trim_right

        if (trim_left > 0L) {
          match_str <- substring(match_str, trim_left + 1L)
        }

        if (trim_right > 0L) {
          match_str <- substring(
            match_str,
            1L,
            nchar(match_str, type = "chars") - trim_right
          )
        }

        return(match_str)
      }

      missing_width <- display_width - match_width
      left_missing <- max(0L, desired_left_pad - lead_spaces)
      right_missing <- max(0L, desired_right_pad - trail_spaces)
      left_pad <- min(left_missing, missing_width)
      right_pad <- min(right_missing, missing_width - left_pad)
      remaining_pad <- missing_width - left_pad - right_pad

      paste0(
        strrep(" ", left_pad + remaining_pad),
        match_str,
        strrep(" ", right_pad)
      )
    }

    base_bounds <- function(x) {
      base_pos <- gregexpr("[A-Za-z-]", x)[[1]]
      if (length(base_pos) == 1L && base_pos[1] == -1L) {
        return(list(left = 0L, right = 0L))
      }

      first_base <- base_pos[1]
      last_base <- base_pos[length(base_pos)]

      list(
        left = first_base - 1L,
        right = nchar(x, type = "chars") - last_base
      )
    }

    q_seq <- strip_label(query_aln, "Query")
    r_seq <- strip_label(ref_aln, "Ref")
    q_core <- extract_core(q_seq)
    r_core <- extract_core(r_seq)

    q_bounds <- base_bounds(q_seq)
    r_bounds <- base_bounds(r_seq)

    left_pad <- max(q_bounds$left, r_bounds$left)
    right_pad <- max(q_bounds$right, r_bounds$right)
    overhang_left <- max(
      edge_lowercase_width(q_core, "left"),
      edge_lowercase_width(r_core, "left")
    )
    overhang_right <- max(
      edge_lowercase_width(q_core, "right"),
      edge_lowercase_width(r_core, "right")
    )
    display_width <- max(
      nchar(q_seq, type = "chars"),
      nchar(r_seq, type = "chars")
    )
    label_prefix_width <- max(
      nchar(query_aln, type = "chars") - nchar(q_seq, type = "chars"),
      nchar(ref_aln, type = "chars") - nchar(r_seq, type = "chars")
    )
    match_str <- normalize_match_string(
      match_str = match_str,
      display_width = display_width,
      label_prefix_width = label_prefix_width,
      desired_left_pad = left_pad + overhang_left,
      desired_right_pad = right_pad + overhang_right
    )

    list(
      query = paste0("tRF:      ", q_seq),
      match = paste0("          ", match_str),
      ref   = paste0("target:   ", r_seq)
    )
  }

  if (nrow(df) > 200000) {
    warning("Static self-contained HTML reports with >200,000 rows will be very large and may still be slow in a browser.")
  }

  report_rows <- vector("list", nrow(df))
  if (nrow(df) > 0) {
    message("Serializing report rows for HTML...")
    pb <- utils::txtProgressBar(min = 0, max = nrow(df), style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_len(nrow(df))) {
    gene <- safe_field(df, i, "gencode_gene_name")
    anticodon <- safe_field(df, i, "tRNA_anticodon", "")
    gc_loc <- safe_field(df, i, "gencode_location")
    st_loc <- safe_field(df, i, "stringtie_location")
    imprint_status <- if ("imprint_status" %in% names(df) && !is.na(df$imprint_status[i])) {
      as.character(df$imprint_status[i])
    } else {
      NA_character_
    }
    ltr_info <- if ("LTR" %in% names(df) && !is.na(df$LTR[i]) && df$LTR[i]) {
      paste0(safe_field(df, i, "LTR_family"), " / ", safe_field(df, i, "LTR_gene_id"))
    } else {
      "None"
    }

    aln <- realign_display(
      safe_field(df, i, "query_alignment", "N/A"),
      safe_field(df, i, "match_string", ""),
      safe_field(df, i, "ref_alignment", "N/A")
    )

    loc_text <- paste0(
      df$seqnames[i], ":",
      format(df$start[i], big.mark = ","),
      "-",
      format(df$end[i], big.mark = ","),
      " (", df$strand[i], ")"
    )

    report_rows[[i]] <- list(
      index = i,
      score = as.numeric(df$alignment_score[i]),
      trf = safe_field(df, i, "tRF", "N/A"),
      anticodon = if (nzchar(anticodon)) anticodon else "Anticodon N/A",
      gene = gene,
      location = loc_text,
      energy = as.character(df$energy[i]),
      has_ltr = identical(ltr_info != "None", TRUE),
      ltr_info = ltr_info,
      imprint_status = imprint_status,
      gencode_location = gc_loc,
      stringtie_location = st_loc,
      query = aln$query,
      match = aln$match,
      ref = aln$ref,
      search_text = tolower(
        paste(
          safe_field(df, i, "tRF", "N/A"),
          if (nzchar(anticodon)) anticodon else "",
          gene,
          loc_text,
          ltr_info,
          gc_loc,
          st_loc,
          if (is.na(imprint_status)) "" else imprint_status,
          sep = " "
        )
      )
    )

    if (nrow(df) > 0 && (i == 1L || i %% 100L == 0L || i == nrow(df))) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  rows_json <- jsonlite::toJSON(
    report_rows,
    auto_unbox = TRUE,
    null = "null",
    na = "null"
  )
  rows_json <- gsub("</", "<\\/", rows_json, fixed = TRUE)

  html <- paste0(
'<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>', escape_html(title), '</title>
<style>
  * { box-sizing: border-box; }
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    max-width: 1200px; margin: 0 auto; padding: 20px; background: #f6f8fa; color: #24292f;
  }
  h1 { margin-bottom: 4px; }
  .summary { color: #57606a; margin-bottom: 16px; }
  .filter-bar, .page-bar {
    margin-bottom: 16px; display: flex; gap: 8px; flex-wrap: wrap; align-items: center;
  }
  .filter-bar input, .filter-bar select, .page-bar select {
    padding: 8px 12px; border: 1px solid #d0d7de; border-radius: 6px; font-size: 14px; background: white;
  }
  .filter-bar input { width: 350px; }
  .page-bar button {
    padding: 8px 12px; border: 1px solid #d0d7de; border-radius: 6px; background: white; cursor: pointer;
    font-size: 14px;
  }
  .page-bar button:disabled { cursor: default; opacity: 0.5; }
  .page-info { color: #57606a; font-size: 14px; }
  .card {
    background: white; border: 1px solid #d0d7de; border-radius: 6px;
    padding: 16px; margin: 8px 0;
  }
  .meta { margin-bottom: 8px; line-height: 1.8; }
  .meta span { display: inline-block; margin-right: 16px; font-size: 13px; color: #57606a; }
  .trf { font-weight: 600; color: #0550ae !important; }
  .anticodon { font-weight: 600; color: #9a6700 !important; }
  .gene { font-weight: 600; color: #1a7f37 !important; }
  .ltr-yes { color: #8250df !important; font-weight: 600; }
  .ltr-no { color: #9a9a9a !important; }
  .imprinted-yes { color: #b42318 !important; font-weight: 600; }
  .imprinted-no { color: #9a9a9a !important; }
  .alignment {
    background: #f6f8fa; padding: 12px; border-radius: 4px;
    font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
    font-size: 13px; line-height: 1.6; overflow-x: auto; margin: 0;
  }
  .wc { color: #1a7f37; font-weight: bold; }
  .gu { color: #bf8700; font-weight: bold; }
  .mm { color: #cf222e; }
  .count { font-weight: 600; }
  .empty {
    background: white; border: 1px dashed #d0d7de; border-radius: 6px;
    padding: 24px; color: #57606a;
  }
</style>
</head>
<body>
<h1>', escape_html(title), '</h1>
<p class="summary"><span class="count" id="visible-count">', nrow(df), '</span> of <span class="count">', nrow(df), '</span> target sites match current filters</p>
<div class="filter-bar">
  <input type="text" id="search" placeholder="Search gene, tRF, coordinates, LTR family...">
  <select id="ltr-filter">
    <option value="all">All sites</option>
    <option value="ltr">LTR-associated only</option>
    <option value="no-ltr">Non-LTR only</option>
  </select>
  <select id="imprinted-filter">
    <option value="all">All genes</option>
    <option value="imprinted">Imprinted genes only</option>
    <option value="non-imprinted">Non-imprinted only</option>
  </select>
  <select id="sort-order">
    <option value="original">Original order</option>
    <option value="score-desc">Alignment score: high to low</option>
    <option value="score-asc">Alignment score: low to high</option>
  </select>
</div>
<div class="page-bar">
  <button id="prev-page" type="button">Previous</button>
  <button id="next-page" type="button">Next</button>
  <select id="page-size">
    <option value="50">50 per page</option>
    <option value="100" selected>100 per page</option>
    <option value="250">250 per page</option>
    <option value="500">500 per page</option>
  </select>
  <span class="page-info" id="page-info"></span>
</div>
<div id="cards"></div>
<script>
const rows = ', rows_json, ';
const state = {
  query: "",
  ltr: "all",
  imprinted: "all",
  sortOrder: "original",
  page: 1,
  pageSize: 100
};

let filteredRows = rows.slice();

function escapeHtml(value) {
  if (value === null || value === undefined) return "";
  return String(value)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

function colorMatchString(matchStr) {
  if (matchStr === null || matchStr === undefined) return "N/A";
  return Array.from(String(matchStr)).map(function(ch) {
    if (ch === "|") return \'<span class="wc">|</span>\';
    if (ch === ":") return \'<span class="gu">:</span>\';
    if (ch === " ") return " ";
    return \'<span class="mm">\' + escapeHtml(ch) + \'</span>\';
  }).join("");
}

function debounce(fn, delayMs) {
  let timer = null;
  return function() {
    const args = arguments;
    clearTimeout(timer);
    timer = setTimeout(function() {
      fn.apply(null, args);
    }, delayMs);
  };
}

function applyFiltersAndSort(resetPage) {
  const query = state.query.trim().toLowerCase();

  filteredRows = rows.filter(function(row) {
    const matchesSearch = !query || row.search_text.indexOf(query) !== -1;
    const matchesLtr = state.ltr === "all" ||
      (state.ltr === "ltr" && row.has_ltr) ||
      (state.ltr === "no-ltr" && !row.has_ltr);

    const hasImprintStatus = row.imprint_status !== null && row.imprint_status !== undefined;
    const isImprinted = row.imprint_status === "Imprinted";
    const matchesImprinted = state.imprinted === "all" ||
      (state.imprinted === "imprinted" && isImprinted) ||
      (state.imprinted === "non-imprinted" && hasImprintStatus && !isImprinted);

    return matchesSearch && matchesLtr && matchesImprinted;
  });

  filteredRows.sort(function(a, b) {
    if (state.sortOrder === "score-desc") {
      return b.score - a.score;
    }
    if (state.sortOrder === "score-asc") {
      return a.score - b.score;
    }
    return a.index - b.index;
  });

  if (resetPage) {
    state.page = 1;
  }

  const maxPage = Math.max(1, Math.ceil(filteredRows.length / state.pageSize));
  if (state.page > maxPage) {
    state.page = maxPage;
  }
}

function renderCards() {
  const cardsContainer = document.getElementById("cards");
  const visibleCount = document.getElementById("visible-count");
  const pageInfo = document.getElementById("page-info");
  const prevButton = document.getElementById("prev-page");
  const nextButton = document.getElementById("next-page");

  const total = filteredRows.length;
  const totalPages = Math.max(1, Math.ceil(total / state.pageSize));
  const startIdx = total === 0 ? 0 : (state.page - 1) * state.pageSize;
  const endIdx = Math.min(startIdx + state.pageSize, total);
  const pageRows = filteredRows.slice(startIdx, endIdx);

  if (pageRows.length === 0) {
    cardsContainer.innerHTML = \'<div class="empty">No target sites match the current filters.</div>\';
  } else {
    cardsContainer.innerHTML = pageRows.map(function(row) {
      const imprintTag = row.imprint_status
        ? \'<span class="imprinted-tag \' + (row.imprint_status === "Imprinted" ? "imprinted-yes" : "imprinted-no") + \'">Imprint: \' + escapeHtml(row.imprint_status) + \'</span>\'
        : "";
      const ltrClass = row.has_ltr ? "ltr-yes" : "ltr-no";

      return \'<div class="card">\' +
        \'<div class="meta">\' +
          \'<span class="trf">\' + escapeHtml(row.trf) + \'</span>\' +
          \'<span class="anticodon">\' + escapeHtml(row.anticodon) + \'</span>\' +
          \'<span class="gene">\' + escapeHtml(row.gene) + \'</span>\' +
          \'<span class="loc">\' + escapeHtml(row.location) + \'</span>\' +
          \'<span class="score">Score: \' + escapeHtml(row.score) + \' | Energy: \' + escapeHtml(row.energy) + \'</span>\' +
          \'<span class="ltr-tag \' + ltrClass + \'">LTR: \' + escapeHtml(row.ltr_info) + \'</span>\' +
          imprintTag +
          \'<span class="region">GENCODE: \' + escapeHtml(row.gencode_location) + \' | StringTie: \' + escapeHtml(row.stringtie_location) + \'</span>\' +
        \'</div>\' +
        \'<pre class="alignment">\' +
          escapeHtml(row.query) + "\\n" +
          colorMatchString(row.match) + "\\n" +
          escapeHtml(row.ref) +
        \'</pre>\' +
      \'</div>\';
    }).join("\\n");
  }

  visibleCount.textContent = total;
  pageInfo.textContent = total === 0
    ? "No results"
    : "Showing " + (startIdx + 1) + "-" + endIdx + " of " + total + " matching target sites (page " + state.page + " of " + totalPages + ")";

  prevButton.disabled = state.page <= 1;
  nextButton.disabled = state.page >= totalPages;
}

function refresh(resetPage) {
  applyFiltersAndSort(resetPage);
  renderCards();
}

document.getElementById("search").addEventListener("input", debounce(function(event) {
  state.query = event.target.value;
  refresh(true);
}, 120));

document.getElementById("ltr-filter").addEventListener("change", function(event) {
  state.ltr = event.target.value;
  refresh(true);
});

document.getElementById("imprinted-filter").addEventListener("change", function(event) {
  state.imprinted = event.target.value;
  refresh(true);
});

document.getElementById("sort-order").addEventListener("change", function(event) {
  state.sortOrder = event.target.value;
  refresh(true);
});

document.getElementById("page-size").addEventListener("change", function(event) {
  state.pageSize = parseInt(event.target.value, 10);
  refresh(true);
});

document.getElementById("prev-page").addEventListener("click", function() {
  if (state.page > 1) {
    state.page -= 1;
    renderCards();
  }
});

document.getElementById("next-page").addEventListener("click", function() {
  const totalPages = Math.max(1, Math.ceil(filteredRows.length / state.pageSize));
  if (state.page < totalPages) {
    state.page += 1;
    renderCards();
  }
});

refresh(true);
</script>
</body>
</html>')

  writeLines(html, output_file)
  message("Alignment report written to: ", output_file)
}

# Generate reports for key subsets
if (all(c("query_alignment", "match_string", "ref_alignment") %in% names(unique_data))) {

  generate_alignment_report(
    unique_data,
    output_file = 'import/miranda/alignment_report_unique.html',
    title       = "tRF Target Site Alignment Browser — Unique Target Sites"
  )
}
