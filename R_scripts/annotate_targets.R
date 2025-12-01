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
library(DESeq2)

data = read_csv('import/miranda/miranda_output_annotated.csv')
unique_data = read_csv('import/miranda/miranda_output_unique_annotated.csv')

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

scan_genes = read_csv('import/annotation_tables/mouse_scan_genes.csv')
other_gag_genes = read_csv('import/annotation_tables/mouse_gag_genes.csv')

data = dplyr::mutate(data, gag_gene = case_when(gencode_gene_id %in% scan_genes$gene_id ~ T,
                                                gencode_gene_name %in% other_gag_genes$gene_name ~ T,
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

# capture the transcript ID in a metadata column
mcols(gr_exon_tx)$tx_id = rep(names(gr_exonsByTx),
                              elementNROWS(gr_exonsByTx))

# 2) map each hit to any overlapping exon → transcript
mcols(gr_hits)$hit_idx = seq_along(gr_hits)
ex_ol = findOverlaps(gr_hits, gr_exon_tx, type="within")

tx_map = tibble(
  hit_idx = mcols(gr_hits)$hit_idx[queryHits(ex_ol)],
  tx_id   = as.integer(mcols(gr_exon_tx)$tx_id[subjectHits(ex_ol)])  
) %>%
  left_join(transcripts_df %>% dplyr::select(tx_id, tx_name), by = "tx_id") %>%
  group_by(hit_idx) %>%
  summarize(transcript = paste(unique(tx_name), collapse = ";"), .groups = "drop")

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

# Extract blast match per hit if any overlapping transcript matches ORF target_id
orf_match_table = open_reading_frames_annotated_filtered %>%
  dplyr::select(target_id, blast_match)

# Join transcript IDs to ORF matches
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
## Annotate DE in Schaefer et al 
################################################################################

library(DESeq2)

# Data import

filenames = list.files('~/RNA-seq/Schaefer_2022/', pattern="*.txt", full.names=TRUE)

raw_counts = filenames %>%
  map(~ read_table(.x) %>% rename_all(tolower)) %>%
  purrr::reduce(full_join, by = "geneid")

samples = tibble(sample = colnames(raw_counts %>% dplyr::select(-geneid)), 
                 genotype = c('wt', 'wt', 'dicer', 'dicer', 'ago', 'ago', 'drosha', 'drosha'))

# DEseq2

input = raw_counts %>%
  dplyr::select(-geneid) %>%
  mutate_all(~replace(., is.na(.), 0))

rownames(input) = raw_counts$geneid

# Create DESeq2 dataset
dds = DESeqDataSetFromMatrix(countData = input, 
                             colData = samples, 
                             design = ~genotype)

# Run DESeq and filter low-count genes
dds = DESeq(dds)
dds = dds[rowSums(counts(dds)) >= 1, ]

# Extract results with contrast

counter = 1
for (condition in unique(samples$genotype)){
  
  if (condition != 'wt'){
    
    print(glue('Comparing wt to {condition}...'))
    
    results = results(dds, contrast = c('genotype', condition, 'wt'), independentFiltering = TRUE)
    
    temp_results_df = as.data.frame(results)
    
    temp_results_df$condition = condition
    temp_results_df$geneid = rownames(temp_results_df)
    
    if (counter == 1){
      
      results_df = temp_results_df
      
    }
    
    else{
      
      results_df = bind_rows(results_df, temp_results_df)
      
    }
    
    counter = counter + 1
    
    rm(temp_results_df)
    
  }
  
}

results_df = dplyr::rename(results_df, genotype = condition, gene = geneid)

data = dplyr::mutate(data, AGO = case_when(gencode_gene_id %in% dplyr::filter(results_df, genotype == 'ago', padj <= 0.1, log2FoldChange > 0.5)$gene ~ T,
                                    T ~ F))

data = dplyr::mutate(data, DROSHA = case_when(gencode_gene_id %in% dplyr::filter(results_df, genotype == 'drosha', padj <= 0.1, log2FoldChange > 0.5)$gene ~ T,
                                           T ~ F))

data = dplyr::mutate(data, DICER = case_when(gencode_gene_id %in% dplyr::filter(results_df, genotype == 'dicer', padj <= 0.1, log2FoldChange > 0.5)$gene ~ T,
                                              T ~ F))

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
#gr_hits$nearest_tx[queryHits(hits)] = first_splice_sites$tx_name[subjectHits(hits)]

unique_data = as.data.frame(gr_hits)

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

## 

query = GRanges(unique_data)

CLIP = read.table('~/tRF_targets_new/GSE139344_AGO-CLIP/narrow_peak.combined.bed') %>%
  dplyr::rename(seqnames = V1, start = V2, end = V3, p = V4, strand = V6)

subject = GRanges(CLIP)

overlap = GenomicAlignments::findOverlaps(query = query, subject = subject, type = 'any', ignore.strand = F)
overlap = as.data.frame(overlap)

unique_data  = unique_data %>%
  dplyr::mutate(AGO_peak = dplyr::row_number() %in% overlap$queryHits)

df = dplyr::filter(unique_data, AGO_peak == T, AGO == T, DROSHA == F, DICER == F, gencode_location == 'Exon - 5\' UTR')
