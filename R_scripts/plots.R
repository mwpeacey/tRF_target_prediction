
LTR = rtracklayer::readGFF(file = 'import/annotation_tables/mm10_rmsk_TE.gtf') %>%
  dplyr::filter(class_id == 'LTR')

query = makeGRangesFromDataFrame(LTR, keep.extra.columns = TRUE)
end(query[strand(query) == '+']) = end(query[strand(query) == '+']) + 200
start(query[strand(query) == '-']) = start(query[strand(query) == '-']) - 200

gr_hits = query

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

LTR <- LTR %>%
  bind_cols(annot %>% dplyr::select(-hit_idx)) %>%
  mutate(
    location = case_when(
      overlaps_5utr ~ "Exon - 5' UTR",
      overlaps_3utr ~ "Exon - 3' UTR",
      overlaps_cds  ~ "Exon - CDS",
      overlaps_exon ~ "Exon - non-coding",
      TRUE          ~ "Intergenic"
    )
  )


input = dplyr::filter(LTR) %>%
  group_by(location) %>%
  summarize(n = n()) %>%
  dplyr::filter(location != 'Intergenic')

size = 0.5

plot = ggplot(input, aes(x = location, y = n)) +
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