
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

expected = input$n

tRF_targets = unique_data %>%
  dplyr::filter(alignment_score >= 75, LTR == TRUE, !is.na(gencode_transcripts), gencode_location != "Intergenic") %>%
  group_by(gencode_location) %>%
  summarize(n = n()) %>%
  ungroup()

observed = tRF_targets$n

names(expected) = input$location
names(observed) = tRF_targets$gencode_location

# Ensure ordering is consistent
expected = expected[order(names(expected))]
observed = observed[order(names(observed))]

# Chi-square test
chisq.test(x = observed, p = expected / sum(expected))



obs_total = sum(observed)
exp_total = sum(expected)

# Loop over each region
results = map2_dfr(names(observed), observed, ~{
  region = .x
  obs_in = .y
  obs_out = obs_total - obs_in
  exp_in = expected[region]
  exp_out = exp_total - exp_in
  
  mat = matrix(c(obs_in, obs_out, exp_in, exp_out), nrow = 2, byrow = TRUE,
               dimnames = list(Group = c("Observed", "Expected"),
                               Region = c("InRegion", "OutRegion")))
  ft = fisher.test(mat)
  
  tibble(
    Region = region,
    log2_OR = log2(ft$estimate),
    log2_CI_low = log2(ft$conf.int[1]),
    log2_CI_high = log2(ft$conf.int[2]),
    p_value = ft$p.value
  )
})

# Adjust p-values
results = results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

print(results)

plot_data = results %>%
  mutate(
    odds_ratio = 2^log2_OR,
    ci_low = 2^log2_CI_low,
    ci_high = 2^log2_CI_high,
    p_label = ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj))
  ) %>%
  mutate(Region = factor(Region,
                           levels = c("Exon - non-coding", "Exon - 5' UTR", "Exon - CDS", "Exon - 3' UTR"),
                           labels = c("lncRNA", "5' UTR", "CDS", "3' UTR")))

plot = ggplot(plot_data, aes(x = Region, y = odds_ratio)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(aes(label = p_label, y = ci_high + 0.15), size = 3, color = "black") +
  scale_y_continuous(name = "Odds ratio", expand = expansion(mult = c(0.05, 0.1))) +
  xlab('')


pdf("~/PhD/Manuscripts/Domesticated_retrotransposon_targets/figure_1F.pdf",
    width = 2, height = 2.25)

plot + theme_bw() + theme(plot.title = element_text(size = 10),
                          axis.text.x = element_text(size = 8, color = 'black', angle = 45, vjust = 0.5, hjust = 0.5),
                          axis.text.y = element_text(size = 8, color = 'black'),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                          axis.line = element_line(size = stroke),
                          panel.border = element_blank(),
                          legend.text = element_text(size = 8),
                          legend.title = element_text(size = 10),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          strip.text.x = element_text(
                            size = 10),
                          strip.background = element_rect(
                            color="white", fill="white", size=1.5, linetype="solid"))

dev.off()
