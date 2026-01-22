library(AnnotationHub)
library(GenomicFeatures)
library(GenomicRanges)

ah   <- AnnotationHub()
txdb <- ah[["AH75036"]]   

gr_utr5 <- GenomicFeatures::fiveUTRsByTranscript(txdb) |> unlist()

GenomeInfoDb::seqlevelsStyle(gr_utr5) <- "UCSC"

saveRDS(gr_utr5, file = "import/annotation_tables/gencode_mm10_5utr_granges.rds")
