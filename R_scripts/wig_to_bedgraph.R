infile  <- glue("~/Downloads/mm10_FANTOM5_totalCounts_Cyp2b23.wig")
cleaned <- glue("~/Downloads/mm10_FANTOM5_totalCounts_Cyp2b23.bedGraph")

# drop UCSC 'track' line and the '#bedGraph section ...' comment lines, plus blanks
x <- readLines(infile)
x <- x[!grepl("^(track|#)", x)]
x <- x[nzchar(x)]
writeLines(x, cleaned)