################################################################################
# miranda_unique_annotated_to_bed.R
# Convert miranda_output_unique_annotated.csv (or similar annotated CSV)
# to BED format.
# Output columns:
#   chrom, chromStart, chromEnd, name, score, strand, tRNA_anticodon
# where:
#   - name is unique tRF ID when available (e.g., tRF_1)
#   - tRNA_anticodon has comma-delimited values converted to "-"
################################################################################

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

input_csv <- ifelse(
  length(args) >= 1,
  args[1],
  "import/miranda/miranda_output_unique_annotated.csv"
)

output_bed <- ifelse(
  length(args) >= 2,
  args[2],
  "import/miranda/miranda_output_unique_annotated.bed"
)

cat("Input CSV :", input_csv, "\n")
cat("Output BED:", output_bed, "\n")

data <- readr::read_csv(input_csv, show_col_types = FALSE)

pick_first_col <- function(df, candidates) {
  hits <- candidates[candidates %in% names(df)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

chrom_col <- pick_first_col(data, c("seqnames", "chrom", "chr"))
start_col <- pick_first_col(data, c("start"))
end_col <- pick_first_col(data, c("end"))
coord_col <- pick_first_col(data, c("coordinates", "coordinate", "coords"))
strand_col <- pick_first_col(data, c("strand"))
anticodon_col <- pick_first_col(data, c("tRNA_anticodon", "anticodon", "tRNA_anticodons"))
trf_col <- pick_first_col(data, c("tRF", "trf", "unique_tRF", "unique_name", "tRF_ID"))
score_col <- pick_first_col(data, c("alignment_score", "score"))

if (is.na(anticodon_col)) {
  stop("Missing anticodon column. Expected one of: tRNA_anticodon, anticodon, tRNA_anticodons")
}

if (is.na(strand_col)) {
  stop("Missing strand column.")
}

if (is.na(chrom_col) || is.na(start_col) || is.na(end_col)) {
  if (is.na(coord_col)) {
    stop(
      "Missing coordinate columns. Provide either seqnames/start/end ",
      "or a coordinate column (coordinates/coordinate/coords)."
    )
  }

  parsed_coords <- data %>%
    transmute(coord_string = .data[[coord_col]]) %>%
    tidyr::extract(
      col = coord_string,
      into = c("chrom", "start", "end"),
      regex = "^([^:]+):(\\d+)-(\\d+)$",
      remove = TRUE
    )

  data <- data %>%
    bind_cols(parsed_coords)

  chrom_col <- "chrom"
  start_col <- "start"
  end_col <- "end"
}

bed_df <- data %>% mutate(
  .chrom = as.character(.data[[chrom_col]]),
  .start = suppressWarnings(as.integer(.data[[start_col]])),
  .end = suppressWarnings(as.integer(.data[[end_col]])),
  .strand = as.character(.data[[strand_col]]),
  .anticodon = as.character(.data[[anticodon_col]]),
  .name = if (!is.na(trf_col)) as.character(.data[[trf_col]]) else as.character(.data[[anticodon_col]]),
  .score = if (!is.na(score_col)) suppressWarnings(as.integer(round(as.numeric(.data[[score_col]])))) else 0L
) %>%
  transmute(
    chrom = .chrom,
    chromStart = pmax(.start - 1L, 0L), # BED uses 0-based start
    chromEnd = .end,
    name = .name,
    score = dplyr::coalesce(.score, 0L),
    strand = .strand,
    tRNA_anticodon = stringr::str_replace_all(.anticodon, ",", "-")
  ) %>%
  filter(!is.na(chrom), !is.na(chromStart), !is.na(chromEnd), !is.na(strand), !is.na(tRNA_anticodon))

readr::write_tsv(bed_df, output_bed, col_names = FALSE)

cat("Wrote", nrow(bed_df), "BED rows\n")
