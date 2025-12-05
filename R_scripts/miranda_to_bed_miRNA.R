################################################################################
# miranda_to_bed.R
# Combines miRanda summary files into a single output with genomic coordinates
# of hit sites. Writes to .csv and .bed format.
################################################################################

library(tidyverse)
library(glue)
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)

args = commandArgs(TRUE)

summary_directory   = args[1]
output_directory    = args[2]
score_cutoff        = as.numeric(args[3])
utr5_annotation_rds = args[4]  # path to saved 5'UTR GRanges RDS

cat("Summary directory:", summary_directory,   "\n")
cat("Output directory :", output_directory,    "\n")
cat("Score cutoff     :", score_cutoff,        "\n")
cat("5'UTR annotation :", utr5_annotation_rds, "\n")

# Load 5'UTR GRanges ONCE
utr5_gr <- readRDS(utr5_annotation_rds)
GenomeInfoDb::seqlevelsStyle(utr5_gr) <- "UCSC"

# Where we'll write the combined outputs
csv_out <- glue("{output_directory}/miranda_output_miRNA_5UTR_unique_{score_cutoff}.csv")
bed_out <- glue("{output_directory}/miranda_output_miRNA_5UTR_unique_{score_cutoff}.bed")

# If re-running, remove old outputs so we don't append to stale files
if (file.exists(csv_out)) file.remove(csv_out)
if (file.exists(bed_out)) file.remove(bed_out)

# List all miRanda summary files
filenames <- list.files(
  summary_directory,
  pattern   = "^summary",
  full.names = TRUE
)

cat("Found", length(filenames), "summary files\n")

process_one_file <- function(file, first_file = FALSE) {
  message("Processing: ", file)

  # Read as character
  temp <- readr::read_tsv(
    file,
    col_names = FALSE,
    col_types = cols(.default = col_character()),
    progress  = FALSE
  )

  if (nrow(temp) == 0) {
    message("  -> empty file, skipping")
    return(invisible(NULL))
  }

  if (!"X11" %in% names(temp)) {
    message("  -> file has fewer than 11 columns (no valid hits?), skipping")
    return(invisible(NULL))
  }

  temp <- temp %>% filter(!is.na(X11))

  # Drop X9, X10 early
  if (all(c("X9", "X10") %in% names(temp))) {
    temp <- temp[, !names(temp) %in% c("X9", "X10")]
  }

  # Rename columns
  df <- temp %>%
    dplyr::rename(
      tRF             = X1,
      coordinates     = X2,
      alignment_score = X3,
      energy          = X4,
      Z_score         = X5,
      miRNA_position  = X6,
      target_position = X7,
      alignment_length= X8,
      strand          = X11
    ) %>%
    mutate(
      alignment_score = suppressWarnings(as.numeric(alignment_score))
    ) %>%
    filter(!is.na(alignment_score),
           alignment_score >= score_cutoff)

  if (nrow(df) == 0) {
    message("  -> no rows passing score cutoff, skipping")
    return(invisible(NULL))
  }

  # Parse target_position into start/end
  df <- df %>%
    tidyr::separate(
      target_position,
      into    = c("start", "end"),
      convert = TRUE,
      remove  = TRUE
    ) %>%
    filter(!is.na(start), !is.na(end))

  # Parse coordinates: "something::chr:start-end"
  df <- df %>%
    tidyr::separate(coordinates, into = c("A", "B"), sep = "::") %>%
    tidyr::separate(B, into = c("seqnames", "C"), sep = ":") %>%
    tidyr::separate(
      C,
      into    = c("window_start", "window_end"),
      sep     = "-",
      convert = TRUE
    )

  # Compute genomic_start / genomic_end
  df <- df %>%
    mutate(
      genomic_start = if_else(
        strand == "+",
        window_start + start,
        window_end   - end   + 1L
      ),
      genomic_end   = if_else(
        strand == "+",
        window_start + end,
        window_end   - start + 1L
      )
    ) %>%
    filter(!is.na(genomic_start), !is.na(genomic_end)) %>%
    distinct(tRF, seqnames, genomic_start, genomic_end, strand, .keep_all = TRUE)

  if (nrow(df) == 0) {
    message("  -> all rows dropped after coordinate filtering")
    return(invisible(NULL))
  }

  # Restrict to canonical chromosomes
  canon <- paste0("chr", c(1:19, "X", "Y"))
  df <- df %>% filter(seqnames %in% canon)

  if (nrow(df) == 0) {
    message("  -> no hits on canonical chromosomes, skipping")
    return(invisible(NULL))
  }

  # Build GRanges for this miRNA's hits
  gr_hits <- GenomicRanges::GRanges(
    seqnames = df$seqnames,
    ranges   = IRanges::IRanges(start = df$genomic_start, end = df$genomic_end),
    strand   = df$strand
  )
  GenomeInfoDb::seqlevelsStyle(gr_hits) <- "UCSC"

  # 5' UTR overlaps
  hits_5utr <- IRanges::overlapsAny(gr_hits, utr5_gr, ignore.strand = FALSE)

  if (!any(hits_5utr)) {
    message("  -> no 5' UTR overlaps, skipping")
    return(invisible(NULL))
  }

  df <- df[hits_5utr, , drop = FALSE]
  gr_hits <- gr_hits[hits_5utr]

  if (nrow(df) == 0) {
    message("  -> all rows dropped after 5' UTR filtering")
    return(invisible(NULL))
  }

  ## ---------- COLLAPSE OVERLAPPING HITS -------------------

  # Reduce overlapping regions
  reduced_gr <- GenomicRanges::reduce(gr_hits)

  # Map each reduced region to original hits
  ov <- GenomicRanges::findOverlaps(reduced_gr, gr_hits)
  hits_df <- tibble::as_tibble(ov) %>%
    dplyr::mutate(
      alignment_score = df$alignment_score[subjectHits]
    )

  # Best-scoring hit per reduced region
  best_hits <- hits_df %>%
    dplyr::group_by(queryHits) %>%
    dplyr::slice_max(order_by = alignment_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  # Count number of raw hits collapsed into each reduced region
  collapsed_counts <- hits_df %>%
    dplyr::count(queryHits, name = "collapsed_hit_count")

  best_hits <- dplyr::left_join(best_hits, collapsed_counts, by = "queryHits")

  # Subset df to those best hits and update coordinates to reduced regions
  df_collapsed <- df[best_hits$subjectHits, , drop = FALSE]

  df_collapsed$genomic_start <- GenomicRanges::start(reduced_gr[best_hits$queryHits])
  df_collapsed$genomic_end   <- GenomicRanges::end  (reduced_gr[best_hits$queryHits])
  df_collapsed$collapsed_hit_count <- best_hits$collapsed_hit_count

  df <- df_collapsed

  if (nrow(df) == 0) {
    message("  -> no rows left after collapsing, skipping")
    return(invisible(NULL))
  }

  ## ---------- CSV / BED output -----

  csv_df <- df %>%
    select(
      tRF, seqnames, alignment_score, energy,
      genomic_start, genomic_end, alignment_length, strand,
      collapsed_hit_count
    )

  readr::write_csv(
    csv_df,
    file      = csv_out,
    append    = !first_file,
    col_names = first_file
  )

  bed_df <- csv_df %>%
    transmute(
      chrom      = seqnames,
      chromStart = genomic_start - 1L,  # 1-based -> 0-based
      chromEnd   = genomic_end,         # half-open end
      name       = tRF,
      score      = as.integer(round(alignment_score)),
      strand     = strand
    )

  readr::write_tsv(
    bed_df,
    file      = bed_out,
    append    = !first_file,
    col_names = FALSE
  )

  rm(temp, df, csv_df, bed_df, gr_hits, reduced_gr, hits_df, best_hits, collapsed_counts)
  gc()

  invisible(NULL)
}

first <- TRUE
for (f in filenames) {
  process_one_file(f, first_file = first)
  if (first) first <- FALSE
}

