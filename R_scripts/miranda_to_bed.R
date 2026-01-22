################################################################################
# miranda_to_bed.R
# Combines miRanda summary files into a single output with genomic coordinates
# of hit sites. Writes to .csv and .bed format.
################################################################################

library(tidyverse)
library(glue)

args = commandArgs(TRUE)

summary_directory = args[1]
output_directory = args[2]
score_cutoff = as.numeric(args[3])

cat("Summary directory:", summary_directory, "\n")
cat("Output directory :", output_directory,  "\n")
cat("Score cutoff     :", score_cutoff,      "\n")

# Where we'll write the combined outputs
csv_out <- glue("{output_directory}/miranda_output_{score_cutoff}.csv")
bed_out <- glue("{output_directory}/miranda_output_{score_cutoff}.bed")

# If re-running, remove old outputs so we don't append to stale files
if (file.exists(csv_out)) file.remove(csv_out)
if (file.exists(bed_out)) file.remove(bed_out)

# List all miRanda summary files
filenames <- list.files(
  summary_directory,
  pattern = "^summary",
  full.names = TRUE
)

cat("Found", length(filenames), "summary files\n")

process_one_file <- function(file, first_file = FALSE) {
  message("Processing: ", file)

  # Read as character; using readr here is typically faster than base read.csv
  temp <- readr::read_tsv(
    file,
    col_names = FALSE,
    col_types = cols(.default = col_character()),
    progress = FALSE
  )
  
  if (nrow(temp) == 0) {
  message("  -> empty file, skipping")
  return(invisible(NULL))}

  # If there are header/summary lines with too few columns, X11 will be NA.
  # Keep only rows that have something in X11.
  
  if (!"X11" %in% names(temp)) {
   message("  -> file has fewer than 11 columns (no valid hits?), skipping")
   return(invisible(NULL))}

  temp <- temp %>% filter(!is.na(X11))

  # Drop V9, V10 early
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
    # Coerce scores to numeric and filter early
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
  # Use default separator behaviour like your original script, then drop incomplete rows
  df <- df %>%
    tidyr::separate(
      target_position,
      into    = c("start", "end"),
      convert = TRUE,
      remove  = TRUE
    ) %>%
    # keep only rows where both start and end are present
    dplyr::filter(!is.na(start), !is.na(end))

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
    # Drop broken rows
    filter(!is.na(genomic_start), !is.na(genomic_end)) %>%
    # Remove duplicates *within this file*
    distinct(tRF, seqnames, genomic_start, genomic_end, strand, .keep_all = TRUE)

  if (nrow(df) == 0) {
    message("  -> all rows dropped after coordinate filtering")
    return(invisible(NULL))
  }

  # ---------------- CSV output ----------------

  csv_df <- df %>%
    select(
      tRF, seqnames, alignment_score, energy,
      genomic_start, genomic_end, alignment_length, strand
    )

  readr::write_csv(
    csv_df,
    file      = csv_out,
    append    = !first_file,
    col_names = first_file
  )

  # ---------------- BED output ----------------

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

  # Clean up
  rm(temp, df, csv_df, bed_df)
  gc()

  invisible(NULL)
}

first <- TRUE
for (f in filenames) {
  process_one_file(f, first_file = first)
  if (first) first <- FALSE
}

