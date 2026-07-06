## Adds tDRnamer names to the tRF CSVs produced by generate_tRF_fasta.R.
##
## Run this AFTER tDRnamer, e.g.:
##   Rscript miranda/generate_tRF_fasta.R
##   tDRnamer --mode seq --seq import/{genome}_tRF3b.fasta \
##            --db .../reference_database/{genome}/{genome} \
##            -o import/tDRnamer/{genome}
##   Rscript miranda/add_tDRnamer_names.R {genome} {tRF_type}
##
## It reads the tDRnamer output fasta (which carries BOTH the original tRF_N
## name and the new tDR- name in each header), builds a tRF_N -> tDR lookup,
## left-joins it onto import/{genome}_{tRF_type}.csv on the unique_name key,
## and rewrites the CSV with a new `tDRnamer` column.

suppressMessages({
  library(Biostrings)
  library(glue)
  library(stringr)
  library(dplyr)
  library(readr)
})

## ---- Arguments
args     <- commandArgs(trailingOnly = TRUE)
genome   <- if (length(args) >= 1) args[[1]] else "mm10"
tRF_type <- if (length(args) >= 2) args[[2]] else "tRF3b"   # or "tRF3a"

csv_in   <- glue("import/{genome}_{tRF_type}.csv")

## tDRnamer (-o import/tDRnamer/{genome}) writes {genome}-tDR.fa, the named
## fasta whose headers carry both the original tRF_N name and the new tDR- name.
tdr_dir  <- glue("import/tDRnamer")
tdr_fa   <- glue("{tdr_dir}/{genome}-tDR.fa")
if (!file.exists(tdr_fa))
  stop(glue("Expected tDRnamer fasta not found: {tdr_fa}"))

message(glue("Reading CSV : {csv_in}"))
message(glue("Reading fasta: {tdr_fa}"))

## ---- Build tRF_N -> tDR lookup from the fasta headers -------------------------
headers <- names(readBStringSet(tdr_fa))          # full header text, ">" stripped

# The two names sit in the same header separated by "_" or whitespace, in either
# order. Pull each out by its own pattern so ordering doesn't matter.
key <- str_extract(headers, "tRF_[0-9]+")         # original unique_name
tdr <- str_extract(headers, "tDR-[^|_[:space:]]+") # tDRnamer name (stop at | _ or space)

lookup <- tibble(unique_name = key, tDRnamer = tdr) %>%
  filter(!is.na(unique_name)) %>%
  distinct(unique_name, .keep_all = TRUE)

## ---- Join and write back -----------------------------------------------------
df <- read_csv(csv_in, show_col_types = FALSE)

df <- df %>%
  select(-any_of("tDRnamer")) %>%      # drop a stale column if re-running
  left_join(lookup, by = "unique_name")

missing <- sum(is.na(df$tDRnamer))
if (missing > 0)
  warning(glue("{missing} of {nrow(df)} tRFs had no tDRnamer match."))

write_csv(df, csv_in)
message(glue("Wrote {csv_in} with tDRnamer column ({nrow(df) - missing}/{nrow(df)} named)."))
