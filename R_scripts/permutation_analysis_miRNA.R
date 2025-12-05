library(tidyverse)
library(glue)
library(regioneR)
library(GenomicRanges)
library(GenomeInfoDb)

args = commandArgs(TRUE)

tRF_file      = args[1]
miRNA_file    = args[2]
five_UTR      = args[3]
output_directory = args[4]

n_cores <- as.integer(Sys.getenv("NSLOTS"))
if (is.na(n_cores)) n_cores <- 1

#############################################################################
## Load data
#############################################################################

# Import 5' UTR universe
print("Importing universe...")

five_UTR_universe = readRDS(five_UTR)

canonical_chrs = paste0("chr", c(1:19, "X", "Y"))
five_UTR_universe = GenomeInfoDb::keepSeqlevels(
  five_UTR_universe,
  canonical_chrs,
  pruning.mode = "coarse"
)

# Import tRF prediction
print("Importing tRF prediction...")

tRF_data = read.csv(file = tRF_file, header = TRUE)
print(colnames(tRF_data))

# Import miRNA prediction
print("Importing miRNA prediction...")

miRNA_data = read.csv(miRNA_file, header = TRUE)

#############################################################################
## Permutation analysis
#############################################################################

results = tibble(
  tRF_cutoff           = numeric(),
  miRNA_cutoff         = numeric(),
  observed_overlaps    = numeric(),
  mean_random_overlaps = numeric(),
  z_score              = numeric(),
  p_value              = numeric()
)

tRF_cutoffs   = c(70, 75, 80, 85, 90)
miRNA_cutoffs = c(150, 200)

set.seed(123)

for (tRF_cutoff in tRF_cutoffs) {

  message(glue("Processing tRF cutoff: {tRF_cutoff}"))

  ## Filter tRFs by score
  tRF_data_subset = tRF_data %>%
    dplyr::filter(alignment_score >= tRF_cutoff)

  if (nrow(tRF_data_subset) == 0) {
    message("No more tRF hits at cutoff ", tRF_cutoff, ". Terminating outer loop.")
    break
  }

  tRF_data_subset_GRanges = makeGRangesFromDataFrame(
    tRF_data_subset,
    keep.extra.columns = TRUE
  )

  ## Restrict tRFs to 5' UTRs once per tRF cutoff
  tRF_5UTR = subsetByOverlaps(tRF_data_subset_GRanges, five_UTR_universe)

  if (length(tRF_5UTR) == 0) {
    message("No tRF hits in 5' UTR at cutoff ", tRF_cutoff,
            ". Terminating outer loop.")
    break
  }

  ## Inner loop: iterate over miRNA score cutoffs
  for (miRNA_cutoff in miRNA_cutoffs) {

    message(glue("  Processing miRNA cutoff: {miRNA_cutoff}"))

    miRNA_data_subset = miRNA_data %>%
      dplyr::filter(alignment_score >= miRNA_cutoff)

    if (nrow(miRNA_data_subset) == 0) {
      message("  No more miRNA hits at cutoff ", miRNA_cutoff,
              " for tRF cutoff ", tRF_cutoff, ". Breaking inner loop.")
      break  # stop trying higher miRNA cutoffs for this tRF cutoff
    }

    miRNA_data_subset_GRanges = makeGRangesFromDataFrame(
      miRNA_data_subset,
      keep.extra.columns = TRUE
    )

    ## Restrict miRNA hits to 5' UTRs
    miRNA_5UTR = subsetByOverlaps(miRNA_data_subset_GRanges, five_UTR_universe)

    if (length(miRNA_5UTR) == 0) {
      message("  No miRNA hits in 5' UTR at cutoff ", miRNA_cutoff,
              " for tRF cutoff ", tRF_cutoff, ". Breaking inner loop.")
      break
    }

    ## Permutation test restricted to 5' UTR universe
    perm = permTest(
      A                  = tRF_5UTR,
      B                  = miRNA_5UTR,
      ntimes             = 100,
      mc.cores           = n_cores,
      alternative        = "greater",  # or "less"/"auto"
      randomize.function = resampleRegions,
      evaluate.function  = numOverlaps,
      universe           = five_UTR_universe
    )

    # Handle permTestResultsList vs permTestResults
    if (inherits(perm, "permTestResultsList")) {
      perm_results = perm$numOverlaps
    } else {
      perm_results = perm
    }

    results = results %>% add_row(
      tRF_cutoff           = tRF_cutoff,
      miRNA_cutoff         = miRNA_cutoff,
      observed_overlaps    = perm_results$observed,
      mean_random_overlaps = mean(perm_results$permuted),
      z_score              = perm_results$zscore,
      p_value              = perm_results$pval
    )

  } # end miRNA_cutoff loop

} # end tRF_cutoff loop

write.csv(
  results,
  glue("{output_directory}/miRNA_enrichment.csv"),
  row.names = FALSE
)

