library(tidyverse)
library(glue)
library(regioneR)
library(GenomicRanges)
library(GenomeInfoDb)

args = commandArgs(TRUE)

tRF_file = args[1]
miRNA_file = args[2]
five_UTR = args[3]
output_directory = args[4]

n_cores <- as.integer(Sys.getenv("NSLOTS"))
if (is.na(n_cores)) n_cores <- 1

#############################################################################
## Load data
#############################################################################

# Import universe

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

miRNA_data = read.csv(miRNA_file, header = T)

#############################################################################
## Permutation analysis
#############################################################################

results = tibble(
  tRF_cutoff = numeric(),
  miRNA_cutoff = numeric(),
  observed_overlaps = integer(),
  mean_random_overlaps = numeric(),
  z_score = numeric(),
  p_value = numeric()
)

tRF_cutoffs = c(70, 75, 80, 85, 90)
miRNA_cutoffs = c(150, 200)

set.seed(123)

for (tRF_cutoff in tRF_cutoffs){

	# Hard coding. Change later to nested loop.
	miRNA_cutoff = 150

	print(glue('Processing tRF cutoff: {tRF_cutoff}'))

	tRF_data_subset = tRF_data %>% dplyr::filter(alignment_score >= tRF_cutoff)

	if (nrow(tRF_data_subset) == 0) {
		message("No more tRF hits at cutoff ", tRF_cutoff, ". Terminating loop.")
		break
	}

        tRF_data_subset_GRanges = makeGRangesFromDataFrame(tRF_data_subset, keep.extra.columns = TRUE)

	miRNA_data_subset = miRNA_data %>% dplyr::filter(alignment_score >= miRNA_cutoff)

	if (nrow(miRNA_data_subset) == 0) {
                message("No more miRNA hits at cutoff ", miRNA_cutoff, ". Terminating loop.")
                break
        }

	miRNA_data_subset_GRanges = makeGRangesFromDataFrame(miRNA_data_subset, keep.extra.columns = TRUE)

	# Restrict both target site sets to 5' UTR

	tRF_5UTR = subsetByOverlaps(tRF_data_subset_GRanges, five_UTR_universe)
  	miRNA_5UTR = subsetByOverlaps(miRNA_data_subset_GRanges, five_UTR_universe)

	if (length(tRF_5UTR) == 0) {
    		message("No tRF hits in 5' UTR at cutoff ", tRF_cutoff, ". Skipping.")
    		break
  	}

  	if (length(miRNA_5UTR) == 0) {
    		message("No miRNA hits in 5' UTR at cutoff ", miRNA_cutoff, ". Skipping.")
    		break
  	}

	perm = permTest(A = tRF_5UTR,
    			B = miRNA_5UTR,
    			ntimes = 100,
    			mc.cores = n_cores,
    			alternative = "greater",
    			randomize.function = resampleRegions,
    			evaluate.function  = numOverlaps,
    			universe = five_UTR_universe)
	perm_results = perm$numOverlaps 

	results = results %>% add_row(
  tRF_cutoff           = tRF_cutoff,
  miRNA_cutoff         = miRNA_cutoff,
  observed_overlaps    = perm_results$observed,
  mean_random_overlaps = mean(perm_results$permuted),
  z_score              = perm_results$zscore,
  p_value              = perm_results$pval
)

}
	
write.csv(results, glue("{output_directory}/miRNA_enrichment.csv"), row.names = FALSE)

