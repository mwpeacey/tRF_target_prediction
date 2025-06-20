library(tidyverse)
library(glue)
library(GenomicFeatures)
library(regioneR)

args = commandArgs(TRUE)

data_file = args[1]
rmsk_file = args[2]
gtf_file = args[3]
min_cutoff = as.numeric(args[4])
output_directory = args[5]

n_cores <- as.integer(Sys.getenv("NSLOTS"))
if (is.na(n_cores)) n_cores <- 1

# Load data

data = read.csv(file = data_file, header = TRUE) %>%
  dplyr::rename(transcript_start = start, transcript_end = end, start = genome_start, end = genome_end) %>%
  dplyr::filter(alignment_score >= min_cutoff)

LTR = rtracklayer::readGFF(file = rmsk_file) %>%
  dplyr::filter(class_id == 'LTR')

subject = makeGRangesFromDataFrame(LTR, keep.extra.columns = TRUE)

GTF = makeTxDbFromGFF(file = gtf_file)

# Permutation analysis 

transcriptome_background = reduce(exons(GTF))

results = tibble()

cutoffs = seq(from = min_cutoff, by = 5, length.out = 10)

set.seed(123)

for (cutoff in cutoffs){

	print(glue('Processing cutoff: {cutoff}'))

	data_sub = data %>% filter(alignment_score >= cutoff)

	if (nrow(data_sub) == 0) {
		message("No more hits at cutoff ", cutoff, ". Terminating loop.")
		break
	}

	query_sub = makeGRangesFromDataFrame(data_sub, keep.extra.columns = TRUE)

	perm = overlapPermTest(mc.cores = n_cores, A = subject, B = query_sub, ntimes = 10, genome = 'mm10', universe=transcriptome_background)

	perm_results = perm$numOverlaps

	results = results %>% add_row(cutoff = cutoff, 
				      observed_overlaps = perm_results$observed, 
				      mean_random_overlaps = mean(perm_results$permuted), 
				      z_score = perm_results$zscore, 
                                      p_value = perm_results$pval)

}
	
write.csv(results, glue("{output_directory}/LTR_enrichment_by_cutoff.csv"), row.names = FALSE)

