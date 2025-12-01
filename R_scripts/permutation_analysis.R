library(tidyverse)
library(glue)
library(regioneR)

args = commandArgs(TRUE)

data_file = args[1]
rmsk_file = args[2]
min_cutoff = as.numeric(args[3])
output_directory = args[4]

n_cores <- as.integer(Sys.getenv("NSLOTS"))
if (is.na(n_cores)) n_cores <- 1

# Load data

print("Importing miRanda output...")

data = read.csv(file = data_file, header = TRUE)
print(colnames(data))

data = dplyr::filter(data, alignment_score >= min_cutoff)



print("Importing RepeatMasker annotation..")

LTR = read.csv(rmsk_file, header = T)
LTR = LTR[!grepl('int', LTR$repName),]

subject = makeGRangesFromDataFrame(LTR, 
                                   keep.extra.columns = TRUE,
                                   start.field = 'genoStart',
                                   end.field = 'genoEnd',
                                   seqnames.field = 'genoName')

end(subject[strand(subject) == '+']) = end(subject[strand(subject) == '+']) + 200
start(subject[strand(subject) == '-']) = start(subject[strand(subject) == '-']) - 200

# Permutation analysis 

results = tibble(
  cutoff = numeric(),
  observed_overlaps = integer(),
  mean_random_overlaps = numeric(),
  z_score = numeric(),
  p_value = numeric()
)

cutoffs = seq(from = min_cutoff, by = 5, length.out = 20)

print(cutoffs)

set.seed(123)

for (cutoff in cutoffs){

	print(glue('Processing cutoff: {cutoff}'))

	data_sub = data %>% dplyr::filter(alignment_score >= cutoff)

	if (nrow(data_sub) == 0) {
		message("No more hits at cutoff ", cutoff, ". Terminating loop.")
		break
	}

	query_sub = makeGRangesFromDataFrame(data_sub, keep.extra.columns = TRUE)

	perm = overlapPermTest(mc.cores = n_cores, alternative = "greater", A = subject, B = query_sub, ntimes = 100, genome = 'mm10')

	perm_results = perm$numOverlaps

	results = results %>% add_row(cutoff = cutoff, 
				      observed_overlaps = perm_results$observed, 
				      mean_random_overlaps = mean(perm_results$permuted), 
				      z_score = perm_results$zscore, 
                                      p_value = perm_results$pval)

}
	
write.csv(results, glue("{output_directory}/LTR_enrichment_by_cutoff.csv"), row.names = FALSE)

