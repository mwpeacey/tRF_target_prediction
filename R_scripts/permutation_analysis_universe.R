################################################################################
## Universe-restricted permutation test: 3'-tRF target sites vs LTRs
##
## Companion to permutation_analysis.R. That script uses
## regioneR::overlapPermTest, which keeps the target sites fixed and shuffles
## the *LTRs* across the whole of mm10, at a series of score cutoffs.
##
## This script instead keeps the LTRs fixed and shuffles the *target sites*
## to random positions inside a user-supplied universe, at a single alignment
## score cutoff (default 80).
##
## Null model
##   Only target sites lying entirely inside the universe are tested. In each
##   permutation, every such site is moved to a uniformly random position in
##   the universe at which it fits completely, keeping its width (equivalent to
##   `bedtools shuffle -incl <universe> -f 1.0`). Placement is genome-wide
##   across the universe, i.e. sites can move between chromosomes.
##
## Statistic
##   Number of target sites overlapping >= 1 LTR (strand ignored). LTRs are
##   defined as in permutation_analysis.R: RepeatMasker entries without "int"
##   in repName, extended 200 bp at their 3' end.
##
## Usage
##   Rscript permutation_analysis_universe.R <data_file> <rmsk_file> \
##       <universe_file> <output_directory> [cutoff=80] [ntimes=1000]
##
##   data_file      miRanda CSV from miranda_to_bed.R (seqnames, genomic_start,
##                  genomic_end, alignment_score)
##   rmsk_file      RepeatMasker table as CSV (genoName, genoStart, genoEnd,
##                  strand, repName), as used by permutation_analysis.R
##   universe_file  .rds (GRanges, e.g. gencode_mm10_5utr_granges.rds) or
##                  .bed (0-based; first three columns used)
##   cutoff         minimum alignment score (>=). Default 80.
##   ntimes         number of permutations. Default 1000.
##
## Outputs (in output_directory; <label> = universe file name without extension)
##   LTR_enrichment_<label>_s<cutoff>.csv           one-row summary
##   LTR_enrichment_<label>_s<cutoff>_permuted.csv  per-permutation counts
################################################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(glue)
  library(parallel)
})

args = commandArgs(TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript permutation_analysis_universe.R <data_file> <rmsk_file> ",
       "<universe_file> <output_directory> [cutoff=80] [ntimes=1000]")
}

data_file        = args[1]
rmsk_file        = args[2]
universe_file    = args[3]
output_directory = args[4]
cutoff           = if (length(args) >= 5) as.numeric(args[5]) else 80
ntimes           = if (length(args) >= 6) as.integer(args[6]) else 1000L

n_cores <- as.integer(Sys.getenv("NSLOTS"))
if (is.na(n_cores)) n_cores <- 1

canonical_chrs = paste0("chr", c(1:19, "X", "Y"))

universe_label = tools::file_path_sans_ext(basename(universe_file))

message(glue("Score cutoff: >= {cutoff} | permutations: {ntimes} | cores: {n_cores}"))
message(glue("Universe: {universe_file}"))

#############################################################################
## Load universe
#############################################################################

message("Importing universe...")

if (grepl("\\.rds$", universe_file, ignore.case = TRUE)) {
  universe = granges(readRDS(universe_file))
} else if (grepl("\\.bed$", universe_file, ignore.case = TRUE)) {
  bed = read.delim(universe_file, header = FALSE, comment.char = "#",
                   colClasses = c("character", "numeric", "numeric"))[, 1:3]
  bed = bed[!grepl("^(track|browser)", bed[[1]]), ]
  universe = GRanges(bed[[1]], IRanges(bed[[2]] + 1, bed[[3]]))  # 0-based -> 1-based
} else {
  stop("Universe file must be .rds (GRanges) or .bed")
}

# Overlapping features (e.g. 5' UTRs of alternative transcripts) are merged so
# no stretch of sequence is weighted more than once.
universe = keepSeqlevels(universe, intersect(canonical_chrs, seqlevels(universe)),
                         pruning.mode = "coarse")
universe = reduce(universe, ignore.strand = TRUE)
strand(universe) = "*"

message(glue("Universe: {length(universe)} intervals, {sum(width(universe))} bp"))

#############################################################################
## Load LTRs (identical definition to permutation_analysis.R)
#############################################################################

message("Importing RepeatMasker annotation...")

LTR = read.csv(rmsk_file, header = TRUE)
LTR = LTR[!grepl("int", LTR$repName), ]

ltr_gr = makeGRangesFromDataFrame(LTR,
                                  keep.extra.columns = FALSE,
                                  start.field = "genoStart",
                                  end.field = "genoEnd",
                                  seqnames.field = "genoName")

end(ltr_gr[strand(ltr_gr) == "+"])   = end(ltr_gr[strand(ltr_gr) == "+"]) + 200
start(ltr_gr[strand(ltr_gr) == "-"]) = start(ltr_gr[strand(ltr_gr) == "-"]) - 200

ltr_gr = keepSeqlevels(ltr_gr, intersect(canonical_chrs, seqlevels(ltr_gr)),
                       pruning.mode = "coarse")
strand(ltr_gr) = "*"
ltr_index = GNCList(reduce(ltr_gr))  # overlapsAny only needs merged LTR territory

ltr_bp_in_universe = sum(width(GenomicRanges::intersect(universe, reduce(ltr_gr))))
ltr_fraction = ltr_bp_in_universe / sum(width(universe))
message(glue("LTRs (+200 bp) cover {round(100 * ltr_fraction, 2)}% of the universe"))

if (ltr_fraction > 0.95) {
  warning("LTRs cover >95% of the universe: almost every shuffled site will ",
          "overlap an LTR, so this test has little power. Is the universe ",
          "itself defined by LTRs?")
}

#############################################################################
## Load target sites
#############################################################################

message("Importing miRanda output...")

# Read only the columns needed; the full CSV includes bulky alignment strings.
needed = c("seqnames", "genomic_start", "genomic_end", "strand", "alignment_score")
header = names(read.csv(data_file, nrows = 1))
col_classes = ifelse(header %in% needed, NA, "NULL")
data = read.csv(data_file, header = TRUE, colClasses = col_classes)

data = data[data$alignment_score >= cutoff, ]
message(glue("{nrow(data)} target sites with alignment_score >= {cutoff}"))

if (nrow(data) == 0) stop("No target sites at this cutoff.")

sites = makeGRangesFromDataFrame(data,
                                 keep.extra.columns = FALSE,
                                 start.field = "genomic_start",
                                 end.field = "genomic_end")
sites = keepSeqlevels(sites, intersect(canonical_chrs, seqlevels(sites)),
                      pruning.mode = "coarse")
strand(sites) = "*"

n_sites_total = nrow(data)

# Only sites wholly inside the universe can be compared with sites shuffled
# within it.
sites_in_universe = subsetByOverlaps(sites, universe, type = "within")
n_sites_in_universe = length(sites_in_universe)

message(glue("{n_sites_in_universe} of {n_sites_total} sites lie entirely within the universe"))

if (n_sites_in_universe == 0) stop("No target sites lie within the universe.")

#############################################################################
## Shuffle function
#############################################################################

u_chr   = as.character(seqnames(universe))
u_start = start(universe)
u_width = width(universe)

site_widths = width(sites_in_universe)
width_groups = split(seq_along(site_widths), site_widths)

# For each site width, a universe interval is chosen with probability
# proportional to the number of positions at which a site of that width fits
# entirely inside it, then a start is drawn uniformly among those positions.
# Together this gives a uniform position across the whole universe.
slots_by_width = lapply(names(width_groups), function(w) {
  pmax(as.numeric(u_width) - as.numeric(w) + 1, 0)
})
names(slots_by_width) = names(width_groups)

for (w in names(slots_by_width)) {
  if (sum(slots_by_width[[w]]) == 0) {
    stop(glue("No universe interval is long enough for a {w} bp site."))
  }
}

shuffle_sites = function() {
  n = length(site_widths)
  pick = integer(n)
  new_start = numeric(n)
  for (w in names(width_groups)) {
    idx = width_groups[[w]]
    slots = slots_by_width[[w]]
    p = sample.int(length(slots), length(idx), replace = TRUE, prob = slots)
    pick[idx] = p
    new_start[idx] = u_start[p] + floor(runif(length(idx)) * slots[p])
  }
  GRanges(u_chr[pick], IRanges(start = new_start, width = site_widths))
}

count_ltr_overlaps = function(gr) sum(overlapsAny(gr, ltr_index))

#############################################################################
## Permutation test
#############################################################################

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Sanity check: shuffled sites must stay inside the universe
check = shuffle_sites()
stopifnot(all(overlapsAny(check, universe, type = "within")))

observed = count_ltr_overlaps(sites_in_universe)

message("Running permutations...")

permuted = unlist(mclapply(seq_len(ntimes),
                           function(i) count_ltr_overlaps(shuffle_sites()),
                           mc.cores = n_cores,
                           mc.set.seed = TRUE))

mean_random = mean(permuted)
sd_random   = sd(permuted)

results = data.frame(
  cutoff                 = cutoff,
  universe               = universe_label,
  universe_bp            = sum(width(universe)),
  ltr_fraction_universe  = ltr_fraction,
  n_sites_total          = n_sites_total,
  n_sites_in_universe    = n_sites_in_universe,
  observed_overlaps      = observed,
  mean_random_overlaps   = mean_random,
  sd_random_overlaps     = sd_random,
  fold_enrichment        = observed / mean_random,
  z_score                = (observed - mean_random) / sd_random,
  p_value_greater        = (sum(permuted >= observed) + 1) / (ntimes + 1),
  p_value_less           = (sum(permuted <= observed) + 1) / (ntimes + 1),
  ntimes                 = ntimes
)

print(results)

dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

summary_file = glue("{output_directory}/LTR_enrichment_{universe_label}_s{cutoff}.csv")
permuted_file = glue("{output_directory}/LTR_enrichment_{universe_label}_s{cutoff}_permuted.csv")

write.csv(results, summary_file, row.names = FALSE)
write.csv(data.frame(permutation = seq_len(ntimes), overlaps = permuted),
          permuted_file, row.names = FALSE)

message(glue("Wrote {summary_file}"))
message(glue("Wrote {permuted_file}"))
