################################################################################
# permutation_zscore.R
#
# Computes Z-score and p-values from the permutation test output and writes
# a CSV ready to plot as a histogram.
#
# Reads:
#   <outdir>/observed_hits.txt   - single integer (observed hit count)
#   <outdir>/shuffled_hits.tsv   - two columns: iteration, hits
#
# Writes:
#   <outdir>/permutation_histogram_data.csv
#     Columns: iteration, shuffled_hits, observed_hits, z_score,
#              empirical_p, analytical_p
#     One row per iteration. Histogram the shuffled_hits column;
#     observed_hits is the vertical line.
#
# Usage:
#   Rscript permutation_zscore.R <score_threshold> <n_iterations>
#
# Example:
#   Rscript permutation_zscore.R 80 1000
################################################################################

library(tidyverse)

# ── Parameters ──────────────────────────────────────────────────────────────

args <- commandArgs(TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript permutation_zscore.R <score_threshold> <n_iterations>")
}

score_threshold <- as.integer(args[1])
n_iter          <- as.integer(args[2])
outdir          <- 'import/miranda/miranda_permutation/'

cat("Score threshold:", score_threshold, "\n")
cat("Iterations requested:", n_iter, "\n")

# ── Load data ────────────────────────────────────────────────────────────────

obs_hits <- as.integer(readLines(file.path(outdir, "observed_hits.txt"))[1])

shuffled <- read_tsv(
  file.path(outdir, "shuffled_hits.tsv"),
  col_types = cols(iteration = col_integer(), hits = col_integer())
)

# Use first n_iter iterations
shuf <- shuffled %>%
  arrange(iteration) %>%
  slice_head(n = n_iter)

actual_n <- nrow(shuf)
if (actual_n < n_iter) {
  warning("Only ", actual_n, " iterations available (requested ", n_iter, ").")
}

shuf_vals <- shuf$hits

# ── Compute statistics ──────────────────────────────────────────────────────

mean_shuf <- mean(shuf_vals)
sd_shuf   <- sd(shuf_vals)
z_score   <- if (sd_shuf > 0) (obs_hits - mean_shuf) / sd_shuf else NA_real_
n_geq     <- sum(shuf_vals >= obs_hits)
emp_p     <- (n_geq + 1) / (actual_n + 1)
ana_p     <- pnorm(z_score, lower.tail = FALSE)

cat("\n── Summary ──\n")
cat("  Observed hits:  ", obs_hits, "\n")
cat("  Shuffled mean:  ", round(mean_shuf, 1), "\n")
cat("  Shuffled SD:    ", round(sd_shuf, 1), "\n")
cat("  Z-score:        ", round(z_score, 2), "\n")
cat("  Empirical p:    ", formatC(emp_p, format = "e", digits = 2), "\n")
cat("  Analytical p:   ", formatC(ana_p, format = "e", digits = 2), "\n")
cat("  N iterations:   ", actual_n, "\n")

# ── Write plot-ready data ────────────────────────────────────────────────────

out_df <- tibble(
  iteration     = shuf$iteration,
  shuffled_hits = shuf_vals,
  observed_hits = obs_hits,
  z_score       = z_score,
  empirical_p   = emp_p,
  analytical_p  = ana_p
)

outfile <- file.path(outdir, "permutation_histogram_data.csv")
write_csv(out_df, outfile)
cat("\nHistogram data written to:", outfile, "\n")
