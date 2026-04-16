################################################################################
# permutation_zscore.R
#
# Computes Z-score and p-values for the genome-shuffling permutation test.
#
# Reads:
#   <outdir>/observed_hits.txt    - single integer: real hit count
#   <outdir>/shuffled_hits.tsv    - two columns: iteration, hits
#
# Writes:
#   <outdir>/permutation_results.csv
#   <outdir>/permutation_histogram.pdf
#
# Usage:
#   Rscript permutation_zscore.R <outdir>
################################################################################

library(tidyverse)

args <- commandArgs(TRUE)
outdir <- args[1]

# ── Load data ────────────────────────────────────────────────────────────────

observed <- as.integer(readLines(file.path(outdir, "observed_hits.txt")))

shuffled <- read_tsv(
  file.path(outdir, "shuffled_hits.tsv"),
  col_types = cols(iteration = col_integer(), hits = col_integer())
)

cat("Observed hits:", observed, "\n")
cat("Shuffled iterations:", nrow(shuffled), "\n")
cat("Shuffled mean:", mean(shuffled$hits), "\n")
cat("Shuffled SD:", sd(shuffled$hits), "\n")

# ── Compute statistics ───────────────────────────────────────────────────────

n_iter <- nrow(shuffled)
mean_shuffled <- mean(shuffled$hits)
sd_shuffled <- sd(shuffled$hits)

# Z-score
z_score <- (observed - mean_shuffled) / sd_shuffled

# Empirical p-value (one-sided, greater)
# +1 in numerator and denominator for conservative estimate (Phipson & Smyth, 2010)
n_geq <- sum(shuffled$hits >= observed)
empirical_p <- (n_geq + 1) / (n_iter + 1)

# Analytical p-value from Z-score (normal approximation)
analytical_p <- pnorm(z_score, lower.tail = FALSE)

cat("\n── Results ──\n")
cat("Z-score:       ", round(z_score, 2), "\n")
cat("Empirical p:   ", format(empirical_p, scientific = TRUE, digits = 3), "\n")
cat("Analytical p:  ", format(analytical_p, scientific = TRUE, digits = 3), "\n")

# ── Write results ────────────────────────────────────────────────────────────

results <- tibble(
  observed_hits = observed,
  n_iterations = n_iter,
  mean_shuffled_hits = round(mean_shuffled, 1),
  sd_shuffled_hits = round(sd_shuffled, 1),
  z_score = round(z_score, 4),
  empirical_p_value = empirical_p,
  analytical_p_value = analytical_p,
  n_shuffled_geq_observed = n_geq
)

write_csv(results, file.path(outdir, "permutation_results.csv"))

# ── Histogram ────────────────────────────────────────────────────────────────

p <- ggplot(shuffled, aes(x = hits)) +
  geom_histogram(bins = 30, fill = "steelblue", colour = "white", alpha = 0.8) +
  geom_vline(xintercept = observed, colour = "red", linewidth = 1, linetype = "dashed") +
  annotate(
    "text",
    x = observed,
    y = Inf,
    label = paste0("Observed = ", format(observed, big.mark = ",")),
    colour = "red",
    hjust = -0.1,
    vjust = 2,
    size = 3.5
  ) +
  labs(
    title = "Genome-shuffling permutation test",
    subtitle = paste0(
      "Z = ", round(z_score, 2),
      ",  empirical p = ", format(empirical_p, scientific = TRUE, digits = 2),
      ",  n = ", n_iter, " iterations"
    ),
    x = "Total miRanda hits (shuffled windows)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(outdir, "permutation_histogram.pdf"),
  p, width = 7, height = 5
)

cat("\nResults written to:", file.path(outdir, "permutation_results.csv"), "\n")
cat("Histogram written to:", file.path(outdir, "permutation_histogram.pdf"), "\n")
