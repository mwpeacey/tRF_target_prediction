################################################################################
# permutation_zscore.R
#
# Computes Z-scores and p-values at each alignment score threshold for the
# genome-shuffling permutation test.
#
# Reads:
#   <outdir>/observed_hits.tsv   - two columns: cutoff, hits
#   <outdir>/shuffled_hits.tsv   - columns: iteration, hits_70, hits_75, ...
#
# Writes:
#   <outdir>/permutation_results.csv
#   <outdir>/permutation_plot.pdf
#
# Usage:
#   Rscript permutation_zscore.R <outdir>
################################################################################

library(tidyverse)

args <- commandArgs(TRUE)
outdir <- args[1]

# ── Load data ────────────────────────────────────────────────────────────────

observed <- read_tsv(
  file.path(outdir, "observed_hits.tsv"),
  col_types = cols(cutoff = col_integer(), hits = col_integer())
)

shuffled <- read_tsv(
  file.path(outdir, "shuffled_hits.tsv"),
  col_types = cols(.default = col_integer())
)

cat("Observed thresholds:", nrow(observed), "\n")
cat("Shuffled iterations:", nrow(shuffled), "\n")

# ── Compute Z-score and p-values at each threshold ──────────────────────────

# Reshape shuffled from wide to long
shuffled_long <- shuffled %>%
  pivot_longer(
    cols = starts_with("hits_"),
    names_to = "cutoff",
    names_prefix = "hits_",
    values_to = "hits"
  ) %>%
  mutate(cutoff = as.integer(cutoff))

# Summarise null distribution per threshold
null_summary <- shuffled_long %>%
  group_by(cutoff) %>%
  summarise(
    mean_shuffled = mean(hits),
    sd_shuffled = sd(hits),
    n_iterations = n(),
    .groups = "drop"
  )

# Join with observed and compute statistics
results <- observed %>%
  rename(observed_hits = hits) %>%
  inner_join(null_summary, by = "cutoff") %>%
  mutate(
    z_score = if_else(sd_shuffled > 0,
      (observed_hits - mean_shuffled) / sd_shuffled,
      NA_real_
    ),
    # Empirical p-value (one-sided, conservative; Phipson & Smyth 2010)
    n_geq = map2_int(cutoff, observed_hits, function(co, obs) {
      vals <- shuffled_long %>% filter(cutoff == co) %>% pull(hits)
      sum(vals >= obs)
    }),
    empirical_p = (n_geq + 1) / (n_iterations + 1),
    analytical_p = pnorm(z_score, lower.tail = FALSE)
  )

cat("\n── Results ──\n")
print(as.data.frame(results), row.names = FALSE)

# ── Write results ────────────────────────────────────────────────────────────

write_csv(results, file.path(outdir, "permutation_results.csv"))

# ── Plot ─────────────────────────────────────────────────────────────────────

p_top <- ggplot(results, aes(x = cutoff, y = observed_hits)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(x = NULL, y = expression(log[10] ~ "(retained sites)")) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5, 10, 0, 10)
  )

p_bottom <- ggplot(results, aes(x = cutoff, y = z_score)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_line(colour = "firebrick", linewidth = 1) +
  labs(x = "Alignment score", y = "Z-score\n(shuffled genome)") +
  theme_classic(base_size = 14) +
  theme(plot.margin = margin(0, 10, 5, 10))

# Stack panels
if (requireNamespace("cowplot", quietly = TRUE)) {
  p_combined <- cowplot::plot_grid(
    p_top, p_bottom, ncol = 1, align = "v", rel_heights = c(1, 1)
  )
} else {
  # Fallback: use patchwork or save separately
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_combined <- p_top / p_bottom
  } else {
    message("Neither cowplot nor patchwork available. Saving panels separately.")
    ggsave(file.path(outdir, "permutation_plot_top.pdf"), p_top, width = 5, height = 3)
    ggsave(file.path(outdir, "permutation_plot_bottom.pdf"), p_bottom, width = 5, height = 3)
    cat("\nResults written to:", file.path(outdir, "permutation_results.csv"), "\n")
    quit(save = "no")
  }
}

ggsave(
  file.path(outdir, "permutation_plot.pdf"),
  p_combined, width = 5, height = 6
)

cat("\nResults written to:", file.path(outdir, "permutation_results.csv"), "\n")
cat("Plot written to:", file.path(outdir, "permutation_plot.pdf"), "\n")
