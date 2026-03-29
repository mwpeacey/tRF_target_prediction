library(tidyverse)
library(glue)

args <- commandArgs(TRUE)

if (length(args) < 5) {
  stop(
    paste(
      "Usage: Rscript plot_gene_background_overlap_ecdf_miRNA.R",
      "<summary_csv> <gene_metrics_csv> <output_pdf> <tRF_cutoff> <miRNA_cutoff> [count|density]"
    )
  )
}

summary_csv <- args[1]
gene_metrics_csv <- args[2]
output_pdf <- args[3]
tRF_cutoff <- as.numeric(args[4])
miRNA_cutoff <- as.numeric(args[5])
metric_type <- if (length(args) >= 6) args[6] else "count"

if (!metric_type %in% c("count", "density")) {
  stop("metric_type must be either 'count' or 'density'.")
}

summary_df <- readr::read_csv(summary_csv, show_col_types = FALSE)
metrics_df <- readr::read_csv(gene_metrics_csv, show_col_types = FALSE)

target_gene <- summary_df$target_gene[[1]]

summary_row <- summary_df %>%
  dplyr::filter(
    tRF_cutoff == !!tRF_cutoff,
    miRNA_cutoff == !!miRNA_cutoff
  )

if (nrow(summary_row) == 0) {
  stop("Requested cutoff combination was not found in the summary table.")
}

background_df <- metrics_df %>%
  dplyr::filter(
    tRF_cutoff == !!tRF_cutoff,
    miRNA_cutoff == !!miRNA_cutoff,
    gene_name != target_gene,
    tRF_site_count > 0,
    miRNA_site_count > 0
  )

target_df <- metrics_df %>%
  dplyr::filter(
    tRF_cutoff == !!tRF_cutoff,
    miRNA_cutoff == !!miRNA_cutoff,
    gene_name == target_gene
  )

if (nrow(target_df) == 0) {
  stop("Target gene row was not found in the gene metrics table.")
}

if (metric_type == "count") {
  metric_col <- "overlap_site_count"
  target_value <- target_df$overlap_site_count[[1]]
  target_percentile <- summary_row$target_overlap_percentile[[1]]
  x_label <- "Overlapping tRF sites per gene"
  title_suffix <- "Overlap Count"
  subtitle_metric <- "ECDF of overlap counts"
  value_label <- glue("{target_value} overlapping sites")
} else {
  metric_col <- "overlap_sites_per_kb"
  target_value <- target_df$overlap_sites_per_kb[[1]]
  target_percentile <- summary_row$target_density_percentile[[1]]
  x_label <- "Overlapping tRF sites per kb of 5' UTR"
  title_suffix <- "Overlap Density"
  subtitle_metric <- "ECDF of overlap density"
  value_label <- glue("{round(target_value, 2)} sites per kb")
}

background_gene_count <- summary_row$background_gene_count[[1]]

annotation_text <- glue(
  "{target_gene}: {value_label}\n",
  "{round(target_percentile, 1)}th percentile\n",
  "n = {background_gene_count} background genes"
)

plot_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.4),
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 9, hjust = 0)
  )

p <- ggplot(background_df, aes(x = .data[[metric_col]])) +
  stat_ecdf(geom = "step", color = "black", linewidth = 0.7) +
  geom_vline(
    xintercept = target_value,
    color = "firebrick",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  geom_label(
    data = tibble(x = target_value, y = 0.15, label = annotation_text),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 0,
    size = 3.6,
    linewidth = 0.25,
    fill = "white"
  ) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.18))
  ) +
  labs(
    title = glue("{target_gene} 5' UTR {title_suffix} Relative To Background"),
    subtitle = glue("{subtitle_metric} across testable genes at tRF >= {tRF_cutoff}, miRNA >= {miRNA_cutoff}"),
    x = x_label,
    y = "Cumulative fraction of background genes",
    caption = "Black curve: background genes with at least one tRF site and one miRNA site. Red dashed line: Peg3."
  ) +
  plot_theme

ggsave(
  filename = output_pdf,
  plot = p,
  width = 8.5,
  height = 5.5,
  units = "in"
)
