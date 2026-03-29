library(tidyverse)
library(glue)
library(grid)

args <- commandArgs(TRUE)

if (length(args) < 3) {
  stop(
    paste(
      "Usage: Rscript plot_gene_background_miRNA.R",
      "<summary_csv> <gene_metrics_csv> <output_prefix> [focal_tRF_cutoff] [focal_miRNA_cutoff]"
    )
  )
}

summary_csv <- args[1]
gene_metrics_csv <- args[2]
output_prefix <- args[3]
focal_tRF_cutoff <- if (length(args) >= 4) as.numeric(args[4]) else NA_real_
focal_miRNA_cutoff <- if (length(args) >= 5) as.numeric(args[5]) else NA_real_

summary_df <- readr::read_csv(summary_csv, show_col_types = FALSE)
metrics_df <- readr::read_csv(gene_metrics_csv, show_col_types = FALSE)

if (nrow(summary_df) == 0) {
  stop("Summary table is empty.")
}

target_gene <- summary_df$target_gene[[1]]

best_row <- summary_df %>%
  dplyr::filter(target_tRF_site_count > 0, target_overlap_site_count > 0)

if (!is.na(focal_tRF_cutoff)) {
  best_row <- best_row %>%
    dplyr::filter(tRF_cutoff == focal_tRF_cutoff)
}

if (!is.na(focal_miRNA_cutoff)) {
  best_row <- best_row %>%
    dplyr::filter(miRNA_cutoff == focal_miRNA_cutoff)
}

best_row <- best_row %>%
  dplyr::arrange(empirical_p_greater_overlap, empirical_p_greater_density)

if (nrow(best_row) == 0) {
  best_row <- summary_df %>%
    dplyr::filter(target_tRF_site_count > 0, target_overlap_site_count > 0) %>%
    dplyr::arrange(empirical_p_greater_overlap, empirical_p_greater_density)
}

if (nrow(best_row) == 0) {
  best_row <- summary_df %>% dplyr::slice(1)
} else {
  best_row <- best_row %>% dplyr::slice(1)
}

best_tRF_cutoff <- best_row$tRF_cutoff[[1]]
best_miRNA_cutoff <- best_row$miRNA_cutoff[[1]]

plot_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.4),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 9, hjust = 0)
  )

summary_long <- summary_df %>%
  dplyr::select(
    tRF_cutoff,
    miRNA_cutoff,
    target_overlap_site_count,
    target_overlap_sites_per_kb,
    target_overlap_percentile,
    target_density_percentile
  ) %>%
  tidyr::pivot_longer(
    cols = -c(tRF_cutoff, miRNA_cutoff),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      target_overlap_site_count = "Peg3 overlap sites",
      target_overlap_sites_per_kb = "Peg3 overlap density (sites/kb)",
      target_overlap_percentile = "Peg3 percentile by overlap count",
      target_density_percentile = "Peg3 percentile by overlap density"
    ),
    metric = factor(
      metric,
      levels = c(
        "Peg3 overlap sites",
        "Peg3 overlap density (sites/kb)",
        "Peg3 percentile by overlap count",
        "Peg3 percentile by overlap density"
      )
    ),
    miRNA_cutoff = factor(miRNA_cutoff)
  )

summary_plot <- ggplot(
  summary_long,
  aes(x = tRF_cutoff, y = value, color = miRNA_cutoff, group = miRNA_cutoff)
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("150" = "black", "200" = "firebrick")) +
  scale_x_continuous(breaks = sort(unique(summary_df$tRF_cutoff))) +
  labs(
    title = glue("{target_gene} 5' UTR Overlap Summary"),
    subtitle = "Points show Peg3 values across score cutoffs; colors indicate miRNA cutoff",
    x = "tRF alignment score cutoff",
    y = NULL,
    color = "miRNA cutoff"
  ) +
  plot_theme

background_df <- metrics_df %>%
  dplyr::filter(
    tRF_cutoff == best_tRF_cutoff,
    miRNA_cutoff == best_miRNA_cutoff,
    gene_name != target_gene,
    tRF_site_count > 0,
    miRNA_site_count > 0
  )

target_metrics <- metrics_df %>%
  dplyr::filter(
    tRF_cutoff == best_tRF_cutoff,
    miRNA_cutoff == best_miRNA_cutoff,
    gene_name == target_gene
  )

dist_long <- background_df %>%
  dplyr::select(gene_name, overlap_site_count, overlap_sites_per_kb) %>%
  tidyr::pivot_longer(
    cols = c(overlap_site_count, overlap_sites_per_kb),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      overlap_site_count = "Background ECDF: overlap sites",
      overlap_sites_per_kb = "Background ECDF: overlap density (sites/kb)"
    )
  )

target_long <- target_metrics %>%
  dplyr::select(gene_name, overlap_site_count, overlap_sites_per_kb) %>%
  tidyr::pivot_longer(
    cols = c(overlap_site_count, overlap_sites_per_kb),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      overlap_site_count = "Background ECDF: overlap sites",
      overlap_sites_per_kb = "Background ECDF: overlap density (sites/kb)"
    )
  )

distribution_plot <- ggplot(dist_long, aes(x = value)) +
  stat_ecdf(geom = "step", color = "black", linewidth = 0.6) +
  geom_vline(
    data = target_long,
    aes(xintercept = value),
    color = "firebrick",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  facet_wrap(~metric, scales = "free_x", ncol = 2) +
  scale_x_continuous(trans = "log1p") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = glue("Background Comparison At tRF >= {best_tRF_cutoff}, miRNA >= {best_miRNA_cutoff}"),
    subtitle = "Black curve: genes with at least one tRF site and one miRNA site; red dashed line: Peg3",
    x = "Value (log1p scale)",
    y = "Cumulative fraction of background genes",
    caption = glue(
      "Representative cutoff chosen by minimum empirical p-value for greater raw overlap.",
      " Background gene count at this cutoff: {nrow(background_df)}"
    )
  ) +
  plot_theme

png_file <- glue("{output_prefix}.png")
pdf_file <- glue("{output_prefix}.pdf")

png(filename = png_file, width = 1800, height = 1600, res = 200)
grid::grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(summary_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(distribution_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

pdf(pdf_file, width = 11, height = 10)
grid::grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(summary_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(distribution_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()
