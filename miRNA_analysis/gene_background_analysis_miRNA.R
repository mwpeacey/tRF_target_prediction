library(tidyverse)
library(glue)
library(GenomicRanges)
library(GenomeInfoDb)
library(parallel)

args <- commandArgs(TRUE)

if (length(args) < 7) {
  stop(
    paste(
      "Usage: Rscript gene_background_analysis_miRNA.R",
      "<tRF_file> <miRNA_file> <five_UTR_rds> <rmsk_file>",
      "<transcript_gtf> <target_gene> <output_directory> [utr_ltr|utr_only]"
    )
  )
}

tRF_file         <- args[1]
miRNA_file       <- args[2]
five_UTR_file    <- args[3]
rmsk_file        <- args[4]
transcript_gtf   <- args[5]
target_gene      <- args[6]
output_directory <- args[7]
universe_mode    <- if (length(args) >= 8) args[8] else "utr_ltr"

if (!universe_mode %in% c("utr_ltr", "utr_only")) {
  stop("universe_mode must be either 'utr_ltr' or 'utr_only'.")
}

tRF_cutoffs   <- c(70, 75, 80, 85, 90)
miRNA_cutoffs <- c(150, 200)
canonical_chrs <- paste0("chr", c(1:19, "X", "Y"))
n_cores <- as.integer(Sys.getenv("NSLOTS"))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1L

strip_version <- function(x) {
  sub("\\..*$", "", x)
}

extract_gtf_attribute <- function(x, key) {
  stringr::str_match(x, paste0(key, ' "([^"]+)"'))[, 2]
}

make_hit_granges <- function(df) {
  makeGRangesFromDataFrame(
    df,
    keep.extra.columns = TRUE
  )
}

count_overlap_sites <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) {
    return(0L)
  }

  as.integer(sum(countOverlaps(a, b, ignore.strand = FALSE) > 0))
}

message("Loading 5' UTR annotation...")

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

five_utr_gr <- readRDS(five_UTR_file)
five_utr_gr <- keepSeqlevels(
  five_utr_gr,
  canonical_chrs,
  pruning.mode = "coarse"
)
names(five_utr_gr) <- strip_version(names(five_utr_gr))

message("Loading transcript-to-gene map...")

gtf_tbl <- readr::read_tsv(
  transcript_gtf,
  comment = "#",
  col_names = FALSE,
  col_types = cols(.default = col_character()),
  progress = FALSE
)

tx_gene_map <- gtf_tbl %>%
  dplyr::filter(X3 == "transcript") %>%
  dplyr::transmute(
    transcript_id = strip_version(extract_gtf_attribute(X9, "transcript_id")),
    gene_name = dplyr::coalesce(
      extract_gtf_attribute(X9, "gene_name"),
      extract_gtf_attribute(X9, "ref_gene_name")
    )
  ) %>%
  dplyr::filter(!is.na(transcript_id), !is.na(gene_name)) %>%
  dplyr::distinct()

mcols(five_utr_gr)$transcript_id <- names(five_utr_gr)
mcols(five_utr_gr)$gene_name <- tx_gene_map$gene_name[
  match(mcols(five_utr_gr)$transcript_id, tx_gene_map$transcript_id)
]

five_utr_gr <- five_utr_gr[!is.na(mcols(five_utr_gr)$gene_name)]

if (universe_mode == "utr_ltr") {
  message("Loading LTR annotation...")

  ltr_df <- read.csv(rmsk_file, header = TRUE)
  ltr_df <- ltr_df[!grepl("int", ltr_df$repName), ]

  ltr_gr <- makeGRangesFromDataFrame(
    ltr_df,
    keep.extra.columns = TRUE,
    start.field = "genoStart",
    end.field = "genoEnd",
    seqnames.field = "genoName"
  )

  end(ltr_gr[strand(ltr_gr) == "+"]) <- end(ltr_gr[strand(ltr_gr) == "+"]) + 200
  start(ltr_gr[strand(ltr_gr) == "-"]) <- start(ltr_gr[strand(ltr_gr) == "-"]) - 200

  ltr_gr <- keepSeqlevels(ltr_gr, canonical_chrs, pruning.mode = "coarse")

  message("Building gene-level 5' UTR/LTR universe...")

  analysis_universe_gr <- GenomicRanges::intersect(
    five_utr_gr,
    ltr_gr,
    ignore.strand = TRUE
  )
} else {
  message("Building gene-level 5' UTR universe...")
  analysis_universe_gr <- five_utr_gr
}

if (length(analysis_universe_gr) == 0) {
  stop(glue("No ranges were available after applying universe_mode = '{universe_mode}'."))
}

utr_back_ov <- findOverlaps(analysis_universe_gr, five_utr_gr, ignore.strand = FALSE)

utr_ltr_gene_map <- tibble(
  queryHits = queryHits(utr_back_ov),
  gene_name = mcols(five_utr_gr)$gene_name[subjectHits(utr_back_ov)]
) %>%
  dplyr::distinct() %>%
  dplyr::group_by(queryHits) %>%
  dplyr::summarize(
    gene_name = dplyr::first(gene_name),
    n_gene_names = dplyr::n_distinct(gene_name),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_gene_names == 1) %>%
  dplyr::select(-n_gene_names)

analysis_universe_gr <- analysis_universe_gr[utr_ltr_gene_map$queryHits]
mcols(analysis_universe_gr)$gene_name <- utr_ltr_gene_map$gene_name

message(glue("Reducing gene-level universes with up to {n_cores} core(s)..."))

gene_universe_split <- split(analysis_universe_gr, mcols(analysis_universe_gr)$gene_name)
gene_universe_list <- as.list(gene_universe_split)

if (.Platform$OS.type == "unix" && n_cores > 1) {
  gene_universe_list <- parallel::mclapply(
    gene_universe_list,
    GenomicRanges::reduce,
    mc.cores = n_cores
  )
} else {
  gene_universe_list <- lapply(gene_universe_list, GenomicRanges::reduce)
}

gene_universe_grl <- GenomicRanges::GRangesList(gene_universe_list)

gene_lengths <- tibble(
  gene_name = names(gene_universe_grl),
  effective_length_bp = vapply(
    gene_universe_grl,
    function(gr) as.numeric(sum(width(gr))),
    numeric(1)
  )
) %>%
  dplyr::filter(effective_length_bp > 0)

if (!(target_gene %in% gene_lengths$gene_name)) {
  stop(
    glue(
      "Target gene '{target_gene}' has no testable universe ",
      "for universe_mode = '{universe_mode}'."
    )
  )
}

gene_universe_flat <- unlist(gene_universe_grl, use.names = FALSE)
mcols(gene_universe_flat)$gene_name <- rep(
  names(gene_universe_grl),
  elementNROWS(gene_universe_grl)
)

message("Loading tRF and miRNA predictions...")

tRF_data <- read.csv(tRF_file, header = TRUE)
miRNA_data <- read.csv(miRNA_file, header = TRUE)

process_cutoff_combo <- function(tRF_cutoff, miRNA_cutoff) {
  message(glue("Processing tRF cutoff {tRF_cutoff}, miRNA cutoff {miRNA_cutoff}"))

  tRF_subset <- tRF_data %>%
    dplyr::filter(alignment_score >= tRF_cutoff)

  miRNA_subset <- miRNA_data %>%
    dplyr::filter(alignment_score >= miRNA_cutoff)

  if (nrow(tRF_subset) == 0 || nrow(miRNA_subset) == 0) {
    return(list(summary = tibble(), by_gene = tibble()))
  }

  tRF_gr <- make_hit_granges(tRF_subset)
  miRNA_gr <- make_hit_granges(miRNA_subset)

  tRF_ov <- findOverlaps(tRF_gr, gene_universe_flat, ignore.strand = FALSE)
  miRNA_ov <- findOverlaps(miRNA_gr, gene_universe_flat, ignore.strand = FALSE)

  tRF_gene_hits <- tibble(
    tRF_hit = queryHits(tRF_ov),
    gene_name = mcols(gene_universe_flat)$gene_name[subjectHits(tRF_ov)]
  ) %>%
    dplyr::distinct()

  miRNA_gene_hits <- tibble(
    miRNA_hit = queryHits(miRNA_ov),
    gene_name = mcols(gene_universe_flat)$gene_name[subjectHits(miRNA_ov)]
  ) %>%
    dplyr::distinct()

  genes_with_signal <- union(tRF_gene_hits$gene_name, miRNA_gene_hits$gene_name)

  tRF_gene_counts <- tRF_gene_hits %>%
    dplyr::count(gene_name, name = "tRF_site_count")

  miRNA_gene_counts <- miRNA_gene_hits %>%
    dplyr::count(gene_name, name = "miRNA_site_count")

  hit_overlaps <- findOverlaps(tRF_gr, miRNA_gr, ignore.strand = FALSE)

  overlap_gene_counts <- tibble(
    tRF_hit = queryHits(hit_overlaps),
    miRNA_hit = subjectHits(hit_overlaps)
  ) %>%
    dplyr::inner_join(tRF_gene_hits, by = "tRF_hit") %>%
    dplyr::inner_join(miRNA_gene_hits, by = c("miRNA_hit", "gene_name")) %>%
    dplyr::distinct(gene_name, tRF_hit) %>%
    dplyr::count(gene_name, name = "overlap_site_count")

  gene_metrics <- gene_lengths %>%
    dplyr::filter(gene_name %in% genes_with_signal | gene_name == target_gene) %>%
    dplyr::left_join(tRF_gene_counts, by = "gene_name") %>%
    dplyr::left_join(miRNA_gene_counts, by = "gene_name") %>%
    dplyr::left_join(overlap_gene_counts, by = "gene_name") %>%
    dplyr::mutate(
      tRF_site_count = dplyr::coalesce(tRF_site_count, 0L),
      miRNA_site_count = dplyr::coalesce(miRNA_site_count, 0L),
      overlap_site_count = dplyr::coalesce(overlap_site_count, 0L),
      overlap_sites_per_kb = dplyr::if_else(
        effective_length_bp > 0,
        overlap_site_count / (effective_length_bp / 1000),
        NA_real_
      ),
      overlap_fraction_of_tRF_sites = dplyr::if_else(
        tRF_site_count > 0,
        overlap_site_count / tRF_site_count,
        NA_real_
      ),
      universe_mode = universe_mode,
      target_gene = target_gene,
      tRF_cutoff = tRF_cutoff,
      miRNA_cutoff = miRNA_cutoff
    )

  target_metrics <- gene_metrics %>%
    dplyr::filter(gene_name == target_gene)

  if (nrow(target_metrics) == 0) {
    return(list(summary = tibble(), by_gene = gene_metrics))
  }

  background_metrics <- gene_metrics %>%
    dplyr::filter(
      gene_name != target_gene,
      tRF_site_count > 0,
      miRNA_site_count > 0
    )

  target_overlap <- target_metrics$overlap_site_count[[1]]
  target_density <- target_metrics$overlap_sites_per_kb[[1]]

  if (nrow(background_metrics) == 0) {
    overlap_percentile <- NA_real_
    density_percentile <- NA_real_
    empirical_p_greater_overlap <- NA_real_
    empirical_p_less_overlap <- NA_real_
    empirical_p_greater_density <- NA_real_
    empirical_p_less_density <- NA_real_
  } else {
    overlap_percentile <- mean(
      background_metrics$overlap_site_count <= target_overlap
    ) * 100

    density_percentile <- mean(
      background_metrics$overlap_sites_per_kb <= target_density
    ) * 100

    empirical_p_greater_overlap <- (
      sum(background_metrics$overlap_site_count >= target_overlap) + 1
    ) / (nrow(background_metrics) + 1)

    empirical_p_less_overlap <- (
      sum(background_metrics$overlap_site_count <= target_overlap) + 1
    ) / (nrow(background_metrics) + 1)

    empirical_p_greater_density <- (
      sum(background_metrics$overlap_sites_per_kb >= target_density) + 1
    ) / (nrow(background_metrics) + 1)

    empirical_p_less_density <- (
      sum(background_metrics$overlap_sites_per_kb <= target_density) + 1
    ) / (nrow(background_metrics) + 1)
  }

  summary_tbl <- tibble(
    target_gene = target_gene,
    universe_mode = universe_mode,
    tRF_cutoff = tRF_cutoff,
    miRNA_cutoff = miRNA_cutoff,
    target_effective_length_bp = target_metrics$effective_length_bp[[1]],
    target_tRF_site_count = target_metrics$tRF_site_count[[1]],
    target_miRNA_site_count = target_metrics$miRNA_site_count[[1]],
    target_overlap_site_count = target_overlap,
    target_overlap_sites_per_kb = target_density,
    target_overlap_fraction_of_tRF_sites = target_metrics$overlap_fraction_of_tRF_sites[[1]],
    background_gene_count = nrow(background_metrics),
    target_overlap_percentile = overlap_percentile,
    target_density_percentile = density_percentile,
    empirical_p_greater_overlap = empirical_p_greater_overlap,
    empirical_p_less_overlap = empirical_p_less_overlap,
    empirical_p_greater_density = empirical_p_greater_density,
    empirical_p_less_density = empirical_p_less_density
  )

  list(summary = summary_tbl, by_gene = gene_metrics)
}

cutoff_grid <- tidyr::crossing(
  tRF_cutoff = tRF_cutoffs,
  miRNA_cutoff = miRNA_cutoffs
)

worker_fun <- function(i) {
  process_cutoff_combo(
    tRF_cutoff = cutoff_grid$tRF_cutoff[[i]],
    miRNA_cutoff = cutoff_grid$miRNA_cutoff[[i]]
  )
}

message(glue("Evaluating {nrow(cutoff_grid)} cutoff combinations with up to {n_cores} core(s)..."))

combo_results <- if (.Platform$OS.type == "unix" && n_cores > 1) {
  parallel::mclapply(seq_len(nrow(cutoff_grid)), worker_fun, mc.cores = n_cores)
} else {
  lapply(seq_len(nrow(cutoff_grid)), worker_fun)
}

results_summary <- dplyr::bind_rows(purrr::map(combo_results, "summary"))
results_by_gene <- dplyr::bind_rows(purrr::map(combo_results, "by_gene"))

write.csv(
  results_summary,
  glue("{output_directory}/{target_gene}_vs_5UTR_background_summary_{universe_mode}.csv"),
  row.names = FALSE
)

write.csv(
  results_by_gene,
  glue("{output_directory}/{target_gene}_vs_5UTR_background_gene_metrics_{universe_mode}.csv"),
  row.names = FALSE
)
