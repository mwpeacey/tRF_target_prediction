################################################################################
# ka_significance.R
#
# Karlin-Altschul per-tRF E-value / FDR for miRanda hits, in R.
# Port of miranda/permutation_test/karlin_altschul_significance.py.
#
# The Python script is the validated reference. Cross-check ONCE before trusting
# this (see "VALIDATION" at the bottom): run both on the same inputs and confirm
# ka_fdr here matches `fdr` from the Python output to numerical tolerance.
#
# Usage inside annotate_targets.R, immediately AFTER the tRF-id parsing block
# (i.e. once `data$tRF` is the unique id "tRF_1"), and before the unique_data
# collapse so the columns propagate to both exports:
#
#   source('R_scripts/ka_significance.R')
#   null_hist <- readr::read_csv(
#     'import/miranda/genome_shuffle_null/shuffled_score_histogram.csv')
#   data <- add_ka_significance(data, null_hist = null_hist,
#                               null_fraction = 0.20, null_iterations = 1000)
#
# Omit null_hist to self-calibrate from the scan (quick cross-check only; the
# shuffle-derived null is the defensible version).
#
# Adds columns: ka_lambda, ka_evalue, ka_pvalue, ka_fdr, ka_significant.
################################################################################

library(dplyr)
library(stringr)

# Normalise a tRF id to the unique key: strip a leading '>' and everything up to
# and including '|'.  "tDR-55:76-Ala-AGC-1|tRF_1" -> "tRF_1"; "tRF_1" -> "tRF_1".
.ka_key <- function(x) str_remove(str_remove(as.character(x), "^>"), "^.*\\|")

# Fit log(N>=s) = A - lambda*s over the chance-dominated bulk of one tRF's score
# distribution. `scores` integer scores, `counts` their multiplicities (1 each in
# self-calibration; histogram counts for the shuffled null). Returns NULL if the
# distribution is too small to fit.
.fit_ev_tail <- function(scores, counts, min_count = 20, skip_low = 1) {
  n <- sum(counts)
  if (n < 3 * min_count) return(NULL)
  smin <- min(scores); smax <- max(scores)
  svals <- smin:smax
  hist <- as.numeric(tapply(counts, factor(scores, levels = svals), sum))
  hist[is.na(hist)] <- 0
  surv <- rev(cumsum(rev(hist)))                       # N(>=s)

  mean_s <- sum(svals * hist) / n
  var_s  <- sum(hist * (svals - mean_s)^2) / n
  sd_s   <- if (var_s > 0) sqrt(var_s) else NA_real_

  idx <- (skip_low + 1):length(svals)                  # skip capture-threshold edge
  idx <- idx[surv[idx] >= min_count]                   # reliable bins only
  if (length(idx) < 3) return(NULL)

  x <- svals[idx]; y <- log(surv[idx]); w <- surv[idx] # weights ~ counts
  fit <- lm(y ~ x, weights = w)
  lambda <- -coef(fit)[["x"]]; logC <- coef(fit)[["(Intercept)"]]
  yhat <- predict(fit); ybar <- sum(w * y) / sum(w)
  r2 <- 1 - sum(w * (y - yhat)^2) / sum(w * (y - ybar)^2)

  list(lambda = lambda, logC = logC, r2 = r2, smax = smax,
       mean = mean_s, sd = sd_s, n = n)
}

add_ka_significance <- function(data,
                                null_hist       = NULL,
                                null_fraction   = 1,
                                null_iterations = 1,
                                fdr_cutoff      = 0.05,
                                min_count = 20, skip_low = 1,
                                score_col = "alignment_score",
                                trf_col   = "tRF") {

  data$.ka_key <- .ka_key(data[[trf_col]])
  data[[score_col]] <- as.integer(round(as.numeric(data[[score_col]])))

  # ── calibration source: external shuffled null, or self ──────────────────
  if (!is.null(null_hist)) {
    nh <- null_hist
    nh$.ka_key   <- .ka_key(nh[[trf_col]])
    nh[[score_col]] <- as.integer(round(as.numeric(nh[[score_col]])))
    if (!"count" %in% names(nh)) nh$count <- 1L
    calib <- nh %>%
      group_by(.ka_key) %>%
      summarise(sc = list(.data[[score_col]]), ct = list(count), .groups = "drop")
    # scale pooled null down to one full-genome scan: E_real = E_null / (frac*iter)
    log_scale <- log(1 / (null_fraction * null_iterations))
  } else {
    calib <- data %>%
      group_by(.ka_key) %>%
      summarise(sc = list(.data[[score_col]]),
                ct = list(rep(1, dplyr::n())), .groups = "drop")
    log_scale <- 0
  }

  # ── fit per tRF ──────────────────────────────────────────────────────────
  fits <- setNames(
    lapply(seq_len(nrow(calib)),
           function(i) .fit_ev_tail(calib$sc[[i]], calib$ct[[i]], min_count, skip_low)),
    calib$.ka_key)

  fitdf <- tibble(
    .ka_key = names(fits),
    ka_lambda = vapply(fits, function(f) if (is.null(f)) NA_real_ else f$lambda, numeric(1)),
    ka_logC   = vapply(fits, function(f) if (is.null(f)) NA_real_ else f$logC,   numeric(1)))

  # ── per-hit E-value / p-value (vectorised) ───────────────────────────────
  data <- data %>% left_join(fitdf, by = ".ka_key")
  data$ka_evalue <- exp(data$ka_logC - data$ka_lambda * data[[score_col]] + log_scale)
  data$ka_pvalue <- 1 - exp(-data$ka_evalue)

  # ── per-tRF FDR = E_fit(s) / N_observed(>=s) ─────────────────────────────
  # Built from the REAL observed (key, score) pairs only, so scores beyond a
  # tRF's observed maximum are never assigned significance. NO Benjamini-Hochberg
  # on top: the E-value already spans the whole genome per tRF.
  fdr_tbl <- data %>%
    count(.ka_key, .data[[score_col]], name = "n") %>%
    arrange(.ka_key, desc(.data[[score_col]])) %>%
    group_by(.ka_key) %>% mutate(Nobs = cumsum(n)) %>% ungroup() %>%
    left_join(fitdf, by = ".ka_key") %>%
    mutate(E = exp(ka_logC - ka_lambda * .data[[score_col]] + log_scale),
           fdr_raw = pmin(1, E / pmax(Nobs, 1))) %>%
    group_by(.ka_key) %>% arrange(.data[[score_col]], .by_group = TRUE) %>%
    mutate(ka_fdr = cummin(fdr_raw)) %>% ungroup() %>%      # monotone as score rises
    select(.ka_key, all_of(score_col), ka_fdr)

  data <- data %>%
    left_join(fdr_tbl, by = c(".ka_key", score_col)) %>%
    mutate(ka_significant = !is.na(ka_fdr) & ka_fdr < fdr_cutoff) %>%
    select(-.ka_key, -ka_logC)

  message(sprintf("[K&A] %d tRFs fitted; %d/%d hits significant (FDR<%.2g); %d tRFs with >=1 sig site",
                  sum(!vapply(fits, is.null, logical(1))),
                  sum(data$ka_significant, na.rm = TRUE), nrow(data), fdr_cutoff,
                  dplyr::n_distinct(data$.ka_key[data$ka_significant])))
  data
}

################################################################################
# VALIDATION (run once):
#   Rscript -e '
#     source("R_scripts/ka_significance.R")
#     d <- read.csv("import/miranda/new/miranda_output_70.csv")
#     d$tRF <- sub("^.*\\|","", sub("^>","", d$tRF))
#     d <- add_ka_significance(d)                       # self-calibrated
#     write.csv(d[,c("tRF","alignment_score","ka_evalue","ka_fdr","ka_significant")],
#               "ka_R_selfcal.csv", row.names = FALSE)'
#   # then compare to Python self-calibrated output on the same file:
#   #   python3 miranda/permutation_test/karlin_altschul_significance.py \
#   #       import/miranda/new/miranda_output_70.csv ka_py_selfcal.csv
#   # join on (tRF, alignment_score) and confirm ka_fdr == fdr within tolerance.
################################################################################
