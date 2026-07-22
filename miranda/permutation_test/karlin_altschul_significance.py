#!/usr/bin/env python3
"""
karlin_altschul_significance.py

Assign a per-prediction statistical significance to miRanda hits, following the
Karlin-Altschul (PNAS 1990, 87:2264-2268) framework for local-alignment scores.

------------------------------------------------------------------------------
WHAT THIS DOES (and why it is set up this way)
------------------------------------------------------------------------------
Karlin-Altschul (K&A) show that, for a local-alignment scoring scheme applied to
random sequence, the number of distinct high-scoring segments with score >= S is
Poisson-distributed with mean

        E(S) = K * m * n * exp(-lambda * S)                        [K&A, p.2265]

i.e. log E(S) is LINEAR in S with slope -lambda. E(S) is the "E-value": the number
of hits scoring >= S expected BY CHANCE over the search space. The per-hit
significance (probability of at least one chance hit that good) is then

        p(S) = 1 - exp(-E(S))                                      (Poisson, m>=1)

K&A give closed-form lambda and K, but ONLY for ungapped, additive, i.i.d.
letter-pair scores. miRanda uses gapped, seed-weighted, wobble-aware scoring, so
those closed forms do NOT apply. K&A explicitly sanction the alternative:
significance "can be based on theoretical models OR on permutation reconstructions
of the observed data." So we estimate lambda and K*m*n EMPIRICALLY, per tRF, and
only borrow the functional form (log-linear tail / Poisson E-value).

Two calibration modes:

  * self-calibration (default): the genome-wide scan is captured at a permissive
    score threshold (here 70), so it is overwhelmingly chance hits. Per tRF, we fit
    log(#hits >= S) = A - lambda*S over the chance-dominated bulk. The intercept A
    already encodes K*m*n for the ACTUAL search space scanned (whole genome, both
    strands), so E(S) = exp(A - lambda*S) needs no further scaling. This mirrors
    K&A's own examples, where the null model is taken "directly from the [sequence]
    at hand." The handful of true targets in the extreme tail do not drive the fit
    because it is weighted toward the high-count bulk.

  * null-file calibration (--null): if you have a shuffled-genome result file (same
    CSV format, produced by the dinucleotide-preserving k=2 permutation), lambda and
    the rate are estimated from THAT instead. This is more rigorous because the null
    then contains no real targets. Rate is scaled by (real search space / null
    search space); pass --null-fraction for the genome fraction the null covered
    (e.g. 0.20) so the E-value refers to the full genome.

Calibration is per-tRF because lambda depends on the length and composition of BOTH
sequences (K&A Eq. 4), which differ from tRF to tRF. This is why a single genome-wide
score cutoff (e.g. 80) is not equivalent across tRFs.

------------------------------------------------------------------------------
OUTPUTS
------------------------------------------------------------------------------
Per-hit CSV (input columns preserved) with added columns:
  evalue        Expected number of chance hits scoring >= this hit's score,
                genome-wide, for this tRF. THIS IS THE K&A STATISTIC.
  pvalue        1 - exp(-evalue). Per-hit significance.
  z_ka          p-value expressed as a standard-normal deviate, qnorm(1 - p).
                A K&A-calibrated "Z" (correct EVD tail), for those who want a Z.
  z_gauss       (score - mean)/sd of this tRF's score distribution. The naive
                Gaussian Z (what miRanda's built-in -shuffle reports). Provided
                ONLY for comparison; it mis-models the extreme-value tail and will
                generally understate significance. Prefer evalue/pvalue/z_ka.
  fdr           Per-tRF, E-value-based FDR at this hit's score:
                E_t(s)/N_obs_t(s) = expected chance hits / observed hits at score
                >= s. This is the correct set-level control: because the E-value
                already spans the whole genome per tRF, a Benjamini-Hochberg pass
                over the p-values would double-count multiplicity and reject
                everything. Monotone non-increasing as score rises.
  significant   True if fdr < --fdr (default 0.05). These are the predictions to
                carry into subsequent steps.
  sig_evalue    True if evalue < --evalue-sig (default 1.0). BLAST-style per-hit
                E-value threshold, as an alternative/stricter flag.

Per-tRF summary CSV (--summary):
  lambda, log_Kmn, fit R^2 (EVD diagnostic), n_hits, score fit range, and the
  minimum score achieving E<1 and FDR<threshold (NA if never reached -- meaning
  that tRF has no sites distinguishable from chance).

------------------------------------------------------------------------------
USAGE
------------------------------------------------------------------------------
  # self-calibrated (quick cross-check; uses the scan as its own null):
  python3 karlin_altschul_significance.py miranda_output_70.csv out.csv \
      --summary out_summary.csv

  # PREFERRED: external shuffled-genome null (k=2 shuffle, 20% genome, 1000 iters).
  # The null contains no real targets by construction, so the "is score 70 low
  # enough" question does not arise -- every hit in it is chance:
  python3 karlin_altschul_significance.py miranda_output_70.csv out.csv \
      --null shuffled_score_histogram.csv --null-fraction 0.20 --null-iterations 1000 \
      --summary out_summary.csv

Requirements: pandas, numpy, scipy
"""

import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_csv", help="miRanda per-hit CSV (tRF,seqnames,alignment_score,...)")
    p.add_argument("output_csv", help="Output per-hit CSV with significance columns")
    p.add_argument("--summary", default=None, help="Optional per-tRF summary CSV path")
    p.add_argument("--null", default=None,
                   help="Optional shuffled-genome result CSV (same format) for calibration")
    p.add_argument("--null-fraction", type=float, default=1.0,
                   help="Genome fraction the --null file covers (e.g. 0.20). Default 1.0")
    p.add_argument("--null-iterations", type=int, default=1,
                   help="Number of shuffle iterations pooled into --null. The null "
                        "rate is divided by (fraction x iterations) to express the "
                        "E-value per single full-genome scan. Default 1")
    p.add_argument("--score-col", default="alignment_score", help="Score column name")
    p.add_argument("--trf-col", default="tRF", help="tRF/query id column name")
    p.add_argument("--fit-min-count", type=int, default=20,
                   help="Only fit score bins with at least this many hits (default 20)")
    p.add_argument("--fit-skip-low", type=int, default=1,
                   help="Skip this many lowest score bins (avoid capture-threshold edge). Default 1")
    p.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff for 'significant' flag")
    p.add_argument("--evalue-sig", type=float, default=1.0,
                   help="E-value used for the 'min score with E<x' summary column")
    return p.parse_args()


def fit_ev_tail(scores, min_count, skip_low, counts=None):
    """
    Fit log(N>=S) = A - lambda*S over the chance-dominated bulk of a score
    distribution.

    scores : integer scores. counts : optional per-score multiplicities (used when
    the input is a histogram, e.g. the pooled shuffled null). If counts is None
    each score counts once.

    Returns dict with lambda, logC (=A, ~log(K*m*n)), r2, fit score range, n,
    mean, sd. Weighted least squares (weights = survival counts, since
    var(log N) ~ 1/N). Returns None if too few usable points.
    """
    scores = np.asarray(scores, dtype=float)
    if counts is None:
        w_scores = np.ones_like(scores)
    else:
        w_scores = np.asarray(counts, dtype=float)
    n = w_scores.sum()
    if n < 3 * min_count:
        return None

    smin, smax = int(scores.min()), int(scores.max())
    # counts per exact score
    svals = np.arange(smin, smax + 1)
    hist = np.bincount((scores.astype(int) - smin), weights=w_scores,
                       minlength=len(svals)).astype(float)
    surv = np.cumsum(hist[::-1])[::-1]  # N(>=S)

    # weighted mean / sd of the score distribution (for the comparison Gaussian Z)
    mean = float(np.sum(svals * hist) / n)
    var = float(np.sum(hist * (svals - mean) ** 2) / n)
    sd = np.sqrt(var) if var > 0 else np.nan

    # choose fit window: skip the lowest `skip_low` bins (capture-threshold edge),
    # and keep bins whose survival count is reliable (>= min_count)
    start = skip_low
    mask = np.zeros_like(svals, dtype=bool)
    for i in range(start, len(svals)):
        if surv[i] >= min_count:
            mask[i] = True
    if mask.sum() < 3:
        return None

    x = svals[mask].astype(float)
    y = np.log(surv[mask])
    w = surv[mask]  # weights ~ counts

    # weighted linear fit y = A + b*x ; lambda = -b
    W = w.sum()
    xm = np.sum(w * x) / W
    ym = np.sum(w * y) / W
    sxx = np.sum(w * (x - xm) ** 2)
    sxy = np.sum(w * (x - xm) * (y - ym))
    b = sxy / sxx
    A = ym - b * xm
    lam = -b

    # weighted R^2
    yhat = A + b * x
    ss_res = np.sum(w * (y - yhat) ** 2)
    ss_tot = np.sum(w * (y - ym) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {
        "lambda": lam,
        "logC": A,
        "r2": r2,
        "fit_lo": float(x.min()),
        "fit_hi": float(x.max()),
        "n_hits": int(n),
        "mean": mean,
        "sd": sd,
    }


def evalue_from_fit(score, fit, log_scale=0.0):
    """E(S) = exp(logC - lambda*S) * scale.  log_scale = log(real/null search space)."""
    return np.exp(fit["logC"] - fit["lambda"] * score + log_scale)


def main():
    args = parse_args()

    print(f"[K&A] Reading {args.input_csv} ...", file=sys.stderr)
    df = pd.read_csv(args.input_csv)
    # strip the leading '>' that miRanda leaves on the tRF id, if present
    df[args.trf_col] = df[args.trf_col].astype(str).str.lstrip(">")
    sc = args.score_col
    tc = args.trf_col
    df[sc] = pd.to_numeric(df[sc], errors="coerce")
    df = df.dropna(subset=[sc])

    # ---- build calibration source (self or external null) ------------------
    # calib_groups: tRF -> (scores_array, counts_array_or_None)
    log_scale = 0.0
    if args.null:
        print(f"[K&A] Calibrating from null file {args.null} "
              f"(fraction={args.null_fraction}, iterations={args.null_iterations}) ...",
              file=sys.stderr)
        null_df = pd.read_csv(args.null)
        null_df[tc] = null_df[tc].astype(str).str.lstrip(">")
        null_df[sc] = pd.to_numeric(null_df[sc], errors="coerce")
        null_df = null_df.dropna(subset=[sc])
        has_counts = "count" in null_df.columns
        calib_groups = {}
        for k, v in null_df.groupby(tc):
            if has_counts:
                calib_groups[k] = (v[sc].values, v["count"].values)
            else:
                calib_groups[k] = (v[sc].values, None)
        # scale the pooled null count down to one full-genome scan:
        #   E_real(S) = E_null(S) / (iterations * fraction)
        denom = args.null_iterations * args.null_fraction
        if denom > 0:
            log_scale = np.log(1.0 / denom)
    else:
        print("[K&A] Self-calibrating from the scan (chance-dominated tail) ...",
              file=sys.stderr)
        calib_groups = {k: (v[sc].values, None) for k, v in df.groupby(tc)}

    # ---- fit per tRF -------------------------------------------------------
    fits = {}
    summary_rows = []
    for trf, (scores, counts) in calib_groups.items():
        fit = fit_ev_tail(scores, args.fit_min_count, args.fit_skip_low, counts=counts)
        fits[trf] = fit
        if fit is None:
            n_hits = int(counts.sum()) if counts is not None else len(scores)
            summary_rows.append({tc: trf, "n_hits": n_hits, "lambda": np.nan,
                                 "log_Kmn": np.nan, "fit_r2": np.nan,
                                 "fit_lo": np.nan, "fit_hi": np.nan,
                                 f"min_score_E<{args.evalue_sig}": np.nan})
            continue
        # minimum integer score achieving E < evalue-sig
        lo, hi = int(fit["fit_lo"]), int(fit["fit_hi"]) + 40
        min_sig = np.nan
        for s in range(lo, hi + 1):
            if evalue_from_fit(s, fit, log_scale) < args.evalue_sig:
                min_sig = s
                break
        summary_rows.append({
            tc: trf, "n_hits": fit["n_hits"], "lambda": round(fit["lambda"], 4),
            "log_Kmn": round(fit["logC"], 3), "fit_r2": round(fit["r2"], 4),
            "fit_lo": fit["fit_lo"], "fit_hi": fit["fit_hi"],
            f"min_score_E<{args.evalue_sig}": min_sig,
        })

    # ---- assign per-hit statistics (fully vectorised via dict maps) --------
    print("[K&A] Assigning per-hit E-values / p-values ...", file=sys.stderr)
    lam_map = {t: (f["lambda"] if f else np.nan) for t, f in fits.items()}
    logC_map = {t: (f["logC"] if f else np.nan) for t, f in fits.items()}
    mean_map = {t: (f["mean"] if f else np.nan) for t, f in fits.items()}
    sd_map = {t: (f["sd"] if f else np.nan) for t, f in fits.items()}

    lam_arr = df[tc].map(lam_map).to_numpy()
    logC_arr = df[tc].map(logC_map).to_numpy()
    mean_arr = df[tc].map(mean_map).to_numpy()
    sd_arr = df[tc].map(sd_map).to_numpy()
    scores_arr = df[sc].to_numpy()

    ev = np.exp(logC_arr - lam_arr * scores_arr + log_scale)
    with np.errstate(invalid="ignore", divide="ignore"):
        zg = np.where(sd_arr > 0, (scores_arr - mean_arr) / sd_arr, np.nan)

    df["evalue"] = ev
    pv = 1.0 - np.exp(-ev)
    # guard tiny E: p ~ E; also clip for qnorm
    pv = np.where(np.isnan(ev), np.nan, pv)
    df["pvalue"] = pv
    # K&A-calibrated Z = qnorm(1 - p); use -expm1 style via survival for tiny p
    with np.errstate(invalid="ignore"):
        df["z_ka"] = norm.isf(np.clip(pv, 1e-300, 1 - 1e-16))
    df["z_gauss"] = zg

    # ---- per-tRF, E-value-based FDR ----------------------------------------
    # The K&A E-value already spans the WHOLE-GENOME search space for each tRF,
    # so a Benjamini-Hochberg pass over the per-hit p-values would count that
    # genome-wide multiplicity a SECOND time and wrongly reject everything. The
    # correct set-level control is, per tRF and per score cutoff s,
    #     FDR(s) = expected chance hits(>=s) / observed hits(>=s)
    #            = E_t(s) / N_obs_t(s)
    # i.e. the fraction of the hits you would keep that are expected to be chance.
    print("[K&A] Per-tRF E-value-based FDR ...", file=sys.stderr)
    cnt = df.groupby([tc, sc]).size().rename("n").reset_index()
    parts = []
    for trf, sub in cnt.groupby(tc):
        fit = fits.get(trf)
        sub = sub.sort_values(sc).copy()
        svals = sub[sc].to_numpy()
        nobs = np.cumsum(sub["n"].to_numpy()[::-1])[::-1]  # observed hits >= s
        if fit is None:
            sub["fdr"] = np.nan
        else:
            eexp = np.exp(fit["logC"] - fit["lambda"] * svals + log_scale)
            f = np.clip(eexp / np.maximum(nobs, 1e-9), 0, 1)
            sub["fdr"] = np.minimum.accumulate(f)  # non-increasing as score rises
        parts.append(sub)
    cnt = pd.concat(parts, ignore_index=True)
    df = df.merge(cnt[[tc, sc, "fdr"]], on=[tc, sc], how="left")
    df["significant"] = df["fdr"] < args.fdr        # keep these downstream
    df["sig_evalue"] = df["evalue"] < args.evalue_sig

    # per-tRF minimum score reaching FDR<cutoff
    summ = pd.DataFrame(summary_rows)
    minfdr = (df[df["significant"]].groupby(tc)[sc].min()
              .rename(f"min_score_FDR<{args.fdr}"))
    summ = summ.merge(minfdr, on=tc, how="left")

    df.to_csv(args.output_csv, index=False)
    print(f"[K&A] Wrote per-hit results: {args.output_csv}", file=sys.stderr)
    if args.summary:
        summ.to_csv(args.summary, index=False)
        print(f"[K&A] Wrote per-tRF summary: {args.summary}", file=sys.stderr)

    # ---- console report ----------------------------------------------------
    n_sig = int(df["significant"].sum())
    n_trf = summ[tc].nunique()
    n_trf_fit = summ["lambda"].notna().sum()
    med_r2 = summ["fit_r2"].median()
    n_trf_sig = df[df["significant"]][tc].nunique()
    print("\n── Summary ──", file=sys.stderr)
    print(f"  tRFs:                       {n_trf} ({n_trf_fit} fitted)", file=sys.stderr)
    print(f"  Median EVD fit R^2:         {med_r2:.4f}", file=sys.stderr)
    print(f"  Total predictions:          {len(df):,}", file=sys.stderr)
    print(f"  Significant (FDR<{args.fdr}):     {n_sig:,}", file=sys.stderr)
    print(f"  tRFs with >=1 sig site:     {n_trf_sig} / {n_trf}", file=sys.stderr)


if __name__ == "__main__":
    main()
