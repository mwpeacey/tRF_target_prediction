#!/usr/bin/env python3
"""
sample_windows_from_fasta.py

Sample a random subset of the pre-generated genomic window targets for the
permutation null, drawing with equal probability from BOTH the plus- and
minus-strand window FASTAs used for the real miRanda scan.

The real scan treats each strand as an independent target: it runs miRanda
against windows.fa (plus) and windows_minus.fa (minus) separately. This mirrors
that exactly — it pools the records of both files and samples uniformly across
the combined set, so a sampled target is equally likely to be a plus- or
minus-strand window. No reverse-complementation is performed anywhere: the two
strands are just separate target sequences, as in the real prediction. The
sampled targets are written to a single FASTA which is then shuffled (k=2) and
scanned as one target set in each permutation iteration.

Sampling restricts to the canonical chromosomes (CANONICAL below), matching the
downstream analysis. The chromosome is parsed from the start of each header, up
to the first ':', whitespace or '('  (handles ">chr1:0-10000", ">chr1:0-10000(+)",
">chr1"). Plus/minus record headers are tagged with |+ / |- to keep names unique
in the pooled file (the target name is irrelevant to the null, which keys only on
the tRF query and score).

Usage:
    python3 sample_windows_from_fasta.py \
        windows_plus.fa windows_minus.fa out_windows.fa \
        --fraction 0.20 --seed 42

Requirements: Python 3.6+
"""

import argparse
import random
import re
import sys

# Canonical mouse chromosomes (mm10 / GRCm38), matching annotate_targets.R.
CANONICAL = ["chr" + c for c in [str(i) for i in range(1, 20)] + ["X", "Y"]]


def iter_fasta(path):
    header, parts = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].rstrip("\n"), []
            else:
                parts.append(line.rstrip("\n"))
    if header is not None:
        yield header, "".join(parts)


def chrom_of(header):
    """
    Chromosome name from a window header.

    Window FASTAs are produced by `bedtools getfasta -name`, giving headers of
    the form "<bed_name>::<chrom>:<start>-<end>" (the same format parsed in
    miranda_to_bed.R). Everything before '::' is the BED name, so it must be
    dropped before reading the chromosome. Plain "<chrom>:<start>-<end>" and
    bare "<chrom>" headers are also handled.
    """
    tok = header.split()[0]
    if "::" in tok:
        tok = tok.split("::", 1)[1]
    return re.split(r"[:\(\t ]", tok)[0]


def write_record(handle, header, seq, width=80):
    handle.write(f">{header}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def eligible_indices(path, chroms):
    """Return the list of record indices whose chromosome is in `chroms`."""
    return [i for i, (h, _s) in enumerate(iter_fasta(path)) if chrom_of(h) in chroms]


def diagnose(path, chroms, n=5):
    """Print the first few headers and the distinct parsed chromosome names."""
    headers, parsed = [], []
    for i, (h, _s) in enumerate(iter_fasta(path)):
        if i < n:
            headers.append(h)
        c = chrom_of(h)
        if c not in parsed:
            parsed.append(c)
        if len(parsed) >= 25 and i > 1000:
            break
    print(f"\n--- diagnostics for {path} ---", file=sys.stderr)
    print("First headers:", file=sys.stderr)
    for h in headers:
        print(f"    >{h}", file=sys.stderr)
    print(f"Parsed chromosome names (first {len(parsed)} distinct): "
          f"{parsed}", file=sys.stderr)
    print(f"Expected any of: {sorted(chroms)}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("windows_plus")
    ap.add_argument("windows_minus")
    ap.add_argument("out_windows")
    ap.add_argument("--fraction", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--chroms", default=None,
                    help="Comma-separated chromosomes to restrict to. "
                         "Default: the canonical set hardcoded in this script.")
    args = ap.parse_args()

    random.seed(args.seed)
    chroms = set(args.chroms.split(",")) if args.chroms else set(CANONICAL)

    # Pass 1: build the combined eligible pool, tagged by strand.
    plus_elig = eligible_indices(args.windows_plus, chroms)
    minus_elig = eligible_indices(args.windows_minus, chroms)
    pool = [("+", i) for i in plus_elig] + [("-", i) for i in minus_elig]
    if not pool:
        print("ERROR: no eligible windows after chromosome filtering.",
              file=sys.stderr)
        diagnose(args.windows_plus, chroms)
        sys.exit(1)

    n_sample = max(1, int(len(pool) * args.fraction))
    selected = set(random.sample(pool, n_sample))
    sel_plus = {i for s, i in selected if s == "+"}
    sel_minus = {i for s, i in selected if s == "-"}
    print(f"Pooled {len(plus_elig)} plus + {len(minus_elig)} minus canonical "
          f"windows; sampled {n_sample} ({args.fraction:.0%}): "
          f"{len(sel_plus)} plus, {len(sel_minus)} minus", file=sys.stderr)

    # Pass 2: write the selected records from each file, strand-tagged.
    with open(args.out_windows, "w") as out:
        for i, (header, seq) in enumerate(iter_fasta(args.windows_plus)):
            if i in sel_plus:
                write_record(out, header.split()[0] + "|+", seq)
        for i, (header, seq) in enumerate(iter_fasta(args.windows_minus)):
            if i in sel_minus:
                write_record(out, header.split()[0] + "|-", seq)

    print(f"Wrote {n_sample} pooled window targets -> {args.out_windows}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
