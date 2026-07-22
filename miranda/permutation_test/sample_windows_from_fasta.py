#!/usr/bin/env python3
"""
sample_windows_from_fasta.py

Sample a random subset of pre-generated genomic windows for the permutation null.

Instead of re-windowing the genome (select_windows.py + bedtools getfasta), this
samples directly from the SAME window FASTAs that were used for the real miRanda
scan. That guarantees the null and the real analysis share an identical search
space: same window boundaries, coordinates, strand handling and N-masking.

Plus- and minus-strand window files are matched by RECORD ORDER (the Nth record
of the plus file corresponds to the Nth record of the minus file), which is how
the minus file is produced (reverse-complement of the plus, order preserved). The
two files must therefore have the same number of records; this is checked.

Optional --chroms restricts sampling to the given chromosomes (e.g. the canonical
set), so the null matches the chromosomes your downstream analysis keeps. The
chromosome is parsed from the start of each window header, up to the first ':',
whitespace or '('  (handles headers like ">chr1:0-10000", ">chr1:0-10000(+)",
">chr1").

Usage:
    python3 sample_windows_from_fasta.py \
        windows_plus.fa windows_minus.fa \
        out_plus.fa out_minus.fa \
        --fraction 0.20 --seed 42 [--chroms chr1,chr2,...,chrX,chrY]

Requirements: Python 3.6+
"""

import argparse
import random
import re
import sys


def iter_fasta(path):
    """Yield (header_line_without_>, sequence_string) for each record."""
    header = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].rstrip("\n")
                parts = []
            else:
                parts.append(line.rstrip("\n"))
    if header is not None:
        yield header, "".join(parts)


def chrom_of(header):
    """Extract chromosome name from a window header."""
    tok = header.split()[0]                 # first whitespace-delimited token
    return re.split(r"[:\(\t ]", tok)[0]    # up to ':' '(' or whitespace


def write_record(handle, header, seq, width=80):
    handle.write(f">{header}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("windows_plus")
    ap.add_argument("windows_minus")
    ap.add_argument("out_plus")
    ap.add_argument("out_minus")
    ap.add_argument("--fraction", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--chroms", default=None,
                    help="Comma-separated chromosomes to restrict to (default: all)")
    args = ap.parse_args()

    random.seed(args.seed)
    chroms = set(args.chroms.split(",")) if args.chroms else None

    # ── Pass 1: index every plus record, note its chromosome ─────────────────
    idx_chrom = []
    for i, (header, _seq) in enumerate(iter_fasta(args.windows_plus)):
        idx_chrom.append((i, chrom_of(header)))
    n_total = len(idx_chrom)

    eligible = [i for i, c in idx_chrom if (chroms is None or c in chroms)]
    if chroms is not None:
        seen = sorted({c for _, c in idx_chrom})
        kept = sorted({c for i, c in idx_chrom if c in chroms})
        print(f"Chromosome filter: kept {kept}", file=sys.stderr)
        missing = chroms - set(seen)
        if missing:
            print(f"WARNING: requested chromosomes not found in windows: "
                  f"{sorted(missing)}", file=sys.stderr)
    if not eligible:
        sys.exit("ERROR: no eligible windows after chromosome filtering — check "
                 "the window header format with `grep '>' <file> | head`.")

    n_sample = max(1, int(len(eligible) * args.fraction))
    selected = set(random.sample(eligible, n_sample))
    print(f"Sampled {n_sample} of {len(eligible)} eligible windows "
          f"({args.fraction:.0%}); {n_total} total records", file=sys.stderr)

    # ── Pass 2: write selected plus records ──────────────────────────────────
    with open(args.out_plus, "w") as out:
        for i, (header, seq) in enumerate(iter_fasta(args.windows_plus)):
            if i in selected:
                write_record(out, header, seq)

    # ── Pass 3: write the SAME indices from the minus file ───────────────────
    n_minus = 0
    with open(args.out_minus, "w") as out:
        for i, (header, seq) in enumerate(iter_fasta(args.windows_minus)):
            n_minus += 1
            if i in selected:
                write_record(out, header, seq)

    if n_minus != n_total:
        sys.exit(f"ERROR: plus ({n_total}) and minus ({n_minus}) window files have "
                 f"different record counts — cannot match by order. Ensure the minus "
                 f"file is the reverse-complement of the plus file in the same order.")

    print(f"Wrote {n_sample} plus + {n_sample} minus windows.", file=sys.stderr)


if __name__ == "__main__":
    main()
