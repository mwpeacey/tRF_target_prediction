#!/usr/bin/env python3
"""
select_windows.py
Select a random subset of non-overlapping 10 kbp windows from a genome,
excluding windows with high N content.

Generates non-overlapping windows (step = window size) to avoid double-counting
hits in the overlap zones used by the full scanning pipeline. The genome FASTA
is read in a single pass to compute per-window N fractions efficiently. Windows
where >50% of bases are N are excluded before sampling. A random fraction of the
remaining windows is selected and written as a BED file for `bedtools getfasta`.

Usage:
    python select_windows.py <genome.fa> <output.bed> \
        [--fraction 0.2] [--window 10000] [--max-n-frac 0.5] [--seed 42]

Requirements:
    Python 3.6+
"""

import argparse
import random
import sys


def parse_fasta_windows(genome_fa, window_size, max_n_frac, chroms=None):
    """
    Read genome FASTA in a single pass. For each chromosome, accumulate
    sequence and evaluate non-overlapping windows on the fly.

    chroms : optional set of chromosome names to restrict sampling to. When
    given, sequences whose name is not in the set are skipped entirely (not
    accumulated), so the null search space matches the real analysis (e.g. the
    canonical chromosomes, excluding unplaced scaffolds / alt contigs / chrM).

    Returns a list of eligible (chrom, start, end) tuples and counts of
    total/excluded windows.
    """
    eligible = []
    n_total = 0
    n_excluded = 0

    current_chrom = None
    keep = True
    seq_buf = []
    buf_len = 0

    def flush_windows(chrom, seq_buf, buf_len):
        """Process all complete windows from the accumulated sequence."""
        nonlocal eligible, n_total, n_excluded
        seq = "".join(seq_buf).upper()
        start = 0
        while start + window_size <= buf_len:
            chunk = seq[start : start + window_size]
            n_frac = chunk.count("N") / window_size
            n_total += 1
            if n_frac <= max_n_frac:
                eligible.append((chrom, start, start + window_size))
            else:
                n_excluded += 1
            start += window_size

    with open(genome_fa) as f:
        for line in f:
            if line.startswith(">"):
                # Process previous chromosome (only if it was kept)
                if current_chrom is not None and keep:
                    flush_windows(current_chrom, seq_buf, buf_len)
                current_chrom = line[1:].split()[0]
                keep = (chroms is None) or (current_chrom in chroms)
                seq_buf = []
                buf_len = 0
            elif keep:
                stripped = line.rstrip("\n")
                seq_buf.append(stripped)
                buf_len += len(stripped)

        # Final chromosome
        if current_chrom is not None and keep:
            flush_windows(current_chrom, seq_buf, buf_len)

    return eligible, n_total, n_excluded


def main():
    parser = argparse.ArgumentParser(
        description="Sample non-overlapping genomic windows, filtering high-N regions."
    )
    parser.add_argument("genome_fa", help="Path to genome FASTA")
    parser.add_argument("output_bed", help="Output BED file of selected windows")
    parser.add_argument(
        "--fraction",
        type=float,
        default=0.2,
        help="Fraction of eligible windows to sample (default: 0.2)",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=10000,
        help="Window size in bp (default: 10000)",
    )
    parser.add_argument(
        "--max-n-frac",
        type=float,
        default=0.5,
        help="Maximum fraction of N bases allowed per window (default: 0.5)",
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed (default: 42)"
    )
    parser.add_argument(
        "--chroms",
        default=None,
        help="Comma-separated chromosomes to restrict sampling to "
        "(e.g. chr1,chr2,...,chrX,chrY). Default: all sequences in the FASTA. "
        "Use this to match the null search space to the real analysis.",
    )
    args = parser.parse_args()

    random.seed(args.seed)

    chroms = set(args.chroms.split(",")) if args.chroms else None
    if chroms:
        print(f"Restricting to {len(chroms)} chromosomes: "
              f"{','.join(sorted(chroms))}", file=sys.stderr)

    print("Reading genome and filtering windows...", file=sys.stderr)

    eligible, n_total, n_excluded = parse_fasta_windows(
        args.genome_fa, args.window, args.max_n_frac, chroms=chroms
    )

    print(
        f"Generated {n_total} non-overlapping {args.window} bp windows, "
        f"excluded {n_excluded} for >{args.max_n_frac * 100:.0f}% N content, "
        f"{len(eligible)} eligible",
        file=sys.stderr,
    )

    # Sample
    n_sample = max(1, int(len(eligible) * args.fraction))
    selected = sorted(random.sample(eligible, n_sample))

    with open(args.output_bed, "w") as out:
        for chrom, start, end in selected:
            name = f"{chrom}_{start + 1}"
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    print(
        f"Selected {n_sample} of {len(eligible)} eligible windows "
        f"({args.fraction * 100:.0f}%)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
