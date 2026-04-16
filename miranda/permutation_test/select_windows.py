#!/usr/bin/env python3
"""
select_windows.py
Select a random subset of non-overlapping 10 kbp windows from a genome.

Generates non-overlapping windows (step = window size) to avoid double-counting
hits in the overlap zones used by the full scanning pipeline. Then randomly
samples a fraction of those windows and writes a BED file. The BED file can
be fed to `bedtools getfasta` to extract sequences.

Usage:
    python select_windows.py <genome.fa.fai> <output.bed> \
        [--fraction 0.2] [--window 10000] [--seed 42]

Requirements:
    Python 3.6+
"""

import argparse
import random
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Sample non-overlapping genomic windows from a FAI index."
    )
    parser.add_argument("fai", help="Path to genome .fai index file")
    parser.add_argument("output_bed", help="Output BED file of selected windows")
    parser.add_argument(
        "--fraction",
        type=float,
        default=0.2,
        help="Fraction of windows to sample (default: 0.2)",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=10000,
        help="Window size in bp (default: 10000)",
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed (default: 42)"
    )
    args = parser.parse_args()

    random.seed(args.seed)
    window = args.window

    # Read chromosome lengths from FAI
    chroms = []
    with open(args.fai) as f:
        for line in f:
            fields = line.strip().split("\t")
            name = fields[0]
            length = int(fields[1])
            chroms.append((name, length))

    # Generate all non-overlapping windows
    all_windows = []
    for chrom, length in chroms:
        start = 0
        while start + window <= length:
            all_windows.append((chrom, start, start + window))
            start += window
        # Include the final partial window only if it's at least half-size
        if start < length and (length - start) >= window // 2:
            all_windows.append((chrom, start, length))

    n_total = len(all_windows)
    n_sample = max(1, int(n_total * args.fraction))

    selected = sorted(random.sample(all_windows, n_sample))

    with open(args.output_bed, "w") as out:
        for chrom, start, end in selected:
            # BED name encodes coordinates for downstream parsing
            name = f"{chrom}_{start + 1}"
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    print(
        f"Selected {n_sample} of {n_total} non-overlapping {window} bp windows "
        f"({args.fraction * 100:.0f}%)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
