#!/usr/bin/env python3
"""
select_windows.py
Select a random subset of non-overlapping 10 kbp windows from a genome,
excluding windows with high N content.

Generates non-overlapping windows (step = window size) to avoid double-counting
hits in the overlap zones used by the full scanning pipeline. Windows where >50%
of bases are N are excluded before sampling. A random fraction of the remaining
windows is selected and written as a BED file for `bedtools getfasta`.

Usage:
    python select_windows.py <genome.fa> <output.bed> \
        [--fraction 0.2] [--window 10000] [--max-n-frac 0.5] [--seed 42]

Requirements:
    Python 3.6+, samtools (for indexing if .fai missing)
"""

import argparse
import os
import random
import subprocess
import sys


def read_fai(fai_path):
    """Read chromosome names and lengths from a .fai index."""
    chroms = []
    with open(fai_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            chroms.append((fields[0], int(fields[1])))
    return chroms


def get_n_fraction(genome_fa, chrom, start, end):
    """
    Use samtools faidx to extract a region and compute its N fraction.
    """
    region = f"{chrom}:{start + 1}-{end}"  # samtools uses 1-based coordinates
    result = subprocess.run(
        ["samtools", "faidx", genome_fa, region],
        capture_output=True,
        text=True,
    )
    seq = "".join(
        line.strip().upper()
        for line in result.stdout.split("\n")
        if not line.startswith(">")
    )
    if len(seq) == 0:
        return 1.0
    return seq.count("N") / len(seq)


def main():
    parser = argparse.ArgumentParser(
        description="Sample non-overlapping genomic windows, filtering high-N regions."
    )
    parser.add_argument("genome_fa", help="Path to genome FASTA (must have .fai index)")
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
    args = parser.parse_args()

    random.seed(args.seed)
    window = args.window
    fai_path = args.genome_fa + ".fai"

    # Ensure genome is indexed
    if not os.path.exists(fai_path):
        print("Indexing genome...", file=sys.stderr)
        subprocess.run(["samtools", "faidx", args.genome_fa], check=True)

    chroms = read_fai(fai_path)

    # Generate all non-overlapping full-size windows
    all_windows = []
    for chrom, length in chroms:
        start = 0
        while start + window <= length:
            all_windows.append((chrom, start, start + window))
            start += window

    print(
        f"Generated {len(all_windows)} non-overlapping {window} bp windows",
        file=sys.stderr,
    )

    # Filter out high-N windows
    print(
        f"Filtering windows with >{args.max_n_frac * 100:.0f}% N content...",
        file=sys.stderr,
    )
    eligible = []
    n_excluded = 0
    for i, (chrom, start, end) in enumerate(all_windows):
        n_frac = get_n_fraction(args.genome_fa, chrom, start, end)
        if n_frac <= args.max_n_frac:
            eligible.append((chrom, start, end))
        else:
            n_excluded += 1

        # Progress reporting
        if (i + 1) % 10000 == 0:
            print(
                f"  Checked {i + 1}/{len(all_windows)} windows, "
                f"excluded {n_excluded} so far...",
                file=sys.stderr,
            )

    print(
        f"Excluded {n_excluded} windows for N content, "
        f"{len(eligible)} eligible remain",
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
