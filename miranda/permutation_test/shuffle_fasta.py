#!/usr/bin/env python3
"""
shuffle_fasta.py
Dinucleotide-preserving shuffle of each sequence in a FASTA file.

Uses the uShuffle algorithm (Jiang et al., 2008) to shuffle sequences while
preserving exact dinucleotide (k=2) composition. This is critical for
generating appropriate null sequences from genomic DNA, where CpG depletion
and other dinucleotide biases would be destroyed by mononucleotide shuffling.

Sequences containing ambiguous bases (N) are handled by splitting on runs of
N, shuffling each unambiguous segment independently, and reassembling. Windows
that are mostly N (>50%) are written through unchanged and flagged on stderr.

Usage:
    python shuffle_fasta.py <input.fa> <output.fa> [--seed 42] [--klet 2]

Requirements:
    ushuffle (conda install -c bioconda ushuffle)
"""

import argparse
import random
import re
import sys

import ushuffle


def parse_fasta(filepath):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_parts = []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def write_fasta(handle, header, sequence, line_width=80):
    """Write a single FASTA record with wrapped sequence lines."""
    handle.write(f">{header}\n")
    for i in range(0, len(sequence), line_width):
        handle.write(sequence[i : i + line_width] + "\n")


def shuffle_segment(seq_str, k=2):
    """Shuffle a contiguous ACGT segment preserving k-mer composition."""
    if len(seq_str) < k + 1:
        return seq_str
    seq_upper = seq_str.upper()
    shuffled = ushuffle.shuffle(seq_upper, len(seq_upper), k)
    # Handle both str and bytes return types across ushuffle versions
    if isinstance(shuffled, bytes):
        return shuffled.decode("ascii")
    return shuffled


def shuffle_with_ns(seq_str, k=2):
    """
    Shuffle a genomic sequence that may contain N runs.

    Splits on contiguous N blocks, shuffles each non-N segment independently
    (preserving dinucleotide composition), and reassembles with the original
    N blocks in place.
    """
    # Split into segments: alternating non-N and N blocks
    parts = re.split(r"(N+)", seq_str.upper())
    shuffled_parts = []
    for part in parts:
        if len(part) == 0:
            continue
        if part[0] == "N":
            shuffled_parts.append(part)
        else:
            shuffled_parts.append(shuffle_segment(part, k))
    return "".join(shuffled_parts)


def main():
    parser = argparse.ArgumentParser(
        description="Dinucleotide-preserving shuffle of FASTA sequences."
    )
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_fasta", help="Output shuffled FASTA file")
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed (default: None)"
    )
    parser.add_argument(
        "--klet",
        type=int,
        default=2,
        help="k-let size to preserve (default: 2 for dinucleotide)",
    )
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    n_skipped = 0
    n_shuffled = 0

    with open(args.output_fasta, "w") as out_handle:
        for header, seq_str in parse_fasta(args.input_fasta):
            seq_str = seq_str.upper()
            n_count = seq_str.count("N")
            n_frac = n_count / len(seq_str) if len(seq_str) > 0 else 1.0

            if n_frac > 0.5:
                # Mostly N — write unchanged
                write_fasta(out_handle, header, seq_str)
                n_skipped += 1
            else:
                shuffled_seq = shuffle_with_ns(seq_str, k=args.klet)
                write_fasta(out_handle, header, shuffled_seq)
                n_shuffled += 1

    print(
        f"Shuffled {n_shuffled} sequences, skipped {n_skipped} (>50% N)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
