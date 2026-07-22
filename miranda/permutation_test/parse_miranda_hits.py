#!/usr/bin/env python3
"""
parse_miranda_hits.py

Turn miRanda -out files into a per-(tRF, score) COUNT histogram.

Why a histogram and not one row per hit: at a permissive score threshold the
shuffled null produces enormous numbers of hits (20% of the genome x hundreds
of tRFs x thousands of iterations = ~10^9 rows). We never need the individual
hits for null calibration -- only the score DISTRIBUTION per tRF. A histogram
(one row per tRF x integer score) is ~10^4 rows regardless of hit count, so it
is trivial to store and sum across iterations.

miRanda hit lines: each detected hit is written on a line beginning with a
single '>' (the per-pair summaries begin with '>>' and are skipped). The fields
are whitespace-separated:

    >Query  Ref  Score  Energy  Qstart Qend  Rstart Rend  AlnLen  %Id  %Sim

so field 0 (minus the '>') is the tRF/query id and field 2 is the alignment
score. Scores are rounded to the nearest integer to match the integer score
grid of the real prediction table.

Both strands are pooled (the search space includes both), so strand is ignored.

Usage:
    python3 parse_miranda_hits.py --out scores_hist.csv result_file1 result_file2 ...
    python3 parse_miranda_hits.py --out scores_hist.csv results_dir/result_*

Output CSV columns: tRF,alignment_score,count
"""

import argparse
import glob
import sys
from collections import Counter


def parse_file(path, counter):
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">") or line.startswith(">>"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            trf = parts[0].lstrip(">")
            try:
                score = int(round(float(parts[2])))
            except ValueError:
                continue
            counter[(trf, score)] += 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="miRanda -out files (globs allowed)")
    ap.add_argument("--out", required=True, help="Output histogram CSV")
    args = ap.parse_args()

    # expand any globs the shell did not
    files = []
    for pat in args.inputs:
        hits = glob.glob(pat)
        files.extend(hits if hits else [pat])

    counter = Counter()
    n_files = 0
    for f in files:
        try:
            parse_file(f, counter)
            n_files += 1
        except FileNotFoundError:
            print(f"WARNING: missing file {f}", file=sys.stderr)

    with open(args.out, "w") as out:
        out.write("tRF,alignment_score,count\n")
        for (trf, score) in sorted(counter):
            out.write(f"{trf},{score},{counter[(trf, score)]}\n")

    total = sum(counter.values())
    print(f"Parsed {n_files} files, {total} hits, "
          f"{len(counter)} (tRF,score) bins -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
