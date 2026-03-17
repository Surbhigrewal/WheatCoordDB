#!/usr/bin/env python3
"""
00_compute_cs_gene_density.py
=============================
Compute CS gene density in a sliding window along each chromosome,
to be used as the expected gene count denominator when computing
anchor_fraction in 02_build_conversion_table.py.

For each 1 kb position (matching the conversion table resolution),
counts the number of CS genes whose midpoint falls within ±window bp.

Output: TSV with columns:
  cs_chr | cs_pos | cs_gene_density

Usage:
  python3 00_compute_cs_gene_density.py \\
      --cs-gff   outputs/renamed/CS/CS.gff3.gz \\
      --cs-chrom-sizes outputs/renamed/CS/CS.chrom.sizes \\
      --output   outputs/cs_gene_density.tsv.gz \\
      --window   5000000 \\
      --resolution 1000

This only needs to be run ONCE. The output is shared across all
24 assembly reruns.

Runtime: ~5 minutes for the full CS annotation (~107k genes).
"""

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd


WHEAT_CHROMOSOMES = {f"Chr{g}{s}" for g in range(1, 8) for s in ["A", "B", "D"]}


def _open(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")


def load_chrom_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def parse_cs_gene_midpoints(gff_path, feature_type="gene"):
    """
    Parse CS GFF3 and return dict: chr -> sorted array of gene midpoints.
    Only retains canonical wheat chromosomes (Chr1A-Chr7D).
    """
    midpoints = {c: [] for c in WHEAT_CHROMOSOMES}
    n_total = 0
    n_kept = 0

    with _open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _, ftype, start, end = parts[0], parts[1], parts[2], parts[3], parts[4]
            if ftype != feature_type:
                continue
            n_total += 1
            if seqid not in WHEAT_CHROMOSOMES:
                continue
            mid = (int(start) + int(end)) // 2
            midpoints[seqid].append(mid)
            n_kept += 1

    print(f"  Parsed {n_total} {feature_type} features, kept {n_kept} on canonical chromosomes",
          file=sys.stderr)

    # Sort each chromosome
    for c in midpoints:
        midpoints[c] = np.array(sorted(midpoints[c]), dtype=np.int64)

    return midpoints


def compute_gene_density(midpoints_arr, query_positions, window):
    """
    For each query position, count gene midpoints within ±window bp.
    Uses searchsorted for efficiency — O(n log n) rather than O(n*m).
    """
    lo = np.searchsorted(midpoints_arr, query_positions - window, side="left")
    hi = np.searchsorted(midpoints_arr, query_positions + window, side="right")
    return (hi - lo).astype(np.int32)


def main():
    parser = argparse.ArgumentParser(
        description="Compute CS gene density in sliding window for each 1kb position"
    )
    parser.add_argument("--cs-gff",          required=True,
                        help="CS v2.1 annotation GFF3 (plain or .gz)")
    parser.add_argument("--cs-chrom-sizes",  required=True,
                        help="CS chromosome sizes file (two-column TSV)")
    parser.add_argument("--output",          required=True,
                        help="Output TSV.gz path")
    parser.add_argument("--window",          type=int, default=5_000_000,
                        help="Window size in bp (default: 5000000 = ±5 Mb)")
    parser.add_argument("--resolution",      type=int, default=1000,
                        help="Position step size in bp (default: 1000 = 1 kb)")
    parser.add_argument("--feature-type",    default="gene",
                        help="GFF3 feature type to count (default: gene)")
    args = parser.parse_args()

    print(f"Parsing CS GFF: {args.cs_gff}", file=sys.stderr)
    midpoints = parse_cs_gene_midpoints(args.cs_gff, feature_type=args.feature_type)

    print(f"Loading chrom sizes: {args.cs_chrom_sizes}", file=sys.stderr)
    cs_sizes = load_chrom_sizes(args.cs_chrom_sizes)

    rows = []
    for cs_chr in sorted(WHEAT_CHROMOSOMES):
        if cs_chr not in cs_sizes:
            print(f"  SKIP {cs_chr}: not in chrom sizes", file=sys.stderr)
            continue

        chr_len = cs_sizes[cs_chr]
        mids    = midpoints.get(cs_chr, np.array([], dtype=np.int64))
        n_genes = len(mids)

        if n_genes == 0:
            print(f"  WARN {cs_chr}: no genes found", file=sys.stderr)
            continue

        query_pos = np.arange(0, chr_len + args.resolution, args.resolution, dtype=np.int64)
        density   = compute_gene_density(mids, query_pos, args.window)

        print(f"  {cs_chr}: {n_genes} genes, {len(query_pos)} positions, "
              f"mean density {density.mean():.1f} genes/window", file=sys.stderr)

        chr_rows = pd.DataFrame({
            "cs_chr":          cs_chr,
            "cs_pos":          query_pos,
            "cs_gene_density": density,
        })
        rows.append(chr_rows)

    print(f"\nWriting output: {args.output}", file=sys.stderr)
    out_df = pd.concat(rows, ignore_index=True)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(f"  Written {len(out_df):,} rows", file=sys.stderr)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
