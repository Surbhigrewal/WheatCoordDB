#!/usr/bin/env python3
"""
02_build_conversion_table.py
============================
Build coordinate conversion tables from anchor pairs.

For each CS chromosome, creates a dense lookup table mapping
CS positions → target assembly positions using monotone cubic
spline interpolation between gene anchors.

If the anchor file contains a 'segment' column (produced by
01_extract_anchors.py with --translocation-map), each CS chromosome
is processed as one or more segments, each mapping to its own
target chromosome. The segments are stitched into a single
conversion table with a heterogeneous tgt_chr column.

Outputs (per chromosome):
  <outdir>/<chr>_conversion.tsv.gz
      cs_pos | tgt_chr | tgt_pos | anchor_density_5Mb | interp_type

  <outdir>/<chr>_anchors.tsv
      The clean anchor pairs used (all segments)

  <outdir>/<chr>_synteny_blocks.tsv + <chr>_synteny.bed
      For visualisation in genome browsers

Also generates QC plots if matplotlib is available.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator


WHEAT_CHROMOSOMES = [f"{g}{s}" for g in range(1, 8) for s in ["A", "B", "D"]]


# ─────────────────────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────────────────────

def load_chrom_sizes(chrom_sizes_file):
    """Load chromosome sizes into a dict. Accepts 2-col TSV or FASTA .fai."""
    sizes = {}
    with open(chrom_sizes_file) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def get_chr_len(cs_chr, cs_sizes):
    """Look up CS chromosome length, trying Chr/chr prefix variants."""
    for key in [cs_chr, f"Chr{cs_chr}", f"chr{cs_chr}",
                cs_chr.replace("Chr", ""), cs_chr.replace("chr", "")]:
        if key in cs_sizes:
            return cs_sizes[key]
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Spline building
# ─────────────────────────────────────────────────────────────────────────────

def build_spline_converter(cs_anchors, tgt_anchors, cs_start, cs_end, tgt_chr_len):
    """
    Build a PCHIP spline interpolator cs_pos -> tgt_pos for one syntenic segment.

    Handles both forward (tgt increases with cs) and inverted (tgt decreases
    with cs) orientations natively — PCHIP works on monotone sequences in
    either direction.

    Parameters
    ----------
    cs_anchors  : array of CS midpoint positions (unsorted ok)
    tgt_anchors : array of corresponding target midpoint positions
    cs_start    : CS coordinate where this segment begins (for boundary anchor)
    cs_end      : CS coordinate where this segment ends   (for boundary anchor)
    tgt_chr_len : length of target chromosome (for clipping)

    Returns
    -------
    interpolator : callable — f(cs_pos_array) -> tgt_pos_array
    mean_gap     : float — mean inter-anchor distance in CS bp
    interp_type  : str
    """
    sort_idx   = np.argsort(cs_anchors)
    cs_sorted  = cs_anchors[sort_idx]
    tgt_sorted = tgt_anchors[sort_idx]

    # Determine orientation from majority of consecutive differences
    diffs    = np.diff(tgt_sorted)
    inverted = np.sum(diffs < 0) > np.sum(diffs >= 0)

    # Extrapolate boundary anchors linearly from the first/last two real anchors
    if len(cs_sorted) >= 2:
        slope_l   = (tgt_sorted[1]  - tgt_sorted[0])  / max(cs_sorted[1]  - cs_sorted[0],  1)
        slope_r   = (tgt_sorted[-1] - tgt_sorted[-2]) / max(cs_sorted[-1] - cs_sorted[-2], 1)
        tgt_left  = tgt_sorted[0]  - slope_l * (cs_sorted[0]  - cs_start)
        tgt_right = tgt_sorted[-1] + slope_r * (cs_end         - cs_sorted[-1])
    else:
        tgt_left  = float(tgt_sorted[0])
        tgt_right = float(tgt_sorted[0])

    tgt_left  = float(np.clip(tgt_left,  0, tgt_chr_len))
    tgt_right = float(np.clip(tgt_right, 0, tgt_chr_len))

    cs_all  = np.concatenate([[cs_start], cs_sorted, [cs_end]])
    tgt_all = np.concatenate([[tgt_left], tgt_sorted, [tgt_right]])

    # Remove duplicate x values (can occur for near-identical assemblies)
    cs_all, unique_idx = np.unique(cs_all, return_index=True)
    tgt_all = tgt_all[unique_idx]
    tgt_all = np.clip(tgt_all, 0, tgt_chr_len)

    # Check final monotonicity for PCHIP
    inc = bool(np.all(np.diff(tgt_all) >= 0))
    dec = bool(np.all(np.diff(tgt_all) <= 0))

    if inc:
        interpolator = PchipInterpolator(cs_all, tgt_all, extrapolate=True)
        interp_type  = "pchip"
    elif dec:
        # PCHIP handles monotone decreasing sequences directly
        interpolator = PchipInterpolator(cs_all, tgt_all, extrapolate=True)
        interp_type  = "pchip_inverted"
    else:
        # Non-monotone (fragmented/noisy) — safe linear fallback
        _cs, _tgt = cs_all.copy(), tgt_all.copy()
        interpolator = lambda x, c=_cs, t=_tgt: np.interp(x, c, t)
        interp_type  = "linear"

    mean_gap = float(np.mean(np.diff(cs_sorted))) if len(cs_sorted) > 1 else float(cs_end - cs_start)
    return interpolator, mean_gap, interp_type


# ─────────────────────────────────────────────────────────────────────────────
# Anchor density
# ─────────────────────────────────────────────────────────────────────────────

def compute_local_anchor_density(cs_pos_query, cs_anchors, window=5_000_000):
    """Number of anchors within ±window bp of each query position."""
    densities = np.zeros(len(cs_pos_query), dtype=int)
    for i, pos in enumerate(cs_pos_query):
        densities[i] = int(np.sum(
            (cs_anchors >= pos - window) & (cs_anchors <= pos + window)
        ))
    return densities


# ─────────────────────────────────────────────────────────────────────────────
# Synteny blocks
# ─────────────────────────────────────────────────────────────────────────────

def build_synteny_blocks(anchors_df, min_block_genes=3, max_gap=2_000_000):
    """Merge adjacent anchors into synteny blocks; respects tgt_chr boundaries."""
    if len(anchors_df) == 0:
        return pd.DataFrame()

    df = anchors_df.sort_values("cs_midpoint").reset_index(drop=True)
    blocks      = []
    block_start = 0

    for i in range(1, len(df)):
        gap         = df.loc[i, "cs_midpoint"] - df.loc[i-1, "cs_midpoint"]
        chr_changed = df.loc[i, "tgt_chr"] != df.loc[i-1, "tgt_chr"]
        last_row    = (i == len(df) - 1)

        if gap > max_gap or chr_changed or last_row:
            block_end  = i if (last_row and gap <= max_gap and not chr_changed) else i - 1
            block_rows = df.loc[block_start:block_end]

            if len(block_rows) >= min_block_genes:
                blocks.append({
                    "cs_chr":        block_rows["cs_chr"].iloc[0],
                    "cs_start":      int(block_rows["cs_midpoint"].min()),
                    "cs_end":        int(block_rows["cs_midpoint"].max()),
                    "tgt_chr":       block_rows["tgt_chr"].mode()[0],
                    "tgt_start":     int(block_rows["tgt_midpoint"].min()),
                    "tgt_end":       int(block_rows["tgt_midpoint"].max()),
                    "n_anchors":     len(block_rows),
                    "mean_identity": float(block_rows["identity"].mean()),
                })

            block_start = i

    return pd.DataFrame(blocks)


# ─────────────────────────────────────────────────────────────────────────────
# QC dotplot
# ─────────────────────────────────────────────────────────────────────────────

def plot_synteny_dots(anchors_df, cs_chr, target_name, output_path):
    """Dot-plot QC figure; each tgt_chr plotted in a distinct colour."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        cmap    = plt.cm.plasma_r
        id_col  = next((c for c in anchors_df.columns if c.endswith("identity")), None)
        tgt_chrs = sorted(anchors_df["tgt_chr"].unique())
        colours  = plt.cm.tab10(np.linspace(0, 0.9, max(len(tgt_chrs), 1)))

        fig, ax = plt.subplots(figsize=(8, 6))

        for colour, tgt_chr in zip(colours, tgt_chrs):
            grp = anchors_df[anchors_df["tgt_chr"] == tgt_chr]
            if id_col:
                ax.scatter(
                    grp["cs_midpoint"] / 1e6,
                    grp["tgt_midpoint"] / 1e6,
                    c=grp[id_col].clip(0.9, 1.0),
                    cmap=cmap, vmin=0.9, vmax=1.0,
                    s=2, alpha=0.6, linewidths=0, rasterized=True,
                    label=tgt_chr,
                )
            else:
                ax.scatter(
                    grp["cs_midpoint"] / 1e6,
                    grp["tgt_midpoint"] / 1e6,
                    color=colour,
                    s=2, alpha=0.6, linewidths=0, rasterized=True,
                    label=tgt_chr,
                )

        ax.set_xlabel(f"CS {cs_chr} position (Mb)", fontsize=11)
        ax.set_ylabel(f"{target_name} position (Mb)", fontsize=11)
        ax.set_title(f"Synteny: CS {cs_chr} → {target_name}", fontsize=12)
        ax.legend(loc="upper left", fontsize=8, markerscale=4)

        if id_col:
            sm = plt.cm.ScalarMappable(cmap=cmap,
                                       norm=mcolors.Normalize(vmin=0.9, vmax=1.0))
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
            cbar.set_label("Sequence identity", fontsize=9)

        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()
    except ImportError:
        pass


# ─────────────────────────────────────────────────────────────────────────────
# Segment definition
# ─────────────────────────────────────────────────────────────────────────────

def _define_segments(grp, cs_chr_len):
    """
    From a group with a 'segment' column, define ordered segment dicts with
    cs_start, cs_end, tgt_chr, tgt_chr_len.

    Segments are ordered by median CS position. CS coordinate space is
    partitioned at the midpoint between adjacent segment anchor clouds so
    the stitched conversion table has no gaps or overlaps.
    """
    seg_info = []
    for seg_label, seg_grp in grp.groupby("segment"):
        tgt_chr = seg_grp["tgt_chr"].mode()[0]
        seg_info.append({
            "label":       seg_label,
            "tgt_chr":     tgt_chr,
            "tgt_chr_len": int(seg_grp["tgt_end"].max() * 1.1),
            "cs_median":   float(seg_grp["cs_midpoint"].median()),
            "cs_min":      float(seg_grp["cs_midpoint"].min()),
            "cs_max":      float(seg_grp["cs_midpoint"].max()),
            "anchors":     seg_grp.copy(),
        })

    seg_info.sort(key=lambda s: s["cs_median"])

    for i, seg in enumerate(seg_info):
        if i == 0:
            seg["cs_start"] = 0
        else:
            gap_mid = int((seg_info[i-1]["cs_max"] + seg["cs_min"]) / 2)
            seg["cs_start"]         = gap_mid
            seg_info[i-1]["cs_end"] = gap_mid

        if i == len(seg_info) - 1:
            seg["cs_end"] = cs_chr_len

    return seg_info


# ─────────────────────────────────────────────────────────────────────────────
# Per-chromosome processing
# ─────────────────────────────────────────────────────────────────────────────

def process_chromosome(cs_chr, grp, cs_sizes, resolution, min_anchors,
                       target_name, outdir, plot_dir):
    """
    Build conversion table for one CS chromosome.

    If 'segment' column is present AND more than one unique segment exists,
    each segment is interpolated independently then stitched into one table.
    """
    cs_chr_len = get_chr_len(cs_chr, cs_sizes)
    if cs_chr_len == 0:
        print(f"  WARN: no size found for {cs_chr}, estimating from anchors",
              file=sys.stderr)
        cs_chr_len = int(grp["cs_end"].max() * 1.05)

    has_segments = ("segment" in grp.columns and grp["segment"].nunique() > 1)

    if has_segments:
        segments = _define_segments(grp, cs_chr_len)
    else:
        tgt_chr     = grp["tgt_chr"].mode()[0]
        seg_grp     = grp[grp["tgt_chr"] == tgt_chr].copy()
        tgt_chr_len = int(seg_grp["tgt_end"].max() * 1.1)
        segments    = [{"label":       "collinear",
                        "cs_start":    0,
                        "cs_end":      cs_chr_len,
                        "tgt_chr":     tgt_chr,
                        "tgt_chr_len": tgt_chr_len,
                        "anchors":     seg_grp}]

    total_anchors = sum(len(s["anchors"]) for s in segments)
    if total_anchors < min_anchors:
        print(f"  SKIP {cs_chr}: {total_anchors} anchors < min {min_anchors}",
              file=sys.stderr)
        return None

    conv_parts   = []
    all_anchors  = []
    summary_segs = []

    for seg in segments:
        seg_anchors = seg["anchors"].sort_values("cs_midpoint").reset_index(drop=True)
        n = len(seg_anchors)

        if n < 3:
            print(f"  SKIP segment '{seg['label']}' ({cs_chr}): only {n} anchors",
                  file=sys.stderr)
            continue

        interp, mean_gap, interp_type = build_spline_converter(
            cs_anchors  = seg_anchors["cs_midpoint"].values,
            tgt_anchors = seg_anchors["tgt_midpoint"].values,
            cs_start    = seg["cs_start"],
            cs_end      = seg["cs_end"],
            tgt_chr_len = seg["tgt_chr_len"],
        )

        print(f"  {cs_chr} [{seg['label']}] -> {seg['tgt_chr']}: "
              f"{n} anchors, gap {mean_gap/1e6:.2f} Mb, {interp_type}",
              file=sys.stderr)

        cs_pos = np.arange(seg["cs_start"], seg["cs_end"] + resolution, resolution)
        cs_pos = cs_pos[cs_pos <= seg["cs_end"]]

        tgt_pos = np.clip(
            np.round(interp(cs_pos)).astype(int), 0, seg["tgt_chr_len"]
        )
        density = compute_local_anchor_density(cs_pos, seg_anchors["cs_midpoint"].values)

        conv_parts.append(pd.DataFrame({
            "cs_chr":             cs_chr,
            "cs_pos":             cs_pos,
            "tgt_chr":            seg["tgt_chr"],
            "tgt_pos":            tgt_pos,
            "anchor_density_5Mb": density,
            "interp_type":        interp_type,
        }))
        all_anchors.append(seg_anchors)
        summary_segs.append({
            "segment":    seg["label"],
            "tgt_chr":    seg["tgt_chr"],
            "n_anchors":  n,
            "mean_gap_Mb": round(mean_gap / 1e6, 3),
            "interp_type": interp_type,
        })

    if not conv_parts:
        return None

    conv_df     = pd.concat(conv_parts, ignore_index=True).sort_values("cs_pos")
    all_anch_df = pd.concat(all_anchors, ignore_index=True)

    # Save conversion table
    conv_file = outdir / f"{cs_chr}_conversion.tsv.gz"
    conv_df.to_csv(conv_file, sep="\t", index=False, compression="gzip")
    print(f"  Saved: {conv_file}", file=sys.stderr)

    # Save anchors (all segments)
    all_anch_df.to_csv(outdir / f"{cs_chr}_anchors.tsv", sep="\t", index=False)

    # Synteny blocks
    blocks = build_synteny_blocks(all_anch_df)
    if len(blocks) > 0:
        blocks.to_csv(outdir / f"{cs_chr}_synteny_blocks.tsv", sep="\t", index=False)
        bed_file = outdir / f"{cs_chr}_synteny.bed"
        with open(bed_file, "w") as bf:
            bf.write("# CS_chr\tCS_start\tCS_end\ttgt_chr\ttgt_start\ttgt_end"
                     "\tn_anchors\tmean_identity\n")
            for _, b in blocks.iterrows():
                bf.write(f"{b['cs_chr']}\t{b['cs_start']}\t{b['cs_end']}\t"
                         f"{b['tgt_chr']}\t{b['tgt_start']}\t{b['tgt_end']}\t"
                         f"{b['n_anchors']}\t{b['mean_identity']:.4f}\n")

    # Dotplot
    if plot_dir:
        plot_synteny_dots(all_anch_df, cs_chr, target_name,
                         Path(plot_dir) / f"{cs_chr}_dotplot.png")

    n_trans = sum(s["n_anchors"] for s in summary_segs if s["segment"] != "collinear")
    return {
        "cs_chr":                  cs_chr,
        "tgt_chr":                 "/".join(s["tgt_chr"] for s in summary_segs),
        "n_anchors":               sum(s["n_anchors"] for s in summary_segs),
        "n_anchors_collinear":     sum(s["n_anchors"] for s in summary_segs
                                       if s["segment"] == "collinear"),
        "n_anchors_translocation": n_trans,
        "has_translocation":       n_trans > 0,
        "mean_anchor_gap_Mb":      round(
            float(np.mean([s["mean_gap_Mb"] for s in summary_segs])), 3),
        "interp_type":             "/".join(s["interp_type"] for s in summary_segs),
        "n_synteny_blocks":        len(blocks) if len(blocks) > 0 else 0,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Build coordinate conversion tables from anchor pairs"
    )
    parser.add_argument("--anchors",             required=True,
                        help="Anchor TSV from 01_extract_anchors.py")
    parser.add_argument("--cs-chrom-sizes",      required=True,
                        help="CS chromosome sizes (2-col TSV or .fai)")
    parser.add_argument("--output-dir",          required=True,
                        help="Output directory")
    parser.add_argument("--target-name",         required=True,
                        help="Target assembly name (used in plot titles)")
    parser.add_argument("--resolution",          type=int, default=1000,
                        help="Output table resolution in bp (default: 1000)")
    parser.add_argument("--min-anchors-per-chr", type=int, default=20,
                        help="Min total anchors to build table (default: 20)")
    parser.add_argument("--plot-dir",            default=None,
                        help="Output directory for QC dotplots")
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.plot_dir:
        Path(args.plot_dir).mkdir(parents=True, exist_ok=True)

    print(f"Loading anchors: {args.anchors}", file=sys.stderr)
    anchors = pd.read_csv(args.anchors, sep="\t")
    print(f"  Total anchors: {len(anchors)}", file=sys.stderr)

    has_seg = "segment" in anchors.columns
    print(f"  Segment column present: {has_seg}", file=sys.stderr)
    if has_seg:
        print(anchors["segment"].value_counts().to_string(), file=sys.stderr)

    cs_sizes = load_chrom_sizes(args.cs_chrom_sizes)

    chr_summary = []
    for cs_chr, grp in anchors.groupby("cs_chr"):
        print(f"\nProcessing {cs_chr}...", file=sys.stderr)
        result = process_chromosome(
            cs_chr      = cs_chr,
            grp         = grp,
            cs_sizes    = cs_sizes,
            resolution  = args.resolution,
            min_anchors = args.min_anchors_per_chr,
            target_name = args.target_name,
            outdir      = outdir,
            plot_dir    = args.plot_dir,
        )
        if result:
            chr_summary.append(result)

    summary_df = pd.DataFrame(chr_summary)
    summary_df.to_csv(outdir / "conversion_summary.tsv", sep="\t", index=False)

    print(f"\n{'='*60}", file=sys.stderr)
    print(f"SUMMARY for {args.target_name}:", file=sys.stderr)
    print(summary_df.to_string(index=False), file=sys.stderr)
    print(f"\nAll outputs in: {outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
