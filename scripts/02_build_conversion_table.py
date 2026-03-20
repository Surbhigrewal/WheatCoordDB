#!/usr/bin/env python3
"""
02_build_conversion_table.py
============================
Build coordinate conversion tables from anchor pairs.

For each CS chromosome, creates a dense lookup table mapping
CS positions → target assembly positions using monotone cubic
spline interpolation between gene anchors.

Outputs (per chromosome):
  <outdir>/<chr>_conversion.tsv.gz
      cs_pos | tgt_chr | tgt_pos | anchor_density_5Mb |
      anchor_fraction | interp_type

  <outdir>/<chr>_anchors.tsv
      Clean anchor pairs used for this chromosome

  <outdir>/<chr>_synteny_blocks.tsv + <chr>_synteny.bed
      Merged synteny blocks for genome browser visualisation

  <outdir>/conversion_summary.tsv
      Per-chromosome summary statistics

Confidence metric — anchor_fraction (0.0–1.0):
  For each query position, within ±5 Mb:
    anchor_fraction = anchors_in_window / cs_genes_in_window

  Where cs_genes_in_window is precomputed from the CS v2.1 GFF by
  00_compute_cs_gene_density.py. This normalises anchor density by
  the expected gene content, so pericentromeric regions (low gene
  density) are not penalised relative to telomeric regions (high
  gene density). Values are capped at 1.0.

  Thresholds: >=0.8 High, 0.5-0.8 Moderate, <0.5 Low confidence.

  anchor_density_5Mb is retained as the raw count (secondary metric).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator


WHEAT_CHROMOSOMES = [f"{g}{s}" for g in range(1, 8) for s in ["A", "B", "D"]]


def load_chrom_sizes(chrom_sizes_file):
    sizes = {}
    with open(chrom_sizes_file) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def load_cs_gene_density(density_file, cs_chr):
    """
    Load precomputed CS gene density for a single chromosome.
    Returns a numpy array indexed by position // resolution.
    """
    df = pd.read_csv(density_file, sep="\t",
                     compression="gzip" if str(density_file).endswith(".gz") else "infer")
    df = df[df["cs_chr"] == cs_chr].reset_index(drop=True)
    if len(df) == 0:
        return None
    # Return as array sorted by cs_pos
    df = df.sort_values("cs_pos").reset_index(drop=True)
    return df["cs_pos"].values, df["cs_gene_density"].values


def compute_anchor_fraction(cs_pos_query, cs_anchors, cs_gene_density_pos,
                             cs_gene_density_vals, window=5_000_000):
    """
    Compute anchor_fraction for each query position:
      anchor_fraction = anchors_in_window / cs_genes_in_window

    Values capped at 1.0 (can't have more anchors than genes).
    Returns NaN where cs_gene_density == 0 (no genes in window).

    This normalises anchor density by the expected gene content,
    so pericentromeric regions (low gene density) are not penalised
    relative to telomeric regions (high gene density).
    """
    # Anchor counts in window (reuse same logic as compute_local_anchor_density)
    anchor_counts = np.zeros(len(cs_pos_query), dtype=float)
    for i, pos in enumerate(cs_pos_query):
        anchor_counts[i] = np.sum(
            (cs_anchors >= pos - window) & (cs_anchors <= pos + window)
        )

    # CS gene counts in window — interpolate from precomputed density array
    # The density array is at 1kb resolution; use nearest-neighbour lookup
    gene_counts = np.interp(cs_pos_query, cs_gene_density_pos,
                            cs_gene_density_vals.astype(float))
    gene_counts = np.round(gene_counts).astype(float)

    # Compute fraction
    fractions = np.full(len(cs_pos_query), np.nan)
    has_genes = gene_counts > 0
    fractions[has_genes] = np.clip(
        anchor_counts[has_genes] / gene_counts[has_genes], 0.0, 1.0
    )
    return fractions


def build_spline_converter(cs_anchors, tgt_anchors, cs_chr_len, tgt_chr_len):
    """
    Build a PCHIP spline interpolator CS_pos -> tgt_pos from anchor arrays.

    Returns
    -------
    interpolator : callable
    mean_anchor_gap : float  (mean inter-anchor distance, bp)
    interp_type : str
    """
    sort_idx  = np.argsort(cs_anchors)
    cs_sorted = cs_anchors[sort_idx]
    tgt_sorted = tgt_anchors[sort_idx]

    if len(cs_sorted) >= 2:
        # Use robust slope from first/last min(10, N//4) anchors via linear regression
        # to avoid a single rogue anchor pair destabilising the boundary extrapolation
        n_edge = max(2, min(10, len(cs_sorted) // 4))

        # Left boundary slope
        cs_l, tgt_l = cs_sorted[:n_edge], tgt_sorted[:n_edge]
        cs_span_l = cs_l[-1] - cs_l[0]
        if cs_span_l > 0:
            slope_left = np.polyfit(cs_l, tgt_l, 1)[0]
        else:
            slope_left = (tgt_sorted[1] - tgt_sorted[0]) / max(cs_sorted[1] - cs_sorted[0], 1)

        # Right boundary slope
        cs_r, tgt_r = cs_sorted[-n_edge:], tgt_sorted[-n_edge:]
        cs_span_r = cs_r[-1] - cs_r[0]
        if cs_span_r > 0:
            slope_right = np.polyfit(cs_r, tgt_r, 1)[0]
        else:
            slope_right = (tgt_sorted[-1] - tgt_sorted[-2]) / max(cs_sorted[-1] - cs_sorted[-2], 1)

        tgt_left  = tgt_sorted[0]  - slope_left  * cs_sorted[0]
        tgt_right = tgt_sorted[-1] + slope_right * (cs_chr_len - cs_sorted[-1])

        # Clamp to valid range — never extrapolate beyond chromosome bounds
        tgt_left  = max(0,           tgt_left)
        tgt_right = min(tgt_chr_len, tgt_right)
    else:
        tgt_left  = 0
        tgt_right = tgt_chr_len

    cs_all  = np.concatenate([[0], cs_sorted, [cs_chr_len]])
    tgt_all = np.concatenate([[tgt_left], tgt_sorted, [tgt_right]])
    tgt_all = np.clip(tgt_all, 0, tgt_chr_len)

    # Remove duplicate x values (PCHIP requires strictly increasing x)
    cs_all, unique_idx = np.unique(cs_all, return_index=True)
    tgt_all = tgt_all[unique_idx]

    if np.all(np.diff(tgt_all) >= 0):
        interpolator = PchipInterpolator(cs_all, tgt_all, extrapolate=True)
        interp_type  = "pchip"
    elif np.all(np.diff(tgt_all) <= 0):
        interpolator = PchipInterpolator(cs_all, tgt_all[::-1], extrapolate=True)
        interp_type  = "pchip_inverted"
    else:
        interpolator = lambda x: np.interp(x, cs_all, tgt_all)
        interp_type  = "linear"

    mean_anchor_gap = np.mean(np.diff(cs_sorted)) if len(cs_sorted) > 1 else cs_chr_len
    return interpolator, mean_anchor_gap, interp_type


def compute_local_anchor_density(cs_pos_query, cs_anchors, window=5_000_000):
    """
    Count anchors within ±window bp of each query position.
    Retained as secondary metric.
    """
    densities = np.zeros(len(cs_pos_query), dtype=int)
    for i, pos in enumerate(cs_pos_query):
        densities[i] = np.sum(
            (cs_anchors >= pos - window) & (cs_anchors <= pos + window)
        )
    return densities


def compute_collinearity_score(cs_pos_query, cs_anchors, tgt_anchors,
                                interpolator, window=5_000_000,
                                residual_threshold=2_000_000):
    """
    Compute local collinearity score for each query position.

    For each query position, within ±window bp:
      1. Find all anchors in that window
      2. Compute each anchor's residual: |true_tgt_pos - spline_predicted_tgt_pos|
      3. Score = fraction of anchors with residual < residual_threshold

    Returns
    -------
    scores : np.ndarray of float (0.0–1.0), same length as cs_pos_query

    Why this works for introgressions:
      In a normal collinear region, local anchors follow the spline closely
      → small residuals → score near 1.0
      In an introgressed region (e.g. Lancer 2B rye segment), anchors that
      mapped to the wrong chromosome were already filtered in step 01, but
      anchors within the introgression that mapped to the correct chromosome
      but at shifted positions will have large residuals from the spline
      fitted to the flanking collinear anchors → score near 0.0
      In sparse pericentromeric regions, few or no anchors exist → score
      is np.nan, flagged separately as low-density.
    """
    predicted_tgt = np.array(interpolator(cs_anchors), dtype=float)
    residuals     = np.abs(tgt_anchors - predicted_tgt)

    scores = np.full(len(cs_pos_query), np.nan)
    for i, pos in enumerate(cs_pos_query):
        in_window = (cs_anchors >= pos - window) & (cs_anchors <= pos + window)
        n = int(np.sum(in_window))
        if n == 0:
            scores[i] = np.nan   # no anchors nearby — sparse region
        else:
            scores[i] = float(np.sum(residuals[in_window] < residual_threshold)) / n
    return scores


def build_synteny_blocks(anchors_df, min_block_genes=3, max_gap=2_000_000):
    if len(anchors_df) == 0:
        return pd.DataFrame()

    df = anchors_df.sort_values("cs_midpoint").reset_index(drop=True)
    blocks = []
    block_start = 0

    for i in range(1, len(df)):
        gap = df.loc[i, "cs_midpoint"] - df.loc[i-1, "cs_midpoint"]
        if gap > max_gap or i == len(df) - 1:
            block_end  = i if gap <= max_gap else i - 1
            block_genes = df.loc[block_start:block_end]
            if len(block_genes) >= min_block_genes:
                blocks.append({
                    "cs_chr":        block_genes["cs_chr"].iloc[0],
                    "cs_start":      int(block_genes["cs_midpoint"].min()),
                    "cs_end":        int(block_genes["cs_midpoint"].max()),
                    "tgt_chr":       block_genes["tgt_chr"].mode()[0],
                    "tgt_start":     int(block_genes["tgt_midpoint"].min()),
                    "tgt_end":       int(block_genes["tgt_midpoint"].max()),
                    "n_anchors":     len(block_genes),
                    "mean_identity": block_genes["identity"].mean(),
                })
            block_start = i

    return pd.DataFrame(blocks)


def plot_synteny_dots(anchors_df, cs_chr, target_name, output_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        import matplotlib.colors as mcolors

        fig, ax = plt.subplots(figsize=(8, 6))
        cmap = plt.cm.plasma_r
        norm = mcolors.Normalize(vmin=0.9, vmax=1.0)

        id_col = next((c for c in anchors_df.columns if c.endswith("_identity")), None)

        for tgt_chr, grp in anchors_df.groupby("tgt_chr"):
            x = grp["cs_midpoint"] / 1e6
            y = grp["tgt_midpoint"] / 1e6
            if id_col and id_col in grp.columns:
                c = grp[id_col].clip(0.9, 1.0)
                sc = ax.scatter(x, y, c=c, cmap=cmap, norm=norm,
                                s=2, alpha=0.6, linewidths=0, rasterized=True)
            else:
                ax.scatter(x, y, s=2, alpha=0.5, color="#0d0887",
                           linewidths=0, rasterized=True, label=tgt_chr)

        if id_col:
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.04)
            cbar.set_label("Sequence identity", fontsize=9)
            cbar.ax.tick_params(labelsize=8)

        ax.set_xlabel(f"CS {cs_chr} position (Mb)", fontsize=10)
        ax.set_ylabel(f"{target_name} position (Mb)", fontsize=10)
        ax.set_title(f"CS {cs_chr} → {target_name}", fontsize=11)
        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close()
    except ImportError:
        pass


def main():
    parser = argparse.ArgumentParser(
        description="Build coordinate conversion tables from anchor pairs"
    )
    parser.add_argument("--anchors",          required=True)
    parser.add_argument("--cs-chrom-sizes",   required=True)
    parser.add_argument("--output-dir",       required=True)
    parser.add_argument("--target-name",      required=True)
    parser.add_argument("--resolution",       type=int,   default=1000)
    parser.add_argument("--min-anchors-per-chr", type=int, default=20)
    parser.add_argument("--cs-gene-density", default=None,
                        help="Precomputed CS gene density TSV.gz from 00_compute_cs_gene_density.py")
    parser.add_argument("--density-window", type=int, default=5_000_000,
                        help="Window (bp) for anchor_density and anchor_fraction (default: 5000000)")
    parser.add_argument("--plot-dir", default=None)
    parser.add_argument("--overview-plot", default=None,
                        help="Path to save 21-chr overview dotplot PNG")
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.plot_dir:
        Path(args.plot_dir).mkdir(parents=True, exist_ok=True)

    print(f"Loading anchors: {args.anchors}", file=sys.stderr)
    anchors = pd.read_csv(args.anchors, sep="\t")
    print(f"  Total anchors: {len(anchors)}", file=sys.stderr)

    # Normalise assembly-prefixed column names to generic names
    # 01_extract_anchors.py writes e.g. Jagger_tgt_chr; strip prefix for internal use
    prefix = args.target_name + "_"
    rename_map = {c: c[len(prefix):] for c in anchors.columns if c.startswith(prefix)}
    if rename_map:
        anchors = anchors.rename(columns=rename_map)
    required = ["tgt_chr", "tgt_start", "tgt_end", "tgt_midpoint"]
    missing = [c for c in required if c not in anchors.columns]
    if missing:
        print(f"ERROR: missing columns {missing}. Got: {list(anchors.columns)}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading CS chrom sizes: {args.cs_chrom_sizes}", file=sys.stderr)
    cs_sizes = load_chrom_sizes(args.cs_chrom_sizes)

    if args.cs_gene_density:
        print(f"Loading CS gene density: {args.cs_gene_density}", file=sys.stderr)
        gene_density_df = pd.read_csv(
            args.cs_gene_density, sep="\t",
            compression="gzip" if args.cs_gene_density.endswith(".gz") else "infer"
        )
        print(f"  Gene density loaded: {len(gene_density_df):,} rows, "
              f"{gene_density_df['cs_chr'].nunique()} chromosomes", file=sys.stderr)
    else:
        gene_density_df = None
        print("  WARNING: --cs-gene-density not provided; anchor_fraction will be NaN",
              file=sys.stderr)

    chr_summary = []

    # Preload gene density for this run if available
    if gene_density_df is not None:
        gd_by_chr = {
            cs_chr: (sub["cs_pos"].values, sub["cs_gene_density"].values)
            for cs_chr, sub in gene_density_df.groupby("cs_chr")
        }
    else:
        gd_by_chr = {}

    for cs_chr, grp in anchors.groupby("cs_chr"):
        print(f"\nProcessing {cs_chr}...", file=sys.stderr)

        if len(grp) < args.min_anchors_per_chr:
            print(f"  SKIP: only {len(grp)} anchors (min={args.min_anchors_per_chr})",
                  file=sys.stderr)
            continue

        cs_chr_len = cs_sizes.get(
            cs_chr,
            cs_sizes.get(f"Chr{cs_chr}", cs_sizes.get(f"chr{cs_chr}", 0))
        )
        if cs_chr_len == 0:
            print(f"  WARN: could not find size for {cs_chr}", file=sys.stderr)
            cs_chr_len = int(grp["cs_end"].max() * 1.05)

        # Gene density arrays for this cs_chr
        gd_pos, gd_vals = gd_by_chr.get(cs_chr, (None, None))

        # Dense query positions for the full chromosome
        cs_positions = np.arange(0, cs_chr_len + args.resolution, args.resolution)

        # ── Detect segments ───────────────────────────────────────────────────
        # For standard chromosomes: one segment (collinear, one tgt_chr).
        # For translocated chromosomes: multiple segments, each covering a
        # different CS coordinate range and mapping to a different tgt_chr.
        segments = grp.groupby(["tgt_chr", "segment"], sort=False)
        seg_list = [(key, sub.sort_values("cs_midpoint").reset_index(drop=True))
                    for key, sub in segments
                    if len(sub) >= args.min_anchors_per_chr]

        if len(seg_list) == 0:
            print(f"  SKIP: no segments with >= {args.min_anchors_per_chr} anchors",
                  file=sys.stderr)
            continue

        # Sort segments by median cs_midpoint so boundaries are ordered
        seg_list.sort(key=lambda x: x[1]["cs_midpoint"].median())

        # Compute assignment boundary between adjacent segments:
        # midpoint between max cs_pos of seg[i] and min cs_pos of seg[i+1]
        seg_boundaries = []
        for i, (_, seg_df) in enumerate(seg_list):
            seg_start = 0 if i == 0 else seg_boundaries[-1][1] + 1
            if i < len(seg_list) - 1:
                next_df  = seg_list[i + 1][1]
                boundary = int((seg_df["cs_midpoint"].max() +
                                next_df["cs_midpoint"].min()) // 2)
                seg_end  = boundary
            else:
                seg_end  = cs_chr_len
            seg_boundaries.append((seg_start, seg_end))

        # All anchors on this cs_chr pooled (for density — avoids boundary artefacts)
        all_seg_anchors = np.concatenate(
            [seg_df["cs_midpoint"].values for _, seg_df in seg_list]
        )

        # Arrays to fill across all query positions
        tgt_chr_arr     = np.full(len(cs_positions), "", dtype=object)
        tgt_pos_arr     = np.zeros(len(cs_positions), dtype=np.int64)
        density_arr     = np.zeros(len(cs_positions), dtype=np.int32)
        fraction_arr    = np.full(len(cs_positions), np.nan)
        interp_type_arr = np.full(len(cs_positions), "", dtype=object)

        summary_segs = []
        for i, ((tgt_chr, seg_label), seg_df) in enumerate(seg_list):
            seg_start, seg_end = seg_boundaries[i]
            mask = (cs_positions >= seg_start) & (cs_positions <= seg_end)
            seg_positions = cs_positions[mask]

            tgt_chr_len = int(seg_df["tgt_end"].max() * 1.1)

            interpolator, mean_gap, interp_type = build_spline_converter(
                cs_anchors=seg_df["cs_midpoint"].values,
                tgt_anchors=seg_df["tgt_midpoint"].values,
                cs_chr_len=seg_end - seg_start,
                tgt_chr_len=tgt_chr_len,
            )
            print(f"  {cs_chr} [{seg_label}] -> {tgt_chr}: {len(seg_df)} anchors, "
                  f"mean gap {mean_gap/1e6:.2f} Mb, interp: {interp_type}",
                  file=sys.stderr)

            tgt_positions = np.clip(
                np.round(interpolator(seg_positions)).astype(np.int64),
                0, tgt_chr_len
            )

            # Use all_seg_anchors for density so transitions don't cause drops
            seg_density = compute_local_anchor_density(
                seg_positions, all_seg_anchors,
                window=args.density_window,
            )

            if gd_pos is not None:
                seg_fraction = compute_anchor_fraction(
                    seg_positions, all_seg_anchors,
                    gd_pos, gd_vals,
                    window=args.density_window,
                )
            else:
                seg_fraction = np.full(len(seg_positions), np.nan)

            tgt_chr_arr[mask]     = tgt_chr
            tgt_pos_arr[mask]     = tgt_positions
            density_arr[mask]     = seg_density
            fraction_arr[mask]    = seg_fraction
            interp_type_arr[mask] = interp_type

            summary_segs.append({
                "tgt_chr":    tgt_chr,
                "segment":    seg_label,
                "n_anchors":  len(seg_df),
                "mean_gap":   mean_gap,
                "interp_type": interp_type,
            })

        # Log anchor_fraction QC
        valid_af = fraction_arr[~np.isnan(fraction_arr)]
        if len(valid_af) > 0:
            print(f"  Anchor fraction: min={valid_af.min():.2f}, "
                  f"mean={valid_af.mean():.2f}, max={valid_af.max():.2f}",
                  file=sys.stderr)
            n_low = int(np.sum(valid_af < 0.5))
            if n_low > 0:
                frac_low = n_low / len(valid_af)
                print(f"  WARNING: {n_low} positions ({frac_low:.1%}) have "
                      f"anchor_fraction < 0.5", file=sys.stderr)

        # Build output DataFrame — prefix all target columns with assembly name
        pfx = args.target_name + "_"
        conv_df = pd.DataFrame({
            "cs_chr":                      cs_chr,
            "cs_pos":                      cs_positions,
            f"{pfx}tgt_chr":              tgt_chr_arr,
            f"{pfx}tgt_pos":              tgt_pos_arr,
            f"{pfx}anchor_density_5Mb":   density_arr,
            f"{pfx}anchor_fraction":       np.round(fraction_arr, 3),
            f"{pfx}interp_type":          interp_type_arr,
        })

        # Save conversion table
        conv_file = outdir / f"{cs_chr}_conversion.tsv.gz"
        conv_df.to_csv(conv_file, sep="\t", index=False, compression="gzip")
        print(f"  Saved: {conv_file}", file=sys.stderr)

        # Save clean anchors — re-prefix tgt columns with assembly name
        grp_all_segs = pd.concat([seg_df for _, seg_df in seg_list], ignore_index=True)
        pfx = args.target_name + "_"
        tgt_generic = ["tgt_chr", "tgt_start", "tgt_end", "tgt_strand",
                       "tgt_midpoint", "coverage", "identity"]
        anchor_rename = {c: f"{pfx}{c}" for c in tgt_generic if c in grp_all_segs.columns}
        grp_out = grp_all_segs.rename(columns=anchor_rename)
        anchor_file = outdir / f"{cs_chr}_anchors.tsv"
        grp_out.to_csv(anchor_file, sep="\t", index=False)

        # Synteny blocks (per segment, combined)
        all_blocks = []
        for _, seg_df in seg_list:
            blocks = build_synteny_blocks(seg_df)
            if len(blocks) > 0:
                all_blocks.append(blocks)
        if all_blocks:
            blocks_df = pd.concat(all_blocks, ignore_index=True)
            # Prefix tgt columns in synteny blocks
            blk_rename = {c: f"{args.target_name}_{c}"
                          for c in ["tgt_chr","tgt_start","tgt_end"]
                          if c in blocks_df.columns}
            blocks_df = blocks_df.rename(columns=blk_rename)
            blocks_df.to_csv(outdir / f"{cs_chr}_synteny_blocks.tsv", sep="\t", index=False)
            bed_file = outdir / f"{cs_chr}_synteny.bed"
            pfx = args.target_name + "_"
            tc, ts, te = f"{pfx}tgt_chr", f"{pfx}tgt_start", f"{pfx}tgt_end"
            with open(bed_file, "w") as bf:
                bf.write(f"# CS_chr\tCS_start\tCS_end\t{tc}\t{ts}\t{te}\tn_anchors\tmean_identity\n")
                for _, b in blocks_df.iterrows():
                    bf.write(f"{b['cs_chr']}\t{b['cs_start']}\t{b['cs_end']}\t"
                             f"{b[tc]}\t{b[ts]}\t{b[te]}\t"
                             f"{b['n_anchors']}\t{b['mean_identity']:.4f}\n")

        # QC plot
        if args.plot_dir:
            plot_synteny_dots(grp, cs_chr, args.target_name,
                              Path(args.plot_dir) / f"{cs_chr}_dotplot.png")

        # Summary: one row per cs_chr; join segment info if multiple
        tgt_chr_summary = "+".join(s["tgt_chr"] for s in summary_segs)
        total_anchors   = sum(s["n_anchors"] for s in summary_segs)
        mean_gap_all    = float(np.mean([s["mean_gap"] for s in summary_segs]))
        pfx = args.target_name + "_"
        chr_summary.append({
            "cs_chr":               cs_chr,
            f"{pfx}tgt_chr":        tgt_chr_summary,
            "n_anchors":            total_anchors,
            f"{pfx}mean_anchor_gap_Mb":   round(mean_gap_all / 1e6, 3),
            f"{pfx}interp_type":          "+".join(s["interp_type"] for s in summary_segs),
            "n_segments":                 len(summary_segs),
            f"{pfx}mean_anchor_fraction": round(float(np.nanmean(fraction_arr)), 3)
                                          if len(fraction_arr) > 0 else None,
        })

    summary_df = pd.DataFrame(chr_summary)
    summary_df.to_csv(outdir / "conversion_summary.tsv", sep="\t", index=False)

    print(f"\n{'='*60}", file=sys.stderr)
    print(f"SUMMARY for {args.target_name}:", file=sys.stderr)
    print(summary_df.to_string(index=False), file=sys.stderr)
    print(f"\nAll outputs in: {outdir}", file=sys.stderr)

    # ── 21-chromosome overview dotplot ────────────────────────────────────────
    if args.overview_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors

            CHROMOSOMES = [f"Chr{g}{s}" for g in range(1, 8) for s in ["A", "B", "D"]]
            NROWS, NCOLS = 7, 3
            cmap = plt.cm.plasma_r
            norm = mcolors.Normalize(vmin=0.9, vmax=1.0)

            fig, axes = plt.subplots(NROWS, NCOLS, figsize=(14, 22), squeeze=False)
            fig.suptitle(
                f"CS RefSeq v2.1  →  {args.target_name}\nSynteny anchors — all 21 chromosomes",
                fontsize=13, fontweight="bold", y=0.998
            )

            has_identity = False
            pfx = args.target_name + "_"

            for idx, chrom in enumerate(CHROMOSOMES):
                row, col = divmod(idx, NCOLS)
                ax = axes[row][col]
                ax.set_title(f"{chrom[3]}{chrom[4]}", fontsize=9, pad=3)

                anchor_file = outdir / f"{chrom}_anchors.tsv"
                if not anchor_file.exists():
                    ax.text(0.5, 0.5, "no data", ha="center", va="center",
                            transform=ax.transAxes, fontsize=8, color="grey")
                    ax.set_xticks([]); ax.set_yticks([])
                    continue

                df_a = pd.read_csv(anchor_file, sep="\t")
                if df_a.empty:
                    ax.text(0.5, 0.5, "no data", ha="center", va="center",
                            transform=ax.transAxes, fontsize=8, color="grey")
                    ax.set_xticks([]); ax.set_yticks([])
                    continue

                tgt_mid_col = f"{pfx}tgt_midpoint"
                id_col      = f"{pfx}identity"

                if "cs_midpoint" not in df_a.columns or tgt_mid_col not in df_a.columns:
                    ax.text(0.5, 0.5, "column\nerror", ha="center", va="center",
                            transform=ax.transAxes, fontsize=7, color="red")
                    continue

                x = df_a["cs_midpoint"] / 1e6
                y = df_a[tgt_mid_col] / 1e6

                if id_col in df_a.columns:
                    c = df_a[id_col].clip(0.9, 1.0)
                    ax.scatter(x, y, c=c, cmap=cmap, norm=norm,
                               s=1.5, alpha=0.6, linewidths=0, rasterized=True)
                    has_identity = True
                else:
                    ax.scatter(x, y, s=1.5, alpha=0.6, color="#3a86ff",
                               linewidths=0, rasterized=True)

                if row == NROWS - 1:
                    ax.set_xlabel("CS (Mb)", fontsize=7)
                else:
                    ax.tick_params(labelbottom=False)
                if col == 0:
                    ax.set_ylabel(f"{args.target_name} (Mb)", fontsize=7)
                else:
                    ax.tick_params(labelleft=False)
                ax.tick_params(axis="both", labelsize=6)
                ax.text(0.97, 0.03, f"n={len(df_a):,}", ha="right", va="bottom",
                        transform=ax.transAxes, fontsize=6, color="dimgrey")

            if has_identity:
                cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.25])
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                cbar = fig.colorbar(sm, cax=cbar_ax)
                cbar.set_label("Sequence identity", fontsize=8)
                cbar.ax.tick_params(labelsize=7)

            plt.subplots_adjust(left=0.07, right=0.90, top=0.97, bottom=0.04,
                                hspace=0.35, wspace=0.15)

            Path(args.overview_plot).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(args.overview_plot, dpi=200, bbox_inches="tight")
            plt.close(fig)
            print(f"  Overview dotplot saved: {args.overview_plot}", file=sys.stderr)

        except Exception as e:
            print(f"  WARNING: overview dotplot failed: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
