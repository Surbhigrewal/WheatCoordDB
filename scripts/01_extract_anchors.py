#!/usr/bin/env python3
"""
01_extract_anchors.py
=====================
Extract orthologous gene anchor pairs between CS v2.1 and a target assembly
from Liftoff output.

For each successfully lifted-off gene, records:
  - CS chromosome, start, end, strand
  - Target chromosome, start, end, strand
  - Coverage, identity scores from Liftoff

Supports known inter-chromosomal translocations (e.g. 5B/7B Robertsonian in
ArinaLrFor and SY_Mattis) via --translocation-map. For translocated segments,
anchors are assigned to the correct target chromosome and processed with
separate collinearity filtering per segment.

Output: TSV with columns:
  gene_id | cs_chr | cs_start | cs_end | cs_strand | cs_midpoint |
  <assembly>_tgt_chr | <assembly>_tgt_start | <assembly>_tgt_end |
  <assembly>_tgt_strand | <assembly>_tgt_midpoint |
  <assembly>_coverage | <assembly>_identity | segment
  (segment: 'collinear' | 'translocation:<cs_chr>_to_<tgt_chr>')

Target columns are prefixed with the assembly name at output time so that
anchor files are self-describing when used outside the pipeline context.
The segment column is not prefixed as it describes the anchor relationship,
not an assembly-specific attribute.
"""

import argparse
import re
import gzip
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd



def _open(path):
    """Open plain or gzipped text file."""
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

def parse_gff3_attributes(attr_string):
    """Parse GFF3 attribute column into a dict."""
    attrs = {}
    for field in attr_string.strip().split(";"):
        field = field.strip()
        if "=" in field:
            key, val = field.split("=", 1)
            attrs[key.strip()] = val.strip()
    return attrs


def parse_gff3(gff_path, feature_type="gene", chromosomes=None):
    """
    Parse a GFF3 file and return a dict: gene_id -> record dict.
    
    Parameters
    ----------
    gff_path : str | Path
    feature_type : str
        Feature type to extract (default: 'gene')
    chromosomes : set | None
        If provided, only parse features on these chromosomes
    
    Returns
    -------
    dict : gene_id -> {chr, start, end, strand, ...}
    """
    records = {}
    skipped = 0
    
    with _open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_str = parts
            
            if ftype != feature_type:
                continue
            if chromosomes and seqid not in chromosomes:
                skipped += 1
                continue
            
            attrs = parse_gff3_attributes(attrs_str)
            gene_id = attrs.get("ID", attrs.get("gene_id", None))
            if not gene_id:
                continue
            
            # Strip version suffixes if present (e.g. TraesCS1A02G000100.1 -> TraesCS1A02G000100)
            gene_id_base = gene_id.split(".")[0] if "." in gene_id else gene_id
            
            records[gene_id_base] = {
                "chr": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "raw_id": gene_id,
            }
    
    return records


def parse_liftoff_gff3(liftoff_gff, feature_type="gene", min_coverage=0.9, min_identity=0.9):
    """
    Parse Liftoff output GFF3.
    Liftoff adds extra attributes: coverage, sequence_ID, extra_copy_number, 
    valid_ORFs, intervening_sequence.
    
    Returns dict: original_gene_id -> {chr, start, end, strand, coverage, identity}
    Only returns entries passing coverage and identity thresholds.
    """
    records = {}
    failed = 0
    
    with _open(liftoff_gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_str = parts
            
            if ftype != feature_type:
                continue
            
            attrs = parse_gff3_attributes(attrs_str)
            
            # Liftoff stores original ID in 'ID', copies are 'ID_1', 'ID_2', etc.
            gene_id = attrs.get("ID", "")
            
            # Skip copy duplicates (keep only primary mapping)
            if re.search(r"_\d+$", gene_id):
                continue
            
            gene_id_base = gene_id.split(".")[0] if "." in gene_id else gene_id
            
            # Parse Liftoff-specific quality metrics
            coverage = float(attrs.get("coverage", 0.0))
            identity = float(attrs.get("sequence_ID", attrs.get("identity", 0.0)))
            
            if coverage < min_coverage or identity < min_identity:
                failed += 1
                continue
            
            records[gene_id_base] = {
                "chr": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "coverage": coverage,
                "identity": identity,
            }
    
    print(f"  Liftoff parse: {len(records)} passed, {failed} failed QC thresholds",
          file=sys.stderr)
    return records


def detect_cs_chromosomes(cs_records):
    """
    Detect canonical CS chromosome names (1A, 1B, ..., 7D).
    Returns set of chromosome names present in the CS GFF.
    """
    chrs = set(r["chr"] for r in cs_records.values())
    return chrs


def assign_homoeologous_group(chr_name):
    """
    Extract chromosome group number from CS chromosome name.
    Handles both Ensembl style (1A) and chr-prefixed style (chr1A).
    """
    match = re.search(r"(\d[A-Da-d])", chr_name)
    return match.group(1).upper() if match else None


def load_translocation_map(tsv_path, assembly_name):
    """
    Load translocation map TSV for a specific assembly.

    TSV format (tab-separated, # comment lines allowed):
      assembly  cs_chr  cs_bp_proximal  cs_bp_distal  tgt_chr_translocated

    cs_bp_proximal / cs_bp_distal define the CS breakpoint interval (bp).
    Genes with cs_midpoint < cs_bp_proximal map collinearly to the CS chr.
    Genes with cs_midpoint > cs_bp_distal map to tgt_chr_translocated.
    Genes between the two breakpoints are in the transition zone and are
    assigned to whichever side their actual tgt_chr matches.

    Example row:
      SY_Mattis  Chr5B  166000000  222000000  Chr7B

    This means:
      CS Chr5B 1–166 Mb   → SY_Mattis Chr5B  (5BS, collinear)
      CS Chr5B 222 Mb–end → SY_Mattis Chr7B  (5BL translocated onto Chr7B)

    Returns list of dicts, one per translocation segment for this assembly.
    Returns [] if tsv_path is None or assembly not present.
    """
    if tsv_path is None:
        return []

    translocations = []
    with open(tsv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                print(f"  WARN: skipping malformed translocation line: {line}",
                      file=sys.stderr)
                continue
            asm, cs_chr, bp_prox, bp_dist, tgt_chr_trans = parts[:5]
            if asm != assembly_name:
                continue
            translocations.append({
                "cs_chr": cs_chr,
                "bp_proximal": int(bp_prox),   # CS pos: proximal breakpoint
                "bp_distal":   int(bp_dist),   # CS pos: distal breakpoint
                "tgt_chr_trans": tgt_chr_trans, # target chr for distal segment
            })

    if translocations:
        print(f"  Loaded {len(translocations)} translocation segment(s) for {assembly_name}:",
              file=sys.stderr)
        for t in translocations:
            print(f"    {t['cs_chr']} bp {t['bp_proximal']:,}–{t['bp_distal']:,} "
                  f"→ distal segment on {t['tgt_chr_trans']}", file=sys.stderr)

    return translocations


def build_anchor_table(cs_records, liftoff_records):
    """
    Inner join CS and Liftoff records on gene_id_base.
    Returns DataFrame of anchor pairs.
    """
    rows = []
    
    for gene_id, cs_rec in cs_records.items():
        if gene_id not in liftoff_records:
            continue
        tgt_rec = liftoff_records[gene_id]
        
        cs_mid = (cs_rec["start"] + cs_rec["end"]) // 2
        tgt_mid = (tgt_rec["start"] + tgt_rec["end"]) // 2
        
        rows.append({
            "gene_id": gene_id,
            "cs_chr": cs_rec["chr"],
            "cs_start": cs_rec["start"],
            "cs_end": cs_rec["end"],
            "cs_strand": cs_rec["strand"],
            "cs_midpoint": cs_mid,
            "tgt_chr": tgt_rec["chr"],
            "tgt_start": tgt_rec["start"],
            "tgt_end": tgt_rec["end"],
            "tgt_strand": tgt_rec["strand"],
            "tgt_midpoint": tgt_mid,
            "coverage": tgt_rec["coverage"],
            "identity": tgt_rec["identity"],
            "segment": "collinear",  # may be overwritten by translocation filter
        })
    
    df = pd.DataFrame(rows)
    return df


def _filter_one_segment(grp, label, max_outlier_distance=5e6):
    """
    Apply collinearity filtering (linear fit + MAD residual removal) to a
    single syntenic segment (one cs_chr → one tgt_chr block).

    Returns cleaned DataFrame with 'segment' column set to label.
    """
    if len(grp) < 5:
        return grp

    grp = grp.sort_values("cs_midpoint").reset_index(drop=True)

    cs_pos  = grp["cs_midpoint"].values
    tgt_pos = grp["tgt_midpoint"].values

    from numpy.polynomial import polynomial as P
    coeffs    = P.polyfit(cs_pos, tgt_pos, deg=1)
    predicted = P.polyval(cs_pos, coeffs)
    residuals = np.abs(tgt_pos - predicted)

    med_res   = np.median(residuals)
    mad       = np.median(np.abs(residuals - med_res))
    threshold = max(med_res + 5 * mad, max_outlier_distance)

    grp = grp[residuals <= threshold].copy()
    grp["segment"] = label

    # Warn if majority strand-flipped (whole-chr orientation difference)
    if len(grp) > 0:
        same_strand = (grp["cs_strand"] == grp["tgt_strand"]).mean()
        if same_strand < 0.5:
            cs_c  = grp["cs_chr"].iloc[0]
            tgt_c = grp["tgt_chr"].iloc[0]
            print(f"  INFO: {cs_c} -> {tgt_c} [{label}]: "
                  f"majority strand-flipped ({1-same_strand:.0%}), "
                  f"likely assembly orientation difference", file=sys.stderr)

    return grp


def filter_collinear(df, translocations=None, max_outlier_distance=5e6):
    """
    Per chromosome pair, filter out non-collinear (translocated) genes.

    If translocations is provided (list of dicts from load_translocation_map),
    chromosomes with known translocations are split at the CS breakpoints and
    each segment is filtered independently against its correct target chr.

    For translocated chromosomes the logic is:
      proximal segment: cs_midpoint <= bp_proximal  → dominant collinear tgt_chr
      distal   segment: cs_midpoint >= bp_distal    → tgt_chr_trans
      breakpoint zone:  bp_proximal < cs_mid < bp_distal → assigned by actual
                        tgt_chr match (whichever side it matches)

    Standard chromosomes (no translocation entry) use the original dominant-
    target-chr approach.
    """
    if translocations is None:
        translocations = []

    # Build lookup: cs_chr → translocation spec
    trans_lookup = {t["cs_chr"]: t for t in translocations}

    clean_rows = []

    for cs_chr, grp in df.groupby("cs_chr"):
        if len(grp) == 0:
            continue

        # ── Translocation-aware path ──────────────────────────────────────────
        if cs_chr in trans_lookup:
            t = trans_lookup[cs_chr]
            bp_prox     = t["bp_proximal"]
            bp_dist     = t["bp_distal"]
            tgt_trans   = t["tgt_chr_trans"]

            # Proximal segment: cs_mid <= bp_proximal
            prox = grp[grp["cs_midpoint"] <= bp_prox].copy()
            # Distal  segment: cs_mid >= bp_distal  (should map to tgt_trans)
            dist = grp[grp["cs_midpoint"] >= bp_dist].copy()
            # Breakpoint zone: between the two breakpoints
            zone = grp[(grp["cs_midpoint"] > bp_prox) &
                       (grp["cs_midpoint"] < bp_dist)].copy()

            # For proximal: find dominant collinear tgt_chr
            if len(prox) > 0:
                dom_prox = prox["tgt_chr"].value_counts().index[0]
                prox = prox[prox["tgt_chr"] == dom_prox]
                prox = _filter_one_segment(prox, "collinear", max_outlier_distance)
                if len(prox) > 0:
                    clean_rows.append(prox)
                    print(f"  {cs_chr} proximal → {dom_prox}: {len(prox)} anchors",
                          file=sys.stderr)

            # For distal: keep only anchors that mapped to tgt_trans
            if len(dist) > 0:
                dist_on_trans = dist[dist["tgt_chr"] == tgt_trans].copy()
                if len(dist_on_trans) >= 5:
                    seg_label = f"translocation:{cs_chr}_to_{tgt_trans}"
                    dist_on_trans = _filter_one_segment(
                        dist_on_trans, seg_label, max_outlier_distance)
                    if len(dist_on_trans) > 0:
                        clean_rows.append(dist_on_trans)
                        print(f"  {cs_chr} distal (translocated) → {tgt_trans}: "
                              f"{len(dist_on_trans)} anchors", file=sys.stderr)
                else:
                    print(f"  WARN: {cs_chr} distal segment has only "
                          f"{len(dist_on_trans)} anchors on {tgt_trans} — skipping",
                          file=sys.stderr)

            # Breakpoint zone: assign by actual tgt_chr
            if len(zone) > 0:
                # Determine which tgt_chr each gene actually landed on
                dom_prox_chr = prox["tgt_chr"].iloc[0] if len(prox) > 0 else None
                for tgt_c, z_grp in zone.groupby("tgt_chr"):
                    if tgt_c == dom_prox_chr:
                        z_grp = _filter_one_segment(
                            z_grp.copy(), "collinear", max_outlier_distance)
                    elif tgt_c == tgt_trans:
                        seg_label = f"translocation:{cs_chr}_to_{tgt_trans}"
                        z_grp = _filter_one_segment(
                            z_grp.copy(), seg_label, max_outlier_distance)
                    else:
                        continue  # off-target mapping, skip
                    if len(z_grp) > 0:
                        clean_rows.append(z_grp)
            continue

        # ── Standard path (no translocation) ─────────────────────────────────
        tgt_chr_counts  = grp["tgt_chr"].value_counts()
        dominant_tgt_chr = tgt_chr_counts.index[0]
        dominant_fraction = tgt_chr_counts.iloc[0] / len(grp)

        if dominant_fraction < 0.7:
            print(f"  WARN: {cs_chr} has fragmented mapping "
                  f"(best tgt: {dominant_tgt_chr} = {dominant_fraction:.1%})",
                  file=sys.stderr)

        grp_clean = grp[grp["tgt_chr"] == dominant_tgt_chr].copy()
        if len(grp_clean) < 5:
            continue

        grp_clean = _filter_one_segment(grp_clean, "collinear", max_outlier_distance)
        if len(grp_clean) > 0:
            clean_rows.append(grp_clean)

    if not clean_rows:
        return pd.DataFrame()

    result = pd.concat(clean_rows, ignore_index=True)
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Extract ortholog anchor pairs from CS GFF3 and Liftoff output"
    )
    parser.add_argument("--cs-gff", required=True, help="CS v2.1 annotation GFF3")
    parser.add_argument("--liftoff-gff", required=True, help="Liftoff output GFF3")
    parser.add_argument("--output", required=True, help="Output TSV anchor file")
    parser.add_argument("--feature-type", default="gene", help="GFF3 feature type (default: gene)")
    parser.add_argument("--min-coverage", type=float, default=0.9)
    parser.add_argument("--min-identity", type=float, default=0.9)
    parser.add_argument("--no-collinearity-filter", action="store_true",
                        help="Skip collinearity outlier removal")
    parser.add_argument("--translocation-map", default=None,
                        help="TSV file of known translocations "
                             "(assembly, cs_chr, bp_proximal, bp_distal, tgt_chr_trans)")
    parser.add_argument("--assembly-name", required=True,
                        help="Name of target assembly (e.g. Jagger). Used to prefix "
                             "target columns in output: tgt_chr -> <assembly>_tgt_chr etc.")
    args = parser.parse_args()

    print(f"Parsing CS GFF3: {args.cs_gff}", file=sys.stderr)
    cs_records = parse_gff3(args.cs_gff, feature_type=args.feature_type)
    print(f"  CS records: {len(cs_records)} {args.feature_type} features", file=sys.stderr)

    print(f"Parsing Liftoff GFF3: {args.liftoff_gff}", file=sys.stderr)
    liftoff_records = parse_liftoff_gff3(
        args.liftoff_gff,
        feature_type=args.feature_type,
        min_coverage=args.min_coverage,
        min_identity=args.min_identity,
    )

    print("Building anchor table...", file=sys.stderr)
    df = build_anchor_table(cs_records, liftoff_records)
    print(f"  Raw anchor pairs: {len(df)}", file=sys.stderr)

    if not args.no_collinearity_filter and len(df) > 0:
        print("Filtering for collinearity...", file=sys.stderr)
        translocations = load_translocation_map(args.translocation_map, args.assembly_name or "")
        df = filter_collinear(df, translocations=translocations)
        print(f"  Collinear anchors: {len(df)}", file=sys.stderr)

    # Print per-chromosome summary
    if len(df) > 0:
        print("\nPer-chromosome anchor counts:", file=sys.stderr)
        for cs_chr, grp in df.groupby("cs_chr"):
            # Show all unique (tgt_chr, segment) combos for this cs_chr
            for (tgt, seg), sg in grp.groupby(["tgt_chr", "segment"]):
                tag = f" [{seg}]" if seg != "collinear" else ""
                print(f"  {cs_chr} -> {tgt}{tag}: {len(sg)} anchors", file=sys.stderr)

    # Prefix target columns with assembly name
    # cs_* columns and gene_id/segment stay unprefixed; tgt_* columns get assembly prefix
    tgt_cols = ["tgt_chr", "tgt_start", "tgt_end", "tgt_strand", "tgt_midpoint",
                "coverage", "identity"]
    rename_map = {c: f"{args.assembly_name}_{c}" for c in tgt_cols if c in df.columns}
    if rename_map:
        df = df.rename(columns=rename_map)
        print(f"  Prefixed columns: {list(rename_map.values())}", file=sys.stderr)

    # Reorder columns: gene_id, cs_*, <asm>_tgt_*, segment
    cs_cols  = ["gene_id", "cs_chr", "cs_start", "cs_end", "cs_strand", "cs_midpoint"]
    tgt_pfx  = [f"{args.assembly_name}_{c}" for c in tgt_cols if c in rename_map or
                f"{args.assembly_name}_{c}" in df.columns]
    other    = [c for c in df.columns if c not in cs_cols and c not in tgt_pfx]
    final_cols = [c for c in cs_cols + tgt_pfx + other if c in df.columns]
    df = df[final_cols]

    # Save output
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved anchor table: {args.output} ({len(df)} pairs)", file=sys.stderr)


if __name__ == "__main__":
    main()
