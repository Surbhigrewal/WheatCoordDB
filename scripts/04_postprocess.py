#!/usr/bin/env python3
"""
04_postprocess.py
=================
Run once after ALL liftoff array jobs have completed.

Rename, dotplots, and per-assembly column prefixing are handled inline
by 2_run_liftoff_array.sh. This script only needs to:

  1. MASTER ANCHORS  — merges all per-target anchor TSVs into one Excel file:
                       - 'All_targets_merged' sheet (CS columns once, all
                         target columns side by side)
                       - One sheet per assembly

  2. MASTER SUMMARY  — merges all conversion_summary.tsv files into one TSV
                       and one Excel sheet with an 'assembly' column prepended.

Usage (from CS_synteny/ directory):
    conda activate liftoff_env
    python3 04_postprocess.py

    # Or specify targets explicitly:
    python3 04_postprocess.py --targets Jagger Lancer ArinaLrFor
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

ALL_TARGETS = [
    "Jagger", "Lancer", "ArinaLrFor", "Stanley", "Spelt",
    "Mace", "SY_Mattis", "Julius", "Landmark", "Norin61",
    "Aikang58", "Chunmai104", "CS_IAAS", "CS_CAU", "Sumai3",
    "JIN50", "MOV", "CS_v1", "Fielder", "Kariega",
    "Attraktion", "Renan_v2", "Paragon_v3", "Cadenza_v2",
]

CHROMOSOMES = [f"Chr{g}{s}" for g in range(1, 8) for s in ["A", "B", "D"]]

CHR_ORDER = {chrom: i for i, chrom in enumerate(CHROMOSOMES)}

CS_ANCHOR_COLS = ["gene_id", "cs_chr", "cs_start", "cs_end",
                  "cs_strand", "cs_midpoint"]


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def sort_by_cs(df: pd.DataFrame) -> pd.DataFrame:
    """Sort a DataFrame by cs_chr (wheat order) then cs_start."""
    if "cs_chr" in df.columns:
        df = df.copy()
        df["_sort"] = df["cs_chr"].map(CHR_ORDER).fillna(999)
        sort_cols = ["_sort"]
        if "cs_start" in df.columns:
            sort_cols.append("cs_start")
        df = df.sort_values(sort_cols).drop(columns=["_sort"]).reset_index(drop=True)
    return df


def autofit_excel(writer, sheet_name: str, df: pd.DataFrame):
    ws = writer.sheets[sheet_name]
    for col_idx, col_name in enumerate(df.columns, start=1):
        sample = df[col_name].astype(str).head(200)
        max_len = max(len(str(col_name)),
                      int(sample.str.len().max()) if len(sample) else 0)
        ws.column_dimensions[
            ws.cell(row=1, column=col_idx).column_letter
        ].width = min(max_len + 2, 30)


# ─────────────────────────────────────────────────────────────────────────────
# Step 1 — Master anchor Excel
# ─────────────────────────────────────────────────────────────────────────────

def make_master_anchors(anchors_dir: Path, output_path: Path, assemblies: list):
    """
    Merge all <assembly>_anchors.tsv files into one Excel workbook.
    Sheets: 'All_targets_merged' + one per assembly.
    """
    per_target = {}
    for assembly in assemblies:
        p = anchors_dir / f"{assembly}_anchors.tsv"
        if p.exists():
            df = sort_by_cs(read_tsv(p))
            # Strip 'gene:' prefix from gene_id for consistent merging
            if "gene_id" in df.columns:
                df["gene_id"] = df["gene_id"].str.replace(r"^gene:", "", regex=True)
            per_target[assembly] = df
            print(f"  {assembly}: {len(df):,} anchors")
        else:
            print(f"  SKIP {assembly}: {p} not found", file=sys.stderr)

    if not per_target:
        print("  ERROR: no anchor files found", file=sys.stderr)
        return

    # Build CS coordinate base from union of ALL targets so no gene is missing
    cs_cols_present = [c for c in CS_ANCHOR_COLS
                       if c in next(iter(per_target.values())).columns]
    merge_key = "gene_id" if "gene_id" in cs_cols_present else ["cs_chr", "cs_start"]
    dedup_key = ["gene_id"] if isinstance(merge_key, str) else merge_key

    cs_base_frames = [df[cs_cols_present].drop_duplicates(subset=dedup_key)
                      for df in per_target.values()]
    merged = pd.concat(cs_base_frames, ignore_index=True) \
               .drop_duplicates(subset=dedup_key) \
               .reset_index(drop=True)
    print(f"  CS gene universe: {len(merged):,} unique genes across all targets")

    for assembly, df in per_target.items():
        # Keep only target columns (already prefixed) + merge key
        tgt_cols = [c for c in df.columns if c.startswith(f"{assembly}_")]
        keep = ([merge_key] if isinstance(merge_key, str) else merge_key) + tgt_cols
        keep = [c for c in keep if c in df.columns]
        merged = merged.merge(df[keep], on=merge_key, how="left")

    merged = sort_by_cs(merged)
    print(f"  Merged: {len(merged):,} rows × {len(merged.columns)} columns")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        # All targets merged sheet first
        merged.to_excel(writer, sheet_name="All_targets_merged", index=False)
        autofit_excel(writer, "All_targets_merged", merged)

        # Individual sheets
        for assembly, df in per_target.items():
            sheet = assembly[:31]
            df.to_excel(writer, sheet_name=sheet, index=False)
            autofit_excel(writer, sheet, df)

    print(f"  Saved: {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Step 2 — Master conversion summary
# ─────────────────────────────────────────────────────────────────────────────

def make_master_summary(conv_dir: Path, output_dir: Path, assemblies: list):
    """
    Concatenate all conversion_summary.tsv files into:
      - outputs/master_conversion_summary.tsv
      - outputs/master_conversion_summary.xlsx
    """
    frames = []
    for assembly in assemblies:
        p = conv_dir / assembly / "conversion_summary.tsv"
        if p.exists():
            df = read_tsv(p)
            df.insert(0, "assembly", assembly)
            frames.append(df)
            print(f"  {assembly}: {len(df)} chromosomes")
        else:
            print(f"  SKIP {assembly}: conversion_summary.tsv not found",
                  file=sys.stderr)

    if not frames:
        print("  ERROR: no summary files found", file=sys.stderr)
        return

    master = pd.concat(frames, ignore_index=True)

    # Sort by assembly order then chromosome order
    assembly_order = {a: i for i, a in enumerate(assemblies)}
    master["_asm_sort"] = master["assembly"].map(assembly_order).fillna(999)
    chr_col = "cs_chr" if "cs_chr" in master.columns else None
    if chr_col:
        master["_chr_sort"] = master[chr_col].map(CHR_ORDER).fillna(999)
        master = master.sort_values(["_asm_sort", "_chr_sort"]) \
                       .drop(columns=["_asm_sort", "_chr_sort"]) \
                       .reset_index(drop=True)
    else:
        master = master.sort_values("_asm_sort") \
                       .drop(columns=["_asm_sort"]) \
                       .reset_index(drop=True)

    output_dir.mkdir(parents=True, exist_ok=True)

    tsv_out = output_dir / "master_conversion_summary.tsv"
    master.to_csv(tsv_out, sep="\t", index=False)
    print(f"  Saved: {tsv_out}")

    xlsx_out = output_dir / "master_conversion_summary.xlsx"
    with pd.ExcelWriter(xlsx_out, engine="openpyxl") as writer:
        master.to_excel(writer, sheet_name="conversion_summary", index=False)
        autofit_excel(writer, "conversion_summary", master)
    print(f"  Saved: {xlsx_out}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate master anchor and summary files across all assemblies"
    )
    parser.add_argument(
        "--targets", nargs="+", default=ALL_TARGETS,
        help="Assembly names to process (default: all 24)"
    )
    parser.add_argument(
        "--outdir", default="outputs",
        help="Base output directory (default: outputs)"
    )
    args = parser.parse_args()

    outdir      = Path(args.outdir)
    conv_dir    = outdir / "conversion_tables"
    anchors_dir = outdir / "anchors"

    # Find which assemblies actually have outputs
    available = []
    for assembly in args.targets:
        if (conv_dir / assembly).exists():
            available.append(assembly)
        else:
            print(f"WARNING: no conversion_tables/{assembly} — skipping",
                  file=sys.stderr)

    if not available:
        print("ERROR: no assemblies found under", conv_dir, file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(available)} assemblies: {available}\n")

    print("=" * 60)
    print("STEP 1: Master anchor Excel")
    print("=" * 60)
    make_master_anchors(
        anchors_dir,
        outdir / "master_anchors.xlsx",
        available
    )

    print()
    print("=" * 60)
    print("STEP 2: Master conversion summary")
    print("=" * 60)
    make_master_summary(conv_dir, outdir, available)

    print("\n" + "=" * 60)
    print("ALL DONE")
    print("=" * 60)
    print(f"  Anchors : {outdir}/master_anchors.xlsx")
    print(f"  Summary : {outdir}/master_conversion_summary.tsv/.xlsx")


if __name__ == "__main__":
    main()
