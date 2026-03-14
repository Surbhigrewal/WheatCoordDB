#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=liftoff_array
#SBATCH --partition=your_partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=160g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your@email.com
#SBATCH --output=outputs/logs/liftoff_%A_%a.out
#SBATCH --error=outputs/logs/liftoff_%A_%a.err
#
# Each array job runs one target assembly through all steps:
#   1. Liftoff
#   2. Extract anchors
#   3. Build conversion tables + per-chr dotplots
#   4. Rename tgt_* columns to <assembly>_tgt_* in all output files
#   5. Generate 21-chr dotplot for this assembly
#

#
# Run from CS_synteny/ directory:
#   sbatch --array=0-9 2_run_liftoff_array.sh      # 10+ panel
#   sbatch --array=0-24 2_run_liftoff_array.sh     # expanded set

set -eo pipefail

source config.sh

CPUS=48  # matches --ntasks-per-node above

# ---- Conda activate ----
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

# ---- This job's target ----
IDX=$SLURM_ARRAY_TASK_ID
NAME="${TARGET_NAMES[$IDX]}"
GENOME="${TARGET_GENOMES[$IDX]}"

# ---- Translocation map (optional — only affects ArinaLrFor and SY_Mattis) ----
TRANSLOCATION_MAP="translocation_map.tsv"

# ---- Output dirs ----
LIFTOFF_OUT="$OUTDIR/liftoff/$NAME"
mkdir -p "$LIFTOFF_OUT/intermediate" \
         "$OUTDIR/logs" \
         "$OUTDIR/anchors" \
         "$OUTDIR/conversion_tables/$NAME" \
         "$OUTDIR/plots/$NAME"

echo "================================================"
echo "  Task    : $IDX"
echo "  Target  : $NAME"
echo "  CPUs    : $CPUS"
echo "  Genome  : $GENOME"
echo "  CS FASTA: $CS_GENOME"
echo "  CS GFF  : $CS_GFF"
echo "  Output  : $LIFTOFF_OUT"
echo "  Started : $(date)"
echo "================================================"

# ---- Validate inputs exist ----
for F in "$CS_GENOME" "$CS_GFF" "$CS_CHROM_SIZES" "$GENOME"; do
    if [[ ! -f "$F" ]]; then
        echo "ERROR: file not found: $F"
        exit 1
    fi
done

# ============================================================
# STEP 1 — Liftoff
# ============================================================
echo "[1/5] Running Liftoff... $(date)"

liftoff \
    -g "$CS_GFF" \
    -o "$LIFTOFF_OUT/${NAME}_liftoff.gff3" \
    -u "$LIFTOFF_OUT/${NAME}_unmapped.txt" \
    -dir "$LIFTOFF_OUT/intermediate" \
    -p "$CPUS" \
    -m minimap2 \
    -a 0.9 \
    -s 0.9 \
    "$GENOME" "$CS_GENOME" \
    2>&1 | tee "$OUTDIR/logs/liftoff_${NAME}.log"

echo "[1/5] Liftoff done: $NAME — $(date)"

# ============================================================
# STEP 2 — Extract anchor pairs
# ============================================================
echo "[2/5] Extracting anchors... $(date)"

python3 01_extract_anchors.py \
    --cs-gff        "$CS_GFF" \
    --liftoff-gff   "$LIFTOFF_OUT/${NAME}_liftoff.gff3" \
    --output        "$OUTDIR/anchors/${NAME}_anchors.tsv" \
    --feature-type  "$ANCHOR_FEATURE" \
    --min-coverage  "$MIN_COVERAGE" \
    --min-identity  "$MIN_IDENTITY" \
    --assembly-name "$NAME" \
    --translocation-map "$TRANSLOCATION_MAP"

echo "[2/5] Anchors done: $NAME — $(date)"

# ============================================================
# STEP 3 — Build conversion tables + per-chr dotplots
# ============================================================
echo "[3/5] Building conversion tables... $(date)"

python3 02_build_conversion_table.py \
    --anchors             "$OUTDIR/anchors/${NAME}_anchors.tsv" \
    --cs-chrom-sizes      "$CS_CHROM_SIZES" \
    --output-dir          "$OUTDIR/conversion_tables/$NAME" \
    --target-name         "$NAME" \
    --resolution          "$CONV_RESOLUTION" \
    --min-anchors-per-chr "$MIN_ANCHORS" \
    --plot-dir            "$OUTDIR/plots/$NAME"

echo "[3/5] Conversion tables done: $NAME — $(date)"

# ============================================================
# STEP 4 — Rename tgt_* columns to <assembly>_tgt_* in all files
# ============================================================
echo "[4/5] Renaming tgt_* columns... $(date)"

python3 - <<PYEOF
import sys
from pathlib import Path
import pandas as pd

assembly = "$NAME"
conv_dir  = Path("$OUTDIR/conversion_tables") / assembly
anchors_dir = Path("$OUTDIR/anchors")

ANCHOR_TGT_COLS = ["tgt_chr", "tgt_start", "tgt_end", "tgt_strand", "tgt_midpoint"]
CONV_TGT_COLS   = ["tgt_chr", "tgt_pos"]
BLOCK_TGT_COLS  = ["tgt_chr", "tgt_start", "tgt_end"]
SUMMARY_TGT_COLS = ["tgt_chr"]

def prefix_cols(df, tgt_cols, asm):
    rename = {c: f"{asm}_{c}" for c in tgt_cols if c in df.columns}
    return df.rename(columns=rename) if rename else df

def already_renamed(df, asm):
    return any(c.startswith(f"{asm}_") for c in df.columns)

# Chr*_anchors.tsv (in conversion_tables/)
for p in sorted(conv_dir.glob("Chr*_anchors.tsv")):
    df = pd.read_csv(p, sep="\t")
    if already_renamed(df, assembly): continue
    df = prefix_cols(df, ANCHOR_TGT_COLS + ["coverage", "identity"], assembly)
    df.to_csv(p, sep="\t", index=False)
    print(f"  renamed: {p.name}")

# Chr*_conversion.tsv.gz
for p in sorted(conv_dir.glob("Chr*_conversion.tsv.gz")):
    df = pd.read_csv(p, sep="\t", compression="gzip")
    if already_renamed(df, assembly): continue
    df = prefix_cols(df, CONV_TGT_COLS, assembly)
    df.to_csv(p, sep="\t", index=False, compression="gzip")
    print(f"  renamed: {p.name}")

# Chr*_synteny_blocks.tsv
for p in sorted(conv_dir.glob("Chr*_synteny_blocks.tsv")):
    df = pd.read_csv(p, sep="\t")
    if already_renamed(df, assembly): continue
    df = prefix_cols(df, BLOCK_TGT_COLS, assembly)
    df.to_csv(p, sep="\t", index=False)
    print(f"  renamed: {p.name}")

# Chr*_synteny.bed (header comment)
for p in sorted(conv_dir.glob("Chr*_synteny.bed")):
    text = p.read_text()
    if f"{assembly}_tgt" in text: continue
    text = text.replace("tgt_chr",   f"{assembly}_tgt_chr") \
               .replace("tgt_start", f"{assembly}_tgt_start") \
               .replace("tgt_end",   f"{assembly}_tgt_end")
    p.write_text(text)
    print(f"  renamed: {p.name}")

# conversion_summary.tsv
summary_p = conv_dir / "conversion_summary.tsv"
if summary_p.exists():
    df = pd.read_csv(summary_p, sep="\t")
    if not already_renamed(df, assembly):
        df = prefix_cols(df, SUMMARY_TGT_COLS, assembly)
        df.to_csv(summary_p, sep="\t", index=False)
        print(f"  renamed: conversion_summary.tsv")

# outputs/anchors/<assembly>_anchors.tsv
anchor_p = anchors_dir / f"{assembly}_anchors.tsv"
if anchor_p.exists():
    df = pd.read_csv(anchor_p, sep="\t")
    if not already_renamed(df, assembly):
        df = prefix_cols(df, ANCHOR_TGT_COLS + ["coverage", "identity"], assembly)
        df.to_csv(anchor_p, sep="\t", index=False)
        print(f"  renamed: anchors/{assembly}_anchors.tsv")

print(f"Renaming complete for {assembly}")
PYEOF

echo "[4/5] Renaming done: $NAME — $(date)"

# ============================================================
# STEP 5 — 21-chromosome dotplot for this assembly
# ============================================================
echo "[5/5] Generating 21-chr dotplot... $(date)"

python3 - <<PYEOF
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd

assembly  = "$NAME"
conv_dir  = Path("$OUTDIR/conversion_tables")
plots_dir = Path("$OUTDIR/plots")

CHROMOSOMES = [f"Chr{g}{s}" for g in range(1, 8) for s in ["A", "B", "D"]]
NROWS, NCOLS = 7, 3

target_dir = conv_dir / assembly
fig, axes = plt.subplots(NROWS, NCOLS, figsize=(14, 22), squeeze=False)
fig.suptitle(
    f"CS RefSeq v2.1  →  {assembly}\nSynteny anchors — all 21 chromosomes",
    fontsize=13, fontweight="bold", y=0.998
)

cmap = plt.cm.plasma_r
norm = mcolors.Normalize(vmin=0.9, vmax=1.0)
has_identity = False

for idx, chrom in enumerate(CHROMOSOMES):
    row, col = divmod(idx, NCOLS)
    ax = axes[row][col]
    ax.set_title(f"{chrom[3]}{chrom[4]}", fontsize=9, pad=3)

    anchor_file = target_dir / f"{chrom}_anchors.tsv"
    if not anchor_file.exists():
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        continue

    df = pd.read_csv(anchor_file, sep="\t")
    if df.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        continue

    cs_mid_col  = "cs_midpoint"
    tgt_mid_col = next((c for c in df.columns if c.endswith("tgt_midpoint")), None)
    id_col      = next((c for c in df.columns if c.endswith("identity")), None)

    if cs_mid_col not in df.columns or tgt_mid_col is None:
        ax.text(0.5, 0.5, "column\nerror", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="red")
        continue

    x = df[cs_mid_col] / 1e6
    y = df[tgt_mid_col] / 1e6

    if id_col:
        c = df[id_col].clip(0.9, 1.0)
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
        ax.set_ylabel(f"{assembly} (Mb)", fontsize=7)
    else:
        ax.tick_params(labelleft=False)
    ax.tick_params(axis="both", labelsize=6)
    ax.text(0.97, 0.03, f"n={len(df):,}", ha="right", va="bottom",
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

plots_dir.mkdir(parents=True, exist_ok=True)
out_path = plots_dir / f"{assembly}_21chr_dotplot.png"
fig.savefig(out_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out_path}")
PYEOF

echo "[5/5] Dotplot done: $NAME — $(date)"
echo "================================================"
echo "ALL DONE: $NAME (task $IDX) — $(date)"
echo "================================================"
echo ""
echo "NOTE: Once ALL array jobs finish, run once from CS_synteny/:"
echo "  python3 04_postprocess.py --skip-rename --skip-plots"
echo "to generate master_anchors.xlsx and master_conversion_summary.xlsx"
