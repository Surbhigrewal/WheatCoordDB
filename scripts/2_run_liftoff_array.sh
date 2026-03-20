#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=liftoff_array
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=160g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@institution.ac.uk
#SBATCH --output=outputs/logs/liftoff_%A_%a.out
#SBATCH --error=outputs/logs/liftoff_%A_%a.err
#
# Each array job runs one target assembly through all steps:
#   1. Liftoff — project CS gene models onto target assembly
#   2. Extract anchors — anchor pairs with assembly-prefixed columns
#   3. Build conversion tables + per-chr dotplots + 21-chr overview dotplot
#
# Target columns are prefixed with the assembly name at step 2
# (e.g. Jagger_tgt_chr) and remain prefixed in all output files.
# No post-extraction column renaming is required.
#
# After ALL array jobs finish, run once manually:
#   python3 04_postprocess.py
# to generate master_anchors.xlsx and master_conversion_summary.xlsx
#
# Run from CS_synteny/ directory:
#   sbatch --array=0-9 2_run_liftoff_array.sh      # 10+ panel
#   sbatch --array=0-23 2_run_liftoff_array.sh     # all 24 assemblies

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
echo "[1/3] Running Liftoff... $(date)"

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

echo "[1/3] Liftoff done: $NAME — $(date)"

# ============================================================
# STEP 2 — Extract anchor pairs
# Target columns written with assembly prefix (e.g. Jagger_tgt_chr)
# ============================================================
echo "[2/3] Extracting anchors... $(date)"

python3 01_extract_anchors.py \
    --cs-gff            "$CS_GFF" \
    --liftoff-gff       "$LIFTOFF_OUT/${NAME}_liftoff.gff3" \
    --output            "$OUTDIR/anchors/${NAME}_anchors.tsv" \
    --feature-type      "$ANCHOR_FEATURE" \
    --min-coverage      "$MIN_COVERAGE" \
    --min-identity      "$MIN_IDENTITY" \
    --assembly-name     "$NAME" \
    --translocation-map "$TRANSLOCATION_MAP"

echo "[2/3] Anchors done: $NAME — $(date)"

# ============================================================
# STEP 3 — Build conversion tables + per-chr dotplots + 21-chr overview
# ============================================================
echo "[3/3] Building conversion tables... $(date)"

python3 02_build_conversion_table.py \
    --anchors             "$OUTDIR/anchors/${NAME}_anchors.tsv" \
    --cs-chrom-sizes      "$CS_CHROM_SIZES" \
    --cs-gene-density     "outputs/cs_gene_density.tsv.gz" \
    --output-dir          "$OUTDIR/conversion_tables/$NAME" \
    --target-name         "$NAME" \
    --resolution          "$CONV_RESOLUTION" \
    --min-anchors-per-chr "$MIN_ANCHORS" \
    --density-window      5000000 \
    --plot-dir            "$OUTDIR/plots/$NAME" \
    --overview-plot       "$OUTDIR/plots/${NAME}_21chr_dotplot.png"

echo "[3/3] Conversion tables done: $NAME — $(date)"
echo "================================================"
echo "ALL DONE: $NAME (task $IDX) — $(date)"
echo "================================================"
echo ""
echo "NOTE: Once ALL array jobs finish, run once from CS_synteny/:"
echo "  python3 04_postprocess.py"
echo "to generate master_anchors.xlsx and master_conversion_summary.xlsx"
