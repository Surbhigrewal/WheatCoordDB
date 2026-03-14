#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=rename
#SBATCH --partition=your_partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your@email.com
#
# Renames chromosomes to Chr1A..Chr7D in all assemblies + CS GFF3.
#
# Index 0        = CS (also renames the GFF)
# Indices 1-10   = original 10-wheat panel
# Indices 11-24  = 14 new assemblies (../Other_wheat)
#
# Submit new assemblies only (CS + 10-wheat already done):
#   sbatch --array=11-24 1_run_rename_all.sh
#
# Or run a single assembly interactively (e.g. CS_IAAS):
#   SLURM_ARRAY_TASK_ID=12 bash 1_run_rename_all.sh

set -euo pipefail
source config.sh
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

mkdir -p "$OUTDIR/renamed"

# ── Assembly list ─────────────────────────────────────────────────────────────
# Index 0      = CS (special: also renames the GFF)
# Indices 1-10 = original 10-wheat panel (skip guard will skip if already done)
# Indices 11-24 = 14 new assemblies from ../Other_wheat

NAMES=(
    "CS"              # 0
    "Jagger"          # 1
    "Lancer"          # 2
    "ArinaLrFor"      # 3
    "Stanley"         # 4
    "Spelt"           # 5
    "Mace"            # 6
    "SY_Mattis"       # 7
    "Julius"          # 8
    "Landmark"        # 9
    "Norin61"         # 10
    "Aikang58"        # 11
    "Chunmai104"      # 12
    "CS_IAAS"         # 13
    "CS_CAU"          # 14
    "Sumai3"          # 15
    "JIN50"           # 16
    "MOV"             # 17
    "CS_v1"           # 18
    "Fielder"         # 19
    "Kariega"         # 20
    "Attraktion"      # 21
    "Renan_v2"        # 22
    "Paragon_v3"      # 23
    "Cadenza_v2"      # 24
)

FASTAS=(
    "$CS_GENOME_RAW"                                                                                   # 0  CS
    "$ASSEMBLIES_DIR/GCA_903993795.1_10wheat_assembly_jagger_genomic.fasta"                            # 1  Jagger
    "$ASSEMBLIES_DIR/GCA_903993975.1_10wheat_assembly_lancer_genomic.fasta"                            # 2  Lancer
    "$ASSEMBLIES_DIR/GCA_903993985.1_10wheat_assembly_arinaLrFor_genomic.fna"                          # 3  ArinaLrFor
    "$ASSEMBLIES_DIR/GCA_903994155.1_10wheat_assembly_stanley_genomic.fna"                             # 4  Stanley
    "$ASSEMBLIES_DIR/GCA_903994165.1_10wheat_assembly_spelt_genomic.fna"                               # 5  Spelt
    "$ASSEMBLIES_DIR/GCA_903994175.1_10wheat_assembly_mace_genomic.fna"                                # 6  Mace
    "$ASSEMBLIES_DIR/GCA_903994185.1_10wheat_assembly_sy_mattis_genomic.fna"                           # 7  SY_Mattis
    "$ASSEMBLIES_DIR/GCA_903994195.1_10wheat_assembly_julius_genomic.fna"                              # 8  Julius
    "$ASSEMBLIES_DIR/GCA_903995565.1_10wheat_assembly_landmark1_genomic.fna"                           # 9  Landmark
    "$ASSEMBLIES_DIR/GCA_904066035.1_10wheat_assembly_norin61_genomic.fna"                             # 10 Norin61
    "$OTHER_DIR/GCA_025895885.1_ASM2589588v1_genomic.fna.gz"                                          # 11 Aikang58
    "$OTHER_DIR/GCA_039655515.1_ASM3965551v1_genomic.fna.gz"                                          # 12 Chunmai104
    "$OTHER_DIR/GCA_040256815.2_ASM4025681v2_genomic.fna.gz"                                          # 13 CS_IAAS
    "$OTHER_DIR/GCA_046250135.1_ASM4625013v1_genomic.fna.gz"                                          # 14 CS_CAU
    "$OTHER_DIR/GCA_051988995.1_ASM5198899v1_genomic.fna.gz"                                          # 15 Sumai3
    "$OTHER_DIR/GCA_052626965.1_ASM5262696v1_genomic.fna.gz"                                          # 16 JIN50
    "$OTHER_DIR/GCA_052924595.1_ASM5292459v1_genomic.fna.gz"                                          # 17 MOV
    "$OTHER_DIR/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna.gz"                                      # 18 CS_v1
    "$OTHER_DIR/GCA_907166925.1_wheat_cv_fielder_v1_assembly_genomic.fasta"                            # 19 Fielder
    "$OTHER_DIR/GCA_910594105.1_Tae_Kariega_v1_genomic.fna.gz"                                        # 20 Kariega
    "$OTHER_DIR/GCA_918797515.1_wheat_cv_attraktion_v1_genomic.fna.gz"                                # 21 Attraktion
    "$OTHER_DIR/GCA_937894285.1_Triticum_aestivum_Renan_v2.1_genomic.fna.gz"                          # 22 Renan_v2
    "$OTHER_DIR/GCA_949126075.2_EI_Triticum_aestivum_paragon_v3_genomic.fna.gz"                       # 23 Paragon_v3
    "$OTHER_DIR/GCA_965120215.1_EI_Triticum_aestivum_cadenza_v2_genomic.fna.gz"                       # 24 Cadenza_v2
)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
NAME="${NAMES[$IDX]}"
FASTA="${FASTAS[$IDX]}"
OUTDIR_ASM="$OUTDIR/renamed/$NAME"

echo "================================================"
echo "  Index : $IDX"
echo "  Name  : $NAME"
echo "  Input : $FASTA"
echo "  Output: $OUTDIR_ASM"
echo "================================================"

# Skip if already done
if [[ -f "$OUTDIR_ASM/${NAME}.fa" ]]; then
    echo "SKIP — $OUTDIR_ASM/${NAME}.fa already exists"
    exit 0
fi

# CS gets GFF renamed too; targets get FASTA only
if [[ "$NAME" == "CS" ]]; then
    python3 rename_chromosomes.py \
        --fasta  "$FASTA" \
        --gff    "$CS_GFF_RAW" \
        --outdir "$OUTDIR_ASM" \
        --name   "$NAME"
else
    python3 rename_chromosomes.py \
        --fasta  "$FASTA" \
        --outdir "$OUTDIR_ASM" \
        --name   "$NAME"
fi

echo "DONE: $NAME"
