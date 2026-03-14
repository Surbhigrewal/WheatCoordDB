#!/bin/bash
# =============================================================================
# config.sh — CS v2.1 vs 10-wheat panel + 14 new assemblies
# Edit the paths marked !! — everything else is pre-filled
# =============================================================================

# =============================================================================
# !! EDIT THESE PATHS !!
# =============================================================================
CS_DIR="/path/to/reference_dir"                  # Dir with both raw Ensembl CS files
ASSEMBLIES_DIR="/path/to/assemblies_dir"          # Dir with original 10-wheat raw .fna/.fasta files
OTHER_DIR="/path/to/additional_assemblies_dir"       # Dir with 14 new assembly raw files
OUTDIR="outputs"                 # All pipeline outputs go here
# =============================================================================

# ---- Raw CS inputs (Ensembl Plants r62, used by rename_all.sh) ----
CS_GENOME_RAW="$CS_DIR/Triticum_aestivum_refseqv2.IWGSC_RefSeq_v2.1.dna.toplevel.fa.gz"
CS_GFF_RAW="$CS_DIR/Triticum_aestivum_refseqv2.IWGSC_RefSeq_v2.1.62.chr.gff3.gz"

# ---- Renamed outputs (Chr1A..Chr7D) — used by liftoff pipeline ----
CS_GENOME="$OUTDIR/renamed/CS/CS.fa"
CS_GFF="$OUTDIR/renamed/CS/CS.gff3"
CS_CHROM_SIZES="$OUTDIR/renamed/CS/CS.chrom.sizes"

# ---- Target assemblies (renamed outputs) ----
# Indices 0-9:  original 10-wheat panel
# Indices 10-23: 14 new assemblies

TARGET_NAMES=(
    "Jagger"           # 0
    "Lancer"           # 1
    "ArinaLrFor"       # 2
    "Stanley"          # 3
    "Spelt"            # 4
    "Mace"             # 5
    "SY_Mattis"        # 6
    "Julius"           # 7
    "Landmark"         # 8
    "Norin61"          # 9
    "Aikang58"         # 10
    "Chunmai104"       # 11
    "CS_IAAS"          # 12
    "CS_CAU"           # 13
    "Sumai3"           # 14
    "JIN50"            # 15
    "MOV"              # 16
    "CS_v1"            # 17
    "Fielder"          # 18
    "Kariega"          # 19
    "Attraktion"       # 20
    "Renan_v2"         # 21
    "Paragon_v3"       # 22
    "Cadenza_v2"       # 23
)

TARGET_GENOMES=(
    "$OUTDIR/renamed/Jagger/Jagger.fa"
    "$OUTDIR/renamed/Lancer/Lancer.fa"
    "$OUTDIR/renamed/ArinaLrFor/ArinaLrFor.fa"
    "$OUTDIR/renamed/Stanley/Stanley.fa"
    "$OUTDIR/renamed/Spelt/Spelt.fa"
    "$OUTDIR/renamed/Mace/Mace.fa"
    "$OUTDIR/renamed/SY_Mattis/SY_Mattis.fa"
    "$OUTDIR/renamed/Julius/Julius.fa"
    "$OUTDIR/renamed/Landmark/Landmark.fa"
    "$OUTDIR/renamed/Norin61/Norin61.fa"
    "$OUTDIR/renamed/Aikang58/Aikang58.fa"
    "$OUTDIR/renamed/Chunmai104/Chunmai104.fa"
    "$OUTDIR/renamed/CS_IAAS/CS_IAAS.fa"
    "$OUTDIR/renamed/CS_CAU/CS_CAU.fa"
    "$OUTDIR/renamed/Sumai3/Sumai3.fa"
    "$OUTDIR/renamed/JIN50/JIN50.fa"
    "$OUTDIR/renamed/MOV/MOV.fa"
    "$OUTDIR/renamed/CS_v1/CS_v1.fa"
    "$OUTDIR/renamed/Fielder/Fielder.fa"
    "$OUTDIR/renamed/Kariega/Kariega.fa"
    "$OUTDIR/renamed/Attraktion/Attraktion.fa"
    "$OUTDIR/renamed/Renan_v2/Renan_v2.fa"
    "$OUTDIR/renamed/Paragon_v3/Paragon_v3.fa"
    "$OUTDIR/renamed/Cadenza_v2/Cadenza_v2.fa"
)

# ---- Liftoff settings ----
ANCHOR_FEATURE="gene"
MIN_COVERAGE=0.9
MIN_IDENTITY=0.9

# ---- Conversion table settings ----
CONV_RESOLUTION=1000
MIN_ANCHORS=20

# ---- Conda environment ----
CONDA_ENV="liftoff_env"

# ---- SLURM settings ----
USE_SLURM="true"
SLURM_CPUS=16
SLURM_MEM="64G"
SLURM_TIME="12:00:00"
SLURM_PARTITION="your_partition"
LOCAL_CPUS=8
