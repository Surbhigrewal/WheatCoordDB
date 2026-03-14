#!/bin/bash
# =============================================================
# setup_conda_env.sh
# Creates conda environment with all required tools
# Run once on your HPC before launching the pipeline
# =============================================================
# Usage: bash setup_conda_env.sh [env_name]

ENV_NAME="${1:-liftoff_env}"

echo "Creating conda environment: $ENV_NAME"

mamba create -n "$ENV_NAME" -y \
    -c conda-forge -c bioconda \
    python=3.10 \
    liftoff \
    minimap2 \
    samtools \
    bedtools \
    pandas \
    numpy \
    scipy \
    matplotlib \
    seaborn

echo "Environment $ENV_NAME created."
echo "Activate with: conda activate $ENV_NAME"
echo ""

# Verify installs
source activate "$ENV_NAME"
echo "Tool versions:"
echo "  liftoff:   $(liftoff --version 2>&1 | head -1)"
echo "  minimap2:  $(minimap2 --version 2>&1 | head -1)"
echo "  samtools:  $(samtools --version | head -1)"
echo "  python:    $(python --version)"
echo "  pandas:    $(python -c 'import pandas; print(pandas.__version__)')"
echo "  scipy:     $(python -c 'import scipy; print(scipy.__version__)')"
