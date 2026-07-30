#!/usr/bin/env bash
# scripts/run_proteincartography.sh
#
# Runs ProteinCartography in "Cluster mode" using pre-predicted PDB structures.
# Output: interactive HTML structural landscape map.
#
# Install:
#   git clone https://github.com/Arcadia-Science/ProteinCartography
#   cd ProteinCartography && conda env create -f envs/dev.yml
#
# Usage:
#   bash run_proteincartography.sh <structures_dir> <outdir> <threads>
#
# Caveats:
#   - Best for proteins <1200aa; long/multi-domain proteins cluster poorly
#   - All-vs-all TM-score is O(N^2); for >1000 structures ensure ample RAM
#   - Linux or macOS only
#
# Cite: Bigge et al. Arcadia Science 2024

set -euo pipefail

STRUCTURES_DIR=$(realpath "$1")
OUTDIR=$(realpath "$2")
THREADS="${3:-8}"
CART_DIR="${4:-ProteinCartography}"

mkdir -p "${OUTDIR}"

conda run -n proteincartography \
    snakemake \
        --snakefile "${CART_DIR}/Snakefile" \
        --configfile "${CART_DIR}/config/config.yaml" \
        --config \
            input_dir="${STRUCTURES_DIR}" \
            output_dir="${OUTDIR}" \
            mode="cluster" \
        --cores "${THREADS}" \
        --use-conda \
        --rerun-incomplete

echo "ProteinCartography complete. Interactive map: ${OUTDIR}/protein_map.html"
