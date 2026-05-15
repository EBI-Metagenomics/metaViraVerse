#!/usr/bin/env bash
# scripts/run_foldseek.sh
#
# Runs Foldseek structural search against BFVD and optionally PDB100.
# BFVD can be streamed via bfvd.foldseek.com or supplied as a local DB path.
#
# Usage:
#   bash run_foldseek.sh <structures_dir> <bfvd_db> <outdir> <threads> <run_pdb> <pdb_db>

STRUCTURES_DIR="$1"
BFVD_DB="$2"
OUTDIR="$3"
THREADS="${4:-8}"
RUN_PDB="${5:-false}"
PDB_DB="${6:-PDB}"

set -euo pipefail
mkdir -p "${OUTDIR}"

FORMAT="query,target,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,alnlen,pident"

echo "Running Foldseek against BFVD..."
foldseek easy-search \
    "${STRUCTURES_DIR}" \
    "${BFVD_DB}" \
    "${OUTDIR}/bfvd_hits.tsv" \
    "${OUTDIR}/tmp_bfvd" \
    --threads "${THREADS}" \
    --format-output "${FORMAT}" \
    --greedy-best-hits \
    -e 0.001
echo "BFVD search complete: ${OUTDIR}/bfvd_hits.tsv"

if [ "${RUN_PDB}" = "true" ]; then
    echo "Running Foldseek against PDB100..."
    # Download once with: foldseek databases PDB pdb_db tmp
    foldseek easy-search \
        "${STRUCTURES_DIR}" \
        "${PDB_DB}" \
        "${OUTDIR}/pdb_hits.tsv" \
        "${OUTDIR}/tmp_pdb" \
        --threads "${THREADS}" \
        --format-output "${FORMAT}" \
        --greedy-best-hits \
        -e 0.001
    echo "PDB search complete: ${OUTDIR}/pdb_hits.tsv"
fi

rm -rf "${OUTDIR}/tmp_bfvd" "${OUTDIR}/tmp_pdb" 2>/dev/null || true
