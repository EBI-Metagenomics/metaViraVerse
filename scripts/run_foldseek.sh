#!/usr/bin/env bash
# scripts/run_foldseek.sh
#
# Foldseek structural homology search against BFVD and optionally PDB100.
#
# BFVD streaming (no local download):
#   Pass 'BFVD' as bfvd_db — Foldseek resolves it via bfvd.foldseek.com.
#   For offline/HPC: download once and pass the local DB directory path.
#
# Output columns:
#   query, target, evalue, bits, alntmscore, qtmscore, ttmscore,
#   lddt, alnlen, pident
#
# Usage:
#   bash run_foldseek.sh <pdb_dir> <bfvd_db> <outdir> <threads> \
#                        <run_pdb: true|false> <pdb_db> <evalue>
#
# Cite:
#   van Kempen et al. Nature Methods 2024 doi:10.1038/s41592-023-02119-x
#   Kim et al. NAR 2024 doi:10.1093/nar/gkae1119

set -euo pipefail

PDB_DIR="$1"
BFVD_DB="$2"
OUTDIR="$3"
THREADS="${4:-8}"
RUN_PDB="${5:-false}"
PDB_DB="${6:-PDB}"
EVALUE="${7:-0.001}"

FORMAT="query,target,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,alnlen,pident"

mkdir -p "${OUTDIR}"

# Count input PDBs — exit cleanly if directory is empty
N_PDB=$(find "${PDB_DIR}" -maxdepth 1 -name "*.pdb" | wc -l)
if [ "${N_PDB}" -eq 0 ]; then
    echo "WARNING: No .pdb files found in ${PDB_DIR} — writing empty output" >&2
    printf "query\ttarget\tevalue\tbits\talntmscore\tqtmscore\tttmscore\tlddt\talnlen\tpident\n" \
        > "${OUTDIR}/bfvd_hits.tsv"
    exit 0
fi

echo "Running Foldseek vs BFVD (${N_PDB} structures, ${THREADS} threads)..."
foldseek easy-search \
    "${PDB_DIR}"                    \
    "${BFVD_DB}"                    \
    "${OUTDIR}/bfvd_hits.tsv"       \
    "${OUTDIR}/tmp_bfvd"            \
    --threads    "${THREADS}"       \
    --format-output "${FORMAT}"     \
    --greedy-best-hits              \
    -e "${EVALUE}"
echo "  BFVD done: ${OUTDIR}/bfvd_hits.tsv"

# ── Optional PDB100 search ────────────────────────────────────
if [ "${RUN_PDB}" = "true" ]; then
    echo "Running Foldseek vs PDB100..."
    foldseek easy-search \
        "${PDB_DIR}"                \
        "${PDB_DB}"                 \
        "${OUTDIR}/pdb_hits.tsv"    \
        "${OUTDIR}/tmp_pdb"         \
        --threads    "${THREADS}"   \
        --format-output "${FORMAT}" \
        --greedy-best-hits          \
        -e "${EVALUE}"
    echo "  PDB100 done: ${OUTDIR}/pdb_hits.tsv"
fi

rm -rf "${OUTDIR}/tmp_bfvd" "${OUTDIR}/tmp_pdb" 2>/dev/null || true
echo "Foldseek search complete."
