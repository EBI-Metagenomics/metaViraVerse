#!/usr/bin/env python3
"""scripts/filter_pdbs_by_plddt.py

Copies only PDB files with mean pLDDT >= cutoff into outdir.
Used to pre-filter input to ProteinCartography (O(N^2) cost control).

Usage:
    python filter_pdbs_by_plddt.py \
        --pdb-dir        structures/ \
        --confidence-tsv structure_confidence.tsv \
        --plddt-cutoff   0.7 \
        --outdir         filtered_pdbs/
"""

import argparse
import shutil
import sys
from pathlib import Path


def load_passing_ids(tsv: Path, cutoff: float) -> set:
    passing = set()
    with open(tsv) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            seq_id, _, plddt_str, status = parts[0], parts[1], parts[2], parts[3]
            if status != "success":
                continue
            try:
                if float(plddt_str) >= cutoff:
                    passing.add(seq_id)
            except ValueError:
                pass
    return passing


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-dir",        required=True)
    parser.add_argument("--confidence-tsv", required=True)
    parser.add_argument("--plddt-cutoff",   type=float, default=0.7)
    parser.add_argument("--outdir",         required=True)
    args = parser.parse_args()

    outdir  = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    passing = load_passing_ids(Path(args.confidence_tsv), args.plddt_cutoff)

    copied, skipped = 0, 0
    for pdb in sorted(Path(args.pdb_dir).glob("*.pdb")):
        if pdb.stem in passing:
            shutil.copy2(pdb, outdir / pdb.name)
            copied += 1
        else:
            skipped += 1

    print(f"  pLDDT filter: {copied} copied, {skipped} skipped "
          f"(cutoff={args.plddt_cutoff})", file=sys.stderr)
    if copied == 0:
        print("  WARNING: No structures passed pLDDT filter — "
              "ProteinCartography map will be empty.", file=sys.stderr)


if __name__ == "__main__":
    main()
