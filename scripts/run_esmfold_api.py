#!/usr/bin/env python3
"""scripts/run_esmfold_api.py

Calls the ESM Metagenomic Atlas fold API for each sequence in a FASTA file.
Writes one .pdb per sequence and a summary TSV with pLDDT scores.

API endpoint: https://api.esmatlas.com/foldSequence/v1/pdb/
  - POST with raw sequence string as body
  - Returns PDB-format text
  - Public, unauthenticated, rate-limited (~1 req/s recommended)
  - Max sequence length: ~400aa via public server

Usage:
    python run_esmfold_api.py \
        --input filtered_reps.faa \
        --outdir structures/ \
        --summary structure_confidence.tsv \
        --delay 1.2
"""

import argparse
import time
import requests
from pathlib import Path

ESM_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"


def parse_fasta(path):
    seqs = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def extract_plddt_from_pdb(pdb_text):
    """ESMFold stores per-residue pLDDT in the B-factor column of ATOM records."""
    bfactors = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            try:
                bfactors.append(float(line[60:66].strip()))
            except ValueError:
                pass
    return sum(bfactors) / len(bfactors) if bfactors else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--delay", type=float, default=1.2)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    seqs = parse_fasta(Path(args.input))

    with open(args.summary, "w") as summary:
        summary.write("seq_id\tlength\tmean_plddt\tstatus\n")
        for seq_id, seq in seqs.items():
            pdb_path = outdir / f"{seq_id}.pdb"
            if pdb_path.exists():
                print(f"  Skipping {seq_id} (already predicted)")
                continue
            print(f"  Predicting {seq_id} (len={len(seq)})...")
            try:
                response = requests.post(
                    ESM_API,
                    data=seq,
                    headers={"Content-Type": "application/x-www-form-urlencoded"},
                    timeout=120,
                )
                response.raise_for_status()
                pdb_text = response.text
                pdb_path.write_text(pdb_text)
                mean_plddt = extract_plddt_from_pdb(pdb_text)
                # pLDDT from ESM API is on 0-100 scale; normalise to 0-1
                normalised = mean_plddt / 100.0 if mean_plddt is not None else None
                summary.write(f"{seq_id}\t{len(seq)}\t{normalised:.3f}\tsuccess\n")
            except requests.RequestException as e:
                print(f"  FAILED {seq_id}: {e}")
                summary.write(f"{seq_id}\t{len(seq)}\tNA\tfailed\n")
            time.sleep(args.delay)


if __name__ == "__main__":
    main()
