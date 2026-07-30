#!/usr/bin/env python3
"""scripts/run_esmfold_api.py

Calls the ESM Metagenomic Atlas fold API for each sequence in a FASTA file.
Writes one .pdb per sequence and a summary TSV (seq_id, length, mean_plddt, status).

API: https://api.esmatlas.com/foldSequence/v1/pdb/
  - POST raw sequence as body; returns PDB text
  - Public, unauthenticated; rate-limited (~1 req/s)
  - Practical ceiling: ~400aa via public server
  - pLDDT stored in B-factor column (0–100); normalised to 0–1 in output TSV

Usage:
    python run_esmfold_api.py \
        --input  filtered_reps.faa \
        --outdir structures/ \
        --summary structure_confidence.tsv \
        --delay  1.2
"""

import argparse
import time
import sys
import requests
from pathlib import Path

ESM_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_RETRIES = 3
BACKOFF = 5.0  # seconds between retries


def parse_fasta(path: Path) -> dict:
    seqs, current = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def mean_plddt_from_pdb(pdb_text: str) -> float | None:
    """Extract mean pLDDT from B-factor column of ATOM records (0–100 scale)."""
    values = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            try:
                values.append(float(line[60:66].strip()))
            except ValueError:
                pass
    return sum(values) / len(values) if values else None


def fold_with_retry(seq: str, seq_id: str, delay: float) -> tuple[str | None, str]:
    """POST to ESM API with retry/backoff. Returns (pdb_text, status)."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.post(
                ESM_API,
                data=seq,
                headers={"Content-Type": "application/x-www-form-urlencoded"},
                timeout=120,
            )
            r.raise_for_status()
            return r.text, "success"
        except requests.RequestException as e:
            print(f"  [{seq_id}] attempt {attempt}/{MAX_RETRIES} failed: {e}",
                  file=sys.stderr)
            if attempt < MAX_RETRIES:
                time.sleep(BACKOFF * attempt)
    return None, "failed"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",   required=True)
    parser.add_argument("--outdir",  required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--delay",   type=float, default=1.2,
                        help="Seconds between API calls (respect rate limits)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    seqs = parse_fasta(Path(args.input))

    if not seqs:
        print("No sequences in input — writing empty summary.", file=sys.stderr)
        Path(args.summary).write_text("seq_id\tlength\tmean_plddt\tstatus\n")
        return

    with open(args.summary, "w") as fh:
        fh.write("seq_id\tlength\tmean_plddt\tstatus\n")
        for seq_id, seq in seqs.items():
            pdb_path = outdir / f"{seq_id}.pdb"
            if pdb_path.exists():
                print(f"  Skipping {seq_id} (already predicted)")
                continue
            print(f"  Predicting {seq_id} (len={len(seq)})...")
            pdb_text, status = fold_with_retry(seq, seq_id, args.delay)
            if pdb_text:
                pdb_path.write_text(pdb_text)
                raw = mean_plddt_from_pdb(pdb_text)
                # Normalise 0–100 → 0–1
                plddt = f"{raw / 100.0:.3f}" if raw is not None else "NA"
            else:
                plddt = "NA"
            fh.write(f"{seq_id}\t{len(seq)}\t{plddt}\t{status}\n")
            time.sleep(args.delay)


if __name__ == "__main__":
    main()
