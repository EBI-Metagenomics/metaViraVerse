#!/usr/bin/env python3
"""scripts/filter_faa_for_esmfold.py

Filters a FASTA protein file for ESMFold compatibility:
  - Removes sequences longer than MAX_LEN residues
  - Removes sequences with >5% ambiguous residues (X, B, Z, U, O)
  - Writes a summary of what was excluded and why

Usage:
    python filter_faa_for_esmfold.py \
        --input viral_sequences_reps.faa \
        --output filtered_reps.faa \
        --max-len 1000 \
        --log filter_log.tsv
"""

import argparse
from pathlib import Path

AMBIGUOUS = set("XBZUO")
MAX_LEN_DEFAULT = 1000
AMBIG_FRAC_MAX = 0.05


def parse_fasta(path: Path):
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-len", type=int, default=MAX_LEN_DEFAULT)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    seqs = parse_fasta(Path(args.input))
    passed, excluded = {}, []

    for seq_id, seq in seqs.items():
        if len(seq) > args.max_len:
            excluded.append((seq_id, len(seq), "too_long"))
            continue
        ambig_frac = sum(1 for aa in seq if aa.upper() in AMBIGUOUS) / max(len(seq), 1)
        if ambig_frac > AMBIG_FRAC_MAX:
            excluded.append((seq_id, len(seq), f"high_ambiguity_{ambig_frac:.2f}"))
            continue
        passed[seq_id] = seq

    with open(args.output, "w") as out:
        for seq_id, seq in passed.items():
            out.write(f">{seq_id}\n{seq}\n")

    with open(args.log, "w") as log:
        log.write("seq_id\tlength\treason_excluded\n")
        for row in excluded:
            log.write("\t".join(str(x) for x in row) + "\n")

    print(f"Passed: {len(passed)} | Excluded: {len(excluded)}")


if __name__ == "__main__":
    main()
