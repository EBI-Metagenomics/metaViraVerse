#!/usr/bin/env python3
"""scripts/merge_struct_annotations.py

Merges structure prediction confidence, Foldseek BFVD hits, and ECOD
annotations into the existing viral_sequences_reps_stats.tsv output.

New columns added:
  struct_predicted, mean_plddt, struct_status,
  bfvd_top_hit, bfvd_tm_score, bfvd_evalue, bfvd_pident,
  ecod_xgroup, ecod_hgroup, ecod_fgroup

Usage:
    python merge_struct_annotations.py \
        --reps-stats viral_sequences_reps_stats.tsv \
        --confidence structure_confidence.tsv \
        --bfvd-hits  bfvd_hits.tsv \
        --ecod       ecod_annotations.tsv \
        --output     viral_sequences_reps_struct_stats.tsv
"""

import argparse
import csv
from pathlib import Path


def load_tsv(path: Path, key_col: str = None) -> dict:
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = row[key_col] if key_col else next(iter(row.values()))
            data[key] = row
    return data


def load_foldseek_best_hits(path: Path) -> dict:
    """Keep only best hit per query (highest alntmscore)."""
    cols = ["query", "target", "evalue", "bits", "alntmscore",
            "qtmscore", "ttmscore", "lddt", "alnlen", "pident"]
    best = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < len(cols):
                continue
            row = dict(zip(cols, parts))
            qid = row["query"].replace(".pdb", "")
            try:
                score = float(row["alntmscore"])
            except ValueError:
                score = 0.0
            prev_score = float(best[qid].get("alntmscore", 0)) if qid in best else 0.0
            if score > prev_score:
                best[qid] = row
    return best


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps-stats", required=True)
    parser.add_argument("--confidence", required=True)
    parser.add_argument("--bfvd-hits",  required=True)
    parser.add_argument("--ecod",       required=True)
    parser.add_argument("--output",     required=True)
    args = parser.parse_args()

    confidence = load_tsv(Path(args.confidence), "seq_id")
    bfvd       = load_foldseek_best_hits(Path(args.bfvd_hits))
    ecod       = load_tsv(Path(args.ecod), "seq_id")

    new_cols = [
        "struct_predicted", "mean_plddt", "struct_status",
        "bfvd_top_hit", "bfvd_tm_score", "bfvd_evalue", "bfvd_pident",
        "ecod_xgroup", "ecod_hgroup", "ecod_fgroup",
    ]

    with open(args.reps_stats) as infile, open(args.output, "w") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        writer = csv.DictWriter(
            outfile,
            fieldnames=list(reader.fieldnames) + new_cols,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in reader:
            seq_id = next(iter(row.values()))
            conf = confidence.get(seq_id, {})
            b    = bfvd.get(seq_id, {})
            e    = ecod.get(seq_id, {})
            row.update({
                "struct_predicted": "yes" if conf.get("status") == "success" else "no",
                "mean_plddt":       conf.get("mean_plddt", "NA"),
                "struct_status":    conf.get("status", "not_run"),
                "bfvd_top_hit":     b.get("target",     "NA"),
                "bfvd_tm_score":    b.get("alntmscore", "NA"),
                "bfvd_evalue":      b.get("evalue",     "NA"),
                "bfvd_pident":      b.get("pident",     "NA"),
                "ecod_xgroup":      e.get("ecod_xgroup","NA"),
                "ecod_hgroup":      e.get("ecod_hgroup","NA"),
                "ecod_fgroup":      e.get("ecod_fgroup","NA"),
            })
            writer.writerow(row)

    print(f"Written: {args.output}")


if __name__ == "__main__":
    main()
