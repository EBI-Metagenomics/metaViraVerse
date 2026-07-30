#!/usr/bin/env python3
"""scripts/merge_struct_annotations.py

Merges ESMFold confidence, Foldseek BFVD hits, and ECOD annotations
into the existing viral_sequences_reps_stats.tsv, appending 10 structural
columns. Every input sequence gets a row regardless of prediction outcome.

New columns:
  struct_predicted  — yes / no
  mean_plddt        — 0-1 float or NA
  struct_status     — success / failed / not_run
  bfvd_top_hit      — target ID of best BFVD hit or NA
  bfvd_tm_score     — alignment TM-score (0-1) or NA
  bfvd_evalue       — E-value or NA
  bfvd_pident       — sequence identity fraction or NA
  ecod_xgroup       — ECOD X-group or NA
  ecod_hgroup       — ECOD H-group or NA
  ecod_fgroup       — ECOD family or NA

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
import sys
from pathlib import Path

NEW_COLS = [
    "struct_predicted", "mean_plddt",   "struct_status",
    "bfvd_top_hit",     "bfvd_tm_score","bfvd_evalue",  "bfvd_pident",
    "ecod_xgroup",      "ecod_hgroup",  "ecod_fgroup",
]


def load_tsv_by_key(path: Path, key_col: str) -> dict:
    """Load TSV into {key: row_dict}. Raises clearly if key_col missing."""
    data = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if key_col not in (reader.fieldnames or []):
            raise ValueError(
                f"Expected column '{key_col}' in {path.name}. "
                f"Found: {reader.fieldnames}"
            )
        for row in reader:
            data[row[key_col]] = row
    return data


def load_foldseek_best_hits(path: Path) -> dict:
    """
    Load Foldseek output; keep only the best hit per query (highest alntmscore).
    Handles header line if present; gracefully skips malformed lines.
    """
    COLS = ["query", "target", "evalue", "bits", "alntmscore",
            "qtmscore", "ttmscore", "lddt", "alnlen", "pident"]
    best = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("query"):
                continue                      # skip header or blank
            parts = line.split("\t")
            if len(parts) < len(COLS):
                continue
            row = dict(zip(COLS, parts))
            # Strip .pdb suffix that Foldseek appends to query names
            qid = row["query"].removesuffix(".pdb")
            try:
                score = float(row["alntmscore"])
            except ValueError:
                score = 0.0
            prev = float(best[qid].get("alntmscore", 0)) if qid in best else 0.0
            if score > prev:
                best[qid] = row
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps-stats", required=True)
    parser.add_argument("--confidence", required=True)
    parser.add_argument("--bfvd-hits",  required=True)
    parser.add_argument("--ecod",       required=True)
    parser.add_argument("--output",     required=True)
    args = parser.parse_args()

    # Load all annotation sources
    try:
        confidence = load_tsv_by_key(Path(args.confidence), "seq_id")
    except ValueError as e:
        print(f"ERROR loading confidence TSV: {e}", file=sys.stderr)
        sys.exit(1)

    bfvd = load_foldseek_best_hits(Path(args.bfvd_hits))

    try:
        ecod = load_tsv_by_key(Path(args.ecod), "seq_id")
    except ValueError as e:
        print(f"ERROR loading ECOD TSV: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"  Loaded: {len(confidence)} confidence | "
          f"{len(bfvd)} BFVD hits | {len(ecod)} ECOD annotations",
          file=sys.stderr)

    with open(args.reps_stats) as infile, open(args.output, "w") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        if not reader.fieldnames:
            print("ERROR: reps_stats is empty or has no header", file=sys.stderr)
            sys.exit(1)

        # First column is always the seq_id in metaViraVerse stats TSVs
        id_col  = reader.fieldnames[0]
        out_fields = list(reader.fieldnames) + NEW_COLS

        writer = csv.DictWriter(
            outfile,
            fieldnames=out_fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()

        n_written = 0
        for row in reader:
            seq_id = row[id_col]
            conf   = confidence.get(seq_id, {})
            b      = bfvd.get(seq_id, {})
            e      = ecod.get(seq_id, {})

            row.update({
                # ESMFold confidence
                "struct_predicted": "yes" if conf.get("status") == "success" else "no",
                "mean_plddt":       conf.get("mean_plddt", "NA"),
                "struct_status":    conf.get("status",     "not_run"),
                # Foldseek BFVD best hit
                "bfvd_top_hit":     b.get("target",     "NA"),
                "bfvd_tm_score":    b.get("alntmscore", "NA"),
                "bfvd_evalue":      b.get("evalue",     "NA"),
                "bfvd_pident":      b.get("pident",     "NA"),
                # ECOD domain classification
                "ecod_xgroup":      e.get("ecod_xgroup", "NA"),
                "ecod_hgroup":      e.get("ecod_hgroup", "NA"),
                "ecod_fgroup":      e.get("ecod_fgroup", "NA"),
            })
            writer.writerow(row)
            n_written += 1

    print(f"  Written {n_written} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
