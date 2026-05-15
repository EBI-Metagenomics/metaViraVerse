#!/usr/bin/env python3
"""scripts/add_taxonomy_overlay.py

Adds taxonomy and pLDDT columns from metaViraVerse outputs to
ProteinCartography's features table for overlay on the structural map.

Usage:
    python add_taxonomy_overlay.py \
        --pc-features proteincartography_out/features.tsv \
        --taxonomy    viral_sequences_vitap_best.tsv \
        --confidence  structure_confidence.tsv \
        --output      proteincartography_out/features_annotated.tsv
"""

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pc-features", required=True)
    parser.add_argument("--taxonomy",    required=True)
    parser.add_argument("--confidence",  required=True)
    parser.add_argument("--output",      required=True)
    args = parser.parse_args()

    tax = {}
    with open(args.taxonomy) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sid = next(iter(row.values()))
            tax[sid] = row.get("taxonomy", "unclassified")

    conf = {}
    with open(args.confidence) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            conf[row["seq_id"]] = row.get("mean_plddt", "NA")

    with open(args.pc_features) as infile, open(args.output, "w") as outfile:
        reader    = csv.DictReader(infile, delimiter="\t")
        fieldnames = list(reader.fieldnames) + ["taxonomy", "mean_plddt"]
        writer    = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            sid = next(iter(row.values()))
            row["taxonomy"]   = tax.get(sid, "unclassified")
            row["mean_plddt"] = conf.get(sid, "NA")
            writer.writerow(row)


if __name__ == "__main__":
    main()
