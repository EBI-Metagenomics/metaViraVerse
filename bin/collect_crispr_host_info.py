#!/usr/bin/env python3
"""
Collects host lineage information for viral sequences/prophages/plasmids
using CRISPR spacer predictions.

Inputs:
  - predictions file (CRISPRCasFinder-style, spacer vs viral seq with e-values)
  - CRISPR spacer metadata TSV (crispr_id, name, parent)
  - genomes metadata TSV (Genome, ..., Lineage, ...)

Output:
  TSV with columns: viral_seq, host_genome, lineage, evalue
"""

import argparse
import csv
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-p", "--predictions", help="Predictions TSV (e.g. predictions_barley_underscore.tsv)")
    p.add_argument("-c", "--crispr-tsv", help="CRISPR spacer info TSV (crispr_id, name, parent)")
    p.add_argument("-m", "--metadata-tsv", help="Genomes metadata TSV with Genome and Lineage columns")
    p.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
    return p.parse_args()


def load_crispr_parents(crispr_tsv):
    """Return dict: spacer_name -> list of genome IDs (prefix before first underscore)."""
    spacer_parents = {}
    with open(crispr_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            name = row["name"].strip()
            raw_parents = row["parent"].strip()
            genome_ids = []
            for parent in raw_parents.split(","):
                parent = parent.strip()
                if parent:
                    genome_ids.append(parent.split("_")[0])
            spacer_parents[name] = genome_ids
    return spacer_parents


def load_genome_lineages(metadata_tsv):
    """Return dict: genome_id -> lineage string."""
    lineages = {}
    with open(metadata_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            lineages[row["Genome"].strip()] = row["Lineage"].strip()
    return lineages


def parse_predictions(predictions_file):
    """
    Parse predictions TSV, keeping only the single best (lowest e-value) hit
    per viral_seq across all spacers.

    Returns list of (spacer_name, viral_seq, evalue_float, evalue_str).
    """
    best = {}  # viral_seq -> (spacer, evalue_float, evalue_str)
    with open(predictions_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            spacer = parts[0].lstrip(">").strip()
            viral_seq = parts[1].strip()
            evalue_str = parts[2].strip()
            try:
                evalue_float = float(evalue_str)
            except ValueError:
                continue
            if viral_seq not in best or evalue_float < best[viral_seq][1]:
                best[viral_seq] = (spacer, evalue_float, evalue_str)

    return [(spacer, viral_seq, v[1], v[2]) for viral_seq, v in best.items()
            for spacer in [v[0]]]


def main():
    args = parse_args()

    spacer_parents = load_crispr_parents(args.crispr_tsv)
    genome_lineages = load_genome_lineages(args.metadata_tsv)
    predictions = parse_predictions(args.predictions)

    out = open(args.output, "w") if args.output != "-" else sys.stdout
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["viral_seq", "host_genome", "lineage", "evalue"])

    for spacer, viral_seq, _, evalue_str in predictions:
        parents = spacer_parents.get(spacer)
        if not parents:
            writer.writerow([viral_seq, "NA", "NA", evalue_str])
            continue
        for genome_id in parents:
            lineage = genome_lineages.get(genome_id, "NA")
            writer.writerow([viral_seq, genome_id, lineage, evalue_str])

    if args.output != "-":
        out.close()


if __name__ == "__main__":
    main()
