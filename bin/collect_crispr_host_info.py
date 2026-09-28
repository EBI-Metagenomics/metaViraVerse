#!/usr/bin/env python3
"""
Collects host lineage information for viral sequences/prophages/plasmids
using CRISPR spacer predictions.

Inputs:
  - predictions file (CRISPRCasFinder-style, spacer vs viral seq with e-values)
  - CRISPR spacer metadata TSV (crispr_id, name, parent)
  - genomes metadata TSV (Genome, ..., Lineage, ...)
  - rename mapping TSV (original, temporary, short, ...) from rename_contigs.py

Output:
  TSV with columns: viral_seq, host_genome, lineage, evalue
"""

import argparse
import csv
import json
import sys

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-p", "--predictions", help="Predictions TSV (e.g. predictions_barley_underscore.tsv)")
    p.add_argument("-c", "--crispr-tsv", help="CRISPR spacer info TSV (crispr_id, name, parent)")
    p.add_argument("-m", "--metadata-tsv", help="Genomes metadata TSV with Genome and Lineage columns")
    p.add_argument("-r", "--rename-map", required=True, help="Rename mapping TSV with temporary and short columns (from rename_contigs.py)")
    p.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
    p.add_argument("-s", "--stats-json", default="lineage_comparison.json",
                   help="Per-rank host vs viral-sequence MAG lineage comparison JSON (default: lineage_comparison.json)")
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


def load_rename_map(rename_map_tsv):
    """Return dict: temporary name -> short (original) name."""
    rename_map = {}
    with open(rename_map_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rename_map[row["temporary"].strip()] = row["short"].strip()
    return rename_map


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


def lineage_part(parts, idx):
    """Return the lineage part at idx, or '' if the lineage is too short to have it."""
    return parts[idx] if len(parts) > idx else ''


def build_stats(host_lineage, viral_lineage):
    """
    Compare host and viral-sequence MAG lineages rank by rank.

    Returns dict: rank -> {"same": bool, "viral": str, "host": str}
    """
    host_parts = host_lineage.split(';')
    viral_parts = viral_lineage.split(';')
    stats = {}
    for idx, rank in enumerate(RANKS):
        viral = lineage_part(viral_parts, idx)
        host = lineage_part(host_parts, idx)
        stats[rank] = {"same": viral == host, "viral": viral, "host": host}
    return stats


def new_summary():
    return {rank: {"same": 0, "different": 0, "differences": {}} for rank in RANKS}


def add_to_summary(summary, stats):
    """Accumulate one build_stats() result into the per-rank summary."""
    for rank, comparison in stats.items():
        if comparison["same"]:
            summary[rank]["same"] += 1
        else:
            summary[rank]["different"] += 1
            key = f'{comparison["viral"]} from {comparison["host"]}'
            differences = summary[rank]["differences"]
            differences[key] = differences.get(key, 0) + 1


def main():
    args = parse_args()

    spacer_parents = load_crispr_parents(args.crispr_tsv)
    genome_lineages = load_genome_lineages(args.metadata_tsv)
    rename_map = load_rename_map(args.rename_map)
    predictions = parse_predictions(args.predictions)

    out = open(args.output, "w") if args.output != "-" else sys.stdout
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["viral_seq", "host_genome", "host_lineage", "evalue", "viral_seq_lineage", "match"])

    summary = new_summary()
    missing_lineage = 0
    for spacer, viral_seq, _, evalue_str in predictions:
        parents = spacer_parents.get(spacer)
        if not parents:
            writer.writerow([viral_seq, "NA", "NA", evalue_str])
            continue
        for genome_id in parents:
            lineage = genome_lineages.get(genome_id, "NA")
            viral_seq_short = rename_map.get(viral_seq, viral_seq)
            viral_seq_mag = viral_seq_short.split('_')[0]
            viral_seq_mag_lineage = genome_lineages.get(viral_seq_mag, "NA")
            if viral_seq_mag_lineage == lineage:
                match = "Yes"
            else:
                match = "No"

            if lineage == "NA" or viral_seq_mag_lineage == "NA":
                missing_lineage += 1
            else:
                add_to_summary(summary, build_stats(lineage, viral_seq_mag_lineage))
            writer.writerow([viral_seq, genome_id, lineage, evalue_str, viral_seq_mag_lineage, match])

    for rank in RANKS:
        print(f'Same {rank} {summary[rank]["same"]}, different {summary[rank]["different"]}')
    with open(args.stats_json, 'w') as file_out:
        json.dump({"missing_lineage": missing_lineage, "ranks": summary}, file_out, indent=2)

    if args.output != "-":
        out.close()


if __name__ == "__main__":
    main()
