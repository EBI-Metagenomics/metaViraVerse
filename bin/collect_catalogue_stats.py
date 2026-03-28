#!/usr/bin/env python3
"""
Collect catalogue statistics and write them to a JSON file.

Inputs
------
--viral-sequences   FASTA  → number of sequences (> headers)
--plasmids          FASTA  → number of sequences
--prophages         FASTA  → number of sequences
--metadata          TSV    → metadata table including "biomes" column
--clusters-viruses  TSV    → number of viral clusters (unique values in column 1)
--clusters-plasmids TSV    → number of plasmid clusters (unique values in column 1)
--proteins-viruses  FASTA  → number of protein records (non-comment lines)
--proteins-plasmids FASTA  → number of protein records

All input files may be plain text or gzip-compressed.

Output
------
JSON file with keys matching the stat names above.
"""

import argparse
import csv
import gzip
import json
import sys


def open_file(path):
    """Return a text-mode file handle for a plain or gzip-compressed file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def count_fasta_sequences(path):
    """Count sequences in a FASTA file (number of header lines starting with '>')."""
    count = 0
    with open_file(path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def count_unique_biomes(path, col="biomes"):
    """Count unique values in the given column of a TSV file."""
    biomes = set()
    with open_file(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if col not in (reader.fieldnames or []):
            print(
                f"Error: column '{col}' not found in {path}. "
                f"Available: {reader.fieldnames}",
                file=sys.stderr,
            )
            sys.exit(1)
        for row in reader:
            val = row.get(col, "").strip()
            if val:
                biomes.add(val)
    return len(biomes)


def count_clusters(path):
    """Count unique cluster IDs from column 1 of a TSV (skips header if present)."""
    clusters = set()
    with open_file(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            clusters.add(line.split("\t")[0])
    return len(clusters)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect catalogue statistics into a JSON file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--viral-sequences", required=True,
                        help="FASTA file of viral sequences (plain or .gz).")
    parser.add_argument("--plasmids", required=True,
                        help="FASTA file of plasmids (plain or .gz).")
    parser.add_argument("--prophages", required=True,
                        help="FASTA file of prophages (plain or .gz).")
    parser.add_argument("--metadata", required=True,
                        help="TSV file with a 'biomes' column (plain or .gz).")
    parser.add_argument("--clusters-viruses", required=True,
                        help="TSV with cluster IDs in column 1 (plain or .gz).")
    parser.add_argument("--clusters-plasmids", required=True,
                        help="TSV with cluster IDs in column 1 (plain or .gz).")
    parser.add_argument("--proteins-viruses", required=True,
                        help="FAA file of viral proteins (plain or .gz).")
    parser.add_argument("--proteins-plasmids", required=True,
                        help="FAA file of plasmid proteins (plain or .gz).")
    parser.add_argument("-o", "--output", required=True,
                        help="Output JSON file path.")
    return parser.parse_args()


def main():
    args = parse_args()

    viral_seqs_count = count_fasta_sequences(args.viral_sequences)
    plasmids_count = count_fasta_sequences(args.plasmids)
    prophages_count = count_fasta_sequences(args.prophages)
    stats = {
        "total_sequences": viral_seqs_count + plasmids_count + prophages_count,
        "viral_sequences": viral_seqs_count,
        "plasmids":        plasmids_count,
        "prophages":       prophages_count,
        "number_of_biomes":   count_unique_biomes(args.metadata),
        "viral_clusters": count_clusters(args.clusters_viruses),
        "plasmid_clusters": count_clusters(args.clusters_plasmids),
        "total_proteins_viruses":  count_fasta_sequences(args.proteins_viruses),
        "total_proteins_plasmids": count_fasta_sequences(args.proteins_plasmids),
    }

    with open(args.output, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Stats written to {args.output}", file=sys.stderr)
    for k, v in stats.items():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
