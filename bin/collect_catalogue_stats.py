#!/usr/bin/env python3
"""
Collect catalogue statistics and write them to a JSON file.

Inputs
------
--initial-json      JSON   → choose_sequences.py's per-category stats.json:
                             number_of_viral_sequences, number_of_viral_sequence_proteins,
                             number_of_prophages, number_of_prophage_proteins,
                             number_of_plasmids, number_of_plasmid_proteins
--clusters-viruses  TSV    → number of viral clusters (unique values in column 1)
--clusters-plasmids TSV    → number of plasmid clusters (unique values in column 1)
--initial-metadata TSV     → metadata table for initial set of sequences, pre-filtered
--metadata          TSV    → metadata table including "biomes" column, filtered
--excluded-metadata TSV    → metadata of excluded poor quality records
--viruses-reps-proteins  FASTA → proteins of viral cluster representatives
                                 (number of records -> total_proteins_viruses_reps)
--plasmids-reps-proteins FASTA → proteins of plasmid cluster representatives
                                 (number of records -> total_proteins_plasmids_reps)

All TSV and FASTA input files may be plain text or gzip-compressed.

Output
------
JSON file with keys matching the stat names above. viral_sequences/plasmids/
prophages and total_proteins_viruses/total_proteins_plasmids are taken
straight from --initial-json (the viral_sequences and prophage protein
counts are summed into total_proteins_viruses, matching the combined
"viruses" grouping used everywhere else in this stats file).
"""

import argparse
import csv
import gzip
import json
import sys


def open_file(path):
    """Return a text-mode file handle for a plain or gzip-compressed file.

    Compression is detected from the gzip magic bytes, so a compressed file
    without a .gz extension (or bgzip output) is read correctly too.
    """
    with open(path, "rb") as f:
        is_gzip = f.read(2) == b"\x1f\x8b"
    if is_gzip:
        return gzip.open(path, "rt")
    return open(path)


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
                biomes.update(val.split(','))  # biomes can be combined, ex. marine,marine_sediment
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


def count_tsv_lines(path):
    count = 0
    with open_file(path) as f:
        for line in f:
            if 'sequence_id' in line:
                continue
            count += 1
    return count


def count_fasta_records(path):
    """Count records (header lines starting with '>') in a plain or gzipped FASTA."""
    count = 0
    with open_file(path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect catalogue statistics into a JSON file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--initial-json", required=True,
                        help="choose_sequences.py's per-category stats.json.")
    parser.add_argument("--metadata", required=True,
                        help="TSV file with a 'biomes' column (plain or .gz).")
    parser.add_argument("--initial-metadata", required=True,
                        help="TSV file for all input sequences")
    parser.add_argument("--excluded-metadata", required=True,
                        help="TSV file for filtered out sequences")
    parser.add_argument("--clusters-viruses", required=True,
                        help="TSV with cluster IDs in column 1 (plain or .gz).")
    parser.add_argument("--clusters-plasmids", required=True,
                        help="TSV with cluster IDs in column 1 (plain or .gz).")
    parser.add_argument("--viruses-reps-proteins", required=True,
                        help="Protein FASTA of viral cluster representatives (plain or .gz).")
    parser.add_argument("--plasmids-reps-proteins", required=True,
                        help="Protein FASTA of plasmid cluster representatives (plain or .gz).")
    parser.add_argument("-o", "--output", required=True,
                        help="Output JSON file path.")
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.initial_json) as f:
        initial_stats = json.load(f)

    stats = {
        "total_sequences": count_tsv_lines(args.initial_metadata),
        "unique_sequences": count_tsv_lines(args.metadata),
        "qc_excluded_sequences": count_tsv_lines(args.excluded_metadata),
        "viral_sequences": initial_stats["number_of_viral_sequences"],
        "plasmids":        initial_stats["number_of_plasmids"],
        "prophages":       initial_stats["number_of_prophages"],
        "number_of_biomes":   count_unique_biomes(args.metadata),
        "viral_clusters": count_clusters(args.clusters_viruses) - 1,
        "plasmid_clusters": count_clusters(args.clusters_plasmids) - 1,
        "total_proteins_viruses": (
            initial_stats["number_of_viral_sequence_proteins"] + initial_stats["number_of_prophage_proteins"]
        ),
        "total_proteins_plasmids": initial_stats["number_of_plasmid_proteins"],
        "total_proteins_viruses_reps": count_fasta_records(args.viruses_reps_proteins),
        "total_proteins_plasmids_reps": count_fasta_records(args.plasmids_reps_proteins),
    }

    with open(args.output, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Stats written to {args.output}", file=sys.stderr)
    for k, v in stats.items():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
