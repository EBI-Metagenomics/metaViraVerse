#!/usr/bin/env python3
"""
Collect catalogue statistics and write them to a JSON file.

Inputs
------
--initial-json      JSON   → choose_sequences.py's per-category stats.json:
                             number_of_viral_sequences, number_of_viral_sequence_proteins,
                             number_of_prophages, number_of_prophage_proteins,
                             number_of_plasmids, number_of_plasmid_proteins
--clusters-viruses  TSV    → number of viral clusters (unique values in column 1; the reps
                             list of extract_reps_stats.py, no header)
--clusters-plasmids TSV    → number of plasmid clusters (same format)
--map-file          TSV    → rename map of rename_contigs.py: one row per input sequence
                             (number of rows -> total_sequences)
--initial-metadata TSV     → metadata of every unique sequence (after deduplication, before QC)
--metadata          TSV    → metadata table including "biomes" column, filtered
--excluded-metadata TSV    → excluded records with their filter_reason (duplicate, rrna,
                             not_determined)
--viruses-reps-proteins  FASTA → proteins of viral cluster representatives
                                 (number of records -> total_proteins_viruses_reps)
--plasmids-reps-proteins FASTA → proteins of plasmid cluster representatives
                                 (number of records -> total_proteins_plasmids_reps)

All TSV and FASTA input files may be plain text or gzip-compressed.

Output
------
JSON file. Sequence counts:
  total_sequences = duplicate_sequences + deduplicated_sequences
  deduplicated_sequences = qc_excluded_sequences + unique_sequences
  qc_excluded_sequences = rrna_excluded_sequences + not_determined_excluded_sequences
Other keys match the stat names above. viral_sequences/plasmids/
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


# first-column names that mark a header line in a clusters / reps TSV
CLUSTER_HEADER_NAMES = {"object", "cluster", "cluster_rep", "rep_id", "representative", "id", "sequence_id"}


def count_clusters(path):
    """Count unique cluster IDs from column 1 of a TSV.

    The first line is a header when its first column is one of CLUSTER_HEADER_NAMES; the header
    and its repeats (concatenated files) are not counted.
    """
    header = None
    ids = set()
    with open_file(path) as f:
        for line in f:
            if not line.strip():
                continue
            first = line.split("\t")[0].strip()
            if header is None:
                header = first if first.lower() in CLUSTER_HEADER_NAMES else ""
            if header and first == header:
                continue
            ids.add(first)
    return len(ids)


def count_tsv_lines(path):
    count = 0
    with open_file(path) as f:
        for line in f:
            if 'sequence_id' in line:
                continue
            count += 1
    return count


def count_map_rows(path):
    """Number of sequences in the rename map (header lines skipped)."""
    count = 0
    with open_file(path) as f:
        for line in f:
            if line.strip() and not line.startswith("original\t"):
                count += 1
    return count


def count_excluded_by_reason(path):
    """{filter_reason: number of rows} of the excluded table (first column)."""
    counts = {}
    with open_file(path) as f:
        for line in f:
            reason = line.split("\t")[0].strip()
            if not reason or reason == "filter_reason":
                continue
            counts[reason] = counts.get(reason, 0) + 1
    return counts


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
    parser.add_argument("--map-file", required=True,
                        help="Rename map (one row per input sequence; plain or .gz).")
    parser.add_argument("--initial-metadata", required=True,
                        help="TSV of every unique sequence (after deduplication, before QC).")
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

    excluded = count_excluded_by_reason(args.excluded_metadata)
    unknown = set(excluded) - {"duplicate", "rrna", "not_determined"}
    if unknown:
        print(f"Warning: unknown filter_reason value(s) in {args.excluded_metadata}: {sorted(unknown)}",
              file=sys.stderr)
    rrna = excluded.get("rrna", 0)
    not_determined = excluded.get("not_determined", 0)

    stats = {
        "total_sequences": count_map_rows(args.map_file),
        "duplicate_sequences": excluded.get("duplicate", 0),
        "deduplicated_sequences": count_tsv_lines(args.initial_metadata),
        "rrna_excluded_sequences": rrna,
        "not_determined_excluded_sequences": not_determined,
        "qc_excluded_sequences": rrna + not_determined,
        "unique_sequences": count_tsv_lines(args.metadata),
        "viral_sequences": initial_stats["number_of_viral_sequences"],
        "plasmids":        initial_stats["number_of_plasmids"],
        "prophages":       initial_stats["number_of_prophages"],
        "number_of_biomes":   count_unique_biomes(args.metadata),
        "viral_clusters": count_clusters(args.clusters_viruses),
        "plasmid_clusters": count_clusters(args.clusters_plasmids),
        "total_proteins_viruses": (
            initial_stats["number_of_viral_sequence_proteins"] + initial_stats["number_of_prophage_proteins"]
        ),
        "total_proteins_plasmids": initial_stats["number_of_plasmid_proteins"],
        "total_proteins_viruses_reps": count_fasta_records(args.viruses_reps_proteins),
        "total_proteins_plasmids_reps": count_fasta_records(args.plasmids_reps_proteins),
    }

    # the counts must add up; a mismatch is reported but does not stop the run
    for label, left, right in [
        ("total_sequences = duplicate_sequences + deduplicated_sequences",
         stats["total_sequences"], stats["duplicate_sequences"] + stats["deduplicated_sequences"]),
        ("deduplicated_sequences = qc_excluded_sequences + unique_sequences",
         stats["deduplicated_sequences"], stats["qc_excluded_sequences"] + stats["unique_sequences"]),
    ]:
        if left != right:
            print(f"Warning: {label} does not hold ({left} != {right})", file=sys.stderr)

    with open(args.output, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Stats written to {args.output}", file=sys.stderr)
    for k, v in stats.items():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
