#!/usr/bin/env python3
"""
Filter FASTA sequences by a header pattern and optionally rename contigs.

Sequences whose description matches the given regex pattern are written to the
output file.  When a mapping file is supplied, each kept sequence is renamed
from its temporary contig ID to the original name before writing.

Usage:
    separate_sequences.py -i input.fna -p "viral_sequence" -o viral.fna
    separate_sequences.py -i input.fna -p "plasmid" -o plasmid.fna --map names.tsv
"""
from __future__ import annotations

import argparse
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace with input (str), pattern (str), output (str),
        and map (list[str] | None).
    """
    parser = argparse.ArgumentParser(
        description="Filter FASTA sequences by header pattern."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument(
        "-p", "--pattern", required=True, help="Regex pattern to match against sequence headers"
    )
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument(
        "--map",
        required=False,
        nargs="+",
        help="TSV mapping file(s) with original and temporary contig names",
    )
    return parser.parse_args()


def read_map(map_files: list[str]) -> dict[str, str]:
    """Read one or more TSV mapping files from temporary to original contig names.

    Each file is expected to have a header row containing 'original' and 'temporary',
    followed by rows with (original_name, temporary_name) tab-separated columns.

    Args:
        map_files: Paths to mapping TSV files.

    Returns:
        Dict mapping temporary contig name -> original contig name.

    Raises:
        SystemExit: If a temporary name appears in more than one mapping entry.
    """
    mapping: dict[str, str] = {}
    for map_file in map_files:
        with open(map_file) as f:
            for line in f:
                if "original" in line and "temporary" in line:
                    continue
                parts = line.strip().split("\t")
                original, temporary = parts[0], parts[1]
                if temporary in mapping:
                    print(f"Mapping already exists {temporary}. Exit")
                    exit(1)
                mapping[temporary] = original
    return mapping


def separate_sequences(
    input_fasta: str,
    pattern: str,
    output_fasta: str,
    mapping: dict[str, str],
) -> None:
    """Filter sequences by header pattern and write matching records to output.

    Sequences whose BioPython description matches ``pattern`` are kept.  If a
    mapping is provided, the record ID and description are replaced with the
    original contig name before writing.

    Args:
        input_fasta: Path to the input FASTA file.
        pattern: Regex pattern tested against each record's full description.
        output_fasta: Path for the output FASTA file.
        mapping: Dict of temporary -> original contig names (may be empty).
    """
    regex = re.compile(pattern)
    kept: list[SeqRecord] = []
    pattern_filtered = 0

    for record in SeqIO.parse(input_fasta, "fasta"):
        if regex.search(record.description):
            if mapping:
                original = mapping.get(record.id, record.id)
                record.id = original
                record.description = original
            kept.append(record)
        else:
            pattern_filtered += 1

    SeqIO.write(kept, output_fasta, "fasta")
    print(f"Sequences kept: {len(kept)}")
    print(f"Sequences excluded (pattern mismatch): {pattern_filtered}")


def main() -> None:
    args = parse_args()
    mapping: dict[str, str] = read_map(args.map) if args.map else {}
    separate_sequences(args.input, args.pattern, args.output, mapping)


if __name__ == "__main__":
    main()
