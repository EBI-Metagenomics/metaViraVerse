#!/usr/bin/env python3
"""
Split a renamed, combined FASTA/GFF/FAA set into virus/prophage/plasmid groups.

Sequences are bucketed by the 'definition' column of the rename map (written
by rename_contigs.py), stripping any 'third_party_' prefix so that, for
example, both 'virus' and 'third_party_virus' land in the same 'virus'
bucket. The matching GFF and (optionally) FAA records are split the same way,
so cross-category deduplication downstream (choose_sequences.py) has a
consistent fna/gff/faa triple to work with per category.

Usage:
    separate_sequences.py \
        --fna combined.fna --gff combined.gff --faa combined.faa \
        --map combined.tsv --category virus --output-prefix viruses
"""
from __future__ import annotations

import argparse
import csv

from Bio import SeqIO

from utils import parse_attributes


CATEGORIES = ('virus', 'prophage', 'plasmid')


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace with fna, gff, faa, map, category, output_prefix.
    """
    parser = argparse.ArgumentParser(
        description="Split a combined FASTA/GFF/FAA set into a virus/prophage/plasmid group."
    )
    parser.add_argument("--fna", required=True, help="Combined (renamed) input FASTA file")
    parser.add_argument("--gff", required=False, help="Combined (renamed) input GFF file")
    parser.add_argument("--faa", required=False, help="Combined input protein FASTA file (original protein IDs)")
    parser.add_argument("--map", required=True, help="Rename map TSV (must include a 'definition' column)")
    parser.add_argument("--category", required=True, choices=CATEGORIES, help="Category to extract")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
    return parser.parse_args()


def read_definitions(map_file: str) -> dict[str, str]:
    """Read the map file and return {temporary_name: base_category}.

    The 'definition' column may carry a 'third_party_' prefix (e.g.
    'third_party_virus'); that prefix is stripped here so third-party and
    MGnify sequences of the same biological category land in the same bucket.

    Args:
        map_file: Path to the rename map TSV (from rename_contigs.py).

    Returns:
        Dict mapping temporary contig name -> base category ('virus',
        'prophage', 'plasmid', or 'NA').
    """
    definitions: dict[str, str] = {}
    with open(map_file) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            temporary = row['temporary']
            definition = row.get('definition') or 'NA'
            definitions[temporary] = definition.replace('third_party_', '')
    return definitions


def split_fasta(fna_file: str, definitions: dict[str, str], category: str, output_fna: str) -> set[str]:
    """Write sequences whose definition matches ``category`` to ``output_fna``.

    Args:
        fna_file: Path to the combined (renamed) input FASTA file.
        definitions: Dict of temporary contig name -> base category (see ``read_definitions``).
        category: Category to keep ('virus', 'prophage', or 'plasmid').
        output_fna: Path for the output FASTA file.

    Returns:
        Set of temporary contig IDs that were kept.
    """
    kept_ids: set[str] = set()
    with open(output_fna, 'w') as out_f:
        for record in SeqIO.parse(fna_file, "fasta"):
            if definitions.get(record.id) == category:
                SeqIO.write(record, out_f, "fasta")
                kept_ids.add(record.id)
    print(f"Wrote {len(kept_ids)} sequences to {output_fna}")
    return kept_ids


def split_gff(gff_file: str, kept_ids: set[str], output_gff: str) -> set[str]:
    """Write GFF records belonging to ``kept_ids`` to ``output_gff``.

    Sequence identity is read straight from column 1 (the GFF seqid).
    rename_contigs.py already rewrites column 1 to the same temporary name
    used in the combined FASTA and the rename map for every line belonging to
    a sequence (including its CDS lines), so no attribute parsing or pattern
    matching is needed here to find which sequence a line belongs to.

    Args:
        gff_file: Path to the combined (renamed) input GFF file.
        kept_ids: Set of temporary contig IDs to keep (see ``split_fasta``).
        output_gff: Path for the output GFF file.

    Returns:
        Set of CDS protein IDs (original naming) belonging to kept sequences.
    """
    protein_ids: set[str] = set()
    seqs_written: set[str] = set()
    with open(gff_file) as gff_in, open(output_gff, 'w') as gff_out:
        gff_out.write("##gff-version 3\n")
        for line in gff_in:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue

            seq_id = parts[0]
            if seq_id not in kept_ids:
                continue

            gff_out.write(line if line.endswith('\n') else line + '\n')
            if parts[2] == 'CDS':
                attrs, _ = parse_attributes(parts[8])
                protein_id = attrs.get('ID')
                if protein_id:
                    protein_ids.add(protein_id)
            else:
                seqs_written.add(seq_id)
    print(f"Wrote {len(seqs_written)} sequences to {output_gff}")
    return protein_ids


def split_faa(faa_file: str, protein_ids: set[str], output_faa: str) -> None:
    """Write protein records whose ID is in ``protein_ids`` to ``output_faa``.

    Args:
        faa_file: Path to the combined protein FASTA file.
        protein_ids: Set of protein IDs to keep (see ``split_gff``).
        output_faa: Path for the output protein FASTA file.
    """
    count = 0
    with open(output_faa, 'w') as out_f:
        for record in SeqIO.parse(faa_file, "fasta"):
            if record.id in protein_ids:
                SeqIO.write(record, out_f, "fasta")
                count += 1
    print(f"Wrote {count} proteins to {output_faa}")


def main() -> None:
    args = parse_args()
    definitions = read_definitions(args.map)

    output_fna = f"{args.output_prefix}.fna"
    kept_ids = split_fasta(args.fna, definitions, args.category, output_fna)

    if args.gff:
        output_gff = f"{args.output_prefix}.gff"
        protein_ids = split_gff(args.gff, kept_ids, output_gff)

        if args.faa:
            output_faa = f"{args.output_prefix}.faa"
            split_faa(args.faa, protein_ids, output_faa)


if __name__ == "__main__":
    main()
