#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from copy import copy
import os
import fileinput
import re
from Bio import SeqIO


NUM_MGYV_DIGITS = 10


def input_args():
    """Multi fasta rename"""
    parser = argparse.ArgumentParser(
        description="Rename multi fasta"
    )
    parser.add_argument(
        "-f", "--fasta", help="indicate input FASTA file", required=True, nargs='+'
    )
    parser.add_argument(
        "-g", "--gff", help="Input GFF file(s), one per FASTA file", required=False, nargs='+'
    )
    parser.add_argument(
        "-m", "--map", help="map file for names", required=False, default="map.txt"
    )
    parser.add_argument(
        "-p", "--prefix", help="Prefix that would be included to header <prefix><digit>", required=False, default="seq"
    )
    parser.add_argument(
        "--start", help="First digit for renaming, ex. prefix1", required=False, default=1, type=int
    )
    parser.add_argument(
        "--end", help="Last digit for renaming, ex. prefix100. Can be skipped", required=False, type=int
    )
    parser.add_argument(
        "--combine", help="Combine all fasta files into one", action='store_true'
    )
    parser.add_argument(
        "-k", "--keep-viral-identifier", help="Keep info after | in new name", action='store_true'
    )
    parser.add_argument(
        "--biome",
        required=False,
        nargs='+',
        help="Biome label for each input FASTA file (same order as --fasta)"
    )
    parser.add_argument(
        "--type",
        required=False,
        nargs='+',
        help="Source type for each input FASTA file (same order as --fasta)"
    )
    args = parser.parse_args()
    if args.biome and len(args.biome) != len(args.fasta):
        parser.error("--biome must have the same number of values as --fasta")
    if args.type and len(args.type) != len(args.fasta):
        parser.error("--type must have the same number of values as --fasta")
    if args.gff and len(args.gff) != len(args.fasta):
        parser.error("--gff must have the same number of values as --fasta")
    return args


def parse_attrs(attrs_str: str) -> tuple[dict[str, str], list[str]]:
    """Parse a GFF3 column-9 attributes string into a dict and an ordered key list.

    :param attrs_str: Semicolon-separated key=value attribute string from GFF column 9.
    :return: Tuple of (attrs dict, list of keys in original order).
    """
    attrs, order = {}, []
    for part in attrs_str.rstrip(";").split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            if k not in attrs:
                order.append(k)
            attrs[k] = v
    return attrs, order


def define_prefix(prefix, num):
    if prefix == "MGYV":
        accession = f"MGYV{num:0{NUM_MGYV_DIGITS}d}"
    else:
        accession = f"{prefix}{num}"
    return accession


def rename_fasta(
    input_fasta: str,
    prefix: str,
    keep_viral_identifier: bool,
    start_accession: int,
    end_accession: int | None,
    biome: str = 'NA',
    source_type: str = 'NA',
) -> tuple[list, dict[str, str], list, int]:
    """Rename FASTA entries with <prefix><counter> and build a mapping table.

    Args:
        input_fasta: Path to the input FASTA file.
        prefix: Prefix for new accession names (e.g. 'MGYV', 'seq').
        keep_viral_identifier: When True, append the ``|<identifier>`` suffix
            from the original name to the new accession.
        start_accession: Starting counter value for accession numbering.
        end_accession: Expected final counter value (used only for logging).
        biome: Biome label to record in the map file for every sequence.
        source_type: Source type label to record in the map file for every sequence.

    Returns:
        Tuple of (fasta_records, map_dir, map_records, count) where:
        - fasta_records: list of renamed SeqRecord objects.
        - map_dir: dict of original_name -> temporary_name (for GFF renaming).
        - map_records: list of rows for the TSV map file.
        - count: number of sequences processed.
    """
    fasta_records = []
    map_records = []
    map_dir = {}
    count = 0
    print("Renaming " + input_fasta)
    for record in SeqIO.parse(str(input_fasta), "fasta"):
        temporary_name = define_prefix(prefix, start_accession + count)
        name = record.id
        viral_identifier = None
        if keep_viral_identifier:
            if '|' in name:
                viral_identifier = name.split('|')[1]
                temporary_name = temporary_name + '|' + viral_identifier
        new_record = copy(record)
        new_record.id = temporary_name
        new_record.description = temporary_name
        fasta_records.append(new_record)
        # populate mapping
        short_name = name.split(' ')[0]
        map_records.append([name, temporary_name, short_name, biome, source_type])
        map_key = name.replace('|' + viral_identifier, '') if viral_identifier else name
        map_dir[map_key] = temporary_name
        count += 1
    print(f"Assigned accessions with {prefix}: {start_accession} - {start_accession+count}")
    if end_accession:
        print(f"Expected end_accession specified: {end_accession}")
    return fasta_records, map_dir, map_records, count


def rename_gff(
    input_gff: str,
    output_gff: str,
    map_dir: dict[str, str],
    append: bool = False,
) -> int:
    """Rename GFF feature IDs using the temporary-name mapping.

    Reads ``input_gff`` (plain or compressed), replaces feature IDs that appear
    in ``map_dir`` with their new accessions, and writes to ``output_gff``.

    Args:
        input_gff: Path to the input GFF file (may be gzip-compressed).
        output_gff: Path for the output GFF file.
        map_dir: Dict of original_name -> temporary_name.
        append: When True, open ``output_gff`` in append mode and skip
            ``##``-prefixed header lines (used when merging multiple GFFs
            into one combined output).

    Returns:
        Number of feature lines written.
    """
    pattern = r"^(MGYG\d+_\d+)\|([\w.]+)-(\d+):(\d+)$"
    mode = 'a' if append else 'w'
    count = 0
    seqs_count = 0
    print(f"Renaming {input_gff} -> {output_gff}")
    with fileinput.hook_compressed(input_gff, "r") as file_in, open(output_gff, mode) as file_out:
        for line in file_in:
            if '##' in line:
                if not append:
                    file_out.write(line)
                continue
            parts = line.strip().split('\t')
            if len(parts) == 9:
                if parts[2] == "CDS":
                    file_out.write('\t'.join(parts) + '\n')
                    count += 1
                    continue
                attrs, _ = parse_attrs(parts[8])
                feat_id = attrs.get("ID", "")
                id_fna = feat_id
                full_line = '\t'.join(parts)
                if re.fullmatch(pattern, feat_id):
                    # handle case where fna has >MGYG000535629_9 viral_sequence|1:3862
                    # but gff has MGYG000535629_9|viral_sequence-1:3862
                    m = re.fullmatch(pattern, feat_id)
                    id_fna = f"{m.group(1)} {m.group(2)}|{m.group(3)}:{m.group(4)}"
                if id_fna in map_dir:
                    file_out.write(full_line.replace(feat_id, map_dir[id_fna]) + '\n')
                else:
                    file_out.write(full_line + '\n')
                count += 1
                seqs_count += 1
            else:
                file_out.write('\t'.join(parts) + '\n')
    print(f"Wrote {seqs_count} sequences to {output_gff}")
    print(f"Wrote {count} feature lines to {output_gff}")
    return count, seqs_count


def write_fasta(outputname, fasta_records):
    with open(outputname, "w") as fasta_out:
        for record in fasta_records:
            SeqIO.write(record, fasta_out, "fasta")
    print(f"Wrote {len(fasta_records)} sequences to {outputname}")


def write_mapfile(mapfilename: str, map_data: list) -> None:
    """Write the contig mapping table to a TSV file.

    Args:
        mapfilename: Output path for the TSV map file.
        map_data: List of rows, each containing
            [original, temporary, short, biome, type].
    """
    with open(mapfilename, "w") as map_tsv:
        tsv_map = csv.writer(map_tsv, delimiter="\t")
        tsv_map.writerow(["original", "temporary", "short", "biome", "type"])
        for line in map_data:
            tsv_map.writerow(line)


def main():
    args = input_args()
    start_accession = args.start
    end_accession = args.end

    all_fasta_records: list = []
    all_map_records: list = []
    combined_map_dir: dict[str, str] = {}
    count_fasta = 0

    biomes = args.biome if args.biome else ['NA'] * len(args.fasta)
    types  = args.type  if args.type  else ['NA'] * len(args.fasta)
    gffs   = args.gff   if args.gff   else [None]  * len(args.fasta)

    for fna, gff, biome, source_type in zip(args.fasta, gffs, biomes, types):
        fasta_records, map_dir, map_records, count_fasta = rename_fasta(
            fna, args.prefix, args.keep_viral_identifier, start_accession, end_accession,
            biome=biome, source_type=source_type,
        )
        if not args.combine:
            basename = os.path.splitext(os.path.basename(fna))[0]
            write_fasta('renamed_' + os.path.basename(fna), fasta_records)
            write_mapfile(basename + '_' + args.map, map_records)
            if gff:
                output_gff = basename + '_renamed.gff'
                count_gff = rename_gff(gff, output_gff, map_dir)
                if count_gff != count_fasta:
                    print(f"Sanity check failed for {gff}: GFF={count_gff}, FASTA={count_fasta}. Exit")
                    exit(1)
                print("Sanity check passed")
        else:
            all_fasta_records.extend(fasta_records)
            all_map_records.extend(map_records)
            combined_map_dir.update(map_dir)
            start_accession += count_fasta

    if args.combine:
        combined_base = os.path.splitext(args.map)[0]
        write_fasta(combined_base + '.fna', all_fasta_records)
        write_mapfile(args.map, all_map_records)

        if args.gff:
            combined_gff = combined_base + '.gff'
            total_gff = 0
            for i, gff in enumerate(args.gff):
                total_gff += rename_gff(gff, combined_gff, combined_map_dir, append=(i > 0))
            if total_gff != len(all_fasta_records):
                print(f"Sanity check failed: GFF={total_gff}, FASTA={len(all_fasta_records)}. Exit")
                exit(1)
            print("Sanity check passed")


if __name__ == "__main__":
    main()
