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
        "-f", "--fasta", help="MGnify input FASTA file(s)", required=False, nargs='+', default=[]
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
        help="Source type for each input FASTA file (same order as --fasta): 'genome' or 'metagenome'"
    )
    parser.add_argument(
        "--viral-sequence-identifier",
        dest="viral_sequence_identifier",
        required=False,
        help="Regex identifying free viral sequences by their original name (MGnify input), e.g. 'viral_sequence'"
    )
    parser.add_argument(
        "--prophage-identifier",
        dest="prophage_identifier",
        required=False,
        help="Regex identifying prophages by their original name (MGnify input), e.g. 'prophage'"
    )
    parser.add_argument(
        "--plasmid-identifier",
        dest="plasmid_identifier",
        required=False,
        help="Regex identifying plasmids by their original name (MGnify input), e.g. 'plasmid'"
    )
    parser.add_argument(
        "--fasta-tp",
        dest="fasta_tp",
        required=False,
        nargs='+',
        default=[],
        help="Third-party input FASTA file(s), renamed after --fasta, continuing the same accession count"
    )
    parser.add_argument(
        "--gff-tp",
        dest="gff_tp",
        required=False,
        nargs='+',
        help="Third-party input GFF file(s), one per --fasta-tp file"
    )
    parser.add_argument(
        "--biome-tp",
        dest="biome_tp",
        required=False,
        nargs='+',
        help="Biome label for each --fasta-tp file (same order as --fasta-tp)"
    )
    parser.add_argument(
        "--source-tp",
        dest="source_tp",
        required=False,
        nargs='+',
        help="Source label for each --fasta-tp file (same order as --fasta-tp), from the third-party "
             "samplesheet's 'source' column: 'genome' or 'metagenome'"
    )
    parser.add_argument(
        "--type-tp",
        dest="type_tp",
        required=False,
        nargs='+',
        choices=['virus', 'plasmid', 'prophage'],
        help="Sequence type for each --fasta-tp file: 'virus', 'plasmid' or 'prophage' (same order as --fasta-tp)"
    )
    args = parser.parse_args()
    if args.biome and len(args.biome) != len(args.fasta):
        parser.error("--biome must have the same number of values as --fasta")
    if args.type and len(args.type) != len(args.fasta):
        parser.error("--type must have the same number of values as --fasta")
    if args.gff and len(args.gff) != len(args.fasta):
        parser.error("--gff must have the same number of values as --fasta")

    if not args.fasta and not args.fasta_tp:
        parser.error("At least one of --fasta or --fasta-tp must be provided")

    if args.fasta_tp:
        if not args.type_tp or len(args.type_tp) != len(args.fasta_tp):
            parser.error("--type-tp is required and must have the same number of values as --fasta-tp")
        if args.biome_tp and len(args.biome_tp) != len(args.fasta_tp):
            parser.error("--biome-tp must have the same number of values as --fasta-tp")
        if args.source_tp and len(args.source_tp) != len(args.fasta_tp):
            parser.error("--source-tp must have the same number of values as --fasta-tp")
        if args.gff_tp and len(args.gff_tp) != len(args.fasta_tp):
            parser.error("--gff-tp must have the same number of values as --fasta-tp")
        if not args.combine:
            parser.error("--fasta-tp is only supported together with --combine")
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


def classify_definition(
    name: str,
    viral_sequence_identifier: str | None,
    prophage_identifier: str | None,
    plasmid_identifier: str | None,
) -> str:
    """Classify an MGnify sequence as virus/prophage/plasmid from its original name.

    Mirrors the pattern matching performed later by ``separate_sequences.py``:
    each identifier is a regex tested (in this priority order) against the
    sequence's original description, and the first pattern that matches wins.

    Args:
        name: Original sequence description (before renaming).
        viral_sequence_identifier: Regex identifying free viral sequences (e.g. 'viral_sequence').
        prophage_identifier: Regex identifying prophages (e.g. 'prophage').
        plasmid_identifier: Regex identifying plasmids (e.g. 'plasmid').

    Returns:
        'virus', 'prophage', 'plasmid', or 'NA' if no pattern is given or none match.
    """
    if viral_sequence_identifier and re.search(viral_sequence_identifier, name):
        return 'virus'
    if prophage_identifier and re.search(prophage_identifier, name):
        return 'prophage'
    if plasmid_identifier and re.search(plasmid_identifier, name):
        return 'plasmid'
    return 'NA'


def rename_fasta(
    input_fasta: str,
    prefix: str,
    keep_viral_identifier: bool,
    start_accession: int,
    end_accession: int | None,
    biome: str = 'NA',
    source_type: str = 'NA',
    definition: str | None = None,
    viral_sequence_identifier: str | None = None,
    prophage_identifier: str | None = None,
    plasmid_identifier: str | None = None,
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
        source_type: Source type label to record in the map file for every sequence
            (MGnify: 'genome'/'metagenome'; third-party: the given --source-tp value).
        definition: Fixed 'virus'/'prophage'/'plasmid' classification to record for
            every sequence from this file (used for third-party input, where the type
            is already known). When ``None``, the classification is instead derived
            per-record from ``viral_sequence_identifier``/``prophage_identifier``/
            ``plasmid_identifier`` (used for MGnify input).
        viral_sequence_identifier: Regex identifying free viral sequences by original name.
        prophage_identifier: Regex identifying prophages by original name.
        plasmid_identifier: Regex identifying plasmids by original name.

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
        name = record.description
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
        record_definition = definition if definition is not None else classify_definition(
            name, viral_sequence_identifier, prophage_identifier, plasmid_identifier
        )
        map_records.append([name, temporary_name, short_name, biome, source_type, record_definition])
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
    pattern_old_map_viral = r"^MGYG\d+_\d+\|viral_sequence$"
    pattern_old_map_plasmid = r"^plasmid_\d+$"
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
                """
                NOTE: old version of MAP had field 'from_mge=' that was not changed during rename
                """
                if re.fullmatch(pattern_old_map_viral, feat_id):
                    # handle cases with results from MAP v1
                    # viral records: MGYG000321096_52	VIRify	viral_sequence	1	16087	.	.	.	ID=MGYG000321096_52|viral_sequence;gbkey= ...
                    id_fna = feat_id.split('|')[0]
                if re.fullmatch(pattern_old_map_plasmid, feat_id):
                    # handle cases with results from MAP v1 where plasmids were predicted with PPR-meta
                    # plasmid records are: MGYG000321096_22	PPR-meta	plasmid	1	22608	.	.	.	ID=plasmid_1;gbkey=mobile_element;mobile_element_type=plasmid
                    id_fna = parts[0]
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
    return seqs_count


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
            [original, temporary, short, biome, type, definition].
    """
    with open(mapfilename, "w") as map_tsv:
        tsv_map = csv.writer(map_tsv, delimiter="\t")
        tsv_map.writerow(["original", "temporary", "short", "biome", "type", "definition"])
        for line in map_data:
            tsv_map.writerow(line)


def main():
    args = input_args()
    start_accession = args.start
    end_accession = args.end

    all_fasta_records: list = []
    all_map_records: list = []
    combined_map_dir: dict[str, str] = {}
    all_gff_files: list = []
    count_fasta = 0

    biomes = args.biome if args.biome else ['NA'] * len(args.fasta)
    types  = args.type  if args.type  else ['NA'] * len(args.fasta)
    gffs   = args.gff   if args.gff   else [None]  * len(args.fasta)

    for fna, gff, biome, source_type in zip(args.fasta, gffs, biomes, types):
        fasta_records, map_dir, map_records, count_fasta = rename_fasta(
            fna, args.prefix, args.keep_viral_identifier, start_accession, end_accession,
            biome=biome, source_type=source_type,
            viral_sequence_identifier=args.viral_sequence_identifier,
            prophage_identifier=args.prophage_identifier,
            plasmid_identifier=args.plasmid_identifier,
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
            if gff:
                all_gff_files.append(gff)
            start_accession += count_fasta

    if args.combine:
        if args.fasta_tp:
            biomes_tp  = args.biome_tp  if args.biome_tp  else ['NA'] * len(args.fasta_tp)
            sources_tp = args.source_tp if args.source_tp else ['NA'] * len(args.fasta_tp)
            gffs_tp    = args.gff_tp    if args.gff_tp    else [None]  * len(args.fasta_tp)

            for fna, gff, biome, source, seq_type in zip(args.fasta_tp, gffs_tp, biomes_tp, sources_tp, args.type_tp):
                fasta_records, map_dir, map_records, count_fasta = rename_fasta(
                    fna, args.prefix, args.keep_viral_identifier, start_accession, end_accession,
                    biome=biome, source_type=source, definition=f"third_party_{seq_type}",
                )
                all_fasta_records.extend(fasta_records)
                all_map_records.extend(map_records)
                combined_map_dir.update(map_dir)
                if gff:
                    all_gff_files.append(gff)
                start_accession += count_fasta

        if all_fasta_records:
            combined_base = os.path.splitext(args.map)[0]
            write_fasta(combined_base + '.fna', all_fasta_records)
            write_mapfile(args.map, all_map_records)

            if all_gff_files:
                combined_gff = combined_base + '.gff'
                total_gff = 0
                for i, gff in enumerate(all_gff_files):
                    total_gff += rename_gff(gff, combined_gff, combined_map_dir, append=(i > 0))
                if total_gff != len(all_fasta_records):
                    print(f"Sanity check failed: GFF={total_gff}, FASTA={len(all_fasta_records)}. Exit")
                    exit(1)
                print("Sanity check passed")
        else:
            print("No input provided (neither MGnify nor third-party), nothing to write")


if __name__ == "__main__":
    main()
