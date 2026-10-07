#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from copy import copy
import os
import re
import gzip
import bz2
from pathlib import Path
from Bio import SeqIO


NUM_MGYV_DIGITS = 10


def input_args():
    """Multi fasta rename"""
    parser = argparse.ArgumentParser(
        description="Rename multi fasta"
    )
    parser.add_argument(
        "-f", "--fasta", help="MGnify input FASTA file(s), plain or gzip/bzip2-compressed", required=False, nargs='+', default=[]
    )
    parser.add_argument(
        "-g", "--gff", help="Input GFF file(s), one per FASTA file, plain or gzip/bzip2-compressed", required=False, nargs='+'
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
        help="Third-party input FASTA file(s), plain or gzip/bzip2-compressed, renamed after --fasta, continuing the same accession count"
    )
    parser.add_argument(
        "--gff-tp",
        dest="gff_tp",
        required=False,
        nargs='+',
        help="Third-party input GFF file(s), one per --fasta-tp file, plain or gzip/bzip2-compressed"
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


def open_file(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    if path.suffix == ".bz2":
        return bz2.open(path, "rt")
    return open(path, "rt")


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


def define_map_key(s):
    # replace all special symbols to _
    s = re.sub(r"[^a-zA-Z0-9_]", "_", s)
    return s


def rename_fasta(
    input_fasta: str,
    prefix: str,
    start_accession: int,
    biome: str = 'NA',
    source_type: str = 'NA',
    definition: str | None = None,
    viral_sequence_identifier: str | None = None,
    prophage_identifier: str | None = None,
    plasmid_identifier: str | None = None,
    third_party: bool = False
) -> tuple[list, dict[str, str], list, int]:
    """Rename FASTA entries with <prefix><counter> and build a mapping table.

    Args:
        input_fasta: Path to the input FASTA file (may be gzip- or bzip2-compressed).
        prefix: Prefix for new accession names (e.g. 'MGYV', 'seq').
        start_accession: Starting counter value for accession numbering.
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
    print("Processing " + input_fasta)
    with open_file(input_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            temporary_name = define_prefix(prefix, start_accession + count)
            name = record.description

            new_record = copy(record)
            new_record.id = temporary_name
            new_record.description = temporary_name
            fasta_records.append(new_record)

            # populate mapping
            short_name = name.split(' ')[0]
            record_definition = definition if definition is not None else classify_definition(
                name, viral_sequence_identifier, prophage_identifier, plasmid_identifier
            )
            if third_party:
                record_definition = f'third_party_{record_definition}'
            map_records.append([name, temporary_name, short_name, biome, source_type, record_definition])
            # Pyrodigal (used for third-party gene calling) writes only the short
            # FASTA id (record.id, the token before the first whitespace) as the GFF
            # seqid, not the full description -- so third-party lookups must key on
            # record.id to match. MGnify GFFs are resolved via define_catalogue_id's
            # pattern reconstruction instead, which needs the full description.
            map_key = define_map_key(record.id if third_party else name)
            map_dir[map_key] = temporary_name
            count += 1
    print(f"Assigned accessions with {prefix}: {start_accession} - {start_accession+count}")
    return fasta_records, map_dir, map_records, count


def define_catalogue_id(feat_id, contig_id):
    pattern = r"^(MGYG\d+_\d+)\|([\w.]+)-(\d+):(\d+)$"
    pattern_old_map_viral = r"^MGYG\d+_\d+\|viral_sequence$"
    pattern_old_map_plasmid = r"^plasmid_\d+$"

    id_fna = feat_id
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
        id_fna = contig_id
    return id_fna


def define_gff_map_key(parts, third_party, current_id):
    proteins = 0
    if parts[2] == 'CDS':
        proteins += 1
    if third_party:
        """
        We predict proteins with pyrodigal for third party data.
        In pyrodigal output there is no GFF line with sequence description.
        example,
        # # Sequence Data: seqnum=1;seqlen=5757;seqhdr="MGYG000517684_26|plasmid-1:6000"
        # Model Data: version=pyrodigal.v3.7.1;run_type=Single;model="Ab initio";gc_cont=63.64;transl_table=11;uses_sd=1
        MGYG000517684_26|plasmid-1:6000 pyrodigal_v3.7.1        CDS     3       1922    264.5   -       0       ID=MGYG000517684_26|plasmid-1:6000_1;...
        """
        return proteins, define_map_key(parts[0])
    else:
        if parts[2] == 'CDS':
            return proteins, current_id
        else:
            attrs, _ = parse_attrs(parts[8])
            id = attrs.get('ID')
            id_checked = define_catalogue_id(id, parts[0])
            return proteins, define_map_key(id_checked)


def rename_gff(
    input_gff: str,
    map_dir: dict[str, str],
    third_party: bool = False
) -> int:
    """Rename GFF feature IDs using the temporary-name mapping.

    Reads ``input_gff`` (plain or compressed), replaces feature IDs that appear
    in ``map_dir`` with their new accessions.

    Args:
        input_gff: Path to the input GFF file (may be gzip- or bzip2-compressed).
        map_dir: Dict of original_name -> temporary_name.
        third_party: Boolean flag defining input was third_party or not

    Returns:
        Number of feature lines written.
    """
    count = 0
    gff_records = []
    count_proteins = 0
    seqs_set = set()
    print(f"Processing {input_gff} ")
    id_fna = ''
    with open_file(input_gff) as file_in:
        for line in file_in:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) == 9:
                protein, id_fna = define_gff_map_key(parts, third_party, id_fna)
                count_proteins += protein
                if id_fna in map_dir:
                    new_parts = [map_dir[id_fna]] + parts[1:]
                    if not protein:
                        # non-CDS (contig-level) record: point its ID= attribute at
                        # the new accession too. CDS lines keep their protein ID.
                        new_parts[8] = re.sub(
                            r'(?<=ID=)[^;]*', map_dir[id_fna], new_parts[8], count=1
                        )
                    full_line = '\t'.join(new_parts) + '\n'
                else:
                    print(f'Warning: no {id_fna} found in mapping')
                    full_line = '\t'.join(parts) + '\n'
                count += 1
                seqs_set.add(id_fna)
                gff_records.append(full_line)
            else:
                print(f'Adding {line} to gff')
                gff_records.append(line)
    seqs_count = len(list(seqs_set))
    print(f"Found {seqs_count} sequences")
    print(f"Found {count_proteins} feature lines")
    return seqs_count, count_proteins, gff_records


def write_fasta(outputname, fasta_records):
    with open(outputname, "w") as fasta_out:
        for record in fasta_records:
            SeqIO.write(record, fasta_out, "fasta")
    print(f"Wrote {len(fasta_records)} sequences to {outputname}")


def write_gff(outputname, records):
    with open(outputname, "w") as file_out:
        file_out.write('##gff-version 3\n')
        file_out.write(''.join(records))


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


def rename_fasta_and_gff(fastas, gffs, biomes, types, prefix, start_accession, third_party,
                           viral_sequence_identifier, prophage_identifier, plasmid_identifier, combine, map, sources):

    all_fasta_records: list = []
    all_gff_records: list = []
    all_map_records: list = []
    gff_records = []

    for fna, gff, biome, source_type, definition in zip(fastas, gffs, biomes, types, sources):
        fasta_records, map_dir, map_records, count_fasta = rename_fasta(
            input_fasta=fna, prefix=prefix, start_accession=start_accession,
            biome=biome, source_type=source_type, definition=definition,
            viral_sequence_identifier=viral_sequence_identifier,
            prophage_identifier=prophage_identifier,
            plasmid_identifier=plasmid_identifier, third_party=third_party
        )
        if gff:
            gff_seqs_count, gff_count_proteins, gff_records = rename_gff(gff, map_dir=map_dir, third_party=third_party)

            if gff_seqs_count != count_fasta:
                print(f'Sanity check failed GFF={gff_seqs_count}, FASTA={count_fasta}. Exit')
                exit(1)

        if combine:
            all_fasta_records.extend(fasta_records)
            all_gff_records.extend(gff_records)
            all_map_records.extend(map_records)
            start_accession += count_fasta
        else:
            basename = os.path.splitext(os.path.basename(fna))[0]
            write_fasta('renamed_' + os.path.basename(fna), fasta_records)
            write_mapfile(basename + '_' + map, map_records)
            if gff:
                output_gff = basename + '_renamed.gff'
                write_gff(output_gff, gff_records)
    return start_accession, all_fasta_records, all_gff_records, all_map_records


def main():
    args = input_args()

    all_fasta_records: list = []
    all_map_records: list = []
    all_gff_records: list = []
    start_accession = args.start

    # rename catalogue input
    if args.fasta:
        biomes = args.biome if args.biome else ['NA'] * len(args.fasta)
        types = args.type if args.type else ['NA'] * len(args.fasta)
        gffs = args.gff if args.gff else [None] * len(args.fasta)
        sources = [None] * len(args.fasta)

        start_accession_upd, all_fasta_records, all_gff_records, all_map_records = \
            rename_fasta_and_gff(fastas=args.fasta, gffs=gffs,
                               biomes=biomes, types=types, prefix=args.prefix, map=args.map,
                               combine=args.combine, third_party=False,
                               start_accession=start_accession,
                               viral_sequence_identifier=args.viral_sequence_identifier,
                               prophage_identifier=args.prophage_identifier,
                               plasmid_identifier=args.plasmid_identifier, sources=sources)
        start_accession = start_accession_upd

    # rename third party input
    if args.fasta_tp:
        biomes_tp = args.biome_tp if args.biome_tp else ['NA'] * len(args.fasta_tp)
        sources_tp = args.source_tp if args.source_tp else ['NA'] * len(args.fasta_tp)
        gffs_tp = args.gff_tp if args.gff_tp else [None] * len(args.fasta_tp)
        types_tp = args.type_tp

        start_accession_tp, all_fasta_records_tp, all_gff_records_tp, all_map_records_tp = \
            rename_fasta_and_gff(fastas=args.fasta_tp, gffs=gffs_tp, sources=types_tp,
                                 biomes=biomes_tp, types=sources_tp, prefix=args.prefix, map=args.map,
                                 combine=args.combine, third_party=True,
                                 start_accession=start_accession,
                                 viral_sequence_identifier=args.viral_sequence_identifier,
                                 prophage_identifier=args.prophage_identifier,
                                 plasmid_identifier=args.plasmid_identifier)
        all_fasta_records.extend(all_fasta_records_tp)
        all_gff_records.extend(all_gff_records_tp)
        all_map_records.extend(all_map_records_tp)

    if args.combine:
        print('Write combined outputs')
        combined_base = os.path.splitext(args.map)[0]
        write_fasta(combined_base + '.fna', all_fasta_records)
        write_mapfile(args.map, all_map_records)
        write_gff(combined_base + '.gff', all_gff_records)


if __name__ == "__main__":
    main()
