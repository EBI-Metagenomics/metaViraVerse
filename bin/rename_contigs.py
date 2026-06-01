#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import fileinput
import re


NUM_MGYV_DIGITS = 10


def input_args():
    """Multi fasta rename"""
    parser = argparse.ArgumentParser(
        description="Rename multi fasta"
    )
    parser.add_argument(
        "-f", "--input", help="indicate input FASTA file", required=True
    )
    parser.add_argument(
        "-g", "--gff", help="indicate input GFF file", required=False
    )
    parser.add_argument(
        "-m", "--map", help="map file for names", required=False, default="map.txt"
    )
    parser.add_argument(
        "-p", "--prefix", help="Prefix that would be included to header <prefix><digit>", required=False, default="seq"
    )
    parser.add_argument(
        "-k", "--keep-viral-identifier", help="Keep info after | in new name", action='store_true'
    )
    parser.add_argument(
        "--start", help="First digit for renaming, ex. prefix1", required=False, default=0, type=int
    )
    parser.add_argument(
        "--end", help="Last digit for renaming, ex. prefix100. Can be skipped", required=False, type=int
    )
    args = parser.parse_args()
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


def rename_fasta(input_fasta, mapfilename, prefix, keep_viral_identifier, start_accession, end_accession):
    """Rename a multi-fasta fasta entries with <name>.<counter> and store the
    mapping between new and old files in tsv
    """
    name = os.path.basename(input_fasta).split('.')[0]
    ext = 'fasta'
    print(ext)
    output = name + '_renamed.' + ext
    map_dir = {}
    print("Renaming " + input_fasta)
    with fileinput.hook_compressed(input_fasta, "r") as fasta_in:
        with open(output, "w") as fasta_out, open(mapfilename, "w") as map_tsv:
            count = 0
            tsv_map = csv.writer(map_tsv, delimiter="\t")
            tsv_map.writerow(["original", "temporary", "short"])
            for line in fasta_in:
                line = str(line)
                if line.startswith(">"):
                    count += 1
                    temporary_name = define_prefix(prefix, start_accession + count)
                    name = line.strip().replace('>', '')
                    viral_identifier = None
                    if keep_viral_identifier:
                        if '|' in name:
                            viral_identifier = name.split('|')[1]
                            temporary_name = temporary_name + '|' + viral_identifier
                    short_name = name.split(' ')[0]
                    fasta_out.write(f">{temporary_name}\n")
                    tsv_map.writerow([name, temporary_name, short_name])
                    map_key = name.replace('|' + viral_identifier, '') if viral_identifier else name
                    map_dir[map_key] = temporary_name
                else:
                    fasta_out.write(line)
    print(f"Wrote {count} sequences to {output}")
    print(f"Assigned accessions: {start_accession} - {start_accession+count}")
    if end_accession:
        print(f"Expected end_accession specified: {end_accession}")
    return map_dir, count


def rename_gff(input_gff, map_dir):
    """Rename a multi-fasta fasta entries with <name>.<counter> and store the
    mapping between new and old files in tsv
    """
    name = os.path.basename(input_gff).split('.')[0]
    ext = 'gff'
    output = name + '_renamed.' + ext
    pattern = r"^(MGYG\d+_\d+)\|([\w.]+)-(\d+):(\d+)$"

    print("Renaming " + input_gff)
    with fileinput.hook_compressed(input_gff, "r") as file_in:
        with open(output, "w") as file_out:
            count = 0
            for line in file_in:
                if '##' in line:
                    file_out.write(line)
                    count += 1
                    continue
                line = line.strip().split('\t')
                if len(line) == 9:
                    if line[2] == "CDS":
                        file_out.write('\t'.join(line) + '\n')
                        continue
                    else:
                        attrs, _ = parse_attrs(line[8])
                        id = attrs.get("ID", "")
                        id_fna = id
                        full_line = '\t'.join(line)
                        if re.fullmatch(pattern, id):
                            # check case when fna has >MGYG000535629_9 viral_sequence|1:3862
                            # but gff has MGYG000535629_9|viral_sequence-1:3862
                            m = re.fullmatch(pattern, id)
                            id_fna = f"{m.group(1)} {m.group(2)}|{m.group(3)}:{m.group(4)}"
                        if id_fna in map_dir:
                            new_line = full_line.replace(id, map_dir[id_fna])
                            file_out.write(new_line + '\n')
                        else:
                            file_out.write(full_line + '\n')
                        count += 1
                else:
                    file_out.write('\t'.join(line) + '\n')
                    continue
    print(f"Renamed {count} sequences to {output}")
    return count


def main():
    args = input_args()
    map_dir, count_fasta = rename_fasta(args.input, args.map, args.prefix, args.keep_viral_identifier, args.start, args.end)
    if args.gff:
        count_gff = rename_gff(args.gff, map_dir)
        if count_gff == count_fasta:
            print("Sanity check passed")
        else:
            print("Number of renamed records doesn't match. Exit")
            exit(1)


if __name__ == "__main__":
    main()
