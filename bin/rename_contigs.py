#!/usr/bin/env python

import argparse
import csv
import os
import fileinput


def input_args():
    """Multi fasta rename"""
    parser = argparse.ArgumentParser(
        description="Rename multi fasta"
    )
    parser.add_argument(
        "-f", "--input", help="indicate input FASTA file", required=True
    )
    parser.add_argument(
        "-g", "--gff", help="indicate input GFF file", required=True
    )
    parser.add_argument(
        "-m", "--map", help="map file for names", required=False, default="map.txt"
    )
    parser.add_argument(
        "-p", "--prefix", help="Prefix that would be included to header <prefix><digit>", required=False
    )
    args = parser.parse_args()
    return args


def rename_fasta(input_fasta, mapfilename, prefix):
    """Rename a multi-fasta fasta entries with <name>.<counter> and store the
    mapping between new and old files in tsv
    """
    name = os.path.basename(input_fasta).split('.')[0]
    ext = os.path.basename(input_fasta).split('.')[1]
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
                    temporary_name = f"{prefix}_{count}"
                    name = line.strip().replace('>', '')
                    if '|' in name:
                        viral_identifier = name.split('|')[1]
                        temporary_name = temporary_name + '|' + viral_identifier
                    short_name = name.split(' ')[0]
                    fasta_out.write(f">{temporary_name}\n")
                    tsv_map.writerow([name, temporary_name, short_name])
                    map_dir[name.replace('|' + viral_identifier, '')] = temporary_name
                else:
                    fasta_out.write(line)
    print(f"Wrote {count} sequences to {output}.")
    return map_dir


def rename_gff(input_gff, map_dir):
    """Rename a multi-fasta fasta entries with <name>.<counter> and store the
    mapping between new and old files in tsv
    """
    name = os.path.basename(input_gff).split('.')[0]
    ext = os.path.basename(input_gff).split('.')[1]
    output = name + '_renamed.' + ext

    print("Renaming " + input_gff)
    with fileinput.hook_compressed(input_gff, "r") as file_in:
        with open(output, "w") as file_out:
            count = 0
            for line in file_in:
                if '##' in line:
                    file_out.write(line)
                    continue
                line = str(line).strip().split('\t')
                name = line[0]
                data = '\t'.join(line[1:])
                temporary_name = map_dir[name]
                file_out.write(f"{temporary_name}\t{data}\n")
    print(f"Wrote {count} sequences to {output}.")


def main():
    args = input_args()
    map_dir = rename_fasta(args.input, args.map, args.prefix)
    rename_gff(args.gff, map_dir)


if __name__ == "__main__":
    main()
