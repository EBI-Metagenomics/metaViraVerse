#!/usr/bin/env python

import argparse
import re
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Define best host from iphop.csv."
    )
    parser.add_argument("-i", "--input", help="Input CSV from iphop")
    parser.add_argument("-o", "--output", default="best.csv", help="Output ")
    parser.add_argument(
        "--itol", default=None,
        help="Write a Krona-style count TSV compatible with plot_itol.py "
             "(count<TAB>tax1<TAB>tax2<TAB>...)",
    )
    return parser.parse_args()


def get_best(input_file):
    best_hits = {}
    with open(input_file, 'r') as file_in:
        for line in file_in:
            if 'Virus' in line:
                continue
            line = line.strip().split(',')
            virus = line[0]
            additional_score = (float(line[5].split(';')[-1]) if len(line[5].split(';')) > 1 else None)
            metadata = {'taxonomy': line[2], 'score': float(line[4]), 'additional': additional_score}
            if virus in best_hits:
                if best_hits[virus]['score'] > float(line[4]):
                    continue
                elif best_hits[virus]['score'] == float(line[4]):
                    if additional_score and best_hits[virus]['additional']:
                        if best_hits[virus]['additional'] > additional_score:
                            continue
                else:
                    best_hits[virus] = metadata
            else:
                best_hits[virus] = metadata
    return best_hits


def generate_itol_input(best_hits, output_file):
    """Write a Krona-style count TSV for use as input to plot_itol.py.

    Each row: count<TAB>tax1<TAB>tax2<TAB>...
    GTDB-style rank prefixes (d__, p__, c__, o__, f__, g__, s__) are stripped.
    Taxonomies are kept in their original forward order (domain → genus).
    """
    tax_counts = Counter(
        v["taxonomy"] for v in best_hits.values() if v.get("taxonomy")
    )
    with open(output_file, "w") as fh:
        for tax, count in sorted(tax_counts.items()):
            parts = [
                re.sub(r"^[a-z]__", "", p.strip())
                for p in tax.split(";")
                if p.strip()
            ]
            if parts:
                fh.write("\t".join([str(count)] + parts[:4]) + "\n")


def main():
    args = parse_args()
    best_hits = get_best(args.input)
    with open(args.output, 'w') as file_out:
        for i in best_hits:
            file_out.write('\t'.join([i, best_hits[i]['taxonomy'], str(best_hits[i]['score'])]) + '\n')
    if args.itol:
        generate_itol_input(best_hits, args.itol)

if __name__ == '__main__':
    main()
