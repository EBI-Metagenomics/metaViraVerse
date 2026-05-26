#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path

_table_headers = [
    "target name",
    "target accession",
    "tlen",
    "query name",
    "query accession",
    "qlen",
    "full sequence E-value",
    "full sequence score",
    "full sequence bias",
    "#",
    "of",
    "c-Evalue",
    "i-Evalue",
    "domain score",
    "domain bias",
    "hmm coord from",
    "hmm coord to",
    "ali coord from",
    "ali coord to",
    "env coord from",
    "env coord to",
    "acc",
    "description of target"
]


def open_file(path):
    """Open a plain or gzip-compressed file for reading as text."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Format hmmer domain hits table.")
    parser.add_argument("-i", dest="input_table",
                        help="hmmer domain hits table")
    parser.add_argument("-o", "--outname",
                        dest="outfile_name", help="Output table name")
    parser.add_argument("-t", "--hmmer-tool",
                        dest="hmmer_tool", help="what hmm tool was used to generate input table", choices=['hmmsearch', 'hmmscan'])
    args = parser.parse_args()

    # hmmscan table format specified in http://eddylab.org/software/hmmer/Userguide.pdf
    # basically:
    #   The domain table has 22 whitespace-delimited fields
    #   followed by a free text target sequence description, as follows
    domain_table = Path(args.input_table)
    if not domain_table.is_file():
        raise Exception(
            "Input domain hits table missing. Path: " + args.input_table)

    with open(args.outfile_name + ".tsv", "w", newline="") as out_table:
        tsv_writer = csv.writer(out_table, delimiter="\t",
                                quoting=csv.QUOTE_MINIMAL)
        tsv_writer.writerow(_table_headers)
        with open_file(domain_table) as dt_reader:
            for line in dt_reader:
                if line.startswith("#"):
                    continue
                # the last column has free text, so get it by itself to
                # properly espace it
                cols = line.split()
                data = cols[:22]
                if args.hmmer_tool == 'hmmsearch':
                    data[0:3], data[3:6] = data[3:6], data[0:3]
                description = " ".join(cols[22:])
                tsv_writer.writerow(data + [description])
