#!/usr/bin/env python3
"""
Summarise HMM annotation hits from a HMMER --tblout file.

Counts how many times each query accession (e.g. VOG ID) appears, then
joins with an optional metadata TSV (e.g. VOG.tsv) using its ID column.
Output is sorted by count descending.

Output columns: accession, count, <all metadata columns in order>

Usage:
  %(prog)s --input ann.tbl --metadata VOG.tsv --output summary.tsv
  %(prog)s --input ann.tbl.gz --metadata VOG.tsv --output summary.tsv
"""

import argparse
import csv
import gzip
import hashlib
import sys
from collections import Counter


def open_file(path):
    """Open a plain or gzip-compressed file for reading as text."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_tblout(path):
    """
    Parse a HMMER --tblout file and yield query names (col index 2).

    Lines starting with '#' are skipped. Fields are whitespace-separated;
    the first 18 columns are fixed, everything after is the free-text
    description of the target.
    """
    with open_file(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 3:
                continue
            yield cols[2]  # query name = HMM profile accession


def load_metadata(path, id_col):
    """
    Load a metadata TSV into a dict keyed by id_col value.

    Returns (meta dict, ordered list of non-key column names).
    """
    meta = {}
    other_cols = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        other_cols = [c for c in reader.fieldnames if c != id_col]
        for row in reader:
            key = row.get(id_col, "").strip()
            if key:
                meta[key] = {c: row.get(c, "").strip() for c in other_cols}
    return meta, other_cols


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarise HMM annotation hits with optional metadata join.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-i", "--input", required=True,
                        help="HMMER --tblout file (plain or .gz).")
    parser.add_argument("-m", "--metadata", default=None,
                        help="Metadata TSV to join on --meta-id-col (e.g. VOG.tsv).")
    parser.add_argument("--meta-id-col", default="ID",
                        help="Column in --metadata used as join key (default: ID).")
    parser.add_argument("-o", "--output", default="-",
                        help="Output TSV path (default: stdout).")
    parser.add_argument("--compress", action="store_true",
                        help="Gzip-compress the output file and write an md5 checksum file.")
    return parser.parse_args()


def main():
    args = parse_args()

    counts = Counter(parse_tblout(args.input))

    meta, meta_cols = {}, []
    if args.metadata:
        meta, meta_cols = load_metadata(args.metadata, args.meta_id_col)

    out_path = args.output
    if args.compress and out_path != "-" and not out_path.endswith(".gz"):
        out_path = out_path + ".gz"

    if out_path == "-":
        out = sys.stdout
    elif out_path.endswith(".gz"):
        out = gzip.open(out_path, "wt", newline="")
    else:
        out = open(out_path, "w", newline="")

    try:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["accession", "count"] + meta_cols)
        for accession, count in counts.most_common():
            row_meta = meta.get(accession, {})
            writer.writerow([accession, count] + [row_meta.get(c, "") for c in meta_cols])
    finally:
        if out_path != "-":
            out.close()

    if args.compress and out_path != "-":
        md5 = hashlib.md5(open(out_path, "rb").read()).hexdigest()
        md5_path = out_path + ".md5"
        with open(md5_path, "w") as f:
            f.write(f"{md5}  {out_path}\n")
        print(f"md5 written to {md5_path}", file=sys.stderr)

    print(
        f"Done. {len(counts)} unique accessions, "
        f"{sum(counts.values())} total hits.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
