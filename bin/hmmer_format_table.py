#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


def open_file(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def parse_hmmer_header(header_lines):
    """Derive column names and start positions from the HMMER comment header.

    Column boundaries are taken from the START of each dash run in the
    separator line.  Using starts (not ends) is required because HMMER
    right-aligns numbers and values can extend one character past the dash
    span into the gap before the next column.

    Returns:
        col_names  : list of str exactly as written in the header name line
        col_starts : list of int start positions; last column runs to end-of-line
    """
    name_line = dash_line = None
    for raw in header_lines:
        body = raw[1:]  # strip leading '#'
        if "target name" in body:
            name_line = body
        elif re.fullmatch(r"[-\s]+", body.strip()):
            dash_line = body

    if dash_line is None or name_line is None:
        raise ValueError("Could not locate column header lines in the HMMER file.")

    col_starts = [m.start() for m in re.finditer(r"-+", dash_line)]

    col_names = []
    for i, start in enumerate(col_starts):
        end = col_starts[i + 1] if i + 1 < len(col_starts) else None
        chunk = (name_line[start:end] if end is not None else name_line[start:]).strip()
        col_names.append(chunk)

    # Prefix E-value, score, bias by occurrence order:
    # first occurrence → "full <name>", second → "domain <name>".
    _to_prefix = {"E-value", "score", "bias"}
    _prefixes = {1: "full", 2: "domain"}
    _counts = defaultdict(int)
    for i, name in enumerate(col_names):
        if name in _to_prefix:
            _counts[name] += 1
            col_names[i] = f"{_prefixes.get(_counts[name], str(_counts[name]))} {name}"
    col_names_unserscore = [i.replace(' ', '_') for i in col_names]
    return col_names_unserscore, col_starts


def parse_data_line(line, col_starts):
    line = line.rstrip("\r\n")
    fields = []
    for i, start in enumerate(col_starts):
        end = col_starts[i + 1] if i + 1 < len(col_starts) else None
        fields.append((line[start:end] if end is not None else line[start:]).strip())
    return fields


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse an HMMER tblout or domtblout file into a TSV.")
    parser.add_argument("-i", dest="input_table", required=True,
                        help="HMMER hits table (.tbl or .tbl.gz)")
    parser.add_argument("-o", "--outname", dest="outfile_name", required=True,
                        help="Output file path without .tsv extension")
    args = parser.parse_args()

    domain_table = Path(args.input_table)
    if not domain_table.is_file():
        raise FileNotFoundError(f"Input table not found: {args.input_table}")

    header_lines = []
    with open_file(domain_table) as fh:
        for line in fh:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                break

    col_names, col_starts = parse_hmmer_header(header_lines)

    with open(args.outfile_name + ".tsv", "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writerow(col_names)
        with open_file(domain_table) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                writer.writerow(parse_data_line(line, col_starts))
