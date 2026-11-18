#!/usr/bin/env python3
"""Merge a file of the existing catalogue with the same file of a --mode update run.

The update annotates only the representatives of new clusters. A merged file is the
existing catalogue's rows, without the rows of representatives that a merge of clusters
dropped, followed by the update's rows.

Representatives are catalogue IDs (MGYV...), but proteins keep their source IDs, so the
rows to drop are given as one ID file (`drop-ids`): the dropped representatives and the
IDs of their proteins (CDS rows on their seqid in the existing reps GFF).

Subcommands
  drop-ids  --dropped-reps FILE --gff REPS_GFF -o drop_ids.txt
  merge     --format {tsv,csv,gff,fasta,tbl} [--id-column NAME | --id-index N]
            --old FILE --new FILE --drop-ids FILE -o OUT
  reps      --clusters-tsv FILE -o reps.tsv    (representative IDs of a vclust cluster TSV)
  concat    --format {gff,fasta} [--unique-ids] -o OUT FILE [FILE ...]
            the whole catalogue's sequence files: the existing catalogue's, then the update's
            (one ##gff-version line; --unique-ids: fail on an ID in two files)

Formats and the ID of a row
  tsv, csv  header line, then rows; ID in --id-column (header names are stripped) or
            --id-index; repeated header lines (concatenated chunks) are skipped
  gff       seqid (column 1); ##gff-version written once; ##sequence-region follows its seqid
  fasta     first word of the header line
  tbl       HMMER --tblout: first whitespace-separated field; '#' lines: the first file's
            leading comment block is kept, other comment lines are dropped

The update must not annotate an existing row again: an ID in both the kept existing rows
and the new rows is an error (for FASTA files this also catches duplicated protein IDs).
Plain or gzip input; the output is gzip-compressed when its name ends in .gz.
"""

import argparse
import csv
import gzip
import sys

MAX_EXAMPLES = 10


def open_text(path):
    with open(path, "rb") as fh:
        is_gzip = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", newline="") if is_gzip else open(path, newline="")


def open_out(path):
    return gzip.open(path, "wt", newline="") if path.endswith(".gz") else open(path, "w", newline="")


def examples(items):
    items = sorted(items)
    return ", ".join(items[:MAX_EXAMPLES]) + (", ..." if len(items) > MAX_EXAMPLES else "")


def read_ids(path):
    with open_text(path) as fh:
        return {line.strip() for line in fh if line.strip()}


# ---------------------------------------------------------------------------
# drop-ids
# ---------------------------------------------------------------------------

def gff_attribute(attributes, key):
    for field in attributes.split(";"):
        name, _, value = field.partition("=")
        if name.strip() == key:
            return value.strip()
    return None


def drop_ids(dropped_reps_path, gff_path):
    """Dropped representatives and the IDs of their CDS features in the reps GFF."""
    reps = read_ids(dropped_reps_path)
    ids = set(reps)
    seen = set()
    with open_text(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[0] not in reps:
                continue
            seen.add(cols[0])
            if cols[2] == "CDS":
                cds_id = gff_attribute(cols[8], "ID")
                if cds_id:
                    ids.add(cds_id)
    missing = reps - seen
    if missing:
        sys.exit(f"Error: {len(missing)} dropped representative(s) are not in {gff_path}, e.g. {examples(missing)}")
    return ids


# ---------------------------------------------------------------------------
# merge: readers yield (kind, text, id); kind is 'header', 'comment' or 'row'
# ---------------------------------------------------------------------------

def table_reader(path, delimiter, id_column, id_index):
    header = None
    index = id_index
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = next(csv.reader([line.rstrip("\r\n")], delimiter=delimiter))
            if header is None:
                header = [f.strip() for f in fields]
                if id_column is not None:
                    if id_column not in header:
                        sys.exit(f"Error: {path} has no column '{id_column}' (columns: {', '.join(header)})")
                    index = header.index(id_column)
                yield "header", line, None
                continue
            if [f.strip() for f in fields] == header:  # repeated header
                continue
            if index >= len(fields):
                sys.exit(f"Error: {path}: row with {len(fields)} field(s) has no ID column {index + 1}: {line.strip()}")
            yield "row", line, fields[index].strip()


def gff_reader(path):
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("##gff-version"):
                yield "header", line, None
            elif line.startswith("##sequence-region"):
                parts = line.split()
                yield "row", line, parts[1] if len(parts) > 1 else None
            elif line.startswith("#"):
                yield "comment", line, None
            else:
                yield "row", line, line.split("\t", 1)[0]


def fasta_reader(path):
    record, record_id = [], None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if record:
                    yield "row", "".join(record), record_id
                record = [line]
                words = line[1:].split()
                record_id = words[0] if words else ""
            elif record:
                record.append(line if line.endswith("\n") else line + "\n")
            elif line.strip():
                sys.exit(f"Error: {path} does not start with a FASTA header")
    if record:
        yield "row", "".join(record), record_id


def tbl_reader(path):
    leading = True
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                yield ("header" if leading else "comment"), line, None
            else:
                leading = False
                yield "row", line, line.split(None, 1)[0]


def reader(fmt, path, id_column, id_index):
    if fmt in ("tsv", "csv"):
        return table_reader(path, "\t" if fmt == "tsv" else ",", id_column, id_index)
    return {"gff": gff_reader, "fasta": fasta_reader, "tbl": tbl_reader}[fmt](path)


def normalised(header_lines, fmt):
    if fmt in ("tsv", "csv"):
        delimiter = "\t" if fmt == "tsv" else ","
        return [[f.strip() for f in next(csv.reader([line.rstrip("\r\n")], delimiter=delimiter))] for line in header_lines]
    return [line.strip() for line in header_lines]


def merge(fmt, old, new, drop, out, id_column=None, id_index=None):
    parts = {}
    for name, path in (("old", old), ("new", new)):
        header, comments, rows = [], [], []
        for kind, text, row_id in reader(fmt, path, id_column, id_index):
            if kind == "header":
                header.append(text)
            elif kind == "comment":
                comments.append(text)
            else:
                rows.append((text, row_id))
        parts[name] = (header, comments, rows)

    old_header, old_comments, old_rows = parts["old"]
    new_header, new_comments, new_rows = parts["new"]
    if fmt in ("tsv", "csv") and old_header and new_header and normalised(old_header, fmt) != normalised(new_header, fmt):
        sys.exit(f"Error: the columns of {old} and {new} differ:\n  existing: {old_header[0].strip()}\n"
                 f"  update:   {new_header[0].strip()}\n"
                 "Was the existing catalogue annotated with another tool or database version?")
    # tbl: the header block of the first file that has one; gff: one ##gff-version line
    header = old_header or new_header
    if fmt == "gff":
        header = header[:1]

    kept = [(text, row_id) for text, row_id in old_rows if row_id not in drop]
    clash = {row_id for _, row_id in kept} & {row_id for _, row_id in new_rows}
    clash.discard(None)
    if clash:
        sys.exit(f"Error: {len(clash)} ID(s) of the update are already in the existing {old}, e.g. {examples(clash)}")

    if fmt == "fasta":
        seen, duplicated = set(), set()
        for _, row_id in kept + new_rows:
            if row_id in seen:
                duplicated.add(row_id)
            seen.add(row_id)
        if duplicated:
            sys.exit(f"Error: {len(duplicated)} sequence ID(s) occur more than once in {new}, e.g. {examples(duplicated)}")

    with open_out(out) as fh:
        fh.writelines(header)
        if fmt == "gff":
            fh.writelines(old_comments + new_comments)
        for text, _ in kept + new_rows:
            fh.write(text if text.endswith("\n") else text + "\n")

    print(f"{out}: {len(kept)} existing row(s) kept, {len(old_rows) - len(kept)} dropped, {len(new_rows)} new")


def reps_list(clusters_tsv, out):
    reps, seen = [], set()
    with open_text(clusters_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:2] != ["object", "cluster"]:
            sys.exit(f"Error: {clusters_tsv} is not a vclust cluster TSV (header object, cluster)")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2 and cols[1] not in seen and cols[1] != "cluster":
                seen.add(cols[1])
                reps.append(cols[1])
    with open_out(out) as fh:
        fh.writelines(f"{rep}\n" for rep in reps)
    print(f"{out}: {len(reps)} representatives")


def concat(fmt, paths, out, unique_ids=False):
    seen, duplicated, n = set(), set(), 0
    with open_out(out) as fh:
        if fmt == "gff":
            fh.write("##gff-version 3\n")
        for path in paths:
            for kind, text, row_id in reader(fmt, path, None, None):
                if kind == "header":
                    continue
                if kind == "row" and (fmt == "fasta" or not text.startswith("#")):
                    if unique_ids and row_id in seen and fmt == "fasta":
                        duplicated.add(row_id)
                    seen.add(row_id)
                    n += 1
                fh.write(text if text.endswith("\n") else text + "\n")
    if duplicated:
        sys.exit(f"Error: {len(duplicated)} ID(s) occur in more than one of {', '.join(paths)}, e.g. {examples(duplicated)}")
    print(f"{out}: {n} {'sequences' if fmt == 'fasta' else 'rows'} from {len(paths)} file(s)")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("drop-ids", help="IDs to drop: dropped representatives and their proteins")
    p.add_argument("--dropped-reps", required=True)
    p.add_argument("--gff", required=True, help="reps GFF of the existing catalogue")
    p.add_argument("-o", "--output", required=True)

    p = sub.add_parser("merge", help="merge an existing file with the update's file")
    p.add_argument("--format", required=True, choices=["tsv", "csv", "gff", "fasta", "tbl"])
    group = p.add_mutually_exclusive_group()
    group.add_argument("--id-column", help="tsv/csv: name of the ID column")
    group.add_argument("--id-index", type=int, help="tsv/csv: 0-based index of the ID column")
    p.add_argument("--old", required=True)
    p.add_argument("--new", required=True)
    p.add_argument("--drop-ids", required=True)
    p.add_argument("-o", "--output", required=True)

    p = sub.add_parser("reps", help="representative IDs of a vclust cluster TSV")
    p.add_argument("--clusters-tsv", required=True)
    p.add_argument("-o", "--output", required=True)

    p = sub.add_parser("concat", help="concatenate the existing and the update's sequence files")
    p.add_argument("--format", required=True, choices=["gff", "fasta"])
    p.add_argument("--unique-ids", action="store_true", help="fasta: fail on a sequence ID in two files")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("files", nargs="+")

    args = parser.parse_args(argv)
    if args.command == "drop-ids":
        ids = drop_ids(args.dropped_reps, args.gff)
        with open_out(args.output) as fh:
            fh.writelines(f"{i}\n" for i in sorted(ids))
        print(f"{args.output}: {len(ids)} ID(s) to drop")
    elif args.command == "merge":
        if args.format in ("tsv", "csv") and args.id_column is None and args.id_index is None:
            parser.error("--format tsv/csv needs --id-column or --id-index")
        merge(args.format, args.old, args.new, read_ids(args.drop_ids), args.output,
              id_column=args.id_column, id_index=args.id_index)
    elif args.command == "concat":
        concat(args.format, args.files, args.output, args.unique_ids)
    else:
        reps_list(args.clusters_tsv, args.output)


if __name__ == "__main__":
    main()
