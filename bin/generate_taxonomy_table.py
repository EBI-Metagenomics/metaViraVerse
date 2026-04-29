#!/usr/bin/env python3
"""
Generate Krona plot input file from a TSV with a taxonomy column.

Supports two taxonomy formats:

  forward (default)
    GFF-style: Realm;Kingdom;Phylum;...  (URL-encoded, top → bottom)

  reversed  (--reversed)
    VIRify/VITAP TSV-style: most-specific → least-specific, with optional
    [Rank]_Value labels and '-' placeholders.

Output format (Krona): count<TAB>rank1<TAB>rank2<TAB>...
"""

import argparse
import csv
import sys
from collections import Counter
from urllib.parse import unquote


# ── Constants ─────────────────────────────────────────────────────────────────

STANDARD_LEVELS_FORWARD = [
    "realm", "kingdom", "phylum", "class",
    "order", "family", "subfamily", "genus", "species",
]

RANK_ORDER_REVERSED = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]


# ── Forward (GFF-style) taxonomy ──────────────────────────────────────────────

def parse_taxonomy_forward(taxonomy_str, standard_levels=None):
    """
    Parse a forward GFF taxonomy string and fill missing levels.

    Missing/empty positions become 'unclassified_<level>_<parent>'.
    Returns a list of capitalised taxonomy terms.
    """
    if not taxonomy_str or taxonomy_str in ("NA", "unclassified"):
        return []

    taxonomy_str = unquote(taxonomy_str)
    parts = [p.strip() for p in taxonomy_str.split(";")]
    while parts and parts[-1] == "":
        parts.pop()

    if standard_levels is None:
        standard_levels = STANDARD_LEVELS_FORWARD

    filled = []
    last_valid = "root"

    for i, part in enumerate(parts):
        level_name = standard_levels[i] if i < len(standard_levels) else f"level_{i}"
        part_lower = part.lower()

        if part_lower and part_lower not in ("", "unclassified"):
            filled.append(part)
            last_valid = part
        else:
            if level_name == 'realm' and last_valid == 'root':
                filled.append('unclassified')
            else:
                filled.append(f"unclassified_{level_name}_{last_valid}")

    return [x.capitalize() for x in filled]


# ── Reversed (VIRify/VITAP TSV-style) taxonomy ───────────────────────────────

def _extract_value(s):
    """Extract taxon label from one reversed-taxonomy token, or None if absent."""
    s = s.strip()
    if not s or s == "-":
        return None
    if s.startswith("[") and "]_" in s:
        return s.split("]_", 1)[1].replace(" ", "_")
    if not s.startswith("["):
        return s.replace(" ", "_")
    return None


def parse_taxonomy_reversed(taxonomy_str, rank_order=None):
    """
    Parse a reversed taxonomy string (most-specific → least-specific).

    Reverses the string so index 0 = Realm, then zips with rank_order to cap
    depth.  Consecutive duplicate names (VIRify assigns the same taxon to
    several ranks when only genus-level info is available) are collapsed.
    Returns a list of taxon labels.
    """
    if rank_order is None:
        rank_order = RANK_ORDER_REVERSED
    parts_reversed = list(reversed(taxonomy_str.split(";")))
    path = []
    for _rank, val_str in zip(rank_order, parts_reversed):
        val = _extract_value(val_str)
        if val and (not path or path[-1] != val):
            path.append(val)
    return path


# ── Metadata table ────────────────────────────────────────────────────────────

def generate_metadata_table(seq_tax, meta_file, meta_id_col, columns, output_file):
    """
    Generate a metadata table: one row per unique taxonomy path, columns are
    unique (;-joined) values of the requested metadata fields.

    Output format (TSV):
        taxonomy<TAB>col1<TAB>col2...

    where taxonomy is the ;-joined path (same as in the Krona output) and each
    metadata column contains all unique values found for sequences with that
    taxonomy, separated by ';'.

    Args:
        seq_tax:      dict mapping seq_id -> taxonomy key (;-joined path string)
        meta_file:    metadata TSV (e.g. combined_meta.tsv)
        meta_id_col:  column in meta_file that matches sequence IDs in seq_tax
        columns:      list of metadata columns to include (e.g. ['type', 'biomes'])
        output_file:  output TSV path
    """
    from collections import defaultdict

    # Load metadata: id -> {col: value}
    meta = {}
    with open(meta_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mid = row.get(meta_id_col, "").strip()
            if mid:
                meta[mid] = {col: row.get(col, "").strip() for col in columns}
    print(f"Got {len(meta)} lines from metadata input ")

    # Group sequences by taxonomy path, collecting metadata values per column
    tax_meta = defaultdict(lambda: {col: set() for col in columns})

    for seq_id, tax_key in seq_tax.items():
        seq_meta = meta.get(seq_id, {})

        for col in columns:
            val = seq_meta.get(col, "").strip()
            if val:
                tax_meta[tax_key][col].update(val.split(','))


    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["taxonomy"] + columns)
        for tax_key in sorted(tax_meta):
            row = [tax_key]
            for col in columns:
                vals = sorted(tax_meta[tax_key][col])
                row.append(";".join(vals) if vals else "")
            writer.writerow(row)

    print(f"Metadata table written to: {output_file}")
    print(f"  Unique taxonomy paths: {len(tax_meta)}")


# ── Krona output ──────────────────────────────────────────────────────────────

def generate_krona_file(taxonomy_paths, output_file):
    """
    Write a Krona-format file from a list of taxonomy path tuples.

    Format: count<TAB>rank1<TAB>rank2<TAB>...
    """
    counts = Counter(taxonomy_paths)

    with open(output_file, "w") as out:
        for path, count in sorted(counts.items(), key=lambda x: -x[1]):
            out.write(f"{count}\t{chr(9).join(path)}\n")

    print(f"Krona file written to: {output_file}")
    print(f"Total unique taxonomy paths: {len(counts)}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a Krona plot input file from a TSV taxonomy column.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Forward GFF-style taxonomy (output of extract_reps_stats.py)
  %(prog)s --input stats.tsv --output krona.txt

  # Reversed VIRify/VITAP TSV-style taxonomy
  %(prog)s --input salmon_vitap_best.tsv --output krona.txt --reversed

  # Custom column name or separator
  %(prog)s --input data.csv --sep , --taxonomy-col lineage --output krona.txt

  # Custom rank levels (e.g. drop subfamily and species)
  %(prog)s --input stats.tsv --output krona.txt --lineage realm,kingdom,phylum,class,order,family,genus
        """,
    )
    parser.add_argument("--input", required=True,
                        help="TSV (or CSV) file with a taxonomy column.")
    parser.add_argument("--output", required=True,
                        help="Output Krona-format file.")
    parser.add_argument("--taxonomy-col", default="taxonomy",
                        help="Name of the taxonomy column (default: taxonomy).")
    parser.add_argument("--sep", default="\t",
                        help="Input field separator (default: TAB).")
    parser.add_argument("--reversed", action="store_true",
                        help="Taxonomy is in reversed order (VIRify/VITAP TSV style).")
    # ── metadata table ────────────────────────────────────────────────────────
    parser.add_argument("--meta", default=None,
                        help="Metadata file to join with the input sequences "
                             "(e.g. combined_meta.tsv).")
    parser.add_argument("--meta-id-col", default="description",
                        help="Column in --meta that matches sequence IDs from "
                             "--input (default: description).")
    parser.add_argument("--id-col", default=None,
                        help="ID column in --input (default: first column).")
    parser.add_argument("--columns", default=None,
                        help="Comma-separated metadata columns to include in the "
                             "output table (e.g. 'type,biomes'). "
                             "Required when --meta is given.")
    parser.add_argument("--meta-output", default=None,
                        help="Output TSV for the metadata table. "
                             "Required when --meta is given.")
    parser.add_argument(
        "--lineage",
        default=None,
        help=(
            "Comma-separated list of rank names, ordered top → bottom "
            "(e.g. 'realm,kingdom,phylum,class,order,family,genus'). "
            "Controls depth of parsing and fills missing-level labels in "
            "forward mode. Defaults to: " + ",".join(STANDARD_LEVELS_FORWARD)
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()

    lineage = [r.strip() for r in args.lineage.split(",")] if args.lineage else STANDARD_LEVELS_FORWARD

    if args.reversed:
        parse_fn = lambda tax: parse_taxonomy_reversed(tax, rank_order=lineage)
    else:
        parse_fn = lambda tax: parse_taxonomy_forward(tax, standard_levels=lineage)

    need_meta = bool(args.meta)
    if need_meta:
        if not args.columns:
            print("Error: --columns is required when --meta is given.", file=sys.stderr)
            sys.exit(1)
        if not args.meta_output:
            print("Error: --meta-output is required when --meta is given.", file=sys.stderr)
            sys.exit(1)
        columns = [c.strip() for c in args.columns.split(",")]

    paths = []
    seq_tax = {}  # seq_id -> ";".join(path), built only when --meta is given
    with open(args.input) as f:
        reader = csv.DictReader(f, delimiter=args.sep)
        if args.taxonomy_col not in (reader.fieldnames or []):
            print(
                f"Error: column '{args.taxonomy_col}' not found. "
                f"Available columns: {reader.fieldnames}",
                file=sys.stderr,
            )
            sys.exit(1)
        _id_col = args.id_col or reader.fieldnames[0]
        for row in reader:
            tax = row[args.taxonomy_col]
            if not tax or tax == "NA":
                continue
            path = parse_fn(tax)
            if path:
                paths.append(tuple(path))
                if need_meta:
                    seq_id = row.get(_id_col, "").strip()
                    if seq_id:
                        seq_tax[seq_id] = ";".join(path)

    if not paths:
        print("Warning: no taxonomy paths found — output will be empty.", file=sys.stderr)

    generate_krona_file(paths, args.output)

    if need_meta:
        generate_metadata_table(
            seq_tax, args.meta, args.meta_id_col, columns, args.meta_output,
        )


if __name__ == "__main__":
    main()
