#!/usr/bin/env python
import argparse
import csv
import os
from collections import Counter

# define rank order (top → bottom)
rank_order = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]

# Color palette for DATASET_BINARY fields
_COLORS = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990",
    "#dcbeff", "#9a6324", "#fffac8", "#800000", "#aaffc3",
]


def extract_value(s):
    """Return taxon label from a taxonomy part, or None if empty/dash."""
    s = s.strip()
    if not s or s == "-":
        return None
    if s.startswith("[") and "]_" in s:
        # [Rank]_Value format — strip the rank prefix
        return s.split("]_", 1)[1].replace(" ", "_")
    if not s.startswith("["):
        return s.replace(" ", "_")
    return None


def parse_lineage(tax):
    """Parse reversed taxonomy string into a deduplicated ordered path.

    The taxonomy column is ordered most-specific → least-specific, so we
    reverse it to get Realm first.  VIRify sometimes assigns the same taxon
    name to several consecutive ranks when only the genus is known; consecutive
    duplicates are collapsed to avoid redundant nodes in the tree.
    """
    parts_reversed = list(reversed(tax.split(";")))
    path = []
    for _rank, val_str in zip(rank_order, parts_reversed):
        val = extract_value(val_str)
        if val and (not path or path[-1] != val):
            path.append(val)
    return path


def add_path(tree, path):
    """Insert path[:-1] as internal nodes; return the node at path[-2]."""
    node = tree
    for name in path[:-1]:
        node = node.setdefault(name, {})
    return node


def build_newick(node):
    """Recursively build a Newick string with internal node names."""
    if not node:
        return ""
    children = []
    for name, subtree in node.items():
        if subtree:
            children.append(f"({build_newick(subtree)}){name}")
        else:
            children.append(name)
    return ",".join(children)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build an iTOL-ready Newick taxonomy tree from a VIRify/VITAP TSV."
    )
    parser.add_argument("-i", "--input", help="input TSV (id, taxonomy, ...)")
    parser.add_argument("-o", "--output", default="tree.nwk", help="output Newick file (default: tree.nwk)")
    parser.add_argument("--labels", default=None, help="iTOL labels file (default: <output>.labels.txt)")
    parser.add_argument("--counts", default=None, help="iTOL counts bar-chart file (default: <output>.counts.txt)")
    parser.add_argument(
        "--format", choices=["auto", "id", "count"], default="auto",
        help=(
            "Input format (default: auto). "
            "'auto': detect from first row — if col 0 is an integer, use 'count', else use 'id'. "
            "'id': id<TAB>taxonomy<TAB>... — taxonomy in col 1, semicolon-separated, "
            "reversed order (VIRify/VITAP style), has header. "
            "'count': count<TAB>tax1<TAB>tax2<TAB>... — Krona-style, "
            "already split and in forward order, no header."
        ),
    )
    parser.add_argument(
        "--meta", default=None,
        help="Metadata table produced by generate_taxonomy_table.py "
             "(columns: taxonomy, col1, col2, ...). "
             "One iTOL DATASET_BINARY file is created per metadata column.",
    )
    return parser.parse_args()


def detect_format(first_row):
    """Return 'count' if col 0 is an integer, 'id' otherwise."""
    try:
        int(first_row[0])
        return "count"
    except (ValueError, IndexError):
        return "id"


def generate_itol_metadata(meta_file, node_ids, base):
    """Read metadata table and write one iTOL DATASET_BINARY file per column.

    meta_file columns: taxonomy (;-joined path) + one or more metadata columns.
    node_ids: dict  lineage_key → (leaf_id, leaf_name, count)

    For each metadata column:
      - collect all unique values across all taxonomy paths
      - write a DATASET_BINARY file: one binary field (circle, shape=1) per unique value
      - each leaf gets 1 for values present in its metadata, 0 otherwise
    """
    with open(meta_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        meta_cols = [c for c in reader.fieldnames if c != "taxonomy"]
        # tax_key → {col: set of values}
        meta_rows = {}
        for row in reader:
            tax_key = row["taxonomy"].strip()
            meta_rows[tax_key] = {col: set(v for v in row.get(col, "").split(";") if v)
                                  for col in meta_cols}

    for col in meta_cols:
        # Collect all unique values for this column
        all_values = sorted({v for d in meta_rows.values() for v in d.get(col, set())})
        if not all_values:
            continue

        out_path = f"{base}.{col}.binary.txt"
        with open(out_path, "w") as out:
            out.write("DATASET_BINARY\nSEPARATOR COMMA\n")
            out.write(f"DATASET_LABEL,{col}\n")

            # One color per unique value, cycling through palette
            colors = ",".join(_COLORS[i % len(_COLORS)] for i in range(len(all_values)))
            out.write(f"COLOR,{_COLORS[0]}\n")
            out.write(f"FIELD_COLORS,{colors}\n")
            out.write(f"FIELD_LABELS,{','.join(all_values)}\n")
            # shape 1 = circle
            out.write(f"FIELD_SHAPES,{','.join('1' for _ in all_values)}\n")
            out.write("\nDATA\n")

            for lineage_key, (leaf_id, _name, _count) in node_ids.items():
                values_present = meta_rows.get(lineage_key, {}).get(col, set())
                bits = ",".join("1" if v in values_present else "0" for v in all_values)
                out.write(f"{leaf_id},{bits}\n")

        print(f"DATASET_BINARY written to {out_path}  ({len(all_values)} fields: {', '.join(all_values)})")


# ── 1. Count sequences per unique lineage path ────────────────────────────────

args = parse_args()
base = os.path.splitext(args.output)[0]
labels_file = args.labels or f"{base}.labels.txt"
counts_file = args.counts or f"{base}.counts.txt"

lineage_counts = Counter()
with open(args.input) as f:
    reader = csv.reader(f, delimiter="\t")
    first_row = next(reader)
    fmt = args.format if args.format != "auto" else detect_format(first_row)

    if fmt == "id":
        # first_row was the header — discard it and process remaining rows
        for row in reader:
            path = parse_lineage(row[1])
            if path:
                lineage_counts[";".join(path)] += 1
    else:
        # first_row is a data row — process it first, then continue
        for row in [first_row] + list(reader):
            count = int(row[0])
            path = [p.strip() for p in row[1:] if p.strip()]
            if path:
                lineage_counts[";".join(path)] += count

# ── 2. Build tree ─────────────────────────────────────────────────────────────
#
# Following the pattern in plot_tree.py:
#   • all taxonomy levels except the last become internal nodes
#   • the last level is represented by an integer leaf ID
#   • a LABELS file maps each integer ID back to the taxon name
#
# This lets the same taxon name appear both as an internal node (when a more
# specific lineage passes through it) and as a labelled leaf (when a lineage
# stops there), without any duplicate-name conflicts.

tree = {}
node_ids = {}   # lineage_key → (leaf_id, leaf_name, count)
counter = 0

for lineage_key in sorted(lineage_counts):
    path = lineage_key.split(";")
    count = lineage_counts[lineage_key]

    # Navigate / create internal nodes for all levels except the last
    parent = add_path(tree, path)

    # Add leaf with integer ID
    leaf_id = str(counter)
    parent[leaf_id] = {}
    node_ids[lineage_key] = (leaf_id, path[-1], count)
    counter += 1

# ── 3. Write Newick ───────────────────────────────────────────────────────────

newick = f"({build_newick(tree)})root;"
with open(args.output, "w") as out:
    out.write(newick)
print(f"Tree written to {args.output}  ({len(node_ids)} leaves, {sum(lineage_counts.values())} sequences)")

# ── 4. iTOL LABELS file — integer leaf ID → deepest taxon name ───────────────

with open(labels_file, "w") as out:
    out.write("LABELS\nSEPARATOR COMMA\n\nDATA\n")
    for _key, (leaf_id, leaf_name, _count) in node_ids.items():
        out.write(f"{leaf_id},{leaf_name}\n")
print(f"Labels written to {labels_file}")

# ── 5. iTOL DATASET_SIMPLEBAR — sequence count per leaf ──────────────────────

with open(counts_file, "w") as out:
    out.write("DATASET_SIMPLEBAR\nSEPARATOR COMMA\n")
    out.write("DATASET_LABEL,Sequence count\nCOLOR,#4a90d9\n\nDATA\n")
    for _key, (leaf_id, _name, count) in node_ids.items():
        out.write(f"{leaf_id},{count}\n")
print(f"Counts written to {counts_file}")

# ── 6. iTOL DATASET_BINARY — one file per metadata column ────────────────────

if args.meta:
    generate_itol_metadata(args.meta, node_ids, base)
