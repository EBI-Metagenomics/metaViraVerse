#!/usr/bin/env python3
"""
Create a Sankey plot from taxonomy data with proper handling of missing taxonomic levels.

Handles missing/empty taxonomy levels by creating "unclassified_[level]_[parent_taxa]" labels
to ensure all links in the Sankey diagram are properly connected.
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


def parse_taxonomy(taxonomy_str, standard_levels=None, lineage_reverted=False):
    """
    Parse taxonomy string and fill in missing levels.

    Args:
        taxonomy_str: Semicolon-separated taxonomy string (e.g., "Viruses;Heunggongvirae;;Caudoviricetes")
        standard_levels: List of standard taxonomic level names (e.g., ['domain', 'phylum', 'class', 'order'])

    Returns:
        List of taxonomy terms with missing levels filled as "unclassified_[level]_[parent]"
    """
    if not taxonomy_str or taxonomy_str == "NA":
        return []

    # Handle completely unclassified sequences at first level (case-insensitive)
    if taxonomy_str.strip().lower() == "unclassified":
        return ["Unclassified"]

    # URL decode and split
    taxonomy_str = unquote(taxonomy_str)
    parts = [p.strip() for p in taxonomy_str.split(';')]

    if lineage_reverted:
        parts = parts[::-1]

    # remove empty annotations
    while parts and parts[-1] == '' or parts[-1] == '-':
        parts.pop()

    # If no standard levels provided, use generic numbering
    if standard_levels is None:
        standard_levels = [f'level_{i}' for i in range(len(parts))]

    filled_taxonomy = []
    last_valid = "root"

    for i, part in enumerate(parts):
        level_name = standard_levels[i] if i < len(standard_levels) else f'level_{i}'

        if part and part.strip():
            # Check if it's "unclassified" (case-insensitive) and normalize it
            if part.strip().lower() == "unclassified":
                filled_taxonomy.append("Unclassified")
                last_valid = "Unclassified"
            else:
                # Valid taxonomy term
                filled_taxonomy.append(part.strip())
                last_valid = part.strip()
        else:
            # Missing level at first position - treat as Unclassified at root level
            if i == 0:
                filled_taxonomy.append("Unclassified")
                last_valid = "Unclassified"
            else:
                # Missing level at intermediate position - create unclassified label
                unclassified_label = f"unclassified_{level_name}_{last_valid}"
                filled_taxonomy.append(unclassified_label)

    return filled_taxonomy


def open_text(path):
    """Open a plain or gzip-compressed text file (detected by magic bytes) for reading."""
    with open(path, 'rb') as fh:
        is_gzip = fh.read(2) == b'\x1f\x8b'
    if is_gzip:
        return gzip.open(path, 'rt')
    return open(path, 'r')


def is_taxid_lineage(taxonomy_str):
    """True for a lineage written as NCBI taxids ('1;131567;2;1224*'), e.g. a CAT_pack lineage that
    could not be translated into names. Such lineages are ignored (CAT_pack marks taxa with '*')."""
    parts = [p.strip().rstrip("*") for p in unquote(taxonomy_str or "").split(";") if p.strip()]
    return bool(parts) and all(p.isdigit() for p in parts)


OTHER_LABEL = "Other"
OTHER_KEY = "__other__"


def truncate_to_stop_level(taxonomy_list_with_counts, levels, stop_level):
    """
    Cut every taxonomy path after the requested stop level.

    Args:
        taxonomy_list_with_counts: List of (count, taxonomy_path) tuples
        levels: Ordered list of taxonomic level names
        stop_level: Level name to stop at (inclusive)

    Returns:
        List of (count, truncated_taxonomy_path) tuples
    """
    depth = levels.index(stop_level) + 1
    return [(count, path[:depth]) for count, path in taxonomy_list_with_counts]


def apply_top_hits(taxonomy_list_with_counts, top_hits):
    """
    Keep only the top N taxa (by count) at each level; merge the rest into one shared
    "Other" node (the same node for every level). Paths merged into "Other" stop there.

    Levels are processed top-down, so ranking at a level only considers paths whose
    ancestors all made it into the top N.

    Args:
        taxonomy_list_with_counts: List of (count, taxonomy_path) tuples
        top_hits: Number of taxa to keep per level

    Returns:
        List of (count, taxonomy_path) tuples, where collapsed paths end with the
        OTHER_KEY marker
    """
    paths = [(count, list(path)) for count, path in taxonomy_list_with_counts]
    collapsed = [False] * len(paths)
    max_depth = max((len(path) for _, path in paths), default=0)

    for depth in range(max_depth):
        prefix_counts = defaultdict(int)
        for i, (count, path) in enumerate(paths):
            if not collapsed[i] and len(path) > depth:
                prefix_counts[tuple(path[:depth + 1])] += count

        ranked = sorted(prefix_counts.items(), key=lambda item: (-item[1], item[0]))
        keep = {prefix for prefix, _ in ranked[:top_hits]}

        for i, (count, path) in enumerate(paths):
            if not collapsed[i] and len(path) > depth and tuple(path[:depth + 1]) not in keep:
                paths[i] = (count, path[:depth] + [OTHER_KEY])
                collapsed[i] = True

    return paths


def build_sankey_data(taxonomy_list_with_counts):
    """
    Build Sankey diagram data from list of taxonomy paths with counts.

    Args:
        taxonomy_list_with_counts: List of (count, taxonomy_path) tuples

    Returns:
        Tuple of (nodes, links) where:
        - nodes is a list of unique taxonomy terms
        - links is a list of (source_idx, target_idx, value) tuples
    """
    # Use full path as node key to handle same names at different levels
    # node_key_to_label maps "path|to|node" -> "node"
    # node_keys stores all unique keys
    node_keys = set()
    node_key_to_label = {}
    edge_counts = defaultdict(int)

    root_key = "root"
    node_keys.add(root_key)
    node_key_to_label[root_key] = "All sequences"

    # Build nodes and edges with unique keys based on full path
    for count, taxonomy_path in taxonomy_list_with_counts:
        if not taxonomy_path:
            continue

        # Build path keys for each node in the taxonomy
        current_path_parts = []
        prev_key = root_key

        for level, taxon in enumerate(taxonomy_path):
            current_path_parts.append(taxon)
            if taxon == OTHER_KEY:
                # A single "Other" node shared by all levels and parents
                node_key = OTHER_KEY
                label = OTHER_LABEL
            else:
                # Node key includes full path to make it unique
                node_key = "|".join(current_path_parts)
                label = taxon

            # Add node
            node_keys.add(node_key)
            node_key_to_label[node_key] = label

            # Add edge from previous level to current level
            edge_counts[(prev_key, node_key)] += count

            prev_key = node_key

    # Convert to sorted lists for consistent indexing
    sorted_keys = sorted(node_keys)
    node_labels = [node_key_to_label[key] for key in sorted_keys]
    key_to_idx = {key: idx for idx, key in enumerate(sorted_keys)}

    # Build links using indices
    links = []
    for (source_key, target_key), count in edge_counts.items():
        links.append({
            'source': key_to_idx[source_key],
            'target': key_to_idx[target_key],
            'value': count
        })

    return node_labels, links


def sankey_links_table(taxonomy_list_with_counts, levels):
    """
    Rows of the table behind the Sankey plot: one row per link (flow between two nodes),
    with the same nodes and counts as build_sankey_data.

    Args:
        taxonomy_list_with_counts: List of (count, taxonomy_path) tuples, as passed to build_sankey_data
        levels: Ordered list of taxonomic level names

    Returns:
        List of dicts with source_level, source, target_level, target, target_lineage and
        sequences (number of sequences flowing through the link), ordered by level, then by
        decreasing count
    """
    edge_counts = defaultdict(int)
    edge_info = {}
    for count, taxonomy_path in taxonomy_list_with_counts:
        prev_key, prev_label, prev_level = "root", "All sequences", "root"
        for depth, taxon in enumerate(taxonomy_path):
            level = levels[depth] if depth < len(levels) else f"level_{depth + 1}"
            if taxon == OTHER_KEY:
                key, label, lineage = OTHER_KEY, OTHER_LABEL, OTHER_LABEL
            else:
                key = "|".join(taxonomy_path[:depth + 1])
                label, lineage = taxon, ";".join(taxonomy_path[:depth + 1])
            edge_counts[(prev_key, key)] += count
            edge_info[(prev_key, key)] = (depth, prev_level, prev_label, level, label, lineage)
            prev_key, prev_label, prev_level = key, label, level

    rows = []
    for edge, count in edge_counts.items():
        depth, source_level, source, target_level, target, lineage = edge_info[edge]
        rows.append({
            "source_level": source_level,
            "source": source,
            "target_level": target_level,
            "target": target,
            "target_lineage": lineage,
            "sequences": count,
            "_depth": depth,
        })
    rows.sort(key=lambda r: (r["_depth"], -r["sequences"], r["target_lineage"]))
    for row in rows:
        del row["_depth"]
    return rows


TABLE_COLUMNS = ["source_level", "source", "target_level", "target", "target_lineage", "sequences"]


def write_links_table(rows, output_file):
    """Write the rows of sankey_links_table as a TSV with a header."""
    with open(output_file, "w") as f:
        f.write("\t".join(TABLE_COLUMNS) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) for c in TABLE_COLUMNS) + "\n")
    print(f"Sankey table saved to: {output_file}")


def create_sankey_plotly(nodes, links, output_file, title="Viral Taxonomy Sankey Diagram"):
    """
    Create Sankey diagram using plotly.

    Args:
        nodes: List of node labels
        links: List of link dictionaries with 'source', 'target', 'value'
        output_file: Output HTML file path
        title: Plot title
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("Error: plotly is required. Install with: pip install plotly", file=sys.stderr)
        sys.exit(1)

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
        ),
        link=dict(
            source=[link['source'] for link in links],
            target=[link['target'] for link in links],
            value=[link['value'] for link in links],
        )
    )])

    fig.update_layout(
        title_text=title,
        font_size=10,
        height=800
    )

    fig.write_html(output_file)
    print(f"Sankey plot saved to: {output_file}")


def read_tsv_taxonomy(tsv_file, taxonomy_column='taxonomy', delimiter='\t',
                      reps_only=False, id_column='ID', rep_column='Cluster_rep'):
    """
    Read taxonomy data from TSV (or CSV) file.

    Supports two formats:
    1. Standard TSV/CSV with header and taxonomy column
    2. Krona format: count\trank1\trank2\t... (no header)

    Args:
        tsv_file: Path to TSV/CSV file
        taxonomy_column: Name of the taxonomy column (for standard TSV/CSV)
        delimiter: Field separator: '\t' for TSV, ',' for CSV
        reps_only: Keep only cluster representatives: rows where id_column equals rep_column
            (standard TSV/CSV only; ignored for Krona format)
        id_column: Name of the sequence ID column (for reps_only)
        rep_column: Name of the cluster representative column (for reps_only)

    Returns:
        List of (count, taxonomy_string) tuples
    """
    import csv

    taxonomy_data = []

    with open_text(tsv_file) as f:
        # Peek at first line to detect format
        first_line = f.readline().strip()
        if not first_line:  # File is completely empty
            print("File is empty")
            exit(0)
        f.seek(0)

        # Check if it's Krona format (first column is a number)
        first_parts = first_line.split(delimiter)
        is_krona_format = False
        if first_parts:
            try:
                int(first_parts[0])
                is_krona_format = True
            except ValueError:
                is_krona_format = False

        if is_krona_format:
            # Krona format: count\trank1\trank2\t...
            print("Detected Krona format (count\ttaxonomy_ranks)")
            for line in f:
                parts = line.strip().split(delimiter)
                if not parts or not parts[0].strip():
                    continue

                try:
                    count = int(parts[0])
                    # Get taxonomy parts (skip count)
                    taxonomy_parts = parts[1:]

                    if not taxonomy_parts or all(not p.strip() for p in taxonomy_parts):
                        # Empty taxonomy - treat as Unclassified at first level
                        taxonomy_str = 'Unclassified'
                    else:
                        # Join taxonomy parts with semicolon, keep empty strings for proper level handling
                        taxonomy_str = ';'.join(taxonomy_parts)

                    taxonomy_data.append((count, taxonomy_str))
                except (ValueError, IndexError):
                    continue

        else:
            # Standard TSV format with header
            print(f"Detected standard {'CSV' if delimiter == ',' else 'TSV'} format with header")
            reader = csv.DictReader(f, delimiter=delimiter)
            skipped_taxids = 0

            required = [taxonomy_column] + ([id_column, rep_column] if reps_only else [])
            for column in required:
                if column not in reader.fieldnames:
                    print(f"Error: Column '{column}' not found in input file.", file=sys.stderr)
                    print(f"Available columns: {', '.join(reader.fieldnames)}", file=sys.stderr)
                    sys.exit(1)

            for row in reader:
                if reps_only and row.get(id_column) != row.get(rep_column):
                    continue
                tax = row.get(taxonomy_column, '')
                if is_taxid_lineage(tax):
                    skipped_taxids += 1
                elif tax and tax != 'NA':
                    taxonomy_data.append((1, tax))
            if skipped_taxids:
                print(f"Ignored {skipped_taxids} lineages written as taxids")

    return taxonomy_data


def read_gff_taxonomy(gff_file):
    """
    Read taxonomy data directly from GFF file.

    Args:
        gff_file: Path to GFF file

    Returns:
        List of (count, taxonomy_string) tuples
    """
    taxonomy_data = []

    with open_text(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            # Parse attributes
            attrs = {}
            for attr in fields[8].split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attrs[key.strip().lower()] = val.strip()

            tax = attrs.get('taxonomy', '')
            if tax and tax != 'NA' and not is_taxid_lineage(tax):
                taxonomy_data.append((1, tax))

    return taxonomy_data


def main():
    parser = argparse.ArgumentParser(
        description='Create Sankey plot from viral taxonomy data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From TSV file (output of extract_reps_stats.py)
  %(prog)s --input stats.tsv --output taxonomy_sankey.html

  # From Krona format file (count\\ttaxonomy_ranks)
  %(prog)s --input krona.txt --output taxonomy_sankey.html

  # From CSV file
  %(prog)s --input stats.csv --input-format csv --output taxonomy_sankey.html

  # From GFF file directly
  %(prog)s --input viral.gff --input-format gff --output taxonomy_sankey.html

  # Specify custom taxonomy levels
  %(prog)s --input stats.tsv --output sankey.html --levels domain phylum class order family genus

  # Plot down to family only, keeping the 10 most abundant taxa per level (rest -> Other)
  %(prog)s --input stats.tsv --output sankey.html --stop-level family --top-hits 10

Standard viral taxonomy levels (default):
  realm, kingdom, phylum, class, order, family, subfamily, genus, species

TSV/CSV Format Support (--input-format tsv or csv):
  1. Standard TSV/CSV with header: taxonomy column contains semicolon-separated ranks
  2. Krona format (auto-detected): count\\trank1\\trank2\\trank3\\t... (no header)
        """
    )

    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        help='Input file (TSV with taxonomy column, Krona format, or GFF); may be gzip-compressed'
    )

    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output HTML file for Sankey plot'
    )

    parser.add_argument(
        '--output-table',
        type=Path,
        default=None,
        help='Also write the table the plot is drawn from (one row per link: source and target node, '
             'their levels, the target lineage and the number of sequences), as TSV'
    )

    parser.add_argument(
        '--input-format',
        choices=['tsv', 'csv', 'gff'],
        default='tsv',
        help='Input file format: tsv or csv (both auto-detect Krona format) or gff (default: tsv)'
    )

    parser.add_argument(
        '--taxonomy-column',
        default='taxonomy',
        help='Name of taxonomy column in TSV file (default: taxonomy)'
    )

    parser.add_argument(
        '--levels',
        nargs='+',
        default=['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'species'],
        help='Taxonomic level names in order'
    )

    parser.add_argument(
        '--stop-level',
        default=None,
        help='Deepest taxonomic level to plot; must be one of --levels (default: plot all levels)'
    )

    parser.add_argument(
        '--top-hits',
        type=int,
        default=None,
        help='Keep only the N most abundant taxa at each level; the rest are merged into "Other" (default: keep all)'
    )

    parser.add_argument(
        '--reps-only',
        action='store_true',
        help='Plot cluster representatives only: keep the rows where the --id-column value equals '
             'the --rep-column value (TSV/CSV with header only)'
    )

    parser.add_argument(
        '--id-column',
        default='ID',
        help='Name of the sequence ID column, used with --reps-only (default: ID)'
    )

    parser.add_argument(
        '--rep-column',
        default='Cluster_rep',
        help='Name of the cluster representative column, used with --reps-only (default: Cluster_rep)'
    )

    parser.add_argument(
        '--title',
        default='Viral Taxonomy Sankey Diagram',
        help='Plot title'
    )


    parser.add_argument(
        '--lineage-reverted',
        action='store_true',
        help='Specify that argument if your lineage is going from species to realm'
    )

    args = parser.parse_args()

    # --levels is user-configurable, so the valid --stop-level choices are validated here
    if args.stop_level is not None and args.stop_level not in args.levels:
        parser.error(f"argument --stop-level: invalid choice: '{args.stop_level}' "
                     f"(choose from {', '.join(args.levels)})")

    if args.top_hits is not None and args.top_hits < 1:
        parser.error("argument --top-hits: must be a positive integer")

    # Validate input file
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Read taxonomy data
    print(f"Reading taxonomy data from {args.input}...")
    if args.input_format == 'tsv':
        taxonomy_data = read_tsv_taxonomy(args.input, args.taxonomy_column, reps_only=args.reps_only,
                                          id_column=args.id_column, rep_column=args.rep_column)
    elif args.input_format == 'csv':
        taxonomy_data = read_tsv_taxonomy(args.input, args.taxonomy_column, delimiter=',',
                                          reps_only=args.reps_only, id_column=args.id_column,
                                          rep_column=args.rep_column)
    else:
        taxonomy_data = read_gff_taxonomy(args.input)

    if not taxonomy_data:
        print("Error: No taxonomy data found in input file", file=sys.stderr)
        sys.exit(0)

    # Calculate total sequences
    total_sequences = sum(count for count, _ in taxonomy_data)
    print(f"Found {len(taxonomy_data)} unique taxonomy paths")
    print(f"Total sequences: {total_sequences}")

    # Parse and fill missing levels
    print("Processing taxonomy paths and filling missing levels...")
    taxonomy_paths_with_counts = []
    for count, tax_str in taxonomy_data:
        path = parse_taxonomy(tax_str, args.levels, args.lineage_reverted)
        if path:
            taxonomy_paths_with_counts.append((count, path))

    print(f"Processed {len(taxonomy_paths_with_counts)} valid taxonomy paths")

    if args.stop_level is not None:
        print(f"Truncating taxonomy paths at level: {args.stop_level}")
        taxonomy_paths_with_counts = truncate_to_stop_level(
            taxonomy_paths_with_counts, args.levels, args.stop_level
        )

    if args.top_hits is not None:
        print(f"Keeping top {args.top_hits} taxa per level, merging the rest into '{OTHER_LABEL}'")
        taxonomy_paths_with_counts = apply_top_hits(taxonomy_paths_with_counts, args.top_hits)

    # Build Sankey data
    print("Building Sankey diagram data...")
    nodes, links = build_sankey_data(taxonomy_paths_with_counts)

    print(f"Sankey diagram: {len(nodes)} nodes, {len(links)} links")

    # Verify root total
    root_total = sum(link['value'] for link in links if nodes[link['source']] == 'All sequences')
    print(f"Root total: {root_total} (should match total sequences: {total_sequences})")

    # Create plot
    print("Creating Sankey plot...")
    create_sankey_plotly(nodes, links, args.output, args.title)

    if args.output_table is not None:
        write_links_table(sankey_links_table(taxonomy_paths_with_counts, args.levels), args.output_table)

    print("Done!")


if __name__ == '__main__':
    main()
