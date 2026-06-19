#!/usr/bin/env python3
"""
Create a Sankey plot from taxonomy data with proper handling of missing taxonomic levels.

Handles missing/empty taxonomy levels by creating "unclassified_[level]_[parent_taxa]" labels
to ensure all links in the Sankey diagram are properly connected.
"""

import argparse
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
            # Node key includes full path to make it unique
            node_key = "|".join(current_path_parts)

            # Add node
            node_keys.add(node_key)
            node_key_to_label[node_key] = taxon

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
        font_size=20,
        height=800
    )

    fig.write_html(output_file)
    print(f"Sankey plot saved to: {output_file}")


def read_tsv_taxonomy(tsv_file, taxonomy_column='taxonomy'):
    """
    Read taxonomy data from TSV file.

    Supports two formats:
    1. Standard TSV with header and taxonomy column
    2. Krona format: count\trank1\trank2\t... (no header)

    Args:
        tsv_file: Path to TSV file
        taxonomy_column: Name of the taxonomy column (for standard TSV)

    Returns:
        List of (count, taxonomy_string) tuples
    """
    import csv

    taxonomy_data = []

    with open(tsv_file, 'r') as f:
        # Peek at first line to detect format
        first_line = f.readline().strip()
        if not first_line:  # File is completely empty
            print("File is empty")
            exit(0)
        f.seek(0)

        # Check if it's Krona format (first column is a number)
        first_parts = first_line.split('\t')
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
                parts = line.strip().split('\t')
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
            print(f"Detected standard TSV format with header")
            reader = csv.DictReader(f, delimiter='\t')

            if taxonomy_column not in reader.fieldnames:
                print(f"Error: Column '{taxonomy_column}' not found in TSV file.", file=sys.stderr)
                print(f"Available columns: {', '.join(reader.fieldnames)}", file=sys.stderr)
                sys.exit(1)

            for row in reader:
                tax = row.get(taxonomy_column, '')
                if tax and tax != 'NA':
                    taxonomy_data.append((1, tax))

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

    with open(gff_file, 'r') as f:
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
            if tax and tax != 'NA':
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

  # From GFF file directly
  %(prog)s --input viral.gff --input-format gff --output taxonomy_sankey.html

  # Specify custom taxonomy levels
  %(prog)s --input stats.tsv --output sankey.html --levels domain phylum class order family genus

Standard viral taxonomy levels (default):
  realm, kingdom, phylum, class, order, family, subfamily, genus, species

TSV Format Support:
  1. Standard TSV with header: taxonomy column contains semicolon-separated ranks
  2. Krona format (auto-detected): count\\trank1\\trank2\\trank3\\t... (no header)
        """
    )

    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        help='Input file (TSV with taxonomy column, Krona format, or GFF)'
    )

    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output HTML file for Sankey plot'
    )

    parser.add_argument(
        '--input-format',
        choices=['tsv', 'gff'],
        default='tsv',
        help='Input file format: tsv (auto-detects Krona format) or gff (default: tsv)'
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

    # Validate input file
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Read taxonomy data
    print(f"Reading taxonomy data from {args.input}...")
    if args.input_format == 'tsv':
        taxonomy_data = read_tsv_taxonomy(args.input, args.taxonomy_column)
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

    print("Done!")


if __name__ == '__main__':
    main()
