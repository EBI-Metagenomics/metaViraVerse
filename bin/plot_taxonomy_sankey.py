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


def parse_taxonomy(taxonomy_str, standard_levels=None):
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

    # URL decode and split
    taxonomy_str = unquote(taxonomy_str)
    parts = [p.strip() for p in taxonomy_str.split(';')]

    # If no standard levels provided, use generic numbering
    if standard_levels is None:
        standard_levels = [f'level_{i}' for i in range(len(parts))]

    filled_taxonomy = []
    last_valid = "root"

    for i, part in enumerate(parts):
        level_name = standard_levels[i] if i < len(standard_levels) else f'level_{i}'

        if part and part.strip():
            # Valid taxonomy term
            filled_taxonomy.append(part.strip())
            last_valid = part.strip()
        else:
            # Missing level - create unclassified label
            unclassified_label = f"unclassified_{level_name}_{last_valid}"
            filled_taxonomy.append(unclassified_label)

    return filled_taxonomy


def build_sankey_data(taxonomy_list):
    """
    Build Sankey diagram data from list of taxonomy paths.

    Args:
        taxonomy_list: List of taxonomy paths (each path is a list of terms)

    Returns:
        Tuple of (nodes, links) where:
        - nodes is a list of unique taxonomy terms
        - links is a list of (source_idx, target_idx, value) tuples
    """
    # Count edges
    edge_counts = defaultdict(int)

    for taxonomy_path in taxonomy_list:
        if len(taxonomy_path) < 2:
            continue

        for i in range(len(taxonomy_path) - 1):
            source = taxonomy_path[i]
            target = taxonomy_path[i + 1]
            edge_counts[(source, target)] += 1

    # Build unique nodes list
    all_nodes = set()
    for taxonomy_path in taxonomy_list:
        all_nodes.update(taxonomy_path)

    nodes = sorted(list(all_nodes))
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}

    # Build links
    links = []
    for (source, target), count in edge_counts.items():
        links.append({
            'source': node_to_idx[source],
            'target': node_to_idx[target],
            'value': count
        })

    return nodes, links


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


def read_tsv_taxonomy(tsv_file, taxonomy_column='taxonomy'):
    """
    Read taxonomy data from TSV file.

    Args:
        tsv_file: Path to TSV file
        taxonomy_column: Name of the taxonomy column

    Returns:
        List of taxonomy strings
    """
    import csv

    taxonomies = []
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')

        if taxonomy_column not in reader.fieldnames:
            print(f"Error: Column '{taxonomy_column}' not found in TSV file.", file=sys.stderr)
            print(f"Available columns: {', '.join(reader.fieldnames)}", file=sys.stderr)
            sys.exit(1)

        for row in reader:
            tax = row.get(taxonomy_column, '')
            if tax and tax != 'NA':
                taxonomies.append(tax)

    return taxonomies


def read_gff_taxonomy(gff_file):
    """
    Read taxonomy data directly from GFF file.

    Args:
        gff_file: Path to GFF file

    Returns:
        List of taxonomy strings
    """
    taxonomies = []

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
                taxonomies.append(tax)

    return taxonomies


def main():
    parser = argparse.ArgumentParser(
        description='Create Sankey plot from viral taxonomy data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From TSV file (output of extract_reps_stats.py)
  %(prog)s --input stats.tsv --output taxonomy_sankey.html

  # From GFF file directly
  %(prog)s --input viral.gff --input-format gff --output taxonomy_sankey.html

  # Specify custom taxonomy levels
  %(prog)s --input stats.tsv --output sankey.html --levels domain phylum class order family genus

Standard viral taxonomy levels (default):
  realm, kingdom, phylum, class, order, family, subfamily, genus, species
        """
    )

    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        help='Input file (TSV or GFF)'
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
        help='Input file format (default: tsv)'
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

    args = parser.parse_args()

    # Validate input file
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Read taxonomy data
    print(f"Reading taxonomy data from {args.input}...")
    if args.input_format == 'tsv':
        taxonomy_strings = read_tsv_taxonomy(args.input, args.taxonomy_column)
    else:
        taxonomy_strings = read_gff_taxonomy(args.input)

    if not taxonomy_strings:
        print("Error: No taxonomy data found in input file", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(taxonomy_strings)} taxonomy records")

    # Parse and fill missing levels
    print("Processing taxonomy paths and filling missing levels...")
    taxonomy_paths = []
    for tax_str in taxonomy_strings:
        path = parse_taxonomy(tax_str, args.levels)
        if path:
            taxonomy_paths.append(path)

    print(f"Processed {len(taxonomy_paths)} valid taxonomy paths")

    # Build Sankey data
    print("Building Sankey diagram data...")
    nodes, links = build_sankey_data(taxonomy_paths)

    print(f"Sankey diagram: {len(nodes)} nodes, {len(links)} links")

    # Create plot
    print("Creating Sankey plot...")
    create_sankey_plotly(nodes, links, args.output, args.title)

    print("Done!")


if __name__ == '__main__':
    main()
