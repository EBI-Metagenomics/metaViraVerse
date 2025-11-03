#!/usr/bin/env python3
import argparse
import os
import sys
from collections import Counter
from urllib.parse import unquote


def parse_attributes(attr_str):
    """Parse the attributes column from a GFF line into a dictionary."""
    attrs = {}
    for part in attr_str.split(";"):
        if "=" in part:
            key, val = part.split("=", 1)
            attrs[key.strip().lower()] = val.strip()
    return attrs


def read_mapfile(mapfile):
    mapping = {}
    with open(mapfile, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')
            tmp_name = line[1]
            original_name = line[0]
            mapping[tmp_name] = original_name
    return mapping


def read_cluster_structure(viral_list_file, mapping=None):
    """
    Read cluster structure from viral list file.

    Format: rep_id\tmember_id (one member per line, or rep only if single member)

    Args:
        viral_list_file: File with cluster representatives and members
        mapping: Optional contig name mapping

    Returns:
        Tuple of (cluster_reps, cluster_members, rep_original_names) where:
        - cluster_reps: list of cluster representative IDs (unique, after mapping)
        - cluster_members: dict mapping rep_id -> list of member IDs (including rep itself)
        - rep_original_names: dict mapping mapped rep_id -> original rep_id
    """
    cluster_reps = []
    cluster_members = {}
    rep_original_names = {}
    seen_reps = set()

    with open(viral_list_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            original_rep_id = parts[0].strip()
            rep_id = original_rep_id

            # Apply mapping if provided
            if mapping and rep_id in mapping:
                rep_id = mapping[rep_id]

            # Track unique reps
            if rep_id not in seen_reps:
                cluster_reps.append(rep_id)
                cluster_members[rep_id] = []
                rep_original_names[rep_id] = original_rep_id
                seen_reps.add(rep_id)

            # Check if there's a second column with member IDs
            if len(parts) > 1 and parts[1].strip():
                member_id = parts[1].strip()
                # Apply mapping if provided
                if mapping and member_id in mapping:
                    member_id = mapping[member_id]
                cluster_members[rep_id].append(member_id)

    return cluster_reps, cluster_members, rep_original_names


def calculate_mean_genes(cluster_members, all_results):
    """
    Calculate mean number of viral genes per cluster.

    Args:
        cluster_members: dict mapping rep_id -> list of member IDs
        all_results: dict mapping seq_id -> data with checkv_viral_genes

    Returns:
        dict mapping rep_id -> mean_genes
    """
    cluster_mean_genes = {}

    for rep_id, members in cluster_members.items():
        gene_counts = []

        # Include the representative itself if it has data
        if rep_id in all_results:
            genes = all_results[rep_id].get("checkv_viral_genes", "NA")
            if genes != "NA" and genes.strip():
                try:
                    gene_counts.append(float(genes))
                except ValueError:
                    pass

        # Include all cluster members
        for member_id in members:
            if member_id in all_results:
                genes = all_results[member_id].get("checkv_viral_genes", "NA")
                if genes != "NA" and genes.strip():
                    try:
                        gene_counts.append(float(genes))
                    except ValueError:
                        pass

        # Calculate mean
        if gene_counts:
            cluster_mean_genes[rep_id] = sum(gene_counts) / len(gene_counts)
        else:
            cluster_mean_genes[rep_id] = None

    return cluster_mean_genes


def parse_taxonomy(taxonomy_str, standard_levels=None):
    """
    Parse taxonomy string and fill in missing levels for Krona plot.

    Args:
        taxonomy_str: Semicolon-separated taxonomy string
        standard_levels: List of standard taxonomic level names

    Returns:
        List of taxonomy terms with missing levels filled as "unclassified_[level]_[parent]"
    """
    if not taxonomy_str or taxonomy_str == "NA":
        return []

    # URL decode and split
    taxonomy_str = unquote(taxonomy_str)
    parts = [p.strip() for p in taxonomy_str.split(';')]

    # Default viral taxonomy levels
    if standard_levels is None:
        standard_levels = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'species']

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


def generate_krona_file(results, viral_names, krona_output_file, standard_levels=None):
    """
    Generate Krona plot input file from taxonomy data.

    Format: count\ttaxonomy_rank1\ttaxonomy_rank2\t...

    Args:
        results: Dictionary of sequence_id -> data with taxonomy
        viral_names: List of all viral sequence names
        krona_output_file: Output file path for Krona format
        standard_levels: List of standard taxonomic level names
    """
    # Collect all taxonomy paths
    taxonomy_paths = []
    for seq in viral_names:
        if seq in results:
            taxonomy_str = results[seq].get("taxonomy", "NA")
            if taxonomy_str and taxonomy_str != "NA":
                path = parse_taxonomy(taxonomy_str, standard_levels)
                if path:
                    taxonomy_paths.append(tuple(path))

    # Count occurrences of each unique taxonomy path
    taxonomy_counts = Counter(taxonomy_paths)

    # Write Krona format
    with open(krona_output_file, "w") as out:
        for taxonomy_path, count in sorted(taxonomy_counts.items(), key=lambda x: -x[1]):
            # Format: count\trank1\trank2\trank3\t...
            taxonomy = '\t'.join(taxonomy_path)
            out.write(f"{count}\t{taxonomy}\n")

    print(f"📊 Krona file written to: {krona_output_file}")
    print(f"   Total unique taxonomy paths: {len(taxonomy_counts)}")


def extract_viral_data(viral_list_file, gff_file, output_file, mapfile):
    """
    Extract taxonomy, checkv_viral_genes, and checkv_quality for viral sequences found in GFF.

    Args:
        viral_list_file: File with viral sequence names (column 1) and cluster members (column 2)
        gff_file: GFF file with annotations
        output_file: Output TSV file
        mapfile: Optional mapping file for renamed contigs

    Returns:
        Tuple of (results dict, cluster_reps list) for further processing
    """
    mapping = None
    if mapfile:
        mapping = read_mapfile(mapfile)

    # Read cluster structure
    cluster_reps, cluster_members, rep_original_names = read_cluster_structure(viral_list_file, mapping)

    # Get all unique sequence IDs (reps + all members)
    all_seq_ids = set(cluster_reps)
    for members in cluster_members.values():
        all_seq_ids.update(members)

    print(f"📋 Found {len(cluster_reps)} cluster representatives")
    print(f"📋 Total sequences (reps + members): {len(all_seq_ids)}")

    # Extract data from GFF for all sequences
    all_results = {}
    found = set()

    with open(gff_file, "r") as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            attr_str = cols[8]
            attrs = parse_attributes(attr_str)
            seq_id = attrs.get("id", "NA")

            if seq_id in all_seq_ids:
                all_results[seq_id] = {
                    "taxonomy": attrs.get("taxonomy", "NA"),
                    "checkv_viral_genes": attrs.get("checkv_viral_genes", "NA"),
                    "checkv_quality": attrs.get("checkv_quality", "NA"),
                    "virify_quality": attrs.get("virify_quality", "NA"),
                }
                found.add(seq_id)

    print(f"✅ Extracted stats for {len(found)} sequences from GFF")

    # Calculate mean genes per cluster
    cluster_mean_genes = calculate_mean_genes(cluster_members, all_results)

    # Write output for cluster representatives only
    with open(output_file, "w") as out:
        out.write("viral_sequence_name\toriginal_name\ttaxonomy\tcheckv_viral_genes\tcheckv_quality\tvirify_quality\tmean_cluster_genes\n")
        for rep_id in cluster_reps:
            data = all_results.get(rep_id, {
                "taxonomy": "NA",
                "checkv_viral_genes": "NA",
                "checkv_quality": "NA",
                "virify_quality": "NA"
            })

            # Get original name
            original_name = rep_original_names.get(rep_id, rep_id)

            # Format mean genes
            mean_genes = cluster_mean_genes.get(rep_id)
            mean_genes_str = f"{mean_genes:.2f}" if mean_genes is not None else "NA"

            out.write(
                f"{rep_id}\t{original_name}\t{data['taxonomy']}\t{data['checkv_viral_genes']}\t"
                f"{data['checkv_quality']}\t{data['virify_quality']}\t{mean_genes_str}\n"
            )

    print(f"📄 Output written to: {output_file}")
    print(f"   Cluster representatives: {len(cluster_reps)}")

    return all_results, cluster_reps


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract cluster reps stats from GFF file and optionally generate Krona plot format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  %(prog)s --viral-list cluster_reps.txt --gff virify.gff --output stats.tsv --krona krona.txt

Input format for --viral-list:
  Column 1: Cluster representative ID
  Column 2: Cluster member ID (optional, one member per line)

  Example:
  rep1\tmember1
  rep1\tmember2
  rep2\tmember3
        """
    )
    parser.add_argument(
        "--viral-list",
        required=True,
        help="Path to file with cluster representatives (column 1) and members (column 2)."
    )
    parser.add_argument(
        "--gff",
        required=True,
        help="Path to GFF file annotated by VIRify or similar tool."
    )
    parser.add_argument(
        "--mapfile",
        required=False,
        help="Map-file as product of renaming contigs step"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file with stats including mean cluster genes."
    )
    parser.add_argument(
        "--krona",
        required=False,
        help="Optional output file for Krona plot format (count\\ttaxonomy_ranks tab-separated)"
    )

    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_arguments()

    # Check file existence
    for f in [args.viral_list, args.gff]:
        if not os.path.exists(f):
            print(f"Error: File not found: {f}", file=sys.stderr)
            sys.exit(1)

    # Extract viral data
    all_results, cluster_reps = extract_viral_data(
        args.viral_list,
        args.gff,
        args.output,
        args.mapfile
    )

    # Generate Krona file if requested
    if args.krona:
        generate_krona_file(all_results, cluster_reps, args.krona)


if __name__ == "__main__":
    main()
