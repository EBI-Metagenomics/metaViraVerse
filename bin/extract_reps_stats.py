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

    Format:
    for blast result: rep_id\tmembers
    for vclust: member_id\trep_id

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
    input_format = None
    with open(viral_list_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Detect format from the first line only
            if input_format is None:
                if 'object' in line and 'cluster' in line:
                    input_format = 'vclust'
                    continue  # skip the header line
                else:
                    input_format = 'blastn'

            if input_format == 'vclust':
                parts = line.split('\t')
                original_rep_id = parts[1].strip()
                rep_id = original_rep_id
                original_member_id = parts[0].strip()
                member_id = original_member_id

                # Apply mapping if provided
                if mapping and rep_id in mapping:
                    rep_id = mapping[rep_id]
                if mapping and member_id in mapping:
                    member_id = mapping[member_id]

                # Track unique reps
                if rep_id not in cluster_members:
                    cluster_reps.append(rep_id)
                    cluster_members[rep_id] = [member_id]
                    rep_original_names[rep_id] = original_rep_id
                else:
                    cluster_members[rep_id].append(member_id)
            else:
                parts = line.split('\t')
                original_rep_id = parts[0].strip()
                rep_id = original_rep_id
                if mapping and rep_id in mapping:
                    rep_id = mapping[rep_id]

                if rep_id not in cluster_members:
                    cluster_reps.append(rep_id)
                    cluster_members[rep_id] = []
                    rep_original_names[rep_id] = original_rep_id

                if len(parts) >= 2:
                    for member_id in parts[1].split(' '):
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



def parse_taxonomy(taxonomy_str, levels=None):
    """Parse a semicolon-separated taxonomy string into a list of level values.

    URL-encoded separators (%3B) are decoded before splitting.  Empty/NA input
    returns an empty list.  When ``levels`` is provided, empty positions are
    filled with ``Unclassified_{level}_{parent.lower()}`` placeholders.
    Without ``levels``, trailing empty tokens are stripped.

    Args:
        taxonomy_str: Semicolon-separated taxonomy string (may be URL-encoded).
        levels: Optional list of level names (e.g. ['realm', 'kingdom', ...]).

    Returns:
        List of taxonomy level strings.
    """
    if not taxonomy_str or taxonomy_str == 'NA':
        return []
    parts = unquote(taxonomy_str).split(';')
    if levels is None:
        return [p for p in parts if p.strip()]
    result = []
    last_filled = ''
    for i, part in enumerate(parts):
        stripped = part.strip()
        if stripped:
            result.append(stripped)
            last_filled = stripped
        else:
            level_name = levels[i] if i < len(levels) else f'level{i}'
            result.append(f'Unclassified_{level_name}_{last_filled.lower()}')
    return result


def generate_krona_file(results, viral_names, output_file):
    """Write a Krona-compatible TSV from taxonomy strings in results.

    Each unique taxonomy path is written as one line:
    ``count<TAB>level1<TAB>level2...``
    Lines are sorted by count descending.  Sequences with empty/NA taxonomy
    are omitted.

    Args:
        results: Dict mapping seq_id -> {'taxonomy': str, ...}.
        viral_names: Sequence IDs to include.
        output_file: Path for the output TSV file.
    """
    counts: Counter = Counter()
    for name in viral_names:
        tax_str = results.get(name, {}).get('taxonomy', 'NA')
        tax_parts = parse_taxonomy(tax_str)
        if tax_parts:
            counts[tuple(tax_parts)] += 1
    with open(output_file, 'w') as f:
        for tax_tuple, count in sorted(counts.items(), key=lambda x: -x[1]):
            f.write(f"{count}\t{chr(9).join(tax_tuple)}\n")


def extract_viral_data(viral_list_file, gff_file, output_file, output_reps, output_gff, output_proteins, mapfile):
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

    # example, MGYG000517142_2|prophage-311953:328572 into MGYG000517142_2 to match with initial GFF
    reps_names_without_viral_region = set([i.split('|')[0] for i in cluster_reps])

    # Get all unique sequence IDs (reps + all members)
    all_seq_ids = set(cluster_reps)
    for members in cluster_members.values():
        all_seq_ids.update(members)

    print(f"📋 Found {len(cluster_reps)} cluster representatives")
    print(f"📋 Total sequences (reps + members): {len(all_seq_ids)}")

    # Extract data from GFF for all sequences
    all_results = {}
    found = set()
    proteins = set()

    with open(gff_file, "r") as gff, open(output_gff, 'w') as reps_gff:
        reps_gff.write('##gff-version 3\n')

        for line in gff:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            if cols[0] in reps_names_without_viral_region:
                reps_gff.write(line)

            attr_str = cols[8]
            attrs = parse_attributes(attr_str)

            if cols[2] == 'CDS':
                if cols[0] in reps_names_without_viral_region:
                    prot_id = attrs.get("id")
                    if prot_id:
                        proteins.add(prot_id)
                continue

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

    # Write a list of representatives
    if output_reps:
        with open(output_reps, "w") as out:
            for rep_id in cluster_reps:
                out.write(f'{rep_id}\n')

    # Write a list of proteins for representatives
    if output_proteins:
        with open(output_proteins, "w") as out:
            for prot_id in proteins:
                out.write(f'{prot_id}\n')

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
        "--output-reps-list",
        required=False,
        help="Output file with representatives line separated"
    )
    parser.add_argument(
        "--output-reps-gff",
        required=True,
        help="Output file with GFF records for cluster representatives"
    )
    parser.add_argument(
        "--output-reps-proteins",
        required=True,
        help="Output FASTA file with proteins for cluster representatives"
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
        args.output_reps_list,
        args.output_reps_gff,
        args.output_reps_proteins,
        args.mapfile
    )

if __name__ == "__main__":
    main()
