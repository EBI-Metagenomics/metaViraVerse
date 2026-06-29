#!/usr/bin/env python3
"""Extract GFF records, protein IDs, and a plain rep-ID list for cluster representatives.

Stats (taxonomy, checkV, cluster-level aggregates) are now produced by
collect_metadata.py which has access to all structured metadata without
needing to parse GFF attributes.
"""
import argparse
import os
import sys
from utils import parse_attributes, read_input_gff


def read_mapfile(mapfile):
    mapping = {}
    with open(mapfile, "r") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                mapping[parts[1]] = parts[0]   # temporary -> original
    return mapping


def read_cluster_reps(viral_list_file, mapping=None):
    """Return (cluster_reps, rep_original_names).

    cluster_reps      : list of mapped rep IDs (in order of first appearance)
    rep_original_names: dict mapped_rep_id -> original_rep_id
    """
    cluster_reps = []
    rep_original_names = {}
    seen = set()
    input_format = None

    with open(viral_list_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if input_format is None:
                if "object" in line and "cluster" in line:
                    input_format = "vclust"
                    continue
                else:
                    input_format = "blastn"

            parts = line.split("\t")
            if input_format == "vclust":
                original_rep = parts[1].strip()
            else:
                original_rep = parts[0].strip()

            rep_id = mapping[original_rep] if (mapping and original_rep in mapping) else original_rep

            if rep_id not in seen:
                seen.add(rep_id)
                cluster_reps.append(rep_id)
                rep_original_names[rep_id] = original_rep

    return cluster_reps, rep_original_names


def extract_proteins(lines):
    protein_ids = set()
    for line in lines:
        cols = line.strip().split("\t")
        if len(cols) >= 9 and cols[2] == "CDS":
            attrs, _ = parse_attributes(cols[8])
            pid = attrs.get("ID")
            if pid:
                protein_ids.add(pid)
    return protein_ids


def extract_gff_and_proteins(viral_list_file, gff_file, output_reps_list,
                              output_reps_gff, output_proteins, mapfile):
    mapping = read_mapfile(mapfile) if mapfile else None
    cluster_reps, rep_original_names = read_cluster_reps(viral_list_file, mapping)

    print(f"Found {len(cluster_reps)} cluster representatives")

    input_gff, attr_id_to_seq_id, _ = read_input_gff([gff_file])

    proteins = set()
    written_gff = set()
    found = 0

    with open(output_reps_gff, "w") as reps_gff:
        reps_gff.write("##gff-version 3\n")
        for rep_id in cluster_reps:
            original_rep = rep_original_names.get(rep_id, rep_id)
            mgyg_id = attr_id_to_seq_id.get(original_rep)
            if mgyg_id and mgyg_id in input_gff:
                lines = input_gff[mgyg_id]
                if mgyg_id not in written_gff:
                    reps_gff.writelines(lines)
                    written_gff.add(mgyg_id)
                proteins.update(extract_proteins(lines))
                found += 1
            else:
                print(f"No GFF data for rep {rep_id} (original: {original_rep})")

    print(f"Extracted GFF records for {found} representatives")

    if output_reps_list:
        with open(output_reps_list, "w") as out:
            for rep_id in cluster_reps:
                out.write(f"{rep_original_names.get(rep_id, rep_id)}\n")

    if output_proteins:
        with open(output_proteins, "w") as out:
            for prot_id in proteins:
                out.write(f"{prot_id}\n")

    print(f"GFF written to: {output_reps_gff}  ({found} sequences)")
    if output_reps_list:
        print(f"Reps list written to: {output_reps_list}  ({len(cluster_reps)} reps)")
    if output_proteins:
        print(f"Proteins list written to: {output_proteins}  ({len(proteins)} proteins)")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Extract GFF records and protein IDs for cluster representatives. "
            "Cluster stats (taxonomy, checkV, aggregates) are produced by collect_metadata.py."
        )
    )
    parser.add_argument("--viral-list", required=True,
                        help="Cluster TSV (vclust: object/cluster columns; blastn: rep\\tmembers)")
    parser.add_argument("--gff", required=True,
                        help="Full GFF file annotated by VIRify or similar tool")
    parser.add_argument("--mapfile", required=False,
                        help="Map file from rename_contigs step (temporary -> original)")
    parser.add_argument("--output-reps-list", required=False,
                        help="Output file with one representative ID per line")
    parser.add_argument("--output-reps-gff", required=True,
                        help="Output GFF containing only records for cluster representatives")
    parser.add_argument("--output-reps-proteins", required=True,
                        help="Output file with protein IDs for cluster representatives")
    return parser.parse_args()


def main():
    args = parse_arguments()

    for f in [args.viral_list, args.gff]:
        if not os.path.exists(f):
            print(f"Error: File not found: {f}", file=sys.stderr)
            sys.exit(1)

    extract_gff_and_proteins(
        args.viral_list,
        args.gff,
        args.output_reps_list,
        args.output_reps_gff,
        args.output_reps_proteins,
        args.mapfile,
    )


if __name__ == "__main__":
    main()
