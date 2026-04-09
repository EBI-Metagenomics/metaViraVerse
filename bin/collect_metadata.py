#!/usr/bin/env python3

import argparse
import csv
import gzip
from urllib.parse import unquote


MISSING = "missing"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect and combine metadata for viral sequences and plasmids."
    )
    parser.add_argument("--combined_meta", required=True, help="Source metadata from choose_sequences step")
    parser.add_argument("--viruses_cluster", required=True, help="Viruses cluster information after vclust cluster step")
    parser.add_argument("--plasmids_cluster", required=True, help="Plasmids cluster information after vclust cluster step")
    parser.add_argument("--additional_metadata", required=True, help="Concatenated catalogues metadata tables with header")
    parser.add_argument("--viruses_vitap", required=True, help="Output of VITAP tool *_vitap_best.tsv")
    parser.add_argument("--gff", required=True, help="Concatenated GFF file containing all input sequences")
    parser.add_argument("--output_viruses", required=True, help="Output TSV for viral sequences and prophages")
    parser.add_argument("--output_plasmids", required=True, help="Output TSV for plasmids")
    parser.add_argument("--compress", action="store_true", help="Compress output files with gzip (.gz appended)")
    return parser.parse_args()


def load_cluster(path):
    """Return dict: description -> cluster_representative."""
    clusters = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            clusters[row["object"]] = row["cluster"]
    return clusters


def load_all_metadata(path):
    """Return dict: genome_id -> {Genome_accession, Source_lineage, Sample_accession, Study_accession}."""
    meta = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            meta[row["Genome"]] = {
                "Genome_accession": row.get("Genome_accession", MISSING) or MISSING,
                "Source_lineage": row.get("Lineage", MISSING) or MISSING,
                "Sample_accession": row.get("Sample_accession", MISSING) or MISSING,
                "Study_accession": row.get("Study_accession", MISSING) or MISSING,
            }
    return meta


def load_vitap(path):
    """Return dict: Genome_ID (description) -> lineage."""
    vitap = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            vitap[row["Genome_ID"]] = row.get("lineage", MISSING) or MISSING
    return vitap


def parse_gff_attributes(attr_str):
    """Parse GFF attribute string into a dict."""
    attrs = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            key, _, value = part.partition("=")
            attrs[key] = unquote(value)
    return attrs


def load_gff_taxonomy(path):
    """Return dict: description (from ID attribute) -> taxonomy."""
    taxonomy = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            source = fields[1]
            if source != "geNomad" and source != 'VIRify':
                continue
            attrs = parse_gff_attributes(fields[8])
            seq_id = attrs.get("ID", "")
            tax = attrs.get("taxonomy", "")
            if seq_id and tax:
                taxonomy[seq_id] = tax
    return taxonomy


def genome_key(description):
    """Extract genome ID from description: 'MGYG000517752_108|plasmid-1:...' -> 'MGYG000517752'."""
    return description.split("|")[0].split("_")[0]


def open_output(path, compress):
    """Open a text-mode writer, optionally gzip-compressed."""
    if compress:
        path = path + ".gz"
        return path, gzip.open(path, "wt", newline="")
    return path, open(path, "w", newline="")


def write_viruses(rows, clusters, all_meta, vitap, gff_tax, out_path, compress):
    fieldnames = [
        "Sequence_ID", "Cluster_rep",
        "vitap_lineage", "viphogs_lineage",
        "Source", "Biome", "Source_accession", "Source_lineage", "Source_sample", "Source_project",
        "Sequence_length", "Sequence_sha256"
    ]
    out_path, fh = open_output(out_path, compress)
    with fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            desc = row["description"]
            gid = genome_key(desc)
            meta = all_meta.get(gid, {})
            writer.writerow({
                "Sequence_ID": desc,
                "Cluster_rep": clusters.get(desc, MISSING),
                "vitap_lineage": vitap.get(desc, MISSING),
                "viphogs_lineage": gff_tax.get(desc, MISSING),
                "Source": row["type"],
                "Biome": row["biomes"],
                "Source_accession": meta.get("Genome_accession", MISSING),
                "Source_lineage": meta.get("Source_lineage", MISSING),
                "Source_sample": meta.get("Sample_accession", MISSING),
                "Source_project": meta.get("Study_accession", MISSING),
                "Sequence_length": row["sequence_length"],
                "Sequence_sha256": row["sequence_sha256"],
            })
    return out_path


def write_plasmids(rows, clusters, all_meta, out_path, compress):
    fieldnames = [
        "Sequence_ID", "Cluster_rep",
        "Source", "Biome", "Source_accession", "Source_lineage", "Source_sample", "Source_project",
        "Sequence_length", "Sequence_sha256"
    ]
    out_path, fh = open_output(out_path, compress)
    with fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            desc = row["description"]
            gid = genome_key(desc)
            meta = all_meta.get(gid, {})
            writer.writerow({
                "Sequence_ID": desc,
                "Cluster_rep": clusters.get(desc, MISSING),
                "Source": row["type"],
                "Biome": row["biomes"],
                "Source_accession": meta.get("Genome_accession", MISSING),
                "Source_lineage": meta.get("Source_lineage", MISSING),
                "Source_sample": meta.get("Sample_accession", MISSING),
                "Source_project": meta.get("Study_accession", MISSING),
                "Sequence_length": row["sequence_length"],
                "Sequence_sha256": row["sequence_sha256"],
            })
    return out_path


def main():
    args = parse_args()

    print("Loading cluster files...")
    virus_clusters = load_cluster(args.viruses_cluster)
    plasmid_clusters = load_cluster(args.plasmids_cluster)

    print("Loading genome metadata...")
    all_meta = load_all_metadata(args.additional_metadata)

    print("Loading ViTAP taxonomy...")
    vitap = load_vitap(args.viruses_vitap)

    print("Loading GFF taxonomy...")
    gff_tax = load_gff_taxonomy(args.gff)

    print("Splitting combined metadata into viruses and plasmids...")
    virus_rows = []
    plasmid_rows = []
    with open(args.combined_meta) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            desc = row["description"]
            if "plasmid" in desc:
                plasmid_rows.append(row)
            elif "viral_sequence" in desc or "prophage" in desc:
                virus_rows.append(row)

    print(f"  Viruses/prophages: {len(virus_rows)}")
    print(f"  Plasmids: {len(plasmid_rows)}")

    print(f"Writing {args.output_viruses}...")
    written = write_viruses(virus_rows, virus_clusters, all_meta, vitap, gff_tax, args.output_viruses, args.compress)
    print(f"  -> {written}")

    print(f"Writing {args.output_plasmids}...")
    written = write_plasmids(plasmid_rows, plasmid_clusters, all_meta, args.output_plasmids, args.compress)
    print(f"  -> {written}")

    print("Done.")


if __name__ == "__main__":
    main()
