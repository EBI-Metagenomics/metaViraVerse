#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import defaultdict


MISSING = "missing"

CHECKV_FIELDS = [
    "checkv_quality",
    "miuvig_quality",
    "completeness",
    "completeness_method",
    "contamination",
    "provirus",
    "proviral_length",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect and combine metadata for viral sequences and plasmids."
    )
    parser.add_argument("--combined_meta", required=True,
                        help="Source metadata from choose_sequences step (includes checkV fields)")
    parser.add_argument("--viruses_cluster", required=True,
                        help="Viruses cluster TSV from vclust (columns: object, cluster)")
    parser.add_argument("--plasmids_cluster", required=True,
                        help="Plasmids cluster TSV from vclust (columns: object, cluster)")
    parser.add_argument("--additional_metadata", required=True,
                        help="Concatenated catalogue metadata tables with header")
    parser.add_argument("--viruses_vitap", required=True,
                        help="VITAP output *_vitap_best.tsv (columns: Genome_ID, lineage, ...)")
    parser.add_argument("--viphogs_taxonomy", required=True,
                        help="Per-contig ViPhOGs taxonomy TSV (columns: contig_ID, superkingdom, ...)")
    parser.add_argument("--genomad", required=True,
                        help="geNomad virus summary TSV (columns: seq_name, taxonomy, ...)")
    parser.add_argument("--map", required=True,
                        help="Rename map TSV from rename_contigs (columns: original, temporary, short, biome, type)")
    parser.add_argument("--output_viruses", required=True,
                        help="Output TSV for viral sequences and prophages")
    parser.add_argument("--output_plasmids", required=True,
                        help="Output TSV for plasmids")
    parser.add_argument("--output_reps", required=True,
                        help="Output TSV for virus cluster representatives with aggregated stats")
    parser.add_argument("--compress", action="store_true",
                        help="Compress output files with gzip (.gz appended)")
    return parser.parse_args()


# ── Loaders ───────────────────────────────────────────────────────────────────

def load_cluster(path):
    """Return (member_to_rep, rep_to_members).

    member_to_rep : dict  description -> rep_description
    rep_to_members: dict  rep_description -> [member_descriptions]  (rep is included)
    """
    member_to_rep = {}
    rep_to_members = defaultdict(list)
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            obj = row["object"]
            rep = row["cluster"]
            member_to_rep[obj] = rep
            rep_to_members[rep].append(obj)
    return member_to_rep, dict(rep_to_members)


def load_all_metadata(path):
    """Return dict: genome_id -> {Genome_accession, Source_lineage, Sample_accession, Study_accession}."""
    meta = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            meta[row["Genome"]] = {
                "Genome_accession": row.get("Genome_accession", MISSING) or MISSING,
                "Source_lineage":   row.get("Lineage", MISSING) or MISSING,
                "Sample_accession": row.get("Sample_accession", MISSING) or MISSING,
                "Study_accession":  row.get("Study_accession", MISSING) or MISSING,
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


def load_viphogs_taxonomy(path):
    """Return dict: contig_ID (MGYV catalogue ID) -> semicolon-joined lineage string."""
    taxonomy = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        reader.fieldnames = [f.strip() for f in reader.fieldnames]
        ranks = [f for f in reader.fieldnames if f != "contig_ID"]
        for row in reader:
            contig_id = row["contig_ID"].strip()
            lineage = ";".join(row.get(r, "").strip() for r in ranks)
            taxonomy[contig_id] = lineage if any(
                row.get(r, "").strip() for r in ranks
            ) else MISSING
    return taxonomy


def load_genomad(path):
    """Return dict: seq_name (MGYV catalogue ID) -> taxonomy string."""
    genomad = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genomad[row["seq_name"]] = row.get("taxonomy", MISSING) or MISSING
    return genomad


def load_map(path):
    """Return two dicts built from the rename map TSV.

    map.original format: "MGYG..._N type|chunk:coords"  (space before type)
    description format:  "MGYG..._N|type-chunk:coords"  (pipe before type, dash before chunk)

    Returns:
        desc_to_mgyv : dict  description -> temporary MGYV ID
        mgyv_to_desc : dict  temporary MGYV ID -> description
    """
    desc_to_mgyv = {}
    mgyv_to_desc = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            original = row["original"]
            mgyv = row["temporary"]
            if " " in original:
                genome, rest = original.split(" ", 1)
                rest_norm = rest.replace("|", "-", 1)
                desc = f"{genome}|{rest_norm}"
            else:
                desc = original
            desc_to_mgyv[desc] = mgyv
            mgyv_to_desc[mgyv] = desc
    return desc_to_mgyv, mgyv_to_desc


def load_combined_meta(path):
    """Return dict: description -> full row dict."""
    rows = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            desc = row.get("description", "").strip()
            if desc:
                rows[desc] = row
    return rows


# ── Helpers ───────────────────────────────────────────────────────────────────

def genome_key(description):
    """'MGYG000517752_108|plasmid-1:...' -> 'MGYG000517752'."""
    return description.split("|")[0].split("_")[0]


def open_output(path, compress):
    if compress:
        path = path + ".gz"
        return path, gzip.open(path, "wt", newline="")
    return path, open(path, "w", newline="")


def _safe_float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


# ── Writers ───────────────────────────────────────────────────────────────────

def write_viruses(combined_meta_rows, member_to_rep, all_meta, vitap,
                  viphogs_tax, genomad, desc_to_mgyv,
                  out_path, compress):
    fieldnames = [
        "Sequence_ID", "Cluster_rep",
        "vitap_lineage", "viphogs_lineage", "genomad_lineage",
        "Source", "Biome", "Source_accession", "Source_lineage", "Source_sample", "Source_project",
        "Sequence_length", "Sequence_sha256",
    ] + CHECKV_FIELDS

    out_path, fh = open_output(out_path, compress)
    with fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for desc, row in combined_meta_rows.items():
            gid = genome_key(desc)
            meta = all_meta.get(gid, {})
            mgyv = desc_to_mgyv.get(desc, "")
            checkv = {f: row.get(f, MISSING) or MISSING for f in CHECKV_FIELDS}
            writer.writerow({
                "Sequence_ID":      desc,
                "Cluster_rep":      member_to_rep.get(desc, MISSING),
                "vitap_lineage":    vitap.get(desc, MISSING),
                "viphogs_lineage":  viphogs_tax.get(mgyv, MISSING) if mgyv else MISSING,
                "genomad_lineage":  genomad.get(mgyv, MISSING) if mgyv else MISSING,
                "Source":           row.get("type", MISSING),
                "Biome":            row.get("biomes", MISSING),
                "Source_accession": meta.get("Genome_accession", MISSING),
                "Source_lineage":   meta.get("Source_lineage", MISSING),
                "Source_sample":    meta.get("Sample_accession", MISSING),
                "Source_project":   meta.get("Study_accession", MISSING),
                "Sequence_length":  row.get("sequence_length", MISSING),
                "Sequence_sha256":  row.get("sequence_sha256", MISSING),
                **checkv,
            })
    return out_path


def write_plasmids(combined_meta_rows, member_to_rep, all_meta,
                   out_path, compress):
    fieldnames = [
        "Sequence_ID", "Cluster_rep",
        "Source", "Biome", "Source_accession", "Source_lineage", "Source_sample", "Source_project",
        "Sequence_length", "Sequence_sha256",
    ] + CHECKV_FIELDS

    out_path, fh = open_output(out_path, compress)
    with fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for desc, row in combined_meta_rows.items():
            gid = genome_key(desc)
            meta = all_meta.get(gid, {})
            checkv = {f: row.get(f, MISSING) or MISSING for f in CHECKV_FIELDS}
            writer.writerow({
                "Sequence_ID":      desc,
                "Cluster_rep":      member_to_rep.get(desc, MISSING),
                "Source":           row.get("type", MISSING),
                "Biome":            row.get("biomes", MISSING),
                "Source_accession": meta.get("Genome_accession", MISSING),
                "Source_lineage":   meta.get("Source_lineage", MISSING),
                "Source_sample":    meta.get("Sample_accession", MISSING),
                "Source_project":   meta.get("Study_accession", MISSING),
                "Sequence_length":  row.get("sequence_length", MISSING),
                "Sequence_sha256":  row.get("sequence_sha256", MISSING),
                **checkv,
            })
    return out_path


def write_virus_reps_stats(virus_rows, rep_to_members, all_meta, vitap,
                           viphogs_tax, genomad, desc_to_mgyv, mgyv_to_desc,
                           out_path, compress):
    """One row per cluster representative with per-rep and per-cluster aggregated stats.

    Per-rep fields: all taxonomy lineages, checkV quality fields, sequence length.
    Per-cluster aggregates:
      - cluster_size          : total members including rep
      - mean_viral_genes      : mean viral_genes across members with data
      - cluster_biomes        : sorted unique biomes across all members
      - cluster_types         : sorted unique sequence types (viral_sequence / prophage)
      - members_with_checkv   : members with a non-missing checkv_quality
      - cluster_completeness_range : min-max completeness across members (numeric values only)
    """
    fieldnames = [
        "Sequence_ID", "Sequence_ID_catalogue",
        "vitap_lineage", "viphogs_lineage", "genomad_lineage",
        "Source", "Biome",
        "Sequence_length",
        "cluster_size", "mean_viral_genes", "members_with_checkv",
        "cluster_biomes", "cluster_types",
        "cluster_completeness_range",
    ] + CHECKV_FIELDS

    out_path, fh = open_output(out_path, compress)
    with fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for rep_desc, members in sorted(rep_to_members.items()):
            rep_mgyv = desc_to_mgyv.get(rep_desc, "")
            rep_row = virus_rows.get(rep_desc, {})

            # ── per-rep fields ────────────────────────────────────────────────
            checkv = {f: rep_row.get(f, MISSING) or MISSING for f in CHECKV_FIELDS}

            # ── cluster aggregates ────────────────────────────────────────────
            cluster_size = len(members)
            viral_genes_vals = []
            biomes = set()
            types = set()
            members_with_checkv = 0
            completeness_vals = []

            for member_desc in members:
                mrow = virus_rows.get(member_desc, {})
                if not mrow:
                    # member may be a MGYV ID if cluster file uses catalogue IDs
                    alt_desc = mgyv_to_desc.get(member_desc, "")
                    mrow = virus_rows.get(alt_desc, {})

                vg = _safe_float(mrow.get("viral_genes"))
                if vg is not None:
                    viral_genes_vals.append(vg)

                biome = mrow.get("biomes", "")
                if biome:
                    biomes.update(b.strip() for b in biome.split(",") if b.strip())

                seq_type = mrow.get("type", "")
                if seq_type:
                    types.add(seq_type.strip())

                cq = mrow.get("checkv_quality", "")
                if cq and cq not in ("", MISSING, "NA"):
                    members_with_checkv += 1

                comp = _safe_float(mrow.get("completeness"))
                if comp is not None:
                    completeness_vals.append(comp)

            mean_viral_genes = (
                f"{sum(viral_genes_vals) / len(viral_genes_vals):.2f}"
                if viral_genes_vals else MISSING
            )
            completeness_range = (
                f"{min(completeness_vals):.1f}-{max(completeness_vals):.1f}"
                if completeness_vals else MISSING
            )

            writer.writerow({
                "Sequence_ID":               rep_desc,
                "Sequence_ID_catalogue":     rep_mgyv or MISSING,
                "vitap_lineage":             vitap.get(rep_desc, MISSING),
                "viphogs_lineage":           viphogs_tax.get(rep_mgyv, MISSING) if rep_mgyv else MISSING,
                "genomad_lineage":           genomad.get(rep_mgyv, MISSING) if rep_mgyv else MISSING,
                "Source":                    rep_row.get("type", MISSING),
                "Biome":                     rep_row.get("biomes", MISSING),
                "Sequence_length":           rep_row.get("sequence_length", MISSING),
                "cluster_size":              cluster_size,
                "mean_viral_genes":          mean_viral_genes,
                "members_with_checkv":       members_with_checkv,
                "cluster_biomes":            ";".join(sorted(biomes)) or MISSING,
                "cluster_types":             ";".join(sorted(types)) or MISSING,
                "cluster_completeness_range": completeness_range,
                **checkv,
            })
    return out_path


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    print("Loading cluster files...")
    virus_member_to_rep, virus_rep_to_members = load_cluster(args.viruses_cluster)
    plasmid_member_to_rep, _ = load_cluster(args.plasmids_cluster)

    print("Loading genome metadata...")
    all_meta = load_all_metadata(args.additional_metadata)

    print("Loading ViTAP taxonomy...")
    vitap = load_vitap(args.viruses_vitap)

    print("Loading ViPhOGs taxonomy...")
    viphogs_tax = load_viphogs_taxonomy(args.viphogs_taxonomy)

    print("Loading geNomad taxonomy...")
    genomad = load_genomad(args.genomad)

    print("Loading rename map...")
    desc_to_mgyv, mgyv_to_desc = load_map(args.map)

    print("Loading combined metadata (with checkV fields)...")
    all_rows = load_combined_meta(args.combined_meta)

    print("Splitting combined metadata into viruses and plasmids...")
    virus_rows = {}
    plasmid_rows = {}
    for desc, row in all_rows.items():
        if "plasmid" in desc:
            plasmid_rows[desc] = row
        elif "viral_sequence" in desc or "prophage" in desc:
            virus_rows[desc] = row

    print(f"  Viruses/prophages: {len(virus_rows)}")
    print(f"  Plasmids: {len(plasmid_rows)}")

    print(f"Writing {args.output_viruses}...")
    written = write_viruses(
        virus_rows, virus_member_to_rep, all_meta, vitap,
        viphogs_tax, genomad, desc_to_mgyv,
        args.output_viruses, args.compress,
    )
    print(f"  -> {written}")

    print(f"Writing {args.output_plasmids}...")
    written = write_plasmids(
        plasmid_rows, plasmid_member_to_rep, all_meta,
        args.output_plasmids, args.compress,
    )
    print(f"  -> {written}")

    print(f"Writing {args.output_reps}...")
    written = write_virus_reps_stats(
        virus_rows, virus_rep_to_members, all_meta, vitap,
        viphogs_tax, genomad, desc_to_mgyv, mgyv_to_desc,
        args.output_reps, args.compress,
    )
    print(f"  -> {written}")

    print("Done.")


if __name__ == "__main__":
    main()
