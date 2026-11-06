#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

MOBILOME_GFF = "mobilome_full.gff.gz"
CDS_FAA = "predicted_cds.faa.gz"
CONTIGS_FASTA = "filtered_contigs.fasta.gz"


def find_download(paths: list[str], suffix: str) -> str | None:
    """Path (relative to results_dir) of the first download ending with suffix"""
    return next((p for p in paths if p.endswith(suffix)), None)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Builds a TSV (id, map_gff, cds_faa, contigs_fasta) from an analyses CSV export, "
        "for the analyses whose downloads list a mobilome_full.gff.gz"
    )
    parser.add_argument("-i", "--input", required=True, help="Analyses CSV (with a results_dir column)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    parser.add_argument(
        "--check-exists",
        action="store_true",
        help="Skip analyses whose mobilome GFF is not on disk and warn about missing CDS/contigs files",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    csv.field_size_limit(sys.maxsize)  # the downloads/metadata JSON columns are large

    kept = no_mobilome = no_results_dir = missing_file = ext_res_dir_exists = 0
    seen = set()
    with open(args.input, newline="") as fin, \
            open(args.output, "w", newline="") as fout, \
            open('no_mobilome.tsv', 'w', newline="") as fmob, \
            open('analyses_metadata.tsv', 'w', newline="") as fmetadata:
        fmob.write('id\tresults_dir\texternal_results_dir\n')
        fmetadata.write('\t'.join(['accession','study_id', 'sample_id', 'biome_id']) + '\n')
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "map_gff", "cds_faa", "contigs_fasta"])
        for row in csv.DictReader(fin):
            paths = [d["path"] for d in json.loads(row["downloads"] or "[]")]
            files = [find_download(paths, suffix) for suffix in (MOBILOME_GFF, CDS_FAA, CONTIGS_FASTA)]
            if files[0] is None:
                no_mobilome += 1
                fmob.write(row["id"] + '\t' + row["results_dir"] + '\t' + row["external_results_dir"] + '\n')
                continue
            results_dir = row["results_dir"].strip().rstrip("/")
            ext_results_dir = row["external_results_dir"].strip().rstrip("/")
            if not results_dir:
                no_results_dir += 1
                if ext_results_dir:
                    ext_res_dir_exists += 1
                continue
            # results_dir ends with .../<ERZ accession>/V6/assembly
            acc = results_dir.split("/")[-3]
            if acc in seen:
                continue
            if None in files:
                print(f"WARNING: {acc}: downloads have no cds FAA or contigs FASTA, skipped", file=sys.stderr)
                missing_file += 1
                continue
            map_gff, cds_faa, contigs_fasta = (os.path.join(results_dir, f) for f in files)
            if args.check_exists:
                if not os.path.isfile(map_gff):
                    print(f"WARNING: {acc}: {map_gff} does not exist, skipped", file=sys.stderr)
                    missing_file += 1
                    continue
                for path in (cds_faa, contigs_fasta):
                    if not os.path.isfile(path):
                        print(f"WARNING: {acc}: {path} does not exist", file=sys.stderr)
            writer.writerow([acc, map_gff.replace('full', 'clean'), cds_faa, contigs_fasta])
            study_accessions = row.get('ena_accessions', 'NA').replace('}', '').replace('{', '').split(',')
            primary_study = next((x for x in study_accessions if x.startswith("PRJ")), study_accessions[0])
            fmetadata.write('\t'.join([acc, primary_study, row.get('sample_id', 'NA'), row.get('biome_id', 'NA')]) + '\n')
            seen.add(acc)
            kept += 1

    print(
        f"{kept} assemblies written; skipped {no_mobilome} rows without a mobilome GFF in downloads, "
        f"{no_results_dir} without results_dir (external_result_dir exists for {ext_res_dir_exists}) and {missing_file} with missing files",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
