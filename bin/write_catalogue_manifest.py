#!/usr/bin/env python3
"""
Write catalogue_manifest.json: what a later `--mode update` run needs to know
about this catalogue.

The manifest records
  - the accession range used by rename_contigs.py (an update continues after `last`)
  - how sequences were deduplicated (an update must hash the same way)
  - how viruses were clustered (an update needs the vjdb-clustering tables)
  - the software versions of the run
  - every published file an update reads: its path relative to --outdir,
    size and SHA-256, so an update can check it reads the right files

Files are given as `--entry KIND KEY TARGET NAME`:
  KIND    'file' (one file per key) or 'list' (several files under one key)
  KEY     manifest key, e.g. viruses.vjdb_clusters_csv
  TARGET  publish path relative to --outdir; a TARGET ending in '/' is a
          directory and the file keeps its name
  NAME    file name in the working directory
"""

import argparse
import gzip
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone

MANIFEST_VERSION = 1

# choose_sequences.py: seq_hash()
DEDUPLICATION = {
    "hash": "sha256",
    "normalisation": "upper().strip()",
    "column": "sequence_sha256",
}

# keys a vjdb update reads (see pipelines/vjdb-clustering, --mode update)
VJDB_UPDATE_KEYS = ("viruses.vjdb_clusters_csv", "viruses.vjdb_reps_fna")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--entry", nargs=4, action="append", default=[], metavar=("KIND", "KEY", "TARGET", "NAME"),
        help="Published file: kind (file|list), manifest key, publish path, local file name",
    )
    parser.add_argument("--rename-accession", required=True, help="Accession prefix of the catalogue, e.g. MGYV")
    parser.add_argument("--viral-clustering", required=True, help="--viral_clustering of the run")
    parser.add_argument("--mode", default="generate", choices=["generate", "update"], help="--mode of the run")
    parser.add_argument(
        "--parent-check", default=None,
        help="--mode update: existing_catalogue_check.json of the catalogue this run updates",
    )
    parser.add_argument("--pipeline-version", default=None)
    parser.add_argument("--run-name", default=None)
    parser.add_argument("--commit", default=None)
    parser.add_argument("-o", "--output", default="catalogue_manifest.json")
    return parser.parse_args(argv)


def open_text(path):
    with open(path, "rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gzipped else open(path)


def sha256sum(path, block_size=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def target_path(target, name):
    """Publish path of a file: a target ending in '/' is a directory."""
    path = f"{target}{name}" if target.endswith("/") else target
    return os.path.normpath(path)


def describe_file(target, name):
    if not os.path.isfile(name):
        sys.exit(f"Error: {name} not found in the working directory")
    return {"path": target_path(target, name), "size": os.path.getsize(name), "sha256": sha256sum(name)}


def collect_files(entries):
    """Turn the --entry arguments into {key: file} and {key: [file, ...]}."""
    files = {}
    for kind, key, target, name in entries:
        if kind == "file":
            if key in files:
                sys.exit(f"Error: manifest key {key} given twice")
            files[key] = describe_file(target, name)
        elif kind == "list":
            if key in files and not isinstance(files[key], list):
                sys.exit(f"Error: manifest key {key} is both a file and a list")
            files.setdefault(key, []).append(describe_file(target, name))
        else:
            sys.exit(f"Error: entry kind must be 'file' or 'list', got '{kind}'")
    for value in files.values():
        if isinstance(value, list):
            value.sort(key=lambda f: f["path"])
    paths = [f["path"] for value in files.values() for f in (value if isinstance(value, list) else [value])]
    clashes = sorted({p for p in paths if paths.count(p) > 1})
    if clashes:
        sys.exit(f"Error: several files are published to the same path: {clashes}")
    return dict(sorted(files.items()))


def accession_range(map_file, rename_accession):
    """First and last accession number in the rename map (column 'temporary').

    Every input sequence gets an accession, including the ones that are later
    deduplicated or filtered out, so `last` is the last number used.
    """
    pattern = re.compile(rf"^{re.escape(rename_accession)}(\d+)$")
    digits = []
    with open_text(map_file) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "temporary" not in header:
            sys.exit(f"Error: {map_file} has no 'temporary' column")
        col = header.index("temporary")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= col or fields == header:
                continue
            match = pattern.match(fields[col])
            if not match:
                sys.exit(f"Error: accession '{fields[col]}' in {map_file} does not match '{rename_accession}<number>'")
            digits.append(match.group(1))
    if not digits:
        return {"prefix": rename_accession, "first": None, "last": None, "count": 0, "zero_padded_digits": None}
    numbers = [int(d) for d in digits]
    if len(set(numbers)) != len(numbers):
        sys.exit(f"Error: {map_file} assigns the same accession number more than once")
    # zero-padded IDs (MGYV0000000001) have one fixed width; an update must use the same
    widths = {len(d) for d in digits}
    padded = any(d.startswith("0") for d in digits) and len(widths) == 1
    return {
        "prefix": rename_accession,
        "first": min(numbers),
        "last": max(numbers),
        "count": len(numbers),
        "zero_padded_digits": widths.pop() if padded else None,
    }


def read_versions(path):
    """Read the collated versions YAML ("PROCESS:" lines, then indented "tool: version")."""
    versions = {}
    process = None
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            key, _, value = line.strip().partition(":")
            key = key.strip().strip("'\"")
            value = value.strip().strip("'\"")
            if not line[0].isspace():
                process = key
                versions.setdefault(process, {})
            elif process is not None:
                versions[process][key] = value
    return dict(sorted(versions.items()))


def build_manifest(args, files):
    if "rename_map" not in files:
        sys.exit("Error: the rename map (key rename_map) is required")
    rename_map = next(name for kind, key, _t, name in args.entry if key == "rename_map")
    versions_name = next((name for kind, key, _t, name in args.entry if key == "software_versions"), None)

    missing_vjdb = [k for k in VJDB_UPDATE_KEYS if k not in files]
    if args.viral_clustering != "vjdb":
        reason = f"viruses were clustered with --viral_clustering {args.viral_clustering}; updates need 'vjdb'"
    elif missing_vjdb:
        reason = f"missing {', '.join(missing_vjdb)}"
    elif args.mode == "update":
        # until the update merges the annotations with those of the existing catalogue,
        # the annotation tables of an updated catalogue cover only its new clusters
        reason = "annotations of an updated catalogue are not merged with the existing ones yet"
    else:
        reason = None
    update_ready = reason is None

    parent = None
    if args.parent_check:
        with open(args.parent_check) as fh:
            check = json.load(fh)
        parent = {"manifest_sha256": check["manifest_sha256"], **check.get("catalogue", {})}
    elif args.mode == "update":
        sys.exit("Error: --mode update needs --parent-check")

    return {
        "manifest_version": MANIFEST_VERSION,
        "catalogue": {
            "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "mode": args.mode,
            "parent": parent,
            "pipeline_version": args.pipeline_version or None,
            "run_name": args.run_name or None,
            "commit": args.commit or None,
        },
        "accessions": accession_range(rename_map, args.rename_accession),
        "deduplication": DEDUPLICATION,
        "viral_clustering": {
            "method": args.viral_clustering,
            "update_ready": update_ready,
            "not_update_ready_reason": reason,
        },
        "software_versions": read_versions(versions_name) if versions_name else {},
        "files": files,
    }


def main(argv=None):
    args = parse_args(argv)
    files = collect_files(args.entry)
    manifest = build_manifest(args, files)
    with open(args.output, "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")
    print(f"Wrote {args.output}: {len(files)} file keys, update_ready={manifest['viral_clustering']['update_ready']}")


if __name__ == "__main__":
    main()
