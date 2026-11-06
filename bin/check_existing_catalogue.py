#!/usr/bin/env python3
"""
Check that the files of an existing catalogue are the ones its
catalogue_manifest.json describes (--mode update).

Every file listed under "files" in the manifest must be in the working
directory under the name it was published with (the last part of its path),
with the size and SHA-256 recorded in the manifest. All problems are listed
before the script exits, so one run shows everything that is wrong.

Writes a JSON report with the manifest's checksum (the new catalogue records it
as its parent), the accession range and the number of files checked.
"""

import argparse
import hashlib
import json
import os
import sys

SUPPORTED_MANIFEST_VERSIONS = {1}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True, help="catalogue_manifest.json of the existing catalogue")
    parser.add_argument("--rename-accession", required=True, help="--rename_accession of this run")
    parser.add_argument("-o", "--output", default="existing_catalogue_check.json")
    return parser.parse_args(argv)


def sha256sum(path, block_size=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_entries(files):
    """(key, entry) for every file in the manifest; list keys give one pair per file."""
    for key, value in files.items():
        for entry in value if isinstance(value, list) else [value]:
            yield key, entry


def check_manifest(manifest, rename_accession):
    """Problems with the manifest itself (format, update readiness, accession prefix)."""
    problems = []
    version = manifest.get("manifest_version")
    if version not in SUPPORTED_MANIFEST_VERSIONS:
        problems.append(f"manifest_version {version} is not supported (supported: {sorted(SUPPORTED_MANIFEST_VERSIONS)})")
    clustering = manifest.get("viral_clustering", {})
    if not clustering.get("update_ready"):
        problems.append(f"the catalogue cannot be updated: {clustering.get('not_update_ready_reason')}")
    prefix = manifest.get("accessions", {}).get("prefix")
    if prefix != rename_accession:
        problems.append(f"the catalogue uses the accession prefix '{prefix}', this run --rename_accession '{rename_accession}'")
    if not isinstance(manifest.get("files"), dict) or not manifest["files"]:
        problems.append("the manifest lists no files")
    return problems


def check_files(files):
    """Problems with the listed files: missing, repeated name, wrong size or checksum."""
    problems = []
    seen = {}
    checked = 0
    for key, entry in manifest_entries(files):
        name = os.path.basename(entry["path"])
        if name in seen:
            problems.append(f"{key}: {entry['path']} and {seen[name]} have the same file name")
            continue
        seen[name] = entry["path"]
        if not os.path.isfile(name):
            problems.append(f"{key}: {entry['path']} is missing")
            continue
        size = os.path.getsize(name)
        if size != entry["size"]:
            problems.append(f"{key}: {entry['path']} has {size} bytes, the manifest says {entry['size']}")
            continue
        if sha256sum(name) != entry["sha256"]:
            problems.append(f"{key}: {entry['path']} has a different SHA-256 than in the manifest")
            continue
        checked += 1
    return problems, checked


def main(argv=None):
    args = parse_args(argv)
    with open(args.manifest) as fh:
        manifest = json.load(fh)

    problems = check_manifest(manifest, args.rename_accession)
    checked = 0
    if isinstance(manifest.get("files"), dict):
        file_problems, checked = check_files(manifest["files"])
        problems += file_problems

    if problems:
        print(f"ERROR: the existing catalogue ({args.manifest}) does not match its manifest:", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
        sys.exit(1)

    report = {
        "status": "ok",
        "manifest_sha256": sha256sum(args.manifest),
        "manifest_version": manifest["manifest_version"],
        "catalogue": manifest.get("catalogue", {}),
        "accessions": manifest["accessions"],
        "files_checked": checked,
    }
    with open(args.output, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")
    print(f"Existing catalogue OK: {checked} files match the manifest; last accession {manifest['accessions'].get('last')}")


if __name__ == "__main__":
    main()
