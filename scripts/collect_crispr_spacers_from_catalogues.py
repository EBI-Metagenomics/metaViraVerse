#!/usr/bin/env python3

import argparse
import csv
import os
import sys


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Script searches for CrisprcasFinder results in catalogue(s) and greps CRISPRspacer information",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-p",
        "--catalogue-path",
        required=True,
        help="Path to NFS location of catalogue(s)",
        nargs='+'
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Path to save results",
        default='.'
    )
    return parser.parse_args()


def check_path_exists(path_str: str) -> None:
    """Verify that a file or directory exists, exit with error if not."""
    if not os.path.exists(path_str):
        print(f"Error: path not found: {path_str}")
        sys.exit(1)


def parse_gff_attributes(attr_str: str) -> dict:
    """Parse GFF9 attribute string into a key→value dict."""
    attrs = {}
    for part in attr_str.split(';'):
        if '=' in part:
            key, _, value = part.partition('=')
            attrs[key.strip()] = value.strip()
    return attrs


def parse_gff(gff: str) -> list[dict]:
    """Return list of CRISPRspacer records from a GFF file.

    Each record is a dict with keys: crispr_id, name, seq, parent.
    """
    spacers = []
    with open(gff) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'CRISPRspacer':
                continue
            attrs = parse_gff_attributes(parts[8])
            spacers.append({
                'crispr_id': attrs.get('ID', 'missing'),
                'name':      attrs.get('Name', 'missing'),
                'seq':       attrs.get('sequence', 'missing'),
                'parent':    attrs.get('Parent', 'missing'),
            })
    print(f'  Found {len(spacers)} CRISPRspacer records in {os.path.basename(gff)}')
    return spacers


def deduplicate_by_seq(spacers: list[dict]) -> list[dict]:
    """Return unique spacers keeping first occurrence per sequence."""
    seen = {}
    for s in spacers:
        seq = s['seq']
        if seq not in seen:
            seen[seq] = s
    return list(seen.values())


def process_catalogue(catalogue_path: str) -> list[dict]:
    """Collect raw CRISPRspacer records from all reps in a catalogue.

    Expected directory structure:
        catalogue_path/
            MGYG.../
                genome/
                    MGYG..._crisprcasfinder.gff

    Returns:
        List of spacer dicts (crispr_id, name, seq, parent). Not deduplicated.
    """
    check_path_exists(catalogue_path)
    reps = sorted(item for item in os.listdir(catalogue_path) if item.startswith('MGYG'))

    spacers = []
    for rep in reps:
        gff = os.path.join(catalogue_path, rep, 'genome', f'{rep}_crisprcasfinder.gff')
        if not os.path.exists(gff):
            continue
        print(f'  Processing: {rep}')
        spacers.extend(parse_gff(gff))
    return spacers


def write_outputs(unique_spacers: list[dict], output_path: str, prefix: str) -> None:
    """Write deduplicated spacers to TSV and FASTA."""
    os.makedirs(output_path, exist_ok=True)

    final_tsv = os.path.join(output_path, prefix + '_crispr.tsv')
    with open(final_tsv, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['crispr_id', 'name', 'parent'], delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(unique_spacers)
    print(f'Written: {final_tsv}')

    final_fasta = os.path.join(output_path, prefix + '_crispr.fasta')
    with open(final_fasta, 'w') as fh:
        for s in unique_spacers:
            fh.write(f'>{s["name"]}\n{s["seq"]}\n')
    print(f'Written: {final_fasta}')


def main() -> None:
    args = parse_arguments()

    all_spacers = []
    for catalogue_path in args.catalogue_path:
        path_parts = catalogue_path.rstrip('/').split('/')
        catalogue_name = '_'.join(path_parts[-2:])
        print(f'Running search for {catalogue_name}')
        all_spacers.extend(process_catalogue(catalogue_path))

    unique_spacers = deduplicate_by_seq(all_spacers)
    print(f'Total: {len(all_spacers)} spacers → {len(unique_spacers)} unique by sequence')

    prefix = '_'.join(
        '_'.join(p.rstrip('/').split('/')[-2:]) for p in args.catalogue_path
    ) if len(args.catalogue_path) > 1 else '_'.join(args.catalogue_path[0].rstrip('/').split('/')[-2:])

    write_outputs(unique_spacers, args.output_path, prefix)


if __name__ == '__main__':
    main()
