#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
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
        help="Path to NFS location of catalogue(s) or catalogue_results/species_catalogue",
        nargs='+'
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Path to save results",
        default='.'
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Output filename",
        default='crispr_results'
    )
    parser.add_argument(
        "--update",
        action="store_true",
        help="Add spacers from --catalogue-path to a previous output of this script "
             "(requires --previous-fasta and --previous-metadata)"
    )
    parser.add_argument(
        "--previous-fasta",
        help="Previous <prefix>_crispr.fasta (plain or .gz); used with --update"
    )
    parser.add_argument(
        "--previous-metadata",
        help="Previous <prefix>_crispr.tsv (crispr_id, name, parent; plain or .gz); used with --update"
    )
    args = parser.parse_args()
    if args.update and not (args.previous_fasta and args.previous_metadata):
        parser.error("--update requires both --previous-fasta and --previous-metadata")
    if not args.update and (args.previous_fasta or args.previous_metadata):
        parser.error("--previous-fasta/--previous-metadata are only used together with --update")
    return args


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


def open_text(path: str):
    """Open a plain or gzip-compressed text file for reading."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path)


def find_gff(genome_dir: str, rep: str) -> str | None:
    """Return the path of rep's CRISPRCasFinder GFF (plain or .gz) in genome_dir, or None."""
    base = os.path.join(genome_dir, f'{rep}_crisprcasfinder.gff')
    for candidate in (base, base + '.gz'):
        if os.path.exists(candidate):
            return candidate
    return None


def parse_gff(gff: str) -> list[dict]:
    """Return list of CRISPRspacer records from a GFF file (plain or .gz).

    Each record is a dict with keys: crispr_id, name, seq, parent.
    """
    spacers = []
    with open_text(gff) as fh:
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


def seq_hash(seq: str) -> str:
    """Return SHA256 hex digest of a sequence for identity comparison."""
    return hashlib.sha256(seq.encode()).hexdigest()


def deduplicate_by_seq(spacers: list[dict]) -> list[dict]:
    """Return unique spacers by sequence content (SHA256), collecting all unique parents."""
    seen = {}  # sha256 -> dict with merged parents
    for s in spacers:
        key = seq_hash(s['seq'])
        if key not in seen:
            seen[key] = {**s, '_parents': [s['parent']]}
        else:
            if s['parent'] not in seen[key]['_parents']:
                seen[key]['_parents'].append(s['parent'])

    result = []
    for entry in seen.values():
        entry['parent'] = ','.join(entry.pop('_parents'))
        result.append(entry)
    return result


def read_fasta(path: str) -> dict[str, str]:
    """Return name -> sequence from a (multi-line, plain or .gz) FASTA file."""
    sequences = {}
    name, chunks = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    sequences[name] = ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        sequences[name] = ''.join(chunks)
    return sequences


def load_previous(fasta: str, metadata: str) -> list[dict]:
    """Return previous spacers (crispr_id, name, seq, parent) in their original order.

    Records are joined by `name`; metadata rows without a sequence in the
    FASTA are skipped with a warning.
    """
    check_path_exists(fasta)
    check_path_exists(metadata)
    sequences = read_fasta(fasta)

    previous = []
    missing_seq = 0
    with open_text(metadata) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            name = (row.get('name') or '').strip()
            if name not in sequences:
                missing_seq += 1
                continue
            previous.append({
                'crispr_id': row.get('crispr_id', 'missing'),
                'name':      name,
                'seq':       sequences[name],
                'parent':    row.get('parent', 'missing'),
            })
    if missing_seq:
        print(f'Warning: {missing_seq} previous metadata rows have no sequence in {fasta}; skipped')
    return previous


def merge_with_previous(previous: list[dict], new_spacers: list[dict]) -> list[dict]:
    """Add deduplicated new spacers to the previous set.

    - A new spacer whose sequence is already in the previous set is not
      added again; any parents it has that are not yet listed are appended
      to the existing record's comma-separated `parent`.
    - Otherwise it is appended as a new record. If its name is already used
      by a previous spacer with a different sequence, a `_new<N>` suffix is
      added so names stay unique.

    Previous records keep their order and come first; new records follow.
    """
    used_names = {p['name'] for p in previous}
    merged = [dict(p) for p in previous]
    by_hash = {seq_hash(p['seq']): p for p in merged}

    added = updated = renamed = 0
    for s in new_spacers:
        existing = by_hash.get(seq_hash(s['seq']))
        if existing is not None:
            parents = [p for p in existing['parent'].split(',') if p]
            new_parents = [p for p in s['parent'].split(',') if p and p not in parents]
            if new_parents:
                existing['parent'] = ','.join(parents + new_parents)
                updated += 1
            continue

        record = dict(s)
        if record['name'] in used_names:
            n = 1
            while f"{s['name']}_new{n}" in used_names:
                n += 1
            record['name'] = f"{s['name']}_new{n}"
            renamed += 1
        used_names.add(record['name'])
        by_hash[seq_hash(record['seq'])] = record
        merged.append(record)
        added += 1

    print(f'Update: {len(previous)} previous spacers; {added} new spacers added, '
          f'{len(new_spacers) - added} already present ({updated} of them got new parents)')
    if renamed:
        print(f'Warning: {renamed} new spacers shared a name with a different previous spacer '
              f'and were renamed with a _new<N> suffix')
    return merged


def list_rep_dirs(catalogue_path: str) -> list[tuple[str, str]]:
    """Return sorted (rep, rep_dir) pairs for every representative in a catalogue.

    Two layouts are supported:

    - species catalogue (``species_catalogue`` in the path): reps are grouped
      in prefix directories, e.g.
          catalogue_path/MGYG0005/MGYG000520529/genome/MGYG000520529_crisprcasfinder.gff[.gz]
      Every MGYG* directory directly under catalogue_path is treated as a
      group, and every MGYG* directory inside it as a rep.
    - otherwise: reps sit directly under catalogue_path, e.g.
          catalogue_path/MGYG000520529/genome/MGYG000520529_crisprcasfinder.gff[.gz]
    """
    def mgyg_dirs(path: str) -> list[str]:
        return sorted(
            item for item in os.listdir(path)
            if item.startswith('MGYG') and os.path.isdir(os.path.join(path, item))
        )

    if 'species_catalogue' in catalogue_path:
        reps = []
        for group in mgyg_dirs(catalogue_path):
            group_dir = os.path.join(catalogue_path, group)
            reps.extend((rep, os.path.join(group_dir, rep)) for rep in mgyg_dirs(group_dir))
        return sorted(reps)

    return [(rep, os.path.join(catalogue_path, rep)) for rep in mgyg_dirs(catalogue_path)]


def process_catalogue(catalogue_path: str) -> list[dict]:
    """Collect raw CRISPRspacer records from all reps in a catalogue.

    See list_rep_dirs() for the supported directory layouts. Inside each rep
    directory the GFF is expected at genome/<rep>_crisprcasfinder.gff, plain
    or gzip-compressed (.gff.gz).

    Returns:
        List of spacer dicts (crispr_id, name, seq, parent). Not deduplicated.
    """
    check_path_exists(catalogue_path)
    if 'species_catalogue' in catalogue_path:
        print('  Using species_catalogue layout: <group>/<rep>/genome/')

    spacers = []
    for rep, rep_dir in list_rep_dirs(catalogue_path):
        gff = find_gff(os.path.join(rep_dir, 'genome'), rep)
        if gff is None:
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

    if args.update:
        # read previous outputs fully before writing, so they can be overwritten in place
        previous = load_previous(args.previous_fasta, args.previous_metadata)
        unique_spacers = merge_with_previous(previous, unique_spacers)
        print(f'Total after update: {len(unique_spacers)} spacers')

    write_outputs(unique_spacers, args.output_path, args.prefix)


if __name__ == '__main__':
    main()
