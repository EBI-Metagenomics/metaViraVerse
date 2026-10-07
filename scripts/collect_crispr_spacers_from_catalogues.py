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


def contig_from_parent(parent: str) -> str:
    """Return the contig of a CRISPR array Parent, e.g.
    'MGYG000296150_19_103_4169' -> 'MGYG000296150_19' (array start/end dropped).
    For a comma-separated list, the first parent is used.
    """
    first = parent.split(',')[0]
    parts = first.rsplit('_', 2)
    if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
        return parts[0]
    return first


def with_contig(contig: str, value: str) -> str:
    """Prefix `value` with `contig_` unless it already starts with it.

    CRISPRCasFinder names spacers by position and length only (e.g.
    'spacer_1131_35', ID 'sp_1131'), so the same name appears in many
    genomes. Prefixing the contig makes names unique across catalogues.
    """
    if not contig or value == 'missing' or value.startswith(f'{contig}_'):
        return value
    return f'{contig}_{value}'


def parse_gff(gff: str) -> list[dict]:
    """Return list of CRISPRspacer records from a GFF file (plain or .gz).

    Each record is a dict with keys: crispr_id, name, seq, parent. `name` and
    `crispr_id` are prefixed with the contig (GFF column 1) so they are unique
    across genomes; `seq` is upper-cased.
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
            contig = parts[0]
            spacers.append({
                'crispr_id': with_contig(contig, attrs.get('ID', 'missing')),
                'name':      with_contig(contig, attrs.get('Name', 'missing')),
                'seq':       attrs.get('sequence', 'missing').upper(),
                'parent':    attrs.get('Parent', 'missing'),
            })
    print(f'  Found {len(spacers)} CRISPRspacer records in {os.path.basename(gff)}')
    return spacers


def seq_hash(seq: str) -> str:
    """Return SHA256 hex digest of a sequence (case-insensitive) for identity comparison."""
    return hashlib.sha256(seq.upper().encode()).hexdigest()


def split_parents(parent: str) -> list[str]:
    return [p for p in parent.split(',') if p]


def deduplicate_by_seq(spacers: list[dict]) -> list[dict]:
    """Return unique spacers by sequence content (SHA256), collecting all unique parents.

    The first record seen for a sequence keeps its name and ID; `parent`
    values (which may already be comma-separated lists) are merged in order.
    """
    seen = {}  # sha256 -> dict with merged parents
    for s in spacers:
        key = seq_hash(s['seq'])
        if key not in seen:
            seen[key] = {**s, '_parents': []}
        for p in split_parents(s['parent']):
            if p not in seen[key]['_parents']:
                seen[key]['_parents'].append(p)

    result = []
    for entry in seen.values():
        entry['parent'] = ','.join(entry.pop('_parents'))
        result.append(entry)
    return result


def read_fasta_records(path: str) -> list[tuple[str, str]]:
    """Return (name, sequence) pairs in file order from a (multi-line, plain or .gz) FASTA.

    Duplicate names are kept as separate records (a dict would silently keep
    only the last sequence).
    """
    records = []
    name, chunks = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    records.append((name, ''.join(chunks)))
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, ''.join(chunks)))
    return records


def load_previous(fasta: str, metadata: str) -> list[dict]:
    """Return previous spacers (crispr_id, name, seq, parent), deduplicated by sequence.

    This script writes the TSV and FASTA in the same order, so rows are
    joined by position (and checked by name). This stays correct even when
    an older output has repeated names. If the files are not in the same
    order, rows are joined by name instead, which requires unique names.

    Names and IDs are brought to the current naming scheme (contig prefix,
    taken from the parent), sequences are upper-cased, and previous records
    sharing a sequence are merged.
    """
    check_path_exists(fasta)
    check_path_exists(metadata)
    fasta_records = read_fasta_records(fasta)
    with open_text(metadata) as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))

    same_order = len(rows) == len(fasta_records) and all(
        (row.get('name') or '').strip() == name for row, (name, _) in zip(rows, fasta_records)
    )
    if same_order:
        pairs = [(row, seq) for row, (_, seq) in zip(rows, fasta_records)]
    else:
        names = [name for name, _ in fasta_records]
        if len(names) != len(set(names)):
            print(f'Error: {fasta} and {metadata} are not in the same order and the FASTA '
                  f'has repeated names, so sequences cannot be matched to metadata reliably. '
                  f'Regenerate the previous output without --update first.')
            sys.exit(1)
        sequences = dict(fasta_records)
        pairs = []
        for row in rows:
            name = (row.get('name') or '').strip()
            if name in sequences:
                pairs.append((row, sequences[name]))
        if len(pairs) < len(rows):
            print(f'Warning: {len(rows) - len(pairs)} previous metadata rows have no sequence in {fasta}; skipped')

    previous = []
    for row, seq in pairs:
        parent = row.get('parent', 'missing')
        contig = contig_from_parent(parent)
        previous.append({
            'crispr_id': with_contig(contig, row.get('crispr_id', 'missing')),
            'name':      with_contig(contig, (row.get('name') or '').strip()),
            'seq':       seq.upper(),
            'parent':    parent,
        })

    deduplicated = deduplicate_by_seq(previous)
    if len(deduplicated) < len(previous):
        print(f'Previous: {len(previous)} records merged into {len(deduplicated)} unique by sequence')
    return deduplicated


def merge_with_previous(previous: list[dict], new_spacers: list[dict]) -> list[dict]:
    """Add deduplicated new spacers to the previous set, matching by sequence.

    - A new spacer whose sequence is already in the previous set is not
      added again; any parents it has that are not yet listed are appended
      to the existing record's comma-separated `parent`.
    - Otherwise it is appended as a new record.

    Previous records keep their order and come first; new records follow.
    """
    merged = [dict(p) for p in previous]
    by_hash = {seq_hash(p['seq']): p for p in merged}

    added = updated = 0
    for s in new_spacers:
        existing = by_hash.get(seq_hash(s['seq']))
        if existing is not None:
            parents = split_parents(existing['parent'])
            new_parents = [p for p in split_parents(s['parent']) if p not in parents]
            if new_parents:
                existing['parent'] = ','.join(parents + new_parents)
                updated += 1
            continue

        record = dict(s)
        by_hash[seq_hash(record['seq'])] = record
        merged.append(record)
        added += 1

    print(f'Update: {len(previous)} previous spacers; {added} new spacers added, '
          f'{len(new_spacers) - added} already present ({updated} of them got new parents)')
    return merged


def ensure_unique_names(spacers: list[dict]) -> list[dict]:
    """Safety net: suffix `_<N>` to any name that is still repeated.

    With contig-prefixed names this should only happen if one contig has two
    different spacers at the same position and length (e.g. annotations from
    two CRISPRCasFinder runs that disagree), which is reported.
    """
    counts: dict[str, int] = {}
    renamed = 0
    for s in spacers:
        n = counts.get(s['name'], 0)
        counts[s['name']] = n + 1
        if n:
            s['name'] = f"{s['name']}_{n + 1}"
            renamed += 1
    if renamed:
        print(f'Warning: {renamed} spacers had a repeated name and got a _<N> suffix')
    return spacers


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

    write_outputs(ensure_unique_names(unique_spacers), args.output_path, args.prefix)


if __name__ == '__main__':
    main()
