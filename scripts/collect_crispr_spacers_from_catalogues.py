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
    parser.add_argument(
        "--catalogue-metadata",
        nargs='+',
        help="Genome metadata table(s) of the catalogue(s) (TSV with header, genome ID in the "
             "first column; plain or .gz). Concatenated into <prefix>_catalogue_metadata.tsv "
             "with one header and one row per genome"
    )
    parser.add_argument(
        "--previous-catalogue-metadata",
        help="Previous <prefix>_catalogue_metadata.tsv (plain or .gz); used with --update. "
             "Genomes also present in --catalogue-metadata are replaced by the new row"
    )
    args = parser.parse_args()
    if args.update and not (args.previous_fasta and args.previous_metadata):
        parser.error("--update requires both --previous-fasta and --previous-metadata")
    if not args.update and (args.previous_fasta or args.previous_metadata or args.previous_catalogue_metadata):
        parser.error("--previous-fasta/--previous-metadata/--previous-catalogue-metadata "
                     "are only used together with --update")
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


def read_metadata_table(path: str) -> tuple[list[str], list[list[str]]]:
    """Return (header, rows) of a TSV metadata table (plain or .gz); empty lines are skipped."""
    check_path_exists(path)
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader, None)
        if header is None:
            print(f'Warning: {path} is empty; skipped')
            return [], []
        rows = [row for row in reader if row and any(cell.strip() for cell in row)]
    return header, rows


def collect_catalogue_metadata(new_paths: list[str], previous_path: str | None) -> tuple[list[str], list[list[str]]]:
    """Concatenate genome metadata tables, keeping one header and one row per genome.

    Rows are unique by the first column (genome ID). Previous rows come
    first, in their original order; a genome that also appears in a new
    table is replaced in place by the new row (the newer metadata wins).
    Genomes only in the new tables are appended in input order. Within the
    new tables, the first occurrence of a genome is kept.

    All tables must have the same header.
    """
    header: list[str] = []
    rows_by_genome: dict[str, list[str]] = {}

    def check_header(path: str, table_header: list[str]) -> None:
        nonlocal header
        if not header:
            header = table_header
        elif table_header != header:
            print(f'Error: header of {path} differs from the first metadata table:\n'
                  f'  {header}\n  {table_header}')
            sys.exit(1)

    n_previous = 0
    if previous_path:
        table_header, rows = read_metadata_table(previous_path)
        if table_header:
            check_header(previous_path, table_header)
            for row in rows:
                rows_by_genome.setdefault(row[0], row)
            n_previous = len(rows_by_genome)

    replaced = added = repeated = 0
    seen_new: set[str] = set()
    for path in new_paths:
        table_header, rows = read_metadata_table(path)
        if not table_header:
            continue
        check_header(path, table_header)
        for row in rows:
            genome = row[0]
            if genome in seen_new:
                repeated += 1
                continue
            seen_new.add(genome)
            if genome in rows_by_genome:
                replaced += 1
            else:
                added += 1
            rows_by_genome[genome] = row  # dicts keep the original position on reassignment

    print(f'Catalogue metadata: {n_previous} previous genomes, {added} added, '
          f'{replaced} updated from new tables, {repeated} repeated rows in new tables skipped; '
          f'{len(rows_by_genome)} genomes in total')
    return header, list(rows_by_genome.values())


def write_catalogue_metadata(header: list[str], rows: list[list[str]], output_path: str, prefix: str) -> None:
    os.makedirs(output_path, exist_ok=True)
    path = os.path.join(output_path, prefix + '_catalogue_metadata.tsv')
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)
    print(f'Written: {path}')


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

    if args.catalogue_metadata or args.previous_catalogue_metadata:
        header, rows = collect_catalogue_metadata(args.catalogue_metadata or [], args.previous_catalogue_metadata)
        if header:
            write_catalogue_metadata(header, rows, args.output_path, args.prefix)


if __name__ == '__main__':
    main()
