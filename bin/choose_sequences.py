#!/usr/bin/env python3
"""
Deduplicate virus/prophage/plasmid sequences across three category inputs.

Sequences are pre-split by category (virus/prophage/plasmid) upstream by
separate_sequences.py, based on the rename map's 'definition' column. This
script still hashes every sequence and deduplicates *across all three
categories together*: the same underlying sequence can legitimately end up
classified as e.g. a prophage in one sample and a free virus in another, and
such cross-category duplicates need to be caught and collapsed just like
same-category ones.

When the same sequence appears more than once, the record from the
highest-priority source is kept (MGnify over third-party; metagenome/assembly
over genome/MAG within either), its category becomes the deduplicated
sequence's final category, and biomes from every contributing source are
merged. Corresponding GFF and protein (FAA) records are filtered the same way
so that dropped duplicates don't leave orphaned annotations behind.

Input:  Up to three (fna, gff, faa) triples -- one per category.
Output: Per-category and combined deduplicated/filtered FNA, GFF and FAA
        files, plus a metadata TSV.

Usage:
    choose_sequences.py \
        --viruses viruses.fna --viruses-gff viruses.gff --viruses-faa viruses.faa \
        --prophages prophages.fna --prophages-gff prophages.gff --prophages-faa prophages.faa \
        --plasmids plasmids.fna --plasmids-gff plasmids.gff --plasmids-faa plasmids.faa \
        --map mapping.tsv \
        --rrna barrnap.gff \
        --quality quality_summary.tsv \
        --output-prefix combined
"""
from __future__ import annotations

import argparse
import csv
import hashlib

from utils import parse_attributes

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


CATEGORIES = ('virus', 'prophage', 'plasmid')
PLURALS = {'virus': 'viruses', 'prophage': 'prophages', 'plasmid': 'plasmids'}


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace with viruses/viruses_gff/viruses_faa (and the
        equivalent prophages_*/plasmids_* options), map (list[str]),
        rrna (list[str] | None), quality (list[str] | None), output_prefix.
    """
    parser = argparse.ArgumentParser(
        description="Deduplicate sequences across virus/prophage/plasmid categories, "
                    "prioritising assembly over MAG and MGnify over third-party."
    )
    for category in CATEGORIES:
        plural = PLURALS[category]
        parser.add_argument(f"--{plural}", required=False, help=f"{plural.capitalize()} FASTA file")
        parser.add_argument(f"--{plural}-gff", required=False, help=f"GFF file for {plural}")
        parser.add_argument(f"--{plural}-faa", required=False, help=f"Protein FASTA file for {plural}")
    parser.add_argument(
        "--map",
        required=False,
        nargs='+',
        help="Map file(s) with original contig names, biome, type and definition"
    )
    parser.add_argument(
        "--rrna",
        required=False,
        nargs='+',
        help="GFF(s) with rrna predictions (covers virus/prophage sequences)"
    )
    parser.add_argument(
        "--quality",
        required=False,
        nargs='+',
        help="CheckV results quality_summary.tsv (covers virus/prophage sequences)"
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix in output files"
    )
    args = parser.parse_args()

    if not any(getattr(args, PLURALS[category]) for category in CATEGORIES):
        parser.error("At least one of --viruses, --prophages, --plasmids must be provided")

    return args


def type_priority(source_type: str, definition: str) -> int:
    """Rank a sequence's source for deduplication priority (lower wins).

    Priority order: MGnify-metagenome (0) > MGnify-genome (1) >
    third-party-metagenome (2) > third-party-genome (3). Unrecognised source
    types get the lowest priority.

    Args:
        source_type: 'genome' or 'metagenome' -- MGnify's own 'type' column,
            or the third-party samplesheet's 'source' column (both written
            into the rename map's 'type' column by rename_contigs.py).
        definition: The rename map's 'definition' column value (e.g. 'virus'
            or 'third_party_virus'). The 'third_party_' prefix is what marks
            a sequence as coming from the third-party samplesheet rather than
            MGnify.

    Returns:
        Integer priority; lower is kept over higher when sequences collide.
    """
    if source_type == 'metagenome':
        rank = 0
    elif source_type == 'genome':
        rank = 1
    else:
        return 99
    is_third_party = bool(definition) and definition.startswith('third_party_')
    return rank + (2 if is_third_party else 0)


def seq_hash(sequence: str) -> str:
    """Compute SHA256 hash of a sequence string (uppercased, stripped).

    Args:
        sequence: Nucleotide sequence string.

    Returns:
        Hex digest of the SHA256 hash.
    """
    return hashlib.sha256(sequence.upper().strip().encode()).hexdigest()


def read_map(map_files: list[str]) -> dict[str, dict[str, str]]:
    """Read one or more TSV mapping files and return per-contig metadata.

    Expects a header row with at least ``original`` and ``temporary``
    columns. The ``biome``, ``type`` and ``definition`` columns are read when
    present (written by rename_contigs.py); absent columns fall back to
    ``'NA'``.

    Args:
        map_files: Paths to mapping TSV files.

    Returns:
        Dict mapping temporary contig name ->
        ``{'original': str, 'biome': str, 'type': str, 'definition': str}``.

    Raises:
        SystemExit: If a temporary name appears in more than one mapping entry.
    """
    mapping: dict[str, dict[str, str]] = {}
    for map_file in map_files:
        with open(map_file, 'r') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                temporary = row['temporary']
                if temporary in mapping:
                    print(f'Mapping already exists {temporary}. Exit')
                    exit(1)
                mapping[temporary] = {
                    'original':   row['original'],
                    'biome':      row.get('biome') or 'NA',
                    'type':       row.get('type') or 'NA',
                    'definition': row.get('definition') or 'NA',
                }
    return mapping


def parse_rna_gff(gff_files: list[str]) -> set[str]:
    """Parse one or more barrnap GFF files and return the set of sequence IDs with RNA features.

    Recognises rRNA, tRNA, and tmRNA feature types.

    Args:
        gff_files: Paths to GFF files produced by barrnap (or compatible tools).

    Returns:
        Set of sequence IDs (column 1) that carry at least one RNA annotation.
    """
    rna_sequences: set[str] = set()
    if not gff_files:
        return rna_sequences

    for gff_file in gff_files:
        with open(gff_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) >= 3:
                    seq_id = fields[0]
                    feature_type = fields[2]
                    if feature_type in ("rRNA", "tRNA", "tmRNA"):
                        rna_sequences.add(seq_id)
    return rna_sequences


QUALITY_COLUMNS = [
    'provirus', 'proviral_length', 'gene_count', 'viral_genes', 'host_genes',
    'checkv_quality', 'miuvig_quality', 'completeness', 'completeness_method',
    'contamination', 'kmer_freq', 'warnings',
]


def read_quality(quality_files: list[str]) -> dict[str, dict[str, str]]:
    """Read one or more CheckV quality_summary.tsv files.

    The key column ``contig_id`` is used for lookup; ``contig_length`` is
    intentionally excluded.  All remaining columns are retained as strings.

    Args:
        quality_files: Paths to CheckV quality_summary.tsv files.

    Returns:
        Dict mapping contig_id -> quality column dict (see QUALITY_COLUMNS).
    """
    quality: dict[str, dict[str, str]] = {}
    for qfile in quality_files:
        with open(qfile) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                contig_id = row['contig_id']
                if contig_id in quality:
                    print(f'Warning: duplicate quality entry for {contig_id}, keeping first')
                    continue
                quality[contig_id] = {col: row.get(col, 'NA') for col in QUALITY_COLUMNS}
    return quality


def filter_reason(entry: dict) -> str | None:
    """Return the exclusion reason for an entry, or None if it passes.

    Args:
        entry: A ``seen`` dict entry as built in ``choose_seqs()``.

    Returns:
        A short reason string if excluded, None if the sequence should be kept.
    """
    if entry['rrna'] == 'Yes' and entry['category'] != 'plasmid':
        return 'rrna'
    q = entry['quality']
    if q is not None and q.get('checkv_quality') == 'Not-determined' and entry['category'] != 'plasmid':
        try:
            viral_genes = float(q.get('viral_genes') or 0)
            kmer_freq = float(q.get('kmer_freq') or 0)
        except (ValueError, TypeError):
            return None
        if viral_genes == 0 and kmer_freq >= 1.0:
            return 'not_determined'
    return None


def passes_filter(entry: dict) -> bool:
    return filter_reason(entry) is None


def choose_seqs(
    category_fastas: list[tuple[str, str | None]],
    mapping: dict[str, dict[str, str]],
    rna_sequences: set[str],
    quality_data: dict[str, dict[str, str]],
) -> dict[str, dict]:
    """Deduplicate sequences across virus/prophage/plasmid inputs, keeping the highest-priority source.

    For each unique sequence (by SHA256), the record from the source with the
    lowest ``type_priority`` value is retained -- this also decides the final
    category for that sequence when it was classified differently in more
    than one input (e.g. a prophage call in one sample vs. a free virus call
    in another). Biome labels from all sources are merged regardless of which
    record wins.

    Args:
        category_fastas: List of (category, fna_path) pairs, e.g.
            [('virus', 'viruses.fna'), ('prophage', 'prophages.fna'), ('plasmid', 'plasmids.fna')].
            A ``None``/empty fna_path is skipped.
        mapping: Dict of temporary -> {'original', 'biome', 'type', 'definition'} (from read_map).
        rna_sequences: Set of sequence IDs that carry rRNA/tRNA/tmRNA annotations.
        quality_data: Dict of contig_id -> CheckV quality column dict.

    Returns:
        Dict keyed by sequence SHA256 hash; each value is a record entry dict.
    """
    seen: dict[str, dict] = {}

    for category, fna_file in category_fastas:
        if not fna_file:
            continue
        for record in SeqIO.parse(fna_file, "fasta"):
            h = seq_hash(str(record.seq))
            map_entry = mapping.get(record.id, {})
            original_name = map_entry.get('original', record.id)
            source_type   = map_entry.get('type',       'NA')
            definition    = map_entry.get('definition',  'NA')
            biome         = map_entry.get('biome',      'NA')
            rrna = ('Yes' if record.id in rna_sequences else 'No') if rna_sequences else 'not-provided'

            if h not in seen:
                seen[h] = {
                    'record': record,
                    'original_name': original_name,
                    'type': source_type,
                    'definition': definition,
                    'category': category,
                    'biomes': {biome},
                    'seq_id': record.id,
                    'description': record.description,
                    'rrna': rrna,
                    'quality': quality_data.get(record.id),
                }
            else:
                entry = seen[h]
                # Always merge biomes, regardless of which source wins
                entry['biomes'].add(biome)

                current_priority = type_priority(entry['type'], entry['definition'])
                new_priority = type_priority(source_type, definition)
                if new_priority < current_priority:
                    entry['record'] = record
                    entry['original_name'] = original_name
                    entry['type'] = source_type
                    entry['definition'] = definition
                    entry['category'] = category
                    entry['seq_id'] = record.id
                    entry['description'] = record.description
                    entry['rrna'] = rrna
                    entry['quality'] = quality_data.get(record.id)
    return seen


def read_input_gff(gffs: list[str]) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Parse one or more GFF3 files and index records by sequence ID.

    Builds two data structures:
    - ``gff_data``: maps each sequence ID to the list of raw GFF lines that
      belong to it (including its CDS lines, tracked via the most recently
      seen non-CDS feature).
    - ``source_map``: maps each sequence ID to the GFF source field
      (column 2) of its non-CDS feature (e.g. ``VIRify``, ``geNomad``).

    Args:
        gffs: Paths to input GFF3 files.

    Returns:
        Tuple of (gff_data, source_map).
    """
    gff_data: dict[str, list[str]] = {}
    source_map: dict[str, str] = {}
    total_lines = 0
    for gff in gffs:
        with open(gff, 'r') as file_in:
            for line in file_in:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                attrs, _ = parse_attributes(parts[8])
                if parts[2] != 'CDS':
                    if attrs.get('ID'):
                        attr_id = attrs['ID']
                    else:
                        print(f'There is no ID found for {line}')
                    gff_data.setdefault(attr_id, [])

                    if attr_id not in source_map:
                        source_map[attr_id] = parts[1]
                gff_data[attr_id].append(line)
                total_lines += 1
    print(f'Total lines in input GFF: {total_lines}')
    return gff_data, source_map


def read_faa(faa_files: list[str]) -> dict[str, SeqRecord]:
    """Read one or more protein FASTA files into a single ID -> record map.

    Args:
        faa_files: Paths to protein FASTA files (original, unrenamed protein IDs).

    Returns:
        Dict mapping protein ID -> Bio.SeqRecord.SeqRecord.
    """
    proteins: dict[str, SeqRecord] = {}
    for faa_file in faa_files:
        for record in SeqIO.parse(faa_file, "fasta"):
            proteins[record.id] = record
    return proteins


def write_final_files(
    metadata: str,
    filtered_fna: str,
    filtered_gff: str,
    filtered_faa: str,
    filtered_tsv: str,
    excluded_tsv: str,
    seen: dict[str, dict],
    gff_data: dict[str, list[str]],
    source_map: dict[str, str],
    faa_records: dict[str, SeqRecord],
    output_prefix: str,
) -> None:
    """Write all output files from the deduplicated record set.

    Writes every unique sequence to ``metadata``, sequences that pass
    ``passes_filter`` to the combined ``filtered_fna``/``filtered_gff``/
    ``filtered_faa`` (and to their category-specific equivalents), and
    excluded sequences with their reason to ``excluded_tsv``. Only proteins
    belonging to kept sequences are written to the FAA outputs, so dropped
    duplicates don't leave orphaned CDS records behind.

    Args:
        metadata: Path for the full metadata TSV (every unique sequence).
        filtered_fna: Path for the combined filtered FASTA (all categories).
        filtered_gff: Path for the combined filtered GFF.
        filtered_faa: Path for the combined filtered protein FASTA.
        filtered_tsv: Path for the combined filtered metadata TSV.
        excluded_tsv: Path for the excluded sequences TSV (with filter_reason column).
        seen: Dict returned by ``choose_seqs``.
        gff_data: Dict of seq_id -> list of raw GFF lines (from ``read_input_gff``).
        source_map: Dict of seq_id -> GFF source field (from ``read_input_gff``).
        faa_records: Dict of protein_id -> SeqRecord (from ``read_faa``).
        output_prefix: Prefix used to derive the per-category output filenames.
    """
    records_total = 0
    records_filtered = 0
    records_excluded = 0
    gff_records_written = 0
    faa_records_written = 0
    gff_found = 0

    category_fna_handles = {c: open(f"{output_prefix}_{c}_filtered.fna", 'w') for c in CATEGORIES}
    category_gff_handles = {c: open(f"{output_prefix}_{c}_filtered.gff", 'w') for c in CATEGORIES}
    category_faa_handles = {c: open(f"{output_prefix}_{c}_filtered.faa", 'w') for c in CATEGORIES}
    for fh in category_gff_handles.values():
        fh.write("##gff-version 3\n")

    with open(filtered_fna, 'w') as out_fna_f, \
         open(filtered_tsv, 'w') as out_tsv_f, \
         open(excluded_tsv, 'w') as out_excl, \
         open(filtered_gff, 'w') as filt_gff, \
         open(filtered_faa, 'w') as filt_faa, \
         open(metadata, 'w') as metadata_file:
        filt_gff.write("##gff-version 3\n")

        quality_header = '\t'.join(QUALITY_COLUMNS)
        header = (
            f"sequence_id\toriginal_name\tdescription\ttype\tsource_of_prediction\tbiomes\t"
            f"sequence_length\trrna\tsequence_sha256\t{quality_header}\tdefinition\n"
        )
        out_tsv_f.write(header)
        metadata_file.write(header)
        out_excl.write(f"filter_reason\t{header}")

        for h, entry in seen.items():
            record = entry['record']
            new_name = record.description.replace('|', '-').replace(' ', '|')
            record.id = new_name
            record.description = new_name

            biomes_str = ','.join(sorted(entry['biomes']))
            q = entry['quality']
            quality_values = '\t'.join(q[col] for col in QUALITY_COLUMNS) if q else '\t'.join(
                'NA' for _ in QUALITY_COLUMNS)
            prediction_source = source_map.get(entry['seq_id'], "NA")

            tsv_row = (
                f"{entry['seq_id']}\t"
                f"{entry['original_name']}\t"
                f"{entry['original_name'].replace('|', '-').replace(' ', '|')}\t"
                f"{entry['type']}\t"
                f"{prediction_source}\t"
                f"{biomes_str}\t"
                f"{len(record.seq)}\t"
                f"{entry['rrna']}\t"
                f"{h}\t"
                f"{quality_values}\t"
                f"{entry['definition']}\n"
            )
            records_total += 1
            metadata_file.write(tsv_row)

            reason = filter_reason(entry)
            if reason is None:
                category = entry['category']

                # fasta: combined + category-specific
                SeqIO.write([record], out_fna_f, "fasta")
                SeqIO.write([record], category_fna_handles[category], "fasta")

                # metadata tsv
                out_tsv_f.write(tsv_row)
                records_filtered += 1

                # gff: combined + category-specific, and the CDS ids it carries
                gff_records = gff_data.get(entry['seq_id'])
                if gff_records:
                    gff_found += 1
                    filt_gff.writelines(gff_records)
                    category_gff_handles[category].writelines(gff_records)
                    gff_records_written += len(gff_records)

                    for line in gff_records:
                        parts = line.rstrip('\n').split('\t')
                        if len(parts) >= 9 and parts[2] == 'CDS':
                            attrs, _ = parse_attributes(parts[8])
                            protein_id = attrs.get('ID')
                            protein_record = faa_records.get(protein_id) if protein_id else None
                            if protein_record is not None:
                                SeqIO.write([protein_record], filt_faa, "fasta")
                                SeqIO.write([protein_record], category_faa_handles[category], "fasta")
                                faa_records_written += 1
                else:
                    print(f"{entry['seq_id']} has no GFF records")
            else:
                out_excl.write(f"{reason}\t{tsv_row}")
                records_excluded += 1

    for fh in (*category_fna_handles.values(), *category_gff_handles.values(), *category_faa_handles.values()):
        fh.close()

    print(f"Total unique sequences processed: {records_total}")
    print(f"Filtered sequences written to {filtered_fna}: {records_filtered}")
    print(f"Excluded sequences written to {excluded_tsv}: {records_excluded}")
    print(f"Total written lines to filtered GFF {gff_records_written}")
    print(f"Proteins written to filtered FAA: {faa_records_written}")
    if records_filtered != gff_found:
        print(f"GFF entries found for {gff_found} sequences but {records_filtered} written to FASTA. Exit")
        exit(1)
    else:
        print('Sanity check passed')


def main() -> None:
    """Deduplicate sequences across the virus/prophage/plasmid inputs and write merged outputs."""
    args = parse_arguments()

    mapping: dict[str, dict[str, str]] = read_map(args.map) if args.map else {}
    rna_sequences: set[str] = parse_rna_gff(args.rrna) if args.rrna else set()
    quality_data: dict[str, dict[str, str]] = read_quality(args.quality) if args.quality else {}

    category_fastas = [(category, getattr(args, PLURALS[category])) for category in CATEGORIES]

    gff_files = [
        getattr(args, f"{PLURALS[category]}_gff") for category in CATEGORIES
        if getattr(args, f"{PLURALS[category]}_gff")
    ]
    input_gff, source_map = read_input_gff(gff_files)

    faa_files = [
        getattr(args, f"{PLURALS[category]}_faa") for category in CATEGORIES
        if getattr(args, f"{PLURALS[category]}_faa")
    ]
    faa_records = read_faa(faa_files) if faa_files else {}

    seen = choose_seqs(category_fastas, mapping, rna_sequences, quality_data)

    metadata     = f"{args.output_prefix}_metadata.tsv"
    filtered_fna = f"{args.output_prefix}_filtered.fna"
    filtered_gff = f"{args.output_prefix}_filtered.gff"
    filtered_faa = f"{args.output_prefix}_filtered.faa"
    filtered_tsv = f"{args.output_prefix}_filtered.tsv"
    excluded_tsv = f"{args.output_prefix}_excluded.tsv"

    write_final_files(
        metadata, filtered_fna, filtered_gff, filtered_faa, filtered_tsv, excluded_tsv,
        seen, input_gff, source_map, faa_records, args.output_prefix,
    )
    processed = [category for category, fna in category_fastas if fna]
    print(f"Categories processed: {processed}")


if __name__ == '__main__':
    main()
