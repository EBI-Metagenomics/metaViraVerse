#!/usr/bin/env python3
"""
Deduplicate viral/plasmid/prophage sequences across multiple FNA files.

When the same sequence (by SHA256 hash) appears in both assembly and MAG sources,
the assembly record is kept but biomes from all sources are merged into the metadata.

Input:  Multiple FNA files with corresponding type and biome values.
Output: A deduplicated FNA file and a TSV metadata table.

Usage:
    choose_sequences.py \
        --fna sample1.fna sample2.fna \
        --type metagenome / genome \
        --biome marine soil \
        --output-fna combined.fna \
        --output-tsv metadata.tsv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

from Bio import SeqIO


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace with fna (list[str]), type (list[str]), biome (list[str]),
        output_fna (str), and output_tsv (str).
    """
    parser = argparse.ArgumentParser(
        description="Deduplicate sequences across FNA files, prioritising assembly over MAG."
    )
    parser.add_argument(
        "--fna",
        required=True,
        nargs='+',
        help="Input FNA file(s)"
    )
    parser.add_argument(
        "--gff",
        required=True,
        nargs='+',
        help="Input GFF file(s)"
    )
    parser.add_argument(
        "--map",
        required=False,
        nargs='+',
        help="Map file with original contig names and new names"
    )
    parser.add_argument(
        "--rrna",
        required=False,
        nargs='+',
        help="GFF with rrna predictions"
    )
    parser.add_argument(
        "--quality",
        required=False,
        nargs='+',
        help="CheckV results quality_summary.tsv"
    )
    parser.add_argument(
        "--output-fna",
        required=True,
        help="Path for output deduplicated FNA file"
    )
    parser.add_argument(
        "--output-gff",
        required=True,
        help="Path for output deduplicated filtered GFF file"
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        help="Path for output metadata TSV file"
    )
    args = parser.parse_args()

    if len(args.gff) != len(args.fna):
        parser.error("--fna and --gff must have the same number of values")

    return args


# Priority: lower number = higher priority
TYPE_PRIORITY = {
    'metagenome': 0,
    'genome': 1,
    'third_party_virus': 2,
    'third_party_plasmid': 3,
}


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

    Expects a header row with at least ``original`` and ``temporary`` columns.
    The ``biome`` and ``type`` columns are read when present (written by
    rename_contigs.py); absent columns fall back to ``'NA'``.

    Args:
        map_files: Paths to mapping TSV files.

    Returns:
        Dict mapping temporary contig name ->
        ``{'original': str, 'biome': str, 'type': str}``.

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
                    'original': row['original'],
                    'biome':    row.get('biome', 'NA'),
                    'type':     row.get('type',  'NA'),
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


def passes_filter(entry: dict) -> bool:
    """Return True if the sequence passes quality filters for the filtered output set.

    A sequence is excluded when either condition holds:
    - It carries an rRNA/tRNA/tmRNA annotation (rrna == 'Yes') and is not a plasmid.
    - CheckV determined the quality as 'Not-determined' (only when quality data exists).

    Args:
        entry: A ``seen`` dict entry as built in ``main()``.

    Returns:
        True if the sequence should be kept in the filtered output.
    """
    if entry['rrna'] == 'Yes' and entry['viral_type'] != 'plasmid':
        return False
    q = entry['quality']
    if q is not None and q.get('checkv_quality') == 'Not-determined' and entry['viral_type'] != 'plasmid':
        return False
    return True


def parse_attrs(attrs_str: str) -> tuple[dict[str, str], list[str]]:
    """Parse a GFF3 column-9 attributes string into a dict and an ordered key list.

    :param attrs_str: Semicolon-separated key=value attribute string from GFF column 9.
    :return: Tuple of (attrs dict, list of keys in original order).
    """
    attrs, order = {}, []
    for part in attrs_str.rstrip(";").split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            if k not in attrs:
                order.append(k)
            attrs[k] = v
    return attrs, order


def read_input_gff(gffs: list[str]) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Parse one or more GFF3 files and index records by sequence ID.

    Builds two data structures:
    - ``gff_data``: maps each sequence ID (GFF column 1) to the list of raw
      GFF lines that belong to it.
    - ``attr_id_to_seq_id``: maps the ``ID`` attribute of non-CDS features to
      their parent sequence ID, for downstream attribute renaming.

    Args:
        gffs: Paths to input GFF3 files.

    Returns:
        Tuple of (gff_data, attr_id_to_seq_id).
    """
    gff_data: dict[str, list[str]] = {}
    attr_id_to_seq_id: dict[str, str] = {}
    for gff in gffs:
        with open(gff, 'r') as file_in:
            for line in file_in:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                seq_id = parts[0]
                attrs, _ = parse_attrs(parts[8])
                gff_data.setdefault(seq_id, [])
                gff_data[seq_id].append(line)
                if parts[2] != 'CDS':
                    if attrs.get('ID'):
                        attr_id = attrs['ID']
                        if attr_id in attr_id_to_seq_id:
                            print(f"{attr_id} already exists in mapping {attr_id_to_seq_id[attr_id]}")
                        attr_id_to_seq_id.setdefault(attr_id, seq_id)
                    else:
                        print(f'There is no ID found for {line}')
    return gff_data, attr_id_to_seq_id


def choose_seqs(
    fnas: list[str],
    mapping: dict[str, dict[str, str]],
    rna_sequences: set[str],
    quality_data: dict[str, dict[str, str]],
) -> dict[str, dict]:
    """Deduplicate sequences across all FNA files, keeping the highest-priority source.

    For each unique sequence (by SHA256), the record from the source with the
    lowest TYPE_PRIORITY value is retained.  Biome labels from all sources are
    merged.  Quality and rRNA annotations are carried from the winning record.
    Biome and type are read from the mapping (populated by rename_contigs.py).

    Args:
        fnas: Paths to input FNA files.
        mapping: Dict of temporary -> {'original', 'biome', 'type'} (from read_map).
        rna_sequences: Set of sequence IDs that carry rRNA/tRNA/tmRNA annotations.
        quality_data: Dict of contig_id -> CheckV quality column dict.

    Returns:
        Dict keyed by sequence SHA256 hash; each value is a record entry dict.
    """
    seen: dict[str, dict] = {}

    for fna_file in fnas:
        for record in SeqIO.parse(fna_file, "fasta"):
            h = seq_hash(str(record.seq))
            map_entry = mapping.get(record.id, {})
            original_name = map_entry.get('original', record.id)
            source_type   = map_entry.get('type',     'NA')
            biome         = map_entry.get('biome',    'NA')
            if source_type not in TYPE_PRIORITY:
                print(f"Warning: unknown type '{source_type}', treating as lowest priority")
            rrna = ('Yes' if record.id in rna_sequences else 'No') if rna_sequences else 'not-provided'

            if h not in seen:
                seen[h] = {
                    'record': record,
                    'original_name': original_name,
                    'type': source_type,
                    'biomes': {biome},
                    'seq_id': record.id,
                    'description': record.description,
                    'rrna': rrna,
                    'quality': quality_data.get(record.id),
                    'viral_type': ('plasmid' if 'plasmid' in original_name or source_type == 'third_party_plasmid' else 'virus')
                }
            else:
                entry = seen[h]
                # Always merge biomes
                entry['biomes'].add(biome)

                # Replace record if the new source has higher priority (lower number)
                current_priority = TYPE_PRIORITY.get(entry['type'], 99)
                new_priority = TYPE_PRIORITY.get(source_type, 99)
                if new_priority < current_priority:
                    entry['record'] = record
                    entry['original_name'] = original_name
                    entry['type'] = source_type
                    entry['seq_id'] = record.id
                    entry['description'] = record.description
                    entry['rrna'] = rrna
                    entry['quality'] = quality_data.get(original_name)
    return seen


def write_final_files(
    output_fna: str,
    output_gff: str,
    output_tsv: str,
    filtered_fna: str,
    filtered_tsv: str,
    seen: dict[str, dict],
    gff_data: dict[str, list[str]],
    attr_id_to_seq_id: dict[str, str],
) -> None:
    """Write all four output files from the deduplicated record set.

    Writes every unique sequence to ``output_fna`` / ``output_tsv``, and
    sequences that pass ``passes_filter`` to ``filtered_fna`` / ``filtered_tsv``.
    GFF records for each filtered sequence are written to ``output_gff``.

    Args:
        output_fna: Path for the full deduplicated FASTA.
        output_gff: Path for the filtered GFF (sequences passing the filter only).
        output_tsv: Path for the full metadata TSV.
        filtered_fna: Path for the filtered FASTA.
        filtered_tsv: Path for the filtered metadata TSV.
        seen: Dict returned by ``choose_seqs``.
        gff_data: Dict of seq_id -> list of raw GFF lines (from ``read_input_gff``).
        attr_id_to_seq_id: Dict of attribute ID -> seq_id (from ``read_input_gff``).
    """
    # Write outputs
    records_written = 0
    records_filtered = 0
    gff_records_written = 0
    gff_seqs_written = 0
    with open(output_fna, 'w') as out_fna, \
            open(output_tsv, 'w') as out_tsv, \
            open(filtered_fna, 'w') as out_fna_f, \
            open(filtered_tsv, 'w') as out_tsv_f, \
            open(output_gff, 'w') as filt_gff:
        quality_header = '\t'.join(QUALITY_COLUMNS)
        header = (
            f"sequence_id\toriginal_name\tdescription\ttype\tbiomes\t"
            f"sequence_length\trrna\tsequence_sha256\t{quality_header}\n"
        )
        out_tsv.write(header)
        out_tsv_f.write(header)

        for h, entry in seen.items():
            record = entry['record']
            new_name = record.description.replace('|', '-').replace(' ', '|')
            record.id = new_name
            record.description = new_name

            biomes_str = ','.join(sorted(entry['biomes']))
            q = entry['quality']
            quality_values = '\t'.join(q[col] for col in QUALITY_COLUMNS) if q else '\t'.join(
                'NA' for _ in QUALITY_COLUMNS)
            tsv_row = (
                f"{entry['seq_id']}\t"
                f"{entry['original_name']}\t"
                f"{entry['original_name'].replace('|', '-').replace(' ', '|')}\t"
                f"{entry['type']}\t"
                f"{biomes_str}\t"
                f"{len(record.seq)}\t"
                f"{entry['rrna']}\t"
                f"{h}\t"
                f"{quality_values}\n"
            )

            SeqIO.write([record], out_fna, "fasta")
            out_tsv.write(tsv_row)
            records_written += 1

            if passes_filter(entry):
                SeqIO.write([record], out_fna_f, "fasta")
                out_tsv_f.write(tsv_row)
                records_filtered += 1

                gff_records = gff_data.get(entry['seq_id'])
                if gff_records:
                    filt_gff.write(''.join(gff_records))
                    gff_records_written += len(gff_records)
                    gff_seqs_written += 1
                else:
                    print(f"{entry['seq_id']} has no GFF records")

    print(f"Total unique sequences written: {records_written}")
    print(f"Filtered sequences written: {records_filtered}")
    print(f"Written sequences {gff_seqs_written} into GFF")
    print(f"Total written lines to filtered GFF {gff_records_written}")
    if records_filtered != gff_seqs_written:
        print("Number of GFF records does not match number of sequences in fasta file. Exit")
        exit(1)

def main() -> None:
    """Deduplicate sequences across all input FNA files and write merged outputs.

    Reads each FNA file together with its source type and biome label.  When the
    same sequence (by SHA256) appears in multiple sources, the record from the
    highest-priority source (see TYPE_PRIORITY) is kept and biome labels from all
    sources are merged.  Optionally renames contigs via a mapping file and flags
    sequences that carry rRNA/tRNA/tmRNA annotations.

    Outputs:
        - Deduplicated FNA file with pipe-delimited FASTA headers.
        - TSV metadata table with columns:
          sequence_id, original_name, description, type, biomes,
          sequence_length, rrna, sequence_sha256,
          provirus, proviral_length, gene_count, viral_genes, host_genes,
          checkv_quality, miuvig_quality, completeness, completeness_method,
          contamination, kmer_freq, warnings.
    """
    args = parse_arguments()

    # Read mapping file(s): temporary_name -> original_name
    mapping: dict[str, str] = read_map(args.map) if args.map else {}

    # Set of sequence IDs that carry RNA annotations
    rna_sequences: set[str] = parse_rna_gff(args.rrna) if args.rrna else set()

    # CheckV quality data keyed by original contig name
    quality_data: dict[str, dict[str, str]] = read_quality(args.quality) if args.quality else {}

    # Read GFFs
    input_gff, attr_id_to_seq_id = read_input_gff(args.gff)

    seen = choose_seqs(args.fna, mapping, rna_sequences, quality_data)

    # Derive filtered output paths by prepending "filtered_" to the filename
    p_fna = Path(args.output_fna)
    p_tsv = Path(args.output_tsv)
    filtered_fna = p_fna.parent / f"filtered_{p_fna.name}"
    filtered_tsv = p_tsv.parent / f"filtered_{p_tsv.name}"

    write_final_files(args.output_fna, args.output_gff, args.output_tsv, filtered_fna, filtered_tsv, seen, input_gff, attr_id_to_seq_id)
    print(f"Sources of FNA processed: {len(args.fna)}")

if __name__ == '__main__':
    main()
