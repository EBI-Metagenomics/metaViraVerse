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
        "--type",
        required=True,
        nargs='+',
        help="Source type for each FNA file: 'metagenome'/'genome'/'third_party_virus'/'third_party_plasmid'"
    )
    parser.add_argument(
        "--biome",
        required=True,
        nargs='+',
        help="Biome label for each FNA file"
    )
    parser.add_argument(
        "--output-fna",
        required=True,
        help="Path for output deduplicated FNA file"
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        help="Path for output metadata TSV file"
    )
    args = parser.parse_args()

    if len(args.fna) != len(args.type) or len(args.fna) != len(args.biome):
        parser.error("--fna, --type, and --biome must have the same number of values")

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


def read_map(map_files: list[str]) -> dict[str, str]:
    """Read one or more TSV mapping files from temporary to original contig names.

    Each file is expected to have a header row containing 'original' and 'temporary',
    followed by rows with (original_name, temporary_name) tab-separated columns.

    Args:
        map_files: Paths to mapping TSV files.

    Returns:
        Dict mapping temporary contig name -> original contig name.

    Raises:
        SystemExit: If a temporary name appears in more than one mapping entry.
    """
    mapping: dict[str, str] = {}
    for map_file in map_files:
        with open(map_file, 'r') as file_in:
            for line in file_in:
                if 'original' in line and 'temporary' in line:
                    continue
                parts = line.strip().split('\t')
                original = parts[0]
                temporary = parts[1]
                if temporary in mapping:
                    print(f'Mapping already exists {temporary}. Exit')
                    exit(1)
                mapping[temporary] = original
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
    if q is not None and q.get('checkv_quality') == 'Not-determined' and  entry['viral_type'] != 'plasmid':
        return False
    return True


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

    # Dict keyed by sequence hash -> {record, type, biomes, seq_id, description, rrna}
    seen: dict[str, dict] = {}

    for fna_file, source_type, biome in zip(args.fna, args.type, args.biome):
        if source_type not in TYPE_PRIORITY:
            print(f"Warning: unknown type '{source_type}', treating as lowest priority")

        for record in SeqIO.parse(fna_file, "fasta"):
            h = seq_hash(str(record.seq))
            original_name = mapping.get(record.id, record.id)
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

    # Derive filtered output paths by prepending "filtered_" to the filename
    p_fna = Path(args.output_fna)
    p_tsv = Path(args.output_tsv)
    filtered_fna = p_fna.parent / f"filtered_{p_fna.name}"
    filtered_tsv = p_tsv.parent / f"filtered_{p_tsv.name}"

    # Write outputs
    records_written = 0
    records_filtered = 0
    with (
        open(args.output_fna, 'w') as out_fna,
        open(args.output_tsv, 'w') as out_tsv,
        open(filtered_fna, 'w') as out_fna_f,
        open(filtered_tsv, 'w') as out_tsv_f,
    ):
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
            quality_values = '\t'.join(q[col] for col in QUALITY_COLUMNS) if q else '\t'.join('NA' for _ in QUALITY_COLUMNS)
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

    print(f"Total unique sequences written: {records_written}")
    print(f"Filtered sequences written: {records_filtered}")
    print(f"Sources processed: {len(args.fna)}")


if __name__ == '__main__':
    main()
