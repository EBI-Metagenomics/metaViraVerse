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
import hashlib

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
        "--type",
        required=True,
        nargs='+',
        help="Source type for each FNA file: 'metagenome' or 'genome'"
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


def main() -> None:
    args = parse_arguments()

    # Dict keyed by sequence hash -> {record, type, biomes, seq_id, description}
    seen = {}

    for fna_file, source_type, biome in zip(args.fna, args.type, args.biome):
        if source_type not in TYPE_PRIORITY:
            print(f"Warning: unknown type '{source_type}', treating as lowest priority")

        for record in SeqIO.parse(fna_file, "fasta"):
            h = seq_hash(str(record.seq))

            if h not in seen:
                # First time seeing this sequence
                seen[h] = {
                    'record': record,
                    'type': source_type,
                    'biomes': {biome},
                    'seq_id': record.id,
                    'description': record.description,
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
                    entry['type'] = source_type
                    entry['seq_id'] = record.id
                    entry['description'] = record.description

    # Write outputs
    records_written = 0
    with open(args.output_fna, 'w') as out_fna, open(args.output_tsv, 'w') as out_tsv:
        # TSV header
        out_tsv.write("sequence_id\tdescription\ttype\tbiomes\tsequence_sha256\tsequence_length\n")

        for h, entry in seen.items():
            record = entry['record']
            new_name = record.description.replace('|', '-').replace(' ', '|')
            record.id = new_name
            record.description = new_name
            SeqIO.write([record], out_fna, "fasta")

            biomes_str = ','.join(sorted(entry['biomes']))
            out_tsv.write(
                f"{entry['seq_id']}\t"
                f"{entry['description'].replace('|', '-').replace(' ', '|')}\t"
                f"{entry['type']}\t"
                f"{biomes_str}\t"
                f"{h}\t"
                f"{len(record.seq)}\n"
            )
            records_written += 1

    print(f"Total unique sequences written: {records_written}")
    print(f"Sources processed: {len(args.fna)}")


if __name__ == '__main__':
    main()
