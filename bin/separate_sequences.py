#!/usr/bin/env python

import argparse
import re


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter FASTA sequences by pattern and exclude RNA-containing sequences"
    )
    parser.add_argument(
        "-i", "--input", help="Input FASTA file", required=True
    )
    parser.add_argument(
        "-p", "--pattern", help="Pattern to match in sequence headers", required=True
    )
    parser.add_argument(
        "-g", "--rna-gff", help="GFF file with RNA annotations from barrnap", required=False
    )
    parser.add_argument(
        "-o", "--output", help="Output FASTA file", required=True
    )
    return parser.parse_args()


def parse_rna_gff(gff_file):
    """Parse barrnap GFF file and return set of sequence IDs containing RNA."""
    rna_sequences = set()
    if not gff_file:
        return rna_sequences

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


def separate_sequences(input_fasta, pattern, rna_sequences, output_fasta):
    """Filter FASTA sequences by pattern and exclude RNA-containing sequences."""
    regex = re.compile(pattern)
    kept_count = 0
    excluded_rna_count = 0
    pattern_filtered_count = 0

    with open(input_fasta, "r") as fin, open(output_fasta, "w") as fout:
        current_header = None
        current_seq_lines = []
        write_current = False

        for line in fin:
            if line.startswith(">"):
                # Write previous sequence if it passed filters
                if write_current and current_header:
                    fout.write(current_header)
                    fout.writelines(current_seq_lines)
                    kept_count += 1

                # Process new header
                current_header = line
                current_seq_lines = []
                header_text = line.strip()

                # Check if header matches pattern
                if regex.search(header_text):
                    # Extract sequence ID (first word after >)
                    seq_id = header_text[1:].split()[0]
                    # Check if sequence contains RNA
                    if seq_id in rna_sequences:
                        write_current = False
                        excluded_rna_count += 1
                    else:
                        write_current = True
                else:
                    write_current = False
                    pattern_filtered_count += 1
            else:
                current_seq_lines.append(line)

        # Write last sequence if it passed filters
        if write_current and current_header:
            fout.write(current_header)
            fout.writelines(current_seq_lines)
            kept_count += 1

    print(f"Sequences kept: {kept_count}")
    print(f"Sequences excluded (pattern mismatch): {pattern_filtered_count}")
    print(f"Sequences excluded (RNA detected): {excluded_rna_count}")


def main():
    args = parse_args()
    rna_sequences = parse_rna_gff(args.rna_gff)
    if args.rna_gff:
        print(f"Found {len(rna_sequences)} sequences with RNA annotations")
    separate_sequences(args.input, args.pattern, rna_sequences, args.output)


if __name__ == "__main__":
    main()
