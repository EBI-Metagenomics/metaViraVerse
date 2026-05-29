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
        "-o", "--output", help="Output FASTA file", required=True
    )
    return parser.parse_args()



def separate_sequences(input_fasta, pattern, output_fasta):
    """Filter FASTA sequences by pattern and exclude RNA-containing sequences."""
    regex = re.compile(pattern)
    kept_count = 0
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


def main():
    args = parse_args()
    separate_sequences(args.input, args.pattern, args.output)


if __name__ == "__main__":
    main()
