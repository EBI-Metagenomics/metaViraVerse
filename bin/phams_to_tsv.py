#!/usr/bin/env python3
"""Script to map pham IDs to gene headers using pham FASTAs as inputs"""

import argparse
import re
import sys
import pathlib
from collections import defaultdict

from phammseqs import Pham


def parse_args():
    """Parse commandline arguments."""
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument("-i", dest="indir", type=pathlib.Path, help="directory containing pham FASTAs")
    p.add_argument("-o", dest="outfile", type=pathlib.Path, help="TSV file to save phamID-to-gene mapping")
    p.add_argument("-m", "--matrix", dest="matrix_file", type=pathlib.Path,
                   help="TSV file to save pham-genome matrix (rows=phams, columns=genomes)")

    return p.parse_args()


def extract_genome_name(geneid):
    """
    Extract genome name from gene ID.
    Takes first part before space and removes trailing _number suffix.

    Example:
    'ERZ2627000.43-NODE-43|prophage-8330:15600_1 description' -> 'ERZ2627000.43-NODE-43|prophage-8330:15600'
    """
    # Take first part before space
    name = geneid.split(' ')[0]
    # Remove trailing _number suffix
    genome = re.sub(r'_\d+$', '', name)
    return genome


def main():
    """Commandline entry point."""
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    args = parse_args()

    # Read in all the phams
    phams = list()
    for pham_fasta in args.indir.iterdir():
        # All phams output by PhaMMseqs will have .faa extension
        if pham_fasta.suffix != ".faa":
            continue
        p = Pham()
        p.load(pham_fasta)
        phams.append(p)

    # Sort phams
    sorted_phams = sorted(phams, reverse=True)

    # Write the TSV (pham-gene mapping)
    with open(args.outfile, "w") as fh:
        for i, pham in enumerate(sorted_phams):
            for geneid in pham.geneids:
                fh.write(f"pham{i + 1}\t{geneid}\n")

    # Write pham-genome matrix if requested
    if args.matrix_file:
        # Collect all genomes and build pham-genome counts
        all_genomes = set()
        pham_genome_counts = []

        for i, pham in enumerate(sorted_phams):
            pham_name = f"pham{i + 1}"
            genome_counts = defaultdict(int)

            for geneid in pham.geneids:
                genome = extract_genome_name(geneid)
                genome_counts[genome] += 1
                all_genomes.add(genome)

            pham_genome_counts.append((pham_name, genome_counts))

        # Sort genomes for consistent column order
        sorted_genomes = sorted(all_genomes)

        # Write matrix
        with open(args.matrix_file, "w") as fh:
            # Header row
            fh.write("pham\t" + "\t".join(sorted_genomes) + "\n")

            # Data rows
            for pham_name, genome_counts in pham_genome_counts:
                counts = [str(genome_counts.get(genome, 0)) for genome in sorted_genomes]
                fh.write(f"{pham_name}\t" + "\t".join(counts) + "\n")

        print(f"Matrix written: {len(pham_genome_counts)} phams x {len(sorted_genomes)} genomes")

    print("Done!")


if __name__ == "__main__":
    main()
