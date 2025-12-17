#!/usr/bin/env python

"""
Combined Gene Caller: Merges gene predictions from Prodigal and PHANOTATE.

Logic:
1. If both callers predict the same ORF (≥overlap_threshold overlap) - keep Prodigal prediction
2. If PHANOTATE-only ORF, keep it if:
   - Minimum length (≥min_length bp)
   - No strong overlap with a longer Prodigal ORF
"""

import argparse
import re
import gzip
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class ORF:
    """Represents an Open Reading Frame."""
    contig: str
    start: int
    end: int
    strand: int
    header: str
    sequence: str
    source: str

    @property
    def length(self):
        return abs(self.end - self.start) + 1

    def coords(self):
        """Return normalized coordinates (start always < end)."""
        return (min(self.start, self.end), max(self.start, self.end))


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Combine gene predictions from Prodigal and PHANOTATE"
    )
    parser.add_argument(
        "-p", "--prodigal", help="Prodigal nucleotide FASTA file (.fna or .fna.gz)", required=True
    )
    parser.add_argument(
        "-t", "--phanotate", help="PHANOTATE FASTA file (.fna or .faa)", required=True
    )
    parser.add_argument(
        "-o", "--output", help="Output combined FASTA file", required=True
    )
    parser.add_argument(
        "--overlap-threshold", type=float, default=0.7,
        help="Minimum overlap fraction to consider ORFs as same (default: 0.7)"
    )
    parser.add_argument(
        "--min-length", type=int, default=90,
        help="Minimum length in bp for PHANOTATE-only ORFs (default: 90)"
    )
    parser.add_argument(
        "--max-overlap-with-longer", type=float, default=0.5,
        help="Maximum overlap allowed with longer Prodigal ORF for PHANOTATE-only (default: 0.5)"
    )
    return parser.parse_args()


def open_file(filepath):
    """Open regular or gzipped file."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def parse_prodigal_header(header):
    """
    Parse Prodigal FASTA header.
    Format: >contig_1 # 3 # 962 # 1 # ID=1_1;partial=00;...

    Example:
    >ERZ2627000.43-NODE-43-length-15729-cov-7.081090|prophage-8330:15600_1 # 3 # 962 # 1 # ID=1_1;partial=10;start_type=Edge;...
    """
    parts = header.split(' # ')
    if len(parts) >= 4:
        contig_with_orf = parts[0][1:].strip()  # Remove '>' and get contig name with ORF number
        start = int(parts[1])
        end = int(parts[2])
        strand = int(parts[3])

        # Extract base contig name by removing the ORF number suffix (_1, _2, etc.)
        # e.g., "contig|prophage-8330:15600_1" -> "contig|prophage-8330:15600"
        contig = re.sub(r'_\d+$', '', contig_with_orf)

        return contig, start, end, strand
    return None


def parse_phanotate_header(header):
    """
    Parse PHANOTATE FASTA header.
    Format: >contig_CDS_[start..end] or >contig_CDS_[complement(start..end)]

    Examples:
    >ERZ2627000.43-NODE-43-length-15729-cov-7.081090|prophage-8330:15600_CDS_[3..962] [note=score:-3.430764E+05]
    >ERZ2627000.43-NODE-43-length-15729-cov-7.081090|prophage-8330:15600_CDS_[complement(1772..2737)] [note=score:-7.436542E+08]
    """
    header_text = header[1:].strip()  # Remove '>'

    # Format: >contig_CDS_[complement(start..end)] [note=...]
    match = re.match(r'(.+?)_CDS_\[complement\((\d+)\.\.(\d+)\)\]', header_text)
    if match:
        contig = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        strand = -1
        return contig, start, end, strand

    # Format: >contig_CDS_[start..end] [note=...]
    match = re.match(r'(.+?)_CDS_\[(\d+)\.\.(\d+)\]', header_text)
    if match:
        contig = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        strand = 1
        return contig, start, end, strand

    return None


def parse_fasta(filepath, source, parser_func):
    """Parse FASTA file and return list of ORFs."""
    orfs = []
    current_header = None
    current_seq_lines = []

    with open_file(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Process previous sequence
                if current_header:
                    parsed = parser_func(current_header)
                    if parsed:
                        contig, start, end, strand = parsed
                        orf = ORF(
                            contig=contig,
                            start=start,
                            end=end,
                            strand=strand,
                            header=current_header,
                            sequence=''.join(current_seq_lines),
                            source=source
                        )
                        orfs.append(orf)
                current_header = line
                current_seq_lines = []
            else:
                current_seq_lines.append(line)

        # Process last sequence
        if current_header:
            parsed = parser_func(current_header)
            if parsed:
                contig, start, end, strand = parsed
                orf = ORF(
                    contig=contig,
                    start=start,
                    end=end,
                    strand=strand,
                    header=current_header,
                    sequence=''.join(current_seq_lines),
                    source=source
                )
                orfs.append(orf)

    return orfs


def calculate_overlap(orf1, orf2):
    """Calculate overlap fraction between two ORFs."""
    start1, end1 = orf1.coords()
    start2, end2 = orf2.coords()

    # Calculate overlap
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    if overlap_start >= overlap_end:
        return 0.0

    overlap_length = overlap_end - overlap_start

    # Return overlap as fraction of the shorter ORF
    min_length = min(orf1.length, orf2.length)
    return overlap_length / min_length


def combine_predictions(prodigal_orfs, phanotate_orfs, overlap_threshold, min_length, max_overlap_with_longer):
    """
    Combine predictions from both callers.

    Logic:
    1. Keep all Prodigal predictions
    2. For each PHANOTATE ORF:
       - If it overlaps ≥overlap_threshold with a Prodigal ORF -> skip (use Prodigal)
       - If no strong overlap and length ≥ min_length -> keep
       - If overlaps < max_overlap_with_longer with longer Prodigal ORF -> keep
    """
    combined = []
    stats = {
        'prodigal_kept': 0,
        'phanotate_kept': 0,
        'phanotate_skipped_overlap': 0,
        'phanotate_skipped_short': 0,
        'phanotate_skipped_longer_overlap': 0
    }

    # Group Prodigal ORFs by contig
    prodigal_by_contig = defaultdict(list)
    for orf in prodigal_orfs:
        prodigal_by_contig[orf.contig].append(orf)

    # Keep all Prodigal predictions
    combined.extend(prodigal_orfs)
    stats['prodigal_kept'] = len(prodigal_orfs)

    # Process PHANOTATE predictions
    for phan_orf in phanotate_orfs:
        # Check length requirement
        if phan_orf.length < min_length:
            stats['phanotate_skipped_short'] += 1
            continue

        # Find potential overlapping Prodigal ORFs on the same contig
        candidate_prodigal = prodigal_by_contig.get(phan_orf.contig, [])

        # Check for overlaps with Prodigal ORFs
        has_strong_overlap = False
        overlaps_longer_orf = False

        for prod_orf in candidate_prodigal:
            # Only compare ORFs on the same strand
            if phan_orf.strand != prod_orf.strand:
                continue

            overlap = calculate_overlap(phan_orf, prod_orf)

            # Check if this is essentially the same ORF
            if overlap >= overlap_threshold:
                has_strong_overlap = True
                break

            # Check if PHANOTATE ORF overlaps significantly with a longer Prodigal ORF
            if prod_orf.length > phan_orf.length and overlap >= max_overlap_with_longer:
                overlaps_longer_orf = True

        if has_strong_overlap:
            stats['phanotate_skipped_overlap'] += 1
        elif overlaps_longer_orf:
            stats['phanotate_skipped_longer_overlap'] += 1
        else:
            combined.append(phan_orf)
            stats['phanotate_kept'] += 1

    return combined, stats


def write_output(orfs, output_file):
    """Write combined ORFs to output FASTA file."""
    with open(output_file, 'w') as f:
        for orf in orfs:
            # Add source tag to header
            header = orf.header
            if not header.endswith(f'[{orf.source}]'):
                header = f"{header} [{orf.source}]"
            f.write(f"{header}\n")
            # Write sequence in lines of 60 characters
            seq = orf.sequence
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")


def main():
    args = parse_args()

    print(f"Parsing Prodigal predictions from {args.prodigal}")
    prodigal_orfs = parse_fasta(args.prodigal, 'prodigal', parse_prodigal_header)
    print(f"  Found {len(prodigal_orfs)} ORFs")

    print(f"Parsing PHANOTATE predictions from {args.phanotate}")
    phanotate_orfs = parse_fasta(args.phanotate, 'phanotate', parse_phanotate_header)
    print(f"  Found {len(phanotate_orfs)} ORFs")

    print(f"\nCombining predictions...")
    print(f"  Overlap threshold: {args.overlap_threshold}")
    print(f"  Minimum length for PHANOTATE-only: {args.min_length} bp")
    print(f"  Max overlap with longer ORF: {args.max_overlap_with_longer}")

    combined, stats = combine_predictions(
        prodigal_orfs,
        phanotate_orfs,
        args.overlap_threshold,
        args.min_length,
        args.max_overlap_with_longer
    )

    print(f"\nResults:")
    print(f"  Prodigal ORFs kept: {stats['prodigal_kept']}")
    print(f"  PHANOTATE ORFs kept (unique): {stats['phanotate_kept']}")
    print(f"  PHANOTATE ORFs skipped (same as Prodigal): {stats['phanotate_skipped_overlap']}")
    print(f"  PHANOTATE ORFs skipped (too short): {stats['phanotate_skipped_short']}")
    print(f"  PHANOTATE ORFs skipped (overlaps longer): {stats['phanotate_skipped_longer_overlap']}")
    print(f"  Total combined ORFs: {len(combined)}")

    write_output(combined, args.output)
    print(f"\nWrote combined predictions to {args.output}")


if __name__ == "__main__":
    main()
