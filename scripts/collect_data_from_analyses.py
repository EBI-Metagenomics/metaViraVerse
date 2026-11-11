#!/usr/bin/env python3
"""
Script searches for viral and plasmid records in MGnify assembly analyses (mobilome-annotation
GFFs) and extracts corresponding nucleotide and protein sequences.

Input: TSV written by collect_mobilome_from_analyses.py (id, map_gff, cds_faa, contigs_fasta)
Output: per analysis — filtered GFF, FNA (nucleotide), and FAA (protein) files
        containing only viral_sequence, plasmid, and prophage records;
        samplesheet.csv for the pipeline and a report with counts per analysis.

Usage:
    python collect_data_from_analyses.py -i marine_mobilome.tsv -o /output/dir --biome marine
"""
import argparse
import csv
import gzip
import os
import sys
from typing import Optional, TextIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

VIRAL_TYPES = ['viral_sequence', 'plasmid', 'prophage']
SAMPLESHEET_COLUMNS = ['id', 'fna', 'gff', 'faa', 'source', 'biome']
REPORT_COLUMNS = ['id'] + VIRAL_TYPES + ['gff_cds', 'fna_sequences', 'faa_proteins']


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments with input, output_path, biome and source.
    """
    parser = argparse.ArgumentParser(
        description="Script searches for viral records in analyses mobilome GFFs and greps corresponding nucleotide and protein sequences.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Script takes as input the TSV from collect_mobilome_from_analyses.py and creates per analysis:
        output: <id>_viral.gff, <id>_viral.fna and <id>_viral.faa
        plus samplesheet.csv and report.tsv
        """
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="TSV with columns id, map_gff, cds_faa, contigs_fasta"
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Path to save results (filtered gff, fna, faa, samplesheet, report)",
        default='.'
    )
    parser.add_argument(
        "--biome",
        help="Biome written to the samplesheet for every analysis",
        default='marine'
    )
    parser.add_argument(
        "--source",
        help="Source written to the samplesheet for every analysis",
        choices=['genome', 'metagenome'],
        default='metagenome'
    )
    return parser.parse_args()


def open_file(path: str) -> TextIO:
    """Open a plain or gzip-compressed text file."""
    with open(path, 'rb') as handle:
        is_gzip = handle.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if is_gzip else open(path, 'r')


def get_attribute_id(annotation_field: str) -> str:
    """Return the ID= value from GFF column 9.

    Examples:
        ID=ERZ12300306_1|viral_sequence-1:26299;mobile_element_type=phage -> ERZ12300306_1|viral_sequence-1:26299
        ID=ERZ12300306_1_9;inference=...                                    -> ERZ12300306_1_9
    """
    for part in annotation_field.split(';'):
        if part.startswith('ID='):
            return part[len('ID='):]
    return ''


def parse_gff(gff: str) -> tuple[list[dict], dict[str, int]]:
    """Parse a mobilome GFF and extract viral/plasmid/prophage records with their CDS lines.

    CDS lines are assigned to a record when they are on the same contig and lie within the
    record's coordinates (so a prophage gets only its own CDSs, not the whole contig's).

    Args:
        gff: Path to the GFF file (plain or gzip).

    Returns:
        A tuple of:
            - regions: list of dicts with contig, start, end, seq_type, id (GFF ID=), line
              and cds (list of (protein_id, line)), in GFF order
            - stats_counts: number of records per type in VIRAL_TYPES
    """
    stats_counts = {vtype: 0 for vtype in VIRAL_TYPES}
    regions = []
    cds_by_contig = {}
    with open_file(gff) as file_in:
        for line in file_in:
            line = line.rstrip('\n')
            if line.startswith('##FASTA'):
                break
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            contig, seq_type = parts[0], parts[2]
            if seq_type in VIRAL_TYPES:
                # example: ERZ12300306_1  VIRify  viral_sequence  1  26299  ...  ID=ERZ12300306_1|viral_sequence-1:26299;...
                stats_counts[seq_type] += 1
                regions.append({
                    'contig': contig,
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'seq_type': seq_type,
                    'id': get_attribute_id(parts[8]),
                    'line': line,
                    'cds': [],
                })
            elif seq_type == 'CDS':
                cds_by_contig.setdefault(contig, []).append((int(parts[3]), int(parts[4]), get_attribute_id(parts[8]), line))

    for region in regions:
        for start, end, protein_id, line in cds_by_contig.get(region['contig'], []):
            if start >= region['start'] and end <= region['end']:
                region['cds'].append((protein_id, line))
    return regions, stats_counts


def grep_regions(regions: list[dict], fasta_file: str) -> list[SeqRecord]:
    """Extract the region of each record from its contig.

    The contig is matched by the first word of the FASTA header. The output header is the GFF
    ID of the record (e.g. ERZ12300306_1|viral_sequence-1:26299): bin/rename_contigs.py matches
    FNA headers to the GFF by this ID and classifies sequences by the type in it.

    Args:
        regions: records from parse_gff().
        fasta_file: contigs FASTA (plain or gzip).

    Returns:
        List of sliced SeqRecords, in the order the contigs appear in the FASTA.
    """
    contig_to_regions = {}
    for region in regions:
        contig_to_regions.setdefault(region['contig'], []).append(region)

    chosen_records = []
    with open_file(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contig = record.description.split(' ')[0]
            for region in contig_to_regions.get(contig, []):
                sub_record = record[region['start'] - 1:region['end']]  # 1-based inclusive to 0-based slice
                sub_record.id = region['id'] or f"{contig}|{region['seq_type']}-{region['start']}:{region['end']}"
                sub_record.description = ''
                chosen_records.append(sub_record)
    return chosen_records


def grep_proteins(protein_ids: set[str], fasta_file: str) -> list[SeqRecord]:
    """Extract proteins whose ID (first word of the header) is in protein_ids.

    Header example: >ERZ12300306_1_9 # 7018 # 21054 # -1 # ID=1_9;partial=00;start_type=ATG;...

    Args:
        protein_ids: CDS IDs from the GFF.
        fasta_file: CDS protein FASTA (plain or gzip).

    Returns:
        List of matching SeqRecords.
    """
    chosen_records = []
    with open_file(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            protein_id = record.description.split(' ')[0]
            if protein_id in protein_ids:
                record.id = protein_id
                record.description = ''
                chosen_records.append(record)
    return chosen_records


def check_path_exists(path_str: str) -> bool:
    """Warn if a file does not exist."""
    if not os.path.exists(path_str):
        print(f"Warning: no file found {path_str}")
        return False
    return True


def process_analysis(row: dict, output_path: str) -> Optional[dict]:
    """Extract viral/plasmid/prophage data of one analysis.

    Output files created in output_path:
        {id}_viral.gff  — record lines, each followed by its CDS lines
        {id}_viral.fna  — nucleotide sequences (sliced to record coordinates)
        {id}_viral.faa  — protein sequences for CDS within records

    Args:
        row: input TSV row with id, map_gff, cds_faa, contigs_fasta.
        output_path: Directory where output files will be written.

    Returns:
        Report row with counts, or None if an input file is missing.
    """
    analysis_id = row['id']
    print(f'Processing analysis: {analysis_id}')
    if not all(check_path_exists(row[col]) for col in ('map_gff', 'cds_faa', 'contigs_fasta')):
        return None

    regions, stats_counts = parse_gff(row['map_gff'])
    protein_ids = {protein_id for region in regions for protein_id, _ in region['cds']}

    final_gff = os.path.join(output_path, f'{analysis_id}_viral.gff')
    final_fna = os.path.join(output_path, f'{analysis_id}_viral.fna')
    final_faa = os.path.join(output_path, f'{analysis_id}_viral.faa')

    # rename_contigs.py assigns CDS lines to the preceding record, so keep them together
    with open(final_gff, 'w') as out_gff:
        for region in regions:
            out_gff.write(region['line'] + '\n')
            for _, line in region['cds']:
                out_gff.write(line + '\n')

    fna_records = grep_regions(regions, row['contigs_fasta']) if regions else []
    SeqIO.write(fna_records, final_fna, "fasta")

    faa_records = grep_proteins(protein_ids, row['cds_faa']) if protein_ids else []
    SeqIO.write(faa_records, final_faa, "fasta")

    if len(fna_records) != len(regions):
        print(f'Warning: {analysis_id}: {len(regions)} records in GFF, {len(fna_records)} sequences found in {row["contigs_fasta"]}')
    if len(faa_records) != len(protein_ids):
        print(f'Warning: {analysis_id}: {len(protein_ids)} CDS in GFF, {len(faa_records)} proteins found in {row["cds_faa"]}')

    stats_str = ', '.join(f'{vtype}: {count}' for vtype, count in stats_counts.items())
    print(f'Stats: {len(regions)} records with {stats_str}; {len(faa_records)} proteins')
    return {
        'id': analysis_id,
        **stats_counts,
        'gff_cds': sum(len(region['cds']) for region in regions),
        'fna_sequences': len(fna_records),
        'faa_proteins': len(faa_records),
    }


def main() -> None:
    """Main entry point. Processes each analysis and writes the samplesheet and report."""
    args = parse_arguments()
    os.makedirs(args.output_path, exist_ok=True)

    with open(args.input, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))

    samplesheet_path = os.path.join(args.output_path, 'samplesheet.csv')
    report_path = os.path.join(args.output_path, 'report.tsv')
    with open(samplesheet_path, 'w', newline='') as out_samplesheet, open(report_path, 'w', newline='') as out_report:
        samplesheet = csv.DictWriter(out_samplesheet, fieldnames=SAMPLESHEET_COLUMNS)
        samplesheet.writeheader()
        report = csv.DictWriter(out_report, fieldnames=REPORT_COLUMNS, delimiter='\t', lineterminator='\n')
        report.writeheader()

        for row in rows:
            counts = process_analysis(row, args.output_path)
            if counts is None:
                continue
            report.writerow(counts)
            if counts['fna_sequences'] == 0:
                print(f"Warning: {row['id']}: no viral sequences, prophages or plasmids, not added to samplesheet")
                continue
            sample_id = row['id']
            samplesheet.writerow({
                'id': sample_id,
                'fna': os.path.abspath(os.path.join(args.output_path, f'{sample_id}_viral.fna')),
                'gff': os.path.abspath(os.path.join(args.output_path, f'{sample_id}_viral.gff')),
                'faa': os.path.abspath(os.path.join(args.output_path, f'{sample_id}_viral.faa')),
                'source': args.source,
                'biome': args.biome,
            })

    print(f'Samplesheet written to {samplesheet_path}')
    print(f'Report written to {report_path}')


if __name__ == '__main__':
    main()
