#!/usr/bin/env python3
"""
Script searches for viral records in MGnify catalogue(s) GFF files and extracts
corresponding nucleotide and protein sequences.

Input: path(s) to catalogue(s) on NFS
Output: per catalogue — filtered GFF, FNA (nucleotide), and FAA (protein) files
        containing only viral_sequence, plasmid, and prophage records.

Usage:
    python collect_data_from_catalogues.py -p /path/to/catalogue1 /path/to/catalogue2 -o /output/dir
"""
import argparse
import csv
import os
import sys
import urllib.error
import urllib.request
from typing import Optional, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

VIRAL_TYPES = ['viral_sequence', 'plasmid', 'prophage']
MGNIFY_GENOMES_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes"
SAMPLESHEET_COLUMNS = ['id', 'fna', 'gff', 'faa', 'source', 'biome']


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments with catalogue_path (list[str]) and output_path (str).
    """
    parser = argparse.ArgumentParser(
        description="Script searches for viral records in catalogue(s) GFFs and greps corresponding nucleotide and protein sequences.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Script Takes as input path(s) to catalogue(s) and creates 3 files with all found viral_sequences and plasmids:
        output: catalogue_name_version.fna and catalogue_name_version.faa
        """
    )
    parser.add_argument(
        "-p",
        "--catalogue-path",
        required=True,
        help="Path to NFS location of catalogue(s)",
        nargs='+'
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=True,
        help="Path to save results (filtered gff, fna, faa)",
        default='.'
    )
    parser.add_argument(
        "--old",
        action="store_true",
        help="For old catalogues with separate _virify.gff"
    )

    return parser.parse_args()


def detect_sequence_id_and_coords(annotation_field: str) -> tuple[str, Optional[str]]:
    """Parse GFF column 9 (attributes) to extract the sequence ID and genomic coordinates.

    The annotation field follows the format: ID=<seq_id>[|<type>-<start>:<end>];...
    When the pipe-delimited region part is present, coordinates are extracted as "start:end".

    Examples of annotation_field values:
        ID=MGYG_X;...                              -> ("MGYG_X", None)
        ID=MGYG_X|plasmid-1:123149;...             -> ("MGYG_X", "1:123149")
        ID=MGYG_X|viral_sequence-1:6186;...        -> ("MGYG_X", "1:6186")
        ID=MGYG_X|prophage-13145:73981;...         -> ("MGYG_X", "13145:73981")

    Args:
        annotation_field: The 9th column (attributes) from a GFF line.

    Returns:
        A tuple of (sequence_id, coords) where coords is "start:end" string
        or None if coordinates could not be parsed.
    """
    protein = annotation_field.split(';')[0].replace('ID=', '')
    protein_id = protein.split('|')[0]
    try:
        coords = protein.split('|')[1].split('-')[1]
    except (IndexError, ValueError):
        # No pipe-delimited region info (e.g. plain "ID=MGYG_X;...")
        coords = None
    return protein_id, coords


def parse_gff(gff: str) -> tuple[list[str], list[tuple[str, Optional[str], str]], list[str]]:
    """Parse a GFF file and extract viral/plasmid/prophage records with their CDS lines.

    Scans the GFF for lines where column 3 (type) matches one of VIRAL_TYPES.
    For each matching record, collects the record line and all subsequent CDS lines
    belonging to the same sequence ID.

    Args:
        gff: Path to the GFF file.

    Returns:
        A tuple of:
            - gff_records: list of raw GFF lines (viral/plasmid/prophage + their CDS lines)
            - sequence_ids: list of (sequence_id, coords, seq_type) tuples for nucleotide lookup.
              coords is "start:end" or None; seq_type is one of VIRAL_TYPES.
            - protein_ids: list of protein sequence IDs (str) extracted from CDS annotation fields.
    """
    stats_counts = {vtype: 0 for vtype in VIRAL_TYPES}
    gff_records = []
    sequence_ids, protein_ids = [], []
    with open(gff, 'r') as file_in:
        record_found = False
        sequence_id = None
        for line in file_in:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            cur_sequence_id = parts[0]
            seq_type = parts[2]
            annotation_field = parts[8]
            if seq_type in VIRAL_TYPES:
                # example of annotation_field: ID=MGYG000516801_121|plasmid-1:24987;mobile_element_type=plasmid
                # ID corresponds to nucleotide sequence region
                stats_counts[seq_type] += 1
                sequence_id = cur_sequence_id
                seq_id, coords = detect_sequence_id_and_coords(annotation_field)
                sequence_ids.append((sequence_id, coords, seq_type))
                record_found = True
                gff_records.append(line)
                continue
            if record_found:
                if cur_sequence_id == sequence_id:
                    if seq_type == 'CDS':
                        # example: MGYG000516801_371	Prodigal:002006	CDS	29	322	.	+	0
                        # ID=MGYG000516801_07313;inference=ab initio prediction:Prodigal:002006;locus_tag=...
                        # MGYG000516801_07313 is protein record
                        protein_id, coords = detect_sequence_id_and_coords(annotation_field)
                        protein_ids.append(protein_id)
                        gff_records.append(line)
                    else:
                        print(f'No CDS in {line}')
                else:
                    record_found = False
    print(f'Stats: found total GFF {len(gff_records)} records')
    stats_str = ', '.join(f'{vtype}: {count}' for vtype, count in stats_counts.items())
    print(f'Stats: {len(sequence_ids)} sequences with {stats_str}')
    print(f'Stats: protein sequences: {len(protein_ids)}')
    return gff_records, sequence_ids, protein_ids


def check_path_exists(path_str: str) -> None:
    """Verify that a file or directory exists, exit with error if not.

    Args:
        path_str: Path to check.
    """
    if not os.path.exists(path_str):
        print(f"Error: No catalogue path found {path_str}")
        sys.exit(1)


def grep_sequences(
    sequence_ids: Union[list[str], list[tuple[str, Optional[str], str]]],
    fasta_file: str
) -> list[SeqRecord]:
    """Extract sequences from a FASTA file, optionally slicing by coordinates.

    Supports two input modes:
      1. Plain string IDs — returns full matching records (used for protein .faa lookup).
      2. Tuples of (id, coords, seq_type) — when coords is "start:end", extracts the
         subsequence at those 1-based inclusive positions. The output record description
         is set to "{record_id} {seq_type}|{start}:{end}". When coords is None, returns
         the full sequence. Used for nucleotide .fna lookup.

    A single sequence ID may appear multiple times with different coordinates (e.g. multiple
    viral regions on the same contig); all regions are extracted.

    Args:
        sequence_ids: Either a list of sequence ID strings, or a list of
            (sequence_id, coords, seq_type) tuples.
        fasta_file: Path to the FASTA file (.fna or .faa) to search.

    Returns:
        List of matching Bio.SeqRecord objects (full or sliced).
    """
    # Detect whether we have tuples or plain strings
    if sequence_ids and isinstance(sequence_ids[0], tuple):
        # Build a dict: seq_id -> list of (coords, seq_type) pairs (one ID can appear multiple times)
        id_to_regions = {}
        for seq_id, coords, seq_type in sequence_ids:
            id_to_regions.setdefault(seq_id, []).append((coords, seq_type))
        lookup_ids = set(id_to_regions)
    else:
        lookup_ids = set(sequence_ids)
        id_to_regions = None

    chosen_records = []
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id not in lookup_ids:
                continue
            if id_to_regions is None:
                chosen_records.append(record)
            else:
                for coords, seq_type in id_to_regions[record.id]:
                    if coords is None:
                        chosen_records.append(record)
                    else:
                        start, end = coords.split(':')
                        start, end = int(start) - 1, int(end)  # 1-based inclusive to 0-based slice
                        sub_record = record[start:end]
                        sub_record.id = record.id
                        sub_record.description = f"{record.id} {seq_type}|{start + 1}:{end}"
                        chosen_records.append(sub_record)
    return chosen_records


def process_catalogue(catalogue_path: str, output_path: str, catalogue_name: str, old: bool) -> None:
    """Process a single MGnify catalogue: extract viral/plasmid/prophage data.

    Iterates over all species representatives (MGYG* directories) in the catalogue.
    For each representative, parses its GFF to find viral records, then extracts
    the corresponding nucleotide and protein sequences from FNA/FAA files.

    Expected directory structure:
        catalogue_path/
            MGYG.../
                genome/
                    MGYG....gff
                    MGYG....fna
                    MGYG....faa

    Output files created in output_path:
        {catalogue_name}_viral.gff  — filtered GFF lines
        {catalogue_name}_viral.fna  — nucleotide sequences (sliced to viral region coordinates)
        {catalogue_name}_viral.faa  — protein sequences for CDS within viral regions

    Args:
        catalogue_path: Path to the catalogue directory containing MGYG* subdirectories.
        output_path: Directory where output files will be written (created if needed).
        catalogue_name: Base name used for output file naming.
        old: flag to use different structure of folders to find viral results
    """
    check_path_exists(catalogue_path)

    os.makedirs(output_path, exist_ok=True)
    final_gff = os.path.join(output_path, catalogue_name + '_viral.gff')
    final_fna = os.path.join(output_path, catalogue_name + '_viral.fna')
    final_faa = os.path.join(output_path, catalogue_name + '_viral.faa')

    if old:
        reps_short = [item for item in os.listdir(catalogue_path) if item.startswith('MGYG')]
        reps = []
        for rep_short in reps_short:
            reps.extend([f'{rep_short}/{item}' for item in os.listdir(os.path.join(catalogue_path, rep_short)) if item.startswith('MGYG')])
        gff_prefix = '_virify.gff'
    else:
        reps = [item for item in os.listdir(catalogue_path) if item.startswith('MGYG')]
        gff_prefix = '.gff'

    with open(final_gff, 'w') as out_gff, open(final_fna, 'w') as out_fna, open(final_faa, 'w') as out_faa:

        for rep in reps:
            print(f'Processing rep: {rep}')
            if '/' in rep:
                genome_name = rep.split('/')[-1]
            else:
                genome_name = rep
            gff = os.path.join(catalogue_path, rep, 'genome', f'{genome_name}{gff_prefix}')
            fna = os.path.join(catalogue_path, rep, 'genome', f'{genome_name}.fna')
            faa = os.path.join(catalogue_path, rep, 'genome', f'{genome_name}.faa')

            check_path_exists(gff)
            check_path_exists(faa)
            check_path_exists(fna)

            # Parse GFF to find viral/plasmid/prophage records and their CDS proteins
            gff_records, sequence_ids, protein_ids = parse_gff(gff)

            # Write filtered GFF lines
            for record in gff_records:
                out_gff.write(record + '\n')

            # Extract and write nucleotide sequences (sliced to viral region coordinates)
            chosen_records = grep_sequences(sequence_ids, fna)
            SeqIO.write(chosen_records, out_fna, "fasta")

            # Extract and write protein sequences for CDS within viral regions
            chosen_records = grep_sequences(protein_ids, faa)
            SeqIO.write(chosen_records, out_faa, "fasta")


def get_catalogue_name_and_version(catalogue_path: str) -> tuple[str, str]:
    """Extract a catalogue's name and version from its path.

    Assumes the catalogue path ends in "<catalogue_name>/<version>", which
    matches both the NFS layout and the corresponding path on the MGnify
    genomes FTP.

    Args:
        catalogue_path: Path to the catalogue directory, e.g. "/nfs/public/.../human-gut/v2.0".

    Returns:
        A tuple of (catalogue_name, version), e.g. ("human-gut", "v2.0").
    """
    path_parts = catalogue_path.rstrip('/').split('/')
    catalogue_name, version = path_parts[-2:]
    return catalogue_name, version


def fetch_metadata(catalogue_path: str, output_path: str) -> str:
    """Download a catalogue's genomes-all_metadata.tsv from the MGnify genomes FTP.

    Builds the download URL from the catalogue's name and version (the last
    two components of catalogue_path), e.g.:
        https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/<name>/<version>/genomes-all_metadata.tsv

    Args:
        catalogue_path: Path to the catalogue directory (used only to derive name/version).
        output_path: Directory to save the downloaded metadata file to.

    Returns:
        Path to the downloaded metadata TSV file.
    """
    catalogue_name, version = get_catalogue_name_and_version(catalogue_path)
    url = f"{MGNIFY_GENOMES_BASE_URL}/{catalogue_name}/{version}/genomes-all_metadata.tsv"
    dest_path = os.path.join(output_path, f'{catalogue_name}_{version}_metadata.tsv')

    print(f'Fetching metadata for {catalogue_name}/{version} from {url}')
    try:
        urllib.request.urlretrieve(url, dest_path)
    except urllib.error.URLError as error:
        print(f"Error: Could not fetch metadata from {url} ({error})")
        sys.exit(1)

    return dest_path


def concatenate_metadata(metadata_files: list[str], output_file: str) -> None:
    """Concatenate multiple catalogue metadata TSVs into one file with a single header.

    Args:
        metadata_files: Paths to per-catalogue metadata TSV files, in order.
        output_file: Path to write the combined TSV to.
    """
    with open(output_file, 'w') as out_handle:
        for i, metadata_file in enumerate(metadata_files):
            with open(metadata_file, 'r') as in_handle:
                header = in_handle.readline()
                if i == 0:
                    out_handle.write(header)
                out_handle.writelines(in_handle.readlines())


def fetch_catalogues_metadata(catalogue_paths: list[str], output_path: str) -> str:
    """Fetch genomes-all_metadata.tsv for each catalogue and concatenate them.

    Args:
        catalogue_paths: Paths to catalogue directories.
        output_path: Directory to save per-catalogue and combined metadata files to.

    Returns:
        Path to the combined metadata TSV file.
    """
    metadata_files = [fetch_metadata(catalogue_path, output_path) for catalogue_path in catalogue_paths]

    combined_metadata = os.path.join(output_path, 'genomes-all_metadata.tsv')
    concatenate_metadata(metadata_files, combined_metadata)
    print(f'Combined metadata for {len(metadata_files)} catalogue(s) written to {combined_metadata}')
    return combined_metadata


def generate_samplesheet(catalogue_paths: list[str], output_path: str) -> str:
    """Generate a pipeline input samplesheet for the processed catalogues.

    Follows the schema in assets/schema_input.json: one row per catalogue,
    pointing at the filtered GFF/FNA/FAA files written by process_catalogue().
    `id` is "<catalogue_name>_<version>", `source` is always "genome" (MAG
    catalogue data), and `biome` is the catalogue name without its version.

    Args:
        catalogue_paths: Paths to catalogue directories.
        output_path: Directory containing the per-catalogue output files;
            the samplesheet is written here too.

    Returns:
        Path to the generated samplesheet.csv.
    """
    samplesheet_path = os.path.join(output_path, 'samplesheet.csv')
    with open(samplesheet_path, 'w', newline='') as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=SAMPLESHEET_COLUMNS)
        writer.writeheader()
        for catalogue_path in catalogue_paths:
            catalogue_name, version = get_catalogue_name_and_version(catalogue_path)
            sample_id = f'{catalogue_name}_{version}'
            writer.writerow({
                'id': sample_id,
                'fna': os.path.abspath(os.path.join(output_path, f'{sample_id}_viral.fna')),
                'gff': os.path.abspath(os.path.join(output_path, f'{sample_id}_viral.gff')),
                'faa': os.path.abspath(os.path.join(output_path, f'{sample_id}_viral.faa')),
                'source': 'genome',
                'biome': catalogue_name,
            })

    print(f'Samplesheet written to {samplesheet_path}')
    return samplesheet_path


def main() -> None:
    """Main entry point. Parses arguments and processes each catalogue path."""
    args = parse_arguments()

    os.makedirs(args.output_path, exist_ok=True)

    for catalogue_path in args.catalogue_path:
        catalogue_name = '_'.join(get_catalogue_name_and_version(catalogue_path))
        print(f'Running search for {catalogue_name}')
        process_catalogue(catalogue_path=catalogue_path, output_path=args.output_path, catalogue_name=catalogue_name, old=args.old)

    fetch_catalogues_metadata(catalogue_paths=args.catalogue_path, output_path=args.output_path)
    generate_samplesheet(catalogue_paths=args.catalogue_path, output_path=args.output_path)


if __name__ == '__main__':
    main()
