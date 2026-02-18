#!/usr/bin/env python3
import argparse
import os
import sys
from Bio import SeqIO

VIRAL_TYPES = ['viral_sequence', 'plasmid', 'prophage']

def parse_arguments():
    """Parse command line arguments."""
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

    return parser.parse_args()


def detect_protein_id(annotation_field):
    """
    Parse the annotation field (GFF column 9) to extract protein ID and coordinates.
    Possible cases:
    ID=MGYG_X|plasmid-1:123149;...
    ID=MGYG_X|viral_sequence-1:6186;...
    ID=MGYG_X|prophage-13145:73981;...
    """
    protein = annotation_field.split(';')[0].replace('ID=', '')
    protein_id = protein.split('|')[0]
    coords = protein.split('|')[1].split('-')[1]
    if ':' not in coords:
        print(f'Sequence coordinates were not detected in {annotation_field}')
        coords = None
    return protein_id, coords


def parse_gff(gff):
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
                stats_counts[seq_type] += 1
                sequence_id = cur_sequence_id
                protein_id, coords = detect_protein_id(annotation_field)
                protein_ids.append(protein_id)
                sequence_ids.append((sequence_id, coords, seq_type))
                record_found = True
                gff_records.append(line)
                continue
            if record_found:
                if cur_sequence_id == sequence_id:
                    if seq_type == 'CDS':
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


def check_path_exists(path_str):
    if not os.path.exists(path_str):
        print(f"Error: No catalogue path found {path_str}")
        sys.exit(1)


def grep_sequences(sequence_ids, fasta_file):
    """
    Extract sequences from a FASTA file by ID.

    sequence_ids can be:
      - a list of plain string IDs (e.g. protein IDs for .faa lookup)
      - a list of (id, coords, seq_type) tuples (for .fna lookup)
        When coords are provided, the subsequence [start:end] is extracted (1-based, inclusive).
        If coords is None, the full sequence is returned.
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


def process_catalogue(catalogue_path, output_path, catalogue_name):
    """
    Assuming structure on NFS: catalogue / species_rep[MGYG] / genome
    Files of interest:
    MGYG.gff - viral records and corresponding proteins
    MGYG.fna - viral nucleotide sequence
    MGYG.faa - viral proteins
    """
    check_path_exists(catalogue_path)
    reps = [item for item in os.listdir(catalogue_path) if item.startswith('MGYG')]

    os.makedirs(output_path, exist_ok=True)

    final_gff = os.path.join(output_path, catalogue_name + '_viral.gff')
    final_fna = os.path.join(output_path, catalogue_name + '_viral.fna')
    final_faa = os.path.join(output_path, catalogue_name + '_viral.faa')
    with open(final_gff, 'w') as out_gff, open(final_fna, 'w') as out_fna, open(final_faa, 'w') as out_faa:

        for rep in reps:
            print(f'Processing rep: {rep}')
            gff = os.path.join(catalogue_path, rep, 'genome', f'{rep}.gff')
            fna = os.path.join(catalogue_path, rep, 'genome', f'{rep}.fna')
            faa = os.path.join(catalogue_path, rep, 'genome', f'{rep}.faa')

            check_path_exists(gff)
            check_path_exists(faa)
            check_path_exists(fna)

            # find viral records in GFFs
            gff_records, sequence_ids, protein_ids = parse_gff(gff)
            # write filtered GFF
            for record in gff_records:
                out_gff.write(record + '\n')

            # find and write FNAs for chosen sequences
            chosen_records = grep_sequences(sequence_ids, fna)
            SeqIO.write(chosen_records, out_fna, "fasta")

            # find and write proteins for chosen sequences
            chosen_records = grep_sequences(protein_ids, faa)
            SeqIO.write(chosen_records, out_faa, "fasta")


def main():
    """Main entry point."""
    args = parse_arguments()

    for catalogue_path in args.catalogue_path:
        path_parts = catalogue_path.rstrip('/').split('/')
        catalogue_name = '_'.join(path_parts[-2:])
        print(f'Running search for {catalogue_name}')
        process_catalogue(catalogue_path=catalogue_path, output_path=args.output_path, catalogue_name=catalogue_name)


if __name__ == '__main__':
    main()
