#!/usr/bin/env python3
import argparse
import os
import sys
from Bio import SeqIO

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


def detect_protein_id(line):
    metadata = line.split('\t')[8]
    protein_id = metadata.split(';')[0].replace('ID=', '')
    return protein_id


def parse_gff(gff):
    gff_records = []
    viral_count, plasmid_count, phage_count = 0, 0, 0
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
            annotation_field = parts[8]
            if 'viral' in annotation_field or 'plasmid' in annotation_field or 'phage' in annotation_field:
                if 'viral' in annotation_field:
                    viral_count += 1
                elif 'plasmid' in annotation_field:
                    plasmid_count += 1
                elif 'phage' in annotation_field:
                    phage_count += 1
                sequence_id = cur_sequence_id
                sequence_ids.append(sequence_id)
                protein_ids.append(detect_protein_id(line))
                record_found = True
                gff_records.append(line)
                continue
            if record_found:
                if cur_sequence_id == sequence_id:
                    if 'CDS' in line:
                        gff_records.append(line)
                    else:
                        print(f'No CDS in {line}')
                else:
                    record_found = False
    print(f'Found total GFF {len(gff_records)} records')
    print(f'{len(sequence_ids)} sequences with viral: {viral_count}, plasmid: {plasmid_count}, phage: {phage_count}')
    print(f'protein sequences: {len(protein_ids)}')
    return gff_records, sequence_ids, protein_ids


def check_path_exists(path_str):
    if not os.path.exists(path_str):
        print(f"Error: No catalogue path found {path_str}")
        sys.exit(1)


def grep_sequences(sequence_ids, fna):
    id_set = set(sequence_ids)
    chosen_records = []
    with open(fna, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in id_set:
                chosen_records.append(record)
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
