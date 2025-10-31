#!/usr/bin/env python3
import argparse
import os


def parse_attributes(attr_str):
    """Parse the attributes column from a GFF line into a dictionary."""
    attrs = {}
    for part in attr_str.split(";"):
        if "=" in part:
            key, val = part.split("=", 1)
            attrs[key.strip().lower()] = val.strip()
    return attrs


def read_mapfile(mapfile):
    mapping = {}
    with open(mapfile, 'r') as file_in:
        for line in file_in:
            line = line.strip.split('t')
            tmp_name = line[1]
            original_name = line[0]
            mapping[tmp_name] = original_name
    return mapping


def extract_viral_taxonomy(viral_list_file, gff_file, output_file, mapfile):
    """
    Extract taxonomy info for viral sequences found in GFF.
    """
    mapping = None
    if mapfile:
        mapping = read_mapfile(mapfile)

    # Read viral sequence names
    viral_names = []
    with open(viral_list_file, "r") as f:
        for line in f:
            line = line.strip()
            name = line.strip().split('\t')[0]
            # if contigs were renamed then name in gff would not match -> returning it
            if mapping:
                name = name.replace(name, mapping[name])
            viral_names.append(name)

    results = {}
    found = set()

    with open(gff_file, "r") as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            seq_id = cols[0].strip()
            attr_str = cols[8]
            attrs = parse_attributes(attr_str)

            if seq_id in viral_names:
                taxonomy = attrs.get("taxonomy", "NA")
                results[seq_id] = taxonomy
                found.add(seq_id)

    # Ensure all viral names appear in output
    with open(output_file, "w") as out:
        out.write("viral_sequence_name\ttaxonomy\n")
        for seq in viral_names:
            tax = results.get(seq, "NA")
            out.write(f"{seq}\t{tax}\n")

    print(f"✅ Extracted taxonomy for {len(found)} of {len(viral_names)} sequences.")
    print(f"📄 Output written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract cluster reps stats from GFF file."
    )
    parser.add_argument(
        "--viral-list",
        required=True,
        help="Path to file with viral sequence names (first column)."
    )
    parser.add_argument(
        "--gff",
        required=True,
        help="Path to GFF file annotated by VIRify or similar tool."
    )
    parser.add_argument(
        "--mapfile",
        required=False,
        help="Map-file as product of renaming contigs step"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file with viral sequence name and taxonomy."
    )

    args = parser.parse_args()

    # Check file existence
    for f in [args.viral_list, args.gff]:
        if not os.path.exists(f):
            parser.error(f"File not found: {f}")

    extract_viral_taxonomy(args.viral_list, args.gff, args.output, args.mapfile)


if __name__ == "__main__":
    main()
