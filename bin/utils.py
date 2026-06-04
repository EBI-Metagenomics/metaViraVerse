#!/usr/bin/env python3
from __future__ import annotations


def parse_attributes(attrs_str: str) -> tuple[dict[str, str], list[str]]:
    """Parse a GFF3 column-9 attributes string into a dict and an ordered key list.

    :param attrs_str: Semicolon-separated key=value attribute string from GFF column 9.
    :return: Tuple of (attrs dict, list of keys in original order).
    """
    attrs, order = {}, []
    for part in attrs_str.rstrip(";").split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            if k not in attrs:
                order.append(k)
            attrs[k] = v
    return attrs, order


def read_input_gff(gffs: list[str]) -> tuple[dict[str, list[str]], dict[str, str], dict[str, str]]:
    """Parse one or more GFF3 files and index records by sequence ID.

    Builds three data structures:
    - ``gff_data``: maps each sequence ID (GFF column 1) to the list of raw
      GFF lines that belong to it.
    - ``attr_id_to_seq_id``: maps the ``ID`` attribute of non-CDS features to
      their parent sequence ID, for downstream attribute renaming.
    - ``source_map``: maps each sequence ID to the GFF source field (column 2)
      of its first non-CDS feature (e.g. ``VIRify``, ``geNomad``).

    Args:
        gffs: Paths to input GFF3 files.

    Returns:
        Tuple of (gff_data, attr_id_to_seq_id, source_map).
    """
    gff_data: dict[str, list[str]] = {}
    attr_id_to_seq_id: dict[str, str] = {}
    source_map: dict[str, str] = {}
    for gff in gffs:
        with open(gff, 'r') as file_in:
            for line in file_in:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                seq_id = parts[0]
                attrs, _ = parse_attributes(parts[8])
                gff_data.setdefault(seq_id, [])
                gff_data[seq_id].append(line)
                if parts[2] != 'CDS':
                    if attrs.get('ID'):
                        attr_id = attrs['ID']
                        if attr_id in attr_id_to_seq_id:
                            print(f"{attr_id} already exists in mapping {attr_id_to_seq_id[attr_id]}")
                        if attr_id not in source_map:
                            source_map[attr_id] = parts[1]
                        attr_id_to_seq_id.setdefault(attr_id, seq_id)
                    else:
                        print(f'There is no ID found for {line}')
    return gff_data, attr_id_to_seq_id, source_map