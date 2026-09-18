#!/usr/bin/env python3
"""
Add plasmid-evidence attributes to a representative GFF.

Three independent evidence sources:

--table (plaSquid's protein_report.tsv: Contig, Protein, RIP_domain,
MOB_group, Inc_group -- one row per protein with at least one plaSquid hit,
see plasquid_rip_extraction.R): matched by `ID=` against each CDS record,
and its non-missing columns added as new column-9 attributes. "NA" cells
(R's write_delim representation of a missing value) mean "no evidence of
that kind" and are skipped, not written out literally.

--biomarker-report (MOB-suite's plasmids_biomarker_report.txt, from
`mob_typer --biomarker_report_file`): a BLAST-based table of replicon/
relaxase/mate-pair-formation/oriT biomarkers found anywhere on the plasmid
sequence, not attributed to any single protein. Its `sseqid` column holds
the plasmid/contig identifier from the same representative assembly used as
mob_typer's input -- i.e. it matches this GFF's seqid (column 1), not a
CDS `ID=` attribute -- so its `biomarker` (category, e.g. "relaxase") and
`qseqid` (the specific reference biomarker sequence matched, e.g. "IncN_1")
columns are added as comma-joined `mob_suite_biomarker` and
`mobsuite_identifier` attributes to that contig's sequence-level (non-CDS)
record(s), the same record level BacPhlip/taxonomy attributes use elsewhere
in this pipeline (see build_final_gff.py) for whole-contig facts.

--amr-gff (AMRINTEGRATOR's integrated_<prefix>.gff): built by augmenting
--gff itself with AMRfinderPlus/RGI/DeepARG hits, so it shares --gff's `ID=`
values and only contains the plasmids/proteins that actually got an AMR
hit -- a subset, not the full representative set. Rather than hardcoding
which of its attribute keys are "AMR fields", this script computes them:
any column-9 key that appears anywhere in --amr-gff but nowhere in --gff is,
by construction, something AMRINTEGRATOR added -- i.e. AMR-specific -- and
gets carried over to the matching --gff record (any feature type, matched
by `ID=`); keys already present in --gff (ID, product, locus_tag, etc.) are
never re-added, so only genuinely new, AMR-related fields are picked up.

In all cases, existing attributes on a record are never overwritten, and
records with no matching evidence pass through completely unchanged.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict

from utils import parse_attributes

PROTEIN_REPORT_COLUMNS = ("RIP_domain", "MOB_group", "Inc_group")
# biomarker report column -> GFF attribute name
BIOMARKER_REPORT_COLUMNS = {
    "biomarker": "mob_suite_biomarker",
    "qseqid": "mobsuite_identifier",
}


def open_file(path: str):
    """Open a plain or gzip-compressed file for reading as text."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def attrs_to_str(attrs: dict[str, str], order: list[str]) -> str:
    """Serialise attrs back to a GFF column-9 string (original key order first)."""
    parts = []
    seen = set()
    for k in order:
        if k in attrs:
            parts.append(f"{k}={attrs[k]}")
            seen.add(k)
    for k, v in attrs.items():
        if k not in seen:
            parts.append(f"{k}={v}")
    return ";".join(parts)


def load_protein_report(path: str) -> dict[str, dict[str, str]]:
    """Return protein_id -> {RIP_domain: ..., MOB_group: ..., Inc_group: ...}.

    Only non-"NA", non-empty cells are kept, so a protein with e.g. just a
    MOB_group hit contributes a single-key dict, not RIP_domain=NA etc.
    """
    result: dict[str, dict[str, str]] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            protein_id = (row.get("Protein") or "").strip()
            if not protein_id:
                continue
            attrs = {}
            for col in PROTEIN_REPORT_COLUMNS:
                value = (row.get(col) or "").strip()
                if value and value != "NA":
                    attrs[col] = value
            if attrs:
                result[protein_id] = attrs
    return result


def load_biomarker_report(path: str) -> dict[str, dict[str, str]]:
    """Return contig_id -> {mob_suite_biomarker: ..., mobsuite_identifier: ...}.

    ``biomarker`` is one of MOB-suite's fixed categories (replicon, relaxase,
    mate-pair-formation, oriT); ``qseqid`` is the specific reference biomarker
    sequence matched (e.g. "IncN_1"). A contig can have more than one row
    (e.g. both a relaxase and a replicon hit, or repeat hits of the same
    category/reference), so both columns are deduplicated and comma-joined
    into one attribute value each.
    """
    values: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            contig_id = (row.get("sseqid") or "").strip()
            if not contig_id:
                continue
            for column, attr in BIOMARKER_REPORT_COLUMNS.items():
                value = (row.get(column) or "").strip()
                if value:
                    values[contig_id][attr].add(value)
    return {
        contig_id: {attr: ",".join(sorted(vals)) for attr, vals in attrs.items()}
        for contig_id, attrs in values.items()
    }


def load_gff_attribute_keys(path: str) -> set[str]:
    """Return the set of every column-9 attribute key used anywhere in a GFF."""
    keys: set[str] = set()
    with open_file(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs, _ = parse_attributes(cols[8])
            keys.update(attrs.keys())
    return keys


def load_amr_gff(path: str, base_attr_keys: set[str]) -> dict[str, dict[str, str]]:
    """Return record_id -> {new_attr_key: value, ...} for AMRINTEGRATOR's output.

    ``record_id`` is a record's `ID=` attribute -- AMRINTEGRATOR augments the
    same GFF passed in via --gff rather than renumbering it, so its `ID=`
    values match directly, unlike the plaSquid/MOB-suite sources above.

    Only attribute keys absent from ``base_attr_keys`` (i.e. --gff's own
    schema) are kept per record: these are exactly the fields AMRINTEGRATOR
    itself introduced (AMR gene name, tool, coordinates, etc.), so nothing
    AMR-unrelated (ID, product, locus_tag, ...) gets carried over even though
    this GFF repeats every field on every record.
    """
    result: dict[str, dict[str, str]] = {}
    with open_file(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs, _ = parse_attributes(cols[8])
            record_id = attrs.get("ID", "").strip()
            if not record_id:
                continue
            new_attrs = {k: v for k, v in attrs.items() if k not in base_attr_keys}
            if new_attrs:
                result[record_id] = new_attrs
    return result


def apply_evidence(attrs: dict[str, str], evidence: dict[str, str] | None) -> bool:
    """Add ``evidence``'s keys to ``attrs`` in place, skipping any already present.

    Returns True if anything was added.
    """
    if not evidence:
        return False
    modified = False
    for k, v in evidence.items():
        if k not in attrs:
            attrs[k] = v
            modified = True
    return modified


def annotate_gff(
    input_gff: str,
    protein_report: dict[str, dict[str, str]],
    biomarker_report: dict[str, dict[str, str]],
    amr_report: dict[str, dict[str, str]],
    out,
) -> tuple[int, int, int]:
    """Write ``input_gff`` to ``out`` with all three evidence sources added.

    Returns (cds_records_annotated, contig_records_annotated, amr_records_annotated).
    """
    cds_annotated = 0
    contig_annotated = 0
    amr_annotated = 0

    with open_file(input_gff) as f:
        for line in f:
            raw = line.rstrip("\n")

            if raw.startswith("#") or not raw.strip():
                out.write(raw + "\n")
                continue

            cols = raw.split("\t")
            if len(cols) < 9:
                out.write(raw + "\n")
                continue

            attrs, order = parse_attributes(cols[8])
            record_id = attrs.get("ID", "").strip()
            modified = False

            if cols[2] == "CDS":
                if apply_evidence(attrs, protein_report.get(record_id)):
                    modified = True
                    cds_annotated += 1
            else:
                if apply_evidence(attrs, biomarker_report.get(cols[0])):
                    modified = True
                    contig_annotated += 1

            # AMR evidence applies to any feature type -- keyed by ID=, which
            # AMRINTEGRATOR's output shares directly with this GFF.
            if apply_evidence(attrs, amr_report.get(record_id)):
                modified = True
                amr_annotated += 1

            if modified:
                cols[8] = attrs_to_str(attrs, order)

            out.write("\t".join(cols) + "\n")

    return cds_annotated, contig_annotated, amr_annotated


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add plaSquid, MOB-suite and AMR plasmid-evidence attributes to a representative GFF.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-g", "--gff", required=True,
                        help="Representative GFF file (plain or .gz).")
    parser.add_argument("-t", "--table", required=True,
                        help="plaSquid protein_report.tsv (Contig, Protein, RIP_domain, MOB_group, Inc_group).")
    parser.add_argument("-b", "--biomarker-report", default=None,
                        help="MOB-suite plasmids_biomarker_report.txt (optional).")
    parser.add_argument("-a", "--amr-gff", default=None,
                        help="AMRINTEGRATOR's integrated_<prefix>.gff (optional).")
    parser.add_argument("-o", "--output", default="plasquid_annotated.gff",
                        help="Output GFF file (default: plasquid_annotated.gff).")
    return parser.parse_args()


def main():
    args = parse_args()

    protein_report = load_protein_report(args.table)
    print(f"Proteins with plaSquid evidence loaded: {len(protein_report)}", file=sys.stderr)

    biomarker_report = load_biomarker_report(args.biomarker_report) if args.biomarker_report else {}
    print(f"Contigs with MOB-suite biomarker evidence loaded: {len(biomarker_report)}", file=sys.stderr)

    if args.amr_gff:
        base_attr_keys = load_gff_attribute_keys(args.gff)
        amr_report = load_amr_gff(args.amr_gff, base_attr_keys)
    else:
        amr_report = {}
    print(f"Records with AMR evidence loaded: {len(amr_report)}", file=sys.stderr)

    with open(args.output, "w") as out:
        cds_annotated, contig_annotated, amr_annotated = annotate_gff(
            args.gff, protein_report, biomarker_report, amr_report, out
        )

    print(f"CDS records annotated (plaSquid): {cds_annotated}", file=sys.stderr)
    print(f"Contig records annotated (MOB-suite): {contig_annotated}", file=sys.stderr)
    print(f"Records annotated (AMR): {amr_annotated}", file=sys.stderr)
    print(f"Output written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
