#!/usr/bin/env python3
"""
Build a final GFF by merging annotations from multiple sources into a base GFF.

Sources
-------
-i / --input            Base GFF (e.g. viruses_reps.gff)
--bacphlip FILE [FILE]  BacPhlip TSV output (one or more)
                          Adds bacphlip_virulent and bacphlip_temperate to
                          non-CDS (sequence-level) records.
--hmmer FILE [FILE]     HMMER --tblout files, plain or .gz (one or more)
                          Adds {prefix}_accession and {prefix}_description to
                          CDS records, where prefix = filename without .tbl/.gz,
                          lowercased.  First (best-scoring) hit per protein is used.
--gff FILE [FILE]       Additional GFF files (one or more)
                          Adds new attributes to CDS records only; existing
                          attributes in the base GFF are never overwritten.

All input files may be plain text or gzip-compressed.
"""

import argparse
import gzip
import hashlib
import os
import sys
from collections import defaultdict


# ── helpers ───────────────────────────────────────────────────────────────────

def open_file(path):
    """Open a plain or gzip-compressed file for reading as text."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_attrs(attrs_str):
    """Return (dict, ordered-key-list) from a GFF column-9 string."""
    attrs, order = {}, []
    for part in attrs_str.rstrip(";").split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            if k not in attrs:
                order.append(k)
            attrs[k] = v
    return attrs, order


def attrs_to_str(attrs, order):
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


# ── loaders ───────────────────────────────────────────────────────────────────

def load_bacphlip(paths):
    """Return seq_id → {bacphlip_virulent: …, bacphlip_temperate: …}."""
    result = {}
    for path in paths:
        with open_file(path) as f:
            for line in f:
                if "Virulent" in line and "Temperate" in line:
                    continue  # skip header
                row = line.rstrip("\n").split("\t")
                if len(row) < 3:
                    continue
                seq_id = row[0].strip()
                if seq_id:
                    result[seq_id] = {
                        "bacphlip_virulent":  round(float(row[1].strip()), 2),
                        "bacphlip_temperate": round(float(row[2].strip()), 2),
                    }
    return result


def _hmmer_prefix(path):
    """PVOG.tbl.gz → pvog,  myhits.tbl → myhits."""
    base = os.path.basename(path)
    for ext in (".gz", ".tbl"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base.lower()


def load_hmmer(paths):
    """Return protein_id → {f'{prefix}_accession': …, f'{prefix}_description': …}.

    Only the first (best-scoring) hit per protein per file is kept.
    Multiple files contribute independent prefixed keys that are all merged.
    """
    result = defaultdict(dict)
    for path in paths:
        prefix = _hmmer_prefix(path)
        seen = set()
        with open_file(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.split()
                if len(cols) < 3:
                    continue
                target = cols[0]
                if target in seen:
                    continue
                seen.add(target)
                result[target][f"{prefix}_accession"]    = cols[2]
                result[target][f"{prefix}_description"] = " ".join(cols[18:]) if len(cols) > 18 else ""
    return result


def load_extra_gffs(paths):
    """Return protein_id → {attr_key: attr_val, …} from CDS lines of extra GFFs."""
    result = defaultdict(dict)
    for path in paths:
        with open_file(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9 or cols[2] != "CDS":
                    continue
                attrs, _ = parse_attrs(cols[8])
                prot_id = attrs.get("ID", "").strip()
                if prot_id:
                    for k, v in attrs.items():
                        if k != "ID":
                            result[prot_id][k] = v
    return result


# ── main ──────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge annotations from multiple sources into a base GFF.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Base GFF file (plain or .gz).")
    parser.add_argument("-o", "--output", default="final.gff",
                        help="Output GFF file (default: final.gff).")
    parser.add_argument("--bacphlip", nargs="+", default=[],
                        help="BacPhlip TSV file(s).")
    parser.add_argument("--hmmer", nargs="+", default=[],
                        help="HMMER --tblout file(s), plain or .gz.")
    parser.add_argument("--gff", nargs="+", default=[],
                        help="Additional GFF file(s) to merge attributes from.")
    parser.add_argument("--compress-output", action="store_true",
                        help="Gzip-compress the output file and write an md5 checksum file.")
    return parser.parse_args()


def main():
    args = parse_args()

    # ── load annotation sources ───────────────────────────────────────────────
    bacphlip = load_bacphlip(args.bacphlip) if args.bacphlip else {}
    hmmer    = load_hmmer(args.hmmer)       if args.hmmer    else {}
    extra    = load_extra_gffs(args.gff)   if args.gff      else {}

    print(f"BacPhlip sequences loaded : {len(bacphlip)}", file=sys.stderr)
    print(f"HMMER protein hits loaded : {len(hmmer)}",    file=sys.stderr)
    print(f"Extra-GFF proteins loaded : {len(extra)}",    file=sys.stderr)

    # ── resolve output path (add .gz suffix if compressing) ──────────────────
    out_path = args.output
    if args.compress_output and not out_path.endswith(".gz"):
        out_path += ".gz"

    if out_path.endswith(".gz"):
        out = gzip.open(out_path, "wt")
    elif out_path == "-":
        out = sys.stdout
    else:
        out = open(out_path, "w")

    # ── process base GFF line by line ─────────────────────────────────────────
    bacphlip_annotated = hmmer_annotated = extra_annotated = 0

    try:
        with open_file(args.input) as f:
            for line in f:
                raw = line.rstrip("\n")

                # pass through comments and directives unchanged
                if raw.startswith("#") or not raw.strip():
                    out.write(raw + "\n")
                    continue

                cols = raw.split("\t")
                if len(cols) < 9:
                    out.write(raw + "\n")
                    continue

                feat_type = cols[2]
                attrs, order = parse_attrs(cols[8])
                record_id = attrs.get("ID", "").strip()
                modified = False

                if feat_type != "CDS":
                    # sequence-level record: add BacPhlip scores
                    bp = bacphlip.get(record_id, {})
                    for k, v in bp.items():
                        if k not in attrs:
                            attrs[k] = v
                            modified = True
                    if bp:
                        bacphlip_annotated += 1

                else:
                    # CDS: add HMMER hits
                    hm = hmmer.get(record_id, {})
                    existing_product = attrs.get("product", "").lower()
                    for k, v in hm.items():
                        if k not in attrs:
                            # skip description when it matches the existing product
                            if k.endswith("_description") and v.lower() == existing_product:
                                continue
                            attrs[k] = v
                            modified = True
                    if hm:
                        hmmer_annotated += 1

                    # CDS: add new attrs from extra GFFs (no overwrite)
                    ex = extra.get(record_id, {})
                    for k, v in ex.items():
                        if k not in attrs:
                            attrs[k] = v
                            modified = True
                    if ex:
                        extra_annotated += 1

                if modified:
                    cols[8] = attrs_to_str(attrs, order)

                out.write("\t".join(cols) + "\n")

    finally:
        if out is not sys.stdout:
            out.close()

    # ── optional md5 ─────────────────────────────────────────────────────────
    if args.compress_output and out_path != "-":
        md5 = hashlib.md5(open(out_path, "rb").read()).hexdigest()
        md5_path = out_path + ".md5"
        with open(md5_path, "w") as f:
            f.write(f"{md5}  {out_path}\n")
        print(f"md5 written to {md5_path}", file=sys.stderr)

    print(f"BacPhlip annotations added : {bacphlip_annotated}", file=sys.stderr)
    print(f"HMMER annotations added    : {hmmer_annotated}",    file=sys.stderr)
    print(f"Extra-GFF annotations added: {extra_annotated}",    file=sys.stderr)
    print(f"Output written to          : {out_path}",           file=sys.stderr)


if __name__ == "__main__":
    main()
