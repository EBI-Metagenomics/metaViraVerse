#!/usr/bin/env python3
"""
Sanity-check the results directory of a metaViraVerse run and write a Markdown report.

Checks
------
1. File presence: every file the pipeline should publish (see the `output` block
   in main.nf) is present and not empty. Files that depend on run options
   (--skip_*, --predict_host_from_custom_spacers, --viral_clustering,
   --save_intermediates, ...) are only expected when those options were enabled.
   Options are read from pipeline_info/params_*.json (or --params); without it,
   the pipeline defaults are assumed and missing optional files are warnings.
2. Missing values in the metadata tables: for every column, the number and
   share of values that are empty, `missing`, `NA` or `not-provided`. Columns
   that docs/output.md documents as never missing, but have missing values,
   are reported as failures; fully empty columns are reported as warnings.
3. Row counts of the metadata tables against catalogue.json.
4. Sequence identifiers and protein counts, for every set of sequence files
   (all viral sequences, prophages, plasmids, and the cluster representatives):
   * FASTA (nucleotide) headers and GFF column 1 are renamed catalogue
     accessions (<prefix><digits>, e.g. MGYV0000000001);
   * the ID= attribute of non-CDS GFF records (the sequence-level record) is a
     renamed accession too. CDS records keep their original protein IDs, so
     their ID= is not checked;
   * the number of CDS records in the GFF equals the number of proteins in the
     FAA, and their IDs match;
   * sequence and protein counts against catalogue.json where it has them.
   The prefix is --accession-prefix, else rename_accession from the params
   JSON, else MGYV.

Usage
-----
    python scripts/check_catalogue_results.py -r results/ -o sanity_report.md

Exit status is 1 when any check fails (useful in automation), 0 otherwise.
Only the Python standard library is used.
"""

import argparse
import csv
import datetime
import glob
import gzip
import json
import os
import re
import sys
from collections import Counter

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUTPUT_DOCS = os.path.join(REPO_DIR, "docs", "output.md")

MISSING_MARKERS = ("missing", "NA", "not-provided")  # plus the empty string
OK, WARN, FAIL = "OK", "WARN", "FAIL"
ICON = {OK: "✅", WARN: "⚠️", FAIL: "❌"}

# Defaults of the params that decide which files are published (main.nf params block).
DEFAULT_PARAMS = {
    "skip_rrna_detection": False,
    "skip_checkv": False,
    "skip_vitap": False,
    "skip_genomad": False,
    "skip_viphogs": False,
    "cluster_vclust": True,
    "viral_clustering": "default",
    "skip_amrfinderplus_for_viruses": True,
    "skip_deeparg_for_viruses": False,
    "skip_rgi_for_viruses": False,
    "skip_amrfinderplus_for_plasmids": True,
    "skip_deeparg_for_plasmids": False,
    "skip_rgi_for_plasmids": False,
    "skip_iphop": False,
    "predict_host_from_custom_spacers": False,
    "skip_phammseqs": False,
    "skip_mmseqs": False,
    "save_intermediates": False,
}

VR = "viruses/cluster_representatives"
PR = "plasmids/cluster_representatives"


def amr_any(p, kind):
    return any(not p[f"skip_{tool}_for_{kind}"] for tool in ("amrfinderplus", "deeparg", "rgi"))


# (group, glob pattern(s), condition on params or None for always expected).
# A tuple of patterns means "any of these"; the condition gets the params dict.
EXPECTED_FILES = [
    # Summary and metadata
    ("Summary", "catalogue.json", None),
    ("Summary", "README.md", None),
    ("Summary", "viruses/viruses-all-metadata.tsv.gz", None),
    ("Summary", "viruses/viruses-cluster-stats.tsv.gz", None),
    ("Summary", "plasmids/plasmids-all-metadata.tsv.gz", None),
    ("Summary", "viruses/viruses_source_sankey.html", None),
    ("Summary", "plasmids/plasmids_source_sankey.html", None),
    # Additional data
    ("Additional data", "additional_data/rename_mapfile.tsv", None),
    ("Additional data", "additional_data/catalogue_filtered_metadata.tsv", None),
    ("Additional data", "additional_data/catalogue_excluded_records_metadata.tsv", None),
    ("Additional data", "additional_data/viruses_checkv_report.tsv",
     lambda p: p["save_intermediates"] and not p["skip_checkv"]),
    ("Additional data", "additional_data/barrnap/*",
     lambda p: p["save_intermediates"] and not p["skip_rrna_detection"]),
    ("Additional data", "additional_data/rename/*", lambda p: p["save_intermediates"]),
    ("Additional data", "additional_data/choose_sequences/*", lambda p: p["save_intermediates"]),
    # All sequences
    *[("Viruses", f"viruses/{name}.{ext}.gz", None)
      for name in ("viral_sequences", "prophages") for ext in ("fna", "faa", "gff")],
    *[("Plasmids", f"plasmids/plasmids.{ext}.gz", None) for ext in ("fna", "faa", "gff")],
    # Viral clustering
    ("Viruses", "viruses/clustering/vjdb1_merged_reps.csv", lambda p: p["viral_clustering"] == "vjdb"),
    ("Viruses", "viruses/clustering/*", lambda p: not p["cluster_vclust"]),
    # Viral representatives
    ("Viruses", (f"{VR}/viruses_reps.fasta.gz", f"{VR}/viruses_reps.fna.gz"), None),
    ("Viruses", f"{VR}/viruses_reps.faa.gz", None),
    ("Viruses", f"{VR}/viruses_reps.gff.gz", None),
    ("Viruses", f"{VR}/viruses_reps_annotated.gff.gz", None),
    ("Viruses", f"{VR}/functional_annotation/lifestyle/*bacphlip*.gz", None),
    ("Viruses", f"{VR}/functional_annotation/hmmer/*.tbl.gz", None),
    ("Viruses", f"{VR}/functional_annotation/hmmer/*_summary.tsv.gz", None),
    ("Viruses", f"{VR}/functional_annotation/amr/viruses_amrfinderplus.tsv",
     lambda p: not p["skip_amrfinderplus_for_viruses"]),
    ("Viruses", f"{VR}/functional_annotation/amr/viruses_deeparg.tsv",
     lambda p: not p["skip_deeparg_for_viruses"]),
    ("Viruses", f"{VR}/functional_annotation/amr/viruses_rgi.tsv",
     lambda p: not p["skip_rgi_for_viruses"]),
    ("Viruses", f"{VR}/functional_annotation/amr/integrated_viruses.gff", lambda p: amr_any(p, "viruses")),
    ("Viruses", f"{VR}/taxonomy/genomad/viruses.tsv", lambda p: not p["skip_genomad"]),
    ("Viruses", f"{VR}/taxonomy/genomad/*_genomad_taxonomy_counts.tsv", lambda p: not p["skip_genomad"]),
    ("Viruses", f"{VR}/taxonomy/viphogs/*_viphogs_taxonomy_counts.tsv", lambda p: not p["skip_viphogs"]),
    ("Viruses", f"{VR}/taxonomy/vitap/*_vitap_taxonomy_counts.tsv", lambda p: not p["skip_vitap"]),
    ("Viruses", f"{VR}/host_detection/iphop/iphop_genome.csv", lambda p: not p["skip_iphop"]),
    ("Viruses", f"{VR}/host_detection/iphop/iphop_genus.csv", lambda p: not p["skip_iphop"]),
    ("Viruses", f"{VR}/host_detection/spacepharer/*spacers_predictions.tsv",
     lambda p: p["predict_host_from_custom_spacers"]),
    ("Viruses", f"{VR}/host_detection/spacepharer/*custom_chosen_host.tsv",
     lambda p: p["predict_host_from_custom_spacers"]),
    ("Viruses", f"{VR}/protein_clusters_phammseqs/viruses_phams.tsv", lambda p: not p["skip_phammseqs"]),
    ("Viruses", f"{VR}/protein_clusters_phammseqs/viruses_clustering_stats.tsv",
     lambda p: not p["skip_phammseqs"]),
    # Plasmid clustering and representatives
    ("Plasmids", "plasmids/clustering/*", lambda p: not p["cluster_vclust"]),
    ("Plasmids", (f"{PR}/plasmids_reps.fasta.gz", f"{PR}/plasmids_reps.fna.gz"), None),
    ("Plasmids", f"{PR}/plasmids_reps.faa.gz", None),
    ("Plasmids", f"{PR}/plasmids_reps.gff.gz", None),
    ("Plasmids", f"{PR}/plasmids_annotated.gff.gz", None),
    ("Plasmids", f"{PR}/functional_annotation/mobility_stats_final.json", None),
    ("Plasmids", f"{PR}/functional_annotation/plasquid/final_report_per_contig.tsv", None),
    ("Plasmids", f"{PR}/functional_annotation/plasquid/final_report_per_protein.tsv", None),
    ("Plasmids", f"{PR}/functional_annotation/plasquid/mobility_classification.tsv", None),
    ("Plasmids", f"{PR}/functional_annotation/plasquid/mobility_stats_base.json", None),
    ("Plasmids", f"{PR}/functional_annotation/plasquid/evidence/*", lambda p: p["save_intermediates"]),
    ("Plasmids", f"{PR}/functional_annotation/mobsuite-typer/plasmids_report.txt", None),
    ("Plasmids", f"{PR}/functional_annotation/mobsuite-typer/*biomarker_report.txt", None),
    ("Plasmids", f"{PR}/functional_annotation/mobsuite-typer/*mge_report.txt", None),
    ("Plasmids", f"{PR}/functional_annotation/amr/plasmids_amrfinderplus.tsv",
     lambda p: not p["skip_amrfinderplus_for_plasmids"]),
    ("Plasmids", f"{PR}/functional_annotation/amr/plasmids_deeparg.tsv",
     lambda p: not p["skip_deeparg_for_plasmids"]),
    ("Plasmids", f"{PR}/functional_annotation/amr/plasmids_rgi.tsv",
     lambda p: not p["skip_rgi_for_plasmids"]),
    ("Plasmids", f"{PR}/functional_annotation/amr/integrated_plasmids.gff", lambda p: amr_any(p, "plasmids")),
    ("Plasmids", f"{PR}/protein_clusters_mmseqs/protein_catalogue-*.tsv", lambda p: not p["skip_mmseqs"]),
    ("Plasmids", f"{PR}/protein_clusters_mmseqs/plasmids_clustering_stats.tsv", lambda p: not p["skip_mmseqs"]),
    # Run reports
    ("Pipeline info", "pipeline_info/*software_mqc_versions.yml", None),
    ("Pipeline info", "pipeline_info/params_*.json", None),
    ("Pipeline info", "pipeline_info/execution_report*.html", None),
    ("Pipeline info", "multiqc/multiqc_report.html", lambda p: p["save_intermediates"]),
]

# Metadata tables whose missing values are counted, with the docs/output.md table that documents them.
METADATA_TABLES = [
    ("viruses/viruses-all-metadata.tsv.gz", "viruses-all-metadata.tsv.gz"),
    ("viruses/viruses-cluster-stats.tsv.gz", "viruses-cluster-stats.tsv.gz"),
    ("plasmids/plasmids-all-metadata.tsv.gz", "plasmids-all-metadata.tsv.gz"),
    ("additional_data/catalogue_filtered_metadata.tsv", None),
    ("additional_data/catalogue_excluded_records_metadata.tsv", None),
]

# Expected row count of a metadata table, as a sum of catalogue.json keys.
ROW_COUNT_CHECKS = [
    ("viruses/viruses-all-metadata.tsv.gz", ("viral_sequences", "prophages")),
    ("viruses/viruses-cluster-stats.tsv.gz", ("viral_clusters",)),
    ("plasmids/plasmids-all-metadata.tsv.gz", ("plasmids",)),
    ("additional_data/catalogue_filtered_metadata.tsv", ("unique_sequences",)),
    ("additional_data/catalogue_excluded_records_metadata.tsv", ("qc_excluded_sequences",)),
]


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def open_text(path):
    """Open a plain or gzip-compressed text file (detected by magic bytes)."""
    with open(path, "rb") as f:
        is_gzip = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", newline="") if is_gzip else open(path, newline="")


def is_empty(path):
    """True for a zero-byte file or a gzip file without content."""
    if os.path.getsize(path) == 0:
        return True
    try:
        with open_text(path) as f:
            return f.read(1) == ""
    except (OSError, UnicodeDecodeError, EOFError):
        return False  # binary (e.g. .png, .tbi) or truncated gzip: size > 0 is enough


def human_size(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def load_params(results_dir, params_path):
    params = dict(DEFAULT_PARAMS)
    if not params_path:
        hits = sorted(glob.glob(os.path.join(results_dir, "pipeline_info", "params_*.json")))
        params_path = hits[-1] if hits else None  # latest run
    if params_path:
        with open(params_path) as f:
            run_params = json.load(f)
        for key in params:
            if key in run_params and run_params[key] is not None:
                params[key] = run_params[key]
    return params, params_path


def parse_documented_columns(docs_path):
    """Return {table name: {column: can_be_missing (bool)}} from the column tables in docs/output.md."""
    if not docs_path or not os.path.exists(docs_path):
        return {}
    with open(docs_path) as f:
        text = f.read()
    tables = {}
    current = None
    for line in text.splitlines():
        heading = re.search(r"Columns of <code>([^<]+)</code>", line)
        if heading:
            current = heading.group(1)
            tables[current] = {}
            continue
        if current and line.startswith("|"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            col = re.fullmatch(r"`([^`]+)`", cells[0]) if cells else None
            if col and len(cells) > 1:
                tables[current][col.group(1)] = not cells[1].lower().startswith("no")
        elif current and line.startswith("</details>"):
            current = None
    return tables


# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #
def check_files(results_dir, params, params_known):
    """Return a list of dicts: group, pattern, status, found, note."""
    rows = []
    for group, patterns, condition in EXPECTED_FILES:
        patterns = (patterns,) if isinstance(patterns, str) else patterns
        expected = condition is None or bool(condition(params))
        hits = sorted({h for p in patterns for h in glob.glob(os.path.join(results_dir, p)) if os.path.isfile(h)})
        label = " or ".join(patterns)
        if hits:
            empty = [h for h in hits if is_empty(h)]
            size = sum(os.path.getsize(h) for h in hits)
            note = f"{len(hits)} file(s), {human_size(size)}"
            if empty:
                status = FAIL
                note += "; empty: " + ", ".join(os.path.relpath(h, results_dir) for h in empty)
            elif not expected:
                status, note = OK, note + " (not expected with these options)"
            else:
                status = OK
            rows.append(dict(group=group, pattern=label, status=status, note=note))
        elif expected:
            # Optional outputs are only warnings when the run options are unknown.
            status = FAIL if (condition is None or params_known) else WARN
            rows.append(dict(group=group, pattern=label, status=status, note="missing"))
        else:
            rows.append(dict(group=group, pattern=label, status="SKIP", note="not expected with these options"))
    return rows


def count_missing(path):
    """Return (n_rows, header, {column: Counter(marker -> count)})."""
    with open_text(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return 0, [], {}
        counts = {col: Counter() for col in header}
        n_rows = 0
        for row in reader:
            if row == header:  # repeated header from concatenated chunks
                continue
            n_rows += 1
            for i, col in enumerate(header):
                value = row[i].strip() if i < len(row) else ""
                if value == "":
                    counts[col]["empty"] += 1
                elif value in MISSING_MARKERS:
                    counts[col][value] += 1
    return n_rows, header, counts


def check_metadata(results_dir, documented):
    reports = []
    for rel_path, doc_name in METADATA_TABLES:
        path = os.path.join(results_dir, rel_path)
        if not os.path.isfile(path):
            reports.append(dict(path=rel_path, present=False))
            continue
        n_rows, header, counts = count_missing(path)
        doc_cols = documented.get(doc_name, {}) if doc_name else {}
        columns = []
        for col in header:
            c = counts[col]
            total = sum(c.values())
            if doc_cols and col in doc_cols and not doc_cols[col] and total:
                status, note = FAIL, "documented as never missing"
            elif n_rows and total == n_rows:
                status, note = WARN, "all values missing"
            else:
                status, note = OK, ""
            columns.append(dict(column=col, missing=total, breakdown=c, status=status, note=note))
        undocumented = [c for c in header if doc_cols and c not in doc_cols]
        not_in_file = [c for c in doc_cols if c not in header]
        reports.append(dict(path=rel_path, present=True, rows=n_rows, columns=columns,
                            undocumented=undocumented, not_in_file=not_in_file))
    return reports


def check_row_counts(results_dir, metadata_reports):
    path = os.path.join(results_dir, "catalogue.json")
    if not os.path.isfile(path):
        return None, []
    with open(path) as f:
        stats = json.load(f)
    rows_by_path = {r["path"]: r["rows"] for r in metadata_reports if r.get("present")}
    checks = []
    for rel_path, keys in ROW_COUNT_CHECKS:
        if rel_path not in rows_by_path:
            continue
        if not all(isinstance(stats.get(k), int) for k in keys):
            checks.append(dict(path=rel_path, keys=keys, expected=None, actual=rows_by_path[rel_path],
                               status=WARN, note="key(s) not in catalogue.json"))
            continue
        expected = sum(stats[k] for k in keys)
        actual = rows_by_path[rel_path]
        checks.append(dict(path=rel_path, keys=keys, expected=expected, actual=actual,
                           status=OK if expected == actual else FAIL, note=""))
    return stats, checks


# Sets of sequence files to check identifiers in: (name, fna patterns, gff, faa,
# catalogue.json keys summed for the sequence count, ... for the protein count).
# The annotated GFFs of the representatives only get the column 1 check.
SEQUENCE_SETS = [
    ("Viral sequences", ("viruses/viral_sequences.fna.gz",), "viruses/viral_sequences.gff.gz",
     "viruses/viral_sequences.faa.gz", ("viral_sequences",), ()),
    ("Prophages", ("viruses/prophages.fna.gz",), "viruses/prophages.gff.gz",
     "viruses/prophages.faa.gz", ("prophages",), ()),
    ("Plasmids", ("plasmids/plasmids.fna.gz",), "plasmids/plasmids.gff.gz",
     "plasmids/plasmids.faa.gz", ("plasmids",), ("total_proteins_plasmids",)),
    ("Viral cluster representatives", (f"{VR}/viruses_reps.fasta.gz", f"{VR}/viruses_reps.fna.gz"),
     f"{VR}/viruses_reps.gff.gz", f"{VR}/viruses_reps.faa.gz", ("viral_clusters",),
     ("total_proteins_viruses_reps",)),
    ("Plasmid cluster representatives", (f"{PR}/plasmids_reps.fasta.gz", f"{PR}/plasmids_reps.fna.gz"),
     f"{PR}/plasmids_reps.gff.gz", f"{PR}/plasmids_reps.faa.gz", ("plasmid_clusters",),
     ("total_proteins_plasmids_reps",)),
]
ANNOTATED_GFFS = [f"{VR}/viruses_reps_annotated.gff.gz", f"{PR}/plasmids_annotated.gff.gz"]
MAX_EXAMPLES = 5


def accession_prefix(cli_prefix, params_path):
    if cli_prefix:
        return cli_prefix
    if params_path:
        with open(params_path) as f:
            prefix = json.load(f).get("rename_accession")
        if prefix:
            return prefix
    return "MGYV"


class BadIds:
    """Count identifiers that are not renamed accessions, keeping a few examples."""

    def __init__(self):
        self.count = 0
        self.examples = []

    def add(self, value):
        self.count += 1
        if len(self.examples) < MAX_EXAMPLES:
            self.examples.append(value)


def fasta_ids(path):
    """Yield the first word of every FASTA header."""
    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].split(None, 1)[0] if line[1:].strip() else ""


def gff_attribute(attributes, key):
    for field in attributes.split(";"):
        if field.startswith(f"{key}="):
            return field[len(key) + 1:]
    return None


def scan_gff(path, accession_re, check_ids=True):
    """Return seqids, bad column-1 values, bad non-CDS IDs and the CDS IDs of a GFF."""
    seqids, bad_seqid, bad_id, cds_ids = set(), BadIds(), BadIds(), []
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, ftype, attributes = parts[0], parts[2], parts[8]
            if seqid not in seqids:
                seqids.add(seqid)
                if not accession_re.fullmatch(seqid):
                    bad_seqid.add(seqid)
            if not check_ids:
                continue
            feature_id = gff_attribute(attributes, "ID")
            if ftype == "CDS":
                cds_ids.append(feature_id)
            elif feature_id is not None and not accession_re.fullmatch(feature_id):
                bad_id.add(feature_id)
    return seqids, bad_seqid, bad_id, cds_ids


def check_sequence_ids(results_dir, prefix, stats):
    accession_re = re.compile(re.escape(prefix) + r"\d+")
    results = []
    for name, fna_patterns, gff_rel, faa_rel, seq_keys, prot_keys in SEQUENCE_SETS:
        fna = next((p for p in fna_patterns if os.path.isfile(os.path.join(results_dir, p))), None)
        gff = gff_rel if os.path.isfile(os.path.join(results_dir, gff_rel)) else None
        faa = faa_rel if os.path.isfile(os.path.join(results_dir, faa_rel)) else None
        r = dict(name=name, fna=fna or " or ".join(fna_patterns), gff=gff_rel, faa=faa_rel,
                 problems=[], status=OK, n_seqs=None, n_gff_seqs=None, n_cds=None, n_prot=None,
                 exp_seqs=None, exp_prot=None)

        if fna:
            bad_fna, n = BadIds(), 0
            fna_ids = set()
            for seq_id in fasta_ids(os.path.join(results_dir, fna)):
                n += 1
                fna_ids.add(seq_id)
                if not accession_re.fullmatch(seq_id):
                    bad_fna.add(seq_id)
            r["n_seqs"] = n
            if bad_fna.count:
                r["problems"].append((FAIL, f"{bad_fna.count:,} FASTA header(s) are not `{prefix}` accessions",
                                      bad_fna.examples))
        else:
            fna_ids = None
            r["problems"].append((WARN, "nucleotide FASTA not found", []))

        cds_ids = None
        if gff:
            seqids, bad_seqid, bad_id, cds_ids = scan_gff(os.path.join(results_dir, gff), accession_re)
            r["n_gff_seqs"], r["n_cds"] = len(seqids), len(cds_ids)
            if bad_seqid.count:
                r["problems"].append((FAIL, f"{bad_seqid.count:,} GFF sequence name(s) (column 1) are not "
                                            f"`{prefix}` accessions", bad_seqid.examples))
            if bad_id.count:
                r["problems"].append((FAIL, f"{bad_id.count:,} non-CDS GFF record ID(s) are not "
                                            f"`{prefix}` accessions", bad_id.examples))
            if fna_ids is not None:
                not_in_fna = sorted(seqids - fna_ids)
                if not_in_fna:
                    r["problems"].append((FAIL, f"{len(not_in_fna):,} GFF sequence(s) not in the FASTA",
                                          not_in_fna[:MAX_EXAMPLES]))
        else:
            r["problems"].append((WARN, "GFF not found", []))

        if faa:
            prot_ids = list(fasta_ids(os.path.join(results_dir, faa)))
            r["n_prot"] = len(prot_ids)
            if len(set(prot_ids)) != len(prot_ids):
                r["problems"].append((WARN, f"{len(prot_ids) - len(set(prot_ids)):,} duplicated protein ID(s) "
                                            "in the FAA", []))
            if cds_ids is not None:
                if len(cds_ids) != len(prot_ids):
                    r["problems"].append((FAIL, f"{len(cds_ids):,} CDS records in the GFF but "
                                                f"{len(prot_ids):,} proteins in the FAA", []))
                cds_set, prot_set = set(cds_ids), set(prot_ids)
                only_faa, only_gff = sorted(prot_set - cds_set), sorted(cds_set - prot_set - {None})
                if only_faa:
                    r["problems"].append((FAIL, f"{len(only_faa):,} FAA protein(s) without a CDS in the GFF",
                                          only_faa[:MAX_EXAMPLES]))
                if only_gff:
                    r["problems"].append((FAIL, f"{len(only_gff):,} GFF CDS ID(s) without a protein in the FAA",
                                          only_gff[:MAX_EXAMPLES]))
        else:
            r["problems"].append((WARN, "FAA not found", []))

        # Counts against catalogue.json
        if stats:
            if seq_keys and all(isinstance(stats.get(k), int) for k in seq_keys):
                r["exp_seqs"] = sum(stats[k] for k in seq_keys)
                if r["n_seqs"] is not None and r["n_seqs"] != r["exp_seqs"]:
                    r["problems"].append((FAIL, f"{r['n_seqs']:,} sequences in the FASTA, catalogue.json has "
                                                f"{r['exp_seqs']:,} ({' + '.join(seq_keys)})", []))
            if prot_keys and all(isinstance(stats.get(k), int) for k in prot_keys):
                r["exp_prot"] = sum(stats[k] for k in prot_keys)
                if r["n_prot"] is not None and r["n_prot"] != r["exp_prot"]:
                    r["problems"].append((FAIL, f"{r['n_prot']:,} proteins in the FAA, catalogue.json has "
                                                f"{r['exp_prot']:,} ({' + '.join(prot_keys)})", []))

        statuses = [st for st, _, _ in r["problems"]]
        r["status"] = FAIL if FAIL in statuses else WARN if WARN in statuses else OK
        results.append(r)

    for gff_rel in ANNOTATED_GFFS:
        path = os.path.join(results_dir, gff_rel)
        r = dict(name="Annotated GFF", gff=gff_rel, problems=[], annotated=True)
        if not os.path.isfile(path):
            r["problems"].append((WARN, "not found", []))
        else:
            seqids, bad_seqid, _, _ = scan_gff(path, accession_re, check_ids=False)
            r["n_gff_seqs"] = len(seqids)
            if bad_seqid.count:
                r["problems"].append((FAIL, f"{bad_seqid.count:,} sequence name(s) (column 1) are not "
                                            f"`{prefix}` accessions", bad_seqid.examples))
        statuses = [st for st, _, _ in r["problems"]]
        r["status"] = FAIL if FAIL in statuses else WARN if WARN in statuses else OK
        results.append(r)
    return results


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def pct(n, total):
    return f"{100 * n / total:.1f}%" if total else "-"


def render_report(results_dir, params, params_path, file_rows, metadata_reports, stats, count_checks,
                  id_checks=(), prefix="MGYV"):
    all_statuses = (
        [r["status"] for r in file_rows]
        + [c["status"] for m in metadata_reports if m.get("present") for c in m["columns"]]
        + [FAIL for m in metadata_reports if not m.get("present")]
        + [c["status"] for c in count_checks]
        + [r["status"] for r in id_checks]
    )
    n_fail, n_warn = all_statuses.count(FAIL), all_statuses.count(WARN)
    overall = FAIL if n_fail else WARN if n_warn else OK

    out = [
        "# Catalogue sanity check report",
        "",
        f"- **Results directory:** `{os.path.abspath(results_dir)}`",
        f"- **Generated:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"- **Run options:** "
        + (f"`{os.path.relpath(params_path, results_dir)}`" if params_path
           else "not found, pipeline defaults assumed (missing optional files are warnings)"),
        f"- **Overall:** {ICON[overall]} **{overall}** ({n_fail} failed, {n_warn} warnings)",
        "",
        "Status: ✅ OK · ⚠️ WARN (worth a look) · ❌ FAIL (needs fixing) · ➖ not expected with the run options",
        "",
    ]

    # Options that change the expected files
    changed = {k: v for k, v in params.items() if v != DEFAULT_PARAMS[k]}
    if changed:
        out += ["Non-default options affecting outputs: " + ", ".join(f"`{k}={v}`" for k, v in changed.items()), ""]

    # 1. Files
    found = sum(1 for r in file_rows if r["status"] == OK)
    expected = sum(1 for r in file_rows if r["status"] != "SKIP")
    out += ["## 1. File presence", "", f"{found} of {expected} expected outputs found.", ""]
    problems = [r for r in file_rows if r["status"] in (FAIL, WARN)]
    if problems:
        out += ["### Problems", "", "| Status | File | Note |", "| --- | --- | --- |"]
        out += [f"| {ICON[r['status']]} | `{r['pattern']}` | {r['note']} |" for r in problems]
        out.append("")
    out += ["<details>", "<summary>All checked files</summary>", ""]
    for group in dict.fromkeys(r["group"] for r in file_rows):
        out += [f"**{group}**", "", "| Status | File | Note |", "| --- | --- | --- |"]
        for r in file_rows:
            if r["group"] == group:
                icon = ICON.get(r["status"], "➖")
                out.append(f"| {icon} | `{r['pattern']}` | {r['note']} |")
        out.append("")
    out += ["</details>", ""]

    # 2. Missing values
    out += [
        "## 2. Missing values in metadata tables",
        "",
        "Counted as missing: empty cells and " + ", ".join(f"`{m}`" for m in MISSING_MARKERS) + ".",
        "",
    ]
    for m in metadata_reports:
        out += [f"### `{m['path']}`", ""]
        if not m["present"]:
            out += [f"{ICON[FAIL]} File not found.", ""]
            continue
        n = m["rows"]
        with_missing = [c for c in m["columns"] if c["missing"]]
        out += [f"{n:,} rows, {len(m['columns'])} columns, {len(with_missing)} column(s) with missing values.", ""]
        if not n:
            out += [f"{ICON[WARN]} Table has no data rows.", ""]
            continue
        out += ["| Status | Column | Missing | % | " + " | ".join(f"`{x}`" for x in MISSING_MARKERS) + " | empty | Note |",
                "| --- | --- | ---: | ---: | " + " | ".join("---:" for _ in MISSING_MARKERS) + " | ---: | --- |"]
        for c in m["columns"]:
            b = c["breakdown"]
            cells = " | ".join(f"{b.get(x, 0):,}" for x in (*MISSING_MARKERS, "empty"))
            out.append(f"| {ICON[c['status']]} | `{c['column']}` | {c['missing']:,} | {pct(c['missing'], n)} "
                       f"| {cells} | {c['note']} |")
        out.append("")
        if m["undocumented"]:
            out += [f"{ICON[WARN]} Columns not described in docs/output.md: "
                    + ", ".join(f"`{c}`" for c in m["undocumented"]), ""]
        if m["not_in_file"]:
            out += [f"{ICON[WARN]} Documented columns absent from the file: "
                    + ", ".join(f"`{c}`" for c in m["not_in_file"]), ""]

    # 3. Row counts
    out += ["## 3. Row counts against `catalogue.json`", ""]
    if stats is None:
        out += [f"{ICON[FAIL]} `catalogue.json` not found.", ""]
    elif not count_checks:
        out += ["No metadata tables to compare.", ""]
    else:
        out += ["| Status | Table | catalogue.json | Expected | Rows |", "| --- | --- | --- | ---: | ---: |"]
        for c in count_checks:
            exp = f"{c['expected']:,}" if c["expected"] is not None else c["note"]
            out.append(f"| {ICON[c['status']]} | `{c['path']}` | {' + '.join(f'`{k}`' for k in c['keys'])} "
                       f"| {exp} | {c['actual']:,} |")
        out.append("")

    # 4. Identifiers and proteins
    def n(value):
        return f"{value:,}" if isinstance(value, int) else "-"

    out += [
        "## 4. Sequence identifiers and protein counts",
        "",
        f"FASTA headers, GFF column 1 and the `ID=` of non-CDS GFF records must be `{prefix}<digits>` accessions. "
        "CDS records keep their original protein IDs, which must match the FAA headers.",
        "",
        "| Status | Set | Sequences (FASTA) | Expected | Sequences (GFF) | CDS (GFF) | Proteins (FAA) | Expected |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in id_checks:
        if r.get("annotated"):
            out.append(f"| {ICON[r['status']]} | `{r['gff']}` (column 1 only) | - | - | "
                       f"{n(r.get('n_gff_seqs'))} | - | - | - |")
        else:
            out.append(f"| {ICON[r['status']]} | {r['name']} | {n(r['n_seqs'])} | {n(r['exp_seqs'])} | "
                       f"{n(r['n_gff_seqs'])} | {n(r['n_cds'])} | {n(r['n_prot'])} | {n(r['exp_prot'])} |")
    out.append("")
    problems = [r for r in id_checks if r["problems"]]
    if problems:
        out += ["### Problems", ""]
        for r in problems:
            label = r["gff"] if r.get("annotated") else r["name"]
            for status, message, examples in r["problems"]:
                line = f"- {ICON[status]} **{label}**: {message}"
                if examples:
                    line += " (e.g. " + ", ".join(f"`{e}`" for e in examples) + ")"
                out.append(line)
        out.append("")

    return "\n".join(out).rstrip() + "\n", overall


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Sanity-check a metaViraVerse results directory and write a Markdown report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-r", "--results-dir", required=True, help="Pipeline results directory (--outdir).")
    parser.add_argument("-p", "--params",
                        help="params JSON of the run (default: latest <results-dir>/pipeline_info/params_*.json).")
    parser.add_argument("--output-docs", default=DEFAULT_OUTPUT_DOCS,
                        help="docs/output.md, for the documented columns (default: this repo's).")
    parser.add_argument("--accession-prefix",
                        help="Prefix of renamed accessions (default: rename_accession from the params JSON, "
                             "else MGYV).")
    parser.add_argument("-o", "--output", default="sanity_report.md", help="Output Markdown report.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if not os.path.isdir(args.results_dir):
        sys.exit(f"Error: results directory not found: {args.results_dir}")

    params, params_path = load_params(args.results_dir, args.params)
    file_rows = check_files(args.results_dir, params, params_known=params_path is not None)
    metadata_reports = check_metadata(args.results_dir, parse_documented_columns(args.output_docs))
    stats, count_checks = check_row_counts(args.results_dir, metadata_reports)
    prefix = accession_prefix(args.accession_prefix, params_path)
    id_checks = check_sequence_ids(args.results_dir, prefix, stats)
    report, overall = render_report(args.results_dir, params, params_path, file_rows,
                                    metadata_reports, stats, count_checks, id_checks, prefix)

    with open(args.output, "w") as f:
        f.write(report)
    print(f"Sanity check: {overall}. Report written to {args.output}", file=sys.stderr)
    return 1 if overall == FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
