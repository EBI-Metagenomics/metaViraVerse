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


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def pct(n, total):
    return f"{100 * n / total:.1f}%" if total else "-"


def render_report(results_dir, params, params_path, file_rows, metadata_reports, stats, count_checks):
    all_statuses = (
        [r["status"] for r in file_rows]
        + [c["status"] for m in metadata_reports if m.get("present") for c in m["columns"]]
        + [FAIL for m in metadata_reports if not m.get("present")]
        + [c["status"] for c in count_checks]
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
    report, overall = render_report(args.results_dir, params, params_path, file_rows,
                                    metadata_reports, stats, count_checks)

    with open(args.output, "w") as f:
        f.write(report)
    print(f"Sanity check: {overall}. Report written to {args.output}", file=sys.stderr)
    return 1 if overall == FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
