#!/usr/bin/env python3
"""
Generate a README.md for a published viral/plasmid catalogue.

The README contains:
  * a summary: the numbers from the pipeline's catalogue.json, the number of
    studies and the biomes the sequences came from, the metaViraVerse version,
    and the tools and annotations that were run (with versions when known),
  * a short description of how the catalogue was built,
  * links to the pipeline docs: databases.md (tool and database versions),
    methods.md, citations.md, output.md (all files) and faq.md (questions about results).

Which tools are listed follows the run's options (--skip_*, --viral_clustering,
--predict_host_from_custom_spacers, ...), read from pipeline_info/params_*.json.

The pipeline runs it at the end of every run (GENERATE_README, published as
README.md in --outdir). It can also be run on a finished results directory,
e.g. to change the title or links at release time:

    python bin/generate_catalogue_readme.py \\
        --results-dir results/ \\
        --catalogue-name "Marine sediment" --catalogue-version v1.0 \\
        --source-catalogue https://www.ebi.ac.uk/metagenomics/genome-catalogues/marine-sediment-v1-0 \\
        -o README.md

With --results-dir the script finds, when present:
  * catalogue.json
  * pipeline_info/*software_mqc_versions.yml  (pipeline and tool versions)
  * pipeline_info/params_*.json               (run options; latest file)
  * viruses/viruses-all-metadata.tsv.gz, plasmids/plasmids-all-metadata.tsv.gz
    (number of studies and biomes, and the "Metadata statistics" section:
    counts by Source_of_prediction, Source and checkv_quality, completeness and
    contamination ranges)
Each of these can also be given explicitly (--catalogue-json, --versions-yml,
--params, --metadata). Without a params file the pipeline defaults are assumed.
The pipeline version is taken from, in order: --pipeline-version, the versions
YAML, the manifest version in nextflow.config (with a warning).

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
import statistics
import sys

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_NEXTFLOW_CONFIG = os.path.join(REPO_DIR, "nextflow.config")
PIPELINE_NAME = "EBI-Metagenomics/metaviraverse"
PIPELINE_URL = "https://github.com/EBI-Metagenomics/metaViraVerse"

MISSING_VALUES = {"", "missing", "na", "not-provided", "none"}

# Markers of a missing value as written in the metadata tables (see docs/output.md#metadata);
# "" is reported as "empty".
MISSING_MARKERS = ["missing", "NA", "not-provided", ""]

# Columns summarised as counts per value, and as min/max/median/mean, in "Metadata statistics".
# Columns absent from a table (e.g. CheckV columns for plasmids) are skipped.
CATEGORY_COLUMNS = ["Source_of_prediction", "Source", "checkv_quality"]
NUMERIC_COLUMNS = ["completeness", "contamination"]

# Defaults of the params that decide which tools run (main.nf params block).
DEFAULT_PARAMS = {
    "third_party_input": None,
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
    "mmseqs_id_threshold": 0.9,
    "mmseqs_cov_threshold": 0.8,
}

# Display name -> (keys in the software versions YAML, URL)
TOOLS = {
    "CheckV": (["checkv"], "https://bitbucket.org/berkeleylab/checkv"),
    "barrnap": (["barrnap"], "https://github.com/tseemann/barrnap"),
    "vclust": (["vclust"], "https://github.com/refresh-bio/vclust"),
    "vjdb-clustering": ([], "https://github.com/maltesie/vjdb-clustering"),
    "BLASTN + anicalc/aniclust": (["blast"], "https://bitbucket.org/berkeleylab/checkv"),
    "geNomad": (["genomad"], "https://github.com/apcamargo/genomad"),
    "VITAP": (["vitap"], "https://github.com/DrKaiyangZheng/VITAP"),
    "ViPhOGs": ([], "https://github.com/EBI-Metagenomics/emg-viral-pipeline"),
    "iPHoP": (["iphop"], "https://bitbucket.org/srouxjgi/iphop"),
    "SpacePHARER": (["spacepharer"], "https://github.com/soedinglab/spacepharer"),
    "BACPHLIP": (["bacphlip"], "https://github.com/adamhockenberry/bacphlip"),
    "HMMER": (["hmmer"], "http://hmmer.org/"),
    "AMRFinderPlus": (["amrfinder", "amrfinderplus"], "https://github.com/ncbi/amr"),
    "DeepARG": (["deeparg"], "https://github.com/gaarangoa/deeparg"),
    "RGI": (["rgi"], "https://github.com/arpcard/rgi"),
    "PhaMMseqs": (["phammseqs"], "https://github.com/chg60/phammseqs"),
    "MMseqs2": (["mmseqs", "mmseqs2"], "https://github.com/soedinglab/MMseqs2"),
    "MOB-suite": (["mobsuite", "mob_suite", "mob_typer"], "https://github.com/phac-nml/mob-suite"),
    "plaSquid": (["plasquid"], "https://github.com/mgimenez720/plaSquid"),
    "Pyrodigal": (["pyrodigal"], "https://github.com/althonos/pyrodigal"),
}


# --------------------------------------------------------------------------- #
# Input helpers
# --------------------------------------------------------------------------- #
def open_text(path):
    """Open a plain or gzip-compressed text file (detected by magic bytes)."""
    with open(path, "rb") as f:
        is_gzip = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gzip else open(path)


def first_existing(*patterns, last=False):
    """Return the first (or last) path matching any of the glob patterns, or None."""
    for pattern in patterns:
        hits = sorted(glob.glob(pattern))
        if hits:
            return hits[-1] if last else hits[0]
    return None


def read_versions_yml(path):
    """Return {key (lower case): version} for every 'key: version' line of the versions YAML.

    The file is written by nf-core's softwareVersionsToYAML, e.g.
        METAVIRAVERSE:CHECKV:
            checkv: 1.0.3
        Workflow:
            EBI-Metagenomics/metaviraverse: v1.0.0-gabc1234
    It is parsed line by line so that PyYAML is not needed. When a tool appears
    with several versions they are joined with ', '.
    """
    versions = {}
    with open_text(path) as f:
        for line in f:
            match = re.match(r"^\s+['\"]?([^'\":]+?)['\"]?\s*:\s*['\"]?([^'\"]+?)['\"]?\s*$", line)
            if match:
                key, version = match.group(1).strip().lower(), match.group(2).strip()
                versions.setdefault(key, [])
                if version not in versions[key]:
                    versions[key].append(version)
    return {k: ", ".join(v) for k, v in versions.items()}


def version_from_nextflow_config(path):
    """Read the manifest version from nextflow.config."""
    if not os.path.exists(path):
        return None
    with open(path) as f:
        text = f.read()
    manifest = re.search(r"manifest\s*\{(.*?)\n\}", text, re.DOTALL)
    if not manifest:
        return None
    match = re.search(r"^\s*version\s*=\s*['\"]([^'\"]+)['\"]", manifest.group(1), re.MULTILINE)
    return match.group(1) if match else None


def normalise_version(version):
    return version if version.startswith("v") else f"v{version}"


def release_tag(version):
    """Release tag of a version, or None for development versions.

    nf-core appends the commit to the version of a run from git ('v1.0.0-gabc1234');
    that is still the v1.0.0 release.
    """
    match = re.fullmatch(r"(v\d+(?:\.\d+)*)(?:-g[0-9a-f]{7,40})?", version)
    return match.group(1) if match else None


def release_url(version):
    """Release URL for a version; None for development versions."""
    tag = release_tag(version)
    return f"{PIPELINE_URL}/releases/tag/{tag}" if tag else None


def docs_url(version, page="output.md"):
    """Link to a page of docs/ in the pipeline repo: at the release tag, or on dev for other versions."""
    ref = release_tag(version) or "dev"
    return f"{PIPELINE_URL}/blob/{ref}/docs/{page}"


def load_params(path):
    params = dict(DEFAULT_PARAMS)
    if path:
        with open(path) as f:
            run_params = json.load(f)
        for key in params:
            if run_params.get(key) is not None:
                params[key] = run_params[key]
    return params


def collect_metadata_values(paths, columns):
    """Collect sorted unique values of the given columns from metadata TSVs.

    Values may be comma- or semicolon-separated lists; missing markers are skipped.
    """
    values = {col: set() for col in columns}
    for path in paths:
        with open_text(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                for col in columns:
                    for value in re.split(r"[,;]", row.get(col) or ""):
                        value = value.strip()
                        if value.lower() not in MISSING_VALUES and value != col:
                            values[col].add(value)
    return {col: sorted(vals) for col, vals in values.items()}


def missing_marker(value):
    """Return the missing marker of a value ('missing', 'NA', 'not-provided', ''), or None."""
    value = (value or "").strip()
    for marker in MISSING_MARKERS:
        if value.lower() == marker.lower():
            return marker
    return None


def table_label(path):
    """viruses/viruses-all-metadata.tsv.gz -> 'viruses'."""
    name = os.path.basename(path)
    return name.split("-all-metadata")[0] if "-all-metadata" in name else name.split(".")[0]


def metadata_table_stats(path):
    """Read one metadata table and summarise it.

    Returns a dict with:
      rows       number of records
      categories {column: {value: count}} for CATEGORY_COLUMNS present in the table
                 (missing values counted under their marker)
      numeric    {column: {"n", "min", "max", "median", "mean"}} for NUMERIC_COLUMNS present
                 (values that are missing or not numbers are left out; None when no value)
    """
    rows = 0
    categories, numbers = {}, {}
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns = reader.fieldnames or []
        categories = {c: {} for c in CATEGORY_COLUMNS if c in columns}
        numbers = {c: [] for c in NUMERIC_COLUMNS if c in columns}
        for row in reader:
            if row.get(columns[0]) == columns[0]:  # repeated header
                continue
            rows += 1
            for col in list(categories) + list(numbers):
                value = (row.get(col) or "").strip()
                marker = missing_marker(value)
                if col in categories:
                    key = value if marker is None else (marker or "empty")
                    categories[col][key] = categories[col].get(key, 0) + 1
                if col in numbers and marker is None:
                    try:
                        numbers[col].append(float(value))
                    except ValueError:
                        pass
    numeric = {}
    for col, values in numbers.items():
        numeric[col] = {
            "n": len(values),
            "min": min(values) if values else None,
            "max": max(values) if values else None,
            "median": statistics.median(values) if values else None,
            "mean": statistics.fmean(values) if values else None,
        }
    return {"rows": rows, "categories": categories, "numeric": numeric}


def collect_metadata_stats(paths):
    """{table label: metadata_table_stats(path)} for each metadata table."""
    return {table_label(path): metadata_table_stats(path) for path in paths}


def pct_of(count, total):
    return f"{100 * count / total:.1f}%" if total else "-"


def render_metadata_stats(metadata_stats):
    """Markdown lines of the "Metadata statistics" section."""
    lines = []
    for label, st in metadata_stats.items():
        total = st["rows"]
        lines += [f"### {label.capitalize()}", "", f"{total:,} records.", ""]
        for col, counts in st["categories"].items():
            if not counts:
                continue
            lines += [f"`{col}`:", "", "| Value | Records | % |", "| --- | ---: | ---: |"]
            for value, count in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
                lines.append(f"| {value} | {count:,} | {pct_of(count, total)} |")
            lines.append("")
        numeric = {c: v for c, v in st["numeric"].items() if v["n"]}
        if numeric:
            lines += ["| Column | With a value | Min | Max | Median | Mean |",
                      "| --- | ---: | ---: | ---: | ---: | ---: |"]
            for col, v in numeric.items():
                lines.append(
                    f"| `{col}` | {v['n']:,} | {v['min']:.2f} | {v['max']:.2f} "
                    f"| {v['median']:.2f} | {v['mean']:.2f} |"
                )
            lines.append("")
    return lines


# --------------------------------------------------------------------------- #
# What the pipeline ran
# --------------------------------------------------------------------------- #
def pct(fraction):
    """0.9 -> '90%' (values above 1 are taken as percentages already)."""
    value = float(fraction)
    value = value * 100 if value <= 1 else value
    return f"{value:g}%"


def amr_tools(p, kind):
    return [name for name, flag in (("AMRFinderPlus", "amrfinderplus"), ("DeepARG", "deeparg"), ("RGI", "rgi"))
            if not p[f"skip_{flag}_for_{kind}"]]


def viral_clustering(p):
    """(tool, description) of the viral clustering method."""
    if p["viral_clustering"] == "vjdb":
        return "vjdb-clustering", ("MMseqs2 pre-clustering (95% identity) followed by vclust "
                                   "(Leiden, 95% ANI, 85% coverage)")
    if p["cluster_vclust"]:
        return "vclust", "95% ANI and 85% coverage"
    return "BLASTN + anicalc/aniclust", "95% ANI and 85% coverage"


def plasmid_clustering(p):
    if p["cluster_vclust"]:
        return "vclust", "pairs with at least 70% ANI and 50% coverage, Leiden clustering at gANI 0.35"
    return "BLASTN + anicalc/aniclust", "ANI-based greedy clustering"


def annotation_steps(p):
    """Ordered list of (step, [tools], details) that were run, following the run options."""
    steps = []
    qc = [t for t, skip in (("CheckV", p["skip_checkv"]), ("barrnap", p["skip_rrna_detection"])) if not skip]
    if qc:
        steps.append(("Quality control of viral sequences and prophages", qc,
                      "completeness, contamination and rRNA detection"))
    tool, details = viral_clustering(p)
    steps.append(("Viral clustering", [tool], details))
    taxonomy = [t for t, flag in (("geNomad", "skip_genomad"), ("VITAP", "skip_vitap"),
                                  ("ViPhOGs", "skip_viphogs")) if not p[flag]]
    if taxonomy:
        steps.append(("Viral taxonomy", taxonomy, "each tool assigns its own lineage"))
    host = (["iPHoP"] if not p["skip_iphop"] else []) + (
        ["SpacePHARER"] if p["predict_host_from_custom_spacers"] else [])
    if host:
        steps.append(("Viral host prediction", host,
                      "genome- and genus-level host" + (" and CRISPR spacer matches" if "SpacePHARER" in host else "")))
    steps.append(("Viral lifestyle", ["BACPHLIP"], "virulent or temperate"))
    steps.append(("Viral protein function", ["HMMER"], "searches against curated HMM profile databases"))
    if amr_tools(p, "viruses"):
        steps.append(("Antimicrobial resistance genes on viruses", amr_tools(p, "viruses"), ""))
    if not p["skip_phammseqs"]:
        steps.append(("Viral protein clustering", ["PhaMMseqs"], "protein families, phams"))
    tool, details = plasmid_clustering(p)
    steps.append(("Plasmid clustering", [tool], details))
    steps.append(("Plasmid typing and mobility", ["MOB-suite", "plaSquid"],
                  "replicon, relaxase, mate-pair formation, incompatibility group and mobility"))
    if amr_tools(p, "plasmids"):
        steps.append(("Antimicrobial resistance genes on plasmids", amr_tools(p, "plasmids"), ""))
    if not p["skip_mmseqs"]:
        steps.append(("Plasmid protein clustering", ["MMseqs2"],
                      f"linclust at {pct(p['mmseqs_id_threshold'])} identity and "
                      f"{pct(p['mmseqs_cov_threshold'])} coverage"))
    if p["third_party_input"]:
        steps.insert(0, ("Gene calling for third-party sequences", ["Pyrodigal"], ""))
    return steps


def tool_label(name, tool_versions):
    keys, url = TOOLS.get(name, ([], None))
    version = next((tool_versions[k] for k in keys if k in tool_versions), None)
    label = f"[{name}]({url})" if url else name
    if version:
        label += f" {version if version.startswith('v') or ',' in version else 'v' + version}"
    return label


# --------------------------------------------------------------------------- #
# README rendering
# --------------------------------------------------------------------------- #
def fmt(value):
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:,.2f}"
    return str(value)


def num(stats, key, default="an unknown number of"):
    value = stats.get(key)
    return fmt(value) if isinstance(value, (int, float)) else default


def summary_bullets(stats):
    """Plain-language summary sentences, skipping any whose numbers are absent."""
    def has(*keys):
        return all(isinstance(stats.get(k), (int, float)) for k in keys)

    bullets = []
    if has("total_sequences", "duplicate_sequences", "deduplicated_sequences", "qc_excluded_sequences",
           "unique_sequences"):
        qc = []
        if has("rrna_excluded_sequences") and stats["rrna_excluded_sequences"]:
            qc.append(f"{fmt(stats['rrna_excluded_sequences'])} carrying rRNA genes")
        if has("not_determined_excluded_sequences") and stats["not_determined_excluded_sequences"]:
            qc.append(f"{fmt(stats['not_determined_excluded_sequences'])} with an undetermined CheckV quality "
                      "and no viral genes")
        bullets.append(
            f"A total of {fmt(stats['total_sequences'])} viral sequences, prophages and plasmids were collected. "
            f"{fmt(stats['duplicate_sequences'])} identical copies were merged, giving "
            f"{fmt(stats['deduplicated_sequences'])} unique sequences; "
            f"{fmt(stats['qc_excluded_sequences'])} of them were removed at the quality-control step"
            + (f" ({join_and(qc)})" if qc else "")
            + f", leaving {fmt(stats['unique_sequences'])} sequences in the catalogue."
        )
    elif has("total_sequences", "unique_sequences", "qc_excluded_sequences"):  # catalogue.json of older runs
        bullets.append(
            f"A total of {fmt(stats['total_sequences'])} viral sequences, prophages and plasmids were collected; "
            f"{fmt(stats['qc_excluded_sequences'])} were removed at the quality-control step, "
            f"leaving {fmt(stats['unique_sequences'])} unique sequences."
        )
    if has("viral_sequences", "prophages", "viral_clusters"):
        bullets.append(
            f"{fmt(stats['viral_sequences'])} viral sequences and {fmt(stats['prophages'])} prophages "
            f"were clustered into {fmt(stats['viral_clusters'])} viral clusters."
        )
    if has("plasmids", "plasmid_clusters"):
        bullets.append(
            f"{fmt(stats['plasmids'])} plasmids were clustered into {fmt(stats['plasmid_clusters'])} plasmid clusters."
        )
    if has("total_proteins_viruses", "total_proteins_viruses_reps"):
        bullets.append(
            f"{fmt(stats['total_proteins_viruses'])} proteins were predicted on viral sequences and prophages, "
            f"{fmt(stats['total_proteins_viruses_reps'])} of them on viral cluster representatives."
        )
    if has("total_proteins_plasmids", "total_proteins_plasmids_reps"):
        bullets.append(
            f"{fmt(stats['total_proteins_plasmids'])} proteins were predicted on plasmids, "
            f"{fmt(stats['total_proteins_plasmids_reps'])} of them on plasmid cluster representatives."
        )
    return bullets


def methods_paragraphs(p, stats, tool_versions, has_source_catalogues):
    """Prose description of what the pipeline did, following the run options."""
    t = lambda name: tool_label(name, tool_versions)  # noqa: E731
    sources = "the source MAG catalogue(s)" if has_source_catalogues else "MGnify MAG catalogues"
    third_party = (" Sequences from third-party sources were added to them, with genes predicted by "
                   f"{t('Pyrodigal')} where no proteins were provided." if p["third_party_input"] else "")
    paragraphs = []

    # Pre-processing
    qc_parts = []
    if not p["skip_checkv"]:
        qc_parts.append(f"assessed for completeness and contamination with {t('CheckV')}")
    if not p["skip_rrna_detection"]:
        qc_parts.append(f"screened for rRNA genes with {t('barrnap')}")
    qc_sentence = (f" Viral sequences and prophages were {' and '.join(qc_parts)}." if qc_parts else "")
    criteria = []
    if not p["skip_rrna_detection"]:
        criteria.append("carrying rRNA genes")
    if not p["skip_checkv"]:
        criteria.append("with an undetermined CheckV quality, no viral genes and a high k-mer frequency")
    exclusion = (f" Viral sequences and prophages {', or '.join(criteria)} were then excluded "
                 f"({num(stats, 'qc_excluded_sequences')} sequences in total)." if criteria else "")
    merged = (f" ({fmt(stats['duplicate_sequences'])} copies)"
              if isinstance(stats.get("duplicate_sequences"), int) else "")
    paragraphs.append(
        "**Collection and quality control.** Viral sequences, prophages and plasmids predicted in "
        f"{sources} (by [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) and the "
        "[mobilome-annotation-pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline)) "
        f"were collected.{third_party} Every sequence received a catalogue identifier.{qc_sentence} "
        f"Identical sequences were merged by SHA-256 checksum{merged}, keeping one record per unique sequence and "
        f"merging the biomes of all copies.{exclusion}"
    )

    # Viruses
    tool, details = viral_clustering(p)
    taxonomy = [n for n, flag in (("geNomad", "skip_genomad"), ("VITAP", "skip_vitap"),
                                  ("ViPhOGs", "skip_viphogs")) if not p[flag]]
    sentences = [
        f"**Viruses.** The {num(stats, 'viral_sequences')} viral sequences and {num(stats, 'prophages')} "
        f"prophages were pooled and clustered with {t(tool)} ({details}) into "
        f"{num(stats, 'viral_clusters')} clusters. Only cluster representatives were annotated."
    ]
    if taxonomy:
        sentences.append(
            f"Taxonomy was assigned independently by {join_and([t(n) for n in taxonomy])}; "
            "the lineages are reported side by side and not merged into a consensus."
        )
    if not p["skip_iphop"]:
        sentences.append(f"Hosts were predicted at genome and genus level with {t('iPHoP')}"
                         + (f", and from matches to a custom CRISPR spacer database with {t('SpacePHARER')}."
                            if p["predict_host_from_custom_spacers"] else "."))
    elif p["predict_host_from_custom_spacers"]:
        sentences.append(f"Hosts were predicted from matches to a custom CRISPR spacer database with "
                         f"{t('SpacePHARER')}.")
    sentences.append(f"Lifestyle (virulent or temperate) was predicted with {t('BACPHLIP')}, and protein "
                     f"functions were annotated with {t('HMMER')} against curated HMM profile databases.")
    if amr_tools(p, "viruses"):
        sentences.append("Antimicrobial resistance genes were detected with "
                         f"{join_and([t(n) for n in amr_tools(p, 'viruses')])}.")
    if not p["skip_phammseqs"]:
        sentences.append(f"Proteins of the representatives were grouped into families (phams) with {t('PhaMMseqs')}.")
    sentences.append("These annotations were merged into one GFF file of the representatives.")
    paragraphs.append(" ".join(sentences))

    # Plasmids
    tool, details = plasmid_clustering(p)
    sentences = [
        f"**Plasmids.** The {num(stats, 'plasmids')} plasmids were clustered with {t(tool)} ({details}) into "
        f"{num(stats, 'plasmid_clusters')} clusters. The cluster representatives were typed with "
        f"{t('MOB-suite')} (`mob_typer`: replicon, relaxase, mate-pair formation and predicted mobility) and "
        f"searched with {t('plaSquid')} for replication initiator proteins, relaxases and incompatibility groups."
    ]
    if amr_tools(p, "plasmids"):
        sentences.append("Antimicrobial resistance genes were detected with "
                         f"{join_and([t(n) for n in amr_tools(p, 'plasmids')])}.")
    if not p["skip_mmseqs"]:
        sentences.append(f"Their proteins were clustered with {t('MMseqs2')} `linclust` at "
                         f"{pct(p['mmseqs_id_threshold'])} identity and {pct(p['mmseqs_cov_threshold'])} coverage.")
    sentences.append("The genome each plasmid was found in is treated as its host.")
    paragraphs.append(" ".join(sentences))

    # Metadata
    paragraphs.append(
        "**Metadata.** Every sequence was linked back to its original name and to the source genome it came "
        "from, with that genome's taxonomic lineage, sample and study. These are combined with the quality, "
        "clustering, taxonomy and host results into one metadata table for viruses, one for plasmids, and a "
        "table of per-cluster statistics for viruses."
    )
    return paragraphs


def join_and(items):
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def render_readme(args, version, stats, metadata_values, params, tool_versions, metadata_stats=None):
    title_parts = [args.catalogue_name or "MGnify viral and plasmid catalogue"]
    if args.catalogue_version:
        title_parts.append(args.catalogue_version)
    lines = [f"# {' '.join(title_parts)}", ""]

    if args.catalogue_version:
        lines += [f"{args.catalogue_version} release", ""]
    if args.website_url:
        lines += [f"Website URL: {args.website_url}", ""]
    if args.release_date:
        lines += [f"Release date: {args.release_date}", ""]

    # Summary
    lines += ["## Summary", ""]
    bullets = summary_bullets(stats)
    if args.source_catalogue:
        bullets.append("Sequences were collected from the following catalogue(s): " + ", ".join(args.source_catalogue))
    studies = metadata_values.get("Source_project") or []
    if studies:
        bullets.append(f"Sequences come from genomes of {len(studies):,} {'study' if len(studies) == 1 else 'studies'}.")
    biomes = metadata_values.get("Biome") or []
    if biomes:
        bullets.append(f"Sequences come from {len(biomes)} biome(s): " + ", ".join(biomes))
    elif isinstance(stats.get("number_of_biomes"), int):
        bullets.append(f"Sequences come from {fmt(stats['number_of_biomes'])} biome(s).")
    url = release_url(version)
    bullets.append(
        f"The catalogue was generated using the metaViraVerse pipeline {version}"
        + (f": {url}" if url else f" ({PIPELINE_URL}).")
    )
    lines += [f"* {b}" for b in bullets]

    lines += ["", "### Tools and annotations", ""]
    for step, tools, details in annotation_steps(params):
        tools_text = ", ".join(tool_label(n, tool_versions) for n in tools)
        lines.append(f"* {step}: {tools_text}" + (f" ({details})" if details else ""))
    lines += [
        "",
        "Versions of the tools and of the reference databases used for the annotation are listed in "
        f"{docs_url(version, 'databases.md')}",
    ]


    # Methods
    lines += ["", "## How the catalogue was built", ""]
    for paragraph in methods_paragraphs(params, stats, tool_versions, bool(args.source_catalogue)):
        lines += [paragraph, ""]
    lines += [
        f"For more details on each step, see {docs_url(version, 'methods.md')}",
        "",
        f"If you use this catalogue, please also cite the tools listed in {docs_url(version, 'citations.md')}",
        "",
    ]

    # Metadata statistics
    if metadata_stats:
        lines += [
            "## Metadata statistics",
            "",
            "Counts are over all records of the metadata tables. `missing`, `NA` and `not-provided` mark "
            f"values that are not available; what they mean for each column is described in "
            f"{docs_url(version, 'output.md')}#metadata",
            "",
        ]
        lines += render_metadata_stats(metadata_stats)

    # Files
    lines += [
        "## Files",
        "",
        "All files of the catalogue, and the columns of the metadata tables, are described in the pipeline "
        f"documentation: {docs_url(version)}",
        "",
        "## Questions",
        "",
        f"If you have questions about the results, check {docs_url(version, 'faq.md')}",
        "",
        "---",
        "",
        f"Generated on {datetime.date.today().isoformat()} by metaViraVerse {version} ({PIPELINE_URL}).",
    ]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate a README.md for a viral/plasmid catalogue release.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-r", "--results-dir",
                        help="Pipeline results directory; used to find the inputs below.")
    parser.add_argument("-j", "--catalogue-json",
                        help="catalogue.json produced by the pipeline (default: <results-dir>/catalogue.json).")
    parser.add_argument("--versions-yml",
                        help="pipeline_info/*software_mqc_versions.yml of the run (pipeline and tool versions).")
    parser.add_argument("-p", "--params",
                        help="params JSON of the run (default: latest <results-dir>/pipeline_info/params_*.json).")
    parser.add_argument("--pipeline-version",
                        help="Pipeline version to report (overrides --versions-yml).")
    parser.add_argument("--metadata", nargs="*", default=None,
                        help="*-all-metadata.tsv(.gz) tables to count studies, list biomes and summarise "
                             "in the Metadata statistics section "
                             "(default: found in --results-dir; pass with no value to disable).")
    parser.add_argument("--catalogue-name",
                        help="Catalogue name for the title, e.g. 'Marine sediment viral and plasmid catalogue'.")
    parser.add_argument("--catalogue-version", help="Catalogue release, e.g. v1.0.")
    parser.add_argument("--website-url", help="Web page of the catalogue.")
    parser.add_argument("--release-date", help="Release date to print (free text).")
    parser.add_argument("--source-catalogue", nargs="+", default=[],
                        help="URL(s) of the source MAG catalogue(s); separate with spaces or commas.")
    parser.add_argument("-o", "--output", default="README.md", help="Output file (default: README.md).")
    return parser.parse_args(argv)


def resolve_inputs(args):
    rd = args.results_dir
    if rd and not os.path.isdir(rd):
        sys.exit(f"Error: results directory not found: {rd}")

    if not args.catalogue_json and rd:
        args.catalogue_json = first_existing(os.path.join(rd, "catalogue.json"))
    if not args.catalogue_json:
        sys.exit("Error: catalogue.json not found; pass --catalogue-json or --results-dir.")

    if not args.versions_yml and rd:
        args.versions_yml = first_existing(os.path.join(rd, "pipeline_info", "*software_mqc_versions.yml"))

    if not args.params and rd:
        args.params = first_existing(os.path.join(rd, "pipeline_info", "params_*.json"), last=True)
    if not args.params:
        print("Warning: no params JSON found; pipeline defaults assumed for the list of tools.", file=sys.stderr)

    if args.metadata is None:
        args.metadata = []
        if rd:
            for pattern in ("viruses/viruses-all-metadata.tsv*", "plasmids/plasmids-all-metadata.tsv*"):
                hit = first_existing(os.path.join(rd, pattern))
                if hit:
                    args.metadata.append(hit)


def resolve_version(args, tool_versions):
    if args.pipeline_version:
        return normalise_version(args.pipeline_version)
    version = tool_versions.get(PIPELINE_NAME.lower())
    if version:
        return normalise_version(version)
    if args.versions_yml:
        print(f"Warning: no '{PIPELINE_NAME}' entry in {args.versions_yml}.", file=sys.stderr)
    version = version_from_nextflow_config(DEFAULT_NEXTFLOW_CONFIG)
    if version:
        print(
            f"Warning: pipeline version taken from nextflow.config ({version}); it may differ from the "
            "version used for the run. Pass --pipeline-version or --versions-yml to be exact.",
            file=sys.stderr,
        )
        return normalise_version(version)
    sys.exit("Error: could not determine the pipeline version; pass --pipeline-version.")


def main(argv=None):
    args = parse_args(argv)
    args.source_catalogue = [u.strip() for value in args.source_catalogue for u in value.split(",") if u.strip()]
    resolve_inputs(args)
    tool_versions = read_versions_yml(args.versions_yml) if args.versions_yml else {}
    version = resolve_version(args, tool_versions)
    params = load_params(args.params)

    with open_text(args.catalogue_json) as f:
        stats = json.load(f)

    metadata_values = (
        collect_metadata_values(args.metadata, ["Source_project", "Biome"]) if args.metadata else {}
    )
    metadata_stats = collect_metadata_stats(args.metadata) if args.metadata else {}
    readme = render_readme(args, version, stats, metadata_values, params, tool_versions, metadata_stats)
    with open(args.output, "w") as f:
        f.write(readme)
    print(f"README written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
