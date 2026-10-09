#!/usr/bin/env python3
"""
Plots and statistics for protein clustering from PhaMMseqs or MMseqs2.

Clustering input (exactly one):
  --phammseqs  phams_to_tsv.py output, one line per protein:
                   <phamID>\t<proteinID> <function>
               The function is everything after the first whitespace of the
               FASTA header (a `product=` value is used if present).
  --mmseqs     MMseqs2 cluster TSV (`mmseqs createtsv`), one line per protein:
                   <cluster_representative>\t<member>
               The representative's ID is used as the cluster ID.

Functions (optional, --fasta): protein FASTA whose headers carry the function,
e.g. `>MGYG000934530_00007 Protein PucC`. Functions from the FASTA take
precedence over the ones in a PhaMMseqs TSV. With --mmseqs and no --fasta,
all proteins are unannotated and only size plots/statistics are produced.

Outputs (in --outdir, file names start with --prefix):
  _cluster_size_distribution.tsv/.png  Clusters and proteins per cluster size; binned bar plot
  _cluster_rank_abundance.png          Cluster size vs rank (log-log)
  _top_functions.tsv                   All informative functions, ranked by number of proteins
  _top_functions.png                   Top --top functions (default 20)
  _cluster_function_purity.tsv/.png    Per cluster: share of annotated proteins with its main function
  _clustering_stats.tsv                Summary statistics and clustering agreement metrics

"Informative" functions exclude empty descriptions and --exclude values
(default: "hypothetical protein", "putative protein"; case-insensitive,
whole-function match).

Clustering agreement metrics compare the clusters with a reference labelling
of the same proteins: by default the function annotation (informative
proteins only), or a --reference TSV (proteinID<TAB>label). They need
scikit-learn; without it they are skipped. Clusters are the "predicted"
labels and the reference the "true" labels, so:
  homogeneity  ~ each cluster contains one function
  completeness ~ each function sits in one cluster (often low: diverged
                 homologs with the same function form separate clusters)
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BAR_COLOR = "#B0692A"
INK = "#0F201C"
INK_2 = "#47625D"
GRID = "#E3ECE9"

SIZE_BINS = [
    ("1", 1, 1),
    ("2", 2, 2),
    ("3–5", 3, 5),
    ("6–10", 6, 10),
    ("11–50", 11, 50),
    ("51–100", 51, 100),
    (">100", 101, math.inf),
]

DEFAULT_EXCLUDE = ["hypothetical protein", "putative protein"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    clustering = p.add_mutually_exclusive_group(required=True)
    clustering.add_argument("--phammseqs", type=Path,
                            help="PhaMMseqs phams TSV (phamID<TAB>proteinID function); plain or .gz")
    clustering.add_argument("--mmseqs", type=Path,
                            help="MMseqs2 cluster TSV (representative<TAB>member); plain or .gz")
    p.add_argument("--fasta", type=Path, default=None,
                   help="Protein FASTA with the function in the header (>proteinID function); plain or .gz")
    p.add_argument("-o", "--outdir", default=Path("."), type=Path, help="Output directory (default: .)")
    p.add_argument("--prefix", default="clusters", help="Output file prefix (default: clusters)")
    p.add_argument("--top", type=int, default=20, help="Number of functions to plot (default: 20)")
    p.add_argument("--exclude", nargs="+", default=DEFAULT_EXCLUDE,
                   help="Functions treated as uninformative (case-insensitive, whole match). "
                        f"Default: {' / '.join(DEFAULT_EXCLUDE)}")
    p.add_argument("--reference", type=Path, default=None,
                   help="Optional TSV (proteinID<TAB>label) to compare clusters against instead of functions")
    return p.parse_args()


# ── Parsing ───────────────────────────────────────────────────────────────────

def open_text(path: Path):
    """Open a plain or gzip-compressed text file for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def extract_function(description: str) -> str:
    """Return the function from a FASTA description ('product=' value if present)."""
    m = re.search(r"(?:^|;)\s*product=([^;]*)", description)
    function = m.group(1) if m else description
    return " ".join(function.split())


def load_phammseqs(path: Path) -> tuple[list[str], list[str], list[str]]:
    """Return parallel lists (cluster, protein_id, function) from a PhaMMseqs TSV."""
    clusters, proteins, functions = [], [], []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            pham, _, rest = line.partition("\t")
            parts = rest.strip().split(maxsplit=1)
            if not parts:
                continue
            clusters.append(pham.strip())
            proteins.append(parts[0])
            functions.append(extract_function(parts[1]) if len(parts) > 1 else "")
    return clusters, proteins, functions


def load_mmseqs(path: Path) -> tuple[list[str], list[str], list[str]]:
    """Return parallel lists (cluster, protein_id, function) from an MMseqs2 cluster TSV.

    The representative (first column) is the cluster ID; each member line is
    one protein (MMseqs2 lists the representative as its own member too).
    Functions are empty here and filled from --fasta.
    """
    clusters, proteins = [], []
    with open_text(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            clusters.append(parts[0])
            proteins.append(parts[1])
    return clusters, proteins, [""] * len(proteins)


def load_fasta_functions(path: Path) -> dict[str, str]:
    """Return proteinID -> function from FASTA headers (>proteinID function)."""
    functions = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].strip().split(maxsplit=1)
            if parts:
                functions[parts[0]] = extract_function(parts[1]) if len(parts) > 1 else ""
    return functions


def load_reference(path: Path) -> dict[str, str]:
    reference = {}
    with open_text(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] and parts[1]:
                reference[parts[0]] = parts[1]
    return reference


# ── Plot helpers ──────────────────────────────────────────────────────────────

def style_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.set_axisbelow(True)


def save(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"Written: {path}")


def label_bars(ax, bars, values, horizontal=False):
    for bar, value in zip(bars, values):
        if horizontal:
            ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2, f" {value:,}",
                    va="center", ha="left", fontsize=8, color=INK_2)
        else:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{value:,}",
                    va="bottom", ha="center", fontsize=8, color=INK_2)


# ── Size distribution ─────────────────────────────────────────────────────────

def size_distribution(cluster_sizes: Counter, unit: str, outdir: Path, prefix: str):
    sizes = np.array(sorted(cluster_sizes.values(), reverse=True))
    total_proteins = sizes.sum()
    Unit = unit.capitalize()

    # exact sizes table
    per_size = Counter(sizes.tolist())
    path = outdir / f"{prefix}_cluster_size_distribution.tsv"
    with open(path, "w") as fh:
        fh.write("cluster_size\tn_clusters\tn_proteins\tcumulative_fraction_proteins\n")
        cumulative = 0
        for size in sorted(per_size):
            n_proteins = size * per_size[size]
            cumulative += n_proteins
            fh.write(f"{size}\t{per_size[size]}\t{n_proteins}\t{cumulative / total_proteins:.4f}\n")
    print(f"Written: {path}")

    # binned bar plot
    labels = [label for label, _, _ in SIZE_BINS]
    counts = [int(((sizes >= lo) & (sizes <= hi)).sum()) for _, lo, hi in SIZE_BINS]
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    bars = ax.bar(labels, counts, color=BAR_COLOR, width=0.7)
    label_bars(ax, bars, counts)
    ax.set_xlabel(f"{Unit} size (proteins per {unit})", color=INK_2)
    ax.set_ylabel(f"Number of {unit}s", color=INK_2)
    ax.set_title(f"{Unit} size distribution ({len(sizes):,} {unit}s, {total_proteins:,} proteins)",
                 color=INK, loc="left", fontsize=11)
    ax.yaxis.grid(True, color=GRID)
    style_axes(ax)
    save(fig, outdir / f"{prefix}_cluster_size_distribution.png")

    # rank-abundance (log-log)
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    ax.plot(np.arange(1, len(sizes) + 1), sizes, color=BAR_COLOR, linewidth=1.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"{Unit} rank", color=INK_2)
    ax.set_ylabel(f"{Unit} size (proteins)", color=INK_2)
    ax.set_title(f"{Unit} rank-abundance", color=INK, loc="left", fontsize=11)
    ax.grid(True, color=GRID, which="major")
    style_axes(ax)
    save(fig, outdir / f"{prefix}_cluster_rank_abundance.png")


# ── Functions ─────────────────────────────────────────────────────────────────

def is_informative(function: str, exclude: set[str]) -> bool:
    return bool(function) and function.lower() not in exclude


def top_functions(clusters, functions, exclude, top_n, unit: str, outdir: Path, prefix: str):
    n_proteins = Counter()
    clusters_per_function = defaultdict(set)
    for cluster, function in zip(clusters, functions):
        if is_informative(function, exclude):
            n_proteins[function] += 1
            clusters_per_function[function].add(cluster)

    ranked = n_proteins.most_common()
    path = outdir / f"{prefix}_top_functions.tsv"
    with open(path, "w") as fh:
        fh.write("rank\tfunction\tn_proteins\tn_clusters\n")
        for rank, (function, count) in enumerate(ranked, 1):
            fh.write(f"{rank}\t{function}\t{count}\t{len(clusters_per_function[function])}\n")
    print(f"Written: {path}")

    top = ranked[:top_n]
    if not top:
        print("No informative functions found; top functions plot skipped")
        return
    labels = [textwrap.shorten(f, width=60, placeholder="…") for f, _ in top][::-1]
    values = [c for _, c in top][::-1]
    fig, ax = plt.subplots(figsize=(8, max(2.5, 0.32 * len(top) + 1)))
    bars = ax.barh(labels, values, color=BAR_COLOR, height=0.7)
    label_bars(ax, bars, values, horizontal=True)
    ax.set_xlabel("Number of proteins", color=INK_2)
    ax.set_title(f"Top {len(top)} protein functions", color=INK, loc="left", fontsize=11)
    ax.xaxis.grid(True, color=GRID)
    ax.set_xlim(0, max(values) * 1.15)
    style_axes(ax)
    save(fig, outdir / f"{prefix}_top_functions.png")


def function_purity(clusters, functions, exclude, unit: str, outdir: Path, prefix: str) -> list[float]:
    """Per cluster: share of informative proteins carrying the most common function."""
    by_cluster = defaultdict(list)
    sizes = Counter(clusters)
    for cluster, function in zip(clusters, functions):
        if is_informative(function, exclude):
            by_cluster[cluster].append(function)

    rows, purities = [], []
    for cluster, funcs in by_cluster.items():
        counts = Counter(funcs)
        dominant, n_dominant = counts.most_common(1)[0]
        purity = n_dominant / len(funcs)
        rows.append((cluster, sizes[cluster], len(funcs), len(counts), dominant, purity))
        if len(funcs) >= 2:
            purities.append(purity)

    path = outdir / f"{prefix}_cluster_function_purity.tsv"
    with open(path, "w") as fh:
        fh.write("cluster\tcluster_size\tn_annotated\tn_functions\tdominant_function\tpurity\n")
        for row in sorted(rows, key=lambda r: -r[1]):
            fh.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]:.3f}\n")
    print(f"Written: {path}")

    if purities:
        fig, ax = plt.subplots(figsize=(6.5, 3.6))
        ax.hist(purities, bins=np.linspace(0, 1, 21), color=BAR_COLOR, edgecolor="white")
        ax.set_xlabel("Function purity (share of annotated proteins with the main function)", color=INK_2)
        ax.set_ylabel(f"Number of {unit}s", color=INK_2)
        ax.set_title(f"Function purity of {unit}s with ≥2 annotated proteins (n = {len(purities):,})",
                     color=INK, loc="left", fontsize=11)
        ax.yaxis.grid(True, color=GRID)
        style_axes(ax)
        save(fig, outdir / f"{prefix}_cluster_function_purity.png")
    return purities


# ── Statistics ────────────────────────────────────────────────────────────────

def size_stats(cluster_sizes: Counter) -> list[tuple[str, object]]:
    sizes = np.array(list(cluster_sizes.values()), dtype=float)
    n = sizes.sum()
    p = sizes / n
    shannon = float(-(p * np.log(p)).sum())
    sorted_sizes = np.sort(sizes)
    k = len(sorted_sizes)
    gini = float((2 * np.arange(1, k + 1) - k - 1).dot(sorted_sizes) / (k * sorted_sizes.sum()))
    return [
        ("n_proteins", int(n)),
        ("n_clusters", k),
        ("n_singleton_clusters", int((sizes == 1).sum())),
        ("singleton_fraction_of_clusters", round(float((sizes == 1).mean()), 4)),
        ("proteins_in_singletons_fraction", round(float(sizes[sizes == 1].sum() / n), 4)),
        ("mean_cluster_size", round(float(sizes.mean()), 3)),
        ("median_cluster_size", float(np.median(sizes))),
        ("max_cluster_size", int(sizes.max())),
        ("shannon_entropy_cluster_sizes", round(shannon, 4)),
        ("effective_number_of_clusters", round(math.exp(shannon), 1)),
        ("gini_cluster_sizes", round(gini, 4)),
    ]


def agreement_metrics(pred: list[str], true: list[str], label: str) -> list[tuple[str, object]]:
    try:
        from sklearn import metrics
    except ImportError:
        print("Warning: scikit-learn is not installed; clustering agreement metrics skipped")
        return [("agreement_metrics", "skipped (scikit-learn not installed)")]
    if len(pred) < 2:
        return [("agreement_metrics", f"skipped (fewer than 2 proteins with {label})")]
    homogeneity, completeness, v_measure = metrics.homogeneity_completeness_v_measure(true, pred)
    return [
        ("reference_labels", label),
        ("n_proteins_compared", len(pred)),
        ("n_reference_labels", len(set(true))),
        ("mutual_information", round(metrics.mutual_info_score(true, pred), 4)),
        ("normalized_mutual_information", round(metrics.normalized_mutual_info_score(true, pred), 4)),
        ("adjusted_mutual_information", round(metrics.adjusted_mutual_info_score(true, pred), 4)),
        ("adjusted_rand_score", round(metrics.adjusted_rand_score(true, pred), 4)),
        ("homogeneity", round(homogeneity, 4)),
        ("completeness", round(completeness, 4)),
        ("v_measure", round(v_measure, 4)),
    ]


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    exclude = {e.lower() for e in args.exclude}

    if args.phammseqs:
        clusters, proteins, functions = load_phammseqs(args.phammseqs)
        tool, unit, source = "PhaMMseqs", "pham", args.phammseqs
    else:
        clusters, proteins, functions = load_mmseqs(args.mmseqs)
        tool, unit, source = "MMseqs2", "cluster", args.mmseqs
    if not clusters:
        sys.exit(f"Error: no records found in {source}")

    if args.fasta:
        fasta_functions = load_fasta_functions(args.fasta)
        found = sum(protein in fasta_functions for protein in proteins)
        functions = [fasta_functions.get(protein, function) for protein, function in zip(proteins, functions)]
        print(f"Functions from {args.fasta.name}: {found:,} of {len(proteins):,} clustered proteins found")
        if found < len(proteins):
            print(f"Warning: {len(proteins) - found:,} clustered proteins are not in {args.fasta.name}")

    cluster_sizes = Counter(clusters)
    print(f"Loaded {len(proteins):,} proteins in {len(cluster_sizes):,} {unit}s ({tool})")

    size_distribution(cluster_sizes, unit, args.outdir, args.prefix)
    top_functions(clusters, functions, exclude, args.top, unit, args.outdir, args.prefix)
    purities = function_purity(clusters, functions, exclude, unit, args.outdir, args.prefix)

    stats = [("clustering_tool", tool)] + size_stats(cluster_sizes)
    n_informative = sum(is_informative(f, exclude) for f in functions)
    stats += [
        ("n_proteins_with_informative_function", n_informative),
        ("informative_function_fraction", round(n_informative / len(functions), 4)),
        ("mean_function_purity", round(float(np.mean(purities)), 4) if purities else "NA"),
    ]

    if args.reference:
        reference = load_reference(args.reference)
        pairs = [(c, reference[protein]) for c, protein in zip(clusters, proteins) if protein in reference]
        label = f"reference ({args.reference.name})"
    else:
        pairs = [(c, f) for c, f in zip(clusters, functions) if is_informative(f, exclude)]
        label = "function (informative only)"
    stats += agreement_metrics([p for p, _ in pairs], [t for _, t in pairs], label)

    path = args.outdir / f"{args.prefix}_clustering_stats.tsv"
    with open(path, "w") as fh:
        fh.write("metric\tvalue\n")
        for name, value in stats:
            fh.write(f"{name}\t{value}\n")
    print(f"Written: {path}")


if __name__ == "__main__":
    main()
