#!/usr/bin/env python3
"""
Recompute plaSquid's mobility_stats.json using MOB-suite's biomarker report
as the missing mating-pair-formation (MPF) evidence source.

plaSquid's own MOBSEARCH only detects the relaxase (MOB) gene, not MPF/T4SS
conjugation machinery, so plasquid_rip_extraction.R's mobility_stats.json can
never place a contig in "conjugative" from plaSquid's evidence alone (see its
own header comment). MOB-suite's mob_typer, run separately on the same
representative plasmids, does search for MPF genes and records any hit as a
"mate-pair-formation" row (keyed by `sseqid`, the plasmid/contig id) in its
--biomarker_report_file output -- see annotate_plasmid_gff.py's docstring for
the same column layout.

This script takes plaSquid's own per-contig classification
(--classification, mobility_classification.tsv: contig, category) and
promotes every contig still marked "mobilizable" to "conjugative" if
MOB-suite's biomarker report also found MPF evidence on it, then re-derives
mobility_stats.json's summary counts/percentages from that updated
classification -- the same JSON shape plaSquid produces (--mobility-json, if
given, is only cross-checked against it for a sanity warning, not used for
the actual recount, since it carries no per-contig detail to update from).

Contigs plaSquid never classified are ignored: this script only ever
promotes mobilizable -> conjugative, it never introduces new contigs or
changes non_mobilizable calls.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys

MOBILIZABLE = "mobilizable"
CONJUGATIVE = "conjugative"
NON_MOBILIZABLE = "non_mobilizable"
CATEGORIES = (CONJUGATIVE, MOBILIZABLE, NON_MOBILIZABLE)


def load_classification(path: str) -> dict[str, str]:
    """Return contig -> category from plaSquid's mobility_classification.tsv."""
    result: dict[str, str] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            contig = (row.get("contig") or "").strip()
            category = (row.get("category") or "").strip()
            if contig and category:
                result[contig] = category
    return result


def load_mpf_contigs(path: str) -> set[str]:
    """Return contig ids with a mate-pair-formation hit in MOB-suite's biomarker report."""
    contigs: set[str] = set()
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (row.get("biomarker") or "").strip() == "mate-pair-formation":
                contig = (row.get("sseqid") or "").strip()
                if contig:
                    contigs.add(contig)
    return contigs


def promote_conjugative(classification: dict[str, str], mpf_contigs: set[str]) -> dict[str, str]:
    """Return an updated contig -> category mapping.

    Promotes "mobilizable" -> "conjugative" wherever MOB-suite found MPF
    evidence on that contig; "non_mobilizable" contigs are never touched --
    the MOB-suite MPF hit alone doesn't establish that a relaxase is present
    too, which "conjugative" requires (MOB + MPF).
    """
    return {
        contig: CONJUGATIVE if category == MOBILIZABLE and contig in mpf_contigs else category
        for contig, category in classification.items()
    }


def write_mobility_stats_json(classification: dict[str, str], mpf_data_available: bool, path: str) -> None:
    """Write the same JSON shape as plasquid_rip_extraction.R's write_mobility_stats_json()."""
    total = len(classification)
    counts = {c: 0 for c in CATEGORIES}
    for category in classification.values():
        if category in counts:
            counts[category] += 1

    def pct(n: int) -> float:
        return round(100 * n / total, 2) if total > 0 else 0

    categories = {c: {"count": counts[c], "percent": pct(counts[c])} for c in CATEGORIES}
    if not mpf_data_available:
        categories[CONJUGATIVE]["note"] = (
            "No MOB-suite mate-pair-formation evidence was supplied; this count is "
            "carried over unchanged from plaSquid's own classification, where "
            "'conjugative' always stays at 0 (plaSquid detects the MOB relaxase gene "
            "only, not MPF/T4SS conjugation machinery)."
        )

    payload = {
        "total_contigs": total,
        "mpf_data_available": mpf_data_available,
        "categories": categories,
    }
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")


def sanity_check_against_plasquid_json(path: str, expected_total: int) -> None:
    """Warn (don't fail) if --mobility-json's total_contigs disagrees with --classification."""
    try:
        with open(path) as f:
            plasquid_json = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        print(f"Warning: could not read --mobility-json {path}: {e}", file=sys.stderr)
        return

    reported_total = plasquid_json.get("total_contigs")
    if reported_total != expected_total:
        print(
            f"Warning: {path} reports total_contigs={reported_total}, but "
            f"--classification has {expected_total} contigs; proceeding with the "
            "per-contig classification as authoritative.",
            file=sys.stderr,
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Recompute plaSquid's mobility_stats.json using MOB-suite's MPF biomarker evidence.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-c", "--classification", required=True,
                        help="plaSquid's mobility_classification.tsv (contig, category).")
    parser.add_argument("-j", "--mobility-json", default=None,
                        help="plaSquid's mobility_stats.json, for a total_contigs sanity check (optional).")
    parser.add_argument("-b", "--biomarker-report", default=None,
                        help="MOB-suite plasmids_biomarker_report.txt (optional).")
    parser.add_argument("-o", "--output", default="mobility_stats_final.json",
                        help="Output JSON file (default: mobility_stats_final.json).")
    return parser.parse_args()


def main():
    args = parse_args()

    classification = load_classification(args.classification)
    print(f"Contigs loaded from plaSquid classification: {len(classification)}", file=sys.stderr)

    if args.mobility_json:
        sanity_check_against_plasquid_json(args.mobility_json, len(classification))

    mpf_contigs = load_mpf_contigs(args.biomarker_report) if args.biomarker_report else set()
    print(f"Contigs with MOB-suite MPF evidence: {len(mpf_contigs)}", file=sys.stderr)

    updated = promote_conjugative(classification, mpf_contigs)
    promoted = sum(
        1 for contig, category in updated.items()
        if category == CONJUGATIVE and classification.get(contig) == MOBILIZABLE
    )
    print(f"Contigs promoted mobilizable -> conjugative: {promoted}", file=sys.stderr)

    write_mobility_stats_json(updated, mpf_data_available=len(mpf_contigs) > 0, path=args.output)
    print(f"Output written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
