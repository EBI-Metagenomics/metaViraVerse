#!/usr/bin/env python3
"""
Combine plaSquid's mobility_stats.json with MOB-suite's predicted mobility
counts into one summary JSON.

--mobility-json (plaSquid's mobility_stats.json, from plasquid_rip_extraction.R):
its keys (plasquid_total_number, plasquid_rip_domain_present, ...) are copied
into the output unchanged.

--mob-report (MOB-suite's mob_typer `--out_file` report, run with `--multi`):
one row per plasmid. Its `predicted_mobility` column (conjugative /
mobilizable / non-mobilizable) is counted per category and added as
mobsuite_total, mobsuite_conjugative, mobsuite_mobilizable and
mobsuite_non_mobilizable. Repeated header lines (from concatenating chunked
reports) are skipped; rows with an empty or unrecognised value count towards
mobsuite_total only.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys

CONJUGATIVE = "conjugative"
MOBILIZABLE = "mobilizable"
NON_MOBILIZABLE = "non_mobilizable"
CATEGORIES = (CONJUGATIVE, MOBILIZABLE, NON_MOBILIZABLE)


def load_plasquid_stats(path: str) -> dict:
    """Return plaSquid's mobility_stats.json as a dict ({} if unreadable)."""
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        print(f"Warning: could not read --mobility-json {path}: {e}", file=sys.stderr)
        return {}


def count_mob_report(path: str) -> dict[str, int]:
    """Return mobsuite_total and per-category counts of `predicted_mobility`.

    MOB-suite writes "non-mobilizable" with a hyphen; it is normalised to
    "non_mobilizable" to match the output key.
    """
    counts = {category: 0 for category in CATEGORIES}
    total = 0
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id or sample_id == "sample_id":
                continue
            total += 1
            mobility = (row.get("predicted_mobility") or "").strip().lower().replace("-", "_")
            if mobility in counts:
                counts[mobility] += 1
    return {
        "mobsuite_total": total,
        "mobsuite_conjugative": counts[CONJUGATIVE],
        "mobsuite_mobilizable": counts[MOBILIZABLE],
        "mobsuite_non_mobilizable": counts[NON_MOBILIZABLE],
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine plaSquid's mobility_stats.json with MOB-suite's predicted mobility counts.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-j", "--mobility-json", default=None,
                        help="plaSquid's mobility_stats.json (optional).")
    parser.add_argument("-m", "--mob-report", default=None,
                        help="MOB-suite mob_typer report (--out_file, run with --multi) (optional).")
    parser.add_argument("-o", "--output", default="mobility_stats_final.json",
                        help="Output JSON file (default: mobility_stats_final.json).")
    return parser.parse_args()


def main():
    args = parse_args()

    stats = load_plasquid_stats(args.mobility_json) if args.mobility_json else {}
    print(f"plaSquid stats loaded: {len(stats)} keys", file=sys.stderr)

    if args.mob_report:
        mobsuite_stats = count_mob_report(args.mob_report)
        print(f"MOB-suite plasmids counted: {mobsuite_stats['mobsuite_total']}", file=sys.stderr)
        stats.update(mobsuite_stats)

    with open(args.output, "w") as f:
        json.dump(stats, f, indent=2)
        f.write("\n")
    print(f"Output written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
