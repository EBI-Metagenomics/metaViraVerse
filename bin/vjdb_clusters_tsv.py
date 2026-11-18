#!/usr/bin/env python3
"""Convert the vjdb-clustering membership CSV to the cluster TSVs the pipeline reads.

The membership CSV (vjdb_id, copy_of, ..., cluster_rep) covers the whole catalogue.
It is written as the vclust cluster TSV (object, cluster).

With --events (the cluster_events.tsv of a vjdb update), two more files describe what
the update changed:

  <prefix>_clusters_to_annotate.tsv  rows of the clusters created by the update ('new').
      Their representatives are new sequences and have no annotations yet. A cluster
      that grew or absorbed others keeps an old representative, which is already
      annotated, so it is left out.
  <prefix>_dropped_reps.txt  old representatives of clusters merged into another one.
      Their annotations no longer belong to a representative.

Without --events (initial clustering) every cluster is annotated and none is dropped.

CHOOSE_SEQUENCES has already removed identical sequences (sha256 of the upper-cased
sequence), which is stricter than the case-sensitive md5 vjdb-clustering uses. So vjdb
must find no copies: every row has copy_of == vjdb_id. A copy here means the two hashes
disagree, and a later update (which reads seq_hash from this CSV) would see it too, so
the script fails instead of passing it on.
"""

import argparse
import csv
import gzip
import sys

CSV_COLUMNS = ("vjdb_id", "copy_of", "cluster_rep")
EVENT_COLUMNS = ("cluster_rep", "event", "old_cluster_reps")
EVENTS = ("new", "grown", "merged")
MAX_EXAMPLES = 10


def open_text(path):
    with open(path, "rb") as fh:
        is_gzip = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", newline="") if is_gzip else open(path, newline="")


def check_columns(fieldnames, required, path):
    missing = [c for c in required if c not in (fieldnames or [])]
    if missing:
        sys.exit(f"Error: {path} is missing column(s) {', '.join(missing)}")


def examples(items):
    items = list(items)
    shown = ", ".join(items[:MAX_EXAMPLES])
    return shown + (", ..." if len(items) > MAX_EXAMPLES else "")


def read_clusters(path):
    """Return [(vjdb_id, cluster_rep)] in file order, after checking the table."""
    rows = []
    seen = set()
    copies = []
    duplicated = []
    with open_text(path) as fh:
        reader = csv.DictReader(fh)
        check_columns(reader.fieldnames, CSV_COLUMNS, path)
        for row in reader:
            if row["vjdb_id"] == "vjdb_id":  # repeated header
                continue
            seq_id, rep = row["vjdb_id"], row["cluster_rep"]
            if row["copy_of"] != seq_id:
                copies.append(f"{seq_id}={row['copy_of']}")
            if seq_id in seen:
                duplicated.append(seq_id)
            seen.add(seq_id)
            rows.append((seq_id, rep))

    if copies:
        sys.exit(f"Error: vjdb-clustering found {len(copies)} identical sequence(s) that CHOOSE_SEQUENCES kept, "
                 f"e.g. {examples(copies)}")
    if duplicated:
        sys.exit(f"Error: {path} lists {len(duplicated)} vjdb_id(s) more than once, e.g. {examples(duplicated)}")
    cluster_of = dict(rows)
    bad_reps = sorted({rep for _, rep in rows if cluster_of.get(rep) != rep})
    if bad_reps:
        sys.exit(f"Error: {len(bad_reps)} cluster_rep(s) in {path} are not members of their own cluster, "
                 f"e.g. {examples(bad_reps)}")
    return rows


def read_events(path, reps):
    """Return (new_reps, dropped_reps) from the cluster_events.tsv of a vjdb update."""
    new_reps = set()
    dropped = set()
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        check_columns(reader.fieldnames, EVENT_COLUMNS, path)
        for row in reader:
            if row["cluster_rep"] == "cluster_rep":  # repeated header
                continue
            rep, event = row["cluster_rep"], row["event"]
            old = [r for r in (row["old_cluster_reps"] or "").split(";") if r]
            if event not in EVENTS:
                sys.exit(f"Error: {path}: unknown event '{event}' for {rep} (expected one of {', '.join(EVENTS)})")
            if rep not in reps:
                sys.exit(f"Error: {path}: {rep} is not a cluster_rep of the membership CSV")
            if event == "new":
                if old:
                    sys.exit(f"Error: {path}: new cluster {rep} lists old representatives {';'.join(old)}")
                new_reps.add(rep)
            elif event == "merged":
                if rep not in old:
                    sys.exit(f"Error: {path}: merged cluster {rep} is not one of its old representatives")
                dropped.update(r for r in old if r != rep)

    still_reps = sorted(dropped & reps)
    if still_reps:
        sys.exit(f"Error: {path}: {len(still_reps)} merged old representative(s) are still cluster_reps, "
                 f"e.g. {examples(still_reps)}")
    return new_reps, dropped


def write_tsv(path, rows):
    with open(path, "w") as out:
        out.write("object\tcluster\n")
        for seq_id, rep in rows:
            out.write(f"{seq_id}\t{rep}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--clusters-csv", required=True, help="vjdb membership CSV (plain or gzip)")
    parser.add_argument("--events", required=False,
                        help="cluster_events.tsv of a vjdb update; without it every cluster is annotated")
    parser.add_argument("--prefix", required=True, help="Prefix of the output files")
    args = parser.parse_args(argv)

    rows = read_clusters(args.clusters_csv)
    reps = {rep for _, rep in rows}

    if args.events:
        new_reps, dropped = read_events(args.events, reps)
        to_annotate = [(seq_id, rep) for seq_id, rep in rows if rep in new_reps]
    else:
        dropped = set()
        to_annotate = rows

    write_tsv(f"{args.prefix}_clusters.tsv", rows)
    write_tsv(f"{args.prefix}_clusters_to_annotate.tsv", to_annotate)
    with open(f"{args.prefix}_dropped_reps.txt", "w") as out:
        for rep in sorted(dropped):
            out.write(f"{rep}\n")

    n_annotate = len({rep for _, rep in to_annotate})
    print(f"{len(rows)} sequences in {len(reps)} clusters; {n_annotate} cluster(s) to annotate; "
          f"{len(dropped)} old representative(s) dropped")
    if args.events and n_annotate == 0:
        print("Warning: the update created no new cluster; there is nothing to annotate", file=sys.stderr)


if __name__ == "__main__":
    main()
