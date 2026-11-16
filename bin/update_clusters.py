#!/usr/bin/env python3
"""Viral clusters of a catalogue update (--mode update) that do not come from vjdb-clustering.

Subcommands
  to-vjdb  --clusters-tsv FILE --fasta FILE [FILE ...] --prefix PREFIX
      Convert the clustering of an existing catalogue that was not clustered with
      vjdb-clustering (--viral_clustering default) into the tables a vjdb-clustering
      update starts from:
        {prefix}_vjdb_clusters.csv   vjdb_id, copy_of, seq_gc, seq_len, seq_hash, cluster_rep
        {prefix}_vjdb_reps.fna.gz    the representative sequences
      seq_hash is the MD5 of the sequence as written (case-sensitive) and seq_gc the GC
      fraction without ambiguous bases, as vjdb-clustering computes them. Identical
      sequences were removed by the pipeline before clustering, so copy_of = vjdb_id.

  compare  --old-clusters FILE --new-clusters FILE --prefix PREFIX
      Compare the clustering of the whole updated catalogue (existing and new sequences,
      clustered again from scratch) with the existing catalogue's:
        {prefix}_clusters_to_annotate.tsv  clusters whose representative was not a
                                           representative before (vclust format)
        {prefix}_dropped_reps.txt          existing representatives that no longer are
        {prefix}_cluster_events.tsv        cluster_rep, event, old_cluster_reps, n_existing, n_new
      Events (unchanged clusters are not listed):
        new       a new representative, no existing members
        new_rep   a new representative of a cluster with existing members (the
                  representative changed, or existing members moved to it)
        merged    an existing representative kept, with members of other existing clusters
        split     an existing representative kept, some of its members moved to other clusters
        grown     an existing cluster with new members only

Cluster tables: vclust (header 'object<TAB>cluster', one row per sequence) or the CheckV
ANI clustering (no header, 'representative<TAB>comma-separated members'). Plain or gzip.
"""

import argparse
import csv
import gzip
import hashlib
import sys
from collections import defaultdict

MAX_EXAMPLES = 10


def open_text(path):
    with open(path, "rb") as fh:
        is_gzip = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", newline="") if is_gzip else open(path, newline="")


def examples(items):
    items = sorted(items)
    return ", ".join(items[:MAX_EXAMPLES]) + (", ..." if len(items) > MAX_EXAMPLES else "")


def read_clusters(path):
    """member -> representative, in file order. Every representative is its own member."""
    member_to_rep = {}

    def add(member, rep):
        if member in member_to_rep and member_to_rep[member] != rep:
            sys.exit(f"Error: {path}: {member} is in the clusters of {member_to_rep[member]} and {rep}")
        member_to_rep[member] = rep

    with open_text(path) as fh:
        lines = [line.rstrip("\r\n") for line in fh if line.strip()]
    if not lines:
        sys.exit(f"Error: {path} is empty")
    vclust = lines[0].split("\t")[:2] == ["object", "cluster"]
    for line in lines[1:] if vclust else lines:
        cols = line.split("\t")
        if len(cols) < 2:
            sys.exit(f"Error: {path}: expected two columns: {line}")
        if vclust:
            if cols[:2] == ["object", "cluster"]:  # repeated header
                continue
            add(cols[0].strip(), cols[1].strip())
        else:
            rep = cols[0].strip()
            add(rep, rep)
            for member in cols[1].split(","):
                if member.strip():
                    add(member.strip(), rep)
    missing = {rep for rep in member_to_rep.values()} - {m for m, r in member_to_rep.items() if m == r}
    if missing:
        sys.exit(f"Error: {path}: {len(missing)} representative(s) are not members of their own cluster, "
                 f"e.g. {examples(missing)}")
    return member_to_rep


def read_fasta(paths):
    """id (first word of the header) -> sequence; fails on a duplicated ID."""
    seqs = {}
    for path in paths:
        seq_id, parts = None, []
        with open_text(path) as fh:
            for line in list(fh) + [">"]:
                if line.startswith(">"):
                    if seq_id is not None:
                        if seq_id in seqs:
                            sys.exit(f"Error: sequence {seq_id} occurs more than once in {', '.join(paths)}")
                        seqs[seq_id] = "".join(parts)
                    words = line[1:].split()
                    seq_id, parts = (words[0] if words else None), []
                elif seq_id is not None:
                    parts.append(line.strip())
    return seqs


def gc_fraction(seq):
    """Bio.SeqUtils.gc_fraction(seq, ambiguous="remove")."""
    gc = sum(seq.count(x) for x in "CGScgs")
    length = sum(seq.count(x) for x in "ACGTSWacgtsw")
    return gc / length if length else 0.0


# ---------------------------------------------------------------------------
# to-vjdb
# ---------------------------------------------------------------------------

def to_vjdb(clusters_tsv, fasta, prefix):
    member_to_rep = read_clusters(clusters_tsv)
    seqs = read_fasta(fasta)
    no_seq = set(member_to_rep) - set(seqs)
    if no_seq:
        sys.exit(f"Error: {len(no_seq)} clustered sequence(s) are not in {', '.join(fasta)}, e.g. {examples(no_seq)}")
    unclustered = set(seqs) - set(member_to_rep)
    if unclustered:
        sys.exit(f"Error: {len(unclustered)} sequence(s) are in no cluster of {clusters_tsv}, e.g. {examples(unclustered)}")

    csv_path = f"{prefix}_vjdb_clusters.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["vjdb_id", "copy_of", "seq_gc", "seq_len", "seq_hash", "cluster_rep"])
        for member, rep in member_to_rep.items():
            seq = seqs[member]
            writer.writerow([member, member, repr(gc_fraction(seq)), len(seq),
                             hashlib.md5(seq.encode()).hexdigest(), rep])

    reps = list(dict.fromkeys(member_to_rep.values()))
    reps_path = f"{prefix}_vjdb_reps.fna.gz"
    with gzip.open(reps_path, "wt") as fh:
        for rep in reps:
            fh.write(f">{rep}\n{seqs[rep]}\n")
    print(f"{csv_path}: {len(member_to_rep)} sequences in {len(reps)} clusters; {reps_path}")


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------

def compare(old_path, new_path, prefix):
    old = read_clusters(old_path)
    new = read_clusters(new_path)
    lost = set(old) - set(new)
    if lost:
        sys.exit(f"Error: {len(lost)} sequence(s) of the existing catalogue are not in the new clustering "
                 f"{new_path}, e.g. {examples(lost)}")

    old_reps = set(old.values())
    old_members = defaultdict(set)
    for member, rep in old.items():
        old_members[rep].add(member)
    new_members = defaultdict(list)
    for member, rep in new.items():
        new_members[rep].append(member)

    events = []
    annotate = []
    for rep, members in new_members.items():
        existing = [m for m in members if m in old]
        n_new = len(members) - len(existing)
        old_clusters = sorted({old[m] for m in existing})
        if rep not in old_reps:
            event = "new_rep" if existing else "new"
            annotate.extend((m, rep) for m in members)
        elif len(old_clusters) > 1:
            event = "merged"
        elif set(existing) != old_members[rep]:
            event = "split"
        elif n_new:
            event = "grown"
        else:
            continue
        events.append([rep, event, ";".join(old_clusters), len(existing), n_new])

    dropped = sorted(old_reps - set(new_members))

    with open(f"{prefix}_clusters_to_annotate.tsv", "w") as fh:
        fh.write("object\tcluster\n")
        fh.writelines(f"{m}\t{r}\n" for m, r in annotate)
    with open(f"{prefix}_dropped_reps.txt", "w") as fh:
        fh.writelines(f"{r}\n" for r in dropped)
    with open(f"{prefix}_cluster_events.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["cluster_rep", "event", "old_cluster_reps", "n_existing", "n_new"])
        writer.writerows(events)

    counts = defaultdict(int)
    for e in events:
        counts[e[1]] += 1
    print(f"{len(new_members)} clusters ({len(old_reps)} before); "
          + ", ".join(f"{k}: {v}" for k, v in sorted(counts.items()))
          + f"; {len({r for _, r in annotate})} to annotate, {len(dropped)} representative(s) dropped")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("to-vjdb", help="vjdb-clustering tables from another clustering")
    p.add_argument("--clusters-tsv", required=True)
    p.add_argument("--fasta", required=True, nargs="+", help="every clustered sequence (plain or gzip)")
    p.add_argument("--prefix", required=True)

    p = sub.add_parser("compare", help="what a new clustering of the whole catalogue changed")
    p.add_argument("--old-clusters", required=True, help="clusters of the existing catalogue")
    p.add_argument("--new-clusters", required=True, help="clusters of the whole updated catalogue")
    p.add_argument("--prefix", required=True)

    args = parser.parse_args(argv)
    if args.command == "to-vjdb":
        to_vjdb(args.clusters_tsv, args.fasta, args.prefix)
    else:
        compare(args.old_clusters, args.new_clusters, args.prefix)


if __name__ == "__main__":
    main()
