#!/usr/bin/env python3
"""Unit tests for update_clusters.py."""

import csv
import gzip
import hashlib
import os
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import update_clusters as uc

SEQS = {"V1": "ACGTACGTGG", "V2": "acgtNNcc", "V3": "GGGGCCCCAT", "V4": "ATATATAT"}


def fasta(ids):
    return "".join(f">{i} some description\n{SEQS[i][:4]}\n{SEQS[i][4:]}\n" for i in ids)


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cwd = os.getcwd()
        os.chdir(self.tmp.name)

    def tearDown(self):
        os.chdir(self.cwd)
        self.tmp.cleanup()

    @staticmethod
    def write(name, text):
        if name.endswith(".gz"):
            with gzip.open(name, "wt") as fh:
                fh.write(text)
        else:
            Path(name).write_text(text)
        return name

    @staticmethod
    def read_tsv(name):
        with open(name) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))


class TestReadClusters(Base):
    def test_vclust_and_checkv_formats(self):
        self.write("v.tsv", "object\tcluster\nV1\tV1\nV2\tV1\nV3\tV3\n")
        self.write("c.tsv", "V1\tV1,V2\nV3\tV3\n")
        self.assertEqual(uc.read_clusters("v.tsv"), {"V1": "V1", "V2": "V1", "V3": "V3"})
        self.assertEqual(uc.read_clusters("c.tsv"), {"V1": "V1", "V2": "V1", "V3": "V3"})

    def test_rep_outside_its_cluster_fails(self):
        self.write("v.tsv", "object\tcluster\nV2\tV1\n")
        with self.assertRaises(SystemExit) as cm:
            uc.read_clusters("v.tsv")
        self.assertIn("not members of their own cluster", str(cm.exception))

    def test_member_in_two_clusters_fails(self):
        self.write("c.tsv", "V1\tV1,V2\nV3\tV3,V2\n")
        with self.assertRaises(SystemExit):
            uc.read_clusters("c.tsv")


class TestToVjdb(Base):
    def setUp(self):
        super().setUp()
        self.write("clusters.tsv", "object\tcluster\nV1\tV1\nV2\tV1\nV3\tV3\n")
        self.write("viral.fna.gz", fasta(["V1", "V3"]))
        self.write("prophages.fna", fasta(["V2"]))

    def test_tables(self):
        uc.main(["to-vjdb", "--clusters-tsv", "clusters.tsv", "--fasta", "viral.fna.gz", "prophages.fna", "--prefix", "viruses"])
        with open("viruses_vjdb_clusters.csv") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual([r["vjdb_id"] for r in rows], ["V1", "V2", "V3"])
        self.assertEqual([r["cluster_rep"] for r in rows], ["V1", "V1", "V3"])
        self.assertTrue(all(r["copy_of"] == r["vjdb_id"] for r in rows))
        v2 = rows[1]
        self.assertEqual(v2["seq_len"], "8")
        self.assertEqual(v2["seq_hash"], hashlib.md5(b"acgtNNcc").hexdigest())   # case-sensitive, as written
        self.assertAlmostEqual(float(v2["seq_gc"]), 4 / 6)                       # N removed
        with gzip.open("viruses_vjdb_reps.fna.gz", "rt") as fh:
            self.assertEqual(fh.read(), f">V1\n{SEQS['V1']}\n>V3\n{SEQS['V3']}\n")

    def test_accepted_by_vjdb_clusters_tsv(self):
        import vjdb_clusters_tsv
        uc.main(["to-vjdb", "--clusters-tsv", "clusters.tsv", "--fasta", "viral.fna.gz", "prophages.fna", "--prefix", "viruses"])
        vjdb_clusters_tsv.main(["--clusters-csv", "viruses_vjdb_clusters.csv", "--prefix", "check"])
        self.assertEqual(self.read_tsv("check_clusters.tsv"),
                         [{"object": "V1", "cluster": "V1"}, {"object": "V2", "cluster": "V1"}, {"object": "V3", "cluster": "V3"}])

    def test_gc_matches_biopython(self):
        try:
            from Bio.SeqUtils import gc_fraction
        except ImportError:
            self.skipTest("Biopython not installed")
        for seq in SEQS.values():
            self.assertAlmostEqual(uc.gc_fraction(seq), gc_fraction(seq))

    def test_sequence_missing_fails(self):
        with self.assertRaises(SystemExit) as cm:
            uc.main(["to-vjdb", "--clusters-tsv", "clusters.tsv", "--fasta", "viral.fna.gz", "--prefix", "viruses"])
        self.assertIn("not in", str(cm.exception))

    def test_unclustered_sequence_fails(self):
        self.write("more.fna", fasta(["V4"]))
        with self.assertRaises(SystemExit) as cm:
            uc.main(["to-vjdb", "--clusters-tsv", "clusters.tsv", "--fasta", "viral.fna.gz", "prophages.fna", "more.fna",
                     "--prefix", "viruses"])
        self.assertIn("in no cluster", str(cm.exception))


class TestCompare(Base):
    # existing: V1 {V1, V2}, V3 {V3}, V4 {V4}, V5 {V5, V6}, V8 {V8, V9}
    OLD = "object\tcluster\nV1\tV1\nV2\tV1\nV3\tV3\nV4\tV4\nV5\tV5\nV6\tV5\nV8\tV8\nV9\tV8\n"

    def compare(self, new):
        self.write("old.tsv", self.OLD)
        self.write("new.tsv", new)
        uc.main(["compare", "--old-clusters", "old.tsv", "--new-clusters", "new.tsv", "--prefix", "viruses"])
        events = {r["cluster_rep"]: r for r in self.read_tsv("viruses_cluster_events.tsv")}
        annotate = self.read_tsv("viruses_clusters_to_annotate.tsv")
        dropped = Path("viruses_dropped_reps.txt").read_text().split()
        return events, annotate, dropped

    def test_events(self):
        # V1 grew (N1), V3 absorbed V4 (merged), V5's member V6 became the representative of
        # a cluster with V9 (new_rep; V8 lost V9: split), N2 is a new singleton
        events, annotate, dropped = self.compare(
            "object\tcluster\nV1\tV1\nV2\tV1\nN1\tV1\nV3\tV3\nV4\tV3\nV5\tV5\nV6\tV6\nV9\tV6\nV8\tV8\nN2\tN2\n"
        )
        self.assertEqual({k: v["event"] for k, v in events.items()},
                         {"V1": "grown", "V3": "merged", "V5": "split", "V6": "new_rep", "V8": "split", "N2": "new"})
        self.assertEqual(events["V3"]["old_cluster_reps"], "V3;V4")
        self.assertEqual(events["V6"]["old_cluster_reps"], "V5;V8")
        self.assertEqual((events["V1"]["n_existing"], events["V1"]["n_new"]), ("2", "1"))
        self.assertEqual(annotate, [{"object": "V6", "cluster": "V6"}, {"object": "V9", "cluster": "V6"},
                                    {"object": "N2", "cluster": "N2"}])
        self.assertEqual(dropped, ["V4"])

    def test_representative_replaced(self):
        # a longer new sequence N1 became the representative of V1's cluster
        events, annotate, dropped = self.compare(
            "object\tcluster\nN1\tN1\nV1\tN1\nV2\tN1\nV3\tV3\nV4\tV4\nV5\tV5\nV6\tV5\nV8\tV8\nV9\tV8\n"
        )
        self.assertEqual(events, {"N1": {"cluster_rep": "N1", "event": "new_rep", "old_cluster_reps": "V1",
                                         "n_existing": "2", "n_new": "1"}})
        self.assertEqual({r["cluster"] for r in annotate}, {"N1"})
        self.assertEqual(dropped, ["V1"])

    def test_unchanged(self):
        events, annotate, dropped = self.compare(self.OLD)
        self.assertEqual((events, annotate, dropped), ({}, [], []))

    def test_existing_sequence_missing_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.compare("object\tcluster\nV1\tV1\n")
        self.assertIn("not in the new clustering", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
