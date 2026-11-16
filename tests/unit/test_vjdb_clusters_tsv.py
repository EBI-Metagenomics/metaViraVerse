#!/usr/bin/env python3
"""Unit tests for vjdb_clusters_tsv.py."""

import gzip
import os
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import vjdb_clusters_tsv as vct

HEADER = "vjdb_id,copy_of,seq_gc,seq_len,seq_hash,cluster_rep\n"
EVENTS_HEADER = "cluster_rep\tevent\told_cluster_reps\tn_existing\tn_new\n"

# Existing catalogue: clusters V1 {V1, V2}, V3 {V3}, V4 {V4}.
# Update with V5..V8: V5 joins V1 (grown), V6 merges V3 and V4 into V3 (merged),
# V7 and V8 form a new cluster with representative V7.
UPDATED_CSV = HEADER + (
    "V1,V1,0.5,100,h1,V1\n"
    "V2,V2,0.5,100,h2,V1\n"
    "V3,V3,0.5,100,h3,V3\n"
    "V4,V4,0.5,100,h4,V3\n"
    "V5,V5,0.5,100,h5,V1\n"
    "V6,V6,0.5,100,h6,V3\n"
    "V7,V7,0.5,100,h7,V7\n"
    "V8,V8,0.5,100,h8,V7\n"
)
EVENTS = EVENTS_HEADER + (
    "V3\tmerged\tV3;V4\t3\t1\n"
    "V1\tgrown\tV1\t2\t1\n"
    "V7\tnew\t\t0\t2\n"
)


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.cwd = os.getcwd()
        os.chdir(self.dir)

    def tearDown(self):
        os.chdir(self.cwd)
        self.tmp.cleanup()

    def write(self, name, text, compress=False):
        path = self.dir / name
        if compress:
            with gzip.open(path, "wt") as fh:
                fh.write(text)
        else:
            path.write_text(text)
        return str(path)

    def tsv_rows(self, name):
        lines = (self.dir / name).read_text().splitlines()
        self.assertEqual(lines[0], "object\tcluster")
        return [tuple(line.split("\t")) for line in lines[1:]]

    def run_main(self, csv_text, events_text=None, compress=False):
        argv = ["--clusters-csv", self.write("clusters.csv", csv_text, compress), "--prefix", "viruses"]
        if events_text is not None:
            argv += ["--events", self.write("events.tsv", events_text, compress)]
        vct.main(argv)


class TestInitial(Base):
    def test_all_clusters_annotated_none_dropped(self):
        self.run_main(HEADER + "V1,V1,0.5,100,h1,V1\nV2,V2,0.5,100,h2,V1\nV3,V3,0.5,100,h3,V3\n")
        expected = [("V1", "V1"), ("V2", "V1"), ("V3", "V3")]
        self.assertEqual(self.tsv_rows("viruses_clusters.tsv"), expected)
        self.assertEqual(self.tsv_rows("viruses_clusters_to_annotate.tsv"), expected)
        self.assertEqual((self.dir / "viruses_dropped_reps.txt").read_text(), "")

    def test_copy_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_main(HEADER + "V1,V1,0.5,100,h1,V1\nV2,V1,0.5,100,h1,V1\n")
        self.assertIn("identical sequence", str(cm.exception))

    def test_missing_column_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_main("vjdb_id,copy_of\nV1,V1\n")
        self.assertIn("cluster_rep", str(cm.exception))

    def test_duplicated_id_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_main(HEADER + "V1,V1,0.5,100,h1,V1\nV1,V1,0.5,100,h1,V1\n")
        self.assertIn("more than once", str(cm.exception))

    def test_rep_outside_its_cluster_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_main(HEADER + "V1,V1,0.5,100,h1,V2\nV2,V2,0.5,100,h2,V1\n")
        self.assertIn("not members of their own cluster", str(cm.exception))

    def test_repeated_header_skipped(self):
        self.run_main(HEADER + "V1,V1,0.5,100,h1,V1\n" + HEADER + "V2,V2,0.5,100,h2,V2\n")
        self.assertEqual(self.tsv_rows("viruses_clusters.tsv"), [("V1", "V1"), ("V2", "V2")])


class TestUpdate(Base):
    def test_split(self):
        self.run_main(UPDATED_CSV, EVENTS)
        self.assertEqual(len(self.tsv_rows("viruses_clusters.tsv")), 8)
        # only the new cluster: its representative is a new sequence
        self.assertEqual(self.tsv_rows("viruses_clusters_to_annotate.tsv"), [("V7", "V7"), ("V8", "V7")])
        # V4 was merged into V3
        self.assertEqual((self.dir / "viruses_dropped_reps.txt").read_text(), "V4\n")

    def test_gzip_inputs(self):
        self.run_main(UPDATED_CSV, EVENTS, compress=True)
        self.assertEqual(self.tsv_rows("viruses_clusters_to_annotate.tsv"), [("V7", "V7"), ("V8", "V7")])

    def test_no_new_cluster(self):
        self.run_main(UPDATED_CSV.replace(",V7\n", ",V1\n"), EVENTS_HEADER + "V1\tgrown\tV1\t2\t3\n")
        self.assertEqual(self.tsv_rows("viruses_clusters_to_annotate.tsv"), [])

    def test_unknown_rep_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_main(UPDATED_CSV, EVENTS_HEADER + "V9\tnew\t\t0\t1\n")
        self.assertIn("not a cluster_rep", str(cm.exception))

    def test_unknown_event_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_main(UPDATED_CSV, EVENTS_HEADER + "V7\tsplit\t\t0\t2\n")
        self.assertIn("unknown event", str(cm.exception))

    def test_new_cluster_with_old_reps_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.run_main(UPDATED_CSV, EVENTS_HEADER + "V7\tnew\tV4\t0\t2\n")
        self.assertIn("lists old representatives", str(cm.exception))

    def test_dropped_rep_still_rep_fails(self):
        csv_text = UPDATED_CSV.replace("V4,V4,0.5,100,h4,V3", "V4,V4,0.5,100,h4,V4")
        with self.assertRaises(SystemExit) as cm:
            self.run_main(csv_text, EVENTS)
        self.assertIn("still cluster_reps", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
