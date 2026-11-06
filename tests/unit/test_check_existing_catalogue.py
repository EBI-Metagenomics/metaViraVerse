#!/usr/bin/env python3
"""
Unit tests for check_existing_catalogue.py, on manifests written by
write_catalogue_manifest.py (the round trip of generate -> update).
"""

import gzip
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import check_existing_catalogue as cec
import write_catalogue_manifest as wcm

MAP = (
    "original\ttemporary\tshort\tbiome\ttype\tdefinition\n"
    "a\tMGYV0000000001\ta\troot\tgenome\tvirus\n"
    "b\tMGYV0000000002\tb\troot\tgenome\tprophage\n"
)


class TestRoundTrip(unittest.TestCase):
    """Write a manifest as WRITE_CATALOGUE_MANIFEST does, then check it as CHECK_EXISTING_CATALOGUE does."""

    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp = tempfile.TemporaryDirectory()
        os.chdir(self._tmp.name)
        Path("rename_map.tsv").write_text(MAP)
        Path("vjdb1_merged_reps.csv").write_text("vjdb_id,copy_of,seq_gc,seq_len,seq_hash,cluster_rep\n")
        with gzip.open("vjdb1_merged_reps.fna.gz", "wt") as fh:
            fh.write(">MGYV0000000001\nACGT\n")
        Path("h1.tbl").write_text("hit\n")
        Path("h2.tbl").write_text("hit2\n")

    def tearDown(self):
        os.chdir(self._cwd)
        self._tmp.cleanup()

    def write_manifest(self, clustering="vjdb"):
        wcm.main([
            "--rename-accession", "MGYV", "--viral-clustering", clustering,
            "--entry", "file", "rename_map", "additional_data/rename_mapfile.tsv", "rename_map.tsv",
            "--entry", "file", "viruses.vjdb_clusters_csv", "viruses/clustering/viruses_vjdb_clusters.csv", "vjdb1_merged_reps.csv",
            "--entry", "file", "viruses.vjdb_reps_fna", "viruses/clustering/viruses_vjdb_reps.fna.gz", "vjdb1_merged_reps.fna.gz",
            "--entry", "list", "viruses.hmmer_tables", "viruses/hmmer/", "h1.tbl",
            "--entry", "list", "viruses.hmmer_tables", "viruses/hmmer/", "h2.tbl",
        ])

    def publish(self):
        """Lay the files out as the next task sees them: under their published names."""
        os.rename("rename_map.tsv", "rename_mapfile.tsv")
        os.rename("vjdb1_merged_reps.csv", "viruses_vjdb_clusters.csv")
        os.rename("vjdb1_merged_reps.fna.gz", "viruses_vjdb_reps.fna.gz")

    def check(self, rename_accession="MGYV"):
        cec.main(["--manifest", "catalogue_manifest.json", "--rename-accession", rename_accession])
        with open("existing_catalogue_check.json") as fh:
            return json.load(fh)

    def test_unchanged_catalogue_passes(self):
        self.write_manifest()
        self.publish()
        report = self.check()
        self.assertEqual(report["status"], "ok")
        self.assertEqual(report["files_checked"], 5)
        self.assertEqual(report["accessions"]["last"], 2)
        self.assertEqual(report["manifest_sha256"], cec.sha256sum("catalogue_manifest.json"))

    def test_changed_file_fails(self):
        self.write_manifest()
        self.publish()
        Path("viruses_vjdb_clusters.csv").write_text("vjdb_id,copy_of,seq_gc,seq_len,seq_hash,cluster_rep\nX,X,0,1,h,X\n")
        with self.assertRaises(SystemExit):
            self.check()

    def test_same_size_different_content_fails(self):
        self.write_manifest()
        self.publish()
        Path("h1.tbl").write_text("HIT\n")
        with self.assertRaises(SystemExit):
            self.check()

    def test_missing_file_fails(self):
        self.write_manifest()
        self.publish()
        os.remove("h2.tbl")
        with self.assertRaises(SystemExit):
            self.check()

    def test_not_update_ready_fails(self):
        self.write_manifest(clustering="default")
        self.publish()
        with self.assertRaises(SystemExit):
            self.check()

    def test_other_accession_prefix_fails(self):
        self.write_manifest()
        self.publish()
        with self.assertRaises(SystemExit):
            self.check(rename_accession="seq")


class TestCheckManifest(unittest.TestCase):
    def manifest(self, **changes):
        m = {
            "manifest_version": 1,
            "viral_clustering": {"update_ready": True, "not_update_ready_reason": None},
            "accessions": {"prefix": "MGYV", "last": 2},
            "files": {"rename_map": {"path": "a.tsv", "size": 1, "sha256": "x"}},
        }
        m.update(changes)
        return m

    def test_valid(self):
        self.assertEqual(cec.check_manifest(self.manifest(), "MGYV"), [])

    def test_unsupported_version(self):
        problems = cec.check_manifest(self.manifest(manifest_version=99), "MGYV")
        self.assertTrue(any("manifest_version" in p for p in problems))

    def test_no_files(self):
        problems = cec.check_manifest(self.manifest(files={}), "MGYV")
        self.assertTrue(any("no files" in p for p in problems))

    def test_all_problems_listed(self):
        m = self.manifest(manifest_version=99, viral_clustering={"update_ready": False, "not_update_ready_reason": "r"})
        self.assertEqual(len(cec.check_manifest(m, "seq")), 3)


class TestCheckFiles(unittest.TestCase):
    def test_repeated_file_name_reported(self):
        files = {
            "a": {"path": "x/same.tsv", "size": 0, "sha256": ""},
            "b": {"path": "y/same.tsv", "size": 0, "sha256": ""},
        }
        problems, _checked = cec.check_files(files)
        self.assertTrue(any("same file name" in p for p in problems))


if __name__ == "__main__":
    unittest.main()
