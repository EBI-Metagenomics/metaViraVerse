#!/usr/bin/env python3
"""
Unit tests for write_catalogue_manifest.py
"""

import gzip
import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import write_catalogue_manifest as wcm

MAP_HEADER = "original\ttemporary\tshort\tbiome\ttype\tdefinition\n"

VERSIONS_YML = """\
"CHECKV_ENDTOEND":
    checkv: 1.0.3
"VCLUST_CLUSTER":
    vclust: 1.3.1
Workflow:
    EBI-Metagenomics/metaviraverse: v1.0.0dev
    Nextflow: 26.10.4
"""


def write(path, text, compress=False):
    opener = gzip.open if compress else open
    with opener(path, "wt") as fh:
        fh.write(text)


def map_rows(ids):
    return MAP_HEADER + "".join(f"orig{i}\t{acc}\torig{i}\troot\tgenome\tvirus\n" for i, acc in enumerate(ids))


class InTempDir(unittest.TestCase):
    """Run each test in an empty directory, like a Nextflow task."""

    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp = tempfile.TemporaryDirectory()
        os.chdir(self._tmp.name)

    def tearDown(self):
        os.chdir(self._cwd)
        self._tmp.cleanup()


class TestTargetPath(unittest.TestCase):
    def test_directory_target_keeps_name(self):
        self.assertEqual(wcm.target_path("viruses/clustering/", "a.csv"), "viruses/clustering/a.csv")

    def test_current_directory_target(self):
        self.assertEqual(wcm.target_path("./", "catalogue.json"), "catalogue.json")

    def test_file_target_renames(self):
        self.assertEqual(wcm.target_path("viruses/x.gff.gz", "a_final.gff.gz"), "viruses/x.gff.gz")


class TestAccessionRange(InTempDir):
    def test_mgyv_zero_padded(self):
        write("map.tsv", map_rows(["MGYV0000000003", "MGYV0000000001", "MGYV0000000012"]))
        acc = wcm.accession_range("map.tsv", "MGYV")
        self.assertEqual(acc["first"], 1)
        self.assertEqual(acc["last"], 12)
        self.assertEqual(acc["count"], 3)
        self.assertEqual(acc["zero_padded_digits"], 10)

    def test_unpadded_prefix(self):
        write("map.tsv", map_rows(["seq1", "seq2", "seq10"]))
        acc = wcm.accession_range("map.tsv", "seq")
        self.assertEqual(acc["last"], 10)
        self.assertIsNone(acc["zero_padded_digits"])

    def test_gzip_map(self):
        write("map.tsv.gz", map_rows(["seq5"]), compress=True)
        self.assertEqual(wcm.accession_range("map.tsv.gz", "seq")["last"], 5)

    def test_repeated_header_skipped(self):
        write("map.tsv", map_rows(["seq1"]) + map_rows(["seq2"]))
        self.assertEqual(wcm.accession_range("map.tsv", "seq")["count"], 2)

    def test_foreign_prefix_fails(self):
        write("map.tsv", map_rows(["seq1", "MGYV0000000002"]))
        with self.assertRaises(SystemExit):
            wcm.accession_range("map.tsv", "seq")

    def test_repeated_accession_fails(self):
        write("map.tsv", map_rows(["seq1", "seq1"]))
        with self.assertRaises(SystemExit):
            wcm.accession_range("map.tsv", "seq")

    def test_empty_map(self):
        write("map.tsv", MAP_HEADER)
        acc = wcm.accession_range("map.tsv", "seq")
        self.assertEqual(acc["count"], 0)
        self.assertIsNone(acc["last"])


class TestReadVersions(InTempDir):
    def test_processes_and_tools(self):
        write("versions.yml", VERSIONS_YML)
        versions = wcm.read_versions("versions.yml")
        self.assertEqual(versions["CHECKV_ENDTOEND"], {"checkv": "1.0.3"})
        self.assertEqual(versions["Workflow"]["Nextflow"], "26.10.4")


class TestCollectFiles(InTempDir):
    def test_file_and_list_entries(self):
        write("a.csv", "x\n")
        write("h1.tbl", "1\n")
        write("h2.tbl", "2\n")
        files = wcm.collect_files([
            ["file", "viruses.vjdb_clusters_csv", "viruses/clustering/clusters.csv", "a.csv"],
            ["list", "viruses.hmmer_tables", "viruses/hmmer/", "h2.tbl"],
            ["list", "viruses.hmmer_tables", "viruses/hmmer/", "h1.tbl"],
        ])
        entry = files["viruses.vjdb_clusters_csv"]
        self.assertEqual(entry["path"], "viruses/clustering/clusters.csv")
        self.assertEqual(entry["size"], 2)
        self.assertEqual(entry["sha256"], hashlib.sha256(b"x\n").hexdigest())
        self.assertEqual([f["path"] for f in files["viruses.hmmer_tables"]], ["viruses/hmmer/h1.tbl", "viruses/hmmer/h2.tbl"])

    def test_duplicate_key_fails(self):
        write("a.csv", "x\n")
        with self.assertRaises(SystemExit):
            wcm.collect_files([["file", "k", "a/", "a.csv"], ["file", "k", "b/", "a.csv"]])

    def test_same_publish_path_fails(self):
        write("a.csv", "x\n")
        write("b.csv", "y\n")
        with self.assertRaises(SystemExit):
            wcm.collect_files([["file", "k1", "out/c.csv", "a.csv"], ["file", "k2", "out/c.csv", "b.csv"]])

    def test_missing_file_fails(self):
        with self.assertRaises(SystemExit):
            wcm.collect_files([["file", "k", "a/", "nope.csv"]])

    def test_unknown_kind_fails(self):
        write("a.csv", "x\n")
        with self.assertRaises(SystemExit):
            wcm.collect_files([["dir", "k", "a/", "a.csv"]])


class TestMain(InTempDir):
    def setUp(self):
        super().setUp()
        write("rename_map.tsv", map_rows(["MGYV0000000001", "MGYV0000000002"]))
        write("versions.yml", VERSIONS_YML)
        write("vjdb1_merged_reps.csv", "vjdb_id,copy_of,seq_gc,seq_len,seq_hash,cluster_rep\n")
        write("vjdb1_merged_reps.fna.gz", ">MGYV0000000001\nACGT\n", compress=True)

    def run_main(self, clustering, with_vjdb=True):
        argv = [
            "--rename-accession", "MGYV", "--viral-clustering", clustering, "--pipeline-version", "1.0.0dev",
            "--entry", "file", "rename_map", "additional_data/rename_mapfile.tsv", "rename_map.tsv",
            "--entry", "file", "software_versions", "pipeline_info/", "versions.yml",
        ]
        if with_vjdb:
            argv += [
                "--entry", "file", "viruses.vjdb_clusters_csv", "viruses/clustering/viruses_vjdb_clusters.csv", "vjdb1_merged_reps.csv",
                "--entry", "file", "viruses.vjdb_reps_fna", "viruses/clustering/viruses_vjdb_reps.fna.gz", "vjdb1_merged_reps.fna.gz",
            ]
        wcm.main(argv + getattr(self, "extra", []))
        with open("catalogue_manifest.json") as fh:
            return json.load(fh)

    def test_vjdb_catalogue_is_update_ready(self):
        manifest = self.run_main("vjdb")
        self.assertEqual(manifest["manifest_version"], wcm.MANIFEST_VERSION)
        self.assertTrue(manifest["viral_clustering"]["update_ready"])
        self.assertIsNone(manifest["viral_clustering"]["not_update_ready_reason"])
        self.assertEqual(manifest["accessions"]["last"], 2)
        self.assertEqual(manifest["deduplication"]["hash"], "sha256")
        self.assertEqual(manifest["files"]["rename_map"]["path"], "additional_data/rename_mapfile.tsv")
        self.assertEqual(manifest["files"]["software_versions"]["path"], "pipeline_info/versions.yml")
        self.assertEqual(manifest["software_versions"]["VCLUST_CLUSTER"]["vclust"], "1.3.1")

    def test_default_clustering_is_not_update_ready(self):
        manifest = self.run_main("default", with_vjdb=False)
        self.assertFalse(manifest["viral_clustering"]["update_ready"])
        self.assertIn("vjdb", manifest["viral_clustering"]["not_update_ready_reason"])

    def test_vjdb_without_tables_is_not_update_ready(self):
        manifest = self.run_main("vjdb", with_vjdb=False)
        self.assertFalse(manifest["viral_clustering"]["update_ready"])
        self.assertIn("viruses.vjdb_clusters_csv", manifest["viral_clustering"]["not_update_ready_reason"])

    def test_generate_has_no_parent(self):
        manifest = self.run_main("vjdb")
        self.assertEqual(manifest["catalogue"]["mode"], "generate")
        self.assertIsNone(manifest["catalogue"]["parent"])

    def test_update_records_parent(self):
        Path("existing_catalogue_check.json").write_text(json.dumps({
            "status": "ok", "manifest_sha256": "abc",
            "catalogue": {"created": "2026-01-01T00:00:00+00:00", "mode": "generate", "pipeline_version": "1.0.0"},
        }))
        self.extra = ["--mode", "update", "--parent-check", "existing_catalogue_check.json"]
        manifest = self.run_main("vjdb")
        self.assertEqual(manifest["catalogue"]["mode"], "update")
        self.assertEqual(manifest["catalogue"]["parent"]["manifest_sha256"], "abc")
        self.assertEqual(manifest["catalogue"]["parent"]["pipeline_version"], "1.0.0")

    def test_update_without_parent_fails(self):
        self.extra = ["--mode", "update"]
        with self.assertRaises(SystemExit):
            self.run_main("vjdb")

    def test_rename_map_required(self):
        with self.assertRaises(SystemExit):
            wcm.main(["--rename-accession", "MGYV", "--viral-clustering", "vjdb"])


if __name__ == "__main__":
    unittest.main()
