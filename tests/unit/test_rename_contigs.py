#!/usr/bin/env python3
"""
Unit tests for rename_contigs.py
"""

from __future__ import annotations

import csv
import os
import sys
import tempfile
import unittest
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures" / "rename_contigs"
BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import rename_contigs as rc


class TestDefinePrefix(unittest.TestCase):
    def test_mgyv_zero_padded(self):
        self.assertEqual(rc.define_prefix("MGYV", 1), "MGYV0000000001")

    def test_mgyv_large_number(self):
        self.assertEqual(rc.define_prefix("MGYV", 1234567890), "MGYV1234567890")

    def test_generic_prefix(self):
        self.assertEqual(rc.define_prefix("seq", 5), "seq5")

    def test_generic_prefix_no_padding(self):
        self.assertEqual(rc.define_prefix("barley", 42), "barley42")


class TestParseAttrs(unittest.TestCase):
    def setUp(self):
        self.attrs_str = "ID=MGYG000535629_9|viral_sequence-1:3862;virify_quality=LC;checkv_quality=Low-quality"
        self.attrs, self.order = rc.parse_attrs(self.attrs_str)

    def test_id_parsed(self):
        self.assertEqual(self.attrs["ID"], "MGYG000535629_9|viral_sequence-1:3862")

    def test_other_keys_parsed(self):
        self.assertEqual(self.attrs["virify_quality"], "LC")
        self.assertEqual(self.attrs["checkv_quality"], "Low-quality")

    def test_order_preserved(self):
        self.assertEqual(self.order, ["ID", "virify_quality", "checkv_quality"])

    def test_trailing_semicolon_handled(self):
        attrs, order = rc.parse_attrs("ID=test;key=val;")
        self.assertEqual(attrs["key"], "val")
        self.assertEqual(len(order), 2)


class TestRenameFasta(unittest.TestCase):
    def setUp(self):
        self.fna = str(FIXTURES / "file1.fna")
        self.records, self.map_dir, self.map_records, self.count = rc.rename_fasta(
            self.fna, "seq", False, 1, None, biome="rhizosphere", source_type="genome"
        )

    def test_count_two_sequences(self):
        self.assertEqual(self.count, 2)

    def test_record_ids_assigned(self):
        self.assertEqual(self.records[0].id, "seq1")
        self.assertEqual(self.records[1].id, "seq2")

    def test_map_dir_keys(self):
        self.assertIn("MGYG000535629_9 viral_sequence|1:3862", self.map_dir)
        self.assertIn("MGYG000535629_15 viral_sequence|1:5155", self.map_dir)

    def test_map_dir_values(self):
        self.assertEqual(self.map_dir["MGYG000535629_9 viral_sequence|1:3862"], "seq1")
        self.assertEqual(self.map_dir["MGYG000535629_15 viral_sequence|1:5155"], "seq2")

    def test_map_records_columns(self):
        row = self.map_records[0]
        # [original, temporary, short, biome, type]
        self.assertEqual(len(row), 5)
        self.assertEqual(row[0], "MGYG000535629_9 viral_sequence|1:3862")
        self.assertEqual(row[1], "seq1")
        self.assertEqual(row[2], "MGYG000535629_9")
        self.assertEqual(row[3], "rhizosphere")
        self.assertEqual(row[4], "genome")

    def test_start_accession_offset(self):
        _, _, records, _ = rc.rename_fasta(
            self.fna, "seq", False, 10, None
        )
        self.assertEqual(records[0][1], "seq10")
        self.assertEqual(records[1][1], "seq11")

    def test_keep_viral_identifier(self):
        records, map_dir, _, _ = rc.rename_fasta(
            self.fna, "seq", True, 1, None
        )
        self.assertIn("|", records[0].id)
        self.assertTrue(records[0].id.startswith("seq1|"))


class TestRenameGff(unittest.TestCase):
    def setUp(self):
        _, self.map_dir, _, _ = rc.rename_fasta(
            str(FIXTURES / "file1.fna"), "seq", False, 1, None
        )
        self.tmp = tempfile.TemporaryDirectory()
        self.out_gff = str(Path(self.tmp.name) / "out.gff")
        self.seqs_count = rc.rename_gff(
            str(FIXTURES / "file1.gff"), self.out_gff, self.map_dir
        )

    def tearDown(self):
        self.tmp.cleanup()

    def test_returns_sequence_count(self):
        # 2 viral_sequence feature lines (one per contig), not CDS
        self.assertEqual(self.seqs_count, 2)

    def test_output_file_created(self):
        self.assertTrue(os.path.exists(self.out_gff))

    def test_ids_renamed(self):
        with open(self.out_gff) as f:
            content = f.read()
        self.assertIn("ID=seq1", content)
        self.assertIn("ID=seq2", content)
        self.assertNotIn("MGYG000535629_9|viral_sequence-1:3862", content)

    def test_header_written(self):
        with open(self.out_gff) as f:
            first_line = f.readline()
        self.assertTrue(first_line.startswith("##gff-version"))

    def test_append_mode_skips_header(self):
        out2 = str(Path(self.tmp.name) / "appended.gff")
        rc.rename_gff(str(FIXTURES / "file1.gff"), out2, self.map_dir, append=False)
        rc.rename_gff(str(FIXTURES / "file1.gff"), out2, self.map_dir, append=True)
        with open(out2) as f:
            headers = [l for l in f if l.startswith("##gff-version")]
        self.assertEqual(len(headers), 1)


class TestWriteMapfile(unittest.TestCase):
    def test_header_and_rows(self):
        data = [
            ["MGYG000535629_9 viral_sequence|1:3862", "seq1", "MGYG000535629_9", "rhizosphere", "genome"],
            ["MGYG000535629_15 viral_sequence|1:5155", "seq2", "MGYG000535629_15", "rhizosphere", "genome"],
        ]
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            path = f.name
        try:
            rc.write_mapfile(path, data)
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                rows = list(reader)
            self.assertEqual(rows[0], ["original", "temporary", "short", "biome", "type"])
            self.assertEqual(rows[1][1], "seq1")
            self.assertEqual(rows[2][3], "rhizosphere")
            self.assertEqual(len(rows), 3)
        finally:
            os.unlink(path)


class TestMainNoCombine(unittest.TestCase):
    """Integration test: single input FASTA + GFF, no --combine."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, extra_args=None):
        orig_dir = os.getcwd()
        os.chdir(self.tmp.name)
        try:
            sys.argv = [
                "rename_contigs.py",
                "--fasta",   str(FIXTURES / "file1.fna"),
                "--gff",     str(FIXTURES / "file1.gff"),
                "--prefix",  "seq",
                "--biome",   "rhizosphere",
                "--type",    "genome",
            ] + (extra_args or [])
            rc.main()
        finally:
            os.chdir(orig_dir)

    def test_renamed_fasta_exists(self):
        self._run()
        self.assertTrue(os.path.exists(str(Path(self.tmp.name) / "renamed_file1.fna")))

    def test_renamed_fasta_has_two_records(self):
        self._run()
        with open(str(Path(self.tmp.name) / "renamed_file1.fna")) as f:
            headers = [l for l in f if l.startswith(">")]
        self.assertEqual(len(headers), 2)

    def test_map_file_exists(self):
        self._run()
        self.assertTrue(os.path.exists(str(Path(self.tmp.name) / "file1_map.txt")))

    def test_map_file_has_biome_and_type(self):
        self._run()
        with open(str(Path(self.tmp.name) / "file1_map.txt")) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)
        self.assertIn("biome", rows[0])
        self.assertIn("type", rows[0])
        self.assertEqual(rows[1][3], "rhizosphere")
        self.assertEqual(rows[1][4], "genome")

    def test_renamed_gff_exists(self):
        self._run()
        self.assertTrue(os.path.exists(str(Path(self.tmp.name) / "file1_renamed.gff")))

    def test_renamed_gff_ids_updated(self):
        self._run()
        with open(str(Path(self.tmp.name) / "file1_renamed.gff")) as f:
            content = f.read()
        self.assertIn("ID=seq1", content)
        self.assertIn("ID=seq2", content)


class TestMainCombine(unittest.TestCase):
    """Integration test: two input FASTAs + GFFs with --combine."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self):
        orig_dir = os.getcwd()
        os.chdir(self.tmp.name)
        try:
            sys.argv = [
                "rename_contigs.py",
                "--fasta",   str(FIXTURES / "file1.fna"), str(FIXTURES / "file2.fna"),
                "--gff",     str(FIXTURES / "file1.gff"), str(FIXTURES / "file2.gff"),
                "--prefix",  "seq",
                "--biome",   "rhizosphere", "soil",
                "--type",    "genome", "genome",
                "--map",     "combined.tsv",
                "--combine",
            ]
            rc.main()
        finally:
            os.chdir(orig_dir)

    def test_combined_fasta_exists(self):
        self._run()
        self.assertTrue(os.path.exists(str(Path(self.tmp.name) / "combined.fna")))

    def test_combined_fasta_has_four_records(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.fna")) as f:
            headers = [l for l in f if l.startswith(">")]
        self.assertEqual(len(headers), 4)

    def test_combined_fasta_ids_sequential(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.fna")) as f:
            headers = [l.strip().lstrip(">") for l in f if l.startswith(">")]
        self.assertEqual(headers, ["seq1", "seq2", "seq3", "seq4"])

    def test_combined_map_exists(self):
        self._run()
        self.assertTrue(os.path.exists(str(Path(self.tmp.name) / "combined.tsv")))

    def test_combined_map_has_four_entries(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.tsv")) as f:
            rows = f.readlines()
        # 1 header + 4 data rows
        self.assertEqual(len(rows), 5)

    def test_combined_map_biomes_per_file(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.tsv")) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)
        biome_idx = rows[0].index("biome")
        biomes = [r[biome_idx] for r in rows[1:]]
        self.assertEqual(biomes[:2], ["rhizosphere", "rhizosphere"])
        self.assertEqual(biomes[2:], ["soil", "soil"])

    def test_combined_gff_exists(self):
        self._run()
        self.assertTrue(os.path.exists(str(Path(self.tmp.name) / "combined.gff")))

    def test_combined_gff_has_one_header(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.gff")) as f:
            headers = [l for l in f if l.startswith("##gff-version")]
        self.assertEqual(len(headers), 1)

    def test_combined_gff_ids_renamed(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.gff")) as f:
            content = f.read()
        for i in range(1, 5):
            self.assertIn(f"ID=seq{i}", content)

    def test_no_per_file_outputs_in_combine_mode(self):
        self._run()
        self.assertFalse(os.path.exists(str(Path(self.tmp.name) / "renamed_file1.fna")))
        self.assertFalse(os.path.exists(str(Path(self.tmp.name) / "renamed_file2.fna")))


if __name__ == "__main__":
    unittest.main()
