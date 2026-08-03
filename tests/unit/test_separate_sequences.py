#!/usr/bin/env python3
"""
Unit tests for separate_sequences.py
"""

import sys
import tempfile
import unittest
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures" / "separate_sequences"
BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import separate_sequences as ss


class TestReadDefinitions(unittest.TestCase):
    def setUp(self):
        self.definitions = ss.read_definitions(str(FIXTURES / "combined_map.tsv"))

    def test_mgnify_definition_unchanged(self):
        self.assertEqual(self.definitions["seq1"], "virus")
        self.assertEqual(self.definitions["seq3"], "plasmid")
        self.assertEqual(self.definitions["seq4"], "prophage")

    def test_third_party_prefix_stripped(self):
        self.assertEqual(self.definitions["seq2"], "prophage")

    def test_missing_key_has_no_default_category(self):
        self.assertNotIn("does-not-exist", self.definitions)


class TestSplitFasta(unittest.TestCase):
    def setUp(self):
        self.definitions = ss.read_definitions(str(FIXTURES / "combined_map.tsv"))
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _split(self, category):
        out_fna = str(Path(self.tmp.name) / f"{category}.fna")
        kept_ids = ss.split_fasta(str(FIXTURES / "combined.fna"), self.definitions, category, out_fna)
        return kept_ids, out_fna

    def test_virus_bucket(self):
        kept_ids, out_fna = self._split("virus")
        self.assertEqual(kept_ids, {"seq1"})
        with open(out_fna) as f:
            self.assertIn(">seq1", f.read())

    def test_prophage_bucket_includes_third_party(self):
        # seq2's definition is 'third_party_prophage' and seq4's is 'prophage';
        # both must land in the 'prophage' bucket
        kept_ids, _ = self._split("prophage")
        self.assertEqual(kept_ids, {"seq2", "seq4"})

    def test_plasmid_bucket(self):
        kept_ids, _ = self._split("plasmid")
        self.assertEqual(kept_ids, {"seq3"})


class TestSplitGff(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_writes_only_kept_sequences_and_returns_protein_ids(self):
        out_gff = str(Path(self.tmp.name) / "out.gff")
        protein_ids = ss.split_gff(str(FIXTURES / "combined.gff"), {"seq1"}, out_gff)

        self.assertEqual(protein_ids, {"seq1_00001", "seq1_00002"})
        with open(out_gff) as f:
            content = f.read()
        self.assertIn("ID=seq1", content)
        self.assertNotIn("ID=seq2", content)
        self.assertTrue(content.startswith("##gff-version 3\n"))

    def test_multiple_kept_sequences(self):
        out_gff = str(Path(self.tmp.name) / "out.gff")
        protein_ids = ss.split_gff(str(FIXTURES / "combined.gff"), {"seq2", "seq4"}, out_gff)
        self.assertEqual(protein_ids, {"seq2_00001", "seq4_00001"})

    def test_no_kept_sequences_returns_empty(self):
        out_gff = str(Path(self.tmp.name) / "out.gff")
        protein_ids = ss.split_gff(str(FIXTURES / "combined.gff"), set(), out_gff)
        self.assertEqual(protein_ids, set())
        with open(out_gff) as f:
            self.assertEqual(f.read(), "##gff-version 3\n")


class TestSplitFaa(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_keeps_only_requested_proteins(self):
        out_faa = str(Path(self.tmp.name) / "out.faa")
        ss.split_faa(str(FIXTURES / "combined.faa"), {"seq1_00001", "seq1_00002"}, out_faa)
        with open(out_faa) as f:
            content = f.read()
        self.assertIn(">seq1_00001", content)
        self.assertIn(">seq1_00002", content)
        self.assertNotIn(">seq2_00001", content)
        self.assertNotIn(">orphan_protein_not_in_gff", content)


class TestMainIntegration(unittest.TestCase):
    """Run main() end-to-end for one category and check all three outputs."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.prefix = str(Path(self.tmp.name) / "prophage")

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, category):
        sys.argv = [
            "separate_sequences.py",
            "--fna", str(FIXTURES / "combined.fna"),
            "--gff", str(FIXTURES / "combined.gff"),
            "--faa", str(FIXTURES / "combined.faa"),
            "--map", str(FIXTURES / "combined_map.tsv"),
            "--category", category,
            "--output-prefix", self.prefix,
        ]
        ss.main()

    def test_prophage_category_end_to_end(self):
        self._run("prophage")
        with open(self.prefix + ".fna") as f:
            fna_ids = {l.strip().lstrip(">") for l in f if l.startswith(">")}
        self.assertEqual(fna_ids, {"seq2", "seq4"})

        with open(self.prefix + ".gff") as f:
            gff_content = f.read()
        self.assertIn("ID=seq2", gff_content)
        self.assertIn("ID=seq4", gff_content)
        self.assertNotIn("ID=seq1", gff_content)

        with open(self.prefix + ".faa") as f:
            faa_content = f.read()
        self.assertIn(">seq2_00001", faa_content)
        self.assertIn(">seq4_00001", faa_content)
        self.assertNotIn(">seq1_00001", faa_content)
        self.assertNotIn(">orphan_protein_not_in_gff", faa_content)

    def test_gff_and_faa_are_optional(self):
        sys.argv = [
            "separate_sequences.py",
            "--fna", str(FIXTURES / "combined.fna"),
            "--map", str(FIXTURES / "combined_map.tsv"),
            "--category", "virus",
            "--output-prefix", self.prefix,
        ]
        ss.main()
        self.assertTrue(Path(self.prefix + ".fna").exists())
        self.assertFalse(Path(self.prefix + ".gff").exists())
        self.assertFalse(Path(self.prefix + ".faa").exists())


if __name__ == "__main__":
    unittest.main()
