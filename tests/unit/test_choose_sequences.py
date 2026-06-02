#!/usr/bin/env python3
"""
Unit tests for choose_sequences.py
"""

import sys
import tempfile
import unittest
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures" / "choose_sequences"
BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import choose_sequences as cs


class TestSeqHash(unittest.TestCase):
    def test_consistent(self):
        self.assertEqual(cs.seq_hash("ACGT"), cs.seq_hash("ACGT"))

    def test_case_insensitive(self):
        self.assertEqual(cs.seq_hash("acgt"), cs.seq_hash("ACGT"))

    def test_different_sequences_differ(self):
        self.assertNotEqual(cs.seq_hash("ACGT"), cs.seq_hash("TGCA"))


class TestReadMap(unittest.TestCase):
    def setUp(self):
        self.mapping = cs.read_map([str(FIXTURES / "barley10.map.tsv")])

    def test_loads_ten_entries(self):
        self.assertEqual(len(self.mapping), 10)

    def test_temporary_to_original(self):
        # barley1 -> MGYG000535629_9 viral_sequence|1:3862
        self.assertEqual(self.mapping["barley1"]["original"], "MGYG000535629_9 viral_sequence|1:3862")

    def test_plasmid_entry_present(self):
        # barley3 maps to a plasmid original name
        self.assertIn("plasmid", self.mapping["barley3"]["original"])

    def test_biome_column_read(self):
        self.assertEqual(self.mapping["barley1"]["biome"], "rhizosphere")

    def test_type_column_read(self):
        self.assertEqual(self.mapping["barley1"]["type"], "genome")

    def test_missing_key_returns_default(self):
        self.assertEqual(self.mapping.get("notakey", "notakey"), "notakey")


class TestParseRnaGff(unittest.TestCase):
    def setUp(self):
        self.rna_seqs = cs.parse_rna_gff([str(FIXTURES / "barley10.gff")])

    def test_finds_barley1(self):
        self.assertIn("barley1", self.rna_seqs)

    def test_only_one_entry(self):
        self.assertEqual(len(self.rna_seqs), 1)

    def test_empty_list_returns_empty_set(self):
        self.assertEqual(cs.parse_rna_gff([]), set())


class TestReadQuality(unittest.TestCase):
    def setUp(self):
        self.quality = cs.read_quality([str(FIXTURES / "barley10_quality.tsv")])

    def test_loads_ten_entries(self):
        self.assertEqual(len(self.quality), 10)

    def test_low_quality_entry(self):
        self.assertEqual(self.quality["barley1"]["checkv_quality"], "Low-quality")

    def test_not_determined_entry(self):
        self.assertEqual(self.quality["barley9"]["checkv_quality"], "Not-determined")

    def test_contig_length_excluded(self):
        self.assertNotIn("contig_length", self.quality["barley1"])

    def test_all_quality_columns_present(self):
        for col in cs.QUALITY_COLUMNS:
            self.assertIn(col, self.quality["barley1"])


class TestPassesFilter(unittest.TestCase):
    def _entry(self, rrna, viral_type, checkv_quality=None):
        q = {"checkv_quality": checkv_quality} if checkv_quality is not None else None
        return {"rrna": rrna, "viral_type": viral_type, "quality": q}

    def test_clean_virus_passes(self):
        self.assertTrue(cs.passes_filter(self._entry("No", "virus", "Low-quality")))

    def test_rrna_virus_filtered(self):
        self.assertFalse(cs.passes_filter(self._entry("Yes", "virus", "Low-quality")))

    def test_rrna_plasmid_passes(self):
        # Plasmids are exempt from the rRNA filter
        self.assertTrue(cs.passes_filter(self._entry("Yes", "plasmid", "Low-quality")))

    def test_not_determined_virus_filtered(self):
        self.assertFalse(cs.passes_filter(self._entry("No", "virus", "Not-determined")))

    def test_not_determined_plasmid_passes(self):
        # Plasmids are exempt from the quality filter too
        self.assertTrue(cs.passes_filter(self._entry("No", "plasmid", "Not-determined")))

    def test_no_quality_data_passes(self):
        # Entry with no CheckV data should not be filtered on quality
        self.assertTrue(cs.passes_filter(self._entry("No", "virus", None)))

    def test_rrna_not_provided_not_filtered(self):
        # 'not-provided' means rRNA tool wasn't run — do not filter
        self.assertTrue(cs.passes_filter(self._entry("not-provided", "virus", "Low-quality")))


class TestMainIntegration(unittest.TestCase):
    """Run main() end-to-end with fixture files and verify output counts."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out_fna = str(Path(self.tmp.name) / "chosen.fna")
        self.out_gff = str(Path(self.tmp.name) / "chosen.gff")
        self.out_tsv = str(Path(self.tmp.name) / "metadata.tsv")

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, extra_args=None):
        argv = [
            "choose_sequences.py",
            "--fna",        str(FIXTURES / "barley10.fasta"),
            "--gff",        str(FIXTURES / "barley10_viral.gff"),
            "--map",        str(FIXTURES / "barley10.map.tsv"),
            "--rrna",       str(FIXTURES / "barley10.gff"),
            "--quality",    str(FIXTURES / "barley10_quality.tsv"),
            "--output-fna", self.out_fna,
            "--output-gff", self.out_gff,
            "--output-tsv", self.out_tsv,
        ]
        if extra_args:
            argv += extra_args
        sys.argv = argv
        cs.main()

    def test_all_sequences_written(self):
        self._run()
        with open(self.out_tsv) as f:
            lines = f.readlines()
        # 1 header + 10 data rows
        self.assertEqual(len(lines), 11)

    def test_fna_has_ten_records(self):
        self._run()
        with open(self.out_fna) as f:
            headers = [l for l in f if l.startswith(">")]
        self.assertEqual(len(headers), 10)

    def test_filtered_tsv_excludes_rrna_and_not_determined(self):
        self._run()
        filtered_tsv = Path(self.tmp.name) / "filtered_metadata.tsv"
        with open(filtered_tsv) as f:
            lines = f.readlines()
        # barley1 (rRNA) and barley9 (Not-determined, virus) are removed
        # barley3 and barley4 (Not-determined, plasmid) are kept → 8 pass
        self.assertEqual(len(lines), 9)  # 1 header + 8 data rows

    def test_filtered_fna_has_eight_records(self):
        self._run()
        filtered_fna = Path(self.tmp.name) / "filtered_chosen.fna"
        with open(filtered_fna) as f:
            headers = [l for l in f if l.startswith(">")]
        self.assertEqual(len(headers), 8)

    def test_filtered_gff_has_eight_sequences(self):
        self._run()
        with open(self.out_gff) as f:
            seq_lines = [l for l in f if '\tviral_sequence\t' in l or '\tplasmid\t' in l]
        self.assertEqual(len(seq_lines), 8)

    def test_tsv_has_correct_columns(self):
        self._run()
        with open(self.out_tsv) as f:
            header = f.readline().rstrip("\n").split("\t")
        expected_start = ["sequence_id", "original_name", "description", "type",
                          "source_of_prediction", "biomes", "sequence_length", "rrna", "sequence_sha256"]
        self.assertEqual(header[:9], expected_start)
        for col in cs.QUALITY_COLUMNS:
            self.assertIn(col, header)

    def test_source_of_prediction_populated(self):
        self._run()
        with open(self.out_tsv) as f:
            rows = [l.rstrip("\n").split("\t") for l in f]
        header = rows[0]
        src_idx = header.index("source_of_prediction")
        seq_id_idx = header.index("sequence_id")
        by_id = {r[seq_id_idx]: r[src_idx] for r in rows[1:]}
        # barley1 and barley2 are VIRify predictions; barley3 is geNomad
        self.assertEqual(by_id["barley1"], "VIRify")
        self.assertEqual(by_id["barley2"], "VIRify")
        self.assertEqual(by_id["barley3"], "geNomad")

    def test_rrna_column_populated(self):
        self._run()
        with open(self.out_tsv) as f:
            rows = [l.rstrip("\n").split("\t") for l in f]
        header = rows[0]
        rrna_idx = header.index("rrna")
        seq_id_idx = header.index("sequence_id")
        by_id = {r[seq_id_idx]: r[rrna_idx] for r in rows[1:]}
        self.assertEqual(by_id["barley1"], "Yes")
        self.assertEqual(by_id["barley2"], "No")

    def test_biome_from_map_file(self):
        self._run()
        with open(self.out_tsv) as f:
            rows = [l.rstrip("\n").split("\t") for l in f]
        header = rows[0]
        biome_idx = header.index("biomes")
        self.assertEqual(rows[1][biome_idx], "rhizosphere")

    def test_type_from_map_file(self):
        self._run()
        with open(self.out_tsv) as f:
            rows = [l.rstrip("\n").split("\t") for l in f]
        header = rows[0]
        type_idx = header.index("type")
        self.assertEqual(rows[1][type_idx], "genome")


if __name__ == "__main__":
    unittest.main()
