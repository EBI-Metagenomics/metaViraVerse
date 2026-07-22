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


class TestFilterReason(unittest.TestCase):
    def _entry(self, rrna, viral_type, checkv_quality=None, viral_genes="0", kmer_freq="1.0"):
        q = {"checkv_quality": checkv_quality, "viral_genes": viral_genes, "kmer_freq": kmer_freq} \
            if checkv_quality is not None else None
        return {"rrna": rrna, "viral_type": viral_type, "quality": q}

    def test_clean_virus_passes(self):
        self.assertIsNone(cs.filter_reason(self._entry("No", "virus", "Low-quality")))

    def test_rrna_virus_reason(self):
        self.assertEqual(cs.filter_reason(self._entry("Yes", "virus", "Low-quality")), "rrna")

    def test_rrna_plasmid_passes(self):
        self.assertIsNone(cs.filter_reason(self._entry("Yes", "plasmid", "Low-quality")))

    def test_none_virus_reason(self):
        self.assertIsNone(cs.filter_reason(self._entry("No", "virus", "Not-determined", viral_genes="3", kmer_freq="1.0")))

    def test_not_determined_virus_no_viral_genes_passes(self):
        self.assertEqual(
            cs.filter_reason(self._entry("No", "virus", "Not-determined", viral_genes="0", kmer_freq="1.0")),
            "not_determined"
        )

    def test_not_determined_virus_high_kmer_passes(self):
        self.assertIsNone(cs.filter_reason(self._entry("No", "virus", "Not-determined", viral_genes="3", kmer_freq="1.5")))

    def test_not_determined_plasmid_passes(self):
        self.assertIsNone(cs.filter_reason(self._entry("No", "plasmid", "Not-determined", viral_genes="3", kmer_freq="1.0")))

    def test_no_quality_data_passes(self):
        self.assertIsNone(cs.filter_reason(self._entry("No", "virus", None)))

    def test_rrna_not_provided_not_filtered(self):
        self.assertIsNone(cs.filter_reason(self._entry("not-provided", "virus", "Low-quality")))


class TestPassesFilter(unittest.TestCase):
    def _entry(self, rrna, viral_type, checkv_quality=None, viral_genes="0", kmer_freq="1.0"):
        q = {"checkv_quality": checkv_quality, "viral_genes": viral_genes, "kmer_freq": kmer_freq} \
            if checkv_quality is not None else None
        return {"rrna": rrna, "viral_type": viral_type, "quality": q}

    def test_clean_virus_passes(self):
        self.assertTrue(cs.passes_filter(self._entry("No", "virus", "Low-quality")))

    def test_rrna_virus_filtered(self):
        self.assertFalse(cs.passes_filter(self._entry("Yes", "virus", "Low-quality")))

    def test_not_determined_virus_filtered(self):
        self.assertFalse(cs.passes_filter(self._entry("No", "virus", "Not-determined", viral_genes="0", kmer_freq="3.0")))

    def test_no_quality_data_passes(self):
        self.assertTrue(cs.passes_filter(self._entry("No", "virus", None)))


class TestMainIntegration(unittest.TestCase):
    """Run main() end-to-end with fixture files and verify output counts."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.prefix = str(Path(self.tmp.name) / "test")
        self.out_fna = self.prefix + "_filtered.fna"
        self.out_gff = self.prefix + "_filtered.gff"
        self.out_tsv = self.prefix + "_metadata.tsv"

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, extra_args=None):
        argv = [
            "choose_sequences.py",
            "--fna",           str(FIXTURES / "barley10.fasta"),
            "--gff",           str(FIXTURES / "barley10_viral.gff"),
            "--map",           str(FIXTURES / "barley10.map.tsv"),
            "--rrna",          str(FIXTURES / "barley10.gff"),
            "--quality",       str(FIXTURES / "barley10_quality.tsv"),
            "--output-prefix", self.prefix,
        ]
        if extra_args:
            argv += extra_args
        sys.argv = argv
        cs.main()

    def test_metadata_has_all_sequences(self):
        self._run()
        with open(self.out_tsv) as f:
            lines = f.readlines()
        # 1 header + 10 data rows
        self.assertEqual(len(lines), 11)

    def test_filtered_fna_has_nine_records(self):
        # barley1 (rrna) is excluded; barley9 passes the updated not_determined filter
        self._run()
        with open(self.out_fna) as f:
            headers = [l for l in f if l.startswith(">")]
        self.assertEqual(len(headers), 9)

    def test_filtered_tsv_excludes_rrna(self):
        self._run()
        filtered_tsv = self.prefix + "_filtered.tsv"
        with open(filtered_tsv) as f:
            lines = f.readlines()
        # barley1 (rRNA) excluded; 9 sequences pass
        self.assertEqual(len(lines), 10)  # 1 header + 9 data rows

    def test_filtered_gff_has_nine_sequences(self):
        self._run()
        with open(self.out_gff) as f:
            seq_lines = [l for l in f if '\tviral_sequence\t' in l or '\tplasmid\t' in l]
        self.assertEqual(len(seq_lines), 9)

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

    def test_excluded_tsv_has_correct_records(self):
        self._run()
        excluded_tsv = self.prefix + "_excluded.tsv"
        with open(excluded_tsv) as f:
            rows = [l.rstrip("\n").split("\t") for l in f]
        header = rows[0]
        # only barley1 (rrna) is excluded; barley9 passes the not_determined filter
        self.assertEqual(len(rows), 2)  # 1 header + 1 excluded
        self.assertEqual(header[0], "filter_reason")
        self.assertEqual(rows[1][0], "rrna")
        self.assertEqual(rows[1][1], "barley1")

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
