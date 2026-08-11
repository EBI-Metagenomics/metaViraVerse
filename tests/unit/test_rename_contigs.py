#!/usr/bin/env python3
"""
Unit tests for rename_contigs.py
"""

from __future__ import annotations

import bz2
import contextlib
import csv
import gzip
import io
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


class TestOpenFile(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.content = ">seq1\nACGT\n"

    def tearDown(self):
        self.tmp.cleanup()

    def test_plain_file(self):
        path = Path(self.tmp.name) / "plain.fna"
        path.write_text(self.content)
        with rc.open_file(str(path)) as f:
            self.assertEqual(f.read(), self.content)

    def test_gzip_file(self):
        path = Path(self.tmp.name) / "gzipped.fna.gz"
        with gzip.open(path, "wt") as f:
            f.write(self.content)
        with rc.open_file(str(path)) as f:
            self.assertEqual(f.read(), self.content)

    def test_bzip2_file(self):
        path = Path(self.tmp.name) / "bzipped.fna.bz2"
        with bz2.open(path, "wt") as f:
            f.write(self.content)
        with rc.open_file(str(path)) as f:
            self.assertEqual(f.read(), self.content)


class TestDefineMapKey(unittest.TestCase):
    def test_replaces_spaces(self):
        self.assertEqual(rc.define_map_key("MGYG1_1 viral_sequence"), "MGYG1_1_viral_sequence")

    def test_replaces_pipe_dash_colon(self):
        self.assertEqual(rc.define_map_key("MGYG1_1|viral_sequence-1:100"), "MGYG1_1_viral_sequence_1_100")

    def test_full_description_and_gff_attribute_collapse_to_the_same_key(self):
        # This is the trick that lets MGnify's two different representations of the
        # same record (FASTA description vs. GFF attribute) match after define_catalogue_id
        # has normalised the GFF side to the "id region|start:end" form.
        fasta_side = rc.define_map_key("MGYG000535629_9 viral_sequence|1:3862")
        gff_side = rc.define_map_key(rc.define_catalogue_id("MGYG000535629_9|viral_sequence-1:3862", "MGYG000535629_9"))
        self.assertEqual(fasta_side, gff_side)

    def test_alnum_and_underscore_untouched(self):
        self.assertEqual(rc.define_map_key("abc_123"), "abc_123")


class TestClassifyDefinition(unittest.TestCase):
    def test_matches_viral_sequence(self):
        self.assertEqual(
            rc.classify_definition("MGYG1_1 viral_sequence|1:100", "viral_sequence", "prophage", "plasmid"),
            "virus",
        )

    def test_matches_prophage(self):
        self.assertEqual(
            rc.classify_definition("MGYG1_1 prophage|1:100", "viral_sequence", "prophage", "plasmid"),
            "prophage",
        )

    def test_matches_plasmid(self):
        self.assertEqual(
            rc.classify_definition("MGYG1_1 plasmid|1:100", "viral_sequence", "prophage", "plasmid"),
            "plasmid",
        )

    def test_no_match_returns_na(self):
        self.assertEqual(
            rc.classify_definition("MGYG1_1 something_else|1:100", "viral_sequence", "prophage", "plasmid"),
            "NA",
        )

    def test_no_identifiers_given_returns_na(self):
        self.assertEqual(
            rc.classify_definition("MGYG1_1 viral_sequence|1:100", None, None, None),
            "NA",
        )


class TestDefineCatalogueId(unittest.TestCase):
    def test_mgyg_pipe_dash_pattern_transformed(self):
        self.assertEqual(
            rc.define_catalogue_id("MGYG000535629_9|viral_sequence-1:3862", "MGYG000535629_9"),
            "MGYG000535629_9 viral_sequence|1:3862",
        )

    def test_old_map_viral_pattern_transformed(self):
        self.assertEqual(
            rc.define_catalogue_id("MGYG000535629_9|viral_sequence", "MGYG000535629_9"),
            "MGYG000535629_9",
        )

    def test_old_map_plasmid_pattern_uses_contig_id(self):
        self.assertEqual(
            rc.define_catalogue_id("plasmid_1", "MGYG000535629_9"),
            "MGYG000535629_9",
        )

    def test_no_pattern_match_returns_feat_id_unchanged(self):
        # Regression test: previously id_fna was only assigned inside the pattern
        # branches, so a feat_id matching none of them raised UnboundLocalError
        # instead of falling back to the original ID.
        self.assertEqual(
            rc.define_catalogue_id("some_plain_id", "some_plain_id"),
            "some_plain_id",
        )


class TestRenameFasta(unittest.TestCase):
    def setUp(self):
        self.fna = str(FIXTURES / "file1.fna")
        self.records, self.map_dir, self.map_records, self.count = rc.rename_fasta(
            self.fna, "seq", 1, biome="rhizosphere", source_type="genome"
        )

    def test_count_two_sequences(self):
        self.assertEqual(self.count, 2)

    def test_record_ids_assigned(self):
        self.assertEqual(self.records[0].id, "seq1")
        self.assertEqual(self.records[1].id, "seq2")

    def test_map_dir_keys(self):
        self.assertIn(rc.define_map_key("MGYG000535629_9 viral_sequence|1:3862"), self.map_dir)
        self.assertIn(rc.define_map_key("MGYG000535629_15 viral_sequence|1:5155"), self.map_dir)

    def test_map_dir_values(self):
        self.assertEqual(self.map_dir[rc.define_map_key("MGYG000535629_9 viral_sequence|1:3862")], "seq1")
        self.assertEqual(self.map_dir[rc.define_map_key("MGYG000535629_15 viral_sequence|1:5155")], "seq2")

    def test_map_records_columns(self):
        row = self.map_records[0]
        # [original, temporary, short, biome, type, definition]
        self.assertEqual(len(row), 6)
        self.assertEqual(row[0], "MGYG000535629_9 viral_sequence|1:3862")
        self.assertEqual(row[1], "seq1")
        self.assertEqual(row[2], "MGYG000535629_9")
        self.assertEqual(row[3], "rhizosphere")
        self.assertEqual(row[4], "genome")

    def test_definition_defaults_to_na_without_identifiers(self):
        # no viral_sequence_identifier/prophage_identifier/plasmid_identifier given
        self.assertEqual(self.map_records[0][5], "NA")

    def test_definition_classified_from_identifier(self):
        _, _, map_records, _ = rc.rename_fasta(
            self.fna, "seq", 1,
            viral_sequence_identifier="viral_sequence",
            prophage_identifier="prophage",
            plasmid_identifier="plasmid",
        )
        self.assertEqual(map_records[0][5], "virus")
        self.assertEqual(map_records[1][5], "virus")

    def test_definition_fixed_value_overrides_identifiers(self):
        _, _, map_records, _ = rc.rename_fasta(
            self.fna, "seq", 1,
            definition="plasmid",
            viral_sequence_identifier="viral_sequence",
        )
        self.assertEqual(map_records[0][5], "plasmid")

    def test_third_party_prefixes_definition(self):
        _, _, map_records, _ = rc.rename_fasta(
            self.fna, "seq", 1, definition="plasmid", third_party=True,
        )
        self.assertEqual(map_records[0][5], "third_party_plasmid")

    def test_start_accession_offset(self):
        _, _, records, _ = rc.rename_fasta(
            self.fna, "seq", 10
        )
        self.assertEqual(records[0][1], "seq10")
        self.assertEqual(records[1][1], "seq11")

    def test_third_party_map_key_uses_short_id_not_full_description(self):
        # Pyrodigal (used for third-party gene calling) only ever writes the short
        # FASTA id (the token before the first whitespace) as the GFF seqid, never
        # the full multi-word description -- so the third-party map_dir key must be
        # built from record.id, not record.description, or GFF lookups will miss.
        _, map_dir, _, _ = rc.rename_fasta(
            str(FIXTURES / "file4.fna"), "seq", 1, third_party=True,
        )
        self.assertIn(rc.define_map_key("TP4"), map_dir)
        self.assertNotIn(rc.define_map_key("TP4 marine sample virus description text"), map_dir)


class TestDefineGffMapKey(unittest.TestCase):
    def test_third_party_uses_column_one_regardless_of_feature_type(self):
        cds_parts = ["seq1", "pyrodigal", "CDS", "1", "10", ".", "+", "0", "ID=seq1_1"]
        proteins, key = rc.define_gff_map_key(cds_parts, third_party=True, current_id="unused")
        self.assertEqual(proteins, 1)
        self.assertEqual(key, "seq1")

    def test_mgnify_cds_line_inherits_current_id(self):
        cds_parts = ["MGYG1_1", "Prodigal", "CDS", "1", "10", ".", "+", "0", "ID=MGYG1_1_00001"]
        proteins, key = rc.define_gff_map_key(cds_parts, third_party=False, current_id="carried_forward_id")
        self.assertEqual(proteins, 1)
        self.assertEqual(key, "carried_forward_id")

    def test_mgnify_non_cds_line_derives_key_from_id_attribute(self):
        parts = ["MGYG000535629_9", "VIRify", "viral_sequence", "1", "3862", ".", ".", ".",
                 "ID=MGYG000535629_9|viral_sequence-1:3862"]
        proteins, key = rc.define_gff_map_key(parts, third_party=False, current_id=None)
        self.assertEqual(proteins, 0)
        self.assertEqual(key, rc.define_map_key("MGYG000535629_9 viral_sequence|1:3862"))


class TestRenameGff(unittest.TestCase):
    def setUp(self):
        _, self.map_dir, _, _ = rc.rename_fasta(
            str(FIXTURES / "file1.fna"), "seq", 1
        )
        self.seqs_count, self.count_proteins, self.gff_records = rc.rename_gff(
            str(FIXTURES / "file1.gff"), self.map_dir
        )

    def test_returns_sequence_count(self):
        self.assertEqual(self.seqs_count, 2)

    def test_returns_protein_count(self):
        # file1.gff has 2 CDS lines per sequence, 2 sequences = 4 CDS features total
        self.assertEqual(self.count_proteins, 4)

    def test_seqid_column_renamed(self):
        seq_lines = [l for l in self.gff_records if l.split('\t')[2] == 'viral_sequence']
        seqids = {l.split('\t')[0] for l in seq_lines}
        self.assertEqual(seqids, {"seq1", "seq2"})

    def test_cds_lines_also_get_seqid_renamed(self):
        cds_lines = [l for l in self.gff_records if l.split('\t')[2] == 'CDS']
        seqids = {l.split('\t')[0] for l in cds_lines}
        self.assertEqual(seqids, {"seq1", "seq2"})

    def test_cds_protein_id_attribute_untouched(self):
        cds_lines = [l for l in self.gff_records if l.split('\t')[2] == 'CDS']
        content = ''.join(cds_lines)
        self.assertIn("ID=MGYG000535629_00028", content)
        self.assertIn("ID=MGYG000535629_00070", content)

    def test_unmatched_id_falls_back_to_original_line_and_warns(self):
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stderr):
            seqs_count, _, gff_records = rc.rename_gff(str(FIXTURES / "file1.gff"), {})
        self.assertIn("Warning: no", stderr.getvalue())
        # original (unrenamed) seqid is kept when nothing matches
        self.assertTrue(any(l.startswith("MGYG000535629_9\t") for l in gff_records))

    def test_third_party_true_does_not_use_mgyg_pattern(self):
        # With third_party=True, file1's MGYG-encoded feat_id is never passed through
        # define_catalogue_id -- only column 1 (the bare contig id) is used as the key,
        # so it won't match map_dir (which is keyed by the full description) and the
        # line is left with its original (unrenamed) seqid.
        _, _, gff_records = rc.rename_gff(str(FIXTURES / "file1.gff"), self.map_dir, third_party=True)
        seq_lines = [l for l in gff_records if l.split('\t')[2] == 'viral_sequence']
        seqids = {l.split('\t')[0] for l in seq_lines}
        self.assertEqual(seqids, {"MGYG000535629_9", "MGYG000535629_15"})


class TestRenameGffThirdParty(unittest.TestCase):
    """third_party=True keys strictly off column 1, matching what Pyrodigal writes."""

    def test_plain_single_token_id_renamed(self):
        _, map_dir, _, _ = rc.rename_fasta(str(FIXTURES / "file3.fna"), "seq", 1, third_party=True)
        seqs_count, _, gff_records = rc.rename_gff(str(FIXTURES / "file3.gff"), map_dir, third_party=True)
        self.assertEqual(seqs_count, 1)
        plasmid_lines = [l for l in gff_records if l.split('\t')[2] == 'plasmid']
        self.assertEqual(plasmid_lines[0].split('\t')[0], "seq1")

    def test_multi_word_description_still_renamed_via_short_id_key(self):
        # Regression test for the record.id-vs-record.description bug: file4's FASTA
        # header has descriptive text beyond the ID, but Pyrodigal's GFF (file4.gff)
        # only ever wrote the short id "TP4" in column 1.
        _, map_dir, _, _ = rc.rename_fasta(str(FIXTURES / "file4.fna"), "seq", 1, third_party=True)
        seqs_count, _, gff_records = rc.rename_gff(str(FIXTURES / "file4.gff"), map_dir, third_party=True)
        self.assertEqual(seqs_count, 1)
        cds_lines = [l for l in gff_records if l.split('\t')[2] == 'CDS']
        self.assertEqual(cds_lines[0].split('\t')[0], "seq1")


class TestWriteGff(unittest.TestCase):
    def test_writes_header_and_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = str(Path(tmp) / "out.gff")
            rc.write_gff(out, ["seq1\tsrc\tCDS\t1\t10\t.\t+\t0\tID=x\n"])
            with open(out) as f:
                content = f.read()
        self.assertTrue(content.startswith("##gff-version 3\n"))
        self.assertIn("seq1\tsrc\tCDS", content)


class TestWriteMapfile(unittest.TestCase):
    def test_header_and_rows(self):
        data = [
            ["MGYG000535629_9 viral_sequence|1:3862", "seq1", "MGYG000535629_9", "rhizosphere", "genome", "virus"],
            ["MGYG000535629_15 viral_sequence|1:5155", "seq2", "MGYG000535629_15", "rhizosphere", "genome", "virus"],
        ]
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            path = f.name
        try:
            rc.write_mapfile(path, data)
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                rows = list(reader)
            self.assertEqual(rows[0], ["original", "temporary", "short", "biome", "type", "definition"])
            self.assertEqual(rows[1][1], "seq1")
            self.assertEqual(rows[2][3], "rhizosphere")
            self.assertEqual(rows[1][5], "virus")
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
            lines = [l for l in f if not l.startswith("#")]
        seq_lines = [l for l in lines if l.split('\t')[2] == 'viral_sequence']
        seqids = {l.split('\t')[0] for l in seq_lines}
        self.assertEqual(seqids, {"seq1", "seq2"})


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
            lines = [l for l in f if not l.startswith("#")]
        seqids = {l.split('\t')[0] for l in lines}
        self.assertEqual(seqids, {"seq1", "seq2", "seq3", "seq4"})

    def test_no_per_file_outputs_in_combine_mode(self):
        self._run()
        self.assertFalse(os.path.exists(str(Path(self.tmp.name) / "renamed_file1.fna")))
        self.assertFalse(os.path.exists(str(Path(self.tmp.name) / "renamed_file2.fna")))


class TestInputArgsValidation(unittest.TestCase):
    """Argument-parsing validation for the third-party options."""

    def _parse(self, argv):
        sys.argv = ["rename_contigs.py"] + argv
        with contextlib.redirect_stderr(io.StringIO()):
            return rc.input_args()

    def test_no_fasta_and_no_fasta_tp_errors(self):
        with self.assertRaises(SystemExit):
            self._parse(["--type", "genome"])

    def test_fasta_tp_without_type_tp_errors(self):
        with self.assertRaises(SystemExit):
            self._parse([
                "--fasta", str(FIXTURES / "file1.fna"),
                "--fasta-tp", str(FIXTURES / "file2.fna"),
                "--combine",
            ])

    def test_fasta_tp_type_tp_length_mismatch_errors(self):
        with self.assertRaises(SystemExit):
            self._parse([
                "--fasta-tp", str(FIXTURES / "file2.fna"), str(FIXTURES / "file3.fna"),
                "--type-tp", "virus",
                "--combine",
            ])

    def test_fasta_tp_invalid_type_choice_errors(self):
        with self.assertRaises(SystemExit):
            self._parse([
                "--fasta-tp", str(FIXTURES / "file2.fna"),
                "--type-tp", "chromosome",
                "--combine",
            ])

    def test_fasta_tp_without_combine_errors(self):
        with self.assertRaises(SystemExit):
            self._parse([
                "--fasta-tp", str(FIXTURES / "file2.fna"),
                "--type-tp", "virus",
            ])

    def test_fasta_tp_alone_is_valid(self):
        args = self._parse([
            "--fasta-tp", str(FIXTURES / "file2.fna"),
            "--type-tp", "virus",
            "--combine",
        ])
        self.assertEqual(args.fasta, [])
        self.assertEqual(args.fasta_tp, [str(FIXTURES / "file2.fna")])


class TestMainThirdPartyOnly(unittest.TestCase):
    """Integration test: only third-party input, no MGnify samplesheet at all."""

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
                "--fasta-tp",  str(FIXTURES / "file3.fna"),
                "--gff-tp",    str(FIXTURES / "file3.gff"),
                "--biome-tp",  "marine",
                "--source-tp", "metagenome",
                "--type-tp",   "plasmid",
                "--prefix",    "seq",
                "--map",       "combined.tsv",
                "--combine",
            ]
            rc.main()
        finally:
            os.chdir(orig_dir)

    def test_combined_output_exists(self):
        self._run()
        self.assertTrue(os.path.exists(str(Path(self.tmp.name) / "combined.fna")))
        self.assertTrue(os.path.exists(str(Path(self.tmp.name) / "combined.tsv")))
        self.assertTrue(os.path.exists(str(Path(self.tmp.name) / "combined.gff")))

    def test_accession_starts_at_one(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.fna")) as f:
            headers = [l.strip().lstrip(">") for l in f if l.startswith(">")]
        self.assertEqual(headers, ["seq1"])

    def test_map_record_type_is_source_and_definition_is_prefixed(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.tsv")) as f:
            reader = csv.reader(f, delimiter="\t")
            rows = list(reader)
        self.assertEqual(rows[1][1], "seq1")
        self.assertEqual(rows[1][3], "marine")
        self.assertEqual(rows[1][4], "metagenome")           # from --source-tp, not prefixed
        self.assertEqual(rows[1][5], "third_party_plasmid")  # from --type-tp, prefixed

    def test_gff_ids_renamed(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.gff")) as f:
            lines = [l for l in f if not l.startswith("#")]
        seqids = {l.split('\t')[0] for l in lines}
        self.assertEqual(seqids, {"seq1"})


class TestMainMgnifyPlusThirdParty(unittest.TestCase):
    """Integration test: MGnify input combined into one output with third-party viruses and plasmids."""

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
                "--fasta",     str(FIXTURES / "file1.fna"),
                "--gff",       str(FIXTURES / "file1.gff"),
                "--biome",     "rhizosphere",
                "--type",      "genome",
                "--viral-sequence-identifier", "viral_sequence",
                "--prophage-identifier",       "prophage",
                "--plasmid-identifier",        "plasmid",
                "--fasta-tp",  str(FIXTURES / "file3.fna"), str(FIXTURES / "file4.fna"),
                "--gff-tp",    str(FIXTURES / "file3.gff"), str(FIXTURES / "file4.gff"),
                "--biome-tp",  "human gut", "marine",
                "--source-tp", "metagenome", "genome",
                "--type-tp",   "virus", "plasmid",
                "--prefix",    "seq",
                "--map",       "combined.tsv",
                "--combine",
            ]
            rc.main()
        finally:
            os.chdir(orig_dir)

    def test_single_combined_fasta_has_all_records(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.fna")) as f:
            headers = [l.strip().lstrip(">") for l in f if l.startswith(">")]
        self.assertEqual(headers, ["seq1", "seq2", "seq3", "seq4"])

    def test_single_combined_gff_has_one_header(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.gff")) as f:
            headers = [l for l in f if l.startswith("##gff-version")]
        self.assertEqual(len(headers), 1)

    def test_combined_gff_ids_renamed_for_both_mgnify_and_third_party(self):
        # Regression test: MGnify (file1, MGYG-encoded attribute) and third-party
        # (file3: plain single-token id; file4: multi-word description) all need to
        # end up correctly renamed in the shared combined GFF, each via its own
        # matching strategy.
        self._run()
        with open(str(Path(self.tmp.name) / "combined.gff")) as f:
            lines = [l for l in f if not l.startswith("#")]
        seqids = {l.split('\t')[0] for l in lines}
        self.assertEqual(seqids, {"seq1", "seq2", "seq3", "seq4"})

    def test_mgnify_rows_have_type_and_derived_definition(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.tsv")) as f:
            rows = list(csv.reader(f, delimiter="\t"))
        self.assertEqual(rows[1][4], "genome")  # type: from --type (MGnify), unaffected by this change
        self.assertEqual(rows[1][5], "virus")   # definition: derived from --viral-sequence-identifier

    def test_third_party_rows_have_source_as_type_and_prefixed_definition(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.tsv")) as f:
            rows = list(csv.reader(f, delimiter="\t"))
        # rows[3] = file3.fna (1 record, type-tp=virus, source-tp=metagenome)
        self.assertEqual(rows[3][4], "metagenome")
        self.assertEqual(rows[3][5], "third_party_virus")
        # rows[4] = file4.fna (1 record, type-tp=plasmid, source-tp=genome)
        self.assertEqual(rows[4][4], "genome")
        self.assertEqual(rows[4][5], "third_party_plasmid")

    def test_accession_continuity_across_sources(self):
        self._run()
        with open(str(Path(self.tmp.name) / "combined.tsv")) as f:
            rows = list(csv.reader(f, delimiter="\t"))
        temporary_names = [r[1] for r in rows[1:]]
        self.assertEqual(temporary_names, ["seq1", "seq2", "seq3", "seq4"])


if __name__ == "__main__":
    unittest.main()
