#!/usr/bin/env python3
"""
Unit tests for choose_sequences.py
"""

import csv
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


class TestTypePriority(unittest.TestCase):
    def test_mgnify_metagenome_is_highest(self):
        self.assertEqual(cs.type_priority("metagenome", "virus"), 0)

    def test_mgnify_genome_beats_third_party(self):
        mgnify_genome = cs.type_priority("genome", "virus")
        third_party_metagenome = cs.type_priority("metagenome", "third_party_virus")
        self.assertLess(mgnify_genome, third_party_metagenome)

    def test_metagenome_beats_genome_within_same_origin(self):
        self.assertLess(
            cs.type_priority("metagenome", "virus"),
            cs.type_priority("genome", "virus"),
        )
        self.assertLess(
            cs.type_priority("metagenome", "third_party_plasmid"),
            cs.type_priority("genome", "third_party_plasmid"),
        )

    def test_third_party_genome_is_lowest(self):
        ranks = [
            cs.type_priority("metagenome", "virus"),
            cs.type_priority("genome", "virus"),
            cs.type_priority("metagenome", "third_party_virus"),
            cs.type_priority("genome", "third_party_virus"),
        ]
        self.assertEqual(ranks, sorted(ranks))
        self.assertEqual(ranks[-1], cs.type_priority("genome", "third_party_virus"))

    def test_unknown_source_type_is_lowest(self):
        self.assertEqual(cs.type_priority("NA", "virus"), 99)

    def test_missing_definition_treated_as_mgnify(self):
        self.assertEqual(cs.type_priority("metagenome", "NA"), 0)
        self.assertEqual(cs.type_priority("metagenome", ""), 0)


class TestReadMap(unittest.TestCase):
    def setUp(self):
        self.mapping = cs.read_map([str(FIXTURES / "barley10.map.tsv")])

    def test_loads_all_entries(self):
        self.assertEqual(len(self.mapping), 13)

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

    def test_definition_column_read(self):
        self.assertEqual(self.mapping["barley1"]["definition"], "virus")
        self.assertEqual(self.mapping["barley3"]["definition"], "plasmid")
        self.assertEqual(self.mapping["barley11"]["definition"], "prophage")

    def test_third_party_definition_read(self):
        self.assertEqual(self.mapping["barley12"]["definition"], "third_party_prophage")
        self.assertEqual(self.mapping["barley12"]["type"], "metagenome")

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


class TestReadFaa(unittest.TestCase):
    def test_reads_across_multiple_files(self):
        proteins = cs.read_faa([str(FIXTURES / "viruses.faa"), str(FIXTURES / "plasmids.faa")])
        self.assertIn("barley1_00001", proteins)
        self.assertIn("barley3_00001", proteins)
        self.assertEqual(len(proteins), 10)

    def test_missing_file_list_is_empty(self):
        self.assertEqual(cs.read_faa([]), {})


class TestFilterReason(unittest.TestCase):
    def _entry(self, rrna, category, checkv_quality=None, viral_genes="0", kmer_freq="1.0"):
        q = {"checkv_quality": checkv_quality, "viral_genes": viral_genes, "kmer_freq": kmer_freq} \
            if checkv_quality is not None else None
        return {"rrna": rrna, "category": category, "quality": q}

    def test_clean_virus_passes(self):
        self.assertIsNone(cs.filter_reason(self._entry("No", "virus", "Low-quality")))

    def test_rrna_virus_reason(self):
        self.assertEqual(cs.filter_reason(self._entry("Yes", "virus", "Low-quality")), "rrna")

    def test_rrna_prophage_reason(self):
        # prophages get the same rRNA/quality filtering as free viruses
        self.assertEqual(cs.filter_reason(self._entry("Yes", "prophage", "Low-quality")), "rrna")

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
    def _entry(self, rrna, category, checkv_quality=None, viral_genes="0", kmer_freq="1.0"):
        q = {"checkv_quality": checkv_quality, "viral_genes": viral_genes, "kmer_freq": kmer_freq} \
            if checkv_quality is not None else None
        return {"rrna": rrna, "category": category, "quality": q}

    def test_clean_virus_passes(self):
        self.assertTrue(cs.passes_filter(self._entry("No", "virus", "Low-quality")))

    def test_rrna_virus_filtered(self):
        self.assertFalse(cs.passes_filter(self._entry("Yes", "virus", "Low-quality")))

    def test_not_determined_virus_filtered(self):
        self.assertFalse(cs.passes_filter(self._entry("No", "virus", "Not-determined", viral_genes="0", kmer_freq="3.0")))

    def test_no_quality_data_passes(self):
        self.assertTrue(cs.passes_filter(self._entry("No", "virus", None)))


class TestChooseSeqsCrossCategory(unittest.TestCase):
    """Direct (non-main()) tests of the cross-category dedup/priority/biome-merge logic.

    Fixture design (see fixtures/choose_sequences/):
    - barley5  (viruses.fasta):   MGnify, genome,     definition=virus
    - barley12 (prophages.fasta): duplicate of barley5's sequence, third-party,
                                   metagenome, definition=third_party_prophage
      -> barley5 must win (MGnify beats third-party, even genome vs metagenome).
    - barley2  (viruses.fasta):   MGnify, genome,     definition=virus
    - barley17 (prophages.fasta): duplicate of barley2's sequence, MGnify,
                                   metagenome, definition=prophage
      -> barley17 must win (metagenome beats genome within the same MGnify origin),
         and the deduplicated sequence's category flips from 'virus' to 'prophage'.
    """

    def setUp(self):
        self.mapping = cs.read_map([str(FIXTURES / "barley10.map.tsv")])
        self.seen, self.duplicates = cs.choose_seqs(
            [
                ("virus", str(FIXTURES / "viruses.fasta")),
                ("prophage", str(FIXTURES / "prophages.fasta")),
                ("plasmid", str(FIXTURES / "plasmids.fasta")),
            ],
            self.mapping,
            rna_sequences=set(),
            quality_data={},
        )

    def _entry_by_seq_id(self, seq_id):
        return next(e for e in self.seen.values() if e["seq_id"] == seq_id)

    def _duplicate_by_seq_id(self, seq_id):
        return next(e for e in self.duplicates if e["seq_id"] == seq_id)

    def test_total_unique_sequences(self):
        # 10 barley1-10 + barley11 unique + (barley12 dup of barley5) + (barley17 dup of barley2)
        # = 13 records in, 11 unique hashes
        self.assertEqual(len(self.seen), 11)

    def test_mgnify_beats_third_party_even_with_lower_genome_rank(self):
        entry = self._entry_by_seq_id("barley5")
        self.assertEqual(entry["category"], "virus")
        self.assertEqual(entry["definition"], "virus")

    def test_third_party_duplicate_is_not_the_winner(self):
        seq_ids = {e["seq_id"] for e in self.seen.values()}
        self.assertNotIn("barley12", seq_ids)

    def test_biomes_merged_across_winner_and_loser(self):
        entry = self._entry_by_seq_id("barley5")
        self.assertEqual(entry["biomes"], {"rhizosphere", "desert"})

    def test_metagenome_beats_genome_within_mgnify(self):
        entry = self._entry_by_seq_id("barley17")
        self.assertEqual(entry["category"], "prophage")
        self.assertEqual(entry["definition"], "prophage")

    def test_category_can_flip_from_original_bucket(self):
        # barley2 (from viruses.fasta, category='virus') loses to barley17
        # (from prophages.fasta, category='prophage'): barley2 must not survive
        seq_ids = {e["seq_id"] for e in self.seen.values()}
        self.assertNotIn("barley2", seq_ids)
        self.assertIn("barley17", seq_ids)

    def test_biomes_merged_for_metagenome_winner(self):
        entry = self._entry_by_seq_id("barley17")
        self.assertEqual(entry["biomes"], {"rhizosphere", "peat"})

    def test_unrelated_sequence_keeps_its_own_category(self):
        entry = self._entry_by_seq_id("barley11")
        self.assertEqual(entry["category"], "prophage")

    def test_plasmid_category_unaffected(self):
        entry = self._entry_by_seq_id("barley3")
        self.assertEqual(entry["category"], "plasmid")

    def test_losing_records_are_returned_as_duplicates(self):
        seq_ids = {e["seq_id"] for e in self.duplicates}
        self.assertEqual(seq_ids, {"barley12", "barley2"})

    def test_duplicate_count_matches_lost_records(self):
        # 13 records in, 11 unique hashes -> 2 records lost a collision
        self.assertEqual(len(self.duplicates), 2)

    def test_duplicate_carries_its_own_single_biome_not_the_merged_set(self):
        # barley12 lost to barley5; the merged set on the winner is
        # {"rhizosphere", "desert"}, but the duplicate itself only carries
        # its own biome ("desert").
        entry = self._duplicate_by_seq_id("barley12")
        self.assertEqual(entry["biomes"], {"desert"})

    def test_duplicate_keeps_its_own_definition_and_category(self):
        entry = self._duplicate_by_seq_id("barley12")
        self.assertEqual(entry["category"], "prophage")
        self.assertEqual(entry["definition"], "third_party_prophage")

    def test_demoted_former_winner_is_also_a_duplicate(self):
        # barley2 was the first record seen for its hash (so it started as the
        # 'seen' winner) but was later demoted when higher-priority barley17
        # showed up -- it must still end up in duplicates, not just vanish.
        entry = self._duplicate_by_seq_id("barley2")
        self.assertEqual(entry["category"], "virus")


class TestReadInputGff(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_indexes_by_column_one_even_when_id_attribute_is_stale(self):
        # Regression test: rename_contigs.py rewrites column 1 (the seqid) to the new
        # temporary name but leaves the sequence-level ID= attribute as the original,
        # un-renamed value. read_input_gff must key gff_data off column 1 -- entry
        # lookups in write_final_files use entry['seq_id'], which is the FASTA
        # record.id (also the temporary name) -- or a genuinely-renamed GFF would
        # never be found.
        gff_with_stale_attribute = (
            "##gff-version 3\n"
            "seq1\tVIRify\tviral_sequence\t1\t61\t.\t.\t.\t"
            "ID=MGYG000535629_9|viral_sequence-1:3862;checkv_quality=Low-quality\n"
            "seq1\tProdigal:002006\tCDS\t1\t30\t.\t+\t0\tID=MGYG000535629_00028;product=hypothetical protein\n"
        )
        gff_path = str(Path(self.tmp.name) / "stale_attr.gff")
        with open(gff_path, "w") as f:
            f.write(gff_with_stale_attribute)

        gff_data, source_map = cs.read_input_gff([gff_path])

        self.assertIn("seq1", gff_data)
        self.assertEqual(len(gff_data["seq1"]), 2)
        self.assertEqual(source_map["seq1"], "VIRify")

    def test_source_map_ignores_cds_lines(self):
        gff_path = str(Path(self.tmp.name) / "test.gff")
        with open(gff_path, "w") as f:
            f.write(
                "##gff-version 3\n"
                "seq1\tgeNomad\tplasmid\t1\t61\t.\t.\t.\tID=seq1;mobile_element_type=plasmid\n"
                "seq1\tProdigal:002006\tCDS\t1\t30\t.\t+\t0\tID=seq1_00001\n"
            )
        _, source_map = cs.read_input_gff([gff_path])
        self.assertEqual(source_map["seq1"], "geNomad")


class TestMainIntegration(unittest.TestCase):
    """Run main() end-to-end with fixture files and verify output counts."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.prefix = str(Path(self.tmp.name) / "test")

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, extra_args=None):
        argv = [
            "choose_sequences.py",
            "--viruses",       str(FIXTURES / "viruses.fasta"),
            "--viruses-gff",   str(FIXTURES / "viruses.gff"),
            "--viruses-faa",   str(FIXTURES / "viruses.faa"),
            "--prophages",     str(FIXTURES / "prophages.fasta"),
            "--prophages-gff", str(FIXTURES / "prophages.gff"),
            "--prophages-faa", str(FIXTURES / "prophages.faa"),
            "--plasmids",      str(FIXTURES / "plasmids.fasta"),
            "--plasmids-gff",  str(FIXTURES / "plasmids.gff"),
            "--plasmids-faa",  str(FIXTURES / "plasmids.faa"),
            "--map",           str(FIXTURES / "barley10.map.tsv"),
            "--rrna",          str(FIXTURES / "barley10.gff"),
            "--quality",       str(FIXTURES / "barley10_quality.tsv"),
            "--output-prefix", self.prefix,
        ]
        if extra_args:
            argv += extra_args
        sys.argv = argv
        cs.main()

    def _read_tsv(self, path):
        with open(path) as f:
            return list(csv.reader(f, delimiter="\t"))

    def test_metadata_has_all_unique_sequences(self):
        self._run()
        rows = self._read_tsv(self.prefix + "_metadata.tsv")
        # 13 input records, 2 are exact duplicates (barley12 of barley5, barley17 of barley2)
        # -> 11 unique sequences + 1 header
        self.assertEqual(len(rows), 12)

    def test_combined_filtered_fna_excludes_rrna_and_duplicates(self):
        self._run()
        with open(self.prefix + "_filtered.fna") as f:
            headers = [l for l in f if l.startswith(">")]
        # 13 input records -> 11 unique hashes (barley12 dupes barley5, barley17 dupes
        # barley2) -> barley1 additionally excluded by the rRNA filter -> 10 remain
        self.assertEqual(len(headers), 10)

    def test_virus_specific_fna_has_correct_members(self):
        self._run()
        with open(self.prefix + "_virus_filtered.fna") as f:
            ids = {l.strip().lstrip(">") for l in f if l.startswith(">")}
        self.assertIn("barley5", ids)
        self.assertNotIn("barley1", ids)   # excluded by rRNA
        self.assertNotIn("barley2", ids)   # lost the dedup to barley17
        self.assertNotIn("barley12", ids)  # never wins (third-party)
        self.assertNotIn("barley17", ids)  # ends up in the prophage bucket, not virus

    def test_prophage_specific_fna_has_correct_members(self):
        self._run()
        with open(self.prefix + "_prophage_filtered.fna") as f:
            ids = {l.strip().lstrip(">") for l in f if l.startswith(">")}
        self.assertEqual(ids, {"barley11", "barley17"})

    def test_plasmid_specific_fna_unaffected(self):
        self._run()
        with open(self.prefix + "_plasmid_filtered.fna") as f:
            ids = {l.strip().lstrip(">") for l in f if l.startswith(">")}
        self.assertEqual(ids, {"barley3", "barley4"})

    def test_faa_drops_proteins_for_deduplicated_losers(self):
        self._run()
        with open(self.prefix + "_filtered.faa") as f:
            content = f.read()
        # barley5 won over barley12 -> barley5's protein kept, barley12's dropped
        self.assertIn(">barley5_00001", content)
        self.assertNotIn(">barley12_00001", content)
        # barley17 won over barley2 -> barley17's protein kept, barley2's dropped
        self.assertIn(">barley17_00001", content)
        self.assertNotIn(">barley2_00001", content)
        # barley1 is excluded by the rRNA filter entirely -> its protein must not appear
        self.assertNotIn(">barley1_00001", content)

    def test_category_specific_faa_matches_category_fna(self):
        self._run()
        with open(self.prefix + "_prophage_filtered.faa") as f:
            content = f.read()
        self.assertIn(">barley11_00001", content)
        self.assertIn(">barley17_00001", content)
        self.assertNotIn(">barley12_00001", content)

    def test_excluded_tsv_has_rrna_record(self):
        self._run()
        rows = self._read_tsv(self.prefix + "_excluded.tsv")
        header = rows[0]
        self.assertEqual(header[0], "filter_reason")
        by_id = {r[1]: r[0] for r in rows[1:]}
        self.assertEqual(by_id["barley1"], "rrna")

    def test_excluded_tsv_has_duplicate_records(self):
        # barley12 lost to barley5, and barley2 (the original winner for its hash)
        # was later demoted by barley17 -- both must show up in excluded.tsv rather
        # than silently vanishing.
        self._run()
        rows = self._read_tsv(self.prefix + "_excluded.tsv")
        by_id = {r[1]: r[0] for r in rows[1:]}
        self.assertEqual(by_id["barley12"], "duplicate")
        self.assertEqual(by_id["barley2"], "duplicate")

    def test_excluded_tsv_covers_every_filtered_out_record(self):
        # Every record that never makes it into the combined filtered FASTA must
        # appear in excluded.tsv -- whether it lost a quality filter or a dedup
        # collision -- so nothing is filtered out without a trace.
        self._run()
        with open(self.prefix + "_filtered.fna") as f:
            kept_ids = {l.strip().lstrip(">") for l in f if l.startswith(">")}
        metadata_rows = self._read_tsv(self.prefix + "_metadata.tsv")[1:]
        excluded_rows = self._read_tsv(self.prefix + "_excluded.tsv")[1:]

        all_seen_ids = {r[0] for r in metadata_rows}
        excluded_ids = {r[1] for r in excluded_rows}
        duplicate_ids = {"barley12", "barley2"}

        self.assertEqual(excluded_ids, (all_seen_ids - kept_ids) | duplicate_ids)

    def test_tsv_has_correct_columns(self):
        self._run()
        with open(self.prefix + "_metadata.tsv") as f:
            header = f.readline().rstrip("\n").split("\t")
        expected_start = ["sequence_id", "original_name", "description", "type",
                          "source_of_prediction", "biomes", "sequence_length", "rrna", "sequence_sha256"]
        self.assertEqual(header[:9], expected_start)
        for col in cs.QUALITY_COLUMNS:
            self.assertIn(col, header)
        self.assertEqual(header[-1], "definition")

    def test_definition_column_reflects_winning_category(self):
        self._run()
        rows = self._read_tsv(self.prefix + "_filtered.tsv")
        header = rows[0]
        seq_id_idx = header.index("sequence_id")
        definition_idx = header.index("definition")
        by_id = {r[seq_id_idx]: r[definition_idx] for r in rows[1:]}
        self.assertEqual(by_id["barley5"], "virus")
        self.assertEqual(by_id["barley17"], "prophage")

    def test_biome_merged_in_output(self):
        self._run()
        rows = self._read_tsv(self.prefix + "_filtered.tsv")
        header = rows[0]
        seq_id_idx = header.index("sequence_id")
        biome_idx = header.index("biomes")
        by_id = {r[seq_id_idx]: r[biome_idx] for r in rows[1:]}
        self.assertEqual(by_id["barley5"], "desert,rhizosphere")
        self.assertEqual(by_id["barley17"], "peat,rhizosphere")

    def test_source_of_prediction_populated(self):
        self._run()
        rows = self._read_tsv(self.prefix + "_metadata.tsv")
        header = rows[0]
        src_idx = header.index("source_of_prediction")
        seq_id_idx = header.index("sequence_id")
        by_id = {r[seq_id_idx]: r[src_idx] for r in rows[1:]}
        self.assertEqual(by_id["barley1"], "VIRify")
        self.assertEqual(by_id["barley3"], "geNomad")

    def test_rrna_column_populated(self):
        self._run()
        rows = self._read_tsv(self.prefix + "_metadata.tsv")
        header = rows[0]
        rrna_idx = header.index("rrna")
        seq_id_idx = header.index("sequence_id")
        by_id = {r[seq_id_idx]: r[rrna_idx] for r in rows[1:]}
        self.assertEqual(by_id["barley1"], "Yes")
        # barley2 loses its dedup to barley17 (same sequence, higher priority) so it
        # never appears as its own row -- barley6 is a genuine, non-duplicated record
        self.assertEqual(by_id["barley6"], "No")


if __name__ == "__main__":
    unittest.main()
