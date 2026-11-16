#!/usr/bin/env python3
"""Unit tests for collect_metadata.py: rows are keyed by catalogue ID."""

import csv
import os
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import collect_metadata as cm

META_HEADER = "sequence_id\toriginal_name\tdescription\ttype\tsource_of_prediction\tbiomes\tsequence_length\trrna\tsequence_sha256\tcheckv_quality\tdefinition\n"

# Whole catalogue after an update: MGYV1 (existing) and MGYV5 (new) have the same
# original name but different sequences; MGYV5 joined the existing cluster of MGYV1.
METADATA = META_HEADER + (
    "MGYV1\tMGYG1_1|viral_sequence\tMGYG1_1|viral_sequence\tgenome\tMGnify\tsoil\t100\tno\th1\tHigh-quality\tvirus\n"
    "MGYV2\tMGYG2_1|prophage-1:69\tMGYG2_1|prophage-1:69\tgenome\tMGnify\tmarine\t69\tno\th2\tmissing\tprophage\n"
    "MGYV5\tMGYG1_1|viral_sequence\tMGYG1_1|viral_sequence\tgenome\tMGnify\tlake\t100\tno\th5\tLow-quality\tvirus\n"
    "MGYV6\tMGYG3_1|plasmid\tMGYG3_1|plasmid\tgenome\tMGnify\tsoil\t50\tno\th6\tmissing\tplasmid\n"
)
MAP = (
    "original\ttemporary\tshort\tbiome\ttype\tdefinition\n"
    "MGYG1_1|viral_sequence\tMGYV1\tMGYG1_1|viral_sequence\tsoil\tgenome\tvirus\n"
    "MGYG2_1|prophage-1:69\tMGYV2\tMGYG2_1|prophage-1:69\tmarine\tgenome\tprophage\n"
    "MGYG1_1|viral_sequence\tMGYV5\tMGYG1_1|viral_sequence\tlake\tgenome\tvirus\n"
    "MGYG3_1|plasmid\tMGYV6\tMGYG3_1|plasmid\tsoil\tgenome\tplasmid\n"
)
VIRUS_CLUSTERS = "object\tcluster\nMGYV1\tMGYV1\nMGYV5\tMGYV1\nMGYV2\tMGYV2\n"
PLASMID_CLUSTERS = "object\tcluster\nMGYV6\tMGYV6\n"
GENOMES = "Genome\tGenome_accession\tLineage\tSample_accession\tStudy_accession\nMGYG1\tGCA1\td__Bacteria\tS1\tP1\n"
GENOMAD = "seq_name\ttaxonomy\nMGYV1\tViruses;Caudoviricetes\nMGYV2\tViruses\n"


class TestCollectMetadata(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cwd = os.getcwd()
        os.chdir(self.tmp.name)
        for name, text in (("meta.tsv", METADATA), ("map.tsv", MAP), ("vc.tsv", VIRUS_CLUSTERS),
                           ("pc.tsv", PLASMID_CLUSTERS), ("genomes.tsv", GENOMES), ("genomad.tsv", GENOMAD)):
            Path(name).write_text(text)
        sys.argv = [
            "collect_metadata.py", "--combined_meta", "meta.tsv", "--viruses_cluster", "vc.tsv",
            "--plasmids_cluster", "pc.tsv", "--additional_metadata", "genomes.tsv", "--genomad", "genomad.tsv",
            "--map", "map.tsv", "--output_viruses", "viruses.tsv", "--output_plasmids", "plasmids.tsv",
            "--output_reps", "reps.tsv",
        ]
        cm.main()

    def tearDown(self):
        os.chdir(self.cwd)
        self.tmp.cleanup()

    def rows(self, name):
        with open(name) as fh:
            return {r["ID"] if "ID" in r else r["Rep_ID"]: r for r in csv.DictReader(fh, delimiter="\t")}

    def test_same_name_two_ids(self):
        viruses = self.rows("viruses.tsv")
        self.assertEqual(sorted(viruses), ["MGYV1", "MGYV2", "MGYV5"])
        self.assertEqual(viruses["MGYV1"]["Biome"], "soil")
        self.assertEqual(viruses["MGYV5"]["Biome"], "lake")
        self.assertEqual(viruses["MGYV5"]["Cluster_rep"], "MGYV1")
        self.assertEqual(viruses["MGYV5"]["Original_ID"], "MGYG1_1|viral_sequence")
        self.assertEqual(viruses["MGYV1"]["Source_accession"], "GCA1")
        self.assertEqual(viruses["MGYV1"]["genomad_lineage"], "Viruses;Caudoviricetes")

    def test_cluster_stats(self):
        reps = self.rows("reps.tsv")
        self.assertEqual(sorted(reps), ["MGYV1", "MGYV2"])
        self.assertEqual(reps["MGYV1"]["cluster_size"], "2")
        self.assertEqual(reps["MGYV1"]["Rep_biome"], "soil")
        self.assertEqual(reps["MGYV1"]["Rep_genomad_lineage"], "Viruses;Caudoviricetes")

    def test_plasmids(self):
        plasmids = self.rows("plasmids.tsv")
        self.assertEqual(list(plasmids), ["MGYV6"])
        self.assertEqual(plasmids["MGYV6"]["Cluster_rep"], "MGYV6")


if __name__ == "__main__":
    unittest.main()
