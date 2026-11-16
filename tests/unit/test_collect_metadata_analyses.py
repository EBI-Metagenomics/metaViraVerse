#!/usr/bin/env python3
"""Unit tests for collect_metadata.py --analyses_metadata (sequences of MGnify analyses, ERZ...)."""

import csv
import gzip
import os
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import collect_metadata as cm

META = (
    "description\ttype\tsource_of_prediction\tbiomes\tsequence_length\trrna\tsequence_sha256\n"
    "MGYG1_1|viral_sequence\tgenome\tMGnify\troot:Host-associated\t100\tno\th1\n"
    "ERZ2772503|viral_sequence-1:500\tmetagenome\tMGnify\troot:Environmental:Aquatic\t500\tno\th2\n"
    "ERZ2772503|prophage-10:90\tmetagenome\tMGnify\t\t80\tno\th3\n"
    "ERZ999|viral_sequence\tmetagenome\tMGnify\troot:Engineered\t70\tno\th4\n"
    "ERZ2772503|plasmid-1:40\tmetagenome\tMGnify\troot:Environmental:Aquatic\t40\tno\th5\n"
)
MAP = (
    "original\ttemporary\tshort\tbiome\ttype\tdefinition\n"
    "MGYG1_1|viral_sequence\tMGYV1\tMGYG1_1|viral_sequence\tx\tgenome\tvirus\n"
    "ERZ2772503|viral_sequence-1:500\tMGYV2\tERZ2772503|viral_sequence-1:500\tx\tmetagenome\tvirus\n"
    "ERZ2772503|prophage-10:90\tMGYV3\tERZ2772503|prophage-10:90\tx\tmetagenome\tprophage\n"
    "ERZ999|viral_sequence\tMGYV4\tERZ999|viral_sequence\tx\tmetagenome\tvirus\n"
    "ERZ2772503|plasmid-1:40\tMGYV5\tERZ2772503|plasmid-1:40\tx\tmetagenome\tplasmid\n"
)
CLUSTERS = "object\tcluster\nMGYV1\tMGYV1\nMGYV2\tMGYV2\nMGYV3\tMGYV2\nMGYV4\tMGYV4\n"
PLASMID_CLUSTERS = "object\tcluster\nMGYV5\tMGYV5\n"
CATALOGUES = (
    "Genome\tGenome_accession\tLineage\tSample_accession\tStudy_accession\n"
    "MGYG1\tGCA_1\td__Bacteria\tSAMEA1\tPRJEB1\n"
)
# NA sample for ERZ999; ERZ2772503 is listed twice (concatenated exports), once with an empty sample
ANALYSES = (
    "accession\tstudy_id\tsample_id\tbiome_id\n"
    "ERZ2772503\tPRJEB2\t\troot:Environmental:Aquatic:Marine\n"
    "accession\tstudy_id\tsample_id\tbiome_id\n"
    "ERZ2772503\tPRJEB2\tSAMEA2\troot:Environmental:Aquatic:Marine\n"
    "ERZ999\tPRJEB3\tNA\troot:Engineered\n"
)


class TestAnalysesMetadata(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cwd = os.getcwd()
        os.chdir(self.tmp.name)
        for name, text in [("meta.tsv", META), ("map.tsv", MAP), ("clusters.tsv", CLUSTERS),
                           ("plasmid_clusters.tsv", PLASMID_CLUSTERS), ("catalogues.tsv", CATALOGUES)]:
            Path(name).write_text(text)
        with gzip.open("analyses.tsv.gz", "wt") as fh:
            fh.write(ANALYSES)

    def tearDown(self):
        os.chdir(self.cwd)
        self.tmp.cleanup()

    def run_script(self, analyses=True):
        sys.argv = [
            "collect_metadata.py", "--combined_meta", "meta.tsv", "--viruses_cluster", "clusters.tsv",
            "--plasmids_cluster", "plasmid_clusters.tsv", "--additional_metadata", "catalogues.tsv",
            "--map", "map.tsv", "--output_viruses", "viruses.tsv", "--output_plasmids", "plasmids.tsv",
            "--output_reps", "reps.tsv",
        ] + (["--analyses_metadata", "analyses.tsv.gz"] if analyses else [])
        cm.main()

    @staticmethod
    def rows(path):
        with open(path) as fh:
            return {r["ID"]: r for r in csv.DictReader(fh, delimiter="\t")}

    def test_analysis_sequences_get_study_and_sample(self):
        self.run_script()
        rows = self.rows("viruses.tsv")
        for mgyv in ("MGYV2", "MGYV3"):
            self.assertEqual(rows[mgyv]["Source_accession"], "ERZ2772503")
            self.assertEqual(rows[mgyv]["Source_project"], "PRJEB2")
            self.assertEqual(rows[mgyv]["Source_sample"], "SAMEA2")   # the non-empty row wins
            self.assertEqual(rows[mgyv]["Source_lineage"], "missing")

    def test_biome_from_analyses_only_when_missing(self):
        self.run_script()
        rows = self.rows("viruses.tsv")
        self.assertEqual(rows["MGYV2"]["Biome"], "root:Environmental:Aquatic")
        self.assertEqual(rows["MGYV3"]["Biome"], "root:Environmental:Aquatic:Marine")

    def test_na_is_missing(self):
        self.run_script()
        rows = self.rows("viruses.tsv")
        self.assertEqual(rows["MGYV4"]["Source_sample"], "missing")
        self.assertEqual(rows["MGYV4"]["Source_project"], "PRJEB3")

    def test_catalogue_genomes_unchanged(self):
        self.run_script()
        row = self.rows("viruses.tsv")["MGYV1"]
        self.assertEqual((row["Source_accession"], row["Source_sample"], row["Source_project"], row["Source_lineage"]),
                         ("GCA_1", "SAMEA1", "PRJEB1", "d__Bacteria"))

    def test_plasmids(self):
        self.run_script()
        row = self.rows("plasmids.tsv")["MGYV5"]
        self.assertEqual((row["Source_sample"], row["Source_project"]), ("SAMEA2", "PRJEB2"))

    def test_without_analyses_metadata(self):
        self.run_script(analyses=False)
        row = self.rows("viruses.tsv")["MGYV2"]
        self.assertEqual((row["Source_sample"], row["Source_project"]), ("missing", "missing"))


if __name__ == "__main__":
    unittest.main()
