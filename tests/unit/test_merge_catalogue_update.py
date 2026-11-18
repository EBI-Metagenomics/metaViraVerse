#!/usr/bin/env python3
"""Unit tests for merge_catalogue_update.py."""

import gzip
import os
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import merge_catalogue_update as mcu

# Existing catalogue: representatives V1, V3 and V4; the update merged V4 into V3 and
# added the new cluster V7. Proteins keep their source IDs (G1_1 is on V1, ...).
OLD_GFF = (
    "##gff-version 3\n"
    "##sequence-region V1 1 100\n"
    "V1\tVIRify\tviral_sequence\t1\t100\t.\t.\t.\tID=V1\n"
    "V1\tProdigal\tCDS\t1\t60\t.\t+\t0\tID=G1_1\n"
    "V3\tVIRify\tprophage\t1\t100\t.\t.\t.\tID=V3\n"
    "V3\tProdigal\tCDS\t1\t60\t.\t+\t0\tID=G3_1\n"
    "V4\tVIRify\tviral_sequence\t1\t100\t.\t.\t.\tID=V4\n"
    "V4\tProdigal\tCDS\t1\t60\t.\t+\t0\tID=G4_1\n"
    "V4\tProdigal\tCDS\t61\t99\t.\t+\t0\tID=G4_2\n"
)
NEW_GFF = (
    "##gff-version 3\n"
    "V7\tVIRify\tviral_sequence\t1\t100\t.\t.\t.\tID=V7\n"
    "V7\tProdigal\tCDS\t1\t60\t.\t+\t0\tID=G7_1\n"
)


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.cwd = os.getcwd()
        os.chdir(self.dir)
        self.write("dropped_reps.txt", "V4\n")
        self.write("old.gff", OLD_GFF)
        mcu.main(["drop-ids", "--dropped-reps", "dropped_reps.txt", "--gff", "old.gff", "-o", "drop_ids.txt"])

    def tearDown(self):
        os.chdir(self.cwd)
        self.tmp.cleanup()

    def write(self, name, text, compress=False):
        path = self.dir / name
        if compress:
            with gzip.open(path, "wt") as fh:
                fh.write(text)
        else:
            path.write_text(text)
        return name

    def read(self, name):
        path = self.dir / name
        with (gzip.open(path, "rt") if name.endswith(".gz") else open(path)) as fh:
            return fh.read()

    def merge(self, fmt, old_text, new_text, out="out.txt", extra=(), compress=False):
        old = self.write("old_" + out, old_text, compress)
        new = self.write("new_" + out, new_text, compress)
        mcu.main(["merge", "--format", fmt, "--old", old, "--new", new, "--drop-ids", "drop_ids.txt", "-o", out, *extra])
        return self.read(out)


class TestDropIds(Base):
    def test_reps_and_their_proteins(self):
        self.assertEqual(self.read("drop_ids.txt").split(), ["G4_1", "G4_2", "V4"])

    def test_rep_missing_from_gff_fails(self):
        self.write("dropped_reps.txt", "V9\n")
        with self.assertRaises(SystemExit) as cm:
            mcu.main(["drop-ids", "--dropped-reps", "dropped_reps.txt", "--gff", "old.gff", "-o", "x.txt"])
        self.assertIn("not in", str(cm.exception))

    def test_nothing_dropped(self):
        self.write("dropped_reps.txt", "")
        mcu.main(["drop-ids", "--dropped-reps", "dropped_reps.txt", "--gff", "old.gff", "-o", "none.txt"])
        self.assertEqual(self.read("none.txt"), "")


class TestTables(Base):
    def test_tsv_by_column_name(self):
        out = self.merge(
            "tsv",
            "seq_name\tlength\ttaxonomy\nV1\t100\tA\nV3\t100\tB\nV4\t100\tC\n",
            "seq_name\tlength\ttaxonomy\nV7\t100\tD\n",
            extra=["--id-column", "seq_name"],
        )
        self.assertEqual(out, "seq_name\tlength\ttaxonomy\nV1\t100\tA\nV3\t100\tB\nV7\t100\tD\n")

    def test_header_names_are_stripped(self):
        # ViPhOGs taxonomy: the ID column header has leading spaces
        out = self.merge("tsv", "    contig_ID\tgenus\nV4\tx\nV1\ty\n", "    contig_ID\tgenus\nV7\tz\n",
                         extra=["--id-column", "contig_ID"])
        self.assertEqual(out, "    contig_ID\tgenus\nV1\ty\nV7\tz\n")

    def test_tsv_by_index_unnamed_column(self):
        # BACPHLIP: the ID column has no name
        out = self.merge("tsv", "\tVirulent\tTemperate\nV1\t0.9\t0.1\nV4\t0.2\t0.8\n",
                         "\tVirulent\tTemperate\nV7\t0.5\t0.5\n", extra=["--id-index", "0"], compress=True)
        self.assertEqual(out, "\tVirulent\tTemperate\nV1\t0.9\t0.1\nV7\t0.5\t0.5\n")

    def test_protein_rows_of_dropped_reps(self):
        # hAMRonization: protein IDs; G4_1 is a protein of the dropped V4
        out = self.merge("tsv", "input_sequence_id\tgene\nG1_1\tblaA\nG4_1\tblaB\n", "input_sequence_id\tgene\nG7_1\tblaC\n",
                         extra=["--id-column", "input_sequence_id"])
        self.assertEqual(out, "input_sequence_id\tgene\nG1_1\tblaA\nG7_1\tblaC\n")

    def test_csv_repeated_headers_and_quotes(self):
        out = self.merge(
            "csv",
            'Virus,Host genus,Confidence score\nV1,"a;b, c",90\nVirus,Host genus,Confidence score\nV4,x,80\n',
            "Virus,Host genus,Confidence score\nV7,y,70\nV7,z,60\n",
            extra=["--id-column", "Virus"],
        )
        self.assertEqual(out, 'Virus,Host genus,Confidence score\nV1,"a;b, c",90\nV7,y,70\nV7,z,60\n')

    def test_empty_old_file_uses_new_header(self):
        out = self.merge("tsv", "", "seq_name\ttaxonomy\nV7\tD\n", extra=["--id-column", "seq_name"])
        self.assertEqual(out, "seq_name\ttaxonomy\nV7\tD\n")

    def test_column_mismatch_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.merge("tsv", "seq_name\ttaxonomy\nV1\tA\n", "seq_name\tscore\ttaxonomy\nV7\t1\tD\n",
                       extra=["--id-column", "seq_name"])
        self.assertIn("columns", str(cm.exception))

    def test_missing_id_column_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.merge("tsv", "a\tb\nV1\t1\n", "a\tb\nV7\t1\n", extra=["--id-column", "seq_name"])
        self.assertIn("no column", str(cm.exception))

    def test_existing_id_annotated_again_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.merge("tsv", "seq_name\tt\nV1\tA\n", "seq_name\tt\nV1\tB\n", extra=["--id-column", "seq_name"])
        self.assertIn("already in the existing", str(cm.exception))


class TestGff(Base):
    def test_merge(self):
        out = self.merge("gff", OLD_GFF, NEW_GFF, out="merged.gff.gz", compress=True)
        lines = out.splitlines()
        self.assertEqual(lines[0], "##gff-version 3")
        self.assertEqual(sum(line.startswith("##gff-version") for line in lines), 1)
        self.assertIn("##sequence-region V1 1 100", lines)
        self.assertFalse(any(line.startswith("V4") for line in lines))
        self.assertEqual([line.split("\t")[0] for line in lines if not line.startswith("#")],
                         ["V1", "V1", "V3", "V3", "V7", "V7"])


class TestFasta(Base):
    def test_merge_proteins(self):
        out = self.merge("fasta", ">G1_1 hypothetical\nMKV\n>G4_1 x\nMAA\n>G4_2\nMCC\n", ">G7_1 new\nMDD\nEE\n")
        self.assertEqual(out, ">G1_1 hypothetical\nMKV\n>G7_1 new\nMDD\nEE\n")

    def test_duplicated_id_fails(self):
        with self.assertRaises(SystemExit) as cm:
            self.merge("fasta", ">V1\nACGT\n", ">V7\nACGT\n>V7\nACGA\n")
        self.assertIn("more than once", str(cm.exception))


class TestTbl(Base):
    def test_merge(self):
        header = "# target name accession query name accession\n#---- ----\n"
        out = self.merge("tbl", header + "G1_1 - PF1 PF1.1\nG4_1 - PF2 PF2.1\n# Program: hmmsearch\n",
                         header + "G7_1 - PF3 PF3.1\n# Program: hmmsearch\n")
        self.assertEqual(out, header + "G1_1 - PF1 PF1.1\nG7_1 - PF3 PF3.1\n")


class TestReps(Base):
    def test_reps(self):
        self.write("clusters.tsv", "object\tcluster\nV1\tV1\nV2\tV1\nV7\tV7\n")
        mcu.main(["reps", "--clusters-tsv", "clusters.tsv", "-o", "reps.tsv"])
        self.assertEqual(self.read("reps.tsv"), "V1\nV7\n")


class TestConcat(Base):
    def test_gff(self):
        self.write("a.gff.gz", OLD_GFF, compress=True)
        self.write("b.gff", NEW_GFF)
        mcu.main(["concat", "--format", "gff", "-o", "all.gff", "a.gff.gz", "b.gff"])
        lines = self.read("all.gff").splitlines()
        self.assertEqual(sum(line.startswith("##gff-version") for line in lines), 1)
        self.assertIn("##sequence-region V1 1 100", lines)
        self.assertEqual([line.split("\t")[0] for line in lines if not line.startswith("#")],
                         ["V1", "V1", "V3", "V3", "V4", "V4", "V4", "V7", "V7"])

    def test_fasta(self):
        self.write("a.fna.gz", ">V1 x\nAC\nGT\n", compress=True)
        self.write("b.fna", ">V7\nGG")
        mcu.main(["concat", "--format", "fasta", "--unique-ids", "-o", "all.fna.gz", "a.fna.gz", "b.fna"])
        self.assertEqual(self.read("all.fna.gz"), ">V1 x\nAC\nGT\n>V7\nGG\n")

    def test_fasta_duplicated_id(self):
        self.write("a.fna", ">V1\nAC\n")
        self.write("b.fna", ">V1\nGG\n")
        with self.assertRaises(SystemExit) as cm:
            mcu.main(["concat", "--format", "fasta", "--unique-ids", "-o", "all.fna", "a.fna", "b.fna"])
        self.assertIn("more than one", str(cm.exception))
        # proteins keep their source IDs: duplicates are allowed without --unique-ids
        mcu.main(["concat", "--format", "fasta", "-o", "all.faa", "a.fna", "b.fna"])


if __name__ == "__main__":
    unittest.main()
