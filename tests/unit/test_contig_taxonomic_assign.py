#!/usr/bin/env python3
"""Unit tests for contig_taxonomic_assign.py, with a small NCBI taxonomy instead of ete3."""

import os
import sys
import tempfile
import unittest
from pathlib import Path

BIN_DIR = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

import contig_taxonomic_assign as cta

# taxid: (name, rank, parent)
TAXONOMY = {
    1: ("root", "no rank", None),
    10: ("Viruses", "superkingdom", 1),
    20: ("Heunggongvirae", "kingdom", 10),
    30: ("Uroviricota", "phylum", 20),
    40: ("Caudoviricetes", "class", 30),
    50: ("Caudovirales", "order", 40),          # deprecated in ViPhOGs version 4
    60: ("Fam1viridae", "family", 50),
    61: ("Fam2viridae", "family", 50),
    70: ("Gen1virus", "genus", 60),
    71: ("Gen2virus", "genus", 60),
    72: ("Gen3virus", "genus", 61),
}

FACTORS = (
    "name,a,avg_cds,std_cds,b,factor\n"
    "Gen1virus,x,10,2,y,1.0\n"
    "Gen2virus,x,10,2,y,1.0\n"
    "Fam1viridae,x,10,2,y,1.0\n"
    "Fam2viridae,x,3,0,y,1.0\n"
)

HEADER = "Contig\tCDS_ID\tStart\tEnd\tDirection\tBest_hit\tAbs_Evalue_exp\tLabel\n"


def rows(contig, labels, no_hits=0):
    out = [f"{contig}\t{contig}_{i}\t1\t2\t+\tViPhOG{i}\t10\t{label}\n" for i, label in enumerate(labels)]
    out += [f"{contig}\t{contig}_n{i}\t1\t2\t+\tNo hit\t\t\n" for i in range(no_hits)]
    return "".join(out)


class FakeNCBI:
    """The ete3 NCBITaxa methods the script uses, counting the lookups."""

    def __init__(self):
        self.calls = 0
        self.by_name = {name: taxid for taxid, (name, _, _) in TAXONOMY.items()}

    def get_name_translator(self, names):
        self.calls += 1
        return {n: [self.by_name[n]] for n in names if n in self.by_name}

    def get_lineage(self, taxid):
        self.calls += 1
        lineage = []
        while taxid is not None:
            lineage.append(taxid)
            taxid = TAXONOMY[taxid][2]
        return lineage[::-1]

    def get_rank(self, taxids):
        self.calls += 1
        return {t: TAXONOMY[t][1] for t in taxids}

    def get_taxid_translator(self, taxids):
        self.calls += 1
        return {t: TAXONOMY[t][0] for t in taxids}


class TestAssign(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cwd = os.getcwd()
        os.chdir(self.tmp.name)
        Path("factors.csv").write_text(FACTORS)

    def tearDown(self):
        os.chdir(self.cwd)
        self.tmp.cleanup()

    def run_script(self, table, *extra):
        Path("annotation.tsv").write_text(HEADER + table)
        self.ncbi = FakeNCBI()
        cta.main(["-i", "annotation.tsv", "-f", "factors.csv", "-o", "out", *extra], ncbi=self.ncbi)
        lines = Path("out/annotation_taxonomy.tsv").read_text().splitlines()
        self.assertEqual(lines[0], "contig_ID\ttaxonomy")
        return dict(line.split("\t") for line in lines[1:]), [line.split("\t")[0] for line in lines[1:]]

    def test_genus_over_threshold(self):
        out, _ = self.run_script(rows("C1", ["Gen1virus"] * 4 + ["Gen2virus"], no_hits=3))
        self.assertEqual(out["C1"], "Viruses;Heunggongvirae;Uroviricota;;Caudoviricetes;Caudovirales;;Fam1viridae;;Gen1virus")

    def test_family_when_genera_disagree(self):
        # genus: 2 x Gen1 / 2 x Gen2, both below 0.6 and tied -> genus empty; family Fam1 100%
        out, _ = self.run_script(rows("C1", ["Gen1virus", "Gen1virus", "Gen2virus", "Gen2virus"]))
        self.assertEqual(out["C1"], "Viruses;Heunggongvirae;Uroviricota;;Caudoviricetes;Caudovirales;;Fam1viridae;;")

    def test_no_hits(self):
        out, _ = self.run_script(rows("C1", [], no_hits=5))
        self.assertEqual(out["C1"], ";" * 9)

    def test_cds_size_filter(self):
        # Gen3virus and Caudovirales have no factors: never assigned over the threshold.
        # Fam2viridae: at most 3 + 2 x 0 proteins
        out, _ = self.run_script(rows("C4", ["Gen3virus"] * 4) + rows("C3", ["Gen3virus"] * 3))
        self.assertEqual(out["C4"], ";" * 9)
        self.assertEqual(out["C3"], "Viruses;Heunggongvirae;Uroviricota;;Caudoviricetes;Caudovirales;;Fam2viridae;;")

    def test_version4_excludes_deprecated_order(self):
        table = rows("C1", ["Gen1virus"] * 3)
        with_order, _ = self.run_script(table)
        without_order, _ = self.run_script(table, "--version4")
        self.assertIn("Caudovirales", with_order["C1"])
        # the deprecated order is not counted, but the lineage of the assigned genus is the NCBI lineage
        self.assertEqual(without_order["C1"], with_order["C1"])

    def test_unknown_label_and_missing_values(self):
        out, _ = self.run_script(rows("C1", ["Gen1virus", "Gen1virus", "Unknownvirus", "NA"]))
        self.assertTrue(out["C1"].endswith("Gen1virus"))

    def test_file_order_and_one_lookup_per_label(self):
        table = "".join(rows(f"C{i}", ["Gen1virus", "Gen2virus", "Gen1virus"]) for i in range(50))
        _, order = self.run_script(table)
        self.assertEqual(order, [f"C{i}" for i in range(50)])
        # two labels and one assigned taxon, each resolved once, whatever the number of contigs
        self.assertLessEqual(self.ncbi.calls, 12)

    def test_factor_file_required(self):
        Path("annotation.tsv").write_text(HEADER + rows("C1", ["Gen1virus"]))
        with self.assertRaises(SystemExit):
            cta.main(["-i", "annotation.tsv", "-o", "out"], ncbi=FakeNCBI())


if __name__ == "__main__":
    unittest.main()
