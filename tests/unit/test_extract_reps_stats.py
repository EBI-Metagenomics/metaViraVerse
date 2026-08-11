#!/usr/bin/env python3
"""
Unit tests for extract_reps_stats.py

extract_reps_stats.py now only extracts GFF records, protein IDs and a plain
rep-ID list for cluster representatives. Taxonomy/checkV/cluster-level stats
(previously tested here as parse_taxonomy/generate_krona_file/calculate_mean_genes)
have moved to generate_taxonomy_table.py, plot_taxonomy_sankey.py and
collect_metadata.py, and attribute parsing has moved to utils.py.
"""

import unittest
import tempfile
import os
import sys
from pathlib import Path

# Add bin directory to path to import the script
bin_dir = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(bin_dir))

import extract_reps_stats as extract_reps_stats


class TestReadMapfile(unittest.TestCase):
    """Test read_mapfile function"""

    def test_read_mapfile(self):
        """Test reading a mapping file"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("original_name1\ttemp_name1\n")
            f.write("original_name2\ttemp_name2\n")
            f.write("original_name3\ttemp_name3\n")
            mapfile = f.name

        try:
            result = extract_reps_stats.read_mapfile(mapfile)
            self.assertEqual(result['temp_name1'], 'original_name1')
            self.assertEqual(result['temp_name2'], 'original_name2')
            self.assertEqual(result['temp_name3'], 'original_name3')
            self.assertEqual(len(result), 3)
        finally:
            os.unlink(mapfile)

    def test_read_empty_mapfile(self):
        """Test reading an empty mapping file"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            mapfile = f.name

        try:
            result = extract_reps_stats.read_mapfile(mapfile)
            self.assertEqual(result, {})
        finally:
            os.unlink(mapfile)


class TestReadClusterReps(unittest.TestCase):
    """Test read_cluster_reps function"""

    def _write(self, lines):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("\n".join(lines) + "\n")
            return f.name

    def test_read_cluster_reps_blastn_format(self):
        """blastn-style cluster file: rep is the first column, no header"""
        cluster_file = self._write([
            "repA\tmemberA1\tmemberA2",
            "repB\tmemberB1",
        ])

        try:
            cluster_reps, rep_original_names = extract_reps_stats.read_cluster_reps(cluster_file)
            self.assertEqual(cluster_reps, ['repA', 'repB'])
            self.assertEqual(rep_original_names, {'repA': 'repA', 'repB': 'repB'})
        finally:
            os.unlink(cluster_file)

    def test_read_cluster_reps_vclust_format(self):
        """vclust-style cluster file: header contains 'object'/'cluster', rep is the second column"""
        cluster_file = self._write([
            "cluster\tobject",
            "clusterA\trepA",
            "clusterA\trepA",
            "clusterB\trepB",
        ])

        try:
            cluster_reps, rep_original_names = extract_reps_stats.read_cluster_reps(cluster_file)
            self.assertEqual(cluster_reps, ['repA', 'repB'])
            self.assertEqual(rep_original_names, {'repA': 'repA', 'repB': 'repB'})
        finally:
            os.unlink(cluster_file)

    def test_read_cluster_reps_with_mapping(self):
        """Original (temporary) rep names are translated through the mapping"""
        cluster_file = self._write(["temp_rep1\tmember1"])
        mapping = {'temp_rep1': 'original_rep1'}

        try:
            cluster_reps, rep_original_names = extract_reps_stats.read_cluster_reps(cluster_file, mapping)
            self.assertEqual(cluster_reps, ['original_rep1'])
            self.assertEqual(rep_original_names['original_rep1'], 'temp_rep1')
        finally:
            os.unlink(cluster_file)

    def test_read_cluster_reps_dedups_repeated_rep(self):
        """The same representative appearing on multiple lines is only kept once"""
        cluster_file = self._write([
            "rep1\tmember1",
            "rep1\tmember2",
            "rep2\tmember3",
        ])

        try:
            cluster_reps, _ = extract_reps_stats.read_cluster_reps(cluster_file)
            self.assertEqual(cluster_reps, ['rep1', 'rep2'])
        finally:
            os.unlink(cluster_file)

    def test_read_cluster_reps_skips_empty_lines(self):
        """Blank lines in the cluster file are ignored"""
        cluster_file = self._write([
            "rep1\tmember1",
            "",
            "rep2\tmember2",
        ])

        try:
            cluster_reps, _ = extract_reps_stats.read_cluster_reps(cluster_file)
            self.assertEqual(cluster_reps, ['rep1', 'rep2'])
        finally:
            os.unlink(cluster_file)


class TestExtractProteins(unittest.TestCase):
    """Test extract_proteins function"""

    def test_extract_proteins_from_cds_lines(self):
        """Only CDS feature lines contribute protein IDs"""
        lines = [
            "seq1\tVIRify\tregion\t1\t100\t.\t+\t.\tID=seq1\n",
            "seq1\tVIRify\tCDS\t1\t50\t.\t+\t.\tID=seq1_1;Name=protA\n",
            "seq1\tVIRify\tCDS\t60\t90\t.\t+\t.\tID=seq1_2\n",
        ]
        result = extract_reps_stats.extract_proteins(lines)
        self.assertEqual(result, {'seq1_1', 'seq1_2'})

    def test_extract_proteins_skips_cds_without_id(self):
        """CDS lines with no ID attribute don't contribute a protein ID"""
        lines = [
            "seq1\tVIRify\tCDS\t1\t50\t.\t+\t.\tName=protA\n",
        ]
        result = extract_reps_stats.extract_proteins(lines)
        self.assertEqual(result, set())

    def test_extract_proteins_skips_short_lines(self):
        """Lines with fewer than 9 columns are ignored"""
        lines = [
            "seq1\tVIRify\tCDS\t1\t50\n",
        ]
        result = extract_reps_stats.extract_proteins(lines)
        self.assertEqual(result, set())

    def test_extract_proteins_empty_input(self):
        """No lines produces no proteins"""
        self.assertEqual(extract_reps_stats.extract_proteins([]), set())


class TestExtractGffAndProteins(unittest.TestCase):
    """Test extract_gff_and_proteins function (integration)"""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _write(self, name, content):
        path = os.path.join(self.tmpdir, name)
        with open(path, 'w') as f:
            f.write(content)
        return path

    def test_extract_gff_and_proteins_writes_expected_outputs(self):
        viral_list_file = self._write("clusters.tsv", "rep1\tmember1\n")
        gff_file = self._write(
            "input.gff",
            "seq1\tVIRify\tregion\t1\t100\t.\t+\t.\tID=rep1\n"
            "seq1\tVIRify\tCDS\t1\t50\t.\t+\t.\tID=rep1_1\n"
            "seq2\tVIRify\tregion\t1\t80\t.\t+\t.\tID=member1\n"
            "seq2\tVIRify\tCDS\t1\t30\t.\t+\t.\tID=member1_1\n",
        )
        output_reps_list = os.path.join(self.tmpdir, "reps_list.tsv")
        output_reps_gff = os.path.join(self.tmpdir, "reps.gff")
        output_proteins = os.path.join(self.tmpdir, "proteins.tsv")

        extract_reps_stats.extract_gff_and_proteins(
            viral_list_file, gff_file, output_reps_list, output_reps_gff, output_proteins, mapfile=None
        )

        with open(output_reps_list) as f:
            self.assertEqual(f.read().splitlines(), ['rep1'])

        with open(output_reps_gff) as f:
            gff_content = f.read()
        self.assertIn("##gff-version 3", gff_content)
        self.assertIn("ID=rep1\n", gff_content)
        self.assertIn("ID=rep1_1\n", gff_content)
        self.assertNotIn("ID=member1", gff_content)

        with open(output_proteins) as f:
            self.assertEqual(f.read().splitlines(), ['rep1_1'])

    def test_extract_gff_and_proteins_skips_rep_with_no_gff_match(self):
        viral_list_file = self._write("clusters.tsv", "unknown_rep\tmember1\n")
        gff_file = self._write(
            "input.gff",
            "seq1\tVIRify\tregion\t1\t100\t.\t+\t.\tID=rep1\n",
        )
        output_reps_list = os.path.join(self.tmpdir, "reps_list.tsv")
        output_reps_gff = os.path.join(self.tmpdir, "reps.gff")
        output_proteins = os.path.join(self.tmpdir, "proteins.tsv")

        extract_reps_stats.extract_gff_and_proteins(
            viral_list_file, gff_file, output_reps_list, output_reps_gff, output_proteins, mapfile=None
        )

        with open(output_reps_gff) as f:
            self.assertEqual(f.read(), "##gff-version 3\n")
        with open(output_proteins) as f:
            self.assertEqual(f.read(), "")

    def test_extract_gff_and_proteins_uses_mapfile(self):
        # mapfile lines are "original\ttemp" (matching rename_contigs' map file), so
        # read_mapfile() indexes it as temp -> original.
        mapfile = self._write("map.tsv", "originalA\ttempA\n")
        # The cluster TSV (blastn format) is keyed by the renamed/temp sequence ID,
        # same as the GFF's ID= attributes below.
        viral_list_file = self._write("clusters.tsv", "tempA\tmember1\n")
        gff_file = self._write(
            "input.gff",
            "seq1\tVIRify\tregion\t1\t100\t.\t+\t.\tID=tempA\n"
            "seq1\tVIRify\tCDS\t1\t50\t.\t+\t.\tID=tempA_1\n",
        )
        output_reps_list = os.path.join(self.tmpdir, "reps_list.tsv")
        output_reps_gff = os.path.join(self.tmpdir, "reps.gff")
        output_proteins = os.path.join(self.tmpdir, "proteins.tsv")

        extract_reps_stats.extract_gff_and_proteins(
            viral_list_file, gff_file, output_reps_list, output_reps_gff, output_proteins, mapfile=mapfile
        )

        # cluster_reps is keyed by the mapped (original) ID, but the reps list is
        # written from rep_original_names, i.e. the raw (temp) cluster-file ID.
        with open(output_reps_list) as f:
            self.assertEqual(f.read().splitlines(), ['tempA'])
        with open(output_proteins) as f:
            self.assertEqual(f.read().splitlines(), ['tempA_1'])


if __name__ == '__main__':
    unittest.main()
