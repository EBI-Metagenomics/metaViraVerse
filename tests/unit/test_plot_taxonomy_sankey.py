#!/usr/bin/env python3
"""
Unit tests for plot_taxonomy_sankey.py
"""

import unittest
import tempfile
import os
import sys
from pathlib import Path

# Add bin directory to path to import the script
bin_dir = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(bin_dir))

import bin.plot_taxonomy_sankey as plot_taxonomy_sankey


class TestParseTaxonomy(unittest.TestCase):
    """Test parse_taxonomy function"""

    def test_parse_complete_taxonomy(self):
        """Test parsing complete taxonomy string"""
        taxonomy_str = "Viruses;Heunggongvirae;Uroviricota;Caudoviricetes"
        levels = ['realm', 'kingdom', 'phylum', 'class']
        result = plot_taxonomy_sankey.parse_taxonomy(taxonomy_str, levels)
        self.assertEqual(result, ['Viruses', 'Heunggongvirae', 'Uroviricota', 'Caudoviricetes'])

    def test_parse_taxonomy_with_missing_levels(self):
        """Test parsing taxonomy with missing levels"""
        taxonomy_str = "Viruses;;Uroviricota;;Caudovirales"
        levels = ['realm', 'kingdom', 'phylum', 'class', 'order']
        result = plot_taxonomy_sankey.parse_taxonomy(taxonomy_str, levels)
        self.assertEqual(result[0], 'Viruses')
        self.assertEqual(result[1], 'Unclassified')
        self.assertEqual(result[2], 'Uroviricota')
        self.assertIn('unclassified_class', result[3])
        self.assertEqual(result[4], 'Caudovirales')

    def test_parse_taxonomy_empty(self):
        """Test parsing empty taxonomy"""
        result = plot_taxonomy_sankey.parse_taxonomy("")
        self.assertEqual(result, [])

    def test_parse_taxonomy_na(self):
        """Test parsing NA taxonomy"""
        result = plot_taxonomy_sankey.parse_taxonomy("NA")
        self.assertEqual(result, [])

    def test_parse_taxonomy_unclassified_only(self):
        """Test parsing taxonomy with only 'unclassified'"""
        result = plot_taxonomy_sankey.parse_taxonomy("unclassified")
        self.assertEqual(result, ['Unclassified'])

    def test_parse_taxonomy_unclassified_case_insensitive(self):
        """Test parsing taxonomy with 'Unclassified' in different cases"""
        result1 = plot_taxonomy_sankey.parse_taxonomy("UNCLASSIFIED")
        result2 = plot_taxonomy_sankey.parse_taxonomy("UnClassified")
        self.assertEqual(result1, ['Unclassified'])
        self.assertEqual(result2, ['Unclassified'])

    def test_parse_taxonomy_url_encoded(self):
        """Test parsing URL-encoded taxonomy"""
        taxonomy_str = "Viruses%3BHeunggongvirae"
        result = plot_taxonomy_sankey.parse_taxonomy(taxonomy_str)
        self.assertIsInstance(result, list)

    def test_parse_taxonomy_trailing_empty(self):
        """Test parsing taxonomy with trailing empty and dash values"""
        taxonomy_str = "Viruses;Heunggongvirae;;-"
        result = plot_taxonomy_sankey.parse_taxonomy(taxonomy_str)
        self.assertEqual(len(result), 2)

    def test_parse_taxonomy_first_level_empty(self):
        """Test parsing taxonomy with empty first level"""
        taxonomy_str = ";Heunggongvirae;Uroviricota"
        levels = ['realm', 'kingdom', 'phylum']
        result = plot_taxonomy_sankey.parse_taxonomy(taxonomy_str, levels)
        self.assertEqual(result[0], 'Unclassified')

    def test_parse_taxonomy_lineage_reverted(self):
        """Test parsing taxonomy with reverted lineage"""
        taxonomy_str = "Caudoviricetes;Uroviricota;Heunggongvirae;Viruses"
        levels = ['realm', 'kingdom', 'phylum', 'class']
        result = plot_taxonomy_sankey.parse_taxonomy(taxonomy_str, levels, lineage_reverted=True)
        self.assertEqual(result[0], 'Viruses')
        self.assertEqual(result[-1], 'Caudoviricetes')


class TestBuildSankeyData(unittest.TestCase):
    """Test build_sankey_data function"""

    def test_build_sankey_simple(self):
        """Test building Sankey data with simple taxonomy"""
        taxonomy_list = [
            (10, ['Viruses', 'Heunggongvirae']),
            (5, ['Viruses', 'Duplodnaviria'])
        ]
        nodes, links = plot_taxonomy_sankey.build_sankey_data(taxonomy_list)

        self.assertIn('All sequences', nodes)
        self.assertIn('Viruses', nodes)
        self.assertIn('Heunggongvirae', nodes)
        self.assertIn('Duplodnaviria', nodes)

        total_from_root = sum(link['value'] for link in links if nodes[link['source']] == 'All sequences')
        self.assertEqual(total_from_root, 15)

    def test_build_sankey_empty_taxonomy(self):
        """Test building Sankey data with empty taxonomy list"""
        taxonomy_list = []
        nodes, links = plot_taxonomy_sankey.build_sankey_data(taxonomy_list)

        self.assertEqual(len(nodes), 1)
        self.assertEqual(nodes[0], 'All sequences')
        self.assertEqual(len(links), 0)

    def test_build_sankey_duplicate_names_different_paths(self):
        """Test building Sankey data with same taxon name in different paths"""
        taxonomy_list = [
            (5, ['Viruses', 'GroupA', 'Family1']),
            (3, ['Viruses', 'GroupB', 'Family1'])
        ]
        nodes, links = plot_taxonomy_sankey.build_sankey_data(taxonomy_list)

        family1_count = sum(1 for node in nodes if 'Family1' in node)
        self.assertEqual(family1_count, 2)

    def test_build_sankey_counts_accumulate(self):
        """Test that counts accumulate correctly"""
        taxonomy_list = [
            (10, ['Viruses', 'Heunggongvirae', 'Uroviricota']),
            (5, ['Viruses', 'Heunggongvirae', 'Uroviricota'])
        ]
        nodes, links = plot_taxonomy_sankey.build_sankey_data(taxonomy_list)

        viruses_idx = nodes.index('Viruses')
        viruses_out = sum(link['value'] for link in links if link['source'] == viruses_idx)
        self.assertEqual(viruses_out, 15)


class TestReadTsvTaxonomy(unittest.TestCase):
    """Test read_tsv_taxonomy function"""

    def test_read_standard_tsv(self):
        """Test reading standard TSV format"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("viral_sequence_name\ttaxonomy\n")
            f.write("seq1\tViruses;Heunggongvirae\n")
            f.write("seq2\tViruses;Duplodnaviria\n")
            tsv_file = f.name

        try:
            result = plot_taxonomy_sankey.read_tsv_taxonomy(tsv_file, 'taxonomy')
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0][0], 1)
            self.assertIn('Viruses', result[0][1])
        finally:
            os.unlink(tsv_file)

    def test_read_krona_format(self):
        """Test reading Krona format"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("10\tViruses\tHeunggongvirae\n")
            f.write("5\tViruses\tDuplodnaviria\n")
            tsv_file = f.name

        try:
            result = plot_taxonomy_sankey.read_tsv_taxonomy(tsv_file)
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0][0], 10)
            self.assertEqual(result[1][0], 5)
        finally:
            os.unlink(tsv_file)

    def test_read_tsv_with_na(self):
        """Test reading TSV with NA values"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("viral_sequence_name\ttaxonomy\n")
            f.write("seq1\tViruses;Heunggongvirae\n")
            f.write("seq2\tNA\n")
            tsv_file = f.name

        try:
            result = plot_taxonomy_sankey.read_tsv_taxonomy(tsv_file, 'taxonomy')
            self.assertEqual(len(result), 1)
        finally:
            os.unlink(tsv_file)

    def test_read_empty_tsv(self):
        """Test reading empty TSV file"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("")
            tsv_file = f.name

        try:
            with self.assertRaises(SystemExit):
                plot_taxonomy_sankey.read_tsv_taxonomy(tsv_file)
        finally:
            os.unlink(tsv_file)

    def test_read_krona_format_with_empty_taxonomy(self):
        """Test reading Krona format with empty taxonomy"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("10\t\t\n")
            f.write("5\tViruses\tHeunggongvirae\n")
            tsv_file = f.name

        try:
            result = plot_taxonomy_sankey.read_tsv_taxonomy(tsv_file)
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0][1], 'Unclassified')
        finally:
            os.unlink(tsv_file)


class TestReadGffTaxonomy(unittest.TestCase):
    """Test read_gff_taxonomy function"""

    def test_read_gff_taxonomy(self):
        """Test reading taxonomy from GFF file"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gff') as f:
            f.write("##gff-version 3\n")
            f.write("seq1\tVIRify\tregion\t1\t1000\t.\t+\t.\tID=seq1;taxonomy=Viruses;Heunggongvirae\n")
            f.write("seq2\tVIRify\tregion\t1\t2000\t.\t+\t.\tID=seq2;taxonomy=Viruses;Duplodnaviria\n")
            gff_file = f.name

        try:
            result = plot_taxonomy_sankey.read_gff_taxonomy(gff_file)
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0][0], 1)
        finally:
            os.unlink(gff_file)

    def test_read_gff_with_na_taxonomy(self):
        """Test reading GFF with NA taxonomy"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gff') as f:
            f.write("##gff-version 3\n")
            f.write("seq1\tVIRify\tregion\t1\t1000\t.\t+\t.\tID=seq1;taxonomy=Viruses;Heunggongvirae\n")
            f.write("seq2\tVIRify\tregion\t1\t2000\t.\t+\t.\tID=seq2;taxonomy=NA\n")
            gff_file = f.name

        try:
            result = plot_taxonomy_sankey.read_gff_taxonomy(gff_file)
            self.assertEqual(len(result), 1)
        finally:
            os.unlink(gff_file)

    def test_read_gff_with_comments(self):
        """Test reading GFF with comment lines"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gff') as f:
            f.write("##gff-version 3\n")
            f.write("# This is a comment\n")
            f.write("seq1\tVIRify\tregion\t1\t1000\t.\t+\t.\tID=seq1;taxonomy=Viruses\n")
            gff_file = f.name

        try:
            result = plot_taxonomy_sankey.read_gff_taxonomy(gff_file)
            self.assertEqual(len(result), 1)
        finally:
            os.unlink(gff_file)


if __name__ == '__main__':
    unittest.main()
