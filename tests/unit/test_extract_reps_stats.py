#!/usr/bin/env python3
"""
Unit tests for extract_reps_stats.py
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


class TestParseAttributes(unittest.TestCase):
    """Test parse_attributes function"""

    def test_parse_simple_attributes(self):
        """Test parsing simple GFF attributes"""
        attr_str = "ID=seq1;Name=viral_sequence;taxonomy=Viruses"
        result = extract_reps_stats.parse_attributes(attr_str)
        self.assertEqual(result['id'], 'seq1')
        self.assertEqual(result['name'], 'viral_sequence')
        self.assertEqual(result['taxonomy'], 'Viruses')

    def test_parse_attributes_with_spaces(self):
        """Test parsing attributes with spaces"""
        attr_str = "ID=seq1; Name=test ; taxonomy=Viruses;Heunggongvirae"
        result = extract_reps_stats.parse_attributes(attr_str)
        self.assertEqual(result['id'], 'seq1')
        self.assertEqual(result['name'], 'test')

    def test_parse_empty_attributes(self):
        """Test parsing empty attributes"""
        attr_str = ""
        result = extract_reps_stats.parse_attributes(attr_str)
        self.assertEqual(result, {})

    def test_parse_attributes_without_equals(self):
        """Test parsing attributes without equals sign"""
        attr_str = "ID=seq1;invalid;Name=test"
        result = extract_reps_stats.parse_attributes(attr_str)
        self.assertEqual(result['id'], 'seq1')
        self.assertEqual(result['name'], 'test')
        self.assertNotIn('invalid', result)


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


class TestReadClusterStructure(unittest.TestCase):
    """Test read_cluster_structure function"""

    def test_read_cluster_with_members(self):
        """Test reading cluster structure with members"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("rep1\tmember1\n")
            f.write("rep1\tmember2\n")
            f.write("rep2\tmember3\n")
            cluster_file = f.name

        try:
            reps, members, originals = extract_reps_stats.read_cluster_structure(cluster_file)
            self.assertEqual(reps, ['rep1', 'rep2'])
            self.assertEqual(members['rep1'], ['member1', 'member2'])
            self.assertEqual(members['rep2'], ['member3'])
            self.assertEqual(originals['rep1'], 'rep1')
            self.assertEqual(originals['rep2'], 'rep2')
        finally:
            os.unlink(cluster_file)

    def test_read_cluster_without_members(self):
        """Test reading cluster structure with only representatives"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("rep1\n")
            f.write("rep2\n")
            cluster_file = f.name

        try:
            reps, members, originals = extract_reps_stats.read_cluster_structure(cluster_file)
            self.assertEqual(reps, ['rep1', 'rep2'])
            self.assertEqual(members['rep1'], [])
            self.assertEqual(members['rep2'], [])
        finally:
            os.unlink(cluster_file)

    def test_read_cluster_with_mapping(self):
        """Test reading cluster structure with name mapping"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("temp_rep1\ttemp_member1\n")
            cluster_file = f.name

        mapping = {'temp_rep1': 'original_rep1', 'temp_member1': 'original_member1'}

        try:
            reps, members, originals = extract_reps_stats.read_cluster_structure(cluster_file, mapping)
            self.assertEqual(reps, ['original_rep1'])
            self.assertEqual(members['original_rep1'], ['original_member1'])
            self.assertEqual(originals['original_rep1'], 'temp_rep1')
        finally:
            os.unlink(cluster_file)

    def test_read_cluster_with_empty_lines(self):
        """Test reading cluster structure with empty lines"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("rep1\tmember1\n")
            f.write("\n")
            f.write("rep2\tmember2\n")
            cluster_file = f.name

        try:
            reps, members, originals = extract_reps_stats.read_cluster_structure(cluster_file)
            self.assertEqual(len(reps), 2)
            self.assertEqual(reps, ['rep1', 'rep2'])
        finally:
            os.unlink(cluster_file)


class TestCalculateMeanGenes(unittest.TestCase):
    """Test calculate_mean_genes function"""

    def test_calculate_mean_genes_simple(self):
        """Test calculating mean genes for simple cluster"""
        cluster_members = {
            'rep1': ['member1', 'member2']
        }
        all_results = {
            'rep1': {'checkv_viral_genes': '10'},
            'member1': {'checkv_viral_genes': '20'},
            'member2': {'checkv_viral_genes': '30'}
        }

        result = extract_reps_stats.calculate_mean_genes(cluster_members, all_results)
        self.assertEqual(result['rep1'], 20.0)

    def test_calculate_mean_genes_with_na(self):
        """Test calculating mean genes with NA values"""
        cluster_members = {
            'rep1': ['member1', 'member2']
        }
        all_results = {
            'rep1': {'checkv_viral_genes': '10'},
            'member1': {'checkv_viral_genes': 'NA'},
            'member2': {'checkv_viral_genes': '30'}
        }

        result = extract_reps_stats.calculate_mean_genes(cluster_members, all_results)
        self.assertEqual(result['rep1'], 20.0)

    def test_calculate_mean_genes_all_na(self):
        """Test calculating mean genes when all values are NA"""
        cluster_members = {
            'rep1': ['member1']
        }
        all_results = {
            'rep1': {'checkv_viral_genes': 'NA'},
            'member1': {'checkv_viral_genes': 'NA'}
        }

        result = extract_reps_stats.calculate_mean_genes(cluster_members, all_results)
        self.assertIsNone(result['rep1'])

    def test_calculate_mean_genes_missing_members(self):
        """Test calculating mean genes when some members are missing from results"""
        cluster_members = {
            'rep1': ['member1', 'member2']
        }
        all_results = {
            'rep1': {'checkv_viral_genes': '10'},
            'member1': {'checkv_viral_genes': '20'}
        }

        result = extract_reps_stats.calculate_mean_genes(cluster_members, all_results)
        self.assertEqual(result['rep1'], 15.0)


class TestParseTaxonomy(unittest.TestCase):
    """Test parse_taxonomy function"""

    def test_parse_complete_taxonomy(self):
        """Test parsing complete taxonomy string"""
        taxonomy_str = "Viruses;Heunggongvirae;Uroviricota;Caudoviricetes"
        result = extract_reps_stats.parse_taxonomy(taxonomy_str)
        self.assertEqual(result[0], 'Viruses')
        self.assertEqual(result[1], 'Heunggongvirae')
        self.assertEqual(result[2], 'Uroviricota')
        self.assertEqual(result[3], 'Caudoviricetes')

    def test_parse_taxonomy_with_missing_levels(self):
        """Test parsing taxonomy with missing levels"""
        taxonomy_str = "Viruses;;Uroviricota;;Caudovirales"
        levels = ['realm', 'kingdom', 'phylum', 'class', 'order']
        result = extract_reps_stats.parse_taxonomy(taxonomy_str, levels)
        self.assertEqual(result[0], 'Viruses')
        self.assertEqual(result[1], 'Unclassified_kingdom_viruses')
        self.assertEqual(result[2], 'Uroviricota')
        self.assertEqual(result[3], 'Unclassified_class_uroviricota')
        self.assertEqual(result[4], 'Caudovirales')

    def test_parse_empty_taxonomy(self):
        """Test parsing empty taxonomy"""
        result = extract_reps_stats.parse_taxonomy("")
        self.assertEqual(result, [])

    def test_parse_na_taxonomy(self):
        """Test parsing NA taxonomy"""
        result = extract_reps_stats.parse_taxonomy("NA")
        self.assertEqual(result, [])

    def test_parse_taxonomy_url_encoded(self):
        """Test parsing URL-encoded taxonomy"""
        taxonomy_str = "Viruses%3BHeunggongvirae"
        result = extract_reps_stats.parse_taxonomy(taxonomy_str)
        self.assertIn('Viruses', result[0])

    def test_parse_taxonomy_trailing_empty(self):
        """Test parsing taxonomy with trailing empty values"""
        taxonomy_str = "Viruses;Heunggongvirae;;;"
        result = extract_reps_stats.parse_taxonomy(taxonomy_str)
        self.assertEqual(len(result), 2)


class TestGenerateKronaFile(unittest.TestCase):
    """Test generate_krona_file function"""

    def test_generate_krona_file(self):
        """Test generating Krona format file"""
        results = {
            'seq1': {'taxonomy': 'Viruses;Heunggongvirae'},
            'seq2': {'taxonomy': 'Viruses;Heunggongvirae'},
            'seq3': {'taxonomy': 'Viruses;Duplodnaviria'}
        }
        viral_names = ['seq1', 'seq2', 'seq3']

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            krona_file = f.name

        try:
            extract_reps_stats.generate_krona_file(results, viral_names, krona_file)

            with open(krona_file, 'r') as f:
                lines = f.readlines()

            self.assertEqual(len(lines), 2)
            self.assertTrue(lines[0].startswith('2\t'))
            self.assertTrue(lines[1].startswith('1\t'))
        finally:
            os.unlink(krona_file)

    def test_generate_krona_file_with_na(self):
        """Test generating Krona file with NA taxonomy"""
        results = {
            'seq1': {'taxonomy': 'Viruses;Heunggongvirae'},
            'seq2': {'taxonomy': 'NA'}
        }
        viral_names = ['seq1', 'seq2']

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            krona_file = f.name

        try:
            extract_reps_stats.generate_krona_file(results, viral_names, krona_file)

            with open(krona_file, 'r') as f:
                lines = f.readlines()

            self.assertEqual(len(lines), 1)
        finally:
            os.unlink(krona_file)


if __name__ == '__main__':
    unittest.main()
