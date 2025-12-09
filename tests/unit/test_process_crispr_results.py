#!/usr/bin/env python3
"""
Unit tests for process_crispr_results.py
"""

import unittest
import tempfile
import os
import sys
from pathlib import Path
from io import StringIO

# Add bin directory to path to import the script
bin_dir = Path(__file__).parent.parent.parent / "bin"
sys.path.insert(0, str(bin_dir))

import bin.process_crispr_results as process_crispr_results


class TestGetCrisprId(unittest.TestCase):
    """Test get_crispr_id function"""

    def test_get_crispr_id_from_crispr_feature(self):
        """Test extracting CRISPR ID from CRISPR feature line"""
        line = "seq1\tCRISPRCasFinder\tCRISPR\t100\t200\t.\t+\t.\tName=seq1_100_200;evidence_level=3\n"
        result = process_crispr_results.get_crispr_id(line)
        self.assertEqual(result, "seq1_100_200")

    def test_get_crispr_id_from_repeat_feature(self):
        """Test extracting CRISPR ID from repeat feature line"""
        line = "seq1\tCRISPRCasFinder\trepeat_unit\t100\t120\t.\t+\t.\tParent=seq1_100_200;ID=repeat1\n"
        result = process_crispr_results.get_crispr_id(line)
        self.assertEqual(result, "seq1_100_200")

    def test_get_crispr_id_empty(self):
        """Test extracting CRISPR ID when not present"""
        line = "seq1\tCRISPRCasFinder\tCRISPR\t100\t200\t.\t+\t.\tID=other\n"
        result = process_crispr_results.get_crispr_id(line)
        self.assertEqual(result, "")


class TestCalcAtPercentage(unittest.TestCase):
    """Test calc_at_percentage function"""

    def test_calc_at_percentage_50_50(self):
        """Test calculating AT percentage for 50/50 AT/GC"""
        from Bio.Seq import Seq
        seq = Seq("ATGC")
        result = process_crispr_results.calc_at_percentage(seq)
        self.assertEqual(result, "50")

    def test_calc_at_percentage_all_at(self):
        """Test calculating AT percentage for all AT"""
        from Bio.Seq import Seq
        seq = Seq("ATATATATAT")
        result = process_crispr_results.calc_at_percentage(seq)
        self.assertEqual(result, "100")

    def test_calc_at_percentage_no_at(self):
        """Test calculating AT percentage for no AT"""
        from Bio.Seq import Seq
        seq = Seq("GCGCGCGC")
        result = process_crispr_results.calc_at_percentage(seq)
        self.assertEqual(result, "0")

    def test_calc_at_percentage_lowercase(self):
        """Test calculating AT percentage with lowercase"""
        from Bio.Seq import Seq
        seq = Seq("atgc")
        result = process_crispr_results.calc_at_percentage(seq)
        self.assertEqual(result, "50")


class TestCheckEndPosition(unittest.TestCase):
    """Test check_end_position function"""

    def test_check_end_position_within_bounds(self):
        """Test checking end position within contig bounds"""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        seq_records = {
            'contig1': SeqRecord(Seq("ATGCATGCATGC"), id='contig1')
        }
        result = process_crispr_results.check_end_position('contig1', 10, seq_records)
        self.assertEqual(result, 10)

    def test_check_end_position_exceeds_bounds(self):
        """Test checking end position that exceeds contig bounds"""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        seq_records = {
            'contig1': SeqRecord(Seq("ATGCATGCATGC"), id='contig1')
        }
        result = process_crispr_results.check_end_position('contig1', 20, seq_records)
        self.assertEqual(result, 12)


class TestProcessTsv(unittest.TestCase):
    """Test process_tsv function"""

    def test_process_tsv_basic(self):
        """Test processing basic TSV report"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tsv_in:
            tsv_in.write("Strain\tName\tSequence\tCRISPR_Id\tCRISPR_Start\tCRISPR_End\tEvidenceLevel\n")
            tsv_in.write("sample1\tname1\tseq1\tcrispr1\t100\t200\t3\n")
            tsv_in.write("sample2\tname2\tseq2\tcrispr2\t300\t400\t4\n")
            tsv_file = tsv_in.name

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tsv_out:
            output_file = tsv_out.name

        try:
            hits, hq_hits, evidence = process_crispr_results.process_tsv(tsv_file, output_file)

            self.assertIn('seq1', hits)
            self.assertIn('seq2', hits)
            self.assertIn('name1_100_200', hq_hits)
            self.assertIn('name2_300_400', hq_hits)
            self.assertEqual(evidence['name1_100_200'], '3')
            self.assertEqual(evidence['name2_300_400'], '4')
        finally:
            os.unlink(tsv_file)
            os.unlink(output_file)

    def test_process_tsv_low_evidence(self):
        """Test processing TSV with low evidence level"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tsv_in:
            tsv_in.write("Strain\tName\tSequence\tCRISPR_Id\tCRISPR_Start\tCRISPR_End\tEvidenceLevel\n")
            tsv_in.write("sample1\tname1\tseq1\tcrispr1\t100\t200\t1\n")
            tsv_file = tsv_in.name

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tsv_out:
            output_file = tsv_out.name

        try:
            hits, hq_hits, evidence = process_crispr_results.process_tsv(tsv_file, output_file)

            self.assertIn('seq1', hits)
            self.assertEqual(len(hq_hits), 0)
        finally:
            os.unlink(tsv_file)
            os.unlink(output_file)

    def test_process_tsv_empty_lines(self):
        """Test processing TSV with empty lines"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tsv_in:
            tsv_in.write("Strain\tName\tSequence\tCRISPR_Id\tCRISPR_Start\tCRISPR_End\tEvidenceLevel\n")
            tsv_in.write("\n")
            tsv_in.write("sample1\tname1\tseq1\tcrispr1\t100\t200\t3\n")
            tsv_file = tsv_in.name

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tsv_out:
            output_file = tsv_out.name

        try:
            hits, hq_hits, evidence = process_crispr_results.process_tsv(tsv_file, output_file)

            self.assertEqual(len(hits), 1)
        finally:
            os.unlink(tsv_file)
            os.unlink(output_file)

    def test_process_tsv_duplicate_sequences(self):
        """Test processing TSV with duplicate sequence names"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tsv_in:
            tsv_in.write("Strain\tName\tSequence\tCRISPR_Id\tCRISPR_Start\tCRISPR_End\tEvidenceLevel\n")
            tsv_in.write("sample1\tname1\tseq1\tcrispr1\t100\t200\t3\n")
            tsv_in.write("sample1\tname1\tseq1\tcrispr2\t300\t400\t4\n")
            tsv_file = tsv_in.name

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tsv_out:
            output_file = tsv_out.name

        try:
            hits, hq_hits, evidence = process_crispr_results.process_tsv(tsv_file, output_file)

            self.assertEqual(len(hits), 1)
            self.assertIn('seq1', hits)
        finally:
            os.unlink(tsv_file)
            os.unlink(output_file)


class TestFixAnnotation(unittest.TestCase):
    """Test fix_annotation function"""

    def test_fix_annotation_basic(self):
        """Test fixing annotation with sequence and AT%"""
        from Bio.Seq import Seq
        feature_seq = Seq("ATGC")
        at_percentage = "50"
        annotation = "ID=test;at%=0;sequence=UNKNOWN;Name=feature"

        result = process_crispr_results.fix_annotation(feature_seq, at_percentage, annotation)

        self.assertIn("at%=50", result)
        self.assertIn("sequence=ATGC", result)
        self.assertNotIn("UNKNOWN", result)

    def test_fix_annotation_no_unknown(self):
        """Test fixing annotation without UNKNOWN sequence"""
        from Bio.Seq import Seq
        feature_seq = Seq("ATGC")
        at_percentage = "50"
        annotation = "ID=test;at%=0;Name=feature"

        result = process_crispr_results.fix_annotation(feature_seq, at_percentage, annotation)

        self.assertIn("at%=50", result)


class TestAddEvidenceLevel(unittest.TestCase):
    """Test add_evidence_level function"""

    def test_add_evidence_level(self):
        """Test adding evidence level to GFF line"""
        line = "seq1\tCRISPRCasFinder\tCRISPR\t100\t200\t.\t+\t.\tName=seq1_100_200\n"
        evidence_levels = {'seq1_100_200': '3'}

        result = process_crispr_results.add_evidence_level(line, evidence_levels)

        self.assertIn("evidence_level=3", result)
        self.assertTrue(result.strip().endswith("evidence_level=3"))

    def test_add_evidence_level_missing_id(self):
        """Test adding evidence level when CRISPR ID is missing"""
        line = "seq1\tCRISPRCasFinder\tCRISPR\t100\t200\t.\t+\t.\tID=other\n"
        evidence_levels = {}

        # Should not raise exception, just log error
        result = process_crispr_results.add_evidence_level(line, evidence_levels)
        self.assertIsInstance(result, str)


class TestFixGffLine(unittest.TestCase):
    """Test fix_gff_line function"""

    def test_fix_gff_line_negative_start(self):
        """Test fixing GFF line with negative start coordinate"""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO

        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fa') as fasta:
            fasta.write(">seq1\n")
            fasta.write("ATGCATGCATGCATGCATGC\n")
            fasta_file = fasta.name

        try:
            line = "seq1\tCRISPRCasFinder\tCRISPR\t-5\t10\t.\t+\t.\tName=test\n"
            result = process_crispr_results.fix_gff_line(line, fasta_file)

            parts = result.strip().split('\t')
            self.assertEqual(parts[3], '1')
        finally:
            os.unlink(fasta_file)

    def test_fix_gff_line_both_negative(self):
        """Test fixing GFF line with both start and end negative"""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fa') as fasta:
            fasta.write(">seq1\n")
            fasta.write("ATGCATGCATGCATGCATGC\n")
            fasta_file = fasta.name

        try:
            line = "seq1\tCRISPRCasFinder\tCRISPR\t-5\t-1\t.\t+\t.\tName=test\n"
            result = process_crispr_results.fix_gff_line(line, fasta_file)

            self.assertIsNone(result)
        finally:
            os.unlink(fasta_file)

    def test_fix_gff_line_unknown_sequence(self):
        """Test fixing GFF line with UNKNOWN sequence"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fa') as fasta:
            fasta.write(">seq1\n")
            fasta.write("ATGCATGCATGCATGCATGC\n")
            fasta_file = fasta.name

        try:
            line = "seq1\tCRISPRCasFinder\tCRISPR\t1\t10\t.\t+\t.\tName=test;sequence=UNKNOWN;at%=0\n"
            result = process_crispr_results.fix_gff_line(line, fasta_file)

            self.assertIsNotNone(result)
            self.assertNotIn("UNKNOWN", result)
        finally:
            os.unlink(fasta_file)


if __name__ == '__main__':
    unittest.main()
