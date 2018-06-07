#######################################################################
# Unit tests for qc/fastq_strand.py
#######################################################################

import unittest
import tempfile

from auto_process_ngs.qc.fastq_strand import Fastqstrand
from auto_process_ngs.qc.fastq_strand import build_fastq_strand_conf

class TestFastqstrand(unittest.TestCase):
    def setUp(self):
        fastq_strand_text = ""
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            self.fastq_strand_txt = fp.name
            fp.write(fastq_strand_text)
    def tearDown(self):
        try:
            os.remove(self.fastq_strand_txt)
        except Exception:
            pass
    def test_fastq_screen_handle_empty_output_file(self):
        """Fastqstrand handles empty output file
        """
        self.assertRaises(Exception,
                          Fastqstrand,
                          self.fastq_strand_txt)

class TestFastqstrand_0_0_1(unittest.TestCase):
    def setUp(self):
        fastq_strand_text = """#fastq_strand version: 0.0.1	#Aligner: STAR	#Reads in subset: 3
#Genome	1st forward	2nd reverse
hg38	13.13	93.21
mm10	15.02	91.18
"""
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            self.fastq_strand_txt = fp.name
            fp.write(fastq_strand_text)
    def tearDown(self):
        try:
            os.remove(self.fastq_strand_txt)
        except Exception:
            pass
    def test_handle_fastq_strand_0_0_1(self):
        """Fastqstrand handles output from version 0.0.1
        """
        fastq_strand = Fastqstrand(self.fastq_strand_txt)
        self.assertEqual(fastq_strand.version,'0.0.1')
        self.assertEqual(fastq_strand.txt,self.fastq_strand_txt)
        self.assertEqual(fastq_strand.genomes,['hg38','mm10'])
        self.assertEqual(fastq_strand.stats['hg38'].forward,13.13)
        self.assertEqual(fastq_strand.stats['hg38'].reverse,93.21)
        self.assertEqual(fastq_strand.stats['mm10'].forward,15.02)
        self.assertEqual(fastq_strand.stats['mm10'].reverse,91.18)

class TestBuildFastqStrandConfFunction(unittest.TestCase):
    indexes = dict(human="/data/to/hg38/STAR",
                   mouse="/data/to/mm10/STAR",)
    def test_build_fastq_strand_conf_no_organism(self):
        """build_fastq_strand_conf: no organisms
        """
        self.assertEqual(build_fastq_strand_conf(None,self.indexes),
                         None)
    def test_build_fastq_strand_conf_one_organism(self):
        """build_fastq_strand_conf: one organism
        """
        self.assertEqual(build_fastq_strand_conf(("human",),self.indexes),
                         "human\t/data/to/hg38/STAR")
    def test_build_fastq_strand_conf_multiple_organisms(self):
        """build_fastq_strand_conf: multiple organisms
        """
        self.assertEqual(build_fastq_strand_conf(("human","mouse"),self.indexes),
                         "human\t/data/to/hg38/STAR\nmouse\t/data/to/mm10/STAR")
    def test_build_fastq_strand_conf_unrecognised_organism(self):
        """build_fastq_strand_conf: unrecognised organism
        """
        self.assertEqual(build_fastq_strand_conf(("zebrafish",),self.indexes),
                         None)
