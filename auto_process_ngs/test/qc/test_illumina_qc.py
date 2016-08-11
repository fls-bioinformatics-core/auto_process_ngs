#######################################################################
# Unit tests for qc/illumina_qc.py
#######################################################################

import unittest

from auto_process_ngs.utils import AnalysisSample
from auto_process_ngs.qc.illumina_qc import get_fastq_pairs

class TestGetFastqPairsFunction(unittest.TestCase):
    def test_get_fastq_pairs_paired_end(self):
        s = AnalysisSample('PB1')
        s.add_fastq('/data/PB/PB1_ATTAGG_L001_R1_001.fastq')
        s.add_fastq('/data/PB/PB1_ATTAGG_L001_R2_001.fastq')
        s.add_fastq('/data/PB/PB1_GCCAAG_L002_R1_001.fastq')
        s.add_fastq('/data/PB/PB1_GCCAAG_L002_R2_001.fastq')
        fq_pair = get_fastq_pairs(s)
        self.assertEqual(len(fq_pair),2)
        self.assertEqual(fq_pair[0].r1,'/data/PB/PB1_ATTAGG_L001_R1_001.fastq')
        self.assertEqual(fq_pair[0].r2,'/data/PB/PB1_ATTAGG_L001_R2_001.fastq')
        self.assertEqual(fq_pair[1].r1,'/data/PB/PB1_GCCAAG_L002_R1_001.fastq')
        self.assertEqual(fq_pair[1].r2,'/data/PB/PB1_GCCAAG_L002_R2_001.fastq')
    def test_get_fastq_pairs_single_end(self):
        s = AnalysisSample('PB1')
        s.add_fastq('/data/PB/PB1_ATTAGG_L001_R1_001.fastq')
        s.add_fastq('/data/PB/PB1_GCCAAG_L002_R1_001.fastq')
        fq_pair = get_fastq_pairs(s)
        self.assertEqual(len(fq_pair),2)
        self.assertEqual(fq_pair[0].r1,'/data/PB/PB1_ATTAGG_L001_R1_001.fastq')
        self.assertEqual(fq_pair[0].r2,None)
        self.assertEqual(fq_pair[1].r1,'/data/PB/PB1_GCCAAG_L002_R1_001.fastq')
        self.assertEqual(fq_pair[1].r2,None)

from auto_process_ngs.qc.illumina_qc import fastq_screen_output
class TestFastqScreenOutputFunction(unittest.TestCase):
    def test_fastq_screen_output(self):
        self.assertEqual(fastq_screen_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq',
                                             'model_organisms'),
                         ('PB1_ATTAGG_L001_R1_001_model_organisms_screen.png',
                          'PB1_ATTAGG_L001_R1_001_model_organisms_screen.txt'))
    def test_fastq_screen_output_fastqgz(self):
        self.assertEqual(fastq_screen_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz',
                                             'model_organisms'),
                         ('PB1_ATTAGG_L001_R1_001_model_organisms_screen.png',
                          'PB1_ATTAGG_L001_R1_001_model_organisms_screen.txt'))

from auto_process_ngs.qc.illumina_qc import fastqc_output
class TestFastqcOutputFunction(unittest.TestCase):
    def test_fastqc_output(self):
        self.assertEqual(fastqc_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq'),
                         ('PB1_ATTAGG_L001_R1_001_fastqc',
                          'PB1_ATTAGG_L001_R1_001_fastqc.html',
                          'PB1_ATTAGG_L001_R1_001_fastqc.zip'))
    def test_fastqc_output_fastqgz(self):
        self.assertEqual(fastqc_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz'),
                         ('PB1_ATTAGG_L001_R1_001_fastqc',
                          'PB1_ATTAGG_L001_R1_001_fastqc.html',
                          'PB1_ATTAGG_L001_R1_001_fastqc.zip'))
    
from auto_process_ngs.qc.illumina_qc import pretty_print_reads
class TestPrettyPrintReadsFunction(unittest.TestCase):
    def test_pretty_print_reads(self):
        self.assertEqual(pretty_print_reads(1),"1")
        self.assertEqual(pretty_print_reads(12),"12")
        self.assertEqual(pretty_print_reads(117),"117")
        self.assertEqual(pretty_print_reads(1024),"1,024")
        self.assertEqual(pretty_print_reads(33385500),"33,385,500")
        self.assertEqual(pretty_print_reads(112839902),"112,839,902")
        self.assertEqual(pretty_print_reads(10212341927),"10,212,341,927")
