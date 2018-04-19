#######################################################################
# Unit tests for qc/illumina_qc.py
#######################################################################

import unittest

from auto_process_ngs.qc.illumina_qc import IlluminaQC
from auto_process_ngs.qc.illumina_qc import fastq_screen_output
from auto_process_ngs.qc.illumina_qc import fastqc_output

class TestIlluminaQC(unittest.TestCase):
    def test_illumina_qc_single_fastq(self):
        """IlluminaQC: generates commands for single Fastq
        """
        illumina_qc = IlluminaQC("/path/to/qc")
        self.assertEqual(illumina_qc.qc_dir,"/path/to/qc")
        cmds = illumina_qc.commands("/path/to/fastqs/test_S1_R1.fastq.gz")
        self.assertEqual(len(cmds),1)
        self.assertEqual(str(cmds[0]),
                         "illumina_qc.sh "
                         "/path/to/fastqs/test_S1_R1.fastq.gz "
                         "--threads 1 "
                         "--qc_dir /path/to/qc")
    def test_illumina_qc_fastq_pair(self):
        """IlluminaQC: generates commands for Fastq pair
        """
        illumina_qc = IlluminaQC("/path/to/qc")
        self.assertEqual(illumina_qc.qc_dir,"/path/to/qc")
        cmds = illumina_qc.commands("/path/to/fastqs/test_S1_R1.fastq.gz",
                                    "/path/to/fastqs/test_S1_R2.fastq.gz")
        self.assertEqual(len(cmds),2)
        self.assertEqual(str(cmds[0]),
                         "illumina_qc.sh "
                         "/path/to/fastqs/test_S1_R1.fastq.gz "
                         "--threads 1 "
                         "--qc_dir /path/to/qc")
        self.assertEqual(str(cmds[1]),
                         "illumina_qc.sh "
                         "/path/to/fastqs/test_S1_R2.fastq.gz "
                         "--threads 1 "
                         "--qc_dir /path/to/qc")

class TestFastqScreenOutputFunction(unittest.TestCase):
    def test_fastq_screen_output(self):
        """fastq_screen_output: handles .fastq file
        """
        self.assertEqual(fastq_screen_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq',
                                             'model_organisms'),
                         ('PB1_ATTAGG_L001_R1_001_model_organisms_screen.png',
                          'PB1_ATTAGG_L001_R1_001_model_organisms_screen.txt'))
    def test_fastq_screen_output_fastqgz(self):
        """fastq_screen_output: handles fastq.gz file
        """
        self.assertEqual(fastq_screen_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz',
                                             'model_organisms'),
                         ('PB1_ATTAGG_L001_R1_001_model_organisms_screen.png',
                          'PB1_ATTAGG_L001_R1_001_model_organisms_screen.txt'))

class TestFastqcOutputFunction(unittest.TestCase):
    def test_fastqc_output(self):
        """fastqc_output: handles .fastq file
        """
        self.assertEqual(fastqc_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq'),
                         ('PB1_ATTAGG_L001_R1_001_fastqc',
                          'PB1_ATTAGG_L001_R1_001_fastqc.html',
                          'PB1_ATTAGG_L001_R1_001_fastqc.zip'))
    def test_fastqc_output_fastqgz(self):
        """fastqc_output: handles fastq.gz file
        """
        self.assertEqual(fastqc_output('/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz'),
                         ('PB1_ATTAGG_L001_R1_001_fastqc',
                          'PB1_ATTAGG_L001_R1_001_fastqc.html',
                          'PB1_ATTAGG_L001_R1_001_fastqc.zip'))
