#######################################################################
# Unit tests for qc/illumina_qc.py
#######################################################################

import unittest
import os
import tempfile
import shutil

from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mockqc import MockQCOutputs
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.utils import AnalysisSample
from auto_process_ngs.qc.illumina_qc import QCReporter
from auto_process_ngs.qc.illumina_qc import get_fastq_pairs

class TestQCReporter(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None
    def tearDown(self):
        # Remove temporary working dir
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def _make_working_dir(self):
        # Create a temporary working directory
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_QCReporter')
    def _make_analysis_project(self,paired_end=True):
        # Create a mock Analysis Project directory
        self._make_working_dir()
        # Generate names for fastq files to add
        if paired_end:
            reads = (1,2)
        else:
            reads = (1,)
        sample_names = ('PJB1','PJB2')
        fastq_names = []
        for i,sname in enumerate(sample_names,start=1):
            for read in reads:
                fq = "%s_S%d_R%d_001.fastq.gz" % (sname,i,read)
                fastq_names.append(fq)
        self.analysis_dir = MockAnalysisProject('PJB',fastq_names)
        # Create the mock directory
        self.analysis_dir.create(top_dir=self.wd)
        # Populate with fake QC products
        qc_dir = os.path.join(self.wd,self.analysis_dir.name,'qc')
        qc_logs = os.path.join(qc_dir,'logs')
        os.mkdir(qc_dir)
        os.mkdir(qc_logs)
        for fq in fastq_names:
            # FastQC
            MockQCOutputs.fastqc_v0_11_2(fq,qc_dir)
            # Fastq_screen
            MockQCOutputs.fastq_screen_v0_9_2(fq,qc_dir,'model_organisms')
            MockQCOutputs.fastq_screen_v0_9_2(fq,qc_dir,'other_organisms')
            MockQCOutputs.fastq_screen_v0_9_2(fq,qc_dir,'rRNA')
        return os.path.join(self.wd,self.analysis_dir.name)
    def test_qcreporter_paired_end(self):
        analysis_dir = self._make_analysis_project(paired_end=True)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertEqual(reporter.name,'PJB')
        self.assertTrue(reporter.paired_end)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_single_end(self):
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertEqual(reporter.name,'PJB')
        self.assertFalse(reporter.paired_end)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.SE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.SE.html')))

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
