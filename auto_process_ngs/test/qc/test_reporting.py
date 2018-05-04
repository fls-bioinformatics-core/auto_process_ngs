#######################################################################
# Unit tests for qc/reporting.py
#######################################################################

import unittest
import os
import tempfile
import shutil

from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mockqc import MockQCOutputs
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.utils import AnalysisSample
from auto_process_ngs.qc.reporting import QCReporter
from auto_process_ngs.qc.reporting import FastqSet
from auto_process_ngs.qc.reporting import pretty_print_reads

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
        """QCReporter: paired-end data
        """
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
        """QCReporter: single-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertEqual(reporter.name,'PJB')
        self.assertFalse(reporter.paired_end)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.SE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.SE.html')))

class TestFastqSet(unittest.TestCase):
    def test_fastqset_PE(self):
        """FastqSet: handles paired-end data (Fastq pair)
        """
        fqset = FastqSet('/data/PB/PB1_ATTAGG_L001_R1_001.fastq',
                         '/data/PB/PB1_ATTAGG_L001_R2_001.fastq')
        # r1/r2 properties
        self.assertEqual(fqset.r1,'/data/PB/PB1_ATTAGG_L001_R1_001.fastq')
        self.assertEqual(fqset.r2,'/data/PB/PB1_ATTAGG_L001_R2_001.fastq')
        # __getitem__ method
        self.assertEqual(fqset[0],'/data/PB/PB1_ATTAGG_L001_R1_001.fastq')
        self.assertEqual(fqset[1],'/data/PB/PB1_ATTAGG_L001_R2_001.fastq')
        # fastqs property
        self.assertEqual(fqset.fastqs,
                         ['/data/PB/PB1_ATTAGG_L001_R1_001.fastq',
                         '/data/PB/PB1_ATTAGG_L001_R2_001.fastq'])
    def test_fastqset_SE(self):
        """FastqSet: handles single-end data (single Fastq)
        """
        fqset = FastqSet('/data/PB/PB1_ATTAGG_L001_R1_001.fastq')
        # r1/r2 properties
        self.assertEqual(fqset.r1,'/data/PB/PB1_ATTAGG_L001_R1_001.fastq')
        self.assertEqual(fqset.r2,None)
        # __getitem__ method
        self.assertEqual(fqset[0],'/data/PB/PB1_ATTAGG_L001_R1_001.fastq')
        try:
            fqset[1]
            self.fail("Attempt to access index 1 should raise IndexError")
        except IndexError:
            pass
        except Exception:
            self.fail("Attempt to access index 1 should raise IndexError")
        # fastqs property
        self.assertEqual(fqset.fastqs,
                         ['/data/PB/PB1_ATTAGG_L001_R1_001.fastq'])

class TestPrettyPrintReadsFunction(unittest.TestCase):
    def test_pretty_print_reads(self):
        """pretty_print_reads: handles different inputs
        """
        self.assertEqual(pretty_print_reads(1),"1")
        self.assertEqual(pretty_print_reads(12),"12")
        self.assertEqual(pretty_print_reads(117),"117")
        self.assertEqual(pretty_print_reads(1024),"1,024")
        self.assertEqual(pretty_print_reads(33385500),"33,385,500")
        self.assertEqual(pretty_print_reads(112839902),"112,839,902")
        self.assertEqual(pretty_print_reads(10212341927),"10,212,341,927")
