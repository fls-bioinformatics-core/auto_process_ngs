#######################################################################
# Unit tests for qc/reporting.py
#######################################################################

import unittest
import os
import tempfile
import shutil
import zipfile
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.mockqc import MockQCOutputs
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.analysis import AnalysisSample
from auto_process_ngs.qc.reporting import QCReporter
from auto_process_ngs.qc.reporting import FastqSet
from auto_process_ngs.qc.reporting import verify
from auto_process_ngs.qc.reporting import report
from auto_process_ngs.qc.reporting import pretty_print_reads

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestQCReporter(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None
    def tearDown(self):
        # Remove temporary working dir
        if not REMOVE_TEST_OUTPUTS:
            return
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def _make_working_dir(self):
        # Create a temporary working directory
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_QCReporter')
    def _make_analysis_project(self,paired_end=True,fastq_dir=None,
                               qc_dir="qc",fastq_names=None):
        # Create a mock Analysis Project directory
        self._make_working_dir()
        # Generate names for fastq files to add
        if paired_end:
            reads = (1,2)
        else:
            reads = (1,)
        sample_names = ('PJB1','PJB2')
        if fastq_names is None:
            fastq_names = []
            for i,sname in enumerate(sample_names,start=1):
                for read in reads:
                    fq = "%s_S%d_R%d_001.fastq.gz" % (sname,i,read)
                    fastq_names.append(fq)
        self.analysis_dir = MockAnalysisProject('PJB',fastq_names)
        # Create the mock directory
        self.analysis_dir.create(top_dir=self.wd)
        # Populate with fake QC products
        qc_dir = os.path.join(self.wd,self.analysis_dir.name,qc_dir)
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
    def test_qcreporter_single_end(self):
        """QCReporter: single-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.SE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.SE.html')))
    def test_qcreporter_paired_end(self):
        """QCReporter: paired-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=True)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_paired_end_with_non_default_fastq_dir(self):
        """QCReporter: paired-end data with non-default fastq dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_dir="fastqs.non_default")
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_paired_end_with_no_fastq_dir(self):
        """QCReporter: paired-end data with no fastq dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_dir=".")
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_paired_end_with_non_default_qc_dir(self):
        """QCReporter: paired-end data with non-default QC dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   qc_dir="qc.non_default")
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify(qc_dir="qc.non_default"))
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'),
                        qc_dir="qc.non_default")
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_paired_end_with_non_canonical_fastq_names(self):
        """QCReporter: paired-end data with non-canonical fastq names
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_names=
                                                   ("PJB1_S1_R1_001_paired.fastq.gz",
                                                    "PJB1_S1_R2_001_paired.fastq.gz",
                                                    "PJB2_S2_R1_001_paired.fastq.gz",
                                                    "PJB2_S2_R2_001_paired.fastq.gz",))
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.non_canonical.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.non_canonical.html')))
    def test_qcreporter_single_end_make_zip_file(self):
        """QCReporter: single-end data: make ZIP file
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,
                                              'PJB',
                                              'report.SE.html'),
                        make_zip=True)
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'PJB','report.SE.html')))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'PJB','report.SE.PJB.zip')))
        contents = zipfile.ZipFile(
            os.path.join(self.wd,'PJB',
                         'report.SE.PJB.zip')).namelist()
        print(contents)
        expected = (
            'report.SE.PJB/report.SE.html',
            'report.SE.PJB/qc/PJB1_S1_R1_001_fastqc.html',
            'report.SE.PJB/qc/PJB1_S1_R1_001_model_organisms_screen.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_model_organisms_screen.txt',
            'report.SE.PJB/qc/PJB1_S1_R1_001_other_organisms_screen.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_other_organisms_screen.txt',
            'report.SE.PJB/qc/PJB1_S1_R1_001_rRNA_screen.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_rRNA_screen.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_fastqc.html',
            'report.SE.PJB/qc/PJB2_S2_R1_001_model_organisms_screen.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_model_organisms_screen.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_other_organisms_screen.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_other_organisms_screen.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_rRNA_screen.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_rRNA_screen.txt')
        for f in expected:
            self.assertTrue(f in contents,"%s is missing from ZIP file" % f)

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

class TestVerifyFunction(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None
    def tearDown(self):
        # Remove temporary working dir
        if not REMOVE_TEST_OUTPUTS:
            return
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def _make_working_dir(self):
        # Create a temporary working directory
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_QCReporter')
    def _make_analysis_project(self,paired_end=True,fastq_dir=None,
                               qc_dir="qc",fastq_names=None):
        # Create a mock Analysis Project directory
        self._make_working_dir()
        # Generate names for fastq files to add
        if paired_end:
            reads = (1,2)
        else:
            reads = (1,)
        sample_names = ('PJB1','PJB2')
        if fastq_names is None:
            fastq_names = []
            for i,sname in enumerate(sample_names,start=1):
                for read in reads:
                    fq = "%s_S%d_R%d_001.fastq.gz" % (sname,i,read)
                    fastq_names.append(fq)
        self.analysis_dir = MockAnalysisProject('PJB',fastq_names)
        # Create the mock directory
        self.analysis_dir.create(top_dir=self.wd)
        # Populate with fake QC products
        project_dir = os.path.join(self.wd,self.analysis_dir.name)
        UpdateAnalysisProject(AnalysisProject(project_dir)).\
            add_qc_outputs(qc_dir=qc_dir,
                           include_report=False,
                           include_zip_file=False,
                           include_multiqc=False)
        return project_dir
    def test_verify_single_end(self):
        """verify: single-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        self.assertTrue(verify(project))
    def test_verify_paired_end(self):
        """verify: paired-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=True)
        project = AnalysisProject('PJB',analysis_dir)
        self.assertTrue(verify(project))
    def test_verify_paired_end_with_non_default_fastq_dir(self):
        """verify: paired-end data with non-default fastq dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_dir="fastqs.non_default")
        project = AnalysisProject('PJB',analysis_dir)
        self.assertTrue(verify(project))
    def test_verify_paired_end_with_no_fastq_dir(self):
        """verify: paired-end data with no fastq dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_dir=".")
        project = AnalysisProject('PJB',analysis_dir)
        self.assertTrue(verify(project))
    def test_verify_paired_end_with_non_default_qc_dir(self):
        """verify: paired-end data with non-default QC dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   qc_dir="qc.non_default")
        project = AnalysisProject('PJB',analysis_dir)
        self.assertTrue(verify(project,qc_dir="qc.non_default"))
    def test_verify_paired_end_with_non_canonical_fastq_names(self):
        """verify: paired-end data with non-canonical fastq names
        """
        analysis_dir = self._make_analysis_project(
            paired_end=True,
            fastq_names=
            ("PJB1_S1_R1_001_paired.fastq.gz",
             "PJB1_S1_R2_001_paired.fastq.gz",
             "PJB2_S2_R1_001_paired.fastq.gz",
             "PJB2_S2_R2_001_paired.fastq.gz",))
        project = AnalysisProject('PJB',analysis_dir)
        self.assertTrue(verify(project))

class TestReportFunction(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None
    def tearDown(self):
        # Remove temporary working dir
        if not REMOVE_TEST_OUTPUTS:
            return
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def _make_working_dir(self):
        # Create a temporary working directory
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_QCReporter')
    def _make_analysis_project(self,paired_end=True,fastq_dir=None,
                               qc_dir="qc",fastq_names=None):
        # Create a mock Analysis Project directory
        self._make_working_dir()
        # Generate names for fastq files to add
        if paired_end:
            reads = (1,2)
        else:
            reads = (1,)
        sample_names = ('PJB1','PJB2')
        if fastq_names is None:
            fastq_names = []
            for i,sname in enumerate(sample_names,start=1):
                for read in reads:
                    fq = "%s_S%d_R%d_001.fastq.gz" % (sname,i,read)
                    fastq_names.append(fq)
        self.analysis_dir = MockAnalysisProject('PJB',fastq_names)
        # Create the mock directory
        self.analysis_dir.create(top_dir=self.wd)
        # Populate with fake QC products
        qc_dir = os.path.join(self.wd,self.analysis_dir.name,qc_dir)
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
    def test_report_single_end(self):
        """report: single-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        report(project,filename=os.path.join(self.wd,'report.SE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.SE.html')))
    def test_report_paired_end(self):
        """report: paired-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=True)
        project = AnalysisProject('PJB',analysis_dir)
        report(project,filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_report_paired_end_with_non_default_fastq_dir(self):
        """report: paired-end data with non-default fastq dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_dir="fastqs.non_default")
        project = AnalysisProject('PJB',analysis_dir)
        report(project,filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_report_paired_end_with_no_fastq_dir(self):
        """report: paired-end data with no fastq dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_dir=".")
        project = AnalysisProject('PJB',analysis_dir)
        report(project,filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_report_paired_end_with_non_default_qc_dir(self):
        """report: paired-end data with non-default QC dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   qc_dir="qc.non_default")
        project = AnalysisProject('PJB',analysis_dir)
        report(project,filename=os.path.join(self.wd,'report.PE.html'),
               qc_dir="qc.non_default")
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_report_paired_end_with_non_canonical_fastq_names(self):
        """report: paired-end data with non-canonical fastq names
        """
        analysis_dir = self._make_analysis_project(
            paired_end=True,
            fastq_names=
            ("PJB1_S1_R1_001_paired.fastq.gz",
             "PJB1_S1_R2_001_paired.fastq.gz",
             "PJB2_S2_R1_001_paired.fastq.gz",
             "PJB2_S2_R2_001_paired.fastq.gz",))
        project = AnalysisProject('PJB',analysis_dir)
        report(project,
               filename=os.path.join(self.wd,'report.non_canonical.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.non_canonical.html')))
    def test_report_single_end_make_zip_file(self):
        """report: single-end data: make ZIP file
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        report(project,filename=os.path.join(self.wd,
                                             'PJB',
                                             'report.SE.html'),
               make_zip=True)
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'PJB','report.SE.html')))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'PJB','report.SE.PJB.zip')))
        contents = zipfile.ZipFile(
            os.path.join(self.wd,'PJB',
                         'report.SE.PJB.zip')).namelist()
        print(contents)
        expected = (
            'report.SE.PJB/report.SE.html',
            'report.SE.PJB/qc/PJB1_S1_R1_001_fastqc.html',
            'report.SE.PJB/qc/PJB1_S1_R1_001_model_organisms_screen.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_model_organisms_screen.txt',
            'report.SE.PJB/qc/PJB1_S1_R1_001_other_organisms_screen.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_other_organisms_screen.txt',
            'report.SE.PJB/qc/PJB1_S1_R1_001_rRNA_screen.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_rRNA_screen.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_fastqc.html',
            'report.SE.PJB/qc/PJB2_S2_R1_001_model_organisms_screen.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_model_organisms_screen.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_other_organisms_screen.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_other_organisms_screen.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_rRNA_screen.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_rRNA_screen.txt')
        for f in expected:
            self.assertTrue(f in contents,"%s is missing from ZIP file" % f)

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
