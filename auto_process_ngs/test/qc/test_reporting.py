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
                               qc_dir="qc",fastq_names=None,
                               include_seqlens=True,
                               include_cellranger_count=False,
                               include_cellranger_multi=False,
                               cellranger_pipelines=('cellranger',),
                               legacy_screens=False,
                               legacy_cellranger_outs=False):
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
            MockQCOutputs.fastq_screen_v0_9_2(fq,qc_dir,'model_organisms',
                                              legacy=legacy_screens)
            MockQCOutputs.fastq_screen_v0_9_2(fq,qc_dir,'other_organisms',
                                              legacy=legacy_screens)
            MockQCOutputs.fastq_screen_v0_9_2(fq,qc_dir,'rRNA',
                                              legacy=legacy_screens)
            # Sequence lengths
            if include_seqlens:
                MockQCOutputs.seqlens(fq,qc_dir)
        # Cellranger count
        if include_cellranger_count:
            for cellranger in cellranger_pipelines:
                if cellranger == "cellranger":
                    version = "6.1.2"
                    refdata = "/data/refdata-cellranger-2020-A"
                elif cellranger == "cellranger-atac":
                    version = "2.0.0"
                    refdata = "/data/refdata-cellranger-atac-2020-A"
                elif cellranger == "cellranger-arc":
                    version = "2.0.0"
                    refdata = "/data/refdata-cellranger-arc-2020-A"
                project_dir = os.path.join(self.wd,
                                           self.analysis_dir.name)
                if not legacy_cellranger_outs:
                    count_dir = os.path.join("cellranger_count",
                                             version,
                                             os.path.basename(refdata))
                else:
                    count_dir = "cellranger_count"
                UpdateAnalysisProject(AnalysisProject(project_dir)).\
                    add_cellranger_count_outputs(
                        reference_data_path=refdata,
                        qc_dir=qc_dir,
                        cellranger=cellranger,
                        prefix=count_dir)
        # Cellranger multi
        if include_cellranger_multi:
            project_dir = os.path.join(self.wd,self.analysis_dir.name)
            # Add the cellranger multi config.csv file
            multi_config = os.path.join(project_dir,"10x_multi_config.csv")
            with open(multi_config,'wt') as fp:
                fastq_dir = os.path.join(self.wd,
                                         "PJB",
                                         "fastqs")
                fp.write("""[gene-expression]
reference,/data/refdata-cellranger-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,%s,any,PJB1,gene expression,
PJB2_MC,%s,any,PJB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PJB_CML1,CMO301,CML1
PBB_CML2,CMO302,CML2
""" % (fastq_dir,fastq_dir))
            UpdateAnalysisProject(AnalysisProject(project_dir)).\
                add_cellranger_multi_outputs(
                    config_csv=multi_config,
                    qc_dir=qc_dir,
                    prefix=os.path.join("cellranger_multi",
                                        "6.1.2",
                                        "refdata-cellranger-2020-A"))
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
    def test_qcreporter_single_end_no_seq_lens(self):
        """QCReporter: single-end data (no sequence lengths)
        """
        analysis_dir = self._make_analysis_project(paired_end=False,
                                                   include_seqlens=False)
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
    def test_qcreporter_paired_end_cellranger_count(self):
        """QCReporter: paired-end data with cellranger 'count'
        """
        analysis_dir = self._make_analysis_project(
            paired_end=True,
            include_cellranger_count=True)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_paired_end_cellranger_multi(self):
        """QCReporter: paired-end data with cellranger 'multi'
        """
        analysis_dir = self._make_analysis_project(
            paired_end=True,
            include_cellranger_multi=True)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_paired_end_cellranger_count_and_multi(self):
        """QCReporter: paired-end data with cellranger 'count' and 'multi'
        """
        analysis_dir = self._make_analysis_project(
            paired_end=True,
            include_cellranger_count=True,
            include_cellranger_multi=True)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_paired_end_cellranger_count_multiome(self):
        """QCReporter: paired-end data with cellranger 'count' (multiome)
        """
        analysis_dir = self._make_analysis_project(
            paired_end=True,
            cellranger_pipelines=('cellranger-arc',),
            include_cellranger_count=True)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_paired_end_cellranger_count_multiome_and_scrnaseq(self):
        """QCReporter: paired-end data with cellranger 'count' (multiome+scRNAseq)
        """
        analysis_dir = self._make_analysis_project(
            paired_end=True,
            cellranger_pipelines=('cellranger',
                                  'cellranger-arc',),
            include_cellranger_count=True)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_paired_end_legacy_cellranger_count(self):
        """QCReporter: paired-end data with cellranger 'count' (legacy)
        """
        analysis_dir = self._make_analysis_project(
            paired_end=True,
            include_cellranger_count=True,
            legacy_cellranger_outs=True)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
    def test_qcreporter_paired_end_no_seq_lens(self):
        """QCReporter: paired-end data (no sequence lengths)
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   include_seqlens=False)
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
    def test_qcreporter_paired_end_with_legacy_screens(self):
        """QCReporter: paired-end data with legacy screen names
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   legacy_screens=True)
        project = AnalysisProject('PJB',analysis_dir)
        reporter = QCReporter(project)
        self.assertTrue(reporter.verify())
        reporter.report(filename=os.path.join(self.wd,'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,'report.PE.html')))
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
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_model_organisms.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_model_organisms.txt',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_other_organisms.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_other_organisms.txt',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_rRNA.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_rRNA.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_fastqc.html',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_model_organisms.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_model_organisms.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_other_organisms.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_other_organisms.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_rRNA.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_rRNA.txt')
        for f in expected:
            self.assertTrue(f in contents,"%s is missing from ZIP file" % f)
    def test_qcreporter_single_end_make_zip_file_legacy_screens(self):
        """QCReporter: single-end data: make ZIP file (legacy screens)
        """
        analysis_dir = self._make_analysis_project(paired_end=False,
                                                   legacy_screens=True)
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
                               qc_dir="qc",fastq_names=None,
                               legacy_screens=False):
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
                           include_multiqc=False,
                           legacy_screens=legacy_screens)
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
        self.assertTrue(verify(project))
    def test_verify_paired_end_with_no_fastq_dir(self):
        """verify: paired-end data with legacy screen names
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   legacy_screens=True)
        project = AnalysisProject('PJB',analysis_dir)
        self.assertTrue(verify(project))

class TestReportFunction(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None
        self.top_dir = None
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
            self.top_dir = os.path.join(self.wd,"Test")
            os.mkdir(self.top_dir)
    def _make_analysis_project(self,name="PJB",paired_end=True,fastq_dir=None,
                               qc_dir="qc",sample_names=None,fastq_names=None,
                               legacy_screens=False):
        # Create a mock Analysis Project directory
        self._make_working_dir()
        # Generate names for fastq files to add
        if paired_end:
            reads = (1,2)
        else:
            reads = (1,)
        if sample_names is None:
            sample_names = ('PJB1','PJB2')
        if fastq_names is None:
            fastq_names = []
            for i,sname in enumerate(sample_names,start=1):
                for read in reads:
                    fq = "%s_S%d_R%d_001.fastq.gz" % (sname,i,read)
                    fastq_names.append(fq)
        analysis_project = MockAnalysisProject(name,fastq_names)
        # Create the mock directory
        analysis_project.create(top_dir=self.top_dir)
        # Populate with fake QC products
        project_dir = os.path.join(self.top_dir,analysis_project.name)
        UpdateAnalysisProject(AnalysisProject(project_dir)).\
            add_qc_outputs(qc_dir=qc_dir,
                           include_report=False,
                           include_zip_file=False,
                           include_multiqc=True,
                           legacy_screens=legacy_screens)
        return project_dir
    def test_report_single_end(self):
        """report: single-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.SE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.SE.html')))
    def test_report_paired_end(self):
        """report: paired-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=True)
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))
    def test_report_paired_end_with_non_default_fastq_dir(self):
        """report: paired-end data with non-default fastq dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_dir="fastqs.non_default")
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))
    def test_report_paired_end_with_no_fastq_dir(self):
        """report: paired-end data with no fastq dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_dir=".")
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))
    def test_report_paired_end_with_non_default_qc_dir(self):
        """report: paired-end data with non-default QC dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   qc_dir="qc.non_default")
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'),
               qc_dir="qc.non_default")
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))
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
        report((project,),
               filename=os.path.join(self.top_dir,
                                     'report.non_canonical.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.non_canonical.html')))
    def test_report_paired_end_with_legacy_screens(self):
        """report: paired-end data with legacy screen names
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   legacy_screens=True)
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))
    def test_report_single_end_multiple_projects(self):
        """report: single-end data: two projects in one report
        """
        analysis_dir = self._make_analysis_project(name="PJB",
                                                   paired_end=False)
        analysis_dir2 = self._make_analysis_project(name="PJB2",
                                                    paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        project2 = AnalysisProject('PJB2',analysis_dir2)
        report((project,project2,),
               title="QC report: PJB & PJB2",
               filename=os.path.join(self.top_dir,
                                     'report.multiple_projects.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.multiple_projects.html')))
    def test_report_single_end_with_data_dir(self):
        """report: single-end data: use data directory
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),
               filename=os.path.join(self.top_dir,
                                     'PJB',
                                     'report.SE.html'),
               use_data_dir=True)
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB','report.SE.html')))
        self.assertTrue(os.path.isdir(
            os.path.join(self.top_dir,
                         'PJB',
                         'report.SE_data',
                         'Test_PJB',
                         'qc')))
        self.assertTrue(os.path.isdir(
            os.path.join(self.top_dir,'PJB','report.SE_data')))
        contents = os.listdir(os.path.join(self.top_dir,
                                           'PJB',
                                           'report.SE_data',
                                           'Test_PJB',
                                           'qc'))
        print(contents)
        expected = (
            'PJB1_S1_R1_001_fastqc.html',
            'PJB1_S1_R1_001_screen_model_organisms.png',
            'PJB1_S1_R1_001_screen_model_organisms.txt',
            'PJB1_S1_R1_001_screen_other_organisms.png',
            'PJB1_S1_R1_001_screen_other_organisms.txt',
            'PJB1_S1_R1_001_screen_rRNA.png',
            'PJB1_S1_R1_001_screen_rRNA.txt',
            'PJB2_S2_R1_001_fastqc.html',
            'PJB2_S2_R1_001_screen_model_organisms.png',
            'PJB2_S2_R1_001_screen_model_organisms.txt',
            'PJB2_S2_R1_001_screen_other_organisms.png',
            'PJB2_S2_R1_001_screen_other_organisms.txt',
            'PJB2_S2_R1_001_screen_rRNA.png',
            'PJB2_S2_R1_001_screen_rRNA.txt')
        for f in expected:
            self.assertTrue(f in contents,"%s is missing from data dir" % f)
        self.assertTrue(os.path.exists(os.path.join(self.top_dir,
                                                    'PJB',
                                                    'report.SE_data',
                                                    'Test_PJB',
                                                    'multiqc_report.html')),
                        "Missing multiqc_report.html")
    def test_report_single_end_make_zip_file(self):
        """report: single-end data: make ZIP file
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'PJB',
                                                'report.SE.html'),
               make_zip=True)
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB','report.SE.html')))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB','report.SE.PJB.zip')))
        contents = zipfile.ZipFile(
            os.path.join(self.top_dir,'PJB',
                         'report.SE.PJB.zip')).namelist()
        print(contents)
        expected = (
            'report.SE.PJB/report.SE.html',
            'report.SE.PJB/multiqc_report.html',
            'report.SE.PJB/qc/PJB1_S1_R1_001_fastqc.html',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_model_organisms.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_model_organisms.txt',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_other_organisms.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_other_organisms.txt',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_rRNA.png',
            'report.SE.PJB/qc/PJB1_S1_R1_001_screen_rRNA.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_fastqc.html',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_model_organisms.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_model_organisms.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_other_organisms.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_other_organisms.txt',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_rRNA.png',
            'report.SE.PJB/qc/PJB2_S2_R1_001_screen_rRNA.txt')
        for f in expected:
            self.assertTrue(f in contents,"%s is missing from ZIP file" % f)
    def test_report_single_end_make_zip_file_legacy_screens(self):
        """report: single-end data: make ZIP file with legacy screens
        """
        analysis_dir = self._make_analysis_project(paired_end=False,
                                                   legacy_screens=True)
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'PJB',
                                                'report.SE.html'),
               make_zip=True)
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB','report.SE.html')))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB','report.SE.PJB.zip')))
        contents = zipfile.ZipFile(
            os.path.join(self.top_dir,'PJB',
                         'report.SE.PJB.zip')).namelist()
        print(contents)
        expected = (
            'report.SE.PJB/report.SE.html',
            'report.SE.PJB/multiqc_report.html',
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
    def test_report_single_end_make_zip_file_with_data_dir(self):
        """report: single-end data: make ZIP file with data directory
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'PJB',
                                                'report.SE.html'),
               use_data_dir=True,
               make_zip=True)
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB','report.SE.html')))
        self.assertTrue(os.path.isdir(
            os.path.join(self.top_dir,'PJB','report.SE_data')))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB','report.SE.PJB.zip')))
        contents = zipfile.ZipFile(
            os.path.join(self.top_dir,'PJB',
                         'report.SE.PJB.zip')).namelist()
        print(contents)
        expected = (
            'report.SE.PJB/report.SE.html',
            'report.SE.PJB/report.SE_data/Test_PJB/multiqc_report.html',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB1_S1_R1_001_fastqc.html',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB1_S1_R1_001_screen_model_organisms.png',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB1_S1_R1_001_screen_model_organisms.txt',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB1_S1_R1_001_screen_other_organisms.png',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB1_S1_R1_001_screen_other_organisms.txt',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB1_S1_R1_001_screen_rRNA.png',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB1_S1_R1_001_screen_rRNA.txt',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB2_S2_R1_001_fastqc.html',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB2_S2_R1_001_screen_model_organisms.png',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB2_S2_R1_001_screen_model_organisms.txt',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB2_S2_R1_001_screen_other_organisms.png',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB2_S2_R1_001_screen_other_organisms.txt',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB2_S2_R1_001_screen_rRNA.png',
            'report.SE.PJB/report.SE_data/Test_PJB/qc/PJB2_S2_R1_001_screen_rRNA.txt')
        for f in expected:
            self.assertTrue(f in contents,"%s is missing from ZIP file" % f)
    def test_report_single_end_multiple_projects_with_zip_file_no_data_dir(self):
        """report: single-end data: fails with two projects in one report (ZIP file/no data directory)
        """
        analysis_dir = self._make_analysis_project(name="PJB",
                                                   sample_names=('PJB1',),
                                                   paired_end=False)
        analysis_dir2 = self._make_analysis_project(name="PJB2",
                                                    sample_names=('PJB2',),
                                                    paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        project2 = AnalysisProject('PJB2',analysis_dir2)
        self.assertRaises(Exception,
                          report,
                          (project,project2,),
                          title="QC report: PJB & PJB2",
                          filename=os.path.join(
                              self.top_dir,
                              'PJB',
                              'report.multiple_projects.html'),
                          make_zip=True)
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB',
                         'report.multiple_projects.html')))
        self.assertFalse(os.path.exists(
            os.path.join(self.top_dir,'PJB',
                         'report.multiple_projects.PJB.zip')))
    def test_report_single_end_multiple_projects_with_zip_file_duplicated_names_no_data_dir(self):
        """report: single-end data: fails with two projects in one report (duplicated names/ZIP file/no data directory)
        """
        analysis_dir = self._make_analysis_project(name="PJB",
                                                   paired_end=False)
        analysis_dir2 = self._make_analysis_project(name="PJB2",
                                                    paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        project2 = AnalysisProject('PJB2',analysis_dir2)
        self.assertRaises(Exception,
                          report,
                          (project,project2,),
                          title="QC report: PJB & PJB2",
                          filename=os.path.join(
                              self.top_dir,
                              'PJB',
                              'report.multiple_projects.html'),
                          make_zip=True)
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB',
                         'report.multiple_projects.html')))
        self.assertFalse(os.path.exists(
            os.path.join(self.top_dir,'PJB',
                         'report.multiple_projects.PJB.zip')))
    def test_report_single_end_multiple_projects_with_zip_file_duplicated_names_with_data_dir(self):
        """report: single-end data: two projects with duplicated names in one report, with ZIP file, with data directory
        """
        analysis_dir = self._make_analysis_project(name="PJB",
                                                   paired_end=False)
        analysis_dir2 = self._make_analysis_project(name="PJB2",
                                                    paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        project2 = AnalysisProject('PJB2',analysis_dir2)
        report((project,project2,),
               title="QC report: PJB & PJB2",
               filename=os.path.join(self.top_dir,'PJB',
                                     'report.multiple_projects.html'),
               use_data_dir=True,
               make_zip=True)
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB',
                         'report.multiple_projects.html')))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'PJB',
                         'report.multiple_projects.PJB.zip')))
        contents = zipfile.ZipFile(
            os.path.join(self.top_dir,'PJB',
                         'report.multiple_projects.PJB.zip')).namelist()
        print(contents)
        expected = (
            'report.multiple_projects.PJB/report.multiple_projects.html',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/multiqc_report.html',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB1_S1_R1_001_fastqc.html',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB1_S1_R1_001_screen_model_organisms.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB1_S1_R1_001_screen_model_organisms.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB1_S1_R1_001_screen_other_organisms.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB1_S1_R1_001_screen_other_organisms.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB1_S1_R1_001_screen_rRNA.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB1_S1_R1_001_screen_rRNA.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB2_S2_R1_001_fastqc.html',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB2_S2_R1_001_screen_model_organisms.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB2_S2_R1_001_screen_model_organisms.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB2_S2_R1_001_screen_other_organisms.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB2_S2_R1_001_screen_other_organisms.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB2_S2_R1_001_screen_rRNA.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB/qc/PJB2_S2_R1_001_screen_rRNA.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/multiqc_report.html',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB1_S1_R1_001_fastqc.html',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB1_S1_R1_001_screen_model_organisms.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB1_S1_R1_001_screen_model_organisms.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB1_S1_R1_001_screen_other_organisms.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB1_S1_R1_001_screen_other_organisms.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB1_S1_R1_001_screen_rRNA.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB1_S1_R1_001_screen_rRNA.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB2_S2_R1_001_fastqc.html',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB2_S2_R1_001_screen_model_organisms.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB2_S2_R1_001_screen_model_organisms.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB2_S2_R1_001_screen_other_organisms.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB2_S2_R1_001_screen_other_organisms.txt',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB2_S2_R1_001_screen_rRNA.png',
            'report.multiple_projects.PJB/report.multiple_projects_data/Test_PJB2/qc/PJB2_S2_R1_001_screen_rRNA.txt')
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
