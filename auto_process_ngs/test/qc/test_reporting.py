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
from auto_process_ngs.qc.reporting import QCOutputs
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

class TestQCOutputs(unittest.TestCase):
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
            self.wd = tempfile.mkdtemp(suffix='.test_QCOutputs')
    def _make_qc_dir(self,qc_dir,fastq_names,
                     screens=('model_organisms','other_organisms','rRNA',),
                     cellranger_pipelines=('cellranger',),
                     cellranger_samples=None,
                     cellranger_multi_samples=None,
                     include_fastqc=True,
                     include_fastq_screen=True,
                     include_strandedness=True,
                     include_seqlens=True,
                     include_multiqc=True,
                     include_cellranger_count=False,
                     include_cellranger_multi=False,
                     legacy_screens=False,
                     legacy_cellranger_outs=False):
        # Create working directory and qc dir
        self._make_working_dir()
        qc_dir = os.path.join(self.wd,qc_dir)
        os.mkdir(qc_dir)
        # Populate with fake QC products
        for fq in fastq_names:
            # FastQC
            if include_fastqc:
                MockQCOutputs.fastqc_v0_11_2(fq,qc_dir)
            # Fastq_screen
            if include_fastq_screen:
                for screen in screens:
                    MockQCOutputs.fastq_screen_v0_9_2(
                        fq,qc_dir,screen)
            # Strandedness
            if include_strandedness:
                MockQCOutputs.fastq_strand_v0_0_4(fq,qc_dir)
            # Sequence lengths
            if include_seqlens:
                MockQCOutputs.seqlens(fq,qc_dir)
        # MultiQC
        if include_multiqc:
            out_file = "multi%s_report.html" % os.path.basename(qc_dir)
            MockQCOutputs.multiqc(self.wd,
                                  multiqc_html=out_file,
                                  version="1.8")
        # Cellranger count
        if include_cellranger_count:
            for cellranger in cellranger_pipelines:
                # Set defaults
                if cellranger == "cellranger":
                    version = "6.1.2"
                    refdata = "/data/refdata-cellranger-2020-A"
                elif cellranger == "cellranger-atac":
                    version = "2.0.0"
                    refdata = "/data/refdata-cellranger-atac-2020-A"
                elif cellranger == "cellranger-arc":
                    version = "2.0.0"
                    refdata = "/data/refdata-cellranger-arc-2020-A"
                # Set top-level output dir
                if not legacy_cellranger_outs:
                    count_dir = os.path.join("cellranger_count",
                                             version,
                                             os.path.basename(refdata))
                else:
                    count_dir = "cellranger_count"
                # Make pipeline outputs
                for sample in cellranger_samples:
                    MockQCOutputs.cellranger_count(
                        sample,
                        qc_dir,
                        cellranger=cellranger,
                        version=version,
                        reference_data_path=refdata,
                        prefix=count_dir)
        # Cellranger multi
        if include_cellranger_multi:
            # Make cellranger multi config.csv file
            multi_config = os.path.join(self.wd,"10x_multi_config.csv")
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
PJB_CML2,CMO302,CML2
""" % (fastq_dir,fastq_dir))
            # Set top-level output dir
            multi_dir = os.path.join("cellranger_multi",
                                     "6.1.2",
                                     "refdata-cellranger-2020-A")
            # Make outputs
            MockQCOutputs.cellranger_multi(cellranger_multi_samples,
                                           qc_dir,
                                           config_csv=multi_config,
                                           prefix=multi_dir)
        return qc_dir

    def test_qcoutputs_single_end(self):
        """
	QCOutputs: single-end data
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'multiqc',
                          'screens_r1',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)

    def test_qcoutputs_single_end_no_fastqc(self):
        """
	QCOutputs: single-end data (no FastQC)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   include_fastqc=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['multiqc',
                          'screens_r1',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)

    def test_qcoutputs_single_end_no_screens(self):
        """
	QCOutputs: single-end data (no screens)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   include_fastq_screen=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'multiqc',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,[])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)

    def test_qcoutputs_single_end_no_seqlens(self):
        """
	QCOutputs: single-end data (no sequence lengths)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   include_seqlens=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'multiqc',
                          'screens_r1',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,None)
        self.assertEqual(qc_outputs.stats.min_sequence_length,None)
        self.assertEqual(qc_outputs.stats.max_sequence_length,None)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         [])

    def test_qcoutputs_single_end_no_standedness(self):
        """
	QCOutputs: single-end data (no strandedness)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   include_strandedness=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'multiqc',
                          'screens_r1',
                          'sequence_lengths'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)

    def test_qcoutputs_single_end_no_multiqc(self):
        """
	QCOutputs: single-end data (no MultiQC)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   include_multiqc=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'screens_r1',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)

    def test_qcoutputs_single_end_legacy_screen_naming(self):
        """
	QCOutputs: single-end data (legacy screen naming)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB2_S2_R1_001',
                                   ),
                                   legacy_screens=True)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'multiqc',
                          'screens_r1',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB2_S2_R1_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)

    def test_qcoutputs_paired_end(self):
        """
	QCOutputs: paired-end data
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_paired_end_no_fastqc(self):
        """
	QCOutputs: paired-end data (no FastQC)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_fastqc=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_paired_end_no_screens(self):
        """
	QCOutputs: paired-end data (no screens)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_fastq_screen=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,[])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_paired_end_no_seqlens(self):
        """
	QCOutputs: paired-end data (no sequence lengths)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_seqlens=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,None)
        self.assertEqual(qc_outputs.stats.min_sequence_length,None)
        self.assertEqual(qc_outputs.stats.max_sequence_length,None)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),[])

    def test_qcoutputs_paired_end_no_strandedness(self):
        """
	QCOutputs: paired-end data (no strandedness)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_strandedness=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_paired_end_no_multiqc(self):
        """
	QCOutputs: paired-end data (no MultiQC)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_multiqc=False)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_paired_end_legacy_screen_naming(self):
        """
	QCOutputs: paired-end data (legacy screen naming)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   legacy_screens=True)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_10x_cellranger_count(self):
        """
        QCOutputs: 10xGenomics data with cellranger 'count'
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '6.1.2' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_10x_cellranger_count_legacy(self):
        """
        QCOutputs: 10xGenomics data with cellranger 'count' (legacy format)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ),
                                   legacy_cellranger_outs=True)
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '?' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_10x_cellranger_atac_count(self):
        """
        QCOutputs: 10xGenomics data with cellranger-atac 'count'
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R3_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R3_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-atac',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger-atac_count',
                          'fastqc_r1',
                          'fastqc_r3',
                          'multiqc',
                          'screens_r1',
                          'screens_r3',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R3_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R3_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-atac-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r3'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger-atac': [ '2.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r3'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r3'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r3'],76)

    def test_qcoutputs_10x_multiome_gex(self):
        """
        QCOutputs: 10xGenomics multiome data (GEX component)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger',
                                                         'cellranger-arc',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger-arc_count',
                          'cellranger_count',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A',
                          '/data/refdata-cellranger-arc-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '6.1.2' ],
                           'cellranger-arc': [ '2.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_10x_multiome_atac(self):
        """
        QCOutputs: 10xGenomics multiome data (ATAC component)
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R3_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R3_001',
                                   ),
                                   include_cellranger_count=True,
                                   cellranger_pipelines=('cellranger-arc',
                                                         'cellranger-atac',),
                                   cellranger_samples=(
                                       'PJB1',
                                       'PJB2',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger-arc_count',
                          'cellranger-atac_count',
                          'fastqc_r1',
                          'fastqc_r3',
                          'multiqc',
                          'screens_r1',
                          'screens_r3',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R3_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R3_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-arc-2020-A',
                          '/data/refdata-cellranger-atac-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r3'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger-arc': [ '2.0.0' ],
                           'cellranger-atac': [ '2.0.0' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r3'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r3'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r3'],76)

    def test_qcoutputs_10x_cellranger_multi(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi'
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_MC',
                                   ),
                                   cellranger_multi_samples=(
                                       'PJB_CML1',
                                       'PJB_CML2',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '6.1.2' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_10x_cellranger_multi_and_count(self):
        """
        QCOutputs: 10xGenomics CellPlex data with 'cellranger multi' and 'count'
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_GEX_S1_R1_001',
                                       'PJB1_GEX_S1_R2_001',
                                       'PJB2_MC_S2_R1_001',
                                       'PJB2_MC_S2_R2_001',
                                   ),
                                   include_cellranger_count=True,
                                   include_cellranger_multi=True,
                                   cellranger_pipelines=('cellranger',),
                                   cellranger_samples=(
                                       'PJB1_GEX',
                                       'PJB2_MC',
                                   ),
                                   cellranger_multi_samples=(
                                       'PJB_CML1',
                                       'PJB_CML2',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['cellranger_count',
                          'cellranger_multi',
                          'fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_GEX_S1_R1_001',
                          'PJB1_GEX_S1_R2_001',
                          'PJB2_MC_S2_R1_001',
                          'PJB2_MC_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1_GEX','PJB2_MC'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,
                         ['/data/refdata-cellranger-2020-A'])
        self.assertEqual(qc_outputs.multiplexed_samples,
                         ['PJB_CML1','PJB_CML2'])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'cellranger': [ '6.1.2' ],
                           'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_with_non_standard_qc_dir(self):
        """
        QCOutputs: non-standard QC directory
        """
        qc_dir = self._make_qc_dir('my_qc_outs',
                                   fastq_names=(
                                       'PJB1_S1_R1_001',
                                       'PJB1_S1_R2_001',
                                       'PJB2_S2_R1_001',
                                       'PJB2_S2_R2_001',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001',
                          'PJB1_S1_R2_001',
                          'PJB2_S2_R1_001',
                          'PJB2_S2_R2_001'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

    def test_qcoutputs_with_non_canonical_fastq_names(self):
        """
        QCOutputs: non-canonical fastq names
        """
        qc_dir = self._make_qc_dir('qc',
                                   fastq_names=(
                                       'PJB1_S1_R1_001_paired',
                                       'PJB1_S1_R2_001_paired',
                                       'PJB2_S2_R1_001_paired',
                                       'PJB2_S2_R2_001_paired',
                                   ))
        qc_outputs = QCOutputs(qc_dir)
        self.assertEqual(qc_outputs.outputs,
                         ['fastqc_r1',
                          'fastqc_r2',
                          'multiqc',
                          'screens_r1',
                          'screens_r2',
                          'sequence_lengths',
                          'strandedness'])
        self.assertEqual(qc_outputs.fastqs,
                         ['PJB1_S1_R1_001_paired',
                          'PJB1_S1_R2_001_paired',
                          'PJB2_S2_R1_001_paired',
                          'PJB2_S2_R2_001_paired'])
        self.assertEqual(qc_outputs.samples,
                         ['PJB1','PJB2'])
        self.assertEqual(qc_outputs.fastq_screens,
                         ['model_organisms',
                          'other_organisms',
                          'rRNA'])
        self.assertEqual(qc_outputs.cellranger_references,[])
        self.assertEqual(qc_outputs.multiplexed_samples,[])
        self.assertEqual(qc_outputs.reads,['r1','r2'])
        self.assertEqual(qc_outputs.software,
                         { 'fastqc': [ '0.11.3' ],
                           'fastq_screen': [ '0.9.2' ],
                           'fastq_strand': [ '0.0.4' ],
                           'multiqc': [ '1.8' ],
                         })
        self.assertEqual(qc_outputs.stats.max_seqs,37285443)
        self.assertEqual(qc_outputs.stats.min_sequence_length,65)
        self.assertEqual(qc_outputs.stats.max_sequence_length,76)
        self.assertEqual(sorted(
            list(qc_outputs.stats.min_sequence_length_read.keys())),
                         ['r1','r2'])
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r1'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r1'],76)
        self.assertEqual(qc_outputs.stats.min_sequence_length_read['r2'],65)
        self.assertEqual(qc_outputs.stats.max_sequence_length_read['r2'],76)

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
