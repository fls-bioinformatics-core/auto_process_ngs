#######################################################################
# Unit tests for qc/reporting.py
#######################################################################

import unittest
import os
import tempfile
import shutil
import zipfile
from auto_process_ngs.mock import make_mock_analysis_project
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.reporting import report
from auto_process_ngs.qc.reporting import pretty_print_reads

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

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

    def _make_analysis_project(self,name="PJB",
                               protocol=None,paired_end=True,
                               fastq_dir="fastqs",qc_dir="qc",
                               fastq_names=None,
                               sample_names=None,
                               screens=('model_organisms',
                                        'other_organisms',
                                        'rRNA',),
                               include_fastqc=True,
                               include_fastq_screen=True,
                               include_strandedness=True,
                               include_seqlens=True,
                               include_multiqc=True,
                               include_cellranger_count=False,
                               include_cellranger_multi=False,
                               cellranger_pipelines=('cellranger',),
                               cellranger_samples=None,
                               cellranger_multi_samples=None,
                               legacy_screens=False,
                               legacy_cellranger_outs=False):
        # Create a mock Analysis Project directory
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_QCReport')
            self.top_dir = os.path.join(self.wd,"Test")
            os.mkdir(self.top_dir)
        return make_mock_analysis_project(
            name=name,
            top_dir=self.top_dir,
            protocol=protocol,
            paired_end=paired_end,
            fastq_dir=fastq_dir,
            qc_dir=qc_dir,
            fastq_names=fastq_names,
            sample_names=sample_names,
            screens=screens,
            include_fastqc=include_fastqc,
            include_fastq_screen=include_fastq_screen,
            include_strandedness=include_strandedness,
            include_seqlens=include_seqlens,
            include_multiqc=include_multiqc,
            include_cellranger_count=include_cellranger_count,
            include_cellranger_multi=include_cellranger_multi,
            cellranger_pipelines=cellranger_pipelines,
            cellranger_samples=cellranger_samples,
            cellranger_multi_samples=cellranger_multi_samples,
            legacy_screens=legacy_screens,
            legacy_cellranger_outs=legacy_cellranger_outs)

    def test_report_single_end(self):
        """
        report: single-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=False)
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.SE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.SE.html')))

    def test_report_single_end_no_seq_lens(self):
        """
        report: single-end data: no sequence lengths
        """
        analysis_dir = self._make_analysis_project(protocol='standardSE',
                                                   include_seqlens=False)
        project = AnalysisProject(analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.SE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.SE.html')))

    def test_report_paired_end(self):
        """
        report: paired-end data
        """
        analysis_dir = self._make_analysis_project(paired_end=True)
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_paired_end_no_seq_lens(self):
        """
        report: paired-end data: no sequence lengths
        """
        analysis_dir = self._make_analysis_project(protocol='standardPE',
                                                   include_seqlens=False)
        project = AnalysisProject(analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_paired_end_cellranger_count(self):
        """
        report: paired-end data with cellranger 'count'
        """
        analysis_dir = self._make_analysis_project(
            protocol='10x_scRNAseq',
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger',),
            cellranger_samples=('PJB1','PJB2',))
        project = AnalysisProject(analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_paired_end_cellranger_multi(self):
        """
        report: paired-end data with cellranger 'multi'
        """
        analysis_dir = self._make_analysis_project(
            protocol='10x_CellPlex',
            include_cellranger_multi=True,
            cellranger_multi_samples=('PJB_CML1','PJB_CML2',))
        project = AnalysisProject(analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_paired_end_cellranger_count_and_multi(self):
        """
        report: paired-end data with cellranger 'count' and 'multi'
        """
        analysis_dir = self._make_analysis_project(
            protocol='10x_CellPlex',
            include_cellranger_multi=True,
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger',),
            # NB only GEX samples
            cellranger_samples=('PJB1_GEX',),
            cellranger_multi_samples=('PJB_CML1','PJB_CML2',))
        project = AnalysisProject(analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_paired_end_cellranger_count_multiome(self):
        """
        report: paired-end data with cellranger 'count' (multiome)
        """
        analysis_dir = self._make_analysis_project(
            protocol='10x_Multiome_GEX',
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger-arc',),
            cellranger_samples=('PJB1','PJB2',))
        project = AnalysisProject(analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_paired_end_cellranger_count_multiome_and_scrnaseq(self):
        """
        report: paired-end data with cellranger 'count' (multiome+scRNAseq)
        """
        analysis_dir = self._make_analysis_project(
            protocol='10x_Multiome_GEX',
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger',
                                  'cellranger-arc',),
            cellranger_samples=('PJB1','PJB2',))
        project = AnalysisProject(analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_paired_end_legacy_cellranger_count(self):
        """
        report: paired-end data with cellranger 'count' (legacy)
        """
        analysis_dir = self._make_analysis_project(
            protocol='10x_scRNAseq',
            paired_end=True,
            include_cellranger_count=True,
            cellranger_pipelines=('cellranger',),
            cellranger_samples=('PJB1','PJB2',),
            legacy_cellranger_outs=True)
        project = AnalysisProject(analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_paired_end_with_non_default_fastq_dir(self):
        """
        report: paired-end data with non-default fastq dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_dir="fastqs.non_default")
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_paired_end_with_no_fastq_dir(self):
        """
        report: paired-end data with no fastq dir
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   fastq_dir=".")
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_paired_end_with_non_default_qc_dir(self):
        """
        report: paired-end data with non-default QC dir
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
        """
        report: paired-end data with non-canonical fastq names
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
        """
        report: paired-end data with legacy screen names
        """
        analysis_dir = self._make_analysis_project(paired_end=True,
                                                   legacy_screens=True)
        project = AnalysisProject('PJB',analysis_dir)
        report((project,),filename=os.path.join(self.top_dir,
                                                'report.PE.html'))
        self.assertTrue(os.path.exists(
            os.path.join(self.top_dir,'report.PE.html')))

    def test_report_single_end_multiple_projects(self):
        """
        report: single-end data: two projects in one report
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
        """
        report: single-end data: use data directory
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
        """
        report: single-end data: make ZIP file
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
        """
        report: single-end data: make ZIP file with legacy screens
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
        """
        report: single-end data: make ZIP file with data directory
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
        """
        report: single-end data: fails with two projects in one report (ZIP file/no data directory)
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
        """
        report: single-end data: fails with two projects in one report (duplicated names/ZIP file/no data directory)
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
        """
        report: single-end data: two projects with duplicated names in one report, with ZIP file, with data directory
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
        """
        pretty_print_reads: handles different inputs
        """
        self.assertEqual(pretty_print_reads(1),"1")
        self.assertEqual(pretty_print_reads(12),"12")
        self.assertEqual(pretty_print_reads(117),"117")
        self.assertEqual(pretty_print_reads(1024),"1,024")
        self.assertEqual(pretty_print_reads(33385500),"33,385,500")
        self.assertEqual(pretty_print_reads(112839902),"112,839,902")
        self.assertEqual(pretty_print_reads(10212341927),"10,212,341,927")
