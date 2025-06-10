#######################################################################
# Unit tests for qc/utils.py
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.utils import mkdirs
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import MockMultiQC
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.mock10xdata import ATAC_SUMMARY
from auto_process_ngs.mock10xdata import ATAC_SUMMARY_2_0_0
from auto_process_ngs.mock10xdata import METRICS_SUMMARY
from auto_process_ngs.mock10xdata import CELLPLEX_METRICS_SUMMARY
from auto_process_ngs.mock10xdata import CELLPLEX_METRICS_SUMMARY_7_1_0
from auto_process_ngs.mock10xdata import CELLPLEX_METRICS_SUMMARY_8_0_0
from auto_process_ngs.mock10xdata import CELLPLEX_METRICS_SUMMARY_9_0_0
from auto_process_ngs.mock10xdata import MULTIOME_SUMMARY
from auto_process_ngs.qc.utils import verify_qc
from auto_process_ngs.qc.utils import report_qc
from auto_process_ngs.qc.utils import get_bam_basename
from auto_process_ngs.qc.utils import get_seq_data_samples
from auto_process_ngs.qc.utils import filter_fastqs
from auto_process_ngs.qc.utils import set_cell_count_for_project
from auto_process_ngs.qc.utils import read_versions_file

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestVerifyQCFunction(unittest.TestCase):
    """
    Tests for verify_qc function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestVerifyQCFunction')
        # Create a temp 'bin' dir
        self.bin = os.path.join(self.wd,"bin")
        os.mkdir(self.bin)
        # Store original location
        self.pwd = os.getcwd()
        # Store original PATH
        self.path = os.environ['PATH']
        # Move to working dir
        os.chdir(self.wd)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Restore PATH
        os.environ['PATH'] = self.path
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_verify_qc_all_outputs(self):
        """verify_qc: project with all QC outputs present
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Add QC outputs
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs()
        # Do verification
        self.assertTrue(verify_qc(project))

    def test_verify_qc_incomplete_outputs(self):
        """verify_qc: project with some QC outputs missing
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Add QC outputs
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs()
        # Remove an output
        os.remove(os.path.join(self.wd,
                               "PJB",
                               "qc",
                               "PJB1_S1_R1_001_fastqc.html"))
        # Do verification
        self.assertFalse(verify_qc(project))

    def test_verify_qc_no_outputs(self):
        """verify_qc: project with no QC outputs
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        # Do verification
        self.assertFalse(verify_qc(project))

    def test_verify_qc_custom_protocol(self):
        """verify_qc: project with custom QC protocol
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Add QC outputs
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs()
        # Update QC metadata with custom protocol name and specification
        qc_info = project.qc_info("qc")
        qc_info['protocol'] = "Custom_QC"
        qc_info['protocol_specification'] = "Custom_QC:'Custom QC':seq_reads=[r1,r2]:index_reads=[]:qc_modules=[fastq_screen,fastqc,sequence_lengths]"
        qc_info.save()
        self.assertTrue(verify_qc(project))

class TestReportQCFunction(unittest.TestCase):
    """
    Tests for report_qc function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestReportQCFunction')
        # Create a temp 'bin' dir
        self.bin = os.path.join(self.wd,"bin")
        os.mkdir(self.bin)
        # Store original location
        self.pwd = os.getcwd()
        # Store original PATH
        self.path = os.environ['PATH']
        # Make mock MultiQC
        MockMultiQC.create(os.path.join(self.bin, "multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin, os.environ['PATH'])
        # Move to working dir
        os.chdir(self.wd)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Restore PATH
        os.environ['PATH'] = self.path
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_report_qc_all_outputs(self):
        """report_qc: project with all QC outputs present
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Add QC outputs
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_multiqc=False)
        # Do reporting
        self.assertEqual(report_qc(project, multiqc=True), 0)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_report_qc_incomplete_outputs(self):
        """report_qc: project with some QC outputs missing
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Add QC outputs
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_multiqc=False)
        # Remove an output
        os.remove(os.path.join(self.wd,
                               "PJB",
                               "qc",
                               "PJB1_S1_R1_001_fastqc.html"))
        # Do reporting
        self.assertEqual(report_qc(project, multiqc=True), 1)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_report_qc_no_outputs(self):
        """report_qc: project with no QC outputs
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        # Do reporting
        self.assertEqual(report_qc(project, multiqc=True), 1)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Found %s (should be missing)" % f)

    def test_report_qc_specify_output_dir_qc_dir_outside_project_dir(self):
        """
        report_qc: specify output directory (QC dir outside project)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Add QC outputs
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_multiqc=False)
        # Make an alternative output directory
        out_dir = os.path.join(self.wd, "outs")
        os.mkdir(out_dir)
        # Move the QC directory into the output directory
        os.rename(os.path.join(self.wd, "PJB", "qc"),
                  os.path.join(out_dir, "qc"))
        # Do reporting
        self.assertEqual(report_qc(project,
                                   qc_dir=os.path.join(out_dir, "qc"),
                                   out_dir=out_dir,
                                   multiqc=True), 0)
        # Check output and reports
        for f in ("qc_report.html",
                  "multiqc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            "Missing %s" % f)

    def test_report_qc_specify_output_directory_no_zip(self):
        """
        report_qc: specify the output directory (no ZIP)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Add QC outputs
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_multiqc=False)
        # Make an alternative output directory
        out_dir = os.path.join(self.wd, "outs")
        os.mkdir(out_dir)
        # Do reporting
        self.assertEqual(report_qc(project,
                                   zip_outputs=False,
                                   out_dir=out_dir,
                                   multiqc=True), 0)
        # Check output and reports
        for f in ("qc_report.html",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(out_dir, f)),
                            "Missing %s" % f)

class TestGetBamBasename(unittest.TestCase):
    """
    Tests for the 'get_bam_basename' function
    """
    def test_get_bam_basename(self):
        """
        get_bam_basename: check correct BAM name is returned
        """
        self.assertEqual(get_bam_basename("SM1_S1_L001_R1_001.fastq.gz"),
                         "SM1_S1_L001_001")
        self.assertEqual(get_bam_basename("SM1_S1_R1_001.fastq.gz"),
                         "SM1_S1_001")

class TestGetSeqDataSamples(unittest.TestCase):
    """
    Tests for the 'get_seq_data_samples' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestSetCellCountForProject')
    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)
    def _make_mock_analysis_project(self,library_type,
                                    fastq_names=None,
                                    single_cell_platform=None,
                                    bio_samples=None):
        # Default Fastqs
        if not fastq_names:
            fastq_names = ("PJB1_S1_L001_R1_001.fastq.gz",
                           "PJB1_S1_L001_R2_001.fastq.gz",
                           "PJB2_S2_L001_R1_001.fastq.gz",
                           "PJB2_S2_L001_R2_001.fastq.gz",)
        # Create a mock AnalysisProject
        m = MockAnalysisProject('PJB',
                                fastq_names=fastq_names,
                                metadata={'Single cell platform':
                                          single_cell_platform,
                                          'Library type': library_type,
                                          'Biological samples': bio_samples})
        m.create(top_dir=self.wd)
        return os.path.join(self.wd,'PJB')

    def test_get_seq_data_samples(self):
        """
        get_seq_data_samples: standard analysis project
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project("RNA-seq")
        # Check sequence data samples
        self.assertEqual(get_seq_data_samples(project_dir),
                         ["PJB1","PJB2"])

    def test_get_seq_data_samples_defined_in_metadata(self):
        """
        get_seq_data_samples: biological samples defined in project metadata
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project("RNA-seq",
                                                       bio_samples="PJB1")
        # Check sequence data samples
        self.assertEqual(get_seq_data_samples(project_dir),
                         ["PJB1"])

    def test_get_seq_data_samples_10x_cellplex(self):
        """
        get_seq_data_samples: 10x CellPlex project
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "CellPlex",
            single_cell_platform="10xGenomics Chromium 3'v3")
        # Make 10x_multi_config.csv file
        with open(os.path.join(project_dir,"10x_multi_config.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1,{fastq_dir},any,PJB1,gene expression,
PJB2,{fastq_dir},any,PJB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBA,CMO301,PBA
PBB,CMO302,PBB
""".format(fastq_dir=os.path.join(project_dir,'fastqs')))
        # Check sequence data samples
        self.assertEqual(get_seq_data_samples(project_dir),
                         ["PJB1"])

    def test_get_seq_data_samples_10x_cellplex_multiple_configs(self):
        """
        get_seq_data_samples: 10x CellPlex project (multiple configs)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "CellPlex",
            single_cell_platform="10xGenomics Chromium 3'v3",
            fastq_names=("PJB1_GEX_S1_L001_R1_001.fastq.gz",
                         "PJB1_GEX_S1_L001_R2_001.fastq.gz",
                         "PJB1_CML_S2_L001_R1_001.fastq.gz",
                         "PJB1_CML_S2_L001_R2_001.fastq.gz",
                         "PJB2_GEX_S3_L001_R1_001.fastq.gz",
                         "PJB2_GEX_S3_L001_R2_001.fastq.gz",
                         "PJB2_CML_S4_L001_R1_001.fastq.gz",
                         "PJB2_CML_S4_L001_R2_001.fastq.gz",))
        # Make 10x_multi_config.csv files
        with open(os.path.join(project_dir,"10x_multi_config.PJB1.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,{fastq_dir},any,PJB1,gene expression,
PJB1_CML,{fastq_dir},any,PJB1,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBA,CMO301,PBA
PBB,CMO302,PBB
""".format(fastq_dir=os.path.join(project_dir,'fastqs')))
        with open(os.path.join(project_dir,"10x_multi_config.PJB2.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB2_GEX,{fastq_dir},any,PJB2,gene expression,
PJB2_CML,{fastq_dir},any,PJB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBC,CMO303,PBC
PBD,CMO304,PBD
""".format(fastq_dir=os.path.join(project_dir,'fastqs')))
        # Check sequence data samples
        self.assertEqual(get_seq_data_samples(project_dir),
                         ["PJB1_GEX", "PJB2_GEX"])

    def test_get_seq_data_samples_10x_flex(self):
        """
        get_seq_data_samples: 10x Flex project
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "Flex",
            single_cell_platform="10xGenomics Chromium 3'v3",
            fastq_names=("PJB1_S1_L001_R1_001.fastq.gz",
                         "PJB1_S1_L001_R2_001.fastq.gz",))
        # Make 10x_multi_config.csv file
        with open(os.path.join(project_dir,"10x_multi_config.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1,{fastq_dir},any,PJB1,gene expression,

[samples]
sample_id,probe_barcode_ids,description
PBA,BC001,PBA
PBB,BC002,PBB
""".format(fastq_dir=os.path.join(project_dir,'fastqs')))
        # Check sequence data samples
        self.assertEqual(get_seq_data_samples(project_dir),
                         ["PJB1"])

    def test_get_seq_data_samples_10x_immune_profiling(self):
        """
        get_seq_data_samples: 10x Single Cell Immune Profiling project
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "Single Cell Immune Profiling",
            single_cell_platform="10xGenomics Chromium 5'",
            fastq_names=("PJB1_GEX_S1_L001_R1_001.fastq.gz",
                         "PJB1_GEX_S1_L001_R2_001.fastq.gz",
                         "PJB1_TCR_S2_L001_R1_001.fastq.gz",
                         "PJB1_TCR_S2_L001_R2_001.fastq.gz",
                         "PJB2_GEX_S3_L001_R1_001.fastq.gz",
                         "PJB2_GEX_S3_L001_R2_001.fastq.gz",
                         "PJB2_TCC_S4_L001_R1_001.fastq.gz",
                         "PJB2_TCR_S4_L001_R2_001.fastq.gz",))
        # Make 10x_multi_config.csv files
        with open(os.path.join(project_dir,"10x_multi_config.PJB1.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,{fastq_dir},any,PJB1,gene expression,
PJB1_TCR,{fastq_dir},any,PJB1,VDJ-T,
""".format(fastq_dir=os.path.join(project_dir,'fastqs')))
        with open(os.path.join(project_dir,"10x_multi_config.PJB2.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB2_GEX,{fastq_dir},any,PJB2,gene expression,
PJB2_TCR,{fastq_dir},any,PJB2,VDJ-T,
""".format(fastq_dir=os.path.join(project_dir,'fastqs')))
        # Check sequence data samples
        self.assertEqual(get_seq_data_samples(project_dir),
                         ["PJB1_GEX","PJB1_TCR","PJB2_GEX","PJB2_TCR"])

    def test_get_seq_data_samples_bad_cellranger_multi_config(self):
        """
        get_seq_data_samples: exception for 'bad' Cellranger multi config file
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "CellPlex",
            single_cell_platform="10xGenomics Chromium 3'v3")
        # Make 10x_multi_config.csv file
        with open(os.path.join(project_dir,"10x_multi_config.csv"),
                  'wt') as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1,{fastq_dir},any,PJB1,[Gene expression|Multiplexing Capture],
PJB2,{fastq_dir},any,PJB2,[Gene expression|Multiplexing Capture],

[samples]
sample_id,cmo_ids,description
MULTIPLEXED_SAMPLE,CMO301|CMO302|...,DESCRIPTION
PBA,CMO301,PBA
PBB,CMO302,PBB
""".format(fastq_dir=os.path.join(project_dir,'fastqs')))
        # Raises exception when fetching sequence data samples
        self.assertRaises(Exception,
                          get_seq_data_samples,
                          project_dir)

class TestFilterFastqs(unittest.TestCase):

    def test_filter_fastqs(self):
        """
        filter_fastqs: check Fastq names are correctly filtered
        """
        fastqs = ("PJB1_S1_R1_001.fastq.gz",
                  "PJB1_S1_R2_001.fastq.gz",
                  "PJB1_S1_R3_001.fastq.gz",
                  "PJB1_S1_I1_001.fastq.gz",
                  "PJB2_S2_R1_001.fastq.gz",
                  "PJB2_S2_R2_001.fastq.gz",
                  "PJB2_S2_R3_001.fastq.gz",
                  "PJB2_S2_I1_001.fastq.gz",)
        # Filter R1
        self.assertEqual(filter_fastqs(['r1'],fastqs),
                         ["PJB1_S1_R1_001",
                          "PJB2_S2_R1_001",])
        # Filter R1 & R3
        self.assertEqual(filter_fastqs(['r1','r3'],fastqs),
                         ["PJB1_S1_R1_001",
                          "PJB1_S1_R3_001",
                          "PJB2_S2_R1_001",
                          "PJB2_S2_R3_001",])
        # Filter I1
        self.assertEqual(filter_fastqs(['i1'],fastqs),
                         ["PJB1_S1_I1_001",
                          "PJB2_S2_I1_001",])
        # Filter R*
        self.assertEqual(filter_fastqs(['r*'],fastqs),
                         ["PJB1_S1_R1_001",
                          "PJB1_S1_R2_001",
                          "PJB1_S1_R3_001",
                          "PJB2_S2_R1_001",
                          "PJB2_S2_R2_001",
                          "PJB2_S2_R3_001",])
        # Filter *
        self.assertEqual(filter_fastqs(['*'],fastqs),
                         ["PJB1_S1_I1_001",
                          "PJB1_S1_R1_001",
                          "PJB1_S1_R2_001",
                          "PJB1_S1_R3_001",
                          "PJB2_S2_I1_001",
                          "PJB2_S2_R1_001",
                          "PJB2_S2_R2_001",
                          "PJB2_S2_R3_001",])
        # Filter everything
        self.assertEqual(filter_fastqs([],fastqs),[])

class TestSetCellCountForProject(unittest.TestCase):
    """
    Tests for the 'set_cell_count_for_project' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestSetCellCountForProject')
    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)
    def _make_mock_analysis_project(self,single_cell_platform,library_type):
        # Create a mock AnalysisProject
        m = MockAnalysisProject('PJB',
                                fastq_names=("PJB1_S1_L001_R1_001.fastq.gz",
                                             "PJB1_S1_L001_R2_001.fastq.gz",),
                                metadata={'Single cell platform':
                                          single_cell_platform,
                                          'Library type': library_type,})
        m.create(top_dir=self.wd)
        return os.path.join(self.wd,'PJB')

    def test_set_cell_count_for_project_chromium_3v3(self):
        """
        set_cell_count_for_project: test for scRNA-seq (10x Chromium 3'v3)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "scRNA-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "8.0.0",
                                  "refdata-gex-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'wt') as fp:
            fp.write(METRICS_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-gex-GRCh38-2020-A
Cellranger version\t8.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)

    def test_set_cell_count_for_project_chromium_gem_x_3v4(self):
        """
        set_cell_count_for_project: test for scRNA-seq (10x Chromium GEM-X 3'v4)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium GEM-X 3'v4",
            "scRNA-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "8.0.0",
                                  "refdata-gex-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'wt') as fp:
            fp.write(METRICS_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-gex-GRCh38-2020-A
Cellranger version\t8.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)

    def test_set_cell_count_for_project_chromium_next_gem(self):
        """
        set_cell_count_for_project: test for scRNA-seq (10x Chromium Next GEM)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium Next GEM",
            "scRNA-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "8.0.0",
                                  "refdata-gex-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'wt') as fp:
            fp.write(METRICS_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-gex-GRCh38-2020-A
Cellranger version\t8.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)

    def test_set_cell_count_for_project_chromium_next_gem_3v31(self):
        """
        set_cell_count_for_project: test for scRNA-seq (10x Chromium Next GEM 3'v3.1)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium Next GEM 3'v3.1",
            "scRNA-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "8.0.0",
                                  "refdata-gex-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'wt') as fp:
            fp.write(METRICS_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-gex-GRCh38-2020-A
Cellranger version\t8.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)

    def test_set_cell_count_for_atac_project(self):
        """
        set_cell_count_for_project: test for scATAC-seq
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Single Cell ATAC",
            "scATAC-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "1.2.0",
                                  "refdata-cellranger-atac-GRCh38-1.2.0",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        summary_file = os.path.join(counts_dir,
                                    "summary.csv")
        with open(summary_file,'w') as fp:
            fp.write(ATAC_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-atac-GRCh38-1.2.0
Cellranger version\t1.2.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         5682)

    def test_set_cell_count_for_single_nuclei_atac_project(self):
        """
        set_cell_count_for_project: test for snATAC-seq
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Single Cell ATAC",
            "snATAC-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "1.2.0",
                                  "refdata-cellranger-atac-GRCh38-1.2.0",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        summary_file = os.path.join(counts_dir,
                                            "summary.csv")
        with open(summary_file,'w') as fp:
            fp.write(ATAC_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-atac-GRCh38-1.2.0
Cellranger version\t1.2.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         5682)

    def test_set_cell_count_for_atac_project_2_0_0(self):
        """
        set_cell_count_for_project: test for scATAC-seq (Cellranger ATAC 2.0.0)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Single Cell ATAC",
            "scATAC-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "2.0.0",
                                  "refdata-cellranger-atac-GRCh38-2020-A-2.0.0",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        summary_file = os.path.join(counts_dir,
                                    "summary.csv")
        with open(summary_file,'w') as fp:
            fp.write(ATAC_SUMMARY_2_0_0)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0
Cellranger version\t2.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         3582)

    def test_set_cell_count_for_multiome_atac_project(self):
        """
        set_cell_count_for_project: test for single cell multiome ATAC
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Single Cell Multiome",
            "ATAC")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "1.0.0",
                                  "refdata-cellranger-arc-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        summary_file = os.path.join(counts_dir,
                                    "summary.csv")
        with open(summary_file,'w') as fp:
            fp.write(MULTIOME_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-arc-GRCh38-2020-A
Cellranger version\t1.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         744)

    def test_set_cell_count_for_multiome_gex_project(self):
        """
        set_cell_count_for_project: test for single cell multiome GEX
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Single Cell Multiome",
            "GEX")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "1.0.0",
                                  "refdata-cellranger-arc-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        summary_file = os.path.join(counts_dir,
                                    "summary.csv")
        with open(summary_file,'w') as fp:
            fp.write(MULTIOME_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-arc-GRCh38-2020-A
Cellranger version\t1.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         744)

    def test_set_cell_count_for_cellplex_project(self):
        """
        set_cell_count_for_project: test for multiplexed data (CellPlex)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "CellPlex")
        # Build mock cellranger multi output directory
        multi_dir = os.path.join(project_dir,
                                 "qc",
                                 "cellranger_multi",
                                 "6.0.0",
                                 "refdata-cellranger-gex-GRCh38-2020-A",
                                 "outs")
        mkdirs(multi_dir)
        for sample in ("PBA","PBB",):
            sample_dir = os.path.join(multi_dir,
                                      "per_sample_outs",
                                      sample)
            mkdirs(sample_dir)
            summary_file = os.path.join(sample_dir,
                                        "metrics_summary.csv")
            with open(summary_file,'wt') as fp:
                fp.write(CELLPLEX_METRICS_SUMMARY)
            web_summary = os.path.join(sample_dir,
                                       "web_summary.html")
            with open(web_summary,'wt') as fp:
                fp.write("Placeholder for web_summary.html\n")
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-gex-GRCh38-2020-A
Cellranger version\t6.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="multi")
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         10350)

    def test_set_cell_count_for_cellplex_project_multiple_physical_samples(self):
        """
        set_cell_count_for_project: test for multiplexed data (CellPlex, multiple physical samples)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "CellPlex")
        # Build mock cellranger multi output directory
        multi_dir = os.path.join(project_dir,
                                 "qc",
                                 "cellranger_multi",
                                 "9.0.0",
                                 "refdata-cellranger-gex-GRCh38-2020-A")
        mkdirs(multi_dir)
        for sample in ("PJB1", "PJB2"):
            if sample == "PJB1":
                multiplexed_samples = ("PBA", "PBB")
            elif sample == "PJB2":
                multiplexed_samples = ("PBC", "PBD")
            for multiplexed_sample in multiplexed_samples:
                sample_dir = os.path.join(multi_dir,
                                          sample,
                                          "outs",
                                          "per_sample_outs",
                                          multiplexed_sample)
                mkdirs(sample_dir)
                summary_file = os.path.join(sample_dir,
                                            "metrics_summary.csv")
                with open(summary_file,'wt') as fp:
                    fp.write(CELLPLEX_METRICS_SUMMARY)
                web_summary = os.path.join(sample_dir,
                                           "web_summary.html")
                with open(web_summary,'wt') as fp:
                    fp.write("Placeholder for web_summary.html\n")
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-gex-GRCh38-2020-A
Cellranger version\t9.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="multi")
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         20700)

    def test_set_cell_count_for_cellplex_project_710(self):
        """
        set_cell_count_for_project: test for multiplexed data (CellPlex 7.1.0)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "CellPlex")
        # Build mock cellranger multi output directory
        multi_dir = os.path.join(project_dir,
                                 "qc",
                                 "cellranger_multi",
                                 "7.1.0",
                                 "refdata-cellranger-gex-GRCh38-2020-A",
                                 "outs")
        mkdirs(multi_dir)
        for sample in ("PBA","PBB",):
            sample_dir = os.path.join(multi_dir,
                                      "per_sample_outs",
                                      sample)
            mkdirs(sample_dir)
            summary_file = os.path.join(sample_dir,
                                        "metrics_summary.csv")
            with open(summary_file,'wt') as fp:
                fp.write(CELLPLEX_METRICS_SUMMARY_7_1_0)
            web_summary = os.path.join(sample_dir,
                                       "web_summary.html")
            with open(web_summary,'wt') as fp:
                fp.write("Placeholder for web_summary.html\n")
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-gex-GRCh38-2020-A
Cellranger version\t7.1.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="multi")
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         3138)

    def test_set_cell_count_for_cellplex_project_800(self):
        """
        set_cell_count_for_project: test for multiplexed data (CellPlex 8.0.0)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "CellPlex")
        # Build mock cellranger multi output directory
        multi_dir = os.path.join(project_dir,
                                 "qc",
                                 "cellranger_multi",
                                 "8.0.0",
                                 "refdata-cellranger-gex-GRCh38-2020-A",
                                 "outs")
        mkdirs(multi_dir)
        for sample in ("PBA","PBB",):
            sample_dir = os.path.join(multi_dir,
                                      "per_sample_outs",
                                      sample)
            mkdirs(sample_dir)
            summary_file = os.path.join(sample_dir,
                                        "metrics_summary.csv")
            with open(summary_file,'wt') as fp:
                fp.write(CELLPLEX_METRICS_SUMMARY_8_0_0)
            web_summary = os.path.join(sample_dir,
                                       "web_summary.html")
            with open(web_summary,'wt') as fp:
                fp.write("Placeholder for web_summary.html\n")
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-gex-GRCh38-2020-A
Cellranger version\t8.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="multi")
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         3164)

    def test_set_cell_count_for_cellplex_project_next_gem(self):
        """
        set_cell_count_for_project: test for multiplexed data (CellPlex with 10x Next GEM)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium Next GEM",
            "CellPlex")
        # Build mock cellranger multi output directory
        multi_dir = os.path.join(project_dir,
                                 "qc",
                                 "cellranger_multi",
                                 "8.0.0",
                                 "refdata-cellranger-gex-GRCh38-2020-A",
                                 "outs")
        mkdirs(multi_dir)
        for sample in ("PBA","PBB",):
            sample_dir = os.path.join(multi_dir,
                                      "per_sample_outs",
                                      sample)
            mkdirs(sample_dir)
            summary_file = os.path.join(sample_dir,
                                        "metrics_summary.csv")
            with open(summary_file,'wt') as fp:
                fp.write(CELLPLEX_METRICS_SUMMARY_8_0_0)
            web_summary = os.path.join(sample_dir,
                                       "web_summary.html")
            with open(web_summary,'wt') as fp:
                fp.write("Placeholder for web_summary.html\n")
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-gex-GRCh38-2020-A
Cellranger version\t8.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="multi")
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         3164)

    def test_set_cell_count_for_cellplex_project_next_gem_3v31(self):
        """
        set_cell_count_for_project: test for multiplexed data (CellPlex with 10x Next GEM 3'v3.1)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium Next GEM 3'v3.1",
            "CellPlex")
        # Build mock cellranger multi output directory
        multi_dir = os.path.join(project_dir,
                                 "qc",
                                 "cellranger_multi",
                                 "8.0.0",
                                 "refdata-cellranger-gex-GRCh38-2020-A",
                                 "outs")
        mkdirs(multi_dir)
        for sample in ("PBA","PBB",):
            sample_dir = os.path.join(multi_dir,
                                      "per_sample_outs",
                                      sample)
            mkdirs(sample_dir)
            summary_file = os.path.join(sample_dir,
                                        "metrics_summary.csv")
            with open(summary_file,'wt') as fp:
                fp.write(CELLPLEX_METRICS_SUMMARY_8_0_0)
            web_summary = os.path.join(sample_dir,
                                       "web_summary.html")
            with open(web_summary,'wt') as fp:
                fp.write("Placeholder for web_summary.html\n")
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-gex-GRCh38-2020-A
Cellranger version\t8.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="multi")
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         3164)

    def test_set_cell_count_for_cellplex_project_900(self):
        """
        set_cell_count_for_project: test for multiplexed data (CellPlex 9.0.0)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "CellPlex")
        # Build mock cellranger multi output directory
        multi_dir = os.path.join(project_dir,
                                 "qc",
                                 "cellranger_multi",
                                 "9.0.0",
                                 "refdata-cellranger-gex-GRCh38-2020-A",
                                 "outs")
        mkdirs(multi_dir)
        for sample in ("PBA","PBB",):
            sample_dir = os.path.join(multi_dir,
                                      "per_sample_outs",
                                      sample)
            mkdirs(sample_dir)
            summary_file = os.path.join(sample_dir,
                                        "metrics_summary.csv")
            with open(summary_file,'wt') as fp:
                fp.write(CELLPLEX_METRICS_SUMMARY_9_0_0)
            web_summary = os.path.join(sample_dir,
                                       "web_summary.html")
            with open(web_summary,'wt') as fp:
                fp.write("Placeholder for web_summary.html\n")
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-gex-GRCh38-2020-A
Cellranger version\t9.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="multi")
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         3142)

    def test_set_cell_count_for_cellplex_project_with_count(self):
        """
        set_cell_count_for_project: test for multiplexed data (CellPlex) with count output
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "CellPlex")
        # Build mock cellranger multi output directory
        multi_dir = os.path.join(project_dir,
                                 "qc",
                                 "cellranger_multi",
                                 "6.0.0",
                                 "refdata-cellranger-gex-GRCh38-2020-A",
                                 "outs")
        mkdirs(multi_dir)
        for sample in ("PBA","PBB",):
            sample_dir = os.path.join(multi_dir,
                                      "per_sample_outs",
                                      sample)
            mkdirs(sample_dir)
            summary_file = os.path.join(sample_dir,
                                        "metrics_summary.csv")
            with open(summary_file,'wt') as fp:
                fp.write(CELLPLEX_METRICS_SUMMARY)
            web_summary = os.path.join(sample_dir,
                                       "web_summary.html")
            with open(web_summary,'wt') as fp:
                fp.write("Placeholder for web_summary.html\n")
        # Build mock cellranger count output directory in parallel
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "6.0.0",
                                  "refdata-cellranger-gex-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'wt') as fp:
            fp.write(METRICS_SUMMARY)
        web_summary = os.path.join(counts_dir,
                                   "web_summary.html")
        with open(web_summary,'wt') as fp:
            fp.write("Placeholder for web_summary.html\n")
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-gex-GRCh38-2020-A
Cellranger version\t6.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject(project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="multi")
        # Check updated cell count
        self.assertEqual(AnalysisProject(project_dir).info.number_of_cells,
                         10350)
        # Update the cell counts from "count" outputs
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="count")
        # Check updated cell count
        self.assertEqual(AnalysisProject(project_dir).info.number_of_cells,
                         2272)

    def test_set_cell_count_for_immune_profiling_project(self):
        """
        set_cell_count_for_project: test for single cell immune profiling data (Chromium 5')
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 5'",
            "Single Cell Immune Profiling")
        # Build mock cellranger multi output directory
        multi_dir = os.path.join(project_dir,
                                 "qc",
                                 "cellranger_multi",
                                 "9.0.0",
                                 "refdata-cellranger-gex-GRCh38-2020-A",
                                 "outs",
                                 "per_sample_outs",
                                 "PJB1")
        mkdirs(multi_dir)
        summary_file = os.path.join(multi_dir, "metrics_summary.csv")
        with open(summary_file,'wt') as fp:
            fp.write(CELLPLEX_METRICS_SUMMARY)
        web_summary = os.path.join(multi_dir, "web_summary.html")
        with open(web_summary,'wt') as fp:
            fp.write("Placeholder for web_summary.html\n")
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-gex-GRCh38-2020-A
Cellranger version\t9.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="multi")
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         5175)

    def test_set_cell_count_for_immune_profiling_project_multiple_physical_samples(self):
        """
        set_cell_count_for_project: test for single cell immune profiling data (Chromium 5', multiple physical samples)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 5'",
            "Single Cell Immune Profiling")
        # Build mock cellranger multi output directory
        multi_dir = os.path.join(project_dir,
                                 "qc",
                                 "cellranger_multi",
                                 "9.0.0",
                                 "refdata-cellranger-gex-GRCh38-2020-A")
        mkdirs(multi_dir)
        for sample in ("PJB1", "PJB2"):
            sub_dir = os.path.join(multi_dir,
                                   sample,
                                   "outs",
                                   "per_sample_outs",
                                   sample)
            mkdirs(sub_dir)
            summary_file = os.path.join(sub_dir, "metrics_summary.csv")
            with open(summary_file,'wt') as fp:
                fp.write(CELLPLEX_METRICS_SUMMARY)
            web_summary = os.path.join(sub_dir, "web_summary.html")
            with open(web_summary,'wt') as fp:
                fp.write("Placeholder for web_summary.html\n")
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-gex-GRCh38-2020-A
Cellranger version\t9.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir,source="multi")
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         10350)

    def test_set_cell_count_project_missing_library_type(self):
        """
        set_cell_count_for_project: test for scRNA-seq when library not set
        """
        # Set up mock project with library type not set
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            None)
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "5.0.1",
                                  "refdata-gex-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'w') as fp:
            fp.write(METRICS_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-gex-GRCh38-2020-A
Cellranger version\t5.0.1
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)

    def test_set_cell_count_for_project_no_subdirs(self):
        """
        set_cell_count_for_project: test for scRNA-seq (old-style output)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "scRNA-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'w') as fp:
            fp.write(METRICS_SUMMARY)
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)

    def test_set_cell_count_project_missing_library_type_no_subdirs(self):
        """
        set_cell_count_for_project: test for scRNA-seq when library not set (old-style output)
        """
        # Set up mock project with library type not set
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            None)
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'w') as fp:
            fp.write(METRICS_SUMMARY)
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)

    def test_set_cell_count_fails_for_project_with_no_metadata(self):
        """
        set_cell_count_for_project: raises exception for project with no metadata
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(None,None)
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "5.0.1",
                                  "refdata-gex-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'wt') as fp:
            fp.write(METRICS_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-gex-GRCh38-2020-A
Cellranger version\t5.0.1
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Attempting to update the cell counts should raise
        # NotImplementedError
        self.assertRaises(NotImplementedError,
                          set_cell_count_for_project,
                          project_dir)
        # Check cell count wasn't updated
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)

class TestReadVersionsFile(unittest.TestCase):
    """
    Tests for the 'read_versions_file' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestReadVersionsFile')
    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)
    def test_read_versions_file(self):
        """
        read_versions_file: get data from a file
        """
        # Make example '_versions' file
        versions_file = os.path.join(self.wd,"_versions")
        with open(versions_file,'wt') as fp:
            fp.write("star\t2.7.7a\nsamtools\t1.15.1\n")
        # Read version information from file
        self.assertEqual(read_versions_file(versions_file),
                         {
                             "star": ['2.7.7a'],
                             "samtools": ['1.15.1']
                         })
    def test_read_versions_file_add_to_existing_data(self):
        """
        read_versions_file: get data from a file (add to existing data)
        """
        # Make example '_versions' file
        versions_file = os.path.join(self.wd,"_versions")
        with open(versions_file,'wt') as fp:
            fp.write("star\t2.7.7a\nsamtools\t1.15.1\n")
        # Existing packages data
        pkgs = {
            "seqtk": ['1.3'],
            "star": ['2.4.2a']
        }
        # Read version information from file
        self.assertEqual(read_versions_file(versions_file,pkgs),
                         {
                             "star": ['2.4.2a',
                                      '2.7.7a'],
                             "samtools": ['1.15.1'],
                             "seqtk": ['1.3']
                         })
