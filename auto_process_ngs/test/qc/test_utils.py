#######################################################################
# Unit tests for qc/utils.py
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.utils import mkdirs
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.mock10xdata import ATAC_SUMMARY
from auto_process_ngs.mock10xdata import ATAC_SUMMARY_2_0_0
from auto_process_ngs.mock10xdata import METRICS_SUMMARY
from auto_process_ngs.mock10xdata import CELLPLEX_METRICS_SUMMARY
from auto_process_ngs.mock10xdata import MULTIOME_SUMMARY
from auto_process_ngs.qc.utils import verify_qc
from auto_process_ngs.qc.utils import report_qc
from auto_process_ngs.qc.utils import determine_qc_protocol
from auto_process_ngs.qc.utils import set_cell_count_for_project

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestDetermineQCProtocolFunction(unittest.TestCase):
    """
    Tests for determine_qc_protocol function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestDetermineQCProtocolFunction')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_determine_qc_protocol_standardPE(self):
        """determine_qc_protocol: standard paired-end run
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "standardPE")

    def test_determine_qc_protocol_standardSE(self):
        """determine_qc_protocol: standard single-end run
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",))
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "standardSE")

    def test_determine_qc_protocol_icell8(self):
        """determine_qc_protocol: single-cell run (ICELL8)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "ICELL8"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "singlecell")

    def test_determine_qc_protocol_10xchromium3(self):
        """determine_qc_protocol: single-cell run (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "singlecell")

    def test_determine_qc_protocol_10xchromium3v2(self):
        """determine_qc_protocol: single-cell run (10xGenomics Chromium 3'v2)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v2"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "singlecell")

    def test_determine_qc_protocol_10xchromium3v3(self):
        """determine_qc_protocol: single-cell run (10xGenomics Chromium 3'v3)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v3"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "singlecell")

    def test_determine_qc_protocol_10xchromium3_rna_seq(self):
        """determine_qc_protocol: single-cell RNA-seq (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_scRNAseq")

    def test_determine_qc_protocol_10xchromium3v2_rna_seq(self):
        """determine_qc_protocol: single-cell RNA-seq (10xGenomics Chromium 3'v2)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v2",
                                          'Library type':
                                          "scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_scRNAseq")

    def test_determine_qc_protocol_10xchromium3v3_rna_seq(self):
        """determine_qc_protocol: single-cell RNA-seq (10xGenomics Chromium 3'v3)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v3",
                                          'Library type':
                                          "scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_scRNAseq")

    def test_determine_qc_protocol_10xchromium3_sn_rna_seq(self):
        """determine_qc_protocol: single-nuclei RNA-seq (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "snRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_snRNAseq")

    def test_determine_qc_protocol_10xchromium3v3_sn_rna_seq(self):
        """determine_qc_protocol: single-nuclei RNA-seq (10xGenomics Chromium 3'v3)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v3",
                                          'Library type':
                                          "snRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_snRNAseq")

    def test_determine_qc_protocol_10xchromium3_cellplex(self):
        """determine_qc_protocol: cell multiplexing CellPlex (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "CellPlex"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_CellPlex")

    def test_determine_qc_protocol_10xchromium3_cellplex_scrnaseq(self):
        """determine_qc_protocol: cell multiplexing CellPlex scRNA-seq (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "CellPlex scRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_CellPlex")

    def test_determine_qc_protocol_10xchromium3_cellplex_snrnaseq(self):
        """determine_qc_protocol: cell multiplexing CellPlex snRNA-seq (10xGenomics Chromium 3')
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'",
                                          'Library type':
                                          "CellPlex snRNA-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_CellPlex")

    def test_determine_qc_protocol_10xchromium3v3_cellplex(self):
        """determine_qc_protocol: cell multiplexing CellPlex (10xGenomics Chromium 3'v3)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Chromium 3'v3",
                                          'Library type':
                                          "CellPlex"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_CellPlex")

    def test_determine_qc_protocol_10x_single_cell_atac_seq(self):
        """determine_qc_protocol: single-cell ATAC-seq (10xGenomics Single Cell ATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Single Cell ATAC",
                                          'Library type':
                                          "scATAC-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_scATAC")

    def test_determine_qc_protocol_10x_single_nuclei_atac_seq(self):
        """determine_qc_protocol: single-nuclei ATAC-seq (10xGenomics Single Cell ATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Single Cell ATAC",
                                          'Library type':
                                          "snATAC-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_scATAC")

    def test_determine_qc_protocol_icell8_atac_seq(self):
        """determine_qc_protocol: single-cell ATAC-seq (ICELL8)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "ICELL8",
                                          'Library type':
                                          "scATAC-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "ICELL8_scATAC")

    def test_determine_qc_protocol_10x_visium(self):
        """determine_qc_protocol: spatial RNA-seq run (10xGenomics Visium)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Visium",
                                          'Library type':
                                          "scATAC-seq"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Visium")

    def test_determine_qc_protocol_10x_multiome_atac(self):
        """determine_qc_protocol: single cell multiome ATAC run (10xGenomics Multiome ATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"
                                       "PJB2_S2_R3_001.fastq.gz",
                                       "PJB2_S2_R3_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB1_S1_I2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Single Cell Multiome",
                                          'Library type':
                                          "ATAC"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Multiome_ATAC")

    def test_determine_qc_protocol_10x_multiome_gex(self):
        """determine_qc_protocol: single cell multiome GEX run (10xGenomics Multiome GEX)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB1_S1_I2_001.fastq.gz"),
                                metadata={'Single cell platform':
                                          "10xGenomics Single Cell Multiome",
                                          'Library type':
                                          "GEX"})
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",
                                  os.path.join(self.wd,"PJB"))
        self.assertEqual(determine_qc_protocol(project),
                         "10x_Multiome_GEX")

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
        UpdateAnalysisProject(project).add_qc_outputs()
        # Do reporting
        self.assertEqual(report_qc(project),0)
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
        UpdateAnalysisProject(project).add_qc_outputs()
        # Remove an output
        os.remove(os.path.join(self.wd,
                               "PJB",
                               "qc",
                               "PJB1_S1_R1_001_fastqc.html"))
        # Do reporting
        self.assertEqual(report_qc(project),1)
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
        self.assertEqual(report_qc(project),1)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Found %s (should be missing)" % f)

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

    def test_set_cell_count_for_project(self):
        """
        set_cell_count_for_project: test for scRNA-seq
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "scRNA-seq")
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
        set_cell_count_for_project(project_dir)
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
