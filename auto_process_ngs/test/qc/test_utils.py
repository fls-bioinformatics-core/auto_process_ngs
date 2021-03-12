#######################################################################
# Unit tests for qc/utils.py
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.utils import verify_qc
from auto_process_ngs.qc.utils import report_qc
from auto_process_ngs.qc.utils import determine_qc_protocol

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

    def test_determine_qc_protocol_10xchromium3v2_atac_seq(self):
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

    def test_determine_qc_protocol_10xchromium3v2_single_nuclei_atac_seq(self):
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
