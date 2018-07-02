#######################################################################
# Unit tests for qc/illumina_qc.py
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.mock import MockIlluminaQcSh
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.illumina_qc import IlluminaQC
from auto_process_ngs.qc.illumina_qc import fastq_screen_output
from auto_process_ngs.qc.illumina_qc import fastqc_output
from auto_process_ngs.qc.illumina_qc import fastq_strand_output
from auto_process_ngs.qc.illumina_qc import determine_qc_protocol

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestIlluminaQC(unittest.TestCase):
    """
    Tests for IlluminaQC class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestIlluminaQC')
        # Create a temp 'bin' dir
        self.bin = os.path.join(self.wd,"bin")
        os.mkdir(self.bin)
        # Store original PATH
        self.path = os.environ['PATH']

    def tearDown(self):
        # Restore PATH
        os.environ['PATH'] = self.path
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_illumina_qc_version(self):
        """IlluminaQC: fetch version for underlying script
        """
        # Make mock illumina_qc.sh
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"),
                                version="1.3.1")
        illumina_qc = IlluminaQC()
        self.assertEqual(illumina_qc.version(),"1.3.1")

    def test_illumina_qc_check_protocols(self):
        """IlluminaQC: check 'protocol' argument
        """
        # No protocol specified
        self.assertEqual(IlluminaQC().protocol,"standardPE")
        # Standard protocol explicitly specified
        self.assertEqual(IlluminaQC(protocol="standardPE").protocol,
                         "standardPE")
        # Single end
        self.assertEqual(IlluminaQC(protocol="standardSE").protocol,
                         "standardSE")
        # Single cell
        self.assertEqual(IlluminaQC(protocol="singlecell").protocol,
                         "singlecell")
        # Unrecognised
        self.assertRaises(Exception,
                          IlluminaQC,
                          protocol="unknown_protocol")

    def test_illumina_qc_commands_for_single_fastq(self):
        """IlluminaQC: generates default commands for single Fastq
        """
        illumina_qc = IlluminaQC()
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",),
                                    "/path/to/qc")
        self.assertEqual(len(cmds),1)
        self.assertEqual(str(cmds[0]),
                         "illumina_qc.sh "
                         "/path/to/fastqs/test_S1_R1.fastq.gz "
                         "--threads 1 "
                         "--qc_dir /path/to/qc")

    def test_illumina_qc_commands_for_fastq_pair(self):
        """IlluminaQC: generates default commands for Fastq pair
        """
        illumina_qc = IlluminaQC()
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",
                                     "/path/to/fastqs/test_S1_R2.fastq.gz"),
                                    "/path/to/qc")
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

    def test_illumina_qc_commands_for_fastq_pair_with_strandedness(self):
        """IlluminaQC: generates commands for Fastq pair with strandedness
        """
        illumina_qc = IlluminaQC(
            fastq_strand_conf="/path/to/fastq_strand.conf")
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",
                                     "/path/to/fastqs/test_S1_R2.fastq.gz"),
                                    "/path/to/qc")
        self.assertEqual(len(cmds),3)
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
        self.assertEqual(str(cmds[2]),
                         "fastq_strand.py "
                         "-n 1 "
                         "--conf /path/to/fastq_strand.conf "
                         "--outdir /path/to/qc "
                         "/path/to/fastqs/test_S1_R1.fastq.gz "
                         "/path/to/fastqs/test_S1_R2.fastq.gz")

    def test_illumina_qc_commands_for_index_read(self):
        """IlluminaQC: generates empty command list for index read Fastq
        """
        illumina_qc = IlluminaQC()
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_I1.fastq.gz",),
                                    "/path/to/qc")
        self.assertEqual(len(cmds),0)

    def test_illumina_qc_command_with_ungzip_fastqs(self):
        """IlluminaQC: generates command line with --ungzip-fastqs
        """
        illumina_qc = IlluminaQC(ungzip_fastqs=True)
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",),
                                    "/path/to/qc")
        self.assertEqual(len(cmds),1)
        self.assertEqual(str(cmds[0]),
                         "illumina_qc.sh "
                         "/path/to/fastqs/test_S1_R1.fastq.gz "
                         "--ungzip-fastqs "
                         "--threads 1 "
                         "--qc_dir /path/to/qc")

    def test_illumina_qc_command_with_non_default_threads(self):
        """IlluminaQC: generates command line with non-default --threads
        """
        illumina_qc = IlluminaQC(nthreads=8)
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",),
                                    "/path/to/qc")
        self.assertEqual(len(cmds),1)
        self.assertEqual(str(cmds[0]),
                         "illumina_qc.sh "
                         "/path/to/fastqs/test_S1_R1.fastq.gz "
                         "--threads 8 "
                         "--qc_dir /path/to/qc")

    def test_illumina_qc_command_with_fastq_screen_subset(self):
        """IlluminaQC: generates command line with --subset
        """
        illumina_qc = IlluminaQC(fastq_screen_subset=8000)
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",),
                                    "/path/to/qc")
        self.assertEqual(len(cmds),1)
        self.assertEqual(str(cmds[0]),
                         "illumina_qc.sh "
                         "/path/to/fastqs/test_S1_R1.fastq.gz "
                         "--threads 1 "
                         "--subset 8000 "
                         "--qc_dir /path/to/qc")

    def test_illumina_qc_commands_for_single_fastq_SE(self):
        """IlluminaQC: generates default commands for single Fastq (SE)
        """
        illumina_qc = IlluminaQC()
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",),
                                    "/path/to/qc")
        self.assertEqual(len(cmds),1)
        self.assertEqual(str(cmds[0]),
                         "illumina_qc.sh "
                         "/path/to/fastqs/test_S1_R1.fastq.gz "
                         "--threads 1 "
                         "--qc_dir /path/to/qc")

    def test_illumina_qc_commands_for_fastq_pair(self):
        """IlluminaQC: generates default commands for Fastq pair
        """
        illumina_qc = IlluminaQC()
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",
                                     "/path/to/fastqs/test_S1_R2.fastq.gz"),
                                    "/path/to/qc")
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

    def test_illumina_qc_commands_for_SE_with_strandedness(self):
        """IlluminaQC: generates commands for Fastq with strandedness (SE)
        """
        illumina_qc = IlluminaQC(
            protocol="standardSE",
            fastq_strand_conf="/path/to/fastq_strand.conf")
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",),
                                    "/path/to/qc")
        self.assertEqual(len(cmds),2)
        self.assertEqual(str(cmds[0]),
                         "illumina_qc.sh "
                         "/path/to/fastqs/test_S1_R1.fastq.gz "
                         "--threads 1 "
                         "--qc_dir /path/to/qc")
        self.assertEqual(str(cmds[1]),
                         "fastq_strand.py "
                         "-n 1 "
                         "--conf /path/to/fastq_strand.conf "
                         "--outdir /path/to/qc "
                         "/path/to/fastqs/test_S1_R1.fastq.gz")

    def test_illumina_qc_commands_for_fastq_pair_single_cell(self):
        """IlluminaQC: generates default commands for Fastq pair (single cell)
        """
        illumina_qc = IlluminaQC(protocol="singlecell")
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",
                                     "/path/to/fastqs/test_S1_R2.fastq.gz"),
                                    "/path/to/qc")
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

    def test_illumina_qc_commands_for_single_cell_with_strandedness(self):
        """IlluminaQC: generates commands for Fastq pair (single cell, with strand)
        """
        illumina_qc = IlluminaQC(
            protocol="singlecell",
            fastq_strand_conf="/path/to/fastq_strand.conf")
        cmds = illumina_qc.commands(("/path/to/fastqs/test_S1_R1.fastq.gz",
                                     "/path/to/fastqs/test_S1_R2.fastq.gz",),
                                    "/path/to/qc")
        self.assertEqual(len(cmds),3)
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
        self.assertEqual(str(cmds[2]),
                         "fastq_strand.py "
                         "-n 1 "
                         "--conf /path/to/fastq_strand.conf "
                         "--outdir /path/to/qc "
                         "/path/to/fastqs/test_S1_R2.fastq.gz")

    def test_illumina_qc_expected_outputs(self):
        """IlluminaQC: generates correct expected outputs for R1 Fastq
        """
        illumina_qc = IlluminaQC()
        expected_outputs = illumina_qc.expected_outputs(
            "/path/to/fastqs/test_S1_R1.fastq.gz",
            "/path/to/qc")
        reference_outputs = ("/path/to/qc/test_S1_R1_fastqc",
                             "/path/to/qc/test_S1_R1_fastqc.html",
                             "/path/to/qc/test_S1_R1_fastqc.zip",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_rRNA_screen.png",
                             "/path/to/qc/test_S1_R1_rRNA_screen.txt",)
        for e in expected_outputs:
            self.assertTrue(e in reference_outputs,
                            "'%s' shouldn't be predicted" % e)
        for r in reference_outputs:
            self.assertTrue(r in expected_outputs,
                            "'%s' should be predicted" % r)

    def test_illumina_qc_expected_outputs_for_input_pair(self):
        """IlluminaQC: generates correct expected outputs for Fastq R1/R2 pair
        """
        illumina_qc = IlluminaQC()
        expected_outputs = illumina_qc.expected_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",
             "/path/to/fastqs/test_S1_R2.fastq.gz"),
            "/path/to/qc")
        reference_outputs = ("/path/to/qc/test_S1_R1_fastqc",
                             "/path/to/qc/test_S1_R1_fastqc.html",
                             "/path/to/qc/test_S1_R1_fastqc.zip",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_rRNA_screen.png",
                             "/path/to/qc/test_S1_R1_rRNA_screen.txt",
                             "/path/to/qc/test_S1_R2_fastqc",
                             "/path/to/qc/test_S1_R2_fastqc.html",
                             "/path/to/qc/test_S1_R2_fastqc.zip",
                             "/path/to/qc/test_S1_R2_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R2_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R2_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R2_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R2_rRNA_screen.png",
                             "/path/to/qc/test_S1_R2_rRNA_screen.txt",)
        for e in expected_outputs:
            self.assertTrue(e in reference_outputs,
                            "'%s' shouldn't be predicted" % e)
        for r in reference_outputs:
            self.assertTrue(r in expected_outputs,
                            "'%s' should be predicted" % r)

    def test_illumina_qc_expected_outputs_for_input_pair_incl_strand(self):
        """IlluminaQC: generates correct expected outputs for Fastq R1/R2 pair (with strand)
        """
        illumina_qc = IlluminaQC(fastq_strand_conf=True)
        expected_outputs = illumina_qc.expected_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",
             "/path/to/fastqs/test_S1_R2.fastq.gz"),
            "/path/to/qc")
        reference_outputs = ("/path/to/qc/test_S1_R1_fastqc",
                             "/path/to/qc/test_S1_R1_fastqc.html",
                             "/path/to/qc/test_S1_R1_fastqc.zip",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_rRNA_screen.png",
                             "/path/to/qc/test_S1_R1_rRNA_screen.txt",
                             "/path/to/qc/test_S1_R2_fastqc",
                             "/path/to/qc/test_S1_R2_fastqc.html",
                             "/path/to/qc/test_S1_R2_fastqc.zip",
                             "/path/to/qc/test_S1_R2_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R2_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R2_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R2_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R2_rRNA_screen.png",
                             "/path/to/qc/test_S1_R2_rRNA_screen.txt",
                             "/path/to/qc/test_S1_R1_fastq_strand.txt",)
        for e in expected_outputs:
            self.assertTrue(e in reference_outputs,
                            "'%s' shouldn't be predicted" % e)
        for r in reference_outputs:
            self.assertTrue(r in expected_outputs,
                            "'%s' should be predicted" % r)

    def test_illumina_qc_expected_outputs_SE(self):
        """IlluminaQC: generates correct expected outputs for R1 Fastq (SE)
        """
        illumina_qc = IlluminaQC(protocol="standardSE")
        expected_outputs = illumina_qc.expected_outputs(
            "/path/to/fastqs/test_S1_R1.fastq.gz",
            "/path/to/qc")
        reference_outputs = ("/path/to/qc/test_S1_R1_fastqc",
                             "/path/to/qc/test_S1_R1_fastqc.html",
                             "/path/to/qc/test_S1_R1_fastqc.zip",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_rRNA_screen.png",
                             "/path/to/qc/test_S1_R1_rRNA_screen.txt",)
        for e in expected_outputs:
            self.assertTrue(e in reference_outputs,
                            "'%s' shouldn't be predicted" % e)
        for r in reference_outputs:
            self.assertTrue(r in expected_outputs,
                            "'%s' should be predicted" % r)

    def test_illumina_qc_expected_outputs_SE_incl_strand(self):
        """IlluminaQC: generates correct expected outputs (SE, with strand)
        """
        illumina_qc = IlluminaQC(protocol="standardSE",
                                 fastq_strand_conf=True)
        expected_outputs = illumina_qc.expected_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",),
            "/path/to/qc")
        reference_outputs = ("/path/to/qc/test_S1_R1_fastqc",
                             "/path/to/qc/test_S1_R1_fastqc.html",
                             "/path/to/qc/test_S1_R1_fastqc.zip",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_rRNA_screen.png",
                             "/path/to/qc/test_S1_R1_rRNA_screen.txt",
                             "/path/to/qc/test_S1_R1_fastq_strand.txt",)
        for e in expected_outputs:
            self.assertTrue(e in reference_outputs,
                            "'%s' shouldn't be predicted" % e)
        for r in reference_outputs:
            self.assertTrue(r in expected_outputs,
                            "'%s' should be predicted" % r)

    def test_illumina_qc_expected_outputs_for_input_pair_single_cell(self):
        """IlluminaQC: generates correct expected outputs for Fastq R1/R2 pair (single cell)
        """
        illumina_qc = IlluminaQC(protocol="singlecell")
        expected_outputs = illumina_qc.expected_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",
             "/path/to/fastqs/test_S1_R2.fastq.gz"),
            "/path/to/qc")
        reference_outputs = ("/path/to/qc/test_S1_R1_fastqc",
                             "/path/to/qc/test_S1_R1_fastqc.html",
                             "/path/to/qc/test_S1_R1_fastqc.zip",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_rRNA_screen.png",
                             "/path/to/qc/test_S1_R1_rRNA_screen.txt",
                             "/path/to/qc/test_S1_R2_fastqc",
                             "/path/to/qc/test_S1_R2_fastqc.html",
                             "/path/to/qc/test_S1_R2_fastqc.zip",
                             "/path/to/qc/test_S1_R2_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R2_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R2_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R2_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R2_rRNA_screen.png",
                             "/path/to/qc/test_S1_R2_rRNA_screen.txt",)
        for e in expected_outputs:
            self.assertTrue(e in reference_outputs,
                            "'%s' shouldn't be predicted" % e)
        for r in reference_outputs:
            self.assertTrue(r in expected_outputs,
                            "'%s' should be predicted" % r)

    def test_illumina_qc_expected_outputs_for_input_pair_incl_strand_single_cell(self):
        """IlluminaQC: generates correct expected outputs for Fastq R1/R2 pair (single cell, with strand)
        """
        illumina_qc = IlluminaQC(protocol="singlecell",
                                 fastq_strand_conf=True)
        expected_outputs = illumina_qc.expected_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",
             "/path/to/fastqs/test_S1_R2.fastq.gz"),
            "/path/to/qc")
        reference_outputs = ("/path/to/qc/test_S1_R1_fastqc",
                             "/path/to/qc/test_S1_R1_fastqc.html",
                             "/path/to/qc/test_S1_R1_fastqc.zip",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R1_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R1_rRNA_screen.png",
                             "/path/to/qc/test_S1_R1_rRNA_screen.txt",
                             "/path/to/qc/test_S1_R2_fastqc",
                             "/path/to/qc/test_S1_R2_fastqc.html",
                             "/path/to/qc/test_S1_R2_fastqc.zip",
                             "/path/to/qc/test_S1_R2_model_organisms_screen.png",
                             "/path/to/qc/test_S1_R2_model_organisms_screen.txt",
                             "/path/to/qc/test_S1_R2_other_organisms_screen.png",
                             "/path/to/qc/test_S1_R2_other_organisms_screen.txt",
                             "/path/to/qc/test_S1_R2_rRNA_screen.png",
                             "/path/to/qc/test_S1_R2_rRNA_screen.txt",
                             "/path/to/qc/test_S1_R2_fastq_strand.txt",)
        for e in expected_outputs:
            self.assertTrue(e in reference_outputs,
                            "'%s' shouldn't be predicted" % e)
        for r in reference_outputs:
            self.assertTrue(r in expected_outputs,
                            "'%s' should be predicted" % r)

    def test_illumina_qc_expected_outputs_index_read(self):
        """IlluminaQC: predicts no outputs for index read Fastq
        """
        illumina_qc = IlluminaQC()
        expected_outputs = illumina_qc.expected_outputs(
            "/path/to/fastqs/test_S1_I1.fastq.gz",
            "/path/to/qc")
        self.assertEqual(len(expected_outputs),0)

    def test_illumina_qc_expected_outputs_index_read_SE(self):
        """IlluminaQC: predicts no outputs for index read Fastq (SE)
        """
        illumina_qc = IlluminaQC(protocol="standardSE")
        expected_outputs = illumina_qc.expected_outputs(
            "/path/to/fastqs/test_S1_I1.fastq.gz",
            "/path/to/qc")
        self.assertEqual(len(expected_outputs),0)

    def test_illumina_qc_expected_outputs_index_read_single_cell(self):
        """IlluminaQC: predicts no outputs for index read Fastq (single cell)
        """
        illumina_qc = IlluminaQC(protocol="singlecell")
        expected_outputs = illumina_qc.expected_outputs(
            "/path/to/fastqs/test_S1_I1.fastq.gz",
            "/path/to/qc")
        self.assertEqual(len(expected_outputs),0)

    def test_illumina_qc_check_outputs_all_present(self):
        """IlluminaQC: check expected outputs when all present
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",
                             "test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC()
        present,missing = illumina_qc.check_outputs(
            "/path/to/fastqs/test_S1_R1.fastq.gz",qc_dir)
        self.assertEqual(len(missing),0)
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)

    def test_illumina_qc_check_outputs_some_missing(self):
        """IlluminaQC: check expected outputs when some are missing
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",)
        reference_missing = ("test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC()
        present,missing = illumina_qc.check_outputs(
            "/path/to/fastqs/test_S1_R1.fastq.gz",qc_dir)
        # Check present
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)
        # Check missing
        for m in missing:
            self.assertTrue(os.path.dirname(m),qc_dir)
            self.assertTrue(os.path.basename(m) in reference_missing,
                            "'%s' shouldn't be found" % m)
        for r in reference_missing:
            self.assertTrue(os.path.join(qc_dir,r) in missing,
                            "'%s' should exist" % r)

    def test_illumina_qc_check_outputs_all_present_for_pair(self):
        """IlluminaQC: check expected outputs when all present for R1/R2 pair
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",
                             "test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",
                             "test_S1_R2_fastqc",
                             "test_S1_R2_fastqc.html",
                             "test_S1_R2_fastqc.zip",
                             "test_S1_R2_model_organisms_screen.png",
                             "test_S1_R2_model_organisms_screen.txt",
                             "test_S1_R2_other_organisms_screen.png",
                             "test_S1_R2_other_organisms_screen.txt",
                             "test_S1_R2_rRNA_screen.png",
                             "test_S1_R2_rRNA_screen.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC()
        present,missing = illumina_qc.check_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",
             "/path/to/fastqs/test_S1_R2.fastq.gz"),qc_dir)
        print missing
        self.assertEqual(len(missing),0)
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)

    def test_illumina_qc_check_outputs_all_present_for_pair_with_strand(self):
        """IlluminaQC: check expected outputs when all present for R1/R2 pair (with strand)
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",
                             "test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",
                             "test_S1_R2_fastqc",
                             "test_S1_R2_fastqc.html",
                             "test_S1_R2_fastqc.zip",
                             "test_S1_R2_model_organisms_screen.png",
                             "test_S1_R2_model_organisms_screen.txt",
                             "test_S1_R2_other_organisms_screen.png",
                             "test_S1_R2_other_organisms_screen.txt",
                             "test_S1_R2_rRNA_screen.png",
                             "test_S1_R2_rRNA_screen.txt",
                             "test_S1_R1_fastq_strand.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC(fastq_strand_conf=True)
        present,missing = illumina_qc.check_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",
             "/path/to/fastqs/test_S1_R2.fastq.gz"),qc_dir)
        print missing
        self.assertEqual(len(missing),0)
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)

    def test_illumina_qc_check_outputs_for_pair_some_missing(self):
        """IlluminaQC: check expected outputs for R1/R2 pair when some are missing
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",
                             "test_S1_R2_fastqc",
                             "test_S1_R2_fastqc.html",
                             "test_S1_R2_fastqc.zip",
                             "test_S1_R2_model_organisms_screen.png",
                             "test_S1_R2_model_organisms_screen.txt",
                             "test_S1_R2_other_organisms_screen.png",
                             "test_S1_R2_other_organisms_screen.txt",
                             "test_S1_R2_rRNA_screen.png",
                             "test_S1_R2_rRNA_screen.txt",)
        reference_missing = ("test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC()
        present,missing = illumina_qc.check_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",
             "/path/to/fastqs/test_S1_R2.fastq.gz"),
            qc_dir)
        # Check present
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)
        # Check missing
        for m in missing:
            self.assertTrue(os.path.dirname(m),qc_dir)
            self.assertTrue(os.path.basename(m) in reference_missing,
                            "'%s' shouldn't be found" % m)
        for r in reference_missing:
            self.assertTrue(os.path.join(qc_dir,r) in missing,
                            "'%s' should exist" % r)

    def test_illumina_qc_check_outputs_all_present_SE(self):
        """IlluminaQC: check expected outputs when all present (SE)
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",
                             "test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC(protocol="standardSE")
        present,missing = illumina_qc.check_outputs(
            "/path/to/fastqs/test_S1_R1.fastq.gz",qc_dir)
        self.assertEqual(len(missing),0)
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)

    def test_illumina_qc_check_outputs_some_missing_SE(self):
        """IlluminaQC: check expected outputs when some are missing (SE)
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",)
        reference_missing = ("test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC(protocol="standardSE")
        present,missing = illumina_qc.check_outputs(
            "/path/to/fastqs/test_S1_R1.fastq.gz",qc_dir)
        # Check present
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)
        # Check missing
        for m in missing:
            self.assertTrue(os.path.dirname(m),qc_dir)
            self.assertTrue(os.path.basename(m) in reference_missing,
                            "'%s' shouldn't be found" % m)
        for r in reference_missing:
            self.assertTrue(os.path.join(qc_dir,r) in missing,
                            "'%s' should exist" % r)

    def test_illumina_qc_check_outputs_all_present_SE_with_strand(self):
        """IlluminaQC: check expected outputs when all for (SE, with strand)
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",
                             "test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",
                             "test_S1_R1_fastq_strand.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC(protocol="standardSE",
                                 fastq_strand_conf=True)
        present,missing = illumina_qc.check_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",),qc_dir)
        print missing
        self.assertEqual(len(missing),0)
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)

    def test_illumina_qc_check_outputs_all_present_for_pair_single_cell(self):
        """IlluminaQC: check expected outputs when all present for R1/R2 pair (single cell)
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",
                             "test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",
                             "test_S1_R2_fastqc",
                             "test_S1_R2_fastqc.html",
                             "test_S1_R2_fastqc.zip",
                             "test_S1_R2_model_organisms_screen.png",
                             "test_S1_R2_model_organisms_screen.txt",
                             "test_S1_R2_other_organisms_screen.png",
                             "test_S1_R2_other_organisms_screen.txt",
                             "test_S1_R2_rRNA_screen.png",
                             "test_S1_R2_rRNA_screen.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC(protocol="singlecell")
        present,missing = illumina_qc.check_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",
             "/path/to/fastqs/test_S1_R2.fastq.gz"),qc_dir)
        print missing
        self.assertEqual(len(missing),0)
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)

    def test_illumina_qc_check_outputs_all_present_pair_with_strand_single_cell(self):
        """IlluminaQC: check expected outputs when all present for R1/R2 pair (single cell, with strand)
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",
                             "test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",
                             "test_S1_R2_fastqc",
                             "test_S1_R2_fastqc.html",
                             "test_S1_R2_fastqc.zip",
                             "test_S1_R2_model_organisms_screen.png",
                             "test_S1_R2_model_organisms_screen.txt",
                             "test_S1_R2_other_organisms_screen.png",
                             "test_S1_R2_other_organisms_screen.txt",
                             "test_S1_R2_rRNA_screen.png",
                             "test_S1_R2_rRNA_screen.txt",
                             "test_S1_R2_fastq_strand.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC(protocol="singlecell",
                                 fastq_strand_conf=True)
        present,missing = illumina_qc.check_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",
             "/path/to/fastqs/test_S1_R2.fastq.gz"),qc_dir)
        print missing
        self.assertEqual(len(missing),0)
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)

    def test_illumina_qc_check_outputs_for_pair_some_missing_single_cell(self):
        """IlluminaQC: check expected outputs for R1/R2 pair when some are missing (single cell)
        """
        # Make QC dir and (empty) reference files
        qc_dir = os.path.join(self.wd,"qc")
        os.mkdir(qc_dir)
        reference_outputs = ("test_S1_R1_fastqc",
                             "test_S1_R1_fastqc.html",
                             "test_S1_R1_fastqc.zip",
                             "test_S1_R1_model_organisms_screen.png",
                             "test_S1_R1_model_organisms_screen.txt",
                             "test_S1_R1_other_organisms_screen.png",
                             "test_S1_R1_other_organisms_screen.txt",
                             "test_S1_R2_fastqc",
                             "test_S1_R2_fastqc.html",
                             "test_S1_R2_fastqc.zip",
                             "test_S1_R2_model_organisms_screen.png",
                             "test_S1_R2_model_organisms_screen.txt",
                             "test_S1_R2_other_organisms_screen.png",
                             "test_S1_R2_other_organisms_screen.txt",
                             "test_S1_R2_rRNA_screen.png",
                             "test_S1_R2_rRNA_screen.txt",)
        reference_missing = ("test_S1_R1_rRNA_screen.png",
                             "test_S1_R1_rRNA_screen.txt",)
        for r in reference_outputs:
            with open(os.path.join(qc_dir,r),'w') as fp:
                fp.write("test")
        # Get lists of present and missing files
        illumina_qc = IlluminaQC(protocol="singlecell")
        present,missing = illumina_qc.check_outputs(
            ("/path/to/fastqs/test_S1_R1.fastq.gz",
             "/path/to/fastqs/test_S1_R2.fastq.gz"),
            qc_dir)
        # Check present
        for p in present:
            self.assertTrue(os.path.dirname(p),qc_dir)
            self.assertTrue(os.path.basename(p) in reference_outputs,
                            "'%s' shouldn't be found" % p)
        for r in reference_outputs:
            self.assertTrue(os.path.join(qc_dir,r) in present,
                            "'%s' should exist" % r)
        # Check missing
        for m in missing:
            self.assertTrue(os.path.dirname(m),qc_dir)
            self.assertTrue(os.path.basename(m) in reference_missing,
                            "'%s' shouldn't be found" % m)
        for r in reference_missing:
            self.assertTrue(os.path.join(qc_dir,r) in missing,
                            "'%s' should exist" % r)

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

class TestFastqStrandOutputFunction(unittest.TestCase):
    def test_fastq_strand_output(self):
        """fastq_strand_output: handles .fastq file
        """
        self.assertEqual(fastq_strand_output(
            '/data/PB/PB1_ATTAGG_L001_R1_001.fastq'),
                         'PB1_ATTAGG_L001_R1_001_fastq_strand.txt')
    def test_fastq_strand_output_fastqgz(self):
        """fastq_strand_output: handles fastq.gz file
        """
        self.assertEqual(fastq_strand_output(
            '/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz'),
                         'PB1_ATTAGG_L001_R1_001_fastq_strand.txt')
