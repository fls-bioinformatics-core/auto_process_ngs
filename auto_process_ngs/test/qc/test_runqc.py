#######################################################################
# Unit tests for qc/runqc.py
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.mock import MockIlluminaQcSh
from auto_process_ngs.mock import MockMultiQC
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.qc.runqc import RunQC

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestRunQC(unittest.TestCase):
    """
    Tests for RunQC class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestRunQC')
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

    def test_run_qc(self):
        """RunQC: standard QC run
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_run_qc_no_multiqc(self):
        """RunQC: standard QC run (no MultiQC)
        """
        # Make mock illumina_qc.sh
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=False,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd)):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
        for f in ("multiqc_report.html",):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

    def test_run_qc_with_missing_fastq_screen_outputs(self):
        """RunQC: standard QC fails for missing FastQScreen outputs
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"),
                                fastq_screen=False,
                                exit_code=1)
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,1)
        for f in ("qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

    def test_run_qc_with_missing_fastqc_outputs(self):
        """RunQC: standard QC fails for missing FastQC outputs
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"),
                                fastqc=False,
                                exit_code=1)
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,1)
        for f in ("qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

    def test_run_qc_with_missing_multiqc_outputs(self):
        """RunQC: standard QC fails for missing MultiQC outputs
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"),
                           no_outputs=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,1)
        for f in ("qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd)):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
        for f in ("multiqc_report.html",):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)
