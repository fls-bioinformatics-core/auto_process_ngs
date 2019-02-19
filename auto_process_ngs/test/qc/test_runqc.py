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
from auto_process_ngs.mock import MockFastqStrandPy
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.runqc import RunQC
from auto_process_ngs.qc.illumina_qc import IlluminaQC

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
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_run_qc_with_strandedness(self):
        """RunQC: standard QC run with strandedness determination
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        illumina_qc=IlluminaQC(fastq_strand_conf="fastq_strand.conf")
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          illumina_qc=illumina_qc)
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_run_qc_with_missing_strandedness(self):
        """RunQC: standard QC fails with missing strandedness outputs
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"),
                                 no_outputs=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        illumina_qc=IlluminaQC(fastq_strand_conf="fastq_strand.conf")
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          illumina_qc=illumina_qc)
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,1)
        self.assertTrue(os.path.exists(os.path.join(self.wd,"PJB","qc")),
                        "Missing 'qc'")
        for f in ("qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

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
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=False,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
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
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,1)
        self.assertTrue(os.path.exists(os.path.join(self.wd,"PJB","qc")),
                        "Missing 'qc'")
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
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,1)
        self.assertTrue(os.path.exists(os.path.join(self.wd,"PJB","qc")),
                        "Missing 'qc'")
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
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,1)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd)):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
        for f in ("multiqc_report.html",):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

    def test_run_qc_non_default_fastq_dir(self):
        """RunQC: standard QC run using non-default Fastq dir
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",
                                fastq_names=("PJB1_S1_R1_001.fastq.gz",
                                             "PJB1_S1_R2_001.fastq.gz"),
                                fastq_dir="fastqs.cells")
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          fastq_dir="fastqs.cells")
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_run_qc_non_default_output_dir(self):
        """RunQC: standard QC run using non-default output dir
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
                                          os.path.join(self.wd,"PJB")),
                          qc_dir="qc.non_default")
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                     "PJB","qc")),
                         "'qc' exists, but shouldn't")
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "PJB",
                                                   "qc.non_default")),
                         "'qc' directory doesn't exist, but should")
        for f in ("qc.non_default_report.html",
                  "qc.non_default_report.PJB.%s.zip" %
                  os.path.basename(self.wd),
                  "multiqc.non_default_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_run_qc_single_end(self):
        """RunQC: standard QC run (single-end data)
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_run_qc_multiple_projects(self):
        """RunQC: standard QC run (multiple projects)
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis projects
        p = MockAnalysisProject("AB",("AB1_S1_R1_001.fastq.gz",
                                      "AB1_S1_R2_001.fastq.gz",
                                      "AB2_S2_R1_001.fastq.gz",
                                      "AB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        p = MockAnalysisProject("CD",("CD3_S3_R1_001.fastq.gz",
                                      "CD3_S3_R2_001.fastq.gz",
                                      "CD4_S4_R1_001.fastq.gz",
                                      "CD4_S4_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        for p in ("AB","CD"):
            runqc.add_project(AnalysisProject(p,
                                              os.path.join(self.wd,p)))
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        for p in ("AB","CD"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (p,os.path.basename(self.wd)),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(self.wd,p,f)),
                                "Missing %s" % f)

    def test_run_qc_with_index_reads(self):
        """RunQC: standard QC run for project with index reads
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",
                                       "PJB2_S2_I1_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_run_qc_with_batching(self):
        """RunQC: standard QC run with batching
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           max_jobs=1,
                           poll_interval=0.5,
                           batch_size=3)
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_run_qc_with_batching_fails_for_missing_outputs(self):
        """RunQC: standard QC run with batching fails for missing outputs
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
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")))
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           max_jobs=1,
                           poll_interval=0.5,
                           batch_size=3)
        # Check output and reports
        self.assertEqual(status,1)
        self.assertTrue(os.path.exists(os.path.join(self.wd,"PJB","qc")),
                        "Missing 'qc'")
        for f in ("qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

    def test_run_qc_non_default_log_dir(self):
        """RunQC: standard QC run using non-default log dir
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
        # Non-default log dir
        log_dir = os.path.join(self.wd,"logs")
        self.assertFalse(os.path.exists(log_dir),
                         "Log dir '%s' already exists" % log_dir)
        # Set up and run the QC
        runqc = RunQC()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          log_dir=log_dir)
        status = runqc.run(multiqc=True,
                           qc_runner=SimpleJobRunner(),
                           verify_runner=SimpleJobRunner(),
                           report_runner=SimpleJobRunner(),
                           poll_interval=0.5,
                           max_jobs=1)
        # Check output and reports
        self.assertEqual(status,0)
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "PJB",
                                                   "qc")),
                         "'qc' directory doesn't exist, but should")
        for f in ("qc_report.html",
                  "qc_report.PJB.%s.zip" %
                  os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
        # Check log directory
        self.assertTrue(os.path.exists(log_dir),
                        "Log dir '%s' not found" % log_dir)
