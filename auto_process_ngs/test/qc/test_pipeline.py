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
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.pipeline import QCPipeline

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestQCPipeline(unittest.TestCase):
    """
    Tests for QCPipeline class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestQCPipeline')
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

    def test_qcpipeline(self):
        """QCPipeline: standard QC run
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_strandedness(self):
        """QCPipeline: standard QC run with strandedness determination
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
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(fastq_strand_indexes=
                           { 'human': '/data/hg38/star_index' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_missing_strandedness(self):
        """QCPipeline: standard QC fails with missing strandedness outputs
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
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(fastq_strand_indexes=
                           { 'human': '/data/hg38/star_index' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
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

    def test_qcpipeline_no_multiqc(self):
        """QCPipeline: standard QC run (no MultiQC)
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=False)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
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

    def test_qcpipelne_with_missing_fastq_screen_outputs(self):
        """QCPipeline: standard QC fails for missing FastQScreen outputs
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
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

    def test_qcpipeline_with_missing_fastqc_outputs(self):
        """QCPipeline: standard QC fails for missing FastQC outputs
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
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

    def test_qcpipeline_with_missing_multiqc_outputs(self):
        """QCPipeline: standard QC fails for missing MultiQC outputs
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
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

    def test_qcpipeline_non_default_fastq_dir(self):
        """QCPipeline: standard QC run using non-default Fastq dir
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          fastq_dir="fastqs.cells",
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_non_default_output_dir(self):
        """QCPipeline: standard QC run using non-default output dir
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          qc_dir="qc.non_default",
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
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

    def test_qcpipeline_single_end(self):
        """QCPipeline: standard QC run (single-end data)
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_multiple_projects(self):
        """QCPipeline: standard QC run (multiple projects)
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
        runqc = QCPipeline()
        for p in ("AB","CD"):
            runqc.add_project(AnalysisProject(p,
                                              os.path.join(self.wd,p)),
                              multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for p in ("AB","CD"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.%s.zip" % (p,os.path.basename(self.wd)),
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(self.wd,p,f)),
                                "Missing %s" % f)

    def test_qcpipeline_with_index_reads(self):
        """QCPipeline: standard QC run for project with index reads
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_batching(self):
        """QCPipeline: standard QC run with batching
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           batch_size=3,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_batching_fails_for_missing_outputs(self):
        """QCPipeline: standard QC run with batching fails for missing outputs
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           batch_size=3,
                           runners={ 'default': SimpleJobRunner(), })
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

    def test_qcpipeline_non_default_log_dir(self):
        """QCPipeline: standard QC run using non-default log dir
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
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True,
                          log_dir=log_dir)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
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

    def test_qcpipeline_with_cellranger_count_scRNA_seq(self):
        """QCPipeline: single cell RNA-seq QC run with 'cellranger count'
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 3\'v2',
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(fastq_strand_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_transcriptomes=
                           { 'human': '/data/hg38/cellranger' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "qc/cellranger_count",
                  "qc/cellranger_count/PJB1/_cmdline",
                  "qc/cellranger_count/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/PJB2/_cmdline",
                  "qc/cellranger_count/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/PJB1/_cmdline",
                  "cellranger_count/PJB1/outs/web_summary.html",
                  "cellranger_count/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/PJB2/_cmdline",
                  "cellranger_count/PJB2/outs/web_summary.html",
                  "cellranger_count/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_snRNA_seq(self):
        """QCPipeline: single nuclei RNA-seq QC run with 'cellranger count'
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 3\'v2',
                                           'Library type': 'snRNA-seq' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(fastq_strand_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_premrna_references=
                           { 'human': '/data/hg38/cellranger' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "qc/cellranger_count",
                  "qc/cellranger_count/PJB1/_cmdline",
                  "qc/cellranger_count/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/PJB2/_cmdline",
                  "qc/cellranger_count/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/PJB1/_cmdline",
                  "cellranger_count/PJB1/outs/web_summary.html",
                  "cellranger_count/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/PJB2/_cmdline",
                  "cellranger_count/PJB2/outs/web_summary.html",
                  "cellranger_count/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_atac_count(self):
        """QCPipeline: single cell ATAC QC run with 'cellranger-atac count'
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB1_S1_R3_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",
                                       "PJB2_S2_R3_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell ATAC',
                                           'Library type': 'scATAC-seq' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(fastq_strand_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_atac_references=
                           { 'human': '/data/hg38/cellranger' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "qc/cellranger_count",
                  "qc/cellranger_count/PJB1/_cmdline",
                  "qc/cellranger_count/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/PJB1/outs/summary.csv",
                  "qc/cellranger_count/PJB2/_cmdline",
                  "qc/cellranger_count/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/PJB2/outs/summary.csv",
                  "cellranger_count",
                  "cellranger_count/PJB1/_cmdline",
                  "cellranger_count/PJB1/outs/web_summary.html",
                  "cellranger_count/PJB1/outs/summary.csv",
                  "cellranger_count/PJB2/_cmdline",
                  "cellranger_count/PJB2/outs/web_summary.html",
                  "cellranger_count/PJB2/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_no_references(self):
        """QCPipeline: single cell QC run with 'cellranger count' (no reference data)
        """
        # Make mock illumina_qc.sh and multiqc
        MockIlluminaQcSh.create(os.path.join(self.bin,
                                             "illumina_qc.sh"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 3\'v2',
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(fastq_strand_indexes=
                           { 'human': '/data/hg38/star_index' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.%s.zip" % os.path.basename(self.wd),
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
