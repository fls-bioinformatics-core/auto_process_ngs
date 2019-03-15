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
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.pipeline import QCPipeline
from auto_process_ngs.qc.pipeline import copy_project
from auto_process_ngs.qc.pipeline import check_illumina_qc_outputs
from auto_process_ngs.qc.pipeline import check_fastq_strand_outputs

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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          fastq_strand_indexes={ 'human': '/data/hg38/star_index' },
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          fastq_strand_indexes={ 'human': '/data/hg38/star_index' },
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=False)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          fastq_dir="fastqs.cells",
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          qc_dir="qc.non_default",
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        for p in ("AB","CD"):
            runqc.add_project(AnalysisProject(p,
                                              os.path.join(self.wd,p)),
                              multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(poll_interval=0.5,
                           max_jobs=1,
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
        runqc = QCPipeline(runners={ 'default': SimpleJobRunner(), })
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True,
                          log_dir=log_dir)
        status = runqc.run(poll_interval=0.5,
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

class TestCopyProject(unittest.TestCase):
    """
    Tests for the 'copy_project' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCopyProject')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_copy_project(self):
        """
        copy_project: copies project instance
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Make initial project
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Make a copy
        project2 = copy_project(project)
        # Check copy
        self.assertEqual(project.name,project2.name)
        self.assertEqual(project.dirn,project2.dirn)
        self.assertEqual(project.fastq_dir,project2.fastq_dir)
        self.assertEqual(project.fastq_dirs,project2.fastq_dirs)
        self.assertEqual(project.fastqs,project2.fastqs)
        self.assertEqual(project.info.organism,project2.info.organism)

class TestCheckIlluminaQcOutputs(unittest.TestCase):
    """
    Tests for the 'check_illumina_qc_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestIlluminaQcOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_illumina_qc_outputs_standardPE_all_missing(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs missing (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardPE"),
                         project.fastqs)

    def test_check_illumina_qc_outputs_standardPE_all_present(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs present (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardPE"),
                         [])

    def test_check_illumina_qc_outputs_standardPE_some_missing(self):
        """
        check_illumina_qc_outputs: some illumina_qc.sh outputs missing (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove some outputs
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R1_001_model_organisms_screen.txt",):
            os.remove(os.path.join(project.qc_dir,f))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardPE"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_illumina_qc_outputs_standardSE_all_missing(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs missing (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardSE"),
                         project.fastqs)

    def test_check_illumina_qc_outputs_standardSE_all_present(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs present (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardSE"),
                         [])

    def test_check_illumina_qc_outputs_standardSE_some_missing(self):
        """
        check_illumina_qc_outputs: some illumina_qc.sh outputs missing (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove some outputs
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R1_001_model_organisms_screen.txt",):
            os.remove(os.path.join(project.qc_dir,f))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="standardSE"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_illumina_qc_outputs_singlecell_all_missing(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="singlecell"),
                         project.fastqs)

    def test_check_illumina_qc_outputs_singlecell_all_present(self):
        """
        check_illumina_qc_outputs: all illumina_qc.sh outputs present (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="singlecell"),
                         [])

    def test_check_illumina_qc_outputs_singlecell_some_missing(self):
        """
        check_illumina_qc_outputs: some illumina_qc.sh outputs missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # Remove some outputs
        for f in ("PJB1_S1_R2_001_fastqc.html",
                  "PJB1_S1_R2_001_model_organisms_screen.txt",):
            os.remove(os.path.join(project.qc_dir,f))
        # Check
        self.assertEqual(check_illumina_qc_outputs(project,
                                                   qc_dir="qc",
                                                   qc_protocol="singlecell"),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R2_001.fastq.gz")])

class TestCheckFastqStrandOutputs(unittest.TestCase):
    """
    Tests for the 'fastq_strand_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestFastqStrandOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_fastq_strand_outputs_standardPE_missing(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Make fastq_strand.conf
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        with open(fastq_strand_conf,'w') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="standardPE"),
                         [(os.path.join(project.fastq_dir,
                                        "PJB1_S1_R1_001.fastq.gz"),
                           os.path.join(project.fastq_dir,
                                        "PJB1_S1_R2_001.fastq.gz")),])

    def test_check_fastq_strand_outputs_standardPE_present(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (standardPE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="standardPE",
            include_fastq_strand=True,
            include_multiqc=False)
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="standardPE"),
                         [])

    def test_check_fastq_strand_outputs_standardSE_missing(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Make fastq_strand.conf
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        with open(fastq_strand_conf,'w') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="standardSE"),
                         [(os.path.join(project.fastq_dir,
                                        "PJB1_S1_R1_001.fastq.gz"),),])

    def test_check_fastq_strand_outputs_standardSE_present(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (standardSE)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="standardSE",
            include_fastq_strand=True,
            include_multiqc=False)
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="standardSE"),
                         [])

    def test_check_fastq_strand_outputs_singlecell_missing(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Make fastq_strand.conf
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        with open(fastq_strand_conf,'w') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="singlecell"),
                         [(os.path.join(project.fastq_dir,
                                        "PJB1_S1_R2_001.fastq.gz"),),])

    def test_check_fastq_strand_outputs_singlecell_present(self):
        """
        check_fastq_strand_outputs: fastq_strand.py output present (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="singlecell",
            include_fastq_strand=True,
            include_multiqc=False)
        fastq_strand_conf = os.path.join(project.dirn,"fastq_strand.conf")
        # Check the outputs
        self.assertEqual(check_fastq_strand_outputs(project,
                                                    "qc",
                                                    fastq_strand_conf,
                                                    qc_protocol="singlecell"),
                         [])
