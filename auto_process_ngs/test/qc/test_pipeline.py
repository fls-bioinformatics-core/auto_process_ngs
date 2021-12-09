#######################################################################
# Unit tests for qc/runqc.py
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.mock import MockFastqScreen
from auto_process_ngs.mock import MockFastQC
from auto_process_ngs.mock import MockFastqStrandPy
from auto_process_ngs.mock import MockMultiQC
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
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
        # Create a temp 'data' dir
        self.data = os.path.join(self.wd,"data")
        os.mkdir(self.data)
        # Add (empty) FastqScreen conf files
        self.fastq_screens = {
            'model_organisms': "fastq_screen_model_organisms.conf",
            'other_organisms': "fastq_screen_other_organisms.conf",
            'rRNA': "fastq_screen_rRNA.conf",
        }
        for screen in self.fastq_screens:
            conf_file = os.path.join(self.data,
                                     self.fastq_screens[screen])
            with open(conf_file,'wt') as fp:
                fp.write("")
            self.fastq_screens[screen] = conf_file
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
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_no_screens(self):
        """QCPipeline: standard QC run with no screens defined
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_strandedness(self):
        """QCPipeline: standard QC run with strandedness determination
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_missing_strandedness(self):
        """QCPipeline: standard QC fails with missing strandedness outputs
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"),
                                 no_outputs=True)
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,1)
        self.assertTrue(os.path.exists(os.path.join(self.wd,"PJB","qc")),
                        "Missing 'qc'")
        for f in ("qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

    def test_qcpipeline_no_multiqc(self):
        """QCPipeline: standard QC run (no MultiQC)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip"):
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
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"),
                                            no_outputs=True,
                                            exit_code=1)
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,1)
        self.assertTrue(os.path.exists(os.path.join(self.wd,"PJB","qc")),
                        "Missing 'qc'")
        for f in ("qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

    def test_qcpipeline_with_missing_fastqc_outputs(self):
        """QCPipeline: standard QC fails for missing FastQC outputs
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"),
                                       no_outputs=True,
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,1)
        self.assertTrue(os.path.exists(os.path.join(self.wd,"PJB","qc")),
                        "Missing 'qc'")
        for f in ("qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

    def test_qcpipeline_with_missing_multiqc_outputs(self):
        """QCPipeline: standard QC fails for missing MultiQC outputs
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,1)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip"):
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
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_non_default_output_dir(self):
        """QCPipeline: standard QC run using non-default output dir
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
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
                  "qc.non_default_report.PJB.zip",
                  "multiqc.non_default_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_single_end(self):
        """QCPipeline: standard QC run (single-end data)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_multiple_projects(self):
        """QCPipeline: standard QC run (multiple projects)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for p in ("AB","CD"):
            for f in ("qc",
                      "qc_report.html",
                      "qc_report.%s.zip" % p,
                      "multiqc_report.html"):
                self.assertTrue(os.path.exists(os.path.join(self.wd,p,f)),
                                "Missing %s" % f)

    def test_qcpipeline_with_index_reads(self):
        """QCPipeline: standard QC run for project with index reads
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_batching(self):
        """QCPipeline: standard QC run with batching
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           batch_size=3,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_batching_fails_for_missing_outputs(self):
        """QCPipeline: standard QC run with batching fails for missing outputs
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"),
                                       no_outputs=True,
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           batch_size=3,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,1)
        self.assertTrue(os.path.exists(os.path.join(self.wd,"PJB","qc")),
                        "Missing 'qc'")
        for f in ("qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

    def test_qcpipeline_non_default_log_dir(self):
        """QCPipeline: standard QC run using non-default log dir
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "PJB",
                                                   "qc")),
                         "'qc' directory doesn't exist, but should")
        for f in ("qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
        # Check log directory
        self.assertTrue(os.path.exists(log_dir),
                        "Log dir '%s' not found" % log_dir)

    def test_qcpipeline_with_cellranger_count_scRNA_seq_310(self):
        """QCPipeline: single cell RNA-seq QC run with 'cellranger count' (v3.1.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="3.1.0",
                                 assert_include_introns=False)
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-cellranger-GRCh38-1.2.0' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB1/_cmdline",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB2/_cmdline",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB1/_cmdline",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB1/outs/web_summary.html",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB2/_cmdline",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB2/outs/web_summary.html",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_scRNA_seq_501(self):
        """QCPipeline: single cell RNA-seq QC run with 'cellranger count' (v5.0.1)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="5.0.1",
                                 assert_include_introns=False)
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_scRNA_seq_600(self):
        """QCPipeline: single cell RNA-seq QC run with 'cellranger count' (v6.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="6.0.0",
                                 assert_include_introns=False)
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_specify_exe(self):
        """QCPipeline: single cell RNA-seq QC run with 'cellranger count' (specify cellranger exe)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock cellranger not on PATH
        cellranger_bin = os.path.join(self.wd,"cellranger")
        os.mkdir(cellranger_bin)
        MockCellrangerExe.create(os.path.join(cellranger_bin,
                                              "cellranger"),
                                 version="5.0.1")
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           cellranger_exe=os.path.join(cellranger_bin,
                                                       "cellranger"),
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_snRNA_seq_310(self):
        """QCPipeline: single nuclei RNA-seq QC run with 'cellranger count' (v3.1.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="3.1.0",
                                 assert_include_introns=False)
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_premrna_references=
                           { 'human': '/data/refdata-cellranger-GRCh38-1.2.0_premrna' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB1/_cmdline",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB2/_cmdline",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB1/_cmdline",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB1/outs/web_summary.html",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB2/_cmdline",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB2/outs/web_summary.html",
                  "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_snRNA_seq_501(self):
        """QCPipeline: single nuclei RNA-seq QC run with 'cellranger count' (v5.0.1)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        # Mock cellranger 5.0.1 with check on --include-introns
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="5.0.1",
                                 assert_include_introns=True)
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_snRNA_seq_600(self):
        """QCPipeline: single nuclei RNA-seq QC run with 'cellranger count' (v6.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        # Mock cellranger 5.0.1 with check on --include-introns
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="6.0.0",
                                 assert_include_introns=True)
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_atac_count_1_2_0(self):
        """QCPipeline: single cell ATAC QC run with 'cellranger-atac count' (1.2.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"),
                                 version="1.2.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_atac_references=
                           { 'human':
                             '/data/refdata-cellranger-atac-GRCh38-1.2.0' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/_cmdline",
                  "qc/cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/outs/summary.csv",
                  "qc/cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/_cmdline",
                  "qc/cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/outs/summary.csv",
                  "cellranger_count",
                  "cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/_cmdline",
                  "cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/outs/web_summary.html",
                  "cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/outs/summary.csv",
                  "cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/_cmdline",
                  "cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/outs/web_summary.html",
                  "cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_atac_count_2_0_0(self):
        """QCPipeline: single cell ATAC QC run with 'cellranger-atac count' (2.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"),
                                 version="2.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_atac_references=
                           { 'human':
                             '/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv",
                  "cellranger_count",
                  "cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                  "cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                  "cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                  "cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                  "cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                  "cellranger_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_multiome_atac(self):
        """QCPipeline: single cell multiome ATAC QC run
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock multiome ATAC analysis project
        p = MockAnalysisProject("PJB_ATAC",("PJB1_ATAC_S1_R1_001.fastq.gz",
                                            "PJB1_ATAC_S1_R2_001.fastq.gz",
                                            "PJB1_ATAC_S1_R3_001.fastq.gz",
                                            "PJB2_ATAC_S2_R1_001.fastq.gz",
                                            "PJB2_ATAC_S2_R2_001.fastq.gz",
                                            "PJB2_ATAC_S2_R3_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell Multiome',
                                           'Library type': 'ATAC' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_ATAC")),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB_ATAC.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_ATAC",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_multiome_gex(self):
        """QCPipeline: single cell multiome GEX QC run
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock multiome GEX analysis project
        p = MockAnalysisProject("PJB_GEX",("PJB1_GEX_S1_R1_001.fastq.gz",
                                           "PJB1_GEX_S1_R2_001.fastq.gz",
                                           "PJB2_GEX_S2_R1_001.fastq.gz",
                                           "PJB2_GEX_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell Multiome',
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_GEX")),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB_GEX.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_GEX",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_multiome_atac_with_cellranger_arc_count_1_0_0(self):
        """QCPipeline: single cell multiome ATAC QC run with 'cellranger-arc count' (Cellranger ARC 1.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="1.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock multiome ATAC analysis project
        p = MockAnalysisProject("PJB_ATAC",("PJB1_ATAC_S1_R1_001.fastq.gz",
                                            "PJB1_ATAC_S1_R2_001.fastq.gz",
                                            "PJB1_ATAC_S1_R3_001.fastq.gz",
                                            "PJB2_ATAC_S2_R1_001.fastq.gz",
                                            "PJB2_ATAC_S2_R2_001.fastq.gz",
                                            "PJB2_ATAC_S2_R3_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell Multiome',
                                           'Library type': 'ATAC' })
        p.create(top_dir=self.wd)
        # Make mock multiome GEX analysis project (with QC outputs)
        p2 = MockAnalysisProject("PJB_GEX",("PJB1_GEX_S1_R1_001.fastq.gz",
                                            "PJB1_GEX_S1_R2_001.fastq.gz",
                                            "PJB2_GEX_S2_R1_001.fastq.gz",
                                            "PJB2_GEX_S2_R2_001.fastq.gz",),
                                 metadata={ 'Organism': 'Human',
                                            'Single cell platform':
                                            '10xGenomics Single Cell Multiome',
                                            'Library type': 'GEX' })
        p2.create(top_dir=self.wd)
        UpdateAnalysisProject(
            AnalysisProject(os.path.join(self.wd,p2.name))).\
            add_qc_outputs(protocol='10x_Multiome_GEX')
        # Add the 10x_multiome_libraries.info file
        with open(os.path.join(self.wd,
                               "PJB_ATAC",
                               "10x_multiome_libraries.info"),'wt') as fp:
            fp.write("{sample1}\t{working_dir}:{project}/{sample2}\n".format(
                sample1="PJB1_ATAC",
                sample2="PJB1_GEX",
                working_dir=self.wd,
                project="PJB_GEX"))
            fp.write("{sample1}\t{working_dir}:{project}/{sample2}\n".format(
                sample1="PJB2_ATAC",
                sample2="PJB2_GEX",
                working_dir=self.wd,
                project="PJB_GEX"))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_ATAC")),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB_ATAC.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",
                  "cellranger_count",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_ATAC",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_multiome_atac_with_cellranger_arc_count_2_0_0(self):
        """QCPipeline: single cell multiome ATAC QC run with 'cellranger-arc count' (Cellranger ARC 2.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="2.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock multiome ATAC analysis project
        p = MockAnalysisProject("PJB_ATAC",("PJB1_ATAC_S1_R1_001.fastq.gz",
                                            "PJB1_ATAC_S1_R2_001.fastq.gz",
                                            "PJB1_ATAC_S1_R3_001.fastq.gz",
                                            "PJB2_ATAC_S2_R1_001.fastq.gz",
                                            "PJB2_ATAC_S2_R2_001.fastq.gz",
                                            "PJB2_ATAC_S2_R3_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell Multiome',
                                           'Library type': 'ATAC' })
        p.create(top_dir=self.wd)
        # Make mock multiome GEX analysis project (with QC outputs)
        p2 = MockAnalysisProject("PJB_GEX",("PJB1_GEX_S1_R1_001.fastq.gz",
                                            "PJB1_GEX_S1_R2_001.fastq.gz",
                                            "PJB2_GEX_S2_R1_001.fastq.gz",
                                            "PJB2_GEX_S2_R2_001.fastq.gz",),
                                 metadata={ 'Organism': 'Human',
                                            'Single cell platform':
                                            '10xGenomics Single Cell Multiome',
                                            'Library type': 'GEX' })
        p2.create(top_dir=self.wd)
        UpdateAnalysisProject(
            AnalysisProject(os.path.join(self.wd,p2.name))).\
            add_qc_outputs(protocol='10x_Multiome_GEX')
        # Add the 10x_multiome_libraries.info file
        with open(os.path.join(self.wd,
                               "PJB_ATAC",
                               "10x_multiome_libraries.info"),'wt') as fp:
            fp.write("{sample1}\t{working_dir}:{project}/{sample2}\n".format(
                sample1="PJB1_ATAC",
                sample2="PJB1_GEX",
                working_dir=self.wd,
                project="PJB_GEX"))
            fp.write("{sample1}\t{working_dir}:{project}/{sample2}\n".format(
                sample1="PJB2_ATAC",
                sample2="PJB2_GEX",
                working_dir=self.wd,
                project="PJB_GEX"))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_ATAC")),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB_ATAC.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",
                  "cellranger_count",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_ATAC",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_multiome_gex_with_cellranger_arc_count_1_0_0(self):
        """QCPipeline: single cell multiome GEX QC run with 'cellranger-arc count' (Cellranger ARC 1.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="1.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock multiome GEX analysis project
        p = MockAnalysisProject("PJB_GEX",("PJB1_GEX_S1_R1_001.fastq.gz",
                                           "PJB1_GEX_S1_R2_001.fastq.gz",
                                           "PJB2_GEX_S2_R1_001.fastq.gz",
                                           "PJB2_GEX_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell Multiome',
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        # Make mock multiome ATAC analysis project (with QC outputs)
        p2 = MockAnalysisProject("PJB_ATAC",("PJB1_ATAC_S1_R1_001.fastq.gz",
                                             "PJB1_ATAC_S1_R2_001.fastq.gz",
                                             "PJB1_ATAC_S1_R3_001.fastq.gz",
                                             "PJB2_ATAC_S2_R1_001.fastq.gz",
                                             "PJB2_ATAC_S2_R2_001.fastq.gz",
                                             "PJB2_ATAC_S2_R3_001.fastq.gz"),
                                 metadata={ 'Organism': 'Human',
                                            'Single cell platform':
                                            '10xGenomics Single Cell Multiome',
                                            'Library type': 'ATAC' })
        p2.create(top_dir=self.wd)
        UpdateAnalysisProject(
            AnalysisProject(os.path.join(self.wd,p2.name))).\
            add_qc_outputs(protocol='10x_Multiome_ATAC')
        # Add the 10x_multiome_libraries.info file
        with open(os.path.join(self.wd,
                               "PJB_GEX",
                               "10x_multiome_libraries.info"),'wt') as fp:
            fp.write("{sample1}\t{working_dir}:{project}/{sample2}\n".format(
                sample1="PJB1_GEX",
                sample2="PJB1_ATAC",
                working_dir=self.wd,
                project="PJB_ATAC"))
            fp.write("{sample1}\t{working_dir}:{project}/{sample2}\n".format(
                sample1="PJB2_GEX",
                sample2="PJB2_ATAC",
                working_dir=self.wd,
                project="PJB_ATAC"))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_GEX")),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB_GEX.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "qc/cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv",
                  "cellranger_count",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_GEX",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_multiome_gex_with_cellranger_arc_count_2_0_0(self):
        """QCPipeline: single cell multiome GEX QC run with 'cellranger-arc count' (Cellranger ARC 2.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="2.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock multiome GEX analysis project
        p = MockAnalysisProject("PJB_GEX",("PJB1_GEX_S1_R1_001.fastq.gz",
                                           "PJB1_GEX_S1_R2_001.fastq.gz",
                                           "PJB2_GEX_S2_R1_001.fastq.gz",
                                           "PJB2_GEX_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Single Cell Multiome',
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        # Make mock multiome ATAC analysis project (with QC outputs)
        p2 = MockAnalysisProject("PJB_ATAC",("PJB1_ATAC_S1_R1_001.fastq.gz",
                                             "PJB1_ATAC_S1_R2_001.fastq.gz",
                                             "PJB1_ATAC_S1_R3_001.fastq.gz",
                                             "PJB2_ATAC_S2_R1_001.fastq.gz",
                                             "PJB2_ATAC_S2_R2_001.fastq.gz",
                                             "PJB2_ATAC_S2_R3_001.fastq.gz"),
                                 metadata={ 'Organism': 'Human',
                                            'Single cell platform':
                                            '10xGenomics Single Cell Multiome',
                                            'Library type': 'ATAC' })
        p2.create(top_dir=self.wd)
        UpdateAnalysisProject(
            AnalysisProject(os.path.join(self.wd,p2.name))).\
            add_qc_outputs(protocol='10x_Multiome_ATAC')
        # Add the 10x_multiome_libraries.info file
        with open(os.path.join(self.wd,
                               "PJB_GEX",
                               "10x_multiome_libraries.info"),'wt') as fp:
            fp.write("{sample1}\t{working_dir}:{project}/{sample2}\n".format(
                sample1="PJB1_GEX",
                sample2="PJB1_ATAC",
                working_dir=self.wd,
                project="PJB_ATAC"))
            fp.write("{sample1}\t{working_dir}:{project}/{sample2}\n".format(
                sample1="PJB2_GEX",
                sample2="PJB2_ATAC",
                working_dir=self.wd,
                project="PJB_ATAC"))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_GEX")),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB_GEX.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "qc/cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv",
                  "cellranger_count",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "cellranger_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_GEX",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_cellplex_with_cellranger_multi(self):
        """QCPipeline: 10xGenomics Cellplex run with 'cellranger multi'
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="6.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x Cellplex analysis project
        p = MockAnalysisProject("PJB",("PJB1_GEX_S1_R1_001.fastq.gz",
                                       "PJB1_GEX_S1_R2_001.fastq.gz",
                                       "PJB2_MC_S2_R1_001.fastq.gz",
                                       "PJB2_MC_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 3\'v3',
                                           'Library type': 'CellPlex' })
        p.create(top_dir=self.wd)
        # Add the cellranger multi config.csv file
        with open(os.path.join(self.wd,
                               "PJB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(self.wd,
                                     "PJB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,%s,any,PJB1,gene expression,
PJB2_MC,%s,any,PJB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBA,CMO301,PBA
PBB,CMO302,PBB
""" % (fastq_dir,fastq_dir))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-gex-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_multi",
                  "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/_cmdline",
                  "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBA/web_summary.html",
                  "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                  "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBB/web_summary.html",
                  "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                  "qc/cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv",
                  "cellranger_multi",
                  "cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/_cmdline",
                  "cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBA/web_summary.html",
                  "cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                  "cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBB/web_summary.html",
                  "cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                  "cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_10x_scRNAseq_no_project_metadata(self):
        """QCPipeline: single cell RNA-seq QC run with no project metadata
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version='5.0.1')
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={})
        p.create(top_dir=self.wd)
        # Remove the README.info file
        os.remove(os.path.join(self.wd,"PJB","README.info"))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          multiqc=True,
                          qc_protocol="10x_scRNAseq",
                          organism="human")
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                  "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_no_references(self):
        """QCPipeline: single cell QC run with 'cellranger count' (no reference data)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="5.0.1")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_rerun_with_protocol_mismatch(self):
        """QCPipeline: handle QC protocol mismatch when rerunning pipeline
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Add existing QC outputs
        UpdateAnalysisProject(
            AnalysisProject("PJB",os.path.join(self.wd,"PJB"))).add_qc_outputs(
                protocol="standardSE",
                include_fastq_strand=False,
                include_multiqc=True)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject("PJB",
                                          os.path.join(self.wd,"PJB")),
                          qc_protocol="standardPE",
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=0.5,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        # Check output and reports
        self.assertEqual(status,0)
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
