#######################################################################
# Unit tests for qc/modules/fastqc.py ('fastqc' QC module)
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.mock import MockFastQC
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.protocols import QCProtocol
from auto_process_ngs.qc.pipeline import QCPipeline
from auto_process_ngs.qc.modules.fastqc import check_fastqc_outputs
from ..protocols import BaseQCPipelineTestCase

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Polling interval for pipeline
POLL_INTERVAL = 0.1

class TestQCPipelineFastqc(BaseQCPipelineTestCase):
    """
    Tests for 'fastqc' QC module
    """
    def test_qcpipeline_qc_modules_fastqc_pe(self):
        """
        QCPipeline: 'fastqc' QC module (PE data)
        """
        # Make mock QC executables
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastqc",
                              description="Fastqc test",
                              seq_data_reads=['r1','r2',],
                              index_reads=None,
                              qc_modules=("fastqc",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R1_001_fastqc.zip",
                  "PJB1_S1_R1_001_fastqc/",
                  "PJB1_S1_R2_001_fastqc.html",
                  "PJB1_S1_R2_001_fastqc.zip",
                  "PJB1_S1_R2_001_fastqc/",
                  "PJB2_S2_R1_001_fastqc.html",
                  "PJB2_S2_R1_001_fastqc.zip",
                  "PJB2_S2_R1_001_fastqc/",
                  "PJB2_S2_R2_001_fastqc.html",
                  "PJB2_S2_R2_001_fastqc.zip",
                  "PJB2_S2_R2_001_fastqc/"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastqc")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_fastqc_pe_with_r1_index_reads(self):
        """
        QCPipeline: 'fastqc' QC module (PE data, R1 as index reads)
        """
        # Make mock QC executables
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastqc",
                              description="Fastqc test",
                              seq_data_reads=['r2',],
                              index_reads=['r1',],
                              qc_modules=("fastqc",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R1_001_fastqc.zip",
                  "PJB1_S1_R1_001_fastqc/",
                  "PJB1_S1_R2_001_fastqc.html",
                  "PJB1_S1_R2_001_fastqc.zip",
                  "PJB1_S1_R2_001_fastqc/",
                  "PJB2_S2_R1_001_fastqc.html",
                  "PJB2_S2_R1_001_fastqc.zip",
                  "PJB2_S2_R1_001_fastqc/",
                  "PJB2_S2_R2_001_fastqc.html",
                  "PJB2_S2_R2_001_fastqc.zip",
                  "PJB2_S2_R2_001_fastqc/"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastqc")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_fastqc_se(self):
        """
        QCPipeline: 'fastqc' QC module (SE data)
        """
        # Make mock QC executables
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastqc",
                              description="Fastqc test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastqc",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R1_001_fastqc.zip",
                  "PJB1_S1_R1_001_fastqc/",
                  "PJB2_S2_R1_001_fastqc.html",
                  "PJB2_S2_R1_001_fastqc.zip",
                  "PJB2_S2_R1_001_fastqc/"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastqc")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_fastqc_se_with_biological_samples(self):
        """
        QCPipeline: 'fastqc' QC module (SE data with biological samples)
        """
        # Make mock QC executables
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Biological samples': 'PJB1' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastqc",
                              description="Fastqc test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastqc",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R1_001_fastqc.zip",
                  "PJB1_S1_R1_001_fastqc/",
                  "PJB2_S2_R1_001_fastqc.html",
                  "PJB2_S2_R1_001_fastqc.zip",
                  "PJB2_S2_R1_001_fastqc/"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastqc")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_fastqc_se_ignore_index_reads(self):
        """
        QCPipeline: 'fastqc' QC module (SE data, ignore I* index reads)
        """
        # Make mock QC executables
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_I1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastqc",
                              description="Fastqc test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastqc",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R1_001_fastqc.zip",
                  "PJB1_S1_R1_001_fastqc/",
                  "PJB2_S2_R1_001_fastqc.html",
                  "PJB2_S2_R1_001_fastqc.zip",
                  "PJB2_S2_R1_001_fastqc/"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        for f in ("PJB1_S1_I1_001_fastqc.html",
                  "PJB1_S1_I1_001_fastqc.zip",
                  "PJB1_S1_I1_001_fastqc/",
                  "PJB2_S2_I1_001_fastqc.html",
                  "PJB2_S2_I1_001_fastqc.zip",
                  "PJB2_S2_I1_001_fastqc/"):
            self.assertFalse(os.path.exists(os.path.join(qc_dir,f)),
                             "%s: present" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastqc")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_I1_001.fastq.gz,"
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_I1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_fastqc_se_split_lanes(self):
        """
        QCPipeline: 'fastqc' QC module (SE data, split by lane)
        """
        # Make mock QC executables
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastqc",
                              description="Fastqc test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastqc",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol,
                          split_fastqs_by_lane=True)
        status = runqc.run(poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_L001_R1_001_fastqc.html",
                  "PJB1_S1_L001_R1_001_fastqc.zip",
                  "PJB1_S1_L001_R1_001_fastqc/",
                  "PJB2_S2_L001_R1_001_fastqc.html",
                  "PJB2_S2_L001_R1_001_fastqc.zip",
                  "PJB2_S2_L001_R1_001_fastqc/"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastqc")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_L001_R1_001.fastq.gz,"
                         "PJB2_S2_L001_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,True)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_fastqc_se_missing_output(self):
        """
        QCPipeline: 'fastqc' QC module (SE data, missing outputs)
        """
        # Make mock QC executables
        MockFastQC.create(os.path.join(self.bin,"fastqc"),
                          no_outputs=True,
                          exit_code=1)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastqc",
                              description="Fastqc test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastqc",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,1)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastqc")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        self.assertTrue(os.path.exists(os.path.join(self.wd,"PJB","qc")),
                        "Missing 'qc'")
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)

class TestCheckFastQCOutputs(unittest.TestCase):
    """
    Tests for the 'check_fastqc_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestFastQCOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_fastqc_outputs_paired_end_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (paired end)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              read_numbers=(1,2)),
                         project.fastqs)

    def test_check_fastqc_outputs_paired_end_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (paired end)
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
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              read_numbers=(1,2)),
                         [])

    def test_check_fastqc_outputs_paired_end_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (paired end)
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
        # Remove a FastQC output
        os.remove(os.path.join(project.qc_dir,
                               "PJB1_S1_R1_001_fastqc.html"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              read_numbers=(1,2)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastqc_outputs_single_end_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (single end)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              read_numbers=(1,)),
                         project.fastqs)

    def test_check_fastqc_outputs_single_end_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (single end)
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
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              read_numbers=(1,)),
                         [])

    def test_check_fastqc_outputs_single_end_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (single end)
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
        # Remove a FastQC output
        os.remove(os.path.join(project.qc_dir,
                               "PJB1_S1_R1_001_fastqc.html"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              read_numbers=(1,)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastqc_outputs_single_cell_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (single cell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              read_numbers=(1,2,)),
                         project.fastqs)

    def test_check_fastqc_outputs_single_cell_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (single cell)
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
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              read_numbers=(1,2,)),
                         [])

    def test_check_fastqc_outputs_single_cell_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (single cell)
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
        # Remove a FastQC output
        os.remove(os.path.join(project.qc_dir,
                               "PJB1_S1_R1_001_fastqc.html"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              read_numbers=(1,2,)),
                         [os.path.join(project.fastq_dir,
                                       "PJB1_S1_R1_001.fastq.gz")])

    def test_check_fastqc_outputs_from_fastq_list_all_missing(self):
        """
        check_fastqc_outputs: all FastQC outputs missing (Fastq list)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # List of Fastqs (subset of all Fastqs)
        fastqs = ["PJB2_S2_R1_001.fastq.gz",
                  "PJB2_S2_R2_001.fastq.gz"]
        # Get the outputs
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              fastqs=fastqs,
                                              read_numbers=(1,2,)),
                         fastqs)

    def test_check_fastqc_outputs_from_fastq_list_all_present(self):
        """
        check_fastqc_outputs: all FastQC outputs present (Fastq list)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # List of Fastqs (subset of all Fastqs)
        fastqs = ["PJB2_S2_R1_001.fastq.gz",
                  "PJB2_S2_R2_001.fastq.gz"]
        # Remove FastQC outputs for Fastqs not in the list
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R2_001_fastqc.html"):
            os.remove(os.path.join(project.qc_dir,f))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              fastqs=fastqs,
                                              read_numbers=(1,2,)),
                         [])

    def test_check_fastqc_outputs_from_fastq_list_some_missing(self):
        """
        check_fastqc_outputs: some FastQC outputs missing (Fastq list)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Add QC artefacts
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            include_fastq_strand=False,
            include_multiqc=False)
        # List of Fastqs (subset of all Fastqs)
        fastqs = ["PJB2_S2_R1_001.fastq.gz",
                  "PJB2_S2_R2_001.fastq.gz"]
        # Remove FastQC outputs for Fastqs not in the list
        for f in ("PJB1_S1_R1_001_fastqc.html",
                  "PJB1_S1_R2_001_fastqc.html"):
            os.remove(os.path.join(project.qc_dir,f))
        # Remove a FastQC output from the list
        os.remove(os.path.join(project.qc_dir,
                               "PJB2_S2_R1_001_fastqc.html"))
        # Check
        self.assertEqual(check_fastqc_outputs(project,
                                              qc_dir="qc",
                                              fastqs=fastqs,
                                              read_numbers=(1,2,)),
                         ["PJB2_S2_R1_001.fastq.gz"])
