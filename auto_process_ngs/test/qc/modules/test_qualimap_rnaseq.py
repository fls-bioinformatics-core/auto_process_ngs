#######################################################################
# Unit tests for qc/pipeline.py ('qualimap_rnaseq' QC module)
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.mock import MockGtf2bed
from auto_process_ngs.mock import MockQualimap
from auto_process_ngs.mock import MockRSeQC
from auto_process_ngs.mock import MockSamtools
from auto_process_ngs.mock import MockStar
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.protocols import QCProtocol
from auto_process_ngs.qc.pipeline import QCPipeline
from ..protocols import BaseQCPipelineTestCase

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Polling interval for pipeline
POLL_INTERVAL = 0.1

class TestQCPipelineQualimapRnaseq(BaseQCPipelineTestCase):
    """
    Tests for 'qualimap_rnaseq' QC module
    """
    def test_qcpipeline_qc_modules_qualimap_rnaseq_pe(self):
        """
        QCPipeline: 'qualimap_rnaseq' QC module (PE data)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        protocol = QCProtocol(name="qualimap_rnaseq",
                              description="Qualimap_rnaseq test",
                              seq_data_reads=['r1','r2'],
                              index_reads=None,
                              qc_modules=("qualimap_rnaseq",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        qualimap_dir = os.path.join(qc_dir,"qualimap-rnaseq","human")
        for f in ("PJB1_S1_001/qualimapReport.html",
                  "PJB1_S1_001/rnaseq_qc_results.txt",
                  "PJB2_S2_001/qualimapReport.html",
                  "PJB2_S2_001/rnaseq_qc_results.txt"):
            self.assertTrue(os.path.exists(os.path.join(qualimap_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"qualimap_rnaseq")
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
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_qualimap_rnaseq_pe_with_r1_index_reads(self):
        """
        QCPipeline: 'qualimap_rnaseq' QC module (PE data, R1 as index reads)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        protocol = QCProtocol(name="qualimap_rnaseq",
                              description="Qualimap_rnaseq test",
                              seq_data_reads=['r2'],
                              index_reads=['r1'],
                              qc_modules=("qualimap_rnaseq",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        qualimap_dir = os.path.join(qc_dir,"qualimap-rnaseq","human")
        for f in ("PJB1_S1_001/qualimapReport.html",
                  "PJB1_S1_001/rnaseq_qc_results.txt",
                  "PJB2_S2_001/qualimapReport.html",
                  "PJB2_S2_001/rnaseq_qc_results.txt"):
            self.assertTrue(os.path.exists(os.path.join(qualimap_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"qualimap_rnaseq")
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
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_qualimap_rnaseq_se(self):
        """
        QCPipeline: 'qualimap_rnaseq' QC module (SE data)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="qualimap_rnaseq",
                              description="Qualimap_rnaseq test",
                              seq_data_reads=['r1'],
                              index_reads=None,
                              qc_modules=("qualimap_rnaseq",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        qualimap_dir = os.path.join(qc_dir,"qualimap-rnaseq","human")
        for f in ("PJB1_S1_001/qualimapReport.html",
                  "PJB1_S1_001/rnaseq_qc_results.txt",
                  "PJB2_S2_001/qualimapReport.html",
                  "PJB2_S2_001/rnaseq_qc_results.txt"):
            self.assertTrue(os.path.exists(os.path.join(qualimap_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"qualimap_rnaseq")
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
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_qualimap_rnaseq_se_with_biological_samples(self):
        """
        QCPipeline: 'qualimap_rnaseq' QC module (SE data with biological samples)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Biological samples': 'PJB1' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="qualimap_rnaseq",
                              description="Qualimap_rnaseq test",
                              seq_data_reads=['r1'],
                              index_reads=None,
                              qc_modules=("qualimap_rnaseq",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        qualimap_dir = os.path.join(qc_dir,"qualimap-rnaseq","human")
        for f in ("PJB1_S1_001/qualimapReport.html",
                  "PJB1_S1_001/rnaseq_qc_results.txt"):
            self.assertTrue(os.path.exists(os.path.join(qualimap_dir,f)),
                            "%s: missing" % f)
        qualimap_dir = os.path.join(qc_dir,"qualimap-rnaseq","human")
        for f in ("PJB2_S2_001/qualimapReport.html",
                  "PJB2_S2_001/rnaseq_qc_results.txt"):
            self.assertFalse(os.path.exists(os.path.join(qualimap_dir,f)),
                            "%s: present" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"qualimap_rnaseq")
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
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_qualimap_rnaseq_se_split_lanes(self):
        """
        QCPipeline: 'qualimap_rnaseq' QC module (SE data, split by lane)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="qualimap_rnaseq",
                              description="Qualimap_rnaseq test",
                              seq_data_reads=['r1'],
                              index_reads=None,
                              qc_modules=("qualimap_rnaseq",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol,
                          split_fastqs_by_lane=True)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        qualimap_dir = os.path.join(qc_dir,"qualimap-rnaseq","human")
        for f in ("PJB1_S1_L001_001/qualimapReport.html",
                  "PJB1_S1_L001_001/rnaseq_qc_results.txt",
                  "PJB2_S2_L001_001/qualimapReport.html",
                  "PJB2_S2_L001_001/rnaseq_qc_results.txt"):
            self.assertTrue(os.path.exists(os.path.join(qualimap_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"qualimap_rnaseq")
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
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_qualimap_rnaseq_missing_star_index(self):
        """
        QCPipeline: 'qualimap_rnaseq' QC module (no STAR index)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="qualimap_rnaseq",
                              description="Qualimap_rnaseq test",
                              seq_data_reads=['r1'],
                              index_reads=None,
                              qc_modules=("qualimap_rnaseq",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        qualimap_dir = os.path.join(qc_dir,"qualimap-rnaseq","human")
        for f in ("PJB1_S1_001/qualimapReport.html",
                  "PJB1_S1_001/rnaseq_qc_results.txt",
                  "PJB2_S2_001/qualimapReport.html",
                  "PJB2_S2_001/rnaseq_qc_results.txt"):
            self.assertFalse(os.path.exists(os.path.join(qualimap_dir,f)),
                            "%s: present" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"qualimap_rnaseq")
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
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_qualimap_rnaseq_missing_gtf_file(self):
        """
        QCPipeline: 'qualimap_rnaseq' QC module (no GTF file)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="qualimap_rnaseq",
                              description="Qualimap_rnaseq test",
                              seq_data_reads=['r1'],
                              index_reads=None,
                              qc_modules=("qualimap_rnaseq",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        qualimap_dir = os.path.join(qc_dir,"qualimap-rnaseq","human")
        for f in ("PJB1_S1_001/qualimapReport.html",
                  "PJB1_S1_001/rnaseq_qc_results.txt",
                  "PJB2_S2_001/qualimapReport.html",
                  "PJB2_S2_001/rnaseq_qc_results.txt"):
            self.assertFalse(os.path.exists(os.path.join(qualimap_dir,f)),
                            "%s: present" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"qualimap_rnaseq")
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
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
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
