#######################################################################
# Unit tests for qc/modules/cellranger_count.py ('cellanger*_count' QC module)
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.metadata import AnalysisProjectInfo
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.protocols import QCProtocol
from auto_process_ngs.qc.pipeline import QCPipeline
from auto_process_ngs.qc.modules.cellranger_count import check_cellranger_count_outputs
from auto_process_ngs.qc.modules.cellranger_count import check_cellranger_atac_count_outputs
from auto_process_ngs.qc.modules.cellranger_count import check_cellranger_arc_count_outputs
from auto_process_ngs.qc.modules.cellranger_count import filter_10x_pipelines
from ..protocols import BaseQCPipelineTestCase

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Polling interval for pipeline
POLL_INTERVAL = 0.1

class TestQCPipelineCellrangerCount(BaseQCPipelineTestCase):
    """
    Tests for 'cellranger_count' QC module
    """
    def test_qcpipeline_qc_modules_cellranger_count_scrnaseq_310(self):
        """
        QCPipeline: 'cellranger_count' QC module (scRNA-seq, Cellranger v3.1.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="3.1.0",
                                 assert_include_introns=False)
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-cellranger-GRCh38-1.2.0' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
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
                "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,4544)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"3.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-GRCh38-1.2.0")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_scrnaseq_501(self):
        """
        QCPipeline: 'cellranger_count' QC module (scRNA-seq, Cellranger v5.0.1)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="5.0.1",
                                 assert_include_introns=False)
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
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
                "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,4544)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"5.0.1")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_scrnaseq_600(self):
        """
        QCPipeline: 'cellranger_count' QC module (scRNA-seq, Cellranger v6.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="6.0.0",
                                 assert_include_introns=False)
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
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
                "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,4544)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"6.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_scrnaseq_700(self):
        """
        QCPipeline: 'cellranger_count' QC module (scRNA-seq, Cellranger v7.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="7.0.0",
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
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,11058)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"7.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_scrnaseq_800(self):
        """
        QCPipeline: 'cellranger_count' QC module (scRNA-seq, Cellranger v8.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="8.0.0",
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
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,11058)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_snrnaseq_310(self):
        """
        QCPipeline: 'cellranger_count' QC module (snRNA-seq, Cellranger v3.1.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="3.1.0",
                                 assert_include_introns=False)
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_premrna_references=
                           { 'human':
                             '/data/refdata-cellranger-GRCh38-1.2.0_premrna' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
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
                "cellranger_count/3.1.0/refdata-cellranger-GRCh38-1.2.0_premrna/PJB2/outs/metrics_summary.csv",):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,4544)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"3.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-GRCh38-1.2.0_premrna")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_snrnaseq_501(self):
        """
        QCPipeline: 'cellranger_count' QC module (snRNA-seq, Cellranger v5.0.1)
        """
        # Make mock QC executables
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
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
                "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,4544)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"5.0.1")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_snrnaseq_600(self):
        """
        QCPipeline: 'cellranger_count' QC module (snRNA-seq, Cellranger v6.0.0)
        """
        # Make mock QC executables
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
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
                "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,4544)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"6.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_snrnaseq_700(self):
        """
        QCPipeline: 'cellranger_count' QC module (snRNA-seq, Cellranger v7.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="7.0.0",
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,11058)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"7.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_snrnaseq_800(self):
        """
        QCPipeline: 'cellranger_count' QC module (snRNA-seq, Cellranger v8.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="8.0.0",
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,11058)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_extra_project(self):
        """
        QCPipeline: 'cellranger_count' QC module (scRNA-seq with extra project)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="7.0.0",
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
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # Make second mock 10x analysis project
        p2 = MockAnalysisProject("PJB2",("PJB3_S3_R1_001.fastq.gz",
                                         "PJB3_S3_R2_001.fastq.gz",
                                         "PJB4_S4_R1_001.fastq.gz",
                                         "PJB4_S4_R2_001.fastq.gz",),
                                 metadata={ 'Organism': 'Human',
                                            'Single cell platform':
                                            '10xGenomics Chromium 3\'v2',
                                            'Library type': 'scRNA-seq' })
        p2.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           cellranger_extra_projects=[
                               AnalysisProject("PJB2",
                                               os.path.join(self.wd,
                                                            "PJB2"))],
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB3/_cmdline",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB3/outs/web_summary.html",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB3/outs/metrics_summary.csv",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB4/_cmdline",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB4/outs/web_summary.html",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB4/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
            "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB3/_cmdline",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB3/outs/web_summary.html",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB3/outs/metrics_summary.csv",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB4/_cmdline",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB4/outs/web_summary.html",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB4/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "Missing %s" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,11058)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"7.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_no_project_metadata(self):
        """
        QCPipeline: 'cellranger_count' QC module (no project metadata)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version='5.0.1')
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol,
                          organism="human")
        status = runqc.run(cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
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
                "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # No metadata file so number of cells cannot be set
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"human")
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
        self.assertEqual(qc_info.cellranger_version,"5.0.1")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_no_references(self):
        """
        QCPipeline: 'cellranger_count' QC module (scRNA-seq, no reference data)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="5.0.1")
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "cellranger_count"):
            self.assertFalse(os.path.exists(os.path.join(project_dir,f)),
                             "%s: present, should be missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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

    def test_qcpipeline_qc_modules_cellranger_count_specify_exe(self):
        """
        QCPipeline: 'cellranger_count' QC module (scRNA-seq, specify cellranger exe)
        """
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           cellranger_exe=os.path.join(cellranger_bin,
                                                       "cellranger"),
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
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
                "cellranger_count/5.0.1/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,4544)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"5.0.1")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_force_cells(self):
        """
        QCPipeline: 'cellranger_count' QC module (scRNA-seq, use --force-cells)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 assert_force_cells=10000)
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           cellranger_force_cells=10000,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "qc/cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        # NB should be 10000 but this test is correct (as stock
        # data are used for test outputs)
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,11058)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"7.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_dont_set_cell_count(self):
        """
        QCPipeline: 'cellranger_count' QC module (don't set cell count)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="8.0.0",
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
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=(
                                  "cellranger_count(set_cell_count=false)",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_dont_set_metadata(self):
        """
        QCPipeline: 'cellranger_count' QC module (don't set metadata)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="8.0.0",
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
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=(
                                  "cellranger_count(set_metadata=false)",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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

    def test_qcpipeline_qc_modules_cellranger_count_specify_chemistry(self):
        """
        QCPipeline: 'cellranger_count' QC module (specify chemistry)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="8.0.0",
                                 assert_chemistry="ARC-v1",
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
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=(
                                  "cellranger_count(chemistry=ARC-v1)",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,11058)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_specify_library(self):
        """
        QCPipeline: 'cellranger_count' QC module (specify library, Cellranger v5.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="5.0.0",
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
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=(
                                  "cellranger_count(library=snRNA-seq)",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "qc/cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "qc/cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "qc/cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "qc/cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "qc/cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB1/_cmdline",
                "cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/web_summary.html",
                "cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB1/outs/metrics_summary.csv",
                "cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB2/_cmdline",
                "cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/web_summary.html",
                "cellranger_count/5.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,4544)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
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
        self.assertEqual(qc_info.cellranger_version,"5.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_multi_config(self):
        """
        QCPipeline: 'cellranger_count' QC module (use cellranger multi config)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="8.0.0")
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
reference,/data/refdata-cellranger-gex-GRCh38-2024-A
create-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,%s,any,PJB1,gene expression,
PJB2_MC,%s,any,PJB2,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBA,CMO301,PBA
PBB,CMO302,PBB
""" % (fastq_dir,fastq_dir))
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=(
                                  "cellranger_count(cellranger_use_multi_config=true)",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2024-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/_cmdline",
                "qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/outs/web_summary.html",
                "qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/outs/metrics_summary.csv",
                "cellranger_count",
                "cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/_cmdline",
                "cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/outs/web_summary.html",
                "cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/outs/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        for f in (
                "qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2_MC",
                "cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2_MC"):
            self.assertFalse(os.path.exists(os.path.join(project_dir,f)),
                             "%s: present" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        ##self.assertEqual(project_metadata.number_of_cells,5529)
        self.assertEqual(project_metadata.number_of_cells,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_GEX")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_GEX_S1_R1_001.fastq.gz,"
                         "PJB1_GEX_S1_R2_001.fastq.gz,"
                         "PJB2_MC_S2_R1_001.fastq.gz,"
                         "PJB2_MC_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-gex-GRCh38-2024-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_count_multi_config_no_file(self):
        """
        QCPipeline: 'cellranger_count' QC module (use cellranger multi config, no file)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="8.0.0")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x Cellplex analysis project without
        # 10x multi config file
        p = MockAnalysisProject("PJB",("PJB1_GEX_S1_R1_001.fastq.gz",
                                       "PJB1_GEX_S1_R2_001.fastq.gz",
                                       "PJB2_MC_S2_R1_001.fastq.gz",
                                       "PJB2_MC_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 3\'v3',
                                           'Library type': 'CellPlex' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger_count",
                              description="Cellranger_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=(
                                  "cellranger_count(cellranger_use_multi_config=true)",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-gex-GRCh38-2024-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_count",
                "qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX",
                "qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2_MC",
                "cellranger_count",
                "cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX",
                "cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2_MC"):
            self.assertFalse(os.path.exists(os.path.join(project_dir,f)),
                             "%s: present" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_count")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_GEX,PJB2_MC")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_GEX_S1_R1_001.fastq.gz,"
                         "PJB1_GEX_S1_R2_001.fastq.gz,"
                         "PJB2_MC_S2_R1_001.fastq.gz,"
                         "PJB2_MC_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2024-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

class TestCheckCellrangerCountOutputs(unittest.TestCase):
    """
    Tests for the 'check_cellranger_count_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCheckCellrangerCountOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_cellranger_count_outputs_singlecell_missing(self):
        """
        check_cellranger_count_outputs: cellranger count output missing (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="singlecell",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),["PJB1",])

    def test_check_cellranger_count_outputs_singlecell_present(self):
        """
        check_cellranger_count_outputs: cellranger count output present (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="singlecell",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs()
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),[])

    def test_check_cellranger_count_outputs_singlecell_present_with_prefix(self):
        """
        check_cellranger_count_outputs: cellranger count output present with prefix (singlecell)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        cellranger_count_prefix = "cellranger_count/5.0.1/refdata-gex-mm10-2020-A"
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="singlecell",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            prefix=cellranger_count_prefix)
        # Check the outputs
        self.assertEqual(
            check_cellranger_count_outputs(project,
                                           prefix=cellranger_count_prefix),[])

    def test_check_cellranger_count_outputs_10x_scRNAseq_missing(self):
        """
        check_cellranger_count_outputs: cellranger count output missing (10x_scRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3", })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),["PJB1",])

    def test_check_cellranger_count_outputs_10x_scRNAseq_present(self):
        """
        check_cellranger_count_outputs: cellranger count output present (10x_scRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs()
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),[])

    def test_check_cellranger_count_outputs_10x_scRNAseq_present_with_prefix(self):
        """
        check_cellranger_count_outputs: cellranger count output present with prefix (10x_scRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        cellranger_count_prefix = "cellranger_count/5.0.1/refdata-gex-mm10-2020-A"
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            prefix=cellranger_count_prefix)
        # Check the outputs
        self.assertEqual(
            check_cellranger_count_outputs(project,
                                           prefix=cellranger_count_prefix),[])

    def test_check_cellranger_count_outputs_10x_snRNAseq_missing(self):
        """
        check_cellranger_count_outputs: cellranger count output missing (10x_snRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3", })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_snRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),["PJB1",])

    def test_check_cellranger_count_outputs_10x_snRNAseq_present(self):
        """
        check_cellranger_count_outputs: cellranger count output present (10x_snRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_snRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs()
        # Check the outputs
        self.assertEqual(check_cellranger_count_outputs(project),[])

    def test_check_cellranger_count_outputs_10x_snRNAseq_present_with_prefix(self):
        """
        check_cellranger_count_outputs: cellranger count output present with prefix (10x_snRNAseq)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        cellranger_count_prefix = "cellranger_count/3.1.0/refdata-cellranger-GRCh38-3.0.0_premRNA_patch"
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_snRNAseq",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            prefix=cellranger_count_prefix)
        # Check the outputs
        self.assertEqual(
            check_cellranger_count_outputs(project,
                                           prefix=cellranger_count_prefix),[])

class TestCheckCellrangerAtacCountOutputs(unittest.TestCase):
    """
    Tests for the 'check_cellranger_atac_count_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCheckCellrangerAtacCountOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_cellranger_atac_count_outputs_10x_scATAC_missing(self):
        """
        check_cellranger_atac_count_outputs: cellranger-atac count output missing (10x_scATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell ATAC" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scATAC",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check the outputs
        self.assertEqual(check_cellranger_atac_count_outputs(project),["PJB1",])

    def test_check_cellranger_atac_count_outputs_10x_scATAC_present(self):
        """
        check_cellranger_atac_count_outputs: cellranger-atac count output present (10x_scATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell ATAC" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scATAC",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            cellranger="cellranger-atac")
        # Check the outputs
        self.assertEqual(check_cellranger_atac_count_outputs(project),[])

    def test_check_cellranger_atac_count_outputs_10x_scATAC_present_with_prefix(self):
        """
        check_cellranger_atac_count_outputs: cellranger-atac count output present with prefix (10x_scATAC)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell ATAC" })
        p.create(top_dir=self.wd)
        project = AnalysisProject("PJB",os.path.join(self.wd,"PJB"))
        cellranger_count_prefix = "cellranger_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0"
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_scATAC",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            cellranger="cellranger-atac",
            prefix=cellranger_count_prefix)
        # Check the outputs
        self.assertEqual(
            check_cellranger_atac_count_outputs(project,
                                                prefix=cellranger_count_prefix),
            [])

class TestCheckCellrangerArcCountOutputs(unittest.TestCase):
    """
    Tests for the 'check_cellranger_arc_count_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCheckCellrangerArcCountOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_check_cellranger_arc_count_outputs_no_libraries_csv(self):
        """
        check_cellranger_arc_count_outputs: cellranger-arc count no libraries.csv (10x_Multiome_*)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell Multiome",
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_Multiome_GEX",
            include_fastq_strand=False,
            include_multiqc=False)
        # Check the outputs
        self.assertEqual(check_cellranger_arc_count_outputs(project),[])

    def test_check_cellranger_arc_count_outputs_10x_multiome_missing(self):
        """
        check_cellranger_arc_count_outputs: cellranger-arc count output missing (10x_Multiome_*)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell Multiome",
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_Multiome_GEX",
            include_fastq_strand=False,
            include_multiqc=False)
        # Add libraries.csv
        with open(os.path.join(project.qc_dir,"libraries.PJB1.csv"),'wt') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(check_cellranger_arc_count_outputs(project),["PJB1",])

    def test_check_cellranger_arc_count_outputs_10x_multiome_present(self):
        """
        check_cellranger_arc_count_outputs: cellranger-arc count output present (10x_Multiome_*)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell Multiome",
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_Multiome_GEX",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            cellranger="cellranger-arc")
        # Add libraries.csv
        with open(os.path.join(project.qc_dir,"libraries.PJB1.csv"),'wt') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(check_cellranger_arc_count_outputs(project),[])

    def test_check_cellranger_arc_count_outputs_10x_multiome_present_with_prefix(self):
        """
        check_cellranger_arc_count_outputs: cellranger-arc count output present with prefix (10x_Multiome_*)
        """
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Single Cell Multiome",
                                           'Library type': 'GEX' })
        p.create(top_dir=self.wd)
        project = AnalysisProject(os.path.join(self.wd,"PJB"))
        cellranger_count_prefix = "cellranger_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A"
        UpdateAnalysisProject(project).add_qc_outputs(
            protocol="10x_Multiome_GEX",
            include_fastq_strand=False,
            include_multiqc=False)
        UpdateAnalysisProject(project).add_cellranger_count_outputs(
            cellranger="cellranger-arc",
            prefix=cellranger_count_prefix)
        # Add libraries.csv
        with open(os.path.join(project.qc_dir,"libraries.PJB1.csv"),'wt') as fp:
            fp.write("")
        # Check the outputs
        self.assertEqual(
            check_cellranger_arc_count_outputs(project,
                                           prefix=cellranger_count_prefix),[])

class TestFilter10xPipelines(unittest.TestCase):

    def test_filter_10x_pipelines(self):
        """
        filter_10x_pipelines: check pipelines are correctly filtered
        """
        pipelines = (("cellranger","5.0.0","refdata-gex-GRCh38-2.0.0"),
                     ("cellranger","6.0.1","refdata-gex-GRCh38-2020"),
                     ("cellranger","6.1.2","refdata-gex-GRCh38-2020"),
                     ("cellranger-arc","2.0.0","refdata-arc-GRCh38-2020"),
                     ("cellranger-arc","2.0.0","refdata-arc-mm10-2020"),)
        # Specific pipeline
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger","6.1.2","refdata-gex-GRCh38-2020"),
                pipelines),
            [("cellranger","6.1.2","refdata-gex-GRCh38-2020"),])
        # Specific package and reference data, any version
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger","*","refdata-gex-GRCh38-2020"),
                pipelines),
            [("cellranger","6.0.1","refdata-gex-GRCh38-2020"),
             ("cellranger","6.1.2","refdata-gex-GRCh38-2020"),])
        # Specific package and version, any reference
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger-arc","2.0.0","*"),
                pipelines),
            [("cellranger-arc","2.0.0","refdata-arc-GRCh38-2020"),
             ("cellranger-arc","2.0.0","refdata-arc-mm10-2020"),])
        # Specific package, any version, any reference
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger","*","*"),
                pipelines),
            [("cellranger","5.0.0","refdata-gex-GRCh38-2.0.0"),
             ("cellranger","6.0.1","refdata-gex-GRCh38-2020"),
             ("cellranger","6.1.2","refdata-gex-GRCh38-2020"),])
        # No matching package
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger-atac","*","*"),
                pipelines),[])
        # No matching version
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger","3.0.1","*"),
                pipelines),[])
        # No matching reference
        self.assertEqual(
            filter_10x_pipelines(
                ("cellranger","*","refdata-gex-mm10-2020"),
                pipelines),[])
