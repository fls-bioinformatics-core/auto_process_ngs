#######################################################################
# Unit tests for qc/pipeline.py ('cellranger-atac_count' QC module)
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
from ..protocols import BaseQCPipelineTestCase

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Polling interval for pipeline
POLL_INTERVAL = 0.1

class TestQCPipelineCellrangerAtacCount(BaseQCPipelineTestCase):
    """
    Tests for 'cellranger-atac_count' QC module
    """
    def test_qcpipeline_qc_modules_cellranger_atac_count_120(self):
        """
        QCPipeline: 'cellranger-atac_count' QC module (Cellranger-atac 1.2.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"),
                                 version="1.2.0")
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
                                           '10x Single Cell ATAC',
                                           'Library type': 'snATAC-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger-atac_count",
                              description="Cellranger-atac_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-atac_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_atac_references=
                           { 'human':
                             '/data/refdata-cellranger-atac-GRCh38-1.2.0' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger-atac_count",
                "qc/cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/_cmdline",
                "qc/cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/outs/web_summary.html",
                "qc/cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/outs/summary.csv",
                "qc/cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/_cmdline",
                "qc/cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/outs/web_summary.html",
                "qc/cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/outs/summary.csv",
                "cellranger-atac_count",
                "cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/_cmdline",
                "cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/outs/web_summary.html",
                "cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB1/outs/summary.csv",
                "cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/_cmdline",
                "cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/outs/web_summary.html",
                "cellranger-atac_count/1.2.0/refdata-cellranger-atac-GRCh38-1.2.0/PJB2/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,11364)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-atac_count")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB1_S1_R3_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz,"
                         "PJB2_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"1.2.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-atac-GRCh38-1.2.0")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_atac_count_200(self):
        """
        QCPipeline: 'cellranger-atac_count' QC module (Cellranger-atac 2.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"),
                                 version="2.0.0")
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
                                           '10x Single Cell ATAC',
                                           'Library type': 'snATAC-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger-atac_count",
                              description="Cellranger-atac_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-atac_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_atac_references=
                           { 'human':
                             '/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger-atac_count",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv",
                "cellranger-atac_count",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,7164)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-atac_count")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB1_S1_R3_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz,"
                         "PJB2_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_atac_count_force_cells(self):
        """
        QCPipeline: 'cellranger-atac_count' QC module (use --force-cells)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"),
                                 assert_force_cells=10000)
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
                                           '10x Single Cell ATAC',
                                           'Library type': 'snATAC-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger-atac_count",
                              description="Cellranger-atac_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-atac_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_atac_references=
                           { 'human':
                             '/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0' },
                           cellranger_force_cells=10000,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger-atac_count",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv",
                "cellranger-atac_count",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,7164)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-atac_count")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB1_S1_R3_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz,"
                         "PJB2_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_atac_count_dont_set_cell_count(self):
        """
        QCPipeline: 'cellranger-atac_count' QC module (don't set cell count)
        """
        # Make mock QC executables
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
                                           '10x Single Cell ATAC',
                                           'Library type': 'snATAC-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger-atac_count",
                              description="Cellranger-atac_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-atac_count(set_cell_count=false)",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_atac_references=
                           { 'human':
                             '/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger-atac_count",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv",
                "cellranger-atac_count",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-atac_count")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB1_S1_R3_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz,"
                         "PJB2_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_atac_count_dont_set_metadata(self):
        """
        QCPipeline: 'cellranger-atac_count' QC module (don't set metadata)
        """
        # Make mock QC executables
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
                                           '10x Single Cell ATAC',
                                           'Library type': 'snATAC-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger-atac_count",
                              description="Cellranger-atac_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-atac_count(set_metadata=false)",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_atac_references=
                           { 'human':
                             '/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger-atac_count",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv",
                "cellranger-atac_count",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-atac_count")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB1_S1_R3_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz,"
                         "PJB2_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_atac_count_specify_chemistry(self):
        """
        QCPipeline: 'cellranger-atac_count' QC module (specify chemistry)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"),
                                 assert_chemistry="ARC-v1")
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
                                           '10x Single Cell ATAC',
                                           'Library type': 'snATAC-seq' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger-atac_count",
                              description="Cellranger-atac_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-atac_count(chemistry=ARC-v1)",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_atac_references=
                           { 'human':
                             '/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger-atac_count",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv",
                "cellranger-atac_count",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/_cmdline",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/web_summary.html",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1/outs/summary.csv",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/_cmdline",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/web_summary.html",
                "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,7164)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-atac_count")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB1_S1_R2_001.fastq.gz,"
                         "PJB1_S1_R3_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz,"
                         "PJB2_S2_R2_001.fastq.gz,"
                         "PJB2_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
