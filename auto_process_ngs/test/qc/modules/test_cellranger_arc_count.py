#######################################################################
# Unit tests for qc/pipeline.py ('cellranger-arc_count' QC module)
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

class TestQCPipelineCellrangerArcCount(BaseQCPipelineTestCase):
    """
    Tests for 'cellranger-arc_count' QC module
    """
    def test_qcpipeline_qc_modules_cellranger_arc_count_multiome_atac_no_linked_gex(self):
        """
        QCPipeline: 'cellranger-arc_count' QC module (SC multiome ATAC, no linked GEX)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"))
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
        # QC protocol
        protocol = QCProtocol(name="cellranger-arc_count_atac",
                              description="Cellranger-arc_count test",
                              seq_data_reads=['r1','r3'],
                              index_reads=[],
                              qc_modules=("cellranger-arc_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_ATAC")),
                          protocol)
        status = runqc.run(cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB_ATAC")
        for f in (
                "qc/cellranger-arc_count",
                "cellranger-arc_count"):
            self.assertFalse(os.path.exists(os.path.join(project_dir,f)),
                             "%s: present" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB_ATAC","README.info"))
        self.assertEqual(project_metadata.number_of_cells,None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_ATAC","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-arc_count_atac")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_ATAC,PJB2_ATAC")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB_ATAC","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_ATAC_S1_R1_001.fastq.gz,"
                         "PJB1_ATAC_S1_R2_001.fastq.gz,"
                         "PJB1_ATAC_S1_R3_001.fastq.gz,"
                         "PJB2_ATAC_S2_R1_001.fastq.gz,"
                         "PJB2_ATAC_S2_R2_001.fastq.gz,"
                         "PJB2_ATAC_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB_ATAC.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_ATAC",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_arc_count_multiome_gex_no_linked_atac(self):
        """
        QCPipeline: 'cellranger-arc_count' QC module (SC multiome GEX, no linked ATAC)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"))
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
        # QC protocol
        protocol = QCProtocol(name="cellranger-arc_count_gex",
                              description="Cellranger-arc_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-arc_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_GEX")),
                          protocol)
        status = runqc.run(cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB_GEX")
        for f in (
                "qc/cellranger-arc_count",
                "cellranger-arc_count"):
            self.assertFalse(os.path.exists(os.path.join(project_dir,f)),
                             "%s: present" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB_GEX","README.info"))
        self.assertEqual(project_metadata.number_of_cells,None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_GEX","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-arc_count_gex")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_GEX,PJB2_GEX")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB_GEX","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_GEX_S1_R1_001.fastq.gz,"
                         "PJB1_GEX_S1_R2_001.fastq.gz,"
                         "PJB2_GEX_S2_R1_001.fastq.gz,"
                         "PJB2_GEX_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB_GEX.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_GEX",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_arc_count_multiome_atac_100(self):
        """
        QCPipeline: 'cellranger-arc_count' QC module (SC multiome ATAC, Cellranger-ARC 1.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="1.0.0")
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
        # QC protocol
        protocol = QCProtocol(name="cellranger-arc_count_atac",
                              description="Cellranger-arc_count test",
                              seq_data_reads=['r1','r3'],
                              index_reads=[],
                              qc_modules=("cellranger-arc_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_ATAC")),
                          protocol)
        status = runqc.run(cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB_ATAC")
        for f in (
                "qc/cellranger-arc_count",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",
                "cellranger-arc_count",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB_ATAC","README.info"))
        self.assertEqual(project_metadata.number_of_cells,1488)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_ATAC","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-arc_count_atac")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_ATAC,PJB2_ATAC")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB_ATAC","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_ATAC_S1_R1_001.fastq.gz,"
                         "PJB1_ATAC_S1_R2_001.fastq.gz,"
                         "PJB1_ATAC_S1_R3_001.fastq.gz,"
                         "PJB2_ATAC_S2_R1_001.fastq.gz,"
                         "PJB2_ATAC_S2_R2_001.fastq.gz,"
                         "PJB2_ATAC_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"1.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB_ATAC.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_ATAC",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_arc_count_multiome_atac_200(self):
        """
        QCPipeline: 'cellranger-arc_count' QC module (SC multiome ATAC, Cellranger-ARC 2.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="2.0.0")
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
        # QC protocol
        protocol = QCProtocol(name="cellranger-arc_count_atac",
                              description="Cellranger-arc_count test",
                              seq_data_reads=['r1','r3'],
                              index_reads=[],
                              qc_modules=("cellranger-arc_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_ATAC")),
                          protocol)
        status = runqc.run(cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB_ATAC")
        for f in (
                "qc/cellranger-arc_count",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",
                "cellranger-arc_count",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB_ATAC","README.info"))
        self.assertEqual(project_metadata.number_of_cells,1570)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_ATAC","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-arc_count_atac")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_ATAC,PJB2_ATAC")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB_ATAC","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_ATAC_S1_R1_001.fastq.gz,"
                         "PJB1_ATAC_S1_R2_001.fastq.gz,"
                         "PJB1_ATAC_S1_R3_001.fastq.gz,"
                         "PJB2_ATAC_S2_R1_001.fastq.gz,"
                         "PJB2_ATAC_S2_R2_001.fastq.gz,"
                         "PJB2_ATAC_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB_ATAC.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_ATAC",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_arc_count_multiome_atac_210(self):
        """
        QCPipeline: 'cellranger-arc_count' QC module (SC multiome ATAC, Cellranger-ARC 2.1.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="2.1.0")
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
        # QC protocol
        protocol = QCProtocol(name="cellranger-arc_count_atac",
                              description="Cellranger-arc_count test",
                              seq_data_reads=['r1','r3'],
                              index_reads=[],
                              qc_modules=("cellranger-arc_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_ATAC")),
                          protocol)
        status = runqc.run(cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB_ATAC")
        for f in (
                "qc/cellranger-arc_count",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",
                "cellranger-arc_count",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB_ATAC","README.info"))
        self.assertEqual(project_metadata.number_of_cells,1452)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_ATAC","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-arc_count_atac")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_ATAC,PJB2_ATAC")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB_ATAC","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_ATAC_S1_R1_001.fastq.gz,"
                         "PJB1_ATAC_S1_R2_001.fastq.gz,"
                         "PJB1_ATAC_S1_R3_001.fastq.gz,"
                         "PJB2_ATAC_S2_R1_001.fastq.gz,"
                         "PJB2_ATAC_S2_R2_001.fastq.gz,"
                         "PJB2_ATAC_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB_ATAC.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_ATAC",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_arc_count_multiome_snatac(self):
        """
        QCPipeline: 'cellranger-arc_count' QC module (SC multiome snATAC)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="2.0.0")
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
                                           'Library type': 'snATAC' })
        p.create(top_dir=self.wd)
        # Make mock multiome GEX analysis project (with QC outputs)
        p2 = MockAnalysisProject("PJB_GEX",("PJB1_GEX_S1_R1_001.fastq.gz",
                                            "PJB1_GEX_S1_R2_001.fastq.gz",
                                            "PJB2_GEX_S2_R1_001.fastq.gz",
                                            "PJB2_GEX_S2_R2_001.fastq.gz",),
                                 metadata={ 'Organism': 'Human',
                                            'Single cell platform':
                                            '10xGenomics Single Cell Multiome',
                                            'Library type': 'snGEX' })
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
        # QC protocol
        protocol = QCProtocol(name="cellranger-arc_count_atac",
                              description="Cellranger-arc_count test",
                              seq_data_reads=['r1','r3'],
                              index_reads=[],
                              qc_modules=("cellranger-arc_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_ATAC")),
                          protocol)
        status = runqc.run(cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB_ATAC")
        for f in (
                "qc/cellranger-arc_count",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",
                "cellranger-arc_count",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB_ATAC","README.info"))
        self.assertEqual(project_metadata.number_of_cells,1570)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_ATAC","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-arc_count_atac")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_ATAC,PJB2_ATAC")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB_ATAC","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_ATAC_S1_R1_001.fastq.gz,"
                         "PJB1_ATAC_S1_R2_001.fastq.gz,"
                         "PJB1_ATAC_S1_R3_001.fastq.gz,"
                         "PJB2_ATAC_S2_R1_001.fastq.gz,"
                         "PJB2_ATAC_S2_R2_001.fastq.gz,"
                         "PJB2_ATAC_S2_R3_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB_ATAC.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_ATAC",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_arc_count_multiome_gex_100(self):
        """
        QCPipeline: 'cellranger-arc_count' QC module (SC multiome GEX, Cellranger-ARC 1.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="1.0.0")
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
        # QC protocol
        protocol = QCProtocol(name="cellranger-arc_count_gex",
                              description="Cellranger-arc_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-arc_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_GEX")),
                          protocol)
        status = runqc.run(cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB_GEX")
        for f in (
                "qc/cellranger-arc_count",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                "qc/cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv",
                "cellranger-arc_count",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                "cellranger-arc_count/1.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB_GEX","README.info"))
        self.assertEqual(project_metadata.number_of_cells,1488)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_GEX","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-arc_count_gex")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_GEX,PJB2_GEX")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB_GEX","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_GEX_S1_R1_001.fastq.gz,"
                         "PJB1_GEX_S1_R2_001.fastq.gz,"
                         "PJB2_GEX_S2_R1_001.fastq.gz,"
                         "PJB2_GEX_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"1.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB_GEX.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_GEX",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_arc_count_multiome_gex_200(self):
        """
        QCPipeline: 'cellranger-arc_count' QC module (SC multiome GEX, Cellranger-ARC 2.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="2.0.0")
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
        # QC protocol
        protocol = QCProtocol(name="cellranger-arc_count_gex",
                              description="Cellranger-arc_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-arc_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_GEX")),
                          protocol)
        status = runqc.run(cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB_GEX")
        for f in (
                "qc/cellranger-arc_count",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv",
                "cellranger-arc_count",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB_GEX","README.info"))
        self.assertEqual(project_metadata.number_of_cells,1570)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_GEX","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-arc_count_gex")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_GEX,PJB2_GEX")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB_GEX","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_GEX_S1_R1_001.fastq.gz,"
                         "PJB1_GEX_S1_R2_001.fastq.gz,"
                         "PJB2_GEX_S2_R1_001.fastq.gz,"
                         "PJB2_GEX_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB_GEX.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_GEX",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_arc_count_multiome_gex_210(self):
        """
        QCPipeline: 'cellranger-arc_count' QC module (SC multiome GEX, Cellranger-ARC 2.1.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="2.1.0")
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
        # QC protocol
        protocol = QCProtocol(name="cellranger-arc_count_gex",
                              description="Cellranger-arc_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-arc_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_GEX")),
                          protocol)
        status = runqc.run(cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB_GEX")
        for f in (
                "qc/cellranger-arc_count",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                "qc/cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv",
                "cellranger-arc_count",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                "cellranger-arc_count/2.1.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB_GEX","README.info"))
        self.assertEqual(project_metadata.number_of_cells,1452)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_GEX","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-arc_count_gex")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_GEX,PJB2_GEX")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB_GEX","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_GEX_S1_R1_001.fastq.gz,"
                         "PJB1_GEX_S1_R2_001.fastq.gz,"
                         "PJB2_GEX_S2_R1_001.fastq.gz,"
                         "PJB2_GEX_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB_GEX.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_GEX",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_arc_count_multiome_sngex(self):
        """
        QCPipeline: 'cellranger-arc_count' QC module (SC multiome snGEX)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="2.0.0")
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
                                           'Library type': 'snGEX' })
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
                                            'Library type': 'snATAC' })
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
        # QC protocol
        protocol = QCProtocol(name="cellranger-arc_count_gex",
                              description="Cellranger-arc_count test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger-arc_count",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_GEX")),
                          protocol)
        status = runqc.run(cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB_GEX")
        for f in (
                "qc/cellranger-arc_count",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv",
                "cellranger-arc_count",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB_GEX","README.info"))
        self.assertEqual(project_metadata.number_of_cells,1570)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_GEX","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger-arc_count_gex")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_GEX,PJB2_GEX")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB_GEX","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_GEX_S1_R1_001.fastq.gz,"
                         "PJB1_GEX_S1_R2_001.fastq.gz,"
                         "PJB2_GEX_S2_R1_001.fastq.gz,"
                         "PJB2_GEX_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB_GEX.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_GEX",f)),
                            "Missing %s" % f)