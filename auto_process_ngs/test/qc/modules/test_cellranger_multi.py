#######################################################################
# Unit tests for qc/pipeline.py ('cellranger_multi' QC module)
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
from auto_process_ngs.qc.modules.cellranger_multi import expected_outputs
from ..protocols import BaseQCPipelineTestCase

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Polling interval for pipeline
POLL_INTERVAL = 0.1

class TestQCPipelineCellrangerMulti(BaseQCPipelineTestCase):
    """
    Tests for 'cellranger_multi' QC module
    """
    def test_qcpipeline_qc_modules_cellranger_multi_cellplex_600(self):
        """
        QCPipeline: 'cellranger_multi' QC module (Cellplex, Cellranger v6.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="cellplex",
                                 version="6.0.0")
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_multi",
                              description="Cellranger_multi test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_multi",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-cellranger-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
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
                "cellranger_multi/6.0.0/refdata-cellranger-gex-GRCh38-2020-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,10350)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_multi")
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
        self.assertEqual(qc_info.cellranger_version,"6.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_multi_cellplex_800(self):
        """
        QCPipeline: 'cellranger_multi' QC module (Cellplex, Cellranger v8.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="cellplex",
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
        protocol = QCProtocol(name="cellranger_multi",
                              description="Cellranger_multi test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_multi",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-cellranger-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_multi",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/web_summary.html",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/web_summary.html",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv",
                "cellranger_multi",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/web_summary.html",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/web_summary.html",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,3164)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_multi")
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

    def test_qcpipeline_qc_modules_cellranger_multi_cellplex_900(self):
        """
        QCPipeline: 'cellranger_multi' QC module (Cellplex, Cellranger v9.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="cellplex",
                                 version="9.0.0")
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
        protocol = QCProtocol(name="cellranger_multi",
                              description="Cellranger_multi test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_multi",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-cellranger-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_multi",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/web_summary.html",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/web_summary.html",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv",
                "cellranger_multi",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/web_summary.html",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/web_summary.html",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,3142)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_multi")
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
        self.assertEqual(qc_info.cellranger_version,"9.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-gex-GRCh38-2024-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_multi_cellplex_no_config_file(self):
        """
        QCPipeline: 'cellranger_multi' QC module (CellPlex, no config file)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="cellplex",
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
        # QC protocol
        protocol = QCProtocol(name="cellranger_multi",
                              description="Cellranger_multi test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_multi",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-cellranger-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_multi",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/web_summary.html",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/web_summary.html",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv",
                "cellranger_multi",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/web_summary.html",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/web_summary.html",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv"):
            self.assertFalse(os.path.exists(os.path.join(project_dir,f)),
                             "%s: present" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_multi")
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
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_multi_cellplex_bad_config(self):
        """
        QCPipeline: 'cellranger_multi' QC module (Cellplex, bad config file)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="cellplex",
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
        # Add "bad" cellranger multi config.csv file
        with open(os.path.join(self.wd,
                               "PJB",
                               "10x_multi_config.csv"),'wt') as fp:
            fastq_dir = os.path.join(self.wd,
                                     "PJB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2024-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,%s,any,PJB1,gene expression,
PJB2_MC,%s,any,PJB2,Multiplex Capture,

[samples]
sample_id,cmo_ids,description
PBA,CMO302,PBA
PBB,CMO302,PBB
""" % (fastq_dir,fastq_dir))
        # QC protocol
        protocol = QCProtocol(name="cellranger_multi",
                              description="Cellranger_multi test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_multi",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-cellranger-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,1)
        # Check outputs
        project_dir = os.path.join(self.wd,"PJB")
        for f in (
                "qc/cellranger_multi",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/web_summary.html",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/web_summary.html",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv",
                "cellranger_multi",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/web_summary.html",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBA/metrics_summary.csv",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/web_summary.html",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PBB/metrics_summary.csv",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv"):
            self.assertFalse(os.path.exists(os.path.join(project_dir,f)),
                             "%s: present (should be missing)" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells, None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol, None)
        self.assertEqual(qc_info.protocol_specification, None)
        self.assertEqual(qc_info.organism, None)
        self.assertEqual(qc_info.seq_data_samples, None)
        self.assertEqual(qc_info.fastq_dir, None)
        self.assertEqual(qc_info.fastqs, None)
        self.assertEqual(qc_info.fastqs_split_by_lane, None)
        self.assertEqual(qc_info.fastq_screens, None)
        self.assertEqual(qc_info.star_index, None)
        self.assertEqual(qc_info.annotation_bed, None)
        self.assertEqual(qc_info.annotation_gtf, None)
        self.assertEqual(qc_info.cellranger_version, None)
        self.assertEqual(qc_info.cellranger_refdata, None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Present (should be missing): %s" % f)

    def test_qcpipeline_qc_modules_cellranger_multi_flex_710(self):
        """
        QCPipeline: 'cellranger_multi' QC module (Flex, Cellranger v7.1.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="flex",
                                 version="7.1.0")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x Flex analysis project
        p = MockAnalysisProject("PJB",("PJB1_Flex_S1_R1_001.fastq.gz",
                                       "PJB1_Flex_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 3\'v3',
                                           'Library type': 'Flex' })
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
probe-set,/data/Probe_Set_v1.0_GRCh38-2020-A.csv
no-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_Flex,{fastq_dir},any,PJB1,Gene Expression,

[samples]
sample_id,probe_barcode_ids,description
PB1,BC001,PB1
PB2,BC002,PB2
""".format(fastq_dir=fastq_dir))
        # QC protocol
        protocol = QCProtocol(name="cellranger_multi",
                              description="Cellranger_multi test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_multi",))
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
                "qc/cellranger_multi",
                "qc/cellranger_multi/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/_cmdline",
                "qc/cellranger_multi/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PB1/web_summary.html",
                "qc/cellranger_multi/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PB1/metrics_summary.csv",
                "qc/cellranger_multi/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PB2/web_summary.html",
                "qc/cellranger_multi/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PB2/metrics_summary.csv",
                "cellranger_multi",
                "cellranger_multi/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/_cmdline",
                "cellranger_multi/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PB1/web_summary.html",
                "cellranger_multi/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PB1/metrics_summary.csv",
                "cellranger_multi/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PB2/web_summary.html",
                "cellranger_multi/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/outs/per_sample_outs/PB2/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,3138)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_multi")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_Flex")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_Flex_S1_R1_001.fastq.gz,"
                         "PJB1_Flex_S1_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"7.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,
                         "/data/Probe_Set_v1.0_GRCh38-2020-A.csv")
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_multi_flex_800(self):
        """
        QCPipeline: 'cellranger_multi' QC module (Flex, Cellranger v8.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="flex",
                                 version="8.0.0")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x Flex analysis project
        p = MockAnalysisProject("PJB",("PJB1_Flex_S1_R1_001.fastq.gz",
                                       "PJB1_Flex_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 3\'v3',
                                           'Library type': 'Flex' })
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
probe-set,/data/Probe_Set_v1.0_GRCh38-2020-A.csv
no-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_Flex,{fastq_dir},any,PJB1,Gene Expression,

[samples]
sample_id,probe_barcode_ids,description
PB1,BC001,PB1
PB2,BC002,PB2
""".format(fastq_dir=fastq_dir))
        # QC protocol
        protocol = QCProtocol(name="cellranger_multi",
                              description="Cellranger_multi test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_multi",))
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
                "qc/cellranger_multi",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB1/web_summary.html",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB1/metrics_summary.csv",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB2/web_summary.html",
                "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB2/metrics_summary.csv",
                "cellranger_multi",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB1/web_summary.html",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB1/metrics_summary.csv",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB2/web_summary.html",
                "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB2/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,3164)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_multi")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_Flex")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_Flex_S1_R1_001.fastq.gz,"
                         "PJB1_Flex_S1_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-gex-GRCh38-2024-A")
        self.assertEqual(qc_info.cellranger_probeset,
                         "/data/Probe_Set_v1.0_GRCh38-2020-A.csv")
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_multi_flex_900(self):
        """
        QCPipeline: 'cellranger_multi' QC module (Flex, Cellranger v9.0.0)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="flex",
                                 version="9.0.0")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x Flex analysis project
        p = MockAnalysisProject("PJB",("PJB1_Flex_S1_R1_001.fastq.gz",
                                       "PJB1_Flex_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 3\'v3',
                                           'Library type': 'Flex' })
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
probe-set,/data/Probe_Set_v1.0_GRCh38-2020-A.csv
no-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_Flex,{fastq_dir},any,PJB1,Gene Expression,

[samples]
sample_id,probe_barcode_ids,description
PB1,BC001,PB1
PB2,BC002,PB2
""".format(fastq_dir=fastq_dir))
        # QC protocol
        protocol = QCProtocol(name="cellranger_multi",
                              description="Cellranger_multi test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_multi",))
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
                "qc/cellranger_multi",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB1/web_summary.html",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB1/metrics_summary.csv",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB2/web_summary.html",
                "qc/cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB2/metrics_summary.csv",
                "cellranger_multi",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/_cmdline",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB1/web_summary.html",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB1/metrics_summary.csv",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB2/web_summary.html",
                "cellranger_multi/9.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB2/metrics_summary.csv"):
            self.assertTrue(os.path.exists(os.path.join(project_dir,f)),
                            "%s: missing" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,3142)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_multi")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_Flex")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_Flex_S1_R1_001.fastq.gz,"
                         "PJB1_Flex_S1_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,"9.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-gex-GRCh38-2024-A")
        self.assertEqual(qc_info.cellranger_probeset,
                         "/data/Probe_Set_v1.0_GRCh38-2020-A.csv")
        # Check output and reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_cellranger_multi_flex_no_config_file(self):
        """
        QCPipeline: 'cellranger_multi' QC module (Flex, no config file)
        """
        # Make mock QC executables
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="flex",
                                 version="8.0.0")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x Flex analysis project
        p = MockAnalysisProject("PJB",("PJB1_Flex_S1_R1_001.fastq.gz",
                                       "PJB1_Flex_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 3\'v3',
                                           'Library type': 'Flex' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="cellranger_multi",
                              description="Cellranger_multi test",
                              seq_data_reads=['r2',],
                              index_reads=['r1'],
                              qc_modules=("cellranger_multi",))
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
                "qc/cellranger_multi",
                "cellranger_multi"):
            self.assertFalse(os.path.exists(os.path.join(project_dir,f)),
                             "%s: present" % f)
        # Check number of cells
        project_metadata = AnalysisProjectInfo(
            os.path.join(self.wd,"PJB","README.info"))
        self.assertEqual(project_metadata.number_of_cells,None)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"cellranger_multi")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_Flex")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_Flex_S1_R1_001.fastq.gz,"
                         "PJB1_Flex_S1_R2_001.fastq.gz")
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

class TestExpectedOutputs(unittest.TestCase):
    """
    Tests for the 'expected_outputs' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestCheckCellrangerCountOutputs')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_expected_outputs(self):
        """
        expected_outputs: 10x multi config file (multiplexed samples)
        """
        # Make a 10x_multi_config.csv file
        cf_file = os.path.join(self.wd, "10x_multi_config.csv")
        with open(cf_file, "wt") as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2024-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/path/to/fastqs,any,PJB1,gene expression,
PJB1_CML,/path/to/fastqs,any,PJB1,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBA,CMO301,PBA
PBB,CMO302,PBB
""")
        # Generate expected outputs
        self.assertEqual(expected_outputs(cf_file),
                         ["_cmdline",
                          "outs/per_sample_outs/PBA/web_summary.html",
                          "outs/per_sample_outs/PBA/metrics_summary.csv",
                          "outs/per_sample_outs/PBB/web_summary.html",
                          "outs/per_sample_outs/PBB/metrics_summary.csv"])

    def test_expected_outputs_with_prefix(self):
        """
        expected_outputs: 10x multi config file (multiplexed samples, with prefix)
        """
        # Make a 10x_multi_config.csv file
        cf_file = os.path.join(self.wd, "10x_multi_config.csv")
        with open(cf_file, "wt") as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2024-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/path/to/fastqs,any,PJB1,gene expression,
PJB1_CML,/path/to/fastqs,any,PJB1,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBA,CMO301,PBA
PBB,CMO302,PBB
""")
        # Generate expected outputs
        self.assertEqual(expected_outputs(cf_file, prefix="PJB1"),
                         ["PJB1/_cmdline",
                          "PJB1/outs/per_sample_outs/PBA/web_summary.html",
                          "PJB1/outs/per_sample_outs/PBA/metrics_summary.csv",
                          "PJB1/outs/per_sample_outs/PBB/web_summary.html",
                          "PJB1/outs/per_sample_outs/PBB/metrics_summary.csv"])

    def test_expected_outputs_no_multiplexed_samples(self):
        """
        expected_outputs: 10x multi config file (no multiplexed samples)
        """
        # Make a 10x_multi_config.csv file
        cf_file = os.path.join(self.wd, "10x_multi_config.csv")
        with open(cf_file, "wt") as fp:
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2024-A

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,/path/to/fastqs,any,PJB1,gene expression,
PJB1_CML,/path/to/fastqs,any,PJB1,Multiplexing Capture,
""")
        # Generate expected outputs
        self.assertEqual(expected_outputs(cf_file),
                         ["_cmdline"])
