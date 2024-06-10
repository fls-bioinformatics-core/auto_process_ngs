#######################################################################
# Unit tests for qc/pipeline.py ('cellranger_count' QC module)
#######################################################################

# All imports declared in __init__.py file
from . import *

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
