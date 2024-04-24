#######################################################################
# Unit tests for qc/pipeline.py (10x single cell RNA-seq data)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipeline(BaseQCPipelineTestCase):
    """
    Tests for 10xGenomics single cell RNA-seq data
    """
    def test_qcpipeline_with_cellranger_count_scRNA_seq_310(self):
        """QCPipeline: single cell RNA-seq QC run with 'cellranger count' (v3.1.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_scRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-cellranger-GRCh38-1.2.0' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"3.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-GRCh38-1.2.0")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
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
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_scRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"5.0.1")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
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
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_scRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"6.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
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

    def test_qcpipeline_with_cellranger_count_scRNA_seq_700(self):
        """QCPipeline: single cell RNA-seq QC run with 'cellranger count' (v7.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="7.0.0",
                                 assert_include_introns=True)
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_scRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"7.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
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
                  "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_scRNA_seq_800(self):
        """QCPipeline: single cell RNA-seq QC run with 'cellranger count' (v8.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="8.0.0",
                                 assert_include_introns=True)
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_scRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2024-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2024-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/_cmdline",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/_cmdline",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/_cmdline",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/outs/web_summary.html",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/_cmdline",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/outs/web_summary.html",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/outs/metrics_summary.csv",
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
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_snRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_premrna_references=
                           { 'human': '/data/refdata-cellranger-GRCh38-1.2.0_premrna' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_snRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_snRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"3.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-GRCh38-1.2.0_premrna")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
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
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_snRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_snRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_snRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"5.0.1")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
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
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_snRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_snRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_snRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"6.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
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

    def test_qcpipeline_with_cellranger_count_snRNA_seq_700(self):
        """QCPipeline: single nuclei RNA-seq QC run with 'cellranger count' (v7.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="7.0.0",
                                 assert_include_introns=True)
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_snRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_snRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_snRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"7.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
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
                  "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_snRNA_seq_800(self):
        """QCPipeline: single nuclei RNA-seq QC run with 'cellranger count' (v8.0.0)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="8.0.0",
                                 assert_include_introns=True)
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_snRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2024-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_snRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_snRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2024-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/_cmdline",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/outs/web_summary.html",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/outs/metrics_summary.csv",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/_cmdline",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/outs/web_summary.html",
                  "qc/cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/_cmdline",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/outs/web_summary.html",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB1/outs/metrics_summary.csv",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/_cmdline",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/outs/web_summary.html",
                  "cellranger_count/8.0.0/refdata-gex-GRCh38-2024-A/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_with_cellranger_count_with_extra_project(self):
        """QCPipeline: single cell RNA-seq QC run with extra project
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="7.0.0",
                                 assert_include_introns=True)
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
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_scRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           cellranger_extra_projects=[
                               AnalysisProject("PJB2",
                                               os.path.join(self.wd,"PJB2"))],
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"7.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
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
                  "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB4/outs/metrics_summary.csv",
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
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_scRNAseq"),
                          multiqc=True,
                          organism="human")
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"5.0.1")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
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
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_scRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
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
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
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
                          fetch_protocol_definition("10x_scRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           cellranger_exe=os.path.join(cellranger_bin,
                                                       "cellranger"),
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"5.0.1")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
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

    def test_qcpipeline_with_cellranger_count_force_cells(self):
        """QCPipeline: single cell RNA-seq QC run with 'cellranger count --force-cells'
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockFastqStrandPy.create(os.path.join(self.bin,"fastq_strand.py"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 assert_force_cells=10000)
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_scRNAseq"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human': '/data/refdata-gex-GRCh38-2020-A' },
                           cellranger_force_cells=10000,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scRNAseq")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scRNAseq")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"7.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-gex-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
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
                  "cellranger_count/7.0.0/refdata-gex-GRCh38-2020-A/PJB2/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
