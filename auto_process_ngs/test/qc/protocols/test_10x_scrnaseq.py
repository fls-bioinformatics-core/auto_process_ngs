#######################################################################
# Unit tests for qc/pipeline.py (10x single cell RNA-seq data)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipeline(BaseQCPipelineTestCase):
    """
    Tests for 10xGenomics single cell RNA-seq data
    """
    def test_qcpipeline_with_10x_sc_rnaseq_data(self):
        """
        QCPipeline: 10xGenomics scRNA-seq run (10x_scRNAseq)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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

    def test_qcpipeline_with_10x_sn_rnaseq_data(self):
        """
        QCPipeline: 10xGenomics snRNA-seq run (10x_snRNAseq)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
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
