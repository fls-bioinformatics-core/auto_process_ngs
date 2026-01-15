#######################################################################
# Unit tests for qc/pipeline.py (10x single cell ATAC data)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipeline(BaseQCPipelineTestCase):
    """
    Tests for 10xGenomics single cell ATAC data
    """
    def test_qcpipeline_with_10x_sc_atac_data(self):
        """
        QCPipeline: 10xGenomics single cell ATAC QC run (10x_scATAC)
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
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"))
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
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_scATAC"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_atac_references=
                           { 'human':
                             '/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_scATAC")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_scATAC")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
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
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
