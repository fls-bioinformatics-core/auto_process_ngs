#######################################################################
# Unit tests for qc/pipeline.py (10x CellPlex data)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipeline(BaseQCPipelineTestCase):
    """
    Tests for 10xGenomics CellPlex data
    """
    #@unittest.skip("Skipped")
    def test_qcpipeline_with_10x_cellplex_data(self):
        """
        QCPipeline: 10xGenomics Cellplex run (10x_CellPlex)
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
                                 multi_outputs="cellplex",
                                 version="8.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_CellPlex"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-cellranger-gex-GRCh38-2024-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_CellPlex")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_CellPlex")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-gex-GRCh38-2024-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
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
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/multi/multiplexing_analysis/tag_calls_summary.csv",
                  "qc/cellranger_count",
                  "qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/_cmdline",
                  "qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/outs/web_summary.html",
                  "qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/_cmdline",
                  "cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/outs/web_summary.html",
                  "cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_GEX/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
        # Verify missing outputs for multiplex capture samples
        for f in ("PJB2_MC_S2_R2_001_fastq_strand.txt",
                  "PJB2_MC_S2_R2_001_screen_model_organisms.txt",
                  "PJB2_MC_S2_R2_001_screen_other_organisms.txt",
                  "PJB2_MC_S2_R2_001_screen_rRNA.txt",
                  "qualimap-rnaseq/human/PJB2_MC_S2_001",):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                         "PJB",
                                                         "qc",
                                                         f)),
                             "Found %s (shouldn't exist)" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_with_cellplex_data_no_multi_config_file(self):
        """
        QCPipeline: 10xGenomics Cellplex run (10x_CellPlex, no multi config)
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
                                 multi_outputs="cellplex",
                                 version="6.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_CellPlex"),
                          multiqc=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_bed_files=
                           { 'human': self.ref_data['hg38']['bed'] },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           cellranger_transcriptomes=
                           { 'human':
                             '/data/refdata-cellranger-gex-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_CellPlex")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_CellPlex")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"6.0.0")
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
        # Check that unexpected outputs are not present
        for f in ("qc/cellranger_multi",
                  "qc/cellranger_count",):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                         "PJB",f)),
                            "Found %s (shouldn't be present)" % f)
