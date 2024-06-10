#######################################################################
# Unit tests for qc/pipeline.py (10x single cell immune profiling data)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipeline(BaseQCPipelineTestCase):
    """
    Tests for 10xGenomics single cell immune profiling data
    """
    #@unittest.skip("Skipped")
    def test_qcpipeline_10x_immune_profiling_data(self):
        """
        QCPipeline: 10xGenomics SC immune profiling run (10x_ImmuneProfiling)
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
                                 version="8.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x single cell immune profiling project
        p = MockAnalysisProject("PJB",("PJB1_GEX_S1_R1_001.fastq.gz",
                                       "PJB1_GEX_S1_R2_001.fastq.gz",
                                       "PJB1_TCR_S2_R1_001.fastq.gz",
                                       "PJB1_TCR_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 5\'',
                                           'Library type':
                                           'Single Cell Immune Profiling' })
        p.create(top_dir=self.wd)
        # Add the cellranger multi config.csv file
        # (This is only used to identify the GEX samples)
        with open(os.path.join(self.wd,
                               "PJB",
                               "10x_multi_config.PJB1.csv"),'wt') as fp:
            fastq_dir = os.path.join(self.wd,
                                     "PJB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2024-A

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,%s,any,PJB1,gene expression,
PJB1_TCR,%s,any,PJB1,VDJ-T,
""" % (fastq_dir,fastq_dir))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_ImmuneProfiling"),
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
        self.assertEqual(qc_info.protocol,"10x_ImmuneProfiling")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_ImmuneProfiling")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_GEX,PJB1_TCR")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_GEX_S1_R1_001.fastq.gz,"
                         "PJB1_GEX_S1_R2_001.fastq.gz,"
                         "PJB1_TCR_S2_R1_001.fastq.gz,"
                         "PJB1_TCR_S2_R2_001.fastq.gz")
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
        # Verify no cellranger count outputs for VDJ-T samples
        for f in ("qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_TCR",
                  "cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_TCR"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Found %s (shouldn't exist)" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_10x_immune_profiling_data_multiple_samples(self):
        """
        QCPipeline: 10xGenomics SC immune profiling run (10x_ImmuneProfiling, multiple samples)
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
                                 version="7.1.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x single cell immune profiling project
        p = MockAnalysisProject("PJB",("PJB1_GEX_S1_R1_001.fastq.gz",
                                       "PJB1_GEX_S1_R2_001.fastq.gz",
                                       "PJB1_TCR_S2_R1_001.fastq.gz",
                                       "PJB1_TCR_S2_R2_001.fastq.gz",
                                       "PJB2_GEX_S3_R1_001.fastq.gz",
                                       "PJB2_GEX_S3_R2_001.fastq.gz",
                                       "PJB2_TCR_S4_R1_001.fastq.gz",
                                       "PJB2_TCR_S4_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10xGenomics Chromium 5\'',
                                           'Library type':
                                           'Single Cell Immune Profiling' })
        p.create(top_dir=self.wd)
        # Add the cellranger multi config.csv files
        # (This is only used to identify the GEX samples)
        with open(os.path.join(self.wd,
                               "PJB",
                               "10x_multi_config.PJB1.csv"),'wt') as fp:
            fastq_dir = os.path.join(self.wd,
                                     "PJB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB1_GEX,%s,any,PJB1,gene expression,
PJB1_TCR,%s,any,PJB1,VDJ-T,
""" % (fastq_dir,fastq_dir))
        with open(os.path.join(self.wd,
                               "PJB",
                               "10x_multi_config.PJB2.csv"),'wt') as fp:
            fastq_dir = os.path.join(self.wd,
                                     "PJB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2020-A

[vdj]
reference,/data/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB2_GEX,%s,any,PJB2,gene expression,
PJB2_TCR,%s,any,PJB2,VDJ-T,
""" % (fastq_dir,fastq_dir))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_ImmuneProfiling"),
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
        self.assertEqual(qc_info.protocol,"10x_ImmuneProfiling")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_ImmuneProfiling")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,
                         "PJB1_GEX,PJB1_TCR,PJB2_GEX,PJB2_TCR")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_GEX_S1_R1_001.fastq.gz,"
                         "PJB1_GEX_S1_R2_001.fastq.gz,"
                         "PJB1_TCR_S2_R1_001.fastq.gz,"
                         "PJB1_TCR_S2_R2_001.fastq.gz,"
                         "PJB2_GEX_S3_R1_001.fastq.gz,"
                         "PJB2_GEX_S3_R2_001.fastq.gz,"
                         "PJB2_TCR_S4_R1_001.fastq.gz,"
                         "PJB2_TCR_S4_R2_001.fastq.gz")
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
                  "qc/cellranger_count",
                  "qc/cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "qc/cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "qc/cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB1_GEX/outs/metrics_summary.csv",
                  "qc/cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "qc/cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "qc/cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB2_GEX/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB1_GEX/outs/metrics_summary.csv",
                  "cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB2_GEX/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)
        # Verify no cellranger counts outputs for VDJ-T samples
        for f in ("qc/cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB1_TCR",
                  "qc/cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB2_TCR",
                  "cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB1_TCR",
                  "cellranger_count/7.1.0/refdata-cellranger-gex-GRCh38-2020-A/PJB2_TCR"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Found %s (shouldn't exist)" % f)
