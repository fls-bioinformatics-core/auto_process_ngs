#######################################################################
# Unit tests for qc/pipeline.py (10x Flex data)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipeline(BaseQCPipelineTestCase):
    """
    Tests for 10xGenomics Flex data
    """
    #@unittest.skip("Skipped")
    def test_qcpipeline_with_10x_flex_data(self):
        """
        QCPipeline: 10xGenomics Flex run (10x_Flex)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockSeqtk.create(os.path.join(self.bin,"seqtk"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="flex",
                                 version="8.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x Flex analysis project
        p = MockAnalysisProject("PJB",("PJB1_Flex_S1_R1_001.fastq.gz",
                                       "PJB1_Flex_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10x Chromium Flex',
                                           'Library type': 'scRNA-seq' })
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
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_Flex"),
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
        self.assertEqual(qc_info.protocol,"10x_Flex")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_Flex")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_Flex")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_Flex_S1_R1_001.fastq.gz,"
                         "PJB1_Flex_S1_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-gex-GRCh38-2024-A")
        self.assertEqual(qc_info.cellranger_probeset,
                         "/data/Probe_Set_v1.0_GRCh38-2020-A.csv")
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
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
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/outs/per_sample_outs/PB2/metrics_summary.csv",
                  #"qc/cellranger_count",
                  #"qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/_cmdline",
                  #"qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/outs/web_summary.html",
                  #"qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/outs/metrics_summary.csv",
                  #"cellranger_count",
                  #"cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/_cmdline",
                  #"cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/outs/web_summary.html",
                  #"cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_with_10x_flex_data_multiple_samples(self):
        """
        QCPipeline: 10xGenomics Flex run (10x_Flex, multiple samples)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockSeqtk.create(os.path.join(self.bin,"seqtk"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 multi_outputs="flex",
                                 version="8.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x Flex analysis project
        p = MockAnalysisProject("PJB",("PJB1_Flex_S1_R1_001.fastq.gz",
                                       "PJB1_Flex_S1_R2_001.fastq.gz",
                                       "PJB2_Flex_S2_R1_001.fastq.gz",
                                       "PJB2_Flex_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10x Chromium Flex',
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # Add the cellranger multi config.csv files
        with open(os.path.join(self.wd,
                               "PJB",
                               "10x_multi_config.PJB1.csv"),'wt') as fp:
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
        with open(os.path.join(self.wd,
                               "PJB",
                               "10x_multi_config.PJB2.csv"),'wt') as fp:
            fastq_dir = os.path.join(self.wd,
                                     "PJB",
                                     "fastqs")
            fp.write("""[gene-expression]
reference,/data/refdata-cellranger-gex-GRCh38-2024-A
probe-set,/data/Probe_Set_v1.0_GRCh38-2020-A.csv
no-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB2_Flex,{fastq_dir},any,PJB2,Gene Expression,

[samples]
sample_id,probe_barcode_ids,description
PB3,BC001,PB3
PB4,BC002,PB4
""".format(fastq_dir=fastq_dir))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_Flex"),
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
        self.assertEqual(qc_info.protocol,"10x_Flex")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_Flex")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_Flex,PJB2_Flex")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_Flex_S1_R1_001.fastq.gz,"
                         "PJB1_Flex_S1_R2_001.fastq.gz,"
                         "PJB2_Flex_S2_R1_001.fastq.gz,"
                         "PJB2_Flex_S2_R2_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"8.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-gex-GRCh38-2024-A")
        self.assertEqual(qc_info.cellranger_probeset,
                         "/data/Probe_Set_v1.0_GRCh38-2020-A.csv")
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB.zip",
                  "qc/cellranger_multi",
                  "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1/_cmdline",
                  "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1/outs/per_sample_outs/PB1/web_summary.html",
                  "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1/outs/per_sample_outs/PB1/metrics_summary.csv",
                  "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1/outs/per_sample_outs/PB2/web_summary.html",
                  "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1/outs/per_sample_outs/PB2/metrics_summary.csv",
                  "cellranger_multi",
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1/_cmdline",
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1/outs/per_sample_outs/PB1/web_summary.html",
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1/outs/per_sample_outs/PB1/metrics_summary.csv",
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1/outs/per_sample_outs/PB2/web_summary.html",
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1/outs/per_sample_outs/PB2/metrics_summary.csv",
                  "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2/_cmdline",
                  "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2/outs/per_sample_outs/PB3/web_summary.html",
                  "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2/outs/per_sample_outs/PB3/metrics_summary.csv",
                  "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2/outs/per_sample_outs/PB4/web_summary.html",
                  "qc/cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2/outs/per_sample_outs/PB4/metrics_summary.csv",
                  "cellranger_multi",
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2/_cmdline",
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2/outs/per_sample_outs/PB3/web_summary.html",
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2/outs/per_sample_outs/PB3/metrics_summary.csv",
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2/outs/per_sample_outs/PB4/web_summary.html",
                  "cellranger_multi/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB2/outs/per_sample_outs/PB4/metrics_summary.csv",
                  #"qc/cellranger_count",
                  #"qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/_cmdline",
                  #"qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/outs/web_summary.html",
                  #"qc/cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/outs/metrics_summary.csv",
                  #"cellranger_count",
                  #"cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/_cmdline",
                  #"cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/outs/web_summary.html",
                  #"cellranger_count/8.0.0/refdata-cellranger-gex-GRCh38-2024-A/PJB1_Flex/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_with_10x_flex_data_no_multi_config_file(self):
        """
        QCPipeline: 10xGenomics Flex run (10x_Flex, no multi config file)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockPicard.create(os.path.join(self.bin,"picard"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockSeqtk.create(os.path.join(self.bin,"seqtk"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        MockRSeQC.create(os.path.join(self.bin,"geneBody_coverage.py"))
        MockQualimap.create(os.path.join(self.bin,"qualimap"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"),
                                 version="7.1.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock 10x Flex analysis project
        p = MockAnalysisProject("PJB",("PJB1_Flex_S1_R1_001.fastq.gz",
                                       "PJB1_Flex_S1_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           '10x Chromium Flex',
                                           'Library type': 'scRNA-seq' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("10x_Flex"),
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
        self.assertEqual(qc_info.protocol,"10x_Flex")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_Flex")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1_Flex")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_Flex_S1_R1_001.fastq.gz,"
                         "PJB1_Flex_S1_R2_001.fastq.gz")
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
        # Check that unexpected outputs are not present
        for f in ("qc/cellranger_multi",
                  "qc/cellranger_count",):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                         "PJB",f)),
                            "Found %s (shouldn't be present)" % f)

