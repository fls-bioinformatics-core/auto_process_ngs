#######################################################################
# Unit tests for qc/pipeline.py (10x single cell multiome data)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipeline(BaseQCPipelineTestCase):
    """
    Tests for 10xGenomics single cell multiome data
    """
    #@unittest.skip("Skipped")
    def test_qcpipeline_10x_multiome_atac_data_no_linked_gex(self):
        """
        QCPipeline: 10xGenomics single cell multiome ATAC run (10x_Multiome_ATAC, no linked GEX)
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
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"),
                                 version="2.0.0",
                                 assert_chemistry="ARC-v1")
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_ATAC")),
                          fetch_protocol_definition("10x_Multiome_ATAC"),
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
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_ATAC","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_Multiome_ATAC")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_Multiome_ATAC")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"2.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB_ATAC.zip",
                  "qc/cellranger-atac_count",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/_cmdline",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/outs/web_summary.html",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/outs/summary.csv",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/_cmdline",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/outs/web_summary.html",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/outs/summary.csv",
                  "cellranger-atac_count",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/_cmdline",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/outs/web_summary.html",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/outs/summary.csv",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/_cmdline",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/outs/web_summary.html",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_ATAC",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_10x_multiome_gex_data_no_linked_atac(self):
        """
        QCPipeline: 10xGenomics single cell multiome GEX run (10x_Multiome_GEX, no linked ATAC)
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
                                 version="6.0.0",
                                 assert_include_introns=True,
                                 assert_chemistry="ARC-v1")
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_GEX")),
                          fetch_protocol_definition("10x_Multiome_GEX"),
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
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_GEX","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_Multiome_GEX")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_Multiome_GEX")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"2.1.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB_GEX.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/outs/metrics_summary.csv",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/outs/metrics_summary.csv",
                  "cellranger_count",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/outs/metrics_summary.csv",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/outs/metrics_summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_GEX",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_10x_multiome_atac_data_with_linked_gex(self):
        """
        QCPipeline: 10xGenomics single cell multiome ATAC run (10x_Multiome_ATAC, linked GEX)
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
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-atac"),
                                 version="2.0.0",
                                 assert_chemistry="ARC-v1")
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="2.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_ATAC")),
                          fetch_protocol_definition("10x_Multiome_ATAC"),
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
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_ATAC","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_Multiome_ATAC")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_Multiome_ATAC")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB_ATAC.zip",
                  "qc/cellranger-atac_count",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/_cmdline",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/outs/web_summary.html",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/outs/summary.csv",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/_cmdline",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/outs/web_summary.html",
                  "qc/cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/outs/summary.csv",
                  "qc/cellranger-arc_count",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",
                  "cellranger-atac_count",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/_cmdline",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/outs/web_summary.html",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB1_ATAC/outs/summary.csv",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/_cmdline",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/outs/web_summary.html",
                  "cellranger-atac_count/2.0.0/refdata-cellranger-atac-GRCh38-2020-A-2.0.0/PJB2_ATAC/outs/summary.csv",
                  "cellranger-arc_count",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/_cmdline",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/web_summary.html",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_ATAC/outs/summary.csv",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/_cmdline",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/web_summary.html",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_ATAC/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_ATAC",f)),
                            "Missing %s" % f)

    #@unittest.skip("Skipped")
    def test_qcpipeline_10x_multiome_gex_data_with_linked_atac(self):
        """
        QCPipeline: 10xGenomics single cell multiome GEX run (10x_Multiome_ATAC, linked ATAC)
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
                                 version="6.0.0",
                                 assert_include_introns=True,
                                 assert_chemistry="ARC-v1")
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger-arc"),
                                 version="2.0.0")
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
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
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB_GEX")),
                          fetch_protocol_definition("10x_Multiome_GEX"),
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
                           cellranger_arc_references=
                           { 'human':
                             '/data/refdata-cellranger-arc-GRCh38-2020-A' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB_GEX","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"10x_Multiome_GEX")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("10x_Multiome_GEX")))
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,self.ref_data['hg38']['bed'])
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,"2.0.0")
        self.assertEqual(qc_info.cellranger_refdata,
                         "/data/refdata-cellranger-arc-GRCh38-2020-A")
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        for f in ("qc",
                  "qc_report.html",
                  "qc_report.PJB_GEX.zip",
                  "qc/cellranger_count",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/outs/metrics_summary.csv",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "qc/cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/outs/metrics_summary.csv",
                  "qc/cellranger-arc_count",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "qc/cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv",
                  "cellranger_count",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB1_GEX/outs/metrics_summary.csv",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "cellranger_count/6.0.0/refdata-gex-GRCh38-2020-A/PJB2_GEX/outs/metrics_summary.csv",
                  "cellranger-arc_count",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/_cmdline",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/web_summary.html",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB1_GEX/outs/summary.csv",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/_cmdline",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/web_summary.html",
                  "cellranger-arc_count/2.0.0/refdata-cellranger-arc-GRCh38-2020-A/PJB2_GEX/outs/summary.csv",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB_GEX",f)),
                            "Missing %s" % f)
