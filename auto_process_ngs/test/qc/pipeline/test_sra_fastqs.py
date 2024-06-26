#######################################################################
# Unit tests for qc/pipeline.py (SRA Fastqs)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipeline(BaseQCPipelineTestCase):
    """
    Tests for SRA Fastqs
    """
    def test_qcpipeline_sra_paired_end(self):
        """
        QCPipeline: SRA paired-end Fastqs
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
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("SRA",("SRR7089001_1.fastq.gz",
                                       "SRR7089001_2.fastq.gz",
                                       "SRR7089002_1.fastq.gz",
                                       "SRR7089002_2.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"SRA")),
                          fetch_protocol_definition("standardPE"),
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
            os.path.join(self.wd,"SRA","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"standardPE")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("standardPE")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,
                         "SRR7089001,SRR7089002")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"SRA","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "SRR7089001_1.fastq.gz,"
                         "SRR7089001_2.fastq.gz,"
                         "SRR7089002_1.fastq.gz,"
                         "SRR7089002_2.fastq.gz")
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
                  "qc_report.SRA.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "SRA",f)),
                            "Missing %s" % f)
        # Check collated Picard insert sizes
        collated_insert_sizes = os.path.join(self.wd,
                                             "SRA",
                                             "qc",
                                             "insert_sizes.human.tsv")
        self.assertTrue(os.path.exists(collated_insert_sizes),
                        "Missing collated insert sizes TSV")
        with open(collated_insert_sizes,'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Bam file	Mean insert size	Standard deviation	Median insert size	Median absolute deviation
SRR7089001.bam	153.754829	69.675347	139	37
SRR7089002.bam	153.754829	69.675347	139	37
""")

    def test_qc_pipeline_sra_single_end(self):
        """
        QCPipeline: SRA single-end Fastqs
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
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("SRA",("SRR7089001.fastq.gz",
                                       "SRR7089002.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"SRA")),
                          fetch_protocol_definition("standardSE"),
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
            os.path.join(self.wd,"SRA","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"standardSE")
        self.assertEqual(qc_info.protocol_specification,
                         str(fetch_protocol_definition("standardSE")))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,
                         "SRR7089001,SRR7089002")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"SRA","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "SRR7089001.fastq.gz,"
                         "SRR7089002.fastq.gz")
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
                  "qc_report.SRA.zip",
                  "multiqc_report.html"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "SRA",f)),
                            "Missing %s" % f)
