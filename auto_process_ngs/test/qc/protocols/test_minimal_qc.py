#######################################################################
# Unit tests for qc/pipeline.py (minimal QC for SE/PE data)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipeline(BaseQCPipelineTestCase):
    """
    Tests for minimal QC for SE and PE data
    """
    def test_qcpipeline_minimal_PE(self):
        """
        QCPipeline: minimal QC run (paired-end data)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("minimal"),
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
        self.assertEqual(qc_info.protocol,"minimal")
        self.assertEqual(qc_info.protocol_specification,
                         "minimal:'Minimal generic QC':"
                         "seq_reads=[r1,r2]:index_reads=[]:"
                         "qc_modules=[fastq_screen,fastqc,"
                         "sequence_lengths]")
        self.assertEqual(qc_info.organism, None)
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
        self.assertEqual(qc_info.star_index, None)
        self.assertEqual(qc_info.annotation_bed, None)
        self.assertEqual(qc_info.annotation_gtf, None)
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

    def test_qcpipeline_minimal_SE(self):
        """
        QCPipeline: minimal QC run (single-end data)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        MockFastQC.create(os.path.join(self.bin,"fastqc"))
        MockMultiQC.create(os.path.join(self.bin,"multiqc"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"))
        p.create(top_dir=self.wd)
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          fetch_protocol_definition("minimal"),
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
        self.assertEqual(qc_info.protocol,"minimal")
        self.assertEqual(qc_info.protocol_specification,
                         "minimal:'Minimal generic QC':"
                         "seq_reads=[r1]:index_reads=[]:"
                         "qc_modules=[fastq_screen,fastqc,"
                         "sequence_lengths]")
        self.assertEqual(qc_info.organism, None)
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index, None)
        self.assertEqual(qc_info.annotation_bed, None)
        self.assertEqual(qc_info.annotation_gtf, None)
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
