#######################################################################
# Unit tests for qc/pipeline.py ('rseqc_infer_experiment' QC module)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipelineRSeQCInferExperiment(BaseQCPipelineTestCase):
    """
    Tests for '' QC module
    """
    def test_qcpipeline_qc_modules_rseqc_infer_experiment_pe(self):
        """
        QCPipeline: 'rseqc_infer_experiment' QC module (PE data)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="rseqc_infer_experiment",
                              description="RSeQC_infer_experiment test",
                              seq_data_reads=['r1','r2'],
                              index_reads=None,
                              qc_modules=("rseqc_infer_experiment",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        infer_expt_dir = os.path.join(qc_dir,"rseqc_infer_experiment","human")
        for f in ("PJB1_S1_001.infer_experiment.log",
                  "PJB2_S2_001.infer_experiment.log"):
            self.assertTrue(os.path.exists(os.path.join(infer_expt_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"rseqc_infer_experiment")
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
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_rseqc_infer_experiment_pe_with_r1_index_reads(self):
        """
        QCPipeline: 'rseqc_infer_experiment' QC module (PE data, R1 as index reads)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="rseqc_infer_experiment",
                              description="RSeQC_infer_experiment test",
                              seq_data_reads=['r1','r2'],
                              index_reads=None,
                              qc_modules=("rseqc_infer_experiment",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        infer_expt_dir = os.path.join(qc_dir,"rseqc_infer_experiment","human")
        for f in ("PJB1_S1_001.infer_experiment.log",
                  "PJB2_S2_001.infer_experiment.log"):
            self.assertTrue(os.path.exists(os.path.join(infer_expt_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"rseqc_infer_experiment")
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
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_rseqc_infer_experiment_se(self):
        """
        QCPipeline: 'rseqc_infer_experiment' QC module (SE data)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="rseqc_infer_experiment",
                              description="RSeQC_infer_experiment test",
                              seq_data_reads=['r1','r2'],
                              index_reads=None,
                              qc_modules=("rseqc_infer_experiment",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        infer_expt_dir = os.path.join(qc_dir,"rseqc_infer_experiment","human")
        for f in ("PJB1_S1_001.infer_experiment.log",
                  "PJB2_S2_001.infer_experiment.log"):
            self.assertTrue(os.path.exists(os.path.join(infer_expt_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"rseqc_infer_experiment")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_rseqc_infer_experiment_se_with_biological_samples(self):
        """
        QCPipeline: 'rseqc_infer_experiment' QC module (SE data with biological samples)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Biological samples': 'PJB1' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="rseqc_infer_experiment",
                              description="RSeQC_infer_experiment test",
                              seq_data_reads=['r1','r2'],
                              index_reads=None,
                              qc_modules=("rseqc_infer_experiment",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        infer_expt_dir = os.path.join(qc_dir,"rseqc_infer_experiment","human")
        for f in ("PJB1_S1_001.infer_experiment.log",):
            self.assertTrue(os.path.exists(os.path.join(infer_expt_dir,f)),
                            "%s: missing" % f)
        infer_expt_dir = os.path.join(qc_dir,"rseqc_infer_experiment","human")
        for f in ("PJB2_S2_001.infer_experiment.log",):
            self.assertFalse(os.path.exists(os.path.join(infer_expt_dir,f)),
                            "%s: present" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"rseqc_infer_experiment")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_rseqc_infer_experiment_se_split_lanes(self):
        """
        QCPipeline: 'rseqc_infer_experiment' QC module (SE data, split by lane)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="rseqc_infer_experiment",
                              description="RSeQC_infer_experiment test",
                              seq_data_reads=['r1','r2'],
                              index_reads=None,
                              qc_modules=("rseqc_infer_experiment",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol,
                          split_fastqs_by_lane=True)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        infer_expt_dir = os.path.join(qc_dir,"rseqc_infer_experiment","human")
        for f in ("PJB1_S1_L001_001.infer_experiment.log",
                  "PJB2_S2_L001_001.infer_experiment.log"):
            self.assertTrue(os.path.exists(os.path.join(infer_expt_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"rseqc_infer_experiment")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_L001_R1_001.fastq.gz,"
                         "PJB2_S2_L001_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,True)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_rseqc_infer_experiment_missing_star_index(self):
        """
        QCPipeline: 'rseqc_infer_experiment' QC module (no STAR index)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="rseqc_infer_experiment",
                              description="RSeQC_infer_experiment test",
                              seq_data_reads=['r1','r2'],
                              index_reads=None,
                              qc_modules=("rseqc_infer_experiment",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(annotation_gtf_files=
                           { 'human': self.ref_data['hg38']['gtf'] },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        infer_expt_dir = os.path.join(qc_dir,"rseqc_infer_experiment","human")
        for f in ("PJB1_S1_001.infer_experiment.log",
                  "PJB2_S2_001.infer_experiment.log"):
            self.assertFalse(os.path.exists(os.path.join(infer_expt_dir,f)),
                            "%s: present" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"rseqc_infer_experiment")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,self.ref_data['hg38']['gtf'])
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check reports
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                            "Missing %s" % f)

    def test_qcpipeline_qc_modules_rseqc_infer_experiment_missing_gtf_file(self):
        """
        QCPipeline: 'rseqc_infer_experiment' QC module (no GTF file)
        """
        # Make mock QC executables
        MockStar.create(os.path.join(self.bin,"STAR"))
        MockSamtools.create(os.path.join(self.bin,"samtools"))
        MockGtf2bed.create(os.path.join(self.bin,"gtf2bed"))
        MockRSeQC.create(os.path.join(self.bin,"infer_experiment.py"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="rseqc_infer_experiment",
                              description="RSeQC_infer_experiment test",
                              seq_data_reads=['r1','r2'],
                              index_reads=None,
                              qc_modules=("rseqc_infer_experiment",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(star_indexes=
                           { 'human': '/data/hg38/star_index' },
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        infer_expt_dir = os.path.join(qc_dir,"rseqc_infer_experiment","human")
        for f in ("PJB1_S1_001.infer_experiment.log",
                  "PJB2_S2_001.infer_experiment.log"):
            self.assertFalse(os.path.exists(os.path.join(infer_expt_dir,f)),
                            "%s: present" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"rseqc_infer_experiment")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,"/data/hg38/star_index")
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
