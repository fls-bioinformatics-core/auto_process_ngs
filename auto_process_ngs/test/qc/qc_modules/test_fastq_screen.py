#######################################################################
# Unit tests for qc/pipeline.py ('fastq_screen' QC module)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestQCPipelineFastqScreen(BaseQCPipelineTestCase):
    """
    Tests for 'fastq_screen' QC module
    """
    def test_qcpipeline_qc_modules_fastq_screen_pe(self):
        """
        QCPipeline: 'fastq_screen' QC module (PE data)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
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
        protocol = QCProtocol(name="fastq_screen",
                              description="Fastq_screen test",
                              seq_data_reads=['r1','r2'],
                              index_reads=None,
                              qc_modules=("fastq_screen",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_screen_model_organisms.png",
                  "PJB1_S1_R1_001_screen_model_organisms.txt",
                  "PJB1_S1_R1_001_screen_other_organisms.png",
                  "PJB1_S1_R1_001_screen_other_organisms.txt",
                  "PJB1_S1_R1_001_screen_rRNA.png",
                  "PJB1_S1_R1_001_screen_rRNA.txt",
                  "PJB1_S1_R2_001_screen_model_organisms.png",
                  "PJB1_S1_R2_001_screen_model_organisms.txt",
                  "PJB1_S1_R2_001_screen_other_organisms.png",
                  "PJB1_S1_R2_001_screen_other_organisms.txt",
                  "PJB1_S1_R2_001_screen_rRNA.png",
                  "PJB1_S1_R2_001_screen_rRNA.txt",
                  "PJB2_S2_R1_001_screen_model_organisms.png",
                  "PJB2_S2_R1_001_screen_model_organisms.txt",
                  "PJB2_S2_R1_001_screen_other_organisms.png",
                  "PJB2_S2_R1_001_screen_other_organisms.txt",
                  "PJB2_S2_R1_001_screen_rRNA.png",
                  "PJB2_S2_R1_001_screen_rRNA.txt",
                  "PJB2_S2_R2_001_screen_model_organisms.png",
                  "PJB2_S2_R2_001_screen_model_organisms.txt",
                  "PJB2_S2_R2_001_screen_other_organisms.png",
                  "PJB2_S2_R2_001_screen_other_organisms.txt",
                  "PJB2_S2_R2_001_screen_rRNA.png",
                  "PJB2_S2_R2_001_screen_rRNA.txt"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastq_screen")
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,None)
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

    def test_qcpipeline_qc_modules_fastq_screen_pe_with_r1_index_reads(self):
        """
        QCPipeline: 'fastq_screen' QC module (PE data, R1 as index reads)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
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
        protocol = QCProtocol(name="fastq_screen",
                              description="Fastq_screen test",
                              seq_data_reads=['r2'],
                              index_reads=['r1'],
                              qc_modules=("fastq_screen",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R2_001_screen_model_organisms.png",
                  "PJB1_S1_R2_001_screen_model_organisms.txt",
                  "PJB1_S1_R2_001_screen_other_organisms.png",
                  "PJB1_S1_R2_001_screen_other_organisms.txt",
                  "PJB1_S1_R2_001_screen_rRNA.png",
                  "PJB1_S1_R2_001_screen_rRNA.txt",
                  "PJB2_S2_R2_001_screen_model_organisms.png",
                  "PJB2_S2_R2_001_screen_model_organisms.txt",
                  "PJB2_S2_R2_001_screen_other_organisms.png",
                  "PJB2_S2_R2_001_screen_other_organisms.txt",
                  "PJB2_S2_R2_001_screen_rRNA.png",
                  "PJB2_S2_R2_001_screen_rRNA.txt"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        for f in ("PJB1_S1_R1_001_screen_model_organisms.png",
                  "PJB1_S1_R1_001_screen_model_organisms.txt",
                  "PJB1_S1_R1_001_screen_other_organisms.png",
                  "PJB1_S1_R1_001_screen_other_organisms.txt",
                  "PJB1_S1_R1_001_screen_rRNA.png",
                  "PJB1_S1_R1_001_screen_rRNA.txt",
                  "PJB2_S2_R1_001_screen_model_organisms.png",
                  "PJB2_S2_R1_001_screen_model_organisms.txt",
                  "PJB2_S2_R1_001_screen_other_organisms.png",
                  "PJB2_S2_R1_001_screen_other_organisms.txt",
                  "PJB2_S2_R1_001_screen_rRNA.png",
                  "PJB2_S2_R1_001_screen_rRNA.txt"):
            self.assertFalse(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: present" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastq_screen")
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,None)
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

    def test_qcpipeline_qc_modules_fastq_screen_se(self):
        """
        QCPipeline: 'fastq_screen' QC module (SE data)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastq_screen",
                              description="Fastq_screen test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastq_screen",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_screen_model_organisms.png",
                  "PJB1_S1_R1_001_screen_model_organisms.txt",
                  "PJB1_S1_R1_001_screen_other_organisms.png",
                  "PJB1_S1_R1_001_screen_other_organisms.txt",
                  "PJB1_S1_R1_001_screen_rRNA.png",
                  "PJB1_S1_R1_001_screen_rRNA.txt",
                  "PJB2_S2_R1_001_screen_model_organisms.png",
                  "PJB2_S2_R1_001_screen_model_organisms.txt",
                  "PJB2_S2_R1_001_screen_other_organisms.png",
                  "PJB2_S2_R1_001_screen_other_organisms.txt",
                  "PJB2_S2_R1_001_screen_rRNA.png",
                  "PJB2_S2_R1_001_screen_rRNA.txt"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastq_screen")
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,None)
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

    def test_qcpipeline_qc_modules_fastq_screen_se_with_biological_samples(self):
        """
        QCPipeline: 'fastq_screen' QC module (SE data with biological samples)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human',
                                           'Biological samples': 'PJB1' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastq_screen",
                              description="Fastq_screen test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastq_screen",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_screen_model_organisms.png",
                  "PJB1_S1_R1_001_screen_model_organisms.txt",
                  "PJB1_S1_R1_001_screen_other_organisms.png",
                  "PJB1_S1_R1_001_screen_other_organisms.txt",
                  "PJB1_S1_R1_001_screen_rRNA.png",
                  "PJB1_S1_R1_001_screen_rRNA.txt"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastq_screen")
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,None)
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

    def test_qcpipeline_qc_modules_fastq_screen_se_ignore_index_reads(self):
        """
        QCPipeline: 'fastq_screen' QC module (SE data, ignore I* index reads)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_I1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_I1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastq_screen",
                              description="Fastq_screen test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastq_screen",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_screen_model_organisms.png",
                  "PJB1_S1_R1_001_screen_model_organisms.txt",
                  "PJB1_S1_R1_001_screen_other_organisms.png",
                  "PJB1_S1_R1_001_screen_other_organisms.txt",
                  "PJB1_S1_R1_001_screen_rRNA.png",
                  "PJB1_S1_R1_001_screen_rRNA.txt",
                  "PJB2_S2_R1_001_screen_model_organisms.png",
                  "PJB2_S2_R1_001_screen_model_organisms.txt",
                  "PJB2_S2_R1_001_screen_other_organisms.png",
                  "PJB2_S2_R1_001_screen_other_organisms.txt",
                  "PJB2_S2_R1_001_screen_rRNA.png",
                  "PJB2_S2_R1_001_screen_rRNA.txt"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        for f in ("PJB1_S1_I1_001_screen_model_organisms.png",
                  "PJB1_S1_I1_001_screen_model_organisms.txt",
                  "PJB1_S1_I1_001_screen_other_organisms.png",
                  "PJB1_S1_I1_001_screen_other_organisms.txt",
                  "PJB1_S1_I1_001_screen_rRNA.png",
                  "PJB1_S1_I1_001_screen_rRNA.txt",
                  "PJB2_S2_I1_001_screen_model_organisms.png",
                  "PJB2_S2_I1_001_screen_model_organisms.txt",
                  "PJB2_S2_I1_001_screen_other_organisms.png",
                  "PJB2_S2_I1_001_screen_other_organisms.txt",
                  "PJB2_S2_I1_001_screen_rRNA.png",
                  "PJB2_S2_I1_001_screen_rRNA.txt"):
            self.assertFalse(os.path.exists(os.path.join(qc_dir,f)),
                             "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastq_screen")
        self.assertEqual(qc_info.protocol_specification,
                         str(protocol))
        self.assertEqual(qc_info.organism,"Human")
        self.assertEqual(qc_info.seq_data_samples,"PJB1,PJB2")
        self.assertEqual(qc_info.fastq_dir,
                         os.path.join(self.wd,"PJB","fastqs"))
        self.assertEqual(qc_info.fastqs,
                         "PJB1_S1_I1_001.fastq.gz,"
                         "PJB1_S1_R1_001.fastq.gz,"
                         "PJB2_S2_I1_001.fastq.gz,"
                         "PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(qc_info.fastqs_split_by_lane,False)
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,None)
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

    def test_qcpipeline_qc_modules_fastq_screen_se_split_lanes(self):
        """
        QCPipeline: 'fastq_screen' QC module (SE data, split by lane)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastq_screen",
                              description="Fastq_screen test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastq_screen",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol,
                          split_fastqs_by_lane=True)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_L001_R1_001_screen_model_organisms.png",
                  "PJB1_S1_L001_R1_001_screen_model_organisms.txt",
                  "PJB1_S1_L001_R1_001_screen_other_organisms.png",
                  "PJB1_S1_L001_R1_001_screen_other_organisms.txt",
                  "PJB1_S1_L001_R1_001_screen_rRNA.png",
                  "PJB1_S1_L001_R1_001_screen_rRNA.txt",
                  "PJB2_S2_L001_R1_001_screen_model_organisms.png",
                  "PJB2_S2_L001_R1_001_screen_model_organisms.txt",
                  "PJB2_S2_L001_R1_001_screen_other_organisms.png",
                  "PJB2_S2_L001_R1_001_screen_other_organisms.txt",
                  "PJB2_S2_L001_R1_001_screen_rRNA.png",
                  "PJB2_S2_L001_R1_001_screen_rRNA.txt"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastq_screen")
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,None)
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

    def test_qcpipeline_qc_modules_fastq_screen_se_legacy_naming(self):
        """
        QCPipeline: 'fastq_screen' QC module (SE data, legacy naming)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastq_screen",
                              description="Fastq_screen test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastq_screen",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           legacy_screens=True,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_model_organisms_screen.png",
                  "PJB1_S1_R1_001_model_organisms_screen.txt",
                  "PJB1_S1_R1_001_other_organisms_screen.png",
                  "PJB1_S1_R1_001_other_organisms_screen.txt",
                  "PJB1_S1_R1_001_rRNA_screen.png",
                  "PJB1_S1_R1_001_rRNA_screen.txt",
                  "PJB2_S2_R1_001_model_organisms_screen.png",
                  "PJB2_S2_R1_001_model_organisms_screen.txt",
                  "PJB2_S2_R1_001_other_organisms_screen.png",
                  "PJB2_S2_R1_001_other_organisms_screen.txt",
                  "PJB2_S2_R1_001_rRNA_screen.png",
                  "PJB2_S2_R1_001_rRNA_screen.txt"):
            self.assertTrue(os.path.exists(os.path.join(qc_dir,f)),
                            "%s: missing" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastq_screen")
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
        # For legacy screen naming, no screen names are recorded
        self.assertEqual(qc_info.fastq_screens,None)
        self.assertEqual(qc_info.star_index,None)
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

    def test_qcpipeline_qc_modules_fastq_screen_se_no_panels(self):
        """
        QCPipeline: 'fastq_screen' QC module (SE data, no panels)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastq_screen",
                              description="Fastq_screen test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastq_screen",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,0)
        # Check outputs
        qc_dir = os.path.join(self.wd,"PJB","qc")
        for f in ("PJB1_S1_R1_001_screen_model_organisms.png",
                  "PJB1_S1_R1_001_screen_model_organisms.txt",
                  "PJB1_S1_R1_001_screen_other_organisms.png",
                  "PJB1_S1_R1_001_screen_other_organisms.txt",
                  "PJB1_S1_R1_001_screen_rRNA.png",
                  "PJB1_S1_R1_001_screen_rRNA.txt",
                  "PJB2_S2_R1_001_screen_model_organisms.png",
                  "PJB2_S2_R1_001_screen_model_organisms.txt",
                  "PJB2_S2_R1_001_screen_other_organisms.png",
                  "PJB2_S2_R1_001_screen_other_organisms.txt",
                  "PJB2_S2_R1_001_screen_rRNA.png",
                  "PJB2_S2_R1_001_screen_rRNA.txt"):
            self.assertFalse(os.path.exists(os.path.join(qc_dir,f)),
                             "%s: found" % f)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastq_screen")
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

    def test_qcpipeline_qc_modules_fastq_screen_se_missing_output(self):
        """
        QCPipeline: 'fastq_screen' QC module (SE data, missing outputs)
        """
        # Make mock QC executables
        MockFastqScreen.create(os.path.join(self.bin,"fastq_screen"),
                               no_outputs=True,
                               exit_code=1)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz"),
                                metadata={ 'Organism': 'Human' })
        p.create(top_dir=self.wd)
        # QC protocol
        protocol = QCProtocol(name="fastq_screen",
                              description="Fastq_screen test",
                              seq_data_reads=['r1',],
                              index_reads=None,
                              qc_modules=("fastq_screen",))
        # Set up and run the QC
        runqc = QCPipeline()
        runqc.add_project(AnalysisProject(os.path.join(self.wd,"PJB")),
                          protocol)
        status = runqc.run(fastq_screens=self.fastq_screens,
                           poll_interval=POLL_INTERVAL,
                           max_jobs=1,
                           runners={ 'default': SimpleJobRunner(), })
        self.assertEqual(status,1)
        # Check QC metadata
        qc_info = AnalysisProjectQCDirInfo(
            os.path.join(self.wd,"PJB","qc","qc.info"))
        self.assertEqual(qc_info.protocol,"fastq_screen")
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
        self.assertEqual(qc_info.fastq_screens,
                         "model_organisms,other_organisms,rRNA")
        self.assertEqual(qc_info.star_index,None)
        self.assertEqual(qc_info.annotation_bed,None)
        self.assertEqual(qc_info.annotation_gtf,None)
        self.assertEqual(qc_info.cellranger_version,None)
        self.assertEqual(qc_info.cellranger_refdata,None)
        self.assertEqual(qc_info.cellranger_probeset,None)
        # Check output and reports
        self.assertTrue(os.path.exists(os.path.join(self.wd,"PJB","qc")),
                        "Missing 'qc'")
        for f in ("qc_report.html",
                  "qc_report.PJB.zip"):
            self.assertFalse(os.path.exists(os.path.join(self.wd,
                                                        "PJB",f)),
                             "Found %s, shouldn't be present" % f)
