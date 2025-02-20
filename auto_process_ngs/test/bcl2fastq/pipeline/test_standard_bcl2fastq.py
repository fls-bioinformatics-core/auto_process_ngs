#######################################################################
# Unit tests for bcl2fastq/pipeline.py (standard protocol with bcl2fastq)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestMakeFastqs(BaseMakeFastqsTestCase):
    """
    Tests for MakeFastqs pipeline (standard protocol with bcl2fastq)
    """
    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bcl2fastq_2_17(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: use v2.17
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 version='2.17.1.14')
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.17.1.14"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bcl2fastq_2_20(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: use v2.20
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 version="2.20.0.422")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_create_fastq_for_index_read(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: create Fastqs for index reads
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        # Check that --create-fastq-for-index-read is set
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_create_fastq_for_index_read=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       create_fastq_for_index_read=True,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_set_adapters(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: set adapter sequences
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        # Check adapter sequences
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_adapter="ACGTACGTACGT",
                                 assert_adapter2="TGCATGCATGCA")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       adapter_sequence="ACGTACGTACGT",
                       adapter_sequence_read2="TGCATGCATGCA")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_set_minimum_trimmed_read_length(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: set minimum trimmed read length
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        # Check adapter sequences
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_minimum_trimmed_read_length=26)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       minimum_trimmed_read_length=26)
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_set_mask_short_adapter_reads(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: set mask short adapter reads
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        # Check adapter sequences
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_mask_short_adapter_reads=10)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       mask_short_adapter_reads=10)
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_disable_adapter_trimming(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: disable adapter trimming
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        # Check adapter trimming is disabled
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_minimum_trimmed_read_length=0,
                                 assert_mask_short_adapter_reads=0,
                                 assert_adapter="",
                                 assert_adapter2="")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       trim_adapters=False)
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
                       ##default_runner=SimpleJobRunner(join_logs=True),
                       ##verbose=True)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_set_bases_mask(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: set bases mask
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_bases_mask="y101,I8,I8,y101")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       bases_mask="y101,I8,I8,y101")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_ignore_missing_bcls(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: ignore missing BCLs
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_ignore_missing_bcls=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard")
        status = p.run(analysis_dir,
                       ignore_missing_bcls=True,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_truncate_read_lengths(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: truncate read lengths
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_bases_mask="y26n75,I8,I8,y76n25")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       r1_length=26,r2_length=76)
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_multiple_lanes(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: multiple lanes
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_SN7001250_00002_AHGXXXX",
            "hiseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.hiseq)
        # Create mock bcl2fastq
        # Check that bases mask is as expected
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_SN7001250_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_multiple_lanes_no_lane_splitting(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: multiple lanes, no lane splitting
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_SN7001250_00002_AHGXXXX",
            "hiseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.hiseq)
        # Create mock bcl2fastq
        # Check that --no-lane-splitting is set
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_no_lane_splitting=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard")
        status = p.run(analysis_dir,
                       no_lane_splitting=True,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_SN7001250_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_multiple_lanes_reduced_subset(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: multiple lanes (reduced subset of lanes)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_SN7001250_00002_AHGXXXX",
            "hiseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.hiseq)
        # Create mock bcl2fastq
        # Check that bases mask is as expected
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       lanes=[1,2,])
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertEqual(IlluminaData(analysis_dir,"bcl2fastq").lanes,
                         [1,2])
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_SN7001250_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_rerun_completed_pipeline(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: rerun completed pipeline
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Set up mock outputs in analysis directory to mimic
        # output from a previous successful pipeline run
        analysis_dir = os.path.join(self.wd,"analysis")
        make_mock_bcl2fastq2_output(os.path.join(analysis_dir,
                                                 "bcl2fastq"),
                                    lanes=(1,),
                                    reads=('R1','R2'),
                                    sample_sheet=sample_sheet)
        # Make a primary data directory
        os.mkdir(os.path.join(analysis_dir,"primary_data"))
        os.symlink(run_dir,
                   os.path.join(analysis_dir,
                                "primary_data",
                                os.path.basename(run_dir)))
        # Make empty stats files
        for f in ("statistics_full.info",
                  "statistics.info",
                  "per_lane_statistics.info",
                  "per_lane_sample_stats.info",
                  "seq_len_statistics.info",
                  "processing_qc.html",):
            with open(os.path.join(analysis_dir,f),'wt') as fp:
                fp.write("")
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       analyse_barcodes=False)
        status = p.run(analysis_dir,
                       ##default_runner=SimpleJobRunner(join_logs=True),
                       ##verbose=True,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_no_demultiplexing(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: no demultiplexing
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with no barcodes
        samplesheet_no_demultiplexing = """[Header]
IEMFileVersion,4
Assay,Nextera XT

[Reads]
76
76

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
AB1,AB1,,,,,AB,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_no_demultiplexing)
        # Create mock bcl2fastq
        # Check that bases mask is as expected
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_bases_mask="y76,nnnnnn,y76")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_NB500968_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_require_bcl2fastq_version(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: specify version
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 version='2.17.1.14')
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl2fastq>=2.17")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.17.1.14"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_require_bcl2fastq_version_not_found(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: specified version not found
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 version='2.17.1.14')
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Pipeline fails for unsatisfied version requirement
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl2fastq=>2.20")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,1)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,None)
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,None)
        self.assertEqual(p.output.stats_full,None)
        self.assertEqual(p.output.per_lane_stats,None)
        self.assertEqual(p.output.per_lane_sample_stats,None)
        self.assertEqual(p.output.seq_len_stats,None)
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for subdir in ("bcl2fastq",
                       "barcode_analysis",):
            self.assertFalse(os.path.exists(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertFalse(os.path.exists(
                os.path.join(analysis_dir,filen)),
                            "Found file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_missing_fastqs_no_placeholders(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: missing fastqs, no placeholders
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 missing_fastqs=(
                                     "Sample1_S1_L001_R1_001.fastq.gz",
                                     "Sample1_S1_L001_R2_001.fastq.gz",))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,1)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,None)
        self.assertEqual(p.output.stats_full,None)
        self.assertEqual(p.output.per_lane_stats,None)
        self.assertEqual(p.output.per_lane_sample_stats,None)
        self.assertEqual(p.output.seq_len_stats,None)
        self.assertEqual(p.output.missing_fastqs,
                         ["Sample1/Sample1_S1_L001_R1_001.fastq.gz",
                          "Sample1/Sample1_S1_L001_R2_001.fastq.gz"])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for subdir in ("bcl2fastq",
                       "barcode_analysis",):
            self.assertFalse(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertFalse(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Found file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_missing_fastqs_with_placeholders(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: missing fastqs with placeholders
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 missing_fastqs=(
                                     "Sample1_S1_L001_R1_001.fastq.gz",
                                     "Sample1_S1_L001_R2_001.fastq.gz",))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       create_empty_fastqs=True,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,
                         ["Sample1/Sample1_S1_L001_R1_001.fastq.gz",
                          "Sample1/Sample1_S1_L001_R2_001.fastq.gz"])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
            
    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_fails_for_missing_bcl2fastq(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: fails for missing bcl2fastq
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,1)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,None)
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,None)
        self.assertEqual(p.output.stats_full,None)
        self.assertEqual(p.output.per_lane_stats,None)
        self.assertEqual(p.output.per_lane_sample_stats,None)
        self.assertEqual(p.output.seq_len_stats,None)
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for subdir in ("bcl2fastq",
                       "barcode_analysis",):
            self.assertFalse(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertFalse(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Found file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_handle_bcl2fastq2_failure(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: handle bcl2fastq2 failure
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq which will fail (i.e.
        # return non-zero exit code)
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 exit_code=1)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,poll_interval=POLL_INTERVAL)
        self.assertEqual(status,1)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,None)
        self.assertEqual(p.output.stats_full,None)
        self.assertEqual(p.output.per_lane_stats,None)
        self.assertEqual(p.output.per_lane_sample_stats,None)
        self.assertEqual(p.output.seq_len_stats,None)
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for subdir in ("bcl2fastq",
                       "barcode_analysis",):
            self.assertFalse(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertFalse(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Found file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_novaseq(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: handle NovaSeq data
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "221206_A7001250_00042_AHGXXXX",
            "novaseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.novaseq)
        # Create mock bcl2fastq
        # Check that bases mask is as expected
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,
                       sample_sheet,
                       protocol="standard",
                       platform="novaseq")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"novaseq")
        self.assertEqual(p.output.flow_cell_mode,"SP")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "221206_A7001250_00042_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "221206_A7001250_00042_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_find_adapters_with_sliding_window(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: use sliding window for adapter trimming
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        # Check that --find-adapters-with-sliding-window is set
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_find_adapters_with_sliding_window=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       find_adapters_with_sliding_window=True,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,
                         os.path.join(analysis_dir,"statistics.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,"statistics_full.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.info"))
        self.assertEqual(p.output.seq_len_stats,
                         os.path.join(analysis_dir,
                                      "seq_len_statistics.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "seq_len_statistics.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_exception_for_standard_protocol_with_10x_barcodes(self):
        """
        MakeFastqs: standard protocol/bcl2fastq: raise exception for 10x barcodes
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with 10xGenomics Chromium SC indices
        samplesheet_chromium_sc_indices = """[Header]
IEMFileVersion,4
Assay,Nextera XT

[Reads]
76
76

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
smpl1,smpl1,,,A001,SI-GA-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-GA-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_sc_indices)
        # Create mock bcl2fastq and cellranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Using a samplesheet with 10x barcodes should raise
        # an exception
        self.assertRaises(Exception,
                          MakeFastqs,
                          run_dir,
                          sample_sheet)
