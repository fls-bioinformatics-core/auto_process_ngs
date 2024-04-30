#######################################################################
# Unit tests for bcl2fastq/pipeline.py (standard protocol with bcl-convert)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestMakeFastqs(BaseMakeFastqsTestCase):
    """
    Tests for MakeFastqs pipeline (standard protocol with bcl-convert)
    """
    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_3_7_5(self):
        """
        MakeFastqs: standard protocol/bcl-convert: use v3.7.5
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
        # Create mock bcl-convert
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 version="3.7.5")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertEqual(IlluminaData(analysis_dir,"bcl2fastq").lanes,
                         [1,])
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_set_adapters(self):
        """
        MakeFastqs: standard protocol/bcl-convert: set adapter sequences
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
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 assert_adapter1="ACGTACGTACGT",
                                 assert_adapter2="TGCATGCATGCA")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert",
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
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_set_missing_adapter_read2_from_sample_sheet(self):
        """
        MakeFastqs: standard protocol/bcl-convert: automatically set missing read2 adapter from sample sheet
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
            fp.write("""[Settings]
Adapter,ACGTACGTACGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTA,,
Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,,""")
        # Create mock bcl2fastq
        # Check adapter sequences
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 assert_adapter1="ACGTACGTACGT",
                                 assert_adapter2="ACGTACGTACGT")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert",
                       adapter_sequence="ACGTACGTACGT")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_dont_reset_adapter_read2_from_sample_sheet_v2(self):
        """
        MakeFastqs: standard protocol/bcl-convert: don't reset read2 adapter for sample sheet V2
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
            fp.write("""[Settings]
AdapterRead1,ACGTACGTACGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTA,,
Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,,""")
        # Create mock bcl2fastq
        # Check adapter sequences
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 assert_adapter1="ACGTACGTACGT",
                                 assert_adapter2="")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert",
                       adapter_sequence="ACGTACGTACGT")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_set_minimum_trimmed_read_length(self):
        """
        MakeFastqs: standard protocol/bcl-convert: set minimum trimmed read length
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
        # Check minimum trimmed read length
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 assert_minimum_trimmed_read_length=26)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert",
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
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_set_mask_short_adapter_reads(self):
        """
        MakeFastqs: standard protocol/bcl-convert: set mask short adapter reads
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
        # Create mock bcl-convert
        # Check short read masking
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 assert_mask_short_reads=10)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert",
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
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_disable_adapter_trimming(self):
        """
        MakeFastqs: standard protocol/bcl-convert: disable adapter trimming
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
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 assert_minimum_trimmed_read_length=0,
                                 assert_mask_short_reads=0,
                                 assert_adapter1="",
                                 assert_adapter2="")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert",
                       trim_adapters=False)
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_set_override_cycles(self):
        """
        MakeFastqs: standard protocol/bcl-convert: set bases mask/override cycles
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
        # Check override cycles
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 assert_override_cycles="Y76;I6N2;I6N2;Y76")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert",
                       bases_mask="y76,I6n2,I6n2,y76")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_truncate_read_lengths(self):
        """
        MakeFastqs: standard protocol/bcl-convert: truncate read lengths
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
        # Check override cycles
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 assert_override_cycles="Y26N75;I8;I8;Y76N25")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert",
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
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_require_bclconvert_version(self):
        """
        MakeFastqs: standard protocol/bcl-convert: specify version
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
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 version='3.7.5')
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert>=3.7")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_require_bclconvert_version_not_found(self):
        """
        MakeFastqs: standard protocol/bcl-convert: specified version not found
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
                                              "bcl-convert"),
                                 version='3.7.5')
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Pipeline fails for unsatisfied version requirement
        p = MakeFastqs(run_dir,sample_sheet,
                       bcl_converter="bcl-convert=>3.8")
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
                      "processing_qc.html"):
            self.assertFalse(os.path.exists(
                os.path.join(analysis_dir,filen)),
                            "Found file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_multiple_lanes(self):
        """
        MakeFastqs: standard protocol/bcl-convert: multiple lanes
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
            samplesheet_data = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,,TAAGGCGA,,AGAGTAGA,AB,
2,TM2,TM2,,,,TAAGGCGA,,GCGTAAGA,TM,
3,CD3,CD3,,,,GTGAAACG,,TCTTTCCC,CD,
4,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
4,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
4,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
5,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
5,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
5,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
6,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
6,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
6,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
7,ML7,ML7,,,,GCCAATAT,,TCTTTCCC,ML,
8,VL8,VL8,,,,GCCAATAT,,TCTTTCCC,VL,"""
            fp.write(samplesheet_data)
        # Create mock bcl-convert
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 version="3.7.5",
                                 assert_no_lane_splitting=False)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       bcl_converter="bcl-convert")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertEqual(IlluminaData(analysis_dir,"bcl2fastq").lanes,
                         [1,2,3,4,5,6,7,8])
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_SN7001250_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1','L2','L3','L4',
                                           'L5','L6','L7','L8'])
        for fq in ("AB1_S1_L001_R1_001.fastq.gz",
                   "AB1_S1_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)
        for fq in ("TM2_S2_L002_R1_001.fastq.gz",
                   "TM2_S2_L002_R2_001.fastq.gz",
                   "Undetermined_S0_L002_R1_001.fastq.gz",
                   "Undetermined_S0_L002_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L2=1)
        for fq in ("CD3_S3_L003_R1_001.fastq.gz",
                   "CD3_S3_L003_R2_001.fastq.gz",
                   "Undetermined_S0_L003_R1_001.fastq.gz",
                   "Undetermined_S0_L003_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L3=1)
        for fq in ("EB4_S4_L004_R1_001.fastq.gz",
                   "EB4_S4_L004_R2_001.fastq.gz",
                   "EB5_S5_L004_R1_001.fastq.gz",
                   "EB5_S5_L004_R2_001.fastq.gz",
                   "EB6_S6_L004_R1_001.fastq.gz",
                   "EB6_S6_L004_R2_001.fastq.gz",
                   "Undetermined_S0_L004_R1_001.fastq.gz",
                   "Undetermined_S0_L004_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L4=1)
        for fq in ("EB4_S4_L005_R1_001.fastq.gz",
                   "EB4_S4_L005_R2_001.fastq.gz",
                   "EB5_S5_L005_R1_001.fastq.gz",
                   "EB5_S5_L005_R2_001.fastq.gz",
                   "EB6_S6_L005_R1_001.fastq.gz",
                   "EB6_S6_L005_R2_001.fastq.gz",
                   "Undetermined_S0_L005_R1_001.fastq.gz",
                   "Undetermined_S0_L005_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L5=1)
        for fq in ("EB4_S4_L006_R1_001.fastq.gz",
                   "EB4_S4_L006_R2_001.fastq.gz",
                   "EB5_S5_L006_R1_001.fastq.gz",
                   "EB5_S5_L006_R2_001.fastq.gz",
                   "EB6_S6_L006_R1_001.fastq.gz",
                   "EB6_S6_L006_R2_001.fastq.gz",
                   "Undetermined_S0_L006_R1_001.fastq.gz",
                   "Undetermined_S0_L006_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L6=1)
        for fq in ("ML7_S7_L007_R1_001.fastq.gz",
                   "ML7_S7_L007_R2_001.fastq.gz",
                   "Undetermined_S0_L007_R1_001.fastq.gz",
                   "Undetermined_S0_L007_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L7=1)
        for fq in ("VL8_S8_L008_R1_001.fastq.gz",
                   "VL8_S8_L008_R2_001.fastq.gz",
                   "Undetermined_S0_L008_R1_001.fastq.gz",
                   "Undetermined_S0_L008_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L8=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_multiple_lanes_no_lane_splitting(self):
        """
        MakeFastqs: standard protocol/bcl-convert: multiple lanes, no lane splitting
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
            samplesheet_data = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,,TAAGGCGA,,AGAGTAGA,AB,
2,TM2,TM2,,,,TAAGGCGA,,GCGTAAGA,TM,
3,CD3,CD3,,,,GTGAAACG,,TCTTTCCC,CD,
4,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
4,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
4,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
5,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
5,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
5,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
6,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
6,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
6,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
7,ML7,ML7,,,,GCCAATAT,,TCTTTCCC,ML,
8,VL8,VL8,,,,GCCAATAT,,TCTTTCCC,VL,"""
            fp.write(samplesheet_data)
        # Create mock bcl-convert
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 version="3.7.5",
                                 assert_no_lane_splitting=False)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       bcl_converter="bcl-convert")
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
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertEqual(IlluminaData(analysis_dir,"bcl2fastq").lanes,
                         [None])
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_SN7001250_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1','L2','L3','L4',
                                           'L5','L6','L7','L8'])
        for fq in ("AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)
        for fq in ("TM2_S2_R1_001.fastq.gz",
                   "TM2_S2_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L2=1)
        for fq in ("CD3_S3_R1_001.fastq.gz",
                   "CD3_S3_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L3=1)
        for fq in ("EB4_S4_R1_001.fastq.gz",
                   "EB4_S4_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,3,L4=1,L5=1,L6=1)
        for fq in ("EB5_S5_R1_001.fastq.gz",
                   "EB5_S5_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,3,L4=1,L5=1,L6=1)
        for fq in ("EB6_S6_R1_001.fastq.gz",
                   "EB6_S6_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,3,L4=1,L5=1,L6=1)
        for fq in ("ML7_S7_R1_001.fastq.gz",
                   "ML7_S7_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L7=1)
        for fq in ("VL8_S8_R1_001.fastq.gz",
                   "VL8_S8_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L8=1)
        for fq in ("Undetermined_S0_R1_001.fastq.gz",
                   "Undetermined_S0_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,8,
                                  L1=1,L2=1,L3=1,L4=1,
                                  L5=1,L6=1,L7=1,L8=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_multiple_lanes_subsets(self):
        """
        MakeFastqs: standard protocol/bcl-convert: multiple lanes with subsets
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
            samplesheet_data = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,,TAAGGCGA,,AGAGTAGA,AB,
2,TM2,TM2,,,,TAAGGCGA,,GCGTAAGA,TM,
3,CD3,CD3,,,,GTGAAACG,,TCTTTCCC,CD,
4,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
4,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
4,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
5,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
5,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
5,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
6,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
6,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
6,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
7,ML7,ML7,,,,GCCAATAT,,TCTTTCCC,ML,
8,VL8,VL8,,,,GCCAATAT,,TCTTTCCC,VL,"""
            fp.write(samplesheet_data)
        # Create mock bcl-convert
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 version="3.7.5",
                                 assert_no_lane_splitting=False)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       bcl_converter="bcl-convert",
                       lane_subsets=(
                           subset(lanes=(1,2,)),
                           subset(lanes=(3,4,5,6,7,8,))))
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertEqual(IlluminaData(analysis_dir,"bcl2fastq").lanes,
                         [1,2,3,4,5,6,7,8])
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_SN7001250_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1','L2','L3','L4',
                                           'L5','L6','L7','L8'])
        for fq in ("AB1_S1_L001_R1_001.fastq.gz",
                   "AB1_S1_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)
        for fq in ("TM2_S2_L002_R1_001.fastq.gz",
                   "TM2_S2_L002_R2_001.fastq.gz",
                   "Undetermined_S0_L002_R1_001.fastq.gz",
                   "Undetermined_S0_L002_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L2=1)
        for fq in ("CD3_S1_L003_R1_001.fastq.gz",
                   "CD3_S1_L003_R2_001.fastq.gz",
                   "Undetermined_S0_L003_R1_001.fastq.gz",
                   "Undetermined_S0_L003_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L3=1)
        for fq in ("EB4_S2_L004_R1_001.fastq.gz",
                   "EB4_S2_L004_R2_001.fastq.gz",
                   "EB5_S3_L004_R1_001.fastq.gz",
                   "EB5_S3_L004_R2_001.fastq.gz",
                   "EB6_S4_L004_R1_001.fastq.gz",
                   "EB6_S4_L004_R2_001.fastq.gz",
                   "Undetermined_S0_L004_R1_001.fastq.gz",
                   "Undetermined_S0_L004_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L4=1)
        for fq in ("EB4_S2_L005_R1_001.fastq.gz",
                   "EB4_S2_L005_R2_001.fastq.gz",
                   "EB5_S3_L005_R1_001.fastq.gz",
                   "EB5_S3_L005_R2_001.fastq.gz",
                   "EB6_S4_L005_R1_001.fastq.gz",
                   "EB6_S4_L005_R2_001.fastq.gz",
                   "Undetermined_S0_L005_R1_001.fastq.gz",
                   "Undetermined_S0_L005_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L5=1)
        for fq in ("EB4_S2_L006_R1_001.fastq.gz",
                   "EB4_S2_L006_R2_001.fastq.gz",
                   "EB5_S3_L006_R1_001.fastq.gz",
                   "EB5_S3_L006_R2_001.fastq.gz",
                   "EB6_S4_L006_R1_001.fastq.gz",
                   "EB6_S4_L006_R2_001.fastq.gz",
                   "Undetermined_S0_L006_R1_001.fastq.gz",
                   "Undetermined_S0_L006_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L6=1)
        for fq in ("ML7_S5_L007_R1_001.fastq.gz",
                   "ML7_S5_L007_R2_001.fastq.gz",
                   "Undetermined_S0_L007_R1_001.fastq.gz",
                   "Undetermined_S0_L007_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L7=1)
        for fq in ("VL8_S6_L008_R1_001.fastq.gz",
                   "VL8_S6_L008_R2_001.fastq.gz",
                   "Undetermined_S0_L008_R1_001.fastq.gz",
                   "Undetermined_S0_L008_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L8=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_multiple_lanes_subsets_no_lane_splitting(self):
        """
        MakeFastqs: standard protocol/bcl-convert: multiple lanes with subsets, no lane splitting
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
            samplesheet_data = """[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,,TAAGGCGA,,AGAGTAGA,AB,
2,TM2,TM2,,,,TAAGGCGA,,GCGTAAGA,TM,
3,CD3,CD3,,,,GTGAAACG,,TCTTTCCC,CD,
4,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
4,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
4,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
5,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
5,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
5,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
6,EB4,EB4,,,,TAAGGCGA,,TAGATCGC,EB,
6,EB5,EB5,,,,AGGCAGAA,,TAGATCTA,EB,
6,EB6,EB6,,,,AGGCAGGG,,ACTGCATA,EB,
7,ML7,ML7,,,,GCCAATAT,,TCTTTCCC,ML,
8,VL8,VL8,,,,GCCAATAT,,TCTTTCCC,VL,"""
            fp.write(samplesheet_data)
        # Create mock bcl-convert
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 assert_no_lane_splitting=False)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       bcl_converter="bcl-convert",
                       lane_subsets=(
                           subset(lanes=(1,2,)),
                           subset(lanes=(3,4,5,6,7,8,))))
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
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1','L2','L3','L4',
                                           'L5','L6','L7','L8'])
        for fq in ("AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)
        for fq in ("TM2_S2_R1_001.fastq.gz",
                   "TM2_S2_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L2=1)
        for fq in ("CD3_S1_R1_001.fastq.gz",
                   "CD3_S1_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L3=1)
        for fq in ("EB4_S2_R1_001.fastq.gz",
                   "EB4_S2_R2_001.fastq.gz",
                   "EB5_S3_R1_001.fastq.gz",
                   "EB5_S3_R2_001.fastq.gz",
                   "EB6_S4_R1_001.fastq.gz",
                   "EB6_S4_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,3,L4=1,L5=1,L6=1)
        for fq in ("ML7_S5_R1_001.fastq.gz",
                   "ML7_S5_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L7=1)
        for fq in ("VL8_S6_R1_001.fastq.gz",
                   "VL8_S6_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L8=1)
        for fq in ("Undetermined_S0_R1_001.fastq.gz",
                   "Undetermined_S0_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,8,
                                  L1=1,L2=1,L3=1,L4=1,
                                  L5=1,L6=1,L7=1,L8=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_no_lanes_no_lane_splitting(self):
        """
        MakeFastqs: standard protocol/bcl-convert: no lanes, no lane splitting
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
        # Create mock bcl-convert
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 version="3.7.5",
                                 assert_no_lane_splitting=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       bcl_converter="bcl-convert")
        status = p.run(analysis_dir,
                       no_lane_splitting=True,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.flow_cell_mode,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertEqual(IlluminaData(analysis_dir,"bcl2fastq").lanes,
                         [None])
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_R1_001.fastq.gz",
                   "Sample1_S1_R2_001.fastq.gz",
                   "Sample2_S2_R1_001.fastq.gz",
                   "Sample2_S2_R2_001.fastq.gz",
                   "Undetermined_S0_R1_001.fastq.gz",
                   "Undetermined_S0_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bclconvert_rerun_completed_pipeline(self):
        """
        MakeFastqs: standard protocol/bcl-convert: rerun completed pipeline
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
        # Create mock bcl-convert
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 version="3.7.5")
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
                  "processing_qc.html",):
            with open(os.path.join(analysis_dir,f),'wt') as fp:
                fp.write("")
        # Do the test
        p = MakeFastqs(run_dir,
                       sample_sheet,
                       protocol="standard",
                       bcl_converter="bcl-convert",
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
        self.assertEqual(p.output.bclconvert_info,
                         (os.path.join(self.bin,"bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)
