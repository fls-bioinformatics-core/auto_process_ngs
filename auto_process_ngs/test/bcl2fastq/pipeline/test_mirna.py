#######################################################################
# Unit tests for bcl2fastq/pipeline.py (mirna protocol)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestMakeFastqs(BaseMakeFastqsTestCase):
    """
    Tests for MakeFastqs pipeline (mirna protocol)
    """
    #@unittest.skip("Skipped")
    def test_makefastqs_mirna_protocol(self):
        """
        MakeFastqs: 'mirna' protocol/bcl2fastq
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y76,I8,I8,y6",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        # Check --minimum-trimmed-read-length and
        # --mask-short-adapter-reads
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_minimum_trimmed_read_length=10,
                                 assert_mask_short_adapter_reads=0)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="mirna")
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
    def test_makefastqs_mirna_protocol_bclconvert(self):
        """
        MakeFastqs: 'mirna' protocol/bcl-convert
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y76,I8,I8,y6",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl-convert
        # Check --minimum-trimmed-read-length and
        # --mask-short-adapter-reads
        MockBclConvertExe.create(os.path.join(self.bin,
                                              "bcl-convert"),
                                 assert_minimum_trimmed_read_length=10,
                                 assert_mask_short_reads=0,
                                 version="3.7.5")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="mirna",
                       bcl_converter="bcl-convert")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
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
        # Check Fastqs
        stats = FastqStatistics(IlluminaData(analysis_dir,"bcl2fastq"))
        self.assertEqual(stats.lane_names,['L1','L2','L3','L4'])
        for fq in ("Sample1_S1_L001_R1_001.fastq.gz",
                   "Sample1_S1_L001_R2_001.fastq.gz",
                   "Sample2_S2_L001_R1_001.fastq.gz",
                   "Sample2_S2_L001_R2_001.fastq.gz",
                   "Undetermined_S0_L001_R1_001.fastq.gz",
                   "Undetermined_S0_L001_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L1=1)
        for fq in ("Sample1_S1_L002_R1_001.fastq.gz",
                   "Sample1_S1_L002_R2_001.fastq.gz",
                   "Sample2_S2_L002_R1_001.fastq.gz",
                   "Sample2_S2_L002_R2_001.fastq.gz",
                   "Undetermined_S0_L002_R1_001.fastq.gz",
                   "Undetermined_S0_L002_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L2=1)
        for fq in ("Sample1_S1_L003_R1_001.fastq.gz",
                   "Sample1_S1_L003_R2_001.fastq.gz",
                   "Sample2_S2_L003_R1_001.fastq.gz",
                   "Sample2_S2_L003_R2_001.fastq.gz",
                   "Undetermined_S0_L003_R1_001.fastq.gz",
                   "Undetermined_S0_L003_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L3=1)
        for fq in ("Sample1_S1_L004_R1_001.fastq.gz",
                   "Sample1_S1_L004_R2_001.fastq.gz",
                   "Sample2_S2_L004_R1_001.fastq.gz",
                   "Sample2_S2_L004_R2_001.fastq.gz",
                   "Undetermined_S0_L004_R1_001.fastq.gz",
                   "Undetermined_S0_L004_R2_001.fastq.gz",):
            self.assertFastqStats(stats,fq,1,L4=1)
