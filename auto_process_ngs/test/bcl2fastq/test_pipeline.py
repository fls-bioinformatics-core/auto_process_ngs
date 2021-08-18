#######################################################################
# Tests for bcl2fastq/pipeline.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.mock import MockIlluminaRun
from bcftbx.mock import SampleSheets
from bcftbx.IlluminaData import IlluminaData
from auto_process_ngs.mock import MockBcl2fastq2Exe
from auto_process_ngs.mock import MockBclConvertExe
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import Mock10xPackageExe
from auto_process_ngs.mock import make_mock_bcl2fastq2_output
from auto_process_ngs.bcl2fastq.pipeline import MakeFastqs
from auto_process_ngs.bcl2fastq.pipeline import subset
from auto_process_ngs.bcl2fastq.utils import make_custom_sample_sheet
from auto_process_ngs.stats import FastqStatistics

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestMakeFastqs(unittest.TestCase):
    """
    Tests for MakeFastqs pipeline
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestMakeFastqs')
        # Create a temp 'bin' dir
        self.bin = os.path.join(self.wd,"bin")
        os.mkdir(self.bin)
        # Store original location
        self.pwd = os.getcwd()
        # Store original PATH
        self.path = os.environ['PATH']
        # Move to working dir
        os.chdir(self.wd)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Restore PATH
        os.environ['PATH'] = self.path
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def assertFastqStats(self,stats,fq,expected_nreads,**lanes):
        """
        Internal: check the statistics on a Fastq file

        'stats' should be a FastqStatistics object
        'fq' should be a Fastq file name (no path)
        'expected_nreads' is the total number of expected reads
        'lanes' can be lane identifiers specifying the expected
        read count for that lane e.g. 'L1=4' means 'expect 4
        reads associated with lane 1'
        """
        # Look up the Fastq
        data = stats.raw.lookup('Fastq',fq)
        if not data:
            raise AssertionError("'%s': Fastq not found" % fq)
        # Check number of reads
        self.assertEqual(data[0]['Nreads'],
                         expected_nreads,
                         "'%s': contains %d reads, expected %d"
                         % (fq,data[0]['Nreads'],expected_nreads))
        # Check reads per lane
        for lane in lanes:
            expected_lane_nreads = int(lanes[lane])
            try:
                nreads = data[0][lane]
                if not nreads:
                    nreads = 0
            except KeyError:
                raise AssertionError("'%s': lane '%s' not found"
                                     % (fq,lane))
            self.assertEqual(nreads,
                             expected_lane_nreads,
                             "'%s': has %d reads for lane %s, "
                             "expected %d"
                             % (fq,nreads,lane,expected_lane_nreads))

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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
    def test_makefastqs_standard_protocol_set_out_dirs(self):
        """
        MakeFastqs: standard protocol: set output directories
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
                       out_dir="fastqs_out",
                       barcode_analysis_dir="barcodes",
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "fastqs_out",
                       "barcodes",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertEqual(p.output.missing_fastqs,[])
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
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_set_stats_files(self):
        """
        MakeFastqs: standard protocol: set custom names for stats files
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
        # Custom names for statistics files
        stats_file = os.path.join(self.wd,"the_stats.txt")
        stats_full = os.path.join(self.wd,"the_stats_full.txt")
        per_lane_stats = os.path.join(analysis_dir,
                                      "the_per_lane_stats.txt")
        per_lane_sample_stats = os.path.join(analysis_dir,
                                             "the_per_lane_sample_stats.txt")
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       stats_file=stats_file,
                       stats_full=stats_full,
                       per_lane_stats=os.path.basename(per_lane_stats),
                       per_lane_sample_stats=\
                       os.path.basename(per_lane_sample_stats),
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.stats_file,stats_file)
        self.assertEqual(p.output.stats_full,stats_full)
        self.assertEqual(p.output.per_lane_stats,per_lane_stats)
        self.assertEqual(p.output.per_lane_sample_stats,
                         per_lane_sample_stats)
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertEqual(p.output.missing_fastqs,[])
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in (stats_file,
                      stats_full,
                      per_lane_stats,
                      per_lane_sample_stats,
                      os.path.join(analysis_dir,"processing_qc.html")):
            self.assertTrue(os.path.isfile(filen),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_set_id(self):
        """
        MakeFastqs: standard protocol: set identifier for outputs
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
                       name="test",
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                         os.path.join(analysis_dir,"statistics.test.info"))
        self.assertEqual(p.output.stats_full,
                         os.path.join(analysis_dir,
                                      "statistics_full.test.info"))
        self.assertEqual(p.output.per_lane_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_statistics.test.info"))
        self.assertEqual(p.output.per_lane_sample_stats,
                         os.path.join(analysis_dir,
                                      "per_lane_sample_stats.test.info"))
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis_test",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertEqual(p.output.missing_fastqs,[])
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        for filen in ("statistics.test.info",
                      "statistics_full.test.info",
                      "per_lane_statistics.test.info",
                      "per_lane_sample_stats.test.info",
                      "processing_qc_test.html"):
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_find_adapters_with_sliding_window(self):
        """
        MakeFastqs: standard protocol: use sliding window for adapter trimming
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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

    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_disable_fastq_statistics(self):
        """
        MakeFastqs: standard protocol: disable Fastq statistics
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
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       fastq_statistics=False)
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertFalse(os.path.exists(
            os.path.join(analysis_dir,"barcode_analysis")),
                         "Found subdir: barcode_analysis")
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
    def test_makefastqs_standard_protocol_disable_barcode_analysis(self):
        """
        MakeFastqs: standard protocol: disable barcode analysis
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
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       analyse_barcodes=False)
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertFalse(os.path.exists(
            os.path.join(analysis_dir,"barcode_analysis")),
                         "Found subdir: barcode_analysis")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
                       ##default_runner=SimpleJobRunner(join_logs=True),
                       ##verbose=True)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       bases_mask="y101,I8,I8,y101")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                                 assert_override_cycles="y76;I6n2;I6n2;y76")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,1)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,1)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                  "processing_qc.html",):
            with open(os.path.join(analysis_dir,f),'wt') as fp:
                fp.write("")
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard",
                       analyse_barcodes=False)
        status = p.run(analysis_dir,
                       ##default_runner=SimpleJobRunner(join_logs=True),
                       ##verbose=True,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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

    #@unittest.skip("Skipped")
    def test_makefastqs_multiple_lanes_multiple_subsets(self):
        """
        MakeFastqs: multiple lanes (explicitly define multiple subsets)
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
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       protocol="standard",
                       lane_subsets=(
                           subset(lanes=(1,2,)),
                           subset(lanes=(3,4,5,6,7,8,))))
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "save.bcl2fastq.L12",
                       "save.bcl2fastq.L345678",
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

    #@unittest.skip("Skipped")
    def test_makefastqs_multiple_subsets_rerun_completed_pipeline(self):
        """
        MakeFastqs: multiple subsets: rerun completed pipeline
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
            fp.write("""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,TAAGGCGA,S504,AGAGTAGA,AB,
2,TM2,TM1,,,N701,TAAGGCGA,S517,GCGTAAGA,TM,
3,CD3,CD3,,,N701,GTGAAACG,S503,TCTTTCCC,CD,
4,EB4,EB4,,A1,N701,TAAGGCGA,S501,TAGATCGC,EB,
5,EB5,EB5,,A3,N703,AGGCAGAA,S501,TAGATCGC,EB,
6,EB6,EB6,,F3,N703,AGGCAGAA,S506,ACTGCATA,EB,
7,ML7,ML7,,,N701,GCCAATAT,S502,TCTTTCCC,ML,
8,VL8,VL8,,,N701,GCCAATAT,S503,TCTTTCCC,VL,
""")
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
                                    lanes=(1,2,3,4,5,6,7,8),
                                    sample_sheet=sample_sheet)
        dir2 = os.path.join(self.wd,"analysis.L3456768")
        make_mock_bcl2fastq2_output(os.path.join(dir2,"bcl2fastq"),
                                    lanes=(3,4,5,6,7,8),
                                    sample_sheet=sample_sheet)
        dir1 = os.path.join(self.wd,"analysis.L12")
        make_mock_bcl2fastq2_output(os.path.join(dir1,"bcl2fastq"),
                                    lanes=(1,2),
                                    sample_sheet=sample_sheet)
        os.rename(os.path.join(dir1,"bcl2fastq"),
                  os.path.join(analysis_dir,"save.bcl2fastq.L12"))
        os.rename(os.path.join(dir2,"bcl2fastq"),
                  os.path.join(analysis_dir,"save.bcl2fastq.L345678"))
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
        p = MakeFastqs(run_dir,sample_sheet,
                       protocol="standard",
                       lane_subsets=(
                           subset(lanes=(1,2,)),
                           subset(lanes=(3,4,5,6,7,8,))),
                       analyse_barcodes=False)
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "save.bcl2fastq.L12",
                       "save.bcl2fastq.L345678",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for subdir in ("bcl2fastq.L12",
                       "bcl2fastq.L345678",):
            self.assertFalse(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
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

    #@unittest.skip("Skipped")
    def test_makefastqs_multiple_subsets_rerun_incomplete_pipeline(self):
        """
        MakeFastqs: multiple subsets: rerun incomplete pipeline
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
            fp.write("""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,TAAGGCGA,S504,AGAGTAGA,AB,
2,TM2,TM1,,,N701,TAAGGCGA,S517,GCGTAAGA,TM,
3,CD3,CD3,,,N701,GTGAAACG,S503,TCTTTCCC,CD,
4,EB4,EB4,,A1,N701,TAAGGCGA,S501,TAGATCGC,EB,
5,EB5,EB5,,A3,N703,AGGCAGAA,S501,TAGATCGC,EB,
6,EB6,EB6,,F3,N703,AGGCAGAA,S506,ACTGCATA,EB,
7,ML7,ML7,,,N701,GCCAATAT,S502,TCTTTCCC,ML,
8,VL8,VL8,,,N701,GCCAATAT,S503,TCTTTCCC,VL,
""")
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Set up mock outputs in analysis directory to mimic
        # output from a previous incomplete pipeline run
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        dir1 = os.path.join(self.wd,"analysis.L12")
        make_mock_bcl2fastq2_output(os.path.join(dir1,"bcl2fastq"),
                                    lanes=(1,2),
                                    sample_sheet=sample_sheet)
        dir2 = os.path.join(self.wd,"analysis.L3456768")
        make_mock_bcl2fastq2_output(os.path.join(dir2,"bcl2fastq"),
                                    lanes=(3,4,5,6,7,8),
                                    sample_sheet=sample_sheet)
        os.rename(os.path.join(dir1,"bcl2fastq"),
                  os.path.join(analysis_dir,"save.bcl2fastq.L12"))
        os.rename(os.path.join(dir2,"bcl2fastq"),
                  os.path.join(analysis_dir,"save.bcl2fastq.L345678"))
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
        p = MakeFastqs(run_dir,sample_sheet,
                       protocol="standard",
                       lane_subsets=(
                           subset(lanes=(1,2,)),
                           subset(lanes=(3,4,5,6,7,8,))),
                       analyse_barcodes=False)
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "save.bcl2fastq.L12",
                       "save.bcl2fastq.L345678",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for subdir in ("bcl2fastq.L12",
                       "bcl2fastq.L345678",):
            self.assertFalse(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
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

    #@unittest.skip("Skipped")
    def test_makefastqs_multiple_lanes_multiple_protocols(self):
        """
        MakeFastqs: multiple lanes (explicitly use multiple protocols)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_SN7001250_00002_AHGXXXX",
            "hiseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        samplesheet_multiple_protocols = """[Header]
IEMFileVersion,4
Date,11/11/2015
Workflow,GenerateFASTQ
Application,HiSeq FASTQ Only
Assay,Nextera XT
Description,
Chemistry,Amplicon

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,SI-GA-A1,,,AB,
2,TM2,TM1,,,N701,SI-GA-B1,,,TM,
3,CD3,CD3,,,N701,GTGAAACG,S503,TCTTTCCC,CD,
4,EB4,EB4,,A1,N701,TAAGGCGA,S501,TAGATCGC,EB,
5,EB5,EB5,,A3,N703,AGGCAGAA,S501,TAGATCGC,EB,
6,EB6,EB6,,F3,N703,AGGCAGAA,S506,ACTGCATA,EB,
7,ML7,ML7,,,N701,GCCAATAT,S502,TCTTTCCC,ML,
8,VL8,VL8,,,N701,GCCAATAT,S503,TCTTTCCC,VL,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_multiple_protocols)
        # Create mock bcl2fastq and cellranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       lane_subsets=(
                           subset(lanes=(1,2,),
                                  protocol='10x_chromium_sc'),
                           subset(lanes=(3,4,5,6,),
                                  protocol='standard'),
                           subset(lanes=(7,8,),
                                  protocol='mirna')))
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "6.0.0"))
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
                       "save.bcl2fastq.L12",
                       "save.bcl2fastq.L3456",
                       "save.bcl2fastq.L78",
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_multiple_protocols_rerun_completed_pipeline(self):
        """
        MakeFastqs: multiple protocols: rerun completed pipeline
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
            fp.write("""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,TAAGGCGA,S504,AGAGTAGA,AB,
2,TM2,TM1,,,N701,TAAGGCGA,S517,GCGTAAGA,TM,
3,CD3,CD3,,,N701,GTGAAACG,S503,TCTTTCCC,CD,
4,EF4,EF4,,A1,N701,TAAGGCGA,S501,TAGATCGC,EF,
5,GH5,GH5,,A3,N703,AGGCAGAA,S501,TAGATCGC,GH,
6,IJB6,IJ6,,F3,N703,AGGCAGAA,S506,ACTGCATA,IJ,
7,KL7,KL7,,,N701,SI-GA-A1,S502,,KL,
8,MN8,MN8,,,N701,SI-GA-A1,S503,,MN,
""")
        # Well list file
        well_list = os.path.join(self.wd,"WellList.txt")
        with open(well_list,'wt') as fp:
            fp.write("""Row	Col	Candidate	For dispense	Sample	Barcode	State	Cells1	Cells2	Signal1	Signal2	Size1	Size2	Integ Signal1	Integ Signal2	Circularity1	Circularity2	Confidence	Confidence1	Confidence2	Dispense tip	Drop index	Global drop index	Source well	Sequencing count	Image1	Image2
0	57	True	True	Smpl1	TAACCAAG+TAAGGCGA	Good	1	0	105		33		3465		1		1	1	1	1	10	10	A1		img1.tif	img2.tif
""")
        # Create mock bcl2fastq and cellranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger-atac"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Set up mock outputs in analysis directory to mimic
        # output from a previous successful pipeline run
        analysis_dir = os.path.join(self.wd,"analysis")
        make_mock_bcl2fastq2_output(os.path.join(analysis_dir,
                                                 "bcl2fastq"),
                                    lanes=(1,2,3,4,5,6,7,8),
                                    sample_sheet=sample_sheet)
        for lanes in ((1,2,3),(4,),(5,),(6,),(7,),(8,)):
            lanes_id = ''.join([str(l) for l in lanes])
            dirn = os.path.join(self.wd,"analysis.L%s" % lanes_id)
            make_mock_bcl2fastq2_output(os.path.join(dirn,"bcl2fastq"),
                                        lanes=lanes,
                                        sample_sheet=sample_sheet)
            os.rename(os.path.join(dirn,"bcl2fastq"),
                      os.path.join(analysis_dir,
                                   "save.bcl2fastq.L%s" % lanes_id))
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
                  "processing_qc.html",
                  "cellranger_qc_summary_7.html",
                  "cellranger-atac_qc_summary_8.html",):
            with open(os.path.join(analysis_dir,f),'wt') as fp:
                fp.write("")
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       protocol="standard",
                       lane_subsets=(
                           subset(lanes=(1,2,3,),protocol="standard"),
                           subset(lanes=(4,),protocol="mirna"),
                           subset(lanes=(5,),protocol="icell8"),
                           subset(lanes=(6,),protocol="icell8_atac",
                                  icell8_well_list=well_list),
                           subset(lanes=(7,),protocol="10x_chromium_sc"),
                           subset(lanes=(8,),protocol="10x_atac")),
                       analyse_barcodes=False)
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "6.0.0"))
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
                       "save.bcl2fastq.L123",
                       "save.bcl2fastq.L4",
                       "save.bcl2fastq.L5",
                       "save.bcl2fastq.L6",
                       "save.bcl2fastq.L7",
                       "save.bcl2fastq.L8",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for subdir in ("bcl2fastq.L123",
                       "bcl2fastq.L4",
                       "bcl2fastq.L4",
                       "bcl2fastq.L6",
                       "bcl2fastq.L7",
                       "bcl2fastq.L8",):
            self.assertFalse(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_SN7001250_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html",
                      "cellranger_qc_summary_7.html",
                      "cellranger-atac_qc_summary_8.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_multiple_protocols_rerun_incomplete_pipeline(self):
        """
        MakeFastqs: multiple protocols: rerun incomplete pipeline
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
            fp.write("""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,TAAGGCGA,S504,AGAGTAGA,AB,
2,TM2,TM1,,,N701,TAAGGCGA,S517,GCGTAAGA,TM,
3,CD3,CD3,,,N701,GTGAAACG,S503,TCTTTCCC,CD,
4,EF4,EF4,,A1,N701,TAAGGCGA,S501,TAGATCGC,EF,
5,GH5,GH5,,A3,N703,AGGCAGAA,S501,TAGATCGC,GH,
6,IJB6,IJ6,,F3,N703,AGGCAGAA,S506,ACTGCATA,IJ,
7,KL7,KL7,,,N701,SI-GA-A1,S502,,KL,
8,MN8,MN8,,,N701,SI-GA-A1,S503,,MN,
""")
        # Well list file
        well_list = os.path.join(self.wd,"WellList.txt")
        with open(well_list,'wt') as fp:
            fp.write("""Row	Col	Candidate	For dispense	Sample	Barcode	State	Cells1	Cells2	Signal1	Signal2	Size1	Size2	Integ Signal1	Integ Signal2	Circularity1	Circularity2	Confidence	Confidence1	Confidence2	Dispense tip	Drop index	Global drop index	Source well	Sequencing count	Image1	Image2
0	57	True	True	Smpl1	TAACCAAG+TAAGGCGA	Good	1	0	105		33		3465		1		1	1	1	1	10	10	A1		img1.tif	img2.tif
""")
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger-atac"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Set up mock outputs in analysis directory to mimic
        # output from a previous incomplete pipeline run
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        for lanes in ((1,2,3,),(4,),(5,),(6,),(7,),(8,)):
            lanes_id = ''.join([str(l) for l in lanes])
            dirn = os.path.join(self.wd,"analysis.L%s" % lanes_id)
            make_mock_bcl2fastq2_output(os.path.join(dirn,"bcl2fastq"),
                                        lanes=lanes,
                                        sample_sheet=sample_sheet)
            os.rename(os.path.join(dirn,"bcl2fastq"),
                      os.path.join(analysis_dir,
                                   "save.bcl2fastq.L%s" % lanes_id))
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
                  "processing_qc.html",
                  "cellranger_qc_summary_7.html",
                  "cellranger-atac_qc_summary_8.html",):
            with open(os.path.join(analysis_dir,f),'wt') as fp:
                fp.write("")
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       protocol="standard",
                       lane_subsets=(
                           subset(lanes=(1,2,3,),protocol="standard"),
                           subset(lanes=(4,),protocol="mirna"),
                           subset(lanes=(5,),protocol="icell8"),
                           subset(lanes=(6,),protocol="icell8_atac",
                                  icell8_well_list=well_list),
                           subset(lanes=(7,),protocol="10x_chromium_sc"),
                           subset(lanes=(8,),protocol="10x_atac")),
                       analyse_barcodes=False)
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "6.0.0"))
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
                       "save.bcl2fastq.L123",
                       "save.bcl2fastq.L4",
                       "save.bcl2fastq.L5",
                       "save.bcl2fastq.L6",
                       "save.bcl2fastq.L7",
                       "save.bcl2fastq.L8",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for subdir in ("bcl2fastq.L12",
                       "bcl2fastq.L345678",):
            self.assertFalse(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_SN7001250_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html",
                      "cellranger_qc_summary_7.html",
                      "cellranger-atac_qc_summary_8.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_multiple_lanes_automatically_split_by_indices(self):
        """
        MakeFastqs: multiple lanes (automatically split by indices)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_SN7001250_00002_AHGXXXX",
            "hiseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        samplesheet_multiple_protocols = """[Header]
IEMFileVersion,4
Date,11/11/2015
Workflow,GenerateFASTQ
Application,HiSeq FASTQ Only
Assay,Nextera XT
Description,
Chemistry,Amplicon

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,TAAGGCGA,,,AB,
2,TM2,TM1,,,N701,GCGTAAGA,,,TM,
3,CD3,CD3,,,N701,GTGAAACG,S503,TCTTTCCC,CD,
4,EB4,EB4,,A1,N701,TAAGGCGA,S501,TAGATCGC,EB,
5,EB5,EB5,,A3,N703,AGGCAGAA,S501,TAGATCGC,EB,
6,EB6,EB6,,F3,N703,AGGCAGAA,S506,ACTGCATA,EB,
7,ML7,ML7,,,N701,GCCAATAT,S502,TCTTTCCC,ML,
8,VL8,VL8,,,N701,GCCAATAT,S503,TCTTTCCC,VL,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_multiple_protocols)
        # Create mock bcl2fastq and cellranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "save.bcl2fastq.L12",
                       "save.bcl2fastq.L345678",
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

    #@unittest.skip("Skipped")
    def test_makefastqs_multiple_lanes_multiple_subsets_set_out_dirs(self):
        """
        MakeFastqs: multiple lanes (multiple subsets, set output directories)
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
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       protocol="standard",
                       lane_subsets=(
                           subset(lanes=(1,2,)),
                           subset(lanes=(3,4,5,6,7,8,))))
        status = p.run(analysis_dir,
                       out_dir="fastqs_out",
                       barcode_analysis_dir="barcodes",
                       poll_interval=0.5)
        # Check outputs
        self.assertEqual(p.output.platform,"hiseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "fastqs_out",
                       "save.fastqs_out.L12",
                       "save.fastqs_out.L345678",
                       "barcodes",):
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
                         "171020_SN7001250_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_icell8_protocol(self):
        """
        MakeFastqs: 'icell8' protocol
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        # Check that bases mask, --minimum-trimmed-read-length and
        # --mask-short-adapter-reads are as expected
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_bases_mask="y25n51,I8,y76",
                                 assert_minimum_trimmed_read_length=21,
                                 assert_mask_short_adapter_reads=0)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="icell8")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_icell8_protocol_rerun_completed_pipeline(self):
        """
        MakeFastqs: 'icell8' protocol: rerun completed pipeline
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        # Check that bases mask, --minimum-trimmed-read-length and
        # --mask-short-adapter-reads are as expected
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_bases_mask="y25n51,I8,y76",
                                 assert_minimum_trimmed_read_length=21,
                                 assert_mask_short_adapter_reads=0)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Set up mock outputs in analysis directory to mimic
        # output from a previous successful pipeline run
        analysis_dir = os.path.join(self.wd,"analysis")
        updated_sample_sheet = os.path.join(self.wd,
                                            "SampleSheet.updated.csv")
        make_custom_sample_sheet(sample_sheet,updated_sample_sheet)
        make_mock_bcl2fastq2_output(os.path.join(analysis_dir,
                                                 "bcl2fastq"),
                                    lanes=(1,2,3,4),
                                    sample_sheet=updated_sample_sheet)
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
        p = MakeFastqs(run_dir,sample_sheet,protocol="icell8",
                       analyse_barcodes=False)
        status = p.run(analysis_dir,
                       #default_runner=SimpleJobRunner(join_logs=True),
                       #verbose=True,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       "bcl2fastq",):
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_icell8_atac_protocol(self):
        """
        MakeFastqs: 'icell8_atac' protocol
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write("""[Header]
IEMFileVersion,4

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTA,ICELL8_ATAC,
Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,ICELL8_ATAC,
""")
        # Well list file
        well_list = os.path.join(self.wd,"WellList.txt")
        with open(well_list,'wt') as fp:
            fp.write("""Row	Col	Candidate	For dispense	Sample	Barcode	State	Cells1	Cells2	Signal1	Signal2	Size1	Size2	Integ Signal1	Integ Signal2	Circularity1	Circularity2	Confidence	Confidence1	Confidence2	Dispense tip	Drop index	Global drop index	Source well	Sequencing count	Image1	Image2
0	57	True	True	Smpl1	TAACCAAG+TAAGGCGA	Good	1	0	105		33		3465		1		1	1	1	1	10	10	A1		img1.tif	img2.tif
""")
        # Create mock bcl2fastq
        # Check that --create-fastq-for-index-read is set
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_create_fastq_for_index_read=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="icell8_atac",
                       icell8_well_list=well_list)
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_icell8_atac_protocol_rerun_completed_pipeline(self):
        """
        MakeFastqs: 'icell8_atac' protocol: rerun completed pipeline
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write("""[Header]
IEMFileVersion,4

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTA,ICELL8_ATAC,
Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,ICELL8_ATAC,
""")
        # Well list file
        well_list = os.path.join(self.wd,"WellList.txt")
        with open(well_list,'wt') as fp:
            fp.write("""Row	Col	Candidate	For dispense	Sample	Barcode	State	Cells1	Cells2	Signal1	Signal2	Size1	Size2	Integ Signal1	Integ Signal2	Circularity1	Circularity2	Confidence	Confidence1	Confidence2	Dispense tip	Drop index	Global drop index	Source well	Sequencing count	Image1	Image2
0	57	True	True	Smpl1	TAACCAAG+TAAGGCGA	Good	1	0	105		33		3465		1		1	1	1	1	10	10	A1		img1.tif	img2.tif
""")
        # Create mock bcl2fastq
        # Check that --create-fastq-for-index-read is set
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_create_fastq_for_index_read=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Set up mock outputs in analysis directory to mimic
        # output from a previous successful pipeline run
        analysis_dir = os.path.join(self.wd,"analysis")
        sample_sheet_icell8_atac = os.path.join(
            self.wd,"SampleSheet.icell8_atac.csv")
        with open(sample_sheet_icell8_atac,'wt') as fp:
            fp.write("""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Unassigned-Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTA,ICELL8_ATAC,
Unassigned-Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,ICELL8_ATAC,
""")
        make_mock_bcl2fastq2_output(os.path.join(analysis_dir,
                                                 "bcl2fastq"),
                                    lanes=(1,2,3,4),
                                    reads=('R1','R2','I1','I2'),
                                    sample_sheet=sample_sheet_icell8_atac)
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
        p = MakeFastqs(run_dir,sample_sheet,protocol="icell8_atac",
                       icell8_well_list=well_list,
                       analyse_barcodes=False)
        status = p.run(analysis_dir,
                       #default_runner=SimpleJobRunner(join_logs=True),
                       #verbose=True,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       "bcl2fastq",):
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
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_501(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol (Cellranger 5.0.1)
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
                                              "cellranger"),
                                 version="5.0.1")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "5.0.1"))
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
                      "processing_qc.html",
                      "cellranger_qc_summary.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_600(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol (Cellranger 6.0.0)
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
                                              "cellranger"),
                                 version="6.0.0")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "6.0.0"))
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_non_standard_dir_name(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol (non-standard directory name)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        # Rename the run directory
        run_dir = os.path.join(self.wd,"10X_ILLUMINA_RUN")
        os.rename(illumina_run.dirn,run_dir)
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
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"illumina")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "6.0.0"))
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
                                    "10X_ILLUMINA_RUN"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "10X_ILLUMINA_RUN")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_rerun_completed_pipeline(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol: rerun completed pipeline
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
        # Set up mock outputs in analysis directory to mimic
        # output from a previous successful pipeline run
        analysis_dir = os.path.join(self.wd,"analysis")
        make_mock_bcl2fastq2_output(os.path.join(analysis_dir,
                                                 "bcl2fastq"),
                                    lanes=(1,2,3,4),
                                    sample_sheet=sample_sheet,
                                    force_sample_dir=True)
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
                  "processing_qc.html",
                  "cellranger_qc_summary.html",):
            with open(os.path.join(analysis_dir,f),'wt') as fp:
                fp.write("")
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
        status = p.run(analysis_dir,
                       #default_runner=SimpleJobRunner(join_logs=True),
                       #verbose=True,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "6.0.0"))
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_atac_protocol_120(self):
        """
        MakeFastqs: '10x_atac' protocol (Cellranger-ATAC 1.2.0)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with 10xGenomics Chromium SC ATAC-seq indices
        samplesheet_chromium_sc_atac_indices = """[Header]
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
smpl1,smpl1,,,A001,SI-NA-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-NA-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_sc_atac_indices)
        # Create mock bcl2fastq and cellranger-atac
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger-atac"),
                                 version="1.2.0")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_atac")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_atac_info,
                         (os.path.join(self.bin,"cellranger-atac"),
                          "cellranger-atac",
                          "1.2.0"))
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
                      "processing_qc.html",
                      "cellranger-atac_qc_summary.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_atac_protocol_200(self):
        """
        MakeFastqs: '10x_atac' protocol (Cellranger-ATAC 2.0.0)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with 10xGenomics Chromium SC ATAC-seq indices
        samplesheet_chromium_sc_atac_indices = """[Header]
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
smpl1,smpl1,,,A001,SI-NA-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-NA-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_sc_atac_indices)
        # Create mock bcl2fastq and cellranger-atac
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger-atac"),
                                 version="2.0.0")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_atac")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_atac_info,
                         (os.path.join(self.bin,"cellranger-atac"),
                          "cellranger-atac",
                          "2.0.0"))
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_atac_protocol_rerun_completed_pipeline(self):
        """
        MakeFastqs: '10x_atac' protocol: rerun completed pipeline
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with 10xGenomics Chromium SC ATAC-seq indices
        samplesheet_chromium_sc_atac_indices = """[Header]
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
smpl1,smpl1,,,A001,SI-NA-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-NA-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_sc_atac_indices)
        # Create mock bcl2fastq and cellranger-atac
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger-atac"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Set up mock outputs in analysis directory to mimic
        # output from a previous successful pipeline run
        analysis_dir = os.path.join(self.wd,"analysis")
        make_mock_bcl2fastq2_output(os.path.join(analysis_dir,
                                                 "bcl2fastq"),
                                    lanes=(1,2,3,4),
                                    reads=('R1','R2','R3','I1'),
                                    sample_sheet=sample_sheet,
                                    force_sample_dir=True)
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
                  "processing_qc.html",
                  "cellranger-atac_qc_summary.html",):
            with open(os.path.join(analysis_dir,f),'wt') as fp:
                fp.write("")
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_atac")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_atac_info,
                         (os.path.join(self.bin,"cellranger-atac"),
                          "cellranger-atac",
                          "2.0.0"))
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_multiome_protocol_100_gex_data(self):
        """
        MakeFastqs: '10x_multiome' protocol (Cellranger-ARC 1.0.0) (GEX data)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with 10xGenomics Chromium SC indices
        samplesheet_visium_indices = """[Header]
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
            fp.write(samplesheet_visium_indices)
        # Create mock bcl2fastq and spaceranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "cellranger-arc"),
                                 version="1.0.0",
                                 multiome_data="GEX")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_multiome")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_arc_info,
                         (os.path.join(self.bin,"cellranger-arc"),
                          "cellranger-arc",
                          "1.0.0"))
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
                      "processing_qc.html",
                      "cellranger-arc_qc_summary.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_multiome_protocol_100_atac_data(self):
        """
        MakeFastqs: '10x_multiome' protocol (Cellranger-ARC 1.0.0) (ATAC data)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with 10xGenomics Chromium SC indices
        samplesheet_visium_indices = """[Header]
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
            fp.write(samplesheet_visium_indices)
        # Create mock bcl2fastq and spaceranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "cellranger-arc"),
                                 version="1.0.0",
                                 multiome_data="ATAC")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_multiome")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_arc_info,
                         (os.path.join(self.bin,"cellranger-arc"),
                          "cellranger-arc",
                          "1.0.0"))
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
                      "processing_qc.html",
                      "cellranger-arc_qc_summary.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_multiome_protocol_200_gex_data(self):
        """
        MakeFastqs: '10x_multiome' protocol (Cellranger-ARC 2.0.0) (GEX data)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with 10xGenomics indices
        samplesheet_10x_indices = """[Header]
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
            fp.write(samplesheet_10x_indices)
        # Create mock bcl2fastq and spaceranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "cellranger-arc"),
                                 version="2.0.0",
                                 multiome_data="GEX")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_multiome")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_arc_info,
                         (os.path.join(self.bin,"cellranger-arc"),
                          "cellranger-arc",
                          "2.0.0"))
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_multiome_protocol_200_atac_data(self):
        """
        MakeFastqs: '10x_multiome' protocol (Cellranger-ARC 2.0.0) (ATAC data)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with 10xGenomics indices
        samplesheet_10x_indices = """[Header]
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
            fp.write(samplesheet_10x_indices)
        # Create mock bcl2fastq and spaceranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "cellranger-arc"),
                                 version="2.0.0",
                                 multiome_data="ATAC")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_multiome")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_arc_info,
                         (os.path.join(self.bin,"cellranger-arc"),
                          "cellranger-arc",
                          "2.0.0"))
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_visium_protocol(self):
        """
        MakeFastqs: '10x_visium' protocol
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with 10xGenomics Chromium SC indices
        samplesheet_visium_indices = """[Header]
IEMFileVersion,4
Assay,Nextera XT

[Reads]
76
76

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
smpl1,smpl1,,,SI-TT-A1,SI-TT-A1,SI-TT-A1,SI-TT-A1,10xGenomics,
smpl2,smpl2,,,SI-TT-B1,SI-TT-B1,SI-TT-B1,SI-TT-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_visium_indices)
        # Create mock bcl2fastq and spaceranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "spaceranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_visium")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.spaceranger_info,
                         (os.path.join(self.bin,"spaceranger"),
                          "spaceranger",
                          "1.1.0"))
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
                      "processing_qc.html",
                      "spaceranger_qc_summary.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"nextseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,1)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                       poll_interval=0.5)
        self.assertEqual(status,1)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
        status = p.run(analysis_dir,poll_interval=0.5)
        self.assertEqual(status,1)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
                      "processing_qc.html"):
            self.assertFalse(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Found file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_unknown_sequencer(self):
        """
        MakeFastqs: unknown sequencer and no platform specified
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_UNKNOWN_00002_AHGXXXX",
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
                                 platform="miseq")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"illumina")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_UNKNOWN_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_UNKNOWN_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_explicitly_specify_platform(self):
        """
        MakeFastqs: explicitly specify the platform
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_UNKNOWN_00002_AHGXXXX",
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
                                 platform="miseq")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,platform="miseq")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"miseq")
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_UNKNOWN_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_UNKNOWN_00002_AHGXXXX")))
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_exception_for_invalid_protocol(self):
        """
        MakeFastqs: raise exception for invalid protocol
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
        # Supplying an unrecognised protocol should raise
        # an exception
        self.assertRaises(Exception,
                          MakeFastqs,
                          run_dir,
                          sample_sheet,
                          protocol="not_unimplemented")

    #@unittest.skip("Skipped")
    def test_makefastqs_exception_for_invalid_barcodes(self):
        """
        MakeFastqs: raise exception for invalid barcodes
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with special character (backspace)
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write("""[Header],,,,,,,,,
IEMFileVersion,4
Date,11/23/2015
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGNN,,
Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,,
""")
        # Using a samplesheet with invalid barcodes
        # should raise an exception
        self.assertRaises(Exception,
                          MakeFastqs,
                          run_dir,
                          sample_sheet)

    #@unittest.skip("Skipped")
    def test_makefastqs_exception_for_samplesheet_with_invalid_characters(self):
        """
        MakeFastqs: raise exception for invalid characters in sample sheet
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with special character (backspace)
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write("""[Header],,,,,,,,,
IEMFileVersion,4
Date,11/23/2015
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTC,,\b
Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,,
""")
        # Using a samplesheet with invalid characters
        # should raise an exception
        self.assertRaises(Exception,
                          MakeFastqs,
                          run_dir,
                          sample_sheet)

    #@unittest.skip("Skipped")
    def test_makefastqs_exception_for_standard_protocol_with_10x_barcodes(self):
        """
        MakeFastqs: raise exception for 10x barcodes with 'standard' protocol
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

    #@unittest.skip("Skipped")
    def test_makefastqs_exception_for_10x_chromium_sc_protocol_with_non_10x_barcodes(self):
        """
        MakeFastqs: raise exception for non-10x barcodes with '10x_chromium_sc' protocol
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with 10xGenomics Chromium SC indices
        samplesheet_chromium_non_sc_indices = """[Header]
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
smpl1,smpl1,,,A001,TACCTGAC,10xGenomics,
smpl2,smpl2,,,A005,TGACTACC,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_non_sc_indices)
        # Create mock bcl2fastq and cellranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Using a samplesheet with non-10x barcodes should raise
        # an exception
        self.assertRaises(Exception,
                          MakeFastqs,
                          run_dir,
                          sample_sheet,
                          protocol="10x_chromium_sc")

    #@unittest.skip("Skipped")
    def test_makefastqs_exception_for_inconsistent_subsets(self):
        """
        MakeFastqs: raise exception for inconsistent subsets
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
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Defining subsets with overlapping lanes
        # should raise an exception
        self.assertRaises(Exception,
                          MakeFastqs,
                          run_dir,
                          sample_sheet,
                          protocol="standard",
                          lane_subsets=(
                              subset(lanes=(1,2,)),
                              subset(lanes=(3,4,5,6,7,8,)),
                              subset(lanes=(7,8,))))

    #@unittest.skip("Skipped")
    def test_makefastqs_exception_if_project_is_split(self):
        """
        MakeFastqs: raise exception if project is split by subsets
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
            fp.write("""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,N701,TAAGGCGA,S504,AGAGTAGA,AB,
2,TM2,TM1,,,N701,TAAGGCGA,S517,GCGTAAGA,TM,
3,CD3,CD3,,,N701,GTGAAACG,S503,TCTTTCCC,CD,
4,EB4,EB4,,A1,N701,TAAGGCGA,S501,TAGATCGC,EB,
5,EB5,EB5,,A3,N703,AGGCAGAA,S501,TAGATCGC,EB,
6,EB6,EB6,,F3,N703,AGGCAGAA,S506,ACTGCATA,EB,
7,ML7,ML7,,,N701,GCCAATAT,S502,TCTTTCCC,ML,
8,VL8,VL8,,,N701,GCCAATAT,S503,TCTTTCCC,VL,
""")
        # Create mock bcl2fastq
        # Check that bases mask is as expected
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Specifying subsets which split a project
        # should raise an exception
        self.assertRaises(Exception,
                          MakeFastqs,
                          run_dir,
                          sample_sheet,
                          lane_subsets=(
                              subset(lanes=(1,2,3,4)),
                              subset(lanes=(5,6,7,8,))
                          ))

    #@unittest.skip("Skipped")
    def test_makefastqs_force_copy_of_primary_data(self):
        """
        MakeFastqs: force copying of primary data
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
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet)
        status = p.run(analysis_dir,
                       force_copy_of_primary_data=True,
                       poll_interval=0.5)
        self.assertEqual(status,0)
        # Check outputs
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertFalse(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_M00879_00002_AHGXXXX")))
        self.assertTrue(os.path.isdir(
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

    #@unittest.skip("Skipped")
    def test_makefastqs_exception_for_unknown_protocol(self):
        """
        MakeFastqs: raise exception for unknown protocol
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
        # Specifying unrecognised protocol should raise an
        # exception
        self.assertRaises(Exception,
                          MakeFastqs,
                          run_dir,
                          sample_sheet,
                          protocol="undefined_protocol")

    #@unittest.skip("Skipped")
    def test_makefastqs_exception_specifying_lanes_no_lanes_in_samplesheet(self):
        """
        MakeFastqs: raise exception when lanes specified but none in samplesheet
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        run_dir = illumina_run.dirn
        # Sample sheet with no lanes
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write("""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGNN,,
Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,,
""")
        # Using a samplesheet with no lanes while specifying lane
        # numbers should raise an exception
        self.assertRaises(Exception,
                          MakeFastqs,
                          run_dir,
                          sample_sheet,
                          lanes=[1,])
        # Using a samplesheet with no lanes while specifying subsets
        # should raise an exception
        self.assertRaises(Exception,
                          MakeFastqs,
                          run_dir,
                          sample_sheet,
                          lane_subsets=(
                              subset(lanes=[1,]),
                              subset(lanes=[2,3,4])
                          ))

class TestSubsetFunction(unittest.TestCase):
    """
    Tests for 'subset' function
    """
    def test_subset_lanes_only(self):
        """
        subset: create subset with lanes only
        """
        # Basic lane list
        s = subset([1,2,7,8])
        self.assertEqual(s['lanes'],[1,2,7,8])
        # Lanes in random order
        s = subset([8,7,1,2])
        self.assertEqual(s['lanes'],[1,2,7,8])
        # Lanes specified as strings
        s = subset(['1','2','7','8'])
        self.assertEqual(s['lanes'],[1,2,7,8])

    def test_subset_set_parameters(self):
        """
        subset: create subset with parameters
        """
        # Single parameter
        s = subset([1,2],protocol="standard")
        self.assertEqual(s['lanes'],[1,2])
        self.assertEqual(s['protocol'],"standard")
        # Multiple parameters
        s = subset([1,2],
                   protocol="icell8",
                   icell8_well_list="/files/well_list.txt")
        self.assertEqual(s['lanes'],[1,2])
        self.assertEqual(s['protocol'],"icell8")
        self.assertEqual(s['icell8_well_list'],
                         "/files/well_list.txt")

    def test_subset_setting_unrecognised_parameter_raises_keyerror(self):
        """
        subset: setting unrecognised parameter raises KeyError
        """
        with self.assertRaises(KeyError) as ex:
            s = subset([1,2],
                       protocol="standard",
                       missing="not_here")

    def test_subset_accessing_unrecognised_parameter_raises_keyerror(self):
        """
        subset: accessing unrecognised parameter raises KeyError
        """
        s = subset([1,2],protocol="standard")
        with self.assertRaises(KeyError) as ex:
            value = s['missing']
        
        
                   
