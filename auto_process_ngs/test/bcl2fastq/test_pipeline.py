#######################################################################
# Tests for bcl2fastq/pipeline.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.mock import MockIlluminaRun
from bcftbx.mock import SampleSheets
from auto_process_ngs.mock import MockBcl2fastq2Exe
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import make_mock_bcl2fastq2_output
from auto_process_ngs.pipeliner import PipelineParam
from auto_process_ngs.bcl2fastq.pipeline import MakeFastqs
from auto_process_ngs.bcl2fastq.pipeline import PathJoinParam
from auto_process_ngs.bcl2fastq.pipeline import PathExistsParam
from auto_process_ngs.bcl2fastq.pipeline import FunctionParam
from auto_process_ngs.bcl2fastq.pipeline import subset
from auto_process_ngs.bcl2fastq_utils import make_custom_sample_sheet

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
            
    #@unittest.skip("Skipped")
    def test_makefastqs_standard_protocol_bcl2fastq_2_17(self):
        """
        MakeFastqs: standard protocol with bcl2fastq v2.17
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
        MakeFastqs: standard protocol with bcl2fastq v2.20
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
    def test_makefastqs_standard_protocol_create_fastq_for_index_read(self):
        """
        MakeFastqs: standard protocol: create Fastqs for index reads
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
        # Check adapter sequences
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
        # Check adapter sequences
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
        MakeFastqs: standard protocol: set adapter sequences
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
    def test_makefastqs_standard_protocol_set_minimum_trimmed_read_length(self):
        """
        MakeFastqs: standard protocol: set minimum trimmed read length
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
    def test_makefastqs_standard_protocol_set_mask_short_adapter_reads(self):
        """
        MakeFastqs: standard protocol: set mask short adapter reads
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
    def test_makefastqs_standard_protocol_disable_adapter_trimming(self):
        """
        MakeFastqs: standard protocol: disable adapter trimming
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
    def test_makefastqs_standard_protocol_multiple_lanes(self):
        """
        MakeFastqs: standard protocol: multiple lanes
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
        MakeFastqs: standard protocol: multiple lanes, no lane splitting
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
    def test_makefastqs_standard_protocol_rerun_completed_pipeline(self):
        """
        MakeFastqs: standard protocol: rerun completed pipeline
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
                          "2.2.0"))
        self.assertTrue(p.output.acquired_primary_data)
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
                      "processing_qc.html"):
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
7,KL7,KL7,,,N701,GCCAATAT,S502,TCTTTCCC,KL,
8,MN8,MN8,,,N701,GCCAATAT,S503,TCTTTCCC,MN,
""")
        # Create mock bcl2fastq and cellranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger-atac"),
                                 reads=('R1','R2','R3','I1',))
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
                  "processing_qc.html",):
            with open(os.path.join(analysis_dir,f),'wt') as fp:
                fp.write("")
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       protocol="standard",
                       lane_subsets=(
                           subset(lanes=(1,2,3,),protocol="standard"),
                           subset(lanes=(4,),protocol="mirna"),
                           subset(lanes=(5,),protocol="icell8"),
                           subset(lanes=(6,),protocol="icell8_atac"),
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
                          "2.2.0"))
        self.assertTrue(p.output.acquired_primary_data)
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
                      "processing_qc.html"):
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
7,KL7,KL7,,,N701,GCCAATAT,S502,TCTTTCCC,KL,
8,MN8,MN8,,,N701,GCCAATAT,S503,TCTTTCCC,MN,
""")
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger-atac"),
                                 reads=('R1','R2','R3','I1',))
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
                  "processing_qc.html",):
            with open(os.path.join(analysis_dir,f),'wt') as fp:
                fp.write("")
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,
                       protocol="standard",
                       lane_subsets=(
                           subset(lanes=(1,2,3,),protocol="standard"),
                           subset(lanes=(4,),protocol="mirna"),
                           subset(lanes=(5,),protocol="icell8"),
                           subset(lanes=(6,),protocol="icell8_atac"),
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
                          "2.2.0"))
        self.assertTrue(p.output.acquired_primary_data)
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
                      "processing_qc.html"):
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
        MakeFastqs: standard protocol: no demultiplexing
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
    def test_makefastqs_standard_protocol_with_bases_mask(self):
        """
        MakeFastqs: standard protocol: set bases mask
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
    def test_makefastqs_10x_chromium_sc_protocol(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol
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
                          "2.2.0"))
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
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
                  "processing_qc.html",):
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
                          "2.2.0"))
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
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
    def test_makefastqs_10x_atac_protocol(self):
        """
        MakeFastqs: '10x_atac' protocol
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
                                 reads=('R1','R2','R3','I1',))
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
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger-atac"),
                          "cellranger-atac",
                          "2.2.0"))
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
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
                                              "cellranger-atac"),
                                 reads=('R1','R2','R3','I1',))
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
                  "processing_qc.html",):
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
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger-atac"),
                          "cellranger-atac",
                          "2.2.0"))
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
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
    def test_makefastqs_mirna_protocol(self):
        """
        MakeFastqs: 'mirna' protocol
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
    def test_makefastqs_standard_protocol_missing_fastqs_no_placeholders(self):
        """
        MakeFastqs: standard protocol: missing fastqs, no placeholders
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
        MakeFastqs: standard protocol: missing fastqs with placeholders
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
        MakeFastqs: standard protocol: fails for missing bcl2fastq
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
        MakeFastqs: standard protocol: handle bcl2fastq2 failure
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
        self.assertEqual(status,1)
        # Check outputs
        self.assertEqual(p.output.platform,None)
        self.assertEqual(p.output.primary_data_dir,
                         os.path.join(analysis_dir,
                                      "primary_data"))
        self.assertEqual(p.output.bcl2fastq_info,
                         (os.path.join(self.bin,"bcl2fastq"),
                          "bcl2fastq",
                          "2.20.0.422"))
        self.assertEqual(p.output.cellranger_info,None)
        self.assertTrue(p.output.acquired_primary_data)
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_UNKNOWN_00002_AHGXXXX"),):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        self.assertTrue(os.path.islink(
            os.path.join(analysis_dir,
                         "primary_data",
                         "171020_UNKNOWN_00002_AHGXXXX")))
        for subdir in ("bcl2fastq",
                       "barcode_analysis",):
            self.assertFalse(os.path.exists(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "processing_qc.html"):
            self.assertFalse(os.path.exists(
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
    def test_makefastqs_standard_protocol_fails_for_chromium_sc_indices(self):
        """
        MakeFastqs: standard protocol: fails for Chromium SC indices
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
smpl3,smpl3,,,A006,SI-GA-C1,10xGenomics,
smpl4,smpl4,,,A007,SI-GA-D1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_sc_indices)
        # Make an (empty) analysis directory
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="standard")
        status = p.run(analysis_dir,
                       poll_interval=0.5)
        self.assertEqual(status,1)

class TestPathJoinParam(unittest.TestCase):
    """
    Tests for the 'PathJoinParam' class
    """
    def test_pathjoinparam_with_strings(self):
        """
        PathJoinParam: handle static strings
        """
        p1 = "/mnt/data"
        p2 = "test.txt"
        pth = PathJoinParam(p1,p2)
        self.assertEqual(pth.value,"/mnt/data/test.txt")

    def test_pathjoinparam_with_pipelineparams(self):
        """
        PathJoinParam: handle PipelineParams
        """
        p1 = PipelineParam(value="/mnt/data")
        p2 = PipelineParam(value="test.txt")
        pth = PathJoinParam(p1,p2)
        self.assertEqual(pth.value,"/mnt/data/test.txt")
        p2.set("updated.txt")
        self.assertEqual(pth.value,"/mnt/data/updated.txt")

    def test_pathjoinparam_with_mixed_input(self):
        """
        PathJoinParam: handle mixture of string and PipelineParam
        """
        p1 = PipelineParam(value="/mnt/data")
        p2 = "test.txt"
        pth = PathJoinParam(p1,p2)
        self.assertEqual(pth.value,"/mnt/data/test.txt")

class TestPathExistsParam(unittest.TestCase):
    """
    Tests for the 'PathExistsParam' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestPathExistsParam')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)
    
    def test_pathexistsparam_with_strings(self):
        """
        PathExistsParam: handle static strings
        """
        dir_path = self.wd
        file_path = os.path.join(self.wd,"file.txt")
        with open(file_path,'wt') as fp:
            fp.write("test\n")
        missing_path = os.path.join(self.wd,"missing.txt")
        self.assertTrue(PathExistsParam(dir_path).value)
        self.assertTrue(PathExistsParam(file_path).value)
        self.assertFalse(PathExistsParam(missing_path).value)

    def test_pathexistsparam_with_pipelineparams(self):
        """
        PathExistsParam: handle PipelineParams
        """
        dir_path = PipelineParam(value=self.wd)
        file_path = PipelineParam(value=os.path.join(self.wd,"file.txt"))
        with open(file_path.value,'wt') as fp:
            fp.write("test\n")
        missing_path = PipelineParam(
            value=os.path.join(self.wd,"missing.txt"))
        self.assertTrue(PathExistsParam(dir_path).value)
        self.assertTrue(PathExistsParam(file_path).value)
        self.assertFalse(PathExistsParam(missing_path).value)

    def test_pathexistsparam_with_pathjoinparams(self):
        """
        PathExistsParam: handle PathJoinParams
        """
        dir_path = PathJoinParam(self.wd)
        file_path = PathJoinParam(self.wd,"file.txt")
        with open(file_path.value,'wt') as fp:
            fp.write("test\n")
        missing_path = PathJoinParam(self.wd,"missing.txt")
        self.assertTrue(PathExistsParam(dir_path).value)
        self.assertTrue(PathExistsParam(file_path).value)
        self.assertFalse(PathExistsParam(missing_path).value)

class TestFunctionParam(unittest.TestCase):
    """
    Tests for the 'FunctionParam' class
    """
    #@unittest.skip("Skipped")
    def test_functionparam_with_static_values(self):
        """
        FunctionParam: handle static values
        """
        # Lambda function
        func_param = FunctionParam(lambda x: x,"hello")
        self.assertEqual(func_param.value,"hello")
        func_param = FunctionParam(lambda x,y: x + y,1,2)
        self.assertEqual(func_param.value,3)
        # Function with keywords
        def func(x,y,z=None):
            return (x,y,z)
        func_param = FunctionParam(func,"hello","goodbye")
        self.assertEqual(func_param.value,("hello","goodbye",None))
        func_param = FunctionParam(func,"hello","goodbye",z="surprise!")
        self.assertEqual(func_param.value,("hello","goodbye","surprise!"))

    def test_functionparam_with_parameters(self):
        """
        FunctionParam: handle pipeline parameters
        """
        # Lambda function
        p = PipelineParam(value="hello")
        func_param = FunctionParam(lambda x: x,p)
        self.assertEqual(func_param.value,"hello")
        p.set("goodbye")
        self.assertEqual(func_param.value,"goodbye")
        px = PipelineParam(value=1)
        py = PipelineParam(value=2)
        func_param = FunctionParam(lambda x,y: x + y,px,py)
        self.assertEqual(func_param.value,3)
        px.set(3)
        py.set(4)
        self.assertEqual(func_param.value,7)
        # Function with keywords
        def func(x,y,z=None):
            return (x,y,z)
        px = PipelineParam(value="hello")
        py = PipelineParam(value="goodbye")
        pz = PipelineParam(value="surprise!")
        func_param = FunctionParam(func,px,py)
        self.assertEqual(func_param.value,("hello","goodbye",None))
        func_param = FunctionParam(func,px,py,z=pz)
        self.assertEqual(func_param.value,("hello","goodbye","surprise!"))
        px.set("goodbye")
        py.set("hello")
        pz.set("backwards")
        self.assertEqual(func_param.value,("goodbye","hello","backwards"))

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
        
        
                   
