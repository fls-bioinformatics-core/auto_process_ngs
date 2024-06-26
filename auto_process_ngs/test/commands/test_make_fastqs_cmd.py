#######################################################################
# Tests for make_fastqs_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.settings import Settings
from auto_process_ngs.auto_processor import AutoProcess
from bcftbx.mock import MockIlluminaRun
from bcftbx.mock import SampleSheets
from auto_process_ngs.mock import MockBcl2fastq2Exe
from auto_process_ngs.mock import MockBclConvertExe
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import Mock10xPackageExe
from auto_process_ngs.commands.make_fastqs_cmd import make_fastqs
from auto_process_ngs.bcl2fastq.pipeline import subset

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestAutoProcessMakeFastqs(unittest.TestCase):
    """
    Tests for AutoProcess.make_fastqs
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestAutoProcessMakeFastqs')
        # Create settings instance
        # This allows us to set the polling interval for the
        # unit tests
        settings_ini = os.path.join(self.wd,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 0.01
""")
        self.settings = Settings(settings_ini)
        # Create a temp 'bin' dir
        self.bin = os.path.join(self.wd,"bin")
        os.mkdir(self.bin)
        # Store original location
        self.pwd = os.getcwd()
        # Store original PATH
        self.path = os.environ['PATH']
        # Move to working dir
        os.chdir(self.wd)
        # Placeholders for test objects
        self.ap = None

    def tearDown(self):
        # Delete autoprocessor object
        if self.ap is not None:
            del(self.ap)
        # Return to original dir
        os.chdir(self.pwd)
        # Restore PATH
        os.environ['PATH'] = self.path
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    #@unittest.skip("Skipped")
    def test_make_fastqs_standard_protocol(self):
        """make_fastqs: standard protocol
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,"bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_M00879_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="standard")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_M00879_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_M00879_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_standard_protocol_bclconvert(self):
        """make_fastqs: standard protocol (bcl-convert)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl-convert
        MockBclConvertExe.create(os.path.join(self.bin,"bcl-convert"),
                                 version="3.7.5")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_M00879_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="standard",
                    bcl_converter="bcl-convert")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_M00879_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_M00879_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_standard_protocol_override_defaults(self):
        """make_fastqs: standard protocol (override defaults)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,"bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_M00879_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="standard",
                    unaligned_dir="fastqs_out",
                    barcode_analysis_dir="barcodes",
                    stats_file="stats.out")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_M00879_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"fastqs_out")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcodes")
        self.assertEqual(ap.params.stats_file,"stats.out")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_M00879_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs"),
                       "fastqs_out",
                       "barcodes",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("stats.out",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_standard_protocol_specify_id(self):
        """make_fastqs: standard protocol (specify id)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,"bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_M00879_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="standard",name="test")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_M00879_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq_test")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis_test")
        self.assertEqual(ap.params.stats_file,"statistics.test.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_M00879_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_test"),
                       "bcl2fastq_test",
                       "barcode_analysis_test",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.test.info",
                      "statistics_full.test.info",
                      "per_lane_statistics.test.info",
                      "per_lane_sample_stats.test.info",
                      "projects.info",
                      "processing_qc_test.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_standard_protocol_exception_for_chromium_sc_indices(self):
        """make_fastqs: standard protocol with Chromium SC indices raises exception
        """
        # Sample sheet with 10xGenomics Chromium SC indices
        samplesheet_chromium_sc_indices = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
smpl1,smpl1,,,A001,SI-GA-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-GA-B1,10xGenomics,
smpl3,smpl3,,,A006,SI-GA-C1,10xGenomics,
smpl4,smpl4,,,A007,SI-GA-D1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_sc_indices)
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,"bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertRaises(Exception,
                          make_fastqs,
                          ap,
                          protocol="standard")

    #@unittest.skip("Skipped")
    def test_make_fastqs_standard_protocol_stores_bases_mask(self):
        """make_fastqs: standard protocol stores supplied bases mask
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_M00879_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        make_fastqs(ap,protocol="standard",bases_mask="y101,I8,I8,y101")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"y101,I8,I8,y101")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_M00879_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_M00879_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_icell8_protocol(self):
        """make_fastqs: icell8 protocol
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        # Sample sheet
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'wt') as fp:
            fp.write(SampleSheets.miseq)
        # Create mock bcl2fastq
        # Check that bases mask is as expected
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_bases_mask="y25n51,I8,y76")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="icell8")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_NB500968_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_NB500968_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_icell8"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_icell8_atac_protocol(self):
        """make_fastqs: icell8_atac protocol
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
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,
                    protocol="icell8_atac",
                    icell8_well_list=well_list)
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_NB500968_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_NB500968_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_icell8_atac"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        # Check the samples listed in the projects.info file
        with open(os.path.join(analysis_dir,'projects.info'),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
ICELL8_ATAC\tUnassigned-Sample1,Unassigned-Sample2\t.\t.\t.\t.\t.\t.
""")

    #@unittest.skip("Skipped")
    def test_make_fastqs_10x_chromium_sc_protocol(self):
        """make_fastqs: 10x_chromium_sc protocol
        """
        # Sample sheet with 10xGenomics Chromium SC indices
        samplesheet_chromium_sc_indices = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
smpl1,smpl1,,,A001,SI-GA-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-GA-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_sc_indices)
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq and cellranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="10x_chromium_sc")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_NB500968_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_NB500968_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_10x_chromium_sc"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_10x_atac_protocol(self):
        """make_fastqs: 10x_atac protocol
        """
        # Sample sheet with 10xGenomics SC ATAC-seq indices
        samplesheet_10x_atac_indices = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
smpl1,smpl1,,,A001,SI-NA-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-NA-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_10x_atac_indices)
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq and cellranger-atac
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger-atac"),
                                 reads=('R1','R2','R3','I1',))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="10x_atac")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_NB500968_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_NB500968_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_10x_atac"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_10x_visium_protocol(self):
        """make_fastqs: 10x_visium protocol
        """
        # Sample sheet with 10xGenomics Visium indices
        samplesheet_10x_atac_indices = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
smpl1,smpl1,,,A001,SI-TT-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-TT-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_10x_atac_indices)
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq and cellranger-atac
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "spaceranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="10x_visium")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_NB500968_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_NB500968_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_10x_visium"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_10x_visium_protocol_forward_strand_workflow_A(self):
        """make_fastqs: 10x_visium protocol (forward strand workflow 'A')
        """
        # Sample sheet with 10xGenomics Visium indices
        samplesheet_10x_atac_indices = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
smpl1,smpl1,,,A001,SI-TT-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-TT-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_10x_atac_indices)
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq and cellranger-atac
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "spaceranger"),
                                 assert_rc_i2_override='false')
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,
                    protocol="10x_visium",
                    spaceranger_rc_i2_override=False)
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_NB500968_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_NB500968_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_10x_visium"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_10x_visium_protocol_reverse_complement_workflow_B(self):
        """make_fastqs: 10x_visium protocol (reverse complement workflow 'B')
        """
        # Sample sheet with 10xGenomics Visium indices
        samplesheet_10x_atac_indices = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
smpl1,smpl1,,,A001,SI-TT-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-TT-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_10x_atac_indices)
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq and cellranger-atac
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "spaceranger"),
                                 assert_rc_i2_override='true')
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,
                    protocol="10x_visium",
                    spaceranger_rc_i2_override=True)
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_NB500968_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_NB500968_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_10x_visium"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_10x_multiome_protocol_gex(self):
        """make_fastqs: 10x_multiome protocol (GEX data)
        """
        # Sample sheet with 10xGenomics indices
        samplesheet_10x_atac_indices = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
smpl1,smpl1,,,A001,SI-TT-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-TT-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_10x_atac_indices)
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq and cellranger-atac
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "cellranger-arc"),
                                 multiome_data="GEX")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="10x_multiome")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_NB500968_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_NB500968_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_10x_multiome"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_10x_multiome_protocol_atac(self):
        """make_fastqs: 10x_multiome protocol (ATAC data)
        """
        # Sample sheet with 10xGenomics indices
        samplesheet_10x_atac_indices = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
smpl1,smpl1,,,A001,SI-TT-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-TT-B1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_10x_atac_indices)
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I8,I8,y101",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq and cellranger-atac
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "cellranger-arc"),
                                 multiome_data="ATAC")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="10x_multiome")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_NB500968_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_NB500968_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_10x_multiome"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_10x_chromium_sc_protocol_exception_for_standard_indices(self):
        """make_fastqs: 10x_chromium_sc protocol without Chromium SC indices raises exception
        """
        # Sample sheet with standard (non-10xGenomics Chromium SC) indices
        samplesheet_standard_indices = """[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
smpl1,smpl1,,,A001,GATCTACT,10xGenomics,
smpl2,smpl2,,,A005,CTTAGGAC,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_standard_indices)
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq and cellranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        self.assertRaises(Exception,
                          make_fastqs,
                          ap,
                          protocol="10x_chromium_sc")

    #unittest.skip("Skipped")
    def test_make_fastqs_parse_evercode_protocol(self):
        """
        make_fastqs: parse_evercode protocol
        """
        # Sample sheet
        samplesheet_parse_evercode = """[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTA,,
Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_parse_evercode)
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y76,I8,I8,y86",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_bases_mask="y76,I8,I8,y86",
                                 assert_adapter='',
                                 assert_adapter2='')
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_NB500968_00002_AHGXXXX"),
                 sample_sheet=sample_sheet)
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,protocol="parse_evercode")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_NB500968_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_NB500968_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_NB500968_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs_parse_evercode"),
                       "bcl2fastq",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_missing_fastqs_no_placeholders(self):
        """make_fastqs: missing fastqs, no placeholders
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 missing_fastqs=(
                                     "Sample1_S1_L001_R1_001.fastq.gz",
                                     "Sample1_S1_L001_R2_001.fastq.gz",))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_M00879_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        self.assertRaises(Exception,
                          make_fastqs,
                          ap,
                          protocol="standard",
                          create_empty_fastqs=False)
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_M00879_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,None)
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_M00879_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs"),):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for subdir in ("bcl2fastq",):
            self.assertFalse(os.path.exists(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertFalse(os.path.exists(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        self.assertTrue(
            os.path.exists(
                os.path.join(analysis_dir,
                             "logs",
                             "002_make_fastqs",
                             "missing_fastqs.log")))

    #@unittest.skip("Skipped")
    def test_make_fastqs_missing_fastqs_with_placeholders(self):
        """make_fastqs: missing fastqs with placeholders
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 missing_fastqs=(
                                     "Sample1_S1_L001_R1_001.fastq.gz",
                                     "Sample1_S1_L001_R2_001.fastq.gz",))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_M00879_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        make_fastqs(ap,
                    protocol="standard",
                    create_empty_fastqs=True)
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_M00879_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_M00879_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info"):
            self.assertTrue(os.path.exists(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
        self.assertTrue(
            os.path.exists(
                os.path.join(analysis_dir,
                             "logs",
                             "002_make_fastqs",
                             "missing_fastqs.log")))

    #@unittest.skip("Skipped")
    def test_make_fastqs_handle_bcl2fastq2_failure(self):
        """make_fastqs: handle bcl2fastq2 failure
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq which will fail (i.e.
        # return non-zero exit code)
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 exit_code=1)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_M00879_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertEqual(ap.params.unaligned_dir,None)
        self.assertEqual(ap.params.barcode_analysis_dir,None)
        self.assertEqual(ap.params.stats_file,None)
        self.assertRaises(Exception,
                          make_fastqs,
                          ap,
                          protocol="standard")
        # Check outputs
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(
                             self.wd,
                             "171020_M00879_00002_AHGXXXX_analysis",
                             "primary_data"))
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,None)
        analysis_dir = os.path.join(
            self.wd,
            "171020_M00879_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs"),):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for subdir in ("bcl2fastq",):
            self.assertFalse(os.path.exists(
                os.path.join(analysis_dir,subdir)),
                            "Found subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertFalse(os.path.exists(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_explicitly_specify_platform(self):
        """make_fastqs: explicitly specify the platform
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_UNKNOWN_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 platform="miseq")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_UNKNOWN_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,None)
        self.assertEqual(ap.params.barcode_analysis_dir,None)
        self.assertEqual(ap.params.stats_file,None)
        make_fastqs(ap,
                    protocol="standard",
                    platform="miseq")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_UNKNOWN_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_UNKNOWN_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_UNKNOWN_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_specify_platform_via_metadata(self):
        """make_fastqs: implicitly specify the platform via metadata
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_UNKNOWN_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 platform="miseq")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_UNKNOWN_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertTrue(ap.metadata.platform is None)
        ap.metadata["platform"] = "miseq"
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,None)
        self.assertEqual(ap.params.barcode_analysis_dir,None)
        self.assertEqual(ap.params.stats_file,None)
        make_fastqs(ap,protocol="standard")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_UNKNOWN_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_UNKNOWN_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_UNKNOWN_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_make_fastqs_unknown_protocol(self):
        """make_fastqs: fails with unknown protocol
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq and cellranger executables
        MockBcl2fastq2Exe.create(os.path.join(self.bin,"bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_M00879_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        self.assertRaises(Exception,
                          make_fastqs,
                          ap,
                          protocol="undefined_protocol")

    #@unittest.skip("Skipped")
    def test_make_fastqs_specify_lane_subsets(self):
        """make_fastqs: specify lane subsets
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
        # Do the test
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_SN7001250_00002_AHGXXXX"))
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,None)
        self.assertEqual(ap.params.barcode_analysis_dir,None)
        self.assertEqual(ap.params.stats_file,None)
        make_fastqs(ap,
                    protocol="standard",
                    lane_subsets=(
                        subset(lanes=(1,2,)),
                        subset(lanes=(3,4,5,6,7,8,))))
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_SN7001250_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_SN7001250_00002_AHGXXXX_analysis")
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
    def test_make_fastqs_existing_projects_info_file(self):
        """make_fastqs: handle existing projects.info file
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_M00879_00002_AHGXXXX",
            "miseq",
            top_dir=self.wd)
        illumina_run.create()
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,"bcl2fastq"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        ap = AutoProcess(settings=self.settings)
        ap.setup(os.path.join(self.wd,
                              "171020_M00879_00002_AHGXXXX"))
        # Make a projects.info file
        with open(os.path.join(ap.analysis_dir,'projects.info'),"wt") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
Proj1\tPR1,PR2\t.\t.\t.\t.\t.\t.
""")
        ap.params['project_metadata'] = 'projects.info'
        # Do the test
        self.assertTrue(ap.params.sample_sheet is not None)
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertTrue(ap.params.primary_data_dir is None)
        self.assertFalse(ap.params.acquired_primary_data)
        self.assertTrue(ap.params.project_metadata is not None)
        make_fastqs(ap,protocol="standard")
        # Check parameters
        self.assertEqual(ap.params.bases_mask,"auto")
        self.assertEqual(ap.params.primary_data_dir,
                         os.path.join(self.wd,
                                      "171020_M00879_00002_AHGXXXX_analysis",
                                      "primary_data"))
        self.assertTrue(ap.params.acquired_primary_data)
        self.assertEqual(ap.params.unaligned_dir,"bcl2fastq")
        self.assertEqual(ap.params.barcode_analysis_dir,"barcode_analysis")
        self.assertEqual(ap.params.stats_file,"statistics.info")
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "171020_M00879_00002_AHGXXXX_analysis")
        for subdir in (os.path.join("primary_data",
                                    "171020_M00879_00002_AHGXXXX"),
                       os.path.join("logs",
                                    "002_make_fastqs"),
                       "bcl2fastq",
                       "barcode_analysis",):
            self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,subdir)),
                            "Missing subdir: %s" % subdir)
        for filen in ("statistics.info",
                      "statistics_full.info",
                      "per_lane_statistics.info",
                      "per_lane_sample_stats.info",
                      "projects.info",
                      "processing_qc.html"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
