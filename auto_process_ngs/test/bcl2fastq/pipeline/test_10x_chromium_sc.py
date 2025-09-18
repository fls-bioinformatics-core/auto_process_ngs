#######################################################################
# Unit tests for bcl2fastq/pipeline.py (10x_chromium_sc protocol)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestMakeFastqs(BaseMakeFastqsTestCase):
    """
    Tests for MakeFastqs pipeline (10x_chromium_sc protocol)
    """
    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_501(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol (Cellranger 5.0.1)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I10,I10,y101",
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
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
smpl1,smpl1,,,A001,SI-GA-A1,A001,SI-GA-A1,10xGenomics,
smpl2,smpl2,,,A005,SI-GA-B1,A001,SI-GA-A1,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_sc_indices)
        # Create mock bcl2fastq and cellranger
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,
                                              "cellranger"),
                                 version="5.0.1",
                                 assert_bases_mask="Y101,I10,I10,Y101")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
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
            bases_mask="y101,I10,I10,y101",
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
                                 version="6.0.0",
                                 assert_bases_mask="Y101,I10,I10,Y101")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_700(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol (Cellranger 7.0.0)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I10,I10,y101",
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
                                 version="7.0.0",
                                 assert_bases_mask="Y101,I10,I10,Y101")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
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
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "7.0.0"))
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_800(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol (Cellranger 8.0.0)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I10,I10,y101",
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
                                 version="8.0.0",
                                 assert_bases_mask="Y101,I10,I10,Y101")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
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
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "8.0.0"))
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_900(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol (Cellranger 9.0.0)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I10,I10,y101",
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
                                 version="9.0.0",
                                 assert_bases_mask="Y101,I10,I10,Y101")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
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
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "9.0.0"))
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_bcl2fastq_illumina_indexes(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol (bcl2fastq/Illumina indexes)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I10,I10,y101",
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
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
TX1,TX1,,,SI-TT-E9,TGTCCCAACG,SI-TT-E9,TGGACATCGA,10xGenomics,
TX2,TX2,,,SI-TT-F9,GTCCCATCAA,SI-TT-F9,GTCACGTTCG,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_sc_indices)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin,
                                              "bcl2fastq"),
                                 assert_bases_mask="Y101,I10,I10,Y101")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_bclconvert_illumina_indexes(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol (bcl-convert/Illumina indexes)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I10,I10,y101",
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
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
TX1,TX1,,,SI-TT-E9,TGTCCCAACG,SI-TT-E9,TGGACATCGA,10xGenomics,
TX2,TX2,,,SI-TT-F9,GTCCCATCAA,SI-TT-F9,GTCACGTTCG,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_chromium_sc_indices)
        # Create mock bcl-convert
        MockBclConvertExe.create(os.path.join(self.bin, "bcl-convert"),
                                 assert_override_cycles="Y101;I10;I10;Y101")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir, sample_sheet,
                       protocol="10x_chromium_sc",
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
                         (os.path.join(self.bin, "bcl-convert"),
                          "BCL Convert",
                          "3.7.5"))
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
            bases_mask="y101,I10,I10,y101",
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
                                              "cellranger"),
                                 assert_bases_mask="Y101,I10,I10,Y101")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
        status = p.run(analysis_dir,
                       poll_interval=POLL_INTERVAL)
        self.assertEqual(status,0)
        # Check outputs
        self.assertEqual(p.output.platform,"illumina")
        self.assertEqual(p.output.flow_cell_mode,None)
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
                          "7.0.0"))
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
                      "seq_len_statistics.info",
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_chromium_sc_protocol_truncate_reads(self):
        """
        MakeFastqs: '10x_chromium_sc' protocol (truncate reads)
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y101,I10,I10,y101",
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
                                 version="8.0.0",
                                 assert_bases_mask="Y28N73,I10,I10,Y90N11")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc",
                       r1_length=28,
                       r2_length=90)
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
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "8.0.0"))
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
            bases_mask="y101,I10,I10,y101",
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
                  "seq_len_statistics.info",
                  "processing_qc.html",
                  "cellranger_qc_summary.html",):
            with open(os.path.join(analysis_dir,f),'wt') as fp:
                fp.write("")
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_chromium_sc")
        status = p.run(analysis_dir,
                       #default_runner=SimpleJobRunner(join_logs=True),
                       #verbose=True,
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
        self.assertEqual(p.output.cellranger_info,
                         (os.path.join(self.bin,"cellranger"),
                          "cellranger",
                          "7.0.0"))
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
                      "processing_qc.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)
