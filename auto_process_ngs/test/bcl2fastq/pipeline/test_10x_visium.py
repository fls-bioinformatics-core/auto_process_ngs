#######################################################################
# Unit tests for bcl2fastq/pipeline.py (10x_visium protocol)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestMakeFastqs(BaseMakeFastqsTestCase):
    """
    Tests for MakeFastqs pipeline (10x_multiome protocol)
    """
    #@unittest.skip("Skipped")
    def test_makefastqs_10x_visium_protocol_110(self):
        """
        MakeFastqs: '10x_visium' protocol (Spaceranger 1.1.0)
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
                                              "spaceranger"),
                                 version='1.1.0',
                                 assert_bases_mask="Y28N73,I10,I10,Y50N51")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_visium")
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
                      "spaceranger_qc_summary.html",):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,filen)),
                            "Missing file: %s" % filen)

    #@unittest.skip("Skipped")
    def test_makefastqs_10x_visium_protocol_131(self):
        """
        MakeFastqs: '10x_visium' protocol (Spaceranger 1.3.1)
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
                                              "spaceranger"),
                                 version='1.3.1',
                                 assert_bases_mask="Y28N73,I10,I10,Y50N51")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_visium")
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
        self.assertEqual(p.output.spaceranger_info,
                         (os.path.join(self.bin,"spaceranger"),
                          "spaceranger",
                          "1.3.1"))
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
    def test_makefastqs_10x_visium_protocol_211(self):
        """
        MakeFastqs: '10x_visium' protocol (Spaceranger 2.1.1)
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
                                              "spaceranger"),
                                 version='2.1.1',
                                 assert_bases_mask="Y28N73,I10,I10,Y50N51")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_visium")
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
        self.assertEqual(p.output.spaceranger_info,
                         (os.path.join(self.bin,"spaceranger"),
                          "spaceranger",
                          "2.1.1"))
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
    def test_makefastqs_10x_visium_protocol_300(self):
        """
        MakeFastqs: '10x_visium' protocol (Spaceranger 3.0.0)
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
                                              "spaceranger"),
                                 version='3.0.0',
                                 assert_bases_mask="Y28N73,I10,I10,Y50N51")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_visium")
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
        self.assertEqual(p.output.spaceranger_info,
                         (os.path.join(self.bin,"spaceranger"),
                          "spaceranger",
                          "3.0.0"))
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
    def test_makefastqs_10x_visium_bcl2fastq_illumina_indexes(self):
        """
        MakeFastqs: '10x_visium' protocol (bcl2fastq/Illumina indexes)
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
smpl1,smpl1,,,SI-TT-A1,TGTCCCAACG,SI-TT-A1,TGGACATCGA,10xGenomics,
smpl2,smpl2,,,SI-TT-B1,GTCCCATCAA,SI-TT-B1,GTCACGTTCG,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_visium_indices)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin, "bcl2fastq"),
                                 assert_bases_mask="Y28N73,I10,I10,Y50N51")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_visium")
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
    def test_makefastqs_10x_visium_bclconvert_illumina_indexes(self):
        """
        MakeFastqs: '10x_visium' protocol (bcl-convert/Illumina indexes)
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
smpl1,smpl1,,,SI-TT-A1,TGTCCCAACG,SI-TT-A1,TGGACATCGA,10xGenomics,
smpl2,smpl2,,,SI-TT-B1,GTCCCATCAA,SI-TT-B1,GTCACGTTCG,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_visium_indices)
        # Create mock bcl-convert
        MockBclConvertExe.create(os.path.join(self.bin, "bcl-convert"),
                                 assert_override_cycles="Y28N73;I10;I10;Y50N51")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir, sample_sheet,
                       bcl_converter="bcl-convert",
                       protocol="10x_visium")
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
    def test_makefastqs_10x_visium_protocol_truncate_reads(self):
        """
        MakeFastqs: '10x_visium' protocol (truncate reads)
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
                                              "spaceranger"),
                                 version='3.0.0',
                                 assert_bases_mask="Y40N61,I10,I10,Y90N11")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_visium",
                       r1_length=40,
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
        self.assertEqual(p.output.spaceranger_info,
                         (os.path.join(self.bin,"spaceranger"),
                          "spaceranger",
                          "3.0.0"))
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
    def test_makefastqs_10x_visium_protocol_forward_strand_workflow_A(self):
        """
        MakeFastqs: '10x_visium' protocol (forward strand workflow 'A')
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
                                              "spaceranger"),
                                 assert_rc_i2_override='false')
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,
                       sample_sheet,
                       protocol="10x_visium",
                       spaceranger_rc_i2_override=False)
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
        self.assertEqual(p.output.spaceranger_info,
                         (os.path.join(self.bin,"spaceranger"),
                          "spaceranger",
                          "3.0.0"))
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
    def test_makefastqs_10x_visium_protocol_reverse_complement_workflow_B(self):
        """
        MakeFastqs: '10x_visium' protocol (reverse complement workflow 'B')
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
                                              "spaceranger"),
                                 assert_rc_i2_override='true')
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,
                       sample_sheet,
                       protocol="10x_visium",
                       spaceranger_rc_i2_override=True)
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
        self.assertEqual(p.output.spaceranger_info,
                         (os.path.join(self.bin,"spaceranger"),
                          "spaceranger",
                          "3.0.0"))
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
    def test_makefastqs_10x_visium_v1_bcl2fastq_illumina_indexes(self):
        """
        MakeFastqs: '10x_visium_v1' protocol (bcl2fastq/Illumina indexes)
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
smpl1,smpl1,,,SI-TT-A1,TGTCCCAACG,SI-TT-A1,TGGACATCGA,10xGenomics,
smpl2,smpl2,,,SI-TT-B1,GTCCCATCAA,SI-TT-B1,GTCACGTTCG,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_visium_indices)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin, "bcl2fastq"),
                                 assert_bases_mask="Y28N73,I10,I10,Y90N11")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_visium_v1")
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
    def test_makefastqs_10x_visium_hd_bcl2fastq_illumina_indexes(self):
        """
        MakeFastqs: '10x_visium_hd' protocol (bcl2fastq/Illumina indexes)
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
smpl1,smpl1,,,SI-TT-A1,TGTCCCAACG,SI-TT-A1,TGGACATCGA,10xGenomics,
smpl2,smpl2,,,SI-TT-B1,GTCCCATCAA,SI-TT-B1,GTCACGTTCG,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_visium_indices)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin, "bcl2fastq"),
                                 assert_bases_mask="Y43N58,I10,I10,Y50N51")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_visium_hd")
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
    def test_makefastqs_10x_visium_hd_3prime_bcl2fastq_illumina_indexes(self):
        """
        MakeFastqs: '10x_visium_hd_3prime' protocol (bcl2fastq/Illumina indexes)
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
smpl1,smpl1,,,SI-TT-A1,TGTCCCAACG,SI-TT-A1,TGGACATCGA,10xGenomics,
smpl2,smpl2,,,SI-TT-B1,GTCCCATCAA,SI-TT-B1,GTCACGTTCG,10xGenomics,
"""
        sample_sheet = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write(samplesheet_visium_indices)
        # Create mock bcl2fastq
        MockBcl2fastq2Exe.create(os.path.join(self.bin, "bcl2fastq"),
                                 assert_bases_mask="Y43N58,I10,I10,Y75N26")
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_visium_hd_3prime")
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