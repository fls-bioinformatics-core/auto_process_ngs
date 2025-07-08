#######################################################################
# Unit tests for bcl2fastq/pipeline.py (10x_multiome* protocols)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestMakeFastqs(BaseMakeFastqsTestCase):
    """
    Tests for MakeFastqs pipeline (10x_multiome* protocols)
    """
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
        # Create mock bcl2fastq and cellranger-arc
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
        # Create mock bcl2fastq and cellranger-arc
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
        # Create mock bcl2fastq and cellranger-arc
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
        # Create mock bcl2fastq and cellranger-arc
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
    def test_makefastqs_10x_multiome_atac(self):
        """
        MakeFastqs: '10x_multiome_atac' protocol
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y50,I10,I24,y90",
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
        # Create mock bcl2fastq and cellranger-arc
        MockBcl2fastq2Exe.create(os.path.join(self.bin,"bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "cellranger-arc"),
                                 version="2.0.0",
                                 multiome_data="ATAC",
                                 assert_bases_mask="Y50,I8n2,Y24,Y90",
                                 assert_filter_single_index=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_multiome_atac")
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
    def test_makefastqs_10x_multiome_gex(self):
        """
        MakeFastqs: '10x_multiome_gex' protocol
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            bases_mask="y50,I10,I24,y90",
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
        # Create mock bcl2fastq and cellranger-arc
        MockBcl2fastq2Exe.create(os.path.join(self.bin,"bcl2fastq"))
        Mock10xPackageExe.create(os.path.join(self.bin,
                                              "cellranger-arc"),
                                 version="2.0.0",
                                 multiome_data="GEX",
                                 assert_bases_mask="Y28n22,I10,I10n14,Y90",
                                 assert_filter_dual_index=True)
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        analysis_dir = os.path.join(self.wd,"analysis")
        os.mkdir(analysis_dir)
        # Do the test
        p = MakeFastqs(run_dir,sample_sheet,protocol="10x_multiome_gex")
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
