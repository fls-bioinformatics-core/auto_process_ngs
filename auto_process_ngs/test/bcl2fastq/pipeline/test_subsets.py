#######################################################################
# Unit tests for bcl2fastq/pipeline.py (subsetting)
#######################################################################

# All imports declared in __init__.py file
from . import *

class TestMakeFastqs(BaseMakeFastqsTestCase):
    """
    Tests for MakeFastqs pipeline (subsetting)
    """
    #@unittest.skip("Skipped")
    def test_makefastqs_multiple_subsets_multiple_lanes(self):
        """
        MakeFastqs: explicitly define subsets (multiple lanes)
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
                       poll_interval=POLL_INTERVAL)
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
    def test_makefastqs_multiple_subsets_multiple_protocols(self):
        """
        MakeFastqs: explicitly define subsets (multiple protocols)
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
                       poll_interval=POLL_INTERVAL)
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
    def test_makefastqs_multiple_subsets_multiple_protocols_rerun_completed_pipeline(self):
        """
        MakeFastqs: explicitly define subsets (multiple protocols, rerun completed pipeline)
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
5,GH5,GH5,,,N701,SI-GA-A2,S503,,GH,
6,IJ6,IJ6,,,N701,SI-GA-A1,S502,,IJ,
7,KL7,KL7,,,N701,SI-GA-A1,S502,,KL,
8,MN8,MN8,,,N701,SI-GA-A1,S503,,MN,
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
        for lanes in ((1,2,3),(4,),(5,6,),(7,8,)):
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
                           subset(lanes=(5,6,),protocol="10x_chromium_sc"),
                           subset(lanes=(7,8,),protocol="10x_atac")),
                       analyse_barcodes=False)
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "save.bcl2fastq.L123",
                       "save.bcl2fastq.L4",
                       "save.bcl2fastq.L56",
                       "save.bcl2fastq.L78",):
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
    def test_makefastqs_multiple_subsets_multiple_protocols_rerun_incomplete_pipeline(self):
        """
        MakeFastqs: explicitly define subsets (multiple protocols, rerun incomplete pipeline)
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
5,GH5,GH5,,,N701,SI-GA-A2,S503,,GH,
6,IJ6,IJ6,,,N701,SI-GA-A1,S502,,IJ,
7,KL7,KL7,,,N701,SI-GA-A1,S502,,KL,
8,MN8,MN8,,,N701,SI-GA-A1,S503,,MN,
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
        for lanes in ((1,2,3,),(4,),(5,6,),(7,8,)):
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
                           subset(lanes=(5,6,),protocol="10x_chromium_sc"),
                           subset(lanes=(7,8,),protocol="10x_atac")),
                       analyse_barcodes=False)
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
        self.assertEqual(p.output.missing_fastqs,[])
        for subdir in (os.path.join("primary_data",
                                    "171020_SN7001250_00002_AHGXXXX"),
                       "bcl2fastq",
                       "save.bcl2fastq.L123",
                       "save.bcl2fastq.L4",
                       "save.bcl2fastq.L56",
                       "save.bcl2fastq.L78",):
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
    def test_makefastqs_exception_for_inconsistent_subsets(self):
        """
        MakeFastqs: raise exception for inconsistent subset definitions
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
    def test_makefastqs_multiple_subsets_lanes_automatically_split_by_indices(self):
        """
        MakeFastqs: automatically define subsets (lanes split by indices)
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
                       poll_interval=POLL_INTERVAL)
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
    def test_makefastqs_multiple_subsets_set_out_dirs(self):
        """
        MakeFastqs: set output directories for multiple subsets
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
                       poll_interval=POLL_INTERVAL)
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
                   protocol="10x_atac",
                   trim_adapters=False)
        self.assertEqual(s['lanes'],[1,2])
        self.assertEqual(s['protocol'],"10x_atac")
        self.assertEqual(s['trim_adapters'], False)

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
