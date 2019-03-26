#######################################################################
# Tests for tenx_genomics_utils.py module
#######################################################################

import unittest
import tempfile
import os
import shutil
from bcftbx.mock import MockIlluminaRun
from auto_process_ngs.mock import MockBcl2fastq2Exe
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.tenx_genomics_utils import *

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestFlowCellId(unittest.TestCase):
    """
    Tests for the 'flow_cell_id' function
    """
    def test_flow_cell_id(self):
        """
        flow_cell_id: returns correct ID
        """
        self.assertEqual(flow_cell_id("170426_K00311_0033_AHJCY7BBXX"),
                         "HJCY7BBXX")

class TestHasChromiumSCIndices(unittest.TestCase):
    """
    Tests for the 'has_chromium_sc_indices' function
    """
    def setUp(self):
        # Sample sheet contents
        self.sample_sheet_with_chromium_indices = """[Header]
IEMFileVersion,4

[Reads]
76
76

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,smpl1,smpl1,,,A001,SI-GA-A1,10xGenomics,
2,smpl2,smpl2,,,A005,SI-GA-B1,10xGenomics,
3,smpl3,smpl3,,,A006,SI-GA-C1,10xGenomics,
4,smpl4,smpl4,,,A007,SI-GA-D1,10xGenomics,
"""
        self.sample_sheet_with_atac_indices = """[Header]
IEMFileVersion,4

[Reads]
76
76

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,smpl1,smpl1,,,A001,SI-NA-A1,10xGenomics,
2,smpl2,smpl2,,,A005,SI-NA-B1,10xGenomics,
3,smpl3,smpl3,,,A006,SI-NA-C1,10xGenomics,
4,smpl4,smpl4,,,A007,SI-NA-D1,10xGenomics,
"""
        self.sample_sheet_standard_indices = """[Header]
IEMFileVersion,4

[Reads]
76
76

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,smpl1,smpl1,,,A001,GGTTTACT,standard,
2,smpl2,smpl2,,,A005,GTAATCTT,standard,
3,smpl3,smpl3,,,A006,CCACTTAT,standard,
4,smpl4,smpl4,,,A007,CACTCGGA,standard,
"""
        self.sample_sheet_mixed_indices = """[Header]
IEMFileVersion,4

[Reads]
76
76

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,smpl1,smpl1,,,A001,GGTTTACT,standard,
2,smpl2,smpl2,,,A005,GTAATCTT,standard,
3,smpl3,smpl3,,,A006,SI-GA-C1,10xGenomics,
4,smpl4,smpl4,,,A007,SI-GA-D1,10xGenomics,
"""
        # Make temporary working dir
        self.wd = tempfile.mkdtemp(suffix="TestHasChromiumSCIndices")
    def tearDown(self):
        # Remove temp dir
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)
    def _make_sample_sheet(self,content):
        # Make a sample sheet with supplied content in
        # the temporary area
        sample_sheet_csv = os.path.join(self.wd,"SampleSheet.csv")
        with open(sample_sheet_csv,'w') as s:
            s.write(content)
        return sample_sheet_csv
    def test_sample_sheet_all_chromium_sc_3_v2_indices(self):
        """
        has_chromium_indices: sample sheet with Chromium SC 3'v2 indices only
        """
        s = self._make_sample_sheet(
            self.sample_sheet_with_chromium_indices)
        self.assertTrue(has_chromium_sc_indices(s))
    def test_sample_sheet_all_sc_atac_indices(self):
        """
        has_chromium_indices: sample sheet with 10x scATAC-seq indices only
        """
        s = self._make_sample_sheet(
            self.sample_sheet_with_atac_indices)
        self.assertTrue(has_chromium_sc_indices(s))
    def test_sample_sheet_no_chromium_indices(self):
        """
        has_chromium_indices: sample sheet with no Chromium SC indices
        """
        s = self._make_sample_sheet(
            self.sample_sheet_standard_indices)
        self.assertFalse(has_chromium_sc_indices(s))
    def test_sample_sheet_some_chromium_sc_3_v2_indices(self):
        """
        has_chromium_indices: sample sheet with some Chromium SC 3'v2 indices
        """
        s = self._make_sample_sheet(
            self.sample_sheet_mixed_indices)
        self.assertTrue(has_chromium_sc_indices(s))

class TestCellrangerInfo(unittest.TestCase):
    """
    Tests for the cellranger_info function
    """
    def setUp(self):
        # Make temporary working dir
        self.wd = tempfile.mkdtemp(suffix="TestCellrangerInfo")
        # Store the original state of PATH env var
        self.original_path = os.environ['PATH']
    def tearDown(self):
        # Restore the PATH env var
        os.environ['PATH'] = self.original_path
        # Remove temp dir
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)
    def _make_mock_cellranger_201(self):
        # Make a fake cellranger 2.0.1 executable
        cellranger_201 = os.path.join(self.wd,"cellranger")
        with open(cellranger_201,'w') as fp:
            fp.write("#!/bin/bash\ncat <<EOF\ncellranger $1  (2.0.1)\nCopyright (c) 2017 10x Genomics, Inc.  All rights reserved.\n-------------------------------------------------------------------------------\n\nUsage:\n    cellranger mkfastq\n\n    cellranger count\n    cellranger aggr\n    cellranger reanalyze\n    cellranger mkloupe\n    cellranger mat2csv\n\n    cellranger mkgtf\n    cellranger mkref\n\n    cellranger vdj\n\n    cellranger mkvdjref\n\n    cellranger testrun\n    cellranger upload\n    cellranger sitecheckEOF")
        os.chmod(cellranger_201,0775)
        return cellranger_201
    def _make_mock_cellranger_atac_101(self):
        # Make a fake cellranger-atac 1.0.1 executable
        cellranger_atac_101 = os.path.join(self.wd,"cellranger-atac")
        with open(cellranger_atac_101,'w') as fp:
            fp.write("#!/bin/bash\ncat <<EOF\ncellranger-atac  (1.0.1)\nCopyright (c) 2018 10x Genomics, Inc.  All rights reserved.\n-------------------------------------------------------------------------------\n\nUsage:\n    cellranger-atac mkfastq\n\n    cellranger-atac count\n\n    cellranger-atac testrun\n    cellranger-atac upload\n    cellranger-atac sitecheckEOF")
        os.chmod(cellranger_atac_101,0775)
        return cellranger_atac_101

    def test_cellranger_201(self):
        """cellranger_info: collect info for cellranger 2.0.1
        """
        cellranger = self._make_mock_cellranger_201()
        self.assertEqual(cellranger_info(path=cellranger),
                         (cellranger,'cellranger','2.0.1'))

    def test_cellranger_201_on_path(self):
        """cellranger_info: collect info for cellranger 2.0.1 from PATH
        """
        os.environ['PATH'] = "%s%s%s" % (os.environ['PATH'],
                                         os.pathsep,
                                         self.wd)
        cellranger = self._make_mock_cellranger_201()
        self.assertEqual(cellranger_info(name='cellranger'),
                         (cellranger,'cellranger','2.0.1'))

    def test_cellranger_atac_101(self):
        """cellranger_info: collect info for cellranger-atac 1.0.1
        """
        cellranger_atac = self._make_mock_cellranger_atac_101()
        self.assertEqual(cellranger_info(path=cellranger_atac),
                         (cellranger_atac,'cellranger-atac','1.0.1'))

    def test_cellranger_atac_101_on_path(self):
        """cellranger_info: collect info for cellranger-atac 1.0.1 from PATH
        """
        os.environ['PATH'] = "%s%s%s" % (os.environ['PATH'],
                                         os.pathsep,
                                         self.wd)
        cellranger_atac = self._make_mock_cellranger_atac_101()
        self.assertEqual(cellranger_info(name='cellranger-atac'),
                         (cellranger_atac,'cellranger-atac','1.0.1'))

class TestMetricsSummary(unittest.TestCase):
    """
    Tests for the 'MetricsSummary' class
    """
    def test_metrics_summary(self):
        """MetricsSummary: check estimated number of cells is extracted
        """
        metrics_summary = """Estimated Number of Cells,Mean Reads per Cell,Median Genes per Cell,Number of Reads,Valid Barcodes,Reads Mapped Confidently to Transcriptome,Reads Mapped Confidently to Exonic Regions,Reads Mapped Confidently to Intronic Regions,Reads Mapped Confidently to Intergenic Regions,Reads Mapped Antisense to Gene,Sequencing Saturation,Q30 Bases in Barcode,Q30 Bases in RNA Read,Q30 Bases in Sample Index,Q30 Bases in UMI,Fraction Reads in Cells,Total Genes Detected,Median UMI Counts per Cell
"2,272","107,875","1,282","245,093,084",98.3%,69.6%,71.9%,6.1%,3.2%,4.4%,51.3%,98.5%,79.2%,93.6%,98.5%,12.0%,"16,437","2,934"
"""
        m = MetricsSummary(metrics_summary)
        self.assertEqual(m.estimated_number_of_cells,2272)

class TestRunCellrangerMkfastq(unittest.TestCase):
    """
    Tests for the 'run_cellranger_mkfastq' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestRunCellrangerMkfastq')
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

    def test_run_cellranger_mkfastq(self):
        """run_cellranger_mkfastq: check cellranger is executed
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_SN7001250_00002_AHGXXXX",
            "hiseq",
            top_dir=self.wd)
        illumina_run.create()
        # Mock sample sheet with chromium indices
        sample_sheet_file = os.path.join(self.wd,"samplesheet.csv")
        with open(sample_sheet_file,'w') as fp:
            fp.write("""[Header]
IEMFileVersion,4

[Reads]
76
76

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,smpl1,smpl1,,,A001,SI-GA-A1,10xGenomics,
2,smpl2,smpl2,,,A005,SI-GA-B1,10xGenomics,
3,smpl3,smpl3,,,A006,SI-GA-C1,10xGenomics,
4,smpl4,smpl4,,,A007,SI-GA-D1,10xGenomics,
""")
        # Create mock bcl2fastq and cellranger executables
        MockBcl2fastq2Exe.create(os.path.join(self.bin,"bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Output dir
        output_dir = "bcl2fastq"
        self.assertFalse(os.path.exists("HGXXXX"))
        self.assertFalse(os.path.exists("cellranger_qc_summary.html"))
        self.assertFalse(os.path.exists(output_dir))
        # Run 'cellranger mkfastq'
        exit_code = run_cellranger_mkfastq(sample_sheet_file,
                                           illumina_run.dirn,
                                           output_dir)
        # Check outputs
        self.assertEqual(exit_code,0)
        self.assertTrue(os.path.isdir("HGXXXX"))
        self.assertTrue(os.path.isfile("cellranger_qc_summary.html"))
        self.assertTrue(os.path.isdir(output_dir))

    def test_run_cellranger_mkfastq_subset_of_lanes(self):
        """run_cellranger_mkfastq: check cellranger is executed for subset of lanes
        """
        # Create mock source data
        illumina_run = MockIlluminaRun(
            "171020_SN7001250_00002_AHGXXXX",
            "hiseq",
            top_dir=self.wd)
        illumina_run.create()
        # Mock sample sheet with chromium indices
        sample_sheet_file = os.path.join(self.wd,"samplesheet.csv")
        with open(sample_sheet_file,'w') as fp:
            fp.write("""[Header]
IEMFileVersion,4

[Reads]
76
76

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,smpl1,smpl1,,,A001,SI-GA-A1,10xGenomics,
2,smpl2,smpl2,,,A005,SI-GA-B1,10xGenomics,
3,smpl3,smpl3,,,A006,SI-GA-C1,10xGenomics,
4,smpl4,smpl4,,,A007,SI-GA-D1,10xGenomics,
""")
        # Create mock bcl2fastq and cellranger executables
        MockBcl2fastq2Exe.create(os.path.join(self.bin,"bcl2fastq"))
        MockCellrangerExe.create(os.path.join(self.bin,"cellranger"))
        os.environ['PATH'] = "%s:%s" % (self.bin,
                                        os.environ['PATH'])
        # Output dir
        output_dir = "bcl2fastq"
        self.assertFalse(os.path.exists("HGXXXX_34"))
        self.assertFalse(os.path.exists("cellranger_qc_summary_34.html"))
        self.assertFalse(os.path.exists(output_dir))
        # Run 'cellranger mkfastq'
        exit_code = run_cellranger_mkfastq(sample_sheet_file,
                                           illumina_run.dirn,
                                           output_dir,
                                           lanes="3,4")
        # Check outputs
        self.assertEqual(exit_code,0)
        self.assertTrue(os.path.isdir("HGXXXX_34"))
        self.assertTrue(os.path.isfile("cellranger_qc_summary_34.html"))
        self.assertTrue(os.path.isdir(output_dir))
