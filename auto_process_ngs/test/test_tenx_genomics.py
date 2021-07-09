#######################################################################
# Tests for tenx_genomics_utils.py module
#######################################################################

import unittest
import tempfile
import os
import shutil
from bcftbx.mock import MockIlluminaRun
from bcftbx.mock import RunInfoXml
from bcftbx.utils import mkdirs
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.mock import MockBcl2fastq2Exe
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import MockAnalysisDir
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock10xdata import METRICS_SUMMARY
from auto_process_ngs.mock10xdata import ATAC_SUMMARY
from auto_process_ngs.mock10xdata import ATAC_SUMMARY_2_0_0
from auto_process_ngs.mock10xdata import MULTIOME_SUMMARY
from auto_process_ngs.mock10xdata import MULTIOME_SUMMARY_2_0_0
from auto_process_ngs.mock10xdata import MULTIOME_LIBRARIES
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

class TestHas10xindices(unittest.TestCase):
    """
    Tests for the 'has_10x_indices' function
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
        self.sample_sheet_with_chromium_indices_A10 = """[Header]
IEMFileVersion,4

[Reads]
76
76

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
1,smpl1,smpl1,,,A001,SI-GA-A10,10xGenomics,
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
        self.sample_sheet_with_visium_indices = """[Header]
IEMFileVersion,4

[Reads]
76
76

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index,index2,Sample_Project,Description
1,smpl1,smpl1,,,A001,SI-TT-A1,,SI-TT-A1,,10xGenomics,
2,smpl2,smpl2,,,A005,SI-TT-B1,,SI-TT-B1,,10xGenomics,
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
        self.sample_sheet_no_indices = """[Header]
IEMFileVersion,4

[Reads]
76
76

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Sample_Project,Description
1,smpl1,smpl1,,GGTTTACT,standard,
2,smpl2,smpl2,,GTAATCTT,standard,
3,smpl3,smpl3,,SI-GA-C1,10xGenomics,
4,smpl4,smpl4,,SI-GA-D1,10xGenomics,
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
        has_10x_indices: sample sheet with Chromium SC 3'v2 indices only
        """
        s = self._make_sample_sheet(
            self.sample_sheet_with_chromium_indices)
        self.assertTrue(has_10x_indices(s))
    def test_sample_sheet_with_star10_chromium_sc_3_v2_index(self):
        """
        has_10x_indices: sample sheet with '*10' Chromium SC 3'v2 index
        """
        s = self._make_sample_sheet(
            self.sample_sheet_with_chromium_indices_A10)
        self.assertTrue(has_10x_indices(s))
    def test_sample_sheet_all_sc_atac_indices(self):
        """
        has_10x_indices: sample sheet with 10x scATAC-seq indices only
        """
        s = self._make_sample_sheet(
            self.sample_sheet_with_atac_indices)
        self.assertTrue(has_10x_indices(s))
    def test_sample_sheet_all_sc_atac_indices(self):
        """
        has_10x_indices: sample sheet with 10x Visium indices only
        """
        s = self._make_sample_sheet(
            self.sample_sheet_with_visium_indices)
        self.assertTrue(has_10x_indices(s))
    def test_sample_sheet_no_chromium_indices(self):
        """
        has_10x_indices: sample sheet with no Chromium SC indices
        """
        s = self._make_sample_sheet(
            self.sample_sheet_standard_indices)
        self.assertFalse(has_10x_indices(s))
    def test_sample_sheet_no_indices(self):
        """
        has_10x_indices: sample sheet with no indices
        """
        s = self._make_sample_sheet(
            self.sample_sheet_no_indices)
        self.assertFalse(has_10x_indices(s))
    def test_sample_sheet_some_chromium_sc_3_v2_indices(self):
        """
        has_10x_indices: sample sheet with some Chromium SC 3'v2 indices
        """
        s = self._make_sample_sheet(
            self.sample_sheet_mixed_indices)
        self.assertTrue(has_10x_indices(s))

class TestGetBasesMask10xAtac(unittest.TestCase):
    """
    Tests for the 'get_bases_mask_10x_atac' function
    """
    def setUp(self):
        # Create a temporary working dir
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        # Remove working dir
        if self.wd is not None:
            shutil.rmtree(self.wd)

    def test_get_bases_mask_10x_atac_I8_I16_dual_indexes(self):
        """get_bases_mask_10x_atac: run with I8,I16 dual indexes
        """
        # Make a single index RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.create("171020_NB500968_00002_AHGXXXX",
                                       "y76,I8,I16,y76",4,12))
        self.assertEqual(get_bases_mask_10x_atac(run_info_xml),
                         "y76,I8,y16,y76")

    def test_get_bases_mask_10x_atac_I16_I16_dual_indexes(self):
        """get_bases_mask_10x_atac: run with I16,I16 dual indexes
        """
        # Make a single index RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.create("171020_NB500968_00002_AHGXXXX",
                                       "y76,I16,I16,y76",4,12))
        self.assertEqual(get_bases_mask_10x_atac(run_info_xml),
                         "y76,I8nnnnnnnn,y16,y76")

    def test_get_bases_mask_10x_atac_too_short_index(self):
        """get_bases_mask_10x_atac: run with too-short index
        """
        # Make a single index RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.create("171020_NB500968_00002_AHGXXXX",
                                       "y76,I6,I16,y76",4,12))
        self.assertRaises(Exception,
                          get_bases_mask_10x_atac,
                          run_info_xml)

    def test_get_bases_mask_10x_atac_wrong_number_of_reads(self):
        """get_bases_mask_10x_atac: run with wrong number of reads
        """
        # Make a single index RunInfo.xml file
        run_info_xml = os.path.join(self.wd,"RunInfo.xml")
        with open(run_info_xml,'w') as fp:
            fp.write(RunInfoXml.create("171020_NB500968_00002_AHGXXXX",
                                       "y76,I8,y76",4,12))
        self.assertRaises(Exception,
                          get_bases_mask_10x_atac,
                          run_info_xml)

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
        os.chmod(cellranger_201,0o775)
        return cellranger_201
    def _make_mock_cellranger_501(self):
        # Make a fake cellranger 5.0.1 executable
        cellranger_501 = os.path.join(self.wd,"cellranger")
        with open(cellranger_501,'w') as fp:
            fp.write("#!/bin/bash\necho -n cellranger cellranger-5.0.1")
        os.chmod(cellranger_501,0o775)
        return cellranger_501
    def _make_mock_cellranger_600(self):
        # Make a fake cellranger 6.0.0 executable
        cellranger_600 = os.path.join(self.wd,"cellranger")
        with open(cellranger_600,'w') as fp:
            fp.write("#!/bin/bash\necho -n cellranger cellranger-6.0.0")
        os.chmod(cellranger_600,0o775)
        return cellranger_600
    def _make_mock_cellranger_atac_101(self):
        # Make a fake cellranger-atac 1.0.1 executable
        cellranger_atac_101 = os.path.join(self.wd,"cellranger-atac")
        with open(cellranger_atac_101,'w') as fp:
            fp.write("#!/bin/bash\ncat <<EOF\ncellranger-atac  (1.0.1)\nCopyright (c) 2018 10x Genomics, Inc.  All rights reserved.\n-------------------------------------------------------------------------------\n\nUsage:\n    cellranger-atac mkfastq\n\n    cellranger-atac count\n\n    cellranger-atac testrun\n    cellranger-atac upload\n    cellranger-atac sitecheckEOF")
        os.chmod(cellranger_atac_101,0o775)
        return cellranger_atac_101
    def _make_mock_cellranger_atac_200(self):
        # Make a fake cellranger-atac 2.0.0 executable
        cellranger_atac_200 = os.path.join(self.wd,"cellranger-atac")
        with open(cellranger_atac_200,'w') as fp:
            fp.write("#!/bin/bash\necho -n cellranger-atac cellranger-atac-2.0.0")
        os.chmod(cellranger_atac_200,0o775)
        return cellranger_atac_200
    def _make_mock_cellranger_arc_100(self):
        # Make a fake cellranger-atac 1.0.0 executable
        cellranger_arc_100 = os.path.join(self.wd,"cellranger-arc")
        with open(cellranger_arc_100,'w') as fp:
            fp.write("#!/bin/bash\necho -n cellranger-arc cellranger-arc-1.0.0")
        os.chmod(cellranger_arc_100,0o775)
        return cellranger_arc_100

    def _make_mock_cellranger_arc_200(self):
        # Make a fake cellranger-atac 1.0.0 executable
        cellranger_arc = os.path.join(self.wd,"cellranger-arc")
        with open(cellranger_arc,'w') as fp:
            fp.write("#!/bin/bash\necho cellranger-arc cellranger-arc-2.0.0")
        os.chmod(cellranger_arc,0o775)
        return cellranger_arc

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

    def test_cellranger_501(self):
        """cellranger_info: collect info for cellranger 5.0.1
        """
        cellranger = self._make_mock_cellranger_501()
        self.assertEqual(cellranger_info(path=cellranger),
                         (cellranger,'cellranger','5.0.1'))

    def test_cellranger_600(self):
        """cellranger_info: collect info for cellranger 6.0.0
        """
        cellranger = self._make_mock_cellranger_600()
        self.assertEqual(cellranger_info(path=cellranger),
                         (cellranger,'cellranger','6.0.0'))

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

    def test_cellranger_atac_200(self):
        """cellranger_info: collect info for cellranger-atac 2.0.0
        """
        cellranger_atac = self._make_mock_cellranger_atac_200()
        self.assertEqual(cellranger_info(path=cellranger_atac),
                         (cellranger_atac,'cellranger-atac','2.0.0'))

    def test_cellranger_atac_200_on_path(self):
        """cellranger_info: collect info for cellranger-atac 2.0.0 from PATH
        """
        os.environ['PATH'] = "%s%s%s" % (os.environ['PATH'],
                                         os.pathsep,
                                         self.wd)
        cellranger_atac = self._make_mock_cellranger_atac_200()
        self.assertEqual(cellranger_info(name='cellranger-atac'),
                         (cellranger_atac,'cellranger-atac','2.0.0'))

    def test_cellranger_arc_100(self):
        """cellranger_info: collect info for cellranger-arc 1.0.0
        """
        cellranger_arc = self._make_mock_cellranger_arc_100()
        self.assertEqual(cellranger_info(path=cellranger_arc),
                         (cellranger_arc,'cellranger-arc','1.0.0'))

    def test_cellranger_arc_100_on_path(self):
        """cellranger_info: collect info for cellranger-arc 1.0.0 from PATH
        """
        os.environ['PATH'] = "%s%s%s" % (os.environ['PATH'],
                                         os.pathsep,
                                         self.wd)
        cellranger_arc = self._make_mock_cellranger_arc_100()
        self.assertEqual(cellranger_info(name='cellranger-arc'),
                         (cellranger_arc,'cellranger-arc','1.0.0'))

    def test_cellranger_arc_200(self):
        """cellranger_info: collect info for cellranger-arc 2.0.0
        """
        cellranger_arc = self._make_mock_cellranger_arc_200()
        self.assertEqual(cellranger_info(path=cellranger_arc),
                         (cellranger_arc,'cellranger-arc','2.0.0'))

    def test_cellranger_arc_200_on_path(self):
        """cellranger_info: collect info for cellranger-arc 2.0.0 from PATH
        """
        os.environ['PATH'] = "%s%s%s" % (os.environ['PATH'],
                                         os.pathsep,
                                         self.wd)
        cellranger_arc = self._make_mock_cellranger_arc_200()
        self.assertEqual(cellranger_info(name='cellranger-arc'),
                         (cellranger_arc,'cellranger-arc','2.0.0'))

class TestSpacerangerInfo(unittest.TestCase):
    """
    Tests for the spaceranger_info function
    """
    def setUp(self):
        # Make temporary working dir
        self.wd = tempfile.mkdtemp(suffix="TestSpacerangerInfo")
        # Store the original state of PATH env var
        self.original_path = os.environ['PATH']
    def tearDown(self):
        # Restore the PATH env var
        os.environ['PATH'] = self.original_path
        # Remove temp dir
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)
    def _make_mock_spaceranger_110(self):
        # Make a fake spaceranger 1.0.1 executable
        spaceranger_110 = os.path.join(self.wd,"spaceranger")
        with open(spaceranger_110,'w') as fp:
            fp.write("#!/bin/bash\necho -n spaceranger 1.1.0")
        os.chmod(spaceranger_110,0o775)
        return spaceranger_110

    def test_spaceranger_110(self):
        """spaceranger_info: collect info for spaceranger 1.1.0
        """
        spaceranger = self._make_mock_spaceranger_110()
        self.assertEqual(spaceranger_info(path=spaceranger),
                         (spaceranger,'spaceranger','1.1.0'))

    def test_spaceranger_110_on_path(self):
        """spaceranger_info: collect info for spaceranger 1.1.0 from PATH
        """
        os.environ['PATH'] = "%s%s%s" % (os.environ['PATH'],
                                         os.pathsep,
                                         self.wd)
        spaceranger = self._make_mock_spaceranger_110()
        self.assertEqual(spaceranger_info(name='spaceranger'),
                         (spaceranger,'spaceranger','1.1.0'))

class TestGexSummary(unittest.TestCase):
    """
    Tests for the 'GexSummary' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestGexSummary')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)
            
    def test_gex_summary(self):
        """GexSummary: check metrics are extracted from CSV file
        """
        summary_csv = os.path.join(self.wd,"metrics_summary.csv")
        with open(summary_csv,'w') as fp:
            fp.write(METRICS_SUMMARY)
        m = GexSummary(summary_csv)
        self.assertEqual(m.estimated_number_of_cells,2272)
        self.assertEqual(m.mean_reads_per_cell,107875)
        self.assertEqual(m.median_genes_per_cell,1282)

class TestAtacSummary(unittest.TestCase):
    """
    Tests for the 'AtacSummary' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestAtacSummary')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_atac_summary_pre_2_0_0(self):
        """AtacSummary: get metrics from summary.csv (Cellranger ATAC pre-2.0.0)
        """
        summary_csv = os.path.join(self.wd,"summary.csv")
        with open(summary_csv,'w') as fp:
            fp.write(ATAC_SUMMARY)
        s = AtacSummary(summary_csv)
        self.assertEqual(s.version,"1.0.1")
        self.assertEqual(s.cells_detected,6748)
        self.assertEqual(s.annotated_cells,5682)
        self.assertEqual(s.median_fragments_per_cell,16119.5)
        self.assertEqual(s.frac_fragments_overlapping_targets,
                         0.575082094792)
        self.assertEqual(s.frac_fragments_overlapping_peaks,
                         0.556428090013)
        self.assertEqual(s.tss_enrichment_score,6.91438390781)
        self.assertRaises(AttributeError,
                          getattr,
                          s,'estimated_number_of_cells')

    def test_atac_summary_2_0_0(self):
        """AtacSummary: get metrics from summary.csv (Cellranger ATAC 2.0.0)
        """
        summary_csv = os.path.join(self.wd,"summary.csv")
        with open(summary_csv,'w') as fp:
            fp.write(ATAC_SUMMARY_2_0_0)
        s = AtacSummary(summary_csv)
        self.assertEqual(s.version,"2.0.0")
        self.assertEqual(s.estimated_number_of_cells,3582)
        self.assertEqual(s.median_fragments_per_cell,51354.5)
        self.assertEqual(s.frac_fragments_overlapping_peaks,0.4856)
        self.assertEqual(s.tss_enrichment_score,6.8333)
        self.assertRaises(AttributeError,
                          getattr,
                          s,'cells_detected')
        self.assertRaises(AttributeError,
                          getattr,
                          s,'annotated_cells')
        self.assertRaises(AttributeError,
                          getattr,
                          s,'frac_fragments_overlapping_targets')

class TestMultiomeSummary(unittest.TestCase):
    """
    Tests for the 'MultiomeSummary' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestMultiomeSummary')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_atac_summary_multiome_arc_1_0_0(self):
        """MultiomeSummary: check metrics are extracted from CSV file (Cellranger ARC 1.0.0)
        """
        summary_csv = os.path.join(self.wd,"summary.csv")
        with open(summary_csv,'wt') as fp:
            fp.write(MULTIOME_SUMMARY)
        s = MultiomeSummary(summary_csv)
        self.assertEqual(s.estimated_number_of_cells,744)
        self.assertEqual(
            s.atac_median_high_quality_fragments_per_cell,8079)
        self.assertEqual(s.gex_median_genes_per_cell,1490)

    def test_atac_summary_multiome_arc_2_0_0(self):
        """MultiomeSummary: check metrics are extracted from CSV file (Cellranger ARC 2.0.0)
        """
        summary_csv = os.path.join(self.wd,"summary.csv")
        with open(summary_csv,'wt') as fp:
            fp.write(MULTIOME_SUMMARY_2_0_0)
        s = MultiomeSummary(summary_csv)
        self.assertEqual(s.estimated_number_of_cells,785)
        self.assertEqual(
            s.atac_median_high_quality_fragments_per_cell,9.0)
        self.assertEqual(s.gex_median_genes_per_cell,15.0)

class TestMultiomeLibraries(unittest.TestCase):
    """
    Tests for the 'MultiomeLibraries' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestMultiomeLibraries')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_multiome_libraries(self):
        """MultiomeLibraries: read in data from a file
        """
        multiome_libraries_info = os.path.join(self.wd,
                                               "10x_multiome_libraries.info")
        with open(multiome_libraries_info,'wt') as fp:
            fp.write(MULTIOME_LIBRARIES)
        m = MultiomeLibraries(multiome_libraries_info)
        # Check list of samples
        self.assertEqual(m.local_samples,['PB1_ATAC','PB2_ATAC'])
        # Check lists of linked samples
        self.assertEqual(m.linked_samples('PB1_ATAC'),
                         ['NEXTSEQ_210111/12:PB_GEX/PB1_GEX',])
        self.assertEqual(m.linked_samples('PB2_ATAC'),
                         ['NEXTSEQ_210111/12:PB_GEX/PB2_GEX',])

    def test_multiome_libraries_linked_projects(self):
        """MultiomeLibraries: return linked projects
        """
        multiome_libraries_info = os.path.join(self.wd,
                                               "10x_multiome_libraries.info")
        with open(multiome_libraries_info,'wt') as fp:
            fp.write(MULTIOME_LIBRARIES)
        m = MultiomeLibraries(multiome_libraries_info)
        # Create run with linked project
        run = MockAnalysisDir("210111_NB01234_00012_ABXXXXX",
                              "nextseq",
                              project_metadata={
                                  'PB_GEX': {
                                      'Library type': 'GEX',
                                  },
                              },
                              lanes=(1,2,3,4),
                              top_dir=self.wd)
        run.add_fastq_batch('PB_GEX','PB1_GEX','PB1_GEX_S1',lanes=(1,2,3,4))
        run.add_fastq_batch('PB_GEX','PB2_GEX','PB2_GEX_S1',lanes=(1,2,3,4))
        run.create()
        # Check list of linked projects
        linked_projects = m.linked_projects()
        self.assertEqual(len(linked_projects),1)
        self.assertEqual(linked_projects[0].dirn,
                         os.path.join(self.wd,
                                      "210111_NB01234_00012_ABXXXXX_analysis",
                                      "PB_GEX"))

    def test_multiome_libraries_write_libraries_csv(self):
        """MultiomeLibraries: write 'libraries.csv' files
        """
        multiome_libraries_info = os.path.join(self.wd,
                                               "10x_multiome_libraries.info")
        with open(multiome_libraries_info,'wt') as fp:
            fp.write(MULTIOME_LIBRARIES)
        m = MultiomeLibraries(multiome_libraries_info)
        # Create run with linked project
        run = MockAnalysisDir("210111_NB01234_00012_ABXXXXX",
                              "nextseq",
                              project_metadata={
                                  'PB_GEX': {
                                      'Library type': 'GEX',
                                  },
                              },
                              lanes=(1,2,3,4),
                              top_dir=self.wd)
        run.add_fastq_batch('PB_GEX','PB1_GEX','PB1_GEX_S1',lanes=(1,2,3,4))
        run.add_fastq_batch('PB_GEX','PB2_GEX','PB2_GEX_S1',lanes=(1,2,3,4))
        run.create()
        # Write libraries.csv file
        libraries_csv = os.path.join(self.wd,"libraries.csv")
        self.assertFalse(os.path.exists(libraries_csv))
        m.write_libraries_csv('PB1_ATAC',
                              '/data/PB_ATAC/fastqs',
                              'ATAC',
                              filen=libraries_csv)
        self.assertTrue(os.path.exists(libraries_csv))
        with open(libraries_csv,'rt') as fp:
            contents = fp.read()
            self.assertEqual(contents,
                             """fastqs,sample,library_type
/data/PB_ATAC/fastqs,PB1_ATAC,Chromatin Accessibility
%s/210111_NB01234_00012_ABXXXXX_analysis/PB_GEX/fastqs,PB1_GEX,Gene Expression
""" % self.wd)

class TestSetCellCountForProject(unittest.TestCase):
    """
    Tests for the 'set_cell_count_for_project' function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestSetCellCountForProject')
    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)
    def _make_mock_analysis_project(self,single_cell_platform,library_type):
        # Create a mock AnalysisProject
        m = MockAnalysisProject('PJB',
                                fastq_names=("PJB1_S1_L001_R1_001.fastq.gz",
                                             "PJB1_S1_L001_R2_001.fastq.gz",),
                                metadata={'Single cell platform':
                                          single_cell_platform,
                                          'Library type': library_type,})
        m.create(top_dir=self.wd)
        return os.path.join(self.wd,'PJB')
    def test_set_cell_count_for_project(self):
        """
        set_cell_count_for_project: test for scRNA-seq
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "scRNA-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "5.0.1",
                                  "refdata-gex-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'wt') as fp:
            fp.write(METRICS_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-gex-GRCh38-2020-A
Cellranger version\t5.0.1
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)
    def test_set_cell_count_for_atac_project(self):
        """
        set_cell_count_for_project: test for scATAC-seq
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Single Cell ATAC",
            "scATAC-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "1.2.0",
                                  "refdata-cellranger-atac-GRCh38-1.2.0",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        summary_file = os.path.join(counts_dir,
                                    "summary.csv")
        with open(summary_file,'w') as fp:
            fp.write(ATAC_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-atac-GRCh38-1.2.0
Cellranger version\t1.2.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         5682)
    def test_set_cell_count_for_single_nuclei_atac_project(self):
        """
        set_cell_count_for_project: test for snATAC-seq
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Single Cell ATAC",
            "snATAC-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "1.2.0",
                                  "refdata-cellranger-atac-GRCh38-1.2.0",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        summary_file = os.path.join(counts_dir,
                                            "summary.csv")
        with open(summary_file,'w') as fp:
            fp.write(ATAC_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-atac-GRCh38-1.2.0
Cellranger version\t1.2.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         5682)
    def test_set_cell_count_for_atac_project_2_0_0(self):
        """
        set_cell_count_for_project: test for scATAC-seq (Cellranger ATAC 2.0.0)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Single Cell ATAC",
            "scATAC-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "2.0.0",
                                  "refdata-cellranger-atac-GRCh38-2020-A-2.0.0",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        summary_file = os.path.join(counts_dir,
                                    "summary.csv")
        with open(summary_file,'w') as fp:
            fp.write(ATAC_SUMMARY_2_0_0)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-atac-GRCh38-2020-A-2.0.0
Cellranger version\t2.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         3582)
    def test_set_cell_count_for_multiome_atac_project(self):
        """
        set_cell_count_for_project: test for single cell multiome ATAC
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Single Cell Multiome",
            "ATAC")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "1.0.0",
                                  "refdata-cellranger-arc-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        summary_file = os.path.join(counts_dir,
                                    "summary.csv")
        with open(summary_file,'w') as fp:
            fp.write(MULTIOME_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-arc-GRCh38-2020-A
Cellranger version\t1.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         744)
    def test_set_cell_count_for_multiome_gex_project(self):
        """
        set_cell_count_for_project: test for single cell multiome GEX
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Single Cell Multiome",
            "GEX")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "1.0.0",
                                  "refdata-cellranger-arc-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        summary_file = os.path.join(counts_dir,
                                    "summary.csv")
        with open(summary_file,'w') as fp:
            fp.write(MULTIOME_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-cellranger-arc-GRCh38-2020-A
Cellranger version\t1.0.0
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         744)
    def test_set_cell_count_project_missing_library_type(self):
        """
        set_cell_count_for_project: test for scRNA-seq when library not set
        """
        # Set up mock project with library type not set
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            None)
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "5.0.1",
                                  "refdata-gex-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'w') as fp:
            fp.write(METRICS_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-gex-GRCh38-2020-A
Cellranger version\t5.0.1
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)
    def test_set_cell_count_for_project_no_subdirs(self):
        """
        set_cell_count_for_project: test for scRNA-seq (old-style output)
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            "scRNA-seq")
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'w') as fp:
            fp.write(METRICS_SUMMARY)
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)
    def test_set_cell_count_project_missing_library_type_no_subdirs(self):
        """
        set_cell_count_for_project: test for scRNA-seq when library not set (old-style output)
        """
        # Set up mock project with library type not set
        project_dir = self._make_mock_analysis_project(
            "10xGenomics Chromium 3'v3",
            None)
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'w') as fp:
            fp.write(METRICS_SUMMARY)
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Update the cell counts
        print("Updating number of cells")
        set_cell_count_for_project(project_dir)
        # Check updated cell count
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         2272)
    def test_set_cell_count_fails_for_project_with_no_metadata(self):
        """
        set_cell_count_for_project: raises exception for project with no metadata
        """
        # Set up mock project
        project_dir = self._make_mock_analysis_project(None,None)
        # Add metrics_summary.csv
        counts_dir = os.path.join(project_dir,
                                  "qc",
                                  "cellranger_count",
                                  "5.0.1",
                                  "refdata-gex-GRCh38-2020-A",
                                  "PJB1",
                                  "outs")
        mkdirs(counts_dir)
        metrics_summary_file = os.path.join(counts_dir,
                                            "metrics_summary.csv")
        with open(metrics_summary_file,'wt') as fp:
            fp.write(METRICS_SUMMARY)
        # Add QC info file
        with open(os.path.join(project_dir,"qc","qc.info"),'wt') as fp:
            fp.write("""Cellranger reference datasets\t/data/refdata-gex-GRCh38-2020-A
Cellranger version\t5.0.1
""")
        # Check initial cell count
        print("Checking number of cells")
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
        # Attempting to update the cell counts should raise
        # NotImplementedError
        self.assertRaises(NotImplementedError,
                          set_cell_count_for_project,
                          project_dir)
        # Check cell count wasn't updated
        self.assertEqual(AnalysisProject("PJB1",
                                         project_dir).info.number_of_cells,
                         None)
