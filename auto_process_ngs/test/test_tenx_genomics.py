#######################################################################
# Tests for tenx_genomics_utils.py module
#######################################################################

import unittest
import tempfile
import os
import shutil
from auto_process_ngs.tenx_genomics_utils import *

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

class TestHasChromiumIndices(unittest.TestCase):
    """
    Tests for the 'has_chromium_indices' function
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
        self.wd = tempfile.mkdtemp(suffix="TestHasChromiumIndices")
    def tearDown(self):
        # Remove temp dir
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
        self.assertTrue(has_chromium_indices(s))
    def test_sample_sheet_no_chromium_indices(self):
        """
        has_chromium_indices: sample sheet with no Chromium indices
        """
        s = self._make_sample_sheet(
            self.sample_sheet_standard_indices)
        self.assertFalse(has_chromium_indices(s))
    def test_sample_sheet_some_chromium_sc_3_v2_indices(self):
        """
        has_chromium_indices: sample sheet with some Chromium SC 3'v2 indices
        """
        s = self._make_sample_sheet(
            self.sample_sheet_mixed_indices)
        self.assertTrue(has_chromium_indices(s))

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
        shutil.rmtree(self.wd)
    def _make_mock_cellranger_201(self):
        # Make a fake cellranger 2.0.1 executable
        cellranger_201 = os.path.join(self.wd,"cellranger")
        with open(cellranger_201,'w') as fp:
            fp.write("#!/bin/bash\ncat <<EOF\ncellranger $1  (2.0.1)\nCopyright (c) 2017 10x Genomics, Inc.  All rights reserved.\n-------------------------------------------------------------------------------\n\nUsage:\n    cellranger mkfastq\n\n    cellranger count\n    cellranger aggr\n    cellranger reanalyze\n    cellranger mkloupe\n    cellranger mat2csv\n\n    cellranger mkgtf\n    cellranger mkref\n\n    cellranger vdj\n\n    cellranger mkvdjref\n\n    cellranger testrun\n    cellranger upload\n    cellranger sitecheckEOF")
        os.chmod(cellranger_201,0775)
        return cellranger_201
            
    def test_cellranger_201(self):
        """cellranger_info: collect info for cellranger 2.0.1
        """
        cellranger = self._make_mock_cellranger_201()
        self.assertEqual(cellranger_info(path=cellranger),
                         (cellranger,'cellranger','2.0.1'))
