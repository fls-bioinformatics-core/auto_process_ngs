#######################################################################
# Tests for tenx_genomics/utils.py module
#######################################################################

import unittest
import tempfile
import os
import shutil
from bcftbx.mock import RunInfoXml
from auto_process_ngs.tenx.utils import *

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
    def _make_mock_cellranger_710(self):
        # Make a fake cellranger 7.1.0 executable
        cellranger_710 = os.path.join(self.wd,"cellranger")
        with open(cellranger_710,'w') as fp:
            fp.write("#!/bin/bash\necho cellranger cellranger-7.1.0")
        os.chmod(cellranger_710,0o775)
        return cellranger_710
    def _make_mock_cellranger_800(self):
        # Make a fake cellranger 8.0.0 executable
        cellranger_800 = os.path.join(self.wd,"cellranger")
        with open(cellranger_800,'w') as fp:
            fp.write("#!/bin/bash\necho cellranger cellranger-8.0.0")
        os.chmod(cellranger_800,0o775)
        return cellranger_800
    def _make_mock_cellranger_900(self):
        # Make a fake cellranger 9.0.0 executable
        cellranger_900 = os.path.join(self.wd,"cellranger")
        with open(cellranger_900,'w') as fp:
            fp.write("#!/bin/bash\necho cellranger cellranger-9.0.0")
        os.chmod(cellranger_900,0o775)
        return cellranger_900
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

    def test_cellranger_710(self):
        """cellranger_info: collect info for cellranger 7.1.0
        """
        cellranger = self._make_mock_cellranger_710()
        self.assertEqual(cellranger_info(path=cellranger),
                         (cellranger,'cellranger','7.1.0'))

    def test_cellranger_800(self):
        """cellranger_info: collect info for cellranger 8.0.0
        """
        cellranger = self._make_mock_cellranger_800()
        self.assertEqual(cellranger_info(path=cellranger),
                         (cellranger,'cellranger','8.0.0'))

    def test_cellranger_900(self):
        """cellranger_info: collect info for cellranger 8.0.0
        """
        cellranger = self._make_mock_cellranger_900()
        self.assertEqual(cellranger_info(path=cellranger),
                         (cellranger,'cellranger','9.0.0'))

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
    def _make_mock_spaceranger_131(self):
        # Make a fake spaceranger 1.3.1 executable
        spaceranger_131 = os.path.join(self.wd,"spaceranger")
        with open(spaceranger_131,'w') as fp:
            fp.write("#!/bin/bash\necho -n spaceranger spaceranger-1.3.1")
        os.chmod(spaceranger_131,0o775)
        return spaceranger_131
    def _make_mock_spaceranger_211(self):
        # Make a fake spaceranger 2.1.1 executable
        spaceranger_211 = os.path.join(self.wd,"spaceranger")
        with open(spaceranger_211,'w') as fp:
            fp.write("#!/bin/bash\necho -n spaceranger spaceranger-2.1.1")
        os.chmod(spaceranger_211,0o775)
        return spaceranger_211
    def _make_mock_spaceranger_300(self):
        # Make a fake spaceranger 2.1.1 executable
        spaceranger_300 = os.path.join(self.wd,"spaceranger")
        with open(spaceranger_300,'w') as fp:
            fp.write("#!/bin/bash\necho -n spaceranger spaceranger-3.0.0")
        os.chmod(spaceranger_300,0o775)
        return spaceranger_300

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

    def test_spaceranger_131(self):
        """spaceranger_info: collect info for spaceranger 1.3.1
        """
        spaceranger = self._make_mock_spaceranger_131()
        self.assertEqual(spaceranger_info(path=spaceranger),
                         (spaceranger,'spaceranger','1.3.1'))

    def test_spaceranger_131_on_path(self):
        """spaceranger_info: collect info for spaceranger 1.3.1 from PATH
        """
        os.environ['PATH'] = "%s%s%s" % (os.environ['PATH'],
                                         os.pathsep,
                                         self.wd)
        spaceranger = self._make_mock_spaceranger_131()
        self.assertEqual(spaceranger_info(name='spaceranger'),
                         (spaceranger,'spaceranger','1.3.1'))

    def test_spaceranger_211(self):
        """spaceranger_info: collect info for spaceranger 2.1.1
        """
        spaceranger = self._make_mock_spaceranger_211()
        self.assertEqual(spaceranger_info(path=spaceranger),
                         (spaceranger,'spaceranger','2.1.1'))

    def test_spaceranger_211_on_path(self):
        """spaceranger_info: collect info for spaceranger 2.1.1 from PATH
        """
        os.environ['PATH'] = "%s%s%s" % (os.environ['PATH'],
                                         os.pathsep,
                                         self.wd)
        spaceranger = self._make_mock_spaceranger_211()
        self.assertEqual(spaceranger_info(name='spaceranger'),
                         (spaceranger,'spaceranger','2.1.1'))

    def test_spaceranger_300(self):
        """spaceranger_info: collect info for spaceranger 3.0.0
        """
        spaceranger = self._make_mock_spaceranger_300()
        self.assertEqual(spaceranger_info(path=spaceranger),
                         (spaceranger,'spaceranger','3.0.0'))

    def test_spaceranger_300_on_path(self):
        """spaceranger_info: collect info for spaceranger 3.0.0 from PATH
        """
        os.environ['PATH'] = "%s%s%s" % (os.environ['PATH'],
                                         os.pathsep,
                                         self.wd)
        spaceranger = self._make_mock_spaceranger_300()
        self.assertEqual(spaceranger_info(name='spaceranger'),
                         (spaceranger,'spaceranger','3.0.0'))

class TestMakeMultiConfigTemplate(unittest.TestCase):
    """
    Tests for the make_multi_config_template function
    """
    def setUp(self):
        # Make temporary working dir
        self.wd = tempfile.mkdtemp(suffix="TestMakeMultiConfigTemplate")

    def tearDown(self):
        # Remove temp dir
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_make_multi_config_template_default(self):
        """
        make_multi_config_template: check default template
        """
        expected_content = """[gene-expression]
reference,/path/to/transcriptome
#force-cells,n
create-bam,true
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate

[samples]
sample_id,cmo_ids,description
MULTIPLEXED_SAMPLE,CMO1|CMO2|...,DESCRIPTION
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_cellplex_710(self):
        """
        make_multi_config_template: check CellPlex template (7.1.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
#force-cells,n
#no-bam,true|false
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_CML,/runs/novaseq_50/fastqs,any,PJB_CML,[Gene Expression|Multiplexing Capture],
PJB_GEX,/runs/novaseq_50/fastqs,any,PJB_GEX,[Gene Expression|Multiplexing Capture],

[samples]
sample_id,cmo_ids,description
MULTIPLEXED_SAMPLE,CMO1|CMO2|...,DESCRIPTION
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_CML","PJB_GEX"),
                                   library_type="CellPlex",
                                   cellranger_version="7.0.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_cellplex_800(self):
        """
        make_multi_config_template: check CellPlex template (8.0.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
#force-cells,n
create-bam,true
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_CML,/runs/novaseq_50/fastqs,any,PJB_CML,[Gene Expression|Multiplexing Capture],
PJB_GEX,/runs/novaseq_50/fastqs,any,PJB_GEX,[Gene Expression|Multiplexing Capture],

[samples]
sample_id,cmo_ids,description
MULTIPLEXED_SAMPLE,CMO1|CMO2|...,DESCRIPTION
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_CML","PJB_GEX"),
                                   library_type="CellPlex",
                                   cellranger_version="8.0.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_cellplex_900(self):
        """
        make_multi_config_template: check CellPlex template (9.0.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
#force-cells,n
create-bam,true
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_CML,/runs/novaseq_50/fastqs,any,PJB_CML,[Gene Expression|Multiplexing Capture],
PJB_GEX,/runs/novaseq_50/fastqs,any,PJB_GEX,[Gene Expression|Multiplexing Capture],

[samples]
sample_id,cmo_ids,description
MULTIPLEXED_SAMPLE,CMO1|CMO2|...,DESCRIPTION
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_CML","PJB_GEX"),
                                   library_type="CellPlex",
                                   cellranger_version="9.0.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_cellplex_scrnaseq_800(self):
        """
        make_multi_config_template: check CellPlex scRNA-seq template (8.0.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
#force-cells,n
create-bam,true
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_CML,/runs/novaseq_50/fastqs,any,PJB_CML,[Gene Expression|Multiplexing Capture],
PJB_GEX,/runs/novaseq_50/fastqs,any,PJB_GEX,[Gene Expression|Multiplexing Capture],

[samples]
sample_id,cmo_ids,description
MULTIPLEXED_SAMPLE,CMO1|CMO2|...,DESCRIPTION
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_CML","PJB_GEX"),
                                   library_type="CellPlex scRNA-seq",
                                   cellranger_version="8.0.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_cellplex_scrnaseq_900(self):
        """
        make_multi_config_template: check CellPlex scRNA-seq template (9.0.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
#force-cells,n
create-bam,true
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_CML,/runs/novaseq_50/fastqs,any,PJB_CML,[Gene Expression|Multiplexing Capture],
PJB_GEX,/runs/novaseq_50/fastqs,any,PJB_GEX,[Gene Expression|Multiplexing Capture],

[samples]
sample_id,cmo_ids,description
MULTIPLEXED_SAMPLE,CMO1|CMO2|...,DESCRIPTION
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_CML","PJB_GEX"),
                                   library_type="CellPlex scRNA-seq",
                                   cellranger_version="9.0.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_flex_700(self):
        """
        make_multi_config_template: check Flex template (7.1.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
probe-set,/data/mm10_probe_set.csv
#force-cells,n
no-bam,true
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_Flex,/runs/novaseq_50/fastqs,any,PJB_Flex,[Gene Expression|Antibody Capture],

[samples]
sample_id,probe_barcode_ids,description
MULTIPLEXED_SAMPLE,BC001|BC002|...,DESCRIPTION
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   probe_set="/data/mm10_probe_set.csv",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_Flex",),
                                   no_bam=True,
                                   library_type="Flex",
                                   cellranger_version="7.1.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_flex_800(self):
        """
        make_multi_config_template: check Flex template (8.0.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
probe-set,/data/mm10_probe_set.csv
#force-cells,n
create-bam,false
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_Flex,/runs/novaseq_50/fastqs,any,PJB_Flex,[Gene Expression|Antibody Capture],

[samples]
sample_id,probe_barcode_ids,description
MULTIPLEXED_SAMPLE,BC001|BC002|...,DESCRIPTION
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   probe_set="/data/mm10_probe_set.csv",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_Flex",),
                                   no_bam=True,
                                   library_type="Flex",
                                   cellranger_version="8.0.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_flex_900(self):
        """
        make_multi_config_template: check Flex template (9.0.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
probe-set,/data/mm10_probe_set.csv
#force-cells,n
create-bam,false
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_Flex,/runs/novaseq_50/fastqs,any,PJB_Flex,[Gene Expression|Antibody Capture],

[samples]
sample_id,probe_barcode_ids,description
MULTIPLEXED_SAMPLE,BC001|BC002|...,DESCRIPTION
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   probe_set="/data/mm10_probe_set.csv",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_Flex",),
                                   no_bam=True,
                                   library_type="Flex",
                                   cellranger_version="9.0.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_immune_profiling_710(self):
        """
        make_multi_config_template: check Single Cell Immune Profiling template (7.0.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
#force-cells,n
#no-bam,true|false
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

#[vdj]
#reference,/path/to/vdj/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_CML,/runs/novaseq_50/fastqs,any,PJB_CML,[Gene Expression|Antibody Capture|VDJ-B|VDJ-T],
PJB_GEX,/runs/novaseq_50/fastqs,any,PJB_GEX,[Gene Expression|Antibody Capture|VDJ-B|VDJ-T],
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_CML","PJB_GEX"),
                                   library_type="Single Cell Immune Profiling",
                                   cellranger_version="7.1.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_immune_profiling_800(self):
        """
        make_multi_config_template: check Single Cell Immune Profiling template (8.0.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
#force-cells,n
create-bam,true
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

#[vdj]
#reference,/path/to/vdj/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_CML,/runs/novaseq_50/fastqs,any,PJB_CML,[Gene Expression|Antibody Capture|VDJ-B|VDJ-T],
PJB_GEX,/runs/novaseq_50/fastqs,any,PJB_GEX,[Gene Expression|Antibody Capture|VDJ-B|VDJ-T],
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_CML","PJB_GEX"),
                                   library_type="Single Cell Immune Profiling",
                                   cellranger_version="8.0.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_immune_profiling_900(self):
        """
        make_multi_config_template: check Single Cell Immune Profiling template (9.0.0)
        """
        expected_content = """[gene-expression]
reference,/data/mm10_transcriptome
#force-cells,n
create-bam,true
#cmo-set,/path/to/custom/cmo/reference

#[feature]
#reference,/path/to/feature/reference

#[vdj]
#reference,/path/to/vdj/reference

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
PJB_CML,/runs/novaseq_50/fastqs,any,PJB_CML,[Gene Expression|Antibody Capture|VDJ-B|VDJ-T],
PJB_GEX,/runs/novaseq_50/fastqs,any,PJB_GEX,[Gene Expression|Antibody Capture|VDJ-B|VDJ-T],
"""
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        make_multi_config_template(out_file,
                                   reference="/data/mm10_transcriptome",
                                   fastq_dir="/runs/novaseq_50/fastqs",
                                   samples=("PJB_CML","PJB_GEX"),
                                   library_type="Single Cell Immune Profiling",
                                   cellranger_version="9.0.0")
        self.assertTrue(os.path.exists(out_file))
        with open(out_file,'rt') as fp:
            actual_content = '\n'.join([line for line in fp.read().split('\n')
                                        if not line.startswith('##')])
        self.assertEqual(expected_content,actual_content)

    def test_make_multi_config_template_unsupported_library(self):
        """
        make_multi_config_template: exception for unsupported library type
        """
        out_file = os.path.join(self.wd,"10x_multi_config.csv")
        self.assertRaises(Exception,
                          make_multi_config_template,
                          out_file,
                          reference="/data/mm10_transcriptome",
                          fastq_dir="/runs/novaseq_50/fastqs",
                          samples=("PJB_CML","PJB_GEX"),
                          library_type="Superduper scRNA-seq",
                          cellranger_version="8.0.0")
        self.assertFalse(os.path.exists(out_file))
