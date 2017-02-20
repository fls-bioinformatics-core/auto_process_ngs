#######################################################################
# Tests for icell8_utils.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.icell8_utils import ICell8WellList

well_list_data = """Row	Col	Candidate	For dispense	Sample	Barcode	State	Cells1	Cells2	Signal1	Signal2	Size1	Size2	Integ Signal1	Integ Signal2	Circularity1	Circularity2	Confidence	Confidence1	Confidence2	Dispense tip	Drop index	Global drop index	Source well	Sequencing count	Image1	Image2
0	4	True	True	ESC2	AACCTTCCTTA	Good	1	0	444		55		24420		0.9805677		1	1	1	1	4	5	A1	Pos0_Hoechst_A01.tif	Pos0_TexasRed_A01.tif
0	6	True	True	ESC2	AACGAACGCTC	Good	1	0	251		21		5271		1		0.8972501	0.8972501	1	1	5	7	A1		Pos1_Hoechst_A02.tif	Pos1_TexasRed_A02.tif
0	20	True	True	d1.2	AACCAATCGTC	Good	1	0	298		36		10728		1		1	1	1	2	9	12	A2		Pos3_Hoechst_A04.tif	Pos3_TexasRed_A04.tif
0	21	True	True	d1.2	AACCAACGCAA	Good	1	0	389		45		17505		1		1	1	1	2	10	13	A2		Pos3_Hoechst_A04.tif	Pos3_TexasRed_A04.tif
0	29	True	True	ESC2	AACGCCAAGAC	Good	1	0	377		67		25259		0.9272023		0.9970149	0.9970149	1	1	6	6A1		Pos4_Hoechst_A05.tif	Pos4_TexasRed_A05.tif
0	36	True	True	d1.2	AACCGCCTAAC	Good	1	0	261		24		6264		1		1	1	1	2	1	4	A2		Pos6_Hoechst_A07.tif	Pos6_TexasRed_A07.tif
0	43	True	True	d1.2	AACCGCGCTCA	Good	1	0	287		38		10906		1		0.97	0.97	1	2	8	11	A2		Pos7_Hoechst_A08.tif	Pos7_TexasRed_A08.tif
0	54	True	True	ESC2	AACCTTGCAAG	Good	1	0	290		35		10150		0.9533346		1	1	1	1	7	7	A1	Pos9_Hoechst_A10.tif	Pos9_TexasRed_A10.tif
0	63	True	True	d1.2	AACCAACCGCA	Good	1	0	247		24		5928		1		1	1	1	2	4	7	A2		Pos10_Hoechst_A11.tif	Pos10_TexasRed_A11.tif
"""

# ICell8WellList
class TestICell8WellList(unittest.TestCase):
    """Tests for the ICell8WellList class
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.ICell8WellList')
        # Test file
        self.well_list = os.path.join(self.wd,'test.fq')
        with open(self.well_list,'w') as fp:
            fp.write(well_list_data)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_icell8welllist_barcodes(self):
        """ICell8WellList: list barcodes from data file
        """
        self.assertEqual(ICell8WellList(self.well_list).barcodes(),
                         ['AACCTTCCTTA',
                          'AACGAACGCTC',
                          'AACCAATCGTC',
                          'AACCAACGCAA',
                          'AACGCCAAGAC',
                          'AACCGCCTAAC',
                          'AACCGCGCTCA',
                          'AACCTTGCAAG',
                          'AACCAACCGCA'])
    def test_icell8welllist_samples(self):
        """ICell8WellList: get sample (=cell type) for barcode
        """
        well_list = ICell8WellList(self.well_list)
        self.assertEqual(well_list.sample('AACCTTCCTTA'),'ESC2')
        self.assertEqual(well_list.sample('AACGAACGCTC'),'ESC2')
        self.assertEqual(well_list.sample('AACCAATCGTC'),'d1.2')
        self.assertEqual(well_list.sample('AACCAACGCAA'),'d1.2')
        self.assertEqual(well_list.sample('AACGCCAAGAC'),'ESC2')
        self.assertEqual(well_list.sample('AACCGCCTAAC'),'d1.2')
        self.assertEqual(well_list.sample('AACCGCGCTCA'),'d1.2')
        self.assertEqual(well_list.sample('AACCTTGCAAG'),'ESC2')
        self.assertEqual(well_list.sample('AACCAACCGCA'),'d1.2')
