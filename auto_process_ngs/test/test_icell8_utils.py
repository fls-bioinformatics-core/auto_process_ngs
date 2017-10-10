#######################################################################
# Tests for icell8_utils.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
from bcftbx.FASTQFile import FastqRead
from auto_process_ngs.icell8_utils import ICell8WellList
from auto_process_ngs.icell8_utils import ICell8Read1
from auto_process_ngs.icell8_utils import ICell8ReadPair
from auto_process_ngs.icell8_utils import ICell8FastqIterator
from auto_process_ngs.icell8_utils import ICell8Stats
from auto_process_ngs.icell8_utils import ICell8FastqIterator
from auto_process_ngs.icell8_utils import normalize_sample_name
from auto_process_ngs.icell8_utils import get_icell8_bases_mask

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

icell8_read_pair = {
    'r1': """@NB500968:70:HCYMKBGX2:1:11101:22672:1659 1:N:0:1
GTTCCTGATTAAGTCAAGTGCTGGGG
+
AAAAAEEEEEEEEEEEEEEEE//6//
""",
    'r2': """@NB500968:70:HCYMKBGX2:1:11101:22672:1659 2:N:0:1
CCCATGAGACTTAAGAATCACACAGACCTTGGACTTTCCTGATTTCACGGGACGCTGCTCTGAGAGTGAAATTGGGCCTTCTGTAAATATGTGAAGTGTGGTTTCTTTTCAAACCTTATATGGCCCTGCA
+
AAAAAEAEEEEEEEE/EEEEEEEEE/EEEEEEEEEEEEEEAAEAEEAEEEE<EEEEAEEEEEEEAEEE<EEEEEAAEEEEEEAAEEEAAEAE<AEEA/</<AEE<AEEEEA/6AAEE/AEE/AAEEA<A<
"""
}

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
        """ICell8WellList: list samples from data file
        """
        self.assertEqual(ICell8WellList(self.well_list).samples(),
                         ['ESC2','d1.2'])
    def test_icell8welllist_sample(self):
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
        self.assertRaises(KeyError,well_list.sample,'NNNNNNNNNNN')

# ICell8Read1
class TestICell8Read1(unittest.TestCase):
    """Tests for the ICell8Read1 class
    """
    def _fastqread(self,s):
        # Return populated FastqRead object given a multiline
        # string
        return FastqRead(*s.rstrip('\n').split('\n'))
    def test_icell_read1_init(self):
        """ICell8Read1: create R1 read
        """
        r1 = ICell8Read1(self._fastqread(icell8_read_pair['r1']))
        self.assertEqual(r1.read,self._fastqread(icell8_read_pair['r1']))
    def test_icell8_read1_barcode(self):
        """ICell8Read1: get barcode
        """
        r1 = ICell8Read1(self._fastqread(icell8_read_pair['r1']))
        self.assertEqual(r1.barcode,"GTTCCTGATTA")
        self.assertEqual(r1.barcode_quality,"AAAAAEEEEEE")
    def test_icell8_read1_umi(self):
        """ICell8Read1: get UMI
        """
        r1 = ICell8Read1(self._fastqread(icell8_read_pair['r1']))
        self.assertEqual(r1.umi,"AGTCAAGTGCTGGG")
        self.assertEqual(r1.umi_quality,"EEEEEEEEEE//6/")
    def test_icell8_read1_min_barcode_quality(self):
        """ICell8Read1: get minimum quality score for barcode
        """
        r1 = ICell8Read1(self._fastqread(icell8_read_pair['r1']))
        self.assertEqual(r1.min_barcode_quality,'A')
    def test_icell8_read1_min_umi_quality(self):
        """ICell8Read1: get minimum quality score for UMI
        """
        r1 = ICell8Read1(self._fastqread(icell8_read_pair['r1']))
        self.assertEqual(r1.min_umi_quality,'/')

# ICell8ReadPair
class TestICell8ReadPair(unittest.TestCase):
    """Tests for the ICell8ReadPair class
    """
    def _fastqread(self,s):
        # Return populated FastqRead object given a multiline
        # string
        return FastqRead(*s.rstrip('\n').split('\n'))
    def test_icell_read_pair_init(self):
        """ICell8ReadPair: create read pair
        """
        pair = ICell8ReadPair(self._fastqread(icell8_read_pair['r1']),
                              self._fastqread(icell8_read_pair['r2']))
        self.assertEqual(pair.r1,self._fastqread(icell8_read_pair['r1']))
        self.assertEqual(pair.r2,self._fastqread(icell8_read_pair['r2']))
    def test_icell_read_pair_init_bad_pair(self):
        """ICell8ReadPair: fail when reads are not a pair
        """
        alt_r1 = """@NB500968:70:HCYMKBGX2:1:11101:24365:2047 1:N:0:1
AGAAGAGTACCTGGAAAATGTTGGCG
+
6AAAAEEEEEEEEEEAEEEEE/<6//
"""
        self.assertRaises(Exception,
                          ICell8ReadPair,
                          self._fastqread(alt_r1),
                          self._fastqread(icell8_read_pair['r2']))
    def test_icell8_read_pair_barcode(self):
        """ICell8ReadPair: get barcode
        """
        pair = ICell8ReadPair(self._fastqread(icell8_read_pair['r1']),
                              self._fastqread(icell8_read_pair['r2']))
        self.assertEqual(pair.barcode,"GTTCCTGATTA")
        self.assertEqual(pair.barcode_quality,"AAAAAEEEEEE")
    def test_icell8_read_pair_umi(self):
        """ICell8ReadPair: get UMI
        """
        pair = ICell8ReadPair(self._fastqread(icell8_read_pair['r1']),
                              self._fastqread(icell8_read_pair['r2']))
        self.assertEqual(pair.umi,"AGTCAAGTGCTGGG")
        self.assertEqual(pair.umi_quality,"EEEEEEEEEE//6/")
    def test_icell8_read_pair_min_barcode_quality(self):
        """ICell8ReadPair: get minimum quality score for barcode
        """
        pair = ICell8ReadPair(self._fastqread(icell8_read_pair['r1']),
                              self._fastqread(icell8_read_pair['r2']))
        self.assertEqual(pair.min_barcode_quality,'A')
    def test_icell8_read_pair_min_umi_quality(self):
        """ICell8ReadPair: get minimum quality score for UMI
        """
        pair = ICell8ReadPair(self._fastqread(icell8_read_pair['r1']),
                              self._fastqread(icell8_read_pair['r2']))
        self.assertEqual(pair.min_umi_quality,'/')

# ICell8FastqIterator
icell8_fastq_r1 = """@NB500968:70:HCYMKBGX2:1:11101:22672:1659 1:N:0:1
GTTCCTGATTAAGTCAAGTGCTGGGG
+
AAAAAEEEEEEEEEEEEEEEE//6//
@NB500968:70:HCYMKBGX2:1:11101:24365:2047 1:N:0:1
AGAAGAGTACCTGGAAAATGTTGGCG
+
6AAAAEEEEEEEEEEAEEEEE/<6//
@NB500968:70:HCYMKBGX2:1:11101:24470:3201 1:N:0:1
GTCTGCAACGCGGAGGCCGGATCGCG
+
AAAAAEEEEEEEEEEEEEEEE/////
"""

icell8_fastq_r2 = """@NB500968:70:HCYMKBGX2:1:11101:22672:1659 2:N:0:1
CCCATGAGACTTAAGAATCACACAGACCTTGGACTTTCCTGATTTCACGGGACGCTGCTCTGAGAGTGAAATTGGGCCTTCTGTAAATATGTGAAGTGTGGTTTCTTTTCAAACCTTATATGGCCCTGCA
+
AAAAAEAEEEEEEEE/EEEEEEEEE/EEEEEEEEEEEEEEAAEAEEAEEEE<EEEEAEEEEEEEAEEE<EEEEEAAEEEEEEAAEEEAAEAE<AEEA/</<AEE<AEEEEA/6AAEE/AEE/AAEEA<A<
@NB500968:70:HCYMKBGX2:1:11101:24365:2047 2:N:0:1
CTACACGACGCGGGAACCGGGCGTGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGCGGATTTCTGAGTTCGAGGCCAGCCTGGTCTACAAAGTGAGTTCCAGGAAAGCCAGGGCTA
+
AAAAAEEEE6AEEEEEEA/EEEEEEEEEAEEE6EEEE/EEEEEEEEAEEAEEEEEAEEEEEE/<<E<EEEEEEAEAE//AEEAA<AE/EEAEEEEEEAEEAA/AEAEAE/EEA<<EE/<A<EE<A//<A/
@NB500968:70:HCYMKBGX2:1:11101:24470:3201 2:N:0:1
TTCCCTACACGACGCGGGGGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAGAGGAAAGAGGAAAGAAAGAGG
+
AAAAAEEEEEEEEEE///EA/EEEEEEEEEAEEEEEEEEEEEEE<EEE/EEEEEEAEEE/E<EEEEEAE6AA<AE//<///////////////////////////////////////////////////6
"""

class TestICell8FastqIterator(unittest.TestCase):
    """Tests for the ICell8FastqIterator class
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.ICell8FastqIterator')
        # Test files
        self.r1 = os.path.join(self.wd,'icell8.r1.fq')
        with open(self.r1,'w') as fp:
            fp.write(icell8_fastq_r1)
        self.r2 = os.path.join(self.wd,'icell8.r2.fq')
        with open(self.r2,'w') as fp:
            fp.write(icell8_fastq_r2)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_icell8fastqiterator_over_pairs(self):
        """ICell8FastqIterator: iterate over read pairs
        """
        fqr1_data = ""
        fqr2_data = ""
        for i,pair in enumerate(ICell8FastqIterator(self.r1,self.r2)):
            self.assertTrue(isinstance(pair.r1,FastqRead))
            self.assertTrue(isinstance(pair.r2,FastqRead))
            n = i*4
            self.assertEqual(str(pair.r1),
                             '\n'.join(icell8_fastq_r1.split('\n')[n:n+4]))
            self.assertEqual(str(pair.r2),
                             '\n'.join(icell8_fastq_r2.split('\n')[n:n+4]))
            fqr1_data += "%s\n" % pair.r1
            fqr2_data += "%s\n" % pair.r2
        self.assertEqual(fqr1_data,icell8_fastq_r1)
        self.assertEqual(fqr2_data,icell8_fastq_r2)

class TestICell8Stats(unittest.TestCase):
    """Tests for the ICell8Stats class
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.ICell8Stats')
        # Test files
        self.r1 = os.path.join(self.wd,'icell8.r1.fq')
        with open(self.r1,'w') as fp:
            fp.write(icell8_fastq_r1)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_icell8stats(self):
        """ICell8Stats: collect stats from Icell8 R1 FASTQ
        """
        fastqs = (self.r1,)
        stats = ICell8Stats(*fastqs)
        self.assertEqual(stats.nreads(),3)
        self.assertEqual(stats.barcodes(),['AGAAGAGTACC',
                                           'GTCTGCAACGC',
                                           'GTTCCTGATTA',])
        self.assertEqual(stats.nreads('GTTCCTGATTA'),1)
        self.assertEqual(stats.nreads('AGAAGAGTACC'),1)
        self.assertEqual(stats.nreads('GTCTGCAACGC'),1)
        self.assertEqual(stats.distinct_umis(),['AGTCAAGTGCTGGG',
                                                'GGAGGCCGGATCGC',
                                                'TGGAAAATGTTGGC'])
        self.assertEqual(stats.distinct_umis('GTTCCTGATTA'),
                         ['AGTCAAGTGCTGGG'])
        self.assertEqual(stats.distinct_umis('AGAAGAGTACC'),
                         ['TGGAAAATGTTGGC'])
        self.assertEqual(stats.distinct_umis('GTCTGCAACGC'),
                         ['GGAGGCCGGATCGC'])
    def test_icell8stats_multicore(self):
        """ICell8Stats: collect stats from Icell8 R1 FASTQ (multicore)
        """
        fastqs = (self.r1,self.r1,)
        stats = ICell8Stats(*fastqs,nprocs=2)
        self.assertEqual(stats.nreads(),6)
        self.assertEqual(stats.barcodes(),['AGAAGAGTACC',
                                           'GTCTGCAACGC',
                                           'GTTCCTGATTA',])
        self.assertEqual(stats.nreads('GTTCCTGATTA'),2)
        self.assertEqual(stats.nreads('AGAAGAGTACC'),2)
        self.assertEqual(stats.nreads('GTCTGCAACGC'),2)
        self.assertEqual(stats.distinct_umis(),['AGTCAAGTGCTGGG',
                                                'GGAGGCCGGATCGC',
                                                'TGGAAAATGTTGGC'])
        self.assertEqual(stats.distinct_umis('GTTCCTGATTA'),
                         ['AGTCAAGTGCTGGG'])
        self.assertEqual(stats.distinct_umis('AGAAGAGTACC'),
                         ['TGGAAAATGTTGGC'])
        self.assertEqual(stats.distinct_umis('GTCTGCAACGC'),
                         ['GGAGGCCGGATCGC'])

class TestNormalizeSampleNameFunction(unittest.TestCase):
    """
    Tests for the normalize_sample_name function
    """
    def test_normalize_permissible_names(self):
        """
        normalize_sample_name: works for permissible names
        """
        self.assertEqual(normalize_sample_name("ESC1"),"ESC1")
        self.assertEqual(normalize_sample_name("d1.2"),"d1.2")
    def test_normalize_names_with_illegal_chars(self):
        """
        normalize_sample_name: works for names with illegal characters
        """
        self.assertEqual(normalize_sample_name("Neg Ctrl"),"Neg_Ctrl")
        self.assertEqual(normalize_sample_name("S1/S2"),"S1_S2")

class TestGetIcell8BasesMaskFunction(unittest.TestCase):
    """
    Tests for the get_icell8_bases_mask function
    """
    def test_get_icell8_bases_mask(self):
        """
        get_icell8_bases_mask: reset bases mask
        """
        self.assertEqual(get_icell8_bases_mask("y101,I7,y101"),
                         "y25n76,I7,y101")
        self.assertEqual(get_icell8_bases_mask("y250,I8,I8,y250"),
                         "y25n225,I8,I8,y250")
        self.assertEqual(get_icell8_bases_mask("y25,I8,I8,y250"),
                         "y25,I8,I8,y250")
