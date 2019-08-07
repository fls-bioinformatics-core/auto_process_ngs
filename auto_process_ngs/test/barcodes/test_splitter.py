#######################################################################
# Tests for barcodes/splitter.py module
#######################################################################

import unittest
import os
import shutil
import tempfile

from auto_process_ngs.barcodes.splitter import HammingMetrics
from auto_process_ngs.barcodes.splitter import BarcodeMatcher
from auto_process_ngs.barcodes.splitter import split_single_end
from auto_process_ngs.barcodes.splitter import split_paired_end

class TestHammingMetrics(unittest.TestCase):
    def test_hamming_exact(self):
        self.assertEqual(HammingMetrics.hamming_distance('AGGTCTA','AGGTCTA'),0)
    def test_hamming_one_mismatch(self):
        self.assertEqual(HammingMetrics.hamming_distance('AGGTCTA','ACGTCTA'),1)
    def test_hamming_two_mismatches(self):
        self.assertEqual(HammingMetrics.hamming_distance('AGGTCTA','ACCTCTA'),2)

class TestHammingWithN(unittest.TestCase):
    def test_hamming_with_N_exact(self):
        self.assertEqual(HammingMetrics.hamming_distance_with_N('AGGTCTA','AGGTCTA'),0)
    def test_hamming_with_N_one_mismatch(self):
        self.assertEqual(HammingMetrics.hamming_distance_with_N('AGGTCTA','ACGTCTA'),1)
        self.assertEqual(HammingMetrics.hamming_distance_with_N('ACNTCTA','ACGTCTA'),1)
        self.assertEqual(HammingMetrics.hamming_distance_with_N('ACNTCTA','ACNTCTA'),1)
    def test_hamming_with_N_two_mismatches(self):
        self.assertEqual(HammingMetrics.hamming_distance_with_N('AGGTCTA','ACCTCTA'),2)
        self.assertEqual(HammingMetrics.hamming_distance_with_N('ANNTCTA','ACCTCTA'),2)
        self.assertEqual(HammingMetrics.hamming_distance_with_N('ANNTCTA','ANNTCTA'),2)
    def test_hamming_different_lengths(self):
        self.assertRaises(ValueError,HammingMetrics.hamming_distance_with_N,
                          'AGGTCTAA','AGGTCTA')

class TestBarcodeMatcher(unittest.TestCase):
    def test_barcodematcher_exact_match(self):
        b = BarcodeMatcher(('TTGCTA','AGGTCT'))
        self.assertEqual(b.match('AGGTCT'),'AGGTCT')
        self.assertEqual(b.match('TTGCTA'),'TTGCTA')
        self.assertEqual(b.match('CGGTCT'),None)
        self.assertEqual(b.match('CGTTCT'),None)
    def test_barcodematcher_one_mismatch(self):
        b = BarcodeMatcher(('TTGCTA','AGGTCT'),max_dist=1)
        self.assertEqual(b.match('AGGTCT'),'AGGTCT')
        self.assertEqual(b.match('TTGCTA'),'TTGCTA')
        self.assertEqual(b.match('CGGTCT'),'AGGTCT')
        self.assertEqual(b.match('CGTTCT'),None)
    def test_barcodematcher_two_mismatches(self):
        b = BarcodeMatcher(('TTGCTA','AGGTCT'),max_dist=2)
        self.assertEqual(b.match('AGGTCT'),'AGGTCT')
        self.assertEqual(b.match('TTGCTA'),'TTGCTA')
        self.assertEqual(b.match('CGGTCT'),'AGGTCT')
        self.assertEqual(b.match('CGTTCT'),'AGGTCT')
    def test_barcodematcher_list_sequences(self):
        b = BarcodeMatcher(('TTGCTA','AGGTCT','GCCTAT'))
        self.assertEqual(b.sequences,['AGGTCT','GCCTAT','TTGCTA'])
    def test_barcodematcher_ambigiuous_sequences(self):
        self.assertRaises(Exception,BarcodeMatcher,('AGGTCTA','AGGTCTA'))
        self.assertRaises(Exception,BarcodeMatcher,('AGGTCTC','AGGTCTA'),max_dist=1)
        self.assertRaises(Exception,BarcodeMatcher,('AGGTCCC','AGGTCTA'),max_dist=2)


fastq_r1 = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
TTTACAACTAGCTTCTCTTTTTCTT
+
>AA?131@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 1:N:0:TAAGGCGA
TCCCCAGTCTCAGCCCTACTCCACT
+
11>>>11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 1:N:0:GCCTTACC
CCACCACGCCTGGCTAATTTTTTTT
+
1>1>>AAAAAFAGGGGBFGGGGGG0
@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:TAAGNNNA
CAGCAATATACACTTCACTCTGCAT
+
111>>1FBFFFFGGGGCBGEG3A3D
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:TAAGGCGA
GATAAAGACAGAGTCTTAATTAAAC
+
11>1>11DFFCFFDGGGB3BF313A
"""

fastq_r2 = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 2:N:0:TAAGGCGA
CTTCTTCTTTTTTTTTTTGTCTATC
+
1>11>13BBD111AAE00/B2B222
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 2:N:0:TAAGGCGA
TTTTTTTTCTTTTTTTTTTCCTTTT
+
11>>1>110B3BF1A00/A01B22B
@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 2:N:0:GCCTTACC
TTTCTATTCTCCTTCTCATAAAAAA
+
111113333B311B1B133331110
@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 2:N:0:TAAGNNNA
TTTTTTTTTCTTCTATTCTCAGATG
+
1>1>111100BB3B3B22B222121
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 2:N:0:TAAGGCGA
TTTTTTTTCCATCTTTTTTAATTTT
+
>1>111>1013B1BF3BFE01BB22
"""

class TestSplitSingleEnd(unittest.TestCase):
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_split_single_end')
        # Test file
        self.fastq = os.path.join(self.wd,'test.fq')
        open(self.fastq,'w').write(fastq_r1)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_split_single_end(self):
        matcher = BarcodeMatcher(('TAAGGCGA','GCCTTACC'))
        split_single_end(matcher,(self.fastq,),output_dir=self.wd)
        self.assertTrue(os.path.exists(os.path.join(self.wd,'TAAGGCGA.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'GCCTTACC.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'undetermined.fastq')))
        self.assertEqual(open(os.path.join(self.wd,'TAAGGCGA.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
TTTACAACTAGCTTCTCTTTTTCTT
+
>AA?131@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 1:N:0:TAAGGCGA
TCCCCAGTCTCAGCCCTACTCCACT
+
11>>>11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:TAAGGCGA
GATAAAGACAGAGTCTTAATTAAAC
+
11>1>11DFFCFFDGGGB3BF313A
""")
        self.assertEqual(open(os.path.join(self.wd,'GCCTTACC.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 1:N:0:GCCTTACC
CCACCACGCCTGGCTAATTTTTTTT
+
1>1>>AAAAAFAGGGGBFGGGGGG0
""")
        self.assertEqual(open(os.path.join(self.wd,'undetermined.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:TAAGNNNA
CAGCAATATACACTTCACTCTGCAT
+
111>>1FBFFFFGGGGCBGEG3A3D
""")

class TestSplitPairedEnd(unittest.TestCase):
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_split_paired_end')
        # Test file
        self.fastq_r1 = os.path.join(self.wd,'test_r1.fq')
        self.fastq_r2 = os.path.join(self.wd,'test_r2.fq')
        open(self.fastq_r1,'w').write(fastq_r1)
        open(self.fastq_r2,'w').write(fastq_r2)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_split_paired_end(self):
        matcher = BarcodeMatcher(('TAAGGCGA','GCCTTACC'))
        split_paired_end(matcher,((self.fastq_r1,self.fastq_r2),),output_dir=self.wd)
        self.assertTrue(os.path.exists(os.path.join(self.wd,'TAAGGCGA_R1.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'TAAGGCGA_R2.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'GCCTTACC_R1.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'GCCTTACC_R2.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'undetermined_R1.fastq')))
        self.assertTrue(os.path.exists(os.path.join(self.wd,'undetermined_R2.fastq')))
        self.assertEqual(open(os.path.join(self.wd,'TAAGGCGA_R1.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
TTTACAACTAGCTTCTCTTTTTCTT
+
>AA?131@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 1:N:0:TAAGGCGA
TCCCCAGTCTCAGCCCTACTCCACT
+
11>>>11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:TAAGGCGA
GATAAAGACAGAGTCTTAATTAAAC
+
11>1>11DFFCFFDGGGB3BF313A
""")
        self.assertEqual(open(os.path.join(self.wd,'TAAGGCGA_R2.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 2:N:0:TAAGGCGA
CTTCTTCTTTTTTTTTTTGTCTATC
+
1>11>13BBD111AAE00/B2B222
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 2:N:0:TAAGGCGA
TTTTTTTTCTTTTTTTTTTCCTTTT
+
11>>1>110B3BF1A00/A01B22B
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 2:N:0:TAAGGCGA
TTTTTTTTCCATCTTTTTTAATTTT
+
>1>111>1013B1BF3BFE01BB22
""")
        self.assertEqual(open(os.path.join(self.wd,'GCCTTACC_R1.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 1:N:0:GCCTTACC
CCACCACGCCTGGCTAATTTTTTTT
+
1>1>>AAAAAFAGGGGBFGGGGGG0
""")
        self.assertEqual(open(os.path.join(self.wd,'GCCTTACC_R2.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 2:N:0:GCCTTACC
TTTCTATTCTCCTTCTCATAAAAAA
+
111113333B311B1B133331110
""")
        self.assertEqual(open(os.path.join(self.wd,'undetermined_R1.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:TAAGNNNA
CAGCAATATACACTTCACTCTGCAT
+
111>>1FBFFFFGGGGCBGEG3A3D
""")
        self.assertEqual(open(os.path.join(self.wd,'undetermined_R2.fastq'),'r').read(),
                         """@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 2:N:0:TAAGNNNA
TTTTTTTTTCTTCTATTCTCAGATG
+
1>1>111100BB3B3B22B222121
""")
