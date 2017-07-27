#######################################################################
# Tests for fastq_utils.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.fastq_utils import assign_barcodes_single_end
from auto_process_ngs.fastq_utils import pair_fastqs

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

fastq_r1_out = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TTTAC
AACTAGCTTCTCTTTTTCTT
+
31@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 1:N:0:TCCCC
AGTCTCAGCCCTACTCCACT
+
11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 1:N:0:CCACC
ACGCCTGGCTAATTTTTTTT
+
AAAAAFAGGGGBFGGGGGG0
@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:CAGCA
ATATACACTTCACTCTGCAT
+
1FBFFFFGGGGCBGEG3A3D
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:GATAA
AGACAGAGTCTTAATTAAAC
+
11DFFCFFDGGGB3BF313A
"""

# assign_barcodes_single_end
class TestAssignBarcodesSingleEnd(unittest.TestCase):
    """Tests for the assign_barcodes_single_end function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_split_single_end')
        # Test file
        self.fastq_in = os.path.join(self.wd,'test.fq')
        open(self.fastq_in,'w').write(fastq_r1)
        # Output file
        self.fastq_out = os.path.join(self.wd,'out.fq')
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_assign_barcodes_single_end(self):
        """assign_barcodes_single_end: extract inline barcodes from first 5 bases
        """
        nreads = assign_barcodes_single_end(self.fastq_in,
                                            self.fastq_out)
        self.assertEqual(nreads,5)
        self.assertEqual(open(self.fastq_out,'r').read(),
                         fastq_r1_out)

# pair_fastqs
fastq1_r1 = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
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
"""
fastq1_r2 = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 2:N:0:TAAGGCGA
TTTACAACTAGCTTCTCTTTTTCTT
+
>AA?131@C1FCGGGG1BFFGF1F3
@MISEQ:34:000000000-A7PHP:1:1101:16449:1793 2:N:0:TAAGGCGA
TCCCCAGTCTCAGCCCTACTCCACT
+
11>>>11CDFFFBCGGG1AFE11AB
@MISEQ:34:000000000-A7PHP:1:1101:15171:1796 2:N:0:GCCTTACC
CCACCACGCCTGGCTAATTTTTTTT
+
1>1>>AAAAAFAGGGGBFGGGGGG0
"""
fastq2_r1 = """@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 1:N:0:TAAGNNNA
CAGCAATATACACTTCACTCTGCAT
+
111>>1FBFFFFGGGGCBGEG3A3D
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 1:N:0:TAAGGCGA
GATAAAGACAGAGTCTTAATTAAAC
+
11>1>11DFFCFFDGGGB3BF313A
"""
fastq2_r2 = """@MISEQ:34:000000000-A7PHP:1:1101:18777:1797 2:N:0:TAAGNNNA
CAGCAATATACACTTCACTCTGCAT
+
111>>1FBFFFFGGGGCBGEG3A3D
@MISEQ:34:000000000-A7PHP:1:1101:18622:1812 2:N:0:TAAGGCGA
GATAAAGACAGAGTCTTAATTAAAC
+
11>1>11DFFCFFDGGGB3BF313A
"""
class TestPairFastqs(unittest.TestCase):
    """Tests for the pair_fastqs function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_split_single_end')
        # Test files
        self.fastq1_r1 = os.path.join(self.wd,'test1_r1.fq')
        with open(self.fastq1_r1,'w') as fp:
            fp.write(fastq1_r1)
        self.fastq1_r2 = os.path.join(self.wd,'test1_r2.fq')
        with open(self.fastq1_r2,'w') as fp:
            fp.write(fastq1_r2)
        self.fastq2_r1 = os.path.join(self.wd,'test2_r1.fq')
        with open(self.fastq2_r1,'w') as fp:
            fp.write(fastq2_r1)
        self.fastq2_r2 = os.path.join(self.wd,'test2_r2.fq')
        with open(self.fastq2_r2,'w') as fp:
            fp.write(fastq2_r2)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_pair_fastqs(self):
        """pair_fastqs: pair up a set of FASTQs
        """
        fastqs = [self.fastq2_r2,
                  self.fastq1_r1,
                  self.fastq2_r1,
                  self.fastq1_r2,]
        self.assertEqual(pair_fastqs(fastqs),
                         ([(self.fastq1_r1,self.fastq1_r2),
                           (self.fastq2_r1,self.fastq2_r2)],
                         []))
    def test_pair_fastqs_unpaired(self):
        """pair_fastqs: handle paired and unpaired FASTQs
        """
        fastqs = [self.fastq2_r2,
                  self.fastq1_r1,
                  self.fastq2_r1,]
        self.assertEqual(pair_fastqs(fastqs),
                         ([(self.fastq2_r1,self.fastq2_r2),],
                         [self.fastq1_r1,]))
    def test_pair_fastqs_no_pairs(self):
        """pair_fastqs: handle set of FASTQ with no pairs
        """
        fastqs = [self.fastq2_r2,
                  self.fastq1_r1,]
        self.assertEqual(pair_fastqs(fastqs),
                         ([],[self.fastq1_r1,self.fastq2_r2]))
