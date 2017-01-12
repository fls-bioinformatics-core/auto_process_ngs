#######################################################################
# Tests for stats.py module
#######################################################################

import os
import unittest
import tempfile
import shutil
import gzip
from auto_process_ngs.stats import FastqReadCounter

# FastqReadCounter
fastq_data = """@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
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
fastq_multi_lane_data = """@NB500968:10:H57NTAFXX:1:11101:4210:1091 1:N:0:CGGCAGAA
CTCCAGTTCTGAGTAACTTCAAGGG
+
AAAAAEEEEEEEEE6/EEEA//EE/
@NB500968:10:H57NTAFXX:1:11101:21113:1097 1:N:0:AGGCAGAA
CCCACACACCTATGTTTATTGCAGC
+
AAAAAEEEEEEEEEEEEEEEEEEEE
@NB500968:10:H57NTAFXX:1:11101:20639:1107 1:N:0:AGGCAGAA
GGATGATGGACCCGGAGCACATAAA
+
AAAAAEEEEEEEEEEEEEEEEEEEE
@NB500968:10:H57NTAFXX:1:11101:11994:1111 1:N:0:AGGCAGAA
GATTGTATGTCTCAAGAGATGCAGG
+
AAAAAEEEEEEEEEEEEEEEEEEEE
@NB500968:10:H57NTAFXX:2:11101:16822:1113 1:N:0:AGGCAGAA
TAAGGTGGAGTGGGTTTGGGGCTAG
+
6AAAAAEEEEEEEEAEEEEEEEAEE
@NB500968:10:H57NTAFXX:3:11101:7366:1119 1:N:0:AGGCAGAA
CCCCTCACCTGTGCTTTCTCGTCCT
+
AAAAAEAEEEEEEEEEEEEEEEEEE
@NB500968:10:H57NTAFXX:3:11101:2547:1129 1:N:0:AGGAAGAA
AGACTGAGCCGAATTGGTAAATAGT
+
AAAAAE/EEE/EEEE6/AAEEEAEE
@NB500968:10:H57NTAFXX:4:11101:12127:1131 1:N:0:AGGCAGAA
GCCCCAACTGTCTGGCGGGGAGGAC
+
AAAAAEEEEEEEEEEEEEEEEEEEE
@NB500968:10:H57NTAFXX:4:11101:11578:1134 1:N:0:AGGCAGAA
CACATCTACAACGTTATCGTCACAG
+
AAAAAEAEE/EEAEEEEEAEE/EEE
@NB500968:10:H57NTAFXX:4:11101:2817:1161 1:N:0:AGGCAGAA
GACCAACCCTGGGGTTAGTATAGCT
+
A/AAAEEEEEEEEEEE<EEEEEEEA
@NB500968:10:H57NTAFXX:1:11101:18892:1167 1:N:0:AGGCAGAA
CTTCTGTGGAACGAGGGTTTATTTT
+
AAAAAEEEEEEEEEEEEEEEEEEEE
@NB500968:10:H57NTAFXX:4:11101:14759:1191 1:N:0:AGGCAGAA
CTTATACACATCTCCGAGCCCACGA
+
AAAAAEEEEEEEEEEEEEEEE/EE/
"""
class TestFastqReadCounter(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None

    def _make_working_dir(self):
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_FastqReadCounter')

    def _make_fastq(self,name,contents):
        # Create a FASTQ file under the working directory
        # called "name" and populated with "contents"
        # Working directory will be created if not already
        # set up
        # If "name" ends with ".gz" then "contents" will be
        # gzipped
        # Returns the path to the file
        self._make_working_dir()
        filen = os.path.join(self.wd,name)
        if filen.endswith(".gz"):
            with gzip.GzipFile(filen,'wb') as fp:
                fp.write(contents)
        else:
            with open(filen,'w') as fp:
                fp.write(contents)
        return filen

    def tearDown(self):
        # Remove temporary working dir
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)

    def test_simple(self):
        readcounter = FastqReadCounter.simple
        fq = self._make_fastq("test_S1_L001_R1_001.fastq",
                              fastq_data)
        self.assertEqual(readcounter(fq),5)
        fq = self._make_fastq("test_S2_R1_001.fastq",
                              fastq_multi_lane_data)
        self.assertEqual(readcounter(fq),12)

    def test_simple_gz(self):
        readcounter = FastqReadCounter.simple
        fq = self._make_fastq("test_S1_L001_R1_001.fastq.gz",
                              fastq_data)
        self.assertEqual(readcounter(fq),5)
        fq = self._make_fastq("test_S2_R1_001.fastq.gz",
                              fastq_multi_lane_data)
        self.assertEqual(readcounter(fq),12)

    def test_fastqiterator(self):
        readcounter = FastqReadCounter.fastqiterator
        fq = self._make_fastq("test_S1_L001_R1_001.fastq",
                              fastq_data)
        self.assertEqual(readcounter(fq),5)
        fq = self._make_fastq("test_S2_R1_001.fastq",
                              fastq_multi_lane_data)
        self.assertEqual(readcounter(fq),12)

    def test_fastqiterator_gz(self):
        readcounter = FastqReadCounter.fastqiterator
        fq = self._make_fastq("test_S1_L001_R1_001.fastq.gz",
                              fastq_data)
        self.assertEqual(readcounter(fq),5)
        fq = self._make_fastq("test_S2_R1_001.fastq.gz",
                              fastq_multi_lane_data)
        self.assertEqual(readcounter(fq),12)

    def test_zcat_wc(self):
        readcounter = FastqReadCounter.zcat_wc
        fq = self._make_fastq("test_S1_L001_R1_001.fastq",
                              fastq_data)
        self.assertEqual(readcounter(fq),5)
        fq = self._make_fastq("test_S2_R1_001.fastq",
                              fastq_multi_lane_data)
        self.assertEqual(readcounter(fq),12)

    def test_zcat_wc_gz(self):
        readcounter = FastqReadCounter.zcat_wc
        fq = self._make_fastq("test_S1_L001_R1_001.fastq.gz",
                              fastq_data)
        self.assertEqual(readcounter(fq),5)
        fq = self._make_fastq("test_S2_R1_001.fastq.gz",
                              fastq_multi_lane_data)
        self.assertEqual(readcounter(fq),12)

    def test_reads_per_lane(self):
        readcounter = FastqReadCounter.reads_per_lane
        fq = self._make_fastq("test_S1_L001_R1_001.fastq",
                              fastq_data)
        self.assertEqual(readcounter(fq),{ 1: 5 })
        fq = self._make_fastq("test_S2_R1_001.fastq",
                              fastq_multi_lane_data)
        self.assertEqual(readcounter(fq),{ 1: 5,
                                           2: 1,
                                           3: 2,
                                           4: 4 })

    def test_reads_per_lane_gz(self):
        readcounter = FastqReadCounter.reads_per_lane
        fq = self._make_fastq("test_S1_L001_R1_001.fastq.gz",
                              fastq_data)
        self.assertEqual(readcounter(fq),{ 1: 5 })
        fq = self._make_fastq("test_S2_R1_001.fastq.gz",
                              fastq_multi_lane_data)
        self.assertEqual(readcounter(fq),{ 1: 5,
                                           2: 1,
                                           3: 2,
                                           4: 4 })
