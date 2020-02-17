#######################################################################
# Tests for fastq_utils.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
import gzip
from auto_process_ngs.fastq_utils import IlluminaFastqAttrs
from auto_process_ngs.fastq_utils import FastqReadCounter
from auto_process_ngs.fastq_utils import assign_barcodes_single_end
from auto_process_ngs.fastq_utils import get_read_number
from auto_process_ngs.fastq_utils import get_read_count
from auto_process_ngs.fastq_utils import pair_fastqs
from auto_process_ngs.fastq_utils import pair_fastqs_by_name
from auto_process_ngs.fastq_utils import group_fastqs_by_name
from auto_process_ngs.fastq_utils import remove_index_fastqs

# Test data
fastq_data = u"""@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
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
fastq_multi_lane_data = u"""@NB500968:10:H57NTAFXX:1:11101:4210:1091 1:N:0:CGGCAGAA
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

# IlluminaFastqAttrs
class TestIlluminaFastqAttrs(unittest.TestCase):
    """Tests for the IlluminaFastqAttrs class
    """
    def test_full_name(self):
        """IlluminaFastqAttrs: full Illumina-style fastq name
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,
                         'NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')

    def test_full_name_dual_index(self):
        """IlluminaFastqAttrs: full Illumina-style fastq name with dual index
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_ACAGTG-GTTCAC_L003_R2_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,
                         'NH1_ChIP-seq_Gli1_ACAGTG-GTTCAC_L003_R2_001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG-GTTCAC')
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG-GTTCAC_L003_R2_001')

    def test_full_name_dual_index_with_plus(self):
        """IlluminaFastqAttrs: full Illumina-style fastq name with dual index using '+'
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_ACAGTG+GTTCAC_L003_R2_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,
                         'NH1_ChIP-seq_Gli1_ACAGTG+GTTCAC_L003_R2_001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG+GTTCAC')
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG+GTTCAC_L003_R2_001')

    def test_full_name_blc2fastq2(self):
        """IlluminaFastqAttrs: Illumina fastq name from bcl2fastq2
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_S4_L003_R2_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,
                         'NH1_ChIP-seq_Gli1_S4_L003_R2_001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,4)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_S4_L003_R2_001')

    def test_index_read_blc2fastq2(self):
        """IlluminaFastqAttrs: Illumina index read fastq name from bcl2fastq2
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_S4_L003_I1_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.sample_number,4)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,1)
        self.assertEqual(fq.set_number,1)
        self.assertTrue(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_S4_L003_I1_001')

    def test_name_no_lane_blc2fastq2(self):
        """IlluminaFastqAttrs: Illumina fastq name from bcl2fastq2 (without lane)
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_S4_R2_001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,
                         'NH1_ChIP-seq_Gli1_S4_R2_001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,4)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_S4_R2_001')

    def test_name_only(self):
        """IlluminaFastqAttrs: reduced fastq name (sample name only)
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,None)
        self.assertEqual(fq.set_number,None)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1')

    def test_name_only_paired_end(self):
        """IlluminaFastqAttrs: reduced fastq name (sample name only, paired end)
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_R2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_R2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_R2')

    def test_name_and_lane(self):
        """IlluminaFastqAttrs: reduced fastq name (sample name and lane)
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_L001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_L001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,1)
        self.assertEqual(fq.read_number,None)
        self.assertEqual(fq.set_number,None)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_L001')

    def test_name_and_lane_paired_end(self):
        """IlluminaFastqAttrs: reduced fastq name (sample name and lane, paired end)
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_L001_R2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_L001_R2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,1)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_L001_R2')

    def test_name_and_tag(self):
        """IlluminaFastqAttrs: reduced fastq name (sample name and barcode)
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_ACAGTG')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_ACAGTG')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,None)
        self.assertEqual(fq.set_number,None)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG')

    def test_name_and_tag_paired_end(self):
        """IlluminaFastqAttrs: reduced fastq name (sample name and barcode, paired end)
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_ACAGTG_R2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_ACAGTG_R2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_R2')

    def test_name_tag_and_lane(self):
        """IlluminaFastqAttrs: reduced fastq name (sample name, barcode and lane)
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_ACAGTG_L001')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_ACAGTG_L001')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,1)
        self.assertEqual(fq.read_number,None)
        self.assertEqual(fq.set_number,None)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_L001')

    def test_name_tag_and_lane_paired_end(self):
        """IlluminaFastqAttrs: reduced fastq name (sample name, barcode and lane, paired end)
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq_Gli1_ACAGTG_L001_R2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_ACAGTG_L001_R2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,1)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_L001_R2')

    def test_AGTC_sample_names(self):
        """IlluminaFastqAttrs: sample names consisting of letters 'A', 'G', 'T' and 'C'
        """
        for name in ('A','G','T','C','AGCT'):
            fq = IlluminaFastqAttrs('%s_R1' % name)
            self.assertEqual(fq.sample_name,name)
            self.assertEqual(fq.sample_number,None)
            self.assertEqual(fq.barcode_sequence,None)
            self.assertEqual(fq.lane_number,None)
            self.assertEqual(fq.read_number,1)
            self.assertEqual(fq.set_number,None)
            self.assertFalse(fq.is_index_read)
            self.assertEqual(str(fq),'%s_R1' % name)

    def test_non_standard_sample_name(self):
        """IlluminaFastqAttrs: non-standard Fastq names with sample name only
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq.r2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq')
        self.assertEqual(fq.basename,'NH1_ChIP-seq.r2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq.r2')

    def test_non_standard_sample_name_with_dots(self):
        """IlluminaFastqAttrs: non-standard Fastq names with sample name containing dots
        """
        fq = IlluminaFastqAttrs('NH1.2.r2')
        self.assertEqual(fq.sample_name,'NH1.2')
        self.assertEqual(fq.basename,'NH1.2.r2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,None)
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,None)
        self.assertEqual(str(fq),'NH1.2.r2')

    def test_non_standard_sample_name_and_barcode(self):
        """IlluminaFastqAttrs: non-standard Fastq names with sample name and barcode
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq.ACAGTG.r2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq')
        self.assertEqual(fq.basename,'NH1_ChIP-seq.ACAGTG.r2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)

    def test_non_standard_sample_name_and_dual_barcode(self):
        """IlluminaFastqAttrs: non-standard Fastq names with sample name and dual barcode
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq.ACAGTG-GTTCAC.r2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq')
        self.assertEqual(fq.basename,'NH1_ChIP-seq.ACAGTG-GTTCAC.r2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG-GTTCAC')
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)

    def test_non_standard_sample_name_and_dual_barcode_with_plus(self):
        """IlluminaFastqAttrs: non-standard Fastq names with sample name and dual barcode using '+'
        """
        fq = IlluminaFastqAttrs('NH1_ChIP-seq.ACAGTG+GTTCAC.r2')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq')
        self.assertEqual(fq.basename,'NH1_ChIP-seq.ACAGTG+GTTCAC.r2')
        self.assertEqual(fq.extension,'')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG+GTTCAC')
        self.assertEqual(fq.lane_number,None)
        self.assertEqual(fq.read_number,2)

    def test_input_is_full_path(self):
        """IlluminaFastqAttrs: input as full path to Fastq file
        """
        fq = IlluminaFastqAttrs('/data/Project_NH/Sample_NH1/NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001.fastq.gz')
        self.assertEqual(fq.sample_name,'NH1_ChIP-seq_Gli1')
        self.assertEqual(fq.basename,'NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')
        self.assertEqual(fq.extension,'.fastq.gz')
        self.assertEqual(fq.sample_number,None)
        self.assertEqual(fq.barcode_sequence,'ACAGTG')
        self.assertEqual(fq.lane_number,3)
        self.assertEqual(fq.read_number,2)
        self.assertEqual(fq.set_number,1)
        self.assertFalse(fq.is_index_read)
        self.assertEqual(str(fq),'NH1_ChIP-seq_Gli1_ACAGTG_L003_R2_001')

# FastqReadCounter
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

# assign_barcodes_single_end
fastq_r1 = u"""@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TAAGGCGA
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

fastq_r1_out = u"""@MISEQ:34:000000000-A7PHP:1:1101:12552:1774 1:N:0:TTTAC
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
class TestAssignBarcodesSingleEnd(unittest.TestCase):
    """Tests for the assign_barcodes_single_end function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_assign_barcodes_single_end')
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
        self.empty_r1 = os.path.join(self.wd,'empty_r1.fq')
        with open(self.empty_r1,'w') as fp:
            fp.write('')
        self.empty_r2 = os.path.join(self.wd,'empty_r2.fq')
        with open(self.empty_r2,'w') as fp:
            fp.write('')
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
    def test_pair_fastqs_empty_files(self):
        """pair_fastqs: handle set of FASTQs with 'empty' pairs
        """
        fastqs = [self.empty_r1,
                  self.empty_r2]
        self.assertEqual(pair_fastqs(fastqs),
                         ([],[self.empty_r1,self.empty_r2]))

# pair_fastqs_by_name
class TestPairFastqsByNameFunction(unittest.TestCase):
    """
    Tests for the pair_fastqs_by_name function
    """
    def test_pair_fastqs_by_name_SE(self):
        """
        pair_fastqs_by_name: single-end fastqs
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(pair_fastqs_by_name(fastqs),
                         [("/data/PJB1_S1_R1_001.fastq.gz",),
                          ("/data/PJB2_S2_R1_001.fastq.gz",)])

    def test_pair_fastqs_by_name_PE(self):
        """
        pair_fastqs_by_name: paired-end fastqs
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB1_S1_R2_001.fastq.gz",
                  "/data/PJB2_S2_R1_001.fastq.gz",
                  "/data/PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(pair_fastqs_by_name(fastqs),
                         [("/data/PJB1_S1_R1_001.fastq.gz",
                           "/data/PJB1_S1_R2_001.fastq.gz"),
                          ("/data/PJB2_S2_R1_001.fastq.gz",
                           "/data/PJB2_S2_R2_001.fastq.gz")])

    def test_pair_fastqs_by_name_PE_with_index_read(self):
        """
        pair_fastqs_by_name: paired-end fastqs with index reads
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB1_S1_R2_001.fastq.gz",
                  "/data/PJB1_S1_I1_001.fastq.gz",
                  "/data/PJB2_S2_R1_001.fastq.gz",
                  "/data/PJB2_S2_R2_001.fastq.gz",
                  "/data/PJB2_S2_I1_001.fastq.gz")
        self.assertEqual(pair_fastqs_by_name(fastqs),
                         [("/data/PJB1_S1_I1_001.fastq.gz",),
                          ("/data/PJB1_S1_R1_001.fastq.gz",
                           "/data/PJB1_S1_R2_001.fastq.gz"),
                          ("/data/PJB2_S2_I1_001.fastq.gz",),
                          ("/data/PJB2_S2_R1_001.fastq.gz",
                           "/data/PJB2_S2_R2_001.fastq.gz")])

    def test_pair_fastqs_by_name_mixed_SE_and_PE(self):
        """
        pair_fastqs_by_name: mixture of single- and paired-end fastqs
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB1_S1_R2_001.fastq.gz",
                  "/data/PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(pair_fastqs_by_name(fastqs),
                         [("/data/PJB1_S1_R1_001.fastq.gz",
                           "/data/PJB1_S1_R2_001.fastq.gz"),
                          ("/data/PJB2_S2_R1_001.fastq.gz",)])

    def test_pair_fastqs_by_name_unpaired_R2(self):
        """
        pair_fastqs_by_name: handle unpaired R2 fastqs
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB1_S1_R2_001.fastq.gz",
                  "/data/PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(pair_fastqs_by_name(fastqs),
                         [("/data/PJB1_S1_R1_001.fastq.gz",
                           "/data/PJB1_S1_R2_001.fastq.gz"),
                          ("/data/PJB2_S2_R2_001.fastq.gz",)])

# group_fastqs_by_name
class TestGroupFastqsByNameFunction(unittest.TestCase):
    """
    Tests for the group_fastqs_by_name function
    """
    def test_group_fastqs_by_name_SE(self):
        """
        group_fastqs_by_name: single-end fastqs
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(group_fastqs_by_name(fastqs),
                         [["/data/PJB1_S1_R1_001.fastq.gz",],
                          ["/data/PJB2_S2_R1_001.fastq.gz",]])

    def test_group_fastqs_by_name_PE(self):
        """
        group_fastqs_by_name: paired-end fastqs
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB1_S1_R2_001.fastq.gz",
                  "/data/PJB2_S2_R1_001.fastq.gz",
                  "/data/PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(group_fastqs_by_name(fastqs),
                         [["/data/PJB1_S1_R1_001.fastq.gz",
                           "/data/PJB1_S1_R2_001.fastq.gz"],
                          ["/data/PJB2_S2_R1_001.fastq.gz",
                           "/data/PJB2_S2_R2_001.fastq.gz"]])

    def test_group_fastqs_by_name_PE_with_index_read(self):
        """
        group_fastqs_by_name: paired-end fastqs with index reads
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB1_S1_R2_001.fastq.gz",
                  "/data/PJB1_S1_I1_001.fastq.gz",
                  "/data/PJB2_S2_R1_001.fastq.gz",
                  "/data/PJB2_S2_R2_001.fastq.gz",
                  "/data/PJB2_S2_I1_001.fastq.gz")
        self.assertEqual(group_fastqs_by_name(fastqs),
                         [["/data/PJB1_S1_R1_001.fastq.gz",
                           "/data/PJB1_S1_R2_001.fastq.gz",
                           "/data/PJB1_S1_I1_001.fastq.gz",],
                          ["/data/PJB2_S2_R1_001.fastq.gz",
                           "/data/PJB2_S2_R2_001.fastq.gz",
                           "/data/PJB2_S2_I1_001.fastq.gz",]])

    def test_group_fastqs_by_name_PE_with_index_read_pair(self):
        """
        group_fastqs_by_name: paired-end fastqs with index read pair
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB1_S1_R2_001.fastq.gz",
                  "/data/PJB1_S1_I1_001.fastq.gz",
                  "/data/PJB1_S1_I2_001.fastq.gz",
                  "/data/PJB2_S2_R1_001.fastq.gz",
                  "/data/PJB2_S2_R2_001.fastq.gz",
                  "/data/PJB2_S2_I1_001.fastq.gz",
                  "/data/PJB2_S2_I2_001.fastq.gz")
        self.assertEqual(group_fastqs_by_name(fastqs),
                         [["/data/PJB1_S1_R1_001.fastq.gz",
                           "/data/PJB1_S1_R2_001.fastq.gz",
                           "/data/PJB1_S1_I1_001.fastq.gz",
                           "/data/PJB1_S1_I2_001.fastq.gz",],
                          ["/data/PJB2_S2_R1_001.fastq.gz",
                           "/data/PJB2_S2_R2_001.fastq.gz",
                           "/data/PJB2_S2_I1_001.fastq.gz",
                           "/data/PJB2_S2_I2_001.fastq.gz",]])

    def test_group_fastqs_by_name_with_R3_and_index_read(self):
        """
        group_fastqs_by_name: R1, R2 and R3 fastqs with index reads
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB1_S1_R2_001.fastq.gz",
                  "/data/PJB1_S1_R3_001.fastq.gz",
                  "/data/PJB1_S1_I1_001.fastq.gz",
                  "/data/PJB2_S2_R1_001.fastq.gz",
                  "/data/PJB2_S2_R2_001.fastq.gz",
                  "/data/PJB2_S2_R3_001.fastq.gz",
                  "/data/PJB2_S2_I1_001.fastq.gz")
        self.assertEqual(group_fastqs_by_name(fastqs),
                         [["/data/PJB1_S1_R1_001.fastq.gz",
                           "/data/PJB1_S1_R2_001.fastq.gz",
                           "/data/PJB1_S1_R3_001.fastq.gz",
                           "/data/PJB1_S1_I1_001.fastq.gz",],
                          ["/data/PJB2_S2_R1_001.fastq.gz",
                           "/data/PJB2_S2_R2_001.fastq.gz",
                           "/data/PJB2_S2_R3_001.fastq.gz",
                           "/data/PJB2_S2_I1_001.fastq.gz",]])

    def test_group_fastqs_by_name_mixed_SE_and_PE(self):
        """
        group_fastqs_by_name: mixture of single- and paired-end fastqs
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB1_S1_R2_001.fastq.gz",
                  "/data/PJB2_S2_R1_001.fastq.gz")
        self.assertEqual(group_fastqs_by_name(fastqs),
                         [["/data/PJB1_S1_R1_001.fastq.gz",
                           "/data/PJB1_S1_R2_001.fastq.gz"],
                          ["/data/PJB2_S2_R1_001.fastq.gz",]])

    def test_group_fastqs_by_name_unpaired_R2(self):
        """
        group_fastqs_by_name: handle unpaired R2 fastqs
        """
        fastqs = ("/data/PJB1_S1_R1_001.fastq.gz",
                  "/data/PJB1_S1_R2_001.fastq.gz",
                  "/data/PJB2_S2_R2_001.fastq.gz")
        self.assertEqual(group_fastqs_by_name(fastqs),
                         [["/data/PJB1_S1_R1_001.fastq.gz",
                           "/data/PJB1_S1_R2_001.fastq.gz"],
                          ["/data/PJB2_S2_R2_001.fastq.gz",]])

# get_read_number
class TestGetReadNumber(unittest.TestCase):
    """
    Tests for the get_read_number function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_get_read_number')
        # Test files
        self.fastq_r1 = os.path.join(self.wd,'test.r1.fq')
        with open(self.fastq_r1,'w') as fp:
            fp.write(fastq1_r1)
        self.fastq_r2 = os.path.join(self.wd,'test.r2.fq')
        with open(self.fastq_r2,'w') as fp:
            fp.write(fastq1_r2)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_get_read_number_r1(self):
        """get_read_number: check read number for R1 Fastq file
        """
        self.assertEqual(get_read_number(self.fastq_r1),1)
    def test_get_read_number_r2(self):
        """get_read_number: check read number for R2 Fastq file
        """
        self.assertEqual(get_read_number(self.fastq_r2),2)

# get_read_number
class TestGetReadNumber(unittest.TestCase):
    """
    Tests for the get_read_number function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_get_read_number')
        # Test files
        self.fastq_r1 = os.path.join(self.wd,'test.r1.fq')
        with open(self.fastq_r1,'w') as fp:
            fp.write(fastq1_r1)
        self.fastq_r2 = os.path.join(self.wd,'test.r2.fq')
        with open(self.fastq_r2,'w') as fp:
            fp.write(fastq1_r2)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_get_read_number_r1(self):
        """get_read_number: check read number for R1 Fastq file
        """
        self.assertEqual(get_read_number(self.fastq_r1),1)
    def test_get_read_number_r2(self):
        """get_read_number: check read number for R2 Fastq file
        """
        self.assertEqual(get_read_number(self.fastq_r2),2)

# get_read_count
class TestGetReadCount(unittest.TestCase):
    """
    Tests for the get_read_count function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.test_get_read_number')
        # Test files
        self.fastq1 = os.path.join(self.wd,'test1.r1.fq')
        with open(self.fastq1,'w') as fp:
            fp.write(fastq1_r1)
        self.fastq2 = os.path.join(self.wd,'test2.r1.fq')
        with open(self.fastq2,'w') as fp:
            fp.write(fastq2_r1)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_get_read_count_single_fastq(self):
        """get_read_count: check read count for a single Fastq
        """
        self.assertEqual(get_read_count((self.fastq1,)),3)
    def test_get_read_count_multiple_fastqs(self):
        """get_read_count: check read count for multiple Fastqs
        """
        self.assertEqual(get_read_count((self.fastq1,self.fastq2)),5)
    def test_get_read_count_no_fastqs(self):
        """get_read_count: check read count for no Fastqs
        """
        self.assertEqual(get_read_count(list()),0)

# remove_index_fastqs
class TestRemoveIndexFastqs(unittest.TestCase):
    """
    Tests for the remove_index_fastqs function
    """
    def test_remove_index_fastqs_no_index_reads(self):
        """remove_index_fastqs: list with no index Fastqs is unchanged
        """
        fastqs = ["PJB1_S1_R1_001.fastq.gz",
                  "PJB1_S1_R2_001.fastq.gz",
                  "PJB2_S2_R1_001.fastq.gz",
                  "PJB2_S2_R2_001.fastq.gz",]
        self.assertEqual(remove_index_fastqs(fastqs),fastqs)
    def test_remove_index_fastqs_i1_removed(self):
        """remove_index_fastqs: I1 index Fastqs are removed
        """
        fastqs_in = ["PJB1_S1_R1_001.fastq.gz",
                     "PJB1_S1_R2_001.fastq.gz",
                     "PJB1_S1_I1_001.fastq.gz",
                     "PJB2_S2_R1_001.fastq.gz",
                     "PJB2_S2_R2_001.fastq.gz",
                     "PJB2_S2_I1_001.fastq.gz",]
        fastqs_out = ["PJB1_S1_R1_001.fastq.gz",
                      "PJB1_S1_R2_001.fastq.gz",
                      "PJB2_S2_R1_001.fastq.gz",
                      "PJB2_S2_R2_001.fastq.gz",]
        self.assertEqual(remove_index_fastqs(fastqs_in),fastqs_out)
    def test_remove_index_fastqs_i1_i2_pairs_removed(self):
        """remove_index_fastqs: I1/I2 index Fastq pairs are removed
        """
        fastqs_in = ["PJB1_S1_R1_001.fastq.gz",
                     "PJB1_S1_R2_001.fastq.gz",
                     "PJB1_S1_I1_001.fastq.gz",
                     "PJB1_S1_I2_001.fastq.gz",
                     "PJB2_S2_R1_001.fastq.gz",
                     "PJB2_S2_R2_001.fastq.gz",
                     "PJB2_S2_I1_001.fastq.gz",
                     "PJB2_S2_I2_001.fastq.gz",]
        fastqs_out = ["PJB1_S1_R1_001.fastq.gz",
                      "PJB1_S1_R2_001.fastq.gz",
                      "PJB2_S2_R1_001.fastq.gz",
                      "PJB2_S2_R2_001.fastq.gz",]
        self.assertEqual(remove_index_fastqs(fastqs_in),fastqs_out)

