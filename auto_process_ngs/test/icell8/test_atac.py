#######################################################################
# Tests for icell8.atac.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
import gzip
from auto_process_ngs.icell8.atac import reverse_complement
from auto_process_ngs.icell8.atac import update_fastq_read_index
from auto_process_ngs.icell8.atac import split_fastq
from auto_process_ngs.icell8.atac import assign_reads
from auto_process_ngs.icell8.atac import concat_fastqs

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = False

fastq_data = {
    'R1': """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:1
TAAACATTCTGGGGGTTGGGGTGAGGTNTNNNNNNNNA
+
AA/AAEEEEEEEEEEEAEEEEAEAEEE#E########E
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 1:N:0:1
AGTCTGAGATGTCTGAATCTGATCTTCNAANNNNNNAT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######EE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 1:N:0:1
CAGCTGTTCTCATCATGATCTTTATAATTTNNNNNNC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######E
@NB500968:115:HWJNYBGX9:1:11101:9835:1054 1:N:0:1
TTTCTGTAGTGTGGCGTGTTGGTGTNGNCNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEE#A#E########E
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 1:N:0:1
AAATATGGCGAGGAAAACTGAAAAAGGNGNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEA#/########E
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 1:N:0:1
CTCCTTCTCTGATTGATCAGATAGCTCNTGNNNNNN
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######
@NB500968:115:HWJNYBGX9:1:11101:24409:1053 1:N:0:1
ACTCTCTTCTAATGGAGGACTGTGTNANANNNNNNNNT
+
AAAAAEEEEEEEEEEEEEEEEEEEE#E#E########E
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 1:N:0:1
GTGTGGTGGTGTGTACTCCTCCAATCCCAGNNNNTNC
+
AAAAAEEEEEEEEEAEEEEEEEEEAEEA/E####<#E
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 1:N:0:1
GACATACTAGGAGACCCAGACAACTACATACCNGCTAA
+
AAAAAEEEEEEE/EEAEEAEAEEAEEEEEAE/#/EEEE
@NB500968:115:HWJNYBGX9:1:11101:21842:1053 1:N:0:1
ACAAAAAATTGCTCCCCTATCAATTNTNANNNNNNNNT
+
AAAAAEEEEEEEEEEEEEEEEE/EE#E#E########E
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 1:N:0:1
CAGCTAAGAGCATCGAGGGGGCGCCGAGAGNNNNNNGG
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######EE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 1:N:0:1
AGATATAGCATTCCCACGAATAAATAATATNANNTNTT
+
AAAAA6EE/A/A///AEEAEEAEEEEEEEE#/##E#AE
""",
    'R2': """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 2:N:0:1
TTTCCACCCAGAAGGATGGGAGCAGATGTTAATAAC
+
AAAAAAEEEEEEAEEEEEAAEEEEE/EEEAEEAAEE
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 2:N:0:1
GAATTCTGTAGTCTGAGCCACCTGTCCAGTGAGCCCGG
+
AAAAAEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 2:N:0:1
ATATGGTGGAGGGCAGCCATGAAGTCATTCTAAATTTG
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@NB500968:115:HWJNYBGX9:1:11101:9835:1054 2:N:0:1
AACGTGCTCCTTCTATCCGGGCAACACCAACACGCCA
+
AAAAAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEA
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 2:N:0:1
TCCTACAGTGGACTTTTCTAAATTTTCCACCTTTTTCA
+
AAAAAEEEEEAEEEEAEEEEEEEEEE//<EAAAAEEAE
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 2:N:0:1
CTGAAGGGACTTGGGCCATGATAGAAAGTGAGAATTTA
+
AAAAAEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEE
@NB500968:115:HWJNYBGX9:1:11101:24409:1053 2:N:0:1
ACCCTATTGTGTAATGTGCATGACATATGGCATACTAA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 2:N:0:1
ATATAGACCAGGGGGCATCTGCTGGGATTGGAGGAGTA
+
AAAAAE6EEEE/EEEEEEEEEEEEEAAEEEEEEE/<EE
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 2:N:0:1
GTGTTTAGTGGATTAGCTGGTATGTAGTTGTCTGGGT
+
AAAAAEAEEEEEEEEE/EE6EEEEEEEEEEEEEEEEE
@NB500968:115:HWJNYBGX9:1:11101:21842:1053 2:N:0:1
ATTGCTAATATTCATCCTATGTGGGCAATTGATGAATA
+
AAAAAEEEE<EEEEEEEEEEEEEEEEEEEEAEEEEEE/
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 2:N:0:1
TCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEAEEEEAEE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 2:N:0:1
CTATTGATGATGCTAGTAGAAGGAGAAATGATGGTGGT
+
6AAAAE/AE/A<AAE//EEAEE/EEEEEEEEEE/EEE/
""",
    'I1': """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:1
AGTAGATT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 1:N:0:1
ATCACTCG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 1:N:0:1
AGTAGATT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:9835:1054 1:N:0:1
TAAGGCGA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 1:N:0:1
GCTACGCT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 1:N:0:1
CATCCTGT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:24409:1053 1:N:0:1
CGAGGCTG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 1:N:0:1
GATCCAAA
+
A/A/////
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 1:N:0:1
TTCCATAT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:21842:1053 1:N:0:1
AAGAGGCA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 1:N:0:1
CAATCTTA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 1:N:0:1
AAGAGGCA
+
AA6A66AA
""",
    'I2': """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 2:N:0:1
GCCACGTC
+
AAA/AEEE
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 2:N:0:1
CCCGCAGG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 2:N:0:1
GCCACGTC
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:9835:1054 2:N:0:1
CTTGGTTA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 2:N:0:1
CGTTGCTG
+
6AAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 2:N:0:1
TGTAGATT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:24409:1053 2:N:0:1
CCCTATCG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 2:N:0:1
TGCACGAA
+
AAAAAEE/
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 2:N:0:1
CCCTATCG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:21842:1053 2:N:0:1
TACTCCTT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 2:N:0:1
TTTCATCA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 2:N:0:1
TACTCCTT
+
AA/AAEEE
""",
    'B000': """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:1
TAAACATTCTGGGGGTTGGGGTGAGGTNTNNNNNNNNA
+
AA/AAEEEEEEEEEEEAEEEEAEAEEE#E########E
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 1:N:0:1
AGTCTGAGATGTCTGAATCTGATCTTCNAANNNNNNAT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######EE
""",
    'B001': """@NB500968:115:HWJNYBGX9:1:11101:23324:1057 1:N:0:1
CAGCTGTTCTCATCATGATCTTTATAATTTNNNNNNC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######E
@NB500968:115:HWJNYBGX9:1:11101:9835:1054 1:N:0:1
TTTCTGTAGTGTGGCGTGTTGGTGTNGNCNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEE#A#E########E
""",
    'B002': """@NB500968:115:HWJNYBGX9:1:11101:4921:1055 1:N:0:1
AAATATGGCGAGGAAAACTGAAAAAGGNGNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEA#/########E
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 1:N:0:1
CTCCTTCTCTGATTGATCAGATAGCTCNTGNNNNNN
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######
""",
}

well_list_data = """Row	Col	Candidate	For dispense	Sample	Barcode	State	Cells1	Cells2	Signal1	Signal2	Size1	Size2	Integ Signal1	Integ Signal2	Circularity1	Circularity2	Confidence	Confidence1	Confidence2	Dispense tip	Drop index	Global drop index	Source well	Sequencing count	Image1	Image2
0	57	True	True	1	TAACCAAG+TAAGGCGA	Good	1	0	105		33		3465		1		1	1	1	1	10	10	A1	Pos9_1-Hoechst_A10.tif	Pos9_3-Texas Red_A10.tif
2	27	True	True	1	CAGCAACG+GCTACGCT	Good	1	0	96		23		2208		1		1	1	1	1	28	34	A1	Pos4_1-Hoechst_A05.tif	Pos4_3-Texas Red_A05.tif
11	34	True	True	2	CGATAGGG+CGAGGCTG	Good	1	0	101		23		2323		1		0.97	0.97	1	3	35	41	B1	Pos17_1-Hoechst_B06.tif	Pos17_3-Texas Red_B06.tif
20	1	True	True	3	AAGGAGTA+AAGAGGCA	Good	1	0	97		18		1746		1		1	1	1	5	26	32	C1	Pos36_1-Hoechst_D01.tif	Pos36_3-Texas Red_D01.tif
24	2	True	True	3	TGATGAAA+CAATCTTA	Good	1	0	83		13		1079		1		1	1	1	5	55	67	C1	Pos48_1-Hoechst_E01.tif	Pos48_3-Texas Red_E01.tif
41	51	True	True	1	AATCTACA+CATCCTGT	Good	1	0	98		24		2352		1		0.98	0.98	1	1	163	206	A1	Pos80_1-Hoechst_G09.tif	Pos80_3-Texas Red_G09.tif
46	34	True	True	2	CGATAGGG+TTCCATAT	Good	1	0	94		18		1692		1		0.97	0.97	1	3	129	159	B1	Pos89_1-Hoechst_H06.tif	Pos89_3-Texas Red_H06.tif
52	27	True	True	2	TTCGTGCA+GATCCAAA	Good	1	0	93		18		1674		1		1	1	1	3	192	236	B1	Pos100_1-Hoechst_I05.tif	Pos100_3-Texas Red_I05.tif
65	66	True	True	10	GACGTGGC+AGTAGATT	Good	1	0	105		41		4305		1		0.8611708	0.8611708	1	8	136	169	D2		Pos131_1-Hoechst_K12.tif	Pos131_3-Texas Red_K12.tif
66	22	True	True	10	CCTGCGGG+ATCACTCG	Good	1	0	195		21		4095		1		0.97	0.97	1	8	153	192	D2	Pos135_1-Hoechst_L04.tif	Pos135_3-Texas Red_L04.tif
"""

# reverse_complement
class TestReverseComplementFunction(unittest.TestCase):
    """
    Tests for the reverse_complement function
    """
    def test_reverse_complement(self):
        """
        reverse_complement: check sequences are reverse complemented
        """
        self.assertEqual(reverse_complement("ATGC"),"GCAT")

    def test_reverse_complement_handle_ns(self):
        """
        reverse_complement: handle Ns in sequence
        """
        self.assertEqual(reverse_complement("ATGCN"),"NGCAT")

# update_fastq_read_index
class TestUpdateFastqReadIndexFunction(unittest.TestCase):
    """
    Tests for the update_fastq_read_index function
    """
    def test_update_fastq_read_index(self):
        """
        update_fastq_read_index: check index sequence is updated
        """
        read_in = ["@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:1",
                   "TAAACATTCTGGGGGTTGGGGTGAGGTNTNNNNNNNNA",
                   "+",
                   "AA/AAEEEEEEEEEEEAEEEEAEAEEE#E########E'"]
        read_out = update_fastq_read_index(read_in,"TAACCAAG+TAAGGCGA")
        self.assertEqual(read_out,
                         ["@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:TAACCAAG+TAAGGCGA",
                          "TAAACATTCTGGGGGTTGGGGTGAGGTNTNNNNNNNNA",
                          "+",
                          "AA/AAEEEEEEEEEEEAEEEEAEAEEE#E########E'"])

# split_fastq
class TestSplitFastqFunction(unittest.TestCase):
    """
    Tests for the split_fastq function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.SplitFastqs')
        # Test file
        self.fastq = os.path.join(self.wd,'test.fq')
        with open(self.fastq,'w') as fp:
            fp.write(fastq_data['R1'])
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            if REMOVE_TEST_OUTPUTS:
                shutil.rmtree(self.wd)
    def test_split_fastq(self):
        """
        split_fastq: splits Fastq file into batches
        """
        # Split the test Fastq into batches of 5 reads
        fastqs = split_fastq((self.fastq,5,self.wd))
        # Check the returned list
        self.assertEqual(fastqs,[os.path.join(self.wd,"test_B000.fastq"),
                                 os.path.join(self.wd,"test_B001.fastq"),
                                 os.path.join(self.wd,"test_B002.fastq")])
        # Check the contents
        self.assertEqual(open(os.path.join(self.wd,"test_B000.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:1
TAAACATTCTGGGGGTTGGGGTGAGGTNTNNNNNNNNA
+
AA/AAEEEEEEEEEEEAEEEEAEAEEE#E########E
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 1:N:0:1
AGTCTGAGATGTCTGAATCTGATCTTCNAANNNNNNAT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######EE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 1:N:0:1
CAGCTGTTCTCATCATGATCTTTATAATTTNNNNNNC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######E
@NB500968:115:HWJNYBGX9:1:11101:9835:1054 1:N:0:1
TTTCTGTAGTGTGGCGTGTTGGTGTNGNCNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEE#A#E########E
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 1:N:0:1
AAATATGGCGAGGAAAACTGAAAAAGGNGNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEA#/########E
""")
        self.assertEqual(open(os.path.join(self.wd,"test_B001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:12850:1056 1:N:0:1
CTCCTTCTCTGATTGATCAGATAGCTCNTGNNNNNN
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######
@NB500968:115:HWJNYBGX9:1:11101:24409:1053 1:N:0:1
ACTCTCTTCTAATGGAGGACTGTGTNANANNNNNNNNT
+
AAAAAEEEEEEEEEEEEEEEEEEEE#E#E########E
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 1:N:0:1
GTGTGGTGGTGTGTACTCCTCCAATCCCAGNNNNTNC
+
AAAAAEEEEEEEEEAEEEEEEEEEAEEA/E####<#E
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 1:N:0:1
GACATACTAGGAGACCCAGACAACTACATACCNGCTAA
+
AAAAAEEEEEEE/EEAEEAEAEEAEEEEEAE/#/EEEE
@NB500968:115:HWJNYBGX9:1:11101:21842:1053 1:N:0:1
ACAAAAAATTGCTCCCCTATCAATTNTNANNNNNNNNT
+
AAAAAEEEEEEEEEEEEEEEEE/EE#E#E########E
""")
        self.assertEqual(open(os.path.join(self.wd,"test_B002.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:12602:1058 1:N:0:1
CAGCTAAGAGCATCGAGGGGGCGCCGAGAGNNNNNNGG
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######EE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 1:N:0:1
AGATATAGCATTCCCACGAATAAATAATATNANNTNTT
+
AAAAA6EE/A/A///AEEAEEAEEEEEEEE#/##E#AE
""")

# assign_reads
class TestAssignReadsFunction(unittest.TestCase):
    """
    Tests for the assign_reads function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.AssignReads')
        # Test Fastqs
        self.fastqs = []
        for read in ('R1','R2','I1','I2'):
            fastq = os.path.join(self.wd,'test_%s_B002.fq' % read)
            with open(fastq,'w') as fp:
                fp.write(fastq_data[read])
            self.fastqs.append(fastq)
        # Well list
        self.well_list = os.path.join(self.wd,'well_list.txt')
        with open(self.well_list,'w') as fp:
            fp.write(well_list_data)
        # Expected barcodes
        self.expected_barcode_counts = {
            'TTCGTGCA+GATCCAAA': 1,
            'AAGGAGTA+AAGAGGCA': 2,
            'CCTGCGGG+ATCACTCG': 1,
            'CGATAGGG+TTCCATAT': 1,
            'CAGCAACG+GCTACGCT': 1,
            'GACGTGGC+AGTAGATT': 2,
            'CGATAGGG+CGAGGCTG': 1,
            'TGATGAAA+CAATCTTA': 1,
            'TAACCAAG+TAAGGCGA': 1,
            'AATCTACA+CATCCTGT': 1,
            'unassigned': 0,
        }
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            if REMOVE_TEST_OUTPUTS:
                shutil.rmtree(self.wd)
    def test_assign_reads_to_samples(self):
        """
        assign_reads: assigns reads to samples
        """
        # NB for these data I1 and I2 have to be flipped and
        # I1 has to be reverse complemented
        result = assign_reads((self.fastqs[0],
                               self.fastqs[1],
                               self.fastqs[2],
                               self.fastqs[3],
                               self.well_list,
                               'samples',
                               True,
                               'i1',
                               False,
                               self.wd,
                               "unassigned",))
        # Check the returned data
        batch_id,barcode_counts,unassigned_barcodes_file = result
        self.assertEqual(batch_id,"B002")
        for barcode in self.expected_barcode_counts:
            self.assertEqual(barcode_counts[barcode],
                             self.expected_barcode_counts[barcode])
        self.assertEqual(unassigned_barcodes_file,
                         os.path.join(self.wd,
                                      "B002",
                                      "unassigned_barcodes.txt"))
        # Check the output files
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "1_S1_R1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:9835:1054 1:N:0:1
TTTCTGTAGTGTGGCGTGTTGGTGTNGNCNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEE#A#E########E
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 1:N:0:1
AAATATGGCGAGGAAAACTGAAAAAGGNGNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEA#/########E
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 1:N:0:1
CTCCTTCTCTGATTGATCAGATAGCTCNTGNNNNNN
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "1_S1_R2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:9835:1054 2:N:0:1
AACGTGCTCCTTCTATCCGGGCAACACCAACACGCCA
+
AAAAAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEA
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 2:N:0:1
TCCTACAGTGGACTTTTCTAAATTTTCCACCTTTTTCA
+
AAAAAEEEEEAEEEEAEEEEEEEEEE//<EAAAAEEAE
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 2:N:0:1
CTGAAGGGACTTGGGCCATGATAGAAAGTGAGAATTTA
+
AAAAAEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "1_S1_I1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:9835:1054 1:N:0:1
TAAGGCGA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 1:N:0:1
GCTACGCT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 1:N:0:1
CATCCTGT
+
AAAAAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "1_S1_I2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:9835:1054 2:N:0:1
CTTGGTTA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 2:N:0:1
CGTTGCTG
+
6AAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 2:N:0:1
TGTAGATT
+
AAAAAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "2_S2_R1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:24409:1053 1:N:0:1
ACTCTCTTCTAATGGAGGACTGTGTNANANNNNNNNNT
+
AAAAAEEEEEEEEEEEEEEEEEEEE#E#E########E
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 1:N:0:1
GTGTGGTGGTGTGTACTCCTCCAATCCCAGNNNNTNC
+
AAAAAEEEEEEEEEAEEEEEEEEEAEEA/E####<#E
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 1:N:0:1
GACATACTAGGAGACCCAGACAACTACATACCNGCTAA
+
AAAAAEEEEEEE/EEAEEAEAEEAEEEEEAE/#/EEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "2_S2_R2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:24409:1053 2:N:0:1
ACCCTATTGTGTAATGTGCATGACATATGGCATACTAA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 2:N:0:1
ATATAGACCAGGGGGCATCTGCTGGGATTGGAGGAGTA
+
AAAAAE6EEEE/EEEEEEEEEEEEEAAEEEEEEE/<EE
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 2:N:0:1
GTGTTTAGTGGATTAGCTGGTATGTAGTTGTCTGGGT
+
AAAAAEAEEEEEEEEE/EE6EEEEEEEEEEEEEEEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "2_S2_I1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:24409:1053 1:N:0:1
CGAGGCTG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 1:N:0:1
GATCCAAA
+
A/A/////
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 1:N:0:1
TTCCATAT
+
AAAAAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "2_S2_I2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:24409:1053 2:N:0:1
CCCTATCG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 2:N:0:1
TGCACGAA
+
AAAAAEE/
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 2:N:0:1
CCCTATCG
+
AAAAAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "3_S3_R1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:21842:1053 1:N:0:1
ACAAAAAATTGCTCCCCTATCAATTNTNANNNNNNNNT
+
AAAAAEEEEEEEEEEEEEEEEE/EE#E#E########E
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 1:N:0:1
CAGCTAAGAGCATCGAGGGGGCGCCGAGAGNNNNNNGG
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######EE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 1:N:0:1
AGATATAGCATTCCCACGAATAAATAATATNANNTNTT
+
AAAAA6EE/A/A///AEEAEEAEEEEEEEE#/##E#AE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "3_S3_R2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:21842:1053 2:N:0:1
ATTGCTAATATTCATCCTATGTGGGCAATTGATGAATA
+
AAAAAEEEE<EEEEEEEEEEEEEEEEEEEEAEEEEEE/
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 2:N:0:1
TCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEAEEEEAEE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 2:N:0:1
CTATTGATGATGCTAGTAGAAGGAGAAATGATGGTGGT
+
6AAAAE/AE/A<AAE//EEAEE/EEEEEEEEEE/EEE/
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "3_S3_I1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:21842:1053 1:N:0:1
AAGAGGCA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 1:N:0:1
CAATCTTA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 1:N:0:1
AAGAGGCA
+
AA6A66AA
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "3_S3_I2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:21842:1053 2:N:0:1
TACTCCTT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 2:N:0:1
TTTCATCA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 2:N:0:1
TACTCCTT
+
AA/AAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "10_S4_R1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:1
TAAACATTCTGGGGGTTGGGGTGAGGTNTNNNNNNNNA
+
AA/AAEEEEEEEEEEEAEEEEAEAEEE#E########E
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 1:N:0:1
AGTCTGAGATGTCTGAATCTGATCTTCNAANNNNNNAT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######EE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 1:N:0:1
CAGCTGTTCTCATCATGATCTTTATAATTTNNNNNNC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######E
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "10_S4_R2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 2:N:0:1
TTTCCACCCAGAAGGATGGGAGCAGATGTTAATAAC
+
AAAAAAEEEEEEAEEEEEAAEEEEE/EEEAEEAAEE
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 2:N:0:1
GAATTCTGTAGTCTGAGCCACCTGTCCAGTGAGCCCGG
+
AAAAAEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 2:N:0:1
ATATGGTGGAGGGCAGCCATGAAGTCATTCTAAATTTG
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "10_S4_I1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:1
AGTAGATT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 1:N:0:1
ATCACTCG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 1:N:0:1
AGTAGATT
+
AAAAAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "10_S4_I2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 2:N:0:1
GCCACGTC
+
AAA/AEEE
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 2:N:0:1
CCCGCAGG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 2:N:0:1
GCCACGTC
+
AAAAAEEE
""")
    def test_assign_reads_to_samples_update_index_sequences(self):
        """
        assign_reads: assigns reads to samples and update index sequences
        """
        # NB for these data I1 and I2 have to be flipped and
        # I1 has to be reverse complemented
        result = assign_reads((self.fastqs[0],
                               self.fastqs[1],
                               self.fastqs[2],
                               self.fastqs[3],
                               self.well_list,
                               'samples',
                               True,
                               'i1',
                               True,
                               self.wd,
                               "unassigned",))
        # Check the returned data
        batch_id,barcode_counts,unassigned_barcodes_file = result
        self.assertEqual(batch_id,"B002")
        for barcode in self.expected_barcode_counts:
            self.assertEqual(barcode_counts[barcode],
                             self.expected_barcode_counts[barcode])
        self.assertEqual(unassigned_barcodes_file,
                         os.path.join(self.wd,
                                      "B002",
                                      "unassigned_barcodes.txt"))
        # Check the output files
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "1_S1_R1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:9835:1054 1:N:0:TAACCAAG+TAAGGCGA
TTTCTGTAGTGTGGCGTGTTGGTGTNGNCNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEE#A#E########E
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 1:N:0:CAGCAACG+GCTACGCT
AAATATGGCGAGGAAAACTGAAAAAGGNGNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEA#/########E
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 1:N:0:AATCTACA+CATCCTGT
CTCCTTCTCTGATTGATCAGATAGCTCNTGNNNNNN
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "1_S1_R2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:9835:1054 2:N:0:TAACCAAG+TAAGGCGA
AACGTGCTCCTTCTATCCGGGCAACACCAACACGCCA
+
AAAAAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEA
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 2:N:0:CAGCAACG+GCTACGCT
TCCTACAGTGGACTTTTCTAAATTTTCCACCTTTTTCA
+
AAAAAEEEEEAEEEEAEEEEEEEEEE//<EAAAAEEAE
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 2:N:0:AATCTACA+CATCCTGT
CTGAAGGGACTTGGGCCATGATAGAAAGTGAGAATTTA
+
AAAAAEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "1_S1_I1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:9835:1054 1:N:0:TAACCAAG+TAAGGCGA
TAAGGCGA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 1:N:0:CAGCAACG+GCTACGCT
GCTACGCT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 1:N:0:AATCTACA+CATCCTGT
CATCCTGT
+
AAAAAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "1_S1_I2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:9835:1054 2:N:0:TAACCAAG+TAAGGCGA
CTTGGTTA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 2:N:0:CAGCAACG+GCTACGCT
CGTTGCTG
+
6AAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 2:N:0:AATCTACA+CATCCTGT
TGTAGATT
+
AAAAAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "2_S2_R1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:24409:1053 1:N:0:CGATAGGG+CGAGGCTG
ACTCTCTTCTAATGGAGGACTGTGTNANANNNNNNNNT
+
AAAAAEEEEEEEEEEEEEEEEEEEE#E#E########E
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 1:N:0:TTCGTGCA+GATCCAAA
GTGTGGTGGTGTGTACTCCTCCAATCCCAGNNNNTNC
+
AAAAAEEEEEEEEEAEEEEEEEEEAEEA/E####<#E
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 1:N:0:CGATAGGG+TTCCATAT
GACATACTAGGAGACCCAGACAACTACATACCNGCTAA
+
AAAAAEEEEEEE/EEAEEAEAEEAEEEEEAE/#/EEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "2_S2_R2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:24409:1053 2:N:0:CGATAGGG+CGAGGCTG
ACCCTATTGTGTAATGTGCATGACATATGGCATACTAA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 2:N:0:TTCGTGCA+GATCCAAA
ATATAGACCAGGGGGCATCTGCTGGGATTGGAGGAGTA
+
AAAAAE6EEEE/EEEEEEEEEEEEEAAEEEEEEE/<EE
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 2:N:0:CGATAGGG+TTCCATAT
GTGTTTAGTGGATTAGCTGGTATGTAGTTGTCTGGGT
+
AAAAAEAEEEEEEEEE/EE6EEEEEEEEEEEEEEEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "2_S2_I1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:24409:1053 1:N:0:CGATAGGG+CGAGGCTG
CGAGGCTG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 1:N:0:TTCGTGCA+GATCCAAA
GATCCAAA
+
A/A/////
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 1:N:0:CGATAGGG+TTCCATAT
TTCCATAT
+
AAAAAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "2_S2_I2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:24409:1053 2:N:0:CGATAGGG+CGAGGCTG
CCCTATCG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:3460:1059 2:N:0:TTCGTGCA+GATCCAAA
TGCACGAA
+
AAAAAEE/
@NB500968:115:HWJNYBGX9:1:11101:22996:1060 2:N:0:CGATAGGG+TTCCATAT
CCCTATCG
+
AAAAAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "3_S3_R1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:21842:1053 1:N:0:AAGGAGTA+AAGAGGCA
ACAAAAAATTGCTCCCCTATCAATTNTNANNNNNNNNT
+
AAAAAEEEEEEEEEEEEEEEEE/EE#E#E########E
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 1:N:0:TGATGAAA+CAATCTTA
CAGCTAAGAGCATCGAGGGGGCGCCGAGAGNNNNNNGG
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######EE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 1:N:0:AAGGAGTA+AAGAGGCA
AGATATAGCATTCCCACGAATAAATAATATNANNTNTT
+
AAAAA6EE/A/A///AEEAEEAEEEEEEEE#/##E#AE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "3_S3_R2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:21842:1053 2:N:0:AAGGAGTA+AAGAGGCA
ATTGCTAATATTCATCCTATGTGGGCAATTGATGAATA
+
AAAAAEEEE<EEEEEEEEEEEEEEEEEEEEAEEEEEE/
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 2:N:0:TGATGAAA+CAATCTTA
TCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEAEEEEAEE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 2:N:0:AAGGAGTA+AAGAGGCA
CTATTGATGATGCTAGTAGAAGGAGAAATGATGGTGGT
+
6AAAAE/AE/A<AAE//EEAEE/EEEEEEEEEE/EEE/
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "3_S3_I1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:21842:1053 1:N:0:AAGGAGTA+AAGAGGCA
AAGAGGCA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 1:N:0:TGATGAAA+CAATCTTA
CAATCTTA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 1:N:0:AAGGAGTA+AAGAGGCA
AAGAGGCA
+
AA6A66AA
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "3_S3_I2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:21842:1053 2:N:0:AAGGAGTA+AAGAGGCA
TACTCCTT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:12602:1058 2:N:0:TGATGAAA+CAATCTTA
TTTCATCA
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:1239:1059 2:N:0:AAGGAGTA+AAGAGGCA
TACTCCTT
+
AA/AAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "10_S4_R1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:GACGTGGC+AGTAGATT
TAAACATTCTGGGGGTTGGGGTGAGGTNTNNNNNNNNA
+
AA/AAEEEEEEEEEEEAEEEEAEAEEE#E########E
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 1:N:0:CCTGCGGG+ATCACTCG
AGTCTGAGATGTCTGAATCTGATCTTCNAANNNNNNAT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######EE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 1:N:0:GACGTGGC+AGTAGATT
CAGCTGTTCTCATCATGATCTTTATAATTTNNNNNNC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######E
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "10_S4_R2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 2:N:0:GACGTGGC+AGTAGATT
TTTCCACCCAGAAGGATGGGAGCAGATGTTAATAAC
+
AAAAAAEEEEEEAEEEEEAAEEEEE/EEEAEEAAEE
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 2:N:0:CCTGCGGG+ATCACTCG
GAATTCTGTAGTCTGAGCCACCTGTCCAGTGAGCCCGG
+
AAAAAEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 2:N:0:GACGTGGC+AGTAGATT
ATATGGTGGAGGGCAGCCATGAAGTCATTCTAAATTTG
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "10_S4_I1_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:GACGTGGC+AGTAGATT
AGTAGATT
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 1:N:0:CCTGCGGG+ATCACTCG
ATCACTCG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 1:N:0:GACGTGGC+AGTAGATT
AGTAGATT
+
AAAAAEEE
""")
        self.assertEqual(open(os.path.join(self.wd,
                                           "B002",
                                           "10_S4_I2_001.fastq")).read(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 2:N:0:GACGTGGC+AGTAGATT
GCCACGTC
+
AAA/AEEE
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 2:N:0:CCTGCGGG+ATCACTCG
CCCGCAGG
+
AAAAAEEE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 2:N:0:GACGTGGC+AGTAGATT
GCCACGTC
+
AAAAAEEE
""")
    def test_assign_reads_to_barcodes(self):
        """
        assign_reads: assigns reads to barcodes
        """
        # NB for these data I1 and I2 have to be flipped and
        # I1 has to be reverse complemented
        result = assign_reads((self.fastqs[0],
                               self.fastqs[1],
                               self.fastqs[2],
                               self.fastqs[3],
                               self.well_list,
                               'barcodes',
                               True,
                               'i1',
                               False,
                               self.wd,
                               "unassigned",))
        # Check the returned data
        batch_id,barcode_counts,unassigned_barcodes_file = result
        self.assertEqual(batch_id,"B002")
        for barcode in self.expected_barcode_counts:
            self.assertEqual(barcode_counts[barcode],
                             self.expected_barcode_counts[barcode])
        self.assertEqual(unassigned_barcodes_file,
                         os.path.join(self.wd,
                                      "B002",
                                      "unassigned_barcodes.txt"))
        # Check the output files
        for f in ("10_S4_CCTGCGGG+ATCACTCG_I1_001.fastq",
                  "10_S4_CCTGCGGG+ATCACTCG_I2_001.fastq",
                  "10_S4_CCTGCGGG+ATCACTCG_R1_001.fastq",
                  "10_S4_CCTGCGGG+ATCACTCG_R2_001.fastq",
                  "10_S4_GACGTGGC+AGTAGATT_I1_001.fastq",
                  "10_S4_GACGTGGC+AGTAGATT_I2_001.fastq",
                  "10_S4_GACGTGGC+AGTAGATT_R1_001.fastq",
                  "10_S4_GACGTGGC+AGTAGATT_R2_001.fastq",
                  "1_S1_AATCTACA+CATCCTGT_I1_001.fastq",
                  "1_S1_AATCTACA+CATCCTGT_I2_001.fastq",
                  "1_S1_AATCTACA+CATCCTGT_R1_001.fastq",
                  "1_S1_AATCTACA+CATCCTGT_R2_001.fastq",
                  "1_S1_CAGCAACG+GCTACGCT_I1_001.fastq",
                  "1_S1_CAGCAACG+GCTACGCT_I2_001.fastq",
                  "1_S1_CAGCAACG+GCTACGCT_R1_001.fastq",
                  "1_S1_CAGCAACG+GCTACGCT_R2_001.fastq",
                  "1_S1_TAACCAAG+TAAGGCGA_I1_001.fastq",
                  "1_S1_TAACCAAG+TAAGGCGA_I2_001.fastq",
                  "1_S1_TAACCAAG+TAAGGCGA_R1_001.fastq",
                  "1_S1_TAACCAAG+TAAGGCGA_R2_001.fastq",
                  "2_S2_CGATAGGG+CGAGGCTG_I1_001.fastq",
                  "2_S2_CGATAGGG+CGAGGCTG_I2_001.fastq",
                  "2_S2_CGATAGGG+CGAGGCTG_R1_001.fastq",
                  "2_S2_CGATAGGG+CGAGGCTG_R2_001.fastq",
                  "2_S2_CGATAGGG+TTCCATAT_I1_001.fastq",
                  "2_S2_CGATAGGG+TTCCATAT_I2_001.fastq",
                  "2_S2_CGATAGGG+TTCCATAT_R1_001.fastq",
                  "2_S2_CGATAGGG+TTCCATAT_R2_001.fastq",
                  "2_S2_TTCGTGCA+GATCCAAA_I1_001.fastq",
                  "2_S2_TTCGTGCA+GATCCAAA_I2_001.fastq",
                  "2_S2_TTCGTGCA+GATCCAAA_R1_001.fastq",
                  "2_S2_TTCGTGCA+GATCCAAA_R2_001.fastq",
                  "3_S3_AAGGAGTA+AAGAGGCA_I1_001.fastq",
                  "3_S3_AAGGAGTA+AAGAGGCA_I2_001.fastq",
                  "3_S3_AAGGAGTA+AAGAGGCA_R1_001.fastq",
                  "3_S3_AAGGAGTA+AAGAGGCA_R2_001.fastq",
                  "3_S3_TGATGAAA+CAATCTTA_I1_001.fastq",
                  "3_S3_TGATGAAA+CAATCTTA_I2_001.fastq",
                  "3_S3_TGATGAAA+CAATCTTA_R1_001.fastq",
                  "3_S3_TGATGAAA+CAATCTTA_R2_001.fastq",):
            self.assertTrue(os.path.exists(os.path.join(self.wd,
                                                        "B002",
                                                        f)))

# concat_fastqs
class TestConcatFastqsFunction(unittest.TestCase):
    """
    Tests for the concat_fastqs function
    """
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.ConcatFastqs')
        # Final directory
        self.final_dir = os.path.join(self.wd,"final")
        os.mkdir(self.final_dir)
    def tearDown(self):
        # Remove temporary working dir
        if os.path.isdir(self.wd):
            if REMOVE_TEST_OUTPUTS:
                shutil.rmtree(self.wd)
    def test_concat_fastqs_for_sample(self):
        """
        concat_fastqs: concatenate Fastq files for sample
        """
        # Make test Fastqs
        self.fastqs = []
        for batch in ('B000','B001','B002',):
            os.mkdir(os.path.join(self.wd,batch))
            fastq = os.path.join(self.wd,batch,'PJB1_S1_R1_001.fastq')
            with open(fastq,'w') as fp:
                fp.write(fastq_data[batch])
            self.fastqs.append(fastq)
        # Concatenate
        fastq = concat_fastqs(('PJB1',
                               1,
                               None,
                               None,
                               'R1',
                               ('B000','B001','B002',),
                               self.wd,
                               self.final_dir))
        # Check output
        self.assertEqual(fastq,
                         os.path.join(self.final_dir,
                                      "PJB1_S1_R1_001.fastq.gz"))
        self.assertEqual(gzip.open(fastq).read().decode(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:1
TAAACATTCTGGGGGTTGGGGTGAGGTNTNNNNNNNNA
+
AA/AAEEEEEEEEEEEAEEEEAEAEEE#E########E
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 1:N:0:1
AGTCTGAGATGTCTGAATCTGATCTTCNAANNNNNNAT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######EE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 1:N:0:1
CAGCTGTTCTCATCATGATCTTTATAATTTNNNNNNC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######E
@NB500968:115:HWJNYBGX9:1:11101:9835:1054 1:N:0:1
TTTCTGTAGTGTGGCGTGTTGGTGTNGNCNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEE#A#E########E
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 1:N:0:1
AAATATGGCGAGGAAAACTGAAAAAGGNGNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEA#/########E
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 1:N:0:1
CTCCTTCTCTGATTGATCAGATAGCTCNTGNNNNNN
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######
""")
        # Check original files were removed
        for batch in ('B000','B001','B002',):
            self.assertFalse(os.path.exists(
                os.path.join(self.wd,
                             batch,
                             "PJB1_S1_R1_001.fastq")))

    def test_concat_fastqs_for_barcode(self):
        """
        concat_fastqs: concatenate Fastq files for barcode
        """
        # Make test Fastqs
        self.fastqs = []
        for batch in ('B000','B001','B002',):
            os.mkdir(os.path.join(self.wd,batch))
            fastq = os.path.join(self.wd,batch,'PJB1_S1_TTCGTGCA+GATCCAAA_R1_001.fastq')
            with open(fastq,'w') as fp:
                fp.write(fastq_data[batch])
            self.fastqs.append(fastq)
        # Concatenate
        fastq = concat_fastqs(('PJB1',
                               1,
                               'TTCGTGCA+GATCCAAA',
                               None,
                               'R1',
                               ('B000','B001','B002',),
                               self.wd,
                               self.final_dir))
        # Check output
        self.assertEqual(fastq,
                         os.path.join(self.final_dir,
                                      "PJB1_S1_TTCGTGCA+GATCCAAA_R1_001.fastq.gz"))
        self.assertEqual(gzip.open(fastq).read().decode(),
                         """@NB500968:115:HWJNYBGX9:1:11101:4820:1056 1:N:0:1
TAAACATTCTGGGGGTTGGGGTGAGGTNTNNNNNNNNA
+
AA/AAEEEEEEEEEEEAEEEEAEAEEE#E########E
@NB500968:115:HWJNYBGX9:1:11101:11115:1057 1:N:0:1
AGTCTGAGATGTCTGAATCTGATCTTCNAANNNNNNAT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######EE
@NB500968:115:HWJNYBGX9:1:11101:23324:1057 1:N:0:1
CAGCTGTTCTCATCATGATCTTTATAATTTNNNNNNC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE######E
@NB500968:115:HWJNYBGX9:1:11101:9835:1054 1:N:0:1
TTTCTGTAGTGTGGCGTGTTGGTGTNGNCNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEE#A#E########E
@NB500968:115:HWJNYBGX9:1:11101:4921:1055 1:N:0:1
AAATATGGCGAGGAAAACTGAAAAAGGNGNNNNNNNNA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEA#/########E
@NB500968:115:HWJNYBGX9:1:11101:12850:1056 1:N:0:1
CTCCTTCTCTGATTGATCAGATAGCTCNTGNNNNNN
+
AAAAAEEEEEEEEEEEEEEEEEEEEEE#EE######
""")
        # Check original files were removed
        for batch in ('B000','B001','B002',):
            self.assertFalse(os.path.exists(
                os.path.join(self.wd,
                             batch,
                             "PJB1_S1_TTCGTGCA+GATCCAAA_R1_001.fastq")))

