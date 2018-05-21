#######################################################################
# Unit tests for qc/fastq_screen.py
#######################################################################

import unittest
import tempfile

from auto_process_ngs.qc.fastq_screen import Fastqscreen

class TestFastqscreen(unittest.TestCase):
    def setUp(self):
        screen_text = ""
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            self.fastq_screen_txt = fp.name
            fp.write(screen_text)
    def tearDown(self):
        try:
            os.remove(self.fastq_screen_txt)
        except Exception:
            pass
    def test_fastq_screen_handle_empty_output_file(self):
        """FastqScreen handles empty output file
        """
        self.assertRaises(Exception,
                          Fastqscreen,
                          self.fastq_screen_txt)

class TestFastqscreen_v0_4_1(unittest.TestCase):
    def setUp(self):
        screen_text = """#Fastq_screen version: 0.4.1
Library	%Unmapped	%One_hit_one_library	%Multiple_hits_one_library	%One_hit_multiple_libraries	%Multiple_hits_multiple_libraries
hg19	98.10	0.02	0.27	0.55	1.06
mm9	35.92	47.46	10.18	3.56	2.88
rn4	93.01	0.18	0.17	3.87	2.77
dm3	99.97	0.00	0.00	0.01	0.02
ws200	99.98	0.00	0.00	0.00	0.02
ecoli	96.52	0.43	3.05	0.00	0.00
saccer	99.97	0.00	0.00	0.00	0.03
PhiX	99.33	0.67	0.00	0.00	0.00
Vectors	99.99	0.00	0.00	0.01	0.00
SpR6	100.00	0.00	0.00	0.00	0.00

%Hit_no_libraries: 30.80
"""
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            self.fastq_screen_txt = fp.name
            fp.write(screen_text)
    def tearDown(self):
        try:
            os.remove(self.fastq_screen_txt)
        except Exception:
            pass
    def test_handle_fastq_screen_v0_4_1(self):
        """FastqScreen handles output from v0.4.1
        """
        screen = Fastqscreen(self.fastq_screen_txt)
        self.assertEqual(screen.version,'0.4.1')
        self.assertEqual(screen.txt,self.fastq_screen_txt)
        self.assertEqual(screen.libraries,['hg19','mm9','rn4','dm3','ws200',
                                           'ecoli','saccer','PhiX','Vectors',
                                           'SpR6'])
        self.assertEqual(screen.no_hits,30.80)

class TestFastqscreen_v0_4_2(unittest.TestCase):
    def setUp(self):
        screen_text = """#Fastq_screen version: 0.4.2	#Reads in subset: 1000000
Library	#Reads_processed	#Unmapped	%Unmapped	#One_hit_one_library	%One_hit_one_library	#Multiple_hits_one_library	%Multiple_hits_one_library	#One_hit_multiple_libraries	%One_hit_multiple_libraries	Multiple_hits_multiple_libraries	%Multiple_hits_multiple_libraries
hg19	89393	89213	99.80	1	0.00	0	0.00	11	0.01	168	0.19
mm9	89393	89157	99.74	11	0.01	5	0.01	2	0.00	218	0.24
rn4	89393	89170	99.75	2	0.00	1	0.00	8	0.01	212	0.24
dm3	89393	89391	100.00	0	0.00	0	0.00	0	0.00	2	0.00
ws200	89393	89391	100.00	0	0.00	0	0.00	1	0.00	1	0.00
ecoli	89393	89393	100.00	0	0.00	0	0.00	0	0.00	0	0.00
saccer	89393	89392	100.00	0	0.00	0	0.00	1	0.00	0	0.00
PhiX	89393	89393	100.00	0	0.00	0	0.00	0	0.00	0	0.00
Vectors	89393	89393	100.00	0	0.00	0	0.00	0	0.00	0	0.00
SpR6	89393	89393	100.00	0	0.00	0	0.00	0	0.00	0	0.00

%Hit_no_libraries: 99.73
"""
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            self.fastq_screen_txt = fp.name
            fp.write(screen_text)
    def tearDown(self):
        try:
            os.remove(self.fastq_screen_txt)
        except Exception:
            pass
    def test_handle_fastq_screen_v0_4_2(self):
        """FastqScreen handles output from v0.4.2
        """
        screen = Fastqscreen(self.fastq_screen_txt)
        self.assertEqual(screen.version,'0.4.2')
        self.assertEqual(screen.txt,self.fastq_screen_txt)
        self.assertEqual(screen.libraries,['hg19','mm9','rn4','dm3','ws200',
                                           'ecoli','saccer','PhiX','Vectors',
                                           'SpR6'])
        self.assertEqual(screen.no_hits,99.73)



class TestFastqscreen_v0_5_2(unittest.TestCase):
    def setUp(self):
        screen_text = """#Fastq_screen version: 0.5.2	#Reads in subset: 1000000
Genome	#Reads_processed	#Unmapped	%Unmapped	#One_hit_one_genome	%One_hit_one_genome	#Multiple_hits_one_genome	%Multiple_hits_one_genome	#One_hit_multiple_genomes	%One_hit_multiple_genomes	Multiple_hits_multiple_genomes	%Multiple_hits_multiple_genomes
hg19	597781	118378	19.80	405856	67.89	72519	12.13	577	0.10	451	0.08
mm9	597781	596881	99.85	10	0.00	1	0.00	534	0.09	355	0.06
rn4	597781	596973	99.86	3	0.00	4	0.00	510	0.09	291	0.05
dm3	597781	597731	99.99	0	0.00	0	0.00	19	0.00	31	0.01
ws200	597781	597748	100.00	0	0.00	0	0.00	6	0.00	27	0.00
ecoli	597781	597780	100.00	0	0.00	1	0.00	0	0.00	0	0.00
saccer	597781	597757	100.00	0	0.00	0	0.00	2	0.00	22	0.00
PhiX	597781	597781	100.00	0	0.00	0	0.00	0	0.00	0	0.00
Vectors	597781	597768	100.00	0	0.00	1	0.00	0	0.00	12	0.00
SpR6	597781	597781	100.00	0	0.00	0	0.00	0	0.00	0	0.00

%Hit_no_genomes: 19.80
"""
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            self.fastq_screen_txt = fp.name
            fp.write(screen_text)
    def tearDown(self):
        try:
            os.remove(self.fastq_screen_txt)
        except Exception:
            pass
    def test_handle_fastq_screen_v0_5_2(self):
        """FastqScreen handles output from v0.5.2
        """
        screen = Fastqscreen(self.fastq_screen_txt)
        self.assertEqual(screen.version,'0.5.2')
        self.assertEqual(screen.txt,self.fastq_screen_txt)
        self.assertEqual(screen.libraries,['hg19','mm9','rn4','dm3','ws200',
                                           'ecoli','saccer','PhiX','Vectors',
                                           'SpR6'])
        self.assertEqual(screen.no_hits,19.80)
