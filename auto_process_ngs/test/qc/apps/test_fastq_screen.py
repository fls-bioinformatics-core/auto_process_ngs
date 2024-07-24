#######################################################################
# Unit tests for qc/fastq_screen.py
#######################################################################

import os
import unittest
import tempfile

from auto_process_ngs.qc.apps.fastq_screen import Fastqscreen
from auto_process_ngs.qc.apps.fastq_screen import fastq_screen_output_files

class TestFastqscreen(unittest.TestCase):
    def setUp(self):
        screen_text = ""
        with tempfile.NamedTemporaryFile(mode='wt',delete=False) as fp:
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
        with tempfile.NamedTemporaryFile(mode='wt',delete=False) as fp:
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
        with tempfile.NamedTemporaryFile(mode='wt',delete=False) as fp:
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
        with tempfile.NamedTemporaryFile(mode='wt',delete=False) as fp:
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

class TestFastqscreen_v0_15_3(unittest.TestCase):
    def setUp(self):
        screen_text = """#Fastq_screen version: 0.15.3	#Aligner: bowtie	#Reads in subset: 100000
Genome	#Reads_processed	#Unmapped	%Unmapped	#One_hit_one_genome	%One_hit_one_genome	#Multiple_hits_one_genome	%Multiple_hits_one_genome	#One_hit_multiple_genomes	%One_hit_multiple_genomes	Multiple_hits_multiple_genomes	%Multiple_hits_multiple_genomes
hg19	94997	43983	46.30	6630	6.98	14774	15.55	3648	3.84	25962	27.33
mm9	94997	72416	76.24	4	0.00	1	0.00	5767	6.07	16809	17.69
rn4	94997	65932	69.41	1	0.00	0	0.00	2501	2.63	26563	27.96
dm3	94997	55447	58.37	26266	27.65	11117	11.70	95	0.10	2072	2.18
ws200	94997	93187	98.09	0	0.00	0	0.00	634	0.67	1176	1.24
ecoli	94997	94993	100.00	0	0.00	4	0.00	0	0.00	0	0.00
saccer	94997	92724	97.60	7	0.01	5	0.01	0	0.00	2261	2.38
PhiX	94997	94997	100.00	0	0.00	0	0.00	0	0.00	0	0.00
Vectors	94997	94990	100.00	1	0.00	4	0.00	2	0.00	0	0.00
SpR6	94997	94994	100.00	0	0.00	3	0.00	0	0.00	0	0.00

%Hit_no_genomes: 6.90
"""
        with tempfile.NamedTemporaryFile(mode='wt',delete=False) as fp:
            self.fastq_screen_txt = fp.name
            fp.write(screen_text)
    def tearDown(self):
        try:
            os.remove(self.fastq_screen_txt)
        except Exception:
            pass
    def test_handle_fastq_screen_v0_15_32(self):
        """FastqScreen handles output from v0.15.3
        """
        screen = Fastqscreen(self.fastq_screen_txt)
        self.assertEqual(screen.version,'0.15.3')
        self.assertEqual(screen.txt,self.fastq_screen_txt)
        self.assertEqual(screen.libraries,['hg19','mm9','rn4','dm3','ws200',
                                           'ecoli','saccer','PhiX','Vectors',
                                           'SpR6'])
        self.assertEqual(screen.no_hits,6.90)

class TestFastqScreenOutputFilesFunction(unittest.TestCase):

    def test_fastq_screen_output_files(self):
        """
        fastq_screen_output_files: handles .fastq file
        """
        self.assertEqual(
            fastq_screen_output_files(
                '/data/PB/PB1_ATTAGG_L001_R1_001.fastq',
                'model_organisms'),
            ('PB1_ATTAGG_L001_R1_001_screen_model_organisms.png',
             'PB1_ATTAGG_L001_R1_001_screen_model_organisms.txt'))

    def test_fastq_screen_output_files_fastqgz(self):
        """
        fastq_screen_output_files: handles fastq.gz file
        """
        self.assertEqual(
            fastq_screen_output_files(
                '/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz',
                'model_organisms'),
            ('PB1_ATTAGG_L001_R1_001_screen_model_organisms.png',
             'PB1_ATTAGG_L001_R1_001_screen_model_organisms.txt'))

    def test_fastq_screen_output_files_legacy_naming(self):
        """
        fastq_screen_output_files: handles legacy naming convention
        """
        self.assertEqual(
            fastq_screen_output_files(
                '/data/PB/PB1_ATTAGG_L001_R1_001.fastq.gz',
                'model_organisms',
                legacy=True),
            ('PB1_ATTAGG_L001_R1_001_model_organisms_screen.png',
             'PB1_ATTAGG_L001_R1_001_model_organisms_screen.txt'))
