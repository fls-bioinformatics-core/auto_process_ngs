#######################################################################
# Unit tests for qc/fastqc.py
#######################################################################

import unittest
import os
import shutil
import tempfile

from auto_process_ngs.qc.fastqc import Fastqc
from auto_process_ngs.qc.fastqc import FastqcSummary
from auto_process_ngs.qc.fastqc import FastqcData

class TestFastqc_v0_11_3(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestFastqc')
        # Make a directory for fake Fastqc output
        self.fastqc_dir = os.path.join(self.dirn,
                                       "PJB_S1_L001_R1_001_fastqc")
        os.mkdir(self.fastqc_dir)
        # Make a fake Fastqc HTML file
        fastqc_html = ""
        self.fastqc_html = "%s.html" % self.fastqc_dir
        with open(self.fastqc_html,'w') as fp:
            fp.write(fastqc_html)
        # Make a fake Fastqc data file
        fastqc_data = """##FastQC	0.11.3
>>Basic Statistics	pass
#Measure	Value
Filename	PJB_S1_L001_R1_001.fastq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	2973921
Sequences flagged as poor quality	0
Sequence length	35-75
%GC	41
>>END_MODULE
>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Per
centile
1	33.653393281126164	34.0	34.0	34.0	34.0	34.0
2	33.693315995952815	34.0	34.0	34.0	34.0	34.0
3	33.74285093652454	34.0	34.0	34.0	34.0	34.0
4	33.75196012267979	34.0	34.0	34.0	34.0	34.0
>>END_MODULE
>>Per tile sequence quality	fail
#Tile	Base	Mean
1101	1	-0.27315747640444243
1101	2	-0.04655524562563329
1101	3	-0.044408092234554886
1101	4	0.02674347384745346
>>END_MODULE
>>Per sequence quality scores	pass
#Quality	Count
2	70.0
3	0.0
4	0.0
>>END_MODULE
>>Per base sequence content	fail
#Base	G	A	T	C
1	38.482761913760974	23.707004823039217	14.204713013530268	23.60552024966954
2	13.34707085190213	22.24691149623838	42.94041631541056	21.465601336448934
3	15.901839063221393	28.581290723711444	31.732625474510996	23.78424473855617
4	11.473305152141112	18.386798800612407	41.42588179434679	28.714014252899688
>>END_MODULE
>>Per sequence GC content	warn
#GC Content	Count
0	70.0
1	35.0
2	0.5
3	1.0
4	1.5
>>END_MODULE
>>Per base N content	pass
#Base	N-Count
1	0.002353794872156994
2	0.002353794872156994
3	0.002353794872156994
4	0.002353794872156994
>>END_MODULE
>>Sequence Length Distribution	warn
#Length	Count
35	8826.0
36	2848.0
37	4666.0
38	4524.0
>>END_MODULE
>>Sequence Duplication Levels	pass
#Total Deduplicated Percentage	93.95730549982757
#Duplication Level	Percentage of deduplicated	Percentage of total
1	97.73363240226425	91.82788757227388
2	1.3763867023900318	2.5864317176472094
3	0.3364050209524425	0.9482312797591361
4	0.16969367993881057	0.6377584370960314
>>END_MODULE
>>Overrepresented sequences	pass
>>END_MODULE
>>Adapter Content	pass
#Position	Illumina Universal Adapter	Illumina Small RNA Adapter	
Nextera Transposase Sequence	SOLID Small RNA Adapter
1	0.0	0.0	0.0	0.0
2	0.0	0.0	0.0	0.0
3	0.0	0.0	0.0	0.0
4	0.0	0.0	3.36256410308142E-5	0.0
>>END_MODULE
>>Kmer Content	fail
#Sequence	Count	PValue	Obs/Exp Max	Max Obs/Exp Position
CTATACT	2310	0.0	22.463886	4
TATACTG	6135	0.0	21.145714	5
TCTATAC	2460	0.0	21.094137	3
AGTACCC	4205	0.0	19.634653	27
>>END_MODULE
"""
        self.fastqc_data = os.path.join(self.fastqc_dir,
                                        "fastqc_data.txt")
        with open(self.fastqc_data,'w') as fp:
            fp.write(fastqc_data)
        # Make a fake Fastqc summary file
        fastqc_summary = """PASS	Basic Statistics	PJB_S1_L001_R1_001.fastq.gz
PASS	Per base sequence quality	PJB_S1_L001_R1_001.fastq.gz
FAIL	Per tile sequence quality	PJB_S1_L001_R1_001.fastq.gz
PASS	Per sequence quality scores	PJB_S1_L001_R1_001.fastq.gz
FAIL	Per base sequence content	PJB_S1_L001_R1_001.fastq.gz
WARN	Per sequence GC content	PJB_S1_L001_R1_001.fastq.gz
PASS	Per base N content	PJB_S1_L001_R1_001.fastq.gz
WARN	Sequence Length Distribution	PJB_S1_L001_R1_001.fastq.gz
PASS	Sequence Duplication Levels	PJB_S1_L001_R1_001.fastq.gz
PASS	Overrepresented sequences	PJB_S1_L001_R1_001.fastq.gz
PASS	Adapter Content	PJB_S1_L001_R1_001.fastq.gz
FAIL	Kmer Content	PJB_S1_L001_R1_001.fastq.gz
"""
        self.fastqc_summary = os.path.join(self.fastqc_dir,
                                           "summary.txt")
        with open(self.fastqc_summary,'w') as fp:
            fp.write(fastqc_summary)
        # Make empty subdirs
        for d in ('Icons','Images',):
            os.mkdir(os.path.join(self.fastqc_dir,d))
        # Make empty files
        for f in ('fastqc.fo','fastqc_report.html',):
            with open(os.path.join(self.fastqc_dir,f),'w') as fp:
                fp.write("")
        # Make empty zip archive
        self.fastqc_zip = "%s.zip" % self.fastqc_dir
        with open(self.fastqc_zip,'w') as fp:
            fp.write("")
    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)
    def test_handle_fastqc_v0_11_3(self):
        """Fastqc handles output from v0.11.3
        """
        # Test the top-level instance
        fastqc = Fastqc(self.fastqc_dir)
        self.assertEqual(fastqc.version,"0.11.3")
        self.assertEqual(fastqc.dir,self.fastqc_dir)
        self.assertEqual(fastqc.html_report,self.fastqc_html)
        self.assertEqual(fastqc.zip,self.fastqc_zip)
        # Test the summary instance
        summary = fastqc.summary
        self.assertTrue(isinstance(summary,FastqcSummary))
        self.assertEqual(summary.path,self.fastqc_summary)
        self.assertEqual(summary.modules,['Basic Statistics',
                                          'Per base sequence quality',
                                          'Per tile sequence quality',
                                          'Per sequence quality scores',
                                          'Per base sequence content',
                                          'Per sequence GC content',
                                          'Per base N content',
                                          'Sequence Length Distribution',
                                          'Sequence Duplication Levels',
                                          'Overrepresented sequences',
                                          'Adapter Content',
                                          'Kmer Content',])
        self.assertEqual(summary.passes,['Basic Statistics',
                                         'Per base sequence quality',
                                         'Per sequence quality scores',
                                         'Per base N content',
                                         'Sequence Duplication Levels',
                                         'Overrepresented sequences',
                                         'Adapter Content',])
        self.assertEqual(summary.warnings,['Per sequence GC content',
                                           'Sequence Length Distribution',])
        self.assertEqual(summary.failures,['Per tile sequence quality',
                                           'Per base sequence content',
                                           'Kmer Content',])
        # Test the data instance
        data = fastqc.data
        self.assertTrue(isinstance(data,FastqcData))
        self.assertEqual(data.version,"0.11.3")
        self.assertEqual(data.path,self.fastqc_data)
        # Data for a module
        self.assertEqual(data.data('Sequence Length Distribution'),
                         ['#Length\tCount',
                          '35\t8826.0',
                          '36\t2848.0',
                          '37\t4666.0',
                          '38\t4524.0'])
        self.assertEqual(data.data('Unknown module'),None)
        # Basic statistics data
        self.assertEqual(data.basic_statistics('Filename'),'PJB_S1_L001_R1_001.fastq.gz')
        self.assertEqual(data.basic_statistics('File type'),'Conventional base calls')
        self.assertEqual(data.basic_statistics('Encoding'),'Sanger / Illumina 1.9')
        self.assertEqual(data.basic_statistics('Total Sequences'),'2973921')
        self.assertEqual(data.basic_statistics('Sequences flagged as poor quality'),'0')
        self.assertEqual(data.basic_statistics('Sequence length'),'35-75')
        self.assertEqual(data.basic_statistics('%GC'),'41')
        self.assertRaises(Exception,
                          data.basic_statistics,
                          "unknown_thing")

class TestFastqc_v0_12_1(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestFastqc')
        # Make a directory for fake Fastqc output
        self.fastqc_dir = os.path.join(self.dirn,
                                       "PJB_S1_L001_R1_001_fastqc")
        os.mkdir(self.fastqc_dir)
        # Make a fake Fastqc HTML file
        fastqc_html = ""
        self.fastqc_html = "%s.html" % self.fastqc_dir
        with open(self.fastqc_html,'w') as fp:
            fp.write(fastqc_html)
        # Make a fake Fastqc data file
        fastqc_data = """##FastQC	0.12.1
>>Basic Statistics	pass
#Measure	Value
Filename	PJB_S1_L001_R1_001.fastq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	379988
Total Bases	28.6 Mbp
Sequences flagged as poor quality	0
Sequence length	35-76
%GC	58
>>END_MODULE
>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile
1	34.265750497384126	37.0	32.0	37.0	32.0	37.0
2	34.7533895807236	37.0	32.0	37.0	32.0	37.0
3	35.59621093297683	37.0	37.0	37.0	32.0	37.0
4	36.022850721601735	37.0	37.0	37.0	37.0	37.0
>>END_MODULE
>>Per tile sequence quality	pass
#Tile	Base	Mean
11101	1	0.014892076641700669
11101	2	0.5392987472679067
11101	3	0.4010864755168484
11101	4	0.2713441615718324
>>END_MODULE
>>Per sequence quality scores	pass
#Quality	Count
2	378.0
3	0.0
4	0.0
>>END_MODULE
>>Per base sequence content	fail
#Base	G	A	T	C
1	61.87297489528727	8.874107636785121	4.119491056610732	25.13342641131688
2	31.41159889885539	9.43875870971141	30.063882193332542	29.085760198100658
3	27.028636160067443	15.680866197739665	22.199214942437894	35.091282699755
4	34.08864581192798	20.47905832266359	19.68212067085026	25.750175194558167
>>END_MODULE
>>Per sequence GC content	fail
#GC Content	Count
0	383.0
1	194.0
2	2.5
3	2.0
4	4.0
5	5.5
6	10.5
7	31.0
8	48.0
9	27.5
10	9.0
11	11.0
12	11.0
13	27.5
14	66.0
15	97.5
16	112.5
17	97.5
18	194.0
19	405.0
20	487.5
21	386.5
22	308.0
23	386.5
24	441.0
25	463.0
26	619.0
27	781.0
28	933.5
29	1116.0
30	1723.0
31	2366.0
32	2576.0
33	2833.0
34	3847.5
35	4400.5
36	3855.0
37	3390.0
38	3451.5
39	3811.5
40	4494.0
41	5475.5
42	6020.5
43	7011.0
44	8165.0
45	10546.0
46	12945.5
47	15506.0
48	15858.0
49	14363.5
50	14618.5
51	16301.5
52	17653.5
53	17592.0
54	17348.5
55	15777.5
56	15497.0
57	16449.0
58	16966.0
59	16200.5
60	14229.5
61	12520.5
62	11457.5
63	10760.5
64	10181.5
65	8558.5
66	6979.5
67	6524.5
68	7462.5
69	7558.0
70	5574.0
71	4362.0
72	4217.5
73	5378.0
74	5561.0
75	5156.0
76	5686.5
77	5401.5
78	4605.5
79	4628.5
80	5793.5
81	5343.0
82	4179.5
83	4405.0
84	6078.5
85	6786.5
86	3908.5
87	1902.5
88	2359.0
89	3164.5
90	2408.5
91	905.0
92	536.5
93	489.5
94	352.5
95	226.0
96	190.0
97	176.0
98	137.0
99	69.0
100	28.0
>>END_MODULE
>>Per base N content	pass
#Base	N-Count
1	0.09947682558396581
2	0.10079265661020874
3	0.10474014968893755
4	0.10789814415192059
>>END_MODULE
>>Sequence Length Distribution	warn
#Length	Count
35	524.0
36	11.0
37	19.0
38	23.0
>>END_MODULE
>>Sequence Duplication Levels	fail
#Total Deduplicated Percentage	39.0958974034817
#Duplication Level	Percentage of total
1	34.413393047069995
2	4.2292688859395104
3	1.9722575226717567
4	1.2680464075212772
5	1.0428983513617631
6	0.9202475480158026
7	0.717802189536058
8	0.7082105897502663
9	0.6733068805951249
>10	15.273708710399417
>50	9.697387969765732
>100	22.353220725266553
>500	4.946018179000487
>1k	1.784232993106257
>5k	0.0
>10k+	0.0
>>END_MODULE
>>Overrepresented sequences	warn
#Sequence	Count	Percentage	Possible Source
GGCGTACGGAAGACCCGCTCCCCGGCGCCGCTCGTGGGGGGCCCAAGTCC	2067	0.5439645462488288	No Hit
GTACGGAAGACCCGCTCCCCGGCGCCGCTCGTGGGGGGCCCAAGTCCTTC	1273	0.335010579281451	No Hit
GTTGGATTGTTCACCCACTAATAGGGAACGTGAGCTGGGTTTAGACCGTC	1201	0.3160626125035528	No Hit
GGGAGTTTGACTGGGGCGGTACACCTGTCAAACGGTAACGCAGGTGTCCT	1158	0.3047464656778635	No Hit
>>END_MODULE
>>Adapter Content	pass
#Position	Illumina Universal Adapter	Illumina Small RNA 3' Adapter	Illumina Small RNA 5' Adapter	Nextera Transposase Sequence	PolyA	PolyG
1	0.0	0.0	0.0	0.0	0.0010526648209943472	0.0
2	0.0	0.0	2.631662052485868E-4	0.0	0.0021053296419886944	0.0
3	0.0	0.0	2.631662052485868E-4	0.0	0.0021053296419886944	5.263324104971736E-4
4	0.0	0.0	2.631662052485868E-4	0.0	0.0028948282577344548	7.894986157457604E-4
>>END_MODULE
"""
        self.fastqc_data = os.path.join(self.fastqc_dir,
                                        "fastqc_data.txt")
        with open(self.fastqc_data,'w') as fp:
            fp.write(fastqc_data)
        # Make a fake Fastqc summary file
        fastqc_summary = """PASS	Basic Statistics	PJB_S1_L001_R1_001.fastq.gz
PASS	Per base sequence quality	PJB_S1_L001_R1_001.fastq.gz
PASS	Per tile sequence quality	PJB_S1_L001_R1_001.fastq.gz
PASS	Per sequence quality scores	PJB_S1_L001_R1_001.fastq.gz
FAIL	Per base sequence content	PJB_S1_L001_R1_001.fastq.gz
FAIL	Per sequence GC content	PJB_S1_L001_R1_001.fastq.gz
PASS	Per base N content	PJB_S1_L001_R1_001.fastq.gz
WARN	Sequence Length Distribution	PJB_S1_L001_R1_001.fastq.gz
FAIL	Sequence Duplication Levels	PJB_S1_L001_R1_001.fastq.gz
WARN	Overrepresented sequences	PJB_S1_L001_R1_001.fastq.gz
PASS	Adapter Content	PJB_S1_L001_R1_001.fastq.gz
"""
        self.fastqc_summary = os.path.join(self.fastqc_dir,
                                           "summary.txt")
        with open(self.fastqc_summary,'w') as fp:
            fp.write(fastqc_summary)
        # Make empty subdirs
        for d in ('Icons','Images',):
            os.mkdir(os.path.join(self.fastqc_dir,d))
        # Make empty files
        for f in ('fastqc.fo','fastqc_report.html',):
            with open(os.path.join(self.fastqc_dir,f),'w') as fp:
                fp.write("")
        # Make empty zip archive
        self.fastqc_zip = "%s.zip" % self.fastqc_dir
        with open(self.fastqc_zip,'w') as fp:
            fp.write("")
    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)
    def test_handle_fastqc_v0_12_1(self):
        """Fastqc handles output from v0.12.1
        """
        # Test the top-level instance
        fastqc = Fastqc(self.fastqc_dir)
        self.assertEqual(fastqc.version,"0.12.1")
        self.assertEqual(fastqc.dir,self.fastqc_dir)
        self.assertEqual(fastqc.html_report,self.fastqc_html)
        self.assertEqual(fastqc.zip,self.fastqc_zip)
        # Test the summary instance
        summary = fastqc.summary
        self.assertTrue(isinstance(summary,FastqcSummary))
        self.assertEqual(summary.path,self.fastqc_summary)
        self.assertEqual(summary.modules,['Basic Statistics',
                                          'Per base sequence quality',
                                          'Per tile sequence quality',
                                          'Per sequence quality scores',
                                          'Per base sequence content',
                                          'Per sequence GC content',
                                          'Per base N content',
                                          'Sequence Length Distribution',
                                          'Sequence Duplication Levels',
                                          'Overrepresented sequences',
                                          'Adapter Content',])
        self.assertEqual(summary.passes,['Basic Statistics',
                                         'Per base sequence quality',
                                         'Per tile sequence quality',
                                         'Per sequence quality scores',
                                         'Per base N content',
                                         'Adapter Content',])
        self.assertEqual(summary.warnings,['Sequence Length Distribution',
                                           'Overrepresented sequences',])
        self.assertEqual(summary.failures,['Per base sequence content',
                                           'Per sequence GC content',
                                           'Sequence Duplication Levels',])
        # Test the data instance
        data = fastqc.data
        self.assertTrue(isinstance(data,FastqcData))
        self.assertEqual(data.version,"0.12.1")
        self.assertEqual(data.path,self.fastqc_data)
        # Data for a module
        self.assertEqual(data.data('Sequence Length Distribution'),
                         ['#Length\tCount',
                          '35\t524.0',
                          '36\t11.0',
                          '37\t19.0',
                          '38\t23.0'])
        self.assertEqual(data.data('Unknown module'),None)
        # Basic statistics data
        self.assertEqual(data.basic_statistics('Filename'),'PJB_S1_L001_R1_001.fastq.gz')
        self.assertEqual(data.basic_statistics('File type'),'Conventional base calls')
        self.assertEqual(data.basic_statistics('Encoding'),'Sanger / Illumina 1.9')
        self.assertEqual(data.basic_statistics('Total Sequences'),'379988')
        self.assertEqual(data.basic_statistics('Sequences flagged as poor quality'),'0')
        self.assertEqual(data.basic_statistics('Sequence length'),'35-76')
        self.assertEqual(data.basic_statistics('%GC'),'58')
        self.assertRaises(Exception,
                          data.basic_statistics,
                          "unknown_thing")
