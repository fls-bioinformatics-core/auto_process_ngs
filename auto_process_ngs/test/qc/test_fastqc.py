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
