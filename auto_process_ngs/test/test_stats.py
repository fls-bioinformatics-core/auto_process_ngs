#######################################################################
# Tests for stats.py module
#######################################################################

import os
import unittest
import tempfile
import shutil
import gzip
import cStringIO
from bcftbx.mock import MockIlluminaData
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.IlluminaData import IlluminaData
from auto_process_ngs.stats import FastqStatistics
from auto_process_ngs.stats import FastqStats
from auto_process_ngs.stats import FastqReadCounter
from auto_process_ngs.stats import collect_fastq_data

# Test data
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
fastq_r1_data = """@HISEQ:1:000000000-A2Y1L:1:1101:19264:2433 1:N:0:AGATCGC
AGATAGCCGA
+
?????BBB@B
@HISEQ:1:000000000-A2Y1L:1:1101:18667:2435 1:N:0:AGATCGC
ATATATTCAT
+
?????BBBDD
@HISEQ:1:000000000-A2Y1L:1:1101:17523:2436 1:N:0:AGATCGC
CATCACTACC
+
?<,<?BBBBB
"""
fastq_r2_data = """@HISEQ:1:000000000-A2Y1L:1:1101:19264:2433 2:N:0:AGATCGC
GCCGATATGC
+
??A??ABBDD
@HISEQ:1:000000000-A2Y1L:1:1101:18667:2435 2:N:0:AGATCGC
GATGACATCA
+
?????BBBDD
@HISEQ:1:000000000-A2Y1L:1:1101:17523:2436 2:N:0:AGATCGC
GAATATAGAA
+
??AAABBBDD
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

# Mock classes
class AugmentedMockIlluminaData(MockIlluminaData):
    """
    Modified  MockIlluminaData to add content to FASTQ files

    This version adds a 'populate_fastq' method that can
    be used to generate or append content to a FASTQ file.
    """
    def __init__(self,name,package,**kws):
        MockIlluminaData.__init__(self,name,package,**kws)
    def populate_fastq(self,path,nreads,lane=None,append=False):
        """
        Write or append content to a FASTQ file

        Writes the specified number of reads to the named
        FASTQ file.

        Read headers will be written to have the appropriate
        read number (extracted from the FASTQ name) and
        lane number.

        The sequence and quality data in the reads are
        identical for every read.

        Arguments:
          path (str): path to the FASTQ file
          nreads (int): number of reads to generate
          lane (int): optional, explicitly specify the lane
            number for read headers (by default will be
            extracted from the FASTQ name, if possible)
          append (boolean): if True then append the reads to
            the FASTQ file (default is to overwrite any
            existing reads)

        Raises:
          Exception if lane is not supplied and can't be
            extracted from the FASTQ name.
        """
        # Extract read number from name
        read_number = IlluminaFastq(path).read_number
        # Extract lane from name, if not explicitly specified
        if lane is None:
            lane = IlluminaFastq(path).lane_number
            if lane is None:
                raise Exception("Couldn't get lane from FASTQ name: %s" %
                                path)
        # Create or append fake reads to a FASTQ file
        fastq = os.path.join(self.unaligned_dir,path)
        if append:
            mode = 'ab'
        else:
            mode = 'wb'
        # Write the required number of reads to the file
        with gzip.GzipFile(fastq,mode) as fq:
            if read_number == 1:
                r = """@HISEQ:1:000000000-A2Y1L:%s:1101:19264:2433 1:N:0:AGATCGC
AGATAGCCGA
+
?????BBB@B
""" % lane
            else:
                r = """@HISEQ:1:000000000-A2Y1L:%s:1101:19264:2433 2:N:0:AGATCGC
GCCGATATGC
+
??A??ABBDD
""" % lane
            for i in xrange(nreads):
                fq.write(r)

# FastqStatistics
class TestFastqStatisticsCasava(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestFastqStats')
    def tearDown(self):
        # Remove the temporary test directory
        #shutil.rmtree(self.dirn)
        pass
    def _setup_casava(self):
        # Create a mock casava dir structure
        mock_data = AugmentedMockIlluminaData(
            '151125_S00879_0001_000000000-ABCDE1_analysis',
            'casava',
            unaligned_dir='bcl2fastq',
            paired_end=True,
            top_dir=self.dirn)
        mock_data.add_fastq_batch('AB','AB1','AB1_GCCAAT',lanes=[1,2])
        mock_data.add_fastq_batch('AB','AB2','AB2_AGTCAA',lanes=[1,2])
        mock_data.add_fastq_batch('CDE','CDE3','CDE3_GCCAAT',lanes=[3,4])
        mock_data.add_fastq_batch('CDE','CDE4','CDE4_AGTCAA',lanes=[3,4])
        mock_data.add_undetermined(lanes=(1,2,3,4))
        # Create on disk
        mock_data.create()
        # Populate the FASTQs with 'fake' reads
        reads = {
            "AB": {
                "AB1": { 1:3, 2:2 },
                "AB2": { 1:5, 2:3 },
            },
            "CDE": {
                "CDE3": { 3:7, 4:2 },
                "CDE4": { 3:8, 4:6 },
            },
            "Undetermined_indices": {
                "lane1": { 1: 2 },
                "lane2": { 2: 1 },
                "lane3": { 3: 4 },
                "lane4": { 4: 3 },
            },
        }
        barcodes = {
            "AB1": "GCCAAT",
            "AB2": "AGTCAA",
            "CDE3": "GCCAAT",
            "CDE4": "AGTCAA"
        }
        for project in reads:
            for sample in reads[project]:
                for lane in reads[project][sample]:
                    for read_number in (1,2):
                        try:
                            barcode = barcodes[sample]
                        except KeyError:
                            barcode = "Undetermined"
                        nreads = reads[project][sample][lane]
                        if project != "Undetermined_indices":
                            project_name = "Project_%s" % project
                        else:
                            project_name = project
                        sample_name = "Sample_%s" % sample
                        fastq = os.path.join(
                            project_name,sample_name,
                            "%s_%s_L%03d_R%d_001.fastq.gz" %
                            (sample,barcode,lane,read_number))
                        mock_data.populate_fastq(fastq,nreads)
        # Store the location of the mock data
        self.illumina_data = mock_data.dirn
    def test_fastqstatistics_casava(self):
        self._setup_casava()
        fqstatistics = FastqStatistics(
            IlluminaData(
                self.illumina_data,
                unaligned_dir="bcl2fastq"))
        self.assertEqual(fqstatistics.lane_names,
                         ['L1','L2','L3','L4'])
        # Check "raw" stored data
        self.assertEqual(fqstatistics.raw.header(),
                         ['Project',
                          'Sample',
                          'Fastq',
                          'Size',
                          'Nreads',
                          'Paired_end',
                          'Read_number',
                          'L1','L2','L3','L4'])
        self.assertEqual(len(fqstatistics.raw),24)
        expected = [
            ['AB','AB1','AB1_GCCAAT_L001_R1_001.fastq.gz',3,{'L1':3}],
            ['AB','AB1','AB1_GCCAAT_L001_R2_001.fastq.gz',3,{'L1':3}],
            ['AB','AB1','AB1_GCCAAT_L002_R1_001.fastq.gz',2,{'L2':2}],
            ['AB','AB1','AB1_GCCAAT_L002_R2_001.fastq.gz',2,{'L2':2}],
            ['AB','AB2','AB2_AGTCAA_L001_R1_001.fastq.gz',5,{'L1':5}],
            ['AB','AB2','AB2_AGTCAA_L001_R2_001.fastq.gz',5,{'L1':5}],
            ['AB','AB2','AB2_AGTCAA_L002_R1_001.fastq.gz',3,{'L2':3}],
            ['AB','AB2','AB2_AGTCAA_L002_R2_001.fastq.gz',3,{'L2':3}],
            ['CDE','CDE3','CDE3_GCCAAT_L003_R1_001.fastq.gz',7,{'L3':7}],
            ['CDE','CDE3','CDE3_GCCAAT_L003_R2_001.fastq.gz',7,{'L3':7}],
            ['CDE','CDE3','CDE3_GCCAAT_L004_R1_001.fastq.gz',2,{'L4':2}],
            ['CDE','CDE3','CDE3_GCCAAT_L004_R2_001.fastq.gz',2,{'L4':2}],
            ['CDE','CDE4','CDE4_AGTCAA_L003_R1_001.fastq.gz',8,{'L3':8}],
            ['CDE','CDE4','CDE4_AGTCAA_L003_R2_001.fastq.gz',8,{'L3':8}],
            ['CDE','CDE4','CDE4_AGTCAA_L004_R1_001.fastq.gz',6,{'L4':6}],
            ['CDE','CDE4','CDE4_AGTCAA_L004_R2_001.fastq.gz',6,{'L4':6}],
            ['Undetermined_indices','lane1',
             'lane1_Undetermined_L001_R1_001.fastq.gz',2,{'L1':2}],
            ['Undetermined_indices','lane1',
             'lane1_Undetermined_L001_R2_001.fastq.gz',2,{'L1':2}],
            ['Undetermined_indices','lane2',
             'lane2_Undetermined_L002_R1_001.fastq.gz',1,{'L2':1}],
            ['Undetermined_indices','lane2',
             'lane2_Undetermined_L002_R2_001.fastq.gz',1,{'L2':1}],
            ['Undetermined_indices','lane3',
             'lane3_Undetermined_L003_R1_001.fastq.gz',4,{'L3':4}],
            ['Undetermined_indices','lane3',
             'lane3_Undetermined_L003_R2_001.fastq.gz',4,{'L3':4}],
            ['Undetermined_indices','lane4',
             'lane4_Undetermined_L004_R1_001.fastq.gz',3,{'L4':3}],
            ['Undetermined_indices','lane4',
             'lane4_Undetermined_L004_R2_001.fastq.gz',3,{'L4':3}],
        ]
        for line,expctd in zip(fqstatistics.raw,expected):
            self.assertEqual(line['Project'],expctd[0])
            self.assertEqual(line['Sample'],expctd[1])
            self.assertEqual(line['Fastq'],expctd[2])
            self.assertEqual(line['Nreads'],expctd[3])
            for lane in ('L1','L2','L3','L4'):
                if lane in expctd[4]:
                    self.assertEqual(line[lane],expctd[4][lane])
                else:
                    self.assertEqual(line[lane],'')
            self.assertEqual(line['Read_number'],
                             IlluminaFastq(expctd[2]).read_number)
    def test_report_per_lane_sample_stats(self):
        fp = cStringIO.StringIO()
        self._setup_casava()
        fqstatistics = FastqStatistics(
            IlluminaData(
                self.illumina_data,
                unaligned_dir="bcl2fastq"))
        fqstatistics.report_per_lane_sample_stats(fp=fp)
        self.assertEqual(fp.getvalue(),"""
Lane 1
Total reads = 10
- AB/AB1	3	30.0%
- AB/AB2	5	50.0%
- Undetermined_indices/lane1	2	20.0%

Lane 2
Total reads = 6
- AB/AB1	2	33.3%
- AB/AB2	3	50.0%
- Undetermined_indices/lane2	1	16.7%

Lane 3
Total reads = 19
- CDE/CDE3	7	36.8%
- CDE/CDE4	8	42.1%
- Undetermined_indices/lane3	4	21.1%

Lane 4
Total reads = 11
- CDE/CDE3	2	18.2%
- CDE/CDE4	6	54.5%
- Undetermined_indices/lane4	3	27.3%
""")
    def test_report_per_lane_summary_stats(self):
        fp = cStringIO.StringIO()
        self._setup_casava()
        fqstatistics = FastqStatistics(
            IlluminaData(
                self.illumina_data,
                unaligned_dir="bcl2fastq"))
        fqstatistics.report_per_lane_summary_stats(fp=fp)
        self.assertEqual(fp.getvalue(),"""#Lane	Total reads	Assigned reads	Unassigned reads	%assigned	%unassigned
Lane 1	10	8	2	80.0	20.0
Lane 2	6	5	1	83.33	16.67
Lane 3	19	15	4	78.95	21.05
Lane 4	11	8	3	72.73	27.27
""")

class TestFastqStatisticsBcl2fastq2(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestFastqStats')
    def tearDown(self):
        # Remove the temporary test directory
        #shutil.rmtree(self.dirn)
        pass
    def _setup_bcl2fastq2(self):
        # Create a mock bcl2fastq dir structure
        mock_data = AugmentedMockIlluminaData(
            '151125_S00879_0001_000000000-ABCDE1_analysis',
            'bcl2fastq2',
            unaligned_dir='bcl2fastq',
            paired_end=True,
            no_lane_splitting=False,
            top_dir=self.dirn)
        mock_data.add_fastq_batch('AB','AB1','AB1_S1',lanes=[1,2])
        mock_data.add_fastq_batch('AB','AB2','AB2_S2',lanes=[1,2])
        mock_data.add_fastq_batch('CDE','CDE3','CDE3_S3',lanes=[3,4])
        mock_data.add_fastq_batch('CDE','CDE4','CDE4_S4',lanes=[3,4])
        mock_data.add_undetermined(lanes=(1,2,3,4))
        # Create on disk
        mock_data.create()
        # Populate the FASTQs with 'fake' reads
        reads = {
            "AB": {
                "AB1": { 1:3, 2:2 },
                "AB2": { 1:5, 2:3 },
            },
            "CDE": {
                "CDE3": { 3:7, 4:2 },
                "CDE4": { 3:8, 4:6 },
            },
            "Undetermined_indices": {
                "Undetermined": { 1: 2, 2: 1, 3: 4, 4: 3 },
            },
        }
        s_indices = {
            "AB1": "S1",
            "AB2": "S2",
            "CDE3": "S3",
            "CDE4": "S4"
        }
        for project in reads:
            for sample in reads[project]:
                for lane in reads[project][sample]:
                    for read_number in (1,2):
                        try:
                            s_index = s_indices[sample]
                        except KeyError:
                            s_index = "S0"
                        nreads = reads[project][sample][lane]
                        if project != "Undetermined_indices":
                            project_name = project
                        else:
                            project_name = ""
                        fastq = os.path.join(
                            project_name,
                            "%s_%s_L%03d_R%d_001.fastq.gz" %
                            (sample,s_index,lane,read_number))
                        mock_data.populate_fastq(fastq,nreads)
        # Store the location of the mock data
        self.illumina_data = mock_data.dirn
    def test_fastqstatistics_bcl2fastq2(self):
        self._setup_bcl2fastq2()
        fqstatistics = FastqStatistics(
            IlluminaData(
                self.illumina_data,
                unaligned_dir="bcl2fastq"))
        self.assertEqual(fqstatistics.lane_names,
                         ['L1','L2','L3','L4'])
        # Check "raw" stored data
        self.assertEqual(fqstatistics.raw.header(),
                         ['Project',
                          'Sample',
                          'Fastq',
                          'Size',
                          'Nreads',
                          'Paired_end',
                          'Read_number',
                          'L1','L2','L3','L4'])
        self.assertEqual(len(fqstatistics.raw),24)
        expected = [
            ['AB','AB1','AB1_S1_L001_R1_001.fastq.gz',3,{'L1':3}],
            ['AB','AB1','AB1_S1_L001_R2_001.fastq.gz',3,{'L1':3}],
            ['AB','AB1','AB1_S1_L002_R1_001.fastq.gz',2,{'L2':2}],
            ['AB','AB1','AB1_S1_L002_R2_001.fastq.gz',2,{'L2':2}],
            ['AB','AB2','AB2_S2_L001_R1_001.fastq.gz',5,{'L1':5}],
            ['AB','AB2','AB2_S2_L001_R2_001.fastq.gz',5,{'L1':5}],
            ['AB','AB2','AB2_S2_L002_R1_001.fastq.gz',3,{'L2':3}],
            ['AB','AB2','AB2_S2_L002_R2_001.fastq.gz',3,{'L2':3}],
            ['CDE','CDE3','CDE3_S3_L003_R1_001.fastq.gz',7,{'L3':7}],
            ['CDE','CDE3','CDE3_S3_L003_R2_001.fastq.gz',7,{'L3':7}],
            ['CDE','CDE3','CDE3_S3_L004_R1_001.fastq.gz',2,{'L4':2}],
            ['CDE','CDE3','CDE3_S3_L004_R2_001.fastq.gz',2,{'L4':2}],
            ['CDE','CDE4','CDE4_S4_L003_R1_001.fastq.gz',8,{'L3':8}],
            ['CDE','CDE4','CDE4_S4_L003_R2_001.fastq.gz',8,{'L3':8}],
            ['CDE','CDE4','CDE4_S4_L004_R1_001.fastq.gz',6,{'L4':6}],
            ['CDE','CDE4','CDE4_S4_L004_R2_001.fastq.gz',6,{'L4':6}],
            ['Undetermined_indices','lane1',
             'Undetermined_S0_L001_R1_001.fastq.gz',2,{'L1':2}],
            ['Undetermined_indices','lane1',
             'Undetermined_S0_L001_R2_001.fastq.gz',2,{'L1':2}],
            ['Undetermined_indices','lane2',
             'Undetermined_S0_L002_R1_001.fastq.gz',1,{'L2':1}],
            ['Undetermined_indices','lane2',
             'Undetermined_S0_L002_R2_001.fastq.gz',1,{'L2':1}],
            ['Undetermined_indices','lane3',
             'Undetermined_S0_L003_R1_001.fastq.gz',4,{'L3':4}],
            ['Undetermined_indices','lane3',
             'Undetermined_S0_L003_R2_001.fastq.gz',4,{'L3':4}],
            ['Undetermined_indices','lane4',
             'Undetermined_S0_L004_R1_001.fastq.gz',3,{'L4':3}],
            ['Undetermined_indices','lane4',
             'Undetermined_S0_L004_R2_001.fastq.gz',3,{'L4':3}],
        ]
        for line,expctd in zip(fqstatistics.raw,expected):
            self.assertEqual(line['Project'],expctd[0])
            self.assertEqual(line['Sample'],expctd[1])
            self.assertEqual(line['Fastq'],expctd[2])
            self.assertEqual(line['Nreads'],expctd[3])
            for lane in ('L1','L2','L3','L4'):
                if lane in expctd[4]:
                    self.assertEqual(line[lane],expctd[4][lane])
                else:
                    self.assertEqual(line[lane],'')
            self.assertEqual(line['Read_number'],
                             IlluminaFastq(expctd[2]).read_number)
    def test_report_per_lane_sample_stats(self):
        fp = cStringIO.StringIO()
        self._setup_bcl2fastq2()
        fqstatistics = FastqStatistics(
            IlluminaData(
                self.illumina_data,
                unaligned_dir="bcl2fastq"))
        fqstatistics.report_per_lane_sample_stats(fp=fp)
        self.assertEqual(fp.getvalue(),"""
Lane 1
Total reads = 10
- AB/AB1	3	30.0%
- AB/AB2	5	50.0%
- Undetermined_indices/lane1	2	20.0%

Lane 2
Total reads = 6
- AB/AB1	2	33.3%
- AB/AB2	3	50.0%
- Undetermined_indices/lane2	1	16.7%

Lane 3
Total reads = 19
- CDE/CDE3	7	36.8%
- CDE/CDE4	8	42.1%
- Undetermined_indices/lane3	4	21.1%

Lane 4
Total reads = 11
- CDE/CDE3	2	18.2%
- CDE/CDE4	6	54.5%
- Undetermined_indices/lane4	3	27.3%
""")
    def test_report_per_lane_summary_stats(self):
        fp = cStringIO.StringIO()
        self._setup_bcl2fastq2()
        fqstatistics = FastqStatistics(
            IlluminaData(
                self.illumina_data,
                unaligned_dir="bcl2fastq"))
        fqstatistics.report_per_lane_summary_stats(fp=fp)
        self.assertEqual(fp.getvalue(),"""#Lane	Total reads	Assigned reads	Unassigned reads	%assigned	%unassigned
Lane 1	10	8	2	80.0	20.0
Lane 2	6	5	1	83.33	16.67
Lane 3	19	15	4	78.95	21.05
Lane 4	11	8	3	72.73	27.27
""")

class TestFastqStatisticsBcl2fastq2NoLaneSplitting(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestFastqStats')
    def tearDown(self):
        # Remove the temporary test directory
        #shutil.rmtree(self.dirn)
        pass
    def _setup_bcl2fastq2_no_lane_splitting(self):
        # Create mock bcl2fastq2 dir structure with no lane splitting
        mock_data = AugmentedMockIlluminaData(
            '151125_S00879_0001_000000000-ABCDE1_analysis',
            'bcl2fastq2',
            unaligned_dir='bcl2fastq',
            paired_end=True,
            no_lane_splitting=True,
            top_dir=self.dirn)
        mock_data.add_fastq_batch('AB','AB1','AB1_S1',lanes=[1,2])
        mock_data.add_fastq_batch('AB','AB2','AB2_S2',lanes=[1,2])
        mock_data.add_fastq_batch('CDE','CDE3','CDE3_S3',lanes=[3,4])
        mock_data.add_fastq_batch('CDE','CDE4','CDE4_S4',lanes=[3,4])
        mock_data.add_undetermined(lanes=(1,2,3,4))
        # Create on disk
        mock_data.create()
        # Populate the FASTQs with 'fake' reads
        reads = {
            "AB": {
                "AB1": { 1:3, 2:2 },
                "AB2": { 1:5, 2:3 },
            },
            "CDE": {
                "CDE3": { 3:7, 4:2 },
                "CDE4": { 3:8, 4:6 },
            },
            "Undetermined_indices": {
                "Undetermined": { 1: 2, 2: 1, 3: 4, 4: 3 },
            },
        }
        s_indices = {
            "AB1": "S1",
            "AB2": "S2",
            "CDE3": "S3",
            "CDE4": "S4"
        }
        for project in reads:
            for sample in reads[project]:
                for lane in reads[project][sample]:
                    for read_number in (1,2):
                        try:
                            s_index = s_indices[sample]
                        except KeyError:
                            s_index = "S0"
                        nreads = reads[project][sample][lane]
                        if project != "Undetermined_indices":
                            project_name = project
                        else:
                            project_name = ""
                        fastq = os.path.join(
                            project_name,
                            "%s_%s_R%d_001.fastq.gz" %
                            (sample,s_index,read_number))
                        mock_data.populate_fastq(fastq,nreads,lane=lane,
                                                 append=True)
        # Store the location of the mock data
        self.illumina_data = mock_data.dirn
    def test_fastqstatistics_bcl2fastq2_no_lane_splitting(self):
        self._setup_bcl2fastq2_no_lane_splitting()
        fqstatistics = FastqStatistics(
            IlluminaData(
                self.illumina_data,
                unaligned_dir="bcl2fastq"))
        self.assertEqual(fqstatistics.lane_names,
                         ['L1','L2','L3','L4'])
        self.assertEqual(fqstatistics.raw.header(),
                         ['Project',
                          'Sample',
                          'Fastq',
                          'Size',
                          'Nreads',
                          'Paired_end',
                          'Read_number',
                          'L1','L2','L3','L4'])
        # Check "raw" stored data
        self.assertEqual(len(fqstatistics.raw),10)
        expected = [
            ['AB','AB1','AB1_S1_R1_001.fastq.gz',5,{'L1':3,'L2':2}],
            ['AB','AB1','AB1_S1_R2_001.fastq.gz',5,{'L1':3,'L2':2}],
            ['AB','AB2','AB2_S2_R1_001.fastq.gz',8,{'L1':5,'L2':3}],
            ['AB','AB2','AB2_S2_R2_001.fastq.gz',8,{'L1':5,'L2':3}],
            ['CDE','CDE3','CDE3_S3_R1_001.fastq.gz',9,{'L3':7,'L4':2}],
            ['CDE','CDE3','CDE3_S3_R2_001.fastq.gz',9,{'L3':7,'L4':2}],
            ['CDE','CDE4','CDE4_S4_R1_001.fastq.gz',14,{'L3':8,'L4':6}],
            ['CDE','CDE4','CDE4_S4_R2_001.fastq.gz',14,{'L3':8,'L4':6}],
            ['Undetermined_indices','undetermined',
             'Undetermined_S0_R1_001.fastq.gz',10,
             {'L1':2,'L2':1,'L3':4,'L4':3}],
            ['Undetermined_indices','undetermined',
             'Undetermined_S0_R2_001.fastq.gz',10,
             {'L1':2,'L2':1,'L3':4,'L4':3}],
        ]
        for line,expctd in zip(fqstatistics.raw,expected):
            self.assertEqual(line['Project'],expctd[0])
            self.assertEqual(line['Sample'],expctd[1])
            self.assertEqual(line['Fastq'],expctd[2])
            self.assertEqual(line['Nreads'],expctd[3])
            for lane in ('L1','L2','L3','L4'):
                if lane in expctd[4]:
                    self.assertEqual(line[lane],expctd[4][lane])
                else:
                    self.assertEqual(line[lane],'')
            self.assertEqual(line['Read_number'],
                             IlluminaFastq(expctd[2]).read_number)
    def test_report_per_lane_sample_stats(self):
        fp = cStringIO.StringIO()
        self._setup_bcl2fastq2_no_lane_splitting()
        fqstatistics = FastqStatistics(
            IlluminaData(
                self.illumina_data,
                unaligned_dir="bcl2fastq"))
        fqstatistics.report_per_lane_sample_stats(fp=fp)
        self.assertEqual(fp.getvalue(),"""
Lane 1
Total reads = 10
- AB/AB1	3	30.0%
- AB/AB2	5	50.0%
- Undetermined_indices/undetermined	2	20.0%

Lane 2
Total reads = 6
- AB/AB1	2	33.3%
- AB/AB2	3	50.0%
- Undetermined_indices/undetermined	1	16.7%

Lane 3
Total reads = 19
- CDE/CDE3	7	36.8%
- CDE/CDE4	8	42.1%
- Undetermined_indices/undetermined	4	21.1%

Lane 4
Total reads = 11
- CDE/CDE3	2	18.2%
- CDE/CDE4	6	54.5%
- Undetermined_indices/undetermined	3	27.3%
""")
    def test_report_per_lane_summary_stats(self):
        fp = cStringIO.StringIO()
        self._setup_bcl2fastq2_no_lane_splitting()
        fqstatistics = FastqStatistics(
            IlluminaData(
                self.illumina_data,
                unaligned_dir="bcl2fastq"))
        fqstatistics.report_per_lane_summary_stats(fp=fp)
        self.assertEqual(fp.getvalue(),"""#Lane	Total reads	Assigned reads	Unassigned reads	%assigned	%unassigned
Lane 1	10	8	2	80.0	20.0
Lane 2	6	5	1	83.33	16.67
Lane 3	19	15	4	78.95	21.05
Lane 4	11	8	3	72.73	27.27
""")

# FastqStats
class TestFastqStats(unittest.TestCase):
    def test_fastqstats_r1(self):
        fqs = FastqStats(
            "/data/fastqs/test_S1_L001_R1_001.fastq",
            "Proj","test")
        self.assertEqual(fqs.fastq,"/data/fastqs/test_S1_L001_R1_001.fastq")
        self.assertEqual(fqs.project,"Proj")
        self.assertEqual(fqs.sample,"test")
        self.assertEqual(fqs.name,"test_S1_L001_R1_001.fastq")
        self.assertEqual(fqs.lanes,[])
        self.assertEqual(fqs.read_number,1)
        self.assertEqual(fqs.nreads,None)
        self.assertEqual(fqs.fsize,None)
        self.assertEqual(fqs.reads_by_lane,{})

    def test_fastqstats_r2(self):
        fqs = FastqStats(
            "/data/fastqs/test_S1_L001_R2_001.fastq",
            "Proj","test")
        self.assertEqual(fqs.fastq,"/data/fastqs/test_S1_L001_R2_001.fastq")
        self.assertEqual(fqs.project,"Proj")
        self.assertEqual(fqs.sample,"test")
        self.assertEqual(fqs.name,"test_S1_L001_R2_001.fastq")
        self.assertEqual(fqs.lanes,[])
        self.assertEqual(fqs.read_number,2)
        self.assertEqual(fqs.nreads,None)
        self.assertEqual(fqs.fsize,None)
        self.assertEqual(fqs.reads_by_lane,{})

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

# collect_fastq_data
class TestCollectFastqData(unittest.TestCase):
    def setUp(self):
        # Temporary working dir (if needed)
        self.wd = None

    def _make_working_dir(self):
        if self.wd is None:
            self.wd = tempfile.mkdtemp(suffix='.test_collect_fastq_data')

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

    def test_collect_fastq_data_r1_lane_in_name(self):
        fastq = self._make_fastq("test_S1_L001_R1_001.fastq",
                                 fastq_r1_data)
        fqs = FastqStats(fastq,"Proj","test")
        self.assertEqual(fqs.fastq,fastq)
        collect_fastq_data(fqs)
        self.assertEqual(fqs.project,"Proj")
        self.assertEqual(fqs.sample,"test")
        self.assertEqual(fqs.name,"test_S1_L001_R1_001.fastq")
        self.assertEqual(fqs.lanes,[1])
        self.assertEqual(fqs.read_number,1)
        self.assertEqual(fqs.nreads,3)
        self.assertEqual(fqs.fsize,os.path.getsize(fastq))
        self.assertEqual(fqs.reads_by_lane,{1: 3})

    def test_collect_fastq_data_r1_no_lane_in_name(self):
        fastq = self._make_fastq("test_S1_R1_001.fastq",
                                 fastq_r1_data)
        fqs = FastqStats(fastq,"Proj","test")
        self.assertEqual(fqs.fastq,fastq)
        collect_fastq_data(fqs)
        self.assertEqual(fqs.project,"Proj")
        self.assertEqual(fqs.sample,"test")
        self.assertEqual(fqs.name,"test_S1_R1_001.fastq")
        self.assertEqual(fqs.lanes,[1])
        self.assertEqual(fqs.read_number,1)
        self.assertEqual(fqs.nreads,3)
        self.assertEqual(fqs.fsize,os.path.getsize(fastq))
        self.assertEqual(fqs.reads_by_lane,{1: 3})

    def test_collect_fastq_data_r2_lane_in_name(self):
        fastq = self._make_fastq("test_S1_L001_R2_001.fastq",
                                 fastq_r2_data)
        fqs = FastqStats(fastq,"Proj","test")
        self.assertEqual(fqs.fastq,fastq)
        collect_fastq_data(fqs)
        self.assertEqual(fqs.project,"Proj")
        self.assertEqual(fqs.sample,"test")
        self.assertEqual(fqs.name,"test_S1_L001_R2_001.fastq")
        self.assertEqual(fqs.lanes,[])
        self.assertEqual(fqs.read_number,2)
        self.assertEqual(fqs.nreads,3)
        self.assertEqual(fqs.fsize,os.path.getsize(fastq))
        self.assertEqual(fqs.reads_by_lane,{})

    def test_collect_fastq_data_r2_no_lane_in_name(self):
        fastq = self._make_fastq("test_S1_R2_001.fastq",
                                 fastq_r2_data)
        fqs = FastqStats(fastq,"Proj","test")
        self.assertEqual(fqs.fastq,fastq)
        collect_fastq_data(fqs)
        self.assertEqual(fqs.project,"Proj")
        self.assertEqual(fqs.sample,"test")
        self.assertEqual(fqs.name,"test_S1_R2_001.fastq")
        self.assertEqual(fqs.lanes,[])
        self.assertEqual(fqs.read_number,2)
        self.assertEqual(fqs.nreads,3)
        self.assertEqual(fqs.fsize,os.path.getsize(fastq))
        self.assertEqual(fqs.reads_by_lane,{})

    def test_collect_fastq_data_multi_lanes(self):
        fastq = self._make_fastq("test_S1_R1_001.fastq",
                                 fastq_multi_lane_data)
        fqs = FastqStats(fastq,"Proj","test")
        self.assertEqual(fqs.fastq,fastq)
        collect_fastq_data(fqs)
        self.assertEqual(fqs.project,"Proj")
        self.assertEqual(fqs.sample,"test")
        self.assertEqual(fqs.name,"test_S1_R1_001.fastq")
        self.assertEqual(fqs.lanes,[1,2,3,4])
        self.assertEqual(fqs.read_number,1)
        self.assertEqual(fqs.nreads,12)
        self.assertEqual(fqs.fsize,os.path.getsize(fastq))
        self.assertEqual(fqs.reads_by_lane,{ 1: 5,
                                             2: 1,
                                             3: 2,
                                             4: 4 })
