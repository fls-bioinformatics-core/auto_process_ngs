#######################################################################
# Unit tests for qc/plots.py
#######################################################################

import unittest
import shutil
import os
import base64
import tempfile

from auto_process_ngs.qc.plots import encode_png
from auto_process_ngs.qc.plots import uscreenplot
from auto_process_ngs.qc.plots import uboxplot
from auto_process_ngs.qc.plots import ufastqcplot
from auto_process_ngs.qc.plots import ustackedbar
from auto_process_ngs.qc.plots import ustrandplot

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestEncodePng(unittest.TestCase):
    """
    Tests for the encode_png function
    """
    def setUp(self):
        # Make a test PNG file
        self.png_base64_data = "iVBORw0KGgoAAAANSUhEUgAAAAQAAAAECAIAAAAmkwkpAAAAIElEQVR4nGNkYPjPwMDw/z8DAwMDEwMSYIGIMTJiyAAAu6IFB//RUJMAAAAASUVORK5CYII="
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            self.png = fp.name
            fp.write(base64.decodestring(self.png_base64_data))

    def tearDown(self):
        try:
            os.remove(self.png)
        except Exception:
            pass

    def test_enode_png(self):
        """
        encode_png: turn PNG file into Base64 encoded string
        """
        self.assertEqual(encode_png(self.png),
                         "data:image/png;base64,%s" %
                         self.png_base64_data)

class TestUScreenPlot(unittest.TestCase):
    """
    Tests for the uscreenplot function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestUBoxplot')
        screen1_data = """#Fastq_screen version: 0.4.2	#Reads in subset: 1000000
Library	#Reads_processed	#Unmapped	%Unmapped	#One_hit_one_library	%One_hit_one_library	#Multiple_hits_one_library	%Multiple_hits_one_library	#One_hit_multiple_libraries	%One_hit_multiple_libraries	Multiple_hits_multiple_libraries	%Multiple_hits_multiple_libraries
hg19	99859	89708	89.83	42	0.04	117	0.12	780	0.78	9212	9.23
mm9	99859	18069	18.09	50157	50.23	13291	13.31	7087	7.10	11255	11.27
rn4	99859	81116	81.23	134	0.13	161	0.16	5678	5.69	12770	12.79

%Hit_no_libraries: 99.73
"""
        screen2_data = """#Fastq_screen version: 0.4.2	#Reads in subset: 1000000
Library	#Reads_processed	#Unmapped	%Unmapped	#One_hit_one_library	%One_hit_one_library	#Multiple_hits_one_library	%Multiple_hits_one_library	#One_hit_multiple_libraries	%One_hit_multiple_libraries	Multiple_hits_multiple_libraries	%Multiple_hits_multiple_libraries
dm3	99859	99836	99.98	0	0.00	0	0.00	1	0.00	22	0.02
ws200	99859	99840	99.98	0	0.00	0	0.00	4	0.00	15	0.02
ecoli	99859	99836	99.98	4	0.00	18	0.02	0	0.00	1	0.00

%Hit_no_libraries: 16.25
"""
        self.screens = list()
        for screen_data,screen_file in zip((screen1_data,
                                            screen2_data),
                                           ("model_organisms_screen.txt",
                                            "other_organisms_screen.txt")):
            screen_file = os.path.join(self.wd,screen_file)
            with open(screen_file,'w') as fp:
                fp.write(screen_data)
            self.screens.append(screen_file)
        # Reference data
        self.png_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGQAAAAUCAIAAAD0og/CAAAAtElEQVR4nO2YwQrEIAxEp4vf2t7iJ5R+QnL0Z91DqctiWXagJRXyTmpAhoGMkklVEfzHBFSgVCzbd2Gt1UcRYGYARMRLwClm9vLWMBJhFkGYRZBUTUSAefWW8nxSzpLzhReW/qhiaeutqzq+JCzRhgRhFkGYRdAC/irm34dDPyN7wJeWwX0AY6gMvpVoQ4IwiyDMIvj84Pf90AF8NwnHSORpnKryndtEGxKEWQTJWwCHb2K8AUVlLWZTGFiYAAAAAElFTkSuQmCC"
        
    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_uscreenplot_to_file(self):
        """uscreenplot: write PNG to file
        """
        outfile = os.path.join(self.wd,"uscreenplot.png")
        self.assertEqual(uscreenplot(self.screens,
                                     outfile=outfile),
                         outfile)
        self.assertEqual(encode_png(outfile),self.png_base64_data)

    def test_uscreenplot_to_base64(self):
        """uscreenplot: write PNG as Base64 encoded string
        """
        self.assertEqual(uscreenplot(self.screens,
                                     inline=True),
                         self.png_base64_data)

class TestUBoxplot(unittest.TestCase):
    """
    Tests for the uboxplot function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestUBoxplot')
        fastqc_data = """##FastQC	0.11.3
>>Basic Statistics	pass
#Measure	Value
Filename	GS10_S92_R1_001.fastq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	28559451
Sequences flagged as poor quality	0
Sequence length	35-76
%GC	48
>>END_MODULE
>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile
1	30.97962481841825	32.0	32.0	32.0	32.0	32.0
2	31.260586276675976	32.0	32.0	32.0	32.0	32.0
3	35.571439450989445	37.0	37.0	37.0	32.0	37.0
4	35.96997564834142	37.0	37.0	37.0	37.0	37.0
5	36.11201426806139	37.0	37.0	37.0	37.0	37.0
6	39.82052498138007	41.0	41.0	41.0	41.0	41.0
7	39.75080921548527	41.0	41.0	41.0	41.0	41.0
8	39.8729139086042	41.0	41.0	41.0	41.0	41.0
9	39.88060719374472	41.0	41.0	41.0	41.0	41.0
10	39.880201513677555	41.0	41.0	41.0	41.0	41.0
>>END_MODULE
"""
        self.fastqc_data_txt = os.path.join(self.wd,"fastqc_data.txt")
        with open(self.fastqc_data_txt,'w') as fp:
            fp.write(fastqc_data)
        # Reference data
        self.png_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAqCAIAAAAF/F3cAAAAnUlEQVR4nO3UsQ3CMBAF0H8RS0CRNVjDFKzhm+V7DYp4HOoUSc0EnyKRMMZAiyKu+Jb1zvK5sZFEUdG93IKkpGEaXpPkDs9lBiAAa9YsIc857EOe83gZVw6H09ovld3dsiQSkuHJHuwezTBMuc1kUn2y4HfVRfflNW1OpEH9uW+ze5RwvAFoZEemT3d/GW3rfG3nL4z2542wVd9eVXdP5kha/kjTTgAAAABJRU5ErkJggg=="

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_uboxplot_to_file(self):
        """uboxplot: write PNG to file
        """
        outfile = os.path.join(self.wd,"uboxplot.png")
        self.assertEqual(uboxplot(self.fastqc_data_txt,
                                  outfile=outfile),
                         outfile)
        self.assertEqual(encode_png(outfile),self.png_base64_data)

    def test_uboxplot_to_base64(self):
        """uboxplot: write PNG as Base64 encoded string
        """
        self.assertEqual(uboxplot(self.fastqc_data_txt,
                                  inline=True),
                         self.png_base64_data)

class TestUFastqcPlot(unittest.TestCase):
    """
    Tests for the ufastqcplot function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestUStackedBar')
        fastqc_summary_data = """PASS	Basic Statistics	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Per base sequence quality	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Per tile sequence quality	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Per sequence quality scores	ES1_GTCCGC_L008_R1_001.fastq.gz
FAIL	Per base sequence content	ES1_GTCCGC_L008_R1_001.fastq.gz
WARN	Per sequence GC content	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Per base N content	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Sequence Length Distribution	ES1_GTCCGC_L008_R1_001.fastq.gz
FAIL	Sequence Duplication Levels	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Overrepresented sequences	ES1_GTCCGC_L008_R1_001.fastq.gz
PASS	Adapter Content	ES1_GTCCGC_L008_R1_001.fastq.gz
FAIL	Kmer Content	ES1_GTCCGC_L008_R1_001.fastq.gz
"""
        self.fastqc_summary_txt = os.path.join(self.wd,"summary.txt")
        with open(self.fastqc_summary_txt,'w') as fp:
            fp.write(fastqc_summary_data)
        # Reference data
        self.png_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB4AAAAwCAIAAACJ9F2zAAAAcUlEQVR4nO2UUQrAIAxD27FT7jJjl9k14+8QZUJppJD8ik94kTgAy8lpZv748Ax36NUjcnkb2uWahpZrHlqueehk16H4uCcDaroOC5lXtYB+Jzavn/5rutaG8NByzUNrr7vU3Wv960+0ITy0NqRLIroBHQJCUS0i89gAAAAASUVORK5CYII="
        
    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_ufastqcplot_to_file(self):
        """ufastqcplot: write PNG to file
        """
        outfile = os.path.join(self.wd,"ufastqcplot.png")
        self.assertEqual(ufastqcplot(self.fastqc_summary_txt,
                                     outfile=outfile),
                         outfile)
        self.assertEqual(encode_png(outfile),self.png_base64_data)

    def test_ufastqcplot_to_base64(self):
        """ufastqcplot: write PNG as Base64 encoded string
        """
        self.assertEqual(ufastqcplot(self.fastqc_summary_txt,
                                     inline=True),
                         self.png_base64_data)

class TestUStackedBar(unittest.TestCase):
    """
    Tests for the ustackedbar function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestUStackedBar')
        # Reference data
        self.png_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGQAAAAUCAIAAAD0og/CAAAAXklEQVR4nO3QQRGAQAwEwQ0S0Hgv/PDCIxY4CzcFeTEtYJOaRMsqSfK03hjn3bq/4jr2lwtVtX3yyk8YCzAWYCzAWICxAGMBxgKMBRgLMBZgLMBYgLEAYwHGAoylHhNXSAQk743TiwAAAABJRU5ErkJggg=="

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_ustackedbar_to_file(self):
        """ustackedbar: write PNG to file
        """
        outfile = os.path.join(self.wd,"ustackedbar.png")
        self.assertEqual(ustackedbar((3,4,2),
                                     outfile=outfile,
                                     colors=((0,0,255),
                                             (100,149,237),
                                             (255,255,255))),
                         outfile)
        self.assertEqual(encode_png(outfile),self.png_base64_data)

    def test_ustackedbar_to_base64(self):
        """ustackedbar: write PNG as Base64 encoded string
        """
        self.assertEqual(ustackedbar((3,4,2),
                                     inline=True,
                                     colors=((0,0,255),
                                             (100,149,237),
                                             (255,255,255))),
                         self.png_base64_data)

class TestUStrandPlot(unittest.TestCase):
    """
    Tests for the ustrandplot function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestUStrandPlot')
        fastq_strand_data = """#fastq_strand version: 0.0.1	#Aligner: STAR	#Reads in subset: 10000
#Genome	1st forward	2nd reverse
human	73.23	86.13
mouse	2.82	100.09
"""
        self.fastq_strand_txt = os.path.join(self.wd,"fastq_strand.txt")
        with open(self.fastq_strand_txt,'w') as fp:
            fp.write(fastq_strand_data)
        # Reference data
        self.png_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADIAAAAZCAIAAAD8NuoTAAAAaElEQVR4nO2VoREAMQgEIfPl0H8jFMSLiBdB5EQ+J25NDDNZwQ5eVcbHuC3QQ6r1zMfdT/+UmTtjEWFmrt0CkBYCqdZ/Ja60bapEHGkhkGp9JW7erNOoRBxpIZBq3byJLTNBlYhAqvUCT7kbJ8xKcU8AAAAASUVORK5CYII="
        self.png_base64_data_dynamic = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADIAAAAyCAIAAACRXR/mAAAAh0lEQVR4nO3XMQ6AIBAFUdZ4HO5/EQ6EBYWFFEyBbDGvsTFxQkx+iN57yec6HTCXNOsej4jY/aXW2sprtdaS9rTMIswizCLC8QHMIv7bxK/pSrqJnFmEWUTSLDeRMItImvVO9eL1cjenmjOLMItImuVUE2YRJ6+vU+NfT3paZhFmEWYRbiLxAKbiG1mx7IKbAAAAAElFTkSuQmCC"

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_ustrandplot_to_file(self):
        """ustrandplot: write PNG to file
        """
        outfile = os.path.join(self.wd,"ustrandplot.png")
        self.assertEqual(ustrandplot(self.fastq_strand_txt,
                                     outfile=outfile),
                         outfile)
        self.assertEqual(encode_png(outfile),self.png_base64_data)

    def test_ustrandplot_to_base64(self):
        """ustrandplot: write PNG as Base64 encoded string
        """
        self.assertEqual(ustrandplot(self.fastq_strand_txt,
                                     inline=True),
                         self.png_base64_data)

    def test_ustrandplot_to_file_dynamic(self):
        """ustrandplot: write PNG to file (dynamic mode)
        """
        outfile = os.path.join(self.wd,"ustrandplot.png")
        self.assertEqual(ustrandplot(self.fastq_strand_txt,
                                     outfile=outfile,
                                     dynamic=True),
                         outfile)
        self.assertEqual(encode_png(outfile),
                         self.png_base64_data_dynamic)

    def test_ustrandplot_to_base64_dynamic(self):
        """ustrandplot: write PNG as Base64 encoded string (dynamic mode)
        """
        self.assertEqual(ustrandplot(self.fastq_strand_txt,
                                     inline=True,
                                     dynamic=True),
                         self.png_base64_data_dynamic)
