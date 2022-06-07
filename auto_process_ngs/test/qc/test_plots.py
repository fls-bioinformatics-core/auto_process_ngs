#######################################################################
# Unit tests for qc/plots.py
#######################################################################

import unittest
import shutil
import os
import base64
import tempfile

from auto_process_ngs.qc.plots import Plot
from auto_process_ngs.qc.plots import encode_png
from auto_process_ngs.qc.plots import uscreenplot
from auto_process_ngs.qc.plots import uboxplot
from auto_process_ngs.qc.plots import ufastqcplot
from auto_process_ngs.qc.plots import ustackedbar
from auto_process_ngs.qc.plots import ustrandplot

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestPlot(unittest.TestCase):
    """
    Tests for the Plot class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestPlot')
        # Define colours
        self.rgb_blue = (0,0,255)
        self.rgb_green = (0,128,0)
        self.rgb_grey = (145,145,145)
        self.rgb_red = (255,0,0)
        self.rgb_yellow = (255,255,0)

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_plot_make_empty_plot(self):
        """
        Plot: make empty plot
        """
        outfile = os.path.join(self.wd,"empty_plot.png")
        empty_plot_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAKCAIAAAA7N+mxAAAAF0lEQVR4nGP8//8/A7mAiWydo5pHjGYAM38DEWQYrPYAAAAASUVORK5CYII="
        p = Plot(20,10)
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        self.assertEqual(encode_png(outfile),
                         empty_plot_base64_data)

    def test_plot_make_empty_plot_with_bbox(self):
        """
        Plot: make empty plot with bounding box
        """
        outfile = os.path.join(self.wd,"empty_plot_bbox.png")
        empty_plot_bbox_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAKCAIAAAA7N+mxAAAAL0lEQVR4nO3NsREAMAzCQCXnedkFJvYOFKmi/k/HNm0DSCpkkltvgY9f4gGSdHgBjBwGgSEOxLEAAAAASUVORK5CYII="
        p = Plot(20,10)
        p.bbox(self.rgb_grey)
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        self.assertEqual(encode_png(outfile),
                         empty_plot_bbox_base64_data)

    def test_plot_make_empty_plot_with_striping(self):
        """
        Plot: make empty plot with striping
        """
        outfile = os.path.join(self.wd,"empty_plot_striping.png")
        empty_plot_striping_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAKCAIAAAA7N+mxAAAAHElEQVR4nGP8z8DA8J+BgZEckomBAjCqeWRoBgAMCgsSRyKjXwAAAABJRU5ErkJggg=="
        p = Plot(20,10)
        p.stripe(self.rgb_red,self.rgb_yellow)
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        self.assertEqual(encode_png(outfile),
                         empty_plot_striping_base64_data)

    def test_plot_hline(self):
        """
        Plot: plot horizontal line
        """
        outfile = os.path.join(self.wd,"hline.png")
        hline_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAKCAIAAAA7N+mxAAAAIklEQVR4nGP8//8/A7mAiWydA6qZkXwfMzCwMIy8AKNIMwD61gYPdVB52wAAAABJRU5ErkJggg=="
        p = Plot(20,10)
        p.hline(5,self.rgb_red)
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        self.assertEqual(encode_png(outfile),
                         hline_base64_data)

    def test_plot_vline(self):
        """
        Plot: plot vertical line
        """
        outfile = os.path.join(self.wd,"vline.png")
        vline_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAKCAIAAAA7N+mxAAAAHklEQVR4nGP8//8/Ax7AyMiAWwETPp2EwKjmkaEZALDBBRHE0mGAAAAAAElFTkSuQmCC"
        p = Plot(20,10)
        p.vline(10,self.rgb_red)
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        self.assertEqual(encode_png(outfile),
                         vline_base64_data)

    def test_block(self):
        """
        Plot: draw solid block
        """
        outfile = os.path.join(self.wd,"block.png")
        block_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAKCAIAAAA7N+mxAAAAK0lEQVR4nGP8//8/A7mAiWydDAwMLFCakZE0ff//U2rzSNQMC22yYpsimwF5/AcTuAyJrwAAAABJRU5ErkJggg=="
        p = Plot(20,10)
        p.block((2,2),(18,8),self.rgb_red)
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        self.assertEqual(encode_png(outfile),block_base64_data)

    def test_block_with_striping(self):
        """
        Plot: draw block with striping
        """
        outfile = os.path.join(self.wd,"block_with_striping.png")
        block_with_striping_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAKCAIAAAA7N+mxAAAALUlEQVR4nGP8//8/A7mAiWydDAwMLFCakZHhPwMDIwOxJMN/Sm0eiZoZByyeAY20Dw9p9T90AAAAAElFTkSuQmCC"
        p = Plot(20,10)
        p.block((2,2),(18,8),self.rgb_red,self.rgb_yellow)
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        print(encode_png(outfile))
        self.assertEqual(encode_png(outfile),
                         block_with_striping_base64_data)

    def test_plot(self):
        """
        Plot: plot dataset
        """
        outfile = os.path.join(self.wd,"plot.png")
        data = { 0: 23,
                 1: 33,
                 2: 30,
                 3: 53,
                 4: 42,
                 5: 48,
                 6: 52,
                 7: 60,
                 8: 74,
                 9: 67,
                 10: 64,
                 11: 87,
                 12: 88,
                 13: 73,
                 14: 99,
                 15: 113,
                 16: 145,
                 17: 152,
                 18: 170,
                 19: 212,
        }
        plot_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAKCAIAAAA7N+mxAAAATElEQVR4nJ3RQQrAIBAEweng/788HjYERSPJ9FWLFRfbyoIrlJLsFIOkljHZ/zFo+KMPz4Z72izfJ9ftqsAiZ7yCsd1G8fH4XAvMUwfKhhgZZ/aOfQAAAABJRU5ErkJggg=="
        p = Plot(20,10)
        p.plot(data,self.rgb_red)
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        print(encode_png(outfile))
        self.assertEqual(encode_png(outfile),plot_base64_data)

    def test_plot_with_fill(self):
        """
        Plot: plot dataset (fill under curve)
        """
        outfile = os.path.join(self.wd,"plot.png")
        data = { 0: 23,
                 1: 33,
                 2: 30,
                 3: 53,
                 4: 42,
                 5: 48,
                 6: 52,
                 7: 60,
                 8: 74,
                 9: 67,
                 10: 64,
                 11: 87,
                 12: 88,
                 13: 73,
                 14: 99,
                 15: 113,
                 16: 145,
                 17: 152,
                 18: 170,
                 19: 212,
        }
        plot_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAIAAAACUFjqAAAANUlEQVR4nGP8//8/A27AhEeOgZERrzQ+3YyMhAzHKc3ICKFZcElgSKNKIEljk4BqwRcoDAwAzDQEIYRsFcgAAAAASUVORK5CYII="
        p = Plot(10,10)
        p.plot(data,self.rgb_red,fill=True)
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        print(encode_png(outfile))
        self.assertEqual(encode_png(outfile),plot_base64_data)

    def test_plot_range(self):
        """
        Plot: plot range defined by pair of datasets
        """
        outfile = os.path.join(self.wd,"plot_range.png")
        data1 = {
            0: 23,
            1: 33,
            2: 30,
            3: 53,
            4: 42,
            5: 48,
            6: 52,
            7: 60,
            8: 74,
            9: 67
        }
        data2 = {
            0: 64,
            1: 87,
            2: 88,
            3: 73,
            4: 99,
            5: 113,
            6: 145,
            7: 152,
            8: 170,
            9: 212
        }
        plot_range_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAIAAAACUFjqAAAAQElEQVR4nI2PQQqAQBDDkmX//+V6EEfUWdkeEwqtSWijwOjdlYVWgGS8kZYDZlsqMG965rnUFFI+L1we+1u+qQ+PZRMTiub8GAAAAABJRU5ErkJggg=="
        p = Plot(10,10)
        p.plot_range(data1,data2,self.rgb_red)
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        print(encode_png(outfile))
        self.assertEqual(encode_png(outfile),
                         plot_range_base64_data)

    def test_stackedbar(self):
        """
        Plot: make stacked bar
        """
        outfile = os.path.join(self.wd,"stacked_bar.png")
        data = [ 30.0, 20.0, 15.0, 10.0, 100.0 ]
        stacked_bar_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAKCAIAAAA7N+mxAAAAO0lEQVR4nGP8//8/A7mAiWydDAwMLFCakZGBgYHhPwMDAwNjIwNDw/+JEyfh0ZaXl0epzSNRM+OAxTMAXhYMDwbgrT8AAAAASUVORK5CYII="
        p = Plot(20,10)
        p.bar(data,
              (2,2),(18,8),
              (self.rgb_red,
               self.rgb_yellow,
               self.rgb_green,
               self.rgb_blue,
               self.rgb_grey))
        p.save(outfile)
        self.assertTrue(os.path.exists(outfile))
        print(encode_png(outfile))
        self.assertEqual(encode_png(outfile),
                         stacked_bar_base64_data)

    def test_plot_rebin_data(self):
        """
        Plot: rebin data along x-axis
        """
        data = { 0: 23,
                 1: 33,
                 2: 30,
                 3: 53,
                 4: 42,
                 5: 48,
                 6: 52,
                 7: 60,
                 8: 74,
                 9: 67,
                 10: 64,
                 11: 87,
                 12: 88,
                 13: 73,
                 14: 99,
                 15: 113,
                 16: 145,
                 17: 152,
                 18: 170,
                 19: 212,
        }
        p = Plot(10,10)
        rebinned_data = p.rebin_data(data)
        self.assertEqual(rebinned_data.mean,
                         {
                             0: 28.0,
                             1: 41.5,
                             2: 45.0,
                             3: 56.0,
                             4: 70.5,
                             5: 75.5,
                             6: 80.5,
                             7: 106.0,
                             8: 148.5,
                             9: 191.0
                         })
        self.assertEqual(rebinned_data.min,
                         {
                             0: 23,
                             1: 30,
                             2: 42,
                             3: 52,
                             4: 67,
                             5: 64,
                             6: 73,
                             7: 99,
                             8: 145,
                             9: 170
                         })
        self.assertEqual(rebinned_data.max,
                         {
                             0: 33,
                             1: 53,
                             2: 48,
                             3: 60,
                             4: 74,
                             5: 87,
                             6: 88,
                             7: 113,
                             8: 152,
                             9: 212
                         })

    def test_plot_normalise_data(self):
        """
        Plot: normalise data along y-axis
        """
        data = {
            0: 25,
            1: 35,
            2: 30,
            3: 55,
            4: 40,
            5: 50,
            6: 50,
            7: 60,
            8: 75,
            9: 60
        }
        p = Plot(10,10)
        self.assertEqual(p.normalise_data(data),
                         {
                             0: 3.0,
                             1: 4.2,
                             2: 3.6,
                             3: 6.6,
                             4: 4.8,
                             5: 6.0,
                             6: 6.0,
                             7: 7.2,
                             8: 9.0,
                             9: 7.2
                         })

class TestEncodePng(unittest.TestCase):
    """
    Tests for the encode_png function
    """
    def setUp(self):
        # Make a test PNG file
        self.png_base64_data = "iVBORw0KGgoAAAANSUhEUgAAAAQAAAAECAIAAAAmkwkpAAAAIElEQVR4nGNkYPjPwMDw/z8DAwMDEwMSYIGIMTJiyAAAu6IFB//RUJMAAAAASUVORK5CYII="
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            self.png = fp.name
            fp.write(base64.decodebytes(self.png_base64_data.encode()))

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
        self.png_base64_data = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAqCAIAAAAF/F3cAAAAnUlEQVR4nO3UsQ0CMQwF0O8TS0Bxa7BGKFgjnuVnDYpkHOor7momMMVFIheMaEHg4kfRcxSniZBEU1G13YKkmeU5PyfJHbYlAiAANXs2Q1lK2IeylOkyVQ6HU+03a7uHdUkkzAQbe7BqFEGei89ksv5kw69qiKrra3xOpMDG8+izanRvrUwmAMcb3Hw32o/z1c9PGO3P38LSfXtd3QFX10hamZN+QgAAAABJRU5ErkJggg=="

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
