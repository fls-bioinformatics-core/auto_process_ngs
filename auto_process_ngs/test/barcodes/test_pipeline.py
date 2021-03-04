#######################################################################
# Tests for barcodes/pipeline.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
import gzip
from bcftbx.mock import MockIlluminaData
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaFastq
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.barcodes.pipeline import AnalyseBarcodes

# Set to False to preserve test outputs
REMOVE_TEST_OUTPUTS = True

class TestAnalyseBarcodes(unittest.TestCase):
    """
    Tests for the AnalyseBarcodes pipeline class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestAnalyseBarcodesPipeline')
        # Store original location
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.wd)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def _insert_fastq_reads(self,dirn):
        # Add read data to fastqs
        illuminadata = IlluminaData(dirn,unaligned_dir="bcl2fastq")
        for p in illuminadata.projects:
            for s in p.samples:
                for fq in s.fastq:
                    fastq = os.path.join(s.dirn,fq)
                    lane = IlluminaFastq(fq).lane_number
                    self.assertTrue(os.path.exists(fastq),
                                    "Missing %s" % fastq)
                    with gzip.open(fastq,'wb') as fp:
                        data = """@ILLUMINA-545855:49:FC61RLR:%d:1:10979:1695 1:N:0:TCCTGA
GCATACTCAGCTTTAGTAATAAGTGTGATTCTGGTA
+
IIIIIHIIIGHHIIDGHIIIIIIHIIIIIIIIIIIH
                            """ % lane
                        fp.write(data.encode())

    def test_analyse_barcodes_with_bcl2fastq_dir_no_samplesheet(self):
        """
        AnalyseBarcodes: bcl2fastq directory as input (no samplesheet)
        """
        # Make a mock bcl2fastq output directory
        datadir = MockIlluminaData(
            os.path.join(self.wd,
                         "200428_M00879_0087_000000000-AGEW9"),
            "bcl2fastq2",
            unaligned_dir="bcl2fastq",
            paired_end=True)
        datadir.add_fastq_batch("AB","AB1","AB1_S1")
        datadir.add_fastq_batch("AB","AB2","AB2_S2")
        datadir.add_fastq_batch("CDE","CDE3","CDE3_S3")
        datadir.add_fastq_batch("CDE","CDE4","CDE4_S4")
        datadir.add_fastq_batch("","Undetermined","Undetermined_S0")
        datadir.create()
        # Add data to Fastq files
        self._insert_fastq_reads(
            os.path.join(self.wd,"200428_M00879_0087_000000000-AGEW9"))
        # Set up and run pipeline
        p = AnalyseBarcodes(bcl2fastq_dir=
                            os.path.join(self.wd,
                                         "200428_M00879_0087_000000000-AGEW9",
                                         "bcl2fastq"))
        exit_code = p.run(os.path.join(self.wd,"barcode_analysis"),
                          working_dir=self.wd,
                          poll_interval=0.5)
        # Check outputs
        self.assertEqual(exit_code,0)
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis")),
                        "Missing dir: barcode_analysis")
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis",
                                                   "counts")),
                            "Missing dir: barcode_analysis/counts")
        for f in (
                "AB.AB1_S1_L001_R1_001.fastq.gz.counts",
                "AB.AB2_S2_L001_R1_001.fastq.gz.counts",
                "CDE.CDE3_S3_L001_R1_001.fastq.gz.counts",
                "CDE.CDE4_S4_L001_R1_001.fastq.gz.counts",
                "__undetermined__.Undetermined_S0_L001_R1_001.fastq.gz.counts"):
            self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                        "barcode_analysis",
                                                        "counts",f)),
                            "Missing file: %s" % f)
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.report")),
                        "Missing file: barcodes.report")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.xls")),
                        "Missing file: barcodes.xls")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.html")),
                        "Missing file: barcodes.html")
        # Check that the report content is non-trivial
        barcodes_report = os.path.join(self.wd,
                                       "barcode_analysis",
                                       "barcodes.report")
        with open(barcodes_report,'rt') as fp:
            contents = fp.read()
            self.assertTrue("Barcode analysis for lane #1" in contents)
            self.assertTrue("#Rank\tIndex\tSample\tN_seqs\tN_reads\t%reads\t(%Total_reads)" in contents)
            # Expect 12 lines of content in total
            self.assertEqual(contents.count('\n'),12)

    def test_analyse_barcodes_with_bcl2fastq_dir_and_samplesheet(self):
        """
        AnalyseBarcodes: bcl2fastq directory as input (with samplesheet)
        """
        # Make a mock bcl2fastq output directory
        datadir = MockIlluminaData(
            os.path.join(self.wd,
                         "200428_M00879_0087_000000000-AGEW9"),
            "bcl2fastq2",
            unaligned_dir="bcl2fastq",
            paired_end=True)
        datadir.add_fastq_batch("AB","AB1","AB1_S1")
        datadir.add_fastq_batch("AB","AB2","AB2_S2")
        datadir.add_fastq_batch("CDE","CDE3","CDE3_S3")
        datadir.add_fastq_batch("CDE","CDE4","CDE4_S4")
        datadir.add_fastq_batch("","Undetermined","Undetermined_S0")
        datadir.create()
        # Add data to Fastq files
        self._insert_fastq_reads(
            os.path.join(self.wd,"200428_M00879_0087_000000000-AGEW9"))
        # Create sample sheet
        sample_sheet = os.path.join(self.wd,"custom_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
AB1,AB1,,,D701,CGTGTAGG,D501,GACCTGAA,AB,
AB2,AB2,,,D702,CGTGTAGG,D501,ATGTAACT,AB,
CDE3,CDE3,,,D701,GACCTGAA,D501,CGTGTAGG,CDE,
CDE4,CDE4,,,D702,ATGTAACT,D501,CGTGTAGG,CDE,
""")
        # Set up and run pipeline
        p = AnalyseBarcodes(bcl2fastq_dir=
                            os.path.join(self.wd,
                                         "200428_M00879_0087_000000000-AGEW9",
                                         "bcl2fastq"))
        exit_code = p.run(os.path.join(self.wd,"barcode_analysis"),
                          sample_sheet=sample_sheet,
                          working_dir=self.wd,
                          poll_interval=0.5)
        # Check outputs
        self.assertEqual(exit_code,0)
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis")),
                        "Missing dir: barcode_analysis")
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis",
                                                   "counts")),
                            "Missing dir: barcode_analysis/counts")
        for f in (
                "AB.AB1_S1_L001_R1_001.fastq.gz.counts",
                "AB.AB2_S2_L001_R1_001.fastq.gz.counts",
                "CDE.CDE3_S3_L001_R1_001.fastq.gz.counts",
                "CDE.CDE4_S4_L001_R1_001.fastq.gz.counts",
                "__undetermined__.Undetermined_S0_L001_R1_001.fastq.gz.counts"):
            self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                        "barcode_analysis",
                                                        "counts",f)),
                            "Missing file: %s" % f)
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.report")),
                        "Missing file: barcodes.report")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.xls")),
                        "Missing file: barcodes.xls")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.html")),
                        "Missing file: barcodes.html")
        # Check that the report content is non-trivial
        barcodes_report = os.path.join(self.wd,
                                       "barcode_analysis",
                                       "barcodes.report")
        with open(barcodes_report,'rt') as fp:
            contents = fp.read()
            self.assertTrue("Barcode analysis for lane #1" in contents)
            self.assertTrue("#Rank\tIndex\tSample\tN_seqs\tN_reads\t%reads\t(%Total_reads)" in contents)
            self.assertTrue("Problems detected:\n * Underrepresented samples" in contents)
            self.assertTrue("   1\tTCCTGA\t\t1\t4\t100.0%\t(100.0%)" in contents)
            self.assertTrue("The following samples are underrepresented:" in contents)
            for line in ("AB1\tCGTGTAGG+GACCTGAA\t\t<0.1%",
	                 "AB2\tCGTGTAGG+ATGTAACT\t\t<0.1%",
	                 "CDE3\tGACCTGAA+CGTGTAGG\t\t<0.1%",
	                 "CDE4\tATGTAACT+CGTGTAGG\t\t<0.1%",):
                self.assertTrue(line in contents)
            # Expect at least 12 lines of content in total
            self.assertTrue(contents.count('\n') >= 12)

    def test_analyse_barcodes_with_bcl2fastq_dir_and_samplesheet_empty_index(self):
        """
        AnalyseBarcodes: bcl2fastq directory as input (with samplesheet, empty index)
        """
        # Make a mock bcl2fastq output directory
        datadir = MockIlluminaData(
            os.path.join(self.wd,
                         "200428_M00879_0087_000000000-AGEW9"),
            "bcl2fastq2",
            unaligned_dir="bcl2fastq",
            paired_end=True)
        datadir.add_fastq_batch("AB","AB1","AB1_S1")
        datadir.create()
        # Add data to Fastq files
        self._insert_fastq_reads(
            os.path.join(self.wd,"200428_M00879_0087_000000000-AGEW9"))
        # Create sample sheet with single empty index
        sample_sheet = os.path.join(self.wd,"custom_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
AB1,AB1,,,,,,,AB,
""")
        # Set up and run pipeline
        p = AnalyseBarcodes(bcl2fastq_dir=
                            os.path.join(self.wd,
                                         "200428_M00879_0087_000000000-AGEW9",
                                         "bcl2fastq"))
        exit_code = p.run(os.path.join(self.wd,"barcode_analysis"),
                          sample_sheet=sample_sheet,
                          working_dir=self.wd,
                          poll_interval=0.5)
        # Check outputs
        self.assertEqual(exit_code,0)
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis")),
                        "Missing dir: barcode_analysis")
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis",
                                                   "counts")),
                        "Missing dir: barcode_analysis/counts")
        self.assertTrue(os.path.isfile(
            os.path.join(
                self.wd,
                "barcode_analysis",
                "counts",
                "AB.AB1_S1_L001_R1_001.fastq.gz.counts")),
                        "Missing file: AB.AB1_S1_L001_R1_001.fastq.gz.counts")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.report")),
                        "Missing file: barcodes.report")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.xls")),
                        "Missing file: barcodes.xls")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.html")),
                        "Missing file: barcodes.html")
        # Check that the report content is non-trivial
        barcodes_report = os.path.join(self.wd,
                                       "barcode_analysis",
                                       "barcodes.report")
        with open(barcodes_report,'rt') as fp:
            contents = fp.read()
            self.assertTrue("Barcode analysis for lane #1" in contents)
            self.assertTrue("#Rank\tIndex\tSample\tN_seqs\tN_reads\t%reads\t(%Total_reads)" in contents)
            self.assertTrue("Problems detected:\n * Underrepresented samples" in contents)
            self.assertTrue("   1\tTCCTGA\t\t1\t1\t100.0%\t(100.0%)" in contents)
            self.assertTrue("The following samples are underrepresented:" in contents)
            self.assertTrue("AB1\t\t\t<0.1%" in contents)
            # Expect at least 12 lines of content in total
            self.assertTrue(contents.count('\n') >= 12)

    def test_analyse_barcodes_with_samplesheet(self):
        """
        AnalyseBarcodes: sample sheet as input
        """
        # Create sample sheet
        sample_sheet = os.path.join(self.wd,"custom_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
AB1,AB1,,,D701,CGTGTAGG,D501,GACCTGAA,AB,
AB2,AB2,,,D702,CGTGTAGG,D501,ATGTAACT,AB,
CDE3,CDE3,,,D701,GACCTGAA,D501,CGTGTAGG,CDE,
CDE4,CDE4,,,D702,ATGTAACT,D501,CGTGTAGG,CDE,
""")
        # Set up pipeline before bcl2fastq directory exists
        p = AnalyseBarcodes(sample_sheet=sample_sheet)
        # Create the bcl2fastq directory before running pipeline
        datadir = MockIlluminaData(
            os.path.join(self.wd,
                         "200428_M00879_0087_000000000-AGEW9"),
            "bcl2fastq2",
            unaligned_dir="bcl2fastq",
            paired_end=True)
        datadir.add_fastq_batch("AB","AB1","AB1_S1")
        datadir.add_fastq_batch("AB","AB2","AB2_S2")
        datadir.add_fastq_batch("CDE","CDE3","CDE3_S3")
        datadir.add_fastq_batch("CDE","CDE4","CDE4_S4")
        datadir.add_fastq_batch("","Undetermined","Undetermined_S0")
        datadir.create()
        # Add data to Fastq files
        self._insert_fastq_reads(
            os.path.join(self.wd,"200428_M00879_0087_000000000-AGEW9"))
        # Run the pipeline
        exit_code = p.run(os.path.join(self.wd,"barcode_analysis"),
                          bcl2fastq_dir=os.path.join(
                              self.wd,
                              "200428_M00879_0087_000000000-AGEW9",
                              "bcl2fastq"),
                          working_dir=self.wd,
                          poll_interval=0.5)
        # Check outputs
        self.assertEqual(exit_code,0)
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis")),
                        "Missing dir: barcode_analysis")
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis",
                                                   "counts")),
                            "Missing dir: barcode_analysis/counts")
        for f in (
                "AB.AB1_S1_L001_R1_001.fastq.gz.counts",
                "AB.AB2_S2_L001_R1_001.fastq.gz.counts",
                "CDE.CDE3_S3_L001_R1_001.fastq.gz.counts",
                "CDE.CDE4_S4_L001_R1_001.fastq.gz.counts",
                "__undetermined__.Undetermined_S0_L001_R1_001.fastq.gz.counts"):
            self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                        "barcode_analysis",
                                                        "counts",f)),
                            "Missing file: %s" % f)
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.report")),
                        "Missing file: barcodes.report")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.xls")),
                        "Missing file: barcodes.xls")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.html")),
                        "Missing file: barcodes.html")
        # Check that the report content is non-trivial
        barcodes_report = os.path.join(self.wd,
                                       "barcode_analysis",
                                       "barcodes.report")
        with open(barcodes_report,'rt') as fp:
            contents = fp.read()
            self.assertTrue("Barcode analysis for lane #1" in contents)
            self.assertTrue("#Rank\tIndex\tSample\tN_seqs\tN_reads\t%reads\t(%Total_reads)" in contents)
            self.assertTrue("Problems detected:\n * Underrepresented samples" in contents)
            self.assertTrue("   1\tTCCTGA\t\t1\t4\t100.0%\t(100.0%)" in contents)
            self.assertTrue("The following samples are underrepresented:" in contents)
            for line in ("AB1\tCGTGTAGG+GACCTGAA\t\t<0.1%",
	                 "AB2\tCGTGTAGG+ATGTAACT\t\t<0.1%",
	                 "CDE3\tGACCTGAA+CGTGTAGG\t\t<0.1%",
	                 "CDE4\tATGTAACT+CGTGTAGG\t\t<0.1%",):
                self.assertTrue(line in contents)
            # Expect at least 12 lines of content in total
            self.assertTrue(contents.count('\n') >= 12)

    def test_analyse_barcodes_with_multi_lane_samplesheet(self):
        """
        AnalyseBarcodes: multi-lane sample sheet as input
        """
        # Create sample sheet
        sample_sheet = os.path.join(self.wd,"custom_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,AB1,AB1,,,D701,CGTGTAGG,D501,GACCTGAA,AB,
1,AB2,AB2,,,D702,CGTGTAGG,D501,ATGTAACT,AB,
2,CDE3,CDE3,,,D701,GACCTGAA,D501,CGTGTAGG,CDE,
2,CDE4,CDE4,,,D702,ATGTAACT,D501,CGTGTAGG,CDE,
""")
        # Set up pipeline before bcl2fastq directory exists
        p = AnalyseBarcodes(sample_sheet=sample_sheet)
        # Create the bcl2fastq directory before running pipeline
        datadir = MockIlluminaData(
            os.path.join(self.wd,
                         "200428_M00879_0087_000000000-AGEW9"),
            "bcl2fastq2",
            unaligned_dir="bcl2fastq",
            paired_end=True)
        datadir.add_fastq_batch("AB","AB1","AB1_S1",lanes=(1,))
        datadir.add_fastq_batch("AB","AB2","AB2_S2",lanes=(1,))
        datadir.add_fastq_batch("CDE","CDE3","CDE3_S3",lanes=(2,))
        datadir.add_fastq_batch("CDE","CDE4","CDE4_S4",lanes=(2,))
        datadir.add_fastq_batch("","Undetermined","Undetermined_S0",
                                lanes=(1,2))
        datadir.create()
        # Add data to Fastq files
        self._insert_fastq_reads(
            os.path.join(self.wd,"200428_M00879_0087_000000000-AGEW9"))
        # Run the pipeline
        exit_code = p.run(os.path.join(self.wd,"barcode_analysis"),
                          bcl2fastq_dir=os.path.join(
                              self.wd,
                              "200428_M00879_0087_000000000-AGEW9",
                              "bcl2fastq"),
                          working_dir=self.wd,
                          poll_interval=0.5)
        # Check outputs
        self.assertEqual(exit_code,0)
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis")),
                        "Missing dir: barcode_analysis")
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis",
                                                   "counts")),
                            "Missing dir: barcode_analysis/counts")
        for f in (
                "AB.AB1_S1_L001_R1_001.fastq.gz.counts",
                "AB.AB2_S2_L001_R1_001.fastq.gz.counts",
                "CDE.CDE3_S3_L002_R1_001.fastq.gz.counts",
                "CDE.CDE4_S4_L002_R1_001.fastq.gz.counts",
                "__undetermined__.Undetermined_S0_L001_R1_001.fastq.gz.counts"):
            self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                        "barcode_analysis",
                                                        "counts",f)),
                            "Missing file: %s" % f)
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.report")),
                        "Missing file: barcodes.report")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.xls")),
                        "Missing file: barcodes.xls")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.html")),
                        "Missing file: barcodes.html")
        # Check that the report content is non-trivial
        barcodes_report = os.path.join(self.wd,
                                       "barcode_analysis",
                                       "barcodes.report")
        with open(barcodes_report,'rt') as fp:
            contents = fp.read()
            self.assertTrue("Barcode analysis for lane #1" in contents)
            self.assertTrue("Barcode analysis for lane #2" in contents)
            self.assertTrue("#Rank\tIndex\tSample\tN_seqs\tN_reads\t%reads\t(%Total_reads)" in contents)
            self.assertTrue("Problems detected:\n * Underrepresented samples" in contents)
            self.assertTrue("   1\tTCCTGA\t\t1\t2\t100.0%\t(100.0%)" in contents)
            self.assertTrue("The following samples are underrepresented:" in contents)
            for line in ("AB1\tCGTGTAGG+GACCTGAA\t\t<0.1%",
	                 "AB2\tCGTGTAGG+ATGTAACT\t\t<0.1%",
	                 "CDE3\tGACCTGAA+CGTGTAGG\t\t<0.1%",
	                 "CDE4\tATGTAACT+CGTGTAGG\t\t<0.1%",):
                self.assertTrue(line in contents)
            # Expect at least 12 lines of content in total
            self.assertTrue(contents.count('\n') >= 12)

    def test_analyse_barcodes_with_samplesheet_and_10x_indices(self):
        """
        AnalyseBarcodes: sample sheet with 10xGenomics indices
        """
        # Create sample sheet
        sample_sheet = os.path.join(self.wd,"custom_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
AB1,AB1,,,D501,SI-GA-A2,AB,
AB2,AB2,,,D501,SI-GA-B2,AB,
CDE3,CDE3,,,D501,SI-GA-C2,CDE,
CDE4,CDE4,,,D501,SI-GA-D2,CDE,
""")
        # Set up pipeline before bcl2fastq directory exists
        p = AnalyseBarcodes(sample_sheet=sample_sheet)
        # Create the bcl2fastq directory before running pipeline
        datadir = MockIlluminaData(
            os.path.join(self.wd,
                         "200428_M00879_0087_000000000-AGEW9"),
            "bcl2fastq2",
            unaligned_dir="bcl2fastq",
            paired_end=True)
        datadir.add_fastq_batch("AB","AB1","AB1_S1")
        datadir.add_fastq_batch("AB","AB2","AB2_S2")
        datadir.add_fastq_batch("CDE","CDE3","CDE3_S3")
        datadir.add_fastq_batch("CDE","CDE4","CDE4_S4")
        datadir.add_fastq_batch("","Undetermined","Undetermined_S0")
        datadir.create()
        # Add data to Fastq files
        self._insert_fastq_reads(
            os.path.join(self.wd,"200428_M00879_0087_000000000-AGEW9"))
        # Run the pipeline
        exit_code = p.run(os.path.join(self.wd,"barcode_analysis"),
                          bcl2fastq_dir=os.path.join(
                              self.wd,
                              "200428_M00879_0087_000000000-AGEW9",
                              "bcl2fastq"),
                          working_dir=self.wd,
                          poll_interval=0.5)
        # Check outputs
        self.assertEqual(exit_code,0)
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis")),
                        "Missing dir: barcode_analysis")
        self.assertTrue(os.path.isdir(os.path.join(self.wd,
                                                   "barcode_analysis",
                                                   "counts")),
                            "Missing dir: barcode_analysis/counts")
        for f in (
                "AB.AB1_S1_L001_R1_001.fastq.gz.counts",
                "AB.AB2_S2_L001_R1_001.fastq.gz.counts",
                "CDE.CDE3_S3_L001_R1_001.fastq.gz.counts",
                "CDE.CDE4_S4_L001_R1_001.fastq.gz.counts",
                "__undetermined__.Undetermined_S0_L001_R1_001.fastq.gz.counts"):
            self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                        "barcode_analysis",
                                                        "counts",f)),
                            "Missing file: %s" % f)
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.report")),
                        "Missing file: barcodes.report")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.xls")),
                        "Missing file: barcodes.xls")
        self.assertTrue(os.path.isfile(os.path.join(self.wd,
                                                    "barcode_analysis",
                                                    "barcodes.html")),
                        "Missing file: barcodes.html")
        # Check that the report content is non-trivial
        barcodes_report = os.path.join(self.wd,
                                       "barcode_analysis",
                                       "barcodes.report")
        with open(barcodes_report,'rt') as fp:
            contents = fp.read()
            self.assertTrue("Barcode analysis for lane #1" in contents)
            self.assertTrue("#Rank\tIndex\tSample\tN_seqs\tN_reads\t%reads\t(%Total_reads)" in contents)
            # Expect 12 lines of content in total
            self.assertEqual(contents.count('\n'),12)

    def test_analyse_barcodes_with_no_inputs(self):
        """
        AnalyseBarcodes: raise exception if no inputs supplied
        """
        self.assertRaises(Exception,
                          AnalyseBarcodes,
                          bcl2fastq_dir=None,
                          sample_sheet=None)

    def test_analyse_barcodes_with_bad_samplesheet(self):
        """
        AnalyseBarcodes: raise exception for 'bad' samplesheet as input
        """
        # Create "bad" sample sheet with mixture of empty and
        # non-empty indices in a lane
        sample_sheet = os.path.join(self.wd,"custom_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
AB1,AB1,,,D701,CGTGTAGG,D501,GACCTGAA,AB,
AB2,AB2,,,D702,CGTGTAGG,D501,ATGTAACT,AB,
CDE3,CDE3,,,,,,,CDE,
CDE4,CDE4,,,,,,,CDE,
""")
        self.assertRaises(Exception,
                          AnalyseBarcodes,
                          bcl2fastq_dir=None,
                          sample_sheet=sample_sheet)

    def test_analyse_barcodes_with_bcl2fastq_dir_and_bad_samplesheet(self):
        """
        AnalyseBarcodes: raise exception for bcl2fastq directory as input using 'bad' samplesheet
        """
        # Make a mock bcl2fastq output directory
        datadir = MockIlluminaData(
            os.path.join(self.wd,
                         "200428_M00879_0087_000000000-AGEW9"),
            "bcl2fastq2",
            unaligned_dir="bcl2fastq",
            paired_end=True)
        datadir.add_fastq_batch("AB","AB1","AB1_S1")
        datadir.add_fastq_batch("AB","AB2","AB2_S2")
        datadir.add_fastq_batch("CDE","CDE3","CDE3_S3")
        datadir.add_fastq_batch("CDE","CDE4","CDE4_S4")
        datadir.add_fastq_batch("","Undetermined","Undetermined_S0")
        datadir.create()
        # Add data to Fastq files
        self._insert_fastq_reads(
            os.path.join(self.wd,"200428_M00879_0087_000000000-AGEW9"))
        # Create "bad" sample sheet with mixture of empty and
        # non-empty indices in a lane
        sample_sheet = os.path.join(self.wd,"custom_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
AB1,AB1,,,D701,CGTGTAGG,D501,GACCTGAA,AB,
AB2,AB2,,,D702,CGTGTAGG,D501,ATGTAACT,AB,
CDE3,CDE3,,,,,,,CDE,
CDE4,CDE4,,,,,,,CDE,
""")
        # Set up and run pipeline
        p = AnalyseBarcodes(bcl2fastq_dir=
                            os.path.join(self.wd,
                                         "200428_M00879_0087_000000000-AGEW9",
                                         "bcl2fastq"))
        self.assertRaises(Exception,
                          AnalyseBarcodes.run,
                          os.path.join(self.wd,"barcode_analysis"),
                          sample_sheet=sample_sheet,
                          working_dir=self.wd,
                          poll_interval=0.5)
