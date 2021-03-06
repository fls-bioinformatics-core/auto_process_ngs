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
                    self.assertTrue(os.path.exists(fastq),
                                    "Missing %s" % fastq)
                    with gzip.open(fastq,'wb') as fp:
                        fp.write(
                            """@ILLUMINA-545855:49:FC61RLR:2:1:10979:1695 1:N:0:TCCTGA
GCATACTCAGCTTTAGTAATAAGTGTGATTCTGGTA
+
IIIIIHIIIGHHIIDGHIIIIIIHIIIIIIIIIIIH
                            """.encode())

    def test_analyse_barcodes_with_bcl2fastq_dir(self):
        """
        AnalyseBarcodes: bcl2fastq directory as input
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

    def test_analyse_barcodes_with_no_inputs(self):
        """
        AnalyseBarcodes: raise exception if no inputs supplied
        """
        self.assertRaises(Exception,
                          AnalyseBarcodes,
                          bcl2fastq_dir=None,
                          sample_sheet=None)
