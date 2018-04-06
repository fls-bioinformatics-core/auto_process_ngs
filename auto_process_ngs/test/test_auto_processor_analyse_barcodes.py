#######################################################################
# Tests for autoprocessor.py module: make_fastqs
#######################################################################

import unittest
import tempfile
import shutil
import os
import gzip
from bcftbx.IlluminaData import IlluminaData
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory

# Set to False to preserve test outputs
REMOVE_TEST_OUTPUTS = True

class TestAutoProcessMakeFastqs(unittest.TestCase):
    """
    Tests for AutoProcess.make_fastqs
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestAutoProcessAnalyseBarcodes')
        # Store original location
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.wd)
        # Placeholders for test objects
        self.ap = None

    def tearDown(self):
        # Delete autoprocessor object
        if self.ap is not None:
            del(self.ap)
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def _insert_fastq_reads(self,dirn):
        # Add read data to fastqs
        illuminadata = IlluminaData(dirn,
                                    unaligned_dir="bcl2fastq")
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
""")

    def test_analyse_barcodes(self):
        """analyse_barcodes: test with defaults
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "160621" },
            top_dir=self.wd)
        mockdir.create(no_project_dirs=True)
        # Add data to Fastq files
        self._insert_fastq_reads(mockdir.dirn)
        # Analyse barcodes
        ap = AutoProcess(mockdir.dirn)
        ap.analyse_barcodes()
        # Check outputs
        analysis_dir = os.path.join(
            self.wd,
            "160621_M00879_0087_000000000-AGEW9_analysis")
        self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,"barcode_analysis")),
                            "Missing dir: barcode_analysis")
        self.assertTrue(os.path.isdir(
                os.path.join(analysis_dir,"barcode_analysis","counts")),
                            "Missing dir: barcode_analysis/counts")
        for f in ("AB.AB1_S1_R1_001.fastq.gz.counts",
                  "AB.AB2_S2_R1_001.fastq.gz.counts",
                  "CDE.CDE3_S3_R1_001.fastq.gz.counts",
                  "CDE.CDE4_S4_R1_001.fastq.gz.counts",
                  "undetermined.Undetermined_S0_R1_001.fastq.gz.counts"):
            self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,"barcode_analysis","counts",f)),
                            "Missing file: %s" % f)
        self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,
                             "barcode_analysis",
                             "barcodes.report")),
                        "Missing file: barcodes.report")
        self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,
                             "barcode_analysis",
                             "barcodes.xls")),
                        "Missing file: barcodes.xls")
        self.assertTrue(os.path.isfile(
                os.path.join(analysis_dir,
                             "barcode_analysis",
                             "barcodes.html")),
                        "Missing file: barcodes.html")

        
        
