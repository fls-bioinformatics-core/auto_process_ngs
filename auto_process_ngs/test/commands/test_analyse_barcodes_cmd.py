#######################################################################
# Tests for analyse_barcodes_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
import gzip
from bcftbx.IlluminaData import IlluminaData
from auto_process_ngs.settings import Settings
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.commands.analyse_barcodes_cmd import analyse_barcodes

# Set to False to preserve test outputs
REMOVE_TEST_OUTPUTS = True

class TestAutoProcessAnalyseBarcodes(unittest.TestCase):
    """
    Tests for AutoProcess.analyse_barcodes
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestAutoProcessAnalyseBarcodes')
        # Create settings instance
        # This allows us to set the polling interval for the
        # unit tests
        settings_ini = os.path.join(self.wd,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 0.5
""")
        self.settings = Settings(settings_ini)
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

    def test_analyse_barcodes_with_stored_bases_mask(self):
        """analyse_barcodes: test with stored bases mask
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            bases_mask='y76,I6,y76',
            metadata={ "instrument_datestamp": "160621" },
            top_dir=self.wd)
        mockdir.create(no_project_dirs=True)
        # Add data to Fastq files
        self._insert_fastq_reads(mockdir.dirn)
        # Populate the samplesheet
        sample_sheet = os.path.join(mockdir.dirn,"custom_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
AB1,AB1,,,D701,CGTGTAGG,D501,GACCTGAA,AB,
AB2,AB2,,,D702,CGTGTAGG,D501,ATGTAACT,AB,
CDE3,CDE3,,,D701,GACCTGAA,D501,CGTGTAGG,CDE,
CDE4,CDE4,,,D702,ATGTAACT,D501,CGTGTAGG,CDE,
""")
        # Analyse barcodes
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        analyse_barcodes(ap)
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

    def test_analyse_barcodes_no_stored_bases_mask(self):
        """analyse_barcodes: test with no stored bases mask
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
        # Populate the samplesheet
        sample_sheet = os.path.join(mockdir.dirn,"custom_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
AB1,AB1,,,D701,CGTGTAGG,D501,GACCTGAA,AB,
AB2,AB2,,,D702,CGTGTAGG,D501,ATGTAACT,AB,
CDE3,CDE3,,,D701,GACCTGAA,D501,CGTGTAGG,CDE,
CDE4,CDE4,,,D702,ATGTAACT,D501,CGTGTAGG,CDE,
""")
        # Analyse barcodes
        ap = AutoProcess(mockdir.dirn,
                         settings=self.settings)
        analyse_barcodes(ap)
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

        
        
