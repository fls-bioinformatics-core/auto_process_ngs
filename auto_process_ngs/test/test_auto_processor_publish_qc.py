#######################################################################
# Tests for 'auto_process publish_qc' command
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.utils import ZipArchive
from auto_process_ngs.qc.illumina_qc import expected_qc_outputs
from auto_process_ngs.mock import MockAnalysisDirFactory

class TestAutoProcessPublishQc(unittest.TestCase):
    """
    Tests for AutoProcess.publish_qc
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessPublishQc')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)
        # Placeholders for test objects
        self.ap = None

    def tearDown(self):
        # Delete autoprocessor object
        if self.ap is not None:
            del(self.ap)
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def _make_file(self,filen,content=None):
        # Make a file
        with open(filen,'w') as fp:
            if content is not None:
                fp.write("%s" % content)
            else:
                fp.write("")

    def _add_processing_report(self,dirn):
        # Add mock processing report
        self._make_file(os.path.join(dirn,"processing_qc.html"))

    def test_publish_qc_metadata_missing(self):
        """publish_qc: raises exception if metadata not set
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn)
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish
        self.assertRaises(Exception,
                          ap.publish_qc,
                          location=publication_dir)

    def test_publish_qc_processing_qc(self):
        """publish_qc: processing QC report only
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        ap = AutoProcess(mockdir.dirn)
        # Add processing report
        self._add_processing_report(ap.analysis_dir)
        # Make a mock publication area
        publication_dir = os.path.join(self.dirn,'QC')
        os.mkdir(publication_dir)
        # Publish QC
        ap.publish_qc(location=publication_dir)
        self.assertTrue(
            os.path.exists(
                os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             "processing_qc.html")))
        # Check outputs
        outputs = ("index.html",
                   "processing_qc.html")
        for item in outputs:
            f = os.path.join(publication_dir,
                             "160621_K00879_0087_000000000-AGEW9_analysis",
                             item)
            self.assertTrue(os.path.exists(f),"Missing %s" % f)
