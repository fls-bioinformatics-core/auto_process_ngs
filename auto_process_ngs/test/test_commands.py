#######################################################################
# Tests for commands.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.commands import archive

# Unit tests

class TestArchiveCommand(unittest.TestCase):
    """
    Tests for the 'archive' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestArchiveCommand')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_archive_to_staging(self):
        """archive: test copying to staging archive dir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=False,
                         final=False)
        self.assertEqual(status,0)
        # Check that staging dir exists
        self.assertTrue(
            os.path.exists(
                os.path.join(
                    archive_dir,
                    "2017",
                    "miseq",
                    "__170901_M00879_0087_000000000-AGEW9_analysis.pending")))
        self.assertEqual(len(os.listdir(final_dir)),1)

    def test_archive_to_final(self):
        """archive: test copying to final archive dir
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Make a mock archive directory
        archive_dir = os.path.join(self.dirn,"archive")
        final_dir = os.path.join(archive_dir,
                                 "2017",
                                 "miseq")
        os.makedirs(final_dir)
        self.assertTrue(os.path.isdir(final_dir))
        self.assertEqual(len(os.listdir(final_dir)),0)
        # Make autoprocess instance and set required metadata
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        ap.set_metadata("source","testing")
        ap.set_metadata("run_number","87")
        # Do archiving op
        status = archive(ap,
                         archive_dir=archive_dir,
                         year='2017',platform='miseq',
                         read_only_fastqs=False,
                         final=True)
        self.assertEqual(status,0)
        # Check that staging dir exists
        self.assertTrue(
            os.path.exists(
                os.path.join(
                    archive_dir,
                    "2017",
                    "miseq",
                    "170901_M00879_0087_000000000-AGEW9_analysis")))
        self.assertEqual(len(os.listdir(final_dir)),1)
