#######################################################################
# Tests for clone_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.commands.clone_cmd import clone

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestAutoProcessClone(unittest.TestCase):
    """
    Tests for AutoProcess.clone
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessClone')
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
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn)

    def test_clone_analysis_dir(self):
        """
        clone: copies an analysis directory using symlinks
        """
        # Make a source analysis dir
        analysis_dir = MockAnalysisDirFactory.bcl2fastq2(
            "190116_M01234_0002_AXYZ123",
            platform="miseq",
            paired_end=True,
            no_lane_splitting=False,
            include_stats_files=True,
            top_dir=self.dirn)
        analysis_dir.create()
        # Make a copy
        clone_dir = os.path.join(self.dirn,"190116_M01234_0002_AXYZ123_copy")
        self.assertFalse(os.path.exists(clone_dir))
        ap = AutoProcess(analysis_dir.dirn)
        clone(ap,clone_dir)
        self.assertTrue(os.path.isdir(clone_dir))
        # Check contents
        for subdir in ('logs','ScriptCode'):
            d = os.path.join(clone_dir,subdir)
            self.assertTrue(os.path.isdir(d),"Missing '%s'" % subdir)
        for filen in ('SampleSheet.orig.csv',
                      'custom_SampleSheet.csv',
                      'auto_process.info',
                      'metadata.info',
                      'statistics.info',):
            f = os.path.join(clone_dir,filen)
            self.assertTrue(os.path.isfile(f),"Missing '%s'" % filen)
        # Check unassigned
        unassigned = os.path.join(clone_dir,'bcl2fastq')
        self.assertTrue(os.path.islink(unassigned))
        # Check projects
        for proj in ('AB','CDE','undetermined'):
            d = os.path.join(clone_dir,proj)
            self.assertTrue(os.path.isdir(d),"Missing '%s'" % proj)

    def test_clone_analysis_dir_copy_fastqs(self):
        """
        clone: copies an analysis directory
        """
        # Make a source analysis dir
        analysis_dir = MockAnalysisDirFactory.bcl2fastq2(
            "190116_M01234_0002_AXYZ123",
            platform="miseq",
            paired_end=True,
            no_lane_splitting=False,
            include_stats_files=True,
            top_dir=self.dirn)
        analysis_dir.create()
        # Make a copy
        clone_dir = os.path.join(self.dirn,"190116_M01234_0002_AXYZ123_copy")
        self.assertFalse(os.path.exists(clone_dir))
        ap = AutoProcess(analysis_dir.dirn)
        clone(ap,clone_dir,copy_fastqs=True)
        self.assertTrue(os.path.isdir(clone_dir))
        # Check contents
        for subdir in ('logs','ScriptCode'):
            d = os.path.join(clone_dir,subdir)
            self.assertTrue(os.path.isdir(d),"Missing '%s'" % subdir)
        for filen in ('SampleSheet.orig.csv',
                      'custom_SampleSheet.csv',
                      'auto_process.info',
                      'metadata.info',
                      'statistics.info',):
            f = os.path.join(clone_dir,filen)
            self.assertTrue(os.path.isfile(f),"Missing '%s'" % filen)
        # Check unassigned
        unassigned = os.path.join(clone_dir,'bcl2fastq')
        self.assertTrue(os.path.isdir(unassigned))
        # Check projects
        for proj in ('AB','CDE','undetermined'):
            d = os.path.join(clone_dir,proj)
            self.assertTrue(os.path.isdir(d),"Missing '%s'" % proj)

    def test_clone_fails_if_target_dir_exists(self):
        """
        clone: raises an exception if target dir already exists 
        """
        # Make a source analysis dir
        analysis_dir = MockAnalysisDirFactory.bcl2fastq2(
            "190116_M01234_0002_AXYZ123",
            platform="miseq",
            paired_end=True,
            no_lane_splitting=False,
            include_stats_files=True,
            top_dir=self.dirn)
        analysis_dir.create()
        # Make target dir
        clone_dir = os.path.join(self.dirn,"190116_M01234_0002_AXYZ123_copy")
        os.mkdir(clone_dir)
        # Try to copy source dir
        ap = AutoProcess(analysis_dir.dirn)
        self.assertRaises(Exception,
                          clone,
                          ap,clone_dir)

        
