#######################################################################
# Tests for update_fastq_stats_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.commands.update_fastq_stats_cmd import update_fastq_stats

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestAutoProcessUpdateFastqStats(unittest.TestCase):
    """
    Tests for AutoProcess.update_fastq_stats
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestAutoProcessMakeFastqs')
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

    def test_update_fastq_stats(self):
        """update_fastq_stats: generates statistics files
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '190104_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "190104" },
            top_dir=self.wd)
        mockdir.create(no_project_dirs=True)
        # Statistics files
        stats_files = (
            "statistics.info",
            "statistics_full.info",
            "per_lane_statistics.info",
            "per_lane_sample_stats.info",
        )
        # Check stats files don't already exist
        for filen in stats_files:
            self.assertFalse(os.path.exists(os.path.join(mockdir.dirn,filen)),
                             "%s: file exists, but shouldn't" %
                             filen)
        # Update (i.e. generate) stats
        ap = AutoProcess(mockdir.dirn)
        update_fastq_stats(ap)
        # Check files now exist
        for filen in stats_files:
            self.assertTrue(os.path.exists(os.path.join(mockdir.dirn,filen)),
                             "%s: missing" % filen)
