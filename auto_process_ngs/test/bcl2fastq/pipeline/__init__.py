#######################################################################
# Base module for unit tests for bcl2fastq/pipeline.py
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.mock import MockIlluminaRun
from bcftbx.mock import SampleSheets
from bcftbx.IlluminaData import IlluminaData
from auto_process_ngs.mock import MockBcl2fastq2Exe
from auto_process_ngs.mock import MockBclConvertExe
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import Mock10xPackageExe
from auto_process_ngs.mock import make_mock_bcl2fastq2_output
from auto_process_ngs.bcl2fastq.pipeline import MakeFastqs
from auto_process_ngs.bcl2fastq.pipeline import subset
from auto_process_ngs.bcl2fastq.utils import make_custom_sample_sheet
from auto_process_ngs.stats import FastqStatistics

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Polling interval for pipeline
POLL_INTERVAL = 0.1

class BaseMakeFastqsTestCase(unittest.TestCase):
    """
    Base class for tests of MakeFastqs class

    Provides setUp() and tearDown() wrappers and common
    environment for tests.

    Following properties are available to subclasses:

    - wd: path to temporary working directory
    - bin: path to a 'bin' directory for mock executables
    - data: path to a temporary 'data' directory with mock
      reference data files

    setUp() moves to the working directory automatically

    Tests must populate the 'bin' directory themselves with
    the required mock executables. The 'bin' directory is
    automatically prepended to the PATH.

    tearDown() moves back to the original directory and
    restores the PATH environment variable. It also
    deletes the working directory and all its contents
    unless the module-wide 'REMOVE_TEST_OUTPUTS' variable
    is set to False.

    Additionally provides a method assertFastqStats()
    which can be used to check the statistics for a Fastq
    file within the tests.
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestMakeFastqs')
        # Create a temp 'bin' dir
        self.bin = os.path.join(self.wd,"bin")
        os.mkdir(self.bin)
        # Store original location
        self.pwd = os.getcwd()
        # Store original PATH
        self.path = os.environ['PATH']
        # Move to working dir
        os.chdir(self.wd)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Restore PATH
        os.environ['PATH'] = self.path
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def assertFastqStats(self,stats,fq,expected_nreads,**lanes):
        """
        Internal: check the statistics on a Fastq file

        'stats' should be a FastqStatistics object
        'fq' should be a Fastq file name (no path)
        'expected_nreads' is the total number of expected reads
        'lanes' can be lane identifiers specifying the expected
        read count for that lane e.g. 'L1=4' means 'expect 4
        reads associated with lane 1'
        """
        # Look up the Fastq
        data = stats.raw.lookup('Fastq',fq)
        if not data:
            raise AssertionError("'%s': Fastq not found" % fq)
        # Check number of reads
        self.assertEqual(data[0]['Nreads'],
                         expected_nreads,
                         "'%s': contains %d reads, expected %d"
                         % (fq,data[0]['Nreads'],expected_nreads))
        # Check reads per lane
        for lane in lanes:
            expected_lane_nreads = int(lanes[lane])
            try:
                nreads = data[0][lane]
                if not nreads:
                    nreads = 0
            except KeyError:
                raise AssertionError("'%s': lane '%s' not found"
                                     % (fq,lane))
            self.assertEqual(nreads,
                             expected_lane_nreads,
                             "'%s': has %d reads for lane %s, "
                             "expected %d"
                             % (fq,nreads,lane,expected_lane_nreads))
