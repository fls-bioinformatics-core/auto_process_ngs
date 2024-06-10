#######################################################################
# Base module for unit tests for qc/pipeline.py
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.metadata import AnalysisProjectInfo
from auto_process_ngs.metadata import AnalysisProjectQCDirInfo
from auto_process_ngs.mock import MockFastqScreen
from auto_process_ngs.mock import MockFastQC
from auto_process_ngs.mock import MockFastqStrandPy
from auto_process_ngs.mock import MockGtf2bed
from auto_process_ngs.mock import MockSeqtk
from auto_process_ngs.mock import MockStar
from auto_process_ngs.mock import MockSamtools
from auto_process_ngs.mock import MockPicard
from auto_process_ngs.mock import MockRSeQC
from auto_process_ngs.mock import MockQualimap
from auto_process_ngs.mock import MockMultiQC
from auto_process_ngs.mock import MockCellrangerExe
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.qc.pipeline import QCPipeline
from auto_process_ngs.qc.protocols import QCProtocol
from auto_process_ngs.qc.protocols import fetch_protocol_definition

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Polling interval for pipeline
POLL_INTERVAL = 0.1

class BaseQCPipelineTestCase(unittest.TestCase):
    """
    Base class for tests of QCPipeline class

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
    """    
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestQCPipeline')
        # Create a temp 'bin' dir
        self.bin = os.path.join(self.wd,"bin")
        os.mkdir(self.bin)
        # Create a temp 'data' dir
        self.data = os.path.join(self.wd,"data")
        os.mkdir(self.data)
        # Add (empty) FastqScreen conf files
        self.fastq_screens = {
            'model_organisms': "fastq_screen_model_organisms.conf",
            'other_organisms': "fastq_screen_other_organisms.conf",
            'rRNA': "fastq_screen_rRNA.conf",
        }
        for screen in self.fastq_screens:
            conf_file = os.path.join(self.data,
                                     self.fastq_screens[screen])
            with open(conf_file,'wt') as fp:
                fp.write("")
            self.fastq_screens[screen] = conf_file
        # Add (empty) reference data files
        self.ref_data = dict()
        for build in ('hg38','mm10',):
            self.ref_data[build] = {}
            build_dir = os.path.join(self.data,build)
            os.mkdir(build_dir)
            for ext in ('bed','gtf'):
                f = os.path.join(build_dir,"%s.%s" % (build,ext))
                with open(f,'wt') as fp:
                    fp.write("")
                self.ref_data[build][ext] = f
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
