#######################################################################
# Tests for cli/auto_process.py "metadata" command
#######################################################################

import unittest
import tempfile
import shutil
import os
from textwrap import dedent
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.cli.auto_process import main as auto_process
from auto_process_ngs.mock import MockAnalysisDirFactory

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Unit tests

class TestMetadataCommand(unittest.TestCase):
    """
    Tests for the 'auto_process.py metadata' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestMetadataCmd')
        # Create empty settings file
        self.local_settings_file = os.path.join(self.dirn, "auto_process.ini")
        with open(self.local_settings_file, "wt") as s:
            s.write(dedent(f"""
            """))
        # Temporary point config to local version
        self.auto_process_conf = os.environ.get('AUTO_PROCESS_CONF')
        os.environ['AUTO_PROCESS_CONF'] = self.local_settings_file

    def tearDown(self):
        # Restore configuration environment variable
        if self.auto_process_conf is not None:
            os.environ['AUTO_PROCESS_CONF'] = self.auto_process_conf
        else:
            del(os.environ['AUTO_PROCESS_CONF'])
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            if os.path.isfile(name):
                os.chmod(os.path.dirname(name),0o755)
                os.chmod(name,0o655)
                os.remove(name)
            elif os.path.isdir(name):
                os.chmod(os.path.dirname(name),0o755)
                os.chmod(name,0o755)
                os.rmdir(name)
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn,onerror=del_rw)

    def test_set_analysis_dir_metadata(self):
        """
        auto_process.py metadata: set analysis directory metadata
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": None,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Set top-level metadata
        self.assertEqual(auto_process(["metadata",
                                       "--set",
                                       "run_number=87",
                                       mockdir.dirn]), 0)
        # Check metadata was updated
        self.assertEqual(AutoProcess(mockdir.dirn).metadata.run_number, 87)

    def test_set_project_metadata(self):
        """
        auto_process.py metadata: set project metadata
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": None,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" }
            },
            top_dir=self.dirn)
        mockdir.create()
        # Set project-level metadata item
        self.assertEqual(auto_process(["metadata",
                                       "--set",
                                       "AB:organism=Mouse",
                                       mockdir.dirn]), 0)
        # Check metadata was updated
        self.assertTrue("organism" in AnalysisProject(os.path.join(mockdir.dirn, "AB")).info)
        self.assertEqual(
            AnalysisProject(os.path.join(mockdir.dirn, "AB")).info.organism,
            "Mouse")

    def test_set_custom_project_metadata(self):
        """
        auto_process.py metadata: set custom metadata for a project
        """
        # Update local settings
        with open(self.local_settings_file, "a") as fp:
            fp.write(dedent("""[metadata]
            custom_project_metadata = order_numbers
            """))
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": None,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq" },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq" },
            },
            top_dir=self.dirn)
        mockdir.create()
        # Set project-level metadata item
        self.assertEqual(auto_process(["metadata",
                                       "--set",
                                       "AB:order_numbers=#00124",
                                       mockdir.dirn]), 0)
        # Check metadata was updated
        self.assertTrue("order_numbers" in AnalysisProject(os.path.join(mockdir.dirn, "AB")).info)
        self.assertEqual(
            AnalysisProject(os.path.join(mockdir.dirn, "AB")).info.order_numbers,
            "#00124")

    def test_update_existing_custom_project_metadata(self):
        """
        auto_process.py metadata: update existing custom metadata for a project
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "source": "testing",
                       "run_number": None,
                       "sequencer_model": "MiSeq" },
            project_metadata={
                "AB": { "User": "Alison Bell",
                        "Library type": "RNA-seq",
                        "Organism": "Human",
                        "PI": "Audrey Bower",
                        "Sequencer model": "MiSeq",
                        "Order numbers": None },
                "CDE": { "User": "Charles David Edwards",
                         "Library type": "ChIP-seq",
                         "Organism": "Mouse",
                         "PI": "Colin Delaney Eccleston",
                         "Sequencer model": "MiSeq",
                         "Order numbers": None },
            },
            top_dir=self.dirn)
        mockdir.create()
        # Set project-level metadata item
        self.assertEqual(auto_process(["metadata",
                                       "--set",
                                       "CDE:order_numbers=#00125",
                                       mockdir.dirn]), 0)
        # Check metadata was updated
        self.assertTrue("order_numbers" in AnalysisProject(os.path.join(mockdir.dirn, "CDE")).info)
        self.assertEqual(
            AnalysisProject(os.path.join(mockdir.dirn, "CDE")).info.order_numbers,
            "#00125")