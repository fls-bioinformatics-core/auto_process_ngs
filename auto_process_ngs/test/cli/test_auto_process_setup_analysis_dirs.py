#######################################################################
# Tests for cli/auto_process.py "setup_analysis_dirs" command
#######################################################################

import unittest
import tempfile
import shutil
import os
from textwrap import dedent
from auto_process_ngs.cli.auto_process import main as auto_process
from auto_process_ngs.mock import MockAnalysisDirFactory

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Unit tests

class TestSetupAnalysisDirsCommand(unittest.TestCase):
    """
    Tests for the 'auto_process.py setup_analysis_dirs' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestSetupAnalysisDirsCmd')
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

    def test_setup_analysis_dirs(self):
        """
        auto_process.py setup_analysis_dirs: create new analysis dirs
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={"instrument_datestamp": "170901"},
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Update 'projects.info' file
        projects_info = os.path.join(mockdir.dirn, "projects.info")
        with open(projects_info, "wt") as fp:
            fp.write("""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tRNA-seq\t.\tHuman\tAudrey Benson\t1% PhiX
CDE\tCDE3,CDE4\tClive David Edwards\tChIP-seq\t.\tMouse\tClaudia Divine Eccleston\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz"],
            "CDE": ["CDE3_S3_R1_001.fastq.gz",
                    "CDE3_S3_R2_001.fastq.gz",
                    "CDE4_S4_R1_001.fastq.gz",
                    "CDE4_S4_R2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",
                             "Undetermined_S0_R1_001.fastq.gz"]
        }
        # Set up the project dirs
        auto_process(["setup_analysis_dirs", mockdir.dirn])
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn, project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir, fq)
                self.assertTrue(os.path.exists(fastq))

    def test_setup_analysis_dirs_with_custom_project_metadata(self):
        """
        auto_process.py setup_analysis_dirs: create new analysis dirs with custom metadata
        """
        # Extend settings file to add custom project metadata items
        with open(self.local_settings_file, "a") as fp:
            fp.write(dedent("""[metadata]
            custom_project_metadata = order_numbers,bioinformatics_analysts,external_analysts
            """))
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={"instrument_datestamp": "170901"},
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Update 'projects.info' file
        projects_info = os.path.join(mockdir.dirn, "projects.info")
        with open(projects_info, "wt") as fp:
            fp.write("""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tRNA-seq\t.\tHuman\tAudrey Benson\t1% PhiX
CDE\tCDE3,CDE4\tClive David Edwards\tChIP-seq\t.\tMouse\tClaudia Divine Eccleston\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz"],
            "CDE": ["CDE3_S3_R1_001.fastq.gz",
                    "CDE3_S3_R2_001.fastq.gz",
                    "CDE4_S4_R1_001.fastq.gz",
                    "CDE4_S4_R2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",
                             "Undetermined_S0_R1_001.fastq.gz"]
        }
        project_metadata_items = {
            "Project name",
            "Run",
            "Platform",
            "Sequencer model",
            "User",
            "PI",
            "Organism",
            "Library type",
            "Single cell platform",
            "Number of cells",
            "Paired_end",
            "Primary fastqs",
            "Samples",
            "Biological samples",
            "Multiplexed samples",
            "Comments",
            "Order numbers",
            "Bioinformatics analysts",
            "External analysts",
        }
        # Set up the project dirs
        auto_process(["setup_analysis_dirs", mockdir.dirn])
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn, project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check Fastqs
            fastqs_dir = os.path.join(project_dir_path,
                                      "fastqs")
            self.assertTrue(os.path.exists(fastqs_dir))
            for fq in projects[project]:
                fastq = os.path.join(fastqs_dir, fq)
                self.assertTrue(os.path.exists(fastq))
            # Check metadata items in README.info
            expected_metadata = project_metadata_items.copy()
            with open(readme_file, "rt") as fp:
                for line in fp:
                    item = line.split("\t")[0]
                    # Check that item is expected
                    self.assertTrue(item in expected_metadata,
                                    f"Found unexpected metadata item: '{item}'")
                    expected_metadata.remove(item)
            # Check no items are missing
            self.assertTrue(len(expected_metadata) == 0,
                            f"Missing project metadata items: '{expected_metadata}'")