#######################################################################
# Tests for setup_analysis_dirs_cmd.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
from textwrap import dedent
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.commands.setup_analysis_dirs_cmd import setup_analysis_dirs
from auto_process_ngs.settings import Settings

# Unit tests

class TestSetupAnalysisDirs(unittest.TestCase):
    """
    Tests for the 'setup_analysis_dirs' command
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestSetupAnalysisDirsCommand')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        def del_rw(action,name,excinfo):
            # Explicitly remove read only files/
            # dirs
            os.chmod(os.path.dirname(name),0o755)
            os.chmod(name,0o655)
            os.remove(name)
        shutil.rmtree(self.dirn,onerror=del_rw)

    def test_setup_analysis_dirs(self):
        """
        setup_analysis_dirs: test create new analysis dirs
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Check project dirs don't exist
        for project in ("AB","CDE"):
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
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
        project_metadata_items = { "Project name",
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
                                   "Comments" }
        # Set up the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
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

    def test_setup_analysis_dirs_with_identifier(self):
        """
        setup_analysis_dirs: test create new analysis dirs with an identifier
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Check project dirs don't exist
        for project in ("AB","CDE"):
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
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
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap,name="test_id")
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,
                                            "%s_%s" % (project,
                                                       "test_id"))
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))

    def test_setup_analysis_dirs_ignore_commented_projects(self):
        """
        setup_analysis_dirs: ignore commented line in 'projects.info'
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Check project dirs don't exist
        for project in ("AB","CDE"):
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
#AB\tAB1,AB2\tAlan Brown\tRNA-seq\t.\tHuman\tAudrey Benson\t1% PhiX
CDE\tCDE3,CDE4\tClive David Edwards\tChIP-seq\t.\tMouse\tClaudia Divine Eccleston\t1% PhiX
""")
        # Expected data
        projects = {
            "CDE": ["CDE3_S3_R1_001.fastq.gz",
                    "CDE3_S3_R2_001.fastq.gz",
                    "CDE4_S4_R1_001.fastq.gz",
                    "CDE4_S4_R2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",
                             "Undetermined_S0_R1_001.fastq.gz"]
        }
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
        # Check that 'AB' and '#AB' projects weren't created
        for project in ('AB','#AB'):
            self.assertFalse(os.path.exists(
                os.path.join(mockdir.dirn,project)))

    def test_setup_analysis_dirs_with_custom_metadata_from_config(self):
        """
        setup_analysis_dirs: test create new analysis dirs with custom metadata from config
        """
        # Make a minimal config file
        local_settings_file = os.path.join(self.dirn, "local_settings.ini")
        with open(local_settings_file, "wt") as fp:
            fp.write(dedent("""[metadata]
            custom_project_metadata = order_numbers,bioinformatics_analysts,external_analysts
            """))
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Check project dirs don't exist
        for project in ("AB","CDE"):
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
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
        project_metadata_items = { "Project name",
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
                                   "External analysts"}
        # Set up the project dirs with additional custom metadata items
        ap = AutoProcess(analysis_dir=mockdir.dirn,
                         settings=Settings(local_settings_file))
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check metadata items
            expected_metadata = project_metadata_items.copy()
            with open(readme_file, "rt") as fp:
                for line in fp:
                    item = line.split("\t")[0]
                    # Check that item is expected
                    self.assertTrue(item in expected_metadata,
                                    f"{project}: found unexpected metadata item: '{item}'")
                    expected_metadata.remove(item)
            # Check no items are missing
            self.assertTrue(len(expected_metadata) == 0,
                            f"{project}: missing project metadata items: '{expected_metadata}'")

    def test_setup_analysis_dirs_with_explicit_custom_metadata(self):
        """
        setup_analysis_dirs: test create new analysis dirs with explicit custom metadata
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Check project dirs don't exist
        for project in ("AB","CDE"):
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
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
        project_metadata_items = { "Project name",
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
                                   "External analysts"}
        # Set up the project dirs with additional custom metadata items
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap, custom_metadata_items=["order_numbers",
                                                       "bioinformatics_analysts",
                                                       "external_analysts"])
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertTrue(os.path.exists(project_dir_path))
            # Check README.info file
            readme_file = os.path.join(project_dir_path,
                                       "README.info")
            self.assertTrue(os.path.exists(readme_file))
            # Check metadata items
            expected_metadata = project_metadata_items.copy()
            with open(readme_file, "rt") as fp:
                for line in fp:
                    item = line.split("\t")[0]
                    # Check that item is expected
                    self.assertTrue(item in expected_metadata,
                                    f"{project}: found unexpected metadata item: '{item}'")
                    expected_metadata.remove(item)
            # Check no items are missing
            self.assertTrue(len(expected_metadata) == 0,
                            f"{project}: missing project metadata items: '{expected_metadata}'")

    def test_setup_analysis_dirs_10x_visium(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Visium
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            reads=('R1','R2','R3','I1'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tFFPE Spatial GEX\t10xGenomics Visium\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_R3_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_R3_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
        # Check that Visium project includes directory
        # for images
        self.assertTrue(os.path.isdir(os.path.join(mockdir.dirn,
                                                   "AB",
                                                   "Visium_images")))

    def test_setup_analysis_dirs_10x_visium_cytassist(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Visium (CytAssist)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            reads=('R1','R2','R3','I1'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tFFPE Spatial GEX\t10xGenomics Visium (CytAssist)\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_R3_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_R3_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
        # Check that Visium project includes directory
        # for images
        self.assertTrue(os.path.isdir(os.path.join(mockdir.dirn,
                                                   "AB",
                                                   "Visium_images")))

    def test_setup_analysis_dirs_10x_cytassist_visium_legacy(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x CytAssist Visium (legacy)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            reads=('R1','R2','R3','I1'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tFFPE Spatial RNA-seq\t10xGenomics CytAssist Visium\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_R3_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_R3_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
        # Check that Visium project includes directory
        # for images
        self.assertTrue(os.path.isdir(os.path.join(mockdir.dirn,
                                                   "AB",
                                                   "Visium_images")))

    def test_setup_analysis_dirs_10x_multiome(self):
        """
        setup_analysis_dirs: test create new analysis dir for 10x Multiome
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            reads=('R1','R2','I1'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tGEX\t10xGenomics Single Cell Multiome\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multiome_libraries.info
            template_libraries_file = os.path.join(
                project_dir_path,
                "10x_multiome_libraries.info.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_libraries_file),
                                "Missing %s" % template_libraries_file)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_libraries_file),
                                 "Found %s" % template_libraries_file)

    def test_setup_analysis_dirs_10x_cellplex_800(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x CellPlex (8.0.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/8.0.0/cellranger",
                        "cellranger",
                        "8.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tscRNA-seq\t10x Chromium 3' CellPlex\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz",
                   "AB2_S2_I2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_cellplex_900(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x CellPlex (9.0.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/9.0.0/cellranger",
                        "cellranger",
                        "9.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tscRNA-seq\t10x Chromium 3' CellPlex\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz",
                   "AB2_S2_I2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_cellplex_710_legacy_metadata(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x CellPlex (7.1.0, legacy metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/7.1.0/cellranger",
                        "cellranger",
                        "7.1.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tCellPlex scRNA-seq\t10xGenomics Chromium 3'v3\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz",
                   "AB2_S2_I2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain '#no-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\n#no-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_cellplex_800_legacy_metadata(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x CellPlex (8.0.0, legacy metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/8.0.0/cellranger",
                        "cellranger",
                        "8.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tCellPlex scRNA-seq\t10xGenomics Chromium 3'v3\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz",
                   "AB2_S2_I2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_cellplex_900_legacy_metadata(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x CellPlex (9.0.0, legacy metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/9.0.0/cellranger",
                        "cellranger",
                        "9.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tCellPlex scRNA-seq\t10xGenomics Chromium 3'v3\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",
                   "AB2_S2_R1_001.fastq.gz",
                   "AB2_S2_R2_001.fastq.gz",
                   "AB2_S2_I1_001.fastq.gz",
                   "AB2_S2_I2_001.fastq.gz"],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_flex_800(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Flex (8.0.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={
                "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/8.0.0/cellranger",
                        "cellranger",
                        "8.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tGEX\t10x Chromium Flex\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_flex_900(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Flex (9.0.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={
                "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/9.0.0/cellranger",
                        "cellranger",
                        "9.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tGEX\t10x Chromium Flex\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_flex_710_legacy_metadata(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Flex (7.1.0, legacy metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={
                "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/7.1.0/cellranger",
                        "cellranger",
                        "7.1.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tFlex\t10xGenomics Chromium 3'v3\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain '#no-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\n#no-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_flex_800_legacy_metadata(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Flex (8.0.0, legacy metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={
                "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/8.0.0/cellranger",
                        "cellranger",
                        "8.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tFlex\t10xGenomics Chromium 3'v3\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_flex_900_legacy_metadata(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Flex (9.0.0, legacy metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={
                "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/9.0.0/cellranger",
                        "cellranger",
                        "9.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tFlex\t10xGenomics Chromium 3'v3\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_gem_x_flex_legacy_metadata(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x GEM-X Flex (legacy metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={
                "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/9.0.0/cellranger",
                        "cellranger",
                        "9.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tFlex\t10xGenomics Chromium GEM-X\tHuman\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_immune_profiling_800(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Single Cell Immune Profiling (8.0.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/8.0.0/cellranger",
                        "cellranger",
                        "8.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tImmune Profiling\t10x Chromium 5'\tMouse\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_immune_profiling_900(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Single Cell Immune Profiling (9.0.0)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/9.0.0/cellranger",
                        "cellranger",
                        "9.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tImmune Profiling\t10x Chromium 5'\tMouse\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_immune_profiling_710_legacy_metadata(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Single Cell Immune Profiling (7.1.0, legacy metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/7.1.0/cellranger",
                        "cellranger",
                        "7.1.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tSingle Cell Immune Profiling\t10xGenomics Chromium 5'\tMouse\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain '#no-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\n#no-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_immune_profiling_800_legacy_metadata(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Single Cell Immune Profiling (8.0.0, legacy metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/8.0.0/cellranger",
                        "cellranger",
                        "8.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tSingle Cell Immune Profiling\t10xGenomics Chromium 5'\tMouse\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_10x_immune_profiling_900_legacy_metadata(self):
        """
        setup_analysis_dirs: create new analysis dir for 10x Single Cell Immune Profiling (9.0.0, legacy metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901",
                "processing_software": {
                    "cellranger" : (
                        "/usr/local/cellranger/9.0.0/cellranger",
                        "cellranger",
                        "9.0.0"
                    )
                }
            },
            reads=('R1','R2','I1','I2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tSingle Cell Immune Profiling\t10xGenomics Chromium 5'\tMouse\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",
                   "AB1_S1_I1_001.fastq.gz",
                   "AB1_S1_I2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))
            # Check template 10x_multi_config.csv
            template_multi_config_file = os.path.join(
                project_dir_path,
                "10x_multi_config.csv.template")
            if project == "AB":
                # Template file should exist
                self.assertTrue(os.path.exists(template_multi_config_file),
                                "Missing %s" % template_multi_config_file)
                # Template should contain 'create-bam'
                with open(template_multi_config_file,'rt') as fp:
                    template = fp.read()
                    self.assertTrue(template.find("\ncreate-bam,") != -1)
            else:
                # No template file
                self.assertFalse(os.path.exists(template_multi_config_file),
                                 "Found %s" % template_multi_config_file)

    def test_setup_analysis_dirs_biorad_ddseq_single_cell_rnaseq(self):
        """
        setup_analysis_dirs: create new analysis dir for Bio-Rad ddSEQ Single Cell 3' RNA-Seq
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '250321_A01899_0263_AHVHGMDRX',
            'novaseq6000',
            metadata={ "instrument_datestamp": "250321",
                "processing_software": {
                    "bcl2fastq" : (
                        "/usr/local/bcl2fastq/2.20/bcl2fastq",
                        "bcl2fastq",
                        "2.20"
                    )
                }
            },
            reads=('R1','R2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tscRNA-seq\tBio-Rad ddSEQ Single Cell 3' RNA-Seq\tMouse\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))

    def test_setup_analysis_dirs_biorad_ddseq_single_cell_atac(self):
        """
        setup_analysis_dirs: create new analysis dir for Bio-Rad ddSEQ Single Cell ATAC
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '250321_A01899_0263_AHVHGMDRX',
            'novaseq6000',
            metadata={ "instrument_datestamp": "250321",
                "processing_software": {
                    "bcl2fastq" : (
                        "/usr/local/bcl2fastq/2.20/bcl2fastq",
                        "bcl2fastq",
                        "2.20"
                    )
                }
            },
            reads=('R1','R2'),
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq")))
        print(os.listdir(os.path.join(mockdir.dirn,"bcl2fastq","AB")))
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1\tAlan Brown\tscATAC-seq\tBio-Rad ddSEQ Single Cell ATAC\tMouse\tAudrey Benson\t1% PhiX
""")
        # Expected data
        projects = {
            "AB": ["AB1_S1_R1_001.fastq.gz",
                   "AB1_S1_R2_001.fastq.gz",],
            "undetermined": ["Undetermined_S0_R1_001.fastq.gz",]
        }
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))

    def test_setup_analysis_dirs_missing_metadata(self):
        """
        setup_analysis_dirs: raise exception if metadata not set
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Attempt to set up the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        self.assertRaises(Exception,
                          setup_analysis_dirs,ap)

    def test_setup_analysis_dirs_ignore_missing_metadata(self):
        """
        setup_analysis_dirs: test create new analysis dirs (ignore missing metadata)
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
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
        # Check project dirs don't exist
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
            self.assertFalse(os.path.exists(project_dir_path))
        # Setup the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        setup_analysis_dirs(ap,ignore_missing_metadata=True)
        # Check project dirs and contents
        for project in projects:
            project_dir_path = os.path.join(mockdir.dirn,project)
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
                fastq = os.path.join(fastqs_dir,fq)
                self.assertTrue(os.path.exists(fastq))

    def test_setup_analysis_dirs_bad_single_cell_platform(self):
        """
        setup_analysis_dirs: raise exception for bad single cell platform
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '170901_M00879_0087_000000000-AGEW9',
            'miseq',
            metadata={ "instrument_datestamp": "170901" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Add required metadata to 'projects.info'
        projects_info = os.path.join(mockdir.dirn,"projects.info")
        with open(projects_info,"w") as fp:
            fp.write(
"""#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\tAlan Brown\tRNA-seq\t.\tHuman\tAudrey Benson\t1% PhiX
CDE\tCDE3,CDE4\tClive David Edwards\tscRNA-seq\t11xGenomics Chromium 5'v0\tMouse\tClaudia Divine Eccleston\t1% PhiX
""")
        # Attempt to set up the project dirs
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        self.assertRaises(Exception,
                          setup_analysis_dirs,ap)
