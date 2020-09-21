#######################################################################
# Tests for autoprocessor.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.mock import MockIlluminaRun
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.analysis import AnalysisProject
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import MockAnalysisProject

# Unit tests

class TestAutoProcess(unittest.TestCase):

    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcess')
        # Store original location
        self.pwd = os.getcwd()
        # Move to working directory
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_analysis_dir_path(self):
        """AutoProcess: analysis dir path is absolute and normalized
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up new AutoProcess instance
        ap = AutoProcess()
        self.assertEqual(ap.analysis_dir,None)
        # Make a mock analysis dir
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Create Autoprocess instances from different
        # forms of path and check stored value
        rel_path = "160621_M00879_0087_000000000-AGEW9_analysis"
        abs_path = os.path.join(self.dirn,rel_path)
        rel_unnormalised = os.path.join("..",
                                        os.path.basename(self.dirn),
                                        rel_path)
        abs_unnormalised = os.path.join(self.dirn,
                                        rel_unnormalised)
        ap = AutoProcess(analysis_dir=abs_path)
        self.assertEqual(ap.analysis_dir,abs_path)
        ap = AutoProcess(analysis_dir=rel_path)
        self.assertEqual(ap.analysis_dir,abs_path)
        ap = AutoProcess(analysis_dir=abs_unnormalised)
        self.assertEqual(ap.analysis_dir,abs_path)
        ap = AutoProcess(analysis_dir=rel_unnormalised)
        self.assertEqual(ap.analysis_dir,abs_path)

    def test_set_log_dir(self):
        """AutoProcess.set_log_dir: sets correct log directory
        """
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Load into AutoProcess object
        ap = AutoProcess(analysis_dir=mockdir.dirn)
        # Check the initial log directory
        self.assertEqual(ap.log_dir,
                         os.path.join(mockdir.dirn,
                                      "logs"))
        # Set a log subdir
        subdir = ap.get_log_subdir("testing")
        self.assertEqual(subdir,"001_testing")
        ap.set_log_dir("001_testing")
        self.assertEqual(ap.log_dir,
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "001_testing"))
        self.assertEqual(ap.log_path("test.log"),
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "001_testing",
                                      "test.log"))
        self.assertEqual(ap.log_path("test_1",
                                     "test1.log"),
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "001_testing",
                                      "test_1",
                                      "test1.log"))
        # Set a second log subdir
        subdir = ap.get_log_subdir("testing")
        self.assertEqual(subdir,"002_testing")
        ap.set_log_dir("002_testing")
        self.assertEqual(ap.log_dir,
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "002_testing"))
        self.assertEqual(ap.log_path("test.log"),
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "002_testing",
                                      "test.log"))
        self.assertEqual(ap.log_path("test_2",
                                     "test2.log"),
                         os.path.join(mockdir.dirn,
                                      "logs",
                                      "002_testing",
                                      "test_2",
                                      "test2.log"))

class TestAutoProcessMakeProjectMetadataFile(unittest.TestCase):
    """
    Test for the 'make_project_metadata_file' method
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcess')
        # Store original location
        self.pwd = os.getcwd()
        # Move to working directory
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_make_project_metadata_file(self):
        """
        AutoProcess.make_project_metadata_file: new 'projects.info' from bcl2fastq output
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the projects.info file
        os.remove(os.path.join(mockdir.dirn,"projects.info"))
        # Create a new projects.info file
        AutoProcess(mockdir.dirn).make_project_metadata_file()
        # Check outputs
        self.assertTrue(os.path.exists(
            os.path.join(mockdir.dirn,"projects.info")))
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
""")

    def test_make_project_metadata_file_no_bcl2fastq_output(self):
        """
        AutoProcess.make_project_metadata_file: new 'projects.info' (no bcl2fastq output)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the projects.info file and the bcl2fastq output dir
        os.remove(os.path.join(mockdir.dirn,"projects.info"))
        shutil.rmtree(os.path.join(mockdir.dirn,"bcl2fastq"))
        # Create a new projects.info file
        AutoProcess(mockdir.dirn).make_project_metadata_file()
        # Check outputs
        self.assertTrue(os.path.exists(
            os.path.join(mockdir.dirn,"projects.info")))
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
""")

    def test_make_project_metadata_file_exception_if_file_exists(self):
        """
        AutoProcess.make_project_metadata_file: raise exception if 'projects.info' already exists
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Raise exception if projects.info file exists already
        self.assertRaises(Exception,
                          AutoProcess(mockdir.dirn).make_project_metadata_file)

class TestAutoProcessUpdateProjectMetadataFile(unittest.TestCase):
    """
    Test for the 'update_project_metadata_file' method
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcess')
        # Store original location
        self.pwd = os.getcwd()
        # Move to working directory
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_update_project_metadata_file_empty_from_bcl2fastq_output(self):
        """
        AutoProcess.update_project_metadata_file: update 'empty' file from bcl2fastq output
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create empty projects.info file
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\n")
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
""")

    def test_update_project_metadata_file_partial_from_bcl2fastq_output(self):
        """
        AutoProcess.update_project_metadata_file: update partial file from bcl2fastq output
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create projects.info file with one project already listed
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\nCDE\tCDE3,CDE4\t.\t.\t.\t.\t.\tKeep me")
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\tKeep me
""")

    def test_update_project_metadata_file_missing_from_bcl2fastq_output(self):
        """
        AutoProcess.update_project_metadata_file: make missing file and populate from bcl2fastq output
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove projects.info file
        os.remove(os.path.join(mockdir.dirn,"projects.info"))
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
""")

    def test_update_project_metadata_file_comment_out_missing_project(self):
        """
        AutoProcess.update_project_metadata_file: missing project is commented out
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create projects.info file with one project already listed
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\nFG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me")
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output - missing project kept but commented out
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            print(fp.read())
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
#FG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me
""")

    def test_update_project_metadata_file_uncomment_existing_project(self):
        """
        AutoProcess.update_project_metadata_file: existing project is uncommented
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create projects.info file with one project already listed
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\n#CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\tKeep me")
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output - missing project kept but commented out
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\tKeep me
""")

    def test_update_project_metadata_file_dont_comment_missing_project_when_dir_is_present(self):
        """
        AutoProcess.update_project_metadata_file: don't comment out 'missing' project when dir is present
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create projects.info file with one project already listed
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\nFG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me")
        # Create the corresponding project
        project = MockAnalysisProject('FG',('FG5_S1_R1_001.fastq.gz',
                                            'FG6_S1_R1_001.fastq.gz'))
        project.create(top_dir=mockdir.dirn)
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output - missing project kept but commented out
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            print(fp.read())
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
FG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me
""")

    def test_update_project_metadata_file_dont_uncomment_missing_project_when_dir_is_present(self):
        """
        AutoProcess.update_project_metadata_file: don't uncomment 'missing' project when dir is present
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Create projects.info file with one project already listed
        with open(os.path.join(mockdir.dirn,"projects.info"),'wt') as fp:
            fp.write("#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments\n#FG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me")
        # Create the corresponding project
        project = MockAnalysisProject('FG',('FG5_S1_R1_001.fastq.gz',
                                            'FG6_S1_R1_001.fastq.gz'))
        project.create(top_dir=mockdir.dirn)
        # Update the projects.info file
        AutoProcess(mockdir.dirn).update_project_metadata_file()
        # Check output - missing project kept but commented out
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            print(fp.read())
        with open(os.path.join(mockdir.dirn,"projects.info"),'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
#FG\tFG5,FG6\t.\t.\t.\t.\t.\tKeep me
""")

class TestAutoProcessGetAnalysisProjectsMethod(unittest.TestCase):
    """
    Tests for the 'get_analysis_projects' method
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcess')
        # Store original location
        self.pwd = os.getcwd()
        # Move to working directory
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_no_projects_dot_info_no_project_dirs(self):
        """AutoProcess.get_analysis_projects: no project dirs (no projects.info)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the projects.info file
        os.remove(os.path.join(mockdir.dirn,"projects.info"))
        print(os.listdir(mockdir.dirn))
        # No projects should be listed
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        self.assertEqual(projects,[])

    def test_projects_dot_info_no_project_dirs(self):
        """AutoProcess.get_analysis_projects: no project dirs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # No projects should be listed
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        self.assertEqual(projects,[])

    def test_with_project_dirs(self):
        """AutoProcess.get_analysis_projects: project dirs exist
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        expected = ('AB','CDE','undetermined')
        self.assertEqual(len(projects),len(expected))
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

    def test_with_project_dirs_no_projects_dot_info(self):
        """AutoProcess.get_analysis_projects: project dirs exist (no projects.info)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # Remove the projects.info file
        os.remove(os.path.join(mockdir.dirn,"projects.info"))
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        expected = ('undetermined',)
        self.assertEqual(len(projects),1)
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

    def test_with_project_dirs_select_subset(self):
        """AutoProcess.get_analysis_projects: select subset of projects
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects("C*")
        expected = ('CDE',)
        self.assertEqual(len(projects),len(expected))
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

class TestAutoProcessGetAnalysisProjectsFromDirsMethod(unittest.TestCase):
    """
    Tests for the 'get_analysis_projects_from_dirs' method
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcess')
        # Store original location
        self.pwd = os.getcwd()
        # Move to working directory
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_no_projects(self):
        """AutoProcess.get_analysis_projects_from_dirs: no project dirs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # No projects should be listed
        projects = AutoProcess(mockdir.dirn).get_analysis_projects_from_dirs()
        self.assertEqual(projects,[])

    def test_with_projects(self):
        """AutoProcess.get_analysis_projects_from_dirs: fetches project dirs
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        expected = ('AB','CDE','undetermined')
        self.assertEqual(len(projects),len(expected))
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

    def test_with_projects_select_subset(self):
        """AutoProcess.get_analysis_projects_from_dirs: selects subset of projects
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            top_dir=self.dirn)
        mockdir.create()
        # List the projects
        projects = AutoProcess(mockdir.dirn).get_analysis_projects_from_dirs("C*")
        expected = ('CDE',)
        self.assertEqual(len(projects),len(expected))
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

class TestAutoProcessPairedEndMethod(unittest.TestCase):
    """
    Tests for the 'paired_end' method
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcess')
        # Store original location
        self.pwd = os.getcwd()
        # Move to working directory
        os.chdir(self.dirn)

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_single_end(self):
        """AutoProcess.paired_end: check for single-ended data
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=False,
            top_dir=self.dirn)
        mockdir.create()
        self.assertFalse(AutoProcess(mockdir.dirn).paired_end)

    def test_paired_end(self):
        """AutoProcess.paired_end: check for paired-ended data
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=True,
            top_dir=self.dirn)
        mockdir.create()
        self.assertTrue(AutoProcess(mockdir.dirn).paired_end)

    def test_single_end_no_projects(self):
        """AutoProcess.paired_end: check for single-ended data (no project dirs)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=False,
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        self.assertFalse(AutoProcess(mockdir.dirn).paired_end)

    def test_paired_end_no_projects(self):
        """AutoProcess.paired_end: check for paired-ended data (no project dirs)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=True,
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        self.assertTrue(AutoProcess(mockdir.dirn).paired_end)

    def test_single_end_no_projects_no_aligned(self):
        """AutoProcess.paired_end: check for single-ended data (no project dirs, no unaligned)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=False,
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the unaligned dir
        shutil.rmtree(os.path.join(mockdir.dirn,"bcl2fastq"))
        self.assertEqual(AutoProcess(mockdir.dirn).paired_end,None)

    def test_paired_end_no_projects_no_unaligned(self):
        """AutoProcess.paired_end: check for paired-ended data (no project dirs, no unaligned)
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_K00879_0087_000000000-AGEW9',
            'hiseq',
            metadata={ "run_number": 87,
                       "source": "local" },
            paired_end=True,
            top_dir=self.dirn)
        mockdir.create(no_project_dirs=True)
        # Remove the unaligned dir
        shutil.rmtree(os.path.join(mockdir.dirn,"bcl2fastq"))
        self.assertEqual(AutoProcess(mockdir.dirn).paired_end,None)
