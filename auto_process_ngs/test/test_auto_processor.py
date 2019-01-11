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

    def test_no_projects_dot_info_no_project_dirs_no_unaligned(self):
        """AutoProcess.get_analysis_projects: no project dirs (no projects.info, no unaligned)
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
        # Remove the unaligned dir
        shutil.rmtree(os.path.join(mockdir.dirn,"bcl2fastq"))
        print os.listdir(mockdir.dirn)
        # No projects should be listed
        projects = AutoProcess(mockdir.dirn).get_analysis_projects()
        self.assertEqual(projects,[])

    def test_no_projects_dot_info_no_project_dirs(self):
        """AutoProcess.get_analysis_projects: no project dirs  (no projects.info)
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
        print os.listdir(mockdir.dirn)
        # Listing the projects should raise an exception
        self.assertRaises(Exception,
                          AutoProcess(mockdir.dirn).get_analysis_projects)

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
        # Listing the projects should raise an exception
        self.assertRaises(Exception,
                          AutoProcess(mockdir.dirn).get_analysis_projects)

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
        expected = ('AB','CDE','undetermined')
        self.assertEqual(len(projects),len(expected))
        for p in projects:
            self.assertTrue(isinstance(p,AnalysisProject))
            self.assertTrue(p.name in expected)
        for p in expected:
            matched_projects = [x for x in projects if x.name == p]
            self.assertEqual(len(matched_projects),1)

    def test_with_project_dirs_no_projects_dot_info_no_unaligned(self):
        """AutoProcess.get_analysis_projects: project dirs exist (no projects.info, no unaligned)
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
        expected = ('AB','CDE','undetermined')
        self.assertEqual(len(projects),len(expected))
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

class TestAutoProcessSetup(unittest.TestCase):

    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessSetup')
        # Store original location so we can get back at the end
        self.pwd = os.getcwd()
        # Move to working dir
        os.chdir(self.dirn)
        # Placeholders for test objects
        self.mock_illumina_run = None
        self.analysis_dir = None

    def tearDown(self):
        # Return to original dir
        os.chdir(self.pwd)
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_autoprocess_setup(self):
        """AutoProcess.setup works for mock MISeq run
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '151125_M00879_0001_000000000-ABCDE1',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up autoprocessor
        ap = AutoProcess()
        ap.setup(mock_illumina_run.dirn)
        analysis_dirn = "%s_analysis" % mock_illumina_run.name
        # Check parameters
        self.assertEqual(ap.analysis_dir,
                         os.path.join(self.dirn,analysis_dirn))
        self.assertEqual(ap.params.data_dir,mock_illumina_run.dirn)
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(self.dirn,analysis_dirn,
                                      'custom_SampleSheet.csv'))
        self.assertEqual(ap.params.bases_mask,'auto')
        # Check metadata
        self.assertEqual(ap.metadata.run_name,
                         "151125_M00879_0001_000000000-ABCDE1")
        self.assertEqual(ap.metadata.run_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"miseq")
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.assay,"TruSeq HT")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.instrument_name,"M00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"151125")
        self.assertEqual(ap.metadata.instrument_run_number,"1")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        # Delete to force write of data to disk
        del(ap)
        # Check directory exists
        self.assertTrue(os.path.isdir(analysis_dirn))
        # Check files exists
        for filen in ('SampleSheet.orig.csv',
                      'custom_SampleSheet.csv',
                      'auto_process.info',
                      'metadata.info',):
            self.assertTrue(os.path.exists(os.path.join(analysis_dirn,
                                                        filen)),
                            "Missing file: %s" % filen)
        # Check subdirs have been created
        for subdirn in ('ScriptCode',
                        'logs',):
            self.assertTrue(os.path.isdir(os.path.join(analysis_dirn,
                                                       subdirn)),
                            "Missing subdir: %s" % subdirn)

    def test_autoprocess_setup_non_canonical_run_name(self):
        """AutoProcess.setup handles run name with non-canonical run name
        """
        # Create mock Illumina run directory with missing flow cell ID
        mock_illumina_run = MockIlluminaRun(
            '151125_M00879_0001',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up autoprocessor
        ap = AutoProcess()
        ap.setup(mock_illumina_run.dirn)
        analysis_dirn = "%s_analysis" % mock_illumina_run.name
        # Check parameters
        self.assertEqual(ap.analysis_dir,
                         os.path.join(self.dirn,analysis_dirn))
        self.assertEqual(ap.params.data_dir,mock_illumina_run.dirn)
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(self.dirn,analysis_dirn,
                                      'custom_SampleSheet.csv'))
        self.assertEqual(ap.params.bases_mask,'auto')
        # Check metadata
        self.assertEqual(ap.metadata.run_name,
                         "151125_M00879_0001")
        self.assertEqual(ap.metadata.run_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"miseq")
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.assay,"TruSeq HT")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.instrument_name,None)
        self.assertEqual(ap.metadata.instrument_datestamp,None)
        self.assertEqual(ap.metadata.instrument_run_number,None)
        self.assertEqual(ap.metadata.instrument_flow_cell_id,None)

    def test_autoprocess_setup_existing_target_dir(self):
        """AutoProcess.setup works when target dir exists
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Make a mock auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Do setup into existing analysis dir
        ap = AutoProcess()
        ap.setup(mock_illumina_run.dirn)
        self.assertTrue(os.path.isdir(
            '160621_M00879_0087_000000000-AGEW9'))

    def test_autoprocess_setup_absolute_paths(self):
        """AutoProcess.setup data dir path is absolute and normalised
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        data_dir_rel = '160621_M00879_0087_000000000-AGEW9'
        data_dir_abs = os.path.join(self.dirn,data_dir_rel)
        data_dir_rel_unnormalised = os.path.join(
            "..",
            os.path.basename(self.dirn),
            data_dir_rel)
        data_dir_abs_unnormalised = os.path.join(
            self.dirn,
            "..",
            os.path.basename(self.dirn),
            data_dir_rel)
        analysis_dir = os.path.join(
            self.dirn,
            "160621_M00879_0087_000000000-AGEW9_analysis")
        # Do setup using absolute path
        ap = AutoProcess()
        ap.setup(data_dir_abs)
        self.assertEqual(ap.params.data_dir,data_dir_abs)
        self.assertEqual(ap.analysis_dir,analysis_dir)
        del(ap)
        shutil.rmtree(analysis_dir)
        # Do setup using relative path
        ap = AutoProcess()
        ap.setup(data_dir_rel)
        self.assertEqual(ap.params.data_dir,data_dir_abs)
        self.assertEqual(ap.analysis_dir,analysis_dir)
        del(ap)
        shutil.rmtree(analysis_dir)
        # Do setup using absolute unnormalized path
        ap = AutoProcess()
        ap.setup(data_dir_abs_unnormalised)
        self.assertEqual(ap.params.data_dir,data_dir_abs)
        self.assertEqual(ap.analysis_dir,analysis_dir)
        del(ap)

    def test_autoprocess_setup_missing_data_directory(self):
        """AutoProcess.setup raises exception if data directory is missing
        """
        # Set up autoprocessor
        ap = AutoProcess()
        self.assertRaises(Exception,
                          ap.setup,
                          os.path.join(
                              self.dirn,
                              '160621_M00879_0087_000000000-AGEW9'))
        self.assertFalse(os.path.exists(
            os.path.join(
                self.dirn,
                '160621_M00879_0087_000000000-AGEW9_analysis')))

    def test_autoprocess_setup_missing_sample_sheet(self):
        """AutoProcess.setup raises exception if sample sheet not found
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '160621_NB00879_0087_000000000-AGEW9',
            'nextseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up autoprocessor
        ap = AutoProcess()
        self.assertRaises(Exception,
                          ap.setup,
                          os.path.join(
                              self.dirn,
                              '160621_NB00879_0087_000000000-AGEW9'))
        self.assertFalse(os.path.exists(
            os.path.join(
                self.dirn,
                '160621_NB00879_0087_000000000-AGEW9_analysis')))

    def test_autoprocess_setup_samplesheet_from_url(self):
        """AutoProcess.setup works when samplesheet is a URL
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '151125_M00879_0001_000000000-ABCDE1',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Copy samplesheet
        sample_sheet = os.path.join(self.dirn,'samplesheet.csv')
        with open(os.path.join(mock_illumina_run.dirn,
                               'Data','Intensities','BaseCalls',
                               'SampleSheet.csv'),'r') as fp1:
            with open(sample_sheet,'w') as fp2:
                fp2.write(fp1.read())
        sample_sheet = "file://%s" % sample_sheet
        print sample_sheet
        # Set up autoprocessor
        ap = AutoProcess()
        ap.setup(mock_illumina_run.dirn,sample_sheet=sample_sheet)
        analysis_dirn = "%s_analysis" % mock_illumina_run.name
        # Check parameters
        self.assertEqual(ap.analysis_dir,
                         os.path.join(self.dirn,analysis_dirn))
        self.assertEqual(ap.params.data_dir,mock_illumina_run.dirn)
        self.assertEqual(ap.params.sample_sheet,
                         os.path.join(self.dirn,analysis_dirn,
                                      'custom_SampleSheet.csv'))
        self.assertEqual(ap.params.bases_mask,'auto')
        # Check metadata
        self.assertEqual(ap.metadata.run_name,
                         "151125_M00879_0001_000000000-ABCDE1")
        self.assertEqual(ap.metadata.run_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"miseq")
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.assay,"TruSeq HT")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.instrument_name,"M00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"151125")
        self.assertEqual(ap.metadata.instrument_run_number,"1")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        # Delete to force write of data to disk
        del(ap)
        # Check directory exists
        self.assertTrue(os.path.isdir(analysis_dirn))
        # Check files exists
        for filen in ('SampleSheet.orig.csv',
                      'custom_SampleSheet.csv',
                      'auto_process.info',
                      'metadata.info',):
            self.assertTrue(os.path.exists(os.path.join(analysis_dirn,
                                                        filen)),
                            "Missing file: %s" % filen)
        # Check subdirs have been created
        for subdirn in ('ScriptCode',
                        'logs',):
            self.assertTrue(os.path.isdir(os.path.join(analysis_dirn,
                                                       subdirn)),
                            "Missing subdir: %s" % subdirn)
