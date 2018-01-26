#######################################################################
# Tests for autoprocessor.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
import gzip
from bcftbx.mock import MockIlluminaRun
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.utils import AnalysisProject
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import MockAnalysisDir
from auto_process_ngs.mock import UpdateAnalysisDir
from auto_process_ngs.mock import UpdateAnalysisProject

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
        self.assertEqual(ap.params.bases_mask,
                         'y101,I8,I8,y101')
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

fastq_reads_r1 = (
    "@HISEQ:1:000000000-A2Y1L:1:1101:19264:2433 1:N:0:AGATCGC",
    "AGATAGCCGA","+","?????BBB@B",
    "@HISEQ:1:000000000-A2Y1L:1:1101:18667:2435 1:N:0:AGATCGC",
    "ATATATTCAT","+","?????BBBDD",
    "@HISEQ:1:000000000-A2Y1L:1:1101:17523:2436 1:N:0:AGATCGC",
    "CATCACTACC","+","?<,<?BBBBB"
)
fastq_reads_r2 = (
    "@HISEQ:1:000000000-A2Y1L:1:1101:19264:2433 2:N:0:AGATCGC",
    "GCCGATATGC","+","??A??ABBDD",
    "@HISEQ:1:000000000-A2Y1L:1:1101:18667:2435 2:N:0:AGATCGC",
    "GATGACATCA","+","?????BBBDD",
    "@HISEQ:1:000000000-A2Y1L:1:1101:17523:2436 2:N:0:AGATCGC",
    "GAATATAGAA","+","??AAABBBDD"
)

class TestAutoProcessMergeFastqDirs(unittest.TestCase):
    """
    Tests for AutoProcess.merge_fastq_dirs
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessMergeFastqDirs')
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
        shutil.rmtree(self.dirn)

    def _setup_bcl2fastq2_no_lane_splitting(self):
        # Create mock bcl2fastq2 dir structure with no lane splitting
        lanes = [1,2,3,4]
        mockdir1 = MockAnalysisDir("161209_K123_0001_BLAH",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.AB',
                                   fmt='bcl2fastq2',
                                   paired_end=True,
                                   no_lane_splitting=True,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir1.add_fastq_batch('AB','AB1','AB1_S1',lanes=lanes)
        mockdir1.add_fastq_batch('AB','AB2','AB2_S2',lanes=lanes)
        m1 = mockdir1.create()
        print m1
        mockdir2 = MockAnalysisDir("161209_K123_0001_BLAH2",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.CDE',
                                   fmt='bcl2fastq2',
                                   paired_end=True,
                                   no_lane_splitting=True,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir2.add_fastq_batch('CDE','CDE3','CDE3_S3',lanes=lanes)
        mockdir2.add_fastq_batch('CDE','CDE4','CDE4_S4',lanes=lanes)
        m2 = mockdir2.create()
        # Move the second 'unaligned' dir into the first analysis dir
        shutil.move(os.path.join(m2,'bcl2fastq.CDE'),m1)
        # Remove unwanted project dirs and files
        shutil.rmtree(os.path.join(m1,'AB'))
        shutil.rmtree(os.path.join(m1,'undetermined'))
        # Add content to the undetermined fastqs
        for n,dirn in enumerate(('bcl2fastq.AB','bcl2fastq.CDE')):
            undetermined_r1 = os.path.join(m1,dirn,'Undetermined_S0_R1_001.fastq.gz')
            undetermined_r2 = os.path.join(m1,dirn,'Undetermined_S0_R2_001.fastq.gz')
            with gzip.GzipFile(undetermined_r1,'wb') as fq:
                for i in xrange(4):
                    fq.write("%s\n" % fastq_reads_r1[n*4+i])
            with gzip.GzipFile(undetermined_r2,'wb') as fq:
                for i in xrange(4):
                    fq.write("%s\n" % fastq_reads_r2[n*4+i])
        print m2
        return m1

    def _setup_bcl2fastq2(self):
        # Create a mock bcl2fastq dir structure
        lanes = [1,2,3,4]
        mockdir1 = MockAnalysisDir("161209_K123_0002_BLAH",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.lanes1-2',
                                   fmt='bcl2fastq2',
                                   paired_end=True,
                                   no_lane_splitting=False,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir1.add_fastq_batch('AB','AB1','AB1_S1',lanes=[1,2])
        mockdir1.add_fastq_batch('AB','AB2','AB2_S2',lanes=[1,2])
        m1 = mockdir1.create()
        print m1
        mockdir2 = MockAnalysisDir("161209_K123_0002_BLAH2",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.lanes3-4',
                                   fmt='bcl2fastq2',
                                   paired_end=True,
                                   no_lane_splitting=False,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir2.add_fastq_batch('CDE','CDE3','CDE3_S3',lanes=[3,4])
        mockdir2.add_fastq_batch('CDE','CDE4','CDE4_S4',lanes=[3,4])
        m2 = mockdir2.create()
        # Move the second 'unaligned' dir into the first analysis dir
        shutil.move(os.path.join(m2,'bcl2fastq.lanes3-4'),m1)
        # Remove unwanted project dirs and files
        shutil.rmtree(os.path.join(m1,'AB'))
        shutil.rmtree(os.path.join(m1,'undetermined'))
        print m2
        return m1

    def _setup_casava(self):
        # Create a mock casava dir structure
        lanes = [1,2,3,4]
        mockdir1 = MockAnalysisDir("161209_K123_0003_BLAH",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.lanes1-2',
                                   fmt='casava',
                                   paired_end=True,
                                   no_lane_splitting=False,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir1.add_fastq_batch('AB','AB1','AB1_GCCAAT',lanes=[1,2])
        mockdir1.add_fastq_batch('AB','AB2','AB2_AGTCAA',lanes=[1,2])
        m1 = mockdir1.create()
        print m1
        mockdir2 = MockAnalysisDir("161209_K123_0003_BLAH2",
                                   "hiseq",
                                   unaligned_dir='bcl2fastq.lanes3-4',
                                   fmt='casava',
                                   paired_end=True,
                                   no_lane_splitting=False,
                                   lanes=lanes,
                                   top_dir=self.dirn)
        mockdir2.add_fastq_batch('CDE','CDE3','CDE3_GCCAAT',lanes=[3,4])
        mockdir2.add_fastq_batch('CDE','CDE4','CDE4_AGTCAA',lanes=[3,4])
        m2 = mockdir2.create()
        # Move the second 'unaligned' dir into the first analysis dir
        shutil.move(os.path.join(m2,'bcl2fastq.lanes3-4'),m1)
        # Remove unwanted project dirs and files
        shutil.rmtree(os.path.join(m1,'AB'))
        shutil.rmtree(os.path.join(m1,'undetermined'))
        print m2
        return m1

    def _assert_dir_exists(self,path):
        self.assertTrue(os.path.isdir(path),
                        "Missing dir '%s'" % path)
        
    def _assert_dir_doesnt_exist(self,path):
        self.assertFalse(os.path.isdir(path),
                         "Dir '%s' shouldn't exist" % path)

    def _assert_file_exists(self,path):
        self.assertTrue(os.path.isfile(path),
                        "Missing file '%s'" % path)

    def _assert_file_doesnt_exist(self,path):
        self.assertFalse(os.path.isfile(path),
                        "File '%s' shouldn't exist" % path)

    def test_bcl2fastq2_no_lane_splitting_new_output_dir(self):
        """
        AutoProcess.merge_fastq_dirs: bcl2fastq v2 output with --no-lane-splitting, new output dir
        """
        analysis_dir = self._setup_bcl2fastq2_no_lane_splitting()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.AB",output_dir="bcl2fastq")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.AB'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.CDE'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.AB'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.CDE'))
        for f in ('AB/AB1_S1_R1_001.fastq.gz',
                  'AB/AB1_S1_R2_001.fastq.gz',
                  'AB/AB2_S2_R1_001.fastq.gz',
                  'AB/AB2_S2_R2_001.fastq.gz',
                  'CDE/CDE3_S3_R1_001.fastq.gz',
                  'CDE/CDE3_S3_R2_001.fastq.gz',
                  'CDE/CDE4_S4_R1_001.fastq.gz',
                  'CDE/CDE4_S4_R2_001.fastq.gz',
                  'Undetermined_S0_R1_001.fastq.gz',
                  'Undetermined_S0_R2_001.fastq.gz',):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq',f))
        # Check merge of undetermined fastqs
        undetermined_r1 = gzip.GzipFile(
            os.path.join(analysis_dir,'bcl2fastq','Undetermined_S0_R1_001.fastq.gz'),
            'rb').read()
        expected_r1 = '\n'.join(fastq_reads_r1[:8])+'\n'
        self.assertEqual(undetermined_r1,expected_r1)
        undetermined_r2 = gzip.GzipFile(
            os.path.join(analysis_dir,'bcl2fastq','Undetermined_S0_R2_001.fastq.gz'),
            'rb').read()
        expected_r2 = '\n'.join(fastq_reads_r2[:8])+'\n'
        self.assertEqual(undetermined_r2,expected_r2)
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Protocol	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2_new_output_dir(self):
        """
        AutoProcess.merge_fastq_dirs: bcl2fastq v2 output, new output dir
        """
        analysis_dir = self._setup_bcl2fastq2()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2",output_dir="bcl2fastq")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
        for f in ('AB/AB1_S1_L001_R1_001.fastq.gz',
                  'AB/AB1_S1_L001_R2_001.fastq.gz',
                  'AB/AB2_S2_L001_R1_001.fastq.gz',
                  'AB/AB2_S2_L001_R2_001.fastq.gz',
                  'AB/AB1_S1_L002_R1_001.fastq.gz',
                  'AB/AB1_S1_L002_R2_001.fastq.gz',
                  'AB/AB2_S2_L002_R1_001.fastq.gz',
                  'AB/AB2_S2_L002_R2_001.fastq.gz',
                  'CDE/CDE3_S3_L003_R1_001.fastq.gz',
                  'CDE/CDE3_S3_L003_R2_001.fastq.gz',
                  'CDE/CDE4_S4_L003_R1_001.fastq.gz',
                  'CDE/CDE4_S4_L003_R2_001.fastq.gz',
                  'CDE/CDE3_S3_L004_R1_001.fastq.gz',
                  'CDE/CDE3_S3_L004_R2_001.fastq.gz',
                  'CDE/CDE4_S4_L004_R1_001.fastq.gz',
                  'CDE/CDE4_S4_L004_R2_001.fastq.gz',
                  'Undetermined_S0_L001_R1_001.fastq.gz',
                  'Undetermined_S0_L001_R2_001.fastq.gz',
                  'Undetermined_S0_L002_R1_001.fastq.gz',
                  'Undetermined_S0_L002_R2_001.fastq.gz',
                  'Undetermined_S0_L003_R1_001.fastq.gz',
                  'Undetermined_S0_L003_R2_001.fastq.gz',
                  'Undetermined_S0_L004_R1_001.fastq.gz',
                  'Undetermined_S0_L004_R2_001.fastq.gz'):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq',f))
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Protocol	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_casava_new_output_dir(self):
        """
        AutoProcess.merge_fastq_dirs: casava/bcl2fastq v1.8.* output, new output dir
        """
        analysis_dir = self._setup_casava()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2",output_dir="bcl2fastq")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
        for f in ('Project_AB/Sample_AB1/AB1_GCCAAT_L001_R1_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L001_R2_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L001_R1_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L001_R2_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R1_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R2_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L002_R1_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L002_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L003_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L003_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L003_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L003_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L004_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L004_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L004_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L004_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R2_001.fastq.gz'):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq',f))
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Protocol	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2_no_lane_splitting(self):
        """
        AutoProcess.merge_fastq_dirs: bcl2fastq v2 output with --no-lane-splitting
        """
        analysis_dir = self._setup_bcl2fastq2_no_lane_splitting()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.AB")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.AB'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.CDE'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.AB'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.CDE'))
        for f in ('AB/AB1_S1_R1_001.fastq.gz',
                  'AB/AB1_S1_R2_001.fastq.gz',
                  'AB/AB2_S2_R1_001.fastq.gz',
                  'AB/AB2_S2_R2_001.fastq.gz',
                  'CDE/CDE3_S3_R1_001.fastq.gz',
                  'CDE/CDE3_S3_R2_001.fastq.gz',
                  'CDE/CDE4_S4_R1_001.fastq.gz',
                  'CDE/CDE4_S4_R2_001.fastq.gz',
                  'Undetermined_S0_R1_001.fastq.gz',
                  'Undetermined_S0_R2_001.fastq.gz',):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq.AB',f))
        # Check merge of undetermined fastqs
        undetermined_r1 = gzip.GzipFile(
            os.path.join(analysis_dir,'bcl2fastq.AB','Undetermined_S0_R1_001.fastq.gz'),
            'rb').read()
        expected_r1 = '\n'.join(fastq_reads_r1[:8])+'\n'
        self.assertEqual(undetermined_r1,expected_r1)
        undetermined_r2 = gzip.GzipFile(
            os.path.join(analysis_dir,'bcl2fastq.AB','Undetermined_S0_R2_001.fastq.gz'),
            'rb').read()
        expected_r2 = '\n'.join(fastq_reads_r2[:8])+'\n'
        self.assertEqual(undetermined_r2,expected_r2)
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Protocol	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2(self):
        """
        AutoProcess.merge_fastq_dirs: bcl2fastq v2 output
        """
        analysis_dir = self._setup_bcl2fastq2()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
        for f in ('AB/AB1_S1_L001_R1_001.fastq.gz',
                  'AB/AB1_S1_L001_R2_001.fastq.gz',
                  'AB/AB2_S2_L001_R1_001.fastq.gz',
                  'AB/AB2_S2_L001_R2_001.fastq.gz',
                  'AB/AB1_S1_L002_R1_001.fastq.gz',
                  'AB/AB1_S1_L002_R2_001.fastq.gz',
                  'AB/AB2_S2_L002_R1_001.fastq.gz',
                  'AB/AB2_S2_L002_R2_001.fastq.gz',
                  'CDE/CDE3_S3_L003_R1_001.fastq.gz',
                  'CDE/CDE3_S3_L003_R2_001.fastq.gz',
                  'CDE/CDE4_S4_L003_R1_001.fastq.gz',
                  'CDE/CDE4_S4_L003_R2_001.fastq.gz',
                  'CDE/CDE3_S3_L004_R1_001.fastq.gz',
                  'CDE/CDE3_S3_L004_R2_001.fastq.gz',
                  'CDE/CDE4_S4_L004_R1_001.fastq.gz',
                  'CDE/CDE4_S4_L004_R2_001.fastq.gz',
                  'Undetermined_S0_L001_R1_001.fastq.gz',
                  'Undetermined_S0_L001_R2_001.fastq.gz',
                  'Undetermined_S0_L002_R1_001.fastq.gz',
                  'Undetermined_S0_L002_R2_001.fastq.gz',
                  'Undetermined_S0_L003_R1_001.fastq.gz',
                  'Undetermined_S0_L003_R2_001.fastq.gz',
                  'Undetermined_S0_L004_R1_001.fastq.gz',
                  'Undetermined_S0_L004_R2_001.fastq.gz'):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2',f))
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Protocol	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_casava(self):
        """
        AutoProcess.merge_fastq_dirs: casava/bcl2fastq v1.8.* output
        """
        analysis_dir = self._setup_casava()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2")
        # Check outputs
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
        for f in ('Project_AB/Sample_AB1/AB1_GCCAAT_L001_R1_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L001_R2_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L001_R1_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L001_R2_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R1_001.fastq.gz',
                  'Project_AB/Sample_AB1/AB1_GCCAAT_L002_R2_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L002_R1_001.fastq.gz',
                  'Project_AB/Sample_AB2/AB2_AGTCAA_L002_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L003_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L003_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L003_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L003_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L004_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE3/CDE3_GCCAAT_L004_R2_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L004_R1_001.fastq.gz',
                  'Project_CDE/Sample_CDE4/CDE4_AGTCAA_L004_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane3/lane3_Undetermined_L003_R2_001.fastq.gz',
                  'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R1_001.fastq.gz',
                  'Undetermined_indices/Sample_lane4/lane4_Undetermined_L004_R2_001.fastq.gz'):
            self._assert_file_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2',f))
        # Check projects.info files
        self._assert_file_exists(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))
        projects_info = open(os.path.join(analysis_dir,'projects.info'),'r').read()
        expected = """#Project	Samples	User	Library	Protocol	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2_no_lane_splitting_dry_run(self):
        """
        AutoProcess.merge_fastq_dirs: dry run on bcl2fastq v2 output with --no-lane-splitting
        """
        analysis_dir = self._setup_bcl2fastq2_no_lane_splitting()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.AB",dry_run=True)
        # Check outputs
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.AB'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.CDE'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.AB'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.CDE'))
        # Check projects.info files
        self._assert_file_doesnt_exist(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))

    def test_bcl2fastq2_dry_run(self):
        """
        AutoProcess.merge_fastq_dirs: dry run on bcl2fastq v2 output
        """
        analysis_dir = self._setup_bcl2fastq2()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2",dry_run=True)
        # Check outputs
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
        # Check projects.info files
        self._assert_file_doesnt_exist(os.path.join(analysis_dir,'save.projects.info'))
        self._assert_file_exists(os.path.join(analysis_dir,'projects.info'))

    def test_casava_dry_run(self):
        """
        AutoProcess.merge_fastq_dirs: dry run on casava/bcl2fastq v1.8.* output
        """
        analysis_dir = self._setup_casava()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        self.ap.merge_fastq_dirs("bcl2fastq.lanes1-2",dry_run=True)
        # Check outputs
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))

class TestAutoProcessImportProject(unittest.TestCase):
    """Tests for AutoProcess.import_project

    """
    def setUp(self):
        self.dirn = tempfile.mkdtemp(suffix='TestAutoProcessImportProject')
        # Make a mock project
        project_dir = os.path.join(self.dirn,'NewProj')
        os.mkdir(project_dir)
        os.mkdir(os.path.join(project_dir,'fastqs'))
        for fq in ('NP01_S1_R1_001.fastq.gz','NP01_S1_R1_001.fastq.gz'):
            open(os.path.join(project_dir,'fastqs',fq),'w').write('')
        open(os.path.join(project_dir,'README.info'),'w').write(
            """Run\t160622_NB5001234_0011_ABCDE5AFXX
Platform\tnextseq
User\tPeter Briggs
PI\tAnne Cleaves
Organism\tHuman
Library type\tRNA-seq
Protocol\tStandard
Paired_end\tY
Samples\t1 sample (NP01)
Comments\t1% PhiX spike in
""")
        self.new_project_dir = project_dir

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_import_project(self):
        """Check AutoProcess.import_project imports a project
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Check that the project is not currently present
        ap = AutoProcess(mockdir.dirn)
        self.assertFalse('NewProj' in [p.name
                                       for p in ap.get_analysis_projects()])
        self.assertFalse('NewProj' in [p.name
                                       for p in ap.get_analysis_projects_from_dirs()])
        self.assertFalse(os.path.exists(os.path.join(ap.analysis_dir,'NewProj')))
        # Import the project
        ap.import_project(self.new_project_dir)
        self.assertTrue('NewProj' in [p.name
                                      for p in ap.get_analysis_projects()])
        self.assertTrue('NewProj' in [p.name
                                      for p in ap.get_analysis_projects_from_dirs()])
        self.assertTrue(os.path.exists(os.path.join(ap.analysis_dir,'NewProj')))
        # Verify via fresh AutoProcess object
        ap2 = AutoProcess(mockdir.dirn)
        self.assertTrue('NewProj' in [p.name
                                      for p in ap2.get_analysis_projects()])
        self.assertTrue('NewProj' in [p.name
                                      for p in ap2.get_analysis_projects_from_dirs()])
        self.assertTrue(os.path.exists(os.path.join(ap2.analysis_dir,'NewProj')))

    def test_import_project_already_in_metadata_file(self):
        """AutoProcess.import_project fails if project exists in projects.info
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Add the project to projects.info
        with open(os.path.join(mockdir.dirn,'projects.info'),'a') as fp:
            fp.write('%s\n' % '\t'.join(('NewProj','NP01',
                                         '.','.','.','.','.')))
        # Import the project
        ap = AutoProcess(mockdir.dirn)
        self.assertRaises(Exception,
                          ap.import_project,self.new_project_dir)


    def test_import_project_directory_already_exists(self):
        """AutoProcess.import_project fails if directory already exists
        """
        # Make an auto-process directory
        mockdir = MockAnalysisDirFactory.bcl2fastq2(
            '160621_M00879_0087_000000000-AGEW9',
            'miseq',
            top_dir=self.dirn)
        mockdir.create()
        # Make an existing subdirectory with same name as target project
        os.mkdir(os.path.join(mockdir.dirn,'NewProj'))
        # Import the project
        ap = AutoProcess(mockdir.dirn)
        self.assertRaises(Exception,
                          ap.import_project,self.new_project_dir)
