#######################################################################
# Tests for clone_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import UpdateAnalysisDir
from auto_process_ngs.metadata import AnalysisDirParameters
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
        ap = AutoProcess(analysis_dir.dirn)
        UpdateAnalysisDir(ap).add_processing_report()
        ap.add_directory("primary_data/190116_M01234_0002_AXYZ123")
        # Make a copy
        clone_dir = os.path.join(self.dirn,"190116_M01234_0002_AXYZ123_copy")
        self.assertFalse(os.path.exists(clone_dir))
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
                      'statistics.info',
                      'statistics_full.info',
                      'per_lane_statistics.info',
                      'per_lane_sample_stats.info',
                      'processing_qc.html',):
            f = os.path.join(clone_dir,filen)
            self.assertTrue(os.path.isfile(f),"Missing '%s'" % filen)
        # Check unaligned
        unaligned = os.path.join(clone_dir,'bcl2fastq')
        self.assertTrue(os.path.islink(unaligned))
        # Check primary data
        primary_data = os.path.join(clone_dir,
                                    'primary_data',
                                    '190116_M01234_0002_AXYZ123')
        self.assertTrue(os.path.islink(primary_data))
        # Check projects
        for proj in ('AB','CDE','undetermined'):
            d = os.path.join(clone_dir,proj)
            self.assertTrue(os.path.isdir(d),"Missing '%s'" % proj)
        # Check parameters
        params = AnalysisDirParameters(filen=os.path.join(
            clone_dir,
            'auto_process.info'))
        self.assertEqual(params.sample_sheet,
                         os.path.join(clone_dir,"custom_SampleSheet.csv"))
        self.assertEqual(params.primary_data_dir,
                         os.path.join(clone_dir,"primary_data"))

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
        ap = AutoProcess(analysis_dir.dirn)
        UpdateAnalysisDir(ap).add_processing_report()
        ap.add_directory("primary_data/190116_M01234_0002_AXYZ123")
        # Make a copy
        clone_dir = os.path.join(self.dirn,"190116_M01234_0002_AXYZ123_copy")
        self.assertFalse(os.path.exists(clone_dir))
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
                      'statistics.info',
                      'statistics_full.info',
                      'per_lane_statistics.info',
                      'per_lane_sample_stats.info',
                      'processing_qc.html',):
            f = os.path.join(clone_dir,filen)
            self.assertTrue(os.path.isfile(f),"Missing '%s'" % filen)
        # Check unaligned
        unaligned = os.path.join(clone_dir,'bcl2fastq')
        self.assertTrue(os.path.isdir(unaligned))
        # Check primary data
        primary_data = os.path.join(clone_dir,
                                    'primary_data',
                                    '190116_M01234_0002_AXYZ123')
        self.assertTrue(os.path.islink(primary_data))
        # Check projects
        for proj in ('AB','CDE','undetermined'):
            d = os.path.join(clone_dir,proj)
            self.assertTrue(os.path.isdir(d),"Missing '%s'" % proj)
        # Check parameters
        params = AnalysisDirParameters(filen=os.path.join(
            clone_dir,
            'auto_process.info'))
        self.assertEqual(params.sample_sheet,
                         os.path.join(clone_dir,"custom_SampleSheet.csv"))
        self.assertEqual(params.primary_data_dir,
                         os.path.join(clone_dir,"primary_data"))

    def test_clone_analysis_dir_no_projects(self):
        """
        clone: copies an analysis directory excluding projects
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
        ap = AutoProcess(analysis_dir.dirn)
        UpdateAnalysisDir(ap).add_processing_report()
        ap.add_directory("primary_data/190116_M01234_0002_AXYZ123")
        # Make a copy
        clone_dir = os.path.join(self.dirn,"190116_M01234_0002_AXYZ123_copy")
        self.assertFalse(os.path.exists(clone_dir))
        clone(ap,clone_dir,exclude_projects=True)
        self.assertTrue(os.path.isdir(clone_dir))
        # Check contents
        for subdir in ('logs','ScriptCode'):
            d = os.path.join(clone_dir,subdir)
            self.assertTrue(os.path.isdir(d),"Missing '%s'" % subdir)
        for filen in ('SampleSheet.orig.csv',
                      'custom_SampleSheet.csv',
                      'auto_process.info',
                      'metadata.info',
                      'statistics.info',
                      'statistics_full.info',
                      'per_lane_statistics.info',
                      'per_lane_sample_stats.info',
                      'processing_qc.html',):
            f = os.path.join(clone_dir,filen)
            self.assertTrue(os.path.isfile(f),"Missing '%s'" % filen)
        # Check unaligned
        unaligned = os.path.join(clone_dir,'bcl2fastq')
        self.assertTrue(os.path.islink(unaligned))
        # Check primary data
        primary_data = os.path.join(clone_dir,
                                    'primary_data',
                                    '190116_M01234_0002_AXYZ123')
        self.assertTrue(os.path.islink(primary_data))
        # Check projects
        for proj in ('AB','CDE','undetermined'):
            d = os.path.join(clone_dir,proj)
            self.assertFalse(os.path.exists(d),"Found '%s'" % proj)
        # Check parameters
        params = AnalysisDirParameters(filen=os.path.join(
            clone_dir,
            'auto_process.info'))
        self.assertEqual(params.sample_sheet,
                         os.path.join(clone_dir,"custom_SampleSheet.csv"))
        self.assertEqual(params.primary_data_dir,
                         os.path.join(clone_dir,"primary_data"))

    def test_clone_analysis_dir_empty_params(self):
        """
        clone: copies an analysis directory when parameter file is empty
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
        ap = AutoProcess(analysis_dir.dirn)
        UpdateAnalysisDir(ap).add_processing_report()
        ap.add_directory("primary_data/190116_M01234_0002_AXYZ123")
        # Remove data from parameter file
        parameter_file = ap.parameter_file
        tmp_parameter_file = os.path.join(self.dirn,'new_params.tmp')
        del(ap)
        with open(parameter_file,'r') as fp:
            with open(tmp_parameter_file,'w') as fpp:
                for line in fp:
                    line = "%s\t." % line.split('\t')[0]
                    fpp.write(line)
        os.remove(parameter_file)
        os.rename(tmp_parameter_file,parameter_file)
        ap = AutoProcess(analysis_dir.dirn)
        # Make a copy
        clone_dir = os.path.join(self.dirn,"190116_M01234_0002_AXYZ123_copy")
        self.assertFalse(os.path.exists(clone_dir))
        clone(ap,clone_dir,exclude_projects=False)
        self.assertTrue(os.path.isdir(clone_dir))
        # Check contents
        for subdir in ('logs','ScriptCode'):
            d = os.path.join(clone_dir,subdir)
            self.assertTrue(os.path.isdir(d),"Missing '%s'" % subdir)
        for filen in ('SampleSheet.orig.csv',
                      'custom_SampleSheet.csv',
                      'auto_process.info',
                      'metadata.info',
                      'statistics.info',
                      'statistics_full.info',
                      'per_lane_statistics.info',
                      'per_lane_sample_stats.info',
                      'processing_qc.html',):
            f = os.path.join(clone_dir,filen)
            self.assertTrue(os.path.isfile(f),"Missing '%s'" % filen)
        # Check unaligned
        unaligned = os.path.join(clone_dir,'bcl2fastq')
        self.assertTrue(os.path.islink(unaligned))
        # Check primary data
        primary_data = os.path.join(clone_dir,
                                    'primary_data',
                                    '190116_M01234_0002_AXYZ123')
        self.assertFalse(os.path.exists(primary_data))
        # Check projects
        for proj in ('AB','CDE','undetermined'):
            d = os.path.join(clone_dir,proj)
            self.assertTrue(os.path.exists(d),"Missing '%s'" % proj)

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
        ap = AutoProcess(analysis_dir.dirn)
        UpdateAnalysisDir(ap).add_processing_report()
        ap.add_directory("primary_data/190116_M01234_0002_AXYZ123")
        # Make target dir
        clone_dir = os.path.join(self.dirn,"190116_M01234_0002_AXYZ123_copy")
        os.mkdir(clone_dir)
        # Try to copy source dir
        self.assertRaises(Exception,
                          clone,
                          ap,clone_dir)

        
