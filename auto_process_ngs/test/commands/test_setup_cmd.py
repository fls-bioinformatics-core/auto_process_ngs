#######################################################################
# Tests for setup_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from bcftbx.mock import MockIlluminaRun
from auto_process_ngs.commands.setup_cmd import setup as setup_

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestAutoProcessSetup(unittest.TestCase):
    """
    Tests for AutoProcess.setup
    """
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
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.dirn)

    def test_autoprocess_setup(self):
        """setup: works for mock MISeq run
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '151125_M00879_0001_000000000-ABCDE1',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up autoprocessor
        ap = AutoProcess()
        setup_(ap,mock_illumina_run.dirn)
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
        """setup: handle run name with non-canonical run name
        """
        # Create mock Illumina run directory with missing flow cell ID
        mock_illumina_run = MockIlluminaRun(
            '151125_M00879_0001',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up autoprocessor
        ap = AutoProcess()
        setup_(ap,mock_illumina_run.dirn)
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
        """setup: works when target dir exists
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
        setup_(ap,mock_illumina_run.dirn)
        self.assertTrue(os.path.isdir(
            '160621_M00879_0087_000000000-AGEW9'))

    def test_autoprocess_setup_absolute_paths(self):
        """setup: works when data dir path is absolute and normalised
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
        setup_(ap,data_dir_abs)
        self.assertEqual(ap.params.data_dir,data_dir_abs)
        self.assertEqual(ap.analysis_dir,analysis_dir)
        del(ap)
        shutil.rmtree(analysis_dir)
        # Do setup using relative path
        ap = AutoProcess()
        setup_(ap,data_dir_rel)
        self.assertEqual(ap.params.data_dir,data_dir_abs)
        self.assertEqual(ap.analysis_dir,analysis_dir)
        del(ap)
        shutil.rmtree(analysis_dir)
        # Do setup using absolute unnormalized path
        ap = AutoProcess()
        setup_(ap,data_dir_abs_unnormalised)
        self.assertEqual(ap.params.data_dir,data_dir_abs)
        self.assertEqual(ap.analysis_dir,analysis_dir)
        del(ap)

    def test_autoprocess_setup_missing_data_directory(self):
        """setup: raises exception if data directory is missing
        """
        # Set up autoprocessor
        ap = AutoProcess()
        self.assertRaises(Exception,
                          setup_,
                          ap,
                          os.path.join(
                              self.dirn,
                              '160621_M00879_0087_000000000-AGEW9'))
        self.assertFalse(os.path.exists(
            os.path.join(
                self.dirn,
                '160621_M00879_0087_000000000-AGEW9_analysis')))

    def test_autoprocess_setup_missing_sample_sheet(self):
        """setup: raises exception if sample sheet not found
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
                          setup_,
                          ap,
                          os.path.join(
                              self.dirn,
                              '160621_NB00879_0087_000000000-AGEW9'))
        self.assertFalse(os.path.exists(
            os.path.join(
                self.dirn,
                '160621_NB00879_0087_000000000-AGEW9_analysis')))

    def test_autoprocess_setup_samplesheet_from_url(self):
        """setup: works when samplesheet is a URL
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
        setup_(ap,mock_illumina_run.dirn,sample_sheet=sample_sheet)
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
