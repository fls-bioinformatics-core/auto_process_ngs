#######################################################################
# Tests for setup_cmd.py module
#######################################################################

import unittest
import tempfile
import shutil
import gzip
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.settings import Settings
from bcftbx.mock import MockIlluminaRun
from bcftbx.mock import MockIlluminaData
from auto_process_ngs.commands.setup_cmd import setup as setup_

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

# Test data
fastq_r1_data = """@HISEQ:1:000000000-A2Y1L:1:1101:19264:2433 1:N:0:AGATCGC
AGATAGCCGA
+
?????BBB@B
@HISEQ:1:000000000-A2Y1L:1:1101:18667:2435 1:N:0:AGATCGC
ATATATTCAT
+
?????BBBDD
@HISEQ:1:000000000-A2Y1L:1:1101:17523:2436 1:N:0:AGATCGC
CATCACTACC
+
?<,<?BBBBB
"""
fastq_r2_data = """@HISEQ:1:000000000-A2Y1L:1:1101:19264:2433 2:N:0:AGATCGC
GCCGATATGC
+
??A??ABBDD
@HISEQ:1:000000000-A2Y1L:1:1101:18667:2435 2:N:0:AGATCGC
GATGACATCA
+
?????BBBDD
@HISEQ:1:000000000-A2Y1L:1:1101:17523:2436 2:N:0:AGATCGC
GAATATAGAA
+
??AAABBBDD
"""

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

    def settings(self,sequencer=None):
        # Create configuration for testing
        # If 'sequencer' is provided then should be a dictionary
        # with keys 'name', 'platform' and 'model' defined
        # e.g. dict(name='SN700215',platform='hiseq',model='HiSeq 2500')
        settings_ini = os.path.join(self.dirn,"auto_process.ini")
        with open(settings_ini,'w') as s:
            s.write("")
            if sequencer:
                s.write("""[sequencer:{name}]
platform = {platform}
model = {model}
""".format(**sequencer))
        return Settings(settings_ini)

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
        ap = AutoProcess(settings=self.settings())
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
        self.assertEqual(ap.metadata.analysis_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"miseq")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"M00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"151125")
        self.assertEqual(ap.metadata.instrument_run_number,"1")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:101bp, I1:8bp, I2:8bp, R2:101bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y101,I8,I8,y101")
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

    def test_autoprocess_setup_novaseq_run(self):
        """setup: works for mock NovaSeq run (set flow cell mode)
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '230711_A00879_0089_000000000-ABCDE1',
            'novaseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Make a samplesheet file
        sample_sheet = os.path.join(self.dirn,
                                    "SampleSheet_NovaSeq_Run_89.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Header]
IEMFileVersion,4

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTA,Project1,
2,Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,Project2,
""")
        # Set up autoprocessor
        ap = AutoProcess(
            settings=self.settings(
                sequencer=dict(name="A00879",
                               platform="novaseq",
                               model="NovaSeq 5000")))
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
                         "230711_A00879_0089_000000000-ABCDE1")
        self.assertEqual(ap.metadata.run_number,None)
        self.assertEqual(ap.metadata.analysis_number,None)
        self.assertEqual(ap.metadata.platform,"novaseq")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"A00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"230711")
        self.assertEqual(ap.metadata.instrument_run_number,"89")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        self.assertEqual(ap.metadata.flow_cell_mode,"SP")
        self.assertEqual(ap.metadata.sequencer_model,"NovaSeq 5000")
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:76bp, I1:10bp, I2:10bp, R2:76bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y76,I10,I10,y76")
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

    def test_autoprocess_setup_specify_facility_run_number(self):
        """setup: specify facility run number
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '151125_M00879_0001_000000000-ABCDE1',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up autoprocessor
        ap = AutoProcess(settings=self.settings())
        setup_(ap,mock_illumina_run.dirn,run_number='1')
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
        self.assertEqual(ap.metadata.run_number,'1')
        self.assertEqual(ap.metadata.analysis_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"miseq")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"M00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"151125")
        self.assertEqual(ap.metadata.instrument_run_number,"1")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:101bp, I1:8bp, I2:8bp, R2:101bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y101,I8,I8,y101")
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

    def test_autoprocess_setup_specify_analysis_number(self):
        """setup: specify analysis number
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '151125_M00879_0001_000000000-ABCDE1',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up autoprocessor
        ap = AutoProcess(settings=self.settings())
        setup_(ap,mock_illumina_run.dirn,analysis_number='2')
        analysis_dirn = "%s_analysis2" % mock_illumina_run.name
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
        self.assertEqual(ap.metadata.analysis_number,"2")
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"miseq")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"M00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"151125")
        self.assertEqual(ap.metadata.instrument_run_number,"1")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:101bp, I1:8bp, I2:8bp, R2:101bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y101,I8,I8,y101")
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
        ap = AutoProcess(settings=self.settings())
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
        self.assertEqual(ap.metadata.analysis_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"miseq")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,None)
        self.assertEqual(ap.metadata.instrument_datestamp,None)
        self.assertEqual(ap.metadata.instrument_run_number,None)
        self.assertEqual(ap.metadata.instrument_flow_cell_id,None)
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:101bp, I1:8bp, I2:8bp, R2:101bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y101,I8,I8,y101")

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
        ap = AutoProcess(settings=self.settings())
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
        ap = AutoProcess(settings=self.settings())
        setup_(ap,data_dir_abs)
        self.assertEqual(ap.params.data_dir,data_dir_abs)
        self.assertEqual(ap.analysis_dir,analysis_dir)
        del(ap)
        shutil.rmtree(analysis_dir)
        # Do setup using relative path
        ap = AutoProcess(settings=self.settings())
        setup_(ap,data_dir_rel)
        self.assertEqual(ap.params.data_dir,data_dir_abs)
        self.assertEqual(ap.analysis_dir,analysis_dir)
        del(ap)
        shutil.rmtree(analysis_dir)
        # Do setup using absolute unnormalized path
        ap = AutoProcess(settings=self.settings())
        setup_(ap,data_dir_abs_unnormalised)
        self.assertEqual(ap.params.data_dir,data_dir_abs)
        self.assertEqual(ap.analysis_dir,analysis_dir)
        del(ap)

    def test_autoprocess_setup_missing_data_directory(self):
        """setup: raises exception if data directory is missing
        """
        # Set up autoprocessor
        ap = AutoProcess(settings=self.settings())
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
        ap = AutoProcess(settings=self.settings())
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

    def test_autoprocess_setup_external_samplesheet(self):
        """setup: specify an external samplesheet file
        """
        # Create mock Illumina run directory with no samplesheet
        mock_illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Make a samplesheet file
        sample_sheet = os.path.join(self.dirn,"external_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Header]
IEMFileVersion,4
Date,01/22/2019
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTA,,
Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,,
""")
        # Set up autoprocessor
        ap = AutoProcess(settings=self.settings())
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
                         "171020_NB500968_00002_AHGXXXX")
        self.assertEqual(ap.metadata.run_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"nextseq")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"NB500968")
        self.assertEqual(ap.metadata.instrument_datestamp,"171020")
        self.assertEqual(ap.metadata.instrument_run_number,"2")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "AHGXXXX")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:76bp, I1:6bp, R2:76bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y76,I6,y76")
        # Delete to force write of data to disk
        del(ap)
        # Check directory exists
        self.assertTrue(os.path.isdir(analysis_dirn))
        # Check files exist
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
        print(sample_sheet)
        # Set up autoprocessor
        ap = AutoProcess(settings=self.settings())
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
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"M00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"151125")
        self.assertEqual(ap.metadata.instrument_run_number,"1")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:101bp, I1:8bp, I2:8bp, R2:101bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y101,I8,I8,y101")
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

    def test_autoprocess_setup_samplesheet_with_lanes(self):
        """setup: specify an external samplesheet file with lanes
        """
        # Create mock Illumina run directory with no samplesheet
        mock_illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Make a samplesheet file
        sample_sheet = os.path.join(self.dirn,"external_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Header]
IEMFileVersion,4
Date,01/22/2019
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTA,,
2,Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,,
""")
        # Set up autoprocessor
        ap = AutoProcess(settings=self.settings())
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
                         "171020_NB500968_00002_AHGXXXX")
        self.assertEqual(ap.metadata.run_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"nextseq")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"NB500968")
        self.assertEqual(ap.metadata.instrument_datestamp,"171020")
        self.assertEqual(ap.metadata.instrument_run_number,"2")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "AHGXXXX")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:76bp, I1:6bp, R2:76bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y76,I6,y76")
        # Delete to force write of data to disk
        del(ap)
        # Check directory exists
        self.assertTrue(os.path.isdir(analysis_dirn))
        # Check files exist
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

    def test_autoprocess_setup_samplesheet_with_lanes_and_trailing_lines(self):
        """setup: specify an external samplesheet file with lanes (trailing lines)
        """
        # Create mock Illumina run directory with no samplesheet
        mock_illumina_run = MockIlluminaRun(
            "171020_NB500968_00002_AHGXXXX",
            "nextseq",
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Make a samplesheet file
        sample_sheet = os.path.join(self.dirn,"external_SampleSheet.csv")
        with open(sample_sheet,'w') as fp:
            fp.write("""[Header]
IEMFileVersion,4
Date,01/22/2019
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,TruSeq HT
Description,
Chemistry,Amplicon

[Reads]
101
101

[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,Sample1,Sample1,,,D701,CGTGTAGG,D501,GACCTGTA,,
2,Sample2,Sample2,,,D702,CGTGTAGG,D501,ATGTAACT,,
,,,,,,,,,,,
,,,,,,,,,,,
,,,,,,,,,,,
""")
        # Set up autoprocessor
        ap = AutoProcess(settings=self.settings())
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
                         "171020_NB500968_00002_AHGXXXX")
        self.assertEqual(ap.metadata.run_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"nextseq")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"NB500968")
        self.assertEqual(ap.metadata.instrument_datestamp,"171020")
        self.assertEqual(ap.metadata.instrument_run_number,"2")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "AHGXXXX")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:76bp, I1:6bp, R2:76bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y76,I6,y76")
        # Delete to force write of data to disk
        del(ap)
        # Check directory exists
        self.assertTrue(os.path.isdir(analysis_dirn))
        # Check files exist
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

    def test_autoprocess_setup_import_extra_files(self):
        """setup: check extra files can be imported
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '151125_M00879_0001_000000000-ABCDE1',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Create additional files
        extra_file1 = os.path.join(self.dirn,'extra_file1.txt')
        with open(extra_file1,'wt') as fp:
            fp.write(u"This is extra file 1\n")
        extra_file2 = os.path.join(self.dirn,'extra_file2.txt')
        with open(extra_file2,'wt') as fp:
            fp.write(u"This is extra file 2\n")
        # Make extra_file2 into a URL
        extra_file2 = "file://%s" % extra_file2
        # Set up autoprocessor
        ap = AutoProcess(settings=self.settings())
        setup_(ap,mock_illumina_run.dirn,extra_files=(extra_file1,
                                                      extra_file2))
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
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"M00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"151125")
        self.assertEqual(ap.metadata.instrument_run_number,"1")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:101bp, I1:8bp, I2:8bp, R2:101bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y101,I8,I8,y101")
        # Delete to force write of data to disk
        del(ap)
        # Check directory exists
        self.assertTrue(os.path.isdir(analysis_dirn))
        # Check files exist
        for filen in ('SampleSheet.orig.csv',
                      'custom_SampleSheet.csv',
                      'auto_process.info',
                      'metadata.info',
                      'extra_file1.txt',
                      'extra_file2.txt',):
            self.assertTrue(os.path.exists(os.path.join(analysis_dirn,
                                                        filen)),
                            "Missing file: %s" % filen)
        # Check subdirs have been created
        for subdirn in ('ScriptCode',
                        'logs',):
            self.assertTrue(os.path.isdir(os.path.join(analysis_dirn,
                                                       subdirn)),
                            "Missing subdir: %s" % subdirn)

    def test_autoprocess_setup_from_casava_outputs(self):
        """setup: get Fastqs from existing CASAVA-style outputs
        """
        # Create mock Illumina run directory
        casava_outputs = MockIlluminaData(
            '151125_M00879_0001_000000000-ABCDE1',
            'casava',
            paired_end=True,
            top_dir=self.dirn)
        casava_outputs.add_fastq('AB','AB1','AB1_GCCAAT_L001_R1_001.fastq.gz')
        casava_outputs.add_fastq('AB','AB1','AB1_GCCAAT_L001_R2_001.fastq.gz')
        casava_outputs.add_fastq('AB','AB2','AB2_AGTCAA_L001_R1_001.fastq.gz')
        casava_outputs.add_fastq('AB','AB2','AB2_AGTCAA_L001_R2_001.fastq.gz')
        casava_outputs.add_fastq('CDE','CDE3','CDE3_GCCAAT_L001_R1_001.fastq.gz')
        casava_outputs.add_fastq('CDE','CDE3','CDE3_GCCAAT_L001_R2_001.fastq.gz')
        casava_outputs.add_fastq('CDE','CDE4','CDE4_AGTCAA_L001_R1_001.fastq.gz')
        casava_outputs.add_fastq('CDE','CDE4','CDE4_AGTCAA_L001_R2_001.fastq.gz')
        casava_outputs.create()
        # Populate Fastqs
        for root,dirs,files in os.walk(casava_outputs.dirn):
            for f in files:
                if f.endswith(".fastq.gz"):
                    with gzip.open(os.path.join(root,f),'wt') as fq:
                        try:
                            f.index("_R1_")
                            fq.write(fastq_r1_data)
                        except ValueError:
                            fq.write(fastq_r2_data)
        # Set up autoprocessor
        run_dir = "151125_M00879_0001_000000000-ABCDE1"
        analysis_dir = "%s_analysis" % run_dir
        ap = AutoProcess(settings=self.settings())
        setup_(ap,run_dir,unaligned_dir=casava_outputs.unaligned_dir)
        # Check parameters
        self.assertEqual(ap.analysis_dir,
                         os.path.join(self.dirn,analysis_dir))
        self.assertEqual(ap.params.unaligned_dir,casava_outputs.unaligned_dir)
        self.assertEqual(ap.params.data_dir,os.path.join(self.dirn,run_dir))
        self.assertEqual(ap.params.sample_sheet,None)
        self.assertEqual(ap.params.bases_mask,"auto")
        # Check metadata
        self.assertEqual(ap.metadata.run_name,run_dir)
        self.assertEqual(ap.metadata.run_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"miseq")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"M00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"151125")
        self.assertEqual(ap.metadata.instrument_run_number,"1")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,None)
        self.assertEqual(ap.metadata.default_bases_mask,None)
        # Delete to force write of data to disk
        del(ap)
        # Check directory exists
        self.assertTrue(os.path.isdir(analysis_dir))
        # Check files exists
        for filen in ('auto_process.info',
                      'metadata.info',
                      'projects.info',):
            self.assertTrue(os.path.exists(os.path.join(analysis_dir,
                                                        filen)),
                            "Missing file: %s" % filen)
        # Check contents of projects.info
        projects_info = os.path.join(analysis_dir,"projects.info")
        with open(projects_info,'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
""")
        # Check subdirs have been created
        for subdirn in ('ScriptCode',
                        'logs',):
            self.assertTrue(os.path.isdir(os.path.join(analysis_dir,
                                                       subdirn)),
                            "Missing subdir: %s" % subdirn)

    def test_autoprocess_setup_from_bcl2fastq2_outputs(self):
        """setup: get Fastqs from existing bcl2fastq2-style outputs
        """
        # Create mock Illumina run directory
        bcl2fastq2_outputs = MockIlluminaData(
            '151125_M00879_0001_000000000-ABCDE1',
            'bcl2fastq2',
            paired_end=True,
            top_dir=self.dirn)
        bcl2fastq2_outputs.add_fastq('AB','AB1','AB1_S1_L001_R1_001.fastq.gz')
        bcl2fastq2_outputs.add_fastq('AB','AB1','AB1_S1_L001_R2_001.fastq.gz')
        bcl2fastq2_outputs.add_fastq('AB','AB2','AB2_S2_L001_R1_001.fastq.gz')
        bcl2fastq2_outputs.add_fastq('AB','AB2','AB2_S2_L001_R2_001.fastq.gz')
        bcl2fastq2_outputs.add_fastq('CDE','CDE3','CDE3_S3_L001_R1_001.fastq.gz')
        bcl2fastq2_outputs.add_fastq('CDE','CDE3','CDE3_S3_L001_R2_001.fastq.gz')
        bcl2fastq2_outputs.add_fastq('CDE','CDE4','CDE4_S4_L001_R1_001.fastq.gz')
        bcl2fastq2_outputs.add_fastq('CDE','CDE4','CDE4_S4_L001_R2_001.fastq.gz')
        bcl2fastq2_outputs.create()
        # Populate Fastqs
        for root,dirs,files in os.walk(bcl2fastq2_outputs.dirn):
            for f in files:
                if f.endswith(".fastq.gz"):
                    with gzip.open(os.path.join(root,f),'wt') as fq:
                        try:
                            f.index("_R1_")
                            fq.write(fastq_r1_data)
                        except ValueError:
                            fq.write(fastq_r2_data)
        # Set up autoprocessor
        run_dir = "151125_M00879_0001_000000000-ABCDE1"
        analysis_dir = "%s_analysis" % run_dir
        ap = AutoProcess(settings=self.settings())
        setup_(ap,run_dir,unaligned_dir=bcl2fastq2_outputs.unaligned_dir)
        # Check parameters
        self.assertEqual(ap.analysis_dir,
                         os.path.join(self.dirn,analysis_dir))
        self.assertEqual(ap.params.unaligned_dir,bcl2fastq2_outputs.unaligned_dir)
        self.assertEqual(ap.params.data_dir,os.path.join(self.dirn,run_dir))
        self.assertEqual(ap.params.sample_sheet,None)
        self.assertEqual(ap.params.bases_mask,"auto")
        # Check metadata
        self.assertEqual(ap.metadata.run_name,run_dir)
        self.assertEqual(ap.metadata.run_number,None)
        self.assertEqual(ap.metadata.source,None)
        self.assertEqual(ap.metadata.platform,"miseq")
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"M00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"151125")
        self.assertEqual(ap.metadata.instrument_run_number,"1")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,None)
        self.assertEqual(ap.metadata.run_configuration,None)
        self.assertEqual(ap.metadata.default_bases_mask,None)
        # Delete to force write of data to disk
        del(ap)
        # Check directory exists
        self.assertTrue(os.path.isdir(analysis_dir))
        # Check files exists
        for filen in ('auto_process.info',
                      'metadata.info',
                      'projects.info',):
            self.assertTrue(os.path.exists(os.path.join(analysis_dir,
                                                        filen)),
                            "Missing file: %s" % filen)
        # Check contents of projects.info
        projects_info = os.path.join(analysis_dir,"projects.info")
        with open(projects_info,'rt') as fp:
            self.assertEqual(fp.read(),
                             """#Project\tSamples\tUser\tLibrary\tSC_Platform\tOrganism\tPI\tComments
AB\tAB1,AB2\t.\t.\t.\t.\t.\t.
CDE\tCDE3,CDE4\t.\t.\t.\t.\t.\t.
""")
        # Check subdirs have been created
        for subdirn in ('ScriptCode',
                        'logs',):
            self.assertTrue(os.path.isdir(os.path.join(analysis_dir,
                                                       subdirn)),
                            "Missing subdir:  %s" % subdirn)

    def test_autoprocess_setup_fails_for_empty_fastq_dir(self):
        """setup: fails for empty unaligned Fastq directory
        """
        # Create empty Fastq directory
        run_dir = os.path.join(self.dirn,"151125_M00879_0001_000000000-ABCDE1")
        unaligned_dir = os.path.join(run_dir,"Unaligned")
        os.mkdir(run_dir)
        os.mkdir(unaligned_dir)
        # Set up autoprocessor
        analysis_dir = "%s_analysis" % run_dir
        ap = AutoProcess(settings=self.settings())
        self.assertRaises(Exception,
                          setup_,
                          ap,
                          run_dir,
                          unaligned_dir=unaligned_dir)
        self.assertFalse(os.path.exists(analysis_dir))

    def test_autoprocess_setup_fails_for_invalid_fastq_dir(self):
        """setup: fails for invalid unaligned Fastq directory
        """
        # Create and populate invalid Fastq directory
        run_dir = os.path.join(self.dirn,"151125_M00879_0001_000000000-ABCDE1")
        unaligned_dir = os.path.join(run_dir,"Unaligned")
        os.mkdir(run_dir)
        os.mkdir(unaligned_dir)
        for f in ('AB1_GCCAAT_L001_R1_001.fastq.gz',
                  'AB1_GCCAAT_L001_R2_001.fastq.gz',
                  'AB2_AGTCAA_L001_R1_001.fastq.gz',
                  'AB2_AGTCAA_L001_R2_001.fastq.gz',):
            if f.endswith(".fastq.gz"):
                with gzip.open(os.path.join(unaligned_dir,f),'wt') as fq:
                    try:
                        f.index("_R1_")
                        fq.write(fastq_r1_data)
                    except ValueError:
                        fq.write(fastq_r2_data)
        # Set up autoprocessor
        analysis_dir = "%s_analysis" % run_dir
        ap = AutoProcess(settings=self.settings())
        self.assertRaises(Exception,
                          setup_,
                          ap,
                          run_dir,
                          unaligned_dir=unaligned_dir)
        self.assertFalse(os.path.exists(analysis_dir))

    def test_autoprocess_setup_stores_sequencer_model(self):
        """setup: stores sequencer model
        """
        # Create mock Illumina run directory
        mock_illumina_run = MockIlluminaRun(
            '151125_M00879_0001_000000000-ABCDE1',
            'miseq',
            top_dir=self.dirn)
        mock_illumina_run.create()
        # Set up autoprocessor
        ap = AutoProcess(settings=self.settings(
            sequencer=dict(name="M00879",
                           platform="miseq",
                           model="MiSeq")))
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
        self.assertEqual(ap.metadata.bcl2fastq_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.cellranger_software,None)
        self.assertEqual(ap.metadata.instrument_name,"M00879")
        self.assertEqual(ap.metadata.instrument_datestamp,"151125")
        self.assertEqual(ap.metadata.instrument_run_number,"1")
        self.assertEqual(ap.metadata.instrument_flow_cell_id,
                         "000000000-ABCDE1")
        self.assertEqual(ap.metadata.flow_cell_mode,None)
        self.assertEqual(ap.metadata.sequencer_model,"MiSeq")
        self.assertEqual(ap.metadata.run_configuration,
                         "R1:101bp, I1:8bp, I2:8bp, R2:101bp")
        self.assertEqual(ap.metadata.default_bases_mask,
                         "y101,I8,I8,y101")
