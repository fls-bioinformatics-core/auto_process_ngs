#######################################################################
# Tests for autoprocessor.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from bcftbx.mock import MockIlluminaRun
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.mock import MockAnalysisDir

# Unit tests

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
        self.analysis_dir = None

    def tearDown(self):
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
        os.remove(os.path.join(m1,'projects.info'))
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
        os.remove(os.path.join(m1,'projects.info'))
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
        os.remove(os.path.join(m1,'projects.info'))
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

    def test_bcl2fastq2_no_lane_splitting(self):
        """
        AutoProcess.merge_fastq_dirs: bcl2fastq v2 output with --no-lane-splitting
        """
        analysis_dir = self._setup_bcl2fastq2_no_lane_splitting()
        # Merge the unaligned dirs
        ap = AutoProcess(analysis_dir)
        ap.merge_fastq_dirs("bcl2fastq.AB")
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

    def test_bcl2fastq2(self):
        """
        AutoProcess.merge_fastq_dirs: bcl2fastq v2 output
        """
        analysis_dir = self._setup_bcl2fastq2()
        # Merge the unaligned dirs
        ap = AutoProcess(analysis_dir)
        ap.merge_fastq_dirs("bcl2fastq.lanes1-2")
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

    def test_casava(self):
        """
        AutoProcess.merge_fastq_dirs: casava/bcl2fastq v1.8.* output
        """
        analysis_dir = self._setup_casava()
        # Merge the unaligned dirs
        ap = AutoProcess(analysis_dir)
        ap.merge_fastq_dirs("bcl2fastq.lanes1-2")
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
