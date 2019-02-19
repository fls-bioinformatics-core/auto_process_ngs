#######################################################################
# Tests for merge_fastq_dirs.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
import gzip
from auto_process_ngs.mock import MockAnalysisDir
from auto_process_ngs.settings import Settings
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.commands.merge_fastq_dirs_cmd import merge_fastq_dirs

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
        # Create settings instance
        # This allows us to set the polling interval for the
        # unit tests
        settings_ini = os.path.join(self.dirn,"settings.ini")
        with open(settings_ini,'w') as s:
            s.write("""[general]
poll_interval = 0.5
""")
        self.settings = Settings(settings_ini)
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
        merge_fastq_dirs: bcl2fastq v2 output with --no-lane-splitting, new output dir
        """
        analysis_dir = self._setup_bcl2fastq2_no_lane_splitting()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir,
                              settings=self.settings)
        merge_fastq_dirs(self.ap,"bcl2fastq.AB",output_dir="bcl2fastq")
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
        expected = """#Project	Samples	User	Library	SC_Platform	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2_new_output_dir(self):
        """
        merge_fastq_dirs: bcl2fastq v2 output, new output dir
        """
        analysis_dir = self._setup_bcl2fastq2()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir,
                              settings=self.settings)
        merge_fastq_dirs(self.ap,
                         "bcl2fastq.lanes1-2",
                         output_dir="bcl2fastq")
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
        expected = """#Project	Samples	User	Library	SC_Platform	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_casava_new_output_dir(self):
        """
        merge_fastq_dirs: casava/bcl2fastq v1.8.* output, new output dir
        """
        analysis_dir = self._setup_casava()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir)
        merge_fastq_dirs(self.ap,
                         "bcl2fastq.lanes1-2",
                         output_dir="bcl2fastq")
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
        expected = """#Project	Samples	User	Library	SC_Platform	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2_no_lane_splitting(self):
        """
        merge_fastq_dirs: bcl2fastq v2 output with --no-lane-splitting
        """
        analysis_dir = self._setup_bcl2fastq2_no_lane_splitting()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir,
                              settings=self.settings)
        merge_fastq_dirs(self.ap,"bcl2fastq.AB")
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
        expected = """#Project	Samples	User	Library	SC_Platform	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2(self):
        """
        merge_fastq_dirs: bcl2fastq v2 output
        """
        analysis_dir = self._setup_bcl2fastq2()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir,
                              settings=self.settings)
        merge_fastq_dirs(self.ap,"bcl2fastq.lanes1-2")
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
        expected = """#Project	Samples	User	Library	SC_Platform	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_casava(self):
        """
        merge_fastq_dirs: casava/bcl2fastq v1.8.* output
        """
        analysis_dir = self._setup_casava()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir,
                              settings=self.settings)
        merge_fastq_dirs(self.ap,"bcl2fastq.lanes1-2")
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
        expected = """#Project	Samples	User	Library	SC_Platform	Organism	PI	Comments
AB	AB1,AB2	.	.	.	.	.	.
CDE	CDE3,CDE4	.	.	.	.	.	.
"""
        self.assertEqual(projects_info,expected)

    def test_bcl2fastq2_no_lane_splitting_dry_run(self):
        """
        merge_fastq_dirs: dry run on bcl2fastq v2 output with --no-lane-splitting
        """
        analysis_dir = self._setup_bcl2fastq2_no_lane_splitting()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir,
                              settings=self.settings)
        merge_fastq_dirs(self.ap,"bcl2fastq.AB",dry_run=True)
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
        merge_fastq_dirs: dry run on bcl2fastq v2 output
        """
        analysis_dir = self._setup_bcl2fastq2()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir,
                              settings=self.settings)
        merge_fastq_dirs(self.ap,"bcl2fastq.lanes1-2",dry_run=True)
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
        merge_fastq_dirs: dry run on casava/bcl2fastq v1.8.* output
        """
        analysis_dir = self._setup_casava()
        # Merge the unaligned dirs
        self.ap = AutoProcess(analysis_dir,
                              settings=self.settings)
        merge_fastq_dirs(self.ap,"bcl2fastq.lanes1-2",dry_run=True)
        # Check outputs
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.lanes1-2'))
        self._assert_dir_doesnt_exist(os.path.join(analysis_dir,'save.bcl2fastq.lanes3-4'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes1-2'))
        self._assert_dir_exists(os.path.join(analysis_dir,'bcl2fastq.lanes3-4'))
