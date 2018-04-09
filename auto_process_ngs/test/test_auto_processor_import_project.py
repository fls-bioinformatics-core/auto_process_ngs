#######################################################################
# Tests for autoprocessor.py module: import_project
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.mock import MockAnalysisDirFactory
from auto_process_ngs.auto_processor import AutoProcess

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
SC_Platform\t.
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
