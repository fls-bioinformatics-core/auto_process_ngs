#######################################################################
# Tests for tenx_genomics/multiome.py module
#######################################################################

import unittest
import os
import shutil
import tempfile
from auto_process_ngs.mock import MockAnalysisDir
from auto_process_ngs.mock10xdata import MULTIOME_LIBRARIES
from auto_process_ngs.mock10xdata import MULTIOME_LIBRARIES_NO_RUN
from auto_process_ngs.tenx.multiome import *

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestMultiomeLibraries(unittest.TestCase):
    """
    Tests for the 'MultiomeLibraries' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestMultiomeLibraries')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_multiome_libraries(self):
        """MultiomeLibraries: read in data from a file
        """
        multiome_libraries_info = os.path.join(self.wd,
                                               "10x_multiome_libraries.info")
        with open(multiome_libraries_info,'wt') as fp:
            fp.write(MULTIOME_LIBRARIES)
        m = MultiomeLibraries(multiome_libraries_info)
        # Check list of samples
        self.assertEqual(m.local_samples,['PB1_ATAC','PB2_ATAC'])
        # Check lists of linked samples
        self.assertEqual(m.linked_samples('PB1_ATAC'),
                         ['NEXTSEQ_210111/12:PB_GEX/PB1_GEX',])
        self.assertEqual(m.linked_samples('PB2_ATAC'),
                         ['NEXTSEQ_210111/12:PB_GEX/PB2_GEX',])

    def test_multiome_libraries_linked_projects(self):
        """MultiomeLibraries: return linked projects
        """
        multiome_libraries_info = os.path.join(self.wd,
                                               "10x_multiome_libraries.info")
        with open(multiome_libraries_info,'wt') as fp:
            fp.write(MULTIOME_LIBRARIES)
        m = MultiomeLibraries(multiome_libraries_info)
        # Create run with linked project
        run = MockAnalysisDir("210111_NB01234_00012_ABXXXXX",
                              "nextseq",
                              project_metadata={
                                  'PB_GEX': {
                                      'Library type': 'GEX',
                                  },
                              },
                              lanes=(1,2,3,4),
                              top_dir=self.wd)
        run.add_fastq_batch('PB_GEX','PB1_GEX','PB1_GEX_S1',lanes=(1,2,3,4))
        run.add_fastq_batch('PB_GEX','PB2_GEX','PB2_GEX_S2',lanes=(1,2,3,4))
        run.create()
        # Check list of linked projects
        linked_projects = m.linked_projects()
        self.assertEqual(len(linked_projects),1)
        self.assertEqual(linked_projects[0].dirn,
                         os.path.join(self.wd,
                                      "210111_NB01234_00012_ABXXXXX_analysis",
                                      "PB_GEX"))

    def test_multiome_libraries_linked_projects_same_run(self):
        """MultiomeLibraries: return linked projects from same run
        """
        # Create run with linked projects
        run = MockAnalysisDir("210111_NB01234_00012_ABXXXXX",
                              "nextseq",
                              project_metadata={
                                  'PB_ATAC': {
                                      'Library type': 'ATAC',
                                  },
                                  'PB_GEX': {
                                      'Library type': 'GEX',
                                  },
                              },
                              lanes=(1,2,3,4),
                              top_dir=self.wd)
        run.add_fastq_batch('PB_ATAC','PB1_ATAC','PB1_ATAC_S1',lanes=(1,2,3,4))
        run.add_fastq_batch('PB_ATAC','PB2_ATAC','PB2_ATAC_S2',lanes=(1,2,3,4))
        run.add_fastq_batch('PB_GEX','PB1_GEX','PB1_GEX_S3',lanes=(1,2,3,4))
        run.add_fastq_batch('PB_GEX','PB2_GEX','PB2_GEX_S4',lanes=(1,2,3,4))
        run.create()
        # Put multiome libraries file into ATAC project
        multiome_libraries_info = os.path.join(run.dirn,
                                               "PB_ATAC",
                                               "10x_multiome_libraries.info")
        with open(multiome_libraries_info,'wt') as fp:
            fp.write(MULTIOME_LIBRARIES_NO_RUN)
        m = MultiomeLibraries(multiome_libraries_info)
        # Check list of linked projects
        linked_projects = m.linked_projects()
        self.assertEqual(len(linked_projects),1)
        self.assertEqual(linked_projects[0].dirn,
                         os.path.join(self.wd,
                                      "210111_NB01234_00012_ABXXXXX_analysis",
                                      "PB_GEX"))

    def test_multiome_libraries_write_libraries_csv(self):
        """MultiomeLibraries: write 'libraries.csv' files
        """
        multiome_libraries_info = os.path.join(self.wd,
                                               "10x_multiome_libraries.info")
        with open(multiome_libraries_info,'wt') as fp:
            fp.write(MULTIOME_LIBRARIES)
        m = MultiomeLibraries(multiome_libraries_info)
        # Create run with linked project
        run = MockAnalysisDir("210111_NB01234_00012_ABXXXXX",
                              "nextseq",
                              project_metadata={
                                  'PB_GEX': {
                                      'Library type': 'GEX',
                                  },
                              },
                              lanes=(1,2,3,4),
                              top_dir=self.wd)
        run.add_fastq_batch('PB_GEX','PB1_GEX','PB1_GEX_S1',lanes=(1,2,3,4))
        run.add_fastq_batch('PB_GEX','PB2_GEX','PB2_GEX_S1',lanes=(1,2,3,4))
        run.create()
        # Write libraries.csv file
        libraries_csv = os.path.join(self.wd,"libraries.csv")
        self.assertFalse(os.path.exists(libraries_csv))
        m.write_libraries_csv('PB1_ATAC',
                              '/data/PB_ATAC/fastqs',
                              'ATAC',
                              filen=libraries_csv)
        self.assertTrue(os.path.exists(libraries_csv))
        with open(libraries_csv,'rt') as fp:
            contents = fp.read()
            self.assertEqual(contents,
                             """fastqs,sample,library_type
/data/PB_ATAC/fastqs,PB1_ATAC,Chromatin Accessibility
%s/210111_NB01234_00012_ABXXXXX_analysis/PB_GEX/fastqs,PB1_GEX,Gene Expression
""" % self.wd)
