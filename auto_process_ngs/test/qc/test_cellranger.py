#######################################################################
# Unit tests for qc/cellranger.py
#######################################################################

import unittest
import os
import shutil
import tempfile
from auto_process_ngs.mock import MockAnalysisProject
from auto_process_ngs.mock import UpdateAnalysisProject
from auto_process_ngs.analysis import AnalysisProject

from auto_process_ngs.qc.cellranger import CellrangerCount

class TestCellrangerCount(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestCellrangerCount')
        # Make mock analysis project
        p = MockAnalysisProject("PJB",("PJB1_S1_R1_001.fastq.gz",
                                       "PJB1_S1_R2_001.fastq.gz",
                                       "PJB2_S2_R1_001.fastq.gz",
                                       "PJB2_S2_R2_001.fastq.gz",),
                                metadata={ 'Organism': 'Human',
                                           'Single cell platform':
                                           "10xGenomics Chromium 3'v3" })
        p.create(top_dir=self.dirn)
        self.project = AnalysisProject("PJB",os.path.join(self.dirn,"PJB"))
    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)
    def test_cellrangercount(self):
        """
        CellrangerCount: check outputs from cellranger count
        """
        # Add cellranger count outputs
        UpdateAnalysisProject(self.project).add_cellranger_count_outputs()
        # Do tests
        count_dir = os.path.join(self.project.dirn,"cellranger_count","PJB1")
        cellranger_count = CellrangerCount(count_dir)
        self.assertEqual(cellranger_count.dir,count_dir)
        self.assertEqual(cellranger_count.sample_name,"PJB1")
        self.assertEqual(cellranger_count.metrics_csv,
                         os.path.join(count_dir,"outs","metrics_summary.csv"))
        self.assertEqual(cellranger_count.web_summary,
                         os.path.join(count_dir,"outs","web_summary.html"))
