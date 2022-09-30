#######################################################################
# Unit tests for qc/picard.py
#######################################################################

import unittest
import os
import tempfile
import shutil

from auto_process_ngs.mockqcdata import PICARD_COLLECT_INSERT_SIZE_METRICS
from auto_process_ngs.qc.picard import CollectInsertSizeMetrics

class TestCollectInsertSizeMetrics(unittest.TestCase):
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestPicard')
        # Make a fake CollectInsertSizeMetrics file
        self.insert_size_metrics_out = os.path.join(
            self.dirn,
            "test.insert_size_metrics.txt")
        with open(self.insert_size_metrics_out,'wt') as fp:
            fp.write(PICARD_COLLECT_INSERT_SIZE_METRICS)

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_collect_insert_size_metrics(self):
        """
        CollectInsertSizeMetrics: read in metrics file
        """
        insert_size = CollectInsertSizeMetrics(self.insert_size_metrics_out)
        self.assertEqual(insert_size.metrics_file,
                         self.insert_size_metrics_out)
        self.assertEqual(insert_size.metrics['MEAN_INSERT_SIZE'],
                         153.754829)
        self.assertEqual(insert_size.histogram[28],1)
        self.assertEqual(insert_size.histogram[29],1)
        self.assertEqual(insert_size.histogram[30],1)
        self.assertEqual(insert_size.histogram[31],1)
        self.assertEqual(insert_size.histogram[32],1)
        self.assertEqual(insert_size.histogram[33],2)
        self.assertEqual(insert_size.histogram[34],3)
