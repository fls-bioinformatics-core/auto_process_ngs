#######################################################################
# Unit tests for qc/picard.py
#######################################################################

import unittest
import os
import tempfile
import shutil

from auto_process_ngs.mockqcdata import PICARD_COLLECT_INSERT_SIZE_METRICS
from auto_process_ngs.qc.apps.picard import CollectInsertSizeMetrics
from auto_process_ngs.qc.apps.picard import picard_collect_insert_size_metrics_output

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

class TestPicardCollectInsertSizeMetricsOutputFunction(unittest.TestCase):
    
    def test_picard_collect_insert_size_metrics_output_fastq(self):
        """
        picard_collect_insert_size_metrics_output: no prefix (Fastq file)
        """
        self.assertEqual(
            picard_collect_insert_size_metrics_output(
                '/data/PB/PB1_ATTAGG_L001_R1_001.fastq'),
            ('PB1_ATTAGG_L001_R1_001.insert_size_metrics.txt',
             'PB1_ATTAGG_L001_R1_001.insert_size_histogram.pdf'))

    def test_picard_collect_insert_size_metrics_output_bam_file(self):
        """
        picard_collect_insert_size_metrics_output: no prefix (BAM file)
        """
        self.assertEqual(
            picard_collect_insert_size_metrics_output(
                '/data/PB/PB1_ATTAGG_L001_R1_001.bam'),
            ('PB1_ATTAGG_L001_R1_001.insert_size_metrics.txt',
             'PB1_ATTAGG_L001_R1_001.insert_size_histogram.pdf'))

    def test_picard_collect_insert_size_metrics_output_with_prefix(self):
        """
        picard_collect_insert_size_metrics_output: with prefix
        """
        self.assertEqual(
            picard_collect_insert_size_metrics_output(
                '/data/PB/PB1_ATTAGG_L001_R1_001.fastq',
                prefix="picard/human"),
            ('picard/human/PB1_ATTAGG_L001_R1_001.insert_size_metrics.txt',
             'picard/human/PB1_ATTAGG_L001_R1_001.insert_size_histogram.pdf'))
