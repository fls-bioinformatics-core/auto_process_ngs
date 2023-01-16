#######################################################################
# Tests for tenx_genomics/metrics.py module
#######################################################################

import unittest
import tempfile
import os
import shutil
from auto_process_ngs.mock10xdata import METRICS_SUMMARY
from auto_process_ngs.mock10xdata import ATAC_SUMMARY
from auto_process_ngs.mock10xdata import ATAC_SUMMARY_2_0_0
from auto_process_ngs.mock10xdata import CELLPLEX_METRICS_SUMMARY
from auto_process_ngs.mock10xdata import MULTIOME_SUMMARY
from auto_process_ngs.mock10xdata import MULTIOME_SUMMARY_2_0_0
from auto_process_ngs.tenx.metrics import *

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestGexSummary(unittest.TestCase):
    """
    Tests for the 'GexSummary' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestGexSummary')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)
            
    def test_gex_summary(self):
        """GexSummary: check metrics are extracted from CSV file
        """
        summary_csv = os.path.join(self.wd,"metrics_summary.csv")
        with open(summary_csv,'w') as fp:
            fp.write(METRICS_SUMMARY)
        m = GexSummary(summary_csv)
        self.assertEqual(m.estimated_number_of_cells,2272)
        self.assertEqual(m.mean_reads_per_cell,107875)
        self.assertEqual(m.median_genes_per_cell,1282)

class TestAtacSummary(unittest.TestCase):
    """
    Tests for the 'AtacSummary' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestAtacSummary')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_atac_summary_pre_2_0_0(self):
        """AtacSummary: get metrics from summary.csv (Cellranger ATAC pre-2.0.0)
        """
        summary_csv = os.path.join(self.wd,"summary.csv")
        with open(summary_csv,'w') as fp:
            fp.write(ATAC_SUMMARY)
        s = AtacSummary(summary_csv)
        self.assertEqual(s.version,"1.0.1")
        self.assertEqual(s.cells_detected,6748)
        self.assertEqual(s.annotated_cells,5682)
        self.assertEqual(s.median_fragments_per_cell,16119.5)
        self.assertEqual(s.frac_fragments_overlapping_targets,
                         0.575082094792)
        self.assertEqual(s.frac_fragments_overlapping_peaks,
                         0.556428090013)
        self.assertEqual(s.tss_enrichment_score,6.91438390781)
        self.assertRaises(AttributeError,
                          getattr,
                          s,'estimated_number_of_cells')

    def test_atac_summary_2_0_0(self):
        """AtacSummary: get metrics from summary.csv (Cellranger ATAC 2.0.0)
        """
        summary_csv = os.path.join(self.wd,"summary.csv")
        with open(summary_csv,'w') as fp:
            fp.write(ATAC_SUMMARY_2_0_0)
        s = AtacSummary(summary_csv)
        self.assertEqual(s.version,"2.0.0")
        self.assertEqual(s.estimated_number_of_cells,3582)
        self.assertEqual(s.median_fragments_per_cell,51354.5)
        self.assertEqual(s.frac_fragments_overlapping_peaks,0.4856)
        self.assertEqual(s.tss_enrichment_score,6.8333)
        self.assertRaises(AttributeError,
                          getattr,
                          s,'cells_detected')
        self.assertRaises(AttributeError,
                          getattr,
                          s,'annotated_cells')
        self.assertRaises(AttributeError,
                          getattr,
                          s,'frac_fragments_overlapping_targets')

class TestMultiomeSummary(unittest.TestCase):
    """
    Tests for the 'MultiomeSummary' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestMultiomeSummary')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_atac_summary_multiome_arc_1_0_0(self):
        """MultiomeSummary: check metrics are extracted from CSV file (Cellranger ARC 1.0.0)
        """
        summary_csv = os.path.join(self.wd,"summary.csv")
        with open(summary_csv,'wt') as fp:
            fp.write(MULTIOME_SUMMARY)
        s = MultiomeSummary(summary_csv)
        self.assertEqual(s.estimated_number_of_cells,744)
        self.assertEqual(
            s.atac_median_high_quality_fragments_per_cell,8079)
        self.assertEqual(s.gex_median_genes_per_cell,1490)

    def test_atac_summary_multiome_arc_2_0_0(self):
        """MultiomeSummary: check metrics are extracted from CSV file (Cellranger ARC 2.0.0)
        """
        summary_csv = os.path.join(self.wd,"summary.csv")
        with open(summary_csv,'wt') as fp:
            fp.write(MULTIOME_SUMMARY_2_0_0)
        s = MultiomeSummary(summary_csv)
        self.assertEqual(s.estimated_number_of_cells,785)
        self.assertEqual(
            s.atac_median_high_quality_fragments_per_cell,9.0)
        self.assertEqual(s.gex_median_genes_per_cell,15.0)

class TestMultiplexSummary(unittest.TestCase):
    """
    Tests for the 'MultiplexSummary' class
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestMultiplexSummary')

    def tearDown(self):
        # Remove the temporary test directory
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_multiplex_summary(self):
        """MultiplexSummary: check metrics are extracted from CSV file
        """
        summary_csv = os.path.join(self.wd,"metrics_summary.csv")
        with open(summary_csv,'w') as fp:
            fp.write(CELLPLEX_METRICS_SUMMARY)
        m = MultiplexSummary(summary_csv)
        self.assertEqual(m.cells,5175)
        self.assertEqual(m.median_reads_per_cell,20052)
        self.assertEqual(m.median_genes_per_cell,3086)
        self.assertEqual(m.total_genes_detected,21260)
        self.assertEqual(m.median_umi_counts_per_cell,10515)
