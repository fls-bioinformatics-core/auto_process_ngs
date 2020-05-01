#######################################################################
# Unit tests for bcl2fastq/reporting.py
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.bcl2fastq.reporting import ProcessingQCReport
from auto_process_ngs.bcl2fastq.reporting import detect_processing_qc_warnings

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = True

class TestProcessingQCReport(unittest.TestCase):
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.TestProcessingQCReport')
    def tearDown(self):
        # Remove temporary working dir
        if not REMOVE_TEST_OUTPUTS:
            return
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_processing_qc_report(self):
        """
        ProcessingQCReport: generate HTML report
        """
        # Make mock stats files
        stats_file = os.path.join(self.wd,"statistics_full.info")
        with open(stats_file,'wt') as fp:
            fp.write(
"""        #Project	Sample	Fastq	Size	Nreads	Paired_end	Read_number	L1
LFB	LFB1	LFB1_S1_L001_R1_001.fastq.gz	16.5M	379988	Y	1	379988
LFB	LFB1	LFB1_S1_L001_R2_001.fastq.gz	16.9M	379988	Y	2	379988
LFB	LFB2	LFB2_S2_L001_R1_001.fastq.gz	14.0M	327235	Y	1	327235
LFB	LFB2	LFB2_S2_L001_R2_001.fastq.gz	13.9M	327235	Y	2	327235
Undetermined_indices	lane1	Undetermined_S0_L001_R1_001.fastq.gz	22.8M	436774	Y	1	436774
Undetermined_indices	lane1	Undetermined_S0_L001_R2_001.fastq.gz	22.9M	436774	Y	2	436774
""")
        per_lane_stats_file = os.path.join(
            self.wd,"per_lane_statistics.info")
        with open(per_lane_stats_file,'wt') as fp:
            fp.write(
"""#Lane	Total reads	Assigned reads	Unassigned reads	%assigned	%unassigned
Lane 1	13132812	12696038	436774	96.67	3.33
""")
        per_lane_sample_stats_file = os.path.join(
            self.wd,"per_lane_sample_stats.info")
        with open(per_lane_sample_stats_file,'wt') as fp:
            fp.write(
"""
Lane 1
Total reads = 13132812
- LFB/LFB1	379988	2.9%
- LFB/LFB2	327235	2.5%
- Undetermined_indices/lane1	436774	3.3%
""")
        # Generate HTML report
        reporter = ProcessingQCReport(self.wd,
                                      stats_file,
                                      per_lane_stats_file,
                                      per_lane_sample_stats_file)
        reporter.write(os.path.join(self.wd,"processing_qc_report.html"))
        self.assertTrue(os.path.exists(
            os.path.join(self.wd,"processing_qc_report.html")))

class TestDetectProcessingQCWarnings(unittest.TestCase):
    def setUp(self):
        # Temporary working dir
        self.wd = tempfile.mkdtemp(suffix='.TestDetectProcessingQCWarnings')
    def tearDown(self):
        # Remove temporary working dir
        if not REMOVE_TEST_OUTPUTS:
            return
        if self.wd is not None and os.path.isdir(self.wd):
            shutil.rmtree(self.wd)
    def test_detect_warnings(self):
        """
        detect_processing_qc_warnings: report contains warnings
        """
        mock_report_file = os.path.join(self.wd,"report.html")
        with open(mock_report_file,'wt') as fp:
            fp.write("<h1>Test</h1>\n<div id='status' class='hide'>\n<p>Status: WARNINGS</p>\n</div>")
        self.assertTrue(detect_processing_qc_warnings(mock_report_file))
    def test_detect_no_warnings(self):
        """
        detect_processing_qc_warnings: report doesn't contain warnings
        """
        mock_report_file = os.path.join(self.wd,"report.html")
        with open(mock_report_file,'wt') as fp:
            fp.write("<h1>Test</h1>\n<div id='status' class='hide'>\n<p>Status: OK</p>\n</div>")
        self.assertFalse(detect_processing_qc_warnings(mock_report_file))
