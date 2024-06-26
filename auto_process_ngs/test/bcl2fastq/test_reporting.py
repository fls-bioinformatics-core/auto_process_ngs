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
    """
    Tests for the ProcessingQCReport class
    """
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
        """ProcessingQCReport: standard report
        """
        # Create test data
        analysis_dir = os.path.join(self.wd,
                                    "180430_K00311_0001_ABCDEFGHXX_analysis")
        os.mkdir(analysis_dir)
        per_lane_sample_stats = os.path.join(analysis_dir,
                                             "per_lane_sample_stats.info")
        with open(per_lane_sample_stats,'w') as fp:
            fp.write("""
Lane 1
Total reads = 136778255
- AB/AB1	25058003	18.3%
- AB/AB2	22330927	16.3%
- AB/AB3	34509382	25.2%
- AB/AB4	27283286	19.9%
- Undetermined_indices/undetermined	27596657	20.2%

Lane 2
Total reads = 136778255
- CDE/CDE1	25058003	18.3%
- CDE/CDE2	22330927	16.3%
- CDE/CDE3	34509382	25.2%
- CDE/CDE4	27283286	19.9%
- Undetermined_indices/undetermined	27596657	20.2%
""")
        per_lane_statistics = os.path.join(analysis_dir,
                                           "per_lane_statistics.info")
        with open(per_lane_statistics,'w') as fp:
            fp.write("""#Lane	Total reads	Assigned reads	Unassigned reads	%assigned	%unassigned
Lane 1	136778255	109181598	27596657	79.8	20.2
Lane 2	136778255	109181598	27596657	79.8	20.2
""")
        statistics_full = os.path.join(analysis_dir,
                                       "statistics_full.info")
        with open(statistics_full,'w') as fp:
            fp.write("""#Project	Sample	Fastq	Size	Nreads	Paired_end	Read_number	L1	L2
AB	AB1	AB1_S1_R1_001.fastq.gz	1.0G	25058003	Y	1	25058003	
AB	AB1	AB1_S1_R2_001.fastq.gz	1.1G	25058003	Y	2	25058003	
AB	AB2	AB2_S2_R1_001.fastq.gz	941.7M	22330927	Y	1	22330927	
AB	AB2	AB2_S2_R2_001.fastq.gz	1.0G	22330927	Y	2	22330927	
AB	AB3	AB3_S3_R1_001.fastq.gz	1.4G	34509382	Y	1	34509382	
AB	AB3	AB3_S3_R2_001.fastq.gz	1.6G	34509382	Y	2	34509382	
AB	AB4	AB4_S4_R1_001.fastq.gz	1.1G	27283286	Y	1	27283286	
AB	AB4	AB4_S4_R2_001.fastq.gz	1.2G	27283286	Y	2	27283286	
CDE	CDE1	CDE1_S5_R1_001.fastq.gz	1.0G	25058003	Y	1		25058003
CDE	CDE1	CDE1_S5_R2_001.fastq.gz	1.1G	25058003	Y	2		25058003
CDE	CDE2	CDE2_S6_R1_001.fastq.gz	941.7M	22330927	Y	1		22330927
CDE	CDE2	CDE2_S6_R2_001.fastq.gz	1.0G	22330927	Y	2		22330927
CDE	CDE3	CDE3_S7_R1_001.fastq.gz	1.4G	34509382	Y	1		34509382
CDE	CDE3	CDE3_S7_R2_001.fastq.gz	1.6G	34509382	Y	2		34509382
CDE	CDE4	CDE4_S8_R1_001.fastq.gz	1.1G	27283286	Y	1		27283286
CDE	CDE4	CDE4_S8_R2_001.fastq.gz	1.2G	27283286	Y	2		27283286
Undetermined_indices	undetermined	Undetermined_S0_R1_001.fastq.gz	22.0G	27596657	Y	1	27596657	27596657
Undetermined_indices	undetermined	Undetermined_S0_R2_001.fastq.gz	24.0G	27596657	Y	2	27596657	27596657
""")
        seq_len_statistics = os.path.join(analysis_dir,
                                          "seq_len_statistics.info")
        with open(seq_len_statistics,'w') as fp:
            fp.write("""#Fastq	min	max	mean	nreads	masked	padded	masked_frac	padded_frac
AB1_S1_R1_001.fastq.gz	35	76	75.26169877916482	25058003	293	85	0.16576242228118512	0.0480880747232107
AB1_S1_R2_001.fastq.gz	35	76	75.26169877916482	25058003	293	85	0.16576242228118512	0.0480880747232107
AB2_S2_R1_001.fastq.gz	35	76	75.26169877916482	22330927	293	85	0.16576242228118512	0.0480880747232107
AB2_S2_R2_001.fastq.gz	35	76	75.26169877916482	22330927	293	85	0.16576242228118512	0.0480880747232107
AB3_S3_R1_001.fastq.gz	35	76	75.26169877916482	34509382	293	85	0.16576242228118512	0.0480880747232107
AB3_S3_R2_001.fastq.gz	35	76	75.26169877916482	34509382	293	85	0.16576242228118512	0.0480880747232107
AB4_S4_R1_001.fastq.gz	35	76	75.26169877916482	27283286	293	85	0.16576242228118512	0.0480880747232107
AB4_S4_R2_001.fastq.gz	35	76	75.26169877916482	27283286	293	85	0.16576242228118512	0.0480880747232107
CDE1_S5_R1_001.fastq.gz	35	76	75.26169877916482	25058003	293	85	0.16576242228118512	0.0480880747232107
CDE1_S5_R2_001.fastq.gz	35	76	75.26169877916482	25058003	293	85	0.16576242228118512	0.0480880747232107
CDE2_S6_R1_001.fastq.gz	35	76	75.26169877916482	22330927	293	85	0.16576242228118512	0.0480880747232107
CDE2_S6_R2_001.fastq.gz	35	76	75.26169877916482	22330927	293	85	0.16576242228118512	0.0480880747232107
CDE3_S7_R1_001.fastq.gz	35	76	75.26169877916482	34509382	293	85	0.16576242228118512	0.0480880747232107
CDE3_S7_R2_001.fastq.gz	35	76	75.26169877916482	34509382	293	85	0.16576242228118512	0.0480880747232107
CDE4_S8_R1_001.fastq.gz	35	76	75.26169877916482	27283286	293	85	0.16576242228118512	0.0480880747232107
CDE4_S8_R2_001.fastq.gz	35	76	75.26169877916482	27283286	293	85	0.16576242228118512	0.0480880747232107
Undetermined_S0_R1_001.fastq.gz	35	76	75.26169877916482	27596657	293	85	0.16576242228118512	0.0480880747232107
Undetermined_S0_R2_001.fastq.gz	35	76	75.26169877916482	27596657	293	85	0.16576242228118512	0.0480880747232107
""")
        # Generate QC report
        output_html = os.path.join(analysis_dir,
                                   "processing_report.html")
        self.assertFalse(os.path.exists(output_html))
        processing_qc_report = ProcessingQCReport(
            analysis_dir,
            statistics_full,
            per_lane_statistics,
            per_lane_sample_stats,
            seq_len_statistics
        )
        processing_qc_report.write(output_html)
        self.assertTrue(os.path.exists(output_html))
        # Check the HTML
        with open(output_html,'rt') as fp:
            html = fp.read()
        # No warnings
        self.assertTrue(html.find("<p>Status: OK</p>") > -1)

    def test_processing_qc_report_no_seq_len_stats(self):
        """ProcessingQCReport: handle missing sequence length stats
        """
        # Create test data
        analysis_dir = os.path.join(self.wd,
                                    "180430_K00311_0001_ABCDEFGHXX_analysis")
        os.mkdir(analysis_dir)
        per_lane_sample_stats = os.path.join(analysis_dir,
                                             "per_lane_sample_stats.info")
        with open(per_lane_sample_stats,'w') as fp:
            fp.write("""
Lane 1
Total reads = 136778255
- AB/AB1	25058003	18.3%
- AB/AB2	22330927	16.3%
- AB/AB3	34509382	25.2%
- AB/AB4	27283286	19.9%
- Undetermined_indices/undetermined	27596657	20.2%

Lane 2
Total reads = 136778255
- CDE/CDE1	25058003	18.3%
- CDE/CDE2	22330927	16.3%
- CDE/CDE3	34509382	25.2%
- CDE/CDE4	27283286	19.9%
- Undetermined_indices/undetermined	27596657	20.2%
""")
        per_lane_statistics = os.path.join(analysis_dir,
                                           "per_lane_statistics.info")
        with open(per_lane_statistics,'w') as fp:
            fp.write("""#Lane	Total reads	Assigned reads	Unassigned reads	%assigned	%unassigned
Lane 1	136778255	109181598	27596657	79.8	20.2
Lane 2	136778255	109181598	27596657	79.8	20.2
""")
        statistics_full = os.path.join(analysis_dir,
                                       "statistics_full.info")
        with open(statistics_full,'w') as fp:
            fp.write("""#Project	Sample	Fastq	Size	Nreads	Paired_end	Read_number	L1	L2
AB	AB1	AB1_S1_R1_001.fastq.gz	1.0G	25058003	Y	1	25058003	
AB	AB1	AB1_S1_R2_001.fastq.gz	1.1G	25058003	Y	2	25058003	
AB	AB2	AB2_S2_R1_001.fastq.gz	941.7M	22330927	Y	1	22330927	
AB	AB2	AB2_S2_R2_001.fastq.gz	1.0G	22330927	Y	2	22330927	
AB	AB3	AB3_S3_R1_001.fastq.gz	1.4G	34509382	Y	1	34509382	
AB	AB3	AB3_S3_R2_001.fastq.gz	1.6G	34509382	Y	2	34509382	
AB	AB4	AB4_S4_R1_001.fastq.gz	1.1G	27283286	Y	1	27283286	
AB	AB4	AB4_S4_R2_001.fastq.gz	1.2G	27283286	Y	2	27283286	
CDE	CDE1	CDE1_S5_R1_001.fastq.gz	1.0G	25058003	Y	1		25058003
CDE	CDE1	CDE1_S5_R2_001.fastq.gz	1.1G	25058003	Y	2		25058003
CDE	CDE2	CDE2_S6_R1_001.fastq.gz	941.7M	22330927	Y	1		22330927
CDE	CDE2	CDE2_S6_R2_001.fastq.gz	1.0G	22330927	Y	2		22330927
CDE	CDE3	CDE3_S7_R1_001.fastq.gz	1.4G	34509382	Y	1		34509382
CDE	CDE3	CDE3_S7_R2_001.fastq.gz	1.6G	34509382	Y	2		34509382
CDE	CDE4	CDE4_S8_R1_001.fastq.gz	1.1G	27283286	Y	1		27283286
CDE	CDE4	CDE4_S8_R2_001.fastq.gz	1.2G	27283286	Y	2		27283286
Undetermined_indices	undetermined	Undetermined_S0_R1_001.fastq.gz	22.0G	27596657	Y	1	27596657	27596657
Undetermined_indices	undetermined	Undetermined_S0_R2_001.fastq.gz	24.0G	27596657	Y	2	27596657	27596657
""")
        # Generate QC report
        output_html = os.path.join(analysis_dir,
                                   "processing_report.html")
        self.assertFalse(os.path.exists(output_html))
        processing_qc_report = ProcessingQCReport(
            analysis_dir,
            statistics_full,
            per_lane_statistics,
            per_lane_sample_stats
        )
        processing_qc_report.write(output_html)
        self.assertTrue(os.path.exists(output_html))
        # Check the HTML
        with open(output_html,'rt') as fp:
            html = fp.read()
        # No warnings
        self.assertTrue(html.find("<p>Status: OK</p>") > -1)

    def test_processing_qc_report_empty_fastqs(self):
        """ProcessingQCReport: handle empty Fastqs
        """
        # Create test data
        analysis_dir = os.path.join(self.wd,
                                    "180430_K00311_0001_ABCDEFGHXX_analysis")
        os.mkdir(analysis_dir)
        per_lane_sample_stats = os.path.join(analysis_dir,
                                             "per_lane_sample_stats.info")
        with open(per_lane_sample_stats,'w') as fp:
            fp.write("""
Lane 1
Total reads = 79937946
- AB/AB1	25058003	31.4%
- AB/AB2	0	0.0%
- AB/AB3	0	0.0%
- AB/AB4	27283286	34.1%
- Undetermined_indices/undetermined	27596657	34.5%

Lane 2
Total reads = 114447328
- CDE/CDE1	25058003	21.9%
- CDE/CDE2	0	0.0%
- CDE/CDE3	34509382	30.2%
- CDE/CDE4	27283286	23.8%
- Undetermined_indices/undetermined	27596657	24.1%
""")
        per_lane_statistics = os.path.join(analysis_dir,
                                           "per_lane_statistics.info")
        with open(per_lane_statistics,'w') as fp:
            fp.write("""#Lane	Total reads	Assigned reads	Unassigned reads	%assigned	%unassigned
Lane 1	79937946	52341289	27596657	65.5	35.5
Lane 2	114447328	86850671	27596657	75.9	24.1
""")
        statistics_full = os.path.join(analysis_dir,
                                       "statistics_full.info")
        with open(statistics_full,'w') as fp:
            fp.write("""#Project	Sample	Fastq	Size	Nreads	Paired_end	Read_number	L1	L2
AB	AB1	AB1_S1_R1_001.fastq.gz	1.0G	25058003	Y	1	25058003	
AB	AB1	AB1_S1_R2_001.fastq.gz	1.1G	25058003	Y	2	25058003	
AB	AB2	AB2_S2_R1_001.fastq.gz	0.0K	0	Y	1		
AB	AB2	AB2_S2_R2_001.fastq.gz	0.0K	0	Y	2		
AB	AB3	AB3_S3_R1_001.fastq.gz	0.0K	0	Y	1		
AB	AB3	AB3_S3_R2_001.fastq.gz	0.0k	0	Y	2		
AB	AB4	AB4_S4_R1_001.fastq.gz	1.1G	27283286	Y	1	27283286	
AB	AB4	AB4_S4_R2_001.fastq.gz	1.2G	27283286	Y	2	27283286	
CDE	CDE1	CDE1_S5_R1_001.fastq.gz	1.0G	25058003	Y	1		25058003
CDE	CDE1	CDE1_S5_R2_001.fastq.gz	1.1G	25058003	Y	2		25058003
CDE	CDE2	CDE2_S6_R1_001.fastq.gz	0.0K	0	Y	1		
CDE	CDE2	CDE2_S6_R2_001.fastq.gz	0.0K	0	Y	2		
CDE	CDE3	CDE3_S7_R1_001.fastq.gz	1.4G	34509382	Y	1		34509382
CDE	CDE3	CDE3_S7_R2_001.fastq.gz	1.6G	34509382	Y	2		34509382
CDE	CDE4	CDE4_S8_R1_001.fastq.gz	1.1G	27283286	Y	1		27283286
CDE	CDE4	CDE4_S8_R2_001.fastq.gz	1.2G	27283286	Y	2		27283286
Undetermined_indices	undetermined	Undetermined_S0_R1_001.fastq.gz	22.0G	27596657	Y	1	27596657	27596657
Undetermined_indices	undetermined	Undetermined_S0_R2_001.fastq.gz	24.0G	27596657	Y	2	27596657	27596657
""")
        # Generate QC report
        output_html = os.path.join(analysis_dir,
                                   "processing_report.html")
        self.assertFalse(os.path.exists(output_html))
        processing_qc_report = ProcessingQCReport(
            analysis_dir,
            statistics_full,
            per_lane_statistics,
            per_lane_sample_stats
        )
        processing_qc_report.write(output_html)
        self.assertTrue(os.path.exists(output_html))
        # Check the HTML
        with open(output_html,'rt') as fp:
            html = fp.read()
        # Warnings
        self.assertTrue(html.find("<p>Status: WARNINGS</p>") > -1)

    def test_processing_qc_report_empty_lane(self):
        """ProcessingQCReport: handle empty lane
        """
        # Create test data
        analysis_dir = os.path.join(self.wd,
                                    "180430_K00311_0001_ABCDEFGHXX_analysis")
        os.mkdir(analysis_dir)
        per_lane_sample_stats = os.path.join(analysis_dir,
                                             "per_lane_sample_stats.info")
        with open(per_lane_sample_stats,'w') as fp:
            fp.write("""
Lane 1
Total reads = 0

Lane 2
Total reads = 0114447328
- CDE/CDE1	25058003	21.9%
- CDE/CDE2	0	0.0%
- CDE/CDE3	34509382	30.2%
- CDE/CDE4	27283286	23.8%
- Undetermined_indices/undetermined	27596657	24.1%
""")
        per_lane_statistics = os.path.join(analysis_dir,
                                           "per_lane_statistics.info")
        with open(per_lane_statistics,'w') as fp:
            fp.write("""#Lane	Total reads	Assigned reads	Unassigned reads	%assigned	%unassigned
Lane 1	0	0	0	0.0	0.0
Lane 2	114447328	86850671	27596657	75.9	24.1
""")
        statistics_full = os.path.join(analysis_dir,
                                       "statistics_full.info")
        with open(statistics_full,'w') as fp:
            fp.write("""#Project	Sample	Fastq	Size	Nreads	Paired_end	Read_number	L1	L2
AB	AB1	AB1_S1_R1_001.fastq.gz	0.0K	0	Y	1		
AB	AB1	AB1_S1_R2_001.fastq.gz	0.0K	0	Y	2		
AB	AB2	AB2_S2_R1_001.fastq.gz	0.0K	0	Y	1		
AB	AB2	AB2_S2_R2_001.fastq.gz	0.0K	0	Y	2		
AB	AB3	AB3_S3_R1_001.fastq.gz	0.0K	0	Y	1		
AB	AB3	AB3_S3_R2_001.fastq.gz	0.0k	0	Y	2		
AB	AB4	AB4_S4_R1_001.fastq.gz	1.1G	0	Y	1		
AB	AB4	AB4_S4_R2_001.fastq.gz	1.2G	0	Y	2		
CDE	CDE1	CDE1_S5_R1_001.fastq.gz	1.0G	0	Y	1      		25058003
CDE	CDE1	CDE1_S5_R2_001.fastq.gz	1.1G	0	Y	2		25058003
CDE	CDE2	CDE2_S6_R1_001.fastq.gz	0.0K	0	Y	1		
CDE	CDE2	CDE2_S6_R2_001.fastq.gz	0.0K	0	Y	2		
CDE	CDE3	CDE3_S7_R1_001.fastq.gz	1.4G	34509382	Y	1		34509382
CDE	CDE3	CDE3_S7_R2_001.fastq.gz	1.6G	34509382	Y	2		34509382
CDE	CDE4	CDE4_S8_R1_001.fastq.gz	1.1G	27283286	Y	1		27283286
CDE	CDE4	CDE4_S8_R2_001.fastq.gz	1.2G	27283286	Y	2		27283286
Undetermined_indices	undetermined	Undetermined_S0_R1_001.fastq.gz	1.0K	0	Y	1	0	
Undetermined_indices	undetermined	Undetermined_S0_R2_001.fastq.gz	1.0K	0	Y	2	0	
""")
        # Generate QC report
        output_html = os.path.join(analysis_dir,
                                   "processing_report.html")
        self.assertFalse(os.path.exists(output_html))
        processing_qc_report = ProcessingQCReport(
            analysis_dir,
            statistics_full,
            per_lane_statistics,
            per_lane_sample_stats
        )
        processing_qc_report.write(output_html)
        self.assertTrue(os.path.exists(output_html))
        # Check the HTML
        with open(output_html,'rt') as fp:
            html = fp.read()
        # Warnings
        self.assertTrue(html.find("<p>Status: WARNINGS</p>") > -1)

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
