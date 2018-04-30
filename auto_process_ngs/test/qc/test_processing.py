#######################################################################
# Tests for qc/processing.py module
#######################################################################

import unittest
import tempfile
import shutil
import os
from auto_process_ngs.auto_processor import AutoProcess
from auto_process_ngs.qc.processing import report_processing_qc

# Set to False to keep test output dirs
REMOVE_TEST_OUTPUTS = False

class TestReportProcessingQC(unittest.TestCase):
    """
    Tests for report_processing_qc function
    """
    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix='TestReportProcessingQC')

    def tearDown(self):
        if REMOVE_TEST_OUTPUTS:
            shutil.rmtree(self.wd)

    def test_report_processing_qc(self):
        """report_processing_qc: standard report
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
        report_processing_qc(AutoProcess(analysis_dir),output_html)
        self.assertTrue(os.path.exists(output_html))

    def test_report_processing_qc_empty_fastqs(self):
        """report_processing_qc: report with empty Fastqs
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
        report_processing_qc(AutoProcess(analysis_dir),output_html)
        self.assertTrue(os.path.exists(output_html))
