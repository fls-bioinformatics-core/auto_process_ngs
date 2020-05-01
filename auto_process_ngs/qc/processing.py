#!/usr/bin/env python
#
# Library for QC reporting of processing stages

import os
from ..bcl2fastq.reporting import ProcessingQCReport

#######################################################################
# Functions
#######################################################################

def get_absolute_file_path(p,base=None):
    """
    Get absolute path for supplied path

    Arguments:
      p (str): path
      base (str): optional, base directory to
        use if p is relative

    Returns:
      String: absolute path for p.
    """
    if not os.path.isabs(p):
        if base is not None:
            p = os.path.join(os.path.abspath(base),p)
        else:
            p = os.path.abspath(p)
    return p

def report_processing_qc(ap,html_file):
    """
    Generate HTML report for processing statistics

    Arguments:
      ap (AutoProcess): AutoProcess instance to report the
        processing from
      html_file (str): destination path and file name for
        HTML report
    """
    # Per-lane statistics
    per_lane_stats_file = ap.params.per_lane_stats_file
    if per_lane_stats_file is None:
        per_lane_stats_file = "per_lane_statistics.info"
    per_lane_stats_file = get_absolute_file_path(per_lane_stats_file,
                                                 base=ap.analysis_dir)
    # Per lane by sample statistics
    per_lane_sample_stats_file = get_absolute_file_path(
        "per_lane_sample_stats.info",
        base=ap.analysis_dir)
    # Per fastq statistics
    stats_file = get_absolute_file_path("statistics_full.info",
                                        base=ap.analysis_dir)
    if not os.path.exists(stats_file):
        if ap.params.stats_file is not None:
            stats_file = ap.params.stats_file
        else:
            stats_file = "statistics.info"
    stats_file = get_absolute_file_path(stats_file,
                                        base=ap.analysis_dir)
    # Generate the report
    ProcessingQCReport(ap.analysis_dir,
                       stats_file,
                       per_lane_stats_file,
                       per_lane_sample_stats_file).write(html_file)

def detect_processing_qc_warnings(html_file):
    """
    Look for warning text in processing_qc.html file

    Arguments:
      html_file (str): path to HTML report file

    Returns:
      Boolean: True if warnings were found, False if not.
    """
    with open(html_file) as fp:
        for line in fp:
            if "Status: WARNINGS" in line:
                return True
                break
    return False
