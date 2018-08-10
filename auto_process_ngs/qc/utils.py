#!/usr/bin/env python
#
#     utils: utility classes and functions for QC
#     Copyright (C) University of Manchester 2018 Peter Briggs
#
"""
Provides utility classes and functions for analysis project QC.

Provides the following functions:

- verify_qc: verify the QC run for a project
- report_qc: generate report for the QC run for a project
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
import uuid
import tempfile
import shutil
from .runqc import ProjectQC
from auto_process_ngs.settings import Settings
from auto_process_ngs.simple_scheduler import SimpleScheduler

# Module-specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Functions
#######################################################################

def verify_qc(project,qc_dir=None,illumina_qc=None,runner=None,
              log_dir=None):
    """
    Verify the QC run for a project

    Arguments:
      project (AnalysisProject): analysis project
        to verify the QC for
      qc_dir (str): optional, specify the subdir with
        the QC outputs being verified
      illumina_qc (IlluminaQC): optional, configured
        IlluminaQC object to use for QC verification
      runner (JobRunner): optional, job runner to use
        for running the verification
      log_dir (str): optional, specify a directory to
        write logs to

    Returns:
      Boolean: True if QC passes verification, otherwise
        False.
    """
    # Sort out runners
    if runner is None:
        runner = Settings().general.default_runner
    # Set up QC project
    project = ProjectQC(project,
                        illumina_qc=illumina_qc,
                        qc_dir=qc_dir,
                        log_dir=log_dir)
    # Set up and start scheduler
    sched = SimpleScheduler()
    sched.start()
    # QC check for project
    project.check_qc(sched,
                     name="verify_qc",
                     runner=runner)
    sched.wait()
    return project.verify()

def report_qc(project,qc_dir=None,illumina_qc=None,
              report_html=None,zip_outputs=True,multiqc=False,
              runner=None,log_dir=None):
    """
    Generate report for the QC run for a project

    Arguments:
      project (AnalysisProject): analysis project
        to report the QC for
      qc_dir (str): optional, specify the subdir with
        the QC outputs being reported
      illumina_qc (IlluminaQC): optional, configured
        IlluminaQC object to use for QC reporting
      report_html (str): optional, path to the name of
        the output QC report
      zip_outputs (bool): if True then also generate ZIP
        archive with the report and QC outputs
      multiqc (bool): if True then also generate MultiQC
        report
      runner (JobRunner): optional, job runner to use
        for running the reporting
      log_dir (str): optional, specify a directory to
        write logs to

    Returns:
      Integer: exit code from reporting job (zero indicates
        success, non-zero indicates a problem).
    """
    # Sort out runners
    if runner is None:
        runner = Settings().general.default_runner
    # Set up QC project
    project = ProjectQC(project,
                        illumina_qc=illumina_qc,
                        qc_dir=qc_dir,
                        log_dir=log_dir)
    # Set up and start scheduler
    sched = SimpleScheduler()
    sched.start()
    # Generate QC for project
    project.report_qc(sched,
                      report_html=report_html,
                      multiqc=multiqc,
                      zip_outputs=zip_outputs,
                      runner=runner)
    sched.wait()
    return project.reporting_status
