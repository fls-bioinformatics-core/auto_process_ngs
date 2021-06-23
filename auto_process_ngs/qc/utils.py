#!/usr/bin/env python
#
#     utils: utility classes and functions for QC
#     Copyright (C) University of Manchester 2018-2021 Peter Briggs
#
"""
Provides utility classes and functions for analysis project QC.

Provides the following functions:

- determine_qc_protocol: get QC protocol for a project
- verify_qc: verify the QC run for a project
- report_qc: generate report for the QC run for a project
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from ..command import Command
from ..settings import Settings
from ..simple_scheduler import SchedulerJob

# Module-specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Functions
#######################################################################

def determine_qc_protocol(project):
    """
    Determine the QC protocol for a project

    Arguments:
      project (AnalysisProject): project instance

    Return:
      String: QC protocol for the project
    """
    # Standard protocols
    if project.info.paired_end:
        protocol = "standardPE"
    else:
        protocol = "standardSE"
    # Single cell protocols
    if project.info.single_cell_platform is not None:
        # Default
        protocol = "singlecell"
        single_cell_platform = project.info.single_cell_platform
        library_type = project.info.library_type
        if single_cell_platform.startswith('10xGenomics Chromium 3\''):
            if library_type == "scRNA-seq":
                # 10xGenomics scATAC-seq
                protocol = "10x_scRNAseq"
            elif library_type == "snRNA-seq":
                # 10xGenomics snRNA-seq
                protocol = "10x_snRNAseq"
        elif library_type in ("scATAC-seq",
                              "snATAC-seq",):
            if single_cell_platform == "10xGenomics Single Cell ATAC":
                # 10xGenomics scATAC-seq
                protocol = "10x_scATAC"
            elif single_cell_platform == "ICELL8":
                # ICELL8 scATAC-seq
                protocol = "ICELL8_scATAC"
    # Spatial RNA-seq
    if project.info.single_cell_platform == "10xGenomics Visium":
        # 10xGenomics Visium spatial transcriptomics
        protocol = "10x_Visium"
    # Multiome ATAC+GEX
    if project.info.single_cell_platform == "10xGenomics Single Cell Multiome":
        if library_type == "ATAC":
            # 10xGenomics single cell Multiome ATAC
            protocol = "10x_Multiome_ATAC"
        elif library_type == "GEX":
            # 10xGenomics single cell Multiome gene expression
            protocol = "10x_Multiome_GEX"
    return protocol

def verify_qc(project,qc_dir=None,fastq_dir=None,qc_protocol=None,
              runner=None,log_dir=None):
    """
    Verify the QC run for a project

    Arguments:
      project (AnalysisProject): analysis project
        to verify the QC for
      qc_dir (str): optional, specify the subdir with
        the QC outputs being verified
      fastq_dir (str): optional, specify a non-default
        directory with Fastq files being verified
      qc_protocol (str): optional, QC protocol to
        verify against
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
    # Construct command for QC verification
    verify_cmd = Command(
        "reportqc.py",
        "--verify")
    if qc_protocol is not None:
        verify_cmd.add_args("--protocol",qc_protocol)
    if qc_dir is not None:
        verify_cmd.add_args("--qc_dir",qc_dir)
    if fastq_dir is not None:
        verify_cmd.add_args("--fastq_dir",fastq_dir)
    verify_cmd.add_args(project.dirn)
    # Run the command
    verify = SchedulerJob(runner,
                          verify_cmd.command_line,
                          name="verify_qc.%s" % project.name,
                          working_dir=project.dirn,
                          log_dir=log_dir)
    verify.start()
    try:
        verify.wait()
    except KeyboardInterrupt as ex:
        logger.warning("Keyboard interrupt, terminating QC verification")
        verify.terminate()
        raise ex
    # Return boolean based on the exit code
    return (verify.exit_code == 0)

def report_qc(project,qc_dir=None,fastq_dir=None,qc_protocol=None,
              report_html=None,zip_outputs=True,multiqc=False,
              force=False,runner=None,log_dir=None):
    """
    Generate report for the QC run for a project

    Arguments:
      project (AnalysisProject): analysis project
        to report the QC for
      qc_dir (str): optional, specify the subdir with
        the QC outputs being reported
      fastq_dir (str): optional, specify a non-default
        directory with Fastq files being verified
      qc_protocol (str): optional, QC protocol to
        verify against
      report_html (str): optional, path to the name of
        the output QC report
      zip_outputs (bool): if True then also generate ZIP
        archive with the report and QC outputs
      multiqc (bool): if True then also generate MultiQC
        report
      force (bool): if True then force generation of
        QC report even if verification fails
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
    # Basename for the outputs
    if qc_dir is None:
        qc_base = os.path.basename(project.qc_dir)
    else:
        qc_base = os.path.basename(qc_dir)
    # Report HTML file name
    if report_html is None:
        out_file = '%s_report.html' % qc_base
    else:
        out_file = report_html
    if not os.path.isabs(out_file):
        out_file = os.path.join(project.dirn,out_file)
    # Report title
    if project.info.run is None:
        title = "%s" % project.name
    else:
        title = "%s/%s" % (project.info.run,
                           project.name)
    if fastq_dir is not None:
        title = "%s (%s)" % (title,fastq_dir)
    title = "%s: QC report" % title
    # Construct command for reporting
    report_cmd = Command(
        "reportqc.py",
        "--filename",out_file,
        "--title",title)
    if qc_protocol is not None:
        report_cmd.add_args("--protocol",qc_protocol)
    if qc_dir is not None:
        report_cmd.add_args("--qc_dir",qc_dir)
    if fastq_dir is not None:
        report_cmd.add_args("--fastq_dir",fastq_dir)
    if multiqc:
        report_cmd.add_args("--multiqc")
    if zip_outputs:
        report_cmd.add_args("--zip")
    if force:
        report_cmd.add_args("--force")
    report_cmd.add_args(project.dirn)
    # Run the command
    report = SchedulerJob(runner,
                          report_cmd.command_line,
                          name="report_qc.%s" % project.name,
                          working_dir=project.dirn,
                          log_dir=log_dir)
    report.start()
    try:
        report.wait()
    except KeyboardInterrupt as ex:
        logger.warning("Keyboard interrupt, terminating QC reporting")
        report.terminate()
        raise ex
    # Return the exit code
    return report.exit_code
