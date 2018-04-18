#!/usr/bin/env python
#
#     run_qc_cmd.py: implement auto process run_qc command
#     Copyright (C) University of Manchester 2018 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import logging
from auto_process_ngs.applications import Command
from auto_process_ngs.qc.runqc import RunQC
from bcftbx.JobRunner import fetch_runner

# Fetch configuration settings
import auto_process_ngs.settings
__settings = auto_process_ngs.settings.Settings()

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Command functions
#######################################################################

def run_qc(ap,projects=None,max_jobs=4,ungzip_fastqs=False,
           fastq_screen_subset=100000,nthreads=1,
           runner=None,fastq_dir=None,qc_dir=None,
           report_html=None,run_multiqc=True):
    """Run QC pipeline script for projects

    Run the illumina_qc.sh script to perform QC on projects.

    Note that if all QC outputs already exist for a project then
    the QC will *not* be run for that project.

    A subset of projects can be selected for QC by setting the
    'projects' argument to a name or pattern, only matching
    projects will be examined.

    Arguments:
      projects (str): specify a pattern to match one or more
        projects to run the QC for (default is to run QC for all
        projects)
      max_jobs (int): maximum number of jobs that will be
        scheduled to run at one time (passed to the scheduler;
        default is 4, set to zero to remove the limit)
      ungzip_fastqs (bool): if True then run the QC script with
        the '--ungzip-fastqs' option to create decompressed
        copies of any fastq.gz inputs (default: False i.e. don't
        decompress the input files)
      fastq_screen_subset (int): subset of reads to use in
        FastQScreen, set to zero or None to use all reads
        (default: 100000)
      nthreads (int): specify number of threads to run the QC jobs
        with (default: 1)
      runner (str): specify a non-default job runner to use for
        the QC jobs
      fastq_dir (str): specify the subdirectory to take the
        Fastq files from; will be used for all projects that are
        processed (default: 'fastqs')
      qc_dir (str): specify a non-standard directory to write the
        QC outputs to; will be used for all projects that are
        processed (default: 'qc')
      report_html (str): specify the name for the output HTML QC
        report (default: '<QC_DIR>_report.html')
      run_multiqc (bool): if True then run MultiQC at the end of
        the QC run (default)

    Returns:
      Integer: UNIX-style integer returncode where 0 = successful
        termination, non-zero indicates an error.
    """
    # Check QC script version
    compatible_versions = ('1.3.0','1.3.1')
    print "Getting QC script information"
    status,qc_script_info = Command('illumina_qc.sh',
                                    '--version').subprocess_check_output()
    print "Using QC script %s" % qc_script_info.strip()
    version = qc_script_info.strip().split()[-1]
    if version not in compatible_versions:
        logger.error("QC script version is %s, needs %s" %
                     (version,'/'.join(compatible_versions)))
        return 1
    # Process project pattern matching
    if projects is None:
        project_pattern = '*'
        sample_pattern = '*'
    else:
        project_pattern = projects.split('/')[0]
        try:
            sample_pattern = projects.split('/')[1]
        except IndexError:
            sample_pattern = '*'
    # Get project dir data
    projects = ap.get_analysis_projects(project_pattern)
    # Check we have projects
    if len(projects) == 0:
        logger.warning("No projects found for QC analysis")
        return 1
    # Set up runners
    default_runner = __settings.general.default_runner
    if runner is not None:
        qc_runner = fetch_runner(runner)
    else:
        qc_runner = __settings.runners.qc
    # Set up the QC for each project
    runqc = RunQC()
    for project in projects:
        runqc.add_project(project,
                          fastq_dir=fastq_dir,
                          sample_pattern=sample_pattern,
                          qc_dir=qc_dir,
                          ungzip_fastqs=ungzip_fastqs)
    # Run the QC
    status = runqc.run(multiqc=True,
                       qc_runner=qc_runner,
                       verify_runner=default_runner,
                       report_runner=default_runner,
                       max_jobs=max_jobs)
    return status
