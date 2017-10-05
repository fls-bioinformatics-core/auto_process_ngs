#!/usr/bin/env python
#
#     tenx_genomics_utils.py: utility functions for handling 10xGenomics data
#     Copyright (C) University of Manchester 2017 Peter Briggs
#

"""
tenx_genomics_utils.py

Utility classes and functions for processing the outputs from 10xGenomics's
Chromium SC 3'v2 system:

- flow_cell_id
- has_chromium_sc_indices
- cellranger_info
- make_qc_summary_html
- build_fastq_path_dir
- run_cellranger_mkfastq
- run_cellranger_count
- run_cellranger_count_for_project
- add_cellranger_args

"""

#######################################################################
# Imports
#######################################################################

import os
import re
import json
import shutil
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.utils import find_program
from bcftbx.utils import mkdirs
from .applications import Command
from .bcl2fastq_utils import available_bcl2fastq_versions
from .bcl2fastq_utils import bcl_to_fastq_info
from .simple_scheduler import SchedulerJob
from .simple_scheduler import SimpleScheduler
from .simple_scheduler import SchedulerReporter
from .docwriter import Document
from .docwriter import List
from .docwriter import Link
from .docwriter import Table
from .utils import get_numbered_subdir
from .utils import AnalysisProject
from .utils import ZipArchive
import css_rules

# Initialise logging
import logging
logger = logging.getLogger(__name__)

#######################################################################
# Functions
#######################################################################

def flow_cell_id(run_name):
    """
    Extract the flow cell ID from a run name

    For example for run name "170426_K00311_0033_AHJCY7BBXX"
    the extracted flow cell ID will be "HJCY7BBXX".

    Arguments:
      run_name (str): path to the run name to extract
        flow cell ID from

    Returns:
      String: the extracted flow cell ID.
    """
    flow_cell_id = os.path.basename(run_name).split("_")[-1]
    return flow_cell_id[1:]

def has_chromium_sc_indices(sample_sheet):
    """
    Check if a sample sheet contains Chromium SC indices

    The Chromium SC indices can be obtained from:

    https://support.10xgenomics.com/permalink/27rGqWvNYYuqkgeS66sksm

    The Chromium SC 3'v2 indices are of the form:

    SI-GA-[A-H][1-12]

    e.g. 'SI-GA-B11'

    Arguments:
      sample_sheet (str): path to the sample sheet CSV
        file to check

    Returns:
      Boolean: True if the sample sheet contains at least
        one Chromium SC index, False if not.
    """
    index_pattern = re.compile(r"SI-GA-[A-H](1[1-2]|[1-9])$")
    s = SampleSheet(sample_sheet)
    for line in s:
        if index_pattern.match(line['index']):
            return True
    return False

def make_qc_summary_html(json_file,html_file):
    """
    Make HTML report for cellranger mkfastqs processing stats

    Arguments:
      json_file (str): path to JSON file output from
        cellranger mkfastq command
      html_file (str): path to output HTML file
    """
    # Load the JSON data
    with open(json_file,'r') as fp:
        data = json.load(fp)
    # Initialise the HTML report
    qc_summary = Document("mkfastq QC report")
    qc_summary.add_css_rule(css_rules.QC_REPORT_CSS_RULES)
    qc_summary.add_css_rule("table { font-size: 80%;\n"
                            "        font-family: sans-serif; }")
    qc_summary.add_css_rule("td { text-align: right; }")
    # Add table of contents
    toc = qc_summary.add_section("Contents",name="toc")
    toc_list = List()
    toc.add(toc_list)
    # General information
    general_info = qc_summary.add_section("General information")
    toc_list.add_item(Link("General information",general_info))
    data_items = ['run_id',
                  'experiment_name',
                  '10x_software_version',
                  'bcl2fastq_version',
                  'bcl2fastq_args',
                  'rta_version',]
    tbl = Table(columns=['Parameter','Value'])
    for item in data_items:
        tbl.add_row(Parameter=item,Value=data[item])
    general_info.add(tbl)
    # Get the sample names
    sample_names = data['sample_qc'].keys()
    # Get names of the associated data items
    sample0 = sample_names[0]
    item_names = data['sample_qc'][sample0]['all'].keys()
    # Report QC for each sample in tables
    for sample in sample_names:
        sample_qc = qc_summary.add_section("Sample: %s" % sample)
        toc_list.add_item(Link("Sample: %s" % sample,sample_qc))
        # Set up the table
        tbl = Table(['items',],items="")
        for item in item_names:
            tbl.add_row(items=item)
        # Lanes
        lanes = data['sample_qc'][sample].keys()
        for lane in lanes:
            column = "%s" % lane 
            tbl.append_columns(column)
            # Add the data
            for i,item in enumerate(item_names):
                tbl.set_value(i,column,data['sample_qc'][sample][lane][item])
        # Add to the document
        sample_qc.add(tbl)
    # Write the report
    qc_summary.write(html_file)

def build_fastq_path_dir(project_dir):
    """
    Create directory mimicking output from cellranger mkfastq

    This function creates and populates a 'cellranger mkfastq'
    style 'fastq_path' directory from an autoprocess analysis
    project, which can then be used as input to 'cellranger
    count'.

    The new directory will be called 'cellranger_fastq_path'
    and will created in the project directory, and will be
    populated by links to the Fastq files in the project.

    Arguments:
      project_dir (str): path to the project directory in
        which to create the 'fastq_path' directory

    Returns:
      String: path to the 'cellranger_fastq_path' directory.
    """
    project = AnalysisProject(os.path.basename(project_dir.rstrip(os.sep)),
                              os.path.abspath(project_dir))
    fastq_path_dir = os.path.join(project.dirn,
                                  "cellranger_fastq_path")
    mkdirs(fastq_path_dir)
    mkdirs(os.path.join(fastq_path_dir,"Reports"))
    mkdirs(os.path.join(fastq_path_dir,"Stats"))
    fq_dir = os.path.join(fastq_path_dir,project.name)
    mkdirs(fq_dir)
    for fastq in project.fastqs:
        print fastq
        link_name = os.path.join(fq_dir,os.path.basename(fastq))
        if os.path.exists(link_name):
            logging.warning("%s: already exists" % link_name)
            continue
        target = os.path.relpath(fastq,fq_dir)
        logging.debug("Linking: %s -> %s" % (link_name,target))
        os.symlink(target,link_name)
    return fastq_path_dir

def cellranger_info(path=None):
    """
    Retrieve information on the cellranger software

    If called without any arguments this will locate the first
    cellranger executable that is available on the user's PATH,
    and attempts to extract the version.

    Alternatively if the path to an executable is supplied then
    the version will be determined from that instead.

    If no version is identified then the script path is still
    returned, but without any version info.

    In all cases the package package name will be returned as
    'cellranger'.

    Returns:
      Tuple: tuple consisting of (PATH,PACKAGE,VERSION) where PATH
        is the full path for the cellranger program, PACKAGE is
        'cellranger', and VERSION is the package version. If any
        value can't be determined then it will be returned as an
        empty string.
    """
    # Initialise
    cellranger_path = ''
    package_name = 'cellranger'
    package_version = ''
    # Locate the core script
    if not path:
        cellranger_path = find_program('cellranger')
    else:
        cellranger_path = os.path.abspath(path)
    # Identify the version
    if os.path.basename(cellranger_path) == 'cellranger':
        # Run the program to get the version
        version_cmd = Command(cellranger_path,'--version')
        output = version_cmd.subprocess_check_output()[1]
        for line in output.split('\n'):
            if line.startswith('cellranger'):
                # Extract version from line of the form
                # cellranger  (2.0.1)
                try:
                    package_version = line.split('(')[-1].strip(')')
                except Exception as ex:
                    logging.warning("Unable to get version from '%s': %s" %
                                    (line,ex))
    else:
        # No package supplied or located
        logging.warning("Unable to identify cellranger package "
                        "from '%s'" % cellranger_path)
    # Return what we found
    return (cellranger_path,package_name,package_version)

def run_cellranger_mkfastq(sample_sheet,
                           primary_data_dir,
                           output_dir,
                           lanes=None,
                           bases_mask=None,
                           cellranger_jobmode=None,
                           cellranger_maxjobs=None,
                           cellranger_mempercore=None,
                           cellranger_jobinterval=None,
                           log_dir=None,
                           dry_run=False):
    """
    Wrapper for running 'cellranger mkfastq'

    Runs the 10xGenomics 'cellranger mkfastq' command to
    generate Fastqs from bcl files for Chromium single-cell
    data.

    Arguments:
      sample_sheet (str): path to input samplesheet with
        10xGenomics barcode indices
      primary_data_dir (str): path to the top-level
        directory holding the sequencing data
      output_dir (str): path to the output directory
      lanes (str): optional, specify the subset of lanes
        to process (default is to process all lanes
        in the run)
      bases_mask (str): optional, specify an alternative
        bases mask setting (default is to let cellranger
        determine the bases mask automatically)
      cellranger_jobmode (str): specify the job mode to
        pass to cellranger (default: None)
      cellranger_maxjobs (int): specify the maximum
        number of jobs to pass to cellranger (default:
        None)
      cellranger_mempercore (int): specify the memory
        per core (in Gb) to pass to cellranger (default:
        None)
      cellranger_jobinterval (int): specify the interval
        between launching jobs (in ms) to pass to
        cellranger (default: None)
      log_dir (str): path to a directory to write logs
        (default: current working directory)
      dry_run (bool): if True then only report actions
        that would be performed but don't run anything

    Returns:
      Integer: exit code from the cellranger command.
    """
    # Check we have cellranger
    cellranger = find_program('cellranger')
    if not cellranger:
        raise Exception("No cellranger package found")
    print "Using cellranger %s: %s" % (cellranger_info(cellranger)[-1],
                                       cellranger)
    # Check we have bcl2fastq
    bcl2fastq = find_program('bcl2fastq')
    if not bcl2fastq:
        raise Exception("No bcl2fastq package found")
    bcl2fastq = available_bcl2fastq_versions(
        paths=(os.path.dirname(bcl2fastq),),
        reqs='>=2.17')
    if not bcl2fastq:
        raise Exception("No appropriate bcl2fastq software "
                        "located")
    bcl2fastq = bcl2fastq[0]
    print "Using bcl2fastq %s: %s" % (
        bcl_to_fastq_info(bcl2fastq)[-1],
        bcl2fastq)
    # Construct the command
    cmd = Command("cellranger","mkfastq",
                  "--samplesheet",sample_sheet,
                  "--run",primary_data_dir,
                  "--output-dir",output_dir)
    if lanes is not None:
        cmd.add_args("--lanes=%s" % lanes)
    if bases_mask is not None:
        cmd.add_args("--use-bases-mask=%s" % bases_mask)
    add_cellranger_args(cmd,
                        jobmode=cellranger_jobmode,
                        mempercore=cellranger_mempercore,
                        maxjobs=cellranger_maxjobs,
                        jobinterval=cellranger_jobinterval)
    # Run the command
    print "Running %s" % cmd
    if not dry_run:
        # Make a log directory
        if log_dir is None:
            log_dir = os.getcwd()
        else:
            log_dir = os.path.abspath(log_dir)
        # Submit the job
        cellranger_mkfastq_job = SchedulerJob(
            SimpleJobRunner(join_logs=True),
            cmd.command_line,
            name='cellranger_mkfastq',
            working_dir=os.getcwd(),
            log_dir=log_dir)
        cellranger_mkfastq_job.start()
        try:
            cellranger_mkfastq_job.wait()
        except KeyboardInterrupt,ex:
            logging.warning("Keyboard interrupt, terminating cellranger")
            cellranger_mkfastq_job.terminate()
            raise ex
        exit_code = cellranger_mkfastq_job.exit_code
        print "cellranger mkfastq completed: exit code %s" % exit_code
        if exit_code != 0:
            logging.error("cellranger mkfastq exited with an error")
            return exit_code
        # Deal with the QC summary report
        flow_cell_dir = flow_cell_id(primary_data_dir)
        if lanes is not None:
            lanes_suffix = "_%s" % lanes.replace(',','')
        else:
            lanes_suffix = ""
        flow_cell_dir = "%s%s" % (flow_cell_dir,lanes_suffix)
        if not os.path.isdir(flow_cell_dir):
            logging.error("No output directory '%s'" % flow_cell_dir)
            return -1
        json_file = os.path.join(flow_cell_dir,
                                 "outs",
                                 "qc_summary.json")
        html_file = "cellranger_qc_summary%s.html" % lanes_suffix
        make_qc_summary_html(json_file,html_file)
        return exit_code

def run_cellranger_count(fastq_dir,
                         transcriptome,
                         cellranger_jobmode='sge',
                         cellranger_maxjobs=None,
                         cellranger_mempercore=None,
                         cellranger_jobinterval=None,
                         max_jobs=4,
                         log_dir=None,
                         dry_run=False,
                         summary_only=True):
    """
    Wrapper for running 'cellranger count'

    Runs the 10xGenomics 'cellranger count' command to
    perform single library analysis on Fastqs from
    Chromium single-cell samples.

    If the supplied 'fastq_dir' is a 'cellranger mkfastq'
    or 'bcl2fastq' output directory then the analysis
    will be run for each of the projects.

    Arguments:
      fastq_dir (str): path of the 'fastq_path' folder
        from 'cellranger mkfastq', or the output folder
        from 'bcl2fastq' (or with a similar structure),
        or any folder containing Fastq files
      transcriptome (str): path to the cellranger
        compatible transcriptome reference data
        directory
      cellranger_jobmode (str): specify the job mode to
        pass to cellranger (default: None)
      cellranger_maxjobs (int): specify the maximum
        number of jobs to pass to cellranger (default:
        None)
      cellranger_mempercore (int): specify the memory
        per core (in Gb) to pass to cellranger (default:
        None)
      cellranger_jobinterval (int): specify the interval
        between launching jobs (in ms) to pass to
        cellranger (default: None)
      max_jobs (int):
      log_dir (str): path to a directory to write logs
        (default: current working directory)
      dry_run (bool): if True then only report actions
        that would be performed but don't run anything
      summary_only (bool): if True then only collect
        the output 'web_summary.html' file, otherwise
        copy all outputs (warning: this can be very
        large)

    Returns:
      Integer: exit code from the cellranger command.
    """
    # Input data
    sample_names = {}
    try:
        illumina_data = IlluminaData(os.getcwd(),
                                     unaligned_dir=fastq_dir)
        for project in illumina_data.projects:
            sample_names[project.name] = []
            for sample in project.samples:
                sample_names[project.name].append(sample.name)
    except IlluminaDataError:
        logging.critical("Couldn't load data from '%s'" %
                         fastq_dir)
        return 1
    print "Samples: %s" % sample_names
    projects = sample_names.keys()

    # Set up a scheduler
    sched_reporter = SchedulerReporter(
        job_start="SCHEDULER: Started  #%(job_number)d: %(job_name)s:\n-- %(command)s",
        job_end=  "SCHEDULER: Finished #%(job_number)d: %(job_name)s"
    )
    sched_reporter = SchedulerReporter()
    sched = SimpleScheduler(max_concurrent=max_jobs,
                            reporter=sched_reporter)
    sched.start()

    # Make a log directory
    if not dry_run:
        if log_dir is None:
            log_dir = os.getcwd()
        log_dir = get_numbered_subdir("cellranger_count",
                                      parent_dir=log_dir,
                                      full_path=True)
        mkdirs(log_dir)

    # Submit the cellranger count jobs
    jobs = []
    for project in projects:
        print "Project: %s" % project
        for sample in sample_names[project]:
            print "Sample: %s" % sample
            # Check if outputs already exist
            count_dir = os.path.abspath(
                os.path.join(project,
                             "cellranger_count",
                             sample,
                             "outs"))
            if os.path.isdir(count_dir):
                print "-- %s: outputs exist, nothing to do" % sample
                continue
            else:
                print "-- %s: setting up cellranger count" % sample
            # Set up job for this sample
            work_dir = os.path.abspath("tmp.cellranger_count.%s.%s" %
                                       (project,sample))
            mkdirs(work_dir)
            cmd = Command("cellranger","count",
                          "--id",sample,
                          "--fastqs",os.path.abspath(fastq_dir),
                          "--sample",sample,
                          "--transcriptome",transcriptome)
            add_cellranger_args(cmd,
                                jobmode=cellranger_jobmode,
                                mempercore=cellranger_mempercore,
                                maxjobs=cellranger_maxjobs,
                                jobinterval=cellranger_jobinterval)
            print "Running: %s" % cmd
            if not dry_run:
                job = sched.submit(cmd,
                                   name="cellranger_count.%s.%s" %
                                   (project,
                                    sample),
                                   log_dir=log_dir,
                                   wd=work_dir)
                jobs.append(job)
    sched.wait()
    sched.stop()

    # If dry run then stop here
    if dry_run:
        return 0

    # Finished, check the exit status
    retval = 0
    for job in jobs:
        retval += job.exit_code
    if retval != 0:
        logging.critical("One or more jobs finished with non-zero "
                         "exit code")
        return retval

    # Handle outputs
    for project in projects:
        print "Project: %s" % project
        for sample in sample_names[project]:
            print "Sample: %s" % sample
            # Destination for count output
            count_dir = os.path.abspath(
                os.path.join(project,
                             "cellranger_count",
                             sample))
            mkdirs(count_dir)
            # Copy the cellranger count outputs
            outs_dir = os.path.join("tmp.cellranger_count.%s.%s"
                                    % (project,sample),
                                    sample,
                                    "outs")
            if not summary_only:
                # Collect all outputs
                print "Copying contents of %s to %s" % (outs_dir,count_dir)
                shutil.copytree(outs_dir,count_dir)
            else:
                # Only collect the web summary
                count_dir = os.path.join(count_dir,"outs")
                mkdirs(count_dir)
                print "Copying web_summary.html from %s to %s" % (outs_dir,count_dir)
                shutil.copy(os.path.join(outs_dir,"web_summary.html"),count_dir)

    # Create a report and zip archive for each project
    pwd = os.getcwd()
    analysis_dir = os.path.basename(pwd)
    for project in projects:
        # Descend into project dir
        os.chdir(project)
        # Set up zip file
        report_zip = os.path.join("cellranger_count_report.%s.%s.zip" %
                                  (project,analysis_dir))
        zip_file = ZipArchive(report_zip,
                              prefix="cellranger_count_report.%s.%s" %
                              (project,analysis_dir))
        # Construct index page
        print "Making report for project %s" % project
        count_report = Document("%s: cellranger count" % project)
        count_report.add_css_rule(css_rules.QC_REPORT_CSS_RULES)
        summaries = count_report.add_section()
        summaries.add("Reports from cellranger count for each sample:")
        summary_links = List()
        for sample in sample_names[project]:
            # Link to summary for sample
            web_summary = os.path.join("cellranger_count",
                                       sample,
                                       "outs",
                                       "web_summary.html")
            print "Adding web summary (%s) for %s" % (web_summary,
                                                      sample)
            summary_links.add_item(Link("%s" % sample,
                                        web_summary))
            # Add to the zip file
            zip_file.add_file(web_summary)
        summaries.add(summary_links)
        # Write the report and add to the zip file
        html_file = "cellranger_count_report.html"
        count_report.write(html_file)
        zip_file.add_file(html_file)
        # Finish
        zip_file.close()
        os.chdir(pwd)
    # Done
    return retval

def run_cellranger_count_for_project(project_dir,
                                     transcriptome,
                                     cellranger_jobmode='sge',
                                     cellranger_maxjobs=None,
                                     cellranger_mempercore=None,
                                     cellranger_jobinterval=None,
                                     max_jobs=4,
                                     log_dir=None,
                                     dry_run=False,
                                     summary_only=True):
    """
    Wrapper to run 'cellranger count' on a project

    Runs the 10xGenomics 'cellranger count' command to
    perform single library analysis on Fastqs from
    Chromium single-cell samples in an autoprocess
    analysis project directory.

    Arguments:
      project_dir (str): path to the analysis project
        containing Fastq files
      transcriptome (str): path to the cellranger
        compatible transcriptome reference data
        directory
      cellranger_jobmode (str): specify the job mode to
        pass to cellranger (default: None)
      cellranger_maxjobs (int): specify the maximum
        number of jobs to pass to cellranger (default:
        None)
      cellranger_mempercore (int): specify the memory
        per core (in Gb) to pass to cellranger (default:
        None)
      cellranger_jobinterval (int): specify the interval
        between launching jobs (in ms) to pass to
        cellranger (default: None)
      max_jobs (int):
      log_dir (str): path to a directory to write logs
        (default: current working directory)
      dry_run (bool): if True then only report actions
        that would be performed but don't run anything
      summary_only (bool): if True then only collect
        the output 'web_summary.html' file, otherwise
        copy all outputs (warning: this can be very
        large)

    Returns:
      Integer: exit code from the cellranger command.
    """
    # Build the 'fastq_path' dir for cellranger
    fastq_dir = build_fastq_path_dir(project_dir)
    # Run the count procedure
    return run_cellranger_count(
        fastq_dir,
        transcriptome,
        cellranger_jobmode=cellranger_jobmode,
        cellranger_maxjobs=cellranger_maxjobs,
        cellranger_mempercore=cellranger_mempercore,
        cellranger_jobinterval=cellranger_jobinterval,
        max_jobs=max_jobs,
        log_dir=log_dir,
        dry_run=dry_run,
        summary_only=summary_only)

def add_cellranger_args(cellranger_cmd,
                        jobmode=None,
                        maxjobs=None,
                        mempercore=None,
                        jobinterval=None):
    """
    Configure options for cellranger

    Given a Command instance for running cellranger,
    add the appropriate options (e.g. --jobmode)
    according to the supplied arguments.

    Arguments:
      cellranger_cmd (Command): Command instance for
        running cellranger
      jobmode (str): if specified, will be passed to the
        --jobmode option
      maxjobs (int): if specified, will be passed to the
        --mempercore option
      mempercore (int): if specified, will be passed to
        the --maxjobs option
      jobinterval (int):  if specified, will be passed to
        the --jobinterval option

    Returns:
      Command: the original command updated with the
        appropriate options.
    """
    if jobmode is not None:
        cellranger_cmd.add_args("--jobmode=%s" % jobmode)
    if mempercore is not None:
        cellranger_cmd.add_args("--mempercore=%s" % mempercore)
    if maxjobs is not None:
        cellranger_cmd.add_args("--maxjobs=%s" % maxjobs)
    if jobinterval is not None:
        cellranger_cmd.add_args("--jobinterval=%s" % jobinterval)
    return cellranger_cmd

import sys
if __name__ == "__main__":
    json_file = sys.argv[1]
    html_file = "qc_summary.html"
    report_cellranger_mkfastqs_qc_summary(json_file,html_file)
