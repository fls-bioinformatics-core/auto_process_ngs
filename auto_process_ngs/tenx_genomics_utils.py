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
- has_chromium_indices
- cellranger_info
- make_qc_summary_html
- run_cellranger_mkfastq
- add_cellranger_args

"""

#######################################################################
# Imports
#######################################################################

import os
import re
import json
from bcftbx.IlluminaData import SampleSheet
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.utils import find_program
from .applications import Command
from .simple_scheduler import SchedulerJob
from .docwriter import Document
from .docwriter import List
from .docwriter import Link
from .docwriter import Table
import css_rules

# Initialise logging
import logging
logger = logging.getLogger(__name__)

#######################################################################
# Functions
#######################################################################

def flow_cell_id(run_name):
    """
    Extract the flow cell ID from the run name
    """
    flow_cell_id = os.path.basename(run_name).split("_")[-1]
    return flow_cell_id[1:]

def has_chromium_indices(sample_sheet):
    """
    Check if a sample sheet contains Chromium indices

    The Chromium indices can be obtained from:

    https://support.10xgenomics.com/permalink/27rGqWvNYYuqkgeS66sksm

    The Chromium SC 3'v2 indices are of the form:

    SI-GA-[A-H][1-12]

    e.g. 'SI-GA-B11'

    Arguments:
      sample_sheet (str): path to the sample sheet CSV
        file to check

    Returns:
      Boolean: True if the sample sheet contains at least
        one Chromium index, False if not.
    """
    index_pattern = re.compile(r"SI-GA-[A-H](1[1-2]|[1-9])$")
    s = SampleSheet(sample_sheet)
    for line in s:
        if index_pattern.match(line['index']):
            return True
    return False

def make_qc_summary_html(json_file,html_file):
    """
    Generate HTML report of processing stats from cellranger mkfastqs
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
    # Construct the command
    cmd = Command("cellranger","mkfastq",
                  "--samplesheet",sample_sheet,
                  "--run",primary_data_dir,
                  "--output-dir",output_dir)
    if lanes is not None:
        cmd.add_args("--lanes=%s" % lanes)
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
