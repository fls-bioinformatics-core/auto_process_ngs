#!/usr/bin/env python
#
#     tenx_genomics_utils.py: utility functions for handling 10xGenomics data
#     Copyright (C) University of Manchester 2017 Peter Briggs
#

"""
tenx_genomics_utils.py

Utility classes and functions for processing the outputs from 10xGenomics's
Chromium SC 3'v2 system:

- make_qc_summary_html
"""

#######################################################################
# Imports
#######################################################################

import os
import json
from .docwriter import Document
from .docwriter import List
from .docwriter import Link
from .docwriter import Table
import css_rules

#######################################################################
# Functions
#######################################################################

def flow_cell_id(run_name):
    """
    Extract the flow cell ID from the run name
    """
    flow_cell_id = os.path.basename(run_name).split("_")[-1]
    return flow_cell_id[1:]

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
