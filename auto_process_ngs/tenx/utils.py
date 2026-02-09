#!/usr/bin/env python
#
#     tenx/utils.py: utility functions for handling 10xGenomics data
#     Copyright (C) University of Manchester 2023-2026 Peter Briggs
#

"""
Utility functions for processing the outputs from 10x Genomics
pipelines:

- flow_cell_id
- has_10x_indices
- has_chromium_sc_indices
- cellranger_info
- spaceranger_info
- make_qc_summary_html
- add_cellranger_args
- make_multi_config_template
"""

#######################################################################
# Imports
#######################################################################

import os
import re
import json
from textwrap import dedent
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import split_run_name_full
from bcftbx.utils import find_program
from ..command import Command
from ..docwriter import Document
from ..docwriter import List
from ..docwriter import Link
from ..docwriter import Table
from .. import css_rules

# Default version of cellranger
from . import DEFAULT_CELLRANGER_VERSION

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
    ds,inst,run,prefix,flow_cell_id = split_run_name_full(run_name)
    return flow_cell_id

def has_10x_indices(sample_sheet):
    """
    Check if a sample sheet contains 10xGenomics-format indices

    The Chromium SC 3'v2 indices are of the form:

    SI-GA-[A-H][1-12]

    e.g. 'SI-GA-B11' (see
    https://support.10xgenomics.com/permalink/27rGqWvNYYuqkgeS66sksm)

    For scATAC-seq the indices are assumed to be of the form:

    SI-NA-[A-H][1-12]

    e.g. 'SI-NA-G9'

    For Visium data the indices are assumed to be of the form:

    SI-(TT|TS)-[A-H][1-12]

    e.g. 'SI-TT-B1'

    Arguments:
      sample_sheet (str): path to the sample sheet CSV
        file to check

    Returns:
      Boolean: True if the sample sheet contains at least
        one 10xGenomics-style index, False if not.
    """
    index_pattern = re.compile(r"SI-(GA|NA|TT|TS)-[A-H](1[0-2]|[1-9])$")
    s = SampleSheet(sample_sheet)
    for line in s:
        try:
            if index_pattern.match(line['index']):
                return True
        except KeyError:
            pass
    return False

def has_chromium_sc_indices(sample_sheet):
    """
    Wrapper for 'has_10x_indices'.

    Maintained for backwards compatibility
    """
    return has_10x_indices(sample_sheet)

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
    sample_names = list(data['sample_qc'].keys())
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

def cellranger_info(path=None,name=None):
    """
    Retrieve information on the cellranger software

    If called without any arguments this will locate the first
    cellranger executable that is available on the user's PATH,
    and attempts to extract the version.

    Alternatively if the path to an executable is supplied then
    the version will be determined from that instead.

    If no version is identified then the script path is still
    returned, but without any version info.

    If a 'path' is supplied then the package name will be taken
    from the basename; otherwise the package name can be supplied
    via the 'name' argument. If neither are supplied then the
    package name defaults to 'cellranger'.

    Returns:
      Tuple: tuple consisting of (PATH,PACKAGE,VERSION) where PATH
        is the full path for the cellranger program, PACKAGE is
        'cellranger', and VERSION is the package version. If any
        value can't be determined then it will be returned as an
        empty string.
    """
    # Initialise
    cellranger_path = ''
    if name is None:
        if path:
            name = os.path.basename(path)
        else:
            name = 'cellranger'
    package_name = name
    package_version = ''
    # Locate the core script
    if not path:
        cellranger_path = find_program(package_name)
    else:
        cellranger_path = os.path.abspath(path)
    # Identify the version
    if os.path.basename(cellranger_path) == package_name:
        # Run the program to get the version
        version_cmd = Command(cellranger_path,'--version')
        output = version_cmd.subprocess_check_output()[1]
        for line in output.split('\n'):
            if package_name in ('cellranger',
                                'cellranger-atac',):
                try:
                    if line.startswith("%s %s-" % (package_name,
                                                   package_name)):
                        # Extract version from line of the form
                        # cellranger cellranger-5.0.1
                        package_version = line.split('-')[-1]
                    elif line.startswith("%s " % package_name) and \
                         line.endswith(")"):
                        # Extract version from line of the form
                        # cellranger ... (2.0.1)
                        package_version = line.split('(')[-1].strip(')')
                    elif line.startswith("%s " % package_name):
                        # Extract version from line of the form
                        # cellranger 10.0.0
                        package_version = line.split(' ')[-1]
                    else:
                        # Raise an exception
                        raise("unrecognised version format")
                except Exception as ex:
                    logger.warning("Unable to get version from '%s': "
                                   "%s" % (line,ex))
            elif package_name == 'cellranger-arc':
                # Extract version from line of the form
                # cellranger-arc cellranger-arc-1.0.0
                try:
                    package_version = line.split('-')[-1]
                except Exception as ex:
                    logger.warning("Unable to get version from '%s': "
                                   "%s" % (line,ex))
            if package_version:
                # Acquired version, stop processing lines
                break
    else:
        # No package supplied or located
        logger.warning("Unable to identify %s package from '%s'" %
                       (name,cellranger_path))
    # Return what we found
    return (cellranger_path,package_name,package_version)

def spaceranger_info(path=None,name=None):
    """
    Retrieve information on the spaceranger software

    If called without any arguments this will locate the first
    spaceranger executable that is available on the user's PATH,
    and attempts to extract the version.

    Alternatively if the path to an executable is supplied then
    the version will be determined from that instead.

    If no version is identified then the script path is still
    returned, but without any version info.

    If a 'path' is supplied then the package name will be taken
    from the basename; otherwise the package name can be supplied
    via the 'name' argument. If neither are supplied then the
    package name defaults to 'cellranger'.

    Returns:
      Tuple: tuple consisting of (PATH,PACKAGE,VERSION) where PATH
        is the full path for the spaceranger program, PACKAGE is
        'spaceranger', and VERSION is the package version. If any
        value can't be determined then it will be returned as an
        empty string.
    """
    # Initialise
    spaceranger_path = ''
    if name is None:
        if path:
            name = os.path.basename(path)
        else:
            name = 'spaceranger'
    package_name = name
    package_version = ''
    # Locate the core script
    if not path:
        spaceranger_path = find_program(package_name)
    else:
        spaceranger_path = os.path.abspath(path)
    # Identify the version
    if os.path.basename(spaceranger_path) == package_name:
        # Run the program to get the version
        version_cmd = Command(spaceranger_path,'--version')
        output = version_cmd.subprocess_check_output()[1]
        for line in output.split('\n'):
            try:
                if line.startswith("%s %s-" % (package_name,
                                               package_name)):
                    # Extract version from line of the form
                    # spaceranger spaceranger-1.3.1
                    package_version = line.split('-')[-1]
                elif line.startswith("%s " % package_name):
                    # Extract version from line of the from
                    # spaceranger 1.1.0
                    package_version = line.split()[-1]
                else:
                    # Raise an exception
                    raise("unrecognised version format")
            except Exception as ex:
                logger.warning("Unable to get version from '%s': %s" %
                               (output,ex))
            if package_version:
                # Acquired version, stop processing lines
                break
    else:
        # No package supplied or located
        logger.warning("Unable to identify spaceranger package "
                       "from '%s'" % spaceranger_path)
    # Return what we found
    return (spaceranger_path,package_name,package_version)

def add_cellranger_args(cellranger_cmd,
                        jobmode=None,
                        maxjobs=None,
                        mempercore=None,
                        jobinterval=None,
                        localcores=None,
                        localmem=None,
                        disable_ui=False):
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
        the --maxjobs option (only if jobmode is not
        "local")
      jobinterval (int):  if specified, will be passed to
        the --jobinterval option
      localcores (int): if specified, will be passed to
        the --localcores option (only if jobmode is
        "local")
      localmem (int): if specified, will be passed to the
        the --localmem option (only if jobmode is
        "local")
      disable_ui (bool): if True, add the --disable-ui
        option (default is not to add it)

    Returns:
      Command: the original command updated with the
        appropriate options.
    """
    if jobmode is not None:
        cellranger_cmd.add_args("--jobmode=%s" % jobmode)
    if jobmode == "local":
        if localcores is not None:
            cellranger_cmd.add_args("--localcores=%s" %
                                    localcores)
        if localmem is not None:
            cellranger_cmd.add_args("--localmem=%s" %
                                    localmem)
    else:
        if mempercore is not None:
            cellranger_cmd.add_args("--mempercore=%s" %
                                    mempercore)
    if maxjobs is not None:
        cellranger_cmd.add_args("--maxjobs=%s" % maxjobs)
    if jobinterval is not None:
        cellranger_cmd.add_args("--jobinterval=%s" % jobinterval)
    if disable_ui:
        cellranger_cmd.add_args("--disable-ui")
    return cellranger_cmd


def make_multi_config_template(f, reference=None, fastq_dir=None,
                               samples=None, multiplexing=None, extensions=None,
                               no_bam=None, include_probe_set=None, probe_set=None,
                               cellranger_version=None):
    """
    Write a template configuration file for 'cellranger multi'

    Generates a template for the 'cellranger multi'
    configuration file. Specific options and sections are
    included or omitted depending on the arguments supplied
    to this function.

    The format and parameters for different data types are
    described in the 10x Genomics 'cellranger' documentation:

    * https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts

    Arguments:
      f (str): path that output template file will be
        written to
      reference (str): path to reference transcriptome
      fastq_dir (str): path to directory with Fastq files
      samples (list): list of sample names
      multiplexing (str): type of multiplexing (one of
        'cellplex', 'flex' or 'ocm', or None for singleplex
        data)
      extensions (list): list of "product extensions" (one
        or more of 'CSP', 'VDJ-T', 'VDJ-B') or None if there
        are no extensions
      no_bam (bool): if set then will be the value of the
        'no-bam' setting
      include_probe_set (bool): if set then the 'probe-set'
        setting will be included in the template (defaults
        to False; if 'probe_set' is defined then will be
        set to True automatically)
      probe_set (str): path to probe set CSV file
      cellranger_version (str): optionally specify the
        target CellRanger version number (or None to use
        the default version)
    """
    if not samples:
        samples = []
    if not extensions:
        extensions = []
    if probe_set:
        include_probe_set = True
    # Target CellRanger version
    if cellranger_version is None:
        cellranger_version = DEFAULT_CELLRANGER_VERSION
    else:
        cellranger_version = str(cellranger_version)
    # Major version
    try:
        cellranger_major_version = int(cellranger_version.split(".")[0])
    except Exception as ex:
        raise Exception(f"'{cellranger_version}': unable to extract "
                        f"Cellranger major version from version string: "
                        f"{ex}")
    if extensions is None:
        extensions = []
    # Write template file
    with open(f,'wt') as fp:
        # Header
        fp.write("## 10x_multi_config.csv\n"
                 "## See: https://www.10xgenomics.com/support/"
                 "software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts\n")
        # Gene expression section
        fp.write("[gene-expression]\n")
        fp.write("reference,%s\n" %
                 (reference if reference
                  else "/path/to/transcriptome"))
        if include_probe_set:
            fp.write("probe-set,%s\n" %
                     (probe_set if probe_set else "/path/to/probe_set"))
        fp.write("#force-cells,n\n")
        if cellranger_major_version == 7:
            # Cellranger 7.* targetted
            if no_bam is not None:
                fp.write("no-bam,%s\n" % str(no_bam).lower())
            else:
                fp.write("#no-bam,true|false\n")
        elif cellranger_major_version >= 8:
            # Cellranger >= 8.* targetted
            if no_bam is not None:
                fp.write("create-bam,%s\n" % str(not no_bam).lower())
            else:
                fp.write("create-bam,true\n")
        fp.write("#cmo-set,/path/to/custom/cmo/reference\n")
        # Feature section
        if "CSP" in extensions:
            fp.write(dedent("""
            #[feature]
            #reference,/path/to/feature_ref
            """))
        # V(D)J section
        if "VDJ-B" in extensions or "VDJ-T" in extensions:
            fp.write(dedent("""
            [vdj]
            #reference,/path/to/vdj_reference
            """))
        # Libraries section
        fp.write(dedent("""
        [libraries]
        fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
        """))
        tenx_library_types = ["Gene Expression",]
        if samples and multiplexing:
            if multiplexing in ("cellplex", "ocm", "flex"):
                tenx_library_types.append("Multiplexing Capture")
        if extensions:
            if "CSP" in extensions:
                tenx_library_types.append("Antibody Capture")
            if "VDJ-B" in extensions or "VDJ-T" in extensions:
                tenx_library_types.extend(["VDJ-B", "VDJ-T"])
        for sample in samples:
            library_type = "[%s]" % "|".join(tenx_library_types)
            fp.write("{sample},{fastqs_dir},any,{sample},{library_type},\n".format(
                sample=sample,
                fastqs_dir=(fastq_dir if fastq_dir else "/path/to/fastqs"),
                library_type=library_type))
        # Multiplexed samples section
        if multiplexing == "cellplex":
            fp.write(dedent("""
            [samples]
            sample_id,cmo_ids,description
            MULTIPLEXED_SAMPLE,CMO1|CMO2|...,DESCRIPTION
            """))
        elif multiplexing == "flex":
            fp.write(dedent("""
            [samples]
            sample_id,probe_barcode_ids,description
            MULTIPLEXED_SAMPLE,BC001|BC002|...,DESCRIPTION
            """))
        elif multiplexing == "ocm":
            fp.write(dedent("""[samples]
            sample_id,ocm_barcode_ids,description
            MULTIPLEXED_SAMPLE,OB1|OB2|...,DESCRIPTION
            """))
